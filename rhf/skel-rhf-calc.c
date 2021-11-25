/*
 * Interface routines to ofmo_*() function
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>

#include "ofmo-parallel.h"
#include "ofmo-prof.h"
#include "ofmo-scf.h"
#include "ofmo-integ.h"
#include "ofmo-twoint.h"
#include "ofmo-mat.h"
#include "skel-rhf-main.h"
#include "skel-rhf-data.h"
#include "skel-rhf-calc.h"

// ----------------

int skel_rhf_init_array(
        skel_rhf_config_t* config, skel_rhf_data_t* data,
        skel_rhf_array_t* array )
{
  int nao = data->nao;
  int nao2 = nao*(nao+1)/2;
  int nat = data->natom;
  double *D, *G, *TEMP, *S, *H, *F, *C, *moe;
  double *aopop, *atpop;
  D = (double*)malloc(sizeof(double) * nao2);
  if (D==NULL) return -1;
  S = (double*)malloc(sizeof(double) * nao2);
  if (S==NULL) return -1;
  H = (double*)malloc(sizeof(double) * nao2);
  if (H==NULL) return -1;
  C = (double*)malloc(sizeof(double) * nao*nao);
  if (C==NULL) return -1;
  moe = (double*)malloc(sizeof(double) * nao);
  if (moe==NULL) return -1;
  aopop = (double*)malloc(sizeof(double) * nao);
  if (aopop==NULL) return -1;
  atpop = (double*)malloc(sizeof(double) * nat);
  if (atpop==NULL) return -1;

  array->D = D;
  array->S = S;
  array->H = H;
  array->C = C;
  array->moe = moe;
  array->aopop = aopop;
  array->atpop = atpop;

  return 0;
}

// ----------------

void skel_rhf_fin_array( skel_rhf_array_t* array )
{
  Free(array->D);
  Free(array->S);
  Free(array->H);
  Free(array->C);
  Free(array->moe);
  Free(array->aopop);
  Free(array->atpop);
}

// ----------------

int skel_rhf_oneint(
        skel_rhf_config_t* config, skel_rhf_data_t* data,
        skel_rhf_array_t* array )
{
  int nao = data->nao;
  int nao2 = nao*(nao+1)/2;
  double *S = array->S;
  double *H = array->H;

  memset( S, '\0', sizeof(double)*nao2 );
  memset( H, '\0', sizeof(double)*nao2 );

#pragma omp parallel
  {
    int nworkers, workerid;
    nworkers = omp_get_num_threads();
    workerid = omp_get_thread_num();
    ofmo_integ_oneint_sorted(
        nworkers, workerid,
        data->maxlqn, data->leading_cs,
        data->shel_tem, data->shel_atm, data->shel_add, data->shel_ini,
        data->atom_x, data->atom_y, data->atom_z,
        data->prim_exp, data->prim_coe,
        data->natom, data->atomic_number, S, H );
  }
  ofmo_scale_diag( nao, 0.5e0, H );

  return 0;
}

// ----------------

int skel_rhf_init_density(
        skel_rhf_config_t* config, skel_rhf_data_t* data,
        skel_rhf_array_t* array )
{
  int nao = data->nao;
  int nao2 = nao*(nao+1)/2;
  double *D = array->D;

  if (strlen(config->density) > 0) {
    if (config->myrank == 0) {
      read_packed_matrix_binary( nao, config->density, D );
    }
    MPI_Bcast( D, nao2, MPI_DOUBLE, 0, config->comm );
  } else {
    ofmo_scf_init_density_ehuckel(
        data->natom, data->ncs, data->nao, data->maxlqn, data->nocc,
        data->atomic_number, data->leading_cs,
        data->shel_atm, data->shel_ini,
        array->S, array->D, array->aopop, array->atpop );
  }
  return 0;
}

// ----------------

int skel_rhf_twoint_first(
        skel_rhf_config_t* config, skel_rhf_data_t* data,
        skel_rhf_array_t* array )
{
#pragma omp parallel
  {
    int incore2 = -1;
    int workerid, nworkers;
    int nthreads, mythread;
    size_t my_buffer_size = 0;
    // for profiling
    double et0, et1;
    nthreads = omp_get_num_threads();
    mythread = omp_get_thread_num();
    workerid = config->myrank * nthreads + mythread;
    nworkers = config->nprocs * nthreads;
    my_buffer_size = config->eribfsz/nthreads;
    int offset = 0;
    ofmo_integ_set_loop_offset( mythread, offset );

    incore2 = ofmo_integ_twoint_first( nworkers, workerid, my_buffer_size,
        data->maxlqn, data->shel_atm, data->shel_ini,
        data->atom_x, data->atom_y, data->atom_z,
        data->leading_cs_pair,
        data->csp_schwarz, data->csp_ics, data->csp_jcs,
        data->csp_leading_ps_pair,
        data->psp_zeta, data->psp_dkps, data->psp_xiza );
  }

  return 0;
}

// ----------------

double skel_rhf_scf(
        skel_rhf_config_t* config, skel_rhf_data_t* data,
        skel_rhf_array_t* array )
{
  double energy, Enuc;
#pragma omp parallel
  {
    int mythread = omp_get_thread_num();
    int offset = 0;
    ofmo_integ_set_loop_offset( mythread, offset );
  }

  // nuclear repulsion
  Enuc = ofmo_calc_nuclear_repulsion( data->natom, data->atomic_number,
      data->atom_x, data->atom_y, data->atom_z);
  if ( config->myrank == 0 ) {
    printf( "Enuc = %17.10f\n", Enuc );
    printf( "scfd scfe eps_ps4 eps_sch eps_eri: %8.2e %8.2e %8.2e %8.2e %8.2e\n",
        config->scfd, config->scfe, config->eps_ps4, config->eps_sch, config->eps_eri);
    fflush(stdout);
  }

  // scf
  ofmo_twoint_eps_ps4(config->eps_ps4);
  ofmo_twoint_eps_eri(config->eps_eri);
  ofmo_twoint_eps_sch(config->eps_sch);
  
  // If you want MO ERI
  const int nao = data->nao;
  const int nao_4 = nao * nao * nao * nao;
  double * mo_tei = (double * )malloc(sizeof(double) * nao_4);

  if (!config->dryrun) {
    ofmo_scf_rhf( config->comm,
        data->maxlqn, Enuc, data->ncs, data->nao,
        data->leading_cs,
        data->shel_atm, data->shel_ini,
        data->atom_x, data->atom_y, data->atom_z,
        data->leading_cs_pair,
        data->csp_schwarz, data->csp_ics, data->csp_jcs,
        data->csp_leading_ps_pair,
        data->psp_zeta, data->psp_dkps, data->psp_xiza,
        data->natom, data->nocc,
        array->S, array->H,
        config->maxscfcyc, config->scfe, config->scfd,
        array->D, array->C, mo_tei, array->moe, &energy );
  }
  const int nao_2 = nao * nao;
  const int nao_3 = nao_2 * nao;
  int imo, imo4, jmo, jmo3, kmo, kmo2, lmo, mo_idx;
  int count_mo = 0;
  printf("===MO ERI output===\n");
  for(imo=0, imo4=0; imo<nao;   imo++, imo4+=nao_3){
  for(jmo=0, jmo3=0; jmo<imo+1; jmo++, jmo3+=nao_2){
  for(kmo=0, kmo2=0; kmo<imo+1; kmo++, kmo2+=nao){
  for(lmo=0; lmo<kmo+1; lmo++, count_mo++){
      /* mo_idx = imo * nao * nao * nao
              +jmo * nao * nao
              +kmo * nao
              +lmo; */
      mo_idx = imo4 + jmo3 + kmo2 + lmo;
      printf("%2d %2d %2d %2d | %.7f\n", imo, jmo, kmo, lmo, mo_tei[mo_idx]);
  }}}}
  fflush(stdout);
  printf("===END MO ERI output===\n");
  free(mo_tei);
  return energy;
}

// ----------------
