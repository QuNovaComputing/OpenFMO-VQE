/*
 * Functions to handle molecular data
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <libgen.h>

#include "ofmo-parallel.h"
#include "ofmo-string.h"
#include "ofmo-basis.h"
#include "ofmo-basis-database.h"
#include "ofmo-mat.h"
#include "ofmo-integ.h"
#include "ofmo-cutoff.h"
#include "skel-rhf-main.h"
#include "skel-rhf-data.h"

#ifdef USE_CUDA
#include "cuda/cudalib.h"
#endif

FILE* fp_debug = NULL;

//#define _MOI_A0_  0.52917724924	// Bohr radius
#define _MOI_A0_  BOHR_RADIUS	// Bohr radius
#define TWO	2.0e0
#define HALF	0.5e0

static int verbose = false;

#define delete_pre_space(s) ofmo_delete_prepost_space(s)
#define delete_prepost_space(s) ofmo_delete_prepost_space(s)

static void dscale(int n, double alpha, double* src) {
    int i;
    for (i=0; i<n; i++) src[i] *= alpha;
}

static int readable_file( const char *filename)
{
  int ret;
  FILE *fp;
  if ((fp=fopen(filename, "r")) != NULL) fclose(fp);
  ret = (fp != NULL);
  if (verbose && !ret) fprintf(stderr, "NOT readable file: %s\n", filename);
  return ret;
}

int read_packed_matrix_binary(const int nao, const char *filename,
	double *ap) {
    int nao2, naod;
    FILE *fp;
    if ( (fp = fopen(filename, "rb")) == NULL) {
	printf("Can not open file (%s)!\n", filename);
	exit(1);
    }
    if (fread(&naod, sizeof(int), 1, fp) < 1) {
	printf("unexpected EOF\n");
	exit(1);
    }
    if (naod != nao) {
	printf( "number of basis function (NAO) "
		"are different between input files (%d vs. %d)\n",
		nao, naod);
	exit(1);
    }
    nao2 = nao * (nao + 1) / 2;
    if (fread(ap, sizeof(double), nao2, fp) != nao2) {
	printf("unexpected EOF\n");
	exit(1);
    }
    fclose(fp);
    return 0;
}

// Read just number of atom from input file
static int get_nat(const char* filename) {
    int nat;
    char ss[MAXSTRLEN];
    char *p;
    FILE *fp;

    if ( (fp=fopen(filename, "r")) == NULL) {
	printf("In get_nat: ERROR\n");
	printf("    File %s can\'t be opened\n", filename);
	return (-1);
    }
    // skip 1st line
    (void)fgets(ss, MAXSTRLEN, fp);

    // read 2nd line
    (void)fgets(ss, MAXSTRLEN, fp);
    (void)delete_pre_space(ss);
    p = strtok(ss, " 	");
    p = strtok(NULL, " 	");
    nat = atoi(p);
    fclose(fp);
    return(nat);
}

static int read_file( const char *filename, const int nat,
	char *basis_name, int *charge, int atomic_number[],
	double atom_x[], double atom_y[], double atom_z[] ) {
    char ss[MAXSTRLEN];
    char *p;
    int i;
    FILE *fp;

    if ( (fp=fopen(filename, "r")) == NULL) {
	printf("In read_file: ERROR\n");
	printf("    File %s can\'t be opened\n", filename);
	return -1;
    }
    // line 1: basis set
    (void)fgets(ss, MAXSTRLEN, fp);
    (void)delete_prepost_space(ss);
    (void)strcpy(basis_name, ss);

    // line 2: charge, number of atom
    (void)fgets(ss, MAXSTRLEN, fp);
    (void)delete_pre_space(ss);
    p = strtok(ss, " 	");
    (*charge) = atoi(p);

    // line 3-: elements, coordinates
    for (i=0; i<nat; i++) {
	if (fgets(ss, MAXSTRLEN, fp) == NULL) {
	    printf("unexpected EOF in read_file\n");
	    return -1;
	}
	(void)delete_pre_space(ss);
	p = strtok(ss, "    ");
	atomic_number[i] = AtomicNumber(p);
	if ( (p = strtok(NULL,"     ")) == NULL) {
	    printf("Illegal x-coordinate of %d-th atom\n", (i+1) );
	    return -1;
	}
	atom_x[i] = atof(p);
	if ( (p = strtok(NULL,"     ")) == NULL) {
	    printf("Illegal y-coordinate of %d-th atom\n", (i+1) );
	    return -1;
	}
	atom_y[i] = atof(p);
	if ( (p = strtok(NULL,"     ")) == NULL) {
	    printf("Illegal z-coordinate of %d-th atom\n", (i+1) );
	    return -1;
	}
	atom_z[i] = atof(p);
    }
    fclose(fp);
    return 0;
}

static void show_atomic_data( const int nat, const int atomic_number[],
	const double atom_x[], const double atom_y[],
	const double atom_z[] ) {
    int iat, n;
    n = (verbose || nat<10)? nat: 10;
    for ( iat=0; iat<n; iat++ ) {
        printf( "%2s  ", AtomicSymbol(atomic_number[iat]) );
	printf("%10.5f  %10.5f  %10.5f\n",
		atom_x[iat], atom_y[iat], atom_z[iat] );
    }
}

static int set_sorted_basis(
	const int nat, const int nsbs, char **basis_name, const int atomic_number[],
	const int atom_basis[],
	int *maxlqn, int *ncs, int *nao, int *nps,
	int **leading_cs, int **shel_tem, int **shel_atm, int **shel_add,
	int **shel_ini, double **prim_exp, double **prim_coe, int **sao2uao ) {
    int *ushel_lqn, *ushel_tem, *ushel_atm, *ushel_add, *ushel_ini;
    double *uprim_exp, *uprim_coe;
    int ierr = 0;
    ierr += ofmo_get_basis_size( nat, nsbs, basis_name, atomic_number, atom_basis,
	    maxlqn, ncs, nao, nps );
    ierr += ofmo_alloc_unsorted_basis( *ncs, *nps,
	    &ushel_lqn, &ushel_tem, &ushel_atm, &ushel_add, &ushel_ini,
	    &uprim_exp, &uprim_coe );
    ierr += ofmo_alloc_sorted_basis( *maxlqn, *ncs, *nao, *nps,
	    leading_cs, shel_tem, shel_atm, shel_add, shel_ini,
	    prim_exp, prim_coe, sao2uao );
    ierr += ofmo_assign_basis( nat, nsbs, basis_name, atomic_number, atom_basis,
	    ushel_lqn, ushel_tem, ushel_atm, ushel_add, ushel_ini,
	    uprim_exp, uprim_coe );
    ierr += ofmo_sort_basis( *maxlqn, *ncs,
	    ushel_lqn, ushel_tem, ushel_atm, ushel_add, ushel_ini,
	    uprim_exp, uprim_coe,
	    *leading_cs, *shel_tem, *shel_atm, *shel_add, *shel_ini,
	    *prim_exp, *prim_coe, *sao2uao );
    Free( ushel_lqn );
    Free( ushel_tem );
    Free( ushel_atm );
    Free( ushel_add );
    Free( ushel_ini );
    Free( uprim_exp );
    Free( uprim_coe );
    return ierr;
}

// ----------------

int skel_rhf_init_data( skel_rhf_config_t* config, skel_rhf_data_t* data )
{
    int ierr;
    int i;
    int *atom_basis;
    char basis_name[MAXSTRLEN];
    char *basis_name_list[1];
    double ia0;
    // basis set
    int *shel_tem, *shel_atm, *shel_add, *shel_ini;
    int *leading_cs, *sao2uao;
    double *prim_exp, *prim_coe;
    int nat, ncs, nao, nps, maxlqn;
    // atomic
    double *atom_x, *atom_y, *atom_z;
    int *atomic_number;
    int charge;
    int nelec, nocc;
    // CS-pair and PS-pair parameters
    int *leading_cs_pair, *csp_leading_ps_pair, *csp_ics, *csp_jcs;
    double *csp_schwarz, *psp_zeta, *psp_dkps, *psp_xiza;
    //
    int ics, ics0, ics1, jcs, ncspair, npspair, lqn;
    char CSTYPE[] = "SPDFGHIJ";
    //
    // MPI
    MPI_Comm comm = config->comm;
    int nprocs = config->nprocs;
    int myrank = config->myrank;

    verbose = config->verbose;

    if ( myrank == 0 ) {
      char *input=config->input;
      char *density=config->density;
        if (!readable_file(input) ||
            (density[0] != '\0' && !readable_file(density))) {
	    MPI_Abort( comm, 1 );
            return -1;
        }
	nat = get_nat( input );
	if (nat < 1) {
	    printf("Illegal number of atoms\n");
	    MPI_Abort( comm, 1 );
	}
	atomic_number = (int*)malloc(sizeof(int) * nat);
	atom_x = (double*)malloc(sizeof(double) * nat);
	atom_y = (double*)malloc(sizeof(double) * nat);
	atom_z = (double*)malloc(sizeof(double) * nat);
	ierr = read_file( input, nat, basis_name, &charge,
		atomic_number, atom_x, atom_y, atom_z );
	if ( ierr < 0 ) {
	    printf("Failure in reading input data\n");
	    MPI_Abort( comm, 1 );
	}
	printf("nat = %d\n", nat );
	printf("=============== INPUT DATA (A) ================\n");
	show_atomic_data( nat, atomic_number, atom_x, atom_y, atom_z );
	printf("===============================================\n");
        // translate Bohr to AU
	ia0 = 1.0 / _MOI_A0_;
	dscale( nat, ia0, atom_x );
	dscale( nat, ia0, atom_y );
	dscale( nat, ia0, atom_z );
	// debug
	fp_debug = stdout;

        // #electron
        for (i=0, nelec=0; i<nat; i++ ) nelec += atomic_number[i];
        nelec -= charge;
        nocc = nelec/2;
        printf( "nelec = %d,  nocc = %d\n", nelec, nocc );
    }
    MPI_Bcast( &nat, 1, MPI_INT, 0, comm );
    MPI_Bcast( &charge, 1, MPI_INT, 0, comm );
    MPI_Bcast( basis_name, MAXSTRLEN, MPI_CHAR, 0, comm );
    if ( myrank != 0 ) {
	atomic_number = (int*)malloc(sizeof(int) * nat);
	atom_x = (double*)malloc(sizeof(double) * nat );
	atom_y = (double*)malloc(sizeof(double) * nat );
	atom_z = (double*)malloc(sizeof(double) * nat );
    }
    MPI_Bcast( atomic_number, nat, MPI_INT, 0, comm );
    MPI_Bcast( atom_x, nat, MPI_DOUBLE, 0, comm );
    MPI_Bcast( atom_y, nat, MPI_DOUBLE, 0, comm );
    MPI_Bcast( atom_z, nat, MPI_DOUBLE, 0, comm );
    MPI_Bcast( &nelec, 1, MPI_INT, 0, comm );
    MPI_Bcast( &nocc, 1, MPI_INT, 0, comm );
    data->natom = nat;
    data->charge = charge;
    data->atomic_number = atomic_number;
    data->atom_x = atom_x;
    data->atom_y = atom_y;
    data->atom_z = atom_z;
    data->nelec = nelec;
    data->nocc = nocc;
    strncpy(data->basis_name, basis_name, MAXSTRLEN);

    // sort basis set
    basis_name_list[0] = basis_name;
    set_sorted_basis(
	    nat, 1, basis_name_list, atomic_number, NULL,
	    &maxlqn, &ncs, &nao, &nps, &leading_cs,
	    &shel_tem, &shel_atm, &shel_add, &shel_ini,
	    &prim_exp, &prim_coe, &sao2uao );
    data->maxlqn = maxlqn;
    data->ncs = ncs;
    data->nao = nao;
    data->nps = nps;
    data->leading_cs = leading_cs;
    data->shel_tem = shel_tem;
    data->shel_atm = shel_atm;
    data->shel_add = shel_add;
    data->shel_ini = shel_ini;
    data->prim_exp = prim_exp;
    data->prim_coe = prim_coe;
    data->sao2uao = sao2uao;
    if ( myrank == 0 ) {
	printf("nat=%d, ncs=%d, nao=%d, nps=%d\n",
		nat, ncs, nao, nps );
	fflush(stdout);
    }

    // Number of PS-pairs
    npspair = 0;
    for (ics=0; ics<ncs; ics++) {
	for (jcs=0; jcs<=ics; jcs++)
	    npspair += (shel_tem[ics]*shel_tem[jcs]);
    }
    data->npspair = npspair;

    // debug
    if (myrank == 0) {
	for (lqn=0; lqn<=maxlqn; lqn++) {
	    ics0 = leading_cs[lqn];
	    ics1 = leading_cs[lqn+1];
	    printf("%c-TYPE: [%3d, %3d)\n", CSTYPE[lqn], ics0, ics1);
            if (verbose) {
	    printf("    ICS: TEM, ATM, ADD, INI\n");
	    for (ics=ics0; ics<ics1; ics++) {
		printf("    %3d: %3d, %3d, %3d, %3d\n", ics,
			shel_tem[ics], shel_atm[ics],
			shel_add[ics], shel_ini[ics] );
	    }
	    printf("\n");
            }
	}
	fflush(stdout);
    }

    return 0;
}

// ----------------

void skel_rhf_fin_data( skel_rhf_data_t* data )
{
    data->natom = 0;
    Free( data->atomic_number);
    Free( data->atom_x );
    Free( data->atom_y );
    Free( data->atom_z );
    Free( data->leading_cs );
    Free( data->shel_tem );
    Free( data->shel_atm );
    Free( data->shel_add );
    Free( data->shel_ini );
    Free( data->sao2uao );
    Free( data->prim_exp );
    Free( data->prim_coe );

    Free( data->leading_cs_pair );
    Free( data->csp_leading_ps_pair );
    Free( data->csp_ics );
    Free( data->csp_jcs );
    Free( data->csp_schwarz );
    Free( data->psp_zeta );
    Free( data->psp_dkps );
    Free( data->psp_xiza );
}

// ----------------

int skel_rhf_cutoff_make_table(
    skel_rhf_config_t* config, skel_rhf_data_t* data )
{
  int rc = 0;
  int master = (config->myrank == 0);
  int ncs = data->ncs;
  int maxlqn = data->maxlqn;
  int npspair = data->npspair;
  int ncspair, maxlqn2;

  ncspair = ncs*(ncs+1)/2;
  maxlqn2 = (maxlqn+1) * (maxlqn+2) / 2;

  data->leading_cs_pair     = (int*)malloc( sizeof(int)*(maxlqn2+2) );
  data->csp_leading_ps_pair = (int*)malloc( sizeof(int)*(ncspair+1) );
  data->csp_ics             = (int*)malloc( sizeof(int)*ncspair );
  data->csp_jcs             = (int*)malloc( sizeof(int)*ncspair );
  data->csp_schwarz         = (double*)malloc( sizeof(double) * ncspair );
  data->psp_zeta = (double*)malloc( sizeof(double) * npspair );
  data->psp_dkps = (double*)malloc( sizeof(double) * npspair );
  data->psp_xiza = (double*)malloc( sizeof(double) * npspair );

  ofmo_integ_init( maxlqn );

  if ( master ) {
    printf("ncspair = %d, npspair = %d \n", ncspair, data->npspair );
    fflush(stdout);
  }
  ofmo_cutoff_make_table( maxlqn, data->leading_cs,
      data->shel_tem, data->shel_atm, data->shel_add,
      data->atom_x, data->atom_y, data->atom_z,
      data->prim_exp, data->prim_coe,
      data->leading_cs_pair, data->csp_schwarz,
      data->csp_ics, data->csp_jcs,
      data->csp_leading_ps_pair, data->psp_zeta,
      data->psp_dkps, data->psp_xiza );

  if ( master ) {
    ofmo_cutoff_show_table( maxlqn, data->leading_cs,
        data->shel_tem, data->leading_cs_pair,
        data->csp_leading_ps_pair);
  }

  return 0;
}

