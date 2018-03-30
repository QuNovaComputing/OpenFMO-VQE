#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <strings.h>
#include <assert.h>

#include <cuda.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "cudalib.h"
#include "cuda-fmt.h"
#include "cuda-fmt-m.h"
#include "cuda-integ.h"
#include "cuda-twoint-direct.h"
#include "cuda-ifc4c-calc.h"
#include "check.h"

#ifdef __cplusplus
extern "C" {
#endif
//#include "integ/ofmo-twoint-direct.h"
#include "integ/ofmo-index.h"
extern FILE* fp_prof; // from common/ofmo-prof.h
#ifdef __cplusplus
}
#endif

static int CUDA_OPTSYNC = FALSE;
//static int CUDA_OPTSYNC = TRUE;

static int SCF_Malloced = FALSE;

static int CUDA_NDEV = 0;
static struct dev_Data *devData = NULL;

static int nL[] = {1, 3, 6, 10};

static int CUDA_ME = 0;
static int CUDA_NP = 1;

static int CUDA_NBLOCKS  = 0;
static int CUDA_NTHREADS = 0;

/* ------------------------------------- */

__device__ __constant__ int ncs;
__device__ __constant__ int nat;
__device__ __constant__ int maxlqn;
__device__ __constant__ int ncspair;
__device__ __constant__ int npspair;
__device__ __constant__ int nao;
__device__ __constant__ int max_num_klcs;
__device__ __constant__ int num_Gkl;
__device__ __constant__ int num_Gi;
// DFT
__device__ __constant__ double x_coef;

//__device__ __constant__ int *p_ijcounter;
__device__ int *ijcounter;

__device__ __constant__ int *shel_atm;
__device__ __constant__ int *shel_ini;
__device__ __constant__ double *atom_x;
__device__ __constant__ double *atom_y;
__device__ __constant__ double *atom_z;
__device__ __constant__ int *leading_cs;
__device__ __constant__ int *leading_cs_pair;
__device__ __constant__ int *csp_leading_ps_pair;
__device__ __constant__ int *csp_ics;
__device__ __constant__ int *csp_jcs;
__device__ __constant__ int *sklcs_b;
__device__ __constant__ double *psp_zeta;
__device__ __constant__ double *psp_dkps;
__device__ __constant__ double *psp_xiza;
__device__ __constant__ double *DFACT;
__device__ __constant__ float *csp_schwarz;
__device__ __constant__ float *Dcs;
__device__ __constant__ int *sorted_csp;
__device__ __constant__ double *Ds;
__device__ __constant__ double *G_d;
__device__ __constant__ double *G_b;
__device__ __constant__ double *Gkl_b;
__device__ __constant__ double *Gi_t;
//__device__ __constant__ double *Gj_t;
__device__ __constant__ double *work_b;

/* ------------------------------------- */
__host__ static void cuda_init_dev_Data(struct dev_Data *dev)
{
  dev->ijcounter = NULL;
  for (int i=0; i<10; i++) dev->fmt_table[i] = NULL;
  for (int i=0; i<9; i++) dev->fmt_m_table[i] = NULL;
  dev->shel_atm = NULL;
  dev->shel_ini = NULL;
  dev->atom_x = NULL;
  dev->atom_y = NULL;
  dev->atom_z = NULL;
  dev->leading_cs = NULL;
  dev->leading_cs_pair = NULL;
  dev->csp_leading_ps_pair = NULL;
  dev->csp_ics = NULL;
  dev->csp_jcs = NULL;
  dev->sklcs_b = NULL;
  dev->psp_zeta = NULL;
  dev->psp_dkps = NULL;
  dev->psp_xiza = NULL;
  dev->DFACT = NULL;
  dev->csp_schwarz = NULL;
  dev->Dcs = NULL;
  dev->sorted_csp = NULL;
  dev->Ds = NULL;
  dev->G_d = NULL;
  dev->G_b = NULL;
  dev->Gkl_b = NULL;
  dev->Gi_t = NULL;
  dev->work_b = NULL;
// for ifc4c
  dev->shel_atm_frg = NULL;
  dev->shel_ini_frg = NULL;
  dev->atom_x_frg = NULL;
  dev->atom_y_frg = NULL;
  dev->atom_z_frg = NULL;
  dev->leading_cs_frg = NULL;
  dev->leading_cs_pair_frg = NULL;
  dev->csp_leading_ps_pair_frg = NULL;
  dev->csp_ics_frg = NULL;
  dev->csp_jcs_frg = NULL;
  dev->psp_zeta_frg = NULL;
  dev->psp_dkps_frg = NULL;
  dev->psp_xiza_frg = NULL;
  dev->csp_schwarz_frg = NULL;
  dev->shel_atm_mon = NULL;
  dev->shel_ini_mon = NULL;
  dev->atom_x_mon = NULL;
  dev->atom_y_mon = NULL;
  dev->atom_z_mon = NULL;
  dev->leading_cs_mon = NULL;
  dev->leading_cs_pair_mon = NULL;
  dev->csp_leading_ps_pair_mon = NULL;
  dev->csp_ics_mon = NULL;
  dev->csp_jcs_mon = NULL;
  dev->psp_zeta_mon = NULL;
  dev->psp_dkps_mon = NULL;
  dev->psp_xiza_mon = NULL;
  dev->csp_schwarz_mon = NULL;
  dev->D_mon = NULL;
  dev->V_frg = NULL;
  dev->V_frgP = NULL;
}



/* ------------------------------------- */
__host__ struct dev_Data *cuda_SCF_get_dev_Data(const int ndev)
{
  int NDEV = cuda_get_numDevice();
  if (ndev<0 || ndev>=NDEV) return NULL;
  return devData+ndev;
}

__host__ int cuda_set_Device(int idev)
{
  int NDEV = cuda_get_numDevice();
  if (NDEV > 1) {
    cudaError_t err;
    err = cudaSetDevice(idev);
    if (err != cudaSuccess) return -1;
  }
  return 0;
}

__host__ int cuda_get_numDevice(void)
{
#if 0
  int ndev = 0;
#ifdef _OPENMP
#pragma omp master
  {
    ndev = CUDA_NDEV;
  }
#else
    ndev = CUDA_NDEV;
#endif
#else
  int ndev = CUDA_NDEV;
#ifdef _OPENMP
  int master = FALSE;
#pragma omp master
  master = TRUE;
  if (CUDA_NDEV>0 && !master) ndev = -1;
#endif
#endif
  return ndev;
}

__host__ int cuda_use_Device(void)
{
  return (CUDA_NDEV>0);
}

/* ------------------------------------- */

__host__ int cuda_get_numBlocks(void)
{
  return CUDA_NBLOCKS;
}

__host__ int cuda_get_numThreads(void)
{
  return CUDA_NTHREADS;
}

__host__ int cuda_set_BT(const int nblk0, const int nthb0, const int maxlqn)
{
  int i;
  int nblk=0;
  int nthb=0;
  int La, Lb, Lc, Ld;
  int Lab, Lcd, Labcd;

  for ( La=0; La<=maxlqn; La++ ) {
    for ( Lb=0; Lb<=La; Lb++ ) {
      Lab = La*(La+1)/2 + Lb;
      for ( Lc=0; Lc<=La; Lc++ ) {
        for ( Ld=0; Ld<=(Lc==La? Lb : Lc ); Ld++ ) {
          Lcd = Lc*(Lc+1)/2 + Ld;
          Labcd = Lab*(Lab+1)/2 + Lcd;
          //if (nblk0>0) dim2e[Labcd][0] = nblk0;
          if (nblk0>0 && dim2e[Labcd][0]!=0) dim2e[Labcd][0] = nblk0;
          if (nthb0>0) dim2e[Labcd][1] = nthb0;
          nblk = MAX2(nblk, dim2e[Labcd][0]);
          nthb = MAX2(nthb, dim2e[Labcd][1]);
          if (dim2e[Labcd][1]%WARP_SIZE != 0) return -1;
        }
      }
    }
  }
  CUDA_NBLOCKS = nblk;
  CUDA_NTHREADS = nthb;
  return nblk*nthb;
}

/* ------------------------------------- */

__host__ int cuda_get_Nprocs(void)
{
  return CUDA_NP;
}

__host__ int cuda_get_myRank(void)
{
  return CUDA_ME;
}
/* ------------------------------------- */
__host__ int cuda_find2eType(char *str)
{
  int i;
  int n = cuda_get_num_types();
  if (str==NULL) return -1;
  for (i=0; i<n; i++)
    if (strncasecmp(str, cuda_s2e[i], 4)==0) return i;
  return -1;
}

__host__ int cuda_set_optCPU(char *optstr)
{
  int n = cuda_get_num_types();
  char *str, *str1, *str2;
  char *savptr1, *savptr2;
  char *stype1, *stype2;
  char *idx;
  int type1, type2;
  int i;

  //printf("optCPU = %s\n", optstr);
  for (str1=optstr;; str1=NULL) {
    str = strtok_r(str1, ", ", &savptr1);
    if (str==NULL) break;
    idx = index(str, '-');
    stype1 = strtok_r(str, "-", &savptr2);
    stype2 = strtok_r(NULL, "- ", &savptr2);
    type1 = cuda_find2eType(stype1);
    type2 = cuda_find2eType(stype2);
    //printf("%d %d\n", type1, type2);
    if (type1<0) return -1;
    if (idx == NULL) dim2e[type1][0] = 0;
    else if (stype2 != NULL)
      for (i=type1;i<=type2;i++) dim2e[i][0] = 0;
    else if (stype1 == str)
      for (i=type1;i<=n;i++) dim2e[i][0] = 0;
    else
      for (i=0;i<=type1;i++) dim2e[i][0] = 0;
  }
  //for (i=0;i<n;i++) printf("(%2d) %s %1d\n",i,s2e[i],dim2e[i][0]);
  return 0;
}

__host__ int cuda_get_optCPU(int Labcd)
{
  return (dim2e[Labcd][0]>0);
}

/* ------------------------------------- */
__host__ int cuda_set_optsync(int optsync)
{
  CUDA_OPTSYNC = optsync;
  return CUDA_OPTSYNC;
}

__host__ int cuda_get_optsync(void)
{
  return CUDA_OPTSYNC;
}

/* ------------------------------------- */
__host__ void cuda_Print_DEFS(void)
{
  printf("DEFS: "
#ifdef SORT_CSP
  "SORT_CSP "
#endif
#ifdef SORT_CSP_SCHWARZ
  "SORT_CSP_SCHWARZ "
#endif
#ifdef SORT_CSP_REV
  "SORT_CSP_REV "
#endif
  );
  if (cuda_get_optsync()) printf("SYNC_TWOINT ");
  printf("\n");

  printf("CUDA_DEFS: "
#ifdef GPU_DLB
  "GPU_DLB "
#endif
#ifdef USE_ATOMIC
  "USE_ATOMIC "
#endif
#ifdef USE_INSTANT_SCHWARZ
  "USE_INSTANT_SCHWARZ "
#endif
#ifdef ADD_FULL_NAO
  "ADD_FULL_NAO "
#endif
#ifdef DLB_KL
  "DLB_KL "
#endif
/*
#ifdef DLB_KL_SSSS
  "DLB_KL_SSSS "
#endif
#ifdef DLB_KL_PSSS
  "DLB_KL_PSSS "
#endif
#ifdef DLB_KL_PSPS
  "DLB_KL_PSPS "
#endif
#ifdef DLB_KL_PPSS
  "DLB_KL_PPSS "
#endif
#ifdef DLB_KL_PPPS
  "DLB_KL_PPPS "
#endif
#ifdef DLB_KL_PPPP
  "DLB_KL_PPPP "
#endif
#ifdef DLB_KL_DSSS
  "DLB_KL_DSSS "
#endif
#ifdef DLB_KL_DSPS
  "DLB_KL_DSPS "
#endif
#ifdef DLB_KL_DSPP
  "DLB_KL_DSPP "
#endif
#ifdef DLB_KL_DSDS
  "DLB_KL_DSDS "
#endif
#ifdef DLB_KL_DPSS
  "DLB_KL_DPSS "
#endif
#ifdef DLB_KL_DPPS
  "DLB_KL_DPPS "
#endif
#ifdef DLB_KL_DPPP
  "DLB_KL_DPPP "
#endif
#ifdef DLB_KL_DPDS
  "DLB_KL_DPDS "
#endif
#ifdef DLB_KL_DPDP
  "DLB_KL_DPDP "
#endif
#ifdef DLB_KL_DDSS
  "DLB_KL_DDSS "
#endif
#ifdef DLB_KL_DDPS
  "DLB_KL_DDPS "
#endif
#ifdef DLB_KL_DDPP
  "DLB_KL_DDPP "
#endif
#ifdef DLB_KL_DDDS
  "DLB_KL_DDDS "
#endif
#ifdef DLB_KL_DDDP
  "DLB_KL_DDDP "
#endif
#ifdef DLB_KL_DDDD
  "DLB_KL_DDDD "
#endif
*/
#ifdef SORT_IJ_SCHWARZ
  "SORT_IJ_SCHWARZ "
#endif
#ifdef SORT_KL_SCHWARZ
  "SORT_KL_SCHWARZ "
#endif
  );
#ifdef CUDA_FMT_M
#ifdef CUDA_FMT_M_SM
  printf("CUDA_FMT_M(%d,%d,shared)",CUDA_FMT_M, CUDA_FMT_M_NEXP);
#else
  printf("CUDA_FMT_M(%d,%d)",CUDA_FMT_M, CUDA_FMT_M_NEXP);
#endif
#endif
  printf("\n");
}

/* ------------------------------------- */

//__host__ int cuda_Init(const int ndev, const int myrank, const int nprocs)
__host__ int cuda_Init_Sub(const int ndev, const int myrank, const int nprocs)
{
  int NDEV = 0;
#ifdef _OPENMP
  int master = FALSE;
#pragma omp master
  {
    master = TRUE;
  }
  if (!master) return 0;
#endif
  CUDA_ME = myrank;
  CUDA_NP = nprocs;

  if (ndev < 0) return 1;
  if (ndev == 0) {
    CUDA_NDEV = 0;
    return 0;
  };
  cudaGetDeviceCount(&NDEV);
  if (NDEV>1) NDEV = 1; // Force single dev.
  //if (ndev > NDEV) return -1;
  if (ndev > NDEV) {
    char hn[16];
    gethostname(hn, 16);
    printf("%s: ndev %d > NDEV %d\n",hn, myrank, ndev, NDEV);
  }
  assert(ndev <= NDEV);

  devData=(struct dev_Data *)malloc(ndev*sizeof(struct dev_Data));
  if (devData==NULL) return -1;
  for (int i=0; i<ndev; i++) {
    struct dev_Data *dev=devData+i;
    cuda_init_dev_Data(dev);
  }

  int i;
  cudaError_t err;
  struct cudaDeviceProp prop;
  struct dev_Data *dev;
  if (ndev==1) {
//    i=myrank/4;
//    cudaSetDevice(i);
    cudaDeviceReset();
    cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
    cudaGetDevice(&i);
    err=cudaGetDeviceProperties(&prop, i);
    if (err!=cudaSuccess) return -2;
    dev = devData;
    dev->numSM = prop.multiProcessorCount;
  } else {
  for (i=0; i<ndev; i++) {
    cuda_set_Device(i);
    cudaDeviceReset();
    cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
    err=cudaGetDeviceProperties(&prop, i);
    if (err!=cudaSuccess) return -2;
    dev = devData + i;
    dev->numSM = prop.multiProcessorCount;
  }
  }
  for (i=0; i<21; i++) {
    dim2e[i][0] *= dev->numSM;
    dim2e[i][1] *= WARP_SIZE;
  }
  {
    int nblk=-1, nthd=-1;
    char *p;
    if ((p=getenv("OFMO_CUDA_IFC4C_BLK"))!=NULL) nblk = atoi(p);
    if ((p=getenv("OFMO_CUDA_IFC4C_THD"))!=NULL) nthd = atoi(p);
    for (i=0; i<12*3; i++) {
      if (dim_ifc4c[i][0]!=0) {
        if (nblk>=0) dim_ifc4c[i][0] = nblk;
        if (nthd>0) dim_ifc4c[i][1] = nthd;
        dim_ifc4c[i][0] *= dev->numSM;
        dim_ifc4c[i][1] *= WARP_SIZE;
      }
    }
  }
  cudaGetDevice(&i);
  if(fp_prof) fprintf(fp_prof, "(%2d) cuda_Init(): %d %d\n",myrank,NDEV,i);
  if(fp_prof) fprintf(fp_prof, "(%2d) cuda_Init(): %d %d [%s (%d)] %d\n",
        myrank, NDEV, i, prop.name , CUDA_ARCH, dev->numSM);
  CUDA_NDEV = ndev;

  return 0;
}

__host__ int cuda_Finalize_Sub(void)
{
  int i;
  int NDEV = cuda_get_numDevice();

  if (NDEV<=0) return 2;
  if (devData==NULL) return -1;
  //if (CUDA_ME==0) cuda_print_wifc4c();
  cuda_SCF_Free();
  cuda_FMT_Finalize();
#ifdef CUDA_FMT_M
  cuda_FMT_m_Finalize();
#endif
  free(devData);
  devData = NULL;
  for (i=0; i<NDEV; i++) {
    cuda_set_Device(i);
    //cudaDeviceReset();
  }
  CUDA_NDEV = 0;

  return 0;
}

#define McudaDevMalloc(name, num, kind) \
  { ret--;\
    if (dev->name!=NULL) {ret-=1000;break;}; \
    if ((err=cudaMalloc((void **)&(dev->name), (num)*sizeof(kind)))!=cudaSuccess) {break;}}

/* ------------------------------------- */
__host__ int cuda_SCF_Malloc(const int ncs, const int nat, const int maxlqn,
    const int ncspair, const int npspair, const int nao,
    const int max_num_klcs, const int num_Gkl, const int num_Gi)
{
  char fn[]="cuda_SCF_Malloc";
  cudaError_t err;
  int i, ret;
  int NDEV = cuda_get_numDevice();
  int maxlqn2 = (maxlqn+1)*(maxlqn+2)/2;
  int nblk = cuda_get_numBlocks();
  int nthb = cuda_get_numThreads();

  if (NDEV<=0) return 2;
  if (devData==NULL) return -1;
  if (SCF_Malloced) return 1;
  if(fp_prof) fprintf(fp_prof, "(%d) %s()\n", CUDA_ME, fn);
  for (i=0; i<NDEV; i++) {
    struct dev_Data *dev = devData + i;
    cuda_set_Device(i);
    ret = -(i+1)*100;
    McudaDevMalloc(Gi_t, num_Gi*nblk*nthb, double);
    McudaDevMalloc(ijcounter, 21, int);
    McudaDevMalloc(shel_atm, ncs, int);
    McudaDevMalloc(shel_ini, ncs+1, int);
    McudaDevMalloc(atom_x, nat, double);
    McudaDevMalloc(atom_y, nat, double);
    McudaDevMalloc(atom_z, nat, double);
    McudaDevMalloc(leading_cs, maxlqn+2, int);
    McudaDevMalloc(leading_cs_pair, maxlqn2+2, int);
    McudaDevMalloc(csp_leading_ps_pair, ncspair+1, int);
    McudaDevMalloc(csp_ics, ncspair, int);
    McudaDevMalloc(csp_jcs, ncspair, int);
    McudaDevMalloc(sklcs_b, max_num_klcs*nblk, int);
    McudaDevMalloc(psp_zeta, npspair, double);
    McudaDevMalloc(psp_dkps, npspair, double);
    McudaDevMalloc(psp_xiza, npspair, double);
    McudaDevMalloc(DFACT, 36, double);
    McudaDevMalloc(csp_schwarz, ncspair, float);
    McudaDevMalloc(Dcs, ncs*ncs, float);
    McudaDevMalloc(sorted_csp, ncspair, int);
    McudaDevMalloc(Ds, nao*nao, double);
    McudaDevMalloc(G_d, nao*nao, double);
    McudaDevMalloc(G_b, nao*nao*nblk, double);
    McudaDevMalloc(Gkl_b, num_Gkl*nblk, double);
    McudaDevMalloc(work_b, WORK_SIZE*nblk*nthb, double);
    ret = 0;
  }
  if (ret<0 && fp_prof) {
    if (ret<-1000) {
      fprintf(fp_prof, "(%d) %s() %d: may not clean pointer!\n", CUDA_ME, fn, ret);
    } else {
      fprintf(fp_prof, "(%d) %s() %d: %s\n", CUDA_ME, fn, ret, cudaGetErrorString(err));
    }
    fprintf(fp_prof, "(%d) %s(%d,%d,%d,%d,%d,%d,\n"
                     "        %d,%d,%d)\n", CUDA_ME, fn,
                     ncs, nat, maxlqn, ncspair, npspair, nao,
                     max_num_klcs, num_Gkl, num_Gi);
    fflush(fp_prof);
    exit(ret);
  }

  SCF_Malloced = TRUE;
  //if(fp_prof) fprintf(fp_prof, "(%d) cuda_SCF_Malloc() done\n", CUDA_ME);
  //if(fp_prof) fflush(fp_prof); 

  return 0;
}
#undef McudaDevMalloc

#define McudaDevFree(name) \
    { ret--; \
      if ((err=cudaFree(dev->name))!=cudaSuccess) {break;};\
      dev->name = NULL;}

__host__ int cuda_SCF_Free(void)
{
  char fn[]="cuda_SCF_Free";
  cudaError_t err;
  int i,ret;
  int NDEV = cuda_get_numDevice();

  if (NDEV<=0) return 2;
  if (devData==NULL) return -1;
//  cuda_FMT_Finalize();
  if (!SCF_Malloced) return 1;
  if(fp_prof) fprintf(fp_prof, "(%d) %s()\n", CUDA_ME, fn);

  for (i=0; i<NDEV; i++) {
    ret = -(i+1)*100;
    struct dev_Data *dev = devData + i;
    cuda_set_Device(i);
    McudaDevFree(Gi_t);
    McudaDevFree(ijcounter);
    McudaDevFree(shel_atm);
    McudaDevFree(shel_ini);
    McudaDevFree(atom_x);
    McudaDevFree(atom_y);
    McudaDevFree(atom_z);
    McudaDevFree(leading_cs);
    McudaDevFree(leading_cs_pair);
    McudaDevFree(csp_leading_ps_pair);
    McudaDevFree(csp_ics);
    McudaDevFree(csp_jcs);
    McudaDevFree(sklcs_b);
    McudaDevFree(psp_zeta);
    McudaDevFree(psp_dkps);
    McudaDevFree(psp_xiza);
    McudaDevFree(DFACT);
    McudaDevFree(csp_schwarz);
    McudaDevFree(Dcs);
    McudaDevFree(sorted_csp);
    McudaDevFree(Ds);
    McudaDevFree(G_d);
    McudaDevFree(G_b);
    McudaDevFree(Gkl_b);
    McudaDevFree(work_b);
    ret = 0;
  }
  if (ret<0 && fp_prof) {
    fprintf(fp_prof, "(%d) %s() %d: %s\n", CUDA_ME, fn, ret, cudaGetErrorString(err));
    fflush(fp_prof);
    exit(ret);
  }
  SCF_Malloced = FALSE;

  return 0;
}
#undef McudaDevFree

/* ------------------------------------- */
#define McudaSetSymbol(a) checkCudaErrors(cudaMemcpyToSymbol(a, &(dev->a), sizeof(dev->a)))
#define McudaMemcpyToSymbol2(a,b,c) checkCudaErrors(cudaMemcpyToSymbol(a, b, c))

__host__ int cuda_SCF_Init_Symbol(
    const int ncs0, const int nat0, const int maxlqn0,
    const int ncspair0, const int npspair0, const int nao0,
    const int max_num_klcs0, const int num_Gkl0, const int num_Gi0,
    const struct dev_Data *dev)
{
    cudaError_t err;

    if (dev==NULL) return -1;
    McudaMemcpyToSymbol2(ncs, &ncs0, sizeof(int));
    McudaMemcpyToSymbol2(nat, &nat0, sizeof(int));
    McudaMemcpyToSymbol2(maxlqn, &maxlqn0, sizeof(int));
    McudaMemcpyToSymbol2(ncspair, &ncspair0, sizeof(int));
    McudaMemcpyToSymbol2(npspair, &npspair0, sizeof(int));
    McudaMemcpyToSymbol2(nao, &nao0, sizeof(int));
    McudaMemcpyToSymbol2(max_num_klcs, &max_num_klcs0, sizeof(int));
    McudaMemcpyToSymbol2(num_Gkl, &num_Gkl0, sizeof(int));
    McudaMemcpyToSymbol2(num_Gi, &num_Gi0, sizeof(int));
//    McudaSetSymbol(p_ijcounter);
    McudaSetSymbol(ijcounter);
    McudaSetSymbol(shel_atm);
    McudaSetSymbol(shel_ini);
    McudaSetSymbol(atom_x);
    McudaSetSymbol(atom_y);
    McudaSetSymbol(atom_z);
    McudaSetSymbol(leading_cs);
    McudaSetSymbol(leading_cs_pair);
    McudaSetSymbol(csp_leading_ps_pair);
    McudaSetSymbol(csp_ics);
    McudaSetSymbol(csp_jcs);
    McudaSetSymbol(sklcs_b);
    McudaSetSymbol(psp_zeta);
    McudaSetSymbol(psp_dkps);
    McudaSetSymbol(psp_xiza);
    McudaSetSymbol(DFACT);
    McudaSetSymbol(csp_schwarz);
    McudaSetSymbol(Dcs);
    McudaSetSymbol(sorted_csp);
    McudaSetSymbol(Ds);
    McudaSetSymbol(G_d);
    McudaSetSymbol(G_b);
    McudaSetSymbol(Gkl_b);
    McudaSetSymbol(Gi_t);
    McudaSetSymbol(work_b);

    return 0;
}
#undef McudaSetSymbol
#undef McudaSetSymbol2

__host__ int cuda_SCF_Finalize(void)
{
  int i;
  int NDEV = cuda_get_numDevice();
  int ret;

  if (NDEV<=0) return 2;
  if (devData==NULL) return -1;
  ret=cuda_SCF_Free();
  if(fp_prof&&ret!=0) fprintf(fp_prof, "(%d) ERROR(%d) in cuda_SCF_Free()\n", CUDA_ME, ret);

  return ret;
}
   
__host__ int cuda_SCF_Init(const int ncs, const int nat, const int maxlqn,
    const int ncspair, const int npspair, const int nao,
    const int *shel_atm, const int *shel_ini,
    const double *atom_x, const double *atom_y, const double *atom_z,
    const int *leading_cs,
    const int *leading_cs_pair, const int *csp_leading_ps_pair,
    const int *csp_ics, const int *csp_jcs,
    const double *psp_zeta, const double *psp_dkps, const double *psp_xiza,
    const float *csp_schwarz,
    const int nblk, const int nthb
    )
{
  cudaError_t err;
  int *sorted_csp;
  double *DFACT=ofmo_getadd_dfact();
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();

  if (NDEV<=0) return 2;
  if (devData==NULL) return -1;

  int maxlqn2 = (maxlqn+1)*(maxlqn+2)/2;
  int max_num_klcs = cuda_max_num_klcs(maxlqn, leading_cs_pair);
//  int num_Gkl = ncspair*nL[maxlqn]*nL[maxlqn];
  int num_Gkl = cuda_num_Gkl(maxlqn, leading_cs_pair);
  num_Gkl = ((num_Gkl-1)/16+1)*16+32; // 16words = 128bytes align;
  int num_Gi = nao * nL[maxlqn] * 2;
//  int num_Gi = cuda_num_Gi(maxlqn, leading_cs_pair, shel_ini);
  num_Gi = ((num_Gi-1)/16+1)*16+32; // 16words = 128bytes align;

  ret = cuda_set_BT(nblk, nthb, maxlqn);
  if (ret<0) return ret;

    checkCudaErrors(cudaDeviceSynchronize());
  if(fp_prof) fprintf(fp_prof, "(%d) cuda_SCF_Init()\n", CUDA_ME);
  if (!SCF_Malloced) {
    ret = cuda_SCF_Malloc(ncs, nat, maxlqn, ncspair, npspair, nao,
        max_num_klcs, num_Gkl, num_Gi);
    if (ret<0) return ret;
  }
#ifdef SORT_INDEX_SCHWARZ
  if ((sorted_csp=(int *)malloc(ncspair*sizeof(int)))==NULL) return -1;
  ret = cuda_sort_csp_schwarz(sorted_csp, ncspair, maxlqn, leading_cs_pair, csp_schwarz);
  if (ret<0) return ret;
#endif
  for (i=0; i<NDEV; i++) {
    struct dev_Data *dev = devData + i;
    cuda_set_Device(i);
    ret = cuda_SCF_Init_Symbol(ncs, nat, maxlqn, ncspair, npspair, nao,
       max_num_klcs, num_Gkl, num_Gi, dev);
    if (ret<0) return ret;
    McudaMemcpyH2D(dev->shel_atm, shel_atm, ncs*sizeof(int));
    McudaMemcpyH2D(dev->shel_ini, shel_ini, ncs*sizeof(int));
    McudaMemcpyH2D(dev->shel_ini+ncs, &nao, sizeof(int));
    McudaMemcpyH2D(dev->atom_x, atom_x, nat*sizeof(double));
    McudaMemcpyH2D(dev->atom_y, atom_y, nat*sizeof(double));
    McudaMemcpyH2D(dev->atom_z, atom_z, nat*sizeof(double));
    McudaMemcpyH2D(dev->leading_cs, leading_cs, (maxlqn+2)*sizeof(int));
    McudaMemcpyH2D(dev->leading_cs_pair, leading_cs_pair, (maxlqn2+2)*sizeof(int));
    McudaMemcpyH2D(dev->csp_leading_ps_pair, csp_leading_ps_pair, (ncspair+1)*sizeof(int));
    McudaMemcpyH2D(dev->csp_ics, csp_ics, ncspair*sizeof(int));
    McudaMemcpyH2D(dev->csp_jcs, csp_jcs, ncspair*sizeof(int));
    McudaMemcpyH2D(dev->psp_zeta, psp_zeta, npspair*sizeof(double));
    McudaMemcpyH2D(dev->psp_dkps, psp_dkps, npspair*sizeof(double));
    McudaMemcpyH2D(dev->psp_xiza, psp_xiza, npspair*sizeof(double));
    McudaMemcpyH2D(dev->DFACT, DFACT, 36*sizeof(double));
    McudaMemcpyH2D(dev->csp_schwarz, csp_schwarz, ncspair*sizeof(float));
#ifdef SORT_INDEX_SCHWARZ
    McudaMemcpyH2D(dev->sorted_csp, sorted_csp, ncspair*sizeof(int));
#endif
  }
  cuda_set_Device(0);
#ifdef SORT_INDEX_SCHWARZ
  free(sorted_csp);
#endif

//  printf("SCF_Ini:(%d/%d) %d\n",CUDA_ME,CUDA_NP,NDEV);
    checkCudaErrors(cudaDeviceSynchronize());
  //if(fp_prof) fprintf(fp_prof, "(%d) cuda_SCF_Init() done\n", CUDA_ME);
  //if (fp_prof) fflush(fp_prof);

  return 0;
}

/* ------------------------------------- */

__host__ int cuda_max_num_klcs(const int maxlqn, const int *leading_cs_pair)
{
  int La, Lb, Lab;
  int mincs, maxcs, nkl = 0;
  for ( La=0; La<=maxlqn; La++ ) {
    for ( Lb=0; Lb<=La; Lb++ ) {
      Lab = La*(La+1)/2 + Lb;
      mincs = leading_cs_pair[Lab];
      maxcs = leading_cs_pair[Lab+1];
      nkl = MAX2(nkl, maxcs-mincs);
    }
  }
  return nkl;
}

__host__ int cuda_num_Gkl(const int maxlqn, const int *leading_cs_pair)
{
  int La, Lb, Lab;
  int Lc, Ld, Lcd;
  int mincs, maxcs, nkl = 0;
  for ( La=0; La<=maxlqn; La++ ) {
    for ( Lb=0; Lb<=La; Lb++ ) {
      Lab = La*(La+1)/2 + Lb;
      for ( Lc=0; Lc<=La; Lc++ ) {
        for ( Ld=0; Ld<=(Lc==La? Lb : Lc ); Ld++ ) {
          Lcd = Lc*(Lc+1)/2 + Ld;
      mincs = leading_cs_pair[Lab];
      maxcs = leading_cs_pair[Lab+1];
      int nkl_t = (maxcs-mincs)*nL[Lc]*nL[Ld];
      nkl = MAX2(nkl, nkl_t);
      mincs = leading_cs_pair[Lcd];
      maxcs = leading_cs_pair[Lcd+1];
      nkl_t = (maxcs-mincs)*nL[La]*nL[Lb];
      nkl = MAX2(nkl, nkl_t);
        }
      }
    }
  }
  return nkl;
}

__host__ int cuda_num_Gi(const int maxlqn, const int *leading_cs_pair,
    const int *shel_ini)
{
  int La, Lb, Lab;
  int mincs, maxcs, ni = 0;
  int num_Gi, maxao;

  La = maxlqn;
  maxcs = leading_cs_pair[La+1];
  maxao = shel_ini[maxcs];
  num_Gi = maxao * nL[maxlqn] * 2;
  num_Gi = ((num_Gi-1)/16+1)*16; // 16words = 128bytes align;
  /*
  for ( La=0; La<=maxlqn; La++ ) {
    for ( Lb=0; Lb<=La; Lb++ ) {
    }
  }
  */
  return num_Gi;
}

/* ------------------------------------- */

const float *fschwarz; 

__host__ static int cuda_sort_csp_schwarz_cmp(const void *a, const void *b)
{
  int ret = 0;
  float f = fschwarz[*(int *)a] - fschwarz[*(int *)b];
  if (f>0) ret = 1;
  else if (f<0) ret = -1;
  return ret;
}

__host__ int cuda_sort_csp_schwarz(int *sorted_csp, const int ncspair,
    const int maxlqn, const int *leading_cs_pair, const float *csp_schwarz)
{
  int i;
  int La, Lb, Lab;
  int mincs, maxcs, ncs;
  int *csp;

  for (i=0;i<ncspair;i++) sorted_csp[i]=i;
  fschwarz = csp_schwarz;

  for ( La=0; La<=maxlqn; La++ ) {
    for ( Lb=0; Lb<=La; Lb++ ) {
      Lab = La*(La+1)/2 + Lb;
      mincs = leading_cs_pair[Lab];
      maxcs = leading_cs_pair[Lab+1];
      ncs = maxcs - mincs;
      qsort(&sorted_csp[mincs],ncs,sizeof(int),cuda_sort_csp_schwarz_cmp);
    }
  }
//  for (i=0;i<ncspair;i++) 
//    printf("%4d %4d %f\n",i, sorted_csp[i], csp_schwarz[sorted_csp[i]]);

  return 0;
}

/* ------------------------------------- */
#ifndef ZERO
#define ZERO 0.0e0
#endif

__global__ void gpu_Gmat_clr(void)
{
  int i;

  double *G   = G_b   + blockIdx.x * nao * nao;
  double *Gkl = Gkl_b + blockIdx.x * num_Gkl;

  __threadfence();
  for (i=threadIdx.x; i<nao*nao; i+=blockDim.x) G[i]=ZERO;
  for (i=threadIdx.x; i<num_Gkl; i+=blockDim.x) Gkl[i]=ZERO;
//  for (i=threadIdx.x+blockIdx.x*blockDim.x; i<nao*nao; i+=blockDim.x*gridDim.x) G_d[i]=ZERO;
  __syncthreads();
  __threadfence();
}

__global__ void gpu_Gmat_add(void)
{
  int i,j;
  __syncthreads();
  __threadfence();
  for (i=threadIdx.x+blockIdx.x*blockDim.x; i<nao*nao; i+=blockDim.x*gridDim.x){
    double *G = G_b + i;
    double Gtmp=ZERO;
    for (j=0; j<gridDim.x; j++,G+=nao*nao) {
      Gtmp += *G;
    }
    G_d[i] = Gtmp;
  }
  __syncthreads();
  __threadfence();
}

__host__ int cuda_genGmat_Init(const int ncs, const int nao,
    const float *Dcs, const double *D, const double x_coef0)
{
  cudaError_t err;
  int i;
  int NDEV = cuda_get_numDevice();
  int nblk = cuda_get_numBlocks();
  int nthb = cuda_get_numThreads();
  int c0[21];
  extern int counter_ini_type[];

  if (NDEV<=0) return 2;
  if (devData==NULL) return -1;
  if (!SCF_Malloced) return -2;
  for (i=0; i<21; i++) {
    c0[i] = 0;
    if      (counter_ini_type[i] == 1) c0[i] = dim2e[i][0];
    else if (counter_ini_type[i] == 2) c0[i] = dim2e[i][0]*NIJCSW;
  }

  for (i=0; i<NDEV; i++) {
    struct dev_Data *dev = devData + i;
    cuda_set_Device(i);
    McudaMemcpyH2D(dev->Dcs, Dcs, ncs*ncs*sizeof(float));
    McudaMemcpyH2D(dev->Ds, D, nao*nao*sizeof(double));
#ifdef GPU_DLB
    McudaMemcpyH2D(dev->ijcounter, c0, 21*sizeof(int));
#endif
#define McudaMemcpyToSymbol2(a,b,c) checkCudaErrors(cudaMemcpyToSymbol(a, b, c))
    McudaMemcpyToSymbol2(x_coef, &x_coef0, sizeof(double));
#undef McudaMemcpyToSymbol2

    gpu_Gmat_clr <<< nblk, nthb >>> ();
    //cudaDeviceSynchronize();
    checkCudaErrors(cudaDeviceSynchronize());
  }
  cuda_set_Device(0);
//  printf("Ini:(%d/%d) %16.10g %16.10g\n",CUDA_ME,CUDA_NP,Dcs[0],D[0]);

  return 0;
}

__host__ int cuda_genGmat_Add(const int nao, double *G)
{
  cudaError_t err;
  int i,j;
  double *Gtmp;
  int NDEV = cuda_get_numDevice();
  int nblk = cuda_get_numBlocks();
  int nthb = cuda_get_numThreads();

  if (NDEV<=0) return 2;
  if (devData==NULL) return -1;
  if (!SCF_Malloced) return -2;

  Gtmp = (double *)malloc(nao*nao*sizeof(double));
  if (Gtmp==NULL) return -3;

  for (i=0; i<NDEV; i++) {
    struct dev_Data *dev = devData + i;
    cuda_set_Device(i);
    //cudaDeviceSynchronize();
    checkCudaErrors(cudaDeviceSynchronize());

    gpu_Gmat_add <<< nblk, nthb >>> ();

    cudaDeviceSynchronize();
    McudaMemcpyD2H(Gtmp, dev->G_d, nao*nao*sizeof(double));
//  printf("Add:(%d/%d) %16.10g\n",CUDA_ME,CUDA_NP,Gtmp[0]);
    for (j=0; j<nao*nao; j++) G[j] += Gtmp[j];
  }
  cuda_set_Device(0);
  free(Gtmp);

  return 0;
}

/* ------------------------------------- */
