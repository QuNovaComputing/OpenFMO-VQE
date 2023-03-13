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

#include "cuda-drv.h"
#include "cudalib.h"
#include "cuda-ifc4c.h"
#include "cuda-ifc4c-os.h"

#ifdef __cplusplus
extern "C" {
#endif
//#include "integ/ofmo-twoint-direct.h"
//#include "integ/ofmo-index.h"
extern FILE* fp_prof; // from common/ofmo-prof.h
#ifdef __cplusplus
}
#endif

#ifndef ZERO
#define ZERO 0.0e0
#endif
static int NCS_F;
static int NCS_M;

static struct {
  int maxlqn;
  int max_num_klcs;
  int nat_f;
  int ncs_f;
  int nao_f;
  int ncspair_f;
  int npspair_f;
  int nat_m;
  int ncs_m;
  int nao_m;
  int ncspair_m;
  int npspair_m;
} ifc4cParam;

/* ------------------------------------- */
__device__ __constant__ int ncs_frg;
__device__ __constant__ int nat_frg;
__device__ __constant__ int nao_frg;
__device__ __constant__ int ncspair_frg;
__device__ __constant__ int npspair_frg;
__device__ __constant__ int ncs_mon;
__device__ __constant__ int nat_mon;
__device__ __constant__ int nao_mon;
__device__ __constant__ int ncspair_mon;
__device__ __constant__ int npspair_mon;

__device__ __constant__ int *shel_atm_frg;
__device__ __constant__ int *shel_ini_frg;
__device__ __constant__ double *atom_x_frg;
__device__ __constant__ double *atom_y_frg;
__device__ __constant__ double *atom_z_frg;
__device__ __constant__ int *leading_cs_pair_frg;
__device__ __constant__ int *leading_cs_frg;
__device__ __constant__ int *csp_leading_ps_pair_frg;
__device__ __constant__ int *csp_ics_frg;
__device__ __constant__ int *csp_jcs_frg;
__device__ __constant__ float *csp_schwarz_frg;
__device__ __constant__ double *psp_zeta_frg;
__device__ __constant__ double *psp_dkps_frg;
__device__ __constant__ double *psp_xiza_frg;
__device__ __constant__ int *shel_atm_mon;
__device__ __constant__ int *shel_ini_mon;
__device__ __constant__ double *atom_x_mon;
__device__ __constant__ double *atom_y_mon;
__device__ __constant__ double *atom_z_mon;
__device__ __constant__ int *leading_cs_pair_mon;
__device__ __constant__ int *leading_cs_mon;
__device__ __constant__ int *csp_leading_ps_pair_mon;
__device__ __constant__ int *csp_ics_mon;
__device__ __constant__ int *csp_jcs_mon;
__device__ __constant__ float *csp_schwarz_mon;
__device__ __constant__ double *psp_zeta_mon;
__device__ __constant__ double *psp_dkps_mon;
__device__ __constant__ double *psp_xiza_mon;
__device__ __constant__ double *D_mon;
__device__ __constant__ double *V_frg;
__device__ __constant__ double *V_frgP;

/* ------------------------------------- */
static int ifc4c_Initialized = FALSE;

static double *Vifc4c = NULL;

#if 0
#define cudaDevMallocAndSetSymbol(name, num, kind) \
  {ret--; \
   if (cudaMalloc((void **)&(dev->name), (num)*sizeof(kind))!=cudaSuccess) break; \
   if (cudaMemcpyToSymbol(name, &(dev->name), sizeof(dev->name))!=cudaSuccess) break;}
#else
#define McudaDevMallocAndSetSymbol(name, num, kind) \
  {ret--;\
   if (dev->name!=NULL) {ret-=1000; break;}; \
   if ((err=cudaMalloc((void **)&(dev->name), (num)*sizeof(kind)))!=cudaSuccess) break; \
   if ((err=cudaMemcpyToSymbol(name, &(dev->name), sizeof(dev->name)))!=cudaSuccess) break;}
#endif

// arg: maximum value to malloc
__host__ int cuda_ifc4c_Init(const int maxlqn, const int max_num_klcs0,
    const int nat_f, const int ncs_f, const int nao_f,
    const int ncspair_f, const int npspair_f,
    const int nat_m, const int ncs_m, const int nao_m,
    const int ncspair_m, const int npspair_m)
{
  char fn[]="cuda_ifc4c_Init";
  cudaError_t err;
  int i,ret;
  int NDEV = cuda_get_numDevice();
  int maxlqn2 = (maxlqn+1)*(maxlqn+2)/2;

  if (NDEV<=0) return 2;
  if (devData==NULL) return -1;
  if (ifc4c_Initialized) return 1;
  if(fp_prof) fprintf(fp_prof, "(%d) %s(%d,%d)\n", CUDA_ME, fn, nao_f, nao_m);

  //ret = cuda_set_BT(nblk, nthb, maxlqn);
  //if (ret<0) return ret;
  //int nblk = cuda_get_numBlocks();
  //int nthb = cuda_get_numThreads();
  int nblk=0, nthb=0;
  for (i=0;i<6*6;i++) nblk=MAX2(nblk, dim_ifc4c[i][0]);

  double *DFACT0=ofmo_getadd_dfact();

  for (i=0; i<NDEV; i++) {
    ret = -(i+1)*100;
    struct dev_Data *dev = devData + i;
    cuda_set_Device(i);
    cudaMemcpyToSymbol(max_num_klcs, &max_num_klcs0, sizeof(int));
    McudaDevMallocAndSetSymbol(sklcs_b, max_num_klcs0*nblk, int);
    McudaDevMallocAndSetSymbol(shel_atm_frg, ncs_f, int);
    McudaDevMallocAndSetSymbol(shel_ini_frg, ncs_f+1, int);
    McudaDevMallocAndSetSymbol(atom_x_frg, nat_f, double);
    McudaDevMallocAndSetSymbol(atom_y_frg, nat_f, double);
    McudaDevMallocAndSetSymbol(atom_z_frg, nat_f, double);
    McudaDevMallocAndSetSymbol(leading_cs_frg, maxlqn+2, int);
    McudaDevMallocAndSetSymbol(leading_cs_pair_frg, maxlqn2+2, int);
    McudaDevMallocAndSetSymbol(csp_leading_ps_pair_frg, ncspair_f+1, int);
    McudaDevMallocAndSetSymbol(csp_ics_frg, ncspair_f, int);
    McudaDevMallocAndSetSymbol(csp_jcs_frg, ncspair_f, int);
    McudaDevMallocAndSetSymbol(psp_zeta_frg, npspair_f, double);
    McudaDevMallocAndSetSymbol(psp_dkps_frg, npspair_f, double);
    McudaDevMallocAndSetSymbol(psp_xiza_frg, npspair_f, double);
    McudaDevMallocAndSetSymbol(csp_schwarz_frg, ncspair_f, float);
    McudaDevMallocAndSetSymbol(shel_atm_mon, ncs_m, int);
    McudaDevMallocAndSetSymbol(shel_ini_mon, ncs_m+1, int);
    McudaDevMallocAndSetSymbol(atom_x_mon, nat_m, double);
    McudaDevMallocAndSetSymbol(atom_y_mon, nat_m, double);
    McudaDevMallocAndSetSymbol(atom_z_mon, nat_m, double);
    McudaDevMallocAndSetSymbol(leading_cs_mon, maxlqn+2, int);
    McudaDevMallocAndSetSymbol(leading_cs_pair_mon, maxlqn2+2, int);
    McudaDevMallocAndSetSymbol(csp_leading_ps_pair_mon, ncspair_m+1, int);
    McudaDevMallocAndSetSymbol(csp_ics_mon, ncspair_m, int);
    McudaDevMallocAndSetSymbol(csp_jcs_mon, ncspair_m, int);
    McudaDevMallocAndSetSymbol(psp_zeta_mon, npspair_m, double);
    McudaDevMallocAndSetSymbol(psp_dkps_mon, npspair_m, double);
    McudaDevMallocAndSetSymbol(psp_xiza_mon, npspair_m, double);
    McudaDevMallocAndSetSymbol(csp_schwarz_mon, ncspair_m, float);
    McudaDevMallocAndSetSymbol(D_mon, (nao_m*nao_m+nao_m)/2, double);
    McudaDevMallocAndSetSymbol(V_frg, (nao_f*nao_f+nao_f)/2, double);
    McudaDevMallocAndSetSymbol(V_frgP, (nao_f*nao_f+nao_f)/2, double);
    McudaDevMallocAndSetSymbol(Dcs, ncs_m*ncs_m, float);
    McudaDevMallocAndSetSymbol(ijcounter, 6*6, int);
#ifdef SORT_INDEX_SCHWARZ
    McudaDevMallocAndSetSymbol(sorted_csp, ncspair_f, int);
#endif
    McudaDevMallocAndSetSymbol(DFACT, 36, double);
    ret --;
    McudaMemcpyH2D(dev->DFACT, DFACT0, 36*sizeof(double));
    ret = 0;
  }
  if (ret<0 && fp_prof) {
    if (ret<-1000) {
      fprintf(fp_prof, "(%d) %s() %d: may not clean pointer!\n", CUDA_ME, fn, ret);
    } else {
      fprintf(fp_prof, "(%d) %s() %d: %s\n", CUDA_ME, fn, ret, cudaGetErrorString(err));
    }
    fprintf(fp_prof, "(%d) %s(%d,%d,\n"
                     "        %d,%d,%d,%d,%d,\n"
                     "        %d,%d,%d,%d,%d)\n", CUDA_ME, fn,
                     maxlqn, max_num_klcs0,
                     nat_f, ncs_f, nao_f, ncspair_f, npspair_f,
                     nat_m, ncs_m, nao_m, ncspair_m, npspair_m);
    fflush(fp_prof);
    exit(ret);
  }
  assert(Vifc4c==NULL);
  if ((Vifc4c=(double *)malloc((nao_f*nao_f+nao_f)/2*sizeof(double)))==NULL) exit(-1);

  NCS_F=ncs_f;
  NCS_M=ncs_m;
  ifc4c_Initialized = TRUE;

  ifc4cParam.maxlqn = maxlqn;
  ifc4cParam.max_num_klcs= max_num_klcs0;
  ifc4cParam.nat_f = nat_f;
  ifc4cParam.ncs_f = ncs_f;
  ifc4cParam.nao_f = nao_f;
  ifc4cParam.ncspair_f = ncspair_f;
  ifc4cParam.npspair_f = npspair_f;
  ifc4cParam.nat_m = nat_m;
  ifc4cParam.ncs_m = ncs_m;
  ifc4cParam.nao_m = nao_m;
  ifc4cParam.ncspair_m = ncspair_m;
  ifc4cParam.npspair_m = npspair_m;

  return 0;
}

#undef McudaDevMallocAndSetSymbol

#define McudaDevFree(name) \
  {ret--; \
    if ((err=cudaFree(dev->name))!=cudaSuccess) break; dev->name = NULL;}

__host__ int cuda_ifc4c_Finalize(void)
{
  char fn[]="cuda_ifc4c_Finalize";
  cudaError_t err;
  int i,ret;
  int NDEV = cuda_get_numDevice();

  if (NDEV<=0) return 2;
  if (devData==NULL) return -1;
  if (!ifc4c_Initialized) return 1;
  if(fp_prof) fprintf(fp_prof, "(%d) %s()\n", CUDA_ME, fn);

  for (i=0; i<NDEV; i++) {
    ret = -(i+1)*100;
    struct dev_Data *dev = devData + i;

    cuda_set_Device(i);
    McudaDevFree(sklcs_b);
    McudaDevFree(shel_atm_frg);
    McudaDevFree(shel_ini_frg);
    McudaDevFree(atom_x_frg);
    McudaDevFree(atom_y_frg);
    McudaDevFree(atom_z_frg);
    McudaDevFree(leading_cs_frg);
    McudaDevFree(leading_cs_pair_frg);
    McudaDevFree(csp_leading_ps_pair_frg);
    McudaDevFree(csp_ics_frg);
    McudaDevFree(csp_jcs_frg);
    McudaDevFree(psp_zeta_frg);
    McudaDevFree(psp_dkps_frg);
    McudaDevFree(psp_xiza_frg);
    McudaDevFree(csp_schwarz_frg);
    McudaDevFree(shel_atm_mon);
    McudaDevFree(shel_ini_mon);
    McudaDevFree(atom_x_mon);
    McudaDevFree(atom_y_mon);
    McudaDevFree(atom_z_mon);
    McudaDevFree(leading_cs_mon);
    McudaDevFree(leading_cs_pair_mon);
    McudaDevFree(csp_leading_ps_pair_mon);
    McudaDevFree(csp_ics_mon);
    McudaDevFree(csp_jcs_mon);
    McudaDevFree(psp_zeta_mon);
    McudaDevFree(psp_dkps_mon);
    McudaDevFree(psp_xiza_mon);
    McudaDevFree(csp_schwarz_mon);
    McudaDevFree(D_mon);
    McudaDevFree(V_frg);
    McudaDevFree(V_frgP);
    McudaDevFree(Dcs);
    McudaDevFree(ijcounter);
#ifdef SORT_INDEX_SCHWARZ
    McudaDevFree(sorted_csp);
#endif
    McudaDevFree(DFACT);
    ret = 0;
  }
  if (ret<0 && fp_prof) {
    fprintf(fp_prof, "(%d) %s(): %d\n", CUDA_ME, fn, ret);
    fflush(fp_prof);
    exit(ret);
  }

  free(Vifc4c); Vifc4c = NULL;
  ifc4c_Initialized = FALSE;

  return 0;
}
#undef McudaDevFree

/* ------------------------------------- */
#define McudaMemcpyToSymbolH2D(dst, src, kind) \
  {ret--; \
   if (cudaMemcpyToSymbol(dst, src, sizeof(kind))!=cudaSuccess) return ret;}

__host__ int cuda_ifc4c_SetData_Symbol(const int iorj, const int maxlqn0,
    const int nat0, const int ncs0, const int nao0,
    const int ncspair0, const int npspair0)
{
  cudaError_t err;
  int ret=-100*(iorj+1);

  if (!ifc4c_Initialized) return -1;

  if (iorj==0) {
    McudaMemcpyToSymbolH2D(ncs_frg, &ncs0, int);
    McudaMemcpyToSymbolH2D(nao_frg, &nao0, int);
    McudaMemcpyToSymbolH2D(nat_frg, &nat0, int);
    McudaMemcpyToSymbolH2D(ncspair_frg, &ncspair0, int);
    McudaMemcpyToSymbolH2D(npspair_frg, &npspair0, int);
    McudaMemcpyToSymbolH2D(maxlqn, &maxlqn0, int);
  } else {
    McudaMemcpyToSymbolH2D(ncs_mon, &ncs0, int);
    McudaMemcpyToSymbolH2D(nao_mon, &nao0, int);
    McudaMemcpyToSymbolH2D(nat_mon, &nat0, int);
    McudaMemcpyToSymbolH2D(ncspair_mon, &ncspair0, int);
    McudaMemcpyToSymbolH2D(npspair_mon, &npspair0, int);
  }

  return 0;
}
#undef McudaMemcpyToSymbolH2D

#if 1
#define McudaDevMemcpyH2D(name, frg, num, kind) \
  {ret--; narg=(num); \
   if (cudaMemcpy(dev->name##_##frg, name, (num)*sizeof(kind),cudaMemcpyHostToDevice)!=cudaSuccess) break;}
#else
#define McudaDevMemcpyH2D(name, frg, num, kind) \
    checkCudaErrors(cudaMemcpy(dev->name##_##frg, name, (num)*sizeof(kind),cudaMemcpyHostToDevice));
#endif

__host__ int cuda_ifc4c_SetData(const int iorj, const int maxlqn,
    const int nat, const int ncs, const int nao,
    const int ncspair, const int npspair,
    const int *shel_atm, const int *shel_ini,
    const double *atom_x, const double *atom_y, const double *atom_z,
    const int *leading_cs_pair, const int *csp_leading_ps_pair,
    const int *csp_ics, const int *csp_jcs,
    const double *psp_zeta, const double *psp_dkps, const double *psp_xiza,
    const float *csp_schwarz, const double *D)
{
  char fn[]="cuda_ifc4c_SetData";
  cudaError_t err;
  int i,ret,narg;
  int NDEV = cuda_get_numDevice();
  int maxlqn2 = (maxlqn+1)*(maxlqn+2)/2;

  if (NDEV<=0) return 2;
  if (devData==NULL) return -1;
  if (!ifc4c_Initialized) return -1;

#ifdef SORT_INDEX_SCHWARZ
  int *sorted_csp;
  if (iorj==0) {
    if ((sorted_csp=(int *)malloc(ncspair*sizeof(int)))==NULL) return -1;
    ret = cuda_sort_csp_schwarz(sorted_csp, ncspair, maxlqn, leading_cs_pair, csp_schwarz);
    if (ret<0) return ret;
  }
#endif

  for (i=0; i<NDEV; i++) {
    struct dev_Data *dev = devData + i;
    cuda_set_Device(i);
    checkCudaErrors(cudaDeviceSynchronize());
    ret =cuda_ifc4c_SetData_Symbol(iorj, maxlqn, nat, ncs, nao,
        ncspair, npspair);
    if (ret<0) return ret;
    ret = -(i+1)*1000 - iorj*100;
    if (iorj==0) {
//      for (int j=0;j<(nao*nao+nao)/2;j++) Vifc4c[j]=0e0;
//      err = cudaMemcpyH2D(dev->V_frg, Vifc4c, (nao*nao+nao)/2*sizeof(double));
      gpu_ifc4c_ClearVfrg <<< dev->numSM, 768 >>> (nao, dev->V_frg);
      gpu_ifc4c_ClearVfrg <<< dev->numSM, 768 >>> (nao, dev->V_frgP);
      McudaDevMemcpyH2D(shel_atm, frg, ncs, int);
      McudaDevMemcpyH2D(shel_ini, frg, ncs+1, int);
      McudaDevMemcpyH2D(atom_x, frg, nat, double);
      McudaDevMemcpyH2D(atom_y, frg, nat, double);
      McudaDevMemcpyH2D(atom_z, frg, nat, double);
      McudaDevMemcpyH2D(leading_cs_pair, frg, maxlqn2+2, int);
      McudaDevMemcpyH2D(csp_leading_ps_pair, frg, ncspair+1, int);
      McudaDevMemcpyH2D(csp_ics, frg, ncspair, int);
      McudaDevMemcpyH2D(csp_jcs, frg, ncspair, int);
      McudaDevMemcpyH2D(psp_zeta, frg, npspair, double);
      McudaDevMemcpyH2D(psp_dkps, frg, npspair, double);
      McudaDevMemcpyH2D(psp_xiza, frg, npspair, double);
      McudaDevMemcpyH2D(csp_schwarz, frg, ncspair, float);
#ifdef SORT_INDEX_SCHWARZ
      ret--;
      McudaMemcpyH2D(dev->sorted_csp, sorted_csp, ncspair*sizeof(int));
#endif
      checkCudaErrors(cudaDeviceSynchronize());
    } else {
      McudaDevMemcpyH2D(shel_atm, mon, ncs, int);
      McudaDevMemcpyH2D(shel_ini, mon, ncs+1, int);
      McudaDevMemcpyH2D(atom_x, mon, nat, double);
      McudaDevMemcpyH2D(atom_y, mon, nat, double);
      McudaDevMemcpyH2D(atom_z, mon, nat, double);
      McudaDevMemcpyH2D(leading_cs_pair, mon, maxlqn2+2, int);
      McudaDevMemcpyH2D(csp_leading_ps_pair, mon, ncspair+1, int);
      McudaDevMemcpyH2D(csp_ics, mon, ncspair, int);
      McudaDevMemcpyH2D(csp_jcs, mon, ncspair, int);
      McudaDevMemcpyH2D(psp_zeta, mon, npspair, double);
      McudaDevMemcpyH2D(psp_dkps, mon, npspair, double);
      McudaDevMemcpyH2D(psp_xiza, mon, npspair, double);
      McudaDevMemcpyH2D(csp_schwarz, mon, ncspair, float);
      McudaDevMemcpyH2D(D, mon, (nao*nao+nao)/2, double);
    }
    ret = 0;
  }
#ifdef SORT_INDEX_SCHWARZ
  if (iorj==0) free(sorted_csp);
#endif
  if (ret<0 && fp_prof) {
    //fprintf(fp_prof, "(%d) %s(): %d\n", CUDA_ME, fn, ret);
    fprintf(fp_prof, "(%d) %s(): %d %d\n", CUDA_ME, fn, ret, narg);
    fprintf(fp_prof, "(%d) %s(): %d %d %d %d %d\n", CUDA_ME, fn, nat, ncs, nao, ncspair, npspair);
    if (iorj==0) fprintf(fp_prof, "(%d) %s(): %d %d %d %d %d\n", CUDA_ME, fn, ifc4cParam.nat_f, ifc4cParam.ncs_f, ifc4cParam.nao_f, ifc4cParam.ncspair_f, ifc4cParam.npspair_f);
    else         fprintf(fp_prof, "(%d) %s(): %d %d %d %d %d\n", CUDA_ME, fn, ifc4cParam.nat_m, ifc4cParam.ncs_m, ifc4cParam.nao_m, ifc4cParam.ncspair_m, ifc4cParam.npspair_m);

    fflush(fp_prof);
    exit(ret);
  }
  //if (iorj==0) assert(ncs<=NCS_F);
  //else         assert(ncs<=NCS_M);
  //if (iorj==0) NCS_F=ncs;
  //else         NCS_M=ncs;
  return 0;
}

#undef McudaDevMemcpyH2D

/* ------------------------------------- */

__global__ void gpu_ifc4c_ClearVfrg(const int nao, double *V)
{
  int tidx = blockIdx.x * blockDim.x + threadIdx.x;
  int nth  = gridDim.x * blockDim.x;
//  int nao = nao_frg;

  __threadfence();
  for (int i=tidx; i<(nao*nao+nao)/2; i+=nth) V[i]=ZERO;
  __syncthreads();
  __threadfence();
}

__global__ void gpu_ifc4c_PsumVfrg(const int nao, double *V, double *VP)
{
  int tidx = blockIdx.x * blockDim.x + threadIdx.x;
  int nth  = gridDim.x * blockDim.x;

  __threadfence();
  for (int i=tidx; i<(nao*nao+nao)/2; i+=nth) {
    VP[i] += V[i];
    V[i]   = ZERO;
  }
  __syncthreads();
  __threadfence();
}


/* ------------------------------------- */

__host__ int cuda_ifc4c_calc_Init(void)
{
  cudaError_t err;
  int NDEV = cuda_get_numDevice();
  int c0[6*6];

  if (NDEV<=0) return 2;
  if (devData==NULL) return -1;
  if (!ifc4c_Initialized) return -1;
  for (int i=0; i<6*6; i++) c0[i] = dim_ifc4c[i][0];

  for (int i=0; i<NDEV; i++) {
    struct dev_Data *dev = devData + i;
    cuda_set_Device(i);
    McudaMemcpyH2D(dev->ijcounter, c0, 6*6*sizeof(int));
  }
  return 0;
}

/* ------------------------------------- */

__host__ int cuda_ifc4c_SetDcs(const int ncs, const float *Dcs)
{
  cudaError_t err;
  int NDEV = cuda_get_numDevice();

  if (NDEV<=0) return 2;
  if (devData==NULL) return -1;
  if (!ifc4c_Initialized) return -1;
  char fn[]="cuda_ifc4c_SetDcs";
  //assert(ncs<=NCS_M);

  for (int i=0; i<NDEV; i++) {
    struct dev_Data *dev = devData + i;
    cuda_set_Device(i);
    McudaMemcpyH2D(dev->Dcs, Dcs, ncs*ncs*sizeof(float));
  }
  return 0;
}

/* ------------------------------------- */

__host__ int cuda_ifc4c_GetVfrg(const int nao, double *V)
{
  char fn[]="cuda_ifc4c_GetVfrg";
  cudaError_t err;
  int NDEV = cuda_get_numDevice();
  int nao2=(nao*nao+nao)/2;

  if (NDEV<=0) return 2;
  if (devData==NULL) return -1;
  if (!ifc4c_Initialized) return -1;

  for (int i=0; i<NDEV; i++) {
    struct dev_Data *dev = devData + i;
    cuda_set_Device(i);
    checkCudaErrors(cudaDeviceSynchronize());
    gpu_ifc4c_PsumVfrg <<< dev->numSM, 768 >>> (nao, dev->V_frgP, dev->V_frg);
    McudaMemcpyD2H(Vifc4c, dev->V_frg, nao2*sizeof(double));
    for (int j=0; j<nao2; j++) V[j]+=Vifc4c[j];
  }
  return 0;
}

/* ------------------------------------- */

__host__ int cuda_ifc4c_PsumVfrg(const int nao)
{
  char fn[]="cuda_ifc4c_PsumVfrg";
  cudaError_t err;
  int NDEV = cuda_get_numDevice();
  int nao2=(nao*nao+nao)/2;

  if (NDEV<=0) return 2;
  if (devData==NULL) return -1;
  if (!ifc4c_Initialized) return -1;

  for (int i=0; i<NDEV; i++) {
    struct dev_Data *dev = devData + i;
    cuda_set_Device(i);
    checkCudaErrors(cudaDeviceSynchronize());
    gpu_ifc4c_PsumVfrg <<< dev->numSM, 768 >>> (nao, dev->V_frg, dev->V_frgP);
  }
  return 0;
}
