/**
 * fmt-m.c
 *
 * */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include <cuda.h>
#include "cudalib.h"
#include "cuda-fmt-m.h"

#ifdef __cplusplus
extern "C" {
#endif
extern FILE* fp_prof; // from common/ofmo-prof.h
#ifdef __cplusplus
}
#endif

static int FMT_m_Initialized = -1;

/* ------------------------------------- */

__device__ __constant__ double *FMT_m_table[9];
__device__ __constant__ double *FMT_m_table0;
__device__ __constant__ double *FMT_m_table1;
__device__ __constant__ double *FMT_m_table2;
__device__ __constant__ double *FMT_m_table3;
__device__ __constant__ double *FMT_m_table4;
__device__ __constant__ double *FMT_m_table5;
__device__ __constant__ double *FMT_m_table6;
__device__ __constant__ double *FMT_m_table7;
__device__ __constant__ double *FMT_m_table8;
__device__ __constant__ size_t FMT_m_size[9];
__device__ __constant__ double FMT_m_delta[9]; // mmax+1
__device__ __constant__ double FMT_m_dinv[11]; // ??
__device__ __constant__ double FMT_m_dinv2[8+10]; // mmax+nexp

__device__ __constant__ double FMT_m_sqrt_pi_2;
__device__ __constant__ double FMT_m_sqrt_pi2;

/* ------------------------------------- */

static size_t fmt_m_size[9];

__host__ size_t cuda_FMT_m_get_size(const int m)
{
  size_t size = 0;
  int mmax = FMT_m_Initialized;
  if (m>=0 && m<=mmax) size = fmt_m_size[m];
  return size;
}

/* ------------------------------------- */
__host__ int cuda_FMT_m_Init(double **dtmp, const size_t *mtmp,
    const int *ndivs, const int mmax)
{
  cudaError_t err;
  int i;

  if (mmax>8) exit(-1);

  int NDEV = cuda_get_numDevice();

  if (NDEV<=0) return 2;
  if (FMT_m_Initialized>mmax) return 1;
  if(fp_prof) fprintf(fp_prof, "cuda_FMT_m_Init(): %d\n",mmax);
  for (i=0; i<NDEV; i++) {
    struct dev_Data *dev;
    cuda_set_Device(i);
    dev = cuda_SCF_get_dev_Data(i);
    dev->fmt_m_mmax = mmax;
    for (int m=0; m<=mmax; m++) {
      err = cudaMalloc((void **)&(dev->fmt_m_table[m]), mtmp[m]*sizeof(double));
      if (err!=cudaSuccess) exit(1);
      McudaMemcpyH2D(dev->fmt_m_table[m], dtmp[m], mtmp[m]*sizeof(double));
    }
    err = cudaMemcpyToSymbol(FMT_m_table, dev->fmt_m_table, (mmax+1)*sizeof(double *));
    if (err!=cudaSuccess) exit(1);
    err = cudaMemcpyToSymbol(FMT_m_table0, &(dev->fmt_m_table[0]), sizeof(double *));
    if (err!=cudaSuccess) exit(1);
    err = cudaMemcpyToSymbol(FMT_m_table1, &(dev->fmt_m_table[1]), sizeof(double *));
    if (err!=cudaSuccess) exit(1);
    err = cudaMemcpyToSymbol(FMT_m_table2, &(dev->fmt_m_table[2]), sizeof(double *));
    if (err!=cudaSuccess) exit(1);
    err = cudaMemcpyToSymbol(FMT_m_table3, &(dev->fmt_m_table[3]), sizeof(double *));
    if (err!=cudaSuccess) exit(1);
    err = cudaMemcpyToSymbol(FMT_m_table4, &(dev->fmt_m_table[4]), sizeof(double *));
    if (err!=cudaSuccess) exit(1);
    err = cudaMemcpyToSymbol(FMT_m_table5, &(dev->fmt_m_table[5]), sizeof(double *));
    if (err!=cudaSuccess) exit(1);
    err = cudaMemcpyToSymbol(FMT_m_table6, &(dev->fmt_m_table[6]), sizeof(double *));
    if (err!=cudaSuccess) exit(1);
    err = cudaMemcpyToSymbol(FMT_m_table7, &(dev->fmt_m_table[7]), sizeof(double *));
    if (err!=cudaSuccess) exit(1);
    err = cudaMemcpyToSymbol(FMT_m_table8, &(dev->fmt_m_table[8]), sizeof(double *));
    if (err!=cudaSuccess) exit(1);
    err = cudaMemcpyToSymbol(FMT_m_size, mtmp, (mmax+1)*sizeof(size_t));
    if (err!=cudaSuccess) exit(1);
    {
      double *delta = (double *)malloc((mmax+1)*sizeof(double));
      for (int m=0; m<=mmax; m++) delta[m] = 1e0/(double)ndivs[m];
      err = cudaMemcpyToSymbol(FMT_m_delta, delta, (mmax+1)*sizeof(double));
      if (err!=cudaSuccess) exit(1);
      free(delta);
    }
    {
      int ndinv=11;
      double *dinv = (double *)malloc(ndinv*sizeof(double));
      dinv[0] = dinv[1] = 1.e0;
      for (int k=2; k<ndinv; k++ ) dinv[k] = 1.e0 / (double)k;
      err = cudaMemcpyToSymbol(FMT_m_dinv, dinv, ndinv*sizeof(double));
      if (err!=cudaSuccess) exit(1);
      free(dinv);
    }
    {
      int ndinv2=mmax+10;
      double *dinv2 = (double *)malloc(ndinv2*sizeof(double));
      dinv2[0] = 1.e0;
      for (int k=1; k<ndinv2; k++ ) dinv2[k] = 1.e0/ (double)(2*k+1);
      err = cudaMemcpyToSymbol(FMT_m_dinv2, dinv2, ndinv2*sizeof(double));
      if (err!=cudaSuccess) exit(1);
      free(dinv2);
    }
    {
      double sqrt_pi_2 = sqrt( atan( 1.e0 ) );
      double sqrt_pi2  = sqrt( 2.e0*atan( 1.e0 ) );
      err = cudaMemcpyToSymbol(FMT_m_sqrt_pi_2, &sqrt_pi_2, sizeof(double));
      if (err!=cudaSuccess) exit(1);
      err = cudaMemcpyToSymbol(FMT_m_sqrt_pi2, &sqrt_pi2, sizeof(double));
      if (err!=cudaSuccess) exit(1);
    }
  }
  for (int m=0; m<=mmax; m++) fmt_m_size[m] = mtmp[m];
  FMT_m_Initialized = mmax;

  return 0;
}

__host__ int cuda_FMT_m_Finalize(void)
{
  cudaError_t err;
  int i;
  int NDEV = cuda_get_numDevice();

  if (NDEV<=0) return 2;
  if (FMT_m_Initialized<0) return 1;
  if(fp_prof) fprintf(fp_prof, "cuda_FMT_m_Finalize():\n");
  for (i=0; i<NDEV; i++) {
    struct dev_Data *dev = cuda_SCF_get_dev_Data(i);
    int mmax = dev->fmt_m_mmax;
    cuda_set_Device(i);
    for (int m=0; m<=mmax; m++) { 
      err = cudaFree(dev->fmt_m_table[m]);
      if (err!=cudaSuccess) exit(1);
    }
  }
  FMT_m_Initialized = -1;

  return 0;
}

/* ------------------------------------- */

#if   CUDA_FMT_M_NEXP == 6
#include "cuda-fmt-method1-8-6-12.cu"
#include "cuda-fmt-method2-8-6-12.cu"
#include "cuda-fmt-method3-8-6-12.cu"
#elif CUDA_FMT_M_NEXP == 8
#include "cuda-fmt-method1-8-8-12.cu"
#include "cuda-fmt-method2-8-8-12.cu"
#include "cuda-fmt-method3-8-8-12.cu"
#elif CUDA_FMT_M_NEXP == 10
#include "cuda-fmt-method1-8-10-12.cu"
#include "cuda-fmt-method2-8-10-12.cu"
#include "cuda-fmt-method3-8-10-12.cu"
#endif

/* ------------------------------------- */

