/**
 * @file fmt.c
 * 各種分子積分計算で用いる誤差関数の計算に関する関数群
 *
 * */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include <cuda.h>
#include "cudalib.h"
#include "cuda-fmt.h"

#ifdef __cplusplus
extern "C" {
#endif
extern FILE* fp_prof; // from common/ofmo-prof.h
#ifdef __cplusplus
}
#endif

static int FMT_Initialized = -1;

/* ------------------------------------- */

__device__ __constant__ double *FMT_fmt_table;
__device__ __constant__ double FMT_fmt_step_size;
__device__ __constant__ double FMT_fmt_inv_step_size;
__device__ __constant__ double FMT_pi_div2;
__device__ __constant__ int FMT_fmt_max_m1;

__device__ __constant__ double *FMT_fmt_table0;
__device__ __constant__ double *FMT_fmt_table1;
__device__ __constant__ double *FMT_fmt_table2;
__device__ __constant__ double *FMT_fmt_table3;
__device__ __constant__ double *FMT_fmt_table4;
__device__ __constant__ double *FMT_fmt_table5;
__device__ __constant__ double *FMT_fmt_table6;
__device__ __constant__ double *FMT_fmt_table7;
__device__ __constant__ double *FMT_fmt_table8;

/* ------------------------------------- */
__host__ int cuda_FMT_Init(double **dtmp, size_t *mtmp,
    const double step_size, const int max_m1)
{
  cudaError_t err;
  int i;
  double inv_step_size = 1.0/step_size;
  double pi_div2 = 2.0 * atan(1.0);
  int NDEV = cuda_get_numDevice();

  if (NDEV<=0) return 2;
  if (FMT_Initialized>=max_m1) return 1;
  if(fp_prof) fprintf(fp_prof, "cuda_FMT_Init(): %d\n",max_m1);
  for (i=0; i<NDEV; i++) {
    struct dev_Data *dev;
    cuda_set_Device(i);
    dev = cuda_SCF_get_dev_Data(i);
    if (FMT_Initialized<0) { // 1st time
      for (int m=0; m<=9; m++) {
        err = cudaMalloc((void **)&(dev->fmt_table[m]), mtmp[m]*sizeof(double));
        checkCudaErrors(err);
        if (err!=cudaSuccess) exit(1);
      }
      err = cudaMemcpyToSymbol(FMT_fmt_table0, &(dev->fmt_table[0]), sizeof(double *));
      checkCudaErrors(err);
      if (err!=cudaSuccess) exit(1);
      err = cudaMemcpyToSymbol(FMT_fmt_table1, &(dev->fmt_table[1]), sizeof(double *));
      checkCudaErrors(err);
      if (err!=cudaSuccess) exit(1);
      err = cudaMemcpyToSymbol(FMT_fmt_table2, &(dev->fmt_table[2]), sizeof(double *));
      checkCudaErrors(err);
      if (err!=cudaSuccess) exit(1);
      err = cudaMemcpyToSymbol(FMT_fmt_table3, &(dev->fmt_table[3]), sizeof(double *));
      checkCudaErrors(err);
      if (err!=cudaSuccess) exit(1);
      err = cudaMemcpyToSymbol(FMT_fmt_table4, &(dev->fmt_table[4]), sizeof(double *));
      checkCudaErrors(err);
      if (err!=cudaSuccess) exit(1);
      err = cudaMemcpyToSymbol(FMT_fmt_table5, &(dev->fmt_table[5]), sizeof(double *));
      checkCudaErrors(err);
      if (err!=cudaSuccess) exit(1);
      err = cudaMemcpyToSymbol(FMT_fmt_table6, &(dev->fmt_table[6]), sizeof(double *));
      checkCudaErrors(err);
      if (err!=cudaSuccess) exit(1);
      err = cudaMemcpyToSymbol(FMT_fmt_table7, &(dev->fmt_table[7]), sizeof(double *));
      checkCudaErrors(err);
      if (err!=cudaSuccess) exit(1);
      err = cudaMemcpyToSymbol(FMT_fmt_table8, &(dev->fmt_table[8]), sizeof(double *));
      checkCudaErrors(err);
      if (err!=cudaSuccess) exit(1);
      err = cudaMemcpyToSymbol(FMT_fmt_table, &(dev->fmt_table[9]), sizeof(double *));
      checkCudaErrors(err);
      if (err!=cudaSuccess) exit(1);

      for (int m=0; m<=9; m++) {
        McudaMemcpyH2D(dev->fmt_table[m], dtmp[m], mtmp[m]*sizeof(double));
      }
    } else { // !(FMT_initialize<0)
    if(fp_prof) fprintf(fp_prof, "cuda_FMT_Init(): %d %d\n",FMT_Initialized, max_m1);
      if (dev->fmt_table[9]!=NULL) cudaFree(dev->fmt_table[9]);
      err = cudaMalloc((void **)&(dev->fmt_table[9]), mtmp[9]*sizeof(double));
      checkCudaErrors(err);
      if (err!=cudaSuccess) exit(1);
      err = cudaMemcpyToSymbol(FMT_fmt_table, &(dev->fmt_table[9]), sizeof(double *));
      checkCudaErrors(err);
      if (err!=cudaSuccess) exit(1);
      McudaMemcpyH2D(dev->fmt_table[9], dtmp[9], mtmp[9]*sizeof(double));
    }
    err = cudaMemcpyToSymbol(FMT_fmt_step_size, &step_size, sizeof(double));
    checkCudaErrors(err);
    if (err!=cudaSuccess) exit(1);
    err = cudaMemcpyToSymbol(FMT_fmt_inv_step_size, &inv_step_size, sizeof(double));
    checkCudaErrors(err);
    if (err!=cudaSuccess) exit(1);
    err = cudaMemcpyToSymbol(FMT_pi_div2, &pi_div2, sizeof(double));
    checkCudaErrors(err);
    if (err!=cudaSuccess) exit(1);
    err = cudaMemcpyToSymbol(FMT_fmt_max_m1, &max_m1, sizeof(int));
    checkCudaErrors(err);
    if (err!=cudaSuccess) exit(1);
  }
  FMT_Initialized = max_m1;

  return 0;
}

__host__ int cuda_FMT_Finalize(void)
{
  cudaError_t err;
  int i;
  int NDEV = cuda_get_numDevice();

  if (NDEV<=0) return 2;
  if (FMT_Initialized<0) return 1;
  if(fp_prof) fprintf(fp_prof, "cuda_FMT_Finalize()\n");
  for (i=0; i<NDEV; i++) {
    struct dev_Data *dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
    for (int m=0; m<9; m++) { 
      err = cudaFree(dev->fmt_table[m]);
      if (err!=cudaSuccess) exit(1);
    }
  }
  FMT_Initialized = -1;

  return 0;
}

/* ------------------------------------- */

/** 誤差関数を計算する関数
 *
 * \c fmt_initialize 関数で作成された誤差関数テーブルを用いるなどして
 * 誤差関数を計算する
 *
 * @attention
 * @li 事前に初期化ルーチンを \c fmt_initialize を呼び出しておく必要がある
 *
 * @param[out] f[] 計算した誤差関数を代入する配列
 *     （\f$ \tt{f[0]}\sim \tt{f[m]} \f$
 * @param[in] m 誤差関数の次数
 * @param[in] t 誤差関数の引数
 * @param[in] coef 誤差関数に掛ける定数
 *
 * @return なし
 *
 * @ingroup integ-fmt
 * */

__device__ void gpu_fmt(double f[], const int m, const double t, const double coef){
    int i;
    const double FMT_inv2 = 1.0 / 2.0;
//    int i,ts;
//    double d,t_inv,nu;
//    int ts0;

    // main
    if (t <= (2*m+36)) {
        const double FMT_inv6 = 1.0 / 6.0;
        double d;
        int ts,ts0;
	ts = 0.5 + t * FMT_fmt_inv_step_size;
	d = ts * FMT_fmt_step_size - t;
	ts0 = ts * FMT_fmt_max_m1;
#pragma unroll
	for (i=0; i<=m; i++){
	    f[i] = coef * (((LDG(FMT_fmt_table[ts0 + i + 3]) * FMT_inv6   * d
			   + LDG(FMT_fmt_table[ts0 + i + 2]) * FMT_inv2 ) * d
			   + LDG(FMT_fmt_table[ts0 + i + 1])            ) * d
			   + LDG(FMT_fmt_table[ts0 + i + 0])            );
	}
    } else {
        double t_inv;
        int nu;
	t_inv = FMT_inv2 / t;
	nu = 1;
	f[0] = coef * sqrt(FMT_pi_div2*t_inv);
#pragma unroll
	for (i=1; i<=m; i++){
	    f[i] = t_inv * f[i-1] * nu;
	    nu += 2;
	}
    }
}

__device__ void gpu_fmt0(double f[], const double t, const double coef){
    int i;
    const int m = 0;
    const double FMT_inv2 = 1.0 / 2.0;
//    int i,ts;
//    double d,t_inv,nu;
//    int ts0;

    // main
    if (t <= 36) {
        const double FMT_inv6 = 1.0 / 6.0;
        double d;
        int ts,ts0;
        int pos;
	ts = 0.5 + t * FMT_fmt_inv_step_size;
	d = ts * FMT_fmt_step_size - t;
        pos = ts * 4;
        f[0] = coef * ((( LDG(FMT_fmt_table0[pos + 3])  * d
                        + LDG(FMT_fmt_table0[pos + 2])) * d
                        + LDG(FMT_fmt_table0[pos + 1])) * d
                        + LDG(FMT_fmt_table0[pos + 0]));
    } else {
        double t_inv;
	t_inv = FMT_inv2 / t;
	f[0] = coef * sqrt(FMT_pi_div2*t_inv);
    }
}

__device__ void gpu_fmt1(double f[], const double t, const double coef){
    int i;
    const int m = 1;
    const double FMT_inv2 = 1.0 / 2.0;
//    int i,ts;
//    double d,t_inv,nu;
//    int ts0;

    // main
    if (t <= (2*m+36)) {
        const double FMT_inv6 = 1.0 / 6.0;
        double d;
        int ts,ts0;
        int pos;
	ts = 0.5 + t * FMT_fmt_inv_step_size;
	d = ts * FMT_fmt_step_size - t;
	pos = ts * (1+4);
#pragma unroll
	for (i=0; i<=m; i++){
	    f[i] = coef * (((LDG(FMT_fmt_table1[pos + i + 3]) * FMT_inv6   * d
			   + LDG(FMT_fmt_table1[pos + i + 2]) * FMT_inv2 ) * d
			   + LDG(FMT_fmt_table1[pos + i + 1])            ) * d
			   + LDG(FMT_fmt_table1[pos + i + 0])            );
	}
    } else {
        double t_inv;
	t_inv = FMT_inv2 / t;
	f[0] = coef * sqrt(FMT_pi_div2*t_inv);
        f[1] = t_inv * f[0];
    }
}

__device__ void gpu_fmt2(double f[], const double t, const double coef){
    int i;
    const int m = 2;
    const double FMT_inv2 = 1.0 / 2.0;
//    int i,ts;
//    double d,t_inv,nu;
//    int ts0;

    // main
    if (t <= (2*m+36)) {
        const double FMT_inv6 = 1.0 / 6.0;
        double d;
        int ts,ts0;
        int pos;
	ts = 0.5 + t * FMT_fmt_inv_step_size;
	d = ts * FMT_fmt_step_size - t;
	pos = ts * (2+4);
#pragma unroll
	for (i=0; i<=m; i++){
	    f[i] = coef * (((LDG(FMT_fmt_table2[pos + i + 3]) * FMT_inv6   * d
			   + LDG(FMT_fmt_table2[pos + i + 2]) * FMT_inv2 ) * d
			   + LDG(FMT_fmt_table2[pos + i + 1])            ) * d
			   + LDG(FMT_fmt_table2[pos + i + 0])            );
	}
    } else {
        double t_inv;
        int nu;
	t_inv = FMT_inv2 / t;
	nu = 1;
	f[0] = coef * sqrt(FMT_pi_div2*t_inv);
#pragma unroll
	for (i=1; i<=m; i++){
	    f[i] = t_inv * f[i-1] * nu;
	    nu += 2;
	}
    }
}

__device__ void gpu_fmt3(double f[], const double t, const double coef){
    int i;
    const int m = 3;
    const double FMT_inv2 = 1.0 / 2.0;
//    int i,ts;
//    double d,t_inv,nu;
//    int ts0;

    // main
    if (t <= (2*m+36)) {
        const double FMT_inv6 = 1.0 / 6.0;
        double d;
        int ts,ts0;
        int pos;
	ts = 0.5 + t * FMT_fmt_inv_step_size;
	d = ts * FMT_fmt_step_size - t;
	pos = ts * (3+4);
#pragma unroll
	for (i=0; i<=m; i++){
	    f[i] = coef * (((LDG(FMT_fmt_table3[pos + i + 3]) * FMT_inv6   * d
			   + LDG(FMT_fmt_table3[pos + i + 2]) * FMT_inv2 ) * d
			   + LDG(FMT_fmt_table3[pos + i + 1])            ) * d
			   + LDG(FMT_fmt_table3[pos + i + 0])            );
	}
    } else {
        double t_inv;
        int nu;
	t_inv = FMT_inv2 / t;
	nu = 1;
	f[0] = coef * sqrt(FMT_pi_div2*t_inv);
#pragma unroll
	for (i=1; i<=m; i++){
	    f[i] = t_inv * f[i-1] * nu;
	    nu += 2;
	}
    }
}

__device__ void gpu_fmt4(double f[], const double t, const double coef){
    int i;
    const int m = 4;
    const double FMT_inv2 = 1.0 / 2.0;
//    int i,ts;
//    double d,t_inv,nu;
//    int ts0;

    // main
    if (t <= (2*m+36)) {
        const double FMT_inv6 = 1.0 / 6.0;
        double d;
        int ts,ts0;
        int pos;
	ts = 0.5 + t * FMT_fmt_inv_step_size;
	d = ts * FMT_fmt_step_size - t;
	pos = ts * (4+4);
#pragma unroll
	for (i=0; i<=m; i++){
	    f[i] = coef * (((LDG(FMT_fmt_table4[pos + i + 3]) * FMT_inv6   * d
			   + LDG(FMT_fmt_table4[pos + i + 2]) * FMT_inv2 ) * d
			   + LDG(FMT_fmt_table4[pos + i + 1])            ) * d
			   + LDG(FMT_fmt_table4[pos + i + 0])            );
	}
    } else {
        double t_inv;
        int nu;
	t_inv = FMT_inv2 / t;
	nu = 1;
	f[0] = coef * sqrt(FMT_pi_div2*t_inv);
#pragma unroll
	for (i=1; i<=m; i++){
	    f[i] = t_inv * f[i-1] * nu;
	    nu += 2;
	}
    }
}

__device__ void gpu_fmt_p0(volatile double f[], const int m, const double t, const double coef,
    const int nidx, const int idx){

    const double FMT_inv2 = 1.0 / 2.0;
        const double FMT_inv6 = 1.0 / 6.0;
        double d;
        int ts,ts0;
	ts = 0.5 + t * FMT_fmt_inv_step_size;
	d = ts * FMT_fmt_step_size - t;
	ts0 = ts * FMT_fmt_max_m1;
#pragma unroll
	for (int i=idx; i<=m; i+=nidx){
	    f[i] = coef * (((LDG(FMT_fmt_table[ts0 + i + 3]) * FMT_inv6   * d
			   + LDG(FMT_fmt_table[ts0 + i + 2]) * FMT_inv2 ) * d
			   + LDG(FMT_fmt_table[ts0 + i + 1])            ) * d
			   + LDG(FMT_fmt_table[ts0 + i + 0])            );
	}
}

__device__ void gpu_fmt_p1(volatile double f[], const int m, const double t, const double coef,
    const int nidx, const int idx){

        const double FMT_inv2 = 1.0 / 2.0;
	double t_inv = FMT_inv2 / t;
	int nu = 1;

	f[0] = coef * sqrt(FMT_pi_div2*t_inv);
#pragma unroll
	for (int i=1; i<=m; i++){
	    f[i] = t_inv * f[i-1] * nu;
	    nu += 2;
	}

/*
        double ft;
	ft = coef * sqrt(FMT_pi_div2*t_inv);
#pragma unroll
	for (int i=1; i<=MIN2(idx,m); i++){
	    ft *= t_inv * nu;
	    nu += 2;
	}
        if (idx<=m) f[idx] = ft;
*/
        
        /*
        double t_inv;
        double t_inv2;
	t_inv = FMT_inv2 / t;
        t_inv2 = t_inv*t_inv;
	f[0] = coef * sqrt(FMT_pi_div2*t_inv);
	f[1] = f[0]*t_inv;
	f[2] = f[1]*t_inv*3;
	for (int i=(3+idx); i<=m; i+=3){
	    f[i] = t_inv2 * f[i-2] * (1+(i-1)*2) * (1+(i-2)*2) * (1+(i-3)*2);
	}
        */
}

__device__ void gpu_fmt_p(double f[], const int m, const double t, const double coef,
    const int nidx, const int idx){

    if (t <= (2*m+36)) {
      volatile double *s = f;
      gpu_fmt_p0(s, m, t, coef, nidx, idx);
    } else {
      volatile double *s = f;
      gpu_fmt_p1(s, m, t, coef, nidx, idx);
    }
}
   
