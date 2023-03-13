#ifndef _CUDA_FMT_H_
#define _CUDA_FMT_H_

#ifdef __cplusplus
extern "C" {
#endif

int cuda_FMT_Init(double **dtmp, size_t *mtmp,
        const double step_size, const int max_m1);
int cuda_FMT_Finalize(void);

#ifdef __cplusplus
}
#endif

#ifdef __CUDACC__
__device__ void gpu_fmt(double f[], const int m, const double t, const double coef);
__device__ void gpu_fmt0(double f[], const double t, const double coef);
__device__ void gpu_fmt1(double f[], const double t, const double coef);
__device__ void gpu_fmt2(double f[], const double t, const double coef);
__device__ void gpu_fmt3(double f[], const double t, const double coef);
__device__ void gpu_fmt4(double f[], const double t, const double coef);
__device__ void gpu_fmt5(double f[], const double t, const double coef);
__device__ void gpu_fmt6(double f[], const double t, const double coef);
__device__ void gpu_fmt7(double f[], const double t, const double coef);
__device__ void gpu_fmt8(double f[], const double t, const double coef);
__device__ void gpu_fmt_p(double f[], const int m, const double t, const double coef,
    const int nidx, const int idx);
#endif /* __CUDACC__ */


#endif /* _CUDA_FMT_H_ */
