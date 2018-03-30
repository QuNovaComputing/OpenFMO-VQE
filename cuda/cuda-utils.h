#ifndef _CUDA_UTILS_H_
#define _CUDA_UTILS_H_

#ifdef __CUDACC__

__device__ double atomicAdd(double* address, double val);
__device__ double atomicSub(double* address, double val);
__device__ void warpReduce(volatile double *sdata, const int tid);
__device__ void grp9Reduce(volatile double *sdata, const int tid);
__device__ void grp6Reduce(volatile double *sdata, const int tid);
#if CUDA_ARCH >= 350
__device__ double warpReduceG(double v, const int tid);
__device__ double grp9ReduceG(double v, const int tid);
__device__ double grp6ReduceG(double v, const int tid);
#endif

#endif /* __CUDACC__ */
#endif /* _CUDA_UTILS_H_ */
