#include "cuda-utils.h"

__device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
      assumed = old;
      old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val + __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}


__device__ double atomicSub(double* address, double val)
{
  /*
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
      assumed = old;
      old = atomicCAS(address_as_ull, assumed,
                      __double_as_longlong(__longlong_as_double(assumed) - val));
    } while (assumed != old);
    return __longlong_as_double(old);
  */
  return atomicAdd(address, -val);
}

/* --------------------------------------------------- */
#if CUDA_ARCH >= 350
__device__   __inline__    double  shfl_xor ( double value,  int const lane )
{
  return  __hiloint2double( __shfl_xor(__double2hiint(value),lane),
                            __shfl_xor(__double2loint(value),lane));
}

__device__   __inline__    double  shfl_down ( double value,  int const delta )
{
  return  __hiloint2double( __shfl_down(__double2hiint(value),delta),
                            __shfl_down(__double2loint(value),delta));
}

__device__ double warpReduceG(double v, const int tid)
{
  double t;
  v += shfl_down(v, 16);
  v += shfl_down(v, 8);
  v += shfl_down(v, 4);
  v += shfl_down(v, 2);
  v += shfl_down(v, 1);
  return v;
}

__device__ double grp9ReduceG(double v, const int tid)
{
  v += shfl_down(v, 6) + shfl_down(v, 3);
  v += shfl_down(v, 2) + shfl_down(v, 1);
  return v;
}

__device__ double grp6ReduceG(double v, const int tid)
{
  v += shfl_down(v, 3);
  v += shfl_down(v, 2) + shfl_down(v, 1);
  return v;
}
#endif

__device__ void warpReduce(volatile double *sdata, const int tid)
{
#if CUDA_ARCH >= 350
  double v = sdata[tid];
  sdata[tid] = warpReduceG(v, tid);
#else
  if (tid<WARP_SIZE/2) {
//  sdata[tid] += sdata[tid+32];
  sdata[tid] += sdata[tid+16];
  sdata[tid] += sdata[tid+8];
  sdata[tid] += sdata[tid+4];
  sdata[tid] += sdata[tid+2];
  sdata[tid] += sdata[tid+1];
  }
#endif
}

__device__ void grp9Reduce(volatile double *sdata, const int tid)
{
//#ifdef FERMI
  if (tid<3)  sdata[tid] += sdata[tid+3] + sdata[tid+6];
  if (tid==0) sdata[tid] += sdata[tid+1] + sdata[tid+2];
/*
#else
  double v = sdata[tid];
  sdata[tid] = grp9ReduceG(v, tid);
#endif
*/
}

__device__ void grp6Reduce(volatile double *sdata, const int tid)
{
//#ifdef FERMI
  if (tid<3)  sdata[tid] += sdata[tid+3];
  if (tid==0) sdata[0] += sdata[1] + sdata[2];
/*
#else
  double v = sdata[tid];
  sdata[tid] = grp6ReduceG(v, tid);
#endif
*/
}

