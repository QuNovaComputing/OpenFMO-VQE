#ifndef _CUDA_ROOT_H_
#define _CUDA_ROOT_H_

#ifdef __cplusplus
extern "C" {
#endif

__device__ void gpu_calc_root( const int nroot,
    const double T, double *U, double *W );

__device__ void gpu_root1( const double T, double *U, double *W );
__device__ void gpu_root2( const double T, double *U, double *W );
__device__ void gpu_root3( const double T, double *U, double *W );
__device__ void gpu_root4( const double T, double *U, double *W );
__device__ void gpu_root5( const double T, double *U, double *W );

#ifdef __cplusplus
}
#endif

#endif /* _CUDA_ROOT_H_ */
