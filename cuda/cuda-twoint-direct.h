/**
 * @file ofmo-twoint-direct.h
 * ２電子積分計算、および、２電子ハミルトン行列計算を
 * 行う関数群のプロトタイプ定義をしているヘッダファイル
 * */
#ifndef _CUDA_TWOINT_DIRECT_H_
#define _CUDA_TWOINT_DIRECT_H_

#define NIJCSW 2


__global__ void gpu_twoint_direct_counter_init(int counter);

__global__ void gpu_twoint_direct_ssss_( const int nwks, const int iwk,
    const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_twoint_direct_psss_( const int nwks, const int iwk,
    const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_twoint_direct_psps_( const int nwks, const int iwk,
    const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_twoint_direct_ppss_( const int nwks, const int iwk,
    const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_twoint_direct_ppps_( const int nwks, const int iwk,
    const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_twoint_direct_pppp_( const int nwks, const int iwk,
    const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_twoint_direct_dsss_( const int nwks, const int iwk,
    const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_twoint_direct_dsps_( const int nwks, const int iwk,
    const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_twoint_direct_dspp_( const int nwks, const int iwk,
    const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_twoint_direct_dsds_( const int nwks, const int iwk,
    const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_twoint_direct_dpss_( const int nwks, const int iwk,
    const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_twoint_direct_dpps_( const int nwks, const int iwk,
    const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_twoint_direct_dppp_( const int nwks, const int iwk,
    const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_twoint_direct_dpds_( const int nwks, const int iwk,
    const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_twoint_direct_dpdp_( const int nwks, const int iwk,
    const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_twoint_direct_ddss_( const int nwks, const int iwk,
    const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_twoint_direct_ddps_( const int nwks, const int iwk,
    const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_twoint_direct_ddpp_( const int nwks, const int iwk,
    const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_twoint_direct_ddds_( const int nwks, const int iwk,
    const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_twoint_direct_dddp_( const int nwks, const int iwk,
    const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_twoint_direct_dddd_( const int nwks, const int iwk,
    const float eps_eri, const float eps_ps4, const float eps_sch );

#define CUDA_EBUF_FULL   1
#define CUDA_EBUF_NOFULL 0

#endif /* _CUDA_TWOINT_DIRECT_H_ */
