#ifndef _CUDA_IFC4C_OS_H_
#define _CUDA_IFC4C_OS_H_

__global__ void gpu_ifc4c_os_ssss( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_ssps( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_sspp( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_ssds( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_ssdp( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_ssdd( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_psss( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_psps( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_pspp( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_psds( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_psdp( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_psdd( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_ppss( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_ppps( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_pppp( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_ppds( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_ppdp( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_ppdd( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_dsss( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_dsps( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_dspp( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_dsds( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_dsdp( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_dsdd( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_dpss( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_dpps( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_dppp( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_dpds( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_dpdp( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_dpdd( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_ddss( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_ddps( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_ddpp( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_ddds( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_dddp( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );
__global__ void gpu_ifc4c_os_dddd( const int nwks, const int iwk,
        const float eps_eri, const float eps_ps4, const float eps_sch );

#endif /* _CUDA_IFC4C_OS_H_ */
