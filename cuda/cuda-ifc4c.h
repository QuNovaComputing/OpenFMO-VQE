#ifndef _CUDA_IFC4C_H_
#define _CUDA_IFC4C_H_

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------- */
int cuda_ifc4c_Init(const int maxlqn, const int max_num_klcs,
    const int nat_f, const int ncs_f, const int nao_f,
    const int ncspair_f, const int npspair_f,
    const int nat_m, const int ncs_m, const int nao_m,
    const int ncspair_m, const int npspair_m);
int cuda_ifc4c_Finalize(void);
int cuda_ifc4c_SetData_Symbol(const int iorj, const int maxlqn0,
    const int nat0, const int ncs0, const int nao0,
    const int ncspair0, const int npspair0);
int cuda_ifc4c_SetData(const int iorj, const int maxlqn,
    const int nat, const int ncs, const int nao,
    const int ncspair, const int npspair,
    const int *shel_atm, const int *shel_ini,
    const double *atom_x, const double *atom_y, const double *atom_z,
    const int *leading_cs_pair, const int *csp_leading_ps_pair,
    const int *csp_ics, const int *csp_jcs,
    const double *psp_zeta, const double *psp_dkps, const double *psp_xiza,
    const float *csp_schwarz, const double *D);
int cuda_ifc4c_calc_Init(void);
int cuda_ifc4c_SetDcs(const int ncs, const float *Dcs);
int cuda_ifc4c_GetVfrg(const int nao, double *V_frg);
int cuda_ifc4c_PsumVfrg(const int nao);
/* ------------------------------------- */
#ifdef __CUDACC__
//__global__ void gpu_ifc4c_ClearVfrg(void);
__global__ void gpu_ifc4c_ClearVfrg(const int nao, double *V);
__global__ void gpu_ifc4c_PsumVfrg(const int nao, double *V, double *VP);
#endif /* __CUDACC__ */

#ifdef __cplusplus
}
#endif

#endif /* _CUDA_IFC4C_H_ */
