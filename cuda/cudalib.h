#ifndef _CUDALIB_H_
#define _CUDALIB_H_

#ifdef __cplusplus
extern "C" {
#endif

#define WARP_SIZE 32

#if defined(SORT_IJ_SCHWARZ) || defined(SORT_KL_SCHWARZ)
#define SORT_INDEX_SCHWARZ
#endif

#if defined(DLB_KL)
#define DLB_KL_SSSS
#define DLB_KL_PSSS
#define DLB_KL_PSPS
#define DLB_KL_PPSS
#define DLB_KL_PPPS
#define DLB_KL_PPPP
#define DLB_KL_DSSS
#define DLB_KL_DSPS
#define DLB_KL_DSPP
#define DLB_KL_DSDS
#define DLB_KL_DPSS
#define DLB_KL_DPPS
#define DLB_KL_DPPP
#define DLB_KL_DPDS
#define DLB_KL_DPDP
#define DLB_KL_DDSS
#define DLB_KL_DDPS
#define DLB_KL_DDPP
#define DLB_KL_DDDS
#define DLB_KL_DDDP
#define DLB_KL_DDDD
#endif

/* ------------------------------------- */

#define WORK_SIZE 30


struct dev_Data {
  int numSM;
  double *fmt_table[10];
  double *fmt_m_table[9]; // (!) mmax + 1
  int fmt_m_mmax;
//  int *p_ijcounter;
  int *ijcounter;
  int *shel_atm;
  int *shel_ini;
  double *atom_x;
  double *atom_y;
  double *atom_z;
  int *leading_cs;
  int *leading_cs_pair;
  int *csp_leading_ps_pair;
  int *csp_ics;
  int *csp_jcs;
  int *sklcs_b;
  double *psp_zeta;
  double *psp_dkps;
  double *psp_xiza;
  double *DFACT;
  float *csp_schwarz;
  float *Dcs;
  int *sorted_csp;
  double *Ds;
  double *G_d;
  double *G_b;
  double *Gkl_b;
  double *Gi_t;
//  double *Gj_t;
  double *work_b;
// for ifc4c
  int *shel_atm_frg;
  int *shel_ini_frg;
  double *atom_x_frg;
  double *atom_y_frg;
  double *atom_z_frg;
  int *leading_cs_frg;
  int *leading_cs_pair_frg;
  int *csp_leading_ps_pair_frg;
  int *csp_ics_frg;
  int *csp_jcs_frg;
  double *psp_zeta_frg;
  double *psp_dkps_frg;
  double *psp_xiza_frg;
  float *csp_schwarz_frg;
  int *shel_atm_mon;
  int *shel_ini_mon;
  double *atom_x_mon;
  double *atom_y_mon;
  double *atom_z_mon;
  int *leading_cs_mon;
  int *leading_cs_pair_mon;
  int *csp_leading_ps_pair_mon;
  int *csp_ics_mon;
  int *csp_jcs_mon;
  double *psp_zeta_mon;
  double *psp_dkps_mon;
  double *psp_xiza_mon;
  float *csp_schwarz_mon;
  double *D_mon;
  double *V_frg;
  double *V_frgP;
};

/* ------------------------------------- */

#ifdef DUMMY
#ifdef __CUDACC__
extern __device__ __constant__ int ncs;
extern __device__ __constant__ int nat;
extern __device__ __constant__ int maxlqn;
extern __device__ __constant__ int ncspair;
extern __device__ __constant__ int npspair;
extern __device__ __constant__ int nao;
extern __device__ __constant__ int max_num_klcs;
extern __device__ __constant__ int num_Gkl;
extern __device__ __constant__ int num_Gi;

extern __device__ __constant__ int *p_ijcounter;
extern __device__ __constant__ int *shel_atm;
extern __device__ __constant__ int *shel_ini;
extern __device__ __constant__ double *atom_x;
extern __device__ __constant__ double *atom_y;
extern __device__ __constant__ double *atom_z;
extern __device__ __constant__ int *leading_cs;
extern __device__ __constant__ int *leading_cs_pair;
extern __device__ __constant__ int *csp_leading_ps_pair;
extern __device__ __constant__ int *csp_ics;
extern __device__ __constant__ int *csp_jcs;
extern __device__ __constant__ int *sklcs_b;
extern __device__ __constant__ double *psp_zeta;
extern __device__ __constant__ double *psp_dkps;
extern __device__ __constant__ double *psp_xiza;
extern __device__ __constant__ double *dfact;
extern __device__ __constant__ float *csp_schwarz;
extern __device__ __constant__ float *Dcs;
extern __device__ __constant__ int *sorted_csp;
extern __device__ __constant__ double *Ds;
extern __device__ __constant__ double *G_d;
extern __device__ __constant__ double *G_b;
extern __device__ __constant__ double *Gkl_b;
extern __device__ __constant__ double *work_b;
#endif /* __CUDACC__ */
#endif /* DUMMY */

/* ------------------------------------- */
#define NNAO(i) (((i)+1)*((i)+2)/2)

#define McudaMemcpyH2D(dest,src,size) \
      checkCudaErrors(cudaMemcpy((dest),(src),(size), cudaMemcpyHostToDevice))
#define McudaMemcpyD2H(dest,src,size) \
      checkCudaErrors(cudaMemcpy((dest),(src),(size), cudaMemcpyDeviceToHost))

#if CUDA_ARCH >= 350
#define LDG(a) __ldg(&(a))
#ifdef CUDA_FMT_M_SM
#define LDM(a) (a)
#else
#define LDM(a) __ldg(&(a))
#endif
#else
#define LDG(a) (a)
#define LDM(a) (a)
#endif

#ifndef TRUE
#define TRUE  (1==1)
#endif
#ifndef FALSE
#define FALSE (1==0)
#endif

#ifndef MIN2
#define MIN2(a, b) (((a)<(b))? (a): (b))
#endif
#ifndef MAX2
#define MAX2(a, b) (((a)>(b))? (a): (b))
#endif

/* ------------------------------------- */
struct dev_Data *cuda_SCF_get_dev_Data(const int ndev);
int cuda_set_Device(int idev);
int cuda_get_numDevice(void); 
int cuda_use_Device(void); // boolean
int cuda_get_numBlocks(void);
int cuda_get_numThreads(void);
int cuda_get_Nprocs(void);
int cuda_get_myRank(void);
/* ------------------------------------- */
int cuda_find2eType(char *str);
int cuda_set_optCPU(char *optstr);
int cuda_get_optCPU(int Labcd);
int cuda_set_optsync(int optsync);
int cuda_get_optsync(void);
void cuda_Print_DEFS(void);
int cuda_Init_Sub(const int ndev, const int myrank, const int nprocs);
int cuda_Finalize_Sub(void);
int cuda_SCF_Init_Symbol(
    const int ncs0, const int nat0, const int maxlqn0,
    const int ncspair0, const int npspair0, const int nao0,
    const int max_num_klcs0, const int num_Gkl0, const int num_Gi0,
    const struct dev_Data *dev);
int cuda_SCF_Malloc(const int ncs, const int nat, const int maxlqn,
        const int ncspair, const int npspair, const int nao,
        const int max_num_klcs, const int num_Gkl, const int num_Gi);
int cuda_SCF_Free(void);
int cuda_SCF_Finalize(void);
int cuda_SCF_Init(const int ncs, const int nat, const int maxlqn,
    const int ncspair, const int npspair, const int nao,
    const int *shel_atm, const int *shel_ini,
    const double *atom_x, const double *atom_y, const double *atom_z,
    const int *leading_cs,
    const int *leading_cs_pair, const int *csp_leading_ps_pair,
    const int *csp_ics, const int *csp_jcs,
    const double *psp_zeta, const double *psp_dkps, const double *psp_xiza,
    const float *csp_schwarz,
    const int nblk, const int nthb
    );
/* ------------------------------------- */
int cuda_max_num_klcs(const int maxlqn, const int *leading_cs_pair);
int cuda_num_Gkl(const int maxlqn, const int *leading_cs_pair);
int cuda_num_Gi(const int maxlqn, const int *leading_cs_pair,
        const int *shel_ini);
/* ------------------------------------- */
int cuda_sort_csp_schwarz(int *sorted_csp, const int ncspair,
        const int maxlqn, const int *leading_cs_pair, const float *csp_schwarz);
/* ------------------------------------- */
int cuda_genGmat_Init(const int ncs, const int nao,
        const float *Dcs, const double *D, const double x_coef);
int cuda_genGmat_Add(const int nao, double *G);
/* ------------------------------------- */

#ifdef __cplusplus
}
#endif

#endif /* _CUDALIB_H_ */
