#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#undef CUDA_IFC4C_BENCH

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif
#include "integ/ofmo-integ.h"
#include "integ/ofmo-ifc4c.h"
  extern FILE* fp_prof; // from common/ofmo-prof.h
#ifdef __cplusplus
}
#endif

#include <cuda.h>
#include "cuda-drv.h"
#include "cudalib.h"
#include "cuda-ifc4c.h"
#include "cuda-ifc4c-calc.h"
#include "cuda-fmt-m.h"

// Default {#block/#SMX, #threads/WARP_SIZE} set for each 2e-type
// This array will be converted to actual {#block, #threads} set
// in cuda_Init_Sub() by multiplying #SMX and WARP_SIZE, respectively.
int dim_ifc4c[][2] = {
#if CUDA_ARCH >= 350
 { 14,  7}, // ssss
 { 14,  8}, // ssps
 { 14,  8}, // sspp
 { 14,  6}, // ssds
 { 13,  8}, // ssdp
 { 14,  6}, // ssdd
 { 11,  8}, // psss
 { 11,  6}, // psps
 {  9,  8}, // pspp
 {  8,  7}, // psds
 {  9,  8}, // psdp
 {  0,  0}, // psdd
 {  6,  8}, // ppss
 {  4,  8}, // ppps
 {  3,  8}, // pppp
 {  3,  7}, // ppds
 {  0,  0}, // ppdp
 {  0,  0}, // ppdd
 {  5,  6}, // dsss
 {  6,  8}, // dsps
 {  5,  8}, // dspp
 {  5,  8}, // dsds
 {  5,  8}, // dsdp
 {  0,  0}, // dsdd
 {  4,  8}, // dpss
 {  3,  8}, // dpps
// {  3,  8}, // dppp
 {  0,  0}, // dppp
 {  3,  7}, // dpds
// {  3,  7}, // dpdp
 {  0,  0}, // dpdp
 {  0,  0}, // dpdd
 {  0,  0}, // ddss
 {  0,  0}, // ddps
 {  0,  0}, // ddpp
 {  0,  0}, // ddds
 {  0,  0}, // dddp
 {  0,  0}, // dddd
#else /* FERMI */
 { 12,  8}, // ssss
 { 12,  8}, // ssps
 { 12,  8}, // sspp
 { 12,  8}, // ssds
 { 12,  8}, // ssdp
 {  8,  6}, // ssdd
 {  6,  8}, // psss
 {  6,  6}, // psps
 {  6,  6}, // pspp
 {  6,  6}, // psds
 {  6,  6}, // psdp
 {  0,  0}, // psdd
 {  3,  7}, // ppss
 {  3,  9}, // ppps
 {  3,  6}, // pppp
// {  0,  0},
 {  3,  6}, // ppds
 {  0,  0}, // ppdp
 {  0,  0}, // ppdd
 {  3,  8}, // dsss
 {  3,  6}, // dsps
 {  7,  6}, // dspp
 {  7,  6}, // dsds
 {  9,  6}, // dsdp
// {  0,  0},
 {  0,  0}, // dsdd
 {  3,  6}, // dpss
 {  3,  6}, // dpps
 {  3,  6}, // dppp
// {  0,  0},
 {  3,  6}, // dpds
// {  0,  0},
 {  3,  7}, // dpdp
// {  0,  0},
 {  0,  0}, // dpdd
 {  0,  0}, // ddss
 {  0,  0}, // ddps
 {  0,  0}, // ddpp
 {  0,  0}, // ddds
 {  0,  0}, // dddp
 {  0,  0}, // dddd
#endif
};

int ifc4c_counter_ini_type[] = { // 0 for 0, 1 for nblk, 2 for nblk*NIJCSW
  1, 1, 1, 1, 1, 1, // ssxx
  1, 1, 1, 1, 1, 1, // psxx
  1, 1, 1, 1, 1, 1, // ppxx
  1, 1, 1, 1, 1, 1, // dsxx
  1, 1, 1, 1, 1, 1, // dpxx
  1, 1, 1, 1, 1, 1, // ddxx
};

/* ------------------------------------- */
/* ---- ssss ---- */
__host__ int cuda_ifc4c_os_ssss(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int Lab = La*(La+1)/2 + Lb;
  int Lcd = Lc*(Lc+1)/2 + Ld;
  int Labcd=Lab*6+Lcd;
  Labcd=0;

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns = 0;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);


//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim_ifc4c[Labcd][0];
  nthb = dim_ifc4c[Labcd][1];

  nwarps = nthb/WARP_SIZE;
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#if 0
#ifdef CUDA_FMT_M_SM
  Ns = cuda_FMT_m_get_size(0)*sizeof(double);
  Ns += nthb * 2 *  sizeof(double);
#else
  Ns = nthb * 3 * sizeof(double);
#endif
#ifdef DLB_KL_SSSS
  Ns += nwarps * sizeof(int);
#endif
#endif
  Ns += nthb *  sizeof(double); // sV
  Ns += nthb *  sizeof(double); // SSSS
#ifdef DLB_KL
  Ns += nwarps * sizeof(int); // klcsw
#endif

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    int c0 = nblk;
#endif
    gpu_ifc4c_os_ssss <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

/* ------------------------------------- */ /* ---- ssps ---- */
__host__ int cuda_ifc4c_os_ssps(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int Lab = La*(La+1)/2 + Lb;
  int Lcd = Lc*(Lc+1)/2 + Ld;
  int Labcd=Lab*6+Lcd;

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns = 0;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);


//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim_ifc4c[Labcd][0];
  nthb = dim_ifc4c[Labcd][1];

  nwarps = nthb/WARP_SIZE;
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns += cuda_FMT_m_get_size(1)*sizeof(double); // tbl
#endif
  Ns += nthb * sizeof(double); // sV
#ifdef DLB_KL
  Ns += nwarps * sizeof(int); // klcsw
#endif

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    int c0 = nblk;
#endif
    gpu_ifc4c_os_ssps <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

/* ------------------------------------- */ /* ---- sspp ---- */
__host__ int cuda_ifc4c_os_sspp(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int Lab = La*(La+1)/2 + Lb;
  int Lcd = Lc*(Lc+1)/2 + Ld;
  int Labcd=Lab*6+Lcd;

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns = 0;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);


//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim_ifc4c[Labcd][0];
  nthb = dim_ifc4c[Labcd][1];

  nwarps = nthb/WARP_SIZE;
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns += cuda_FMT_m_get_size(2)*sizeof(double); // tbl
#endif
  Ns += nthb * sizeof(double); // sV
#ifdef DLB_KL
  Ns += nwarps * sizeof(int); // klcsw
#endif

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    int c0 = nblk;
#endif
    gpu_ifc4c_os_sspp <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

/* ------------------------------------- */ /* ---- ssds ---- */
__host__ int cuda_ifc4c_os_ssds(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int Lab = La*(La+1)/2 + Lb;
  int Lcd = Lc*(Lc+1)/2 + Ld;
  int Labcd=Lab*6+Lcd;

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns = 0;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);


//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim_ifc4c[Labcd][0];
  nthb = dim_ifc4c[Labcd][1];

  nwarps = nthb/WARP_SIZE;
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns += cuda_FMT_m_get_size(2)*sizeof(double); // tbl
#endif
  Ns += nthb * sizeof(double); // sV
#ifdef DLB_KL
  Ns += nwarps * sizeof(int); // klcsw
#endif

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    int c0 = nblk;
#endif
    gpu_ifc4c_os_ssds <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

/* ------------------------------------- */ /* ---- ssdp ---- */
__host__ int cuda_ifc4c_os_ssdp(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int Lab = La*(La+1)/2 + Lb;
  int Lcd = Lc*(Lc+1)/2 + Ld;
  int Labcd=Lab*6+Lcd;

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns = 0;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);


//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim_ifc4c[Labcd][0];
  nthb = dim_ifc4c[Labcd][1];

  nwarps = nthb/WARP_SIZE;
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns += cuda_FMT_m_get_size(3)*sizeof(double); // tbl
#endif
  Ns += nthb * sizeof(double); // sV
#ifdef DLB_KL
  Ns += nwarps * sizeof(int); // klcsw
#endif

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    int c0 = nblk;
#endif
    gpu_ifc4c_os_ssdp <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

/* ------------------------------------- */ /* ---- ssdd ---- */
__host__ int cuda_ifc4c_os_ssdd(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int Lab = La*(La+1)/2 + Lb;
  int Lcd = Lc*(Lc+1)/2 + Ld;
  int Labcd=Lab*6+Lcd;

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns = 0;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);


//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim_ifc4c[Labcd][0];
  nthb = dim_ifc4c[Labcd][1];

  nwarps = nthb/WARP_SIZE;
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns += cuda_FMT_m_get_size(4)*sizeof(double); // tbl
#endif
  Ns += nthb * sizeof(double); // sV
#ifdef DLB_KL
  Ns += nwarps * sizeof(int); // klcsw
#endif

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    int c0 = nblk;
#endif
    gpu_ifc4c_os_ssdd <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

/* ------------------------------------- */ /* ---- psss ---- */
__host__ int cuda_ifc4c_os_psss(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int Lab = La*(La+1)/2 + Lb;
  int Lcd = Lc*(Lc+1)/2 + Ld;
  int Labcd=Lab*6+Lcd;
  int Nab = NNAO(La)*NNAO(Lb);

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns = 0;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);


//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim_ifc4c[Labcd][0];
  nthb = dim_ifc4c[Labcd][1];

  nwarps = nthb/WARP_SIZE;
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns += cuda_FMT_m_get_size(1)*sizeof(double); // tbl
#endif
  Ns += nthb * sizeof(double); // sV
#ifdef DLB_KL
  Ns += nwarps * sizeof(int); // klcsw
#endif
  Ns += nwarps * Nab * sizeof(double); // sVw

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    int c0 = nblk;
#endif
    gpu_ifc4c_os_psss <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

/* ------------------------------------- */ /* ---- psps ---- */
__host__ int cuda_ifc4c_os_psps(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int Lab = La*(La+1)/2 + Lb;
  int Lcd = Lc*(Lc+1)/2 + Ld;
  int Labcd=Lab*6+Lcd;
  int Nab = NNAO(La)*NNAO(Lb);

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns = 0;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);


//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim_ifc4c[Labcd][0];
  nthb = dim_ifc4c[Labcd][1];

  nwarps = nthb/WARP_SIZE;
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns += cuda_FMT_m_get_size(2)*sizeof(double); // tbl
#endif
  Ns += nthb * sizeof(double); // sV
#ifdef DLB_KL
  Ns += nwarps * sizeof(int); // klcsw
#endif
  Ns += nwarps * Nab * sizeof(double); // sVw

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    int c0 = nblk;
#endif
    gpu_ifc4c_os_psps <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

/* ------------------------------------- */ /* ---- pspp ---- */
__host__ int cuda_ifc4c_os_pspp(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int Lab = La*(La+1)/2 + Lb;
  int Lcd = Lc*(Lc+1)/2 + Ld;
  int Labcd=Lab*6+Lcd;
  int Nab = NNAO(La)*NNAO(Lb);

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns = 0;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);


//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim_ifc4c[Labcd][0];
  nthb = dim_ifc4c[Labcd][1];

  nwarps = nthb/WARP_SIZE;
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns += cuda_FMT_m_get_size(3)*sizeof(double); // tbl
#endif
  Ns += nthb * sizeof(double); // sV
#ifdef DLB_KL
  Ns += nwarps * sizeof(int); // klcsw
#endif
  Ns += nwarps * Nab * sizeof(double); // sVw

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    int c0 = nblk;
#endif
    gpu_ifc4c_os_pspp <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

/* ------------------------------------- */ /* ---- psds ---- */
__host__ int cuda_ifc4c_os_psds(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int Lab = La*(La+1)/2 + Lb;
  int Lcd = Lc*(Lc+1)/2 + Ld;
  int Labcd=Lab*6+Lcd;
  int Nab = NNAO(La)*NNAO(Lb);

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns = 0;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);


//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim_ifc4c[Labcd][0];
  nthb = dim_ifc4c[Labcd][1];

  nwarps = nthb/WARP_SIZE;
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns += cuda_FMT_m_get_size(3)*sizeof(double); // tbl
#endif
  Ns += nthb * sizeof(double); // sV
#ifdef DLB_KL
  Ns += nwarps * sizeof(int); // klcsw
#endif
  Ns += nwarps * Nab * sizeof(double); // sVw

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    int c0 = nblk;
#endif
    gpu_ifc4c_os_psds <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

/* ------------------------------------- */ /* ---- psdp ---- */
__host__ int cuda_ifc4c_os_psdp(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int Lab = La*(La+1)/2 + Lb;
  int Lcd = Lc*(Lc+1)/2 + Ld;
  int Labcd=Lab*6+Lcd;
  int Nab = NNAO(La)*NNAO(Lb);

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns = 0;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);


//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim_ifc4c[Labcd][0];
  nthb = dim_ifc4c[Labcd][1];

  nwarps = nthb/WARP_SIZE;
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns += cuda_FMT_m_get_size(4)*sizeof(double); // tbl
#endif
  Ns += nthb * sizeof(double); // sV
#ifdef DLB_KL
  Ns += nwarps * sizeof(int); // klcsw
#endif
  Ns += nwarps * Nab * sizeof(double); // sVw

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    int c0 = nblk;
#endif
    gpu_ifc4c_os_psdp <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

/* ------------------------------------- */ /* ---- ppss ---- */
__host__ int cuda_ifc4c_os_ppss(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int Lab = La*(La+1)/2 + Lb;
  int Lcd = Lc*(Lc+1)/2 + Ld;
  int Labcd=Lab*6+Lcd;
  int Nab = NNAO(La)*NNAO(Lb);

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns = 0;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);


//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim_ifc4c[Labcd][0];
  nthb = dim_ifc4c[Labcd][1];

  nwarps = nthb/WARP_SIZE;
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns += cuda_FMT_m_get_size(2)*sizeof(double); // tbl
#endif
  Ns += nthb * sizeof(double); // sV
#ifdef DLB_KL
  Ns += nwarps * sizeof(int); // klcsw
#endif
  Ns += Nab * sizeof(double); // sVw

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    int c0 = nblk;
#endif
    gpu_ifc4c_os_ppss <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

/* ------------------------------------- */ /* ---- ppps ---- */
__host__ int cuda_ifc4c_os_ppps(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int Lab = La*(La+1)/2 + Lb;
  int Lcd = Lc*(Lc+1)/2 + Ld;
  int Labcd=Lab*6+Lcd;
  int Nab = NNAO(La)*NNAO(Lb);

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns = 0;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);


//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim_ifc4c[Labcd][0];
  nthb = dim_ifc4c[Labcd][1];

  nwarps = nthb/WARP_SIZE;
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns += cuda_FMT_m_get_size(3)*sizeof(double); // tbl
#endif
  Ns += nthb * sizeof(double); // sV
#ifdef DLB_KL
  Ns += nwarps * sizeof(int); // klcsw
#endif
  Ns += Nab * sizeof(double); // sVw

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    int c0 = nblk;
#endif
    gpu_ifc4c_os_ppps <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

/* ------------------------------------- */ /* ---- pppp ---- */
__host__ int cuda_ifc4c_os_pppp(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int Lab = La*(La+1)/2 + Lb;
  int Lcd = Lc*(Lc+1)/2 + Ld;
  int Labcd=Lab*6+Lcd;
  int Nab = NNAO(La)*NNAO(Lb);

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns = 0;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);


//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim_ifc4c[Labcd][0];
  nthb = dim_ifc4c[Labcd][1];

  nwarps = nthb/WARP_SIZE;
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns += cuda_FMT_m_get_size(4)*sizeof(double); // tbl
#endif
  Ns += nthb * sizeof(double); // sV
#ifdef DLB_KL
  Ns += nwarps * sizeof(int); // klcsw
#endif
  Ns += Nab * sizeof(double); // sVw

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    int c0 = nblk;
#endif
    gpu_ifc4c_os_pppp <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

/* ------------------------------------- */ /* ---- ppds ---- */
__host__ int cuda_ifc4c_os_ppds(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int Lab = La*(La+1)/2 + Lb;
  int Lcd = Lc*(Lc+1)/2 + Ld;
  int Labcd=Lab*6+Lcd;
  int Nab = NNAO(La)*NNAO(Lb);

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns = 0;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);


//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim_ifc4c[Labcd][0];
  nthb = dim_ifc4c[Labcd][1];

  nwarps = nthb/WARP_SIZE;
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns += cuda_FMT_m_get_size(4)*sizeof(double); // tbl
#endif
  Ns += nthb * sizeof(double); // sV
#ifdef DLB_KL
  Ns += nwarps * sizeof(int); // klcsw
#endif
  Ns += Nab * sizeof(double); // sVw

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    int c0 = nblk;
#endif
    gpu_ifc4c_os_ppds <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

/* ------------------------------------- */ /* ---- dsss ---- */
__host__ int cuda_ifc4c_os_dsss(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int Lab = La*(La+1)/2 + Lb;
  int Lcd = Lc*(Lc+1)/2 + Ld;
  int Labcd=Lab*6+Lcd;
  int Nab = NNAO(La)*NNAO(Lb);

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns = 0;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);


//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim_ifc4c[Labcd][0];
  nthb = dim_ifc4c[Labcd][1];

  nwarps = nthb/WARP_SIZE;
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns += cuda_FMT_m_get_size(2)*sizeof(double); // tbl
#endif
  Ns += nthb * sizeof(double); // sV
#ifdef DLB_KL
  Ns += nwarps * sizeof(int); // klcsw
#endif
  Ns += nwarps * Nab * sizeof(double); // sVw

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    int c0 = nblk;
#endif
    gpu_ifc4c_os_dsss <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

/* ------------------------------------- */ /* ---- dsps ---- */
__host__ int cuda_ifc4c_os_dsps(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int Lab = La*(La+1)/2 + Lb;
  int Lcd = Lc*(Lc+1)/2 + Ld;
  int Labcd=Lab*6+Lcd;
  int Nab = NNAO(La)*NNAO(Lb);

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns = 0;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);


//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim_ifc4c[Labcd][0];
  nthb = dim_ifc4c[Labcd][1];

  nwarps = nthb/WARP_SIZE;
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns += cuda_FMT_m_get_size(3)*sizeof(double); // tbl
#endif
  Ns += nthb * sizeof(double); // sV
#ifdef DLB_KL
  Ns += nwarps * sizeof(int); // klcsw
#endif
  //Ns += nwarps * Nab * sizeof(double); // sVw
  Ns += Nab * sizeof(double); // sVw

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    int c0 = nblk;
#endif
    gpu_ifc4c_os_dsps <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

/* ------------------------------------- */ /* ---- dspp ---- */
__host__ int cuda_ifc4c_os_dspp(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int Lab = La*(La+1)/2 + Lb;
  int Lcd = Lc*(Lc+1)/2 + Ld;
  int Labcd=Lab*6+Lcd;
  int Nab = NNAO(La)*NNAO(Lb);

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns = 0;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);


//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim_ifc4c[Labcd][0];
  nthb = dim_ifc4c[Labcd][1];

  nwarps = nthb/WARP_SIZE;
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns += cuda_FMT_m_get_size(4)*sizeof(double); // tbl
#endif
  Ns += nthb * sizeof(double); // sV
#ifdef DLB_KL
  Ns += nwarps * sizeof(int); // klcsw
#endif
  //Ns += nwarps * Nab * sizeof(double); // sVw
  Ns += Nab * sizeof(double); // sVw

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    int c0 = nblk;
#endif
    gpu_ifc4c_os_dspp <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

/* ------------------------------------- */ /* ---- dsds ---- */
__host__ int cuda_ifc4c_os_dsds(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int Lab = La*(La+1)/2 + Lb;
  int Lcd = Lc*(Lc+1)/2 + Ld;
  int Labcd=Lab*6+Lcd;
  int Nab = NNAO(La)*NNAO(Lb);

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns = 0;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);


//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim_ifc4c[Labcd][0];
  nthb = dim_ifc4c[Labcd][1];

  nwarps = nthb/WARP_SIZE;
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns += cuda_FMT_m_get_size(4)*sizeof(double); // tbl
#endif
  Ns += nthb * sizeof(double); // sV
#ifdef DLB_KL
  Ns += nwarps * sizeof(int); // klcsw
#endif
  //Ns += nwarps * Nab * sizeof(double); // sVw
  Ns += Nab * sizeof(double); // sVw

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    int c0 = nblk;
#endif
    gpu_ifc4c_os_dsds <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

/* ------------------------------------- */ /* ---- dsdp ---- */
__host__ int cuda_ifc4c_os_dsdp(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int Lab = La*(La+1)/2 + Lb;
  int Lcd = Lc*(Lc+1)/2 + Ld;
  int Labcd=Lab*6+Lcd;
  int Nab = NNAO(La)*NNAO(Lb);

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns = 0;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);


//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim_ifc4c[Labcd][0];
  nthb = dim_ifc4c[Labcd][1];

  nwarps = nthb/WARP_SIZE;
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns += cuda_FMT_m_get_size(5)*sizeof(double); // tbl
#endif
  Ns += nthb * sizeof(double); // sV
#ifdef DLB_KL
  Ns += nwarps * sizeof(int); // klcsw
#endif
  //Ns += nwarps * Nab * sizeof(double); // sVw
  Ns += Nab * sizeof(double); // sVw

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    int c0 = nblk;
#endif
    gpu_ifc4c_os_dsdp <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

/* ------------------------------------- */ /* ---- dpss ---- */
__host__ int cuda_ifc4c_os_dpss(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int Lab = La*(La+1)/2 + Lb;
  int Lcd = Lc*(Lc+1)/2 + Ld;
  int Labcd=Lab*6+Lcd;
  int Nab = NNAO(La)*NNAO(Lb);

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns = 0;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);


//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim_ifc4c[Labcd][0];
  nthb = dim_ifc4c[Labcd][1];

  nwarps = nthb/WARP_SIZE;
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns += cuda_FMT_m_get_size(3)*sizeof(double); // tbl
#endif
  Ns += nthb * sizeof(double); // sV
#ifdef DLB_KL
  Ns += nwarps * sizeof(int); // klcsw
#endif
  //Ns += nwarps * Nab * sizeof(double); // sVw
  Ns += Nab * sizeof(double); // sVw

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    int c0 = nblk;
#endif
    gpu_ifc4c_os_dpss <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

/* ------------------------------------- */ /* ---- dpps ---- */
__host__ int cuda_ifc4c_os_dpps(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int Lab = La*(La+1)/2 + Lb;
  int Lcd = Lc*(Lc+1)/2 + Ld;
  int Labcd=Lab*6+Lcd;
  int Nab = NNAO(La)*NNAO(Lb);

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns = 0;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);


//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim_ifc4c[Labcd][0];
  nthb = dim_ifc4c[Labcd][1];

  nwarps = nthb/WARP_SIZE;
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns += cuda_FMT_m_get_size(4)*sizeof(double); // tbl
#endif
  Ns += nthb * sizeof(double); // sV
#ifdef DLB_KL
  Ns += nwarps * sizeof(int); // klcsw
#endif
  //Ns += nwarps * Nab * sizeof(double); // sVw
  Ns += Nab * sizeof(double); // sVw

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    int c0 = nblk;
#endif
    gpu_ifc4c_os_dpps <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

/* ------------------------------------- */ /* ---- dppp ---- */
__host__ int cuda_ifc4c_os_dppp(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int Lab = La*(La+1)/2 + Lb;
  int Lcd = Lc*(Lc+1)/2 + Ld;
  int Labcd=Lab*6+Lcd;
  int Nab = NNAO(La)*NNAO(Lb);

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns = 0;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);


//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim_ifc4c[Labcd][0];
  nthb = dim_ifc4c[Labcd][1];

  nwarps = nthb/WARP_SIZE;
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns += cuda_FMT_m_get_size(5)*sizeof(double); // tbl
#endif
  Ns += nthb * sizeof(double); // sV
#ifdef DLB_KL
  Ns += nwarps * sizeof(int); // klcsw
#endif
  //Ns += nwarps * Nab * sizeof(double); // sVw
  Ns += Nab * sizeof(double); // sVw

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    int c0 = nblk;
#endif
    gpu_ifc4c_os_dppp <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

/* ------------------------------------- */ /* ---- dpds ---- */
__host__ int cuda_ifc4c_os_dpds(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int Lab = La*(La+1)/2 + Lb;
  int Lcd = Lc*(Lc+1)/2 + Ld;
  int Labcd=Lab*6+Lcd;
  int Nab = NNAO(La)*NNAO(Lb);

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns = 0;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);


//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim_ifc4c[Labcd][0];
  nthb = dim_ifc4c[Labcd][1];

  nwarps = nthb/WARP_SIZE;
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns += cuda_FMT_m_get_size(5)*sizeof(double); // tbl
#endif
  Ns += nthb * sizeof(double); // sV
#ifdef DLB_KL
  Ns += nwarps * sizeof(int); // klcsw
#endif
  //Ns += nwarps * Nab * sizeof(double); // sVw
  Ns += Nab * sizeof(double); // sVw

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    int c0 = nblk;
#endif
    gpu_ifc4c_os_dpds <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

/* ------------------------------------- */ /* ---- dpdp ---- */
__host__ int cuda_ifc4c_os_dpdp(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int Lab = La*(La+1)/2 + Lb;
  int Lcd = Lc*(Lc+1)/2 + Ld;
  int Labcd=Lab*6+Lcd;
  int Nab = NNAO(La)*NNAO(Lb);

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns = 0;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);


//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim_ifc4c[Labcd][0];
  nthb = dim_ifc4c[Labcd][1];

  nwarps = nthb/WARP_SIZE;
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns += cuda_FMT_m_get_size(6)*sizeof(double); // tbl
#endif
  Ns += nthb * sizeof(double); // sV
#ifdef DLB_KL
  Ns += nwarps * sizeof(int); // klcsw
#endif
  //Ns += nwarps * Nab * sizeof(double); // sVw
  Ns += Nab * sizeof(double); // sVw

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    int c0 = nblk;
#endif
    gpu_ifc4c_os_dpdp <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}
/* ------------------------------------- */


static int (*cuda_host_ifc4c_calc_a[])(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) = {
    // original
    /*ofmo_ifc4c_ssss__, ofmo_ifc4c_ssps__, ofmo_ifc4c_sspp__,
    ofmo_ifc4c_ssds__, ofmo_ifc4c_ssdp__, ofmo_ifc4c_ssdd__,
    ofmo_ifc4c_psss__, ofmo_ifc4c_psps__, ofmo_ifc4c_pspp__,
    ofmo_ifc4c_psds__, ofmo_ifc4c_psdp__, ofmo_ifc4c_psdd__,
    ofmo_ifc4c_ppss__, ofmo_ifc4c_ppps__, ofmo_ifc4c_pppp__,
    ofmo_ifc4c_ppds__, ofmo_ifc4c_ppdp__, ofmo_ifc4c_ppdd__,
    ofmo_ifc4c_dsss__, ofmo_ifc4c_dsps__, ofmo_ifc4c_dspp__,
    ofmo_ifc4c_dsds__, ofmo_ifc4c_dsdp__, ofmo_ifc4c_dsdd__,
    ofmo_ifc4c_dpss__, ofmo_ifc4c_dpps__, ofmo_ifc4c_dppp__,
    ofmo_ifc4c_dpds__, ofmo_ifc4c_dpdp__, ofmo_ifc4c_dpdd__,
    ofmo_ifc4c_ddss__, ofmo_ifc4c_ddps__, ofmo_ifc4c_ddpp__,
    ofmo_ifc4c_ddds__, ofmo_ifc4c_dddp__, ofmo_ifc4c_dddd__,*/
    // OS
    ofmo_ifc4c_os_ssss, ofmo_ifc4c_os_ssps, ofmo_ifc4c_os_sspp,
    ofmo_ifc4c_os_ssds, ofmo_ifc4c_os_ssdp, ofmo_ifc4c_os_ssdd,
    ofmo_ifc4c_os_psss, ofmo_ifc4c_os_psps, ofmo_ifc4c_os_pspp,
    ofmo_ifc4c_os_psds, ofmo_ifc4c_os_psdp, ofmo_ifc4c_os_psdd,
    ofmo_ifc4c_os_ppss, ofmo_ifc4c_os_ppps, ofmo_ifc4c_os_pppp,
    ofmo_ifc4c_os_ppds, ofmo_ifc4c_os_ppdp, ofmo_ifc4c_os_ppdd,
    ofmo_ifc4c_os_dsss, ofmo_ifc4c_os_dsps, ofmo_ifc4c_os_dspp,
    ofmo_ifc4c_os_dsds, ofmo_ifc4c_os_dsdp, ofmo_ifc4c_os_dsdd,
    ofmo_ifc4c_os_dpss, ofmo_ifc4c_os_dpps, ofmo_ifc4c_os_dppp,
    ofmo_ifc4c_os_dpds, ofmo_ifc4c_os_dpdp, ofmo_ifc4c_os_dpdd,
    ofmo_ifc4c_os_ddss, ofmo_ifc4c_os_ddps, ofmo_ifc4c_os_ddpp,
    ofmo_ifc4c_os_ddds, ofmo_ifc4c_os_dddp, ofmo_ifc4c_os_dddd,
    // Rys
    /*ofmo_ifc4c_rys_ssss, ofmo_ifc4c_rys_ssps, ofmo_ifc4c_rys_sspp,
    ofmo_ifc4c_rys_ssds, ofmo_ifc4c_rys_ssdp, ofmo_ifc4c_rys_ssdd,
    ofmo_ifc4c_rys_psss, ofmo_ifc4c_rys_psps, ofmo_ifc4c_rys_pspp,
    ofmo_ifc4c_rys_psds, ofmo_ifc4c_rys_psdp, ofmo_ifc4c_rys_psdd,
    ofmo_ifc4c_rys_ppss, ofmo_ifc4c_rys_ppps, ofmo_ifc4c_rys_pppp,
    ofmo_ifc4c_rys_ppds, ofmo_ifc4c_rys_ppdp, ofmo_ifc4c_rys_ppdd,
    ofmo_ifc4c_rys_dsss, ofmo_ifc4c_rys_dsps, ofmo_ifc4c_rys_dspp,
    ofmo_ifc4c_rys_dsds, ofmo_ifc4c_rys_dsdp, ofmo_ifc4c_rys_dsdd,
    ofmo_ifc4c_rys_dpss, ofmo_ifc4c_rys_dpps, ofmo_ifc4c_rys_dppp,
    ofmo_ifc4c_rys_dpds, ofmo_ifc4c_rys_dpdp, ofmo_ifc4c_rys_dpdd,
    ofmo_ifc4c_rys_ddss, ofmo_ifc4c_rys_ddps, ofmo_ifc4c_rys_ddpp,
    ofmo_ifc4c_rys_ddds, ofmo_ifc4c_rys_dddp, ofmo_ifc4c_rys_dddd,*/
};

static int (*cuda_ifc4c_calc_a[])(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) = {
    // OS
    cuda_ifc4c_os_ssss, cuda_ifc4c_os_ssps, cuda_ifc4c_os_sspp,
    cuda_ifc4c_os_ssds, cuda_ifc4c_os_ssdp, cuda_ifc4c_os_ssdd,
    cuda_ifc4c_os_psss, cuda_ifc4c_os_psps, cuda_ifc4c_os_pspp,
    cuda_ifc4c_os_psds, cuda_ifc4c_os_psdp, ofmo_ifc4c_os_psdd,
    cuda_ifc4c_os_ppss, cuda_ifc4c_os_ppps, cuda_ifc4c_os_pppp,
    cuda_ifc4c_os_ppds, ofmo_ifc4c_os_ppdp, ofmo_ifc4c_os_ppdd,
    cuda_ifc4c_os_dsss, cuda_ifc4c_os_dsps, cuda_ifc4c_os_dspp,
    cuda_ifc4c_os_dsds, cuda_ifc4c_os_dsdp, ofmo_ifc4c_os_dsdd,
    cuda_ifc4c_os_dpss, cuda_ifc4c_os_dpps, cuda_ifc4c_os_dppp,
    cuda_ifc4c_os_dpds, cuda_ifc4c_os_dpdp, ofmo_ifc4c_os_dpdd,
    ofmo_ifc4c_os_ddss, ofmo_ifc4c_os_ddps, ofmo_ifc4c_os_ddpp,
    ofmo_ifc4c_os_ddds, ofmo_ifc4c_os_dddp, ofmo_ifc4c_os_dddd,
};

static char *sifc4c[] = {
#if 0
  "ssss", "ssps", "sspp", "ssds", "ssdp", "ssdd",
  "psss", "psps", "pspp", "psds", "psdp", "psdd",
  "ppss", "ppps", "pppp", "ppds", "ppdp", "ppdd",
  "dsss", "dsps", "dspp", "dsds", "dsdp", "dsdd",
  "dpss", "dpps", "dppp", "dpds", "dpdp", "dpdd",
  "ddss", "ddps", "ddpp", "ddds", "dddp", "dddd",
#else
  "(ss,ss)", "(ss,ps)", "(ss,pp)", "(ss,ds)", "(ss,dp)", "(ss,dd)",
  "(ps,ss)", "(ps,ps)", "(ps,pp)", "(ps,ds)", "(ps,dp)", "(ps,dd)",
  "(pp,ss)", "(pp,ps)", "(pp,pp)", "(pp,ds)", "(pp,dp)", "(pp,dd)",
  "(ds,ss)", "(ds,ps)", "(ds,pp)", "(ds,ds)", "(ds,dp)", "(ds,dd)",
  "(dp,ss)", "(dp,ps)", "(dp,pp)", "(dp,ds)", "(dp,dp)", "(dp,dd)",
  "(dd,ss)", "(dd,ps)", "(dd,pp)", "(dd,ds)", "(dd,dp)", "(dd,dd)",
#endif
};
static double wifc4c[] = {
  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,
  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,
  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,
  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,
  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,
  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,
};

void cuda_print_wifc4c(void) {
  int i;
  int n = 6*6;

  int NDEV = cuda_get_numDevice();
//  if (NDEV<=0) return;
  if (CUDA_ME!=0) return;
#ifdef CUDA_IFC4C_BENCH
#ifdef _OPENMP
#pragma omp master
#endif
  {
    printf("--- wifc4c ---\n");
    for (i=0; i<n; i++) {
      if (wifc4c[i]>=0) {
        int nb=0, nt=0;
        if (NDEV>0&&dim_ifc4c[i][0]!=0) {
          nb = dim_ifc4c[i][0];
          nt = dim_ifc4c[i][1];
        }
//        printf("%4s %8.4f (%3d,%3d)\n",sifc4c[i],wifc4c[i],nb,nt);
        printf("%7s %8.4f (%3d,%3d)\n",sifc4c[i],wifc4c[i],nb,nt);
      }
      wifc4c[i]=-1.0e0;
    }
    printf("-----------\n");
    fflush(stdout);
  }
#endif /* CUDA_IFC4C_BENCH */
}

__host__ int cuda_ifc4c_calc( const int idev,
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] )
{
  int NDEV = cuda_get_numDevice();
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int Lab = La*(La+1)/2 + Lb;
  int Lcd = Lc*(Lc+1)/2 + Ld;
  int Labcd = Lab*6 + Lcd;

#if 0
#pragma omp master
  {
    int ijcs0 = leading_cs_pair_frg[Lab];
    int ijcs1 = leading_cs_pair_frg[Lab+1];
    int klcs0 = leading_cs_pair_mon[Lcd];
    int klcs1 = leading_cs_pair_mon[Lcd+1];
    int nijcs = ijcs1-ijcs0+1;
    int nklcs = klcs1-klcs0+1;
    if (fp_prof) fprintf(fp_prof, "%s: %d-%d\n", sifc4c[Labcd], nijcs, nklcs);
  }
#endif

#ifdef CUDA_IFC4C_BENCH
#pragma omp master
  checkCudaErrors(cudaDeviceSynchronize());
#pragma omp barrier
  double w0,w1;
#pragma omp master
  w0 = cuda_Wtime();
#endif /* CUDA_IFC4C_BENCH */
  if (NDEV==0||dim_ifc4c[Labcd][0]==0) {
    if (idev==0)
    cuda_host_ifc4c_calc_a[Labcd](
        pnworkers, pworkerid,
        pLa, pLb, pLc, pLd,
        shel_atm_frg, shel_ini_frg,
        atom_x_frg, atom_y_frg, atom_z_frg,
        leading_cs_pair_frg,
        csp_schwarz_frg, csp_ics_frg, csp_jcs_frg,
        csp_leading_ps_pair_frg,
        psp_zeta_frg, psp_dkps_frg, psp_xiza_frg,
        shel_atm_mon, shel_ini_mon,
        atom_x_mon, atom_y_mon, atom_z_mon,
        leading_cs_pair_mon,
        csp_schwarz_mon, csp_ics_mon, csp_jcs_mon,
        csp_leading_ps_pair_mon,
        psp_zeta_mon, psp_dkps_mon, psp_xiza_mon,
        D_mon, V_frg );
  } else {
    if (idev==1)
    cuda_ifc4c_calc_a[Labcd](
        pnworkers, pworkerid,
        pLa, pLb, pLc, pLd,
        shel_atm_frg, shel_ini_frg,
        atom_x_frg, atom_y_frg, atom_z_frg,
        leading_cs_pair_frg,
        csp_schwarz_frg, csp_ics_frg, csp_jcs_frg,
        csp_leading_ps_pair_frg,
        psp_zeta_frg, psp_dkps_frg, psp_xiza_frg,
        shel_atm_mon, shel_ini_mon,
        atom_x_mon, atom_y_mon, atom_z_mon,
        leading_cs_pair_mon,
        csp_schwarz_mon, csp_ics_mon, csp_jcs_mon,
        csp_leading_ps_pair_mon,
        psp_zeta_mon, psp_dkps_mon, psp_xiza_mon,
        D_mon, V_frg );
#ifndef CUDA_IFC4C_BENCH
  }
#else
    checkCudaErrors(cudaDeviceSynchronize());
  }
#pragma omp barrier
#pragma omp master
  {
  w1 = cuda_Wtime();
  if (wifc4c[Labcd]<0) wifc4c[Labcd] = 0.0;
  wifc4c[Labcd] += w1-w0;
  }
#endif /* CUDA_IFC4C_BENCH */

  return 0;
}

