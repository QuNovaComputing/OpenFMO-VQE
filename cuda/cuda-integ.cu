#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif
#include "integ/ofmo-twoint.h"
#include "integ/ofmo-twoint-direct.h"
#include "integ/ofmo-os-xxxx.h"
#include "integ/ofmo-rys-xxxx.h"
  extern FILE* fp_prof; // from common/ofmo-prof.h
#ifdef __cplusplus
}
#endif

#include <cuda.h>
#include "cuda-drv.h"
#include "cudalib.h"
#include "cuda-twoint-direct.h"
#include "cuda-integ.h"
#include "cuda-fmt-m.h"

// Default {#block/#SMX, #threads/WARP_SIZE} set for each 2e-type
// This array will be converted to actual {#block, #threads} set
// in cuda_Init_Sub() by multiplying #SMX and WARP_SIZE, respectively.
int dim2e[][2] = {
#if CUDA_ARCH >= 350
  {  9,  2}, // ssss
  {  8,  2}, // psss
  {  8,  2}, // psps
  {  8,  2}, // ppss
  {  8,  2}, // ppps
  {  8,  1}, // pppp
//  { 10,  1}, // dsss
  {  9,  1}, // dsss
  {  8,  1}, // dsps
  {  8,  1}, // dspp
  {  0,  0}, // dsds
  {  0,  0}, // dpss
  {  8,  1}, // dpps
  {  0,  0}, // dppp
  {  0,  0}, // dpds
  {  0,  0}, // dpdp
  {  0,  0}, // ddss
  {  0,  0}, // ddps
  {  0,  0}, // ddpp
  {  0,  0}, // ddds
  {  0,  0}, // dddp
  {  0,  0}, // dddd
#else /* FERMI */
  {  5,  2}, // ssss
  {  4,  2}, // psss
  {  7,  1}, // psps
  {  6,  1}, // ppss
  {  5,  1}, // ppps
  {  5,  1}, // pppp
  {  5,  1}, // dsss
  {  3,  1}, // dsps
  {  5,  1}, // dspp
  {  0,  0}, // dsds
  {  0,  0}, // dpss
  {  4,  1}, // dpps
  {  0,  0}, // dppp
  {  0,  0}, // dpds
  {  0,  0}, // dpdp
  {  0,  0}, // ddss
  {  0,  0}, // ddps
  {  0,  0}, // ddpp
  {  0,  0}, // ddds
  {  0,  0}, // dddp
  {  0,  0}, // dddd
#endif
  {  0,  0}  // buff
};

int counter_ini_type[] = { // 0 for 0, 1 for nblk, 2 for nblk*NIJCSW
  1,                // sxxx
  2, 1, 2, 2, 1,    // pxxx             
  2, 2, 2, 0,       // dsxx
  0, 2, 0, 0, 0,    // dpxx
  0, 0, 0, 0, 0, 0  // ddxx
};

/* ------------------------------------- */
__host__ int cuda_twoint_direct_ssss_(
        // paralleization
        const int *pnworkers, const int *pworkerid,
        // integral type data
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        // basis set & cutoff table data
        const int shel_atm[], const int shel_ini[],
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // concerned about buffered direct method
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int nao = *pnao;

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
//  int nwks = np * NDEV * NBLOCKS;
//  int iwk0 = me * NDEV * NBLOCKS;
//  nwkblk = nworkers * ndev * gridDim.x;
//  iwkblk = workerid * ndev * gridDim.x + idev * gridDim.x + blockIdx.x;
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns;
  int Labcd = 0;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);

//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = 96; // best
  nthb = 64;
  nblk = 64; // better for nblk <= 64
  nthb = 96;
  nblk = dim2e[Labcd][0];
  nthb = dim2e[Labcd][1];

  nwarps = nthb/WARP_SIZE;
//  if (nthb%WARP_SIZE!=0) exit(1);
//  if (nwarps>NTHREADS) exit(1);
//  if (nblk>NBLOCKS) exit(1);
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns = cuda_FMT_m_get_size(0)*sizeof(double);
  Ns += nthb * 2 *  sizeof(double);
#else
  Ns = nthb * 3 * sizeof(double);
#endif
#ifdef DLB_KL_SSSS
  Ns += nwarps * sizeof(int);
#endif

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
//    gpu_twoint_direct_counter_init <<< NBLOCKS, NTHREADS >>> (NBLOCKS);
    int c0 = nblk;
    //gpu_twoint_direct_counter_init <<< dimGrid, dimBlock >>> (Labcd, c0);
    //cudaMemcpyH2D(dev->ijcounter + Labcd, &c0, sizeof(int));
#endif
//    gpu_twoint_direct_ssss_ <<< NBLOCKS, NTHREADS >>> (nwks, iwk,
    gpu_twoint_direct_ssss_ <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

__host__ int cuda_twoint_direct_psss_(
        // paralleization
        const int *pnworkers, const int *pworkerid,
        // integral type data
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        // basis set & cutoff table data
        const int shel_atm[], const int shel_ini[],
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // concerned about buffered direct method
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int nao = *pnao;

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
//  int nwks = np * NDEV * NBLOCKS;
//  int iwk0 = me * NDEV * NBLOCKS;
//  nwkblk = nworkers * ndev * gridDim.x;
//  iwkblk = workerid * ndev * gridDim.x + idev * gridDim.x + blockIdx.x;
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns;
  int Labcd = 1;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);

//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = 128; // best
  nthb = 32;
  nblk = 64;  // better for nblk <= 64
  nthb = 64;
  nblk = dim2e[Labcd][0];
  nthb = dim2e[Labcd][1];

  nwarps = nthb/WARP_SIZE;
//  if (NTHREADS%WARP_SIZE!=0) exit(1);
//  if (nwarps>NTHREADS) exit(1);
//  if (nblk>NBLOCKS) exit(1);
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
  Ns = 0;
#ifdef CUDA_FMT_M_SM
  Ns += cuda_FMT_m_get_size(1)*sizeof(double);
  Ns += nthb * 3 * sizeof(double);
#else
  Ns += nthb * 3 * sizeof(double);
#endif
#ifdef DLB_KL_PSSS
  Ns += nwarps * sizeof(int);
#endif

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
//    gpu_twoint_direct_counter_init <<< NBLOCKS, NTHREADS >>> (NBLOCKS);
//    gpu_twoint_direct_counter_init <<< dimGrid, dimBlock >>> (Labcd, nblk*NIJCSW);
    int c0 = nblk*NIJCSW;
    //gpu_twoint_direct_counter_init <<< dimGrid, dimBlock >>> (Labcd, c0);
    //cudaMemcpyH2D(dev->ijcounter + Labcd, &c0, sizeof(int));
#endif
    gpu_twoint_direct_psss_ <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

__host__ int cuda_twoint_direct_psps_(
        // paralleization
        const int *pnworkers, const int *pworkerid,
        // integral type data
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        // basis set & cutoff table data
        const int shel_atm[], const int shel_ini[],
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // concerned about buffered direct method
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int nao = *pnao;

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
//  int nwks = np * NDEV * NBLOCKS;
//  int iwk0 = me * NDEV * NBLOCKS;
//  nwkblk = nworkers * ndev * gridDim.x;
//  iwkblk = workerid * ndev * gridDim.x + idev * gridDim.x + blockIdx.x;
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns;
  int Labcd = 2;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);

//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = 96; // best
  nthb = 32;
  nblk = 48; // better for nblk <= 64
  nthb = 64;
  nblk = dim2e[Labcd][0];
  nthb = dim2e[Labcd][1];

  nwarps = nthb/WARP_SIZE;
//  if (NTHREADS%WARP_SIZE!=0) exit(1);
//  if (nwarps>NTHREADS) exit(1);
//  if (nblk>NBLOCKS) exit(1);
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns = cuda_FMT_m_get_size(2)*sizeof(double);
  Ns += nthb * sizeof(double);
#else
  Ns = nthb * 3 * sizeof(double);
#endif
#ifdef DLB_KL_PSPS
  Ns += nwarps * sizeof(int);
#endif


  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    //gpu_twoint_direct_counter_init <<< dimGrid, dimBlock >>> (Labcd, nblk);
    int c0 = nblk;
    //gpu_twoint_direct_counter_init <<< dimGrid, dimBlock >>> (Labcd, c0);
    //cudaMemcpyH2D(dev->ijcounter + Labcd, &c0, sizeof(int));
#endif
    gpu_twoint_direct_psps_ <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

__host__ int cuda_twoint_direct_ppss_(
        // paralleization
        const int *pnworkers, const int *pworkerid,
        // integral type data
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        // basis set & cutoff table data
        const int shel_atm[], const int shel_ini[],
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // concerned about buffered direct method
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int nao = *pnao;
  int last_ijcs = *plast_ijcs, last_klcs = *plast_klcs;

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
//  int nwks = np * NDEV * NBLOCKS;
//  int iwk0 = me * NDEV * NBLOCKS;
//  nwkblk = nworkers * ndev * gridDim.x;
//  iwkblk = workerid * ndev * gridDim.x + idev * gridDim.x + blockIdx.x;
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns;
  int Labcd = 3;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);

//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim2e[Labcd][0];
  nthb = dim2e[Labcd][1];
  if (last_ijcs!=-1||last_klcs!=-1) exit(1);

  nwarps = nthb/WARP_SIZE;
//  if (NTHREADS%WARP_SIZE!=0) exit(1);
//  if (nwarps>NTHREADS) exit(1);
//  if (nblk>NBLOCKS) exit(1);
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns = cuda_FMT_m_get_size(2)*sizeof(double);
  Ns += nthb * sizeof(double);
#else
  Ns = nthb * 3 * sizeof(double);
#endif
#ifdef DLB_KL_PPSS
  Ns += nwarps * sizeof(int);
#endif

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
//    gpu_twoint_direct_counter_init <<< NBLOCKS, NTHREADS >>> (NBLOCKS);
//    gpu_twoint_direct_counter_init <<< dimGrid, dimBlock >>> (Labcd, nblk*NIJCSW);
    int c0 = nblk*NIJCSW;
    //gpu_twoint_direct_counter_init <<< dimGrid, dimBlock >>> (Labcd, c0);
    //cudaMemcpyH2D(dev->ijcounter + Labcd, &c0, sizeof(int));
#endif
    gpu_twoint_direct_ppss_ <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

__host__ int cuda_twoint_direct_ppps_(
        // paralleization
        const int *pnworkers, const int *pworkerid,
        // integral type data
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        // basis set & cutoff table data
        const int shel_atm[], const int shel_ini[],
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // concerned about buffered direct method
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int nao = *pnao;
  int last_ijcs = *plast_ijcs, last_klcs = *plast_klcs;

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
//  int nwks = np * NDEV * NBLOCKS;
//  int iwk0 = me * NDEV * NBLOCKS;
//  nwkblk = nworkers * ndev * gridDim.x;
//  iwkblk = workerid * ndev * gridDim.x + idev * gridDim.x + blockIdx.x;
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns;
  int Labcd = 4;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);

//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim2e[Labcd][0];
  nthb = dim2e[Labcd][1];
  if (last_ijcs!=-1||last_klcs!=-1) exit(1);

  nwarps = nthb/WARP_SIZE;
//  if (NTHREADS%WARP_SIZE!=0) exit(1);
//  if (nwarps>NTHREADS) exit(1);
//  if (nblk>NBLOCKS) exit(1);
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns = cuda_FMT_m_get_size(3)*sizeof(double);
  Ns += nthb * sizeof(double);
#else
  Ns = nthb * 3 * sizeof(double);
#endif
//  Ns = nwarps*(WARP_SIZE/9) * (3*3+4) * sizeof(double);
#ifdef DLB_KL_PPPS
  Ns += nwarps * sizeof(int);
#endif
//  cudaFuncSetCacheConfig(gpu_twoint_direct_ppps_, cudaFuncCachePreferShared);

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
//    gpu_twoint_direct_counter_init <<< NBLOCKS, NTHREADS >>> (NBLOCKS);
//    gpu_twoint_direct_counter_init <<< dimGrid, dimBlock >>> (Labcd, nblk*NIJCSW);
    int c0 = nblk*NIJCSW;
    //gpu_twoint_direct_counter_init <<< dimGrid, dimBlock >>> (Labcd, c0);
    //cudaMemcpyH2D(dev->ijcounter + Labcd, &c0, sizeof(int));
#endif
    gpu_twoint_direct_ppps_ <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

__host__ int cuda_twoint_direct_pppp_(
        // paralleization
        const int *pnworkers, const int *pworkerid,
        // integral type data
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        // basis set & cutoff table data
        const int shel_atm[], const int shel_ini[],
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // concerned about buffered direct method
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int nao = *pnao;
  int last_ijcs = *plast_ijcs, last_klcs = *plast_klcs;

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
//  int nwks = np * NDEV * NBLOCKS;
//  int iwk0 = me * NDEV * NBLOCKS;
//  nwkblk = nworkers * ndev * gridDim.x;
//  iwkblk = workerid * ndev * gridDim.x + idev * gridDim.x + blockIdx.x;
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns;
  int Labcd = 5;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);

//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim2e[Labcd][0];
  nthb = dim2e[Labcd][1];

  nwarps = nthb/WARP_SIZE;
//  if (NTHREADS%WARP_SIZE!=0) exit(1);
//  if (nwarps>NTHREADS) exit(1);
//  if (nblk>NBLOCKS) exit(1);
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns = cuda_FMT_m_get_size(4)*sizeof(double);
  Ns += nthb * 1 * sizeof(double);
#else
  Ns = nthb * 3 * sizeof(double);
#endif
#ifdef DLB_KL_PPPP
  Ns += nwarps * sizeof(int);
#endif


  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    //gpu_twoint_direct_counter_init <<< dimGrid, dimBlock >>> (Labcd, nblk);
    int c0 = nblk;
    //gpu_twoint_direct_counter_init <<< dimGrid, dimBlock >>> (Labcd, c0);
    //cudaMemcpyH2D(dev->ijcounter + Labcd, &c0, sizeof(int));
#endif
    gpu_twoint_direct_pppp_ <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

__host__ int cuda_twoint_direct_dsss_(
        // paralleization
        const int *pnworkers, const int *pworkerid,
        // integral type data
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        // basis set & cutoff table data
        const int shel_atm[], const int shel_ini[],
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // concerned about buffered direct method
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int nao = *pnao;
  int last_ijcs = *plast_ijcs, last_klcs = *plast_klcs;

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
//  int nwks = np * NDEV * NBLOCKS;
//  int iwk0 = me * NDEV * NBLOCKS;
//  nwkblk = nworkers * ndev * gridDim.x;
//  iwkblk = workerid * ndev * gridDim.x + idev * gridDim.x + blockIdx.x;
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns;
  int Labcd = 6;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);

//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim2e[Labcd][0];
  nthb = dim2e[Labcd][1];
  if (last_ijcs!=-1||last_klcs!=-1) exit(1);

  nwarps = nthb/WARP_SIZE;
//  if (NTHREADS%WARP_SIZE!=0) exit(1);
//  if (nwarps>NTHREADS) exit(1);
//  if (nblk>NBLOCKS) exit(1);
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns = cuda_FMT_m_get_size(2)*sizeof(double);
  Ns += nthb * sizeof(double);
#else
  Ns = nthb * 3 * sizeof(double);
#endif
#ifdef DLB_KL_DSSS
  Ns += nwarps * sizeof(int);
#endif

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    //gpu_twoint_direct_counter_init <<< dimGrid, dimBlock >>> (Labcd, nblk*NIJCSW);
    int c0 = nblk*NIJCSW;
    //gpu_twoint_direct_counter_init <<< dimGrid, dimBlock >>> (Labcd, c0);
    //cudaMemcpyH2D(dev->ijcounter + Labcd, &c0, sizeof(int));
#endif
    gpu_twoint_direct_dsss_ <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

__host__ int cuda_twoint_direct_dsps_(
        // paralleization
        const int *pnworkers, const int *pworkerid,
        // integral type data
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        // basis set & cutoff table data
        const int shel_atm[], const int shel_ini[],
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // concerned about buffered direct method
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int nao = *pnao;
  int last_ijcs = *plast_ijcs, last_klcs = *plast_klcs;

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
//  int nwks = np * NDEV * NBLOCKS;
//  int iwk0 = me * NDEV * NBLOCKS;
//  nwkblk = nworkers * ndev * gridDim.x;
//  iwkblk = workerid * ndev * gridDim.x + idev * gridDim.x + blockIdx.x;
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns;
  int Labcd = 7;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);

//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim2e[Labcd][0];
  nthb = dim2e[Labcd][1];
  if (last_ijcs!=-1||last_klcs!=-1) exit(1);

  nwarps = nthb/WARP_SIZE;
//  if (NTHREADS%WARP_SIZE!=0) exit(1);
//  if (nwarps>NTHREADS) exit(1);
//  if (nblk>NBLOCKS) exit(1);
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns = cuda_FMT_m_get_size(3)*sizeof(double);
  Ns += nthb * sizeof(double);
#else
  Ns = nthb * 3 * sizeof(double);
#endif
#ifdef DLB_KL_DSPS
  Ns += nwarps * sizeof(int);
#endif

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    //gpu_twoint_direct_counter_init <<< dimGrid, dimBlock >>> (Labcd, nblk*NIJCSW);
    int c0 = nblk*NIJCSW;
    //gpu_twoint_direct_counter_init <<< dimGrid, dimBlock >>> (Labcd, c0);
    //cudaMemcpyH2D(dev->ijcounter + Labcd, &c0, sizeof(int));
#endif
    gpu_twoint_direct_dsps_ <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}

__host__ int cuda_twoint_direct_dspp_(
        // paralleization
        const int *pnworkers, const int *pworkerid,
        // integral type data
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        // basis set & cutoff table data
        const int shel_atm[], const int shel_ini[],
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // concerned about buffered direct method
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int nao = *pnao;
  int last_ijcs = *plast_ijcs, last_klcs = *plast_klcs;

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
//  int nwks = np * NDEV * NBLOCKS;
//  int iwk0 = me * NDEV * NBLOCKS;
//  nwkblk = nworkers * ndev * gridDim.x;
//  iwkblk = workerid * ndev * gridDim.x + idev * gridDim.x + blockIdx.x;
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns;
  int Labcd = 8;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);

//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim2e[Labcd][0];
  nthb = dim2e[Labcd][1];
  if (last_ijcs!=-1||last_klcs!=-1) exit(1);

  nwarps = nthb/WARP_SIZE;
//  if (NTHREADS%WARP_SIZE!=0) exit(1);
//  if (nwarps>NTHREADS) exit(1);
//  if (nblk>NBLOCKS) exit(1);
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns = cuda_FMT_m_get_size(4)*sizeof(double);
  Ns += nthb * 1 * sizeof(double);
#else
  Ns = nthb * 3 * sizeof(double);
#endif
#ifdef DLB_KL_DSPP
  Ns += nwarps * sizeof(int);
#endif

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    //gpu_twoint_direct_counter_init <<< dimGrid, dimBlock >>> (Labcd, nblk*NIJCSW);
    int c0 = nblk*NIJCSW;
    //gpu_twoint_direct_counter_init <<< dimGrid, dimBlock >>> (Labcd, c0);
    //cudaMemcpyH2D(dev->ijcounter + Labcd, &c0, sizeof(int));
#endif
    gpu_twoint_direct_dspp_ <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}


__host__ int cuda_twoint_direct_dpps_(
        // paralleization
        const int *pnworkers, const int *pworkerid,
        // integral type data
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        // basis set & cutoff table data
        const int shel_atm[], const int shel_ini[],
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // concerned about buffered direct method
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] )
{
  int nworkers = *pnworkers;
  int workerid = *pworkerid;
  int La = *pLa, Lb = *pLb, Lc = *pLc, Ld = *pLd;
  int nao = *pnao;
  int last_ijcs = *plast_ijcs, last_klcs = *plast_klcs;

  cudaError_t err;
  int ret = 0;
  int i;
  int NDEV = cuda_get_numDevice();
  struct dev_Data *dev;
  int np = cuda_get_Nprocs();
  int me = cuda_get_myRank();
//  int nwks = np * NDEV * NBLOCKS;
//  int iwk0 = me * NDEV * NBLOCKS;
//  nwkblk = nworkers * ndev * gridDim.x;
//  iwkblk = workerid * ndev * gridDim.x + idev * gridDim.x + blockIdx.x;
  int nwks = np * NDEV;
  int iwk0 = me * NDEV;
  int nwarps;
  int nblk,nthb;
  size_t Ns;
  int Labcd = 11;
  float eps_eri = ofmo_twoint_eps_eri(0);
  float eps_ps4 = ofmo_twoint_eps_ps4(0);
  float eps_sch = ofmo_twoint_eps_sch(0);

//  nblk = NBLOCKS;
//  nwarps = NTHREADS/WARP_SIZE;
  nblk = dim2e[Labcd][0];
  nthb = dim2e[Labcd][1];
  if (last_ijcs!=-1||last_klcs!=-1) exit(1);

  nwarps = nthb/WARP_SIZE;
//  if (NTHREADS%WARP_SIZE!=0) exit(1);
//  if (nwarps>NTHREADS) exit(1);
//  if (nblk>NBLOCKS) exit(1);
  dim3 dimBlock(WARP_SIZE, nwarps);
  dim3 dimGrid(nblk);
#ifdef CUDA_FMT_M_SM
  Ns = cuda_FMT_m_get_size(4)*sizeof(double);
  Ns += nthb * 1 * sizeof(double);
#else
  Ns = nthb * 3 * sizeof(double);
#endif
#ifdef DLB_KL_DPPS
  Ns += nwarps * sizeof(int);
#endif

  if (NDEV<=0) return 2;
  for (i=0; i<NDEV; i++) {
//    int iwk = iwk0 + i * NBLOCKS;
    int iwk = iwk0 + i;
    dev = cuda_SCF_get_dev_Data(i);
    cuda_set_Device(i);
#ifdef GPU_DLB
    //gpu_twoint_direct_counter_init <<< dimGrid, dimBlock >>> (Labcd, nblk*NIJCSW);
    int c0 = nblk*NIJCSW;
    //gpu_twoint_direct_counter_init <<< dimGrid, dimBlock >>> (Labcd, c0);
    //cudaMemcpyH2D(dev->ijcounter + Labcd, &c0, sizeof(int));
#endif
    gpu_twoint_direct_dpps_ <<< dimGrid, dimBlock, Ns >>> (nwks, iwk,
        eps_eri, eps_ps4, eps_sch);
  }
  cuda_set_Device(0);

  return 0;
}
/* ------------------------------------- */

static inttype TWOINT_INTTYPE = INTTYPE_OS;

inttype cuda_twoint_inttype(inttype type)
{
    if (type>=0) TWOINT_INTTYPE = type;
      return TWOINT_INTTYPE;
}

static int (*cuda_host_twoint_direct[])(
        // paralleization
        const int *pnworkers, const int *pworkerid,
        // integral type data
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        // basis set & cutoff table data
        const int shel_atm[], const int shel_ini[],
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // concerned about buffered direct method
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] ) = {
#if 0
  ofmo_twoint_direct_ssss__,
  // Obara-Saika式（一般式、C言語）
  ofmo_twoint_direct_xxxx, ofmo_twoint_direct_xxxx,
  ofmo_twoint_direct_xxxx, ofmo_twoint_direct_xxxx,
  ofmo_twoint_direct_xxxx, ofmo_twoint_direct_xxxx,
  ofmo_twoint_direct_xxxx, ofmo_twoint_direct_xxxx,
  ofmo_twoint_direct_xxxx, ofmo_twoint_direct_xxxx,
  ofmo_twoint_direct_xxxx, ofmo_twoint_direct_xxxx,
  ofmo_twoint_direct_xxxx, ofmo_twoint_direct_xxxx,
  ofmo_twoint_direct_xxxx, ofmo_twoint_direct_xxxx,
  ofmo_twoint_direct_xxxx, ofmo_twoint_direct_xxxx,
  ofmo_twoint_direct_xxxx, ofmo_twoint_direct_xxxx,
#endif
  ofmo_twoint_direct_ssss__,
  // Obara-Saika式（個別、C言語）
  ofmo_twoint_direct_os_psss, ofmo_twoint_direct_os_psps,
  ofmo_twoint_direct_os_ppss, ofmo_twoint_direct_os_ppps,
  ofmo_twoint_direct_os_pppp, ofmo_twoint_direct_os_dsss,
  ofmo_twoint_direct_os_dsps, ofmo_twoint_direct_os_dspp,
  ofmo_twoint_direct_os_dsds, ofmo_twoint_direct_os_dpss,
  ofmo_twoint_direct_os_dpps, ofmo_twoint_direct_os_dppp,
  ofmo_twoint_direct_os_dpds, ofmo_twoint_direct_os_dpdp,
  ofmo_twoint_direct_os_ddss, ofmo_twoint_direct_os_ddps,
  ofmo_twoint_direct_os_ddpp, ofmo_twoint_direct_os_ddds,
  ofmo_twoint_direct_os_dddp, ofmo_twoint_direct_os_dddd,
#if 0
  ofmo_twoint_direct_ssss__,
  // Rys求積法（一般式、C言語）
  ofmo_twoint_direct_rys_xxxx, ofmo_twoint_direct_rys_xxxx,
  ofmo_twoint_direct_rys_xxxx, ofmo_twoint_direct_rys_xxxx,
  ofmo_twoint_direct_rys_xxxx, ofmo_twoint_direct_rys_xxxx,
  ofmo_twoint_direct_rys_xxxx, ofmo_twoint_direct_rys_xxxx,
  ofmo_twoint_direct_rys_xxxx, ofmo_twoint_direct_rys_xxxx,
  ofmo_twoint_direct_rys_xxxx, ofmo_twoint_direct_rys_xxxx,
  ofmo_twoint_direct_rys_xxxx, ofmo_twoint_direct_rys_xxxx,
  ofmo_twoint_direct_rys_xxxx, ofmo_twoint_direct_rys_xxxx,
  ofmo_twoint_direct_rys_xxxx, ofmo_twoint_direct_rys_xxxx,
  ofmo_twoint_direct_rys_xxxx, ofmo_twoint_direct_rys_xxxx,
  ofmo_twoint_direct_ssss__,
  // Rys求積法（個別、C言語）
  ofmo_twoint_direct_rys_psss, ofmo_twoint_direct_rys_psps,
  ofmo_twoint_direct_rys_ppss, ofmo_twoint_direct_rys_ppps,
  ofmo_twoint_direct_rys_pppp, ofmo_twoint_direct_rys_dsss,
  ofmo_twoint_direct_rys_dsps, ofmo_twoint_direct_rys_dspp,
  ofmo_twoint_direct_rys_dsds, ofmo_twoint_direct_rys_dpss,
  ofmo_twoint_direct_rys_dpps, ofmo_twoint_direct_rys_dppp,
  ofmo_twoint_direct_rys_dpds, ofmo_twoint_direct_rys_dpdp,
  ofmo_twoint_direct_rys_ddss, ofmo_twoint_direct_rys_ddps,
  ofmo_twoint_direct_rys_ddpp, ofmo_twoint_direct_rys_ddds,
  ofmo_twoint_direct_rys_dddp, ofmo_twoint_direct_rys_dddd,
#endif
};


static int (*cuda_twoint_direct[])(
        // paralleization
        const int *pnworkers, const int *pworkerid,
        // integral type data
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        // basis set & cutoff table data
        const int shel_atm[], const int shel_ini[],
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // concerned about buffered direct method
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] ) = {
  //ofmo_twoint_direct_ssss__,
  cuda_twoint_direct_ssss_,
  cuda_twoint_direct_psss_, cuda_twoint_direct_psps_,
  cuda_twoint_direct_ppss_, cuda_twoint_direct_ppps_,
  cuda_twoint_direct_pppp_, cuda_twoint_direct_dsss_,
  cuda_twoint_direct_dsps_, cuda_twoint_direct_dspp_,
  ofmo_twoint_direct_os_dsds, ofmo_twoint_direct_os_dpss,
//  ofmo_twoint_direct_os_dpps, ofmo_twoint_direct_os_dppp,
  cuda_twoint_direct_dpps_, ofmo_twoint_direct_os_dppp,
  ofmo_twoint_direct_os_dpds, ofmo_twoint_direct_os_dpdp,
  ofmo_twoint_direct_os_ddss, ofmo_twoint_direct_os_ddps,
  ofmo_twoint_direct_os_ddpp, ofmo_twoint_direct_os_ddds,
  ofmo_twoint_direct_os_dddp, ofmo_twoint_direct_os_dddd,
};

char *cuda_s2e[] = {
  "ssss", "psss", "psps", "ppss", "ppps", "pppp",
  "dsss", "dsps", "dspp", "dsds", "dpss", "dpps",
  "dppp", "dpds", "dpdp", "ddss", "ddps", "ddpp",
  "ddds", "dddp", "dddd",
  "buff",
};
int cuda_get_num_types(void)
{
  return sizeof(dim2e)/(sizeof(int)*2)-1;
}


#if 0
static double w2e[] = {
   -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,
   -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,
   -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,
   -1.0e0,  -1.0e0,  -1.0e0,
   -1.0e0, // buffered
};
#ifdef _OPENMP
#pragma omp threadprivate(w2e)
#endif

int cuda_get_num_types(void)
{
  return sizeof(w2e)/sizeof(double);
}

void cuda_print_w2e(void) {
  int i;
  int n = cuda_get_num_types();
#ifdef _OPENMP
#pragma omp master
#endif
  {
    printf("--- w2e ---\n");
    for (i=0; i<n; i++) {
      if (w2e[i]>=0) {
        int nb=0, nt=0;
        if (dim2e[i][0]!=0) {
          nb = dim2e[i][0];
          nt = dim2e[i][1];
        }
        printf("%4s: %8.4f (%3d,%3d)\n",s2e[i],w2e[i],nb,nt);
      }
    }
    printf("-----------\n");
  }
}
#endif


__host__ int cuda_calc_twoint_direct(
        const int Labcd,
        // paralleization for CPU
        const int nworkers, const int workerid,
        // integral type data
        const int La, const int Lb, const int Lc, const int Ld,
        // basis set & cutoff table data
        const int shel_atm[], const int shel_ini[],
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // concerned about buffered direct method
        const long max_nzeri, long *p_nzeri,
        double etmp_val[], short int etmp_ind4[],
        const int last_ijcs, const int last_klcs,
        // density matrix & G-matrix data
        const int nao, const double Ds[], double G[])
{
  int NDEV = cuda_get_numDevice();
  double w0,w1;
  int master=TRUE;
  int optsync=cuda_get_optsync();
  //int buffered = (Labcd <= ofmo_twoint_get_global_last_eri_type());
  int gpu = (cuda_use_Device() && cuda_get_optCPU(Labcd)!=0);
  int type = cuda_twoint_inttype(INTTYPE_QUERY);

//  if (optsync) cuda_Barrier();

//  w0 = cuda_Wtime();
//  if (NDEV<0) return 2;
//  if (cuda_use_Device() && dim2e[Labcd][0]!=0 && !buffered) {
  if (gpu) {
    if (NDEV>0&&optsync) checkCudaErrors(cudaDeviceSynchronize());
    //fprintf(stderr,"Labcd: %d\n",Labcd);
        cuda_twoint_direct[Labcd](
                            &nworkers, &workerid,
                            &La, &Lb, &Lc, &Ld,
                            shel_atm, shel_ini,
                            atom_x, atom_y, atom_z,
                            leading_cs_pair,
                            csp_schwarz, csp_ics, csp_jcs,
                            csp_leading_ps_pair,
                            psp_zeta, psp_dkps, psp_xiza,
                            &max_nzeri, p_nzeri,
                            etmp_val, etmp_ind4,
                            &last_ijcs, &last_klcs,
                            &nao, Ds, G );
    if (NDEV>0&&optsync) checkCudaErrors(cudaDeviceSynchronize());
  } else {
    int L = Labcd + type*21;
        cuda_host_twoint_direct[Labcd](
                            &nworkers, &workerid,
                            &La, &Lb, &Lc, &Ld,
                            shel_atm, shel_ini,
                            atom_x, atom_y, atom_z,
                            leading_cs_pair,
                            csp_schwarz, csp_ics, csp_jcs,
                            csp_leading_ps_pair,
                            psp_zeta, psp_dkps, psp_xiza,
                            &max_nzeri, p_nzeri,
                            etmp_val, etmp_ind4,
                            &last_ijcs, &last_klcs,
                            &nao, Ds, G );
  }
//  if (optsync) cuda_Barrier();
//  w1 = cuda_Wtime();
//  if (w2e[Labcd]<0) w2e[Labcd]=0.0e0;
//  w2e[Labcd] += w1-w0;

  return 0;
}
