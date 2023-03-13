//
// cuda driver routines in C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <strings.h>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USE_MPI
#include <mpi.h>
#else
#include "mpi-dummy.h"
#endif

#include "cuda-drv.h"
#include "cudalib.h"

extern FILE* fp_prof; // from common/ofmo-prof.h


static MPI_Comm CUDA_MPI_COMM = MPI_COMM_NULL;


/* ------------------------------------- */

MPI_Comm cuda_Mpi_Comm(void)
{
  return CUDA_MPI_COMM;
}

/* ------------------------------------- */

int cuda_Init(const int ndev, const int myrank, const int nprocs,
    const MPI_Comm comm)
{
  int ret = 0;

#ifdef _OPENMP
#pragma omp master
#endif
  {
  CUDA_MPI_COMM = comm;

 ret = cuda_Init_Sub(ndev, myrank, nprocs);
  }
#ifdef _OPENMP
#pragma omp barrier
#endif
 return ret;
}

int cuda_Reconfig(const int myrank, const int nprocs,
    const MPI_Comm comm)
{
  int ret = 0;
  int ndev = -1;

#ifdef _OPENMP
#pragma omp master
#endif
  {
  CUDA_MPI_COMM = comm;

 ret = cuda_Init_Sub(ndev, myrank, nprocs);
  }
#ifdef _OPENMP
#pragma omp barrier
#endif
 return ret;
}

int cuda_Finalize(void)
{
  int ret = 0;
#ifdef _OPENMP
#pragma omp master
#endif
  {
  ret = cuda_Finalize_Sub();

  CUDA_MPI_COMM = MPI_COMM_WORLD;
  }
#ifdef _OPENMP
#pragma omp barrier
#endif
  return ret;
}

/* ------------------------------------- */

void cuda_Barrier(void)
{
  MPI_Comm comm = cuda_Mpi_Comm();
  int master = TRUE;

#ifdef _OPENMP
  master = (omp_get_thread_num() == 0);
#pragma omp barrier
#endif
#ifdef USE_MPI
  if(master) MPI_Barrier(comm);
#endif
#ifdef _OPENMP
#pragma omp barrier
#endif
}

double cuda_Wtime(void) {
#ifdef USE_MPI
  return MPI_Wtime();
#else
  struct timeval TV;
  gettimeofday(&TV, NULL);
  return(TV.tv_sec + (double)TV.tv_usec*1.0e-6);
#endif
}

/* ------------------------------------- */
