#include <stdio.h>
#include <stdlib.h>

#include "ofmo-parallel.h"
#include "ofmo-def.h"
#include "skel-w2e.h"

#ifdef USE_CUDA
#include "cuda/cuda-integ.h"
#include "cuda/cudalib.h"
#endif

static double w2e[] = {
   -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,
   -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,
   -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,  -1.0e0,
   -1.0e0,  -1.0e0,  -1.0e0,
   -1.0e0,  // buffered
};
#ifdef _OPENMP
#pragma omp threadprivate(w2e)
#endif

static char *s2e[] = {
  "ssss", "psss", "psps", "ppss", "ppps", "pppp",
  "dsss", "dsps", "dspp", "dsds", "dpss", "dpps",
  "dppp", "dpds", "dpdp", "ddss", "ddps", "ddpp",
  "ddds", "dddp", "dddd",
  "buff"
};


static MPI_Comm comm = MPI_COMM_WORLD;
static int sync = false;

void setup_w2e(const MPI_Comm wcomm, const int optsync)
{
  comm = wcomm;
  sync = optsync;
}

static void w2e_Barrier(void)
{
#pragma omp barrier
#ifdef USE_MPI
#pragma omp master
    {
      MPI_Barrier(comm);
    }
#pragma omp barrier
#endif
}

static double t0=0.0;
void start_w2e(void)
{
  if (sync) w2e_Barrier();
  t0 = MPI_Wtime();
}

void set_w2e(int Labcd)
{
  double t1;
  t1 = MPI_Wtime();
  if (sync) w2e_Barrier();
  if (Labcd<0) {Labcd = sizeof(w2e)/sizeof(double) - 1;}
  if (w2e[Labcd]<0) w2e[Labcd] = 0.0;
  w2e[Labcd] += t1-t0;
}

void print_w2e(void)
{
  int i;
  int n = sizeof(w2e)/sizeof(double);
#pragma omp master
  {
    printf("--- w2e ---\n");
    for (i=0; i<n; i++) {
//      if (w2e[i]>=0) {
#ifndef USE_CUDA
        printf("%4s: %8.4f\n",s2e[i],w2e[i]);
#else
        int nb=0, nt=0;
        if (dim2e[i][0]!=0&&cuda_get_numBlocks()!=0) {
          nb = dim2e[i][0];
          nt = dim2e[i][1];
        }
        printf("%4s: %8.4f (%3d,%3d)\n",s2e[i],w2e[i],nb,nt);
#endif
//      }
    }
    printf("-----------\n");
  }
}


