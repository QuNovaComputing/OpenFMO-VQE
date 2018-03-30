#ifndef _SKEL_RHF_MAIN_H_
#define _SKEL_RHF_MAIN_H_

#include "ofmo-def.h"
#include "ofmo-parallel.h"

struct skel_rhf_config_st {
  MPI_Comm comm;
  int nprocs;
  int myrank;
  int nthreads;
  long eribfsz;
  int maxscfcyc;
  double scfe;
  double scfd;
  float eps_ps4;
  float eps_eri;
  float eps_sch;
  int verbose;
  int dryrun;
  int optsync;
  char input[MAXSTRLEN];
  char density[MAXSTRLEN];
  int ndev;
#ifdef USE_CUDA
  int nblk;
  int nthb;
  int type;
#endif
};
typedef struct skel_rhf_config_st skel_rhf_config_t;

#endif /* _SKEL_RHF_MAIN_H_ */

