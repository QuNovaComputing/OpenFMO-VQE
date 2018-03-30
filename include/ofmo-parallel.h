#ifndef _OFMO_PARALLEL_H_
#define _OFMO_PARALLEL_H_

#ifdef _OPENMP
#include <omp.h>
#else
#include "omp-dummy.h"
#endif

#ifdef USE_MPI
#include <mpi.h>
#else
#include "mpi-dummy.h"
#endif

#endif /* _OFMO_PARALLEL_H_ */
