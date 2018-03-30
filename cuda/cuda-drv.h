#ifndef _CUDA_DRV_H_
#define _CUDA_DRV_H_

#ifdef __cplusplus
extern "C" {
#endif

#ifndef __CUDACC__
/* ------------------------------------- */

#ifdef USE_MPI
#include <mpi.h>
#else
#include "mpi-dummy.h"
#endif

MPI_Comm cuda_Mpi_Comm(void);

int cuda_Init(const int ndev, const int myrank, const int nprocs,
    const MPI_Comm comm);
int cuda_Reconfig(const int myrank, const int nprocs,
        const MPI_Comm comm);
int cuda_Finalize(void);

#else /* __CUDACC__ */
/* ------------------------------------- */

void cuda_Barrier(void);
double cuda_Wtime(void);

/* ------------------------------------- */
#endif /* __CUDACC__ */


#ifdef __cplusplus
}
#endif

#endif /* _CUDALIB_H_ */
