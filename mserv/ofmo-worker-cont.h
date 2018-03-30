#ifndef _OFMO_WORKER_CONT_H_
#define _OFMO_WORKER_CONT_H_

extern int ofmo_alloc_worker_container( int maxworkers );
extern int ofmo_set_workerid( int workerid );
extern int ofmo_set_comm_worker( int id, MPI_Comm comm_worker );

extern MPI_Comm* ofmo_getadd_comm_workers();
extern int* ofmo_getadd_workerids();

#endif
