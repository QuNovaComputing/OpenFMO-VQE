#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "ofmo-def.h"

static int MAXWORKERS           = 0;
static int *WORKERIDS           = NULL;
static MPI_Comm *COMM_WORKERS   = NULL;

static void ofmo_dealloc_worker_container() {
    Free( WORKERIDS );
    Free( COMM_WORKERS );
    MAXWORKERS = 0;
}

int ofmo_alloc_worker_container( int maxworkers ) {
    static int called = false;
    if ( called ) return 0;
    MAXWORKERS = maxworkers;
    WORKERIDS  = (int*)malloc( sizeof(int) * maxworkers );
    COMM_WORKERS = (MPI_Comm*)malloc( sizeof(MPI_Comm) * maxworkers );
    for ( int i=0; i<maxworkers; i++ ) WORKERIDS[i] = -1;
    atexit( ofmo_dealloc_worker_container );
    called = true;
    return 0;
}

static int get_id() {
    int id;
    for ( id=0; id<MAXWORKERS; id++ )
	if ( WORKERIDS[id] == -1 ) return id;
    return -1;
}

int ofmo_set_workerid( int workerid ) {
    int id;
    if ( (id=get_id() ) < 0 ) return -1;
    WORKERIDS[ id ] = workerid;
    return id;
}

int ofmo_set_comm_worker( int id, MPI_Comm comm_worker ) {
    COMM_WORKERS[id] = comm_worker;
    return 0;
}

MPI_Comm* ofmo_getadd_comm_workers() {
    return COMM_WORKERS;
}

int* ofmo_getadd_workerids() {
    return WORKERIDS;
}
