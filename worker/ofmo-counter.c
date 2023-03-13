#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include <string.h>
#include <math.h>
#include <limits.h>

#include "ofmo-def.h"
#include "ofmo-prof.h"

#ifdef _OPENMP
#include <omp.h>
#else
#include "omp-dummy.h"
#endif

/* -----------------------------------------------------------------
 * ハイブリッド並列時のワーカー内global counterに関する関数群
 * ----------------------------------------------------------------- */
/* global counter
 * The master thread of master process(rank=0) is special thread for
 * global counter and doesn't perform any jobs.
 * This function must be called from thread-parallel region.
 * MPI_THREAD_SERIALIZED support is needed.
 * */
#define N_INC_STATE 4
#define NCNTR	3

static int nused_counter = 0;

static struct {
    MPI_Comm comm;
    int myrank;
    int nprocs;
    int init_val;
    int current_val;
    int last_val;
    int node_val;
    int finish_flag;
    int maxthreads;
    int inc_level;
    int incs[N_INC_STATE];
    int limit[N_INC_STATE];
} gc[NCNTR];

int ofmo_gc_init(
	const int id,
	MPI_Comm comm, const int init_val,
	const int njobs ) {
    int provided;
    MPI_Query_thread( &provided );
    if ( provided < MPI_THREAD_SERIALIZED ) {
	if ( fp_prof )
	    fdbg( fp_prof,
		    "ERROR: MPI_THREAD_SERIALIZED is not supported\n");
	return -1;
    }
    if ( id<0 || id>=NCNTR ) return -1;
    MPI_Comm_rank( comm, &(gc[id].myrank) );
    MPI_Comm_size( comm, &(gc[id].nprocs) );
    gc[id].comm        = comm;
    gc[id].init_val    = init_val;
    gc[id].current_val = init_val;	// 現在の値
    gc[id].last_val    = init_val + njobs;	// 全体の最後の値
    gc[id].node_val    = init_val;	// 小区間の最後の値
    gc[id].finish_flag = false;
    gc[id].maxthreads  = omp_get_max_threads();

    if ( id == 0 ) {	// for IFC3C
	gc[id].inc_level       = 0;
	gc[id].limit[0]        = gc[id].last_val - (njobs>>1);
	gc[id].limit[1]        = gc[id].last_val - (njobs>>2);
	gc[id].limit[2]        = gc[id].last_val - (njobs>>3);
	gc[id].limit[3]        = INT_MAX;
	gc[id].incs[0]         = gc[id].maxthreads*2;
	gc[id].incs[1]         = gc[id].maxthreads;
	gc[id].incs[2]         = gc[id].maxthreads>>1;
	gc[id].incs[3]         = gc[id].maxthreads>>2;
	if ( gc[id].incs[2] < 2 ) gc[id].incs[2] = 2;
    } else if ( id == 1 ) {	// for IFC2C
	gc[id].inc_level       = 0;
	gc[id].limit[0]        = gc[id].last_val - (njobs>>1);
	gc[id].limit[1]        = gc[id].last_val - (njobs>>2);
	gc[id].limit[2]        = gc[id].last_val - (njobs>>3);
	gc[id].limit[3]        = INT_MAX;
	gc[id].incs[0]         = gc[id].maxthreads*3;
	gc[id].incs[1]         = gc[id].maxthreads*2;
	gc[id].incs[2]         = gc[id].maxthreads;
	gc[id].incs[3]         = gc[id].maxthreads>>1;
    } else if ( id == 2 ) {	// for IFC4C
	gc[id].limit[0]        = gc[id].last_val - (njobs>>1);
	gc[id].limit[1]        = gc[id].last_val - (njobs>>2);
	gc[id].limit[2]        = gc[id].last_val - (njobs>>3);
	gc[id].limit[3]        = INT_MAX;
	gc[id].incs[0]         = gc[id].maxthreads*4;
	gc[id].incs[1]         = gc[id].maxthreads*2;
	gc[id].incs[2]         = gc[id].maxthreads;
	gc[id].incs[3]         = gc[id].maxthreads>>1;
    }
    if ( gc[id].incs[3] < 1 ) gc[id].incs[3] = 1;

    // added
    for ( int i=0; i<N_INC_STATE; i++ ) {
	if ( (init_val+gc[id].incs[i]) < gc[id].limit[i] ) {
	    gc[id].inc_level = i;
	    break;
	}
    }

    nused_counter++;

    return 0;
}

/*
 * カウンタの値を取得する関数
 * スレッド並列領域内で呼び出す
 * */
int ofmo_gc_nxtval( const int id ) {
    int myrank, tag=15, val, mythread;
    int inc, master;
    int buf[2];
    MPI_Status status;
    MPI_Comm comm;

    if ( id<0 || id>=NCNTR ) return -1;
    myrank = gc[id].myrank;
    comm   = gc[id].comm;
    master = 0;
    if ( gc[id].nprocs == 1 ) {
#pragma omp critical
	{
	    val = gc[id].current_val;
	    gc[id].current_val++;
	}
    } else {
	if ( myrank == master ) {
	    mythread = omp_get_thread_num();
	    if ( mythread == 0 ) {	// master-masterの処理
		int nfinished[NCNTR], next, ID;
		if ( nused_counter == 0 ) return INT_MAX;
		for ( int i=0; i<NCNTR; i++ ) nfinished[i] = 0;
		while (1) {
		    MPI_Recv(&ID, 1, MPI_INT, MPI_ANY_SOURCE, tag,
			    comm, &status);
		    inc = gc[ID].incs[ gc[ID].inc_level ];
#pragma omp critical
		    {
			val = gc[ID].current_val;
			gc[ID].current_val += inc;
		    }
		    buf[0] = val; buf[1] = inc;
		    MPI_Send( buf, 2, MPI_INT, status.MPI_SOURCE, tag,
			    comm );

		    next = val + inc;
		    if ( next >= gc[ID].last_val ) nfinished[ID]++;
		    if ( nfinished[ID] >= ( gc[ID].nprocs - 1 ) ) {
			nused_counter--;
			if ( nused_counter == 0 ) break;
		    }
		    if ( next >= gc[ID].limit[gc[ID].inc_level] )
			gc[ID].inc_level++;
		}	// while (1)
		val = INT_MAX;
	    } else {
#pragma omp critical
		{
		    val = gc[id].current_val;
		    gc[id].current_val++;
		}
	    }	// if ( mythread == 0 )
	} else {	// マスタープロセス以外
#pragma omp critical
	    {
		// ローカルにジョブが余っていなかったら
		if ( gc[id].current_val == gc[id].node_val
			&& !gc[id].finish_flag ) {
		    buf[0] = id;
		    MPI_Send( buf, 1, MPI_INT, 0, tag, comm );
		    MPI_Recv( buf, 2, MPI_INT, 0, tag, comm, &status );
		    gc[id].current_val = buf[0];
		    inc                = buf[1];
		    gc[id].node_val = gc[id].current_val + inc;
		    if ( gc[id].node_val >= gc[id].last_val
			    && !gc[id].finish_flag )
			gc[id].finish_flag = true;
		}
		val = gc[id].current_val;
		gc[id].current_val++;
	    }
	}	// if ( myrank == master )
    }	// if ( gc[id].nprocs == 1 )
    return val;
}

