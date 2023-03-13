#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <time.h>

#include <mpi.h>

#include "ofmo-def.h"
#include "ofmo-prof.h"
#include "ofmo-misc.h"

#include "ofmo-worker-cont.h"

#ifdef DEBUG_MODE
static FILE *fp_debug;
#endif

/* */
// 接続を受け付ける最大ワーカー数
// 接続を受け付けたワーカーのコミュニケータ
// 接続を受け付けたワーカーのID
// 0以上なら接続されている、0未満なら接続されていない
// 各memory serverプロセスが保持しているデータの要素数
// 実数データしか考えていない
// comm_workersとworker_idsは１対1で対応している

//static int NDATA          = 0;
static int NFRAG          = 0;
static int **DATA_NELEMS    = NULL;
static int **DATA_TARGET    = NULL;
static int **DATA_OFFSET    = NULL;

static double **DATA_ENTITY = NULL;

static char port_name[MPI_MAX_PORT_NAME];
static char service_name[MAXSTRLEN];

static MPI_Comm COMM_MASTER;

/* accumulate用 */
static double* DTMP = NULL;
static int*    ITMP = NULL;


static void dealloc_for_acc() {
    Free( DTMP );
    Free( ITMP );
}

/* ACCUMULATE操作用のメモリ確保 
 この関数が呼ばれた時点で、AO pop, Atomic popに
 関するデータが作成されているものとする */
static int alloc_for_acc() {
    static int called = false;
    if ( called ) return 0;
    int nao, natom, n;
    nao    = DATA_NELEMS[OFMO_TOTAL_AOPOP][0];
    natom  = DATA_NELEMS[OFMO_TOTAL_ATPOP][0];
    n      = ( nao>= natom ? nao : natom );
    // このデータサイズは取りすぎだが・・・
    DTMP   = (double*)malloc( sizeof(double) * n );
    ITMP   = (int*)malloc( sizeof(int) * n );
    atexit( dealloc_for_acc );
    called = true;
    return 0;
}

static void dealloc() {
    ofmo_free_imatrix( DATA_TARGET );
    ofmo_free_imatrix( DATA_OFFSET );
    ofmo_free_imatrix( DATA_NELEMS );
    //NDATA = 0;
    NFRAG = 0;
}

static int alloc( int ndata, int nfrag ) {
    static int called = false;
    if ( called ) return 0;
    DATA_TARGET = ofmo_alloc_imatrix( ndata, nfrag );
    DATA_NELEMS = ofmo_alloc_imatrix( ndata, nfrag );
    DATA_OFFSET = ofmo_alloc_imatrix( ndata, nfrag );
    NFRAG       = nfrag;
    atexit( dealloc );
    called = true;
    return 0;
}

/* PUT, GET関数の親玉
TODO: エラーチェックを行うこと
    1. targetと自分のランクの整合性
    2. data_idの値の範囲の整合性
    3. ifragの範囲
    4. workeridの範囲
    5. 受け取った要素数と受信すべき要素数の整合性
   */
static int ofmo_mserv_put( int data_id, int ifrag, int src, int id ) {
    int offset, count, nelem, tfrag;
    int tag_ret=OFMO_TAG_RET, tag_data=OFMO_TAG_DAT;
    double *p;
    MPI_Status status;
    MPI_Comm comm;
    static int called = false, nfao_total, nfatom_total;
    static MPI_Comm *comm_workers;
    if ( ! called ) {
	nfao_total   = ofmo_isum2( NFRAG, DATA_NELEMS[OFMO_AOPOP1] );
	nfatom_total = ofmo_isum2( NFRAG, DATA_NELEMS[OFMO_ATPOP1] );
	comm_workers = ofmo_getadd_comm_workers();
	called = true;
    }

    if ( data_id == OFMO_TOTAL_AOPOP || data_id == OFMO_TOTAL_ATPOP ) {
	tfrag = 0;
	nelem = DATA_NELEMS[data_id][0];
    } else if ( ifrag >= 0 ) {
	tfrag = ifrag;
	nelem = DATA_NELEMS[data_id][ifrag];
    } else {
	tfrag = 0;
	if      ( data_id == OFMO_AOPOP1 || data_id == OFMO_AOPOP2 )
	    nelem = nfao_total;
	else if ( data_id == OFMO_ATPOP1 || data_id == OFMO_ATPOP2 )
	    nelem = nfatom_total;
	else if ( data_id == OFMO_ENERGY || data_id == OFMO_ENERGY0 )
	    nelem = NFRAG;
    }

    offset = DATA_OFFSET[data_id][tfrag];
    p      = &DATA_ENTITY[data_id][offset];
    comm     = ( id<0 ? COMM_MASTER : comm_workers[id] );
    MPI_Recv( p, nelem, MPI_DOUBLE, src, tag_data, comm, &status );
    MPI_Get_count( &status, MPI_DOUBLE, &count );
    if ( count != nelem ) {
#ifdef DEBUG_MODE
	fdbg( fp_debug, "Inconsistency between %d vs %d\n", nelem, count );
	fflush( fp_debug );
#endif
	// debug
	if ( fp_prof ) {
	    fprintf( fp_prof, "%d vs %d\n", count, nelem );
	    fflush( fp_prof );
	}
	MPI_Abort( MPI_COMM_WORLD, 17 );
    }
    MPI_Send( MPI_BOTTOM, 0, MPI_BYTE, src, tag_ret, comm );
    return 0;
}

static int ofmo_mserv_get( int data_id, int ifrag, int dst, int id ) {
    int offset, nelem, tag_data=OFMO_TAG_DAT, tfrag;
    double *p;
    MPI_Comm comm;
    static int called = false, nfao_total, nfatom_total;
    static MPI_Comm *comm_workers;
    if ( ! called ) {
	nfao_total   = ofmo_isum2( NFRAG, DATA_NELEMS[OFMO_AOPOP1] );
	nfatom_total = ofmo_isum2( NFRAG, DATA_NELEMS[OFMO_ATPOP1] );
	comm_workers = ofmo_getadd_comm_workers();
	if ( fp_prof ) {
	    fprintf( fp_prof, "nfao_total = %d, nfatom_total = %d\n",
		    nfao_total, nfatom_total );
	    fflush( fp_prof );
	}
	called = true;
    }
    if ( data_id == OFMO_TOTAL_AOPOP || data_id == OFMO_TOTAL_ATPOP ) {
	tfrag = 0;
	nelem = DATA_NELEMS[data_id][0];
    } else if ( ifrag >= 0 ) {
	tfrag = ifrag;
	nelem = DATA_NELEMS[data_id][ifrag];
    } else {
	tfrag = 0;
	if      ( data_id == OFMO_AOPOP1 || data_id == OFMO_AOPOP2 )
	    nelem = nfao_total;
	else if ( data_id == OFMO_ATPOP1 || data_id == OFMO_ATPOP2 )
	    nelem = nfatom_total;
	else if ( data_id == OFMO_ENERGY || data_id == OFMO_ENERGY0 )
	    nelem = NFRAG;
    }
    offset   = DATA_OFFSET[data_id][tfrag];
    p        = &DATA_ENTITY[data_id][offset];
    comm     = ( id<0 ? COMM_MASTER : comm_workers[id] );
    MPI_Send( p, nelem, MPI_DOUBLE, dst, tag_data, comm );

    return 0;
}

/* 受け取ったデータを要素に加算して更新する処理 */
/* 全体のpopulationデータ (OFMO_TOTAL_AOPOP, OFMO_TOTAL_ATPOP）
 だけが更新の対象である */
static int ofmo_mserv_acc( int data_id, int nelem, int src,
	int id ) {
    static int called = false;
    static MPI_Comm *comm_workers;
    if ( ! called ) {
	alloc_for_acc();	// 一時配列の確保
	comm_workers = ofmo_getadd_comm_workers();
	called = true;
    }
    // 全体population以外は受け付けない
    if ( data_id != OFMO_TOTAL_AOPOP && data_id != OFMO_TOTAL_ATPOP ) {
	if ( fp_prof ) {
	    fdbg( fp_prof, "ERROR: Illegal DATA ID (%d)\n", data_id );
	    fflush( fp_prof );
	}
	return -1;
    }
    if ( id < 0 ) {
	if ( fp_prof ) {
	    fdbg( fp_prof, "ERROR: Illgel worker ID (%d)\n", id );
	    fflush( fp_prof );
	}
	return -1;
    }
    int *iconv   = ITMP, i;
    int tag_ret = OFMO_TAG_RET, tag_data = OFMO_TAG_DAT;
    double *data = DTMP, *p;
    MPI_Comm comm;
    MPI_Status status;
    comm = comm_workers[id];
    p    = &DATA_ENTITY[data_id][0];
    MPI_Recv( data, nelem, MPI_DOUBLE, src, tag_data, comm, &status );
    MPI_Recv( iconv, nelem, MPI_INT, src, tag_data, comm, &status );
    for ( i=0; i<nelem; i++ ) p[ iconv[i] ] += data[i];
    MPI_Send( MPI_BOTTOM, 0, MPI_BYTE, src, tag_ret, comm );
    return 0;
}

/* memory server 本体 */
/* 起動時の引数は、以下の二つ
 * mserv_id =atoi( argv[1] ) : memory server ID
 * maxworkers = atoi( argv[2] ) : 接続を受け付ける最大ワーカー数
 * */
int main( int argc, char* argv[] ) {
    int myrank, nprocs, is_root, root=0;
    int mserv_id, maxworkers, flag;
    //char *p;
    MPI_Status status;
    MPI_Request req_master;
    //int *requested;
    int imsg[OFMO_IMSG_SZ], imsg_sz=OFMO_IMSG_SZ;
    int nfrag, ifrag;
    int ndata;
    MPI_Comm comm_master;
    // プロファイル
    char header[MAXSTRLEN], prof_name[MAXSTRLEN];
    time_t result;
    // debug
    int resultlen;
    char hostname[MPI_MAX_PROCESSOR_NAME];

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
    is_root = ( myrank == root );
    MPI_Comm_get_parent( &comm_master );
    COMM_MASTER = comm_master;
    // ---- 引数をすべてのプロセスにBcast ----
    if ( is_root ) {
	mserv_id   = atoi( argv[1] );
	maxworkers = atoi( argv[2] );
	ndata      = atoi( argv[3] );
	nfrag      = atoi( argv[4] );
	strcpy( service_name, argv[5] );
	if ( mserv_id < 0 ) {
	    //dbg("Illegal number of ID (%d)\n", mserv_id );
	    MPI_Abort( MPI_COMM_WORLD, 1 );
	}
	if ( maxworkers < 1 ) {
	    //dbg("Illegal max. # of workers (%d)\n", maxworkers );
	    MPI_Abort( MPI_COMM_WORLD, 2 );
	}
	if ( ndata < 1 ) {
	    //dbg("Illegal # of data(%d)\n", ndata );
	    MPI_Abort( MPI_COMM_WORLD, 3 );
	}
	if ( nfrag < 1 ) {
	    //dbg("Illegal # of fragments (%d)\n", nfrag );
	    MPI_Abort( MPI_COMM_WORLD, 4 );
	}
	MPI_Open_port( MPI_INFO_NULL, port_name );
	MPI_Publish_name( service_name, MPI_INFO_NULL, port_name );
	// debug
	MPI_Get_processor_name( hostname, &resultlen );
	printf("****** LIST OF HOST NAME (MSERV(ID=%d))******\n", mserv_id );
	printf("rank=%3d : %s\n", root, hostname );
	for ( int irank=0; irank<nprocs; irank++ ) {
	    if ( irank == root ) continue;
	    MPI_Recv( hostname, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, irank,
		    0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
	    printf("rank=%3d : %s\n", irank, hostname );
	}
	printf("********************************\n");
	fflush(stdout);
    } else {
	MPI_Get_processor_name( hostname, &resultlen );
	MPI_Send( hostname, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, root,
		0, MPI_COMM_WORLD );
    }
    MPI_Bcast( &mserv_id, 1, MPI_INT, root, MPI_COMM_WORLD );
    MPI_Bcast( &maxworkers, 1, MPI_INT, root, MPI_COMM_WORLD );
    MPI_Bcast( &ndata, 1, MPI_INT, root, MPI_COMM_WORLD );
    MPI_Bcast( &nfrag, 1, MPI_INT, root, MPI_COMM_WORLD );
    // ---- 初期化処理1 memory serverで管理するデータの保存領域作成 ----
    int master_root = 0, total, irank;
    alloc( ndata, nfrag );
    total = ndata * nfrag;
    MPI_Bcast( DATA_NELEMS[0], total, MPI_INT, master_root, comm_master );
    MPI_Bcast( DATA_TARGET[0], total, MPI_INT, master_root, comm_master );
    MPI_Bcast( DATA_OFFSET[0], total, MPI_INT, master_root, comm_master );
    MPI_Send( MPI_BOTTOM, 0, MPI_BYTE, master_root, 0, comm_master );
    //== ここまでで、初期化に伴うマスターからのデータ受信は終了


    // 各プロセスでのデータサイズ取得
    int **subtotal_nelems, *total_nelems;
    subtotal_nelems = ofmo_alloc_imatrix( nprocs, ndata );
    total_nelems    = (int*)malloc( sizeof(int) * nprocs );
    for ( int i=0; i<ndata*nprocs; i++ ) subtotal_nelems[0][i] = 0;
    for ( int i=0; i<nprocs; i++ )          total_nelems[i] = 0;
    // 各データの各プロセス毎の要素数の総和を求める
    for ( int data_id=0; data_id<ndata; data_id++ ) {
	for ( ifrag=0; ifrag<nfrag; ifrag++ ) {
	    irank = DATA_TARGET[data_id][ifrag];
	    subtotal_nelems[irank][data_id] += DATA_NELEMS[data_id][ifrag];
	}
	for ( irank=0; irank<nprocs; irank++ )
	    total_nelems[irank] += subtotal_nelems[irank][data_id];
    }
    DATA_ENTITY = ofmo_alloc_dmatrixv( ndata, subtotal_nelems[myrank] );
    memset( DATA_ENTITY[0], '\0', sizeof(double)*total_nelems[myrank] );
    ofmo_free_imatrix( subtotal_nelems );
    free( total_nelems );

    // ---- memory server用プロファイルを出力するファイルの初期化 ----
    strcpy( header, "mserv-prof");	// 暫定処置
    ofmo_create_mserv_prof_name( header, mserv_id, prof_name, MAXSTRLEN );
    ofmo_prof_init( prof_name, MPI_COMM_WORLD );
#ifdef DEBUG_MODE
    // debug用ファイルオープン（デバッグが終了したら除去）
    char filename[MAXSTRLEN];
    sprintf( filename, "ppp-%d-%d.log", mserv_id, myrank );
    if ( (fp_debug=fopen(filename, "w")) == NULL ) {
	return -1;
    }
#endif
    if ( fp_prof ) {
	result = time( NULL );
	fprintf( fp_prof, "memory server is created at %s\n",
		asctime( localtime(&result) ) );
	fprintf( fp_prof, "ID of memory server = %d\n", mserv_id );
	fprintf( fp_prof, "max. # of workers   = %d\n", maxworkers );
	fprintf( fp_prof, "port name = %s\n", port_name );
	fprintf( fp_prof, "Initialization of memory server is finished\n");
	fflush( fp_prof );
    }
    // ワーカー情報の保存先を確保する
    ofmo_alloc_worker_container( maxworkers );
    // 各種処理要求を待って、データ処理を行う
    int **imsgs, tag_ret=OFMO_TAG_RET;
    MPI_Request *req_workers;
    MPI_Comm *comm_workers;
    int *workerids, worker_root=0;
    req_workers = (MPI_Request*)malloc( sizeof(MPI_Request) * maxworkers );
    imsgs = ofmo_alloc_imatrix( maxworkers, OFMO_IMSG_SZ );
    comm_workers = ofmo_getadd_comm_workers();
    workerids    = ofmo_getadd_workerids();
    MPI_Irecv( imsg, imsg_sz, MPI_INT, master_root, tag_ret, comm_master,
	    &req_master );

    while (1) {
	// マスタープロセスからのシグナル
	MPI_Test( &req_master, &flag, &status );
	if ( flag ) {
	    if        ( imsg[OFMO_I_CMD] == OFMO_ACCEPT ) {
		int id, workerid;
		MPI_Comm comm_worker;
		workerid = imsg[OFMO_I_METHOD];
		MPI_Barrier( MPI_COMM_WORLD );
		if ( is_root )
		    MPI_Send( MPI_BOTTOM, 0, MPI_BYTE, master_root, tag_ret,
			    comm_master );

		MPI_Comm_accept( port_name, MPI_INFO_NULL, root,
			MPI_COMM_WORLD, &comm_worker );
		id = ofmo_set_workerid( workerid );
		ofmo_set_comm_worker( id, comm_worker );
		// マスタープロセスからの、次の要求を待つ
		MPI_Irecv( imsg, imsg_sz, MPI_INT, master_root, tag_ret,
			comm_master, &req_master );
		// 新たに接続したワーカーからの要求待ちを開始
		MPI_Irecv( imsgs[id], imsg_sz, MPI_INT, MPI_ANY_SOURCE,
			tag_ret, comm_worker, &req_workers[id] );
	    } else if ( imsg[OFMO_I_CMD] == OFMO_REJECT ) {
		int id = -1, workerid;
		workerid = imsg[OFMO_I_METHOD];
		MPI_Barrier( MPI_COMM_WORLD );
		for ( int i=0; i<maxworkers; i++ ) {
		    if ( workerids[i] == workerid ) {
			id = i;
			break;
		    }
		}
		if ( id < 0 ) MPI_Abort( MPI_COMM_WORLD, 112 );
		MPI_Cancel( &req_workers[id] );
		MPI_Request_free( &req_workers[id] );
		MPI_Comm_disconnect( &comm_workers[id] );
		workerids[id] = -1;
		MPI_Barrier( MPI_COMM_WORLD );
		if ( is_root ) {
		    MPI_Send( MPI_BOTTOM, 0, MPI_BYTE, master_root, tag_ret,
			    comm_master );
		}
		MPI_Irecv( imsg, imsg_sz, MPI_INT, master_root, tag_ret,
			comm_master, &req_master );
	    } else if ( imsg[OFMO_I_CMD] == OFMO_BARRIER ) {
		MPI_Barrier( MPI_COMM_WORLD );
		if ( is_root ) {
		    MPI_Send( MPI_BOTTOM, 0, MPI_BYTE, master_root, tag_ret,
			    comm_master );
		}
		MPI_Irecv( imsg, imsg_sz, MPI_INT, master_root, tag_ret,
			comm_master, &req_master );
	    } else if ( imsg[OFMO_I_CMD] == OFMO_GET ) {
		int data_id, dst, ifrag;
		data_id  = imsg[OFMO_I_METHOD];
		ifrag    = imsg[OFMO_I_SCC];
		dst      = status.MPI_SOURCE;
		ofmo_mserv_get( data_id, ifrag, dst, -1 );
		MPI_Irecv( imsg, imsg_sz, MPI_INT, master_root, tag_ret,
			comm_master, &req_master );
	    } else if ( imsg[OFMO_I_CMD] == OFMO_PUT ) {
		int data_id, src, ifrag;
		data_id  = imsg[OFMO_I_METHOD];
		ifrag    = imsg[OFMO_I_SCC];
		src      = status.MPI_SOURCE;
		ofmo_mserv_put( data_id, ifrag, src, -1 );
		MPI_Irecv( imsg, imsg_sz, MPI_INT, master_root, tag_ret,
			comm_master, &req_master );
	    } else if ( imsg[OFMO_I_CMD] == OFMO_FINALIZE ) {
		MPI_Barrier( MPI_COMM_WORLD );
		if ( is_root ) 
		    MPI_Send( MPI_BOTTOM, 0, MPI_BYTE, master_root, tag_ret,
			    comm_master );
		break;
	    }
	}
	// ワーカープロセスからの要求
	for ( int id=0; id<maxworkers; id++ ) {
	    if ( workerids[id] < 0 ) continue;
	    MPI_Test( &req_workers[id], &flag, &status );
	    if ( flag ) {
		// データ更新処理(PUT)
		if ( imsgs[id][OFMO_I_CMD] == OFMO_PUT ) {
		    int data_id, src, ifrag;
		    data_id = imsgs[id][OFMO_I_METHOD];
		    ifrag   = imsgs[id][OFMO_I_SCC];
		    src     = status.MPI_SOURCE;
		    ofmo_mserv_put( data_id, ifrag, src, id );
		// データ取得処理(GET)
		} else if ( imsgs[id][OFMO_I_CMD] == OFMO_GET ) {
		    int data_id, dst, ifrag;
		    data_id  = imsgs[id][OFMO_I_METHOD];
		    ifrag    = imsgs[id][OFMO_I_SCC];
		    dst      = status.MPI_SOURCE;
		    ofmo_mserv_get( data_id, ifrag, dst, id );
		// データ加算更新処理(ACCUMULATE)
		} else if ( imsgs[id][OFMO_I_CMD] == OFMO_ACC ) {
		    int data_id, src, nelem;
		    data_id  = imsgs[id][OFMO_I_METHOD];
		    nelem    = imsgs[id][OFMO_I_SCC];
		    src      = status.MPI_SOURCE;
		    ofmo_mserv_acc( data_id, nelem, src, id );
		} else if ( imsgs[id][OFMO_I_CMD] == OFMO_BARRIER ) {
		    MPI_Barrier( MPI_COMM_WORLD );
		    if ( is_root ) {
			MPI_Send( MPI_BOTTOM, 0, MPI_BYTE, worker_root, tag_ret,
				comm_workers[id] );
		    }
		// ゼロクリア処理
		/*} else if ( imsgs[id][OFMO_I_CMD] == OFMO_ZCR ) {*/
		} else {
		    MPI_Abort( MPI_COMM_WORLD, 100 );
		}	// if ( bufs[workerid][0] );
		MPI_Irecv( imsgs[id], imsg_sz, MPI_INT, MPI_ANY_SOURCE,
			tag_ret, comm_workers[id], &req_workers[id] );
	    }	// if (flag )
	}	// for ( workerid );
    }	// while(1)
    if ( req_master != MPI_REQUEST_NULL ) {
	MPI_Cancel( &req_master );
	MPI_Request_free( &req_master );
    }
    MPI_Barrier( MPI_COMM_WORLD );
    //
    if ( fp_prof ) {
	fprintf( fp_prof, "== recv end signal ==\n");
	fflush( fp_prof );
    }
    
    if ( is_root ) {
	MPI_Unpublish_name( service_name, MPI_INFO_NULL, port_name );
	MPI_Close_port( port_name );
    }
#ifdef DEBUG_MODE
    fprintf( fp_debug, "recv. finish signal\n");
    fflush( fp_debug );
#endif
    for ( int id=0; id<maxworkers; id++ ) {
	if ( workerids[id] < 0 ) continue;
	if ( req_workers[id] != MPI_REQUEST_NULL ) {
	    MPI_Cancel( &req_workers[id] );
	    MPI_Request_free( &req_workers[id] );
	}
    }
    //
    if ( fp_prof ) {
	fprintf( fp_prof, "== all requests have been cancelled ==\n");
	fflush( fp_prof );
    }
    ofmo_free_dmatrix( DATA_ENTITY );
    if ( fp_prof ) {
	result = time( NULL );
	fprintf( fp_prof, "=======================================\n");
	fprintf( fp_prof, "memory server (ID=%d) is finalized at %s\n",
		mserv_id, asctime( localtime(&result) ) );
	fflush( fp_prof );
    }
#ifdef DEBUG_MODE
    fclose ( fp_debug );
#endif
    MPI_Finalize();
    return 0;
}
