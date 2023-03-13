/** ワーカーからメモリサーバーへのアクセスを
  行うための関数群

  タイムアウト設定を行うようにした
  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "ofmo-def.h"
#include "ofmo-misc.h"
#include "ofmo-data.h"
#include "ofmo-prof.h"

static int DEF_MSERV_ID      = -1;
static int NMSERVS           = 0;
static char **PORT_NAMES     = NULL;
static char **SERVICE_NAMES  = NULL;
static MPI_Comm *COMM_MSERVS = NULL;
static int **DATA_TARGET     = NULL;
static int *USED_MSERVS      = NULL;
static int CALLED            = false;

static void ofmo_mserv_container_finalize() {
    ofmo_free_imatrix( DATA_TARGET );
    ofmo_free_cmatrix( PORT_NAMES );
    ofmo_free_cmatrix( SERVICE_NAMES );
    Free( COMM_MSERVS );
    Free( USED_MSERVS );
    DEF_MSERV_ID = -1;
    NMSERVS      = 0;
    CALLED       = false;
}

int ofmo_mserv_container_init(
	int def_mserv_id, int ndata, int nfrag,
	int** data_target,
	int nmservs, char* service_name_header ) {
    if ( CALLED ) return 0;
    DEF_MSERV_ID  = def_mserv_id;
    DATA_TARGET   = data_target;
    NMSERVS       = nmservs;
    PORT_NAMES    = ofmo_alloc_cmatrix( nmservs, MPI_MAX_PORT_NAME );
    SERVICE_NAMES = ofmo_alloc_cmatrix( nmservs, MPI_MAX_PORT_NAME );
    COMM_MSERVS   = (MPI_Comm*)malloc( sizeof(MPI_Comm) * nmservs );
    USED_MSERVS   = (int*)malloc( sizeof(int) * nmservs );
    for ( int mserv_id=0; mserv_id<nmservs; mserv_id++ ) {
	sprintf( SERVICE_NAMES[mserv_id], "%s-%d",
		service_name_header, mserv_id );
	USED_MSERVS[mserv_id] = false;
    }
    atexit( ofmo_mserv_container_finalize );
    CALLED = true;
    return 0;
}

static int invalid_mserv_id( int mserv_id ) {
    if ( mserv_id < 0 || mserv_id >= NMSERVS ) return true;
    return false;
}

char* ofmo_getadd_service_name( int mserv_id ) {
    if ( invalid_mserv_id( mserv_id ) ) return NULL;
    return SERVICE_NAMES[mserv_id];
}

char* ofmo_getadd_port_name( int mserv_id ) {
    if ( invalid_mserv_id( mserv_id ) ) return NULL;
    return PORT_NAMES[mserv_id];
}

int ofmo_set_comm_mserv( int mserv_id, MPI_Comm comm_mserv ) {
    if ( invalid_mserv_id( mserv_id ) ) return -1;
    COMM_MSERVS[mserv_id] = comm_mserv;
    return 0;
}

int ofmo_set_default_mserv( int mserv_id ) {
    if ( invalid_mserv_id( mserv_id ) ) return -1;
    DEF_MSERV_ID = mserv_id;
    return 0;
}

int ofmo_connect_mserv( int mserv_id, MPI_Comm comm_worker ) {
    char *service_name, *port_name;
    int myrank, worker_root=0, is_root;
    MPI_Comm_rank( comm_worker, &myrank );
    is_root = ( myrank == worker_root );
    service_name = SERVICE_NAMES[mserv_id];
    port_name    = PORT_NAMES[mserv_id];
    if ( is_root )
	MPI_Lookup_name( service_name, MPI_INFO_NULL, port_name );
    MPI_Comm_connect( port_name, MPI_INFO_NULL, worker_root,
	    comm_worker, &COMM_MSERVS[mserv_id] );
    USED_MSERVS[mserv_id] = true;
    return 0;
}

int ofmo_disconnect_mserv( int mserv_id ) {
    if ( invalid_mserv_id( mserv_id ) ) return -1;
    //MPI_Comm_disconnect( &COMM_MSERVS[mserv_id] );
    USED_MSERVS[mserv_id] = false;
    return 0;
}

/** タイムアウト付きのSend/Recv関数 */

/** 与えられたmserv_id以外の有効なmserv_idを返す
  */
static int get_valid_mserv_id() {
    int i;
    for ( i=0; i<NMSERVS; i++ ) {
	if ( USED_MSERVS[i] ) return i;
    }
    return -1;
}

int ofmo_send_with_timeout( void *buf, int count,
	MPI_Datatype datatype, int dest, int tag,
	MPI_Comm comm_mserv, int sec, int usec ) {
    MPI_Status status;
    MPI_Request req;
    int flag, retval;
    double tout, t0, t1;

    tout = (double)sec + 1.e-6*usec;
    t0 = MPI_Wtime();
    MPI_Isend( buf, count, datatype, dest, tag,
	    comm_mserv, &req );
    while (1) {
	MPI_Test( &req, &flag, &status );
	if ( flag ) {
	    retval = 0;
	    break;
	} else {
	    t1 = MPI_Wtime();
	    if ( (t1-t0) > tout ) {
		if ( fp_prof ) {
		    fdbg( fp_prof,
			    "time(%6.3f) has come (%6.3f)\n",
			    tout, (t1-t0));
		    fflush( fp_prof );
		}
		MPI_Cancel( &req );
		retval = -1;
		break;
	    }
	}
    }
    return retval;
}


int ofmo_recv_with_timeout( void *buf, int count,
	MPI_Datatype datatype, int src, int tag,
	MPI_Comm comm_mserv, MPI_Status *status,
	int sec, int usec ) {
    MPI_Request req;
    int flag, retval;
    double tout, t0, t1;

    tout = (double)sec + 1.e-6*usec;
    t0 = MPI_Wtime();
    MPI_Irecv( buf, count, datatype, src, tag,
	    comm_mserv, &req );
    while (1) {
	MPI_Test( &req, &flag, status );
	if ( flag ) {
	    retval = 0;
	    break;
	} else {
	    t1 = MPI_Wtime();
	    if ( (t1-t0) > tout ) {
		if ( fp_prof ) {
		    fdbg( fp_prof,
			    "time(%6.3f) has come (%6.3f)\n",
			    tout, (t1-t0) );
		    fflush( fp_prof );
		}
		MPI_Cancel( &req );
		retval = -1;
		break;
	    }
	}
    }
    return retval;
}

/** メモリサーバーにデータを書き込む関数

  初期化時に決められたメモリサーバーを使用する
*/
int ofmo_worker_put( int data_id, int ifrag, double *src ) {
    static int nfrag, *nfao, *nfatom, nfao_total, nfatom_total;
    MPI_Status status;
    int nao, natom;
    static int called = false;
    if ( !called ) {
	ofmo_data_get_vals("nfrag nfao nfatom nao natom",
		&nfrag, &nfao, &nfatom, &nao, &natom );
	nfao_total   = ofmo_isum2( nfrag, nfao );
	nfatom_total = ofmo_isum2( nfrag, nfatom );
	called = true;
    }
    int nelem, target, mserv_id;
    int imsg[OFMO_IMSG_SZ], imsg_sz=OFMO_IMSG_SZ;
    int tag_ret=OFMO_TAG_RET, tag_data=OFMO_TAG_DAT;
    if ( ifrag >= nfrag ) {
	if ( fp_prof ) {
	    fdbg( fp_prof, "ERROR: Illegal monomer number (%d)\n", ifrag);
	    fflush( fp_prof );
	}
	return -1;
    }
    nelem = 0;
    if ( data_id == OFMO_DENS1 || data_id == OFMO_DENS2 ) {
	if ( ifrag >= 0 ) {
	    int n;
	    n      = nfao[ifrag];
	    nelem  = (n*n+n)>>1;
	    target = DATA_TARGET[data_id][ifrag];
	}
    } else if ( data_id == OFMO_AOPOP1 || data_id == OFMO_AOPOP2 ) {
	nelem = ( ifrag<0 ? nfao_total : nfao[ifrag] );
	target = DATA_TARGET[data_id][0];
    } else if ( data_id == OFMO_ATPOP1 || data_id == OFMO_ATPOP2 ) {
	nelem = ( ifrag<0 ? nfatom_total : nfatom[ifrag] );
	target = DATA_TARGET[data_id][0];
    } else if ( data_id == OFMO_DISTA ) {
	nelem = nfrag;
	target = DATA_TARGET[data_id][ifrag];
    } else if ( data_id == OFMO_ENERGY || data_id == OFMO_ENERGY0 ) {
	nelem = ( ifrag<0 ? nfrag : 1 );
	target = DATA_TARGET[data_id][0];
    } else if ( data_id == OFMO_TOTAL_AOPOP ) {
	nelem  = nao;
	target = DATA_TARGET[data_id][0];
    } else if ( data_id == OFMO_TOTAL_ATPOP ) {
	nelem  = natom;
	target = DATA_TARGET[data_id][0];
    }
    if ( nelem <= 0 ) {
	if ( fp_prof ) {
	    fdbg( fp_prof,
		    "ERROR: Illegal parameter (ifrag=%d, data ID=%d)\n",
		    ifrag, data_id );
	    fflush( fp_prof );
	}
	return -1;
    }
    imsg[OFMO_I_CMD]    = OFMO_PUT;
    imsg[OFMO_I_METHOD] = data_id;
    imsg[OFMO_I_SCC]    = ifrag;

    int ncomp = 0, ierr, nerr;
    for ( mserv_id=0; mserv_id<NMSERVS; mserv_id++ ) {
	if ( USED_MSERVS[mserv_id] == false ) continue;
	nerr = 0;
	ierr = ofmo_send_with_timeout(
		imsg, imsg_sz, MPI_INT, target, tag_ret,
		COMM_MSERVS[mserv_id], 2, 0 );
	if ( ierr != 0 ) nerr++;
	ierr = ofmo_send_with_timeout(
		src, nelem, MPI_DOUBLE, target, tag_data,
		COMM_MSERVS[mserv_id], 2, 0 );
	if ( ierr != 0 ) nerr++;
	ierr = ofmo_recv_with_timeout(
		MPI_BOTTOM, 0, MPI_BYTE, target, tag_ret,
		COMM_MSERVS[mserv_id], &status, 2, 0 );
	if ( ierr != 0 ) nerr++;

	if ( nerr == 0 ) ncomp++;
	else             USED_MSERVS[mserv_id] = false;
	/*MPI_Send( imsg, imsg_sz, MPI_INT, target, tag_ret,
		COMM_MSERVS[mserv_id] );
	MPI_Send( src, nelem, MPI_DOUBLE, target, tag_data,
		COMM_MSERVS[mserv_id] );
	MPI_Recv( MPI_BOTTOM, 0, MPI_BYTE, target, tag_ret,
		COMM_MSERVS[mserv_id], MPI_STATUS_IGNORE );*/
    }
    if ( ncomp == 0 ) {
	if ( fp_prof ) {
	    fdbg( fp_prof, "No mserv is alive\n");
	    fflush( fp_prof );
	}
	MPI_Abort( MPI_COMM_WORLD, 3001 );
    }
    return 0;
}

/** メモリサーバーからデータを読み込む関数
  */
int ofmo_worker_get( int data_id, int ifrag, double *dest ) {
    static int nfrag, *nfao, *nfatom, nfao_total, nfatom_total;
    static int nao, natom;
    static int called = false;
    if ( !called ) {
	int ierr;
	ierr = ofmo_data_get_vals("nfrag nfao nfatom nao natom",
		&nfrag, &nfao, &nfatom, &nao, &natom );
	if ( ierr != 0 ) return -1;
	nfao_total   = ofmo_isum2( nfrag, nfao );
	nfatom_total = ofmo_isum2( nfrag, nfatom );
	called = true;
    }

    int nelem, target, count, ierr = 0;
    int imsg[OFMO_IMSG_SZ], imsg_sz=OFMO_IMSG_SZ;
    int tag_ret=OFMO_TAG_RET, tag_data=OFMO_TAG_DAT;
    MPI_Status status;
    if ( ifrag >= nfrag ) {
	if ( fp_prof ) {
	    fdbg( fp_prof, "ERROR: Illegal monomer number (%d)\n", ifrag);
	    fflush( fp_prof );
	}
	return -1;
    }
    nelem = 0;
    if ( data_id == OFMO_DENS1 || data_id == OFMO_DENS2 ) {
	if ( ifrag >= 0 ) {
	    int n;
	    n = nfao[ifrag];
	    nelem = (n*n+n)>>1;
	    target = DATA_TARGET[data_id][ifrag];
	}
    } else if ( data_id == OFMO_AOPOP1 || data_id == OFMO_AOPOP2 ) {
	nelem = ( ifrag<0 ? nfao_total : nfao[ifrag] );
	target = DATA_TARGET[data_id][0];
    } else if ( data_id == OFMO_ATPOP1 || data_id == OFMO_ATPOP2 ) {
	nelem = ( ifrag<0 ? nfatom_total : nfatom[ifrag] );
	target = DATA_TARGET[data_id][0];
    } else if ( data_id == OFMO_DISTA ) {
	if ( ifrag >= 0 ) {
	    nelem = nfrag;
	    target = DATA_TARGET[data_id][ifrag];
	}
    } else if ( data_id == OFMO_ENERGY || data_id == OFMO_ENERGY0 ) {
	nelem = ( ifrag<0 ? nfrag : 1 );
	target = DATA_TARGET[data_id][0];
    } else if ( data_id == OFMO_TOTAL_AOPOP ) {
	nelem  = nao;
	target = DATA_TARGET[data_id][0];
    } else if ( data_id == OFMO_TOTAL_ATPOP ) {
	nelem  = natom;
	target = DATA_TARGET[data_id][0];
    }
    if ( nelem <= 0 ) {
	if ( fp_prof ) {
	    fdbg( fp_prof,
		    "ERROR: Illegal parameter (ifrag=%d, data ID=%d)\n",
		    ifrag, data_id );
	    fflush( fp_prof );
	}
	return -1;
    }
    imsg[OFMO_I_CMD]    = OFMO_GET;
    imsg[OFMO_I_METHOD] = data_id;
    imsg[OFMO_I_SCC]    = ifrag;
    /*MPI_Send( imsg, imsg_sz, MPI_INT, target, tag_ret,
	    COMM_MSERVS[DEF_MSERV_ID] );
    MPI_Recv( dest, nelem, MPI_DOUBLE, target, tag_data,
	    COMM_MSERVS[DEF_MSERV_ID], &status );*/
    int nerr;
    while (1) {
	nerr = 0;
	ierr = ofmo_send_with_timeout(
		imsg, imsg_sz, MPI_INT, target, tag_ret,
		COMM_MSERVS[DEF_MSERV_ID], 2, 0 );
	if ( ierr != 0 ) nerr++;
	ierr = ofmo_recv_with_timeout(
		dest, nelem, MPI_DOUBLE, target, tag_data,
		COMM_MSERVS[DEF_MSERV_ID], &status, 2, 0 );
	if ( ierr != 0 ) nerr++;
	if ( nerr == 0 ) {
	    ierr = 0;
	    break;
	} else {
	    USED_MSERVS[DEF_MSERV_ID] = false;
	    DEF_MSERV_ID = get_valid_mserv_id();
	    if ( DEF_MSERV_ID < 0 ) {
		ierr = -1;
		if ( fp_prof ) {
		    fdbg( fp_prof, "No mserv is alive\n");
		    fflush( fp_prof );
		}
		MPI_Abort( MPI_COMM_WORLD, 3002 );
	    }
	}
    }
    MPI_Get_count( &status, MPI_DOUBLE, &count );
    if ( count != nelem ) {
	if ( fp_prof ) {
	    fprintf( fp_prof,
		    "ERROR: there are incomsistency\n"
		    "  dataid = %d, ifrag = %d\n"
		    "  %d elems must be recv, but %d recv, mserv_id=%d\n",
		    data_id, ifrag, nelem, count, DEF_MSERV_ID );
	    fflush( fp_prof );
	}
	ierr = -1;
    }
    return ierr;
}

// メモリサーバー上のデータに送信したデータを加算する
// 引数rankで示されたプロセスのみで実行する
int ofmo_worker_acc( int data_id, int nelem, double *data, int *iconv ) {
    int target, mserv_id;
    int imsg[OFMO_IMSG_SZ], imsg_sz=OFMO_IMSG_SZ;
    int tag_ret=OFMO_TAG_RET, tag_data=OFMO_TAG_DAT;
    MPI_Status status;
    if ( data_id != OFMO_TOTAL_AOPOP && data_id != OFMO_TOTAL_ATPOP ) {
	if ( fp_prof ) {
	    fdbg( fp_prof, "ERROR: Illegal data ID (%d)\n", data_id );
	    fflush( fp_prof );
	}
	return -1;
    }
    imsg[OFMO_I_CMD]    = OFMO_ACC;
    imsg[OFMO_I_METHOD] = data_id;
    imsg[OFMO_I_SCC]    = nelem;

    target = DATA_TARGET[data_id][0];

    int ncomp=0, ierr, nerr;
    for ( mserv_id=0; mserv_id<NMSERVS; mserv_id++ ) {
	if ( USED_MSERVS[mserv_id] == false ) continue;
	nerr = 0;
	ierr = ofmo_send_with_timeout(
		imsg, imsg_sz, MPI_INT, target, tag_ret,
		COMM_MSERVS[mserv_id], 2, 0 );
	if ( ierr != 0 ) nerr++;
	ierr = ofmo_send_with_timeout(
		data, nelem, MPI_DOUBLE, target, tag_data,
		COMM_MSERVS[mserv_id], 2, 0 );
	if ( ierr != 0 ) nerr++;
	ierr = ofmo_send_with_timeout(
		iconv, nelem, MPI_INT, target, tag_data,
		COMM_MSERVS[mserv_id], 2, 0 );
	if ( ierr != 0 ) nerr++;
	ierr = ofmo_recv_with_timeout(
		MPI_BOTTOM, 0, MPI_BYTE, target, tag_ret,
		COMM_MSERVS[mserv_id], &status, 2, 0 );
	if ( ierr != 0 ) nerr++;
	if ( nerr == 0 ) ncomp++;
	else             USED_MSERVS[mserv_id] = false;
	/*MPI_Send( imsg, imsg_sz, MPI_INT, target, tag_ret,
		COMM_MSERVS[mserv_id] );
	MPI_Send( data, nelem, MPI_DOUBLE, target, tag_data,
		COMM_MSERVS[mserv_id] );
	MPI_Send( iconv, nelem, MPI_INT, target, tag_data,
		COMM_MSERVS[mserv_id] );
	MPI_Recv( MPI_BOTTOM, 0, MPI_BYTE, target, tag_ret,
		COMM_MSERVS[mserv_id], MPI_STATUS_IGNORE );*/
    }
    if ( ncomp == 0 ) {
	if ( fp_prof ) {
	    fdbg( fp_prof, "No mserv is alive\n");
	    fflush( fp_prof );
	}
	MPI_Abort( MPI_COMM_WORLD, 3003 );
    }
    return 0;
}

int ofmo_worker_barrierp() {
    int mserv_id, nprocs, irank, root=0, ncomp, nerr, ierr;
    int imsg[OFMO_IMSG_SZ], imsg_sz=OFMO_IMSG_SZ;
    int tag_ret=OFMO_TAG_RET, tag_data=OFMO_TAG_DAT;
    MPI_Status status;
    imsg[OFMO_I_CMD]    = OFMO_BARRIER;
    ncomp = 0;
    for ( mserv_id=0; mserv_id<NMSERVS; mserv_id++ ) {
	nerr = 0;
	if ( USED_MSERVS[mserv_id] == false ) continue;
	MPI_Comm_remote_size( COMM_MSERVS[mserv_id], &nprocs );
	for ( irank=0; irank<nprocs; irank++ ) {
	    ierr = ofmo_send_with_timeout(
		    imsg, imsg_sz, MPI_INT, irank, tag_ret,
		    COMM_MSERVS[mserv_id], 1, 0 );
	    if ( ierr != 0 ) nerr++;
	    /*MPI_Send( imsg, imsg_sz, MPI_INT, irank, tag_ret,
		    COMM_MSERVS[mserv_id] );*/
	}
	ierr = ofmo_recv_with_timeout(
		MPI_BOTTOM, 0, MPI_BYTE, root, tag_ret,
		COMM_MSERVS[mserv_id], &status, 1, 0 );
	if ( ierr != 0 ) nerr++;

	if ( nerr == 0 ) ncomp++;
	else USED_MSERVS[mserv_id] = false;
	/*MPI_Recv( MPI_BOTTOM, 0, MPI_BYTE, root, tag_ret,
		COMM_MSERVS[mserv_id], MPI_STATUS_IGNORE );*/
    }
    if ( ncomp == 0 ) {
	if ( fp_prof ) {
	    fdbg( fp_prof, "No mserv is alive\n");
	    fflush( fp_prof );
	}
	MPI_Abort( MPI_COMM_WORLD, 3004 );
    }
    return 0;
}
