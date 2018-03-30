/* マスタープロセス
   １プロセスであることを前提とする
   RPCを意識した呼び出しにする
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <libgen.h>

#include <math.h>

#include <mpi.h>

#ifdef _OPENMP
#include <omp.h>
#else
#include "omp-dummy.h"
#endif

#include "ofmo-def.h"
#include "ofmo-data.h"

#define MAXTOKLEN MAXSTRLEN

#include "ofmo-string.h"
#include "ofmo-misc.h"

#include "ofmo-data-struct.h"
#include "ofmo-put-get-master.h"

/* 外部で定義された関数 */
extern int ofmo_init( const char* filename, MPI_Comm comm );

extern unsigned char * ofmo_alloc_compressed_data( int* );

extern int ofmo_compress_data( const char filename[],
	unsigned char comped_data[], const int max_comped_data_size );

#define ERROR_IN_MPI 99

/* -----------------------------------------------------
 * モノマーをAO数の大きい順に並び替えたリストの作成
 * ----------------------------------------------------- */
/** ２つの整数を比較する関数
 * */
static int comp2( const void *p1, const void *p2 ) {
    return ( (*(int*)p2) - (*(int*)p1) );	// good
}

static int *frag_order = NULL;

static void dealloc_frag_order() {
    if ( frag_order != NULL ) free( frag_order );
    frag_order = NULL;
}

static int ofmo_init_frag_order() {
    static int called = false;
    if ( called ) return 0;
    int nfrag, i, i2, *nfao;
    if ( ofmo_data_get_vals("nfrag nfao", &nfrag, &nfao) != 0 ) {
	dbg("error\n");
	return -1;
    }
    frag_order = (int*)malloc( sizeof(int) * nfrag *2 );
    if ( frag_order == NULL ) return -1;
    for ( i=0, i2=0; i<nfrag; i++, i2+=2 ) {
	frag_order[i2+0] = nfao[i];
	frag_order[i2+1] = i;
    }
    qsort( frag_order, nfrag, sizeof(int)*2, comp2 );
    for ( i=0, i2=0; i<nfrag; i++, i2+=2 )
	frag_order[i] = frag_order[i2+1];
    atexit( dealloc_frag_order );
    called = true;
    return 0;
}
// -----------------------------------------------------------------
static void ofmo_get_service_name_header( char *header ) {
    char *p;
    if ( (p=getenv("OFMO_MSERV_NAME")) != NULL ) {
	strncpy( header, p, MAXSTRLEN );
    } else {
	snprintf( header, MAXSTRLEN, "mserv-%d", (int)getpid() );
    }
}

// ----- SCC収束条件を決める関数 -----
static int ofmo_get_monomer_scfconv( double dE ) {
    static int scfconv_now = 2, called=false, scfconv;
    static double sccconvd;
    static int nit=0;
    int scfconv_old;
    double de;
    if ( !called ) {
        int sccconv;
        double scfconvd;
        ofmo_data_get_vals("scfconv sccconv", &scfconv, &sccconv );
	called = true;
        scfconvd = pow(10,-scfconv);
        sccconvd = pow(10,-sccconv);
        if (scfconvd>sccconvd) 
            printf("Warning: CONV in $SCF (%8.2e) > CONV in $FMOPRP (%8.2e)!\n",
                   scfconvd, sccconvd);
    }
    de = fabs(dE);
    scfconv_old = scfconv_now;
    nit++;
#if 0
    //if      ( scfconv_now <= 4 && de < 1.e-4 ) scfconv_now = scfconv;
    if      ( scfconv_now <= 8 && de < 1.e-8 ) scfconv_now = scfconv;
    else if ( scfconv_now <= 7 && de < 1.e-7 ) scfconv_now = 8;
    else if ( scfconv_now <= 6 && de < 1.e-6 ) scfconv_now = 7;
    else if ( scfconv_now <= 5 && de < 1.e-5 ) scfconv_now = 6;
    else if ( scfconv_now <= 4 && de < 1.e-4 ) scfconv_now = 5;
    else if ( scfconv_now <= 3 && de < 1.e-3 ) scfconv_now = 4;
    else if ( scfconv_now <= 2 && de < 1.e-1 ) scfconv_now = 3;
#else
    if      ( scfconv_now == 2 && (nit>5 || de < 1.e-1) ) scfconv_now+=2;
    else if ( scfconv_now >= 3 && (nit>3 || de < pow(10,-scfconv_now))) scfconv_now+=2;
#endif
    if (scfconv_now > scfconv || de < 10*sccconvd) scfconv_now = scfconv;
    if ( scfconv_now != scfconv_old ) nit=1;
    return scfconv_now;
}

/* 配列に対するaccumulate処理を行う関数
 */
static void ofmo_accumulate( const int n, const double data[],
	const int index[], double target[] ) {
    int i;
    for ( i=0; i<n; i++ ) target[ index[i] ] += data[i];
}

/* ********************* メモリサーバー関連 *************** */
/** メモリサーバーを起動し、初期化する関数

  @retval 正常終了時 メモリサーバーID(非負整数)
  @retval 異常終了時 -1
  */
static int ofmo_create_mserv(
	int mserv_size, int maxworkers, int ndata, int nfrag,
	int *nelems, int *target, int* offset,
        char *bindir,
	char *service_name_header, int mserv_id, MPI_Comm *comm_mserv ) {
    int root=0, total, ierr;
    char *args[6], smserv_id[MAXSTRLEN], smaxworkers[MAXSTRLEN];
    char sndata[MAXSTRLEN], service_name[MAXSTRLEN], snfrag[MAXSTRLEN];
    char binfname[MAXSTRLEN];

    snprintf( smaxworkers, MAXSTRLEN, "%d", maxworkers );
    snprintf( smserv_id, MAXSTRLEN, "%d", mserv_id );
    snprintf( sndata, MAXSTRLEN, "%d", ndata );
    snprintf( snfrag, MAXSTRLEN, "%d", nfrag );
    snprintf( service_name, MAXSTRLEN, "%s-%d", service_name_header, mserv_id );
    args[0] = smserv_id;
    args[1] = smaxworkers;
    args[2] = sndata;
    args[3] = snfrag;
    args[4] = service_name;
    args[5] = NULL;

    snprintf(binfname, MAXSTRLEN, "%s/ofmo-mserv", bindir);
    //ierr = MPI_Comm_spawn("./ofmo-mserv", args, mserv_size, MPI_INFO_NULL,
    ierr = MPI_Comm_spawn(binfname, args, mserv_size, MPI_INFO_NULL,
	    root, MPI_COMM_WORLD, comm_mserv, MPI_ERRCODES_IGNORE );
    if ( ierr != MPI_SUCCESS ) {
	dbg("Failure in invoking %d-th memory server\n", mserv_id );
	return -1;
    }
    total = nfrag * ndata;
    ierr = MPI_Bcast( nelems, total, MPI_INT, MPI_ROOT, *comm_mserv );
    if ( ierr != MPI_SUCCESS ) { dbg("nelems\n");}
    ierr = MPI_Bcast( target, total, MPI_INT, MPI_ROOT, *comm_mserv );
    if ( ierr != MPI_SUCCESS ) { dbg("target\n");}
    ierr = MPI_Bcast( offset, total, MPI_INT, MPI_ROOT, *comm_mserv );
    if ( ierr != MPI_SUCCESS ) { dbg("offset\n");}

    ierr = MPI_Recv( MPI_BOTTOM, 0, MPI_BYTE, 0, 0, *comm_mserv,
	    MPI_STATUS_IGNORE );
    if ( ierr != MPI_SUCCESS ) { dbg("recv\n");}
    return 0;
}

// 同期でメモリサーバーにジョブを渡す
static int ofmo_call_mserv( MPI_Comm comm_mserv,
	int imsg_sz, int imsg[] ) {
    int mserv_sz, irank, tag_ret=OFMO_TAG_RET, mserv_root=0;
    int rc;
    rc = MPI_Comm_remote_size( comm_mserv, &mserv_sz );
    if (rc != MPI_SUCCESS) { MPI_Abort( comm_mserv, ERROR_IN_MPI ); };
    for ( irank=0; irank<mserv_sz; irank++ ) {
	rc = MPI_Send( imsg, imsg_sz, MPI_INT, irank, tag_ret, comm_mserv );
        if (rc != MPI_SUCCESS) { MPI_Abort( comm_mserv, ERROR_IN_MPI ); };
    }

    rc = MPI_Recv( MPI_BOTTOM, 0, MPI_BYTE, mserv_root, tag_ret,
	    comm_mserv, MPI_STATUS_IGNORE );
    if (rc != MPI_SUCCESS) { MPI_Abort( comm_mserv, ERROR_IN_MPI ); };

    return 0;
}

static int ofmo_destroy_mserv( MPI_Comm comm_mserv ) {
    int imsg[OFMO_IMSG_SZ];
    imsg[OFMO_I_CMD] = OFMO_FINALIZE;
    ofmo_call_mserv( comm_mserv, OFMO_IMSG_SZ, imsg );
    return 0;
}

static int ofmo_mserv_barrier( int nmservs,
	MPI_Comm comm_mservs[] ) {
    int mserv_id, imsg[OFMO_IMSG_SZ];
    imsg[OFMO_I_CMD] = OFMO_BARRIER;
    for ( mserv_id=0; mserv_id<nmservs; mserv_id++ ) {
	ofmo_call_mserv( comm_mservs[mserv_id], OFMO_IMSG_SZ, imsg );
    }
    return 0;
}

// ********************* ワーカー関連 **********************
/** ワーカーの生成+初期化関数
  */
static int ofmo_create_worker( int worker_size, int workerid,
    long eribfsz, int ndev,
	int comped_data_size, unsigned char* comped_data,
	char *input, char* header, char* local_input_dir,
	int nmservs, int ndata, int nfrag, int *target, 
        char *bindir,
	char *service_name_header, int def_mserv_id,
	MPI_Comm *comm_worker ) {
    int ierr, root=0;
    char *args[13];
    char sworkerid[MAXSTRLEN], snmservs[MAXSTRLEN];
    char scomped_data_size[MAXSTRLEN], sndata[MAXSTRLEN];
    char sdef_mserv_id[MAXSTRLEN], snfrag[MAXSTRLEN];
    char binfname[MAXSTRLEN];
    char seribfsz[MAXSTRLEN], sndev[MAXSTRLEN];

    snprintf( sworkerid, MAXSTRLEN, "%d", workerid );
    snprintf( snmservs , MAXSTRLEN, "%d", nmservs );
    snprintf( scomped_data_size, MAXSTRLEN, "%d", comped_data_size );
    snprintf( sndata, MAXSTRLEN, "%d", ndata );
    snprintf( snfrag, MAXSTRLEN, "%d", nfrag );
    snprintf( sdef_mserv_id, MAXSTRLEN, "%d", def_mserv_id );
    snprintf( seribfsz, MAXSTRLEN, "%ld", eribfsz );
    snprintf( sndev, MAXSTRLEN, "%d", ndev );
    args[0] = sworkerid;
    args[1] = snmservs;
    args[2] = scomped_data_size;
    args[3] = sndata;
    args[4] = snfrag;
    args[5] = sdef_mserv_id;
    args[6] = input;
    args[7] = header;
    args[8] = local_input_dir;
    args[9] = service_name_header;
    args[10] = seribfsz;
    args[11] = sndev;
    args[12] = NULL;

    snprintf(binfname, MAXSTRLEN, "%s/ofmo-worker", bindir);
    //ierr = MPI_Comm_spawn( "./ofmo-worker", args, worker_size,
    ierr = MPI_Comm_spawn( binfname, args, worker_size,
	    MPI_INFO_NULL, root, MPI_COMM_WORLD, comm_worker,
	    MPI_ERRCODES_IGNORE );
    if ( ierr != MPI_SUCCESS ) {
	dbg("Failure in invoking %d-th worker\n", workerid );
	return -1;
    }
    ierr = MPI_Bcast( comped_data, comped_data_size, MPI_UNSIGNED_CHAR,
	    MPI_ROOT, *comm_worker );
    if ( ierr != MPI_SUCCESS ) { return -1;}
    ierr = MPI_Bcast( target, ndata*nfrag, MPI_INT,
	    MPI_ROOT, *comm_worker );
    if ( ierr != MPI_SUCCESS ) { return -1;}
    return 0;
}

/* ************ OmniRPCの関数のまね **************** */
// 非同期で、ワーカーにジョブを渡す
static MPI_Request ofmo_call_async_worker( MPI_Comm comm_worker,
	int imsg_sz, int imsg[],
	int dmsg_sz, double dmsg[] ) {
    int tag_ret = OFMO_TAG_RET, worker_root = 0;
    int rc;
    MPI_Request req;
    rc = MPI_Bcast( imsg, imsg_sz, MPI_INT, MPI_ROOT, comm_worker );
    if (rc != MPI_SUCCESS) { MPI_Abort( comm_worker, ERROR_IN_MPI ); };
    rc = MPI_Irecv( dmsg, dmsg_sz, MPI_DOUBLE, worker_root, tag_ret,
	    comm_worker, &req );
    if (rc != MPI_SUCCESS) { MPI_Abort( comm_worker, ERROR_IN_MPI ); };
    return req;
}

// すべてのワーカーからのジョブを待つ
static void ofmo_waitall_worker( int nworkers, MPI_Request reqs[] ) {
    int count = 0, workerid;
    MPI_Status status;
    int rc;
    while (1) {
	rc = MPI_Waitany( nworkers, reqs, &workerid, &status );
        if (rc != MPI_SUCCESS) { MPI_Abort( MPI_COMM_WORLD, ERROR_IN_MPI ); };
//        if (status.MPI_ERROR == MPI_ERR_IN_STATUS) { MPI_Abort( MPI_COMM_WORLD, ERROR_IN_MPI ); };
	count++;
	if ( count == nworkers ) break;
    }
}

static int ofmo_waitany_worker( int nworkers, MPI_Request reqs[] ) {
    MPI_Status status;
    int workerid;
    int rc;
    rc = MPI_Waitany( nworkers, reqs, &workerid, &status );
    if (rc != MPI_SUCCESS) { MPI_Abort( MPI_COMM_WORLD, ERROR_IN_MPI ); };
//    if (status.MPI_ERROR == MPI_ERR_IN_STATUS) { MPI_Abort( MPI_COMM_WORLD, ERROR_IN_MPI ); };
    return workerid;
}

static void ofmo_wait_worker( MPI_Request req ) {
    MPI_Status status;
    int rc;
    rc = MPI_Wait( &req, &status );
    if (rc != MPI_SUCCESS) { MPI_Abort( MPI_COMM_WORLD, ERROR_IN_MPI ); };
//    if (status.MPI_ERROR == MPI_ERR_IN_STATUS) { MPI_Abort( MPI_COMM_WORLD, ERROR_IN_MPI ); };
}


/* ******** 上記の関数の組み合わせ ************* */
// ワーカーにジョブを依頼して、それが終わるまで待つ
static int ofmo_call_worker( MPI_Comm comm_worker,
	int imsg_sz, int imsg[],
	int dmsg_sz, double dmsg[] ) {
    MPI_Request req;
    req = ofmo_call_async_worker( comm_worker,
	    imsg_sz, imsg, dmsg_sz, dmsg );
    ofmo_wait_worker( req );
    return 0;
}

static void ofmo_destroy_worker( MPI_Comm comm_worker ) {
    int imsg[OFMO_IMSG_SZ], imsg_sz=OFMO_IMSG_SZ, dmsg_sz=OFMO_DMSG_SZ;
    double dmsg[OFMO_DMSG_SZ];
    imsg[OFMO_I_CMD] = OFMO_FINALIZE;
    ofmo_call_worker( comm_worker, imsg_sz, imsg, dmsg_sz, dmsg );
}

// すべてのワーカーがジョブ待ち状態になるのを待つ
static void ofmo_worker_barrier( int nworkers, MPI_Comm comm_workers[] ) {
    int workerid, imsg[OFMO_IMSG_SZ];
    double dmsg[OFMO_DMSG_SZ];
    MPI_Request req;
    imsg[OFMO_I_CMD] = OFMO_BARRIER;
    for ( workerid=0; workerid<nworkers; workerid++ ) {
	req = ofmo_call_async_worker( comm_workers[workerid],
		OFMO_IMSG_SZ, imsg, OFMO_DMSG_SZ, dmsg );
	ofmo_wait_worker( req );
    }
}

static int ofmo_update_monomer_data_worker(
	int nworkers, MPI_Comm comm_workers[] ) {
    int workerid, imsg[OFMO_IMSG_SZ];
    double dmsg[OFMO_DMSG_SZ];
    MPI_Request req;
    imsg[OFMO_I_CMD] = OFMO_UPDATE_DATA;
    for ( workerid=0; workerid<nworkers; workerid++ ) {
	req = ofmo_call_async_worker( comm_workers[workerid],
		OFMO_IMSG_SZ, imsg, OFMO_DMSG_SZ, dmsg );
	ofmo_wait_worker( req );
    }
    return 0;
}

static int ofmo_reset_monomer_data_worker(
	int nworkers, MPI_Comm comm_workers[] ) {
    int workerid, imsg[OFMO_IMSG_SZ];
    double dmsg[OFMO_DMSG_SZ];
    MPI_Request req;
    imsg[OFMO_I_CMD] = OFMO_RESET_MON_DATA;
    for ( workerid=0; workerid<nworkers; workerid++ ) {
	req = ofmo_call_async_worker( comm_workers[workerid],
		OFMO_IMSG_SZ, imsg, OFMO_DMSG_SZ, dmsg );
	ofmo_wait_worker( req );
    }
    return 0;
}

/* ワーカから返される値が特別な配列である場合の処理
   (populationデータをワーカから返してもらうとき）
   他の関数では、ジョブをワーカに発行すると、dmsgという配列に
   決まった量のデータが返される、というパターンをとっている。
   しかし、ワーカから直前の計算で得られたpopulationデータを返して
   もらう場合には、データ長がまちまちであり、また、実数、整数データが
   入り混じっているなど、返される値の著しく異なる。
   そのような理由から、別関数を用意している。
*/
static int ofmo_get_pop_data(
	MPI_Comm comm_worker, int mode,
	int bufsz_aop, int *fnao, double aopop_frg[], int fsao2tuao[],
	int bufsz_atp, int *fnat, double atpop_frg[], int fatom2tatom[] ) {
    int imsg[OFMO_IMSG_SZ], sz1, sz2, tag_ret=OFMO_TAG_RET, worker_root=0;
    int ierr = 0;
    MPI_Status status;
    imsg[OFMO_I_CMD] = OFMO_GET_POP_DATA;
    imsg[OFMO_I_METHOD] = mode;
    *fnao = *fnat = 0;
    MPI_Bcast( imsg, OFMO_IMSG_SZ, MPI_INT, MPI_ROOT, comm_worker );
    if ( mode == OFMO_BOTH_POPS || mode == OFMO_AOPOP_ONLY ) {
	MPI_Recv( aopop_frg, bufsz_aop, MPI_DOUBLE, worker_root, tag_ret,
		comm_worker, &status );
	MPI_Get_count( &status, MPI_DOUBLE, &sz1 );
	MPI_Recv( fsao2tuao, bufsz_aop, MPI_INT, worker_root, tag_ret,
		comm_worker, &status );
	MPI_Get_count( &status, MPI_INT, &sz2 );
	if ( sz1 != sz2 ) {
	    dbg("error: size inconsistency between fragment AO population data"
		    " (%d vs. %d)\n", sz1, sz2 );
	    fflush(stdout);
	    ierr = -1;
	} else {
	    *fnao = sz1;
	}
    }
    if ( mode == OFMO_BOTH_POPS || mode == OFMO_ATPOP_ONLY ) {
	MPI_Recv( atpop_frg, bufsz_atp, MPI_DOUBLE, worker_root, tag_ret,
		comm_worker, &status );
	MPI_Get_count( &status, MPI_DOUBLE, &sz1 );
	MPI_Recv( fatom2tatom, bufsz_atp, MPI_INT, worker_root, tag_ret,
		comm_worker, &status );
	MPI_Get_count( &status, MPI_INT, &sz2 );
	if ( sz1 != sz2 ) {
	    dbg("error: size inconsistency between fragment atomic population data"
		    " (%d vs. %d)\n", sz1, sz2 );
	    fflush(stdout);
	    ierr = -1;
	} else {
	    *fnat = sz1;
	}
    }
    return ierr;
}

void show_help(const char *myname, const int verbose)
{
  fprintf(stderr,"Usage: %s [options] [input [InitDens]]\n", myname);
  fprintf(stderr,"  -ng #: # groups\n");
  fprintf(stderr,"  -np #: # total MPI procs\n");
  fprintf(stderr,"  -B #: buffer size / proc (MB, default: %d)\n", 512);
  fprintf(stderr,"  -v: verbose\n");
  fprintf(stderr,"  -h: show this help\n");
#ifdef USE_CUDA
  fprintf(stderr," Options for GPGPU:\n");
  fprintf(stderr,"  -d #: # devices (default:1)\n");
#endif
}

/** マスターのメインプログラム

  １プロセスで起動することを前提としている
  （複数プロセス起動しても、rank=0のプロセス以外は何もしない）
   */
#ifdef FJ_MAIN
int MAIN__( int argc, char *argv[] ) {
#else
int main( int argc, char *argv[] ) {
#endif
    int nprocs;
    double et0, et1, ET0;
    char *p, input[MAXSTRLEN], header[MAXSTRLEN];
//    char port_name_prefix[MAXSTRLEN];
    char local_input_dir[MAXSTRLEN];
    char bindir[MAXSTRLEN];
    char dens[MAXSTRLEN];
    int ngroup=-1, nioprocs=-1, niogroup=-1, ierr, nmaxprocs=-1;
    int ngroupA, nioprocsA, niogroupA, iopos, master;
    int group_size;
    int imsg[OFMO_IMSG_SZ], imsg_sz=OFMO_IMSG_SZ, dmsg_sz=OFMO_DMSG_SZ;
    double dmsg[OFMO_DMSG_SZ];
    int ndev;
    long eribfsz;
    // for debug
    int resultlen;
    char hostname[MPI_MAX_PROCESSOR_NAME];

    // ****************** マスタープロセスの初期化 *********************
    MPI_Init( &argc, &argv );
    ET0 = et0 = MPI_Wtime();
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
    if ( nprocs > 1 ) {
	printf("ERROR: Illegal number of procs (%d)\n", nprocs );
	//MPI_Abort( MPI_COMM_WORLD, 0 );
        MPI_Finalize();
        return 0;
    }
    // 入力データ、環境変数などの読み込み
    // マスタープロセスだけが、直接、これらのデータにアクセスできると仮定
    // 入力ファイル名の取得
    /*
    if ( (p=getenv("OFMO_INPUT")) == NULL ) {
	dbg("OFMO_INPUT variable is not set\n");
	MPI_Abort( MPI_COMM_WORLD, 2);
    }
    strncpy( input, p, MAXSTRLEN );
    input[MAXSTRLEN-1] = '\0';
    */
    input[0] = '\0';
    if ( (p=getenv("OFMO_INPUT")) != NULL ) strncpy( input, p, MAXSTRLEN );
    dens[0] = '\0';
    if ( (p=getenv("OFMO_FILE_NAME")) != NULL ) strncpy( dens, p, MAXSTRLEN );

    // 起動可能プロセス数の取得
    /*
    if ( (p=getenv("OFMO_NPROCS")) == NULL ) {
	dbg("OFMO_NPROCS variable is not set\n");
	MPI_Abort( MPI_COMM_WORLD, 3 );
    }
    if ( (p=getenv("OFMO_NPROCS")) != NULL ) {
    nmaxprocs = atoi( p );
    if ( nmaxprocs < 1 ) {
	 //dbg("Illegal value of OFMO_NPROCS(%d)\n", nmaxprocs );
	 //MPI_Abort( MPI_COMM_WORLD, 4 );
    }
    */
    if ( (p=getenv("OFMO_NPROCS")) != NULL ) nmaxprocs = atoi( p );

    // options
    {
      char myname[MAXSTRLEN];
      char *endptr;
      extern char *optarg;
      extern int optind, opterr;
      int option_index = 0;
      static struct option long_options[] = {
        {"ngroup",    1, 0, 0},
        {"ng",        1, 0, 0},
        {"niogroup",  1, 0, 1},
        {"nioprocs",  1, 0, 2},
        {"master",    1, 0, 3},
        {"iopos",     1, 0, 4},
        {"nmaxprocs", 1, 0, 'n'},
        {"np",        1, 0, 'n'},
        {"buffer",    1, 0, 'B'},
        {"bindir",    1, 0, 5},
        {"scrdir",    1, 0, 6},
        {"ndev",      1, 0, 'd'},
        {"help",      0, 0, 'h'},
        {"verbose",   0, 0, 'v'},
        {0, 0, 0, 0}
      };

      int c;
      int verbose = 0;
      int help = 0;
      int np = -1;
      ngroup = -1;
      niogroup = -1;
      nioprocs = -1;
      ngroupA = -1;
      niogroupA = -1;
      nioprocsA = -1;
      master = -1;
      iopos = -2;

      eribfsz = -1;
      ndev = 0;

      strncpy(myname, basename(argv[0]), MAXSTRLEN);
      bindir[0]='\0';
      local_input_dir[0]='\0';

      while ((c=getopt_long_only(argc, argv, "B:d:n:hv", long_options, &option_index))!=-1) {
        switch(c) {
          case 0:
            ngroup = strtol(optarg, &endptr, 10);
            if (endptr == optarg) help = -1;
            if (ngroup<=0) help = -1;
            ngroupA = ngroup;
            break;
          case 1:
            niogroup = strtol(optarg, &endptr, 10);
            if (endptr == optarg) help = -1;
            if (niogroup<=0) help = -1;
            niogroupA = niogroup;
            break;
          case 2:
            nioprocs = strtol(optarg, &endptr, 10);
            if (endptr == optarg) help = -1;
            if (nioprocs<=0) help = -1;
            nioprocsA = nioprocs;
            break;
          case 3:
            master = strtol(optarg, &endptr, 10);
            if (endptr == optarg) help = -1;
            if (master<0) help = -1;
            break;
          case 4:
            iopos = strtol(optarg, &endptr, 10);
            if (endptr == optarg) help = -1;
            if (iopos<-1) help = -1;
            break;
          case 5:
            strncpy(bindir, optarg, MAXSTRLEN);
            if (bindir[0]=='\0') help = -1;
            break;
          case 6:
            strncpy(local_input_dir, optarg, MAXSTRLEN);
            if (local_input_dir[0]=='\0') help = -1;
            break;
          case 'n':
            np = strtol(optarg, &endptr, 10);
            if (endptr == optarg) help = -1;
            if (np<=0) help = -1;
            break;
          case 'B':
            eribfsz = strtol(optarg, &endptr, 10);
            if (endptr == optarg) help = -1;
            if (eribfsz<0) help = -1;
            break;
          case 'd':
            ndev = strtol(optarg, &endptr, 10);
#ifndef USE_CUDA
            ndev = 0;
#endif
            if (endptr == optarg) help = -1;
            if (ndev<0) help = -1;
            break;
          case 'v':
            verbose = 1;
            break;
          case 'h':
          default:
            help = 1;
            break;
        }
      }
      argc -= optind;
      argv += optind;

      if (help!=0) {
        show_help(myname, verbose);
        MPI_Finalize();
        return 0;
      }

      if (argc > 0) strncpy(input, argv[0], MAXSTRLEN);
      if (argc > 1) strncpy(dens, argv[1], MAXSTRLEN);
      if ( input[0]=='\0' ) {
        dbg("Input file is not specified.\n");
        show_help(myname, verbose);
        MPI_Finalize();
        return 2;
      }
      input[MAXSTRLEN-1] = '\0';
      dens[MAXSTRLEN-1] = '\0';

      if (np>0) nmaxprocs = np;
      if ( nmaxprocs < 1 ) {
        dbg("Illegal number of procs. (%d)\n", nmaxprocs );
        show_help(myname, verbose);
        MPI_Finalize();
        return 4;
      }
    } // for getopt

#if 0
    // メモリサーバーサービス名のprefix名取得
    if ( (p=getenv("OFMO_PORT_NAME_PREFIX")) == NULL ) {
	dbg("OFMO_PORT_NAME_PREFIX is not set\n");
	MPI_Abort( MPI_COMM_WORLD, 5 );
    }
    strncpy( port_name_prefix, p, MAXSTRLEN);
    port_name_prefix[MAXSTRLEN-1] = '\0';
#endif
    // プロファイルデータファイル名のヘッダ取得
    if ( (p=getenv("OFMO_HEADER")) != NULL ) {
	strncpy( header, p, MAXSTRLEN );
	header[MAXSTRLEN-1] = '\0';
    } else {
	int pos;
	pos = strcspn( input, "." );
	strncpy( header, input, pos );
	header[pos] = '\0';
    }
    // path name to mserv and worker 
    if (bindir[0]=='\0') {
    if ( (p=getenv("OFMO_BINDIR")) != NULL )
	strncpy( bindir, p, MAXSTRLEN );
    else strncpy( bindir, "./", MAXSTRLEN );
    }
    bindir[MAXSTRLEN-1] = '\0';
    // ローカルな入力ファイルの置き場所
    if (local_input_dir[0]=='\0') {
    if ( (p=getenv("OFMO_LOCAL_INPUT_DIR")) != NULL )
	strncpy( local_input_dir, p, MAXSTRLEN );
    else strncpy( local_input_dir, "./", MAXSTRLEN );
    }
    local_input_dir[MAXSTRLEN-1] = '\0';
    // 入力データを読み込む
    if ( ofmo_init( input, MPI_COMM_WORLD ) != 0 ) {
	dbg("failure in input from file (%s)\n", input );
	//MPI_Abort( MPI_COMM_WORLD, 6 );
        MPI_Finalize();
        return 6;
    }
    ierr = ofmo_data_get_vals("ngroup niogroup nioprocs",
	    &ngroup, &niogroup, &nioprocs );
    if ( ierr != 0 ) {
	dbg("failure in get from database\n");
	//MPI_Abort( MPI_COMM_WORLD, 66 );
        MPI_Finalize();
        return 66;
    }
    if (ngroupA>0) ngroup = ngroupA;
    if (niogroupA>0) niogroup = niogroupA;
    if (nioprocsA>0) nioprocs = nioprocsA;
    ierr = ofmo_data_put_vals("ngroup niogroup nioprocs",
	    &ngroup, &niogroup, &nioprocs );
    if (master>=0) ofmo_data_put_vals("master", &master);
    if (iopos>=-1) ofmo_data_put_vals("iopos", &iopos);

    if (eribfsz<0) {
      ierr = ofmo_data_get_vals("nintic", &eribfsz );
      if (ierr!=0) eribfsz = 0;
    }

    // 各ワーカーのプロセス数（グループサイズ）の計算
    group_size = ( nmaxprocs - (nioprocs*niogroup) - 1 ) / ngroup;

    printf("ngroup       = %d\n", ngroup );
    printf("nioprocs     = %d\n", nioprocs );
    printf("niogroup     = %d\n", niogroup );
    printf("nmaxprocs    = %d\n", nmaxprocs );
    printf("group size   = %d\n", group_size );
    printf("# of threads = %d\n", omp_get_max_threads() );
    printf("nintic       = %ld\n", eribfsz );
#ifdef USE_CUDA
    printf("# GPGPU      = %d\n", ndev );
#endif
    printf("input file name       = %s\n", input );
    printf("prefix of output file = %s\n", header );
//    printf("port name prefix      = %s\n", port_name_prefix );
    printf("local input dir.      = %s\n", local_input_dir );
    printf("# of total invoked procs = %d (=%diog*%diop+1+%dgr*%dpr)\n",
	    (nioprocs*niogroup + 1 + ngroup*group_size),
	    niogroup, nioprocs, ngroup, group_size);

    if ( group_size < 1 ) {
	dbg("Illgel group size (%d)\n", group_size );
	//MPI_Abort( MPI_COMM_WORLD, 7 );
        MPI_Finalize();
        return 7;
    }

    // debug
    MPI_Get_processor_name( hostname, &resultlen );
    printf("master host name = %s\n", hostname );
    printf("-----------------------------------------------------\n");

    // FMO計算で用いるデータ構造の決定
    int retval, maxworkers=ngroup, nfrag;
    int ndata=11, nmservs=niogroup, mserv_size=nioprocs, mserv_id;
    int *data_target, *data_offset, *data_nelems;
    char service_name_header[MAXSTRLEN];
    ofmo_data_get_vals("nfrag", &nfrag );
    ofmo_make_data_struct( ndata, nfrag, mserv_size );
    printf(" # of data = %d\n", ofmo_get_ndata() );
    printf("nfrag = %d\n", nfrag );
    et1 = MPI_Wtime();
    printf("##T %38s : lap time= %10.6f ( total etime= %10.6f )\n",
	    "Finish Initialization", (et1-et0), (et1-ET0) );
    fflush(stdout);

    // populationデータ作成のための初期化
    double *atpop_total, *aopop_total, *aopop_frg, *atpop_frg;
    int *fsao2tuao, *fatom2tatom, fnao, fnat, bufsz_aop, bufsz_atp;
    int nao, natom, maxnfao, maxnfatom, nbody;
    ofmo_data_get_vals("nao natom maxnfao maxnfatom nbody",
	    &nao, &natom, &maxnfao, &maxnfatom, &nbody );
    bufsz_aop   = maxnfao  * nbody;
    bufsz_atp   = maxnfatom * nbody;
    atpop_total = (double*)malloc( sizeof(double) * natom );
    aopop_total = (double*)malloc( sizeof(double) * nao );
    aopop_frg   = (double*)malloc( sizeof(double) * bufsz_aop );
    atpop_frg   = (double*)malloc( sizeof(double) * bufsz_atp );
    fsao2tuao   = (int*)malloc( sizeof(int) * bufsz_aop );
    fatom2tatom = (int*)malloc( sizeof(int) * bufsz_atp );
    memset( atpop_total, '\0', sizeof(double)*natom );
    memset( aopop_total, '\0', sizeof(double)*nao );

    // *************** メモリサーバーの起動 *****************
    MPI_Comm *comm_mservs;
    int *used_mservs;
    et0 = MPI_Wtime();
    comm_mservs = (MPI_Comm*)malloc( sizeof(MPI_Comm) * nmservs );
    used_mservs = (int*)malloc( sizeof(int) * nmservs );
    for ( mserv_id=0; mserv_id<nmservs; mserv_id++ )
	used_mservs[mserv_id] = false;
    data_target = ofmo_getadd_data_target( 0 );
    data_nelems = ofmo_getadd_data_nelems( 0 );
    data_offset = ofmo_getadd_data_offset( 0 );
    ofmo_get_service_name_header( service_name_header );
    for ( mserv_id=0; mserv_id<nmservs; mserv_id++ ) {
	retval = ofmo_create_mserv( mserv_size, maxworkers, ndata,
		nfrag, data_nelems, data_target, data_offset,
                bindir,
		service_name_header, mserv_id,
		&comm_mservs[mserv_id] );
	if ( retval < 0 ) {
	    dbg("Failure in ofmo_mserv_create (%d)\n", mserv_id );
	    MPI_Abort( MPI_COMM_WORLD, 11 );
	}
	used_mservs[mserv_id] = true;
    }

    et1 = MPI_Wtime();
    printf("##T %38s : lap time= %10.6f ( total etime= %10.6f )\n",
	    "Finish to Create Memory Servers", (et1-et0), (et1-ET0) );
    fflush(stdout);

    // ******************** workerの起動 *********************
    int workerid, nworkers=ngroup, worker_size=group_size;
    int comped_data_size, max_comped_data_size;
    int limit, ndiv, nmod, def_mserv_id;
    int *def_mserv_ids;
    unsigned char *comped_data;
    MPI_Comm *comm_workers;
    et0 = MPI_Wtime();
    comm_workers = (MPI_Comm*)malloc( sizeof(MPI_Comm) * nworkers );
    def_mserv_ids = (int*)malloc( sizeof(int) * nworkers );
    comped_data = ofmo_alloc_compressed_data( &max_comped_data_size );
    comped_data_size = ofmo_compress_data( input, comped_data,
	    max_comped_data_size );
    ndiv = nworkers/nmservs;
    nmod = nworkers%nmservs;
    limit =  nmod * (ndiv+1);
    if ( comped_data_size < 1 ) {
	dbg("Illegal size of compressed data(%d)\n", comped_data_size);
	MPI_Abort( MPI_COMM_WORLD, 20 );
    }
    for ( workerid=0; workerid<nworkers; workerid++ ) {
	def_mserv_id =
	    ( workerid<limit ? workerid/(ndiv+1) : (workerid-limit)/ndiv );
	ierr = ofmo_create_worker( worker_size, workerid,
            eribfsz, ndev,
		comped_data_size, comped_data,
		input, header, local_input_dir,
		nmservs, ndata, nfrag, data_target,
                bindir,
		service_name_header, def_mserv_id,
		&comm_workers[workerid] );
	def_mserv_ids[workerid] = def_mserv_id;
	if ( ierr < 0 ) {
	    dbg("Failure in ofmo_worker_create(%d)\n", workerid );
	    MPI_Abort( MPI_COMM_WORLD, 21 );
	}
    }
    // debug
    et1 = MPI_Wtime();
    printf("##T %38s : lap time= %10.6f ( total etime= %10.6f )\n",
	    "Finish to Create Workers", (et1-et0), (et1-ET0) );
    fflush(stdout);

    // *************** メモリサーバーとワーカーとの接続要求 *************
    et0 = MPI_Wtime();
    imsg[OFMO_I_CMD] = OFMO_ACCEPT;
    for ( workerid=0; workerid<nworkers; workerid++ ) {
	for ( mserv_id=0; mserv_id<nmservs; mserv_id++ ) {
	    imsg[OFMO_I_METHOD] = workerid;
	    ofmo_call_mserv( comm_mservs[mserv_id], imsg_sz, imsg );
	    imsg[OFMO_I_METHOD] = mserv_id;
	    ofmo_call_worker( comm_workers[workerid], imsg_sz, imsg,
		    dmsg_sz, dmsg );
	}
    }

    et1 = MPI_Wtime();
    printf("##T %38s : lap time= %10.6f ( total etime= %10.6f )\n",
	    "Finish to Connect Mservs & Workers", (et1-et0), (et1-ET0) );
    fflush(stdout);
    // ワーカーのバリア
    ofmo_worker_barrier( nworkers, comm_workers );

    int ifrag;
    double **dmsgs;
    MPI_Request *req_workers;
    int **issued_job_list;

    /* ワーカーに指示したジョブのリスト */
    issued_job_list = ofmo_alloc_imatrix( nworkers, (MAXNJOB*2) );
    for ( workerid=0; workerid<nworkers; workerid++ ) {
	for ( int i=0; i<(MAXNJOB*2); i++ )
	    issued_job_list[workerid][i]=-1;
    }
    ofmo_init_frag_order();	// モノマーをAO数の大きい順にソート

    req_workers = (MPI_Request*)malloc( sizeof(MPI_Request) * ngroup );
    dmsgs = ofmo_alloc_dmatrix( ngroup, OFMO_DMSG_SZ );

    // ******************* 初期密度行列計算 *********************
    int calc_init_dens = true;
    et0 = MPI_Wtime();
    // ファイルからデータを読み込む
    //if ( (p=getenv("OFMO_FILE_NAME")) != NULL ) {
    if ( dens[0] != '\0' ) {
        p = dens;
	if ( ofmo_read_and_put_all( nmservs, comm_mservs, p ) == 0 ) {
	    calc_init_dens = false;
	    ofmo_reset_monomer_data_worker( nworkers, comm_workers );
	}
    }
    if ( calc_init_dens ) {
	ifrag = 0;
	imsg[OFMO_I_CMD]  = OFMO_INIT_DENS;
	imsg[OFMO_I_NMON] = 1;
	// まずはじめに、すべてのワーカーに一通り、ジョブを投げる
	for ( workerid=0; workerid<nworkers; workerid++ ) {
	    imsg[OFMO_I_MON1] = ifrag;
	    req_workers[workerid] =
		ofmo_call_async_worker( comm_workers[workerid],
			OFMO_IMSG_SZ, imsg, OFMO_DMSG_SZ, dmsgs[workerid] );
	    ifrag += group_size;
	}
	while (1) {
	    if ( ifrag >= nfrag ) break;
	    workerid = ofmo_waitany_worker( nworkers, req_workers );
	    imsg[OFMO_I_MON1] = ifrag;
	    req_workers[workerid] =
		ofmo_call_async_worker( comm_workers[workerid],
			OFMO_IMSG_SZ, imsg, OFMO_DMSG_SZ, dmsgs[workerid] );
	    ifrag += group_size;
	}
	// 最初に戻ってきた１つのワーカーに対して、モノマー間距離の
	// 計算を依頼する
	workerid = ofmo_waitany_worker( nworkers, req_workers );
	imsg[OFMO_I_CMD] = OFMO_DISTANCE;
	req_workers[workerid] =
	    ofmo_call_async_worker( comm_workers[workerid],
		    OFMO_IMSG_SZ, imsg, OFMO_DMSG_SZ, dmsgs[workerid] );
	printf("calc. inter-frag distance list = %d\n", workerid );
	fflush( stdout );
	// すべてのワーカープロセスでの密度行列計算が終了するのを待つ
	ofmo_waitall_worker( nworkers, req_workers );
	// ここでバリア同期をとる
	ofmo_worker_barrier( nworkers, comm_workers );
	// モノマーデータの交換（更新）を行う
	//ofmo_mserv_barrier( nmservs, comm_mservs );
	ofmo_update_monomer_data_worker( nworkers, comm_workers );
	ofmo_update_monomer_data_master();
    }
    et1 = MPI_Wtime();
    printf("##T %38s : lap time= %10.6f ( total etime= %10.6f )\n",
	    "Finish to Calculate Initial Density", (et1-et0), (et1-ET0) );
    fflush(stdout);
    if ( nbody == 0 ) goto final;
    // ============= モノマーSCC計算 =================
    double etscc0, etscc1;
    int itera, maxscc, sccconv, scfconv_now;
    double total_energy, total_energy0;
    double total_energy0_prev=0.e0, dE=HUGE_VAL;
    double scc_tol;
    imsg[OFMO_I_CMD]    = OFMO_SCF;
    imsg[OFMO_I_METHOD] = OFMO_RHF;
    imsg[OFMO_I_NMON]   = 1;
    ofmo_data_get_vals("maxscc sccconv", &maxscc, &sccconv );
    scc_tol = pow( 0.1e0, (double)sccconv );
    if ( calc_init_dens == false ) {	// ファイルから読み込んだ場合
	double *e0s;
	e0s = (double*)malloc( sizeof(double) * nfrag );
	ofmo_master_get( comm_mservs[0], OFMO_ENERGY0, -1, e0s );
	total_energy0_prev=0.e0;
	for ( ifrag=0; ifrag<nfrag; ifrag++ ) total_energy0_prev += e0s[ifrag];
	free( e0s );
    }
    printf("---- start SCC (maxscc=%d, sccconv=%d) ----\n",
	    maxscc, sccconv );
    fflush(stdout);
    etscc0 = MPI_Wtime();
    for ( itera=1; itera<=maxscc; itera++ ) {
	et0 = MPI_Wtime();
	imsg[OFMO_I_SCC]  = itera;
	imsg[OFMO_I_CONV] = scfconv_now = ofmo_get_monomer_scfconv( dE );
	total_energy = total_energy0 = 0.e0;
	// すべてのモノマー計算を各グループに割り振る
	for ( int i=0; i<nfrag; i++ ) {
	    imsg[OFMO_I_MON1] = frag_order[i];
	    //imsg[OFMO_I_MON1] = i;
	    if ( i < nworkers ) {
		workerid = i;
	    } else {
		workerid = ofmo_waitany_worker( nworkers, req_workers );
		total_energy  += dmsgs[workerid][OFMO_D_ENERGY];
		total_energy0 += dmsgs[workerid][OFMO_D_ENERGY0];
		issued_job_list[workerid][0] = -1;
	    }
	    req_workers[workerid] = ofmo_call_async_worker(
		    comm_workers[workerid],
		    imsg_sz, imsg, dmsg_sz, dmsgs[workerid] );
	    issued_job_list[workerid][0] = imsg[OFMO_I_MON1];
	}
	// すべてのワーカーの終了を待つ
	ofmo_waitall_worker( nworkers, req_workers );

	for ( workerid=0; workerid<nworkers; workerid++ ) {
	    total_energy  += dmsgs[workerid][OFMO_D_ENERGY];
	    total_energy0 += dmsgs[workerid][OFMO_D_ENERGY0];
	    issued_job_list[workerid][0] = -1;
	}
	et1 = MPI_Wtime();
	printf("##SCC itera=%2d  conv=%d ", itera, scfconv_now );
	printf(" energy = %16.8f  energy0 = %16.8f  ( %10.6f )\n",
		total_energy, total_energy0, (et1-et0) );
	fflush( stdout );
	// モノマーデータの交換（更新）を行う
	ofmo_worker_barrier( nworkers, comm_workers );
	//ofmo_mserv_barrier( nmservs, comm_mservs );
	ofmo_update_monomer_data_worker( nworkers, comm_workers );
	ofmo_update_monomer_data_master();
	//
	dE = total_energy0 - total_energy0_prev;
	if ( fabs(dE) < scc_tol ) {
	    printf("==== SCC converged ====\n");
	    break;
	}
	total_energy0_prev = total_energy0;
    }
    etscc1 = MPI_Wtime();
    printf("##T %38s : lap time= %10.6f ( total etime= %10.6f )\n",
	    "---- Finish SCC procedure ----", (etscc1-etscc0),
	    (etscc1-ET0) );
    fflush(stdout);
    // === モノマーのデータ（密度行列、AO populationなど）をファイルに出力
    if ( (p=getenv("OFMO_FILE_NAME")) != NULL ) {
	ofmo_get_and_write( comm_mservs[0], p );
    }
    // debug
    //goto final;

    /*// テスト
    // すべてのメモリサーバーを一旦killして、再起動する
    // そのあと、ファイルからデータを読み込む
    // まず、ワーカーがメモリサーバーとの接続を解除する
    ofmo_worker_barrier( nworkers, comm_workers );
    // すべてのメモリサーバーをkillする
    for ( mserv_id=0; mserv_id<nmservs; mserv_id++ ) {
	ofmo_destroy_mserv( comm_mservs[mserv_id] );
    }
    printf("-- destroyed all mservs --\n");
    fflush( stdout );
    // 新規にメモリサーバーを起動する
    et0 = MPI_Wtime();
    data_target = ofmo_getadd_data_target( 0 );
    data_nelems = ofmo_getadd_data_nelems( 0 );
    data_offset = ofmo_getadd_data_offset( 0 );
    ofmo_get_service_name_header( service_name_header );
    for ( mserv_id=0; mserv_id<nmservs; mserv_id++ ) {
	retval = ofmo_create_mserv( mserv_size, maxworkers, ndata,
		nfrag, data_nelems, data_target, data_offset,
                bindir,
		service_name_header, mserv_id,
		&comm_mservs[mserv_id] );
	if ( retval < 0 ) {
	    dbg("Failure in ofmo_mserv_create (%d)\n", mserv_id );
	    MPI_Abort( MPI_COMM_WORLD, 11 );
	}
    }
    et1 = MPI_Wtime();
    printf("##T %38s : lap time= %10.6f ( total etime= %10.6f )\n",
	    "Finish to Create Memory Servers2", (et1-et0), (et1-ET0) );
    fflush(stdout);
    // 起動されたメモリサーバーにファイルから読み込んだデータを送る
    if ( (p=getenv("OFMO_FILE_NAME")) != NULL ) {
	ofmo_read_and_put_all( nmservs, comm_mservs, p );
    }
    // ワーカーとメモリサーバーを接続する
    et0 = MPI_Wtime();
    imsg[OFMO_I_CMD] = OFMO_ACCEPT;
    for ( workerid=0; workerid<nworkers; workerid++ ) {
	for ( mserv_id=0; mserv_id<nmservs; mserv_id++ ) {
	    imsg[OFMO_I_METHOD] = workerid;
	    ofmo_call_mserv( comm_mservs[mserv_id], imsg_sz, imsg );
	    imsg[OFMO_I_METHOD] = mserv_id;
	    ofmo_call_worker( comm_workers[workerid], imsg_sz, imsg,
		    dmsg_sz, dmsg );
	}
    }
    et1 = MPI_Wtime();
    printf("##T %38s : lap time= %10.6f ( total etime= %10.6f )\n",
	    "Finish to Connect Mservs & Workers", (et1-et0), (et1-ET0) );
    fflush(stdout);
    //ofmo_worker_barrier( nworkers, comm_workers );
    ofmo_reset_monomer_data_worker( nworkers, comm_workers );*/
    /*// テスト２
    // メモリサーバーを１つ落として、それを参照しているワーカーの
    // 参照するメモリサーバーを変更する
    if ( nmservs > 1 ) {
	ofmo_destroy_mserv( comm_mservs[nmservs-1] );
	used_mservs[nmservs-1] = false;
    }
    printf("memory server (id=%d) has been destroyed\n",
	    (nmservs-1) );
    fflush(stdout);
    // ワーカーのデフォルトメモリサーバーを変更
    et0 = MPI_Wtime();
    int def_mserv = 0;
    for ( workerid=0; workerid<nworkers; workerid++ ) {
	imsg[OFMO_I_CMD]    = OFMO_DISCONNECT;
	imsg[OFMO_I_METHOD] = nmservs-1;
	ofmo_call_worker( comm_workers[workerid], imsg_sz, imsg,
		dmsg_sz, dmsg );
	if ( def_mserv_ids[workerid] == (nmservs-1) ) {
	    imsg[OFMO_I_CMD] = OFMO_CHANGE_MSERV;
	    imsg[OFMO_I_METHOD] = def_mserv;
	    ofmo_call_worker( comm_workers[workerid], imsg_sz, imsg,
		    dmsg_sz, dmsg );
	    def_mserv_ids[workerid] = def_mserv;
	}
    }
    et1 = MPI_Wtime();
    printf("##T %38s : lap time= %10.6f ( total etime= %10.6f )\n",
	    "Finish to change default mserv", (et1-et0), (et1-ET0) );
    fflush(stdout);*/
    /*// テスト３
    // 1. メモリサーバを複数起動してSCCを行う
    // 2. SCC終了後、１つのメモリサーバを落とす
    // 3. １つのメモリサーバを新規に起動する
    // 4. 新規メモリサーバにデータをコピーする
    // まず、１つのメモリサーバーをkillする
    if ( nmservs > 1 ) {
	ofmo_destroy_mserv( comm_mservs[nmservs-1] );
	used_mservs[nmservs-1] = false;
    }
    for ( workerid=0; workerid<nworkers; workerid++ ) {
	imsg[OFMO_I_CMD]    = OFMO_DISCONNECT;
	imsg[OFMO_I_METHOD] = nmservs-1;
	ofmo_call_worker( comm_workers[workerid], imsg_sz, imsg,
		dmsg_sz, dmsg );
    }
    printf("---- finish to destroy a mserv ----\n");
    fflush( stdout );
    // 新しいメモリサーバーを起動する
    data_target = ofmo_getadd_data_target( 0 );
    data_nelems = ofmo_getadd_data_nelems( 0 );
    data_offset = ofmo_getadd_data_offset( 0 );
    ofmo_get_service_name_header( service_name_header );
    mserv_id = nmservs-1;
    retval = ofmo_create_mserv( mserv_size, maxworkers, ndata,
	    nfrag, data_nelems, data_target, data_offset,
                bindir,
	    service_name_header, mserv_id,
	    &comm_mservs[mserv_id] );
    if ( retval < 0 ) {
	dbg("Failure in ofmo_mserv_create (%d)\n", mserv_id );
	MPI_Abort( MPI_COMM_WORLD, 11 );
    }
    used_mservs[mserv_id] = true;
    printf("----finish to create a mserv ----\n");
    fflush( stdout );
    // 作成したメモリサーバーにデータをコピーする
    ofmo_copy_mserv_data( comm_mservs[0], comm_mservs[mserv_id] );
    printf("---- finish to copy mserv data ----\n");
    fflush(stdout);
    // コネクションを張る
    imsg[OFMO_I_CMD] = OFMO_ACCEPT;
    for ( workerid=0; workerid<nworkers; workerid++ ) {
	imsg[OFMO_I_METHOD] = workerid;
	ofmo_call_mserv( comm_mservs[mserv_id], imsg_sz, imsg );
	imsg[OFMO_I_METHOD] = mserv_id;
	ofmo_call_worker( comm_workers[workerid], imsg_sz, imsg,
		dmsg_sz, dmsg );
    }
    ofmo_worker_barrier( nworkers, comm_workers );
    printf("---- finish to connect to new mserv ----\n");*/
    /*// テスト４：１つのメモリサーバーをkillしたあと、
    // ワーカーに何も知らせなくても、実行が継続できるか？
    if ( nmservs > 1 ) {
	ofmo_destroy_mserv( comm_mservs[nmservs-1] );
	used_mservs[nmservs-1] = false;
    }
    printf("memory server (id=%d) has been destroyed\n",
	    (nmservs-1) );
    fflush(stdout);*/
    // =================== dimerの計算 ====================
    if ( nbody >= 2 ) {
	int jfrag, nscf_dimer=0;
	double de0scf;
	double *dist, ldim;
	double *menergy0;
	int scfconv;
	et0 = MPI_Wtime();
	ofmo_data_get_vals("ldim scfconv", &ldim, &scfconv );
	dist = (double*)malloc( sizeof(double) * nfrag );
	menergy0 = (double*)malloc( sizeof(double) * nfrag );
	ofmo_master_get( comm_mservs[0], OFMO_ENERGY0, -1, menergy0 );
	imsg[OFMO_I_CMD]  = OFMO_SCF;
	imsg[OFMO_I_METHOD] = OFMO_RHF;
	imsg[OFMO_I_SCC]  = maxscc + 100;
	imsg[OFMO_I_CONV] = scfconv;
	imsg[OFMO_I_NMON] = 2;
	de0scf = 0.e0;
	for ( int i=1; i<nfrag; i++ ) {
	    ifrag = frag_order[i];
	    imsg[OFMO_I_MON1] = ifrag;
	    ofmo_master_get( comm_mservs[0], OFMO_DISTA, ifrag, dist );
	    for ( int j=0; j<i; j++ ) {
		jfrag = frag_order[j];
		if ( dist[jfrag] < ldim ) {
		    imsg[OFMO_I_MON2] = jfrag;
		    if ( nscf_dimer < ngroup ) {
			workerid = nscf_dimer;
		    } else {
			workerid = ofmo_waitany_worker( nworkers,
				req_workers );
			de0scf +=
			    ( dmsgs[workerid][OFMO_D_ENERGY0]
			     -menergy0[issued_job_list[workerid][0]]
			     -menergy0[issued_job_list[workerid][1]] );
			de0scf += dmsgs[workerid][OFMO_D_DDV];
			ofmo_get_pop_data( comm_workers[workerid], OFMO_BOTH_POPS,
				bufsz_aop, &fnao, aopop_frg, fsao2tuao,
				bufsz_atp, &fnat, atpop_frg, fatom2tatom );
			ofmo_accumulate( fnao, aopop_frg, fsao2tuao, aopop_total );
			ofmo_accumulate( fnat, atpop_frg, fatom2tatom, atpop_total );
			issued_job_list[workerid][0] = -1;
			issued_job_list[workerid][1] = -1;
		    }
		    req_workers[workerid] = ofmo_call_async_worker(
			    comm_workers[workerid],
			    imsg_sz, imsg, dmsg_sz, dmsgs[workerid] );
		    issued_job_list[workerid][0] = ifrag;
		    issued_job_list[workerid][1] = jfrag;
		    nscf_dimer++;
		}
	    }
	}
	// SCF dimer計算の数が少ない場合の処理
	if ( nscf_dimer < nworkers ) {
	    imsg[OFMO_I_NMON] = 0;
	    for ( workerid=nscf_dimer; workerid<nworkers; workerid++ ) {
		req_workers[workerid] = ofmo_call_async_worker(
			comm_workers[workerid],
			imsg_sz, imsg, dmsg_sz, dmsgs[workerid] );
	    }
	}
	// すべてのワーカーの終了を待つ
	ofmo_waitall_worker( nworkers, req_workers );
	for ( workerid=0; workerid<nworkers; workerid++ ) {
	    de0scf +=
		( dmsgs[workerid][OFMO_D_ENERGY0]
		 -menergy0[issued_job_list[workerid][0]]
		 -menergy0[issued_job_list[workerid][1]] );
	    de0scf += dmsgs[workerid][OFMO_D_DDV];
	    ofmo_get_pop_data( comm_workers[workerid], OFMO_BOTH_POPS,
		    bufsz_aop, &fnao, aopop_frg, fsao2tuao,
		    bufsz_atp, &fnat, atpop_frg, fatom2tatom );
	    ofmo_accumulate( fnao, aopop_frg, fsao2tuao, aopop_total );
	    ofmo_accumulate( fnat, atpop_frg, fatom2tatom, atpop_total );
	    issued_job_list[workerid][0] = -1;
	    issued_job_list[workerid][1] = -1;
	}
	free( dist );
	et1 = MPI_Wtime();
	printf("##T %38s : lap time= %10.6f ( total etime= %10.6f )\n",
		"---- Finish to Calc. Dimer SCF ----",
		(et1-et0), (et1-ET0) );
	printf( "       ( # of scf dimer = %d )\n", nscf_dimer );
	printf("        dE0(SCF dimer) = %16.8f\n", de0scf );
	fflush(stdout);
	int maxnjob=MAXNJOB, nline=NLINE, nes_dimer=0;
	int njob, *joblist, *ifrags;
	double **dists;
	double total_esdimer0;
	int nsend_esjob=0;
	et0 = MPI_Wtime();
	joblist = &imsg[OFMO_I_MON1];
	dists = ofmo_alloc_dmatrix( nline, nfrag );
	ifrags  = (int*)malloc( sizeof(int) * nline );
	imsg[OFMO_I_CMD] = OFMO_APPROX;
	njob = 0;
	total_esdimer0=0.e0;
	for ( int ii=0; ii<nfrag; ii+=nline ) {
	    int ie;
	    ie = ii + nline;
	    if ( ie > nfrag ) ie = nfrag;
	    // フラグメント番号と距離の保存
	    for ( int i=ii, k=0; i<ie; i++, k++ ) {
		ifrags[k] = frag_order[i];
		ofmo_master_get( comm_mservs[0], OFMO_DISTA,
			ifrags[k], dists[k] );
	    }
	    for ( int j=0; j<ie; j++ ) {
		jfrag = frag_order[j];
		for ( int i=ii, k=0; i<ie; i++, k++ ) {
		    if ( j>=i ) continue;
		    if ( dists[k][jfrag] >= ldim ) {
			if ( ifrags[k] > jfrag ) {
			    joblist[njob*2+0] = ifrags[k];
			    joblist[njob*2+1] = jfrag;
			} else {
			    joblist[njob*2+0] = jfrag;
			    joblist[njob*2+1] = ifrags[k];
			}
			njob++;
			nes_dimer++;
			if ( njob == maxnjob ) {
			    if ( nsend_esjob<nworkers ) {
				workerid = nsend_esjob;
			    } else {
				workerid = ofmo_waitany_worker(nworkers,
					req_workers);
				total_esdimer0 +=
				    dmsgs[workerid][OFMO_D_ENERGY0];
			    }
			    imsg[OFMO_I_NMON] = njob;
			    memcpy( issued_job_list[workerid], joblist,
				    sizeof(int)*njob*2 );
			    req_workers[workerid] = ofmo_call_async_worker(
				    comm_workers[workerid],
				    imsg_sz, imsg, dmsg_sz,
				    dmsgs[workerid] );
			    njob=0;
			    nsend_esjob++;
			}	// if ( njob > (maxnjob-nline) )
		    }	// if ( dists[k][jfrag] >= ldim )
		}	// for ( int i=ii)
	    }	// for ( int j )
	}	// for ( int ii );
	// バッファに残った分の転送
	if ( njob != 0 ) {
	    if ( nsend_esjob<nworkers ) {
		workerid = nsend_esjob;
	    } else {
		workerid = ofmo_waitany_worker( nworkers, req_workers );
		total_esdimer0 += dmsgs[workerid][OFMO_D_ENERGY0];
	    }
	    imsg[OFMO_I_NMON] = njob;
	    req_workers[workerid] = ofmo_call_async_worker(
		    comm_workers[workerid],
		    imsg_sz, imsg, dmsg_sz, dmsgs[workerid] );
	    njob=0;
	    nsend_esjob++;
	}
	// ES dimerのジョブを発行していないワーカーが存在するとき
	if ( nsend_esjob < nworkers ) {
	    imsg[OFMO_I_NMON] = 0;
	    for ( workerid=nsend_esjob; workerid<nworkers; workerid++ ) {
		req_workers[workerid] = ofmo_call_async_worker(
			comm_workers[workerid],
			imsg_sz, imsg, dmsg_sz, dmsgs[workerid] );
	    }
	}

	ofmo_waitall_worker( nworkers, req_workers );
	for ( workerid=0; workerid<nworkers; workerid++ ) {
	    total_esdimer0 += dmsgs[workerid][OFMO_D_ENERGY0];
	}

	double total_fmo_energy;
	et1 = MPI_Wtime();
	printf("##T %38s : lap time= %10.6f ( total etime= %10.6f )\n",
		"---- Finish to Calc. ES Dimer ----", (et1-et0),
		(et1-ET0) );
	printf("    (# of ES dimer = %d)\n", nes_dimer );
	printf("    dE0(ES dimer) = %16.8f\n", total_esdimer0 );
	total_fmo_energy = total_energy0 + de0scf + total_esdimer0;
	printf("    total energy(FMO) = %16.8f\n", total_fmo_energy );
	fflush(stdout);

	ofmo_free_dmatrix( dists );
	free ( ifrags );
    }

final:

    // workerの終了
    et0 = MPI_Wtime();
    for ( workerid=0; workerid<nworkers; workerid++ ) {
	ofmo_destroy_worker( comm_workers[workerid] );
	if ( req_workers[workerid] != MPI_REQUEST_NULL ) {
	    MPI_Cancel( &req_workers[workerid] );
	    MPI_Request_free( &req_workers[workerid] );
	}
    }
    et1 = MPI_Wtime();
    printf("##T %38s : lap time= %10.6f ( total etime= %10.6f )\n",
	    "-- All workers have been destroyed --", (et1-et0), (et1-ET0) );
    fflush(stdout);

    if ( nbody > 0 ) {
	// 全体のAO population、ならびに、atomic populationの計算
	int *atomic_number, **msao2tuao, *nfao, *nfatom, **ifatom;
	int *atm_lcs, *ushel_lqn, *ushel_ini;
	int NNAO[] = { 1, 3, 6, 10, 15 };
	char *CTP = "SPDFG";
	double charge;
	et0 = MPI_Wtime();
	ofmo_data_get_vals(
		"msao2tuao atn atm_lcs ushel_lqn ushel_ini nfao"
		" nfatom ifatom",
		&msao2tuao, &atomic_number, &atm_lcs, &ushel_lqn,
		&ushel_ini, &nfao, &nfatom, &ifatom );
	for ( ifrag=0; ifrag<nfrag; ifrag++ ) {
	    ofmo_get_monomer_aopop_master(comm_mservs[0],ifrag,aopop_frg);
	    ofmo_get_monomer_atpop_master(comm_mservs[0],ifrag,atpop_frg);
	    for ( int iao=0; iao<nfao[ifrag]; iao++ )
		aopop_total[ msao2tuao[ifrag][iao] ] += aopop_frg[iao];
	    for ( int iat=0; iat<nfatom[ifrag]; iat++ )
		atpop_total[ ifatom[ifrag][iat] ] += atpop_frg[iat];
	}
	printf("---- ATOMIC POPULATION DATA -----\n");
	printf("   SN AN  ATOMIC POP  NET CHARGE\n");
	printf("---------------------------------\n");
	charge=0.e0;
	for ( int iat=0; iat<natom; iat++ ) {
	    charge += ((double)atomic_number[iat]-atpop_total[iat]);
	    printf(" %4d %2d  %10.3f  %10.3f\n",
		    (iat+1), atomic_number[iat], atpop_total[iat],
		    ((double)atomic_number[iat]-atpop_total[iat]) );
	}
	printf("---------------------------------\n");
	printf("total charge = %10.3f\n", charge );
	/*printf("---------------------------------\n");
	printf(" ATSN   ICS    AO C   AO POP\n");
	printf("---------------------------------\n");
	int ics0, ics1, ics, iao0, iao1, iao, lqn;
	for ( int iat=0; iat<natom; iat++ ) {
	    ics0 = atm_lcs[iat];
	    ics1 = atm_lcs[iat+1];
	    for ( ics=ics0; ics<ics1; ics++ ) {
		lqn  = ushel_lqn[ics];
		iao0 = ushel_ini[ics];
		iao1 = iao0 + NNAO[lqn];
		for ( iao=iao0; iao<iao1; iao++ )
		    printf(" %4d %5d %5d %c %10.3f\n",
			    (iat+1), (ics+1), (iao+1),
			    CTP[lqn], aopop_tot[iao] );
	    }
	}*/
	printf("---------------------------------\n");
	et1 = MPI_Wtime();
	printf("##T %38s : lap time= %10.6f ( total etime= %10.6f )\n",
		"-- calc AO and atomic population --",
		(et1-et0), (et1-ET0) );
	fflush( stdout );
    }

    // memory serverの終了処理
    et0 = MPI_Wtime();
    for ( mserv_id=0; mserv_id<nmservs; mserv_id++ ) {
	if ( used_mservs[mserv_id] )
	    ofmo_destroy_mserv( comm_mservs[mserv_id] );
    }
    et1 = MPI_Wtime();
    printf("##T %38s : lap time= %10.6f ( total etime= %10.6f )\n",
	    "-- All mservs have been destroyed --", (et1-et0), (et1-ET0) );
    fflush(stdout);

    free( def_mserv_ids );
    free( comm_mservs );
    free( used_mservs );
    free( aopop_total );
    free( atpop_total );
    free( aopop_frg );
    free( atpop_frg );
    free( fsao2tuao );
    free( fatom2tatom );
    
    MPI_Finalize();
    return 0;
}
