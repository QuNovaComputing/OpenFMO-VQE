#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <math.h>

#include <mpi.h>

#include "ofmo-def.h"
#include "ofmo-data.h"
#include "ofmo-prof.h"
#include "ofmo-debug.h"
#include "ofmo-misc.h"

#include "ofmo-integ.h"
#include "ofmo-scf.h"

#include "ofmo-mserv-cont.h"
#include "ofmo-monomer-data.h"

#ifdef USE_CUDA
#include "cuda/cuda-drv.h"
#include "cuda/cuda-integ.h"
#endif

#define MAXTOKLEN MAXSTRLEN

#ifdef DEBUG_MODE
FILE *fp_debug;
#endif

/* 外部関数のプロトタイプ宣言 */
extern int ofmo_decompress_data( const char* filename,
	unsigned char* comped_data, int comped_data_size );

extern int ofmo_init( const char *filename, MPI_Comm comm );

// ofmo-calc-frag.cで定義
extern int ofmo_frag_init();
extern int ofmo_calc_fragment_electronic_state(
	MPI_Comm comm, int nmonomer, int monomer_list[], int level,
	double tolscf,
	double *energy, double *energy0, double *ddv,
	int *fnao, double daopop[], int fsau2tuao[],
	int *fnat, double datpop[], int fatom2tatom[] );
extern void ofmo_set_new_scc_flag();
extern int ofmo_monomer_init_density( const int *imsg, MPI_Comm comm );
extern int ofmo_calc_es_dimer( MPI_Comm comm, int njob, int joblist[],
	double *energy0 );

extern int ofmo_make_inter_frag_distance_list( MPI_Comm comm );

//extern int ofmo_twoint_init();

/* */
#ifdef FJ_MAIN
int MAIN__( int argc, char *argv[] ) {
#else
int main( int argc, char *argv[] ) {
#endif
    int required=MPI_THREAD_SERIALIZED, provided;
    int myrank, nprocs, is_root, root=0, master_root=0;
    int comped_data_size, ndata, nfrag, def_mserv_id;
    char prof_name[MAXSTRLEN];
    char input[MAXSTRLEN], header[MAXSTRLEN];
    char local_input_dir[MAXSTRLEN], service_name_header[MAXSTRLEN];
    int ndev;
    long eribfsz;
    //
    MPI_Comm comm_master;
    int workerid, ierr, nmservs;
    // debug
    int resultlen;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    int rc;

    MPI_Init_thread( &argc, &argv, required, &provided );
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
    is_root = ( myrank == root );
    if ( required != provided ) {
	if ( is_root ) {
	    dbg("MPI_THREAD_SERIALIZED is not supported\n");
	    MPI_Abort( MPI_COMM_WORLD, 10001 );
	}
    }
    MPI_Comm_get_parent( &comm_master );
    // ---- プログラム引数を全プロセスに通知 ----
    if ( is_root ) {
	workerid         = atoi( argv[1] );
	nmservs          = atoi( argv[2] );
	comped_data_size = atoi( argv[3] );
	ndata            = atoi( argv[4] );
	nfrag            = atoi( argv[5] );
	def_mserv_id     = atoi( argv[6] );
	strcpy( input,               argv[ 7] );
	strcpy( header,              argv[ 8] );
	strcpy( local_input_dir,     argv[ 9] );
	strcpy( service_name_header, argv[10] );
        eribfsz            = atol( argv[11] );
        ndev     = atoi( argv[12] );
    }
    rc = MPI_Bcast( &workerid, 1, MPI_INT, root, MPI_COMM_WORLD );
    if (rc!=MPI_SUCCESS) { MPI_Abort( MPI_COMM_WORLD, 10099 ); }
    MPI_Bcast( &nmservs, 1, MPI_INT, root, MPI_COMM_WORLD );
    MPI_Bcast( &comped_data_size, 1, MPI_INT, root, MPI_COMM_WORLD );
    MPI_Bcast( &ndata, 1, MPI_INT, root, MPI_COMM_WORLD );
    MPI_Bcast( &nfrag, 1, MPI_INT, root, MPI_COMM_WORLD );
    MPI_Bcast( &def_mserv_id, 1, MPI_INT, root, MPI_COMM_WORLD );
    MPI_Bcast( input, MAXSTRLEN, MPI_CHAR, root, MPI_COMM_WORLD );
    MPI_Bcast( header, MAXSTRLEN, MPI_CHAR, root, MPI_COMM_WORLD );
    MPI_Bcast( local_input_dir, MAXSTRLEN, MPI_CHAR, root, MPI_COMM_WORLD );
    MPI_Bcast( service_name_header, MAXSTRLEN, MPI_CHAR, root,
	    MPI_COMM_WORLD );
    MPI_Bcast( &eribfsz, 1, MPI_LONG, root, MPI_COMM_WORLD );
    MPI_Bcast( &ndev, 1, MPI_INT, root, MPI_COMM_WORLD );
    //WORKERID = workerid;

    // debug
    if ( is_root ) {
	MPI_Get_processor_name( hostname, &resultlen );
	printf("****** LIST OF HOST NAME (WORKER(ID=%d)) ******\n", workerid);
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

    // ---- マスターから各種データを取得 ----
    unsigned char *comped_data;
    int **data_target;
    comped_data = (unsigned char*)malloc( comped_data_size );
    data_target = ofmo_alloc_imatrix( ndata, nfrag );
    MPI_Bcast( comped_data, comped_data_size, MPI_UNSIGNED_CHAR,
	    master_root, comm_master );
    MPI_Bcast( data_target[0], ndata*nfrag, MPI_INT,
	    master_root, comm_master );
    // ** ここまでで、初期化に伴うマスターからのデータ受信が終了

#ifdef DEBUG_MODE
    {
	char filename[MAXSTRLEN];
	sprintf( filename, "ooo-%d-%d.log", workerid, myrank );
	if ( (fp_debug=fopen( filename, "w")) == NULL ) {
	    return -1;
	}
    }
#endif
    // ---- プロファイルデータを出力するための準備 ----
    ofmo_create_worker_prof_name( input, header, workerid, prof_name,
	    MAXSTRLEN );
    ofmo_prof_init( prof_name, MPI_COMM_WORLD );

    // メモリサーバーへの接続初期化
    ofmo_mserv_container_init( def_mserv_id, ndata, nfrag, data_target,
	    nmservs, service_name_header );

    // ---- 入力ファイルからのデータ読み込み ----
    // 圧縮ファイルを展開
    char filename[MAXSTRLEN];
    sprintf( filename,
	    "%s/%s-%d-%d.inp", local_input_dir, header, workerid, myrank );
    ierr = ofmo_decompress_data( filename, comped_data, comped_data_size );
    if ( ierr != 0 ) {
	if ( is_root ) {
	    dbg("Failure in decompress data\n");
	    MPI_Abort( MPI_COMM_WORLD, 10002 );
	}
    }
    if ( ofmo_init( filename, MPI_COMM_WORLD ) != 0 ) {
	if ( is_root ) dbg("Failure in ofmo_init\n");
	MPI_Abort( MPI_COMM_WORLD, 10003 );
    }
    ofmo_data_put_vals("nintic", eribfsz );
    //ofmo_show_input_data();

#ifdef USE_CUDA
//    int ndev = 0;
    if ( fp_prof ) fprintf( fp_prof,"ndev = %d\n", ndev);
    ierr = cuda_Init(ndev, myrank, nprocs, MPI_COMM_WORLD);
    if (ierr<0) {
      dbg("Failure in cuda_Init\n");
      MPI_Abort( MPI_COMM_WORLD, 10100 );
    }
#endif

    // メモリ確保
    int maxnfao, nbody, maxnfatom, maxlqn;
    double *daopop, *datpop;	// ダイマーでのpopulation変位情報
    int *fsao2tuao, *fatom2tatom;
    int fnao, fnat;
    ofmo_data_get_vals("nbody maxnfao maxnfatom maxlqn",
	    &nbody, &maxnfao, &maxnfatom, &maxlqn );
    ofmo_frag_init();
    ofmo_integ_init( maxlqn );
    //ofmo_twoint_init();
    ofmo_scf_init( nbody * maxnfao );
    daopop      = (double*)malloc( sizeof(double) * nbody * maxnfao );
    datpop      = (double*)malloc( sizeof(double) * nbody * maxnfatom );
    fsao2tuao   = (int*)malloc( sizeof(int) * nbody * maxnfao );
    fatom2tatom = (int*)malloc( sizeof(int) * nbody * maxnfatom );
    // ---- マスタープロセスからの計算依頼待ち ----
    int imsg[OFMO_IMSG_SZ], dmsg_sz=OFMO_DMSG_SZ;
    int tag_ret=OFMO_TAG_RET;
    double dmsg[OFMO_DMSG_SZ];
    double *es;
    es = (double*)malloc( sizeof(double) * nfrag );
    //MPI_Status status;
    // imsgの内容（予定）
    //     0 : 計算の種類
    //         OFMO_INIT_DENS = 初期密度行列計算
    //         OFMO_SCF       = 近似なし電子状態計算
    //         OFMO_APPROX    = 近似あり電子状態計算
    //         OFMO_FINALIZE  = すべての計算終了のシグナル
    //     1: 計算手法
    //         OFMO_RHF       = RHF
    //         OFMO_RIMP2     = RIMP2
    //         OFMO_DFT       = DFT
    //     2: SCC繰り返し回数
    //     3: 収束条件
    //     4: 予約1
    //     5: 関係するモノマー数
    //        OFMO_APPROXの場合は、計算する近似ダイマー数
    //     6: モノマー１
    //     7: モノマー２
    //     8: ...
    int *joblist, njob;
    int itol, iscc, iscc_prev = -1;
    int nmonomer, *monomer_list;
    joblist = &imsg[OFMO_I_MON1];
    while (1) {
	memset( dmsg, '\0', sizeof(double)*dmsg_sz );
	// ジョブ依頼のシグナルをマスタープロセスから受け取る
	MPI_Bcast( imsg, OFMO_IMSG_SZ, MPI_INT, master_root,
		comm_master );
	if ( imsg[OFMO_I_CMD] == OFMO_ACCEPT ) {
	    int mserv_id;
	    mserv_id = imsg[1];
	    ofmo_connect_mserv( mserv_id, MPI_COMM_WORLD );
	    if ( is_root ) {
		MPI_Send( dmsg, dmsg_sz, MPI_DOUBLE, master_root, tag_ret,
			comm_master );
	    } 
	} else if ( imsg[OFMO_I_CMD] == OFMO_DISCONNECT ) {
	    int mserv_id;
	    mserv_id = imsg[OFMO_I_METHOD];
	    ofmo_disconnect_mserv( mserv_id );
	    if ( is_root ) {
		MPI_Send( dmsg, dmsg_sz, MPI_DOUBLE, master_root, tag_ret,
			comm_master );
	    } 
	} else if ( imsg[OFMO_I_CMD] == OFMO_RESET_MON_DATA ) {
	    ofmo_reset_monomer_data_workerp();
	    if ( is_root ) {
		MPI_Send( dmsg, dmsg_sz, MPI_DOUBLE, master_root,
			tag_ret, comm_master );
	    }
	} else if ( imsg[OFMO_I_CMD] == OFMO_DISTANCE ) {
	    ofmo_make_inter_frag_distance_list( MPI_COMM_WORLD );
	    if ( is_root ) {
		MPI_Send( dmsg, dmsg_sz, MPI_DOUBLE, master_root,
			tag_ret, comm_master );
	    }
	} else if ( imsg[OFMO_I_CMD] == OFMO_BARRIER ) {
	    if ( is_root ) ofmo_worker_barrierp();
	    MPI_Barrier( MPI_COMM_WORLD );
	    if ( is_root ) {
		MPI_Send( dmsg, dmsg_sz, MPI_DOUBLE, master_root,
			tag_ret, comm_master );
	    }
	} else if ( imsg[OFMO_I_CMD] == OFMO_INIT_DENS ) {
	    ofmo_monomer_init_density( imsg, MPI_COMM_WORLD );
	    if ( is_root ) {
		MPI_Send( dmsg, dmsg_sz, MPI_DOUBLE, master_root,
			tag_ret, comm_master );
	    }
	} else if ( imsg[OFMO_I_CMD] == OFMO_UPDATE_DATA ) {
	    ofmo_update_monomer_data();
	    if ( is_root ) {
		MPI_Send( dmsg, dmsg_sz, MPI_DOUBLE, master_root,
			tag_ret, comm_master );
	    }
	} else if ( imsg[OFMO_I_CMD] == OFMO_CHANGE_MSERV ) {
	    int mserv_id;
	    mserv_id = imsg[OFMO_I_METHOD];
	    ofmo_set_default_mserv( mserv_id );
	    if ( fp_prof ) {
		fprintf( fp_prof,
			"== default memory server has been changed to %d\n",
			mserv_id );
		fflush( fp_prof );
	    }
	    if ( is_root ) {
		MPI_Send( dmsg, dmsg_sz, MPI_DOUBLE, master_root,
			tag_ret, comm_master );
	    }
	} else if ( imsg[OFMO_I_CMD] == OFMO_SCF ) {
	    double energy, energy0, ddv, tolscf;
	    nmonomer     = imsg[OFMO_I_NMON];
	    if ( nmonomer > nbody || nmonomer == 0 ) {
		energy = 0.e0;
		energy0 = 0.e0;
		ddv = 0.e0;
	    } else {
		monomer_list = &imsg[OFMO_I_MON1];
		iscc = imsg[OFMO_I_SCC];
		itol = imsg[OFMO_I_CONV];
		if ( iscc != iscc_prev ) {
		    if ( is_root )
			ofmo_get_monomer_energy( -1, es );
		    MPI_Bcast( es, nfrag, MPI_DOUBLE, root, MPI_COMM_WORLD );
		    ofmo_set_new_scc_flag();
		    iscc_prev = iscc;
		    if ( nmonomer == 1 && fp_prof ) {
			fprintf( fp_prof,
				"==== SCC itera= %d ====\n", iscc );
			fflush( fp_prof );
		    }
		}
		tolscf = pow( 0.1e0, (double)itol );
		if ( nmonomer == 1 ) energy = es[ monomer_list[0] ];
		ofmo_calc_fragment_electronic_state(
			MPI_COMM_WORLD, nmonomer, monomer_list,
			imsg[OFMO_I_METHOD], tolscf,
			&energy, &energy0, &ddv,
			&fnao, daopop, fsao2tuao,
			&fnat, datpop, fatom2tatom );
	    }
	    if ( is_root ) {
		dmsg[OFMO_D_ENERGY]  = energy;
		dmsg[OFMO_D_ENERGY0] = energy0;
		dmsg[OFMO_D_DDV]     = ddv;
		MPI_Send( dmsg, dmsg_sz, MPI_DOUBLE, master_root,
			tag_ret, comm_master );
	    }
	} else if ( imsg[OFMO_I_CMD] == OFMO_GET_POP_DATA ) {
	    int mode;
	    mode = imsg[OFMO_I_METHOD];
	    if ( is_root ) {
		if ( mode == OFMO_BOTH_POPS || mode == OFMO_AOPOP_ONLY ) {
		    MPI_Send( daopop, fnao, MPI_DOUBLE, master_root,
			    tag_ret, comm_master );
		    MPI_Send( fsao2tuao, fnao, MPI_INT, master_root,
			    tag_ret, comm_master );
		}
		if ( mode == OFMO_BOTH_POPS || mode == OFMO_ATPOP_ONLY ) {
		    MPI_Send( datpop, fnat, MPI_DOUBLE, master_root,
			    tag_ret, comm_master );
		    MPI_Send( fatom2tatom, fnat, MPI_INT, master_root,
			    tag_ret, comm_master );
		}
	    }
	} else if ( imsg[OFMO_I_CMD] == OFMO_APPROX ) {
	    double eesdimer0;
	    if ( nbody < 2 ) {
		eesdimer0 = 0.e0;
	    } else {
		njob = imsg[OFMO_I_NMON];
		if ( njob == 0 ) {
		    eesdimer0 = 0.e0;
		} else {
		    ofmo_calc_es_dimer( MPI_COMM_WORLD, njob, joblist,
			    &eesdimer0 );
		}
	    }
	    if ( is_root ) {
		dmsg[OFMO_D_ENERGY0] = eesdimer0;
		MPI_Send( dmsg, dmsg_sz, MPI_DOUBLE, master_root,
			tag_ret, comm_master );
	    }
	} else if ( imsg[OFMO_I_CMD] == OFMO_FINALIZE ) {
	    MPI_Barrier( MPI_COMM_WORLD );
	    if ( is_root ) {
		MPI_Send( dmsg, dmsg_sz, MPI_DOUBLE, master_root,
			tag_ret, comm_master );
	    }
	    break;
	}
    }
    ofmo_show_cache_prof();
    if ( fp_prof ) {
	fprintf( fp_prof, "== finish all calculation ==\n");
	fflush( fp_prof );
    }
    Free( es );
    Free( daopop );
    Free( datpop );
    Free( fsao2tuao );
    Free( fatom2tatom );

#ifdef USE_CUDA
    ierr = cuda_Finalize();
    MPI_Barrier( MPI_COMM_WORLD );
#endif

    MPI_Finalize();
    return 0;
}
