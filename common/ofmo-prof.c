/** プロファイルファイル出力のための関数群

  各ワーカー、メモリサーバーなど、同一のコミュニケータを持つプロセスで
  取得したプロファイル情報を出力するための関数群。
  以下の仮定をしている。

  - 各ワーカー、メモリサーバーで、１つのプロファイル情報を取得する
  - rank=0のプロセスがファイルへの出力を行う
  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "ofmo-parallel.h"

#include "ofmo-def.h"
#include "ofmo-misc.h"

FILE* fp_prof = NULL;

#ifdef OFMO_PROF

static MPI_Comm COMM = MPI_COMM_NULL;
static int initialized = false;

/* 小規模電子状態計算時におけるスレッド毎のプロファイル取得に関する関数 */

#define MAX_THRD_CNTR 12	// カウンタ種類の最大値

static int MAXTHREADS = 1;
static double **thrd_timer_val  = NULL;
static double **thrd_timer_val0 = NULL;
static char **thrd_timer_name  = NULL;

static int *thrd_timer_used    = NULL;	// 使用されているかどうか
static int *thrd_timer_attr   = NULL;	// 時間か、％値か
static double *thrd_timer_tmp    = NULL;

static int timer_val_size     = 0;
static int timer_val_elem     = 0;

static void dealloc_thread_timer() {
    ofmo_free_dmatrix( thrd_timer_val );
    ofmo_free_dmatrix( thrd_timer_val0 );
    ofmo_free_cmatrix( thrd_timer_name );
    Free( thrd_timer_used );
    Free( thrd_timer_attr );
    Free( thrd_timer_tmp );
    timer_val_size = 0;
    timer_val_elem = 0;
}

static int alloc_thread_timer() {
    int maxthreads, i;
    static int called = false;
    if ( !called ) {
	maxthreads = omp_get_max_threads();

	thrd_timer_val  = ofmo_alloc_dmatrix( MAX_THRD_CNTR, maxthreads );
	thrd_timer_val0 = ofmo_alloc_dmatrix( MAX_THRD_CNTR, maxthreads );
	thrd_timer_name = ofmo_alloc_cmatrix( MAX_THRD_CNTR, MAXSTRLEN );

	thrd_timer_used = (int*)malloc( sizeof(int) * MAX_THRD_CNTR );
	thrd_timer_attr = (int*)malloc( sizeof(int) * MAX_THRD_CNTR );
	thrd_timer_tmp =
	    (double*)malloc( sizeof(double) * MAX_THRD_CNTR * maxthreads );

	timer_val_size = MAX_THRD_CNTR * maxthreads * sizeof(double);
	timer_val_elem = MAX_THRD_CNTR * maxthreads;
	atexit( dealloc_thread_timer );
	MAXTHREADS = maxthreads;
	called = true;
    }
    // initialize
    for ( i=0; i<MAX_THRD_CNTR; i++ ) {
	thrd_timer_used[i] = false;
	thrd_timer_attr[i] = 0;	// 0=時間, 1=%値
	thrd_timer_name[i][0] = '\0';
    }
    return 0;
}

static void reset_thread_timer() {
    int i;
    for ( i=0; i<MAX_THRD_CNTR; i++ ) {
	thrd_timer_used[i] = false;
	thrd_timer_name[i][0] = '\0';
    }
}


/* プロセス毎のタイマー */
#define MAX_PROC_CNTR 25
static double proc_timer_val[MAX_PROC_CNTR];
static double proc_timer_val0[MAX_PROC_CNTR];
static char   proc_timer_name[MAX_PROC_CNTR][MAXSTRLEN];
static int    proc_timer_used[MAX_PROC_CNTR];
//static double proc_timer_tmp[MAX_PROC_CNTR];
static int    proc_timer_level[MAX_PROC_CNTR];
static int proc_timer_size      = 0;
static int proc_timer_elem      = 0;
static double *proc_all_val      = NULL;

static void dealloc_proc_timer() {
    Free( proc_all_val );
    proc_all_val    = NULL;
    proc_timer_size = 0;
    proc_timer_elem = 0;
}

static int alloc_proc_timer() {
    static int called = false;
    int nprocs;
    MPI_Comm comm = COMM;
    MPI_Comm_size( comm, &nprocs );
    dealloc_proc_timer();
    proc_all_val = (double*)malloc( sizeof(double)*nprocs*MAX_PROC_CNTR);
    for ( int i=0; i<MAX_PROC_CNTR; i++ ) {
	proc_timer_used[i]  = false;
	proc_timer_level[i] = -1;
    }
    proc_timer_size = sizeof(double) * MAX_PROC_CNTR;
    proc_timer_elem = MAX_PROC_CNTR;
    if ( called ) {
	atexit( dealloc_proc_timer );
	called = true;
    }
    return 0;
}

static void reset_proc_timer() {
    int i;
    for ( i=0; i<MAX_PROC_CNTR; i++ ) {
	proc_timer_used[i]  = false;
    }
}

static void ofmo_prof_finalize() {
    if ( fp_prof != NULL ) {
        fclose( fp_prof );
        fp_prof = NULL;
        COMM = MPI_COMM_NULL;
    }
}

static void ofmo_open_profile( const char *prof_name, MPI_Comm comm ) {
    int myrank, root=0;
    static char sprof[MAXSTRLEN]="";
    assert(fp_prof==NULL);
    MPI_Comm_rank( comm, &myrank );
    if ( myrank == root ) {
      if (prof_name == NULL) {
        fp_prof = stdout;
      } else {
        char *mode="w";
        if (strncmp(sprof, prof_name, MAXSTRLEN)==0) {mode="a";}
	//if ( (fp_prof=fopen( prof_name, "w" )) == NULL ) {
	if ( (fp_prof=fopen( prof_name, mode )) == NULL ) {
	    dbg("failure in file open (%s)\n", prof_name );
	    MPI_Abort( comm, 11 );
	}
        strncpy(sprof, prof_name, MAXSTRLEN);
      }
    } else {
	fp_prof = NULL;
    }
}
/** プロファイルファイルの初期化関数

  引数で渡されているコミュニケータが、既に登録されているものと同一の
  場合には、何もしない（現在使用しているプロファイル用ファイルを使用する）

  異なるコミュニケータを引数として、同一ファイル名を指定すると、
  上書きされてしまう
 * */
int ofmo_prof_init( const char *prof_name, MPI_Comm comm ) {
    int result = MPI_UNEQUAL;
    static int called = false;
    if ( called ) {
        //if (COMM != NULL && COMM != MPI_COMM_NULL) {
        if (COMM != MPI_COMM_NULL) {
            MPI_Comm_compare( comm, COMM, &result );
        }
	if ( result != MPI_IDENT ) {
	    ofmo_open_profile( prof_name, comm );
	    COMM = comm;
	    alloc_proc_timer();
	    alloc_thread_timer();
	}
    } else {
	ofmo_open_profile( prof_name, comm );
	COMM = comm;
	alloc_proc_timer();
	alloc_thread_timer();
	atexit( ofmo_prof_finalize );
	called = true;
    }
    initialized = true;
    return 0;
}

int ofmo_prof_reinit( const char *prof_name, MPI_Comm comm ) {
    ofmo_prof_finalize();
    COMM = MPI_COMM_NULL;
    return ofmo_prof_init(prof_name, comm);
}

/* ======= スレッド毎のタイマーに関連する関数群の定義 ====== */
static int get_thread_timer_id() {
    int id;
    if ( ! initialized ) return -2;
    for ( id=0; id<MAX_THRD_CNTR; id++ ) {
	if ( !thrd_timer_used[id] ) {
	    thrd_timer_used[id] = true;
	    return id;
	}
    }
    return -1;
}

/** スレッド毎のタイマーカウンタを作成する
  * 
  */
int ofmo_create_thread_timer( const char *str, const int attr ) {
    int id;
    id = get_thread_timer_id();
    if ( id >= 0 ) {
	strcpy( thrd_timer_name[id], str );
	thrd_timer_attr[id] = attr;
	memset( thrd_timer_val[id], '\0', sizeof(double)*MAXTHREADS );
    }
    return id;
}

/** すべてのカウンタの値をリセットする */
void ofmo_reset_all_thread_timer() {
    if ( ! initialized ) return;
    memset( thrd_timer_val[0], '\0', timer_val_size );
}

/** タイマーの測定を開始する（経過時間0の時刻を計測、記憶する) */
void ofmo_start_thread_timer( const int id, const int mythread ) {
    thrd_timer_val0[id][mythread] = omp_get_wtime();
}

/** タイマーの経過時間に、直前のofmo_start_thread_timer関数呼び出し
 からの経過時間を加算する */
void ofmo_acc_thread_timer( const int id, const int mythread ) {
    double t0, t1;
    t0 = thrd_timer_val0[id][mythread];
    t1 = omp_get_wtime();
    thrd_timer_val[id][mythread] += (t1 - t0);
}

/** タイマーの値をセットする */
void ofmo_set_thread_timer( const int id, const int mythread,
	const double val ) {
    thrd_timer_val[id][mythread] = val;
}

/** すべてのカウンタの値を表示する */
void ofmo_show_thread_timer_all() {
    int myrank, ithrd, irank, nthreads, i, nprocs, tag=123, root=0;
    MPI_Comm comm = COMM;
    MPI_Comm_rank( comm, &myrank );
    MPI_Comm_size( comm, &nprocs );
    nthreads = omp_get_max_threads();
    if ( myrank == root ) {
	MPI_Status status;
	for ( irank=0; irank<nprocs; irank++ ) {
	    if ( irank == root ) {
		memcpy(thrd_timer_tmp, thrd_timer_val[0], timer_val_size);
	    } else {
		MPI_Recv( thrd_timer_tmp, timer_val_elem, MPI_DOUBLE,
			irank, tag, comm, &status );
	    }
	    for ( ithrd=0; ithrd<nthreads; ithrd++ ) {
		fprintf(fp_prof,"## rank= %3d thrd= %2d", irank, ithrd );
		for ( i=0; i<MAX_THRD_CNTR; i++ ) {
		    if ( !thrd_timer_used[i] ) continue;
		    if ( thrd_timer_attr[i] == 0 ) {
			fprintf(fp_prof, " %s= %9.6f",
				thrd_timer_name[i],
				thrd_timer_tmp[i*nthreads + ithrd] );
		    } else if ( thrd_timer_attr[i] == 1 ) {
			fprintf(fp_prof, " %s= %5.1f",
				thrd_timer_name[i],
				thrd_timer_tmp[i*nthreads + ithrd] );
		    }
		}
		fprintf(fp_prof, "\n");
	    }
	}
    } else {
	MPI_Send( thrd_timer_val[0], timer_val_elem, MPI_DOUBLE,
		root, tag, comm );
    }
    fflush(fp_prof);
    reset_thread_timer();
}

/** 指定されたIDを持つカウンタの値を出力する */
void ofmo_show_thread_timers( int n, int list[] ) {
    int myrank, ithrd, irank, nthreads, i, nprocs, tag=123, root=0;
    MPI_Comm comm = COMM;
    MPI_Comm_rank( comm, &myrank );
    MPI_Comm_size( comm, &nprocs );
    nthreads = omp_get_max_threads();
    if ( myrank == root ) {
	MPI_Status status;
	int id;
	for ( irank=0; irank<nprocs; irank++ ) {
	    if ( irank == root ) {
		memcpy(thrd_timer_tmp, thrd_timer_val[0], timer_val_size);
	    } else {
		MPI_Recv( thrd_timer_tmp, timer_val_elem, MPI_DOUBLE,
			irank, tag, comm, &status );
	    }
	    for ( ithrd=0; ithrd<nthreads; ithrd++ ) {
		fprintf(fp_prof,"## rank= %3d thrd= %2d", irank, ithrd );
		for ( i=0; i<n; i++ ) {
		    id = list[i];
		    if ( !thrd_timer_used[id] ) continue;
		    if ( thrd_timer_attr[id] == 0 ) {
			fprintf(fp_prof, " %s= %9.6f",
				thrd_timer_name[id],
				thrd_timer_tmp[id*nthreads + ithrd] );
		    } else if ( thrd_timer_attr[id] == 1 ) {
			fprintf(fp_prof, " %s= %5.1f",
				thrd_timer_name[id],
				thrd_timer_tmp[id*nthreads + ithrd] );
		    }
		}
		fprintf(fp_prof, "\n");
	    }
	}
    } else {
	MPI_Send( thrd_timer_val[0], timer_val_elem, MPI_DOUBLE,
		root, tag, comm );
    }
    fflush(fp_prof);
    reset_thread_timer();
}

static int get_proc_timer_id() {
    if ( ! initialized ) return -2;
    for ( int id=0; id<MAX_PROC_CNTR; id++ ) {
	if ( !proc_timer_used[id] ) {
	    proc_timer_used[id] = true;
	    return id;
	}
    }
    return -1;
}

int ofmo_create_proc_timer( const char *str, const int level ) {
    int id;
    if ( !initialized ) return -2;
    id = get_proc_timer_id();
    if ( id >= 0 ) {
	strcpy( proc_timer_name[id], str );
	proc_timer_level[id] = level;
	proc_timer_val[id]   = 0.e0;
    }
    return id;
}

void ofmo_reset_all_proc_timer() {
    memset( proc_timer_val, '\0', proc_timer_size );
}

void ofmo_start_proc_timer( const int id ) {
    proc_timer_val0[id] = MPI_Wtime();
}

void ofmo_acc_proc_timer( const int id ) {
    double t;
    t = MPI_Wtime();
    proc_timer_val[id] += ( t - proc_timer_val0[id] );
}

void ofmo_show_proc_timer_all() {
    int myrank, nprocs, irank, id, iat, root=0;
    int max_level, nl;
    MPI_Comm comm = COMM;
    double *val;
    MPI_Comm_rank( comm, &myrank );
    MPI_Comm_size( comm, &nprocs );

    max_level = proc_timer_level[0];
    for ( id=1; id<proc_timer_elem; id++ ) {
	if ( proc_timer_used[id] && proc_timer_level[id] > max_level )
	    max_level = proc_timer_level[id];
    }
    MPI_Gather( proc_timer_val, proc_timer_elem, MPI_DOUBLE,
	    proc_all_val, proc_timer_elem, MPI_DOUBLE, root, comm );
    if ( myrank == root ) {
	for ( iat=0; iat<=max_level; iat++ ) {
	    nl = 0;
	    for ( id=0; id<proc_timer_elem; id++ ) {
		if ( proc_timer_used[id] && proc_timer_level[id] == iat )
		    nl++;
	    }
	    if ( nl == 0 ) continue;
	    for ( irank=0, val=proc_all_val; irank<nprocs;
		    irank++, val+=proc_timer_elem ) {
		fprintf(fp_prof, "##%d rank= %3d", iat, irank );
		for ( id=0; id<proc_timer_elem; id++ ) {
		    if ( proc_timer_level[id] == iat ) {
			fprintf(fp_prof, " %s= %9.6f",
				proc_timer_name[id], val[id] );
		    }
		}
		fprintf(fp_prof, "\n");
	    }
	}
    }
    fflush(fp_prof);
    reset_proc_timer();
}

void ofmo_show_proc_timers( int n, int list[] ) {
    int myrank, nprocs, irank, id, iat, root=0, i;
    int max_level, nl;
    MPI_Comm comm = COMM;
    double *val;
    MPI_Comm_rank( comm, &myrank );
    MPI_Comm_size( comm, &nprocs );

    max_level = proc_timer_level[0];
    for ( id=1; id<proc_timer_elem; id++ ) {
	if ( proc_timer_level[id] > max_level )
	    max_level = proc_timer_level[id];
    }
    MPI_Gather( proc_timer_val, proc_timer_elem, MPI_DOUBLE,
	    proc_all_val, proc_timer_elem, MPI_DOUBLE, root, comm );
    if ( myrank == root ) {
	for ( iat=0; iat<=max_level; iat++ ) {
	    nl=0;
	    for ( i=0; i<n; i++ ) {
		id = list[i];
		if ( proc_timer_level[id] == iat ) nl++;
	    }
	    if ( nl==0 ) continue;
	    for ( irank=0, val=proc_all_val; irank<nprocs;
		    irank++, val+=proc_timer_elem ) {
		fprintf(fp_prof, "##%d rank= %3d", iat, irank );
		for ( i=0; i<n; i++ ) {
		    id = list[i];
		    if ( proc_timer_level[id] == iat ) {
			fprintf(fp_prof, " %s= %9.6f",
				proc_timer_name[id], val[id] );
		    }
		}
		fprintf(fp_prof, "\n");
	    }
	}
    }
    fflush(fp_prof);
    reset_proc_timer();
}

/* ワーカー用プロファイルのファイル名を作成する */
int ofmo_create_worker_prof_name( const char *input, const char *header,
	const int group_id, char *prof_name, size_t maxstrlen ) {
    char prefix[MAXSTRLEN];
    int pos;
    if ( header[0] == '\0' ) {
	pos = strcspn( input, "." );
	strncpy( prefix, input, pos );
	prefix[pos] = '\0';
    } else {
	strncpy( prefix, header, MAXSTRLEN );
	prefix[MAXSTRLEN-1] = '\0';
    }
    snprintf( prof_name, maxstrlen, "%s-worker-prof-%04d.log",
	    prefix, group_id );
    prof_name[maxstrlen-1] = '\0';
    return 0;
}

/* memory server用プロファイルのファイル名を作成する */
int ofmo_create_mserv_prof_name( const char *header, const int mserv_id,
	char *prof_name, size_t maxstrlen ) {
    char prefix[MAXSTRLEN];
    if ( header[0] == '\0' ) strncpy( prefix, "profile", MAXSTRLEN );
    else                     strncpy( prefix, header, MAXSTRLEN );
    prefix[MAXSTRLEN-1] = '\0';
    snprintf( prof_name, maxstrlen, "%s-memserv-prof-%02d.log",
	    prefix, mserv_id );
    prof_name[maxstrlen-1] = '\0';
    return 0;
}

#else

int ofmo_prof_init( const char *prof_name, const int root, MPI_Comm comm ){
    return 0;
}

int ofmo_create_worker_prof_name( const char *input, const char *header,
	const int group_id, char *prof_name, size_t maxstrlen ) {
    return 0;
}

int ofmo_create_mserv_prof_name( const char *header, const int mserv_id,
	char *prof_name, size_t maxstrlen ) {
    return 0;
}

// スレッド毎のタイマー
int ofmo_create_thread_timer( const char *str, const int attr ) {
    return 0; }
void ofmo_reset_all_thread_timer() {};
void ofmo_start_thread_timer( const int id, const int mythread ) {};
void ofmo_acc_thread_timer( const int id, const int mythread ) {};
void ofmo_set_thread_timer( const int id, const int mythread,
	const double val ) {};
void ofmo_show_thread_timer_all() {};
void ofmo_show_thread_timers( int n, int list[] ) {};

// プロセス毎のタイマー
int ofmo_create_proc_timer( const char *str, const int attr ) { return 0;}
void ofmo_reset_all_proc_timer() {};
void ofmo_start_proc_timer( const int id ) {};
void ofmo_acc_proc_timer( const int id ) {};
void ofmo_show_proc_timer_all() {};
void ofmo_show_proc_timers( int n, int list[] ) {};

#endif /* OFMO_PROF */
