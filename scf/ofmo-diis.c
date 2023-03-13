/**
 * @file ofmo-diis.c
 * A file that defines functions for accelerating SCF convergence.
 *
 * @attention
 * - この関数群からは、\c ofmo-mat.h で定義されている関数を呼び出している
 *
 * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ofmo-def.h"
#include "ofmo-mat.h"
#include "ofmo-prof.h"

#ifndef ZERO
#define ZERO	0.e0
#endif
#ifndef TWO
#define TWO	2.e0
#endif
#ifndef HALF
#define HALF	0.5e0
#endif
#ifndef ONE
#define ONE	1.e0
#endif

#ifndef OFMO_MAXDIIS
#define OFMO_MAXDIIS 5
#endif

// DIIS関連
static double *E_PROF_MASTER;
static double *F_PROF_MASTER;
//static double *D_PROF_MASTER;
static double **Eprof;
static double **Fprof;
//static double **Dprof;
static double *BDIIS;	// B行列
static double *B_TMP;	// B行列の作業用配列（連立方程式を解くとき）
static double *bDIIS;	// 連立方程式の解の行列
static double *S12;	// S^{-1/2}を代入
static double *S;
// temporary memory
static double *T1;
static double *T2;
static double *T3;
static int *IPIV;

static int NAONOW = 0;
static int MAXDIIS = 0;
static int MAXNAO = 0;
static int diis_now = 0;	// 現在保存されているDIISの履歴数を表す。
				// 履歴がある場合は、1以上の値をとる
//static size_t diis_total_mem = 0;

//#define Free(p) if ((p)!=NULL) free( (p) )


// debug
/*static int debug_flag = false;
void set_debug_flag( const int flag ) {
    debug_flag = flag;
}*/

/** エラーベクトルの計算
 *
 * DIISで用いるエラーベクトルを計算する関数。次式で計算している。
 * \f[
 * \bold{\epsilon} =
 * \bold{U}^{-1}\left(\bold{FDS} -
 * \bold{SDF}\right) \bold{U},
 * \ \ \bold{U} = \bold{S}^{1/2}
 * \f]
 * 
 * @param[in] n 行列サイズ（通常はAO数）
 * @param[in] Fp[] Fock行列（圧縮U形式）
 * @param[in] Dp[] 密度行列（圧縮U形式）
 * @param[out] Ep[] 計算されたエラーベクトル。エラーベクトルは反対称
 *     行列（\f$ a_{ji} = -a_{ij},\ a_{ii}=0 \f$ ）であるが、その
 *     上三角行列部分が計算・代入される。
 *
 * @return なし
 * @ingroup ofmo-diis
 * */
static double diis_error_vector( const int n, const double Fp[],
	const double Dp[], double Ep[] ) {
    double *D=T1, *F=T2, *DS=T3, *FDS=T1, *E=T2, t, *ES=T3, *SES=T1;
    double alpha = 1.e0, beta = 0.e0;
    int i, j, ij, n2;
    ofmo_unpack_matrix( n, Dp, D );
    ofmo_unpack_matrix( n, Fp, F );
    ofmo_dgemm( n, "T", "N", alpha, D, S, beta, DS );
    ofmo_dgemm( n, "T", "N", alpha, F, DS, beta, FDS );
    for ( i=0; i<n; i++ ) {
	for ( j=0; j<=i; j++ ) {
	    t = E[i*n+j] = FDS[i*n+j] - FDS[j*n+i];
	    E[j*n+i] = -t;
	}
    }
    ofmo_dgemm( n, "T", "N", alpha, E, S12, beta, ES );
    ofmo_dgemm( n, "T", "N", alpha, ES, S12, beta, SES );
    ofmo_pack_matrix( n, SES, Ep );
    n2 = n*(n+1)/2;
    t = 0.e0;
    for ( ij=0; ij<n2; ij++ ) if ( fabs(Ep[ij]) > t ) t = fabs(Ep[ij]);
    return t;
}

/** DIIS関数群で利用しているメモリ領域を開放する関数
 *
 * DIIS関数群で利用しているメモリ領域を開放する関数。実行終了時に
 * 自動的に呼ばれるので、ユーザが露に呼ぶ必要はない。
 * */
static void diis_dealloc() {
    Free(bDIIS);	bDIIS = NULL;
    Free(BDIIS);	BDIIS = NULL;
    Free(B_TMP);	B_TMP = NULL;
    Free(Fprof);	Fprof = NULL;
    Free(Eprof);	Eprof = NULL;
    Free(F_PROF_MASTER);	F_PROF_MASTER = NULL;
    Free(E_PROF_MASTER);	E_PROF_MASTER = NULL;
    Free( S ); S = NULL;

    Free(S12);	S12 = NULL;
    Free(T1);	T1  = NULL;
    Free(T2);	T2  = NULL;
    Free(T3);	T3  = NULL;
    Free(IPIV);	IPIV = NULL;

    //diis_total_mem = 0;
    diis_now = 0;
    MAXDIIS  = 0;
    MAXNAO   = 0;
}

/** DIIS関数の内部で使用する各種配列を確保する関数
 *
 * DIIS関数の内部で使用する各種配列を確保する関数.
 *
 * @param[in] maxdiis DIISで保存する履歴の最大数
 * @param[in] maxnao （最大）AO数
 *
 * @retval  0 正常終了（必要なすべての領域を確保した）
 * @retval -1 異常終了（メモリの確保に失敗した）
 *
 * */
static int diis_alloc( const int maxdiis, const int maxnao) {
    int nao2, i, nnao;
    size_t total_mem;
    MAXDIIS = maxdiis;
    MAXNAO  = maxnao;
    nao2 = maxnao * (maxnao + 1) / 2;
    nnao = maxnao * maxnao;
    total_mem = 0;
    E_PROF_MASTER = (double*)malloc( sizeof(double) * nao2 * maxdiis);
    F_PROF_MASTER = (double*)malloc( sizeof(double) * nao2 * maxdiis);
    if ( E_PROF_MASTER == NULL || F_PROF_MASTER == NULL ) {
	diis_dealloc();
	return -1;
    }
    total_mem = 2 * sizeof(double) * nao2 * maxdiis;
    // 先頭アドレスを保存するポインタ配列確保
    Eprof = (double**)malloc( sizeof(double*) * maxdiis);
    Fprof = (double**)malloc( sizeof(double*) * maxdiis);
    total_mem += 2 * sizeof(double*) * maxdiis;
    if ( Eprof == NULL || Fprof == NULL ) {
	diis_dealloc();
	return -1;
    }
    Eprof[0] = E_PROF_MASTER;
    Fprof[0] = F_PROF_MASTER;
    for (i=1; i<maxdiis; i++) {
	Eprof[i] = Eprof[i-1] + nao2;
	Fprof[i] = Fprof[i-1] + nao2;
    }
    // B行列
    BDIIS = (double*)malloc( sizeof(double) * (maxdiis+1)*(maxdiis+1) );
    B_TMP = (double*)malloc( sizeof(double) * (maxdiis+1)*(maxdiis+1) );
    bDIIS = (double*)malloc( sizeof(double) * (maxdiis+1) );
    if ( BDIIS == NULL || B_TMP == NULL || bDIIS == NULL ) {
	diis_dealloc();
	return -1;
    }
    total_mem += 2 * sizeof(double) * (maxdiis+1)*(maxdiis+1);
    total_mem += sizeof(double) * (maxdiis+1);
    // S^(1/2) & temporary memory
    S   = (double*)malloc( sizeof(double) * nnao );
    S12 = (double*)malloc( sizeof(double) * nnao );
    T1  = (double*)malloc( sizeof(double) * nnao );
    T2  = (double*)malloc( sizeof(double) * nnao );
    T3  = (double*)malloc( sizeof(double) * nnao );
    IPIV = (int*)malloc( sizeof(int) * maxnao );
    if ( S12 == NULL || T1 == NULL || T2 == NULL || T3 == NULL
	    || IPIV == NULL ) {
	diis_dealloc();
	return -1;
    }
    total_mem += 4 * sizeof(double) * nnao + sizeof(int) * maxnao;

    diis_now = 0;
    //diis_total_mem = total_mem;
    if ( fp_prof ) {
	fprintf( fp_prof, "Allocated memory in ofmo_diis = %8.3f MB\n",
		(double)((double)total_mem/(double)(1024*1024) ) );
	fflush( fp_prof );
    }
    return 0;
}

/** DIIS関数の内部で用いる各種配列を確保する関数
 *
 * DIIS関数の内部で用いる各種配列を確保する関数。
 * 実行一回につき１度しかSCF計算を行わないのであれば、露に呼び出す必要は
 * ない（他の外部関数から呼ばれるため）。しかし、FMO計算のように、
 * 実行一回につき複数回のSCF計算を行う場合には、遭遇しうる最大AO数を
 * 引数としてこの関数を呼び出したほうが、無駄な\c malloc を呼び出す
 * 必要がなく、性能面で好ましい。
 *
 * @param[in] maxdiis DIISで保存する履歴の最大数
 * @param[in] maxnao （最大）AO数
 * 
 * @retval  0 正常終了（必要なすべての領域を確保した）
 * @retval -1 異常終了（メモリの確保に失敗した）
 *
 * @ingroup ofmo-diis
 * */
int ofmo_diis_alloc( const int maxdiis, const int maxnao ) {
    int _maxdiis_, _maxnao_;
    static int called = false;
    _maxdiis_ = maxdiis;
    _maxnao_  = maxnao;
    if ( _maxdiis_ > MAXDIIS || _maxnao_ > MAXNAO ) {
	if ( _maxdiis_ < MAXDIIS ) _maxdiis_ = MAXDIIS;
	if ( _maxnao_ < MAXNAO   ) _maxnao_  = MAXNAO;
	diis_dealloc();
	if ( diis_alloc( _maxdiis_, _maxnao_ ) != 0 ) return -1;
	if ( called == false ) {
	    atexit( diis_dealloc );
	    called = true;
	}
    }
    return 0;
}


/** DIIS用のもっとも古い履歴を廃棄する
 *
 * DIIS用のもっとも古い履歴を廃棄する
 *
 * */
static void diis_shift_profile() {
    double *et, *ft, *B = BDIIS;
    int n = MAXDIIS, i, j;
    et = Eprof[0];
    ft = Fprof[0];
    for ( i=1; i<MAXDIIS; i++) {
	Eprof[i-1] = Eprof[i];
	Fprof[i-1] = Fprof[i];
    }
    Eprof[MAXDIIS-1] = et;
    Fprof[MAXDIIS-1] = ft;
    for ( i=1; i<diis_now; i++) {
	for ( j=1; j<diis_now; j++) {
	    B[(i-1)*n + (j-1)] = B[i*n + j];
	}
    }
    diis_now--;
}


/** DIIS用の履歴を登録する
 *
 * DIIS用の履歴を登録する。
 * この関数内部では以下のことを行っている。
 *
 * - 必要な場合には、B行列のシフト
 * - エラーベクトルの作成
 * - 作成したエラーベクトルが関係するB行列要素の計算
 * - 履歴の数 \c diis_now の更新（この値が０であれば、履歴が全く保存
 *   されていないことを意味する。
 *
 * @param[in] nao AO数
 * @param[in] Dp[] 密度行列（圧縮U形式）
 * @param[in] Fp[] Fock行列（圧縮U形式）
 *
 * @return なし
 *
 * @ingroup ofmo-diis
 */
double ofmo_diis_profiling(
	const int nao, const double Dp[], const double Fp[]) {
    int nao2, i;
    double *B = BDIIS, e_max;
    double diis_thr = 1.e-2;
    nao2 = nao * (nao+1) / 2;
    if (diis_now >= MAXDIIS) diis_shift_profile();
    ofmo_dcopy( nao2, Fp, Fprof[diis_now] );
    e_max = diis_error_vector( nao, Fp, Dp, Eprof[diis_now] );
    /*// debug
    if ( debug_flag ) {
	printf("max(e) = %12.6f\n", 2.e0*e_max );
	fflush(stdout);
    }*/
    if ( e_max > diis_thr ) {
	diis_now = 0;
	return -1;
    }
    // B行列の要素を計算する
    for ( i=0; i<=diis_now; i++) {
	B[i*MAXDIIS + diis_now] = B[diis_now*MAXDIIS + i]
	    = 8.e0 * ofmo_dot_product( nao2, Eprof[diis_now], Eprof[i] );
    }
    diis_now++;
    return e_max;
}

/** DIISを用いたFock行列の更新を行う
 * 
 * DIISを用いたFock行列の更新を行う.
 *
 * @param[in] nao AO数
 * @param[in] Dp[] 密度行列（圧縮U形式）
 * @param[in,out] Fp[] Fock行列（圧縮U形式）。入力時にはDIISでの更新前の
 *     素のFock行列。出力時にはDIISでの更新を行ったFock行列。
 *
 * @retval  0 正常終了（正常にDIIS更新が出来た）
 * @retval -1 異常終了（DIIS更新が出来なかった（連立方程式が
 *     解けなかった）
 * */
static int diis_update(
	const int nao, const double Dp[], const double Fp[],
	double Fdiis[] ) {
    double *B_SRC = BDIIS, *B = B_TMP, *c = bDIIS;
    int n, flag = false, nao2, i, j;
    nao2 = nao * (nao+1) / 2;

    // 連立方程式を解くための行列を作成する
    n = diis_now + 1;
    for ( i=0; i<diis_now; i++ ) {
	B[diis_now+i*n] = B[diis_now*n+i] = -1.e0;
	c[i] = 0.e0;
    }
    B[diis_now*n+diis_now] = 0.e0;
    c[diis_now] = -1.e0;
    for ( i=0; i<diis_now; i++ ) {
	for ( j=0; j<diis_now; j++ ) {
	    B[i*n + j] = B_SRC[i*MAXDIIS + j];
	}
    }
    /*// debug
    if ( debug_flag ) {
	printf("\n=================== B matrix =============\n");
	for ( i=0; i<=diis_now; i++ ) {
	    for ( j=0; j<=diis_now; j++ ) printf(" %9.6f", B[i*n+j] );
	    printf("\n");
	}
    }*/
    while (diis_now > 0) {
	if ( ofmo_solv_leq(n, B, c) == 0 ) {
	    flag = true;
	    /*// debug
	    if ( debug_flag ) {
		printf("b=");
		for ( int ii=0; ii<n; ii++ )
		    printf(" %10.5f", c[ii] );
		printf("\n");
		fflush(stdout);
	    }*/
	    break;
	} else {
	    diis_shift_profile();
	    n = diis_now + 1;
	    for ( i=0; i<diis_now; i++ ) {
		B[diis_now+i*n] = B[diis_now*n+i] = -1.e0;
		c[i] = 0.e0;
	    }
	    B[diis_now*n+diis_now] = 0.e0;
	    c[diis_now] = -1.e0;
	    for ( i=0; i<diis_now; i++ ) {
		for ( j=0; j<diis_now; j++ ) {
		    B[i*n + j] = B_SRC[i*MAXDIIS + j];
		}
	    }
	}
    }
    if (flag == false) return -1;	// 連立方程式がどうしても
					// 解けなかった場合
    // Fock行列を再構成する
    memset( Fdiis, '\0', sizeof(double)*nao2 );
    for ( i=0; i<(diis_now-1); i++ )
	ofmo_daxpy( nao2, c[i], Fprof[i], Fdiis );
    ofmo_daxpy( nao2, c[diis_now-1], Fp, Fdiis );

    /*// debug
    if ( debug_flag ) {
	printf("@@ check sum of F = %16.10f\n", check_sum( nao2, Fp ) );
	fflush(stdout);
    }*/

    return 0;
}

/** DIISの初期化関数
 *
 * DIISの初期化関数。
 * 重なり行列をCholesky分解して、その結果を保持する。
 * メモリ確保などの処理も行う。
 *
 * @attention
 * - 新たにSCF計算を始める際に、一回呼び出す必要がある。FMO計算のように
 *   複数回SCF計算を行う場合には、SCF計算を行う度に呼び出す必要がある。
 *
 * @param[in] nao AO数
 * @param[in] Sp[] 重なり行列（圧縮U形式）
 *
 * @retval  0 正常終了（メモリ確保、Cholesky分解が正常に出来た）
 * @retval -1 異常終了（メモリ確保を失敗した、あるいは、与えられた
 *     重なり行列のCholesky分解に失敗した）
 *
 * @ingroup ofmo-diis
 * */
int ofmo_diis_init( const int nao, const double Sp[] ) {
    int maxdiis = -1, maxnao = -1, ival;
    char *p;
    static int called = false;
#ifdef OFMO_MAXDIIS
    if ( called == false ) maxdiis = OFMO_MAXDIIS;
#endif
#ifdef OFMO_MAXNAO
    if ( called == false ) maxnao = OFMO_MAXNAO;
#endif
    if ( called == false && (p=getenv("OFMO_MAX_DIIS")) != NULL ) {
	ival = atoi( p );
	maxdiis = ( ival > maxdiis ? ival : maxdiis );
    }
    if ( called == false && (p=getenv("OFMO_MAX_NAO"))  != NULL ) {
	ival = atoi( p );
	maxnao = ( ival > maxnao ? ival : maxnao );
    }
    maxnao = ( nao > maxnao ? nao : maxnao );

    if ( called == false ) {
	if ( maxnao < 1 ) {
	    dbg("ERROR: Illegal maxnao (%d)\n", maxnao );
	    return -1;
	}
	if ( maxdiis < 1 ) {
	    dbg("ERROR: Illegal maxdiis (%d)\n", maxdiis );
	    return -1;
	}
    }
    ofmo_diis_alloc( maxdiis, maxnao );
    called = true;

    NAONOW = nao;
    double *U=T1, *ev=T2;
    ofmo_unpack_matrix( nao, Sp, U );
    if ( ofmo_solv_SEP( nao, U, ev ) != 0 ) {
	dbg("error\n"); fflush(stdout);
	return -1;
    }
    int loc=-1;
    for ( int i=0; i<nao; i++ ) {
	if ( ev[i] < 1.e-6 ) {
	    loc = i;
	    break;
	} else {
	    ev[i] = sqrt( 1.e0 / ev[i] );
	}
    }
    if ( loc >= 0 ) {
	dbg("error\n"); fflush(stdout);
	return -1;
    }
    for ( int i=0, k=0; i<nao; i++, k+=nao ) {
	for ( int j=0; j<nao; j++ ) {
	    S12[k+j] = U[k+j] * ev[i];
	}
    }
    ofmo_unpack_matrix( nao, Sp, S );
    /*ofmo_unpack_matrix( nao, Sp, S12 );
    info = ofmo_chodec( nao, S12 );*/
    diis_now = 0;
    //if ( info != 0 ) dbg("Failure in Cholesky decomposition\n");
    return 0;
}

/** Fock行列に対してDIISを適用する
 *
 * Fock行列に対してDIISを適用する。与えるフラグによっては、
 * 履歴を保存するだけである。
 *
 * @param[in] nao AO数
 * @param[in] Dp[] 密度行列（圧縮U形式）
 * @param[in,out] Fp[] Fock行列（圧縮U形式）。入力時にはDIISでの更新前の
 *     素のFock行列。出力時にはDIISでの更新を行ったFock行列。
 *     ただし、\c dodiis の値が0の場合には、入力と出力は同じもの。
 * @param[in] dodiis DIIS更新を行うかどうかを表すフラグ。0の場合には、
 *     与えられた密度行列とFock行列を用いた履歴情報を保存するだけで、
 *     Fock行列のDIIS更新は行わない。非0の場合には、履歴情報の保存に
 *     加えて、Fock行列のDIIS更新を行う。
 *     
 * @attention
 *   入力される密度行列、ならびに、Fock行列は、対角要素に
 *   細工されていたらダメ
 *   （電子状態計算プログラムでは、密度行列やFock行列の対角要素を
 *   2倍（あるいは1/2倍）してあることがよくあるので、注意すること）
 *
 * @ingroup ofmo-diis
 * */
int ofmo_diis_update( const int nao, double Dp[], const double Fp[],
	double Fdiis[], const int dodiis ) {
    int info = 0;
    if ( nao != NAONOW ) {
	dbg("Inconsistent in # of AO ( %d vs. %d )\n",
		NAONOW, nao );
	return -1;
    }
    //diis_profiling( nao, Dp, Fp );
    if ( dodiis ) info = diis_update( nao, Dp, Fp, Fdiis );
    return info;
}

