/**
 * @defgroup basis 基底関数の代入などに関わる関数群
 *
 * このクラス（C++のクラスとは違うけど・・・）では、与えられた
 * 組成情報（原子数、原子番号）と基底関数名を元に、基底関数データを
 * 作成する関数を定義している。また、分子積分プログラムの書きやすさを
 * 考慮して、通常の並びの基底関数データだけでなく、軌道量子数の大きさで
 * 並べ替えた基底関数データも扱う。
 *
 * @section unsort_params 非ソート基底関数データ
 * 非ソート基底関数データでは、以下のデータを扱う
 * 
 * @arg \c ushel_lqn[ics] CS番号 \c ics のCSの軌道量子数
 * @arg \c ushel_tem[ics] CS番号 \c ics のCSの縮約長
 * @arg \c ushel_atm[ics] CS番号 \c ics のCSが属する原子の番号
 * @arg \c ushel_add[ics] CS番号 \c ics のCSに含まれるPSの先頭PS番号。
 *     すなわち、CS番号 \c ics のCSに含まれるPSの先頭PSのデータが代入
 *     されている場所（オフセット値）。
 * @arg \c ushel_ini[ics] CS番号 \c ics のCSに含まれるAOの先頭AO番号。
 *     内部ではpure Cartessianとして計算してある。
 * @arg \c uprim_exp[ips] PS番号\c ips のPSの軌道指数
 * @arg \c uprim_coe[ips] PS番号\c ips のPSの規格化定数込みの縮約係数
 *
 * @subsection unsort_exam 非ソート使用例
 * 非ソート基底関数を用いた
 * １電子積分計算を行う場合のループの回し方（概要）
 * @include unsort-oneint.c
 * @section sort_params ソート基底関数データ
 * ソート基底関数データでは、以下のデータを扱う
 *
 * @arg \c maxlqn 最大軌道量子数
 * @arg \c leading_cs[lqn] 軌道量子数 \c lqn の先頭CS番号
 * @arg \c shel_tem[ics] CS番号 \c ics のCSの縮約長
 * @arg \c shel_atm[ics] CS番号 \c ics のCSが属する原子の番号
 * @arg \c shel_add[ics] CS番号 \c ics のCSに含まれるPSの先頭PS番号
 * @arg \c shel_ini[ics] CS番号 \c ics のCSに含まれるAOの先頭AO番号。
 *     内部ではpure Cartessianとして計算してある。
 * @arg \c prim_exp[ips] PS番号 \c ips のPSの軌道指数
 * @arg \c prim_coe[ips] PS番号 \c ips のPSの規格化定数込みの縮約係数
 *
 * @subsection sort_exam ソート基底関数の使用例
 * ソート基底関数を用いた
 * １電子積分計算を行う場合のループの回し方（概要）
 * @include sort-oneint.c
 *
 * @ingroup ofmo
 * */
/**
 * @file ofmo-basis.c
 *
 * 基底関数データの代入を行う関数群を定義しているファイル。
 *
 * 通常並びの基底関数データ、および、軌道量子数の大きさで並び替えた
 * 基底関数データを扱う。
 * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "ofmo-def.h"
#include "ofmo-basis-database.h"

/** すべて小文字に変換する関数
 *
 * \note
 * 文字列の長さのチェックは行っていない
 * */
static int Tolower( char str[] ) {
    int len, i;
    len = strlen( str );
    for ( i=0; i<len; i++ ) str[i] = tolower( str[i] );
    return 0;
}

/** 同じ基底関数名かチェックする関数
 *
 * 与えられた２つの基底関数名が同じか、違うかを判別する関数。
 * 内部ですべて小文字に変換して比較しているだけで、意味の比較
 * （例えば、6-31G**と6-31G(d,p)が等価であるなど）は行っていない。
 *
 * @param[in] bs1[] 基底関数名１
 * @param[in] bs2[] 基底関数名２
 *
 * @retval 0 1つの基底関数名が異なる
 * @retval 1 1つの基底関数名が一致
 *
 * */
static int is_same_basis_name( char bs1[], char bs2[] ) {
    int len1, len2;
    char BS1[MAXSTRLEN], BS2[MAXSTRLEN];
    len1 = strlen( bs1 );
    len2 = strlen( bs2 );
    if ( len1 != len2 ) return 0;
    if ( len1 >= MAXSTRLEN ) {
	dbg("ERROR: Too long string (%d)\n", len1 );
	return 0;	// errorだけど0を返している
    }
    strcpy( BS1, bs1 );
    strcpy( BS2, bs2 );
    Tolower( BS1 );
    Tolower( BS2 );
    return  (strcmp( BS1, BS2)==0 ? 1 : 0 ) ;
}

static int *NNAO = NULL;

static void ofmo_finalize_nnao() {
    if( NNAO != NULL ) free( NNAO );
    NNAO = NULL;
}

/** CSに含まれるAO数の設定関数
 *
 * CSに含まれるAO数を設定している。2011/06/09現在では、pure Cartessian
 * ガウス型関数を想定している
 * */
static int ofmo_init_nnao() {
    int i;
    static int called = false;
    if ( called ) return 0;
    NNAO = (int*)malloc( sizeof(int) * (MAXLQN+1) );
    for ( i=0; i<=MAXLQN; i++ ) NNAO[i] = (i+1)*(i+2)/2;// Cartessian
    atexit( ofmo_finalize_nnao );
    called = true;
    return 0;
}

/**
 * \var NSBS
 *      この計算で用いている基底関数の種類
 * \var *BASIS_NAME[ibs]
 *     \c ibs 番目の基底関数名
 * \var EXIST_BASIS[ibs][atn-1]
 *     基底関数番号\c ibs 、原子番号\c atn の原子種番号 \c ats 。
 *     \c ATOM_ADD[ibs][] と\c ATOM_NCS[ibs][] のパラメタが代入されている
 *     場所を示す。
 * \var ATOM_ADD[ibs][ats]
 *     基底関数番号\c ibs 、原子種番号\c ats の先頭CSパラメタの
 *     オフセット値\c ics0
 * \var ATOM_NCS[ibs][ats]
 *     基底関数番号\c ibs 、原子種番号\c ats のCS数
 * \var SHEL_LQN[ibs][ics]
 *     基底関数番号\c ibs 、CS番号\c ics の軌道量子数
 * \var SHEL_TEM[ibs][ics]
 *     基底関数番号\c ibs 、CS番号\c ics の縮約長
 * \var NNAO[lqn]
 *     軌道量子数\c lqn のCSに含まれるAO数
 * \var SHEL_ADD[ibs][ics]
 *     基底関数番号\c ibs 、CS番号\c ics に
 *     対応する先頭PSパラメターのオフセット値\c ips0
 * \var PRIM_EXP[ibs][ips]
 *     基底関数番号\c ibs 、PS番号\c ips の軌道指数
 * \var PRIM_COE[ibs][ips]
 *     基底関数番号\c ibs 、PS番号\c ips の規格化定数込みの縮約係数 
 *
 * */
static int CALLED        = false;
static int NSBS          = 0;
static char **BASIS_NAME = NULL;
static int **EXIST_BASIS = NULL;
static int **ATOM_ADD    = NULL;
static int **ATOM_NCS    = NULL;
static int **SHEL_LQN    = NULL;
static int **SHEL_TEM    = NULL;
static int **SHEL_ADD    = NULL;
static double **PRIM_EXP = NULL;
static double **PRIM_COE = NULL;

static void ofmo_finalize_database() {
    if ( CALLED ) {
	if ( BASIS_NAME != NULL ) {
	    if ( BASIS_NAME[0] != NULL ) free( BASIS_NAME[0] );
	    Free( BASIS_NAME );
	}
	Free( EXIST_BASIS );
	Free( ATOM_ADD );
	Free( ATOM_NCS );
	Free( SHEL_LQN );
	Free( SHEL_TEM );
	Free( SHEL_ADD );
	Free( PRIM_EXP );
	Free( PRIM_COE );
	NSBS = 0;
	BASIS_NAME  = NULL;
	EXIST_BASIS = NULL;
	ATOM_ADD    = NULL;
	ATOM_NCS    = NULL;
	SHEL_LQN    = NULL;
	SHEL_TEM    = NULL;
	SHEL_ADD    = NULL;
	PRIM_EXP    = NULL;
	PRIM_COE    = NULL;
    }
}

/** 計算で使う予定の基底関数データベースを作成する
 * */
static int ofmo_init_database(
	const int nsbs, char *basis_name[] ) {
    int ibs, ierr;
    if ( CALLED ) {
	if ( nsbs != NSBS ) {
	    dbg("ERROR: Inconsistency in # of kind of basis sets"
		    " %d vs %d\n", nsbs, NSBS );
	    return -1;
	}
	for ( ibs=0; ibs<nsbs; ibs++ ) {
	    if ( !is_same_basis_name(BASIS_NAME[ibs], basis_name[ibs]) ) {
		dbg("ERROR: Inconsistency in basis name (%s vs %s)\n",
			BASIS_NAME[ibs], basis_name[ibs] );
		return -1;
	    }
	}
	return 0;
    } else {
	NSBS = nsbs;
	// basis name
	BASIS_NAME = (char**)malloc( sizeof(char*) * nsbs );
	BASIS_NAME[0] = (char*)malloc( sizeof(char) * MAXSTRLEN * nsbs );
	for ( ibs=1; ibs<nsbs; ibs++ )
	    BASIS_NAME[ibs] = BASIS_NAME[ibs-1] + MAXSTRLEN;
	for ( ibs=0; ibs<nsbs; ibs++ ) {
	    if ( strlen(basis_name[ibs]) >= MAXSTRLEN ) {
		dbg("ERROR: Too long basis name (%s)\n", basis_name[ibs] );
		return -1;
	    }
	    strcpy( BASIS_NAME[ibs], basis_name[ibs] );
	}
	//
	EXIST_BASIS = (int**)malloc( sizeof(int*) * nsbs );
	ATOM_ADD    = (int**)malloc( sizeof(int*) * nsbs );
	ATOM_NCS    = (int**)malloc( sizeof(int*) * nsbs );
	SHEL_LQN    = (int**)malloc( sizeof(int*) * nsbs );
	SHEL_TEM    = (int**)malloc( sizeof(int*) * nsbs );
	SHEL_ADD    = (int**)malloc( sizeof(int*) * nsbs );
	PRIM_EXP    = (double**)malloc( sizeof(double*) * nsbs );
	PRIM_COE    = (double**)malloc( sizeof(double*) * nsbs );
	for ( ibs=0; ibs<nsbs; ibs++ ) {
	    ierr = ofmo_get_basis_data( basis_name[ibs],
		    "exist atom_add atom_ncs shel_lqn shel_tem shel_add "
		    "prim_exp prim_coe",
		    &EXIST_BASIS[ibs], &ATOM_ADD[ibs], &ATOM_NCS[ibs],
		    &SHEL_LQN[ibs], &SHEL_TEM[ibs], &SHEL_ADD[ibs],
		    &PRIM_EXP[ibs], &PRIM_COE[ibs] );
	    if ( ierr != 0 ) {
		dbg("ERROR: No parameter is found for basis (%s)\n",
			basis_name[ibs] );
		return -1;
	    }
	}
	// 規格化
	double pi, cons, pow2[MAXLQN+1], rfact2[MAXLQN+1];
	int lqn, atn, ics0, ics1, ics, ats, ips, ips0, ips1;
	double coef, lqn2;
	pow2[0]  = 1.e0;
	rfact2[0] = 1.e0;
	for ( lqn=1; lqn<=MAXLQN; lqn++ ) {
	    pow2[lqn]   = pow2[lqn-1]*2.e0;
	    rfact2[lqn] = rfact2[lqn-1] *sqrt( 1.e0/(double)(2*lqn-1) );
	}
	pi   = 4.e0 * atan( 1.e0 );
	cons = pow( 2.e0/pi, 0.75e0 );
	for ( ibs=0; ibs<nsbs; ibs++ ) {
	    for ( atn=1; atn<MAXATOMICNUMBER; atn++ ) {
		if ( (ats=EXIST_BASIS[ibs][atn-1]) >= 0 ) {
		    ics0 = ATOM_ADD[ibs][ats];
		    ics1 = ics0 + ATOM_NCS[ibs][ats];
		    for ( ics=ics0; ics<ics1; ics++ ) {
			lqn = SHEL_LQN[ibs][ics];
			lqn2 = ( (double)lqn + 1.5e0) * .5e0;
			coef = cons * pow2[lqn] * rfact2[lqn];
			ips0 = SHEL_ADD[ibs][ics];
			ips1 = ips0 + SHEL_TEM[ibs][ics];
			for ( ips=ips0; ips<ips1; ips++ ) 
			    PRIM_COE[ibs][ips] *=
				( coef * pow(PRIM_EXP[ibs][ips], lqn2 ) );
		    }	// for ( ics; );
		}
	    }
	}
	atexit( ofmo_finalize_database );
	CALLED = true;
    }
    return 0;
}

static int ofmo_basis_init(
	const int nsbs, char *basis_name_list[] ) {
    int ierr;
    ierr  = ofmo_init_nnao();
    ierr += ofmo_init_database( nsbs, basis_name_list );
    return ierr;
}

/** 基底関数データのサイズを取得する関数
 *
 * 与えられた原子数、組成（原子番号リスト）、基底関数リストなどを元に
 * 基底関数データのサイズ（最大軌道量子数、CS数、AO数、PS数）を得る。
 *
 * @param[in] nat 原子数
 * @param[in] nsbs 用いる基底関数の種類
 * @param[in] *basis_name[] 基底関数名のリスト
 * @param[in] atomic_number[] 各原子の原子番号
 * @param[in] atom_basis[] 各原子の基底関数の種類。
 *            \f$ ( 0\le \tt{nat} < \tt{nsbs} ) \f$\n
 *            ただし、NULLポインタの場合には、basis_name[0]で示される
 *            基底関数1種類しか用いないと仮定する。
 * @param[out] *maxlqn 最大軌道量子数
 * @param[out] *ncs CS数
 * @param[out] *nao AO数
 * @param[out] *nps PS数
 *
 * @retval  0 正常終了
 * @retval -1 異常終了（対応する基底関数データがない場合など）
 *
 * \note
 *
 * \b この関数の外部で定義されていて、この関数内で利用している変数\n
 * \arg \c EXIST_BASIS[ibs][atn-1] 基底関数番号\c ibs 、原子番号\c atn の
 *       CS種番号。後述の\c ATOM_ADD[ibs][] と\c ATOM_NCS[ibs][] の
 *       パラメタが代入されている場所\c add 。
 * \arg \c ATOM_ADD[ibs][add] 基底関数番号\c ibs 、CS種番号\c add の
 *       先頭CSパラメタのオフセット値\c ics0
 * \arg \c ATOM_NCS[ibs][add] 基底関数番号\c ibs 、CS種番号\c add のCS数
 * \arg \c SHEL_LQN[ibs][ics] 基底関数番号\c ibs 、CS番号\c ics のCS数
 * \arg \c SHEL_TEM[ibs][ics] 基底関数番号\c ibs 、CS番号\c ics の縮約長
 * \arg \c NNAO[lqn] 軌道量子数\c lqn のCSに含まれるAO数
 *
 * \b この関数から呼ばれる内部関数
 *
 * ofmo_basis_init
 *
 * @ingroup basis
 * */
int ofmo_get_basis_size( const int nat, const int nsbs,
	char *basis_name_list[], const int atomic_number[],
	const int atom_basis[],
	int *maxlqn, int *ncs, int *nao, int *nps ) {
    int mlqn, cs, ao, ps;
    int iat, atn, bsn, add, ics0, ics1, ics, lqn;
    /*
    // debug
    printf("basis_name = %s\n", basis_name[0] );
    fflush(stdout);
    */

    if ( ofmo_basis_init(nsbs, basis_name_list) != 0 ) return -1;
    if ( nat < 1 ) {
	dbg("ERROR: Illegal number of atom (%d)\n", nat );
	fflush(stdout);
	return -1;
    }
    mlqn = -1;
    cs = ao = ps = 0;
    for ( iat=0; iat<nat; iat++ ) {
	atn = atomic_number[iat];
	bsn = ( atom_basis == NULL ? 0 : atom_basis[iat] );
	if ( (add = EXIST_BASIS[bsn][atn-1]) < 0 ) {
	    dbg("ERROR: No basis data (%s, atn=%d)\n",
		    basis_name_list[bsn], atn );
	    fflush(stdout);
	    return -1;
	}
	ics0 = ATOM_ADD[bsn][add];
	ics1 = ics0 + ATOM_NCS[bsn][add];
	cs  += ATOM_NCS[bsn][add];
	for ( ics=ics0; ics<ics1; ics++ ) {
	    lqn = SHEL_LQN[bsn][ics];
	    if ( lqn > mlqn ) mlqn = lqn;
	    ao += NNAO[lqn];
	    ps += SHEL_TEM[bsn][ics];
	}
    }
    *maxlqn = mlqn;
    *ncs    = cs;
    *nao    = ao;
    *nps    = ps;
    return 0;

}
/** 非ソート基底関数データの代入
 *
 * 与えられた原子数、原子組成、基底関数名などのデータを元に、
 * 非ソート基底関数データを対応する配列に代入する。
 * 代入先の配列の領域は、呼び出し時点で確保されている必要がある。
 *
 * @param[in] nat 原子数
 * @param[in] nsbs 用いる基底関数の種類
 * @param[in] *basis_name[] 基底関数名のリスト
 * @param[in] atomic_number[] 各原子の原子番号
 * @param[in] atom_basis[] 各原子の基底関数の種類。
 *            \f$ ( 0\le \tt{nat} < \tt{nsbs} ) \f$\n
 *            ただし、NULLポインタの場合には、basis_name[0]で示される
 *            基底関数1種類しか用いないと仮定する。
 * 
 * @param[out] ushel_lqn[] 軌道量子数
 * @param[out] ushel_tem[] 縮約長
 * @param[out] ushel_atm[] 属する原子の番号
 * @param[out] ushel_add[] 対応する先頭PS番号
 * @param[out] ushel_ini[] 先頭AO番号
 * @param[out] uprim_exp[] 軌道指数
 * @param[out] uprim_coe[] 規格化定数込みの縮約係数
 *
 * @retval  0 正常終了
 * @retval -1 異常終了（対応する基底関数パラメタがない場合など）
 *
 * \note
 *
 * \b この関数の外部で定義されていて、この関数内で利用している変数\n
 * \arg \c EXIST_BASIS[ibs][atn-1] 基底関数番号\c ibs 、原子番号\c atn の
 *       CS種番号。後述の\c ATOM_ADD[ibs][] と\c ATOM_NCS[ibs][] の
 *       パラメタが代入されている場所\c add 。
 * \arg \c ATOM_ADD[ibs][add] 基底関数番号\c ibs 、CS種番号\c add の
 *       先頭CSパラメタのオフセット値\c ics0
 * \arg \c ATOM_NCS[ibs][add] 基底関数番号\c ibs 、CS種番号\c add のCS数
 * \arg \c SHEL_LQN[ibs][ics] 基底関数番号\c ibs 、CS番号\c ics のCS数
 * \arg \c SHEL_TEM[ibs][ics] 基底関数番号\c ibs 、CS番号\c ics の縮約長
 * \arg \c SHEL_ADD[ibs][ics] 基底関数番号\c ibs 、CS番号\c ics に
 * 対応する先頭PSパラメターのオフセット値\c ips0 。
 * \arg \c PRIM_EXP[ibs][ips] 基底関数番号\c ibs 、PS番号\c ips の
 * 軌道指数
 * \arg \c PRIM_COE[ibs][ips] 基底関数番号\c ibs 、PS番号\c ips の
 * 規格化定数込みの縮約係数 
 * \arg \c NNAO[lqn] 軌道量子数\c lqn のCSに含まれるAO数
 *
 * \b この関数から呼ばれる内部関数
 *
 * ofmo_basis_init
 *
 * @ingroup basis
 * */
int ofmo_assign_basis(
	const int nat, const int nsbs, char *basis_name[],
	const int atomic_number[], const int atom_basis[],
	int ushel_lqn[], int ushel_tem[], int ushel_atm[],
	int ushel_add[], int ushel_ini[],
	double uprim_exp[], double uprim_coe[] ) {
    int ncs, nao, nps;
    int ics, ics0, ics1, ips, ips0, ips1;
    int iat, add, tem, lqn, atn, bsn;
    if ( ofmo_basis_init(nsbs, basis_name) != 0 ) return -1;
    if ( nat < 1 ) {
	dbg("ERROR: Illegal number of atom (%d)\n", nat );
	fflush(stdout);
	return -1;
    }
    ncs = nao = nps = 0;
    for ( iat=0; iat<nat; iat++ ) {
	atn = atomic_number[iat];
	bsn = ( atom_basis == NULL ? 0 : atom_basis[iat] );
	if ( (add = EXIST_BASIS[bsn][atn-1]) < 0 ) {
	    dbg("ERROR: No basis data (%s, atn=%d)\n",
		    basis_name[bsn], atn );
	    fflush(stdout);
	    return -1;
	}
	ics0 = ATOM_ADD[bsn][add];
	ics1 = ics0 + ATOM_NCS[bsn][add];
	for ( ics=ics0; ics<ics1; ics++ ) {
	    ushel_lqn[ncs] = lqn = SHEL_LQN[bsn][ics];
	    ushel_tem[ncs] = tem = SHEL_TEM[bsn][ics];
	    ushel_atm[ncs] = iat;
	    ushel_add[ncs] = nps;
	    ushel_ini[ncs] = nao;
	    ips0 = SHEL_ADD[bsn][ics];
	    ips1 = ips0 + tem;
	    for ( ips=ips0; ips<ips1; ips++ ) {
		uprim_exp[nps] = PRIM_EXP[bsn][ips];
		uprim_coe[nps] = PRIM_COE[bsn][ips];
		nps++;
	    }
	    nao += NNAO[lqn];
	    ncs++;
	}
    }
    return 0;
}

/** 基底関数データの並べ替え
 *
 * 非ソート基底関数をソート基底関数に並び替えるための関数。
 * 非ソート基底関数データは既に代入されている必要がある。また、ソート
 * 基底関数データを代入する配列は、この関数呼び出し時に、すでに確保
 * されている必要がある。
 *
 * @param[in] maxlqn 最大軌道量子数
 * @param[in] ncs CS数
 * @param[in] ushel_lqn[] 軌道量子数
 * @param[in] ushel_tem[] 縮約長
 * @param[in] ushel_atm[] 属する原子の番号
 * @param[in] ushel_add[] 対応する先頭PS番号
 * @param[in] ushel_ini[] 先頭AO番号
 * @param[in] uprim_exp[] 軌道指数
 * @param[in] uprim_coe[] 規格化定数込みの縮約係数
 *
 * @param[out] leading_cs[] 各軌道量子数を持つ先頭CS番号
 * @param[out] shel_tem[] 縮約長（ソート済み）
 * @param[out] shel_atm[] 属する原子の番号（ソート済み）
 * @param[out] shel_add[] 対応する先頭PS番号（ソート済み）
 * @param[out] shel_ini[] 先頭AO番号（ソート済み）
 * @param[out] prim_exp[] 軌道指数（ソート済み）
 * @param[out] prim_coe[] 規格化定数込みの縮約係数（ソート済み）
 * @param[out] s2u[] ソート後のAOの、ソート前AO番号
 *
 * \note 外部で定義されていて、関数内部で使っている変数
 *
 * \c NNAO[] 各軌道量子数を持つCSに含まれるAO数
 *
 * \b この関数から呼ばれる内部関数
 *
 * ofmo_basis_init
 *
 * @ingroup basis
 * */
int ofmo_sort_basis(
	const int maxlqn, const int ncs,
	const int ushel_lqn[], const int ushel_tem[],
	const int ushel_atm[], const int ushel_add[],
	const int ushel_ini[],
	const double uprim_exp[], const double uprim_coe[],
	int leading_cs[],
	int shel_tem[], int shel_atm[], int shel_add[], int shel_ini[],
	double prim_exp[], double prim_coe[], int s2u[]
	) {
    int lqn, ics, ips, ips0, ips1, tem, iao0, iao1, iao;
    int cs, ao, ps;
    int nnao;
    ofmo_init_nnao();
    cs = ao = ps = 0;
    for ( lqn=0; lqn<=maxlqn; lqn++ ) {
	leading_cs[lqn] = cs;
	nnao = NNAO[lqn];
	for ( ics=0; ics<ncs; ics++ ) {
	    if ( ushel_lqn[ics] == lqn ) {
		shel_tem[cs] = tem = ushel_tem[ics];
		shel_atm[cs] = ushel_atm[ics];
		shel_add[cs] = ps;
		shel_ini[cs] = ao;
		iao0         = ushel_ini[ics];
		iao1         = iao0 + nnao;
		ips0         = ushel_add[ics];
		ips1         = ips0 + tem;
		for ( iao=iao0; iao<iao1; iao++ ) {
		    s2u[ao] = iao;
		    ao++;
		}
		for ( ips=ips0; ips<ips1; ips++ ) {
		    prim_exp[ps] = uprim_exp[ips];
		    prim_coe[ps] = uprim_coe[ips];
		    ps++;
		}
		cs++;
	    }
	}
    }
    leading_cs[maxlqn+1] = cs;
    return 0;
}

/** 非ソート基底関数データを保存するための領域を確保する
 *
 * CS数 \c ncs 、PS数 \c nps
 * である非ソート基底関数のデータを保存するための領域を確保する関数。
 * 必要なすべての配列が確保できたら、正常終了して呼び出し元に戻る。
 * １つでも確保できない配列があった場合には、確保できた配列用のメモリを
 * すべて開放して、異常終了として呼び出し元に戻る。
 *
 * @param[in] ncs CS数
 * @param[in] nps PS数
 * @param[out] *(*ushel_lqn) 軌道量子数を保存する配列の先頭アドレス
 * @param[out] *(*ushel_tem) 縮約長を保存する配列の先頭アドレス
 * @param[out] *(*ushel_atm) 属する原子を保存する配列の先頭アドレス
 * @param[out] *(*ushel_add) 対応する先頭PS番号を保存する配列の先頭アドレス
 * @param[out] *(*ushel_ini) 先頭AO番号を保存する配列の先頭アドレス
 * @param[out] *(*uprim_exp) 軌道指数を保存する配列の先頭アドレス
 * @param[out] *(*uprim_coe) 縮約係数を保存する配列の先頭アドレス
 * 
 * @retval  0 正常終了（すべての配列を確保できた場合）
 * @retval -1 異常終了（配列の確保ができなかった場合）
 *
 * @ingroup basis
 * */
int ofmo_alloc_unsorted_basis( const int ncs, const int nps,
	int **ushel_lqn, int **ushel_tem, int **ushel_atm,
	int **ushel_add, int **ushel_ini,
	double **uprim_exp, double **uprim_coe ) {
    *ushel_lqn = (int*)malloc( sizeof(int) * ncs );
    *ushel_tem = (int*)malloc( sizeof(int) * ncs );
    *ushel_atm = (int*)malloc( sizeof(int) * ncs );
    *ushel_add = (int*)malloc( sizeof(int) * ncs );
    *ushel_ini = (int*)malloc( sizeof(int) * ncs );
    *uprim_exp = (double*)malloc( sizeof(double) * nps );
    *uprim_coe = (double*)malloc( sizeof(double) * nps );
    if (    *ushel_lqn == NULL || *ushel_tem == NULL ||
	    *ushel_atm == NULL || *ushel_add == NULL ||
	    *ushel_ini == NULL || *uprim_exp == NULL ||
	    *uprim_coe == NULL ) {
	if ( *ushel_lqn != NULL ) free ( *ushel_lqn );
	if ( *ushel_tem != NULL ) free ( *ushel_tem );
	if ( *ushel_atm != NULL ) free ( *ushel_atm );
	if ( *ushel_add != NULL ) free ( *ushel_add );
	if ( *ushel_ini != NULL ) free ( *ushel_ini );
	if ( *uprim_exp != NULL ) free ( *uprim_exp );
	if ( *uprim_coe != NULL ) free ( *uprim_coe );
	return -1;
    }
    return 0;
}

/** * ソート基底関数データを保存するための領域を確保する
 *
 * 最大軌道量子数 \c maxlqn 、CS数 \c ncs 、AO数 \c nao 、PS数 \c nps
 * であるソート基底関数のデータを保存するための領域を確保する関数。
 * 必要なすべての配列が確保できたら、正常終了して呼び出し元に戻る。
 * １つでも確保できない配列があった場合には、確保できた配列用のメモリを
 * すべて開放して、異常終了として呼び出し元に戻る。
 *
 * @param[in] maxlqn 最大軌道量子数
 * @param[in] ncs CS数
 * @param[in] nao AO数
 * @param[in] nps PS数
 * @param[out] *(*leading_cs) 先頭CS番号を記憶する配列の先頭アドレス
 * @param[out] *(*shel_tem) 縮約長を保存する配列の先頭アドレス
 * @param[out] *(*shel_atm) 属する原子を保存する配列の先頭アドレス
 * @param[out] *(*shel_add) 対応する先頭PS番号を保存する配列の先頭アドレス
 * @param[out] *(*shel_ini) 先頭AO番号を保存する配列の先頭アドレス
 * @param[out] *(*prim_exp) 軌道指数を保存する配列の先頭アドレス
 * @param[out] *(*prim_coe) 縮約係数を保存する配列の先頭アドレス
 * @param[out] *(*s2u) 元の並びのAO番号への変換パラメタ用の配列の
 *                     先頭アドレス
 * 
 * @retval  0 正常終了（すべての配列を確保できた場合）
 * @retval -1 異常終了（配列の確保ができなかった場合）
 *
 * @ingroup basis
 * */
int ofmo_alloc_sorted_basis(
	const int maxlqn, const int ncs, const int nao, const int nps,
	int **leading_cs, int **shel_tem, int **shel_atm, int **shel_add,
	int **shel_ini, double **prim_exp, double **prim_coe,
	int **s2u ) {
    *leading_cs = (int*)malloc( sizeof(int) * (maxlqn + 1 + 1 ) );
    *shel_tem = (int*)malloc( sizeof(int) * ncs );
    *shel_atm = (int*)malloc( sizeof(int) * ncs );
    *shel_add = (int*)malloc( sizeof(int) * ncs );
    *shel_ini = (int*)malloc( sizeof(int) * ncs );
    *prim_exp = (double*)malloc( sizeof(double) * nps );
    *prim_coe = (double*)malloc( sizeof(double) * nps );
    *s2u      = (int*)malloc( sizeof(int) * nao );
    if (    *leading_cs == NULL || *s2u == NULL ||
	    *shel_tem == NULL || *shel_atm == NULL ||
	    *shel_add == NULL || *shel_ini == NULL ||
	    *prim_exp == NULL || *prim_coe == NULL ) {
	if ( *leading_cs != NULL ) free( *leading_cs );
	if ( *shel_tem != NULL ) free( *shel_tem );
	if ( *shel_atm != NULL ) free( *shel_atm );
	if ( *shel_add != NULL ) free( *shel_add );
	if ( *shel_ini != NULL ) free( *shel_ini );
	if ( *prim_exp != NULL ) free( *prim_exp );
	if ( *prim_coe != NULL ) free( *prim_coe );
	if ( *s2u      != NULL ) free( *s2u );
	return -1;
    }
    return 0;
}

/** 非ソート基底関数データ出力関数（デバッグ出力用？）
 *
 * 非ソート基底関数データを出力する関数。デバッグ用かな？
 *
 * @param[in] *fp 出力先のファイル記述子
 * @param[in] ncs CS数
 * @param[in] ushel_lqn[] 軌道量子数
 * @param[in] ushel_tem[] 縮約長
 * @param[in] ushel_atm[] 属する原子の番号
 * @param[in] ushel_add[] 対応するPSの先頭PS番号
 * @param[in] ushel_ini[] 先頭AO番号
 * @param[in] uprim_exp[] 軌道指数
 * @param[in] uprim_coe[] 規格化定数込みの縮約係数
 *
 * @return なし
 *
 * @ingroup basis
 *
 * */
void ofmo_show_unsorted_basis_params( FILE *fp, const int ncs,
	const int ushel_lqn[], const int ushel_tem[],
	const int ushel_atm[], const int ushel_add[],
	const int ushel_ini[],
	const double uprim_exp[], const double uprim_coe[] ) {
    int ics, ips, ips0, ips1;
    fprintf(fp, "======= UNSORTED BASIS SET PARAMETERS ========\n");
    fprintf(fp, "ICS LQN TEM ATM ADD INI: EXP0 COE0 EXP1 COE1 ...\n");
    for ( ics=0; ics<ncs; ics++ ) {
	fprintf(fp,"%3d %3d %3d %3d %3d %3d:", ics,
		ushel_lqn[ics], ushel_tem[ics], ushel_atm[ics],
		ushel_add[ics], ushel_ini[ics] );
	ips0 = ushel_add[ics];
	ips1 = ips0 + ushel_tem[ics];
	for ( ips=ips0; ips<ips1; ips++ )
	    fprintf(fp," %10.3e %10.2e", uprim_exp[ips], uprim_coe[ips] );
	fprintf(fp, "\n");
    }
    fprintf(fp, "==============================================\n");
    fflush( fp );
}

/** ソート基底関数データ出力関数（デバッグ出力用？）
 *
 * @param[in] *fp 出力先のファイル記述子
 * @param[in] maxlqn 最大軌道量子数
 * @param[in] nao AO数
 * @param[in] leading_cs[lqn] 軌道量子数 \c lqn の先頭CS番号
 * @param[in] shel_tem[] 縮約長（ソート済み）
 * @param[in] shel_atm[] 属する原子の番号（ソート済み）
 * @param[in] shel_add[] 対応するPSの先頭PS番号（ソート済み）
 * @param[in] shel_ini[] 先頭AO番号（ソート済み）
 * @param[in] prim_exp[] 軌道指数（ソート済み）
 * @param[in] prim_coe[] 規格化定数込みの縮約係数（ソート済み）
 * @param[in] s2u[] 元の並びに於けるAO番号
 *
 * @return なし
 *
 * @ingroup basis
 * */
void ofmo_show_sorted_basis_params( FILE *fp,
	const int maxlqn, const int nao, const int leading_cs[],
	const int shel_tem[], const int shel_atm[],
	const int shel_add[], const int shel_ini[],
	const double prim_exp[], const double prim_coe[],
	const int s2u[] ) {
    int lqn, ics, ics0, ics1, ips, ips0, ips1, iao;
    fprintf(fp, "======= SORTED BASIS SET PARAMETERS ========\n");
    fprintf(fp, "ICS TEM ATM ADD INI: EXP0 COE0 EXP1 COE1 ...\n");
    for ( lqn=0; lqn<=maxlqn; lqn++ ) {
	ics0 = leading_cs[lqn];
	ics1 = leading_cs[lqn+1];
	fprintf(fp, "lqn=%d\n", lqn );
	for ( ics=ics0; ics<ics1; ics++ ) {
	    fprintf(fp, "%3d %3d %3d %3d %3d:", ics,
		    shel_tem[ics], shel_atm[ics], shel_add[ics],
		    shel_ini[ics] );
	    ips0 = shel_add[ics];
	    ips1 = ips0 + shel_tem[ics];
	    for ( ips=ips0; ips<ips1; ips++ )
		fprintf(fp," %10.3e %10.2e",
			prim_exp[ips], prim_coe[ips] );
	    fprintf(fp, "\n");
	}
    }
    fprintf(fp, "-- AO conversion ( sorted -> unsorted ) --");
    for ( iao=0; iao<nao; iao++ )
	fprintf(fp, "%s%3d -> %3d",
		( (iao%5) == 0 ? "\n  " : ", "), iao, s2u[iao] );
    if ( (nao%5) != 0 ) fprintf(fp, "\n");
    fprintf(fp, "============================================\n");
    fflush( fp );
}

/** 圧縮U形式の対称行列の要素を並べ替える
 *
 * 圧縮U形式で保存された対称行列の要素の並びを並び替える。
 *
 * 分子積分コードの単純化のために、基底関数を軌道量子数の大きさで
 * 並び替えて、その並び替えた基底関数データ（ソート基底関数データ）を
 * 用いて、各種分子積分を行っている。この基底関数データの並び替えを
 * 行うと、AOの並びもソートされるようになっている。
 * しかし、並び替えられていない密度行列データを外部から読み込んで
 * 計算を行うなどの場合には、ソートされていない（非ソート）データを
 * ソートしたデータに変換する必要がある。また、ソートされた状態で
 * 計算された密度行列データを、他のプログラムで利用するために出力する、
 * といった場合には、ソートされたデータを非ソートデータに並び替える
 * 必要がある。
 *
 * この関数は、ソートされたの圧縮U形式対称行列の要素を
 * 非ソート（元の並び）の圧縮U形式対称行列に並び替える関数である。
 *
 * @param[in] nao AO数（＝行列サイズ）
 * @param[in] s2u[iao] ソート基底関数データにおける \c iao 番目のAOの
 *     非ソート基底関数データにおけるAO番号。基底関数の並べ替え関数
 *     \c ofmo_sort_basis を呼び出すと得られる。
 * @param[in] src[nao*(nao+1)/2] ソートされた
 *     圧縮U形式対称行列要素が代入されている配列
 * @param[out] dest[nao*(nao+1)/2] ソートされていない（元の並びの）
 *     圧縮U形式対称行列要素が代入される配列
 *
 * @return なし
 *
 * @ingroup basis
 *
 * */
void ofmo_unsort_packed_matrix(const int nao, const int s2u[],
	const double src[], double dist[]) {
    int i, j, c, I,J,IJ;
    c = 0;
    for ( i=0; i<nao; i++) {
	I = s2u[i];
	for ( j=0; j<=i; j++) {
	    J = s2u[j];
	    IJ = (I>=J ? I*(I+1)/2+J : J*(J+1)/2+I);
	    dist[IJ] = src[c];
	    c++;
	}
    }
}

/** 圧縮U形式の対称行列の要素を並べ替える
 *
 * 圧縮U形式で保存された対称行列の要素の並びを並び替える。
 *
 * 分子積分コードの単純化のために、基底関数を軌道量子数の大きさで
 * 並び替えて、その並び替えた基底関数データ（ソート基底関数データ）を
 * 用いて、各種分子積分を行っている。この基底関数データの並び替えを
 * 行うと、AOの並びもソートされるようになっている。
 * しかし、並び替えられていない密度行列データを外部から読み込んで
 * 計算を行うなどの場合には、ソートされていない（非ソート）データを
 * ソートしたデータに変換する必要がある。また、ソートされた状態で
 * 計算された密度行列データを、他のプログラムで利用するために出力する、
 * といった場合には、ソートされたデータを非ソートデータに並び替える
 * 必要がある。
 *
 * この関数は、非ソート状態（元の並び）の圧縮U形式対称行列の要素を
 * ソートされた圧縮U形式対称行列に並び替える関数である。
 *
 * @param[in] nao AO数（＝行列サイズ）
 * @param[in] s2u[iao] ソート基底関数データにおける \c iao 番目のAOの
 *     非ソート基底関数データにおけるAO番号。基底関数の並べ替え関数
 *     \c ofmo_sort_basis を呼び出すと得られる。
 * @param[in] src[nao*(nao+1)/2] 元の並びの（並び替えを行っていない）
 *     圧縮U形式対称行列要素が代入されている配列
 * @param[out] dest[nao*(nao+1)/2] ソートされた圧縮U形式対称行列要素が
 *     代入される配列
 *
 * @return なし
 *
 * @ingroup basis
 *
 * */
void ofmo_sort_packed_matrix(const int nao, const int s2u[],
	const double src[], double dest[]) {
    int i, j, c, I,J,IJ, I2;
    c = 0;
    for ( i=0; i<nao; i++) {
	I = s2u[i];
	I2 = I*(I+1)/2;
	for ( j=0; j<=i; j++) {
	    J = s2u[j];
	    IJ = (I>=J ? I2+J : J*(J+1)/2+I);
	    dest[c] = src[IJ];
	    c++;
	}
    }
}
/** ベクトルをソートする
 *
 * 非ソートのベクトル（たとえば、MO係数行列のある列ベクトル）をソートする
 *
 * @param[in] nao AO数（＝ベクトルサイズ）
 * @param[in] sao2uao[iao] ソート基底関数データにおける \c iao 番目のAOの
 *     非ソート基底関数データにおけるAO番号。基底関数の並べ替え関数
 *     \c ofmo_sort_basis を呼び出すと得られる。
 * @param[in] v0[nao] 元の並びのベクトル
 * @param[out] vs[nao] ソートしたベクトル
 *
 * @return なし
 *
 * @ingroup basis
 *
 * */
void ofmo_sort_vector( const int nao, const int sao2uao[],
	const double v0[], double vs[] ) {
    int i, I;
    for ( i=0; i<nao; i++ ) {
	I = sao2uao[i];
	vs[i] = v0[I];
    }
}

/* \mainpage 基底関数代入に関する関数群
 *
 * \section コンパイル方法
 *
 * コンパイルには、
 * \li 実行可能なテストプログラムを作成する
 * \li ライブラリを作成する
 *
 * という、２つのモードがある。このモードは、makefileの
 * 先頭付近にある \c TEST_MAIN 変数で切り替えることができる。
 *
 * \c TEST_MAIN=YES にすると、実行可能なテストプログラム \c test-basis が
 * 作成される。これは、このコードの開発時に用いたもので、通常、この
 * モードでのコンパイルは必要ない。
 *
 * \c TEST_MAIN=NO にすると、ライブラリファイル \c libofmo-basis.a が
 * 作成される。これは、他のオブジェクトとリンクして使用することを
 * 目的としている。
 *
 * 
 * */

#ifdef TEST_MAIN
// テスト用メイン関数
int main(int argc, char** argv) {
    // input
    int iat, nat = 8;
    int atomic_number[8] = { 1, 6, 7, 8, 1, 6, 8, 7 };
    int atom_basis[8]    = { 0, 0, 0, 0, 1, 1, 1, 1 };
    int ibs, nsbs = 2;
    char *basis_name[] = {"STO-3G", "6-31G*" };
    // output
    int ncs, nao, nps, maxlqn;
    int *ushel_lqn, *ushel_tem, *ushel_atm;
    int *ushel_add, *ushel_ini;
    double *uprim_exp, *uprim_coe;

    printf("------------------- INPUT DATA ------------------------\n");
    printf("# of kind of basis set = %d\n", nsbs );
    printf("basis name list = ");
    for ( ibs=0; ibs<nsbs; ibs++ ) printf(" %s", basis_name[ibs] );
    printf("\n");
    printf("# of atom = %d\n", nat);
    printf("atomic_number[] =");
    for (iat=0; iat<nat; iat++ )
	printf("%s%2d", ( (iat%10)==0 ? "\n  " : " " ),
		atomic_number[iat] );
    if ( (nat%10) != 0 ) printf("\n");
    printf("atom_basis[] =");
    for (iat=0; iat<nat; iat++ )
	printf("%s%2d", ( (iat%10)==0 ? "\n  " : " " ),
		atom_basis[iat] );
    if ( (nat%10) != 0 ) printf("\n");
    printf("-------------------------------------------------------\n");

    if ( ofmo_get_basis_size(
		nat, nsbs, basis_name, atomic_number, atom_basis,
		&maxlqn, &ncs, &nao, &nps ) != 0 ) {
	dbg("error: Failure in ofmo_get_basis_size\n");
	return -1;
    }
    if ( ofmo_alloc_unsorted_basis(ncs, nps, &ushel_lqn, &ushel_tem,
		&ushel_atm, &ushel_add, &ushel_ini, &uprim_exp, &uprim_coe)
	    != 0 ) {
	dbg("error: Failure in memory allocation for unsorted basis\n");
	return -1;
    }
    if ( ofmo_assign_basis(
		nat, nsbs, basis_name, atomic_number, atom_basis,
		ushel_lqn, ushel_tem, ushel_atm, ushel_add, ushel_ini,
		uprim_exp, uprim_coe ) != 0 ) {
	dbg("error: Failure in assignment basis parameters\n");
	return -1;
    }
    ofmo_show_unsorted_basis_params( stdout, ncs,
	    ushel_lqn, ushel_tem, ushel_atm, ushel_add, ushel_ini,
	    uprim_exp, uprim_coe );
    int *shel_tem, *shel_atm, *shel_add, *shel_ini, *leading_cs, *s2u;
    double *prim_exp, *prim_coe;
    if ( ofmo_alloc_sorted_basis(maxlqn, ncs, nao, nps,
		&leading_cs, &shel_tem, &shel_atm, &shel_add, &shel_ini,
		&prim_exp, &prim_coe, &s2u ) != 0 ) {
	dbg("error: Failure in memory allocation for sorted basis\n");
	return -1;
    }
    if ( ofmo_sort_basis( maxlqn, ncs,
		ushel_lqn, ushel_tem, ushel_atm, ushel_add, ushel_ini,
		uprim_exp, uprim_coe,
		leading_cs, shel_tem, shel_atm, shel_add, shel_ini,
		prim_exp, prim_coe, s2u ) != 0 ) {
	dbg("error: failure in sort basis parameters\n");
	return -1;
    }
    ofmo_show_sorted_basis_params( stdout, maxlqn, nao,
	    leading_cs, shel_tem, shel_atm, shel_add, shel_ini,
	    prim_exp, prim_coe, s2u );
    return 0;
}
#endif
