/**
 * @file ofmo-fragment-basis.c
 *
 * フラグメントのソート基底関数を代入するための関数を定義した
 * ファイル
 *
 * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ofmo-def.h"
#include "ofmo-data.h"
#include "ofmo-basis.h"

static int *fatm_lcs = NULL;

static int *fushel_lqn    = NULL;
static int *fushel_tem    = NULL;
static int *fushel_atm    = NULL;
static int *fushel_add    = NULL;
static int *fushel_ini    = NULL;
static double *fuprim_exp = NULL;
static double *fuprim_coe = NULL;

/* 一時データ領域の開放 */
static void dealloc() {
    Free( fatm_lcs );
    Free( fushel_lqn );
    Free( fushel_tem );
    Free( fushel_atm );
    Free( fushel_add );
    Free( fushel_ini );
    Free( fuprim_exp );
    Free( fuprim_coe );
    fatm_lcs     = NULL;
    fushel_lqn   = NULL;
    fushel_tem   = NULL;
    fushel_atm   = NULL;
    fushel_add   = NULL;
    fushel_ini   = NULL;
    fuprim_exp   = NULL;
    fuprim_coe   = NULL;
}

/* 一時データ領域の確保 */
static int alloc() {
    static int called = false;
    if ( called ) return 0;
    int ierr, nbody, maxnfatom, maxnfcs, maxnfps;
    int NAT, NCS, NPS;
    ierr = ofmo_data_get_vals("nbody maxnfatom maxnfcs maxnfps",
	    &nbody, &maxnfatom, &maxnfcs, &maxnfps );
    if ( ierr != 0 ) return -1;
    NAT = nbody * maxnfatom;
    NCS = nbody * maxnfcs;
    NPS = nbody * maxnfps;
    fatm_lcs = (int*)malloc( sizeof(int) * (NAT+1) );
    fushel_lqn = (int*)malloc( sizeof(int) * NCS );
    fushel_tem = (int*)malloc( sizeof(int) * NCS );
    fushel_atm = (int*)malloc( sizeof(int) * NCS );
    fushel_add = (int*)malloc( sizeof(int) * NCS );
    fushel_ini = (int*)malloc( sizeof(int) * NCS );
    fuprim_exp = (double*)malloc( sizeof(double) * NPS );
    fuprim_coe = (double*)malloc( sizeof(double) * NPS );
    atexit( dealloc );
    called = true;
    return 0;
}

/** フラグメントのソート基底関数データを代入する
 *
 * フラグメント（モノマー、ダイマーなど）のソート基底関数データを
 * 代入する。また、フラグメント内のAO（ソート、非ソート）番号と
 * 分子全体の非ソートAO番号の対応表も作成する。
 *
 * @attention
 * - モノマーデータの代入で用いてもよいが、少しだけ
 * （ofmo_get_basis_size関数呼び出しの分）オーバーヘッドがかかる
 * - この関数を呼び出す前に、分子全体の基底関数情報代入関数
 *   \c ofmo_set_basis_for_entire_molecule が呼び出されて、以下の
 *   変数が設定されている必要がある
 *     -# \c atm_lcs[iatm] 分子全体で \c iatm 番目の原子の
 *         （非ソート基底関数における）先頭CS番号
 *     -# \c ushel_ini[ics] 分子全体 \c ics 番目のCSの
 *         （非ソート基底関数における）先頭AO番号
 * - 出力用の配列変数（\c feading_cs や \c fshel_tem など）の領域は
 *   事前に確保してある必要がある（この関数内部では確保しない）
 *
 * @param[in] nat フラグメントの原子数
 * @param[in] nsbs 用いる基底関数の種類
 * @param[in] *basis_name_list[] 基底関数名のリスト
 * @param[in] atomic_number[] フラグメントを構成する原子の元の原子番号
 * @param[in] atom_basis[] 各原子の基底関数の種類
 *         \f$ ( 0\le \tt{nat} < \tt{nsbs} ) \f$。\n
 *         ただし、NULLポインタの場合には、basis_name[0]で示される
 *         基底関数1種類しか用いないと仮定する。
 * @param[in] fatom2tatom[ifat] フラグメント内の \c ifat 番目の原子の
 *     分子全体でのシリアル番号
 *
 * @param[out] *fncs CS数
 * @param[out] *fnao AO数
 * @param[out] *fnps PS数
 * @param[out] *fnpspair PSペア数
 * @param[out] fleading_cs[] 各軌道量子数を持つ先頭CS番号
 * @param[out] fshel_tem[] 縮約長（ソート済み）
 * @param[out] fshel_atm[] 属する原子の番号（ソート済み）
 * @param[out] fshel_add[] 対応する先頭PS番号（ソート済み）
 * @param[out] fshel_ini[] 先頭AO番号（ソート済み）
 * @param[out] fprim_exp[] 軌道指数（ソート済み）
 * @param[out] fprim_coe[] 規格化定数込みの縮約係数（ソート済み）
 * @param[out] fsao2fuao[ifao] フラグメントのソートAO番号 \c ifao のAOの
 *     フラグメントの非ソートAO番号
 * @param[out] fuao2tuao[ifao] フラグメントの非ソートAO番号 \c ifao のAOの
 *     分子全体における非ソートAO番号
 * @param[out] fsao2tuao[ifao] フラグメントのソートAO番号 \c ifao のAOの
 *     分子全体における非ソートAO番号
 *
 * @retval  0 正常終了
 * @retval -1 異常終了（一時配列確保に失敗したなど）
 *
 * @note
 * 各軌道量子数のCSに含まれるAO数 \c NNAO を内部で定義しているが、
 * どこか他の場所(例えば、\c ofmo-basis )で定義した方がよい
 *
 * @ingroup ofmo-input
 * */
int ofmo_assign_fragment_basis(
	// 入力引数
	const int nat, const int nsbs, char *basis_name_list[],
	const int atomic_number[], const int atom_basis[],
	const int fatom2tatom[],
	// 出力用引数
	int *fncs, int *fnao, int *fnps, int *fnpspair,
	int fleading_cs[], int fshel_tem[], int fshel_atm[],
	int fshel_add[], int fshel_ini[],
	double fprim_exp[], double fprim_coe[],
	int fsao2fuao[], int fuao2tuao[], int fsao2tuao[]
	) {
    int maxlqn, ncs, nao, nps, npspair;
    //
    static int called = false;
    static int *atm_lcs, *ushel_ini, NNAO[MAXLQN+1];
    if ( alloc() != 0 ) return -1;
    /* この関数が初めて呼ばれたときの初期処理 */
    if ( !called ) {
	int ierr;
	ierr = ofmo_data_get_vals("atm_lcs ushel_ini",
		&atm_lcs, &ushel_ini);
	if ( ierr != 0 ) return -1;
	for ( int lqn=0; lqn<=MAXLQN; lqn++ ) NNAO[lqn]=(lqn+1)*(lqn+2)/2;
	called = true;
    }
    /* 基底関数のサイズの取得 */
    ofmo_get_basis_size(
	    nat, nsbs, basis_name_list, atomic_number, atom_basis,
	    &maxlqn, &ncs, &nao, &nps );
    /* 非ソート基底関数パラメタの代入 */
    ofmo_assign_basis(
	    nat, nsbs, basis_name_list, atomic_number, atom_basis,
	    fushel_lqn, fushel_tem, fushel_atm, fushel_add, fushel_ini,
	    fuprim_exp, fuprim_coe );
    /* ソート機定款数パラメタの代入 */
    ofmo_sort_basis( maxlqn, ncs,
	    fushel_lqn, fushel_tem, fushel_atm, fushel_add, fushel_ini,
	    fuprim_exp, fuprim_coe,
	    fleading_cs, fshel_tem, fshel_atm, fshel_add, fshel_ini,
	    fprim_exp, fprim_coe, fsao2fuao );
    /* 非ソートモノマーAOから分子全体へのAOへ対応表の作成 */
    // まず、フラグメントの各原子の先頭CS番号を求める
    int atm_now, atm;
    atm_now     = 0;
    fatm_lcs[0] = 0;
    for ( int ics=1; ics<ncs; ics++ ) {
	atm = fushel_atm[ics];
	if ( atm != atm_now ) {
	    fatm_lcs[atm] = ics;
	    atm_now       = atm;
	}
    }
    fatm_lcs[nat] = ncs;
    // 非ソートフラグメントAO番号を非ソート全体AO番号へ変換する表の作成
    for ( int iat=0; iat<nat; iat++ ) {
	int ics, ics0, iao, iao0;
	int ifcs, ifcs0, ifcs1, ifao, ifao0, ifao1;
	atm  = fatom2tatom[iat];
	ics0 = atm_lcs[atm];
	ifcs0 = fatm_lcs[iat];
	ifcs1 = fatm_lcs[iat+1];
	for ( ifcs=ifcs0, ics=ics0; ifcs<ifcs1; ifcs++, ics++ ) {
	    iao0  =  ushel_ini[ics];
	    ifao0 = fushel_ini[ifcs];
	    ifao1 = ifao0 + NNAO[ fushel_lqn[ifcs] ];
	    for ( ifao=ifao0, iao=iao0; ifao<ifao1; ifao++, iao++ ) {
		fuao2tuao[ifao] = iao;
	    }
	}
    }
    // ソートフラグメントAO番号を非ソート全体AO番号へ変換する表の作成
    for ( int iao=0; iao<nao; iao++ )
	fsao2tuao[iao]=fuao2tuao[ fsao2fuao[iao] ];
    /* PSペア数のカウント */
    npspair=0;
    for ( int ics=0; ics<ncs; ics++ ) {
	for ( int jcs=0; jcs<=ics; jcs++ )
	    npspair += fushel_tem[ics]*fushel_tem[jcs];
    }

    *fncs = ncs;
    *fnao = nao;
    *fnps = nps;
    *fnpspair = npspair;
    return 0;
}
