/**
 * @file ofmo-fragment.c
 *
 * 入力データを元にして、各モノマーの情報を作成する関数を定義したファイル
 *
 * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ofmo-def.h"
#include "ofmo-data.h"

#include "ofmo-prof.h"

extern int ofmo_assign_fragment_basis(
	const int nat, const int nsbs, char *basis_name_list[],
	const int atomic_number[], const int atom_basis[],
	const int fatom2tatom[],
	int *fncs, int *fnao, int *fnps, int *fnpspair,
	int fleading_cs[], int fshel_tem[], int fshel_atm[],
	int fshel_add[], int fshel_ini[],
	double fprim_exp[], double fprim_coe[],
	int fsao2fuao[], int fuao2tuao[], int fsao2tuao[] );

/* ---------------------------------------------------------------
 * フラグメント（ダイマー、trimer, ... etc.)のデータを作成する関数たち
 * --------------------------------------------------------------- */
// フラグメントの原子データ、結合情報作成で使用する変数
static int *iamf          = NULL;
static int *frg_mk_tbl1   = NULL;
static int *oatomic_number = NULL; /* 元の原子番号 */
// 非ソート基底関数データ関連
static int *fushel_lqn    = NULL;
static int *fushel_tem    = NULL;
static int *fushel_atm    = NULL;
static int *fushel_add    = NULL;
static int *fushel_ini    = NULL;
static double *fuprim_exp = NULL;
static double *fuprim_coe = NULL;

/** フラグメントの情報を作成する際に用いる一時配列の開放を行う関数
 * */
static void dealloc() {
    Free( iamf );
    Free( frg_mk_tbl1 );
    Free( oatomic_number );

    Free( fushel_lqn );
    Free( fushel_tem );
    Free( fushel_atm );
    Free( fushel_add );
    Free( fushel_ini );
    Free( fuprim_exp );
    Free( fuprim_coe );

    iamf         = NULL;
    frg_mk_tbl1  = NULL;
    oatomic_number = NULL;

    fushel_lqn   = NULL;
    fushel_tem   = NULL;
    fushel_atm   = NULL;
    fushel_add   = NULL;
    fushel_ini   = NULL;
    fuprim_exp   = NULL;
    fuprim_coe   = NULL;
}

/** フラグメントの情報を作成する際に用いる一時配列の確保を行う関数
 * */
static int alloc() {    
    static int called = false;
    if ( called ) return 0;
    int ierr;
    int nbody, nao, maxlqn, maxnfbond;
    int maxnfatom, maxnfcs, maxnfps;
    int NAT, NCS, NPS, NBND;
    ierr = ofmo_data_get_vals(
	    "nbody maxnfcs nao maxlqn maxnfbond maxnfatom maxnfps",
	    &nbody, &maxnfcs, &nao, &maxlqn,
	    &maxnfbond, &maxnfatom, &maxnfps );
    if ( ierr != 0 ) return -1;
    NAT = nbody * maxnfatom;
    NCS = nbody * maxnfcs;
    NPS = nbody * maxnfps;
    NBND = nbody * maxnfbond;
    iamf        = (int*)malloc( sizeof(int) * NAT );
    frg_mk_tbl1 = (int*)malloc( sizeof(int) * NBND );
    oatomic_number = (int*)malloc( sizeof(int) * NAT );
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

/** フラグメント（ダイマー、trimer, ... etc.)の原子データを作成する関数
 *
 * @attention
 * - この関数では以下の変数を代入する
 *     - \c fatomic_number[iatm] フラグメントの \c iatm 番目の原子の核電荷
 *     - \c fatom2tatom [iatm] フラグメントの \c iatm 番目の原子の
 *         分子全体でのシリアル番号
 *     - \c fatom_x[iatm] フラグメントの \c iatm 番目の原子のx座標
 *         （au単位）
 *     - \c fatom_y[iatm] フラグメントの \c iatm 番目の原子のy座標
 *         （au単位）
 *     - \c fatom_z[iatm] フラグメントの \c iatm 番目の原子のz座標
 *         （au単位）
 *
 * @ingroup ofmo-input
 * */
static int ofmo_make_fragment_atomic_data(
	const int nmonomer, const int monomer_list[] ) {
    int imon, jmon, iatm, jatm, ibnd, jbnd, ifrag, jfrag;
    int i, j, ibda, jbda;
    int ix;
    int ierr;
    // フラグメントの原子データ
    static int *fatomic_number = NULL;
    static int *fatom2tatom    = NULL;
    static double *fatom_x     = NULL;
    static double *fatom_y     = NULL;
    static double *fatom_z     = NULL;
    int fnat;
    //
    static int called = false;
    static int *nfatom, *nfbond, **fbondsn1, *woelec, **fbda;
    static int **ifatom, **matomic_number, *atomic_number;
    static double *atom_x, *atom_y, *atom_z;
    static double **matom_x, **matom_y, **matom_z;
    if ( alloc() != 0 ) return -1;
    /* 初めて呼ばれたときの初期処理 */
    if ( !called ) {
	int nbody, maxnfatom;
	int NAT;
	ierr = ofmo_data_get_vals(
		"nfatom nfbond fbsn1 fbda woelec ifatom matn atn "
		"nbody maxnfatom atx aty atz "
		"matx maty matz",
		&nfatom, &nfbond, &fbondsn1, &fbda, &woelec,
		&ifatom, &matomic_number, &atomic_number,
		&nbody, &maxnfatom,
		&atom_x, &atom_y, &atom_z,
		&matom_x, &matom_y, &matom_z );
	if ( ierr != 0 ) return -1;
	NAT = nbody * maxnfatom;
	fatomic_number = (int*)malloc( sizeof(int) * NAT );
	fatom2tatom    = (int*)malloc( sizeof(int) * NAT );
	fatom_x        = (double*)malloc( sizeof(double) * NAT );
	fatom_y        = (double*)malloc( sizeof(double) * NAT );
	fatom_z        = (double*)malloc( sizeof(double) * NAT );
	ierr = ofmo_data_put_vals(
		"fatn fat2tat fatx faty fatz",
		fatomic_number, fatom2tatom, fatom_x, fatom_y, fatom_z );
	if ( ierr != 0 ) return -1;
	called = true;
    } else {
	ierr = ofmo_data_get_vals("fatn fat2tat fatx faty fatz",
		&fatomic_number, &fatom2tatom,
		&fatom_x, &fatom_y, &fatom_z );
	if ( ierr != 0 ) return -1;
    }
    // モノマーの場合の処理（データのコピー）
    if ( nmonomer == 1 ) {
	ifrag = monomer_list[0];
	fnat = nfatom[ifrag];
	memcpy( fatomic_number, matomic_number[ifrag], sizeof(int)*fnat );
	memcpy( fatom2tatom, ifatom[ifrag], sizeof(int)*fnat );
	memcpy( fatom_x, matom_x[ifrag], sizeof(double)*fnat );
	memcpy( fatom_y, matom_y[ifrag], sizeof(double)*fnat );
	memcpy( fatom_z, matom_z[ifrag], sizeof(double)*fnat );
	ofmo_data_put_vals("fnatom", fnat );
	return 0;
    }
    /* 複数のモノマーで構成されるフラグメント（ダイマーなど）のデータを
     * 作成する場合、あるモノマーに含まれる原子が、結合原子(bond detouched
     * atom, BDA)として、別のモノマーに含まれていることがある。
     * また、あるモノマーのBDAが、別のモノマーのBDAとして含まれている場合
     * もある。
     * ここでは、フラグメントを構成するあるモノマーの原子が、他のモノマー
     * に含まれているか、さらに、あるモノマーのBDAがフラグメントを構成する
     * 別のモノマーのBDAとして含まれているか、のチェックを行う。
     * */
    // 初期化
    for ( imon=0, ix=0; imon<nmonomer; imon++ ) {
	for ( iatm=0; iatm<nfatom[ monomer_list[imon] ]; iatm++, ix++ )
	    iamf[ix] = -1;
    }
    for ( imon=0, ix=0; imon<nmonomer; imon++ ) {
	for ( ibnd=0; ibnd<nfbond[ monomer_list[imon] ]; ibnd++, ix++ )
	    frg_mk_tbl1[ix] = 0;
    }
    // データ作成
    int katm, kbnd, latm, lbnd;
    katm = kbnd = 0;
    for ( imon=0; imon<nmonomer; imon++ ) {
	ifrag = monomer_list[imon];
	latm = lbnd = 0;
	for ( jmon=0; jmon<imon; jmon++ ) {
	    jfrag = monomer_list[jmon];
	    for ( ibnd=0; ibnd<nfbond[ifrag]; ibnd++ ) {
		i = fbondsn1[ifrag][ibnd]; /* 結合原子ペアの全体でのSN */
		i = ( i<0 ? -i : i );
		ibda = woelec[i-1];
		for ( jbnd=0; jbnd<nfbond[jfrag]; jbnd++ ) {
		    j = fbondsn1[jfrag][jbnd];
		    j = ( j<0 ? -j : j );
		    jbda = woelec[j-1];
		    /* 電子を与える結合原子が２つのモノマーで同じ場合 */
		    if ( ibda == jbda ) {
			iatm = fbda[ifrag][ibnd]; /* モノマー内でのSN */
			jatm = fbda[jfrag][jbnd];
			if ( iamf[latm+jatm] == -1 )
			    iamf[katm+iatm] = latm+jatm;
			else
			    iamf[katm+iatm] = iamf[latm+jatm];
			if ( i == j ) { /* 同じ結合原子ペアの場合 */
			    frg_mk_tbl1[kbnd+ibnd] = -1;
			    frg_mk_tbl1[lbnd+jbnd] = -1;
			}
			break;
		    }
		}	// for ( jbnd );
	    }		// for ( ibnd );
	    lbnd += nfbond[jfrag];
	    latm += nfatom[jfrag];
	}	// for ( jmon );
	kbnd += nfbond[ifrag];
	katm += nfatom[ifrag];
    }	// for ( imon );
    /* ここまでの結果を用いて、フラグメントの原子データを作成する
     * fatomic_number[iat] フラグメントのiat番目の原子の原子番号（核電荷）
     * oatomic_number[iat] フラグメントのiat番目の原子の元の原子番号
     * fatom2tatom[iat] フラグメントのiat番目の原子の全体に於けるSN
     *
     * */
    int atm, itmp;
    fnat = katm = 0;
    for ( imon=0; imon<nmonomer; imon++ ) {
	ifrag = monomer_list[imon];
	for ( iatm=0; iatm<nfatom[ifrag]; iatm++ ) {
	    if ( iamf[katm+iatm] == -1 ) {
		atm = ifatom[ifrag][iatm];
		iamf[katm+iatm] = fnat;/* フラグメント内での原子のSN */
		fatomic_number[fnat] = matomic_number[ifrag][iatm];
		oatomic_number[fnat] = atomic_number[ atm ];
		fatom2tatom[fnat]    = atm;
		fatom_x[fnat]        = atom_x[atm];
		fatom_y[fnat]        = atom_y[atm];
		fatom_z[fnat]        = atom_z[atm];
		fnat++;
	    } else {
		itmp = iamf[katm+iatm]; /*どの原子のデータと同一か? */
		iamf[katm+iatm] = iamf[itmp];
		fatomic_number[ iamf[katm+iatm] ] +=
		    matomic_number[ifrag][iatm];
	    }
	}
	katm += nfatom[ifrag];
    }
    // 原子数の登録
    ierr = ofmo_data_put_vals("fnatom", fnat );
    if ( ierr != 0 ) return -1;
    return 0;
}

/* ==============================================================
   ダイマー、トリマーなどの基底関数データに関連する関数
   ============================================================== */
static int *fatom_basis = NULL;
static int *ofatomic_number = NULL;
static void dealloc2() {
    Free( fatom_basis );
    Free( ofatomic_number );
    fatom_basis = NULL;
    ofatomic_number = NULL;
}

static int alloc2() {
    static int called = false;
    if ( called ) return 0;
    int ierr, nbody, maxnfatom, *atom_basis;
    ierr = ofmo_data_get_vals("nbody maxnfatom atombs",
	    &nbody, &maxnfatom, &atom_basis );
    if ( ierr != 0 ) return -1;
    if ( atom_basis != NULL ) {
	fatom_basis = (int*)malloc( sizeof(int) * nbody * maxnfatom );
	if ( fatom_basis == NULL ) return -1;
    } else fatom_basis = NULL;

    ofatomic_number = (int*)malloc( sizeof(int) * nbody * maxnfatom );
    if ( ofatomic_number == NULL ) return -1;
    atexit( dealloc2 );
    called = true;
    return 0;
}

/* この関数内部で領域確保する配列 */
// ソート基底関数データ関連
static int *fleading_cs   = NULL;
static int *fshel_tem     = NULL;
static int *fshel_atm     = NULL;
static int *fshel_add     = NULL;
static int *fshel_ini     = NULL;
static double *fprim_exp  = NULL;
static double *fprim_coe  = NULL;
// AO変換テーブル関連
static int *fsao2fuao     = NULL;
static int *fsao2tuao     = NULL;
static int *fuao2tuao     = NULL;
/** モノマー以外のフラグメントの基底関数データのテーブルを作成する
 *
 * @note この関数は、ダイマーなどを構成する原子の情報（原子数、原子番号
 * など）を元にして、基底関数データを作成する。
 * したがって、フラグメントの情報が定義されていない時点で呼び出すと、
 * 結果は不定（あるいは、異常終了）となる。
 *
 * @ingroup ofmo-input
 *
 * */
static int ofmo_make_fragment_basis_data(
	const int nmonomer, const int monomer_list[] ) {
    int ierr, nat, ncs, nao, nps, npspair;
    /* 外部で作成された配列 */
    static int *fatomic_number;
    static int *fatom2tatom;
    static double *fatom_x;
    static double *fatom_y;
    static double *fatom_z;
    static int nsbs, maxlqn;
    static char **basis_name_list;
    static int *atomic_number;
    static int *atom_basis;
    /* モノマーの基底関数データ */
    static int **mshel_tem, **mshel_atm, **mshel_add, **mshel_ini;
    static int **mlcs, *nfcs, *nfao, *nfps;
    static double **mprim_exp, **mprim_coe;
    static int **msao2muao, **msao2tuao, **muao2tuao;
    //
    static int called = false;
    if ( alloc2() != 0 ) return -1;
    if ( !called ) {
	int maxnfatom, maxnfcs, maxnfao, maxnfps, nbody;
	int NCS, NAO, NPS;
	ierr = ofmo_data_get_vals(
		"maxlqn maxnfatom maxnfcs maxnfao maxnfps nbody "
		"atombs nsbs bslst atn",
		&maxlqn, &maxnfatom, &maxnfcs, &maxnfao, &maxnfps, &nbody,
		&atom_basis, &nsbs, &basis_name_list,
		&atomic_number );
	ierr += ofmo_data_get_vals("fatn fat2tat fatx faty fatz",
		&fatomic_number, &fatom2tatom,
		&fatom_x, &fatom_y, &fatom_z );
	ierr += ofmo_data_get_vals("mshel_tem mshel_atm mshel_add mshel_ini "
		"mlcs mprim_exp mprim_coe nfcs nfao nfps "
		"msao2muao msao2tuao muao2tuao",
		&mshel_tem, &mshel_atm, &mshel_add, &mshel_ini,
		&mlcs, &mprim_exp, &mprim_coe, &nfcs, &nfao, &nfps,
		&msao2muao, &msao2tuao, &muao2tuao );

	if ( ierr != 0 ) return -1;
	NCS = nbody * maxnfcs;
	NAO = nbody * maxnfao;
	NPS = nbody * maxnfps;
	fleading_cs    = (int*)malloc( sizeof(int) * (maxlqn+1+1) );
	fshel_tem      = (int*)malloc( sizeof(int) * NCS );
	fshel_atm      = (int*)malloc( sizeof(int) * NCS );
	fshel_add      = (int*)malloc( sizeof(int) * NCS );
	fshel_ini      = (int*)malloc( sizeof(int) * NCS );
	fprim_exp      = (double*)malloc(sizeof(double) * NPS );
	fprim_coe      = (double*)malloc( sizeof(double) * NPS );
	fsao2fuao      = (int*)malloc( sizeof(int) * NAO );
	fsao2tuao      = (int*)malloc( sizeof(int) * NAO );
	fuao2tuao      = (int*)malloc( sizeof(int) * NAO );
	ierr = ofmo_data_put_vals(
		"flcs fshel_tem fshel_atm fshel_add fshel_ini fprim_exp "
		"fprim_coe fsao2fuao fsao2tuao fuao2tuao",
		fleading_cs, fshel_tem, fshel_atm, fshel_add, fshel_ini,
		fprim_exp, fprim_coe, fsao2fuao, fsao2tuao, fuao2tuao );
	if ( ierr != 0 ) return -1;
	called = true;
    }

    // モノマーの場合の処理（データのコピー）
    if ( nmonomer == 1 ) {
	int ifrag;
	ifrag = monomer_list[0];
	ncs = nfcs[ifrag];
	nao = nfao[ifrag];
	nps = nfps[ifrag];
	memcpy( fleading_cs, mlcs[ifrag], sizeof(int)*(maxlqn+1+1) );
	memcpy( fshel_tem, mshel_tem[ifrag], sizeof(int) * ncs );
	memcpy( fshel_atm, mshel_atm[ifrag], sizeof(int) * ncs );
	memcpy( fshel_add, mshel_add[ifrag], sizeof(int) * ncs );
	memcpy( fshel_ini, mshel_ini[ifrag], sizeof(int) * ncs );
	memcpy( fprim_exp, mprim_exp[ifrag], sizeof(double) * nps );
	memcpy( fprim_coe, mprim_coe[ifrag], sizeof(double) * nps );
	memcpy( fsao2fuao, msao2muao[ifrag], sizeof(int) * nao );
	memcpy( fsao2tuao, msao2tuao[ifrag], sizeof(int) * nao );
	memcpy( fuao2tuao, muao2tuao[ifrag], sizeof(int) * nao );
	ofmo_data_put_vals("fncs fnao fnps", ncs, nao, nps );
	return 0;
    }

    ierr = ofmo_data_get_vals("fnatom", &nat );
    for ( int iat=0; iat<nat; iat++ )
	ofatomic_number[iat] = atomic_number[ fatom2tatom[iat] ];
    if ( atom_basis != NULL ) {
	for ( int iat=0; iat<nat; iat++ )
	    fatom_basis[iat] = atom_basis[ fatom2tatom[iat] ];
    }
    ierr = ofmo_assign_fragment_basis(
	    nat, nsbs, basis_name_list, ofatomic_number, fatom_basis,
	    fatom2tatom,
	    &ncs, &nao, &nps, &npspair,
	    fleading_cs, fshel_tem, fshel_atm, fshel_add, fshel_ini,
	    fprim_exp, fprim_coe, fsao2fuao, fuao2tuao, fsao2tuao );
    if ( ierr != 0 ) return -1;
    ierr = ofmo_data_put_vals("fncs fnao fnps", ncs, nao, nps );
    if ( ierr != 0 ) return -1;
    return 0;
}

int ofmo_fragment_init(
	const int nmonomer, const int monomer_list[] ) {
    int ierr;
    ierr  = ofmo_make_fragment_atomic_data( nmonomer, monomer_list );
    ierr += ofmo_make_fragment_basis_data( nmonomer, monomer_list );
    return ierr;
}
