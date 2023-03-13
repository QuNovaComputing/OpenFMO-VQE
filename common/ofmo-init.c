#include <stdio.h>
#include <stdlib.h>

#include "ofmo-parallel.h"

#include "ofmo-def.h"
#include "ofmo-data.h"
#include "ofmo-basis.h"
#include "ofmo-string.h"
#include "ofmo-misc.h"

extern int ofmo_read_gamess_input( const char *filename );
extern int ofmo_assign_fragment_basis(
	// 入力引数
	const int nat, const int nsbs, char *basis_name_list[],
	const int atomic_number[], const int atom_basis[],
	const int fatom2tatom[],
	// 出力引数
	int *fncs, int *fnao, int *fnps, int *fnpspair,
	int fleading_cs[], int fshel_tem[], int fshel_atm[],
	int fshel_add[], int fshel_ini[],
	double fprim_exp[], double fprim_coe[],
	int fsao2fuao[], int fuao2tuao[], int fsao2tuao[]);

/** GAMESSの入力で得られた基底関数情報から、基底関数名を得る
 * */
static int ofmo_get_basis_name( char *basis_name ) {
    int ngauss, ndfunc, nffunc, npfunc;
    char gbasis[MAXSTRLEN], *gb, *lmobs;
    // local variable
    int flag, ierr;
    ierr = ofmo_data_get_vals(
	    "ngauss ndfunc nffunc npfunc gbasis lmobs",
	    &ngauss, &ndfunc, &nffunc, &npfunc, &gb, &lmobs );
    if ( ierr != 0 ) return -1;
    strcpy( gbasis, gb );
    ofmo_toupper( gbasis );
    flag = false;
    if ( strcmp(gbasis, "STO") == 0 ) {
	if ( ngauss == 3 ) {
	    strcpy( basis_name, "STO-3G" );
	    flag = true;
	}
    } else if ( strcmp(gbasis, "N31") == 0 ) {
	if ( ngauss == 6 ) {
	    if ( ndfunc == 1 ) {
		if ( npfunc == 1 ) strcpy(basis_name, "6-31G**");
		else               strcpy(basis_name, "6-31G*");
		flag = true;
	    } else if ( ndfunc < 1 ) {
		strcpy(basis_name, "6-31G");
		flag = true;
	    } else if ( ndfunc == 3 ) {
		strcpy( basis_name, "6-31G(3df,3pd)");
		flag = true;
	    }
	}
    }
    // もし$basis の部分の内容で基底関数が決められなかった場合の処理
    if ( flag == false ) {
	if ( lmobs[0] == '\0' ) {
	    dbg("Can\'t determine basis set name ("
		    "ngauss=%d, ndfunc=%d, nffunc=%d, npfunc=%d,"
		    " gbasis=%s\n",
		    ngauss, ndfunc, nffunc, npfunc,
		    ( (gbasis==NULL)? " " : gbasis) );
	    return -1;
	}
	strcpy( basis_name, lmobs );
    } else {
	if ( strcmp( basis_name, lmobs ) != 0 ) {
	    dbg("WARNING:basis set names are different between %s vs %s\n",
		    basis_name, lmobs );
	}
    }
    return 0;
}

/** モノマーの情報を作成する関数
 *
 * 入力データに基づいて、モノマー情報を作成する関数
 *
 * - この関数で作成されるモノマー情報
 *     - \c nfatom[ifrag] モノマー \c ifrag の原子数（ghost atomを含む）
 *     - \c nfbond[ifrag] モノマー \c ifrag に含まれる結合原子ペアの数
 *     - \c maxnfatom モノマーに含まれる原子数の最大値
 *         (\f$ \max(\tt{nfatom})\f$)
 *     - \c maxnfbond モノマーに含まれる結合原子ペア数の最大値
 *         (\f$ \max(\tt{nfbond})\f$)
 *     - \c ifatom[ifrag][iatm] モノマー \c ifrag の \c iatm 番目の
 *         原子の分子全体における原子のシリアル番号
 *     - \c matomic_number[ifrag][iatm] モノマー \c ifrag の
 *         \c iatm 番目の原子の核電荷。結合原子ペア部分は、元の原子番号
 *         と異なる。
 *     - \c fbondsn1[ifrag][idb] モノマー \c ifrag の \c idb 番目の
 *         結合原子ペアの分子全体における結合原子ペア番号+1。
 *         正の値の場合は、該当モノマーに電子をもらう結合原子が含まれて
 *         いることを表す。逆に、負の値の場合は、該当モノマーに電子を
 *         与える結合原子が含まれていることを表す。
 *     - \c fbda[ifrag][ibd] 電子を与える結合原子のモノマー \c ifrag に
 *         おけるシリアル番号。電子を与える結合原子は、少なくとも
 *         ghost atomとして、必ず、モノマーに含まれる。
 *     - \c matom_x[ifrag][iatm] モノマー \c ifrag の \c iatm 番目の
 *         原子のx座標（au単位）
 *     - \c matom_y[ifrag][iatm] モノマー \c ifrag の \c iatm 番目の
 *         原子のy座標（au単位）
 *     - \c matom_z[ifrag][iatm] モノマー \c ifrag の \c iatm 番目の
 *         原子のz座標（au単位）
 *
 *
 * - この関数で参照されるデータ（要するに、この関数呼び出し時に
 *   以下のデータが代入されていないと、この関数の結果はおかしくなる）。
 *   入力データとして与えられるデータである。
 *     - \c nfrag モノマー数
 *     - \c welec[ibd] 分子全体における結合原子ペア番号 \c idb において、
 *         電子をもらう結合原子の分子全体における原子のシリアル番号
 *     - \c woelec[idb] 分子全体における結合原子ペア番号 \c idb において、
 *         電子を与える結合原子の分子全体における原子のシリアル番号
 *     - \c nbond 分子全体での結合原子ペア数
 *     - \c at2frg[iatm] 分子全体における \c iatm 番目の原子が属する
 *         モノマー番号
 *     - \c icharg[ifrag] モノマー \c ifrag の電荷
 *     - \c atn[iatm] 分子全体における \c iatm 番目の原子の原子番号
 *     - \c atx[iatm] 分子全体における \c iatm 番目の原子のx座標
 *     - \c aty[iatm] 分子全体における \c iatm 番目の原子のy座標
 *     - \c atz[iatm] 分子全体における \c iatm 番目の原子のz座標
 *     - \c natom 分子全体に含まれる原子数
 *
 * @ingroup ofmo-input
 * */
static int ofmo_make_monomer_atomic_data() {
    /* この関数で使用する、外部で定義された変数 */
    int nfrag, *welec, *woelec, nbond, *at2frg, *icharg;
    int *atomic_number, natom;
    double *atx, *aty, *atz;
    /* この関数で定義する変数 */
    int *nfatom, *nfbond, **ifatom, **matomic_number;
    double **matom_x, **matom_y, **matom_z;
    int maxnfatom, maxnfbond;
    /* ローカル変数 */
    int i, ifrag, jfrag, iat, ibsn1, ibda;

    ofmo_data_get_vals(
	    "nfrag welec woelec nbond at2frg icharg atn atx aty atz natom",
	    &nfrag, &welec, &woelec, &nbond, &at2frg, &icharg,
	    &atomic_number, &atx, &aty, &atz, &natom );
    nfatom = (int*)malloc(sizeof(int) * nfrag);
    nfbond = (int*)malloc(sizeof(int) * nfrag);
    /* ************* モノマーの原子数のカウント **************
     * ここで決まる変数
     * nfatom[], maxnfatom
     * nfbond[], maxnfbond
     * */
    for ( i=0; i<nfrag; i++ ) nfatom[i] = nfbond[i] = 0;
    for ( i=0; i<natom; i++ ) nfatom[ at2frg[i] ]++;
    /* count the # of BDAs and ghost H atoms in each monomer */
    for ( i=0; i<nbond; i++ ) {
	ifrag = at2frg[ welec[i] ];
	jfrag = at2frg[ woelec[i] ];
	nfatom[ifrag]++;	/* add ghost H atom */
	nfbond[ifrag]++;	/* add BDA for monomer */
	nfbond[jfrag]++;	/* add BDA for monomer */
    }
    maxnfatom = ofmo_imax( nfrag, nfatom );
    maxnfbond = ofmo_imax( nfrag, nfbond );
    /* *********** 各モノマーのi番目の原子の全体におけるSNなど
     * ここで決まる変数
     * ifatom[ifrag][iatm]
     * matomic_number[ifrag][iatm]
     * */
    /* subtract atomic number from BDA without electron (temporary) */
    for (i=0; i<nbond; i++) atomic_number[ woelec[i] ]--;
    /* make ifatom and matomic_number */
    ifatom         = ofmo_alloc_imatrixv( nfrag, nfatom );
    matomic_number = ofmo_alloc_imatrixv( nfrag, nfatom );
    for ( i=0; i<nfrag; i++ ) nfatom[i] = 0;	/* clear nfatom */
    for ( i=0; i<natom; i++ ) {
	ifrag = at2frg[i];
	ifatom[ ifrag ][ nfatom[ifrag] ] = i;
	matomic_number[ ifrag ][ nfatom[ifrag] ] = atomic_number[i];
	nfatom[ifrag]++;
    }
    for ( i=0; i<nbond; i++ ) {	/* add ghost H atom */
	ifrag = at2frg[ welec[i] ];
	ifatom[ ifrag ][ nfatom[ifrag] ] = woelec[i];
	matomic_number[ ifrag ][ nfatom[ifrag] ] = 1;
	nfatom[ifrag]++;
    }
    /* make fbondsn1 and fbda                      */
    /* fbondsn1[ifrag][idb] = BDAのシリアル番号+1  */
    /*     正数＝電子をもらうモノマーの場合        */
    /*     負数＝電子を与えるモノマーの場合        */
    /* fbda[ifrag][ibd] = 電子を与えるBDAのモノマー内でのシリアル番号 */
    int **fbondsn1, **fbda;
    fbondsn1 = ofmo_alloc_imatrixv( nfrag, nfbond );
    fbda     = ofmo_alloc_imatrixv( nfrag, nfbond );
    for ( i=0; i<nfrag; i++ ) nfbond[i] = 0;	/* clear nfbond */
    for ( i=0; i<nbond; i++ ) {
	ifrag = at2frg[ welec[i] ];
	jfrag = at2frg[ woelec[i] ];
	fbondsn1[ifrag][ nfbond[ifrag] ] =  (i+1);	/* 1-offset */
	fbondsn1[jfrag][ nfbond[jfrag] ] = -(i+1);	/* 1-offset */
	nfbond[ifrag]++;
	nfbond[jfrag]++;
    }
    for ( ifrag=0; ifrag<nfrag; ifrag++ ) {
	for ( i=0; i<nfbond[ifrag]; i++ ) {
	    ibsn1 = abs(fbondsn1[ifrag][i]);
	    ibda  = woelec[ ibsn1 - 1 ];
	    for ( iat=0; iat<nfatom[ifrag]; iat++ )
		if ( ifatom[ifrag][iat] == ibda ) break;
	    if ( iat >= nfatom[ifrag] ) {
		dbg("ERROR: BDA %d in bond %d in monomer %5d\n",
			ibda, ibsn1, ifrag );
		return -1;
	    }
	    fbda[ifrag][i] = iat;
	}
    }

    /* **** モノマーの原子の中心座標の定義
     * ここで定義される変数
     * matom_x[ifrag][iatm]
     * matom_y[ifrag][iatm]
     * matom_z[ifrag][iatm]
     * */
    matom_x = ofmo_alloc_dmatrixv( nfrag, nfatom );
    matom_y = ofmo_alloc_dmatrixv( nfrag, nfatom );
    matom_z = ofmo_alloc_dmatrixv( nfrag, nfatom );
    for ( ifrag=0; ifrag<nfrag; ifrag++ ) {
	for ( i=0; i<nfatom[ifrag]; i++ ) {
	    iat = ifatom[ifrag][i];
	    matom_x[ifrag][i] = atx[iat];
	    matom_y[ifrag][i] = aty[iat];
	    matom_z[ifrag][i] = atz[iat];
	}
    }
    /* restore atomic number  */
    for (i=0; i<nbond; i++) atomic_number[ woelec[i] ]++;
    // 設定内容の保存
    ofmo_data_put_vals("nfatom nfbond maxnfatom maxnfbond ifatom "
	    "matn fbsn1 fbda matx maty matz",
	    nfatom, nfbond, maxnfatom, maxnfbond, ifatom,
	    matomic_number, fbondsn1, fbda, matom_x, matom_y, matom_z );

    return 0;
}

/** モノマーのソート基底関数データ、および、モノマー内のAOと
 * 全体のAOの変換テーブルを作成する
 *
 * すべてのモノマーのソート基底関数データ、および、モノマー内のAOと
 * 全体のAOの変換テーブルを作成する。
 * この関数では、以下の変数を設定する
 *
 * - この関数で設定される変数
 *     - \c nfcs[ifrag] モノマー \c ifrag のCS数
 *     - \c nfao[ifrag] モノマー \c ifrag のAO数
 *     - \c nfps[ifrag] モノマー \c ifrag のPS数
 *     - \c maxnfcs モノマーに含まれるCS数の最大値
 *     - \c maxnfao モノマーに含まれるAO数の最大値
 *     - \c maxnfps モノマーに含まれるPS数の最大値
 *     - \c mleading_cs[ifrag][lqn] モノマー \c ifrag の軌道量子数 \c lqn
 *         の先頭CS番号
 *     - \c mshel_tem[ifrag][ics] モノマー \c ifrag の \c ics 番目のCSの
 *         縮約長
 *     - \c mshel_atm[ifrag][ics] モノマー \c ifrag の \c ics 番目のCSが
 *         属する原子のモノマー内でのシリアル番号
 *     - \c mshel_add[ifrag][ics] モノマー \c ifrag の \c ics 番目のCSに
 *         属するPSの先頭PS番号
 *     - \c mshel_ini[ifrag][ics] モノマー \c ifrag の \c ics 番目のCSに
 *         属するAOの先頭AO番号
 *     - \c mprim_exp[ifrga][ips] モノマー \c ifrag の \c ips 番目のPSの
 *         軌道指数
 *     - \c mprim_coe[ifrga][ips] モノマー \c ifrag の \c ips 番目のPSの
 *         規格化定数込みの縮約係数
 *
 *     - \c msao2muao[ifrag][iao] モノマー \c ifrag のソートAO番号 \c iao
 *         のAOのモノマーにおける非ソートAO番号
 *     - \c msao2tuao[ifrag][iao] モノマー \c ifrag のソートAO番号 \c iao
 *         のAOの分子全体における非ソートAO番号
 *     - \c muao2tuao[ifrag][iao] モノマー \c ifrag の非ソートAO番号 \c iao
 *         のAOの分子全体における非ソートAO番号
 *     - \c maxnpspair モノマーのPSペア数の最大値
 *     .
 * - この関数内部で用いている変数
 *     - \c nfrag モノマー数
 *     - \c maxnfatom モノマーに含まれる原子数(ghost atomを含む)の最大値
 *     - \c ifatom[ifrag][iatm] モノマー \c ifrag の\c iatm 番目の原子の
 *         分子全体におけるシリアル番号
 *     - \c atn[iatm] 分子全体における \c iatm 番目の原子の原子番号
 *     - \c nfatom[ifrag] モノマー \c ifrag に含まれる原子数
 *         (ghost atomを含む)
 *     - \c maxlqn 最大軌道量子数
 *     - \c nsbs 基底関数の種類の数
 *     - \c bslst[ibs] 基底関数番号 \c ibs の基底関数の名前
 *     - \c atombs[iatm] 分子全体のおける \c iatm 番目の原子の基底関数種
 *         の番号。基底関数が１種類の場合には、NULLでも構わない。というか、
 *         NULLだと、基底関数が１種類だとみなす。
 *
 * @attention
 * - モノマーの非ソート基底関数データは作成しない
 *
 * @ingroup ofmo-input
 * */
static int ofmo_make_monomer_basis_data( int myrank ) {
    // 既に代入されていて、この関数内で用いるデータ
    int nfrag, maxnfatom, nsbs, maxlqn;
    int **ifatom, *atomic_number, *nfatom, *atom_basis;
    char **basis_name_list;
    // この関数内部で用いる一時配列、一時変数
    int *matomic_number, *matom_basis;
    int ao, cs, ps, mlqn, ifrag, ierr;
    int tncs, tnao, tnps;
    // この関数で求めたい値を表す変数
    int maxnfcs, maxnfao, maxnfps, maxnpspair;
    int *nfao, *nfcs, *nfps;
    // ソート基底関数
    int **mleading_cs, **mshel_tem, **mshel_atm, **mshel_add, **mshel_ini;
    double **mprim_exp, **mprim_coe;
    // 変換テーブル
    int **msao2tuao, **muao2tuao, **msao2muao;

    // debug
    //flag_dbg = ( myrank == 0 ? 1 : 0 );

    ierr = ofmo_data_get_vals("nfrag maxnfatom", &nfrag, &maxnfatom );
    ierr += ofmo_data_get_vals("ifatom atn nfatom",
	    &ifatom, &atomic_number, &nfatom );
    ierr += ofmo_data_get_vals("maxlqn nsbs bslst atombs",
	    &maxlqn, &nsbs, &basis_name_list, &atom_basis );
    if ( ierr != 0 ) return -1;
    /* 初期設定 */
    // 一時配列（モノマーに含まれる原子の元の原子番号を格納する配列）確保
    matomic_number = (int*)malloc( sizeof(int) * maxnfatom );
    if ( atom_basis != NULL )
	matom_basis = (int*)malloc( sizeof(int) * maxnfatom );
    else
	matom_basis = NULL;

    /* 各モノマーのCS数、AO数、PS数、および、その和と、最大値を求める */
    nfcs = (int*)malloc( sizeof(int) * nfrag );
    nfao = (int*)malloc( sizeof(int) * nfrag );
    nfps = (int*)malloc( sizeof(int) * nfrag );
    maxnfcs = maxnfao = maxnfps = 0;
    tncs = tnao = tnps = 0;
    for ( ifrag=0; ifrag<nfrag; ifrag++ ) {
	for ( int iat=0; iat<nfatom[ifrag]; iat++ )
	    matomic_number[iat] = atomic_number[ ifatom[ifrag][iat] ];
	if ( atom_basis != NULL ) {
	    for ( int iat=0; iat<nfatom[ifrag]; iat++ )
		matom_basis[iat] = atom_basis[ ifatom[ifrag][iat] ];
	}
	ofmo_get_basis_size(
		nfatom[ifrag], nsbs, basis_name_list,
		matomic_number, matom_basis,
		&mlqn, &cs, &ao, &ps );
	nfcs[ifrag] = cs;
	nfao[ifrag] = ao;
	nfps[ifrag] = ps;
	tncs += cs;
	tnao += ao;
	tnps += ps;
	if ( cs > maxnfcs ) maxnfcs = cs;
	if ( ao > maxnfao ) maxnfao = ao;
	if ( ps > maxnfps ) maxnfps = ps;
    }
    ofmo_data_put_vals("maxnfcs maxnfao maxnfps nfcs nfao nfps",
	    maxnfcs, maxnfao, maxnfps, nfcs, nfao, nfps );
    /* モノマーのソート基底関数の代入、および、モノマーのAOから
     * 全体のAOへの変換テーブルの作成                              */
    // メモリ確保
    mleading_cs = (int**)malloc( sizeof(int*) * nfrag );
    mshel_tem = (int**)malloc( sizeof(int*) * nfrag );
    mshel_atm = (int**)malloc( sizeof(int*) * nfrag );
    mshel_add = (int**)malloc( sizeof(int*) * nfrag );
    mshel_ini = (int**)malloc( sizeof(int*) * nfrag );
    mprim_exp = (double**)malloc( sizeof(double*) * nfrag );
    mprim_coe = (double**)malloc( sizeof(double*) * nfrag );
    msao2muao = (int**)malloc( sizeof(int*) * nfrag );
    msao2tuao = (int**)malloc( sizeof(int*) * nfrag );
    muao2tuao = (int**)malloc( sizeof(int*) * nfrag );
    mleading_cs[0] = (int*)malloc( sizeof(int) * nfrag * (maxlqn+1+1) );
    mshel_tem[0] = (int*)malloc( sizeof(int) * tncs );
    mshel_atm[0] = (int*)malloc( sizeof(int) * tncs );
    mshel_add[0] = (int*)malloc( sizeof(int) * tncs );
    mshel_ini[0] = (int*)malloc( sizeof(int) * tncs );
    mprim_exp[0] = (double*)malloc( sizeof(double) * tnps );
    mprim_coe[0] = (double*)malloc( sizeof(double) * tnps );
    msao2muao[0] = (int*)malloc( sizeof(int) * tnao );
    msao2tuao[0] = (int*)malloc( sizeof(int) * tnao );
    muao2tuao[0] = (int*)malloc( sizeof(int) * tnao );
    for ( ifrag=1; ifrag<nfrag; ifrag++ ) {
	mleading_cs[ifrag] = mleading_cs[ifrag-1] + (maxlqn+1+1);
	mshel_tem[ifrag] = mshel_tem[ifrag-1] + nfcs[ifrag-1];
	mshel_atm[ifrag] = mshel_atm[ifrag-1] + nfcs[ifrag-1];
	mshel_add[ifrag] = mshel_add[ifrag-1] + nfcs[ifrag-1];
	mshel_ini[ifrag] = mshel_ini[ifrag-1] + nfcs[ifrag-1];
	mprim_exp[ifrag] = mprim_exp[ifrag-1] + nfps[ifrag-1];
	mprim_coe[ifrag] = mprim_coe[ifrag-1] + nfps[ifrag-1];
	msao2muao[ifrag] = msao2muao[ifrag-1] + nfao[ifrag-1];
	msao2tuao[ifrag] = msao2tuao[ifrag-1] + nfao[ifrag-1];
	muao2tuao[ifrag] = muao2tuao[ifrag-1] + nfao[ifrag-1];
    }
    /* モノマーに対するソート基底関数の代入 */
    maxnpspair = 0;
    for ( ifrag=0; ifrag<nfrag; ifrag++ ) {
	int fncs, fnao, fnps, npspair;
	// モノマーを構成する原子の元の原子番号を代入
	for ( int iat=0; iat<nfatom[ifrag]; iat++ )
	    matomic_number[iat] = atomic_number[ ifatom[ifrag][iat] ];

	ierr = ofmo_assign_fragment_basis(
		nfatom[ifrag], nsbs, basis_name_list, matomic_number, NULL,
		ifatom[ifrag],
		&fncs, &fnao, &fnps, &npspair,
		mleading_cs[ifrag], mshel_tem[ifrag], mshel_atm[ifrag],
		mshel_add[ifrag], mshel_ini[ifrag],
		mprim_exp[ifrag], mprim_coe[ifrag],
		msao2muao[ifrag], muao2tuao[ifrag], msao2tuao[ifrag] );

	if ( npspair > maxnpspair ) maxnpspair = npspair;
    }	// ifrag
    ofmo_data_put_vals("mlcs mshel_tem mshel_atm mshel_add mshel_ini",
	    mleading_cs, mshel_tem, mshel_atm, mshel_add, mshel_ini );
    ofmo_data_put_vals("mprim_exp mprim_coe msao2muao msao2tuao muao2tuao",
	    mprim_exp, mprim_coe, msao2muao, msao2tuao, muao2tuao );
    ofmo_data_put_vals("maxnpspair", maxnpspair );

    free( matomic_number );
    return 0;
}

/** 分子全体の基底関数情報を代入する関数
 *
 * - この関数では、以下のデータが代入される
 *     - \c ncs 分子全体のCS数
 *     - \c nao 分子全体のAO数
 *     - \c nps 分子全体のPS数
 *     - \c maxlqn 最大軌道量子数
 *     - \c ushel_ini[ics] 全体でのCS番号 \c ics 番目のCSに含まれる
 *         AOの先頭AO番号
 *     - \c atm_lcs[iatm] 全体での \c iatm 番目の原子の先頭CS番号
 *         （非ソート基底関数における）
 *
 * - この関数内で参照する変数
 *     - \c nsbs 基底関数の種類の数
 *     - \c bslst[ibs] 基底関数種の番号\c ibs の基底関数名
 *     - \c natom 分子全体の原子数
 *     - \c atn[iatm] 全体での \c iatm 番目の原子の原子番号
 *
 * - この関数で利用している外部関数
 *     - \c ofmo_get_basis_size
 *     - \c ofmo_alloc_unsorted_basis
 *     - \c ofmo_assign_basis
 *     - \c ofmo_data_get_vals
 *     - \c ofmo_data_put_vals
 *
 * @ingroup ofmo-input
 * */
static int ofmo_assign_entire_molecule_basis() {
    int nsbs, natom, *atomic_number;
    char **basis_name_list;
    int ierr;
    int maxlqn, ncs, nao, nps;
    int *ushel_lqn, *ushel_tem, *ushel_atm, *ushel_add, *ushel_ini;
    double *uprim_exp, *uprim_coe;
    ierr = ofmo_data_get_vals("nsbs bslst natom atn",
	    &nsbs, &basis_name_list, &natom, &atomic_number );
    if ( ierr != 0 ) return -1;
    ierr += ofmo_get_basis_size(
	    natom, nsbs, basis_name_list, atomic_number, NULL,
	    &maxlqn, &ncs, &nao, &nps );

    ofmo_alloc_unsorted_basis( ncs, nps,
	    &ushel_lqn, &ushel_tem, &ushel_atm, &ushel_add, &ushel_ini,
	    &uprim_exp, &uprim_coe );
    ofmo_assign_basis(
	    natom, nsbs, basis_name_list, atomic_number, NULL,
	    ushel_lqn, ushel_tem, ushel_atm, ushel_add, ushel_ini,
	    uprim_exp, uprim_coe );
    /* 各原子の先頭CS番号 atm_lcs[] を作成 */
    int *atm_lcs;
    int atm_now, atm;
    atm_lcs = (int*)malloc( sizeof(int) * (natom+1) );
    atm_lcs[0] = 0;
    atm_now    = 0;
    for ( int ics=1; ics<ncs; ics++ ) {
	atm = ushel_atm[ics];
	if ( atm != atm_now ) {
	    atm_lcs[atm]     = ics;
	    atm_now          = atm;
	}
    }
    atm_lcs[natom] = ncs;
    // 登録
    ofmo_data_put_vals("ncs nao nps maxlqn ushel_lqn ushel_ini atm_lcs",
	    ncs, nao, nps, maxlqn, ushel_lqn, ushel_ini, atm_lcs );
    free( ushel_tem );
    free( ushel_atm );
    free( ushel_add );
    free( uprim_exp );
    free( uprim_coe );
    return 0;
}

/** FMO calculation initialization routine
 *
 * FMO計算の初期化ルーチン。この関数呼び出しを行うことで、
 * Initialization routine for FMO calculation. By calling this function,
 * all the information of each fragment is assigned.
 * The current input uses only one type of basis function.
 *
 * @param[in] filename Input file name
 * @param[in] comm A communicator that shares input data.
 * Normally, you can specify \c MPI_COMM_WORLD.
 *
 * @ingroup ofmo-input
 * */
int ofmo_init( const char *filename, MPI_Comm comm ) {
    char **basis_name_list;
    int nsbs=1;
    /* read input file */
    if ( ofmo_read_gamess_input( filename ) != 0 ) {
	dbg("ERROR: in ofmo_read_input_data\n");
	return -1;
    };
    /*  determine basis set name */
    basis_name_list = (char**)malloc( sizeof(char*) );
    basis_name_list[0] = (char*)malloc( sizeof(char) * MAXSTRLEN );
    basis_name_list[0][0] = '\0';
    if ( ofmo_get_basis_name( basis_name_list[0] ) != 0 ) {
	dbg("ERROR: in ofmo_get_basis_name\n");
	return -1;
    };
    ofmo_data_put_vals("nsbs bslst", nsbs, basis_name_list );
    /* set basis set for entire molecule */
    ofmo_assign_entire_molecule_basis();
    /* set monomer data */
    ofmo_make_monomer_atomic_data();
    ofmo_make_monomer_basis_data( 0 );
    return 0;
}
