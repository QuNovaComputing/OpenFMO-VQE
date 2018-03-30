/**
 * @file ofmo-scf.c
 *
 * SCF計算プログラムのサンプルコード。
 * 通信回数を減らすために、MPI_Allreduce + MPI_Bcastでの
 * 実装を行った。
 *
 * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ofmo-def.h"
#include "ofmo-mat.h"
#include "ofmo-integ.h"
#include "ofmo-diis.h"
#include "ofmo-prof.h"
#include "ofmo-twoint.h"
#include "ofmo-scf.h"

#ifdef USE_MPI
#include <mpi.h>
#else
#include "mpi-dummy.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#else
#include "omp-dummy.h"
#endif

#define ZERO	0.e0
#define TWO	2.e0
#define HALF	.5e0

#define EPS_PS4 1.e-30
#define EPS_ERI 1.e-15
#define EPS_SCH 1.e-12

#ifdef USE_CUDA
#include "cuda/cudalib.h"
#include "cuda/cuda-integ.h"
#endif

static enum ofmo_scf_convType convType = scf;

#ifndef MINIMUM_SCF
#define MINIMUM_SCF 0
#endif
static int minscfcyc = MINIMUM_SCF;

void ofmo_scf_set_convType(enum ofmo_scf_convType type)
{
  convType = type;
}
void ofmo_scf_set_minimum_scf(int mincyc)
{
  minscfcyc = mincyc;
}

// debug
//#include "ofmo-twoint.h"

// debug
static double check_sum2( const int n, const double a[] ) {
    double limit=1.e2, sum;
    int i;
    sum = 0.e0;
    for ( i=0; i<n; i++ ) {
	if ( fabs(a[i]) < limit ) sum += a[i];
    }
    return sum;
}

static int RESERVED_NAO = 0;
static double *_U_    = NULL;
static double *_Dold_ = NULL;
static double *_G_    = NULL;
static double *_F_    = NULL;
static double *_dD_   = NULL;
static double *_dG_   = NULL;
static double *_T_    = NULL;

static double *_C_    = NULL;
static double *_ev_   = NULL;
static int    *_iv_   = NULL;

static double *_Foda_   = NULL;
static double *_Doda_   = NULL;

/** SCF関数内部で利用しているメモリの開放を行う関数
 *
 * SCF関数内部で利用しているメモリの開放を行う関数。
 * 実行終了時に呼ばれるので、ユーザが露に呼ぶ必要はない。
 *
 * */
static void ofmo_scf_dealloc() {
    Free( _U_ );
    Free( _Dold_ );
    Free( _G_ );
    Free( _dD_ );
    Free( _dG_ );
    Free( _T_ );
    Free( _F_ );
    Free( _Foda_ );
    Free( _Doda_ );
    Free( _C_ );
    Free( _ev_ );
    Free( _iv_ );

    RESERVED_NAO = 0;

    _U_    = NULL;
    _Dold_ = NULL;
    _G_    = NULL;
    _dD_   = NULL;
    _dG_   = NULL;
    _T_    = NULL;
    _F_    = NULL;
    _Foda_ = NULL;
    _Doda_ = NULL;
    _C_    = NULL;
    _ev_   = NULL;
    _iv_   = NULL;
}

/** SCF関数内部で利用する配列を確保する関数
 *
 * SCF関数内部で利用する配列を確保する関数。
 * 以下の配列を確保している。
 *
 * @arg \c _U_ 重なり行列をCholesky分解した結果（正方行列）
 * @arg \c _Dold_ SCFステップの１つ前の密度行列（圧縮U形式）
 * @arg \c _G_ ２電子ハミルトン行列を代入する配列（圧縮U形式）
 * @arg \c _dD_ SCFステップの１つ前との密度行列差（圧縮U形式）
 * @arg \c _dG_ 密度行列差を元に計算した２電子ハミルトン行列（圧縮U形式）。
 *     incrementalなFock行列生成時に用いる
 * @arg \c _T_ \c MPI_Allreduce 関数などの利用時に用いる一時配列(正方行列）
 * @arg _C_ MO係数行列が代入される配列（正方行列）
 * @arg _ev_ MOエネルギー（固有値）が代入される配列（ベクトル）
 * @arg _iv_ 固有値計算などで用いる整数ベクトル
 *
 * @param[in] nao AO数
 *
 * @retval  0 正常終了（領域の確保が出来た）
 * @retval -1 異常終了（領域確保に失敗した）
 *
 * @ingroup ofmo-rhf
 * */
int ofmo_scf_init( const int nao ) {
    static int called = false;
    if ( nao > RESERVED_NAO ) {
	int nao2, nnao;
	nao2 = nao * (nao+1) / 2;
	nnao = nao * nao;
	ofmo_scf_dealloc();
	_U_    = (double*)malloc( sizeof(double) * nnao );
	_Dold_ = (double*)malloc( sizeof(double) * nao2 );
	_G_    = (double*)malloc( sizeof(double) * nao2 );
	_dD_   = (double*)malloc( sizeof(double) * nao2 );
	_dG_   = (double*)malloc( sizeof(double) * nao2 );
	_T_    = (double*)malloc( sizeof(double) * nnao );
	_F_    = (double*)malloc( sizeof(double) * nao2 );
	_Foda_ = (double*)malloc( sizeof(double) * nao2 );
	_Doda_ = (double*)malloc( sizeof(double) * nao2 );

	_C_    = (double*)malloc( sizeof(double) * nnao );
	_ev_   = (double*)malloc( sizeof(double) * nao );
	_iv_   = (int*)malloc( sizeof(int) * nao );
	if ( _U_ == NULL || _Dold_ == NULL || _G_ == NULL ||
		_dD_ == NULL || _dG_ == NULL || _T_ == NULL ) {
	    dbg("ERROR : Failure in memory allocation\n");
	    ofmo_scf_dealloc();
	    return -1;
	}
	RESERVED_NAO = nao;
	// information of memory size
	if ( fp_prof != NULL ) {
	    double dsize;
	    dsize = (double)(nao2*7 + nnao*3 + nao )*sizeof(double);
	    dsize += (double)( nao * sizeof(int) );
	    dsize /= (double)(1024*1024);	// MB
	    fprintf(fp_prof,
		    "== allocd memory size in ofmo-scf.c = %10.3f MB\n",
		    dsize );
	}

	if ( !called ) {
	    atexit( ofmo_scf_dealloc );
	    called = true;
	}
    }
    return 0;
}

/** 核間反発エネルギーを求める
 *
 * 核間反発エネルギーを計算する。単位はhartreeである。
 *
 * @param[in] nat 原子数
 * @param[in] atomic_number[] 各原子の原子番号
 * @param[in] atom_x[] 原子のx座標（単位 bohr）
 * @param[in] atom_y[] 原子のy座標（単位 bohr）
 * @param[in] atom_z[] 原子のz座標（単位 bohr）
 *
 * @return 計算された核間反発エネルギー
 *
 * @ingroup ofmo-rhf
 * */
double ofmo_calc_nuclear_repulsion(
	const int nat, const int atomic_number[],
	const double atom_x[], const double atom_y[],
	const double atom_z[] ) {
    double *datomic_number, enuc, rab2;
    int i, j;
    ofmo_scf_init( nat );
    /* Data type conversion of atomic number from int to double.
     * Supposing nao > nat */
    datomic_number = _ev_;
    for ( i=0; i<nat; i++ ) datomic_number[i] = (double)atomic_number[i];
    /* calc. nuclear repulsion energy in au.*/
    enuc = 0.e0;
    for ( i=0; i<nat; i++ ) {
	for ( j=0; j<i; j++ ) {
	    rab2 = (atom_x[i]-atom_x[j]) * (atom_x[i]-atom_x[j])
		 + (atom_y[i]-atom_y[j]) * (atom_y[i]-atom_y[j])
		 + (atom_z[i]-atom_z[j]) * (atom_z[i]-atom_z[j]);
	    enuc += (datomic_number[i]*datomic_number[j] * sqrt(1.e0/rab2));
	}
    }
    return enuc;
}

/** 閉殻系の電子密度行列をMO係数行列から求める
 *
 * 閉殻系の電子密度行列を、与えられたMO係数行列を用いて計算する
 *
 * @attention
 * @li 内部で利用している関数 \c ofmo_transpose_matrix ,
 *     \c ofmo_dot_product
 * @li 内部で利用している一時配列 \c _T_
 *
 * @param[in] nao AO数
 * @param[in] nocc 閉殻軌道数。電子数の半分である。
 * @param[in] C[nao*nao] MO係数行列が代入された配列
 * @param[out] D[nao*(nao+1)/2] 計算された密度行列要素が代入される
 *     配列（圧縮U形式）
 *
 * @ingroup ofmo-rhf
 * */
int ofmo_scf_make_rhf_density(const int nao, const int nocc,
	const double C[], double D[]) {
    double *Ct, *ci, *cj;
    int i, j, ij;
    ofmo_scf_init( nao );
    Ct = _T_;
    ofmo_transpose_matrix(nao, C, Ct);
    ij = 0;
    for (i=0, ci=Ct; i<nao; i++, ci+=nao) {
	for (j=0, cj=Ct; j<=i; j++, cj+=nao) {
	    D[ij++] = ofmo_dot_product(nocc, ci, cj);
	}
    }
    return 0;
}

/** Fock行列を計算して、電子エネルギーを求める
 *
 * 与えられた１電子ハミルトン行列と２電子ハミルトン行列をもとに、
 * Fock行列を計算する。また、それらと密度行列を使って、電子エネルギー
 * を求める。
 *
 * @attention
 * @li 与えられる配列は、すべて、圧縮U形式である
 *
 * @param[in] nao AO数（＝行列サイズ）
 * @param[in] D[] 密度行列
 * @param[in] H[] 一電子ハミルトン行列。\c const 形式ではないが、
 *     実質、この関数呼び出しでは値が変化しない
 * @param[out] F[] Fock行列（H+G）
 *
 * @return 計算された電子エネルギーの値
 *
 * @ingroup ofmo-rhf
 * */
double ofmo_scf_rhf_energy(const int nao, const double D[],
	const double H[], const double F[]) {
    int nao2, i;
    double energy;
    nao2 = nao * (nao+1) /2;
    // E = tr{D(H+F)}
    energy = 0.e0;
    for (i=0; i<nao2; i++) energy += D[i] * (H[i] + F[i]);
    energy *= 2.0;
    return energy;
}


/** Mulliken population解析を行う
 *
 * @param[in] nat 原子数
 * @param[in] nao AO数
 * @param[in] maxlqn 最大軌道量子数
 * @param[in] leading_cs[lqn] 軌道量子数 \c lqn の先頭CS番号
 * @param[in] shel_atm[ics] CS番号 \c ics のCSが属する原子の番号
 * @param[in] shel_ini[ics] CS番号 \c ics のCSに含まれるAOの先頭AO番号
 * @param[in] SP[] 重なり行列（圧縮U形式）
 * @param[in] DP[] 密度行列（圧縮U形式）
 * @param[out] aopop[] 計算されたAO population
 * @param[out] atpop[] 計算されたatomic population
 *
 * @ingroup ofmo-rhf
 * */
int ofmo_scf_mulliken_population(
	const int nat, const int nao, const int maxlqn,
	const int leading_cs[], const int shel_atm[], const int shel_ini[],
	const double SP[], const double DP[],
	double aopop[], double atpop[] ) {
    double *D, *S, *DS;
    int lqn, ics, ics0, ics1, iao, iao0, iao1, atm, nnao;
    int NNAO[] = { 1, 3, 6, 10, 15};
    ofmo_scf_init( nao );
    D  = _U_;
    S  = _C_;
    DS = _T_;
    ofmo_unpack_matrix( nao, SP, S );
    ofmo_unpack_matrix( nao, DP, D );
    ofmo_dgemm( nao, "T", "N", 2.e0, D, S, 0.e0, DS );
    memset( aopop, '\0', sizeof(double) * nao );
    memset( atpop, '\0', sizeof(double) * nat );
    for ( lqn=0; lqn<=maxlqn; lqn++ ) {
	nnao = NNAO[lqn];
	ics0 = leading_cs[lqn];
	ics1 = leading_cs[lqn+1];
	for ( ics=ics0; ics<ics1; ics++ ) {
	    atm  = shel_atm[ics];
	    iao0 = shel_ini[ics];
	    iao1 = iao0 + nnao;
	    for ( iao=iao0; iao<iao1; iao++ ) {
		aopop[iao]  = DS[iao*nao + iao];
		atpop[atm] += aopop[iao];
	    }
	}
    }
    return 0;
}

/** 初期電子密度行列作成ルーチン
 *
 * ソート基底関数の並びで、拡張Huckel法を用いた初期密度行列作成を行う。
 *
 * @param[in] nat 原子数
 * @param[in] ncs CS数
 * @param[in] nao AO数
 * @param[in] maxlqn 最大軌道量子数
 * @param[in] nocc 二重占有軌道数（閉殻系では電子数の半分）
 * @param[in] atomic_number[iat] \c iat 番目の原子の原子番号 (
 *     \f$ 0\le \tt{iat} < \tt{nat} \f$)
 * @param[in] shel_atm[ics] \c ics 番目のCSが属する原子の番号
 * @param[in] shel_ini[ics] \c ics 番目のCSの属するAOの先頭AO番号
 * @param[in] leading_cs[lqn] \c lqn の軌道量子数を持つCSの先頭CS番号
 * @param[in] SP[] ソート済み圧縮U形式の重なり行列
 * @param[out] DP[] 計算されたソート済み圧縮U形式の初期密度行列
 *
 * @note
 * デカルト基底関数系のみに対応している
 *
 * @ingroup ofmo-rhf
 * */
int ofmo_scf_init_density_ehuckel(
	const int nat, const int ncs, const int nao, const int maxlqn,
	const int nocc, const int atomic_number[],
	const int leading_cs[], const int shel_atm[], const int shel_ini[],
	const double SP[], double DP[], double aopop[], double atpop[] ) {
    int lqn;
    int *INI;
    int NNAO[] = {1, 3, 6, 10, 15};
    double *F;
    int atm, cs_tmp=0, iatm, iatn;
    int ics, ics0, ics1, iao, iao0, iao1,  iini, ilqn, i2;
    int jcs, jcs0, jcs1, jao, jao0, jao1,  jini, jlqn, jcs_max, jao_max;
    int ii, ij, jj;
    double val, coe;
    double POTE[] = { 0.0,	/* ダミー要素 */
	0.50,0.90,0.20,0.34,0.30,0.41,0.53,0.50,
	0.64,0.79,0.19,0.28,0.22,0.30,0.39,0.38,0.48,0.58,0.16,0.22,
	0.24,0.25,0.25,0.25,0.27,0.29,0.29,0.28,
	0.28,0.35,0.22,0.29,0.36,0.36,0.43,0.51,0.15,0.21,
	0.23,0.15,0.25,0.26,0.27,0.27,0.27,0.31,
	0.28,0.33,0.21,0.27,0.32,0.33,0.38,0.45,
    };

    ofmo_scf_init( nao );

    INI = _iv_;
    F   = _G_;

    // INIの初期化
    // 対応するCSが、内殻(0)、原子価殻(1)、その他(2以上)のどれかを表す
    for ( lqn=0; lqn<=maxlqn; lqn++ ) {
	ics0 = leading_cs[lqn];
	ics1 = leading_cs[lqn+1];
	atm = -1;
	for ( ics=ics0; ics<ics1; ics++ ) {
	    iatm = shel_atm[ics];
	    if ( iatm != atm ) {
		cs_tmp = 0;
		atm = iatm;
	    }
	    iatn = atomic_number[iatm];
	    if        ( iatn <=2  ) {	// H, He
		INI[ics] = ( lqn==0 ? 1 : 2 );
	    } else if ( iatn <=10 ) {	// Li - Ne
		if ( lqn==0 ) INI[ics] = ( cs_tmp==0 ? 0 : 1 );
		else          INI[ics] = ( lqn==1 ? 1 : 2 );
	    } else if ( iatn <=18 ) {	// Na - Ar
		if      ( lqn==0 ) INI[ics] = ( cs_tmp<=1 ? 0 : 1 );
		else if ( lqn==1 ) INI[ics] = ( cs_tmp==0 ? 0 : 1 );
		else               INI[ics] = 2;
	    } else if ( iatn <= 36 ) {	// K - Kr
		if      ( lqn==0 ) INI[ics] = ( cs_tmp<=2 ? 0 : 1 );
		else if ( lqn==1 ) INI[ics] = ( cs_tmp<=1 ? 0 : 1 );
		else               INI[ics] = 2;
	    } else {
		dbg("atomic number (%d) is not supported, yet\n", iatn );
		fflush(stdout);
		return -1;
	    }
	    cs_tmp++;
	}
    }

    // Fock行列の作成（菊池の方法、拡張Huckelの変形)
    for ( lqn=0; lqn<=maxlqn; lqn++ ) {
	ics0 = leading_cs[lqn];
	ics1 = leading_cs[lqn+1];
	/* 対角成分 */
	for ( ics=ics0; ics<ics1; ics++ ) {
	    iao0 = shel_ini[ics];
	    iao1 = iao0 + NNAO[lqn];
	    iatm = shel_atm[ics];
	    iatn = atomic_number[iatm];
	    iini = INI[ics];
	    if      ( iini == 0 ) val = -2.e0;
	    else if ( iini == 1 ) val = -POTE[iatn];
	    else                  val = 1.e3;
	    for ( iao=iao0; iao<iao1; iao++ )
		F[((iao*iao+iao)>>1) + iao ] = val;
	}
    }
    /* 非対角成分 */
    for ( ilqn=0; ilqn<=maxlqn; ilqn++ ) {
	ics0 = leading_cs[ilqn];
	ics1 = leading_cs[ilqn+1];
	for ( jlqn=0; jlqn<=ilqn; jlqn++ ) {
	    jcs0 = leading_cs[jlqn];
	    jcs1 = leading_cs[jlqn+1];
	    for ( ics=ics0; ics<ics1; ics++ ) {
		iao0 = shel_ini[ics];
		iao1 = iao0 + NNAO[ilqn];
		iini = INI[ics];
		jcs_max = ( jlqn == ilqn ? (ics+1) : jcs1 );
		for ( jcs=jcs0; jcs<jcs_max; jcs++ ) {
		    jao0 = shel_ini[jcs];
		    jao1 = jao0 + NNAO[jlqn];
		    jini = INI[jcs];
		    if      (iini>=2 || jini>=2) coe = 0.0;
		    else if (iini==1 && jini==1) coe = 2.0;
		    else if (iini==0 && jini==0) coe = 0.1;
		    else                         coe = 0.5;
		    for ( iao=iao0; iao<iao1; iao++ ) {
			i2 = (iao*iao+iao)>>1;
			ii = i2 + iao;
			jao_max = (jcs==ics ? (iao+1) : jao1);
			for ( jao=jao0; jao<jao_max; jao++ ) {
			    if ( iao == jao ) continue;
			    ij = i2 + jao;
			    jj = ((jao*jao+jao)>>1) + jao ;
			    F[ij] = coe * SP[ij] * ( F[ii] + F[jj] );
			}
		    }
		}
	    }
	}
    }

    // 一般化固有値問題求値+density作成
    double *U, *EV, *C;
    EV = _ev_;
    C  = _C_;
    U  = _U_;
    ofmo_unpack_matrix(nao, SP, U);
    ofmo_chodec(nao, U);
    ofmo_unpack_matrix(nao, F, C);
    ofmo_solv_GSEP(nao, U, C, EV);
    ofmo_scf_make_rhf_density(nao, nocc, C, DP);
    ofmo_scf_mulliken_population( nat, nao, maxlqn, leading_cs,
	    shel_atm, shel_ini, SP, DP, aopop, atpop );

    return 0;
}

/** テスト用RHFルーチン
 *
 * テスト用のRHF-SCFルーチン
 *
 * @param[in] comm MPIのコミュニケータ
 * @param[in] maxlqn 最大軌道量子数
 * @param[in] ncs CS数
 * @param[in] nao AO数
 * @param[in] leading_cs[lqn] 軌道量子数 \c lqn の先頭CS番号
 * @param[in] shel_atm[ics] CS番号 \c ics のCSが属する原子の番号
 * @param[in] shel_ini[ics] CS番号 \c ics のCSに含まれるAOの先頭AO番号
 * @param[in] atom_x[iat] 原子の番号 \c iat のx座標（au単位）
 * @param[in] atom_y[iat] 原子の番号 \c iat のy座標（au単位）
 * @param[in] atom_z[iat] 原子の番号 \c iat のz座標（au単位）
 * @param[in] leading_cs_pair[itype] CSペアタイプ番号 \c itype の
 *     先頭CSペア番号
 * @param[in] csp_schwarz[icsp] CSペア番号 \c icsp のSchwarz積分
 * @param[in] csp_ics[icsp] CSペア番号 \c icsp の1つ目のCS番号
 * @param[in] csp_jcs[icsp] CSペア番号 \c icsp の2つめのCS番号。ただし、
 *     \f$ \tt{csp\_ics[icsp]} \ge \tt{csp\_jcs[icsp]} \f$ である。
 * @param[in] csp_leading_ps_pair[icsp]  CSペア番号 \c icsp に含まれる
 *     PSペアの先頭PSペア番号
 * @param[in] psp_zeta[ipsp] PSペア番号 \c ipsp の軌道指数和
 *     \f$ \zeta = \zeta_a + \zeta_b \f$
 * @param[in] psp_dkps[ipsp] PSペア番号 \c ipsp の線型結合定数
 *     \f[ K_{ab} = \sqrt2 \pi^{5/4} \frac1{\zeta_a+\zeta_b}
 *     \exp\left[ -\frac{\zeta_a \zeta_b}{\zeta_a + \zeta_b}
 *     ( \boldmath A \unboldmath - \boldmath B \unboldmath )^2
 *     \right]\f]
 * @param[in] psp_xiza[ipsp] PSペア番号 \c ipsp の
 *     \f$ \frac{\xi}{\zeta_a} = \frac{\zeta_b}{\zeta_a+\zeta_b} \f$
 *
 * @param[in] nat 原子数
 * @param[in] atomic_number[iat] 原子の番号 \c iat の原子番号
 * @param[in] ncharge 分子の電荷
 *
 * @param[in] S[] 重なり行列（圧縮U形式）
 * @param[in] H[] 一電子ハミルトン行列（圧縮U形式）
 *
 * @param[in] maxscfcyc 最大SCF繰り返し回数
 * @param[in] scfe エネルギーの収束条件
 * @param[in] scfd 密度行列の収束条件
 * 
 * @param[in,out] （入力時）初期密度行列／（出力時）SCF収束時の密度行列
 *     （圧縮U形式）
 * @param[out] F[] SCF収束時のソート済みFock行列（圧縮U形式）
 * @param[out] C[] SCF収束時のソート済みMO係数行列（正方行列）
 * @param[out] moe[] SCF収束時のMOエネルギー（ベクトル）
 * @param[out] *Eelec SCF収束時の電子エネルギー（hartree）
 *
 * @retval  0 正常終了時（SCFが収束した）
 * @retval -1 異常終了時（SCFが収束しなかった）
 *
 * @ingroup ofmo-rhf
 * */
int ofmo_scf_rhf(
	MPI_Comm comm, const int maxlqn, const double Enuc,
	const int ncs, const int nao,
	const int leading_cs[],
	// 基底関数データ
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	// カットオフテーブルデータ
	const int leading_cs_pair[], const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// 分子データ
	const int nat,
	const int nocc,
	// 積分データ
	double S[], double H[],
	// 制御用データ
	const int maxscfcyc, const double scfe, const double scfd,
	// 結果代入用データ
	double D[], double C[], double moe[], double *Etot) {
    int nao2;
    double *U, *Dold, *dG, *dD, *TMP, *F;
    double Eold, Enew, deltaD = 1.e100, deltaE;
    int itera, myrank, nprocs;
    double errdiis=0.0e0;
    // 様々な閾値関連変数
    double tol_diis = 1.e-2;
    // 様々なフラグ
    int dodiis=false, flag_scf;
    //int doshift=false;
    int ierr;
    // ODA
    int koda;
    double *Foda, *Doda, Eoda, lambda1;
    // シフト演算子関連
    //double shifto = 1.0, shiftv =1.0;
    int cid_gmat;
    int tid_init, tid_gmat, tid_allred, tid_diis, tid_diag, tid_bcast;
    int tid_tot;
    // Screening threshold
#if 0
    float eps_ps4 = EPS_PS4;
    float eps_eri = EPS_ERI;
    float eps_sch = EPS_SCH;
    float eps_fac = (scfe>1e-8)? scfe/1e-8: 1.0;
    eps_ps4 *= eps_fac;
    eps_eri *= eps_fac;
    eps_sch *= eps_fac;
    if(fp_prof) fprintf(fp_prof, "scfd scfe eps_sch: %e %e %e \n", scfd, scfe, eps_sch);
#else
    float eps0_ps4=ofmo_twoint_eps_ps4(0);
    float eps0_eri=ofmo_twoint_eps_eri(0);
    float eps0_sch=ofmo_twoint_eps_sch(0);
    float eps_ps4 = eps0_ps4;
    float eps_eri = eps0_eri;
    float eps_sch = eps0_sch;
#endif
    int maxitera = maxscfcyc;
    int minitera = minscfcyc;

    cid_gmat = ofmo_create_thread_timer( "GMAT", 0 );
    tid_init = ofmo_create_proc_timer( "INIT", 2 );
    tid_gmat = ofmo_create_proc_timer( "GMAT", 2 );
    tid_allred = ofmo_create_proc_timer( "ALLRED", 2 );
    tid_diis   = ofmo_create_proc_timer( "DIIS", 2 );
    tid_diag   = ofmo_create_proc_timer( "DIAG", 2 );
    tid_bcast  = ofmo_create_proc_timer( "BCAST", 2 );
    tid_tot    = ofmo_create_proc_timer( "TOT", 2 );

    MPI_Comm_rank( comm, &myrank );
    MPI_Comm_size( comm, &nprocs );

    ofmo_start_proc_timer( tid_init );
    ofmo_start_proc_timer( tid_tot );

    ofmo_scf_init( nao );

    nao2 = nao * (nao+1) / 2;
    U    = _U_;
    Dold = _Dold_;
    //G    = _G_;
    dG   = _dG_;
    dD   = _dD_;
    TMP  = _T_;
    F    = _F_;
    Foda = _Foda_;
    Doda = _Doda_;
    // 重なり行列のCholesky分解
    ofmo_unpack_matrix(nao, S, U);

    ierr = ofmo_chodec(nao, U);

    if (ierr != 0) {
	if (myrank == 0)
	    dbg("ERROR: Failure in Cholesky decomposition of S\n");
	return -1;
    }
    // DIIS用にS=U'UのUを保存
    ofmo_diis_init( nao, S );
    // 初期設定
    //Eold = ZERO;
    Eold = *Etot;
    memset( Dold, '\0', sizeof(double)*nao2 );

    ofmo_dcopy( nao2, D, dD );	/* for incremental Fock */

    ofmo_twoint_alloc_Dcs(ncs); /* for Schwarz screening with Density */

    int incore = -1;
    int last_incore = -1;
#pragma omp parallel
    {
      int mythread = omp_get_thread_num();
      int incore2 = ofmo_twoint_get_last_eri_type( mythread );
#pragma omp critical
      if (incore2 > incore) incore = incore2;
    }
    MPI_Allreduce(&incore, &last_incore, 1, MPI_INT, MPI_MAX, comm);
    ofmo_twoint_set_global_last_eri_type(last_incore);
    if(fp_prof) fprintf(fp_prof, "g_incore: %d\n", last_incore);

#ifdef USE_CUDA
    {
      int nblk = -1, nthb = -1;
      int ret = 0;
      int Lab = maxlqn*(maxlqn+1)/2 + maxlqn;
      int ncspair = leading_cs_pair[Lab+1];
      int npspair = csp_leading_ps_pair[ncspair];
      float *csp_schwarz_f;
      csp_schwarz_f = (float*)malloc( sizeof(float) * ncspair );
      for (int i=0; i<ncspair; i++) csp_schwarz_f[i]=(float)csp_schwarz[i];
      ret = cuda_SCF_Init(ncs, nat, maxlqn, ncspair,npspair, nao,
          shel_atm, shel_ini, atom_x, atom_y, atom_z,
          leading_cs,
          leading_cs_pair, csp_leading_ps_pair, csp_ics, csp_jcs,
          psp_zeta, psp_dkps, psp_xiza, csp_schwarz_f,
          nblk, nthb);
      if (ret<0) exit(ret);
      free(csp_schwarz_f);
    }
#endif

    ofmo_acc_proc_timer( tid_init );

    // ODA
    koda = 0;
    //for (itera=1, flag_scf=false; itera<=maxscfcyc; itera++) {
    for (itera=1, flag_scf=false; itera<=maxitera; itera++) {
	ofmo_start_proc_timer( tid_gmat );
        ofmo_twoint_eps_ps4(eps_ps4);
        ofmo_twoint_eps_eri(eps_eri);
        ofmo_twoint_eps_sch(eps_sch);
#pragma omp parallel
	{
	    int nworkers, workerid;
	    int nthreads, mythread;
	    nthreads = omp_get_num_threads();
	    mythread = omp_get_thread_num();
	    nworkers = nthreads * nprocs;
	    workerid = myrank * nthreads + mythread;
	    ofmo_start_thread_timer( cid_gmat, mythread );
	    ofmo_integ_gen_gmat( nworkers, workerid,
		    maxlqn, shel_atm, shel_ini, atom_x, atom_y, atom_z,
		    leading_cs_pair, leading_cs,
		    csp_schwarz, csp_ics, csp_jcs, csp_leading_ps_pair,
		    psp_zeta, psp_dkps, psp_xiza, nao, dD, dG );
	    ofmo_acc_thread_timer( cid_gmat, mythread );
	}
	ofmo_acc_proc_timer( tid_gmat );
	ofmo_start_proc_timer( tid_allred );
	MPI_Allreduce( dG, TMP, nao2, MPI_DOUBLE, MPI_SUM, comm );
	memcpy( dG, TMP, sizeof(double)*nao2 );
	ofmo_acc_proc_timer( tid_allred );
	/*// debug
	if ( itera== 1 && nao > 80 ) {
	    long nzeri, nzeri0;
	    nzeri0 = ofmo_get_nonzero_eri();
	    MPI_Reduce( &nzeri0, &nzeri, 1, MPI_LONG, MPI_SUM, 0, comm );
	    if ( fp_prof ) {
		fprintf( fp_prof, "nzeri= %12ld\n", nzeri );
	    }
	}*/

	/* form total Fock matrix */
	if ( itera == 1 )
	    for ( int i=0; i<nao2; i++ ) F[i] = H[i] + dG[i];
	else
	    for ( int i=0; i<nao2; i++ ) F[i] += dG[i];
	// Fock行列の計算、エネルギーの計算
	Enew = ofmo_scf_rhf_energy(nao, D, H, F) + Enuc;

	/* DIISをやるかどうかの判定 */
	ofmo_start_proc_timer( tid_diis );
	ofmo_scale_diag( nao, 2.e0, F );
	if ( (errdiis=ofmo_diis_profiling( nao, D, F )) > 0 ) {
	    if ( dodiis == false && fabs(Enew-Eold)<tol_diis ) {
		dodiis = true;
	    }
	} else {
	    dodiis = false;
	}
	double *FE = _T_;
	if ( dodiis ) ofmo_diis_update( nao, D, F, FE, dodiis );
	/* ODA */
	if ( dodiis ) {
	    koda = 0;
	} else {
	    if ( koda == 0 ) {
		memcpy( Foda, F, sizeof(double)*nao2 );
		Eoda = Enew;
	    } else {
		double s2, c, lambda;
		int i;
		ofmo_scale_diag( nao, 0.5e0, Doda );

		s2 = c = 0.e0;
		for ( i=0; i<nao2; i++ ) {
		    s2 += Foda[i] * Doda[i];
		    c  += (F[i]-Foda[i]) * Doda[i];
		}
		c *= 2.e0;
		ofmo_scale_diag( nao, 2.e0, Doda );
		lambda = ( c <= -s2 ? 1.e0 : (-s2/c) );
		lambda1 = 1.e0 - lambda;
		Eoda += lambda*(s2*2.e0+lambda*c);
		for ( i=0; i<nao2; i++ )
		    Foda[i]=lambda1*Foda[i]+lambda*F[i];
	    }
	    koda++;
	}

	double *Fd;
	Fd = ( dodiis ? FE : Foda );
	ofmo_acc_proc_timer( tid_diis );

	ofmo_unpack_matrix(nao, Fd, C);/* solve generalized symetric
					   eigenvalue problem */
	ofmo_scale_diag( nao, 0.5e0, F );

	ofmo_start_proc_timer( tid_diag );
	ierr = ofmo_solv_GSEP( nao, U, C, moe );
	ofmo_acc_proc_timer( tid_diag );
	if ( ierr != 0 ) {
	    dbg("error in solving GSEP\n");
	    MPI_Abort( MPI_COMM_WORLD, -1 );
	}
	/* update density */
	ofmo_dcopy( nao2, D, Dold ); /* store previous density matrix */
	ofmo_scf_make_rhf_density(nao, nocc, C, D);
	// 収束判定を行う
	deltaE = Enew - Eold;
	deltaD = ofmo_max_diff( nao2, D, Dold );
	if ( fp_prof ) {
	    fprintf(fp_prof, "%5d : %17.10f ( %17.10f ) (%10.7f)\n",
		    itera, Enew, deltaE, deltaD);
	}
        if (itera>=minitera) {
          if (convType==scc) {
            if (deltaE<scfe && deltaD<scfd) {
              flag_scf = true;
            }
          } else if (convType==scf) {
            // from rhfuhf in GAMESS
            int cvdens = false;
            int cvengy = false;
            int cvdiis = false;
            double diistol = scfe*1.0e-2;
            cvdens = (deltaD<scfd && fabs(deltaE)<scfe*10)
              || (deltaD<scfd*0.2e0);
            cvengy = (fabs(deltaE)<scfe && deltaD<scfd*2);
            errdiis = fabs(errdiis);
            cvdiis = (dodiis && errdiis<diistol && deltaD<scfd*2);
            if (cvdens || cvengy || cvdiis) { flag_scf = true; };
          } else {
            if (fabs(deltaE)<scfe && deltaD<scfd) {
              flag_scf = true;
            }
          }
        }

	ofmo_start_proc_timer( tid_bcast );
	MPI_Bcast( &flag_scf, 1, MPI_INT, 0, comm );
	ofmo_acc_proc_timer( tid_bcast );
	if ( flag_scf == true ) break;

	// for incremental Fock generation
	for ( int i=0; i<nao2; i++ ) dD[i] = D[i] - Dold[i];
	Eold = Enew;
	/* ODA */
	if ( koda == 1 ) {
	    memcpy( Doda, dD, sizeof(double)*nao2 );
	} else if ( koda > 1 ) {
	    for ( int i=0; i<nao2; i++ )
		Doda[i] = dD[i] + lambda1*Doda[i];
	}
    }	// end for (itera)
    MPI_Bcast( &Enew, 1, MPI_DOUBLE, 0, comm );
    *Etot = Enew;

    ofmo_twoint_free_Dcs();

#ifdef USE_CUDA
    {
      int ret = 0;
      ret = cuda_SCF_Finalize();
      if (ret<0) exit(ret);
    }
#endif
    ofmo_acc_proc_timer( tid_tot );

    ofmo_twoint_eps_ps4(eps0_ps4);
    ofmo_twoint_eps_eri(eps0_eri);
    ofmo_twoint_eps_sch(eps0_sch);

    // 最終処理
    if (flag_scf == true) {
	//if (myrank == 0) printf("==== Allready SCF ====\n");
	return 0;
    } else {
	if (myrank == 0 && fp_prof != NULL)
	    fprintf(fp_prof, "==== SCF not converged ====\n");
	return -1;
    }
}
