/**
 * @file ofmo-cutoff.c
 * ２電子積分、４中心クーロン積分、および、３中心クーロン積分の
 * 計算の際に用いるカットオフテーブルの計算を行う
 * */
/**
 * @defgroup integ-cutoff カットオフテーブル計算を行う関数群
 * この関数では、各CSペアタイプのカットオフテーブル計算ルーチンが
 * 定義されている。
 *
 * 初期化関数を除いて、すべて同じ引数を持つ。以下に、引数とその内容を記す。
 *
 * @param[in] La CSペアのうち、１つ目のCSの軌道量子数
 * @param[in] Lb CSペアのうち、 ２つ目のCSの軌道量子数
 *     （\f$ \tt{La} \ge \tt{Lb} \ge 0 \f$）
 * @param[in] leading_cs[lqn] 軌道量子数 /c lqn の先頭CS番号
 * @param[in] shel_tem[ics] CS番号 \c ics のCSの縮約長
 * @param[in] shel_atm[ics] CS番号 \c ics のCSが属する原子の番号
 * @param[in] shel_add[ics] CS番号 \c ics のCSに含まれるPSの先頭PS番号
 * @param[in] atom_x[iat] 原子の番号 \c iat のx座標（au単位）
 * @param[in] atom_y[iat] 原子の番号 \c iat のy座標（au単位）
 * @param[in] atom_z[iat] 原子の番号 \c iat のz座標（au単位）
 * @param[in] prim_exp[ips] PS番号 \c ips のPSの軌道指数
 * @param[in] prim_coe[ips] PS番号 \c ips のPSの規格化定数込みの縮約係数
 *
 * @param[out] leading_cs_pair[itype] CSペアタイプ番号 \c itype の
 *     先頭CSペア番号
 * @param[out] csp_schwarz[icsp] CSペア番号 \c icsp のSchwarz積分
 * @param[out] csp_ics[icsp] CSペア番号 \c icsp の1つ目のCS番号
 * @param[out] csp_jcs[icsp] CSペア番号 \c icsp の2つめのCS番号。ただし、
 *     \f$ \tt{csp\_ics[icsp]} \ge \tt{csp\_jcs[icsp]} \f$ である。
 * @param[out] csp_leading_ps_pair[icsp]  CSペア番号 \c icsp に含まれる
 *     PSペアの先頭PSペア番号
 * @param[out] psp_zeta[ipsp] PSペア番号 \c ipsp の軌道指数和
 *     \f$ \zeta = \zeta_a + \zeta_b \f$
 * @param[out] psp_dkps[ipsp] PSペア番号 \c ipsp の線型結合定数
 *     \f[ K_{ab} = \sqrt2 \pi^{5/4} \frac1{\zeta_a+\zeta_b}
 *     \exp\left[ -\frac{\zeta_a \zeta_b}{\zeta_a + \zeta_b}
 *     ( \boldmath A \unboldmath - \boldmath B \unboldmath )^2
 *     \right]\f]
 * @param[out] psp_xiza[ipsp] PSペア番号 \c ipsp の
 *     \f$ \frac{\xi}{\zeta_a} = \frac{\zeta_b}{\zeta_a+\zeta_b} \f$
 * 
 * @ingroup integ-med
 * */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <sys/time.h>

#include "ofmo-def.h"
#include "fmt.h"

#define MAXPSPAIR	100	// １つのCSペアあたりの最大PSペア数
				// MAXTEMの２乗でもいいけど、staticに
				// 確保したいので、このように定義する
#define EPS_CS_PAIR	1.e-15
#define EPS_CS_PAIR2	1.e-30
#define EPS_PS_PAIR	1.e-32
#define MAXTEM2		MAXPSPAIR
#define ZERO	0.e0

#include "ofmo-cutoff-core.h"
#include "ofmo-cutoff.h"


static double _CK_ = ZERO;

/** カットオフテーブル作成プログラムの初期化ルーチン
 *
 * @ingroup integ-cutoff
 * */
int ofmo_cutoff_init() {
    static int called = false;
    if ( called ) return 0;
    double pi, t;
#ifdef M_PI
    pi = M_PI;
#else
    pi = 4.e0 * atan( 1.e0 );
#endif
    t = 2.e0 * pi * pi * sqrt( pi );
    _CK_ = sqrt( t );
    called = true;
    return 0;
}

static int schwarz_calc_ps_pair_params(
	const double prim_exp[], const double prim_coe[],
	const int ips0, const int ips1,
	const int jps0, const int jps1,
	const double AB2,
	double vzeta[], double vdkps[], double vxiza[] ) {
    double zeta_a, coef_a, zeta_b, coef_b;
    double sqrz, rz, xiza, Kab, zeta, coef;
    int npps, ips, jps;
    npps = 0;
    for ( ips=ips0; ips<ips1; ips++ ) {
	zeta_a = prim_exp[ips];
	coef_a = prim_coe[ips];
	for ( jps=jps0; jps<jps1; jps++ ) {
	    zeta_b = prim_exp[jps];
	    coef_b = prim_coe[jps];

	    zeta = zeta_a + zeta_b;
	    coef = coef_a * coef_b;

	    sqrz = 1.e0 / sqrt( zeta );
	    rz   = sqrz * sqrz;
	    xiza = zeta_b * rz;
	    Kab  = _CK_ * coef * sqrz * rz * exp( -xiza * zeta_a * AB2 );
	    if ( fabs(Kab) > EPS_PS_PAIR ) {
		vzeta[npps] = rz;
		vdkps[npps] = Kab;
		vxiza[npps] = xiza;
		npps++;
	    }
	}
    }
    return npps;
}

/** SSタイプのカットオフテーブルを作成する
 * @ingroup integ-cutoff
 * */
int ofmo_cutoff_ss_(
	// input arguments
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[], const int shel_add[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	// output arguments
	int leading_cs_pair[],
	double csp_schwarz[], int csp_ics[], int csp_jcs[],
	int csp_leading_ps_pair[],
	double psp_zeta[], double psp_dkps[], double psp_xiza[] ) {
    int ics, ics0, ics1, ips0, ips1, iat;
    int jcs, jcs0,       jps0, jps1, jat;
    int npps, ipps, ncs_pair, nps_pair;
    int i, Lab;
    double BA[3], A[3], B[3], AB2, max_eri;
    double vzeta[MAXTEM2], vdkps[MAXTEM2], vxiza[MAXTEM2];
    int La=*pLa, Lb=*pLb;

    Lab       = La*(La+1)/2 + Lb;
    ncs_pair  = leading_cs_pair[Lab];
    nps_pair  = csp_leading_ps_pair[ncs_pair];

    ics0 = leading_cs[La];	// La=0
    ics1 = leading_cs[La+1];
    jcs0 = leading_cs[Lb];	// Lb=0
    for ( ics=ics0; ics<ics1; ics++ ) {
	ips0 = shel_add[ics];
	ips1 = ips0 + shel_tem[ics];
	iat  = shel_atm[ics];
	A[0]=atom_x[ iat ]; A[1]=atom_y[ iat ]; A[2]=atom_z[ iat ];
	for ( jcs=jcs0; jcs<=ics; jcs++ ) {
	    jps0 = shel_add[jcs];
	    jps1 = jps0 + shel_tem[jcs];
	    jat  = shel_atm[jcs];
	    B[0]=atom_x[ jat ]; B[1]=atom_y[ jat ]; B[2]=atom_z[ jat ];

	    AB2  = ZERO;
	    for ( i=0; i<3; i++ ) {
		BA[i] = B[i] - A[i];
		AB2  += BA[i]*BA[i];
	    }
	    
	    npps = schwarz_calc_ps_pair_params(
		    prim_exp, prim_coe, ips0, ips1, jps0, jps1,AB2,
		    vzeta, vdkps, vxiza );
	    if ( npps == 0 ) continue;
	    
	    max_eri = schwarz_core_ssss_(
		    npps, vzeta, vdkps, vxiza, BA, AB2 );

	    if ( max_eri > EPS_CS_PAIR2 ) {
		csp_schwarz[ncs_pair] = sqrt( max_eri );
		csp_ics[ncs_pair] = ics;
		csp_jcs[ncs_pair] = jcs;
		for ( ipps=0; ipps<npps; ipps++ ) {
		    psp_zeta[nps_pair] = vzeta[ipps];
		    psp_dkps[nps_pair] = vdkps[ipps];
		    psp_xiza[nps_pair] = vxiza[ipps];
		    nps_pair++;
		}
		csp_leading_ps_pair[ncs_pair+1] = nps_pair;
		ncs_pair++;
	    }
	}	// jcs
    }		// ics
    leading_cs_pair[Lab+1] = ncs_pair;
    /*
    // debug
    printf("(SS) : ncs_pair=%d, nps_pair=%d\n", ncs_pair, nps_pair );
    fflush(stdout);
    */
    /*
    // debug
    printf("(s|s) : nps_pair = %d\n", nps_pair );
    */
    return 0;
}

/** PSタイプのカットオフテーブルを作成する
 * @ingroup integ-cutoff
 * */
int ofmo_cutoff_ps_(
	// input arguments
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[], const int shel_add[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	// output arguments
	int leading_cs_pair[],
	double csp_schwarz[], int csp_ics[], int csp_jcs[],
	int csp_leading_ps_pair[],
	double psp_zeta[], double psp_dkps[], double psp_xiza[] ) {
    int ics, ics0, ics1, ips0, ips1, iat;
    int jcs, jcs0, jcs1, jps0, jps1, jat;
    int npps, ipps, ncs_pair, nps_pair;
    int i, Lab;
    double BA[3], A[3], B[3], AB2;
    double max_eri;
    double vzeta[MAXTEM2], vdkps[MAXTEM2], vxiza[MAXTEM2];
    int La=*pLa, Lb=*pLb;

    Lab       = La*(La+1)/2 + Lb;
    ncs_pair  = leading_cs_pair[Lab];
    nps_pair  = csp_leading_ps_pair[ncs_pair];

    ics0 = leading_cs[La];	// La=1
    ics1 = leading_cs[La+1];
    jcs0 = leading_cs[Lb];	// Lb=0
    jcs1 = leading_cs[Lb+1];
    for ( ics=ics0; ics<ics1; ics++ ) {
	ips0 = shel_add[ics];
	ips1 = ips0 + shel_tem[ics];
	iat  = shel_atm[ics];
	A[0]=atom_x[ iat ]; A[1]=atom_y[ iat ]; A[2]=atom_z[ iat ];
	for ( jcs=jcs0; jcs<jcs1; jcs++ ) {
	    jps0 = shel_add[jcs];
	    jps1 = jps0 + shel_tem[jcs];
	    jat  = shel_atm[jcs];
	    B[0]=atom_x[ jat ]; B[1]=atom_y[ jat ]; B[2]=atom_z[ jat ];

	    AB2 = ZERO;
	    for ( i=0; i<3; i++ ) {
		BA[i] = B[i] - A[i];
		AB2  += BA[i]*BA[i];
	    }
	    
	    npps = schwarz_calc_ps_pair_params(
		    prim_exp, prim_coe, ips0, ips1, jps0, jps1, AB2,
		    vzeta, vdkps, vxiza );
	    if ( npps == 0 ) continue;
	    
	    max_eri = schwarz_core_psps_(
		    npps, vzeta, vdkps, vxiza, BA, AB2 );

	    if ( max_eri > EPS_CS_PAIR2 ) {
		csp_schwarz[ncs_pair] = sqrt( max_eri );
		csp_ics[ncs_pair] = ics;
		csp_jcs[ncs_pair] = jcs;
		for ( ipps=0; ipps<npps; ipps++ ) {
		    psp_zeta[nps_pair] = vzeta[ipps];
		    psp_dkps[nps_pair] = vdkps[ipps];
		    psp_xiza[nps_pair] = vxiza[ipps];
		    nps_pair++;
		}
		csp_leading_ps_pair[ncs_pair+1] = nps_pair;
		ncs_pair++;
	    }
	}	// jcs
    }		// ics
    leading_cs_pair[Lab+1] = ncs_pair;
    /*
    // debug
    printf("(PS) : ncs_pair=%d, nps_pair=%d\n", ncs_pair, nps_pair );
    fflush(stdout);
    */
    return 0;
}

/** PPタイプのカットオフテーブルを作成する
 * @ingroup integ-cutoff
 * */
int ofmo_cutoff_pp_(
	// input arguments
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[], const int shel_add[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	// output arguments
	int leading_cs_pair[],
	double csp_schwarz[], int csp_ics[], int csp_jcs[],
	int csp_leading_ps_pair[],
	double psp_zeta[], double psp_dkps[], double psp_xiza[] ) {
    int ics, ics0, ics1, ips0, ips1, iat;
    int jcs, jcs0,       jps0, jps1, jat;
    int npps, ipps, ncs_pair, nps_pair;
    int i, Lab;
    double BA[3], A[3], B[3], AB2;
    double max_eri;
    double vzeta[MAXTEM2], vdkps[MAXTEM2], vxiza[MAXTEM2];
    int La=*pLa, Lb=*pLb;

    Lab       = La*(La+1)/2 + Lb;
    ncs_pair  = leading_cs_pair[Lab];
    nps_pair  = csp_leading_ps_pair[ncs_pair];

    ics0 = leading_cs[La];	// La=1
    ics1 = leading_cs[La+1];
    jcs0 = leading_cs[Lb];	// Lb=1
    for ( ics=ics0; ics<ics1; ics++ ) {
	ips0 = shel_add[ics];
	ips1 = ips0 + shel_tem[ics];
	iat  = shel_atm[ics];
	A[0]=atom_x[ iat ]; A[1]=atom_y[ iat ]; A[2]=atom_z[ iat ];
	for ( jcs=jcs0; jcs<=ics; jcs++ ) {
	    jps0 = shel_add[jcs];
	    jps1 = jps0 + shel_tem[jcs];
	    jat  = shel_atm[jcs];
	    B[0]=atom_x[ jat ]; B[1]=atom_y[ jat ]; B[2]=atom_z[ jat ];

	    AB2  = ZERO;
	    for ( i=0; i<3; i++ ) {
		BA[i] = B[i] - A[i];
		AB2  += BA[i]*BA[i];
	    }
	    
	    npps = schwarz_calc_ps_pair_params(
		    prim_exp, prim_coe, ips0, ips1, jps0, jps1, AB2,
		    vzeta, vdkps, vxiza );
	    if ( npps == 0 ) continue;

	    max_eri = schwarz_core_pppp_(
		    npps, vzeta, vdkps, vxiza, BA, AB2 );
	    
	    if ( max_eri > EPS_CS_PAIR2 ) {
		csp_schwarz[ncs_pair] = sqrt( max_eri );
		csp_ics[ncs_pair] = ics;
		csp_jcs[ncs_pair] = jcs;
		for ( ipps=0; ipps<npps; ipps++ ) {
		    psp_zeta[nps_pair] = vzeta[ipps];
		    psp_dkps[nps_pair] = vdkps[ipps];
		    psp_xiza[nps_pair] = vxiza[ipps];
		    nps_pair++;
		}
		csp_leading_ps_pair[ncs_pair+1] = nps_pair;
		ncs_pair++;
	    }
	}	// jcs
    }		// ics
    leading_cs_pair[Lab+1] = ncs_pair;
    /*
    // debug
    printf("(PP) : ncs_pair=%d, nps_pair=%d\n", ncs_pair, nps_pair );
    fflush(stdout);
    */
    return 0;
}

/** DSタイプのカットオフテーブルを作成する
 * @ingroup integ-cutoff
 * */
int ofmo_cutoff_ds_(
	// input arguments
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[], const int shel_add[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	// output arguments
	int leading_cs_pair[],
	double csp_schwarz[], int csp_ics[], int csp_jcs[],
	int csp_leading_ps_pair[],
	double psp_zeta[], double psp_dkps[], double psp_xiza[] ) {
    int ics, ics0, ics1, ips0, ips1, iat;
    int jcs, jcs0, jcs1, jps0, jps1, jat;
    int npps, ipps, ncs_pair, nps_pair;
    int i, Lab;
    double BA[3], A[3], B[3], AB2;
    double max_eri;
    double vzeta[MAXPSPAIR], vdkps[MAXPSPAIR], vxiza[MAXPSPAIR];
    int La=*pLa, Lb=*pLb;

    Lab       = La*(La+1)/2 + Lb;
    ncs_pair  = leading_cs_pair[Lab];
    nps_pair  = csp_leading_ps_pair[ncs_pair];

    ics0 = leading_cs[La];	// La=2
    ics1 = leading_cs[La+1];
    jcs0 = leading_cs[Lb];	// Lb=0
    jcs1 = leading_cs[Lb+1];
    for ( ics=ics0; ics<ics1; ics++ ) {
	ips0 = shel_add[ics];
	ips1 = ips0 + shel_tem[ics];
	iat  = shel_atm[ics];
	A[0]=atom_x[ iat ]; A[1]=atom_y[ iat ]; A[2]=atom_z[ iat ];
	for ( jcs=jcs0; jcs<jcs1; jcs++ ) {
	    jps0 = shel_add[jcs];
	    jps1 = jps0 + shel_tem[jcs];
	    jat  = shel_atm[jcs];
	    B[0]=atom_x[ jat ]; B[1]=atom_y[ jat ]; B[2]=atom_z[ jat ];

	    AB2  = ZERO;
	    for ( i=0; i<3; i++ ) {
		BA[i] = B[i] - A[i];
		AB2  += BA[i]*BA[i];
	    }
	    
	    npps = schwarz_calc_ps_pair_params(
		    prim_exp, prim_coe, ips0, ips1, jps0, jps1, AB2,
		    vzeta, vdkps, vxiza );
	    if ( npps == 0 ) continue;
	    
	    max_eri = schwarz_core_dsds_(
		    npps, vzeta, vdkps, vxiza, BA, AB2 );

	    if ( max_eri > EPS_CS_PAIR2 ) {
		csp_schwarz[ncs_pair] = sqrt( max_eri );
		csp_ics[ncs_pair] = ics;
		csp_jcs[ncs_pair] = jcs;
		for ( ipps=0; ipps<npps; ipps++ ) {
		    psp_zeta[nps_pair] = vzeta[ipps];
		    psp_dkps[nps_pair] = vdkps[ipps];
		    psp_xiza[nps_pair] = vxiza[ipps];
		    nps_pair++;
		}
		csp_leading_ps_pair[ncs_pair+1] = nps_pair;
		ncs_pair++;
	    }
	}	// jcs
    }		// ics
    leading_cs_pair[Lab+1] = ncs_pair;
    return 0;
}

/** DPタイプのカットオフテーブルを作成する
 * @ingroup integ-cutoff
 * */
int ofmo_cutoff_dp_(
	// input arguments
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[], const int shel_add[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	// output arguments
	int leading_cs_pair[],
	double csp_schwarz[], int csp_ics[], int csp_jcs[],
	int csp_leading_ps_pair[],
	double psp_zeta[], double psp_dkps[], double psp_xiza[] ) {
    int ics, ics0, ics1, ips0, ips1, iat;
    int jcs, jcs0, jcs1, jps0, jps1, jat;
    int npps, ipps, ncs_pair, nps_pair;
    int i, Lab;
    double BA[3], A[3], B[3], AB2;
    double max_eri;
    double vzeta[MAXPSPAIR], vdkps[MAXPSPAIR], vxiza[MAXPSPAIR];
    int La=*pLa, Lb=*pLb;

    Lab       = La*(La+1)/2 + Lb;
    ncs_pair  = leading_cs_pair[Lab];
    nps_pair  = csp_leading_ps_pair[ncs_pair];

    ics0 = leading_cs[La];	// La=2
    ics1 = leading_cs[La+1];
    jcs0 = leading_cs[Lb];	// Lb=1
    jcs1 = leading_cs[Lb+1];
    for ( ics=ics0; ics<ics1; ics++ ) {
	ips0 = shel_add[ics];
	ips1 = ips0 + shel_tem[ics];
	iat  = shel_atm[ics];
	A[0]=atom_x[ iat ]; A[1]=atom_y[ iat ]; A[2]=atom_z[ iat ];
	for ( jcs=jcs0; jcs<jcs1; jcs++ ) {
	    jps0 = shel_add[jcs];
	    jps1 = jps0 + shel_tem[jcs];
	    jat  = shel_atm[jcs];
	    B[0]=atom_x[ jat ]; B[1]=atom_y[ jat ]; B[2]=atom_z[ jat ];

	    AB2  = ZERO;
	    for ( i=0; i<3; i++ ) {
		BA[i] = B[i] - A[i];
		AB2  += BA[i]*BA[i];
	    }
	    
	    npps = schwarz_calc_ps_pair_params(
		    prim_exp, prim_coe, ips0, ips1, jps0, jps1, AB2,
		    vzeta, vdkps, vxiza );
	    if ( npps == 0 ) continue;
	    

	    max_eri = schwarz_core_dpdp_(
		    npps, vzeta, vdkps, vxiza, BA, AB2 );

	    if ( max_eri > EPS_CS_PAIR2 ) {
		csp_schwarz[ncs_pair] = sqrt( max_eri );
		csp_ics[ncs_pair] = ics;
		csp_jcs[ncs_pair] = jcs;
		for ( ipps=0; ipps<npps; ipps++ ) {
		    psp_zeta[nps_pair] = vzeta[ipps];
		    psp_dkps[nps_pair] = vdkps[ipps];
		    psp_xiza[nps_pair] = vxiza[ipps];
		    nps_pair++;
		}
		csp_leading_ps_pair[ncs_pair+1] = nps_pair;
		ncs_pair++;
	    }
	}	// jcs
    }		// ics
    leading_cs_pair[Lab+1] = ncs_pair;
    return 0;
}

/** DDタイプのカットオフテーブルを作成する
 * @ingroup integ-cutoff
 * */
int ofmo_cutoff_dd_(
	// input arguments
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[], const int shel_add[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	// output arguments
	int leading_cs_pair[],
	double csp_schwarz[], int csp_ics[], int csp_jcs[],
	int csp_leading_ps_pair[],
	double psp_zeta[], double psp_dkps[], double psp_xiza[] ) {
    int ics, ics0, ics1, ips0, ips1, iat;
    int jcs, jcs0,       jps0, jps1, jat;
    int npps, ipps, ncs_pair, nps_pair;
    int i, Lab;
    double BA[3], A[3], B[3], AB2;
    double max_eri;
    double vzeta[MAXPSPAIR], vdkps[MAXPSPAIR], vxiza[MAXPSPAIR];
    int La=*pLa, Lb=*pLb;

    Lab       = La*(La+1)/2 + Lb;
    ncs_pair  = leading_cs_pair[Lab];
    nps_pair  = csp_leading_ps_pair[ncs_pair];

    ics0 = leading_cs[La];	// La=2
    ics1 = leading_cs[La+1];
    jcs0 = leading_cs[Lb];	// Lb=2
    for ( ics=ics0; ics<ics1; ics++ ) {
	ips0 = shel_add[ics];
	ips1 = ips0 + shel_tem[ics];
	iat  = shel_atm[ics];
	A[0]=atom_x[ iat ]; A[1]=atom_y[ iat ]; A[2]=atom_z[ iat ];
	for ( jcs=jcs0; jcs<=ics; jcs++ ) {
	    jps0 = shel_add[jcs];
	    jps1 = jps0 + shel_tem[jcs];
	    jat  = shel_atm[jcs];
	    B[0]=atom_x[ jat ]; B[1]=atom_y[ jat ]; B[2]=atom_z[ jat ];

	    AB2  = ZERO;
	    for ( i=0; i<3; i++ ) {
		BA[i] = B[i] - A[i];
		AB2  += BA[i]*BA[i];
	    }
	    
	    npps = schwarz_calc_ps_pair_params(
		    prim_exp, prim_coe, ips0, ips1, jps0, jps1, AB2,
		    vzeta, vdkps, vxiza );
	    if ( npps == 0 ) continue;
	    
	    max_eri = schwarz_core_dddd_(
		    npps, vzeta, vdkps, vxiza, BA, AB2 );

	    if ( max_eri > EPS_CS_PAIR2 ) {
		csp_schwarz[ncs_pair] = sqrt( max_eri );
		csp_ics[ncs_pair] = ics;
		csp_jcs[ncs_pair] = jcs;
		for ( ipps=0; ipps<npps; ipps++ ) {
		    psp_zeta[nps_pair] = vzeta[ipps];
		    psp_dkps[nps_pair] = vdkps[ipps];
		    psp_xiza[nps_pair] = vxiza[ipps];
		    nps_pair++;
		}
		csp_leading_ps_pair[ncs_pair+1] = nps_pair;
		ncs_pair++;
	    }
	}	// jcs
    }		// ics
    leading_cs_pair[Lab+1] = ncs_pair;
    return 0;
}

// デバッグ用コード
#ifndef MAXSTRLEN
#define MAXSTRLEN 256
#endif
/** カットオフテーブルに関する統計データを出力する
 * */
int ofmo_cutoff_show_table( const int maxlqn,
	const int leading_cs[], const int shel_tem[],
	const int leading_cs_pair[],
	const int csp_leading_ps_pair[] ) {
    int Lab, La, Lb;
    int ics, jcs;
    int ncs1, ncs2, ics1, ics2, jcs1, jcs2, ijcs1, ijcs2;
    int ncspair_t, ncspair_s, npspair_t, npspair_s;
    double r_ncspair, r_npspair;
    char cabcd[MAXSTRLEN];
    char cstype[] = {'S', 'P', 'D', 'F', 'G', 'H', 'I', 'J'};
    printf("******************** cutoff results **********************\n");
    printf("%6s: %21s   %26s\n", " ", "NCS PAIR         ",
	    "NPS PAIR           ");
    printf("%6s: %6s/%6s (%6s)   %8s/%8s (%6s)\n",
	    "cstype", "surved", "total", "rate",
	    "surved", "total", "rate");
    for ( La=0; La<=maxlqn; La++ ) {
	ics1 = leading_cs[La];
	ics2 = leading_cs[La+1];
	ncs1 = ics2 - ics1;
	for ( Lb=0; Lb<=La; Lb++ ) {
	    jcs1 = leading_cs[Lb];
	    jcs2 = leading_cs[Lb+1];
	    ncs2 = jcs2 - jcs1;
	    ncspair_t = (La==Lb ? (ncs1*(ncs1+1)/2) : (ncs1*ncs2) );
	    Lab = La*(La+1)/2 + Lb;
	    ijcs1 = leading_cs_pair[Lab];
	    ijcs2 = leading_cs_pair[Lab+1];
	    ncspair_s = ijcs2 - ijcs1;
	    // ここで、このCSペアタイプに属するPSペア数を数える
	    // （この後に及んで・・・）
	    npspair_t = 0;
	    if (La==Lb) {
		for ( ics=ics1; ics<ics2; ics++ ) {
		    for ( jcs=ics1; jcs<=ics; jcs++) {
			npspair_t += shel_tem[ics] * shel_tem[jcs];
		    }
		}
	    } else {
		for ( ics=ics1; ics<ics2; ics++ ) {
		    for ( jcs=jcs1; jcs<jcs2; jcs++) {
			npspair_t += shel_tem[ics] * shel_tem[jcs];
		    }
		}
	    }

	    npspair_s = csp_leading_ps_pair[ijcs2]
		- csp_leading_ps_pair[ijcs1];
	    r_ncspair = (double)ncspair_s / (double)ncspair_t * 100.0;
	    r_npspair = (double)npspair_s / (double)npspair_t * 100.0;
	    sprintf(cabcd, "(%c|%c)", cstype[La], cstype[Lb]);
	    printf("%6s: %6d/%6d (%6.2f)   %8d/%8d (%6.2f)\n", cabcd,
		    ncspair_s, ncspair_t, r_ncspair,
		    npspair_s, npspair_t, r_npspair);
	}
    }
    printf("**********************************************************\n");
    return 0;
}

/**
 * sort cutoff table
 * */
#ifdef SORT_CSP

static int ofmo_cutoff_cmpnum_(const void *a, const void *b)
{
  int ret = 0;
  double t;
  struct ofmo_cutoff_sort_t *a0 = (struct ofmo_cutoff_sort_t *)a;
  struct ofmo_cutoff_sort_t *b0 = (struct ofmo_cutoff_sort_t *)b;
  ret = (a0->iv - b0->iv);
  /*
  if (ret==0) {
    t = (a0->vab - b0->vab);
    if (t>0) ret = 1;
    else if (t<0) ret = -1;
  }
  */
#ifdef SORT_CSP_SCHWARZ
  t = (a0->vab - b0->vab);
  if (t>0) ret = 1;
  else if (t<0) ret = -1;
//  if (fabs(t)<.1) ret = (a0->iv - b0->iv);
#endif
#ifdef SORT_CSP_REV
  ret = -ret;
#endif
  return ret;
}

int ofmo_cutoff_sort_(
        const int La, const int Lb, const int leading_cs[],
        const int shel_tem[], const int shel_atm[], const int shel_add[],
        const double atom_x[], const double atom_y[],
        const double atom_z[],
        const double prim_exp[], const double prim_coe[],
        int leading_cs_pair[],
        double csp_schwarz[], int csp_ics[], int csp_jcs[],
        int csp_leading_ps_pair[],
        double psp_zeta[], double psp_dkps[], double psp_xiza[] )
{
    int Lab, csp0, csp1, ncsp, psp0, psp1, npsp, ipsp;
    struct ofmo_cutoff_sort_t *icsp;
    int *ticsp;
    double *tdcsp, *tdpsp;
    char cstype[] = {'S', 'P', 'D', 'F', 'G', 'H', 'I', 'J'};
    char sLab[]="SS";

    sLab[0]=cstype[La];
    sLab[1]=cstype[Lb];

    Lab  = La*(La+1)/2 + Lb;
    csp0  = leading_cs_pair[Lab];
    csp1  = leading_cs_pair[Lab+1];
    ncsp = csp1 - csp0;
    psp0 = csp_leading_ps_pair[csp0];
    psp1 = csp_leading_ps_pair[csp1];
    npsp = psp1 - psp0;

    icsp = (struct ofmo_cutoff_sort_t *)malloc((ncsp+1) * sizeof(struct ofmo_cutoff_sort_t));
    if (icsp==NULL) exit(1);
    tdcsp = (double *)malloc((ncsp+npsp)*sizeof(double));
    if (tdcsp==NULL) exit(1);
    ticsp = (int *)tdcsp;
    tdpsp = tdcsp;

    for (int i=0; i<ncsp; i++) {
      icsp[i].ic = csp0 + i;
      icsp[i].iv = csp_leading_ps_pair[csp0+i+1] - csp_leading_ps_pair[csp0+i];
      icsp[i].vab = fabs(csp_schwarz[csp0+i]);
//      printf("icsp0:%4d %4d %4d %g\n", i, icsp[i].ic, icsp[i].iv, icsp[i].vab);
    }
    qsort(icsp, ncsp, sizeof(struct ofmo_cutoff_sort_t), ofmo_cutoff_cmpnum_);
//    for (int i=0; i<ncsp; i++) printf("icsp1:%4d %4d %4d %g\n", i, icsp[i].ic, icsp[i].iv, icsp[i].vab);
#ifdef DEBUG_SORT
    {
      int tiv=1, niv=0;
      printf("--- ncsp[%2s] ---\n", sLab);
      for (int i=0; i<ncsp; i++) {
        if (tiv == icsp[i].iv) {
          niv++;
        } else {
          printf("%1d: %d %d\n",Lab, tiv, niv);
          tiv = icsp[i].iv;
          niv = 1;
        }
      }
      printf("%1d: %d %d\n",Lab, tiv, niv);
    }
#endif

    for (int i=0; i<ncsp; i++) tdcsp[i] = csp_schwarz[icsp[i].ic];
    memcpy(&csp_schwarz[csp0], tdcsp, ncsp*sizeof(double));

    for (int i=0; i<ncsp; i++) ticsp[i] = csp_ics[icsp[i].ic];
    memcpy(&csp_ics[csp0], ticsp, ncsp*sizeof(int));

    for (int i=0; i<ncsp; i++) ticsp[i] = csp_jcs[icsp[i].ic];
    memcpy(&csp_jcs[csp0], ticsp, ncsp*sizeof(int));

    ipsp = 0;
    for (int i=0; i<ncsp; i++) {
      memcpy(&tdpsp[ipsp], &psp_zeta[csp_leading_ps_pair[icsp[i].ic]],
          icsp[i].iv * sizeof(double));
      ipsp += icsp[i].iv;
    }
    memcpy(&psp_zeta[psp0], tdpsp, npsp*sizeof(double));

    ipsp = 0;
    for (int i=0; i<ncsp; i++) {
      memcpy(&tdpsp[ipsp], &psp_dkps[csp_leading_ps_pair[icsp[i].ic]],
          icsp[i].iv * sizeof(double));
      ipsp += icsp[i].iv;
    }
    memcpy(&psp_dkps[psp0], tdpsp, npsp*sizeof(double));

    ipsp = 0;
    for (int i=0; i<ncsp; i++) {
      memcpy(&tdpsp[ipsp], &psp_xiza[csp_leading_ps_pair[icsp[i].ic]],
          icsp[i].iv * sizeof(double));
      ipsp += icsp[i].iv;
    }
    memcpy(&psp_xiza[psp0], tdpsp, npsp*sizeof(double));

    ipsp = 0;
    for (int i=0; i<ncsp; i++) {
      csp_leading_ps_pair[csp0 + i] = psp0 + ipsp;
      ipsp += icsp[i].iv;
    }

    free(tdcsp);
    free(icsp);
    return 0;
}
#endif /* SORT_CSP */
