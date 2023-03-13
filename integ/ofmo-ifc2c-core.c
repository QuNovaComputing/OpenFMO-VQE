/**
 * @file ofmo-ifc2c-core.c
 *
 * １つのCSペアに対する２中心クーロン積分を計算する関数 。
 * */
/**
 * @defgroup core-ifc2c ２中心クーロン積分を行う関数群
 *
 * ２中心クーロン相互作用項の計算を行う下位関数。
 * 2011/06/16現在は、(s,s)から(d,d)までの積分を行える
 * ようになっている。
 *
 * 各種積分を行う関数 \c oneint_ifc2c_xx_ はすべて同じ
 * 引数をとる。その引数の説明を記す。
 *
 * @param[in] *ips0 １つ目のCSに属するPSの先頭PS番号
 * @param[in] *nps_i １つ目のCSの縮約長（PS数）
 * @param[in] A[3] １つ目のCSの中心座標(au 単位)
 * @param[in] *ips0 ２つ目のCSに属するPSの先頭PS番号
 * @param[in] *nps_j ２つ目のCSの縮約長（PS数）
 * @param[in] B[3] ２つ目のCSの中心座標(au 単位)
 * @param[in] prim_exp[ips] PS番号 \c ips 番目のPSの軌道指数
 * @param[in] prim_coe[ips] PS番号 \c ips 番目のPSの規格化定数込みの
 *     縮約係数
 * @param[in] *nat 原子数
 * @param[in] atom_x[iat] 原子の番号 \c iat の原子のx座標
 * @param[in] atom_y[iat] 原子の番号 \c iat の原子のy座標
 * @param[in] atom_z[iat] 原子の番号 \c iat の原子のz座標
 * @param[in] atm_pop[iat] 原子の番号 \c iat の原子のatomic population
 *
 * @param[out] XX[] 計算した２中心クーロン積分を格納する配列。
 *     \c XX = \c SS , \c PS , \c PP , \c DS , \c DP , \c DD
 *
 * @ingroup integ-core
 * */
#include <stdio.h>
#include <math.h>
#include "fmt.h"
#include "ofmo-def.h"

#define ZERO 0.e0
#define HALF .5e0
#define ONE  1.e0

static double _2PI_;
static double _SQR3_;

/** 初期化関数
 * ２中心クーロン積分を計算する際に用いる定数を設定する
 * @ingroup core-ifc2c
 * */
int ofmo_ifc2c_init() {
    static int called = false;
    if ( called  ) return 0;
    _2PI_ = 8.e0 * atan(1.e0);
    _SQR3_ = sqrt(3.e0);
    called = true;
    return 0;
}
/** SSタイプの１つの縮約２中心クーロン積分を計算する関数
 * @ingroup core-ifc2c
 * */
int ifc2c_core_ss__(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const double atm_pop[],
	double SS[1] ) {
    int i;
    int ips, jps, kat, ips1, jps1;
    double BA[3], AB2, P[3], PC[3], PC2, pi2;
    double zeta_a, coef_a, zeta_b, coef_b;
    double zeta, coef, xiza, xi, zi, Qk, css, U;
    double ss[1], sum;

    pi2 = _2PI_;

    for ( i=0, AB2=ZERO; i<3; i++ ) {
	BA[i] = B[i] - A[i];
	AB2  += BA[i]*BA[i];
    }
    ips1 = (*ips0) + (*nps_i);
    jps1 = (*jps0) + (*nps_j);
    SS[0] = ZERO;
    for ( ips=(*ips0); ips<ips1; ips++ ) {
	zeta_a = prim_exp[ips];
	coef_a = prim_coe[ips];
	for ( jps=(*jps0); jps<jps1; jps++ ) {
	    zeta_b = prim_exp[jps];
	    coef_b = prim_coe[jps];

	    zeta   = zeta_a + zeta_b;
	    zi     = 1.e0   / zeta ;
	    xiza   = zeta_b * zi;
	    coef   = coef_a * coef_b;
	    xi     = zeta_a * xiza;
	    for ( i=0; i<3; i++ ) P[i]  = A[i] + xiza * BA[i];
	    css = pi2 * coef * zi * exp( -xi*AB2 );
	    sum = ZERO;
	    for ( kat=0; kat<(*nat); kat++ ) {
		PC[0] = P[0] - atom_x[kat];
		PC[1] = P[1] - atom_y[kat];
		PC[2] = P[2] - atom_z[kat];
		PC2 = PC[0]*PC[0] + PC[1]*PC[1] + PC[2]*PC[2];
		Qk = -atm_pop[kat];
		U = zeta * PC2;
		fmt( ss, 0, U, css );
		sum += Qk*ss[0];
	    }	// kat
	    SS[0] += sum;
	}	// jps
    }	// ips
    return 0;
}

/** PSタイプの１つの縮約２中心クーロン積分を計算する関数
 * @ingroup core-ifc2c
 * */
int ifc2c_core_ps__(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const double atm_pop[],
	double PS[3] ) {
    int i;
    int ips, jps, kat, ips1, jps1;
    double BA[3], AB2, P[3], PC[3], PA[3], PC2, pi2;
    double zeta_a, coef_a, zeta_b, coef_b;
    double zeta, coef, xiza, xi, zi, Qk, css, U;
    double ss[1+1];

    pi2 = _2PI_;
    for ( i=0, AB2=ZERO; i<3; i++ ) {
	BA[i] = B[i] - A[i];
	AB2  += BA[i]*BA[i];
    }
    ips1 = (*ips0) + (*nps_i);
    jps1 = (*jps0) + (*nps_j);
    for ( i=0; i<3; i++ ) PS[i] = ZERO;
    for ( ips=(*ips0); ips<ips1; ips++ ) {
	zeta_a = prim_exp[ips];
	coef_a = prim_coe[ips];
	for ( jps=(*jps0); jps<jps1; jps++ ) {
	    zeta_b = prim_exp[jps];
	    coef_b = prim_coe[jps];

	    zeta   = zeta_a + zeta_b;
	    zi     = 1.e0   / zeta ;
	    xiza   = zeta_b * zi;
	    coef   = coef_a * coef_b;
	    xi     = zeta_a * xiza;
	    for ( i=0; i<3; i++ ) {
		P[i]  = A[i] + xiza * BA[i];
		PA[i] = xiza * BA[i];
	    }
	    css = pi2 * coef * zi * exp( -xi*AB2 );
	    for ( kat=0; kat<(*nat); kat++ ) {
		PC[0] = P[0] - atom_x[kat];
		PC[1] = P[1] - atom_y[kat];
		PC[2] = P[2] - atom_z[kat];
		PC2 = PC[0]*PC[0] + PC[1]*PC[1] + PC[2]*PC[2];
		Qk = -atm_pop[kat];
		U = zeta * PC2;
		fmt( ss, 1, U, Qk*css );
		for ( i=0; i<3; i++ ) PS[i] += PA[i]*ss[0] - PC[i]*ss[1];
	    }	// kat
	}	// jps
    }	// ips
    return 0;
}

/** PPタイプの１つの縮約２中心クーロン積分を計算する関数
 * @ingroup core-ifc2c
 * */
int ifc2c_core_pp__(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const double atm_pop[],
	double PP[3*3] ) {
    int i, m, m1;
    int ips, jps, kat, ips1, jps1;
    double BA[3], AB2, P[3], PC[3], PA[3], PC2, pi2;
    double zeta_a, coef_a, zeta_b, coef_b;
    double zeta, coef, xiza, xi, zi, Qk, css, U, zeta2;
    double tmpss;
    double ss[2+1], ps[1+1][3];
    double PS[3], DS[6];

    pi2 = _2PI_;
    for ( i=0, AB2=ZERO; i<3; i++ ) {
	BA[i] = B[i] - A[i];
	AB2  += BA[i]*BA[i];
    }
    ips1 = (*ips0) + (*nps_i);
    jps1 = (*jps0) + (*nps_j);
    for ( i=0; i<3; i++ ) PS[i] = ZERO;
    for ( i=0; i<6; i++ ) DS[i] = ZERO;
    for ( ips=(*ips0); ips<ips1; ips++ ) {
	zeta_a = prim_exp[ips];
	coef_a = prim_coe[ips];
	for ( jps=(*jps0); jps<jps1; jps++ ) {
	    zeta_b = prim_exp[jps];
	    coef_b = prim_coe[jps];

	    zeta   = zeta_a + zeta_b;
	    zi     = 1.e0   / zeta ;
	    xiza   = zeta_b * zi;
	    coef   = coef_a * coef_b;
	    xi     = zeta_a * xiza;
	    zeta2  = HALF   * zi;
	    for ( i=0; i<3; i++ ) {
		P[i]  = A[i] + xiza * BA[i];
		PA[i] = xiza * BA[i];
	    }
	    css = pi2 *coef * zi * exp( -xi*AB2 );
	    for ( kat=0; kat<(*nat); kat++ ) {
		PC[0] = P[0] - atom_x[kat];
		PC[1] = P[1] - atom_y[kat];
		PC[2] = P[2] - atom_z[kat];
		PC2 = PC[0]*PC[0] + PC[1]*PC[1] + PC[2]*PC[2];
		Qk = -atm_pop[kat];
		U = zeta * PC2;
		fmt( ss, 2, U, Qk*css );
		// (p|A|s) (m=0,1)
		for ( m=0, m1=1; m<=1; m++, m1++ ) {
		    for ( i=0; i<3; i++ )
			ps[m][i] = PA[i]*ss[m] - PC[i]*ss[m1];
		}
		for ( i=0; i<3; i++ ) PS[i] += ps[0][i];
		// (d|A|s) m=0
		tmpss = ss[0]-ss[1];
		DS[0] +=  PA[0]*ps[0][0]-PC[0]*ps[1][0]+zeta2*tmpss;
		DS[1] +=  PA[1]*ps[0][1]-PC[1]*ps[1][1]+zeta2*tmpss;
		DS[2] +=  PA[2]*ps[0][2]-PC[2]*ps[1][2]+zeta2*tmpss;
		DS[3] +=  PA[1]*ps[0][0]-PC[1]*ps[1][0];
		DS[4] +=  PA[2]*ps[0][1]-PC[2]*ps[1][1];
		DS[5] +=  PA[0]*ps[0][2]-PC[0]*ps[1][2];
	    }	// kat
	}	// jps
    }	// ips
    // HRR
    PP[0*3+0] = DS[0] - BA[0]*PS[0];
    PP[0*3+1] = DS[3] - BA[1]*PS[0];
    PP[0*3+2] = DS[5] - BA[2]*PS[0];
    PP[1*3+0] = DS[3] - BA[0]*PS[1];
    PP[1*3+1] = DS[1] - BA[1]*PS[1];
    PP[1*3+2] = DS[4] - BA[2]*PS[1];
    PP[2*3+0] = DS[5] - BA[0]*PS[2];
    PP[2*3+1] = DS[4] - BA[1]*PS[2];
    PP[2*3+2] = DS[2] - BA[2]*PS[2];
    return 0;
}

/** DSタイプの１つの縮約２中心クーロン積分を計算する関数
 * @ingroup core-ifc2c
 * */
int ifc2c_core_ds__(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const double atm_pop[],
	double DS[6] ) {
    int i, m, m1;
    int ips, jps, kat, ips1, jps1;
    double BA[3], AB2, P[3], PC[3], PA[3], PC2, pi2, sqr3;
    double zeta_a, coef_a, zeta_b, coef_b;
    double zeta, coef, xiza, xi, zi, Qk, css, U, zeta2;
    double tmpss;
    double ss[2+1], ps[1+1][3];

    pi2  = _2PI_;
    sqr3 = _SQR3_;
    for ( i=0, AB2=ZERO; i<3; i++ ) {
	BA[i] = B[i] - A[i];
	AB2  += BA[i]*BA[i];
    }
    ips1 = (*ips0) + (*nps_i);
    jps1 = (*jps0) + (*nps_j);
    for ( i=0; i<6; i++ ) DS[i] = ZERO;
    for ( ips=(*ips0); ips<ips1; ips++ ) {
	zeta_a = prim_exp[ips];
	coef_a = prim_coe[ips];
	for ( jps=(*jps0); jps<jps1; jps++ ) {
	    zeta_b = prim_exp[jps];
	    coef_b = prim_coe[jps];

	    zeta   = zeta_a + zeta_b;
	    zi     = 1.e0   / zeta ;
	    xiza   = zeta_b * zi;
	    coef   = coef_a * coef_b;
	    xi     = zeta_a * xiza;
	    zeta2  = HALF   * zi;
	    for ( i=0; i<3; i++ ) {
		P[i]  = A[i] + xiza * BA[i];
		PA[i] = xiza * BA[i];
	    }
	    css = pi2 * coef * zi * exp( -xi*AB2 );
	    for ( kat=0; kat<(*nat); kat++ ) {
		PC[0] = P[0] - atom_x[kat];
		PC[1] = P[1] - atom_y[kat];
		PC[2] = P[2] - atom_z[kat];
		PC2 = PC[0]*PC[0] + PC[1]*PC[1] + PC[2]*PC[2];
		Qk = -atm_pop[kat];
		U = zeta * PC2;
		fmt( ss, 2, U, Qk*css );
		// (p|A|s) (m=0,1)
		for ( m=0,  m1=1; m<=1; m++, m1++ ) {
		    for ( i=0; i<3; i++ )
			ps[m][i] = PA[i]*ss[m] - PC[i]*ss[m1];
		}
		// (d|A|s) m=0
		tmpss = ss[0]-ss[1];
		DS[0] +=  PA[0]*ps[0][0]-PC[0]*ps[1][0]+zeta2*tmpss;
		DS[1] +=  PA[1]*ps[0][1]-PC[1]*ps[1][1]+zeta2*tmpss;
		DS[2] +=  PA[2]*ps[0][2]-PC[2]*ps[1][2]+zeta2*tmpss;
		DS[3] +=  PA[1]*ps[0][0]-PC[1]*ps[1][0];
		DS[4] +=  PA[0]*ps[0][2]-PC[0]*ps[1][2];
		DS[5] +=  PA[2]*ps[0][1]-PC[2]*ps[1][1];
	    }	// kat
	}	// jps
    }	// ips
    for ( i=3; i<6; i++ ) DS[i] *= sqr3;
    return 0;
}

/** DPタイプの１つの縮約２中心クーロン積分を計算する関数
 * @ingroup core-ifc2c
 * */
int ifc2c_core_dp__(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const double atm_pop[],
	double DP[6*3] ) {
    int i, m, m1;
    int ips, jps, kat, ips1, jps1;
    double BA[3], AB2, P[3], PC[3], PA[3], PC2, pi2, sqr3;
    double zeta_a, coef_a, zeta_b, coef_b;
    double zeta, coef, xiza, xi, zi, Qk, css, U, zeta2;
    double tmpss;
    double ss[3+1], ps[2+1][3], ds[1+1][6];
    double DS[6], FS[10];

    pi2  = _2PI_;
    sqr3 = _SQR3_;
    for ( i=0, AB2=ZERO; i<3; i++ ) {
	BA[i] = B[i] - A[i];
	AB2  += BA[i]*BA[i];
    }
    ips1 = (*ips0) + (*nps_i);
    jps1 = (*jps0) + (*nps_j);
    for ( i=0; i< 6; i++ ) DS[i] = ZERO;
    for ( i=0; i<10; i++ ) FS[i] = ZERO;
    for ( ips=(*ips0); ips<ips1; ips++ ) {
	zeta_a = prim_exp[ips];
	coef_a = prim_coe[ips];
	for ( jps=(*jps0); jps<jps1; jps++ ) {
	    zeta_b = prim_exp[jps];
	    coef_b = prim_coe[jps];

	    zeta   = zeta_a + zeta_b;
	    zi     = 1.e0   / zeta ;
	    xiza   = zeta_b * zi;
	    coef   = coef_a * coef_b;
	    xi     = zeta_a * xiza;
	    zeta2  = HALF   * zi;
	    for ( i=0; i<3; i++ ) {
		P[i]  = A[i] + xiza * BA[i];
		PA[i] = xiza * BA[i];
	    }
	    css = pi2 * coef * zi * exp( -xi*AB2 );
	    for ( kat=0; kat<(*nat); kat++ ) {
		PC[0] = P[0] - atom_x[kat];
		PC[1] = P[1] - atom_y[kat];
		PC[2] = P[2] - atom_z[kat];
		PC2 = PC[0]*PC[0] + PC[1]*PC[1] + PC[2]*PC[2];
		Qk = -atm_pop[kat];
		U = zeta * PC2;
		fmt( ss, 3, U, Qk*css );
		// (p|A|s) (m=0,2)
		for ( m=0, m1=1; m<=2; m++, m1++ ) {
		    for ( i=0; i<3; i++ )
			ps[m][i] = PA[i]*ss[m] - PC[i]*ss[m1];
		}
		// (d|A|s) (m=0,1)
		for ( m=0, m1=1; m<=1; m++, m1++ ) {
		    tmpss = ss[m]-ss[m1];
		    ds[m][0] = PA[0]*ps[m][0]-PC[0]*ps[m1][0]+zeta2*tmpss;
		    ds[m][1] = PA[1]*ps[m][1]-PC[1]*ps[m1][1]+zeta2*tmpss;
		    ds[m][2] = PA[2]*ps[m][2]-PC[2]*ps[m1][2]+zeta2*tmpss;
		    ds[m][3] = PA[1]*ps[m][0]-PC[1]*ps[m1][0];
		    ds[m][4] = PA[2]*ps[m][1]-PC[2]*ps[m1][1];
		    ds[m][5] = PA[0]*ps[m][2]-PC[0]*ps[m1][2];
		}
		for ( i=0; i<6; i++ ) DS[i] += ds[0][i];
		// (f|A|s) m=0
		FS[0]+=PA[0]*ds[0][0]-PC[0]*ds[1][0]
		    +zi*(ps[0][0]-ps[1][0]);
		FS[1]+=PA[1]*ds[0][1]-PC[1]*ds[1][1]
		    +zi*(ps[0][1]-ps[1][1]);
		FS[2]+=PA[2]*ds[0][2]-PC[2]*ds[1][2]
		    +zi*(ps[0][2]-ps[1][2]);
		FS[3]+=PA[1]*ds[0][0]-PC[1]*ds[1][0];
		FS[4]+=PA[2]*ds[0][1]-PC[2]*ds[1][1];
		FS[5]+=PA[0]*ds[0][2]-PC[0]*ds[1][2];
		FS[6]+=PA[0]*ds[0][1]-PC[0]*ds[1][1];
		FS[7]+=PA[1]*ds[0][2]-PC[1]*ds[1][2];
		FS[8]+=PA[2]*ds[0][0]-PC[2]*ds[1][0];
		FS[9]+=PA[2]*ds[0][3]-PC[2]*ds[1][3];
	    }	// kat
	}	// jps
    }	// ips
    // HRR
    DP[0*3+0] = FS[0] - BA[0]*DS[0];
    DP[0*3+1] = FS[3] - BA[1]*DS[0];
    DP[0*3+2] = FS[8] - BA[2]*DS[0];

    DP[1*3+0] = FS[6] - BA[0]*DS[1];
    DP[1*3+1] = FS[1] - BA[1]*DS[1];
    DP[1*3+2] = FS[4] - BA[2]*DS[1];

    DP[2*3+0] = FS[5] - BA[0]*DS[2];
    DP[2*3+1] = FS[7] - BA[1]*DS[2];
    DP[2*3+2] = FS[2] - BA[2]*DS[2];

    DP[3*3+0] = FS[3] - BA[0]*DS[3];
    DP[3*3+1] = FS[6] - BA[1]*DS[3];
    DP[3*3+2] = FS[9] - BA[2]*DS[3];

    DP[4*3+0] = FS[8] - BA[0]*DS[5];
    DP[4*3+1] = FS[9] - BA[1]*DS[5];
    DP[4*3+2] = FS[5] - BA[2]*DS[5];

    DP[5*3+0] = FS[9] - BA[0]*DS[4];
    DP[5*3+1] = FS[4] - BA[1]*DS[4];
    DP[5*3+2] = FS[7] - BA[2]*DS[4];

    for ( i=3*3; i<6*3; i++ ) DP[i] *= sqr3;
    return 0;
}

/** DDタイプの１つの縮約２中心クーロン積分を計算する関数
 * @ingroup core-ifc2c
 * */
int ifc2c_core_dd__(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const double atm_pop[],
	double DD[6*6] ) {
    int i, j, ix, m, m1;
    int ips, jps, kat, ips1, jps1;
    double BA[3], AB2, P[3], PC[3], PA[3], PC2, pi2, sqr3;
    double zeta_a, coef_a, zeta_b, coef_b;
    double zeta, coef, xiza, xi, zi, Qk, css, U, zeta2, zeta23;
    double tmpss, tmpds[3];
    double ss[4+1], ps[3+1][3], ds[2+1][6], fs[1+1][10];
    double DS[6], FS[10], GS[15];
    double DP[6*3], FP[10*3];
    double coe0, coe;

    pi2  = _2PI_;
    sqr3 = _SQR3_;
    for ( i=0, AB2=ZERO; i<3; i++ ) {
	BA[i] = B[i] - A[i];
	AB2  += BA[i]*BA[i];
    }
    ips1 = (*ips0) + (*nps_i);
    jps1 = (*jps0) + (*nps_j);
    for ( i=0; i< 6; i++ ) DS[i] = ZERO;
    for ( i=0; i<10; i++ ) FS[i] = ZERO;
    for ( i=0; i<15; i++ ) GS[i] = ZERO;
    for ( ips=(*ips0); ips<ips1; ips++ ) {
	zeta_a = prim_exp[ips];
	coef_a = prim_coe[ips];
	for ( jps=(*jps0); jps<jps1; jps++ ) {
	    zeta_b = prim_exp[jps];
	    coef_b = prim_coe[jps];

	    zeta   = zeta_a + zeta_b;
	    zi     = 1.e0   / zeta ;
	    xiza   = zeta_b * zi;
	    coef   = coef_a * coef_b;
	    xi     = zeta_a * xiza;
	    zeta2  = HALF   * zi;
	    zeta23 = 1.5e0  * zi;
	    for ( i=0; i<3; i++ ) {
		P[i]  = A[i] + xiza * BA[i];
		PA[i] = xiza * BA[i];
	    }
	    css = pi2 * coef * zi * exp( -xi*AB2 );
	    for ( kat=0; kat<(*nat); kat++ ) {
		PC[0] = P[0] - atom_x[kat];
		PC[1] = P[1] - atom_y[kat];
		PC[2] = P[2] - atom_z[kat];
		PC2 = PC[0]*PC[0] + PC[1]*PC[1] + PC[2]*PC[2];
		Qk = -atm_pop[kat];
		U = zeta * PC2;
		fmt( ss, 4, U, Qk*css );
		// (p|A|s) (m=0,3)
		for ( m=0, m1=1; m<=3; m++, m1++ ) {
		    for ( i=0; i<3; i++ )
			ps[m][i] = PA[i]*ss[m] - PC[i]*ss[m1];
		}
		// (d|A|s) (m=0,2)
		for ( m=0, m1=1; m<=2; m++, m1++ ) {
		    tmpss = ss[m]-ss[m1];
		    ds[m][0] = PA[0]*ps[m][0]-PC[0]*ps[m1][0]+zeta2*tmpss;
		    ds[m][1] = PA[1]*ps[m][1]-PC[1]*ps[m1][1]+zeta2*tmpss;
		    ds[m][2] = PA[2]*ps[m][2]-PC[2]*ps[m1][2]+zeta2*tmpss;
		    ds[m][3] = PA[1]*ps[m][0]-PC[1]*ps[m1][0];
		    ds[m][4] = PA[2]*ps[m][1]-PC[2]*ps[m1][1];
		    ds[m][5] = PA[0]*ps[m][2]-PC[0]*ps[m1][2];
		}
		for ( i=0; i<6; i++ ) DS[i] += ds[0][i];
		// (f|A|s) (m=0,1)
		for ( m=0, m1=1; m<=1; m++, m1++ ) {
		    fs[m][0]=PA[0]*ds[m][0]-PC[0]*ds[m1][0]
			+zi*(ps[m][0]-ps[m1][0]);
		    fs[m][1]=PA[1]*ds[m][1]-PC[1]*ds[m1][1]
			+zi*(ps[m][1]-ps[m1][1]);
		    fs[m][2]=PA[2]*ds[m][2]-PC[2]*ds[m1][2]
			+zi*(ps[m][2]-ps[m1][2]);
		    fs[m][3]=PA[1]*ds[m][0]-PC[1]*ds[m1][0];
		    fs[m][4]=PA[2]*ds[m][1]-PC[2]*ds[m1][1];
		    fs[m][5]=PA[0]*ds[m][2]-PC[0]*ds[m1][2];
		    fs[m][6]=PA[0]*ds[m][1]-PC[0]*ds[m1][1];
		    fs[m][7]=PA[1]*ds[m][2]-PC[1]*ds[m1][2];
		    fs[m][8]=PA[2]*ds[m][0]-PC[2]*ds[m1][0];
		    fs[m][9]=PA[2]*ds[m][3]-PC[2]*ds[m1][3];
		}
		for ( i=0; i<10; i++ ) FS[i] += fs[0][i];
		// (g|A|s)
		m =0;
		m1=1;
		tmpds[0] = ds[m][0] - ds[m1][0];
		tmpds[1] = ds[m][1] - ds[m1][1];
		tmpds[2] = ds[m][2] - ds[m1][2];
		GS[ 0]+=PA[0]*fs[m][0] - PC[0]*fs[m1][0] + zeta23*tmpds[0];
		GS[ 1]+=PA[1]*fs[m][1] - PC[1]*fs[m1][1] + zeta23*tmpds[1];
		GS[ 2]+=PA[2]*fs[m][2] - PC[2]*fs[m1][2] + zeta23*tmpds[2];
		GS[ 3]+=PA[1]*fs[m][0] - PC[1]*fs[m1][0];
		GS[ 4]+=PA[2]*fs[m][1] - PC[2]*fs[m1][1];
		GS[ 5]+=PA[0]*fs[m][2] - PC[0]*fs[m1][2];
		GS[ 6]+=PA[1]*fs[m][3] - PC[1]*fs[m1][3] + zeta2*tmpds[0];
		GS[ 7]+=PA[2]*fs[m][4] - PC[2]*fs[m1][4] + zeta2*tmpds[1];
		GS[ 8]+=PA[0]*fs[m][5] - PC[0]*fs[m1][5] + zeta2*tmpds[2];
		GS[ 9]+=PA[0]*fs[m][1] - PC[0]*fs[m1][1];
		GS[10]+=PA[1]*fs[m][2] - PC[1]*fs[m1][2];
		GS[11]+=PA[2]*fs[m][0] - PC[2]*fs[m1][0];
		GS[12]+=PA[2]*fs[m][3] - PC[2]*fs[m1][3];
		GS[13]+=PA[0]*fs[m][4] - PC[0]*fs[m1][4];
		GS[14]+=PA[1]*fs[m][5] - PC[1]*fs[m1][5];
	    }	// kat
	}	// jps
    }	// ips
    // HRR
    DP[0*3+0] = FS[0] - BA[0]*DS[0];
    DP[0*3+1] = FS[3] - BA[1]*DS[0];
    DP[0*3+2] = FS[8] - BA[2]*DS[0];
    DP[1*3+0] = FS[6] - BA[0]*DS[1];
    DP[1*3+1] = FS[1] - BA[1]*DS[1];
    DP[1*3+2] = FS[4] - BA[2]*DS[1];
    DP[2*3+0] = FS[5] - BA[0]*DS[2];
    DP[2*3+1] = FS[7] - BA[1]*DS[2];
    DP[2*3+2] = FS[2] - BA[2]*DS[2];
    DP[3*3+0] = FS[3] - BA[0]*DS[3];
    DP[3*3+1] = FS[6] - BA[1]*DS[3];
    DP[3*3+2] = FS[9] - BA[2]*DS[3];
    DP[4*3+0] = FS[9] - BA[0]*DS[4];
    DP[4*3+1] = FS[4] - BA[1]*DS[4];
    DP[4*3+2] = FS[7] - BA[2]*DS[4];
    DP[5*3+0] = FS[8] - BA[0]*DS[5];
    DP[5*3+1] = FS[9] - BA[1]*DS[5];
    DP[5*3+2] = FS[5] - BA[2]*DS[5];
    // (F|V|P) = (G||S) - (B-A)*(F||S)
    FP[0*3+0] = GS[ 0] - BA[0]*FS[0];
    FP[0*3+1] = GS[ 3] - BA[1]*FS[0];
    FP[0*3+2] = GS[11] - BA[2]*FS[0];
    FP[1*3+0] = GS[ 9] - BA[0]*FS[1];
    FP[1*3+1] = GS[ 1] - BA[1]*FS[1];
    FP[1*3+2] = GS[ 4] - BA[2]*FS[1];
    FP[2*3+0] = GS[ 5] - BA[0]*FS[2];
    FP[2*3+1] = GS[10] - BA[1]*FS[2];
    FP[2*3+2] = GS[ 2] - BA[2]*FS[2];
    FP[3*3+0] = GS[ 3] - BA[0]*FS[3];
    FP[3*3+1] = GS[ 6] - BA[1]*FS[3];
    FP[3*3+2] = GS[12] - BA[2]*FS[3];
    FP[4*3+0] = GS[13] - BA[0]*FS[4];
    FP[4*3+1] = GS[ 4] - BA[1]*FS[4];
    FP[4*3+2] = GS[ 7] - BA[2]*FS[4];
    FP[5*3+0] = GS[ 8] - BA[0]*FS[5];
    FP[5*3+1] = GS[14] - BA[1]*FS[5];
    FP[5*3+2] = GS[ 5] - BA[2]*FS[5];
    FP[6*3+0] = GS[ 6] - BA[0]*FS[6];
    FP[6*3+1] = GS[ 9] - BA[1]*FS[6];
    FP[6*3+2] = GS[13] - BA[2]*FS[6];
    FP[7*3+0] = GS[14] - BA[0]*FS[7];
    FP[7*3+1] = GS[ 7] - BA[1]*FS[7];
    FP[7*3+2] = GS[10] - BA[2]*FS[7];
    FP[8*3+0] = GS[11] - BA[0]*FS[8];
    FP[8*3+1] = GS[12] - BA[1]*FS[8];
    FP[8*3+2] = GS[ 8] - BA[2]*FS[8];
    FP[9*3+0] = GS[12] - BA[0]*FS[9];
    FP[9*3+1] = GS[13] - BA[1]*FS[9];
    FP[9*3+2] = GS[14] - BA[2]*FS[9];
    // (D|V|D) = (F||P) - (B-A)*(D||P)
    DD[0*6+0] = FP[0*3+0] - BA[0]*DP[0*3+0];
    DD[0*6+1] = FP[3*3+1] - BA[1]*DP[0*3+1];
    DD[0*6+2] = FP[8*3+2] - BA[2]*DP[0*3+2];
    DD[0*6+3] = FP[0*3+1] - BA[0]*DP[0*3+1];
    DD[0*6+4] = FP[8*3+0] - BA[2]*DP[0*3+0];
    DD[0*6+5] = FP[3*3+2] - BA[1]*DP[0*3+2];

    DD[1*6+0] = FP[6*3+0] - BA[0]*DP[1*3+0];
    DD[1*6+1] = FP[1*3+1] - BA[1]*DP[1*3+1];
    DD[1*6+2] = FP[4*3+2] - BA[2]*DP[1*3+2];
    DD[1*6+3] = FP[6*3+1] - BA[0]*DP[1*3+1];
    DD[1*6+4] = FP[4*3+0] - BA[2]*DP[1*3+0];
    DD[1*6+5] = FP[1*3+2] - BA[1]*DP[1*3+2];

    DD[2*6+0] = FP[5*3+0] - BA[0]*DP[2*3+0];
    DD[2*6+1] = FP[7*3+1] - BA[1]*DP[2*3+1];
    DD[2*6+2] = FP[2*3+2] - BA[2]*DP[2*3+2];
    DD[2*6+3] = FP[5*3+1] - BA[0]*DP[2*3+1];
    DD[2*6+4] = FP[2*3+0] - BA[2]*DP[2*3+0];
    DD[2*6+5] = FP[7*3+2] - BA[1]*DP[2*3+2];

    DD[3*6+0] = FP[3*3+0] - BA[0]*DP[3*3+0];
    DD[3*6+1] = FP[6*3+1] - BA[1]*DP[3*3+1];
    DD[3*6+2] = FP[9*3+2] - BA[2]*DP[3*3+2];
    DD[3*6+3] = FP[3*3+1] - BA[0]*DP[3*3+1];
    DD[3*6+4] = FP[9*3+0] - BA[2]*DP[3*3+0];
    DD[3*6+5] = FP[6*3+2] - BA[1]*DP[3*3+2];

    DD[4*6+0] = FP[8*3+0] - BA[0]*DP[5*3+0];
    DD[4*6+1] = FP[9*3+1] - BA[1]*DP[5*3+1];
    DD[4*6+2] = FP[5*3+2] - BA[2]*DP[5*3+2];
    DD[4*6+3] = FP[8*3+1] - BA[0]*DP[5*3+1];
    DD[4*6+4] = FP[5*3+0] - BA[2]*DP[5*3+0];
    DD[4*6+5] = FP[9*3+2] - BA[1]*DP[5*3+2];

    DD[5*6+0] = FP[9*3+0] - BA[0]*DP[4*3+0];
    DD[5*6+1] = FP[4*3+1] - BA[1]*DP[4*3+1];
    DD[5*6+2] = FP[7*3+2] - BA[2]*DP[4*3+2];
    DD[5*6+3] = FP[9*3+1] - BA[0]*DP[4*3+1];
    DD[5*6+4] = FP[7*3+0] - BA[2]*DP[4*3+0];
    DD[5*6+5] = FP[4*3+2] - BA[1]*DP[4*3+2];

    //
    for ( i=0, ix=0; i<6; i++ ) {
	coe0 = ( i<3? ONE : sqr3 );
	for ( j=0; j<6; j++, ix++ ) {
	    coe = coe0 * (j<3? ONE : sqr3 );
	    DD[ix] *= coe;
	}
    }
    return 0;
}
