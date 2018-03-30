/**
 * @file ofmo-oneint-core.c
 * １つのCSペアに対する１電子積分を計算する関数群
 * */
/**
 * @defgroup core-oneint 1電子積分を行う関数群
 *
 * １つのCSペアに対する１電子積分を計算する関数群。
 * 2011/06/16現在は、(s,s)から(d,d)までの積分を行える
 * ようになっている。
 * この関数は粒度が小さいので、CとFortranで共通のAPIにするため、
 * すべての引数が参照渡しになっている。
 *
 * 各種積分を行う関数 \c oneint_core_xx_ はすべて同じ
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
 * @param[in] atomic_number[iat] 原子の番号 \c iat の原子の原子番号
 *
 * @param[out] OVI[] 計算した重なり積分を代入する配列
 * @param[out] HCORE[] 計算した１電子ハミルトン行列要素
 *     （運動エネルギー積分＋核-引力積分）を代入する配列
 *
 * @ingroup integ-core
 * */
#include <stdio.h>
#include <math.h>
#include "fmt.h"
#include "ofmo-def.h"

#define ZERO 0.e0
#define HALF .5e0

#define EPS_PS_PAIR 1.e-32

static double EPS_PS_PAIR_NAI;

static double _PI_32_;
static double _2PI_;
static double _inv2_;
static double _inv3_;
static double _spi2_;

/** １電子積分関数の初期化を行う関数
 *
 * １電子積分計算で必要な定数などの初期化を行う
 * @ingroup core-oneint
 * */
int ofmo_oneint_init( ) {
    double pi;
    static int called = false;
    if ( called ) return 0;
    pi = 4.e0 * atan( 1.e0 );
    _PI_32_ = pi * sqrt(pi);
    _2PI_   = 2.e0 * pi;
    _inv2_ = 1.e0/2.e0;
    _inv3_ = 1.e0/3.e0;
    _spi2_ = sqrt( 0.5e0*pi );
    EPS_PS_PAIR_NAI = EPS_PS_PAIR;
    called = true;
    return 0;
}

int ofmo_oneint_set_sum_atomic_numbers( const int sum_atomic_number ) {
    if ( sum_atomic_number < 1 ) {
	EPS_PS_PAIR_NAI = EPS_PS_PAIR;
    } else {
	EPS_PS_PAIR_NAI = EPS_PS_PAIR / (double)sum_atomic_number;
    }
    return 0;
}

extern double *FMT_fmt_table0;
extern double *FMT_fmt_table1;
extern double *FMT_fmt_table2;
extern double *FMT_fmt_table3;
extern double *FMT_fmt_table4;

extern double FMT_fmt_step_size;
extern double FMT_fmt_inv_step_size;
extern double FMT_pi_div2;

/** １つのCS対のSSタイプの１電子積分を計算する
 * @ingroup core-oneint
 * */
void oneint_core_ss__(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const int atomic_number[],
	double OVI[], double HCORE[] ) {
    int i, ips, jps, kat, ips1, jps1;
    double zeta_a, zeta_b, zeta, coef_a, coef_b, coef;
    double sqrzi, zi, xiza, xi, xiab2, exp_ab, css;
    double oviss, keiss, naiss;
    double sum, PC2, dq, U;
    double BA[3], P[3], AB2;
    double _pi32_ = _PI_32_, _2pi_ = _2PI_;
    OVI[0]   = ZERO;
    HCORE[0] = ZERO;
    AB2 = ZERO;
    for ( i=0; i<3; i++ ) {
	BA[i] = B[i]-A[i];
	AB2 += BA[i]*BA[i];
    }
    ips1 = *ips0 + (*nps_i);
    jps1 = *jps0 + (*nps_j);
    /*// debug
//#pragma omp master
    {
	printf("FMT_fmt_table0 = %p\n", FMT_fmt_table0 );
	printf("ips1, jps1 = %d, %d\n", ips1, jps1 );
	fflush(stdout);
    }*/

    for ( ips=(*ips0); ips<ips1; ips++ ) {
	zeta_a = prim_exp[ips];
	coef_a = prim_coe[ips];
	for ( jps=(*jps0); jps<jps1; jps++ ) {
	    zeta_b = prim_exp[jps];
	    coef_b = prim_coe[jps];

	    zeta   = zeta_a + zeta_b;
	    sqrzi  = sqrt( 1.e0 / zeta );
	    zi     = sqrzi  * sqrzi;
	    coef   = coef_a * coef_b;

	    xiza  = zeta_b * zi;
	    xi    = xiza * zeta_a;
	    xiab2 = xi * AB2;

	    exp_ab = coef * exp( -xiab2 );
	    oviss = _pi32_ * zi * sqrzi * exp_ab;
	    keiss = xi * ( 3.e0 - 2.e0 * xiab2 ) * oviss;
	    OVI[0]   += oviss;
	    HCORE[0] += keiss;
	    css   = _2pi_ * zi * exp_ab;
	    if ( fabs(css) < EPS_PS_PAIR_NAI ) continue;
	    for ( i=0; i<3; i++ ) P[i] = xiza*BA[i] + A[i];

	    sum = ZERO;
	    for ( kat=0; kat<(*nat); kat++ ) {
		dq   = (double)atomic_number[kat];
		PC2  = (P[0]-atom_x[kat])*(P[0]-atom_x[kat]);
		PC2 += (P[1]-atom_y[kat])*(P[1]-atom_y[kat]);
		PC2 += (P[2]-atom_z[kat])*(P[2]-atom_z[kat]);
		U = zeta * PC2;
		{
		    int it0, pos;
		    double dT;
		    if ( U < 36e0 ) {
			it0 = (int)(0.5e0 + U * FMT_fmt_inv_step_size);
			dT  = it0 * FMT_fmt_step_size - U;
			pos = it0 * 4;
			naiss = css * ((( FMT_fmt_table0[pos+3]   * dT
					+ FMT_fmt_table0[pos+2] ) * dT
					+ FMT_fmt_table0[pos+1] ) * dT
					+ FMT_fmt_table0[pos+0] );
		    } else {
			naiss = css * _spi2_ * sqrt( 1.e0/(U+U) );
		    }
		}
		//fmt( &naiss, 0, U, css );
		sum += (dq*naiss);
	    }
	    HCORE[0] -= sum;
	}
    }
}

/** １つのCS対のPSタイプの１電子積分を計算する
 * @ingroup core-oneint
 * */
void oneint_core_ps__(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const int atomic_number[],
	double OVI[], double HCORE[] ) {
    int i, ips, jps, kat, ips1, jps1;
    double zeta_a, zeta_b, zeta, coef_a, coef_b, coef;
    double sqrzi, zi, xiza, xi, xiab2, xi2, exp_ab, css, cnass;
    double oviss, keiss;
    double PC2, dq, U;
    double BA[3], P[3], PA[3], PC[3], AB2;
    double keips[3], ovips[3], naips[3];
    double naiss[1+1];
    double _pi32_ = _PI_32_, _2pi_ = _2PI_;

    for ( i=0; i<3; i++ ) OVI[i] = HCORE[i] = ZERO;
    AB2 = ZERO;
    for ( i=0; i<3; i++ ) {
	BA[i] = B[i]-A[i];
	AB2 += BA[i]*BA[i];
    }
    ips1 = *ips0 + (*nps_i);
    jps1 = *jps0 + (*nps_j);

    for ( ips=(*ips0); ips<ips1; ips++ ) {
	zeta_a = prim_exp[ips];
	coef_a = prim_coe[ips];
	for ( jps=(*jps0); jps<jps1; jps++ ) {
	    zeta_b = prim_exp[jps];
	    coef_b = prim_coe[jps];

	    zeta   = zeta_a + zeta_b;
	    sqrzi  = sqrt( 1.e0 / zeta );
	    zi     = sqrzi  * sqrzi;
	    coef   = coef_a * coef_b;

	    xiza  = zeta_b * zi;
	    xi    = xiza * zeta_a;
	    xiab2 = xi * AB2;
	    xi2   = 2.e0 * xi;

	    exp_ab = coef * exp( -xiab2 );
	    oviss = _pi32_ * zi * sqrzi * exp_ab;
	    keiss = xi * ( 3.e0 - 2.e0 * xiab2 ) * oviss;
	    css   = _2pi_ * zi * exp_ab;

	    if ( fabs(css) < EPS_PS_PAIR_NAI ) continue;

	    for ( i=0; i<3; i++ ) {
		P[i]  = xiza*BA[i] + A[i];
		PA[i] = xiza*BA[i];
	    }
	    // overlap integral
	    for ( i=0; i<3; i++ ) ovips[i] = PA[i]*oviss;
	    // kinetic energy integral
	    for ( i=0; i<3; i++ ) keips[i] = PA[i]*keiss + xi2*ovips[i];

	    for ( i=0; i<3; i++ ) naips[i] = ZERO;

	    for ( kat=0; kat<(*nat); kat++ ) {
		dq   = (double)atomic_number[kat];
		PC[0] = P[0]-atom_x[kat];
		PC[1] = P[1]-atom_y[kat];
		PC[2] = P[2]-atom_z[kat];
		PC2 = PC[0]*PC[0] + PC[1]*PC[1] + PC[2]*PC[2];
		cnass = dq * css;
		U = zeta * PC2;
		{
		    int it0, pos;
		    double dT, t_inv, st_inv, dT2, dT3;
		    if ( U < 38e0 ) {
			it0 = (int)(0.5e0 + U * FMT_fmt_inv_step_size);
			dT  = it0 * FMT_fmt_step_size - U;
			dT2 = dT * _inv2_;
			dT3 = dT * _inv3_;
			pos = it0 * (1+4);
			naiss[0] = cnass*(((FMT_fmt_table1[pos+3] * dT3
					+ FMT_fmt_table1[pos+2] ) * dT2
					+ FMT_fmt_table1[pos+1] ) * dT
					+ FMT_fmt_table1[pos+0] );
			naiss[1] = cnass*(((FMT_fmt_table1[pos+4] * dT3
					+ FMT_fmt_table1[pos+3] ) * dT2
					+ FMT_fmt_table1[pos+2] ) * dT
					+ FMT_fmt_table1[pos+1] );
		    } else {
			st_inv = sqrt( 0.5e0 / U );
			t_inv  = st_inv * st_inv;
			naiss[0] = cnass * _spi2_ * st_inv;
			naiss[1] = t_inv * naiss[0];
		    }
		}
		//fmt( naiss, 1, U, cnass );
		for ( i=0; i<3; i++ )
		    naips[i] += PA[i]*naiss[0]-PC[i]*naiss[1];
	    }
	    // contraction
	    for ( i=0; i<3; i++ ) {
		OVI[i] += ovips[i];
		HCORE[i] += (keips[i] - naips[i]);
	    }
	}
    }
}

/** １つのCS対のPPタイプの１電子積分を計算する
 * @ingroup core-oneint
 * */
void oneint_core_pp__(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const int atomic_number[],
	double OVI[], double HCORE[] ) {
    int i, m, ips, jps, kat, ips1, jps1;
    double zeta_a, zeta_b, zeta, coef_a, coef_b, coef;
    double sqrzi, zi, xiza, xi, xiab2, xi2, exp_ab, css, zeta2;
    double oviss, keiss, tmp;
    double BA[3], P[3], PA[3], PC[3], PB[3], AB2;
    double keips[3], keipp[3*3], ovips[3], ovipp[3*3];
    double PC2, dq, U;
    double naiss[2+1], naips[2][3], naipp[3*3];
    double _pi32_ = _PI_32_, _2pi_ = _2PI_;

    for ( i=0; i<3*3; i++ ) OVI[i] = HCORE[i] = ZERO;
    AB2 = ZERO;
    for ( i=0; i<3; i++ ) {
	BA[i] = B[i]-A[i];
	AB2 += BA[i]*BA[i];
    }
    ips1 = *ips0 + (*nps_i);
    jps1 = *jps0 + (*nps_j);

    for ( ips=(*ips0); ips<ips1; ips++ ) {
	zeta_a = prim_exp[ips];
	coef_a = prim_coe[ips];
	for ( jps=(*jps0); jps<jps1; jps++ ) {
	    zeta_b = prim_exp[jps];
	    coef_b = prim_coe[jps];

	    zeta   = zeta_a + zeta_b;
	    sqrzi  = sqrt( 1.e0 / zeta );
	    zi     = sqrzi  * sqrzi;
	    coef   = coef_a * coef_b;

	    xiza  = zeta_b * zi;
	    xi    = xiza * zeta_a;
	    xiab2 = xi * AB2;
	    xi2   = 2.e0 * xi;
	    zeta2 = HALF * zi;

	    exp_ab = coef * exp( -xiab2 );
	    oviss = _pi32_ * zi * sqrzi * exp_ab;
	    keiss = xi * ( 3.e0 - 2.e0 * xiab2 ) * oviss;
	    css   = _2pi_ * zi * exp_ab;
	    if ( fabs(css) < EPS_PS_PAIR_NAI ) continue;
	    for ( i=0; i<3; i++ ) {
		P[i]  = xiza*BA[i] + A[i];
		PA[i] = xiza*BA[i];
		PB[i] = xiza*BA[i] - BA[i];
	    }
	    // overlap integral
	    // (p|s)
	    for ( i=0; i<3; i++ ) ovips[i] = PA[i]*oviss;
	    // (p|p)
	    ovipp[0*3+0] = PB[0]*ovips[0] + zeta2*oviss;
	    ovipp[0*3+1] = PB[1]*ovips[0];
	    ovipp[0*3+2] = PB[2]*ovips[0];
	    ovipp[1*3+0] = PB[0]*ovips[1];
	    ovipp[1*3+1] = PB[1]*ovips[1] + zeta2*oviss;
	    ovipp[1*3+2] = PB[2]*ovips[1];
	    ovipp[2*3+0] = PB[0]*ovips[2];
	    ovipp[2*3+1] = PB[1]*ovips[2];
	    ovipp[2*3+2] = PB[2]*ovips[2] + zeta2*oviss;
	    // kinetic energy integral
	    // (p|T|s)
	    for ( i=0; i<3; i++ ) keips[i] = PA[i]*keiss + xi2*ovips[i];
	    // (p|T|p)
	    keipp[0*3+0] = PB[0]*keips[0] + xi2*ovipp[0*3+0] + zeta2*keiss;
	    keipp[0*3+1] = PB[1]*keips[0] + xi2*ovipp[0*3+1];
	    keipp[0*3+2] = PB[2]*keips[0] + xi2*ovipp[0*3+2];
	    keipp[1*3+0] = PB[0]*keips[1] + xi2*ovipp[1*3+0];
	    keipp[1*3+1] = PB[1]*keips[1] + xi2*ovipp[1*3+1] + zeta2*keiss;
	    keipp[1*3+2] = PB[2]*keips[1] + xi2*ovipp[1*3+2];
	    keipp[2*3+0] = PB[0]*keips[2] + xi2*ovipp[2*3+0];
	    keipp[2*3+1] = PB[1]*keips[2] + xi2*ovipp[2*3+1];
	    keipp[2*3+2] = PB[2]*keips[2] + xi2*ovipp[2*3+2] + zeta2*keiss;

	    for ( i=0; i<3*3; i++ ) naipp[i] = ZERO;
	    for ( kat=0; kat<(*nat); kat++ ) {
		dq   = (double)atomic_number[kat];
		PC[0] = P[0]-atom_x[kat];
		PC[1] = P[1]-atom_y[kat];
		PC[2] = P[2]-atom_z[kat];
		PC2 = PC[0]*PC[0] + PC[1]*PC[1] + PC[2]*PC[2];
		U = zeta * PC2;
		{
		    int it0, pos;
		    double dT, t_inv, st_inv, dT2, dT3, cnass;
		    cnass = dq*css;
		    if ( U < 40e0 ) {
			it0 = (int)(0.5e0 + U * FMT_fmt_inv_step_size);
			dT  = it0 * FMT_fmt_step_size - U;
			dT2 = dT * _inv2_;
			dT3 = dT * _inv3_;
			pos = it0 * (2+4);
			naiss[0] = cnass*(((FMT_fmt_table2[pos+3] * dT3
					+ FMT_fmt_table2[pos+2] ) * dT2
					+ FMT_fmt_table2[pos+1] ) * dT
					+ FMT_fmt_table2[pos+0] );
			naiss[1] = cnass*(((FMT_fmt_table2[pos+4] * dT3
					+ FMT_fmt_table2[pos+3] ) * dT2
					+ FMT_fmt_table2[pos+2] ) * dT
					+ FMT_fmt_table2[pos+1] );
			naiss[2] = cnass*(((FMT_fmt_table2[pos+5] * dT3
					+ FMT_fmt_table2[pos+4] ) * dT2
					+ FMT_fmt_table2[pos+3] ) * dT
					+ FMT_fmt_table2[pos+2] );
		    } else {
			st_inv = sqrt( 0.5e0 / U );
			t_inv  = st_inv * st_inv;
			naiss[0] = cnass * _spi2_ * st_inv;
			naiss[1] = t_inv * naiss[0];
			naiss[2] = 3.e0 * t_inv * naiss[1];
		    }
		}
		//fmt( naiss, 2, U, dq*css );
		// -- (P|A(0)|S) --
		for ( m=0; m<=1; m++) {
		    naips[m][0] = PA[0]*naiss[m] - PC[0]*naiss[m+1];
		    naips[m][1] = PA[1]*naiss[m] - PC[1]*naiss[m+1];
		    naips[m][2] = PA[2]*naiss[m] - PC[2]*naiss[m+1];
		}
		// -- (P|A|P) --
		tmp = zeta2*( naiss[0] - naiss[1] );

		naipp[0*3+0] += PB[0]*naips[0][0]-PC[0]*naips[1][0]+tmp;
		naipp[0*3+1] += PB[1]*naips[0][0]-PC[1]*naips[1][0];
		naipp[0*3+2] += PB[2]*naips[0][0]-PC[2]*naips[1][0];
		naipp[1*3+0] += PB[0]*naips[0][1]-PC[0]*naips[1][1];
		naipp[1*3+1] += PB[1]*naips[0][1]-PC[1]*naips[1][1]+tmp;
		naipp[1*3+2] += PB[2]*naips[0][1]-PC[2]*naips[1][1];
		naipp[2*3+0] += PB[0]*naips[0][2]-PC[0]*naips[1][2];
		naipp[2*3+1] += PB[1]*naips[0][2]-PC[1]*naips[1][2];
		naipp[2*3+2] += PB[2]*naips[0][2]-PC[2]*naips[1][2]+tmp;
	    }
	    // contraction
	    for ( i=0; i<3*3; i++ ) {
		OVI[i]   += ovipp[i];
		HCORE[i] += (keipp[i] - naipp[i]);
	    }
	}	// jps
    }		// ips
}

/** １つのCS対のDSタイプの１電子積分を計算する
 * @ingroup core-oneint
 * */
void oneint_core_ds__(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const int atomic_number[],
	double OVI[], double HCORE[] ) {
    int i, m, ips, jps, kat, ips1, jps1;
    double zeta_a, zeta_b, zeta, coef_a, coef_b, coef;
    double sqrzi, zi, xiza, xi, xiab2, xi2, exp_ab, css, zeta2;
    double oviss, keiss, tmp;
    double BA[3], P[3], PA[3], PC[3], AB2;
    double keips[3], keids[6], ovips[3], ovids[6];
    double PC2, dq, U;
    double naiss[2+1], naips[2][3], naids[6];
    double _pi32_ = _PI_32_, _2pi_ = _2PI_;
    double sqr3;

    sqr3 = sqrt(3.e0);

    for ( i=0; i<6; i++ ) OVI[i] = HCORE[i] = ZERO;
    AB2 = ZERO;
    for ( i=0; i<3; i++ ) {
	BA[i] = B[i]-A[i];
	AB2 += BA[i]*BA[i];
    }
    ips1 = *ips0 + (*nps_i);
    jps1 = *jps0 + (*nps_j);

    for ( ips=(*ips0); ips<ips1; ips++ ) {
	zeta_a = prim_exp[ips];
	coef_a = prim_coe[ips];
	for ( jps=(*jps0); jps<jps1; jps++ ) {
	    zeta_b = prim_exp[jps];
	    coef_b = prim_coe[jps];

	    zeta   = zeta_a + zeta_b;
	    sqrzi  = sqrt( 1.e0 / zeta );
	    zi     = sqrzi  * sqrzi;
	    coef   = coef_a * coef_b;

	    xiza   = zeta_b * zi;
	    xi     = xiza * zeta_a;
	    xiab2  = xi * AB2;
	    xi2    = 2.e0 * xi;
	    zeta2  = HALF * zi;

	    exp_ab = coef * exp( -xiab2 );
	    oviss  = _pi32_ * zi * sqrzi * exp_ab;
	    keiss  = xi * ( 3.e0 - 2.e0 * xiab2 ) * oviss;
	    css    = _2pi_ * zi * exp_ab;
	    if ( fabs(css) < EPS_PS_PAIR_NAI ) continue;
	    for ( i=0; i<3; i++ ) {
		P[i]  = xiza*BA[i] + A[i];
		PA[i] = xiza*BA[i];
	    }
	    // overlap integral
	    // -- (P||S) --
	    ovips[0] = PA[0]*oviss;
	    ovips[1] = PA[1]*oviss;
	    ovips[2] = PA[2]*oviss;
	    // -- (D||S) --
	    ovids[0] = PA[0]*ovips[0] + zeta2*oviss;
	    ovids[1] = PA[1]*ovips[1] + zeta2*oviss;
	    ovids[2] = PA[2]*ovips[2] + zeta2*oviss;
	    ovids[3] = PA[0]*ovips[1];
	    ovids[4] = PA[1]*ovips[2];
	    ovids[5] = PA[2]*ovips[0];
	    // kinetic energy integral
	    // -- (P|T|S) --
	    keips[0] = PA[0]*keiss + xi2*ovips[0];
	    keips[1] = PA[1]*keiss + xi2*ovips[1];
	    keips[2] = PA[2]*keiss + xi2*ovips[2];
	    // -- (D|T|S) --
	    tmp = zeta2*keiss - xiza*oviss;
	    keids[0] = PA[0]*keips[0] + xi2*ovids[0] + tmp;
	    keids[1] = PA[1]*keips[1] + xi2*ovids[1] + tmp;
	    keids[2] = PA[2]*keips[2] + xi2*ovids[2] + tmp;
	    keids[3] = PA[0]*keips[1] + xi2*ovids[3];
	    keids[4] = PA[1]*keips[2] + xi2*ovids[4];
	    keids[5] = PA[2]*keips[0] + xi2*ovids[5];

	    for ( i=0; i<6; i++ ) naids[i] = ZERO;
	    for ( kat=0; kat<(*nat); kat++ ) {
		dq   = (double)atomic_number[kat];
		PC[0] = P[0]-atom_x[kat];
		PC[1] = P[1]-atom_y[kat];
		PC[2] = P[2]-atom_z[kat];
		PC2 = PC[0]*PC[0] + PC[1]*PC[1] + PC[2]*PC[2];
		U = zeta * PC2;
		{
		    int it0, pos;
		    double dT, t_inv, st_inv, dT2, dT3, cnass;
		    cnass = dq*css;
		    if ( U < 40e0 ) {
			it0 = (int)(0.5e0 + U * FMT_fmt_inv_step_size);
			dT  = it0 * FMT_fmt_step_size - U;
			dT2 = dT * _inv2_;
			dT3 = dT * _inv3_;
			pos = it0 * (2+4);
			naiss[0] = cnass*(((FMT_fmt_table2[pos+3] * dT3
					+ FMT_fmt_table2[pos+2] ) * dT2
					+ FMT_fmt_table2[pos+1] ) * dT
					+ FMT_fmt_table2[pos+0] );
			naiss[1] = cnass*(((FMT_fmt_table2[pos+4] * dT3
					+ FMT_fmt_table2[pos+3] ) * dT2
					+ FMT_fmt_table2[pos+2] ) * dT
					+ FMT_fmt_table2[pos+1] );
			naiss[2] = cnass*(((FMT_fmt_table2[pos+5] * dT3
					+ FMT_fmt_table2[pos+4] ) * dT2
					+ FMT_fmt_table2[pos+3] ) * dT
					+ FMT_fmt_table2[pos+2] );
		    } else {
			st_inv = sqrt( 0.5e0 / U );
			t_inv  = st_inv * st_inv;
			naiss[0] = cnass * _spi2_ * st_inv;
			naiss[1] = t_inv * naiss[0];
			naiss[2] = 3.e0 * t_inv * naiss[1];
		    }
		}
		//fmt( naiss, 2, U, dq*css );
		// -- (P|A(0)|S) --
		for ( m=0; m<=1; m++) {
		    naips[m][0] = PA[0]*naiss[m] - PC[0]*naiss[m+1];
		    naips[m][1] = PA[1]*naiss[m] - PC[1]*naiss[m+1];
		    naips[m][2] = PA[2]*naiss[m] - PC[2]*naiss[m+1];
		}
		// -- (D|A|S) --
		tmp = zeta2*( naiss[0] - naiss[1] );
		naids[0] += PA[0]*naips[0][0] - PC[0]*naips[1][0] + tmp;
		naids[1] += PA[1]*naips[0][1] - PC[1]*naips[1][1] + tmp;
		naids[2] += PA[2]*naips[0][2] - PC[2]*naips[1][2] + tmp;
		naids[3] += PA[0]*naips[0][1] - PC[0]*naips[1][1];
		naids[4] += PA[1]*naips[0][2] - PC[1]*naips[1][2];
		naids[5] += PA[2]*naips[0][0] - PC[2]*naips[1][0];
	    }
	    // contraction
	    for ( i=0; i<6; i++ ) {
		OVI[i]   += ovids[i];
		HCORE[i] += (keids[i] - naids[i]);
	    }
	}	// jps
    }		// ips
    for ( i=3; i<6; i++ ) {
	OVI[i]   *= sqr3;
	HCORE[i] *= sqr3;
    }
}

/** １つのCS対のDPタイプの１電子積分を計算する
 * @ingroup core-oneint
 * */
void oneint_core_dp__(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const int atomic_number[],
	double OVI[], double HCORE[] ) {
    int i, j, ij, m, ips, jps, kat, ips1, jps1;
    double zeta_a, zeta_b, zeta, coef_a, coef_b, coef;
    double sqrzi, zi, xiza, xi, xiab2, xi2, exp_ab, css, zeta2;
    double oviss, keiss, tmp, tmpps[3];
    double BA[3], P[3], PA[3], PC[3], PB[3], AB2;
    double keips[3], keids[6], keidp[6*3];
    double ovips[3], ovids[6], ovidp[6*3];
    double PC2, dq, U;
    double naiss[3+1], naips[2+1][3], naids[1+1][6], naidp[6*3];
    double _pi32_ = _PI_32_, _2pi_ = _2PI_;
    double sqr3;

    sqr3 = sqrt(3.e0);

    for ( i=0; i<6*3; i++ ) OVI[i] = HCORE[i] = ZERO;
    AB2 = ZERO;
    for ( i=0; i<3; i++ ) {
	BA[i] = B[i]-A[i];
	AB2 += BA[i]*BA[i];
    }
    ips1 = *ips0 + (*nps_i);
    jps1 = *jps0 + (*nps_j);

    for ( ips=(*ips0); ips<ips1; ips++ ) {
	zeta_a = prim_exp[ips];
	coef_a = prim_coe[ips];
	for ( jps=(*jps0); jps<jps1; jps++ ) {
	    zeta_b = prim_exp[jps];
	    coef_b = prim_coe[jps];

	    zeta   = zeta_a + zeta_b;
	    sqrzi  = sqrt( 1.e0 / zeta );
	    zi     = sqrzi  * sqrzi;
	    coef   = coef_a * coef_b;

	    xiza  = zeta_b * zi;
	    xi    = xiza * zeta_a;
	    xiab2 = xi * AB2;
	    xi2   = 2.e0 * xi;
	    zeta2 = HALF * zi;

	    exp_ab = coef * exp( -xiab2 );
	    oviss = _pi32_ * zi * sqrzi * exp_ab;
	    keiss = xi * ( 3.e0 - 2.e0 * xiab2 ) * oviss;
	    css   = _2pi_ * zi * exp_ab;
	    if ( fabs(css) < EPS_PS_PAIR_NAI ) continue;
	    for ( i=0; i<3; i++ ) {
		P[i]  = xiza*BA[i] + A[i];
		PA[i] = xiza*BA[i];
		PB[i] = xiza*BA[i] - BA[i];
	    }
	    // overlap integral
	    // -- (P||S) --
	    ovips[0] = PA[0]*oviss;
	    ovips[1] = PA[1]*oviss;
	    ovips[2] = PA[2]*oviss;
	    // -- (D||S) --
	    ovids[0] = PA[0]*ovips[0] + zeta2*oviss;
	    ovids[1] = PA[1]*ovips[1] + zeta2*oviss;
	    ovids[2] = PA[2]*ovips[2] + zeta2*oviss;
	    ovids[3] = PA[0]*ovips[1];
	    ovids[4] = PA[1]*ovips[2];
	    ovids[5] = PA[2]*ovips[0];
	    // -- (D||P) --
	    ovidp[0*3+0] = PB[0]*ovids[0] + zi   *ovips[0];
	    ovidp[0*3+1] = PB[1]*ovids[0];
	    ovidp[0*3+2] = PB[2]*ovids[0];
	    ovidp[1*3+0] = PB[0]*ovids[1];
	    ovidp[1*3+1] = PB[1]*ovids[1] + zi   *ovips[1];
	    ovidp[1*3+2] = PB[2]*ovids[1];
	    ovidp[2*3+0] = PB[0]*ovids[2];
	    ovidp[2*3+1] = PB[1]*ovids[2];
	    ovidp[2*3+2] = PB[2]*ovids[2] + zi   *ovips[2];
	    ovidp[3*3+0] = PB[0]*ovids[3] + zeta2*ovips[1];
	    ovidp[3*3+1] = PB[1]*ovids[3] + zeta2*ovips[0];
	    ovidp[3*3+2] = PB[2]*ovids[3];
	    ovidp[4*3+0] = PB[0]*ovids[4];
	    ovidp[4*3+1] = PB[1]*ovids[4] + zeta2*ovips[2];
	    ovidp[4*3+2] = PB[2]*ovids[4] + zeta2*ovips[1];
	    ovidp[5*3+0] = PB[0]*ovids[5] + zeta2*ovips[2];
	    ovidp[5*3+1] = PB[1]*ovids[5];
	    ovidp[5*3+2] = PB[2]*ovids[5] + zeta2*ovips[0];
	    // kinetic energy integral
	    // -- (P|T|S) --
	    keips[0] = PA[0]*keiss + xi2*ovips[0];
	    keips[1] = PA[1]*keiss + xi2*ovips[1];
	    keips[2] = PA[2]*keiss + xi2*ovips[2];
	    // -- (D|T|S) --
	    tmp = zeta2*keiss - xiza*oviss;

	    keids[0] = PA[0]*keips[0] + xi2*ovids[0] + tmp;
	    keids[1] = PA[1]*keips[1] + xi2*ovids[1] + tmp;
	    keids[2] = PA[2]*keips[2] + xi2*ovids[2] + tmp;
	    keids[3] = PA[0]*keips[1] + xi2*ovids[3];
	    keids[4] = PA[1]*keips[2] + xi2*ovids[4];
	    keids[5] = PA[2]*keips[0] + xi2*ovids[5];
	    // -- (D|T|P) --
	    keidp[0*3+0] = PB[0]*keids[0] + xi2*ovidp[0*3+0] + zi*keips[0];
	    keidp[0*3+1] = PB[1]*keids[0] + xi2*ovidp[0*3+1];
	    keidp[0*3+2] = PB[2]*keids[0] + xi2*ovidp[0*3+2];
	    keidp[1*3+0] = PB[0]*keids[1] + xi2*ovidp[1*3+0];
	    keidp[1*3+1] = PB[1]*keids[1] + xi2*ovidp[1*3+1] + zi*keips[1];
	    keidp[1*3+2] = PB[2]*keids[1] + xi2*ovidp[1*3+2];
	    keidp[2*3+0] = PB[0]*keids[2] + xi2*ovidp[2*3+0];
	    keidp[2*3+1] = PB[1]*keids[2] + xi2*ovidp[2*3+1];
	    keidp[2*3+2] = PB[2]*keids[2] + xi2*ovidp[2*3+2] + zi*keips[2];
	    keidp[3*3+0] = PB[0]*keids[3] + xi2*ovidp[3*3+0]
		         + zeta2*keips[1];
	    keidp[3*3+1] = PB[1]*keids[3] + xi2*ovidp[3*3+1]
			 + zeta2*keips[0];
	    keidp[3*3+2] = PB[2]*keids[3] + xi2*ovidp[3*3+2];
	    keidp[4*3+0] = PB[0]*keids[4] + xi2*ovidp[4*3+0];
	    keidp[4*3+1] = PB[1]*keids[4] + xi2*ovidp[4*3+1]
			 + zeta2*keips[2];
	    keidp[4*3+2] = PB[2]*keids[4] + xi2*ovidp[4*3+2]
			 + zeta2*keips[1];
	    keidp[5*3+0] = PB[0]*keids[5] + xi2*ovidp[5*3+0]
			 + zeta2*keips[2];
	    keidp[5*3+1] = PB[1]*keids[5] + xi2*ovidp[5*3+1];
	    keidp[5*3+2] = PB[2]*keids[5] + xi2*ovidp[5*3+2]
			 + zeta2*keips[0];

	    for ( i=0; i<6*3; i++ ) naidp[i] = ZERO;
	    for ( kat=0; kat<(*nat); kat++ ) {
		dq   = (double)atomic_number[kat];
		PC[0] = P[0]-atom_x[kat];
		PC[1] = P[1]-atom_y[kat];
		PC[2] = P[2]-atom_z[kat];
		PC2 = PC[0]*PC[0] + PC[1]*PC[1] + PC[2]*PC[2];
		U = zeta * PC2;
		{
		    int it0, pos;
		    double dT, t_inv, st_inv, dT2, dT3, cnass;
		    cnass = dq*css;
		    if ( U < 42e0 ) {
			it0 = (int)(0.5e0 + U * FMT_fmt_inv_step_size);
			dT  = it0 * FMT_fmt_step_size - U;
			dT2 = dT * _inv2_;
			dT3 = dT * _inv3_;
			pos = it0 * (3+4);
			naiss[0] = cnass*(((FMT_fmt_table3[pos+3] * dT3
					+ FMT_fmt_table3[pos+2] ) * dT2
					+ FMT_fmt_table3[pos+1] ) * dT
					+ FMT_fmt_table3[pos+0] );
			naiss[1] = cnass*(((FMT_fmt_table3[pos+4] * dT3
					+ FMT_fmt_table3[pos+3] ) * dT2
					+ FMT_fmt_table3[pos+2] ) * dT
					+ FMT_fmt_table3[pos+1] );
			naiss[2] = cnass*(((FMT_fmt_table3[pos+5] * dT3
					+ FMT_fmt_table3[pos+4] ) * dT2
					+ FMT_fmt_table3[pos+3] ) * dT
					+ FMT_fmt_table3[pos+2] );
			naiss[3] = cnass*(((FMT_fmt_table3[pos+6] * dT3
					+ FMT_fmt_table3[pos+5] ) * dT2
					+ FMT_fmt_table3[pos+4] ) * dT
					+ FMT_fmt_table3[pos+3] );
		    } else {
			st_inv = sqrt( 0.5e0 / U );
			t_inv  = st_inv * st_inv;
			naiss[0] = cnass * _spi2_ * st_inv;
			naiss[1] = t_inv * naiss[0];
			naiss[2] = 3.e0 * t_inv * naiss[1];
			naiss[3] = 5.e0 * t_inv * naiss[2];
		    }
		}
		//fmt( naiss, 3, U, dq*css );
		// -- (P|A(0)|S) --
		for ( m=0; m<=2; m++) {
		    naips[m][0] = PA[0]*naiss[m] - PC[0]*naiss[m+1];
		    naips[m][1] = PA[1]*naiss[m] - PC[1]*naiss[m+1];
		    naips[m][2] = PA[2]*naiss[m] - PC[2]*naiss[m+1];
		}
		// -- (D|A(0)|S) --
		for (m=0; m<=1; m++) {
		    tmp = naiss[m] - naiss[m+1];
		    naids[m][0] = PA[0]*naips[m][0]-PC[0]*naips[m+1][0]
				+ zeta2*tmp;
		    naids[m][1] = PA[1]*naips[m][1]-PC[1]*naips[m+1][1]
				+ zeta2*tmp;
		    naids[m][2] = PA[2]*naips[m][2]-PC[2]*naips[m+1][2]
				+ zeta2*tmp;
		    naids[m][3] = PA[0]*naips[m][1]-PC[0]*naips[m+1][1];
		    naids[m][4] = PA[1]*naips[m][2]-PC[1]*naips[m+1][2];
		    naids[m][5] = PA[2]*naips[m][0]-PC[2]*naips[m+1][0];
		}
		// -- (D|A|P) --
		tmpps[0] = naips[0][0] - naips[1][0];
		tmpps[1] = naips[0][1] - naips[1][1];
		tmpps[2] = naips[0][2] - naips[1][2];
		naidp[0*3+0] += PB[0]*naids[0][0] - PC[0]*naids[1][0]
			     + zi*tmpps[0];
		naidp[0*3+1] += PB[1]*naids[0][0] - PC[1]*naids[1][0];
		naidp[0*3+2] += PB[2]*naids[0][0] - PC[2]*naids[1][0];
		naidp[1*3+0] += PB[0]*naids[0][1] - PC[0]*naids[1][1];
		naidp[1*3+1] += PB[1]*naids[0][1] - PC[1]*naids[1][1]
			     + zi*tmpps[1];
		naidp[1*3+2] += PB[2]*naids[0][1] - PC[2]*naids[1][1];
		naidp[2*3+0] += PB[0]*naids[0][2] - PC[0]*naids[1][2];
		naidp[2*3+1] += PB[1]*naids[0][2] - PC[1]*naids[1][2];
		naidp[2*3+2] += PB[2]*naids[0][2] - PC[2]*naids[1][2]
			     + zi*tmpps[2];
		naidp[3*3+0] += PB[0]*naids[0][3] - PC[0]*naids[1][3]
			     + zeta2*tmpps[1];
		naidp[3*3+1] += PB[1]*naids[0][3] - PC[1]*naids[1][3]
			     + zeta2*tmpps[0];
		naidp[3*3+2] += PB[2]*naids[0][3] - PC[2]*naids[1][3];
		naidp[4*3+0] += PB[0]*naids[0][4] - PC[0]*naids[1][4];
		naidp[4*3+1] += PB[1]*naids[0][4] - PC[1]*naids[1][4]
			     + zeta2*tmpps[2];
		naidp[4*3+2] += PB[2]*naids[0][4] - PC[2]*naids[1][4]
			     + zeta2*tmpps[1];
		naidp[5*3+0] += PB[0]*naids[0][5] - PC[0]*naids[1][5]
			     + zeta2*tmpps[2];
		naidp[5*3+1] += PB[1]*naids[0][5] - PC[1]*naids[1][5];
		naidp[5*3+2] += PB[2]*naids[0][5] - PC[2]*naids[1][5]
			     + zeta2*tmpps[0];
	    }
	    // contraction
	    for ( i=0; i<6*3; i++ ) {
		OVI[i]   += ovidp[i];
		HCORE[i] += (keidp[i] - naidp[i]);
	    }
	}	// jps
    }		// ips
    for ( i=3; i<6; i++ ) {
	for ( j=0; j<3; j++ ) {
	    ij = i*3+j;
	    OVI[ij]   *= sqr3;
	    HCORE[ij] *= sqr3;
	}
    }
}

/** １つのCS対のDDタイプの１電子積分を計算する
 * @ingroup core-oneint
 * */
void oneint_core_dd__(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const int atomic_number[],
	double OVI[], double HCORE[] ) {
    int i, j, ij, m, ips, jps, kat, ips1, jps1;
    double zeta_a, zeta_b, zeta, coef_a, coef_b, coef;
    double sqrzi, zi, xiza, xizb, xi, xiab2, xi2, exp_ab, css, zeta2;
    double oviss, keiss, tmp, tmpps[3], tmppp[3*3], tmpds[6];
    double BA[3], P[3], PA[3], PC[3], PB[3], AB2;
    double keips[3], keids[6], keipp[3*3], keidp[6*3], keidd[6*6];
    double ovips[3], ovids[6], ovipp[3*3], ovidp[6*3], ovidd[6*6];
    double PC2, dq, U;
    double naiss[4+1], naips[3+1][3], naids[2+1][6], naidp[1+1][6*3];
    double naipp[1+1][3*3], naidd[6*6];
    double _pi32_ = _PI_32_, _2pi_ = _2PI_;
    double sqr3, coe_a, coe;

    sqr3  = sqrt(3.e0);

    for ( i=0; i<6*6; i++ ) OVI[i] = HCORE[i] = ZERO;
    AB2 = ZERO;
    for ( i=0; i<3; i++ ) {
	BA[i] = B[i]-A[i];
	AB2 += BA[i]*BA[i];
    }
    ips1 = *ips0 + (*nps_i);
    jps1 = *jps0 + (*nps_j);

    for ( ips=(*ips0); ips<ips1; ips++ ) {
	zeta_a = prim_exp[ips];
	coef_a = prim_coe[ips];
	for ( jps=(*jps0); jps<jps1; jps++ ) {
	    zeta_b = prim_exp[jps];
	    coef_b = prim_coe[jps];

	    zeta   = zeta_a + zeta_b;
	    sqrzi  = sqrt( 1.e0 / zeta );
	    zi     = sqrzi  * sqrzi;
	    coef   = coef_a * coef_b;

	    xiza  = zeta_b * zi;
	    xizb  = zeta_a * zi;
	    xi    = xiza * zeta_a;
	    xiab2 = xi * AB2;
	    xi2   = 2.e0 * xi;
	    zeta2 = HALF * zi;

	    exp_ab = coef * exp( -xiab2 );
	    oviss = _pi32_ * zi * sqrzi * exp_ab;
	    keiss = xi * ( 3.e0 - 2.e0 * xiab2 ) * oviss;
	    css   = _2pi_ * zi * exp_ab;

	    if ( fabs(css) < EPS_PS_PAIR_NAI ) continue;
	    for ( i=0; i<3; i++ ) {
		P[i]  = xiza*BA[i] + A[i];
		PA[i] = xiza*BA[i];
		PB[i] = xiza*BA[i] - BA[i];
	    }
	    // overlap integral
	    // -- (P||S) --
	    ovips[0] = PA[0]*oviss;
	    ovips[1] = PA[1]*oviss;
	    ovips[2] = PA[2]*oviss;
	    // -- (D||S) --
	    ovids[0] = PA[0]*ovips[0] + zeta2*oviss;
	    ovids[1] = PA[1]*ovips[1] + zeta2*oviss;
	    ovids[2] = PA[2]*ovips[2] + zeta2*oviss;
	    ovids[3] = PA[0]*ovips[1];
	    ovids[4] = PA[1]*ovips[2];
	    ovids[5] = PA[2]*ovips[0];
	    // -- (P||P) --
	    ovipp[0*3+0] = PB[0]*ovips[0] + zeta2*oviss;
	    ovipp[0*3+1] = PB[1]*ovips[0];
	    ovipp[0*3+2] = PB[2]*ovips[0];
	    ovipp[1*3+0] = PB[0]*ovips[1];
	    ovipp[1*3+1] = PB[1]*ovips[1] + zeta2*oviss;
	    ovipp[1*3+2] = PB[2]*ovips[1];
	    ovipp[2*3+0] = PB[0]*ovips[2];
	    ovipp[2*3+1] = PB[1]*ovips[2];
	    ovipp[2*3+2] = PB[2]*ovips[2] + zeta2*oviss;
	    // -- (D||P) --
	    ovidp[0*3+0] = PB[0]*ovids[0] + zi   *ovips[0];
	    ovidp[0*3+1] = PB[1]*ovids[0];
	    ovidp[0*3+2] = PB[2]*ovids[0];
	    ovidp[1*3+0] = PB[0]*ovids[1];
	    ovidp[1*3+1] = PB[1]*ovids[1] + zi   *ovips[1];
	    ovidp[1*3+2] = PB[2]*ovids[1];
	    ovidp[2*3+0] = PB[0]*ovids[2];
	    ovidp[2*3+1] = PB[1]*ovids[2];
	    ovidp[2*3+2] = PB[2]*ovids[2] + zi   *ovips[2];
	    ovidp[3*3+0] = PB[0]*ovids[3] + zeta2*ovips[1];
	    ovidp[3*3+1] = PB[1]*ovids[3] + zeta2*ovips[0];
	    ovidp[3*3+2] = PB[2]*ovids[3];
	    ovidp[4*3+0] = PB[0]*ovids[4];
	    ovidp[4*3+1] = PB[1]*ovids[4] + zeta2*ovips[2];
	    ovidp[4*3+2] = PB[2]*ovids[4] + zeta2*ovips[1];
	    ovidp[5*3+0] = PB[0]*ovids[5] + zeta2*ovips[2];
	    ovidp[5*3+1] = PB[1]*ovids[5];
	    ovidp[5*3+2] = PB[2]*ovids[5] + zeta2*ovips[0];
	    // -- (D||D) --
	    ovidd[0*6+0] = PB[0]*ovidp[0*3+0] + zeta2*ovids[0]
			 + zi*ovipp[0*3+0];
	    ovidd[0*6+1] = PB[1]*ovidp[0*3+1] + zeta2*ovids[0];
	    ovidd[0*6+2] = PB[2]*ovidp[0*3+2] + zeta2*ovids[0];
	    ovidd[0*6+3] = PB[0]*ovidp[0*3+1] + zi*ovipp[0*3+1];
	    ovidd[0*6+4] = PB[1]*ovidp[0*3+2];
	    ovidd[0*6+5] = PB[2]*ovidp[0*3+0];
	    ovidd[1*6+0] = PB[0]*ovidp[1*3+0] + zeta2*ovids[1];
	    ovidd[1*6+1] = PB[1]*ovidp[1*3+1] + zeta2*ovids[1]
			 + zi*ovipp[1*3+1];
	    ovidd[1*6+2] = PB[2]*ovidp[1*3+2] + zeta2*ovids[1];
	    ovidd[1*6+3] = PB[0]*ovidp[1*3+1];
	    ovidd[1*6+4] = PB[1]*ovidp[1*3+2] + zi*ovipp[1*3+2];
	    ovidd[1*6+5] = PB[2]*ovidp[1*3+0];
	    ovidd[2*6+0] = PB[0]*ovidp[2*3+0] + zeta2*ovids[2];
	    ovidd[2*6+1] = PB[1]*ovidp[2*3+1] + zeta2*ovids[2];
	    ovidd[2*6+2] = PB[2]*ovidp[2*3+2] + zeta2*ovids[2]
			 + zi*ovipp[2*3+2];
	    ovidd[2*6+3] = PB[0]*ovidp[2*3+1];
	    ovidd[2*6+4] = PB[1]*ovidp[2*3+2];
	    ovidd[2*6+5] = PB[2]*ovidp[2*3+0] + zi*ovipp[2*3+0];
	    ovidd[3*6+0] = PB[0]*ovidp[3*3+0] + zeta2*ovids[3]
			 + zeta2*ovipp[1*3+0];
	    ovidd[3*6+1] = PB[1]*ovidp[3*3+1] + zeta2*ovids[3]
			 + zeta2*ovipp[0*3+1];
	    ovidd[3*6+2] = PB[2]*ovidp[3*3+2] + zeta2*ovids[3];
	    ovidd[3*6+3] = PB[0]*ovidp[3*3+1] + zeta2*ovipp[1*3+1];
	    ovidd[3*6+4] = PB[1]*ovidp[3*3+2] + zeta2*ovipp[0*3+2];
	    ovidd[3*6+5] = PB[2]*ovidp[3*3+0];
	    ovidd[4*6+0] = PB[0]*ovidp[4*3+0] + zeta2*ovids[4];
	    ovidd[4*6+1] = PB[1]*ovidp[4*3+1] + zeta2*ovids[4]
			 + zeta2*ovipp[2*3+1];
	    ovidd[4*6+2] = PB[2]*ovidp[4*3+2] + zeta2*ovids[4]
			 + zeta2*ovipp[1*3+2];
	    ovidd[4*6+3] = PB[0]*ovidp[4*3+1];
	    ovidd[4*6+4] = PB[1]*ovidp[4*3+2] + zeta2*ovipp[2*3+2];
	    ovidd[4*6+5] = PB[2]*ovidp[4*3+0] + zeta2*ovipp[1*3+0];
	    ovidd[5*6+0] = PB[0]*ovidp[5*3+0] + zeta2*ovids[5]
			 + zeta2*ovipp[2*3+0];
	    ovidd[5*6+1] = PB[1]*ovidp[5*3+1] + zeta2*ovids[5];
	    ovidd[5*6+2] = PB[2]*ovidp[5*3+2] + zeta2*ovids[5]
			 + zeta2*ovipp[0*3+2];
	    ovidd[5*6+3] = PB[0]*ovidp[5*3+1] + zeta2*ovipp[2*3+1];
	    ovidd[5*6+4] = PB[1]*ovidp[5*3+2];
	    ovidd[5*6+5] = PB[2]*ovidp[5*3+0] + zeta2*ovipp[0*3+0];
	    // kinetic energy integral
	    // -- (P|T|S) --
	    keips[0] = PA[0]*keiss + xi2*ovips[0];
	    keips[1] = PA[1]*keiss + xi2*ovips[1];
	    keips[2] = PA[2]*keiss + xi2*ovips[2];
	    // -- (D|T|S) --
	    tmp = zeta2*keiss - xiza*oviss;
	    keids[0] = PA[0]*keips[0] + xi2*ovids[0] + tmp;
	    keids[1] = PA[1]*keips[1] + xi2*ovids[1] + tmp;
	    keids[2] = PA[2]*keips[2] + xi2*ovids[2] + tmp;
	    keids[3] = PA[0]*keips[1] + xi2*ovids[3];
	    keids[4] = PA[1]*keips[2] + xi2*ovids[4];
	    keids[5] = PA[2]*keips[0] + xi2*ovids[5];
	    // -- (P|T|P) --
	    keipp[0*3+0]=PB[0]*keips[0]+xi2*ovipp[0*3+0]+zeta2*keiss;
	    keipp[0*3+1]=PB[1]*keips[0]+xi2*ovipp[0*3+1];
	    keipp[0*3+2]=PB[2]*keips[0]+xi2*ovipp[0*3+2];
	    keipp[1*3+0]=PB[0]*keips[1]+xi2*ovipp[1*3+0];
	    keipp[1*3+1]=PB[1]*keips[1]+xi2*ovipp[1*3+1]+zeta2*keiss;
	    keipp[1*3+2]=PB[2]*keips[1]+xi2*ovipp[1*3+2];
	    keipp[2*3+0]=PB[0]*keips[2]+xi2*ovipp[2*3+0];
	    keipp[2*3+1]=PB[1]*keips[2]+xi2*ovipp[2*3+1];
	    keipp[2*3+2]=PB[2]*keips[2]+xi2*ovipp[2*3+2]+zeta2*keiss;
	    // -- (D|T|P) --
	    keidp[0*3+0]=PB[0]*keids[0]+xi2*ovidp[0*3+0]+zi*keips[0];
	    keidp[0*3+1]=PB[1]*keids[0]+xi2*ovidp[0*3+1];
	    keidp[0*3+2]=PB[2]*keids[0]+xi2*ovidp[0*3+2];
	    keidp[1*3+0]=PB[0]*keids[1]+xi2*ovidp[1*3+0];
	    keidp[1*3+1]=PB[1]*keids[1]+xi2*ovidp[1*3+1]+zi*keips[1];
	    keidp[1*3+2]=PB[2]*keids[1]+xi2*ovidp[1*3+2];
	    keidp[2*3+0]=PB[0]*keids[2]+xi2*ovidp[2*3+0];
	    keidp[2*3+1]=PB[1]*keids[2]+xi2*ovidp[2*3+1];
	    keidp[2*3+2]=PB[2]*keids[2]+xi2*ovidp[2*3+2]+zi*keips[2];
	    keidp[3*3+0]=PB[0]*keids[3]+xi2*ovidp[3*3+0]
		        +zeta2*keips[1];
	    keidp[3*3+1]=PB[1]*keids[3]+xi2*ovidp[3*3+1]
			+zeta2*keips[0];
	    keidp[3*3+2]=PB[2]*keids[3]+xi2*ovidp[3*3+2];
	    keidp[4*3+0]=PB[0]*keids[4]+xi2*ovidp[4*3+0];
	    keidp[4*3+1]=PB[1]*keids[4]+xi2*ovidp[4*3+1]
			+zeta2*keips[2];
	    keidp[4*3+2]=PB[2]*keids[4]+xi2*ovidp[4*3+2]
			+zeta2*keips[1];
	    keidp[5*3+0]=PB[0]*keids[5]+xi2*ovidp[5*3+0]
			+zeta2*keips[2];
	    keidp[5*3+1]=PB[1]*keids[5]+xi2*ovidp[5*3+1];
	    keidp[5*3+2]=PB[2]*keids[5]+xi2*ovidp[5*3+2]
			+zeta2*keips[0];
	    // -- (D|T|D) --
	    tmpds[0] = zeta2*keids[0] - xizb*ovids[0];
	    tmpds[1] = zeta2*keids[1] - xizb*ovids[1];
	    tmpds[2] = zeta2*keids[2] - xizb*ovids[2];
	    tmpds[3] = zeta2*keids[3] - xizb*ovids[3];
	    tmpds[4] = zeta2*keids[4] - xizb*ovids[4];
	    tmpds[5] = zeta2*keids[5] - xizb*ovids[5];

	    keidd[0*6+0] = PB[0]*keidp[0*3+0]+xi2*ovidd[0*6+0]
			+  tmpds[0]+zi*keipp[0*3+0];
	    keidd[0*6+1] = PB[1]*keidp[0*3+1]+xi2*ovidd[0*6+1]+tmpds[0];
	    keidd[0*6+2] = PB[2]*keidp[0*3+2]+xi2*ovidd[0*6+2]+tmpds[0];
	    keidd[0*6+3] = PB[0]*keidp[0*3+1]+xi2*ovidd[0*6+3]
			+zi*keipp[0*3+1];
	    keidd[0*6+4] = PB[1]*keidp[0*3+2]+xi2*ovidd[0*6+4];
	    keidd[0*6+5] = PB[2]*keidp[0*3+0]+xi2*ovidd[0*6+5];
	    keidd[1*6+0] = PB[0]*keidp[1*3+0]+xi2*ovidd[1*6+0]+tmpds[1];
	    keidd[1*6+1] = PB[1]*keidp[1*3+1]+xi2*ovidd[1*6+1]
			+  tmpds[1]+zi*keipp[1*3+1];
	    keidd[1*6+2] = PB[2]*keidp[1*3+2]+xi2*ovidd[1*6+2]+tmpds[1];
	    keidd[1*6+3] = PB[0]*keidp[1*3+1]+xi2*ovidd[1*6+3];
	    keidd[1*6+4] = PB[1]*keidp[1*3+2]+xi2*ovidd[1*6+4]
			+zi*keipp[1*3+2];
	    keidd[1*6+5] = PB[2]*keidp[1*3+0]+xi2*ovidd[1*6+5];
	    keidd[2*6+0] = PB[0]*keidp[2*3+0]+xi2*ovidd[2*6+0]+tmpds[2];
	    keidd[2*6+1] = PB[1]*keidp[2*3+1]+xi2*ovidd[2*6+1]+tmpds[2];
	    keidd[2*6+2] = PB[2]*keidp[2*3+2]+xi2*ovidd[2*6+2]
			+  tmpds[2]+zi*keipp[2*3+2];
	    keidd[2*6+3] = PB[0]*keidp[2*3+1]+xi2*ovidd[2*6+3];
	    keidd[2*6+4] = PB[1]*keidp[2*3+2]+xi2*ovidd[2*6+4];
	    keidd[2*6+5] = PB[2]*keidp[2*3+0]+xi2*ovidd[2*6+5]
			+zi*keipp[2*3+0];
	    keidd[3*6+0] = PB[0]*keidp[3*3+0]+xi2*ovidd[3*6+0]
			+  tmpds[3]+zeta2*keipp[1*3+0];
	    keidd[3*6+1] = PB[1]*keidp[3*3+1]+xi2*ovidd[3*6+1]
			+  tmpds[3]+zeta2*keipp[0*3+1];
	    keidd[3*6+2] = PB[2]*keidp[3*3+2]+xi2*ovidd[3*6+2]+tmpds[3];
	    keidd[3*6+3] = PB[0]*keidp[3*3+1]+xi2*ovidd[3*6+3]
			+zeta2*keipp[1*3+1];
	    keidd[3*6+4] = PB[1]*keidp[3*3+2]+xi2*ovidd[3*6+4]
			+zeta2*keipp[0*3+2];
	    keidd[3*6+5] = PB[2]*keidp[3*3+0]+xi2*ovidd[3*6+5];
	    keidd[4*6+0] = PB[0]*keidp[4*3+0]+xi2*ovidd[4*6+0]+tmpds[4];
	    keidd[4*6+1] = PB[1]*keidp[4*3+1]+xi2*ovidd[4*6+1]
			+  tmpds[4]+zeta2*keipp[2*3+1];
	    keidd[4*6+2] = PB[2]*keidp[4*3+2]+xi2*ovidd[4*6+2]
			+  tmpds[4]+zeta2*keipp[1*3+2];
	    keidd[4*6+3] = PB[0]*keidp[4*3+1]+xi2*ovidd[4*6+3];
	    keidd[4*6+4] = PB[1]*keidp[4*3+2]+xi2*ovidd[4*6+4]
			+zeta2*keipp[2*3+2];
	    keidd[4*6+5] = PB[2]*keidp[4*3+0]+xi2*ovidd[4*6+5]
			+zeta2*keipp[1*3+0];
	    keidd[5*6+0] = PB[0]*keidp[5*3+0]+xi2*ovidd[5*6+0]
			+  tmpds[5]+zeta2*keipp[2*3+0];
	    keidd[5*6+1] = PB[1]*keidp[5*3+1]+xi2*ovidd[5*6+1]+tmpds[5];
	    keidd[5*6+2] = PB[2]*keidp[5*3+2]+xi2*ovidd[5*6+2]
			+  tmpds[5]+zeta2*keipp[0*3+2];
	    keidd[5*6+3] = PB[0]*keidp[5*3+1]+xi2*ovidd[5*6+3]
			+zeta2*keipp[2*3+1];
	    keidd[5*6+4] = PB[1]*keidp[5*3+2]+xi2*ovidd[5*6+4];
	    keidd[5*6+5] = PB[2]*keidp[5*3+0]+xi2*ovidd[5*6+5]
			+zeta2*keipp[0*3+0];

	    for ( i=0; i<6*6; i++ ) naidd[i] = ZERO;
	    for ( kat=0; kat<(*nat); kat++ ) {
		dq   = (double)atomic_number[kat];
		PC[0] = P[0]-atom_x[kat];
		PC[1] = P[1]-atom_y[kat];
		PC[2] = P[2]-atom_z[kat];
		PC2 = PC[0]*PC[0] + PC[1]*PC[1] + PC[2]*PC[2];
		U = zeta * PC2;
		{
		    int it0, pos;
		    double dT, t_inv, st_inv, dT2, dT3, cnass;
		    cnass = dq*css;
		    if ( U < 44e0 ) {
			it0 = (int)(0.5e0 + U * FMT_fmt_inv_step_size);
			dT  = it0 * FMT_fmt_step_size - U;
			dT2 = dT * _inv2_;
			dT3 = dT * _inv3_;
			pos = it0 * (4+4);
			naiss[0] = cnass*(((FMT_fmt_table4[pos+3] * dT3
					+ FMT_fmt_table4[pos+2] ) * dT2
					+ FMT_fmt_table4[pos+1] ) * dT
					+ FMT_fmt_table4[pos+0] );
			naiss[1] = cnass*(((FMT_fmt_table4[pos+4] * dT3
					+ FMT_fmt_table4[pos+3] ) * dT2
					+ FMT_fmt_table4[pos+2] ) * dT
					+ FMT_fmt_table4[pos+1] );
			naiss[2] = cnass*(((FMT_fmt_table4[pos+5] * dT3
					+ FMT_fmt_table4[pos+4] ) * dT2
					+ FMT_fmt_table4[pos+3] ) * dT
					+ FMT_fmt_table4[pos+2] );
			naiss[3] = cnass*(((FMT_fmt_table4[pos+6] * dT3
					+ FMT_fmt_table4[pos+5] ) * dT2
					+ FMT_fmt_table4[pos+4] ) * dT
					+ FMT_fmt_table4[pos+3] );
			naiss[4] = cnass*(((FMT_fmt_table4[pos+7] * dT3
					+ FMT_fmt_table4[pos+6] ) * dT2
					+ FMT_fmt_table4[pos+5] ) * dT
					+ FMT_fmt_table4[pos+4] );
		    } else {
			st_inv = sqrt( 0.5e0 / U );
			t_inv  = st_inv * st_inv;
			naiss[0] = cnass * _spi2_ * st_inv;
			naiss[1] = t_inv * naiss[0];
			naiss[2] = 3.e0 * t_inv * naiss[1];
			naiss[3] = 5.e0 * t_inv * naiss[2];
			naiss[4] = 7.e0 * t_inv * naiss[3];
		    }
		}
		//fmt( naiss, 4, U, dq*css );
		// -- (P|A(0)|S) --
		for ( m=0; m<=3; m++) {
		    naips[m][0] = PA[0]*naiss[m] - PC[0]*naiss[m+1];
		    naips[m][1] = PA[1]*naiss[m] - PC[1]*naiss[m+1];
		    naips[m][2] = PA[2]*naiss[m] - PC[2]*naiss[m+1];
		}
		// -- (D|A(0)|S) --
		for (m=0; m<=2; m++) {
		    tmp = naiss[m] - naiss[m+1];
		    naids[m][0] = PA[0]*naips[m][0]-PC[0]*naips[m+1][0]
				+ zeta2*tmp;
		    naids[m][1] = PA[1]*naips[m][1]-PC[1]*naips[m+1][1]
				+ zeta2*tmp;
		    naids[m][2] = PA[2]*naips[m][2]-PC[2]*naips[m+1][2]
				+ zeta2*tmp;
		    naids[m][3] = PA[0]*naips[m][1]-PC[0]*naips[m+1][1];
		    naids[m][4] = PA[1]*naips[m][2]-PC[1]*naips[m+1][2];
		    naids[m][5] = PA[2]*naips[m][0]-PC[2]*naips[m+1][0];
		}
		// -- (P|A|P) --
		for (m=0; m<=1; m++) {
		    tmp = naiss[m] - naiss[m+1];
		    naipp[m][0*3+0]=PB[0]*naips[m][0]-PC[0]*naips[m+1][0]
				   +zeta2*tmp;
		    naipp[m][0*3+1]=PB[1]*naips[m][0]-PC[1]*naips[m+1][0];
		    naipp[m][0*3+2]=PB[2]*naips[m][0]-PC[2]*naips[m+1][0];
		    naipp[m][1*3+0]=PB[0]*naips[m][1]-PC[0]*naips[m+1][1];
		    naipp[m][1*3+1]=PB[1]*naips[m][1]-PC[1]*naips[m+1][1]
				   +zeta2*tmp;
		    naipp[m][1*3+2]=PB[2]*naips[m][1]-PC[2]*naips[m+1][1];
		    naipp[m][2*3+0]=PB[0]*naips[m][2]-PC[0]*naips[m+1][2];
		    naipp[m][2*3+1]=PB[1]*naips[m][2]-PC[1]*naips[m+1][2];
		    naipp[m][2*3+2]=PB[2]*naips[m][2]-PC[2]*naips[m+1][2]
				   +zeta2*tmp;
		}
		// -- (D|A|P) --
		for (m=0; m<=1; m++) {
		    tmpps[0] = naips[m][0]-naips[m+1][0];
		    tmpps[1] = naips[m][1]-naips[m+1][1];
		    tmpps[2] = naips[m][2]-naips[m+1][2];

		    naidp[m][0*3+0]=PB[0]*naids[m][0]-PC[0]*naids[m+1][0]
				   +zi*tmpps[0];
		    naidp[m][0*3+1]=PB[1]*naids[m][0]-PC[1]*naids[m+1][0];
		    naidp[m][0*3+2]=PB[2]*naids[m][0]-PC[2]*naids[m+1][0];
		    naidp[m][1*3+0]=PB[0]*naids[m][1]-PC[0]*naids[m+1][1];
		    naidp[m][1*3+1]=PB[1]*naids[m][1]-PC[1]*naids[m+1][1]
				   +zi*tmpps[1];
		    naidp[m][1*3+2]=PB[2]*naids[m][1]-PC[2]*naids[m+1][1];
		    naidp[m][2*3+0]=PB[0]*naids[m][2]-PC[0]*naids[m+1][2];
		    naidp[m][2*3+1]=PB[1]*naids[m][2]-PC[1]*naids[m+1][2];
		    naidp[m][2*3+2]=PB[2]*naids[m][2]-PC[2]*naids[m+1][2]
				   +zi*tmpps[2];
		    naidp[m][3*3+0]=PB[0]*naids[m][3]-PC[0]*naids[m+1][3]
				   + zeta2*tmpps[1];
		    naidp[m][3*3+1]=PB[1]*naids[m][3]-PC[1]*naids[m+1][3]
				   +zeta2*tmpps[0];
		    naidp[m][3*3+2]=PB[2]*naids[m][3]-PC[2]*naids[m+1][3];
		    naidp[m][4*3+0]=PB[0]*naids[m][4]-PC[0]*naids[m+1][4];
		    naidp[m][4*3+1]=PB[1]*naids[m][4]-PC[1]*naids[m+1][4]
				   +zeta2*tmpps[2];
		    naidp[m][4*3+2]=PB[2]*naids[m][4]-PC[2]*naids[m+1][4]
				   +zeta2*tmpps[1];
		    naidp[m][5*3+0]=PB[0]*naids[m][5]-PC[0]*naids[m+1][5]
				   +zeta2*tmpps[2];
		    naidp[m][5*3+1]=PB[1]*naids[m][5]-PC[1]*naids[m+1][5];
		    naidp[m][5*3+2]=PB[2]*naids[m][5]-PC[2]*naids[m+1][5]
				   +zeta2*tmpps[0];
		}
		// -- (D|A|D) --
		tmpds[0]=naids[0][0]-naids[1][0];
		tmpds[1]=naids[0][1]-naids[1][1];
		tmpds[2]=naids[0][2]-naids[1][2];
		tmpds[3]=naids[0][3]-naids[1][3];
		tmpds[4]=naids[0][4]-naids[1][4];
		tmpds[5]=naids[0][5]-naids[1][5];

		tmppp[0*3+0]=naipp[0][0*3+0]-naipp[1][0*3+0];
		tmppp[0*3+1]=naipp[0][0*3+1]-naipp[1][0*3+1];
		tmppp[0*3+2]=naipp[0][0*3+2]-naipp[1][0*3+2];
		tmppp[1*3+0]=naipp[0][1*3+0]-naipp[1][1*3+0];
		tmppp[1*3+1]=naipp[0][1*3+1]-naipp[1][1*3+1];
		tmppp[1*3+2]=naipp[0][1*3+2]-naipp[1][1*3+2];
		tmppp[2*3+0]=naipp[0][2*3+0]-naipp[1][2*3+0];
		tmppp[2*3+1]=naipp[0][2*3+1]-naipp[1][2*3+1];
		tmppp[2*3+2]=naipp[0][2*3+2]-naipp[1][2*3+2];

		naidd[0*6+0]+=PB[0]*naidp[0][0*3+0]-PC[0]*naidp[1][0*3+0]
			    + zeta2*tmpds[0]+zi*tmppp[0*3+0];
		naidd[0*6+1]+=PB[1]*naidp[0][0*3+1]-PC[1]*naidp[1][0*3+1]
			    + zeta2*tmpds[0];
		naidd[0*6+2]+=PB[2]*naidp[0][0*3+2]-PC[2]*naidp[1][0*3+2]
			    + zeta2*tmpds[0];
		naidd[0*6+3]+=PB[0]*naidp[0][0*3+1]-PC[0]*naidp[1][0*3+1]
			    + zi*tmppp[0*3+1];
		naidd[0*6+4]+=PB[1]*naidp[0][0*3+2]-PC[1]*naidp[1][0*3+2];
		naidd[0*6+5]+=PB[2]*naidp[0][0*3+0]-PC[2]*naidp[1][0*3+0];
		naidd[1*6+0]+=PB[0]*naidp[0][1*3+0]-PC[0]*naidp[1][1*3+0]
			    + zeta2*tmpds[1];
		naidd[1*6+1]+=PB[1]*naidp[0][1*3+1]-PC[1]*naidp[1][1*3+1]
			    + zeta2*tmpds[1]+zi*tmppp[1*3+1];
		naidd[1*6+2]+=PB[2]*naidp[0][1*3+2]-PC[2]*naidp[1][1*3+2]
			    + zeta2*tmpds[1];
		naidd[1*6+3]+=PB[0]*naidp[0][1*3+1]-PC[0]*naidp[1][1*3+1];
		naidd[1*6+4]+=PB[1]*naidp[0][1*3+2]-PC[1]*naidp[1][1*3+2]
			    + zi*tmppp[1*3+2];
		naidd[1*6+5]+=PB[2]*naidp[0][1*3+0]-PC[2]*naidp[1][1*3+0];
		naidd[2*6+0]+=PB[0]*naidp[0][2*3+0]-PC[0]*naidp[1][2*3+0]
			    + zeta2*tmpds[2];
		naidd[2*6+1]+=PB[1]*naidp[0][2*3+1]-PC[1]*naidp[1][2*3+1]
			    + zeta2*tmpds[2];
		naidd[2*6+2]+=PB[2]*naidp[0][2*3+2]-PC[2]*naidp[1][2*3+2]
			    + zeta2*tmpds[2]+zi*tmppp[2*3+2];
		naidd[2*6+3]+=PB[0]*naidp[0][2*3+1]-PC[0]*naidp[1][2*3+1];
		naidd[2*6+4]+=PB[1]*naidp[0][2*3+2]-PC[1]*naidp[1][2*3+2];
		naidd[2*6+5]+=PB[2]*naidp[0][2*3+0]-PC[2]*naidp[1][2*3+0]
			    + zi*tmppp[2*3+0];
		naidd[3*6+0]+=PB[0]*naidp[0][3*3+0]-PC[0]*naidp[1][3*3+0]
			    + zeta2*tmpds[3]+ zeta2*tmppp[1*3+0];
		naidd[3*6+1]+=PB[1]*naidp[0][3*3+1]-PC[1]*naidp[1][3*3+1]
			    + zeta2*tmpds[3]+ zeta2*tmppp[0*3+1];
		naidd[3*6+2]+=PB[2]*naidp[0][3*3+2]-PC[2]*naidp[1][3*3+2]
			    + zeta2*tmpds[3];
		naidd[3*6+3]+=PB[0]*naidp[0][3*3+1]-PC[0]*naidp[1][3*3+1]
			    + zeta2*tmppp[1*3+1];
		naidd[3*6+4]+=PB[1]*naidp[0][3*3+2]-PC[1]*naidp[1][3*3+2]
			    + zeta2*tmppp[0*3+2];
		naidd[3*6+5]+=PB[2]*naidp[0][3*3+0]-PC[2]*naidp[1][3*3+0];
		naidd[4*6+0]+=PB[0]*naidp[0][4*3+0]-PC[0]*naidp[1][4*3+0]
			    + zeta2*tmpds[4];
		naidd[4*6+1]+=PB[1]*naidp[0][4*3+1]-PC[1]*naidp[1][4*3+1]
			    + zeta2*tmpds[4]+ zeta2*tmppp[2*3+1];
		naidd[4*6+2]+=PB[2]*naidp[0][4*3+2]-PC[2]*naidp[1][4*3+2]
			    + zeta2*tmpds[4]+ zeta2*tmppp[1*3+2];
		naidd[4*6+3]+=PB[0]*naidp[0][4*3+1]-PC[0]*naidp[1][4*3+1];
		naidd[4*6+4]+=PB[1]*naidp[0][4*3+2]-PC[1]*naidp[1][4*3+2]
			    + zeta2*tmppp[2*3+2];
		naidd[4*6+5]+=PB[2]*naidp[0][4*3+0]-PC[2]*naidp[1][4*3+0]
			    + zeta2*tmppp[1*3+0];
		naidd[5*6+0]+=PB[0]*naidp[0][5*3+0]-PC[0]*naidp[1][5*3+0]
			    + zeta2*tmpds[5]+ zeta2*tmppp[2*3+0];
		naidd[5*6+1]+=PB[1]*naidp[0][5*3+1]-PC[1]*naidp[1][5*3+1]
			    + zeta2*tmpds[5];
		naidd[5*6+2]+=PB[2]*naidp[0][5*3+2]-PC[2]*naidp[1][5*3+2]
			    + zeta2*tmpds[5]+ zeta2*tmppp[0*3+2];
		naidd[5*6+3]+=PB[0]*naidp[0][5*3+1]-PC[0]*naidp[1][5*3+1]
			    + zeta2*tmppp[2*3+1];
		naidd[5*6+4]+=PB[1]*naidp[0][5*3+2]-PC[1]*naidp[1][5*3+2];
		naidd[5*6+5]+=PB[2]*naidp[0][5*3+0]-PC[2]*naidp[1][5*3+0]
			    + zeta2*tmppp[0*3+0];
	    }
	    // contraction
	    for ( i=0; i<6*6; i++ ) {
		OVI[i]   += ovidd[i];
		HCORE[i] += (keidd[i]-naidd[i]);
	    }
	}	// jps
    }		// ips
    ij = 0;
    for ( i=0; i<6; i++ ) {
	coe_a = (i<3? 1.e0 : sqr3 );
	for ( j=0; j<6; j++ ) {
	    coe = coe_a * (j<3? 1.e0 : sqr3 );
	    OVI[ij]   *= coe;
	    HCORE[ij] *= coe;
	    ij++;
	}
    }
}
