/**
 * @file ofmo-twoint-core-ssss-dspp.c
 * １つのCS４重対に対する２電子積分を計算する関数群。
 * 2011/06/16現在、(ss,ss)～(dd,dd)までの２１種類の２電子積分
 * 計算をサポートしている。
 * このファイルには、(ss,ss)～(ds,pp)タイプの２電子積分計算を
 * 行う関数が記述されている。
 * */
/**
 * @defgroup core-twoint 縮約２電子積分を計算する関数群
 * (ss,ss)～(dd,dd)の１つの縮約２電子積分を計算する関数群
 *
 * Fortranとの親和性を考慮して、すべてポインタ変数になっている。
 *
 * 各関数で、共通の引数を持っているので、それについての説明を記す。
 *
 * @param[in] *nijps １つ目のCSペアに含まれるPSペア数
 * @param[in] vzeta[ijps] １つ目のCSペア中の \c ijps 番目のPSペアの
 *      軌道指数和\f$ \zeta = \zeta_a + \zeta_b \f$
 * @param[in] vdkab[ijps] １つ目のCSペア中の \c ijps 番目のPSペアの
 *      結合係数
 *     \f[ K_{ab} = \sqrt2 \pi^{5/4} \frac1{\zeta_a+\zeta_b}
 *     \exp\left[ -\frac{\zeta_a \zeta_b}{\zeta_a + \zeta_b}
 *     ( \boldmath A \unboldmath - \boldmath B \unboldmath )^2
 *     \right]\f]
 * @param[in] vxiza[ijps] １つ目のCSペア中の \c ijps 番目のPSペアの
 *     \f$\frac{\xi}{\zeta_a}=\frac{\zeta_b}{\zeta_a+\zeta_b}\f$の値
 * @param[in] BA[3] ２つの中心の差ベクトル\f$B-A\f$の各成分の値
 *
 * @param[in] *nklps ２つ目のCSペアに含まれるPSペア数
 * @param[in] veta[klps] ２つ目のCSペア中の \c klps 番目のPSペアの
 *      軌道指数和\f$ \eta = \zeta_c + \zeta_d \f$
 * @param[in] vdkcd[klps] ２つ目のCSペア中の \c klps 番目のPSペアの
 *      結合係数
 *     \f[ K_{cd} = \sqrt2 \pi^{5/4} \frac1{\zeta_c+\zeta_d}
 *     \exp\left[ -\frac{\zeta_c \zeta_d}{\zeta_c + \zeta_d}
 *     ( \boldmath C \unboldmath - \boldmath D \unboldmath )^2
 *     \right]\f]
 * @param[in] vxizc[klps] ２つ目のCSペア中の \c klps 番目のPSペアの
 *     \f$\frac{\xi}{\zeta_c}=\frac{\zeta_d}{\zeta_c+\zeta_d}\f$の値
 * @param[in] DC[3] ２つの中心の差ベクトル\f$D-C\f$の各成分の値
 *
 * @param[out] XXXX[] 計算した縮約２電子積分を格納するための関数。
 *     積分タイプによって変数名は異なる。
 *
 * @ingroup integ-core
 * */
#include <stdio.h>
#include <math.h>

#include "ofmo-def.h"
#include "fmt.h"

#define ZERO 0.e0
#define HALF .5e0

extern double _twoint_inv2_;
extern double _twoint_inv3_;
extern double _twoint_spi2_;

extern double *FMT_fmt_table0;
extern double *FMT_fmt_table1;
extern double *FMT_fmt_table2;
extern double *FMT_fmt_table3;
extern double *FMT_fmt_table4;

extern double FMT_fmt_step_size;
extern double FMT_fmt_inv_step_size;
extern double FMT_pi_div2;
/** １つのCS４重対に対して(ss,ss)タイプの２電子積分を計算する関数
 * @ingroup core-twoint
 * */
void twoint_core_ssss__(
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double SSSS[1] ) {
    int ijps, klps, i;
    double ssss, cssss, zeta, dkab, xiza, eta, xizc, dk;
    double sqrho, rho, PQ2, T;
    // loop unrolling
    int nres;
    double rrho0, rrho1, rrho2, rrho3, dk0, dk1, dk2, dk3;
    double PQ20, PQ21, PQ22, PQ23;
    double sqrho0, sqrho1, sqrho2, sqrho3;
    double xizc0, xizc1, xizc2, xizc3;
    double rho0, rho1, rho2, rho3;
    double cssss0, cssss1, cssss2, cssss3, T0, T1, T2, T3;
    double t00, t01, t02, t03, dt0, dt1, dt2, dt3;
    double F0, F1, F2, F3;
    int it0, it1, it2, it3;
    double PC[3], PQ[3], DC2, PC2, AC2, BA2;
    double BADC2, DCAC2, BAAC2, PCDC2;
    nres= ((*nklps)& 0x03);	// = (*nklps) % 4
    BA2   = BA[0]*BA[0] + BA[1]*BA[1] + BA[2]*BA[2];
    DC2   = DC[0]*DC[0] + DC[1]*DC[1] + DC[2]*DC[2];
    AC2   = AC[0]*AC[0] + AC[1]*AC[1] + AC[2]*AC[2];
    BADC2 = 2.e0*( BA[0]*DC[0] + BA[1]*DC[1] + BA[2]*DC[2] );
    BAAC2 = 2.e0*( BA[0]*AC[0] + BA[1]*AC[1] + BA[2]*AC[2] );
    DCAC2 = 2.e0*( DC[0]*AC[0] + DC[1]*AC[1] + DC[2]*AC[2] );
    SSSS[0] = ZERO;
    for ( ijps=0; ijps<(*nijps); ijps++ ) {
	zeta = vzeta[ijps];
	dkab = vdkab[ijps];
	xiza = vxiza[ijps];

	PC2   = AC2 + xiza*( BAAC2 + xiza*BA2 );
	PCDC2 = DCAC2 + xiza*BADC2;

	for ( klps=0; klps<nres; klps++ ) {
	//for ( klps=0; klps<(*nklps); klps++ ) {
	    eta   = veta[klps];
	    dk    = dkab * vdkcd[klps];
	    xizc  = vxizc[klps];
	    sqrho = sqrt(1.e0/(zeta+eta));
	    PQ2   = PC2 - xizc*( PCDC2 - xizc*DC2 );
	    rho   = sqrho*sqrho;
	    cssss = sqrho * dk;
	    T     = rho * PQ2;
	    {
		int pos;
		double dT;
		if ( T < 36e0 ) {
		    it0 = (int)(0.5e0 + T * FMT_fmt_inv_step_size);
		    dT  = it0 * FMT_fmt_step_size - T;
		    pos = it0 * 4;
		    ssss = (( FMT_fmt_table0[pos+3]   * dT
			    + FMT_fmt_table0[pos+2] ) * dT
			    + FMT_fmt_table0[pos+1] ) * dT
			    + FMT_fmt_table0[pos+0];
		} else {
		    ssss = _twoint_spi2_ * sqrt( 1.e0/(T+T) );
		}
	    }
	    SSSS[0] += cssss * ssss;
	}
	for ( klps=nres; klps<(*nklps); klps+=4 ) {
	    sqrho0 = sqrt( 1.e0 / (zeta+veta[klps+0]) );
	    sqrho1 = sqrt( 1.e0 / (zeta+veta[klps+1]) );
	    sqrho2 = sqrt( 1.e0 / (zeta+veta[klps+2]) );
	    sqrho3 = sqrt( 1.e0 / (zeta+veta[klps+3]) );
	    dk0  = dkab * vdkcd[klps+0];
	    dk1  = dkab * vdkcd[klps+1];
	    dk2  = dkab * vdkcd[klps+2];
	    dk3  = dkab * vdkcd[klps+3];
	    xizc0 = vxizc[klps+0];
	    xizc1 = vxizc[klps+1];
	    xizc2 = vxizc[klps+2];
	    xizc3 = vxizc[klps+3];
	    PQ20   = PC2 - xizc0*( PCDC2 - xizc0*DC2 );
	    PQ21   = PC2 - xizc1*( PCDC2 - xizc1*DC2 );
	    PQ22   = PC2 - xizc2*( PCDC2 - xizc2*DC2 );
	    PQ23   = PC2 - xizc3*( PCDC2 - xizc3*DC2 );
	    rho0   = sqrho0 * sqrho0;
	    rho1   = sqrho1 * sqrho1;
	    rho2   = sqrho2 * sqrho2;
	    rho3   = sqrho3 * sqrho3;
            cssss0 = sqrho0 * dk0;
            cssss1 = sqrho1 * dk1;
            cssss2 = sqrho2 * dk2;
            cssss3 = sqrho3 * dk3;
            T0     = rho0 * PQ20;
            T1     = rho1 * PQ21;
            T2     = rho2 * PQ22;
            T3     = rho3 * PQ23;
	    {
		int pos0, pos1, pos2, pos3;
		double dT0, dT1, dT2, dT3;
		// 0
		if ( T0 < 36e0 ) {
		    it0  = (int)(0.5e0 + T0 * FMT_fmt_inv_step_size);
		    dT0  = it0 * FMT_fmt_step_size - T0;
		    pos0 = it0 * 4;
		    F0 = (( FMT_fmt_table0[pos0+3]   * dT0
			  + FMT_fmt_table0[pos0+2] ) * dT0
			  + FMT_fmt_table0[pos0+1] ) * dT0
			  + FMT_fmt_table0[pos0+0];
		} else {
		    F0 = _twoint_spi2_ * sqrt( 1.e0/(T0+T0) );
		}
		// 1
		if ( T1 < 36e0 ) {
		    it1  = (int)(0.5e0 + T1 * FMT_fmt_inv_step_size);
		    dT1  = it1 * FMT_fmt_step_size - T1;
		    pos1 = it1 * 4;
		    F1 = (( FMT_fmt_table0[pos1+3]   * dT1
			  + FMT_fmt_table0[pos1+2] ) * dT1
			  + FMT_fmt_table0[pos1+1] ) * dT1
			  + FMT_fmt_table0[pos1+0];
		} else {
		    F1 = _twoint_spi2_ * sqrt( 1.e0/(T1+T1) );
		}
		// 2
		if ( T2 < 36e0 ) {
		    it2  = (int)(0.5e0 + T2 * FMT_fmt_inv_step_size);
		    dT2  = it2 * FMT_fmt_step_size - T2;
		    pos2 = it2 * 4;
		    F2 = (( FMT_fmt_table0[pos2+3]   * dT2
			  + FMT_fmt_table0[pos2+2] ) * dT2
			  + FMT_fmt_table0[pos2+1] ) * dT2
			  + FMT_fmt_table0[pos2+0];
		} else {
		    F2 = _twoint_spi2_ * sqrt( 1.e0/(T2+T2) );
		}
		// 3
		if ( T3 < 36e0 ) {
		    it3  = (int)(0.5e0 + T3 * FMT_fmt_inv_step_size);
		    dT3  = it3 * FMT_fmt_step_size - T3;
		    pos3 = it3 * 4;
		    F3 = (( FMT_fmt_table0[pos3+3]   * dT3
			  + FMT_fmt_table0[pos3+2] ) * dT3
			  + FMT_fmt_table0[pos3+1] ) * dT3
			  + FMT_fmt_table0[pos3+0];
		} else {
		    F3 = _twoint_spi2_ * sqrt( 1.e0/(T3+T3) );
		}
	    }	// calc. error functions
	    SSSS[0] += cssss0*F0 + cssss1*F1 + cssss2*F2 + cssss3*F3;
	}	// klps
    }	// ijps
}

/** １つのCS４重対に対して(ps,ss)タイプの２電子積分を計算する関数
 * @ingroup core-twoint
 * */
void twoint_core_psss__(
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double PSSS[3] ) {
    int ijps, klps, i;
    double ssss[1+1], cssss, zeta, dkab, xiza, eta, xizc, dk;
    double rz;
    double sqrho, rho, PQ2, T, F0, F1;
    double PC[3], PA[3], WP[3], PQ[3];
    //
    double VDK[1296], VSQRHO[1296], VPQ2[1296], VRZ[1296];
    double VF0[1296], VF1[1296], VT[1296], VCSSSS[1296];
    double VPQx[1296], VPQy[1296], VPQz[1296];
    double VTXIZA[1296], VTZETA[1296];
    int nprm, iprm;

    for ( i=0; i<3; i++ ) PSSS[i] = ZERO;
    nprm = 0;
    for ( ijps=0; ijps<(*nijps); ijps++ ) {
	zeta = vzeta[ijps];
	dkab = vdkab[ijps];
	xiza = vxiza[ijps];
	for ( i=0; i<3; i++ ) PC[i] = AC[i] + xiza*BA[i];
	for ( klps=0; klps<(*nklps); klps++ ) {
	    eta  = veta[klps];
	    xizc = vxizc[klps];
	    for ( i=0; i<3; i++ ) PQ[i] = PC[i] - xizc*DC[i];

	    VSQRHO[nprm] = sqrt(1.e0/(zeta+eta));
	    VDK[nprm]    = dkab * vdkcd[klps];
	    VPQ2[nprm]   = PQ[0]*PQ[0] + PQ[1]*PQ[1] + PQ[2]*PQ[2];
	    VPQx[nprm]   = PQ[0];
	    VPQy[nprm]   = PQ[1];
	    VPQz[nprm]   = PQ[2];
	    VTXIZA[nprm] = xiza;
	    VTZETA[nprm] = zeta;
	    nprm++;
	}
    }
    //
    for ( iprm=0; iprm<nprm; iprm++ ) {
	sqrho   = VSQRHO[iprm];
	rho     = sqrho * sqrho;
	VCSSSS[iprm] = sqrho * VDK[iprm];
	VT[iprm]     = rho * VPQ2[iprm];
	VRZ[iprm]    = rho * VTZETA[iprm];
    }
    //
    for ( iprm=0; iprm<nprm; iprm++ ) {
	cssss = VCSSSS[iprm];
	T     = VT[iprm];
	rz    = VRZ[iprm];
	xiza  = VTXIZA[iprm];
	{
	    int it0, pos;
	    double dT, t_inv, st_inv, dT2, dT3;
	    if ( T < 38e0 ) {
		it0 = (int)(0.5e0 + T * FMT_fmt_inv_step_size);
		dT  = it0 * FMT_fmt_step_size - T;
		dT2 = dT * _twoint_inv2_;
		dT3 = dT * _twoint_inv3_;
		pos = it0 * (1+4);
		ssss[0] = cssss*(((FMT_fmt_table1[pos+3]   * dT3
				 + FMT_fmt_table1[pos+2] ) * dT2
				 + FMT_fmt_table1[pos+1] ) * dT
				 + FMT_fmt_table1[pos+0] );
		ssss[1] = cssss*(((FMT_fmt_table1[pos+4]   * dT3
				 + FMT_fmt_table1[pos+3] ) * dT2
				 + FMT_fmt_table1[pos+2] ) * dT
				 + FMT_fmt_table1[pos+1] );
	    } else {
		st_inv = sqrt( 1.e0 / (T+T) );
		t_inv  = st_inv * st_inv;
		ssss[0] = cssss * _twoint_spi2_ * st_inv;
		ssss[1] = t_inv * ssss[0];
	    }
	    VF0[iprm] = xiza * ssss[0];
	    VF1[iprm] = rz   * ssss[1];
	}
    }
    //
    for ( iprm=0; iprm<nprm; iprm++ ) {
	F0 = VF0[iprm];
	F1 = VF1[iprm];
	PQ[0] = VPQx[iprm];
	PQ[1] = VPQy[iprm];
	PQ[2] = VPQz[iprm];
	PSSS[0] += BA[0]*F0 - PQ[0]*F1;
	PSSS[1] += BA[1]*F0 - PQ[1]*F1;
	PSSS[2] += BA[2]*F0 - PQ[2]*F1;
    }
}

/** １つのCS４重対に対して(ps,ps)タイプの２電子積分を計算する関数
 * @ingroup core-twoint
 * */
void twoint_core_psps__(
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double PSPS[3*3] ) {
    int ijps, klps, i;
    double cssss, zeta, dkab, xiza, eta, xizc, dk;
    double rz, ze;
    double sqrho, rho, PQ2, T;
    double PC[3], PA[3], WP[3], PQ[3], QC[3], WQ[3];
    double psss0[3], psss1[3];
    //
    double VDK[1296], VSQRHO[1296], VPQ2[1296];
    double VRHO[1296], VF0[1296], VF1[1296], VF2[1296], VT[1296];
    double VCSSSS[1296], VPQx[1296], VPQy[1296], VPQz[1296];
    double VPAx[1296], VPAy[1296], VPAz[1296], VTZETA[1296], VTETA[1296];
    int nprm, iprm;
    //
    double F0, F1, F2, rho2, tmp;
    nprm = 0;
    for ( ijps=0; ijps<(*nijps); ijps++ ) {
	zeta = vzeta[ijps];
	dkab = vdkab[ijps];
	xiza = vxiza[ijps];
	for (i=0; i<3; i++) PC[i] = AC[i] + xiza*BA[i];
	for (i=0; i<3; i++) PA[i] = xiza * BA[i];
	for ( klps=0; klps<(*nklps); klps++ ) {
	    eta  = veta[klps];
	    xizc = vxizc[klps];
	    for ( i=0; i<3; i++ ) PQ[i] = PC[i] - xizc*DC[i];
	    VDK[nprm]  = dkab * vdkcd[klps];
	    VPQ2[nprm] = PQ[0]*PQ[0] + PQ[1]*PQ[1] + PQ[2]*PQ[2];
	    VSQRHO[nprm] = sqrt( 1.e0/(zeta+eta) );
	    VPQx[nprm] = PQ[0];
	    VPQy[nprm] = PQ[1];
	    VPQz[nprm] = PQ[2];
	    //
	    VPAx[nprm] = PA[0];
	    VPAy[nprm] = PA[1];
	    VPAz[nprm] = PA[2];
	    VTZETA[nprm] = zeta;
	    VTETA[nprm]  = eta;
	    nprm++;
	}
    }
    //
    for ( iprm=0; iprm<nprm; iprm++ ) {
	sqrho = VSQRHO[iprm];
	rho   = sqrho * sqrho;
	VCSSSS[iprm] = sqrho * VDK[iprm];
	VT[iprm]     = rho * VPQ2[iprm];
	VRHO[iprm]   = rho;
    }
    //
    for ( iprm=0; iprm<nprm; iprm++ ) {
	cssss = VCSSSS[iprm];
	T     = VT[iprm];
	rho   = VRHO[iprm];
	rho2  = rho * rho;
	{
	    int it0, pos;
	    double dT, t_inv, st_inv, dT2, dT3;
	    if ( T < 40e0 ) {
		it0 = (int)(0.5e0 + T * FMT_fmt_inv_step_size);
		dT  = it0 * FMT_fmt_step_size - T;
		dT2 = dT * _twoint_inv2_;
		dT3 = dT * _twoint_inv3_;
		pos = it0 * (2+4);
		F0 = cssss*(((FMT_fmt_table2[pos+3]   * dT3
			    + FMT_fmt_table2[pos+2] ) * dT2
			    + FMT_fmt_table2[pos+1] ) * dT
			    + FMT_fmt_table2[pos+0] );
		F1 = cssss*(((FMT_fmt_table2[pos+4]   * dT3
			    + FMT_fmt_table2[pos+3] ) * dT2
			    + FMT_fmt_table2[pos+2] ) * dT
			    + FMT_fmt_table2[pos+1] );
		F2 = cssss*(((FMT_fmt_table2[pos+5]   * dT3
			    + FMT_fmt_table2[pos+4] ) * dT2
			    + FMT_fmt_table2[pos+3] ) * dT
			    + FMT_fmt_table2[pos+2] );
	    } else {
		st_inv  = sqrt( 1.e0 / (T+T) );
		t_inv   = st_inv * st_inv;
		F0 = cssss * _twoint_spi2_ * st_inv;
		F1 =        t_inv * F0;
		F2 = 3.e0 * t_inv * F1;
	    }
	    VF0[iprm] =        F0;
	    VF1[iprm] = rho  * F1;
	    VF2[iprm] = rho2 * F2;
	}
    }
    for ( i=0; i<3*3; i++ ) PSPS[i] = ZERO;
    tmp = 0.e0;
    for ( iprm=0; iprm<nprm; iprm++ ) {
	zeta  = VTZETA[iprm];
	eta   = VTETA[iprm];
	PQ[0] = VPQx[iprm];
	PQ[1] = VPQy[iprm];
	PQ[2] = VPQz[iprm];
	PA[0] = VPAx[iprm];
	PA[1] = VPAy[iprm];
	PA[2] = VPAz[iprm];
	F0 = VF0[iprm];
	F1 = VF1[iprm];
	F2 = VF2[iprm];
	//
	ze = zeta*eta;
	for (i=0; i<3; i++) WP[i] = zeta*PQ[i];
	for (i=0; i<3; i++) WQ[i] =  eta*PQ[i];
	for (i=0; i<3; i++) QC[i] = AC[i]+PA[i]-PQ[i];
	// psss(m=0,1)
	for (i=0; i<3; i++) psss0[i] = PA[i]*F0 - WP[i]*F1;
	for (i=0; i<3; i++) psss1[i] = PA[i]*F1 - WP[i]*F2;
	tmp += ze*F1;
	// psps
	PSPS[0*3+0] += QC[0]*psss0[0]+WQ[0]*psss1[0];
	PSPS[0*3+1] += QC[1]*psss0[0]+WQ[1]*psss1[0];
	PSPS[0*3+2] += QC[2]*psss0[0]+WQ[2]*psss1[0];
	PSPS[1*3+0] += QC[0]*psss0[1]+WQ[0]*psss1[1];
	PSPS[1*3+1] += QC[1]*psss0[1]+WQ[1]*psss1[1];
	PSPS[1*3+2] += QC[2]*psss0[1]+WQ[2]*psss1[1];
	PSPS[2*3+0] += QC[0]*psss0[2]+WQ[0]*psss1[2];
	PSPS[2*3+1] += QC[1]*psss0[2]+WQ[1]*psss1[2];
	PSPS[2*3+2] += QC[2]*psss0[2]+WQ[2]*psss1[2];
    }
    PSPS[0*3+0] += 0.5e0*tmp;
    PSPS[1*3+1] += 0.5e0*tmp;
    PSPS[2*3+2] += 0.5e0*tmp;
}

/** １つのCS４重対に対して(pp,ss)タイプの２電子積分を計算する関数
 * @ingroup core-twoint
 * */
void twoint_core_ppss__(
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double PPSS[3*3] ) {
    int ijps, klps, i, m;
    double cssss, zeta, dkab, xiza, eta, xizc, dk;
    double rz, zeta2, tmp;
    double sqrho, rho, PQ2, T;
    double PC[3], PA[3], WP[3], QP[3];
    double ssss[1+2], psss[1+1][3];
    double PSSS[3], DSSS[6];
    for ( i=0; i<3; i++ ) PSSS[i] = ZERO;
    for ( i=0; i<6; i++ ) DSSS[i] = ZERO;
    for ( ijps=0; ijps<(*nijps); ijps++ ) {
	zeta = vzeta[ijps];
	dkab = vdkab[ijps];
	xiza = vxiza[ijps];
	zeta2 = HALF * zeta;
	for ( i=0; i<3; i++ ) {
	    PC[i] = AC[i] + xiza*BA[i];
	    PA[i] = xiza * BA[i];
	}
	for ( klps=0; klps<(*nklps); klps++ ) {
	    eta  = veta[klps];
	    dk   = dkab * vdkcd[klps];
	    xizc = vxizc[klps];
	    PQ2  = ZERO;
	    for ( i=0; i<3; i++ ) {
		QP[i] = xizc*DC[i] - PC[i];
		PQ2  += QP[i]*QP[i];
	    }
	    sqrho = sqrt(1.e0/(zeta+eta));
	    rho   = sqrho*sqrho;
	    rz    = rho * zeta;
	    for ( i=0; i<3; i++ ) WP[i]=rz*QP[i];
	    T     = rho * PQ2;
	    cssss = sqrho * dk;
	    {
		int it0, pos;
		double dT, t_inv, st_inv, dT2, dT3;
		if ( T < 40e0 ) {
		    it0 = (int)(0.5e0 + T * FMT_fmt_inv_step_size);
		    dT  = it0 * FMT_fmt_step_size - T;
		    dT2 = dT * _twoint_inv2_;
		    dT3 = dT * _twoint_inv3_;
		    pos = it0 * (2+4);
		    ssss[0] = cssss*(((FMT_fmt_table2[pos+3] * dT3
				    + FMT_fmt_table2[pos+2] ) * dT2
				    + FMT_fmt_table2[pos+1] ) * dT
				    + FMT_fmt_table2[pos+0] );
		    ssss[1] = cssss*(((FMT_fmt_table2[pos+4] * dT3
				    + FMT_fmt_table2[pos+3] ) * dT2
				    + FMT_fmt_table2[pos+2] ) * dT
				    + FMT_fmt_table2[pos+1] );
		    ssss[2] = cssss*(((FMT_fmt_table2[pos+5] * dT3
				    + FMT_fmt_table2[pos+4] ) * dT2
				    + FMT_fmt_table2[pos+3] ) * dT
				    + FMT_fmt_table2[pos+2] );
		} else {
		    st_inv = sqrt( 1.e0 / (T+T) );
		    t_inv  = st_inv * st_inv;
		    ssss[0] = cssss * _twoint_spi2_ * st_inv;
		    ssss[1] = t_inv * ssss[0];
		    ssss[2] = 3.e0 * t_inv * ssss[1];
		}
	    }
	    //fmt( ssss, 2, T, cssss );
	    // psss
	    for ( m=0; m<=1; m++ ) {
		for ( i=0; i<3; i++ )
		    psss[m][i] = PA[i]*ssss[m]+WP[i]*ssss[m+1];
	    }
	    for (i=0; i<3; i++ ) PSSS[i] += psss[0][i];
	    // dsss
	    tmp = ssss[0] - rz*ssss[1];
	    DSSS[0] += PA[0]*psss[0][0] + WP[0]*psss[1][0] + zeta2*tmp;
	    DSSS[1] += PA[1]*psss[0][1] + WP[1]*psss[1][1] + zeta2*tmp;
	    DSSS[2] += PA[2]*psss[0][2] + WP[2]*psss[1][2] + zeta2*tmp;
	    DSSS[3] += PA[0]*psss[0][1] + WP[0]*psss[1][1];
	    DSSS[4] += PA[1]*psss[0][2] + WP[1]*psss[1][2];
	    DSSS[5] += PA[2]*psss[0][0] + WP[2]*psss[1][0];
	}	// klps
    }	// ijps
    // (P,P|S,S)
    PPSS[0*3+0] = DSSS[0] - BA[0]*PSSS[0];
    PPSS[0*3+1] = DSSS[3] - BA[1]*PSSS[0];
    PPSS[0*3+2] = DSSS[5] - BA[2]*PSSS[0];
    PPSS[1*3+0] = DSSS[3] - BA[0]*PSSS[1];
    PPSS[1*3+1] = DSSS[1] - BA[1]*PSSS[1];
    PPSS[1*3+2] = DSSS[4] - BA[2]*PSSS[1];
    PPSS[2*3+0] = DSSS[5] - BA[0]*PSSS[2];
    PPSS[2*3+1] = DSSS[4] - BA[1]*PSSS[2];
    PPSS[2*3+2] = DSSS[2] - BA[2]*PSSS[2];
}

/** １つのCS４重対に対して(pp,ps)タイプの２電子積分を計算する関数
 * @ingroup core-twoint
 * */
void twoint_core_ppps__(
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double PPPS[3*3*3] ) {
    int ijps, klps, i, m, m1, c0;
    double cssss, zeta, dkab, xiza, eta, xizc, dk;
    double rz, re, ze2, ze22, zeta2, tmp;
    double sqrho, rho, PQ2, T;
    double PC[3], PA[3], WP[3], QP[3], QC[3], WQ[3];
    double ssss[3+1], psss[2+1][3], dsss[1+1][6];
    double PSPS[3*3], DSPS[6*3];
    for ( i=0; i<3*3; i++ ) PSPS[i] = ZERO;
    for ( i=0; i<6*3; i++ ) DSPS[i] = ZERO;
    for ( ijps=0; ijps<(*nijps); ijps++ ) {
	zeta = vzeta[ijps];
	dkab = vdkab[ijps];
	xiza = vxiza[ijps];
	zeta2 = HALF * zeta;
	for ( i=0; i<3; i++ ) {
	    PC[i] = AC[i] + xiza*BA[i];
	    PA[i] = xiza * BA[i];
	}
	for ( klps=0; klps<(*nklps); klps++ ) {
	    eta  = veta[klps];
	    dk   = dkab * vdkcd[klps];
	    xizc = vxizc[klps];
	    PQ2  = ZERO;
	    for ( i=0; i<3; i++ ) {
		QC[i] = xizc*DC[i];
		QP[i] = xizc*DC[i] - PC[i];
		PQ2  += QP[i]*QP[i];
	    }
	    sqrho = sqrt(1.e0/(zeta+eta));
	    rho   = sqrho*sqrho;
	    rz    = rho * zeta;
	    re    = rho * eta;
	    ze2   = re  * zeta2;
	    ze22  = re  * zeta;
	    for ( i=0; i<3; i++ ) {
		WP[i] = rz*QP[i];
		WQ[i] = rz*QP[i] - QP[i];
	    }
	    T     = rho * PQ2;
	    cssss = sqrho * dk;
	    {
		int it0, pos;
		double dT, t_inv, st_inv, dT2, dT3;
		if ( T < 42e0 ) {
		    it0 = (int)(0.5e0 + T * FMT_fmt_inv_step_size);
		    dT  = it0 * FMT_fmt_step_size - T;
		    dT2 = dT * _twoint_inv2_;
		    dT3 = dT * _twoint_inv3_;
		    pos = it0 * (3+4);
		    ssss[0] = cssss*(((FMT_fmt_table3[pos+3] * dT3
				    + FMT_fmt_table3[pos+2] ) * dT2
				    + FMT_fmt_table3[pos+1] ) * dT
				    + FMT_fmt_table3[pos+0] );
		    ssss[1] = cssss*(((FMT_fmt_table3[pos+4] * dT3
				    + FMT_fmt_table3[pos+3] ) * dT2
				    + FMT_fmt_table3[pos+2] ) * dT
				    + FMT_fmt_table3[pos+1] );
		    ssss[2] = cssss*(((FMT_fmt_table3[pos+5] * dT3
				    + FMT_fmt_table3[pos+4] ) * dT2
				    + FMT_fmt_table3[pos+3] ) * dT
				    + FMT_fmt_table3[pos+2] );
		    ssss[3] = cssss*(((FMT_fmt_table3[pos+6] * dT3
				    + FMT_fmt_table3[pos+5] ) * dT2
				    + FMT_fmt_table3[pos+4] ) * dT
				    + FMT_fmt_table3[pos+3] );
		} else {
		    st_inv = sqrt( 1.e0 / (T+T) );
		    t_inv  = st_inv * st_inv;
		    ssss[0] = cssss * _twoint_spi2_ * st_inv;
		    ssss[1] = t_inv * ssss[0];
		    ssss[2] = 3.e0 * t_inv * ssss[1];
		    ssss[3] = 5.e0 * t_inv * ssss[2];
		}
	    }
	    //fmt( ssss, 3, T, cssss );
	    // psss
	    for ( m=0; m<=2; m++ ) {
		m1 = m+1;
		for ( i=0; i<3; i++ )
		    psss[m][i] = PA[i]*ssss[m]+WP[i]*ssss[m1];
	    }
	    // dsss
	    for ( m=0; m<=1; m++) {
		m1 = m+1;
		tmp = ssss[m] - rz*ssss[m1];

		dsss[m][0] = PA[0]*psss[m][0] + WP[0]*psss[m1][0]
		    + zeta2*tmp;
		dsss[m][1] = PA[1]*psss[m][1] + WP[1]*psss[m1][1]
			   + zeta2*tmp;
		dsss[m][2] = PA[2]*psss[m][2] + WP[2]*psss[m1][2]
			   + zeta2*tmp;
		dsss[m][3] = PA[0]*psss[m][1] + WP[0]*psss[m1][1];
		dsss[m][4] = PA[1]*psss[m][2] + WP[1]*psss[m1][2];
		dsss[m][5] = PA[2]*psss[m][0] + WP[2]*psss[m1][0];
	    }

	    // psps
	    PSPS[0*3+0] += QC[0]*psss[0][0] + WQ[0]*psss[1][0]
			+ ze2*ssss[1];
	    PSPS[0*3+1] += QC[1]*psss[0][0] + WQ[1]*psss[1][0];
	    PSPS[0*3+2] += QC[2]*psss[0][0] + WQ[2]*psss[1][0];
	    PSPS[1*3+0] += QC[0]*psss[0][1] + WQ[0]*psss[1][1];
	    PSPS[1*3+1] += QC[1]*psss[0][1] + WQ[1]*psss[1][1]
			+ ze2*ssss[1];
	    PSPS[1*3+2] += QC[2]*psss[0][1] + WQ[2]*psss[1][1];
	    PSPS[2*3+0] += QC[0]*psss[0][2] + WQ[0]*psss[1][2];
	    PSPS[2*3+1] += QC[1]*psss[0][2] + WQ[1]*psss[1][2];
	    PSPS[2*3+2] += QC[2]*psss[0][2] + WQ[2]*psss[1][2]
			+ ze2*ssss[1];
	    // dsps
	    DSPS[0*3+0] += QC[0]*dsss[0][0] + WQ[0]*dsss[1][0]
			+ ze22*psss[1][0];
	    DSPS[0*3+1] += QC[1]*dsss[0][0] + WQ[1]*dsss[1][0];
	    DSPS[0*3+2] += QC[2]*dsss[0][0] + WQ[2]*dsss[1][0];
	    DSPS[1*3+0] += QC[0]*dsss[0][1] + WQ[0]*dsss[1][1];
	    DSPS[1*3+1] += QC[1]*dsss[0][1] + WQ[1]*dsss[1][1]
			+ ze22*psss[1][1];
	    DSPS[1*3+2] += QC[2]*dsss[0][1] + WQ[2]*dsss[1][1];
	    DSPS[2*3+0] += QC[0]*dsss[0][2] + WQ[0]*dsss[1][2];
	    DSPS[2*3+1] += QC[1]*dsss[0][2] + WQ[1]*dsss[1][2];
	    DSPS[2*3+2] += QC[2]*dsss[0][2] + WQ[2]*dsss[1][2]
			+ ze22*psss[1][2];
	    DSPS[3*3+0] += QC[0]*dsss[0][3] + WQ[0]*dsss[1][3]
			+ ze2*psss[1][1];
	    DSPS[3*3+1] += QC[1]*dsss[0][3] + WQ[1]*dsss[1][3]
			+ ze2*psss[1][0];
	    DSPS[3*3+2] += QC[2]*dsss[0][3] + WQ[2]*dsss[1][3];
	    DSPS[4*3+0] += QC[0]*dsss[0][4] + WQ[0]*dsss[1][4];
	    DSPS[4*3+1] += QC[1]*dsss[0][4] + WQ[1]*dsss[1][4]
			+ ze2*psss[1][2];
	    DSPS[4*3+2] += QC[2]*dsss[0][4] + WQ[2]*dsss[1][4]
			+ ze2*psss[1][1];
	    DSPS[5*3+0] += QC[0]*dsss[0][5] + WQ[0]*dsss[1][5]
			+ ze2*psss[1][2];
	    DSPS[5*3+1] += QC[1]*dsss[0][5] + WQ[1]*dsss[1][5];
	    DSPS[5*3+2] += QC[2]*dsss[0][5] + WQ[2]*dsss[1][5]
			+ ze2*psss[1][0];
	}	// klps
    }	// ijps
    // ppps
    for ( c0=0; c0<3; c0++)
        PPPS[0*3*3+0*3+c0] = DSPS[0*3+c0] - BA[0]*PSPS[0*3+c0];
    for ( c0=0; c0<3; c0++)
        PPPS[0*3*3+1*3+c0] = DSPS[3*3+c0] - BA[1]*PSPS[0*3+c0];
    for ( c0=0; c0<3; c0++)
        PPPS[0*3*3+2*3+c0] = DSPS[5*3+c0] - BA[2]*PSPS[0*3+c0];
    for ( c0=0; c0<3; c0++)
        PPPS[1*3*3+0*3+c0] = DSPS[3*3+c0] - BA[0]*PSPS[1*3+c0];
    for ( c0=0; c0<3; c0++)
        PPPS[1*3*3+1*3+c0] = DSPS[1*3+c0] - BA[1]*PSPS[1*3+c0];
    for ( c0=0; c0<3; c0++)
        PPPS[1*3*3+2*3+c0] = DSPS[4*3+c0] - BA[2]*PSPS[1*3+c0];
    for ( c0=0; c0<3; c0++)
        PPPS[2*3*3+0*3+c0] = DSPS[5*3+c0] - BA[0]*PSPS[2*3+c0];
    for ( c0=0; c0<3; c0++)
        PPPS[2*3*3+1*3+c0] = DSPS[4*3+c0] - BA[1]*PSPS[2*3+c0];
    for ( c0=0; c0<3; c0++)
        PPPS[2*3*3+2*3+c0] = DSPS[2*3+c0] - BA[2]*PSPS[2*3+c0];
}

/** １つのCS４重対に対して(pp,pp)タイプの２電子積分を計算する関数
 * @ingroup core-twoint
 * */
void twoint_core_pppp__(
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double PPPP[3*3*3*3] ) {
    int ijps, klps, i, m, m1, c0;
    int ab, ab01, ab10, ab00;
    double cssss, zeta, dkab, xiza, eta, xizc, dk;
    double rz, re, ze2, ze22, zeta2, eta2;
    double tmp, tmp3[3], tmp6[6];
    double sqrho, rho, PQ2, T;
    double PC[3], PA[3], WP[3], QP[3], QC[3], WQ[3];
    double ssss[4+1], psss[3+1][3], dsss[2+1][6];
    double ssps[3], psps[1+1][3*3], dsps[1+1][6*3];
    double PSPS[3*3], DSPS[6*3], PSDS[3*6], DSDS[6*6];
    double PPPS[3*3*3], PPDS[3*3*6];
    for ( i=0; i<3*3; i++ ) PSPS[i] = ZERO;
    for ( i=0; i<6*3; i++ ) DSPS[i] = ZERO;
    for ( i=0; i<3*6; i++ ) PSDS[i] = ZERO;
    for ( i=0; i<6*6; i++ ) DSDS[i] = ZERO;
    for ( ijps=0; ijps<(*nijps); ijps++ ) {
	zeta = vzeta[ijps];
	dkab = vdkab[ijps];
	xiza = vxiza[ijps];
	zeta2 = HALF * zeta;
	for ( i=0; i<3; i++ ) {
	    PC[i] = AC[i] + xiza*BA[i];
	    PA[i] = xiza * BA[i];
	}
	for ( klps=0; klps<(*nklps); klps++ ) {
	    eta  = veta[klps];
	    dk   = dkab * vdkcd[klps];
	    xizc = vxizc[klps];
	    eta2 = eta * HALF;
	    PQ2  = ZERO;
	    for ( i=0; i<3; i++ ) {
		QC[i] = xizc*DC[i];
		QP[i] = xizc*DC[i] - PC[i];
		PQ2  += QP[i]*QP[i];
	    }
	    sqrho = sqrt(1.e0/(zeta+eta));
	    rho   = sqrho*sqrho;
	    rz    = rho * zeta;
	    re    = rho * eta;
	    ze2   = re  * zeta2;
	    ze22  = re  * zeta;
	    for ( i=0; i<3; i++ ) {
		WP[i] = rz*QP[i];
		WQ[i] = rz*QP[i] - QP[i];
	    }
	    T     = rho * PQ2;
	    cssss = sqrho * dk;
	    {
		int it0, pos;
		double dT, t_inv, st_inv, dT2, dT3;
		if ( T < 44e0 ) {
		    it0 = (int)(0.5e0 + T * FMT_fmt_inv_step_size);
		    dT  = it0 * FMT_fmt_step_size - T;
		    dT2 = dT * _twoint_inv2_;
		    dT3 = dT * _twoint_inv3_;
		    pos = it0 * (4+4);
		    ssss[0] = cssss*(((FMT_fmt_table4[pos+3] * dT3
				    + FMT_fmt_table4[pos+2] ) * dT2
				    + FMT_fmt_table4[pos+1] ) * dT
				    + FMT_fmt_table4[pos+0] );
		    ssss[1] = cssss*(((FMT_fmt_table4[pos+4] * dT3
				    + FMT_fmt_table4[pos+3] ) * dT2
				    + FMT_fmt_table4[pos+2] ) * dT
				    + FMT_fmt_table4[pos+1] );
		    ssss[2] = cssss*(((FMT_fmt_table4[pos+5] * dT3
				    + FMT_fmt_table4[pos+4] ) * dT2
				    + FMT_fmt_table4[pos+3] ) * dT
				    + FMT_fmt_table4[pos+2] );
		    ssss[3] = cssss*(((FMT_fmt_table4[pos+6] * dT3
				    + FMT_fmt_table4[pos+5] ) * dT2
				    + FMT_fmt_table4[pos+4] ) * dT
				    + FMT_fmt_table4[pos+3] );
		    ssss[4] = cssss*(((FMT_fmt_table4[pos+7] * dT3
				    + FMT_fmt_table4[pos+6] ) * dT2
				    + FMT_fmt_table4[pos+5] ) * dT
				    + FMT_fmt_table4[pos+4] );
		} else {
		    st_inv = sqrt( 1.e0 / (T+T) );
		    t_inv  = st_inv * st_inv;
		    ssss[0] = cssss * _twoint_spi2_ * st_inv;
		    ssss[1] = t_inv * ssss[0];
		    ssss[2] = 3.e0 * t_inv * ssss[1];
		    ssss[3] = 5.e0 * t_inv * ssss[2];
		    ssss[4] = 7.e0 * t_inv * ssss[3];
		}
	    }
	    //fmt( ssss, 4, T, cssss );
	    // (P,S|S,S)
	    for ( m=0, m1=1; m<=3; m++, m1++) {
		psss[m][0] = PA[0]*ssss[m] + WP[0]*ssss[m1];
		psss[m][1] = PA[1]*ssss[m] + WP[1]*ssss[m1];
		psss[m][2] = PA[2]*ssss[m] + WP[2]*ssss[m1];
	    }	// end m loop
	    // (D,S|S,S)
	    for ( m=0, m1=1; m<=2; m++, m1++) {
		tmp = ssss[m] - rz*ssss[m1];
		dsss[m][0] = PA[0]*psss[m][0] + WP[0]*psss[m1][0]
			   + zeta2*tmp;
		dsss[m][1] = PA[1]*psss[m][1] + WP[1]*psss[m1][1]
			   + zeta2*tmp;
		dsss[m][2] = PA[2]*psss[m][2] + WP[2]*psss[m1][2]
			   + zeta2*tmp;
		dsss[m][3] = PA[0]*psss[m][1] + WP[0]*psss[m1][1];
		dsss[m][4] = PA[1]*psss[m][2] + WP[1]*psss[m1][2];
		dsss[m][5] = PA[2]*psss[m][0] + WP[2]*psss[m1][0];
	    }	// end m loop
	    // (S,S|P,S)
	    ssps[0] = QC[0]*ssss[1] + WQ[0]*ssss[2];
	    ssps[1] = QC[1]*ssss[1] + WQ[1]*ssss[2];
	    ssps[2] = QC[2]*ssss[1] + WQ[2]*ssss[2];
	    // (P,S|P,S)
	    for ( m=0, m1=1; m<=1; m++, m1++) {
		psps[m][0*3+0] = QC[0]*psss[m][0] + WQ[0]*psss[m1][0]
			       + ze2*ssss[m1];
		psps[m][0*3+1] = QC[1]*psss[m][0] + WQ[1]*psss[m1][0];
		psps[m][0*3+2] = QC[2]*psss[m][0] + WQ[2]*psss[m1][0];
		psps[m][1*3+0] = QC[0]*psss[m][1] + WQ[0]*psss[m1][1];
		psps[m][1*3+1] = QC[1]*psss[m][1] + WQ[1]*psss[m1][1]
			       + ze2*ssss[m1];
		psps[m][1*3+2] = QC[2]*psss[m][1] + WQ[2]*psss[m1][1];
		psps[m][2*3+0] = QC[0]*psss[m][2] + WQ[0]*psss[m1][2];
		psps[m][2*3+1] = QC[1]*psss[m][2] + WQ[1]*psss[m1][2];
		psps[m][2*3+2] = QC[2]*psss[m][2] + WQ[2]*psss[m1][2]
			       + ze2*ssss[m1];
	    }	// end m loop
	    // (D,S|P,S)
	    for ( m=0; m<=1; m++) {
		m1 = m+1;
		dsps[m][0*3+0] = QC[0]*dsss[m][0] + WQ[0]*dsss[m1][0]
			       + ze22*psss[m1][0];
		dsps[m][0*3+1] = QC[1]*dsss[m][0] + WQ[1]*dsss[m1][0];
		dsps[m][0*3+2] = QC[2]*dsss[m][0] + WQ[2]*dsss[m1][0];
		dsps[m][1*3+0] = QC[0]*dsss[m][1] + WQ[0]*dsss[m1][1];
		dsps[m][1*3+1] = QC[1]*dsss[m][1] + WQ[1]*dsss[m1][1]
			       + ze22*psss[m1][1];
		dsps[m][1*3+2] = QC[2]*dsss[m][1] + WQ[2]*dsss[m1][1];
		dsps[m][2*3+0] = QC[0]*dsss[m][2] + WQ[0]*dsss[m1][2];
		dsps[m][2*3+1] = QC[1]*dsss[m][2] + WQ[1]*dsss[m1][2];
		dsps[m][2*3+2] = QC[2]*dsss[m][2] + WQ[2]*dsss[m1][2]
			       + ze22*psss[m1][2];
		dsps[m][3*3+0] = QC[0]*dsss[m][3] + WQ[0]*dsss[m1][3]
			       + ze2*psss[m1][1];
		dsps[m][3*3+1] = QC[1]*dsss[m][3] + WQ[1]*dsss[m1][3]
			       + ze2*psss[m1][0];
		dsps[m][3*3+2] = QC[2]*dsss[m][3] + WQ[2]*dsss[m1][3];
		dsps[m][4*3+0] = QC[0]*dsss[m][4] + WQ[0]*dsss[m1][4];
		dsps[m][4*3+1] = QC[1]*dsss[m][4] + WQ[1]*dsss[m1][4]
			       + ze2*psss[m1][2];
		dsps[m][4*3+2] = QC[2]*dsss[m][4] + WQ[2]*dsss[m1][4]
			       + ze2*psss[m1][1];
		dsps[m][5*3+0] = QC[0]*dsss[m][5] + WQ[0]*dsss[m1][5]
			       + ze2*psss[m1][2];
		dsps[m][5*3+1] = QC[1]*dsss[m][5] + WQ[1]*dsss[m1][5];
		dsps[m][5*3+2] = QC[2]*dsss[m][5] + WQ[2]*dsss[m1][5]
			       + ze2*psss[m1][0];
	    }	// end m loop
	    // (P,S|D,S)
	    for ( i=0; i<3; i++) tmp3[i] = psss[0][i] - re*psss[1][i];

	    PSDS[0*6+0] += QC[0]*psps[0][0*3+0] + WQ[0]*psps[1][0*3+0]
			+   eta2*tmp3[0] + ze2*ssps[0];
	    PSDS[0*6+1] += QC[1]*psps[0][0*3+1] + WQ[1]*psps[1][0*3+1]
			+   eta2*tmp3[0];
	    PSDS[0*6+2] += QC[2]*psps[0][0*3+2] + WQ[2]*psps[1][0*3+2]
			+   eta2*tmp3[0];
	    PSDS[0*6+3] += QC[0]*psps[0][0*3+1] + WQ[0]*psps[1][0*3+1]
			+ ze2*ssps[1];
	    PSDS[0*6+4] += QC[1]*psps[0][0*3+2] + WQ[1]*psps[1][0*3+2];
	    PSDS[0*6+5] += QC[2]*psps[0][0*3+0] + WQ[2]*psps[1][0*3+0];
	    PSDS[1*6+0] += QC[0]*psps[0][1*3+0] + WQ[0]*psps[1][1*3+0]
			+   eta2*tmp3[1];
	    PSDS[1*6+1] += QC[1]*psps[0][1*3+1] + WQ[1]*psps[1][1*3+1]
			+   eta2*tmp3[1] + ze2*ssps[1];
	    PSDS[1*6+2] += QC[2]*psps[0][1*3+2] + WQ[2]*psps[1][1*3+2]
			+   eta2*tmp3[1];
	    PSDS[1*6+3] += QC[0]*psps[0][1*3+1] + WQ[0]*psps[1][1*3+1];
	    PSDS[1*6+4] += QC[1]*psps[0][1*3+2] + WQ[1]*psps[1][1*3+2]
			+ ze2*ssps[2];
	    PSDS[1*6+5] += QC[2]*psps[0][1*3+0] + WQ[2]*psps[1][1*3+0];
	    PSDS[2*6+0] += QC[0]*psps[0][2*3+0] + WQ[0]*psps[1][2*3+0]
			+   eta2*tmp3[2];
	    PSDS[2*6+1] += QC[1]*psps[0][2*3+1] + WQ[1]*psps[1][2*3+1]
			+   eta2*tmp3[2];
	    PSDS[2*6+2] += QC[2]*psps[0][2*3+2] + WQ[2]*psps[1][2*3+2]
			+   eta2*tmp3[2] + ze2*ssps[2];
	    PSDS[2*6+3] += QC[0]*psps[0][2*3+1] + WQ[0]*psps[1][2*3+1];
	    PSDS[2*6+4] += QC[1]*psps[0][2*3+2] + WQ[1]*psps[1][2*3+2];
	    PSDS[2*6+5] += QC[2]*psps[0][2*3+0] + WQ[2]*psps[1][2*3+0]
			+ ze2*ssps[0];
	    // (D,S|D,S)
	    for ( i=0; i<6; i++) tmp6[i] = dsss[0][i] - re*dsss[1][i];

	    DSDS[0*6+0] += QC[0]*dsps[0][0*3+0] + WQ[0]*dsps[1][0*3+0]
			+   eta2*tmp6[0] + ze22*psps[1][0*3+0];
	    DSDS[0*6+1] += QC[1]*dsps[0][0*3+1] + WQ[1]*dsps[1][0*3+1]
			+   eta2*tmp6[0];
	    DSDS[0*6+2] += QC[2]*dsps[0][0*3+2] + WQ[2]*dsps[1][0*3+2]
			+   eta2*tmp6[0];
	    DSDS[0*6+3] += QC[0]*dsps[0][0*3+1] + WQ[0]*dsps[1][0*3+1]
			+ ze22*psps[1][0*3+1];
	    DSDS[0*6+4] += QC[1]*dsps[0][0*3+2] + WQ[1]*dsps[1][0*3+2];
	    DSDS[0*6+5] += QC[2]*dsps[0][0*3+0] + WQ[2]*dsps[1][0*3+0];
	    DSDS[1*6+0] += QC[0]*dsps[0][1*3+0] + WQ[0]*dsps[1][1*3+0]
			+   eta2*tmp6[1];
	    DSDS[1*6+1] += QC[1]*dsps[0][1*3+1] + WQ[1]*dsps[1][1*3+1]
			+   eta2*tmp6[1] + ze22*psps[1][1*3+1];
	    DSDS[1*6+2] += QC[2]*dsps[0][1*3+2] + WQ[2]*dsps[1][1*3+2]
			+   eta2*tmp6[1];
	    DSDS[1*6+3] += QC[0]*dsps[0][1*3+1] + WQ[0]*dsps[1][1*3+1];
	    DSDS[1*6+4] += QC[1]*dsps[0][1*3+2] + WQ[1]*dsps[1][1*3+2]
			+ ze22*psps[1][1*3+2];
	    DSDS[1*6+5] += QC[2]*dsps[0][1*3+0] + WQ[2]*dsps[1][1*3+0];
	    DSDS[2*6+0] += QC[0]*dsps[0][2*3+0] + WQ[0]*dsps[1][2*3+0]
			+   eta2*tmp6[2];
	    DSDS[2*6+1] += QC[1]*dsps[0][2*3+1] + WQ[1]*dsps[1][2*3+1]
			+   eta2*tmp6[2];
	    DSDS[2*6+2] += QC[2]*dsps[0][2*3+2] + WQ[2]*dsps[1][2*3+2]
			+   eta2*tmp6[2] + ze22*psps[1][2*3+2];
	    DSDS[2*6+3] += QC[0]*dsps[0][2*3+1] + WQ[0]*dsps[1][2*3+1];
	    DSDS[2*6+4] += QC[1]*dsps[0][2*3+2] + WQ[1]*dsps[1][2*3+2];
	    DSDS[2*6+5] += QC[2]*dsps[0][2*3+0] + WQ[2]*dsps[1][2*3+0]
			+ ze22*psps[1][2*3+0];
	    DSDS[3*6+0] += QC[0]*dsps[0][3*3+0] + WQ[0]*dsps[1][3*3+0]
			+   eta2*tmp6[3] + ze2*psps[1][1*3+0];
	    DSDS[3*6+1] += QC[1]*dsps[0][3*3+1] + WQ[1]*dsps[1][3*3+1]
			+   eta2*tmp6[3] + ze2*psps[1][0*3+1];
	    DSDS[3*6+2] += QC[2]*dsps[0][3*3+2] + WQ[2]*dsps[1][3*3+2]
			+   eta2*tmp6[3];
	    DSDS[3*6+3] += QC[0]*dsps[0][3*3+1] + WQ[0]*dsps[1][3*3+1]
			+ ze2*psps[1][1*3+1];
	    DSDS[3*6+4] += QC[1]*dsps[0][3*3+2] + WQ[1]*dsps[1][3*3+2]
			+ ze2*psps[1][0*3+2];
	    DSDS[3*6+5] += QC[2]*dsps[0][3*3+0] + WQ[2]*dsps[1][3*3+0];
	    DSDS[4*6+0] += QC[0]*dsps[0][4*3+0] + WQ[0]*dsps[1][4*3+0]
			+   eta2*tmp6[4];
	    DSDS[4*6+1] += QC[1]*dsps[0][4*3+1] + WQ[1]*dsps[1][4*3+1]
			+   eta2*tmp6[4] + ze2*psps[1][2*3+1];
	    DSDS[4*6+2] += QC[2]*dsps[0][4*3+2] + WQ[2]*dsps[1][4*3+2]
			+   eta2*tmp6[4] + ze2*psps[1][1*3+2];
	    DSDS[4*6+3] += QC[0]*dsps[0][4*3+1] + WQ[0]*dsps[1][4*3+1];
	    DSDS[4*6+4] += QC[1]*dsps[0][4*3+2] + WQ[1]*dsps[1][4*3+2]
			+ ze2*psps[1][2*3+2];
	    DSDS[4*6+5] += QC[2]*dsps[0][4*3+0] + WQ[2]*dsps[1][4*3+0]
			+ ze2*psps[1][1*3+0];
	    DSDS[5*6+0] += QC[0]*dsps[0][5*3+0] + WQ[0]*dsps[1][5*3+0]
			+   eta2*tmp6[5] + ze2*psps[1][2*3+0];
	    DSDS[5*6+1] += QC[1]*dsps[0][5*3+1] + WQ[1]*dsps[1][5*3+1]
			+   eta2*tmp6[5];
	    DSDS[5*6+2] += QC[2]*dsps[0][5*3+2] + WQ[2]*dsps[1][5*3+2]
			+   eta2*tmp6[5] + ze2*psps[1][0*3+2];
	    DSDS[5*6+3] += QC[0]*dsps[0][5*3+1] + WQ[0]*dsps[1][5*3+1]
			+ ze2*psps[1][2*3+1];
	    DSDS[5*6+4] += QC[1]*dsps[0][5*3+2] + WQ[1]*dsps[1][5*3+2];
	    DSDS[5*6+5] += QC[2]*dsps[0][5*3+0] + WQ[2]*dsps[1][5*3+0]
			+ ze2*psps[1][0*3+0];
	    for ( i=0; i<3*3; i++) PSPS[i] += psps[0][i];
	    for ( i=0; i<6*3; i++) DSPS[i] += dsps[0][i];
	}	// klps
    }	// ijps
    // (P,P|P,S)
    for ( c0=0; c0<3; c0++)
        PPPS[0*3*3+0*3+c0] = DSPS[0*3+c0] - BA[0]*PSPS[0*3+c0];
    for ( c0=0; c0<3; c0++)
        PPPS[0*3*3+1*3+c0] = DSPS[3*3+c0] - BA[1]*PSPS[0*3+c0];
    for ( c0=0; c0<3; c0++)
        PPPS[0*3*3+2*3+c0] = DSPS[5*3+c0] - BA[2]*PSPS[0*3+c0];
    for ( c0=0; c0<3; c0++)
        PPPS[1*3*3+0*3+c0] = DSPS[3*3+c0] - BA[0]*PSPS[1*3+c0];
    for ( c0=0; c0<3; c0++)
        PPPS[1*3*3+1*3+c0] = DSPS[1*3+c0] - BA[1]*PSPS[1*3+c0];
    for ( c0=0; c0<3; c0++)
        PPPS[1*3*3+2*3+c0] = DSPS[4*3+c0] - BA[2]*PSPS[1*3+c0];
    for ( c0=0; c0<3; c0++)
        PPPS[2*3*3+0*3+c0] = DSPS[5*3+c0] - BA[0]*PSPS[2*3+c0];
    for ( c0=0; c0<3; c0++)
        PPPS[2*3*3+1*3+c0] = DSPS[4*3+c0] - BA[1]*PSPS[2*3+c0];
    for ( c0=0; c0<3; c0++)
        PPPS[2*3*3+2*3+c0] = DSPS[2*3+c0] - BA[2]*PSPS[2*3+c0];
    // (P,P|D,S)
    for ( c0=0; c0<6; c0++)
        PPDS[0*3*6+0*6+c0] = DSDS[0*6+c0] - BA[0]*PSDS[0*6+c0];
    for ( c0=0; c0<6; c0++)
        PPDS[0*3*6+1*6+c0] = DSDS[3*6+c0] - BA[1]*PSDS[0*6+c0];
    for ( c0=0; c0<6; c0++)
        PPDS[0*3*6+2*6+c0] = DSDS[5*6+c0] - BA[2]*PSDS[0*6+c0];
    for ( c0=0; c0<6; c0++)
        PPDS[1*3*6+0*6+c0] = DSDS[3*6+c0] - BA[0]*PSDS[1*6+c0];
    for ( c0=0; c0<6; c0++)
        PPDS[1*3*6+1*6+c0] = DSDS[1*6+c0] - BA[1]*PSDS[1*6+c0];
    for ( c0=0; c0<6; c0++)
        PPDS[1*3*6+2*6+c0] = DSDS[4*6+c0] - BA[2]*PSDS[1*6+c0];
    for ( c0=0; c0<6; c0++)
        PPDS[2*3*6+0*6+c0] = DSDS[5*6+c0] - BA[0]*PSDS[2*6+c0];
    for ( c0=0; c0<6; c0++)
        PPDS[2*3*6+1*6+c0] = DSDS[4*6+c0] - BA[1]*PSDS[2*6+c0];
    for ( c0=0; c0<6; c0++)
        PPDS[2*3*6+2*6+c0] = DSDS[2*6+c0] - BA[2]*PSDS[2*6+c0];
    // (a,b|c,d+1) = (a,b|c+1,d) - DC(a,b|c,d)
    // (P,P|P,P)
    for ( ab=0, ab01=0, ab10=0, ab00=0;
            ab<3*3; ab++, ab01 += 3*3, ab10 += 6, ab00 += 3) {
        PPPP[ab01+0*3+0] = PPDS[ab10+0] - DC[0]*PPPS[ab00+0];
        PPPP[ab01+0*3+1] = PPDS[ab10+3] - DC[1]*PPPS[ab00+0];
        PPPP[ab01+0*3+2] = PPDS[ab10+5] - DC[2]*PPPS[ab00+0];
        PPPP[ab01+1*3+0] = PPDS[ab10+3] - DC[0]*PPPS[ab00+1];
        PPPP[ab01+1*3+1] = PPDS[ab10+1] - DC[1]*PPPS[ab00+1];
        PPPP[ab01+1*3+2] = PPDS[ab10+4] - DC[2]*PPPS[ab00+1];
        PPPP[ab01+2*3+0] = PPDS[ab10+5] - DC[0]*PPPS[ab00+2];
        PPPP[ab01+2*3+1] = PPDS[ab10+4] - DC[1]*PPPS[ab00+2];
        PPPP[ab01+2*3+2] = PPDS[ab10+2] - DC[2]*PPPS[ab00+2];
    }
}

/** １つのCS４重対に対して(ds,ss)タイプの２電子積分を計算する関数
 * @ingroup core-twoint
 * */
void twoint_core_dsss__(
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double DSSS[6] ) {
    int ijps, klps, i, m;
    double cssss, zeta, dkab, xiza, eta, xizc, dk;
    double rz, zeta2, tmp;
    double sqrho, rho, PQ2, T;
    double PC[3], PA[3], WP[3], QP[3];
    double ssss[1+2], psss[1+1][3];
    double sqr3;
    sqr3 = sqrt(3.e0);
    for ( i=0; i<6; i++ ) DSSS[i] = ZERO;
    for ( ijps=0; ijps<(*nijps); ijps++ ) {
	zeta = vzeta[ijps];
	dkab = vdkab[ijps];
	xiza = vxiza[ijps];
	zeta2 = HALF * zeta;
	for ( i=0; i<3; i++ ) {
	    PC[i] = AC[i] + xiza*BA[i];
	    PA[i] = xiza * BA[i];
	}
	for ( klps=0; klps<(*nklps); klps++ ) {
	    eta  = veta[klps];
	    dk   = dkab * vdkcd[klps];
	    xizc = vxizc[klps];
	    PQ2  = ZERO;
	    for ( i=0; i<3; i++ ) {
		QP[i] = xizc*DC[i] - PC[i];
		PQ2  += QP[i]*QP[i];
	    }
	    sqrho = sqrt(1.e0/(zeta+eta));
	    rho   = sqrho*sqrho;
	    rz    = rho * zeta;
	    for ( i=0; i<3; i++ ) WP[i]=rz*QP[i];
	    T     = rho * PQ2;
	    cssss = sqrho * dk;
	    {
		int it0, pos;
		double dT, t_inv, st_inv, dT2, dT3;
		if ( T < 40e0 ) {
		    it0 = (int)(0.5e0 + T * FMT_fmt_inv_step_size);
		    dT  = it0 * FMT_fmt_step_size - T;
		    dT2 = dT * _twoint_inv2_;
		    dT3 = dT * _twoint_inv3_;
		    pos = it0 * (2+4);
		    ssss[0] = cssss*(((FMT_fmt_table2[pos+3] * dT3
				    + FMT_fmt_table2[pos+2] ) * dT2
				    + FMT_fmt_table2[pos+1] ) * dT
				    + FMT_fmt_table2[pos+0] );
		    ssss[1] = cssss*(((FMT_fmt_table2[pos+4] * dT3
				    + FMT_fmt_table2[pos+3] ) * dT2
				    + FMT_fmt_table2[pos+2] ) * dT
				    + FMT_fmt_table2[pos+1] );
		    ssss[2] = cssss*(((FMT_fmt_table2[pos+5] * dT3
				    + FMT_fmt_table2[pos+4] ) * dT2
				    + FMT_fmt_table2[pos+3] ) * dT
				    + FMT_fmt_table2[pos+2] );
		} else {
		    st_inv = sqrt( 1.e0 / (T+T) );
		    t_inv  = st_inv * st_inv;
		    ssss[0] = cssss * _twoint_spi2_ * st_inv;
		    ssss[1] = t_inv * ssss[0];
		    ssss[2] = 3.e0 * t_inv * ssss[1];
		}
	    }
	    //fmt( ssss, 2, T, cssss );
	    // psss
	    for ( m=0; m<=1; m++ ) {
		for ( i=0; i<3; i++ )
		    psss[m][i] = PA[i]*ssss[m]+WP[i]*ssss[m+1];
	    }
	    // dsss
	    tmp = ssss[0] - rz*ssss[1];
	    DSSS[0] += PA[0]*psss[0][0] + WP[0]*psss[1][0] + zeta2*tmp;
	    DSSS[1] += PA[1]*psss[0][1] + WP[1]*psss[1][1] + zeta2*tmp;
	    DSSS[2] += PA[2]*psss[0][2] + WP[2]*psss[1][2] + zeta2*tmp;
	    DSSS[3] += PA[0]*psss[0][1] + WP[0]*psss[1][1];
	    DSSS[4] += PA[1]*psss[0][2] + WP[1]*psss[1][2];
	    DSSS[5] += PA[2]*psss[0][0] + WP[2]*psss[1][0];
	}	// klps
    }	// ijps
    for (i=3; i<6; i++ ) DSSS[i] *= sqr3;
}

/** １つのCS４重対に対して(ds,ps)タイプの２電子積分を計算する関数
 * @ingroup core-twoint
 * */
void twoint_core_dsps__(
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double DSPS[6*3] ) {
    int ijps, klps, i, m, m1;
    double cssss, zeta, dkab, xiza, eta, xizc, dk;
    double rz, re, ze2, ze22, zeta2, tmp;
    double sqrho, rho, PQ2, T;
    double PC[3], PA[3], WP[3], QP[3], QC[3], WQ[3];
    double ssss[3+1], psss[2+1][3], dsss[1+1][6];
    double sqr3;
    sqr3 = sqrt(3.e0);
    for ( i=0; i<6*3; i++ ) DSPS[i] = ZERO;
    for ( ijps=0; ijps<(*nijps); ijps++ ) {
	zeta = vzeta[ijps];
	dkab = vdkab[ijps];
	xiza = vxiza[ijps];
	zeta2 = HALF * zeta;
	for ( i=0; i<3; i++ ) {
	    PC[i] = AC[i] + xiza*BA[i];
	    PA[i] = xiza * BA[i];
	}
	for ( klps=0; klps<(*nklps); klps++ ) {
	    eta  = veta[klps];
	    dk   = dkab * vdkcd[klps];
	    xizc = vxizc[klps];
	    PQ2  = ZERO;
	    for ( i=0; i<3; i++ ) {
		QC[i] = xizc*DC[i];
		QP[i] = xizc*DC[i] - PC[i];
		PQ2  += QP[i]*QP[i];
	    }
	    sqrho = sqrt(1.e0/(zeta+eta));
	    rho   = sqrho*sqrho;
	    rz    = rho * zeta;
	    re    = rho * eta;
	    ze2   = re  * zeta2;
	    ze22  = re  * zeta;
	    for ( i=0; i<3; i++ ) {
		WP[i] = rz*QP[i];
		WQ[i] = rz*QP[i] - QP[i];
	    }
	    T     = rho * PQ2;
	    cssss = sqrho * dk;
	    {
		int it0, pos;
		double dT, t_inv, st_inv, dT2, dT3;
		if ( T < 42e0 ) {
		    it0 = (int)(0.5e0 + T * FMT_fmt_inv_step_size);
		    dT  = it0 * FMT_fmt_step_size - T;
		    dT2 = dT * _twoint_inv2_;
		    dT3 = dT * _twoint_inv3_;
		    pos = it0 * (3+4);
		    ssss[0] = cssss*(((FMT_fmt_table3[pos+3] * dT3
				    + FMT_fmt_table3[pos+2] ) * dT2
				    + FMT_fmt_table3[pos+1] ) * dT
				    + FMT_fmt_table3[pos+0] );
		    ssss[1] = cssss*(((FMT_fmt_table3[pos+4] * dT3
				    + FMT_fmt_table3[pos+3] ) * dT2
				    + FMT_fmt_table3[pos+2] ) * dT
				    + FMT_fmt_table3[pos+1] );
		    ssss[2] = cssss*(((FMT_fmt_table3[pos+5] * dT3
				    + FMT_fmt_table3[pos+4] ) * dT2
				    + FMT_fmt_table3[pos+3] ) * dT
				    + FMT_fmt_table3[pos+2] );
		    ssss[3] = cssss*(((FMT_fmt_table3[pos+6] * dT3
				    + FMT_fmt_table3[pos+5] ) * dT2
				    + FMT_fmt_table3[pos+4] ) * dT
				    + FMT_fmt_table3[pos+3] );
		} else {
		    st_inv = sqrt( 1.e0 / (T+T) );
		    t_inv  = st_inv * st_inv;
		    ssss[0] = cssss * _twoint_spi2_ * st_inv;
		    ssss[1] = t_inv * ssss[0];
		    ssss[2] = 3.e0 * t_inv * ssss[1];
		    ssss[3] = 5.e0 * t_inv * ssss[2];
		}
	    }
	    //fmt( ssss, 3, T, cssss );
	    // psss
	    for ( m=0; m<=2; m++ ) {
		m1 = m+1;
		for ( i=0; i<3; i++ )
		    psss[m][i] = PA[i]*ssss[m]+WP[i]*ssss[m1];
	    }
	    // dsss
	    for ( m=0; m<=1; m++) {
		m1 = m+1;
		tmp = ssss[m] - rz*ssss[m1];

		dsss[m][0] = PA[0]*psss[m][0] + WP[0]*psss[m1][0]
		    + zeta2*tmp;
		dsss[m][1] = PA[1]*psss[m][1] + WP[1]*psss[m1][1]
			   + zeta2*tmp;
		dsss[m][2] = PA[2]*psss[m][2] + WP[2]*psss[m1][2]
			   + zeta2*tmp;
		dsss[m][3] = PA[0]*psss[m][1] + WP[0]*psss[m1][1];
		dsss[m][4] = PA[1]*psss[m][2] + WP[1]*psss[m1][2];
		dsss[m][5] = PA[2]*psss[m][0] + WP[2]*psss[m1][0];
	    }
	    // dsps
	    DSPS[0*3+0] += QC[0]*dsss[0][0] + WQ[0]*dsss[1][0]
			+ ze22*psss[1][0];
	    DSPS[0*3+1] += QC[1]*dsss[0][0] + WQ[1]*dsss[1][0];
	    DSPS[0*3+2] += QC[2]*dsss[0][0] + WQ[2]*dsss[1][0];
	    DSPS[1*3+0] += QC[0]*dsss[0][1] + WQ[0]*dsss[1][1];
	    DSPS[1*3+1] += QC[1]*dsss[0][1] + WQ[1]*dsss[1][1]
			+ ze22*psss[1][1];
	    DSPS[1*3+2] += QC[2]*dsss[0][1] + WQ[2]*dsss[1][1];
	    DSPS[2*3+0] += QC[0]*dsss[0][2] + WQ[0]*dsss[1][2];
	    DSPS[2*3+1] += QC[1]*dsss[0][2] + WQ[1]*dsss[1][2];
	    DSPS[2*3+2] += QC[2]*dsss[0][2] + WQ[2]*dsss[1][2]
			+ ze22*psss[1][2];
	    DSPS[3*3+0] += QC[0]*dsss[0][3] + WQ[0]*dsss[1][3]
			+ ze2*psss[1][1];
	    DSPS[3*3+1] += QC[1]*dsss[0][3] + WQ[1]*dsss[1][3]
			+ ze2*psss[1][0];
	    DSPS[3*3+2] += QC[2]*dsss[0][3] + WQ[2]*dsss[1][3];
	    DSPS[4*3+0] += QC[0]*dsss[0][4] + WQ[0]*dsss[1][4];
	    DSPS[4*3+1] += QC[1]*dsss[0][4] + WQ[1]*dsss[1][4]
			+ ze2*psss[1][2];
	    DSPS[4*3+2] += QC[2]*dsss[0][4] + WQ[2]*dsss[1][4]
			+ ze2*psss[1][1];
	    DSPS[5*3+0] += QC[0]*dsss[0][5] + WQ[0]*dsss[1][5]
			+ ze2*psss[1][2];
	    DSPS[5*3+1] += QC[1]*dsss[0][5] + WQ[1]*dsss[1][5];
	    DSPS[5*3+2] += QC[2]*dsss[0][5] + WQ[2]*dsss[1][5]
			+ ze2*psss[1][0];
	}	// klps
    }	// ijps
    for ( i=3*3; i<6*3; i++) DSPS[i] *= sqr3;
}

/** １つのCS４重対に対して(ds,pp)タイプの２電子積分を計算する関数
 * @ingroup core-twoint
 * */
void twoint_core_dspp__(
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double DSPP[6*3*3] ) {
    int ijps, klps, i, m, m1;
    int ab, ab01, ab10, ab00;
    double cssss, zeta, dkab, xiza, eta, xizc, dk;
    double rz, re, ze2, ze22, zeta2, eta2;
    double tmp, tmp6[6];
    double sqrho, rho, PQ2, T;
    double PC[3], PA[3], WP[3], QP[3], QC[3], WQ[3];
    double ssss[4+1], psss[3+1][3], dsss[2+1][6];
    double psps[3*3], dsps[1+1][6*3];
    double DSPS[6*3], DSDS[6*6];
    double sqr3;
    sqr3 = sqrt(3.e0);
    for ( i=0; i<6*3; i++ ) DSPS[i] = ZERO;
    for ( i=0; i<6*6; i++ ) DSDS[i] = ZERO;
    for ( ijps=0; ijps<(*nijps); ijps++ ) {
	zeta = vzeta[ijps];
	dkab = vdkab[ijps];
	xiza = vxiza[ijps];
	zeta2 = HALF * zeta;
	for ( i=0; i<3; i++ ) {
	    PC[i] = AC[i] + xiza*BA[i];
	    PA[i] = xiza * BA[i];
	}
	for ( klps=0; klps<(*nklps); klps++ ) {
	    eta  = veta[klps];
	    dk   = dkab * vdkcd[klps];
	    xizc = vxizc[klps];
	    eta2 = HALF * eta;
	    PQ2  = ZERO;
	    for ( i=0; i<3; i++ ) {
		QC[i] = xizc*DC[i];
		QP[i] = xizc*DC[i] - PC[i];
		PQ2  += QP[i]*QP[i];
	    }
	    sqrho = sqrt(1.e0/(zeta+eta));
	    rho   = sqrho*sqrho;
	    rz    = rho * zeta;
	    re    = rho * eta;
	    ze2   = rz  * eta2;
	    ze22  = rz  * eta;
	    for ( i=0; i<3; i++ ) {
		WP[i] = rz*QP[i];
		WQ[i] = rz*QP[i] - QP[i];
	    }
	    T     = rho * PQ2;
	    cssss = sqrho * dk;
	    {
		int it0, pos;
		double dT, t_inv, st_inv, dT2, dT3;
		if ( T < 44e0 ) {
		    it0 = (int)(0.5e0 + T * FMT_fmt_inv_step_size);
		    dT  = it0 * FMT_fmt_step_size - T;
		    dT2 = dT * _twoint_inv2_;
		    dT3 = dT * _twoint_inv3_;
		    pos = it0 * (4+4);
		    ssss[0] = cssss*(((FMT_fmt_table4[pos+3] * dT3
				    + FMT_fmt_table4[pos+2] ) * dT2
				    + FMT_fmt_table4[pos+1] ) * dT
				    + FMT_fmt_table4[pos+0] );
		    ssss[1] = cssss*(((FMT_fmt_table4[pos+4] * dT3
				    + FMT_fmt_table4[pos+3] ) * dT2
				    + FMT_fmt_table4[pos+2] ) * dT
				    + FMT_fmt_table4[pos+1] );
		    ssss[2] = cssss*(((FMT_fmt_table4[pos+5] * dT3
				    + FMT_fmt_table4[pos+4] ) * dT2
				    + FMT_fmt_table4[pos+3] ) * dT
				    + FMT_fmt_table4[pos+2] );
		    ssss[3] = cssss*(((FMT_fmt_table4[pos+6] * dT3
				    + FMT_fmt_table4[pos+5] ) * dT2
				    + FMT_fmt_table4[pos+4] ) * dT
				    + FMT_fmt_table4[pos+3] );
		    ssss[4] = cssss*(((FMT_fmt_table4[pos+7] * dT3
				    + FMT_fmt_table4[pos+6] ) * dT2
				    + FMT_fmt_table4[pos+5] ) * dT
				    + FMT_fmt_table4[pos+4] );
		} else {
		    st_inv = sqrt( 1.e0 / (T+T) );
		    t_inv  = st_inv * st_inv;
		    ssss[0] = cssss * _twoint_spi2_ * st_inv;
		    ssss[1] = t_inv * ssss[0];
		    ssss[2] = 3.e0 * t_inv * ssss[1];
		    ssss[3] = 5.e0 * t_inv * ssss[2];
		    ssss[4] = 7.e0 * t_inv * ssss[3];
		}
	    }
	    //fmt( ssss, 4, T, cssss );
	    // (P,S|S,S)
	    for ( m=0; m<=3; m++) {
		m1 = m+1;
		psss[m][0] = PA[0]*ssss[m] + WP[0]*ssss[m1];
		psss[m][1] = PA[1]*ssss[m] + WP[1]*ssss[m1];
		psss[m][2] = PA[2]*ssss[m] + WP[2]*ssss[m1];
	    }	// end m loop
	    // (D,S|S,S)
	    for ( m=0; m<=2; m++) {
		m1 = m+1;
		tmp = ssss[m] - rz*ssss[m1];

		dsss[m][0] = PA[0]*psss[m][0] + WP[0]*psss[m1][0]
			   + zeta2*tmp;
		dsss[m][1] = PA[1]*psss[m][1] + WP[1]*psss[m1][1]
			   + zeta2*tmp;
		dsss[m][2] = PA[2]*psss[m][2] + WP[2]*psss[m1][2]
			   + zeta2*tmp;
		dsss[m][3] = PA[0]*psss[m][1] + WP[0]*psss[m1][1];
		dsss[m][4] = PA[1]*psss[m][2] + WP[1]*psss[m1][2];
		dsss[m][5] = PA[2]*psss[m][0] + WP[2]*psss[m1][0];
	    }	// end m loop
	    // (P,S|P,S)
	    psps[0*3+0] = QC[0]*psss[1][0] + WQ[0]*psss[2][0]
			+ ze2*ssss[2];
	    psps[0*3+1] = QC[1]*psss[1][0] + WQ[1]*psss[2][0];
	    psps[0*3+2] = QC[2]*psss[1][0] + WQ[2]*psss[2][0];
	    psps[1*3+0] = QC[0]*psss[1][1] + WQ[0]*psss[2][1];
	    psps[1*3+1] = QC[1]*psss[1][1] + WQ[1]*psss[2][1]
			+ ze2*ssss[2];
	    psps[1*3+2] = QC[2]*psss[1][1] + WQ[2]*psss[2][1];
	    psps[2*3+0] = QC[0]*psss[1][2] + WQ[0]*psss[2][2];
	    psps[2*3+1] = QC[1]*psss[1][2] + WQ[1]*psss[2][2];
	    psps[2*3+2] = QC[2]*psss[1][2] + WQ[2]*psss[2][2]
			+ ze2*ssss[2];
	    // (D,S|P,S)
	    for ( m=0; m<=1; m++) {
		m1 = m+1;
		dsps[m][0*3+0] = QC[0]*dsss[m][0] + WQ[0]*dsss[m1][0]
			       + ze22*psss[m1][0];
		dsps[m][0*3+1] = QC[1]*dsss[m][0] + WQ[1]*dsss[m1][0];
		dsps[m][0*3+2] = QC[2]*dsss[m][0] + WQ[2]*dsss[m1][0];
		dsps[m][1*3+0] = QC[0]*dsss[m][1] + WQ[0]*dsss[m1][1];
		dsps[m][1*3+1] = QC[1]*dsss[m][1] + WQ[1]*dsss[m1][1]
			       + ze22*psss[m1][1];
		dsps[m][1*3+2] = QC[2]*dsss[m][1] + WQ[2]*dsss[m1][1];
		dsps[m][2*3+0] = QC[0]*dsss[m][2] + WQ[0]*dsss[m1][2];
		dsps[m][2*3+1] = QC[1]*dsss[m][2] + WQ[1]*dsss[m1][2];
		dsps[m][2*3+2] = QC[2]*dsss[m][2] + WQ[2]*dsss[m1][2]
			       + ze22*psss[m1][2];
		dsps[m][3*3+0] = QC[0]*dsss[m][3] + WQ[0]*dsss[m1][3]
			       + ze2*psss[m1][1];
		dsps[m][3*3+1] = QC[1]*dsss[m][3] + WQ[1]*dsss[m1][3]
			       + ze2*psss[m1][0];
		dsps[m][3*3+2] = QC[2]*dsss[m][3] + WQ[2]*dsss[m1][3];
		dsps[m][4*3+0] = QC[0]*dsss[m][4] + WQ[0]*dsss[m1][4];
		dsps[m][4*3+1] = QC[1]*dsss[m][4] + WQ[1]*dsss[m1][4]
			       + ze2*psss[m1][2];
		dsps[m][4*3+2] = QC[2]*dsss[m][4] + WQ[2]*dsss[m1][4]
			       + ze2*psss[m1][1];
		dsps[m][5*3+0] = QC[0]*dsss[m][5] + WQ[0]*dsss[m1][5]
			       + ze2*psss[m1][2];
		dsps[m][5*3+1] = QC[1]*dsss[m][5] + WQ[1]*dsss[m1][5];
		dsps[m][5*3+2] = QC[2]*dsss[m][5] + WQ[2]*dsss[m1][5]
			       + ze2*psss[m1][0];
	    }	// end m loop
	    for ( i=0; i<6*3; i++ ) DSPS[i] += dsps[0][i];
	    // (D,S|D,S)
	    for ( i=0; i<6; i++) tmp6[i] = dsss[0][i] - re*dsss[1][i];

	    DSDS[0*6+0] += QC[0]*dsps[0][0*3+0] + WQ[0]*dsps[1][0*3+0]
			+   eta2*tmp6[0] + ze22*psps[0*3+0];
	    DSDS[0*6+1] += QC[1]*dsps[0][0*3+1] + WQ[1]*dsps[1][0*3+1]
			+   eta2*tmp6[0];
	    DSDS[0*6+2] += QC[2]*dsps[0][0*3+2] + WQ[2]*dsps[1][0*3+2]
			+   eta2*tmp6[0];
	    DSDS[0*6+3] += QC[0]*dsps[0][0*3+1] + WQ[0]*dsps[1][0*3+1]
			+ ze22*psps[0*3+1];
	    DSDS[0*6+4] += QC[1]*dsps[0][0*3+2] + WQ[1]*dsps[1][0*3+2];
	    DSDS[0*6+5] += QC[2]*dsps[0][0*3+0] + WQ[2]*dsps[1][0*3+0];
	    DSDS[1*6+0] += QC[0]*dsps[0][1*3+0] + WQ[0]*dsps[1][1*3+0]
			+   eta2*tmp6[1];
	    DSDS[1*6+1] += QC[1]*dsps[0][1*3+1] + WQ[1]*dsps[1][1*3+1]
			+   eta2*tmp6[1] + ze22*psps[1*3+1];
	    DSDS[1*6+2] += QC[2]*dsps[0][1*3+2] + WQ[2]*dsps[1][1*3+2]
			+   eta2*tmp6[1];
	    DSDS[1*6+3] += QC[0]*dsps[0][1*3+1] + WQ[0]*dsps[1][1*3+1];
	    DSDS[1*6+4] += QC[1]*dsps[0][1*3+2] + WQ[1]*dsps[1][1*3+2]
			+ ze22*psps[1*3+2];
	    DSDS[1*6+5] += QC[2]*dsps[0][1*3+0] + WQ[2]*dsps[1][1*3+0];
	    DSDS[2*6+0] += QC[0]*dsps[0][2*3+0] + WQ[0]*dsps[1][2*3+0]
			+   eta2*tmp6[2];
	    DSDS[2*6+1] += QC[1]*dsps[0][2*3+1] + WQ[1]*dsps[1][2*3+1]
			+   eta2*tmp6[2];
	    DSDS[2*6+2] += QC[2]*dsps[0][2*3+2] + WQ[2]*dsps[1][2*3+2]
			+   eta2*tmp6[2] + ze22*psps[2*3+2];
	    DSDS[2*6+3] += QC[0]*dsps[0][2*3+1] + WQ[0]*dsps[1][2*3+1];
	    DSDS[2*6+4] += QC[1]*dsps[0][2*3+2] + WQ[1]*dsps[1][2*3+2];
	    DSDS[2*6+5] += QC[2]*dsps[0][2*3+0] + WQ[2]*dsps[1][2*3+0]
			+ ze22*psps[2*3+0];
	    DSDS[3*6+0] += QC[0]*dsps[0][3*3+0] + WQ[0]*dsps[1][3*3+0]
			+   eta2*tmp6[3] + ze2*psps[1*3+0];
	    DSDS[3*6+1] += QC[1]*dsps[0][3*3+1] + WQ[1]*dsps[1][3*3+1]
			+   eta2*tmp6[3] + ze2*psps[0*3+1];
	    DSDS[3*6+2] += QC[2]*dsps[0][3*3+2] + WQ[2]*dsps[1][3*3+2]
			+   eta2*tmp6[3];
	    DSDS[3*6+3] += QC[0]*dsps[0][3*3+1] + WQ[0]*dsps[1][3*3+1]
			+ ze2*psps[1*3+1];
	    DSDS[3*6+4] += QC[1]*dsps[0][3*3+2] + WQ[1]*dsps[1][3*3+2]
			+ ze2*psps[0*3+2];
	    DSDS[3*6+5] += QC[2]*dsps[0][3*3+0] + WQ[2]*dsps[1][3*3+0];
	    DSDS[4*6+0] += QC[0]*dsps[0][4*3+0] + WQ[0]*dsps[1][4*3+0]
			+   eta2*tmp6[4];
	    DSDS[4*6+1] += QC[1]*dsps[0][4*3+1] + WQ[1]*dsps[1][4*3+1]
			+   eta2*tmp6[4] + ze2*psps[2*3+1];
	    DSDS[4*6+2] += QC[2]*dsps[0][4*3+2] + WQ[2]*dsps[1][4*3+2]
			+   eta2*tmp6[4] + ze2*psps[1*3+2];
	    DSDS[4*6+3] += QC[0]*dsps[0][4*3+1] + WQ[0]*dsps[1][4*3+1];
	    DSDS[4*6+4] += QC[1]*dsps[0][4*3+2] + WQ[1]*dsps[1][4*3+2]
			+ ze2*psps[2*3+2];
	    DSDS[4*6+5] += QC[2]*dsps[0][4*3+0] + WQ[2]*dsps[1][4*3+0]
			+ ze2*psps[1*3+0];
	    DSDS[5*6+0] += QC[0]*dsps[0][5*3+0] + WQ[0]*dsps[1][5*3+0]
			+   eta2*tmp6[5] + ze2*psps[2*3+0];
	    DSDS[5*6+1] += QC[1]*dsps[0][5*3+1] + WQ[1]*dsps[1][5*3+1]
			+   eta2*tmp6[5];
	    DSDS[5*6+2] += QC[2]*dsps[0][5*3+2] + WQ[2]*dsps[1][5*3+2]
			+   eta2*tmp6[5] + ze2*psps[0*3+2];
	    DSDS[5*6+3] += QC[0]*dsps[0][5*3+1] + WQ[0]*dsps[1][5*3+1]
			+ ze2*psps[2*3+1];
	    DSDS[5*6+4] += QC[1]*dsps[0][5*3+2] + WQ[1]*dsps[1][5*3+2];
	    DSDS[5*6+5] += QC[2]*dsps[0][5*3+0] + WQ[2]*dsps[1][5*3+0]
			+ ze2*psps[0*3+0];
	}	// klps
    }	// ijps
    // (D,S|P,P)
    for ( ab=0, ab01=0, ab10=0, ab00=0;
            ab<6; ab++, ab01 += 3*3, ab10 += 6, ab00 += 3) {
        DSPP[ab01+0*3+0] = DSDS[ab10+0] - DC[0]*DSPS[ab00+0];
        DSPP[ab01+0*3+1] = DSDS[ab10+3] - DC[1]*DSPS[ab00+0];
        DSPP[ab01+0*3+2] = DSDS[ab10+5] - DC[2]*DSPS[ab00+0];
        DSPP[ab01+1*3+0] = DSDS[ab10+3] - DC[0]*DSPS[ab00+1];
        DSPP[ab01+1*3+1] = DSDS[ab10+1] - DC[1]*DSPS[ab00+1];
        DSPP[ab01+1*3+2] = DSDS[ab10+4] - DC[2]*DSPS[ab00+1];
        DSPP[ab01+2*3+0] = DSDS[ab10+5] - DC[0]*DSPS[ab00+2];
        DSPP[ab01+2*3+1] = DSDS[ab10+4] - DC[1]*DSPS[ab00+2];
        DSPP[ab01+2*3+2] = DSDS[ab10+2] - DC[2]*DSPS[ab00+2];
    }
    for ( i=3*3*3; i<6*3*3; i++ ) DSPP[i] *= sqr3;
}
