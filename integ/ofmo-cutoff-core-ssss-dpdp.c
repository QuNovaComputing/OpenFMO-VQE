/**
 * @file ofmo-cutoff-core-ssss-dpdp.c
 * １つのSchwarz積分計算関数群（１）
 *
 * １つのSchwarz積分を計算する関数群のうち、(ss,ss)?(dp,dp)タイプの
 * 計算のための関数が定義してあるファイル。
 * (ss,ss)タイプ以外では、計算された複数のSchwarz積分の最大値を返す。
 * */

/**
 * @defgroup core-cutoff Schwarz積分を計算する関数群
 * １つのSchwarz積分を計算する関数
 *
 *
 * すべての関数で同じ引数を持っている。以下に、引数の説明を記す。
 *
 * @param[in] npps PSペアの数
 * @param[in] vzeta[ipps] PSペア番号 \c ipps のPSペアの軌道指数和
 *     \f$ \zeta = \zeta_a + \zeta_b \f$
 * @param[in] vdkps[ipps] PSペア番号 \c ipps のPSペアの結合定数
 *     \f[ K_{ab} = \sqrt2 \pi^{5/4} \frac1{\zeta_a+\zeta_b}
 *     \exp\left[ -\frac{\zeta_a \zeta_b}{\zeta_a + \zeta_b}
 *     ( \boldmath A \unboldmath - \boldmath B \unboldmath )^2
 *     \right]\f]
 * @param[in] vxiza[ipps] PSペア番号 \c ipps のPSペアの
 *     \f$ \frac{\xi}{\zeta_a} = \frac{\zeta_b}{\zeta_a+\zeta_b} \f$
 * @param[in] BA[3] \f$ \boldmath B-A \f$の各成分
 * @param[in] AB2 \f$ (A-B)^2 \f$の値
 *
 * @return Schwarz積分。(ss,ss)タイプ以外では、
 *     計算された複数のSchwarz積分の最大値
 *
 * @ingroup integ-core
 * */
#include <stdio.h>
#include <math.h>

#define ZERO 0.e0
#define HALF .5e0

#define EPS_PS_PAIR 1.e-32

#include "ofmo-def.h"
#include "fmt.h"

/** (ss,ss)タイプのSchwarz積分計算関数
 * @ingroup core-cutoff
 * */
double schwarz_core_ssss_(
	const int npps, const double vzeta[], const double vdkps[],
	const double vxiza[], const double BA[3], const double AB2 ) {
    int ijps, klps;
    double zeta, Kab, xiza;
    double eta,  Kcd, xizc;
    double Kabcd, sqrho, dxi, T, cssss, ssss;
    double SSSS;
    SSSS = ZERO;
    for ( ijps=0; ijps<npps; ijps++ ) {
	zeta = vzeta[ijps];
	Kab  = vdkps[ijps];
	xiza = vxiza[ijps];
	for ( klps=0; klps<npps; klps++ ) {
	    eta  = vzeta[klps];
	    Kcd  = vdkps[klps];
	    xizc = vxiza[klps];

	    Kabcd = Kab * Kcd;
	    sqrho = 1.e0 / sqrt( zeta + eta );
	    dxi   = xizc - xiza;
	    T     = sqrho * sqrho * dxi * dxi * AB2;
	    cssss = sqrho * Kabcd;
	    fmt(&ssss, 0, T, cssss);

	    SSSS += ssss;
	}
    }
    return fabs(SSSS);
}

/** (ps,ps)タイプのSchwarz積分計算関数
 * @ingroup core-cutoff
 * */
double schwarz_core_psps_(
	const int npps, const double vzeta[], const double vdkps[],
	const double vxiza[], const double BA[3], const double AB2 ) {
    int i;
    int ijps, klps;
    double zeta, Kab, xiza;
    double eta,  Kcd, xizc;
    double Kabcd, sqrho, dxi, rho, T, cssss, rz, ze2;
    double PA[3], QC[3], QPi, WP[3], WQ[3];
    double ssss[2+1], psss[(1+1)*3], PSPS[3];
    double maxval;

    for ( i=0; i<3; i++ ) PSPS[i] = ZERO;

    for ( ijps=0; ijps<npps; ijps++ ) {
	zeta = vzeta[ijps];
	Kab  = vdkps[ijps];
	xiza = vxiza[ijps];
	for ( i=0; i<3; i++ ) PA[i] = xiza*BA[i];

	for ( klps=0; klps<npps; klps++ ) {
	    eta  = vzeta[klps];
	    Kcd  = vdkps[klps];
	    xizc = vxiza[klps];

	    Kabcd = Kab * Kcd;
	    sqrho = 1.e0 / sqrt( zeta + eta );
	    dxi   = xizc - xiza;
	    rho   = sqrho * sqrho;
	    T     = rho * dxi * dxi * AB2;
	    cssss = sqrho * Kabcd;
	    rz    = zeta  * rho;
	    ze2   = HALF * eta * rz;
	    for ( i=0; i<3; i++ ) {
		QC[i] = xizc*BA[i];
		QPi   = QC[i] - PA[i];
		WP[i] = rz * QPi;
		WQ[i] = rz * QPi - QPi;
	    }
	    fmt(ssss, 2, T, cssss);

	    // (ps,ss) m=0,1
	    // m=0
	    psss[0*3+0]=PA[0]*ssss[0]+WP[0]*ssss[1];
	    psss[0*3+1]=PA[1]*ssss[0]+WP[1]*ssss[1];
	    psss[0*3+2]=PA[2]*ssss[0]+WP[2]*ssss[1];
	    // m=1
	    psss[1*3+0]=PA[0]*ssss[1]+WP[0]*ssss[2];
	    psss[1*3+1]=PA[1]*ssss[1]+WP[1]*ssss[2];
	    psss[1*3+2]=PA[2]*ssss[1]+WP[2]*ssss[2];
	    // (ps,ps) m=0 ( only diagonal term )
	    PSPS[0]+=QC[0]*psss[0*3+0]+WQ[0]*psss[1*3+0]
		+ze2*ssss[1];
	    PSPS[1]+=QC[1]*psss[0*3+1]+WQ[1]*psss[1*3+1]
		+ze2*ssss[1];
	    PSPS[2]+=QC[2]*psss[0*3+2]+WQ[2]*psss[1*3+2]
		+ze2*ssss[1];
	}
    }
    for ( i=0; i<3; i++ ) PSPS[i]=fabs(PSPS[i]);
    maxval = PSPS[0];
    for ( i=1; i<3; i++ ) if ( PSPS[i]>maxval ) maxval = PSPS[i];

    return maxval;
}

/** (pp,pp)タイプのSchwarz積分計算関数
 * @ingroup core-cutoff
 * */
double schwarz_core_pppp_(
	const int npps, const double vzeta[], const double vdkps[],
	const double vxiza[], const double BA[3], const double AB2 ) {
    int i, m;
    int ijps, klps;
    double zeta, Kab, xiza, zeta2;
    double  eta, Kcd, xizc,   eta2;
    double Kabcd, sqrho, dxi, rho, T, cssss, rz, re, ze2, ze22, tmp1;
    double PA[3], QC[3], QPi, WP[3], WQ[3];
    double ssss[4+1], psss[(3+1)*3], dsss[(2+1)*6];
    double ssps[3], psps[(1+1)*3], dsps[(1+1)*9];
    double PSDS[9], PSPS[3], DSPS[9], DSDS[6];
    double PPPS[9], PPDS[9], PPPP[9];
    double maxval;

    double *psss_m, *psss_m1, *psps_m, *psps_m1,  *dsss_m, *dsss_m1;
    double *dsps_m, *dsps_m1;

    for (i=0; i<3; i++) PSPS[i] = ZERO;
    for (i=0; i<9; i++) DSPS[i] = ZERO;
    for (i=0; i<9; i++) PSDS[i] = ZERO;
    for (i=0; i<6; i++) DSDS[i] = ZERO;

    for ( ijps=0; ijps<npps; ijps++ ) {
	zeta  = vzeta[ijps];
	Kab   = vdkps[ijps];
	xiza  = vxiza[ijps];
	zeta2 = HALF*zeta;
	for ( i=0; i<3; i++ ) PA[i] = xiza*BA[i];
	for ( klps=0; klps<npps; klps++ ) {
	    eta  = vzeta[klps];
	    Kcd  = vdkps[klps];
	    xizc = vxiza[klps];
	    eta2 = HALF*eta;

	    Kabcd = Kab * Kcd;
	    sqrho = 1.e0 / sqrt( zeta + eta );
	    dxi   = xizc - xiza;
	    rho   = sqrho * sqrho;
	    T     = rho * dxi * dxi * AB2;
	    cssss = sqrho * Kabcd;
	    rz    = zeta * rho;
	    re    = eta  * rho;
	    ze2   = eta2 * rz;
	    ze22  = eta  * rz;
	    for ( i=0; i<3; i++ ) {
		QC[i] = xizc*BA[i];
		QPi   = QC[i] - PA[i];
		WP[i] = rz * QPi;
		WQ[i] = rz * QPi - QPi;
	    }
	    fmt(ssss, 4, T, cssss);
	    // (ps,ss) m=0-3
	    //for ( m=0, psss_m=&psss[0]; m<=3; m++, psss_m+=3 ) {
	    for ( m=0, psss_m=psss; m<=3; m++, psss_m+=3 ) {
		psss_m[0]=PA[0]*ssss[m]+WP[0]*ssss[m+1];
		psss_m[1]=PA[1]*ssss[m]+WP[1]*ssss[m+1];
		psss_m[2]=PA[2]*ssss[m]+WP[2]*ssss[m+1];
	    }
	    // (ds,ss) m=0-2
	    for ( m=0, dsss_m=dsss, psss_m=psss, psss_m1=psss+3;
		    m<=2;
		    m++, dsss_m+=6, psss_m+=3, psss_m1+=3) {
		tmp1 = ssss[m]-rz*ssss[m+1];
		dsss_m[0]=PA[0]*psss_m[0]+WP[0]*psss_m1[0]+zeta2*tmp1;
		dsss_m[1]=PA[1]*psss_m[1]+WP[1]*psss_m1[1]+zeta2*tmp1;
		dsss_m[2]=PA[2]*psss_m[2]+WP[2]*psss_m1[2]+zeta2*tmp1;
		dsss_m[3]=PA[0]*psss_m[1]+WP[0]*psss_m1[1];
		dsss_m[4]=PA[1]*psss_m[2]+WP[1]*psss_m1[2];
		dsss_m[5]=PA[2]*psss_m[0]+WP[2]*psss_m1[0];
	    }
	    // (ss,ps) m=1
	    m=1;
	    ssps[0]=QC[0]*ssss[m]+WQ[0]*ssss[m+1];
	    ssps[1]=QC[1]*ssss[m]+WQ[1]*ssss[m+1];
	    ssps[2]=QC[2]*ssss[m]+WQ[2]*ssss[m+1];
	    // (ps,ps) m=0,1 ( only for diagonal term )
	    for ( m=0, psps_m=psps, psss_m=psss, psss_m1=psss+3;
		    m<=1;
		    m++, psps_m+=3, psss_m+=3, psss_m1+=3 ) {
		psps_m[0]=QC[0]*psss_m[0]+WQ[0]*psss_m1[0]+ze2*ssss[m+1];
		psps_m[1]=QC[1]*psss_m[1]+WQ[1]*psss_m1[1]+ze2*ssss[m+1];
		psps_m[2]=QC[2]*psss_m[2]+WQ[2]*psss_m1[2]+ze2*ssss[m+1];
	    }
	    for ( i=0; i<3; i++ ) PSPS[i] += psps[i];
	    // (ds,ps) m=0,1 ( only for diagonal term )
	    for ( m=0, dsps_m=dsps, dsss_m=dsss, dsss_m1=dsss+6,
		    psss_m1=psss+3;
		    m<=1;
		    m++, dsps_m+=9, dsss_m+=6, dsss_m1+=6, psss_m1+=3 ) {
		// (xx,x), (xy,x), (zx,x)
		dsps_m[0]=QC[0]*dsss_m[0]+WQ[0]*dsss_m1[0]+ze22*psss_m1[0];
		dsps_m[1]=QC[0]*dsss_m[3]+WQ[0]*dsss_m1[3]+ze22*psss_m1[1];
		dsps_m[2]=QC[0]*dsss_m[5]+WQ[0]*dsss_m1[5]+ze22*psss_m1[2];
		// (xy,y), (yy,y), (yz,y)
		dsps_m[3]=QC[1]*dsss_m[3]+WQ[1]*dsss_m1[3]+ze22*psss_m1[0];
		dsps_m[4]=QC[1]*dsss_m[1]+WQ[1]*dsss_m1[1]+ze22*psss_m1[1];
		dsps_m[5]=QC[1]*dsss_m[4]+WQ[1]*dsss_m1[4]+ze22*psss_m1[2];
		// (zx,z), (yz,z), (zz,z)
		dsps_m[6]=QC[2]*dsss_m[5]+WQ[2]*dsss_m1[5]+ze22*psss_m1[0];
		dsps_m[7]=QC[2]*dsss_m[4]+WQ[2]*dsss_m1[4]+ze22*psss_m1[1];
		dsps_m[8]=QC[2]*dsss_m[2]+WQ[2]*dsss_m1[2]+ze22*psss_m1[2];
	    }
	    for (i=0; i<9; i++ ) DSPS[i] += dsps[i];
	    // (ps,ds) m=0
	    psps_m  = psps;
	    psps_m1 = psps+3;
	    psss_m  = psss;
	    psss_m1 = psss+3;
		// (x,xx), (x,xy), (x,zx)
	    PSDS[0]+=QC[0]*psps_m[0]+WQ[0]*psps_m1[0]
		+eta2*(psss_m[0]-re*psss_m1[0])+ze2*ssps[0];
	    PSDS[1]+=QC[1]*psps_m[0]+WQ[1]*psps_m1[0];
	    PSDS[2]+=QC[2]*psps_m[0]+WQ[2]*psps_m1[0];
		// (y,xy), (y,yy), (y,yz)
	    PSDS[3]+=QC[0]*psps_m[1]+WQ[0]*psps_m1[1];
	    PSDS[4]+=QC[1]*psps_m[1]+WQ[1]*psps_m1[1]
		+eta2*(psss_m[1]-re*psss_m1[1])+ze2*ssps[1];
	    PSDS[5]+=QC[2]*psps_m[1]+WQ[2]*psps_m1[1];
		// (z,zx), (z,yz), (z,zz)
	    PSDS[6]+=QC[0]*psps_m[2]+WQ[0]*psps_m1[2];
	    PSDS[7]+=QC[1]*psps_m[2]+WQ[1]*psps_m1[2];
	    PSDS[8]+=QC[2]*psps_m[2]+WQ[2]*psps_m1[2]
		+eta2*(psss_m[2]-re*psss_m1[2])+ze2*ssps[2];
	    // (ds,ds) m=0
	    dsps_m  = dsps;
	    dsps_m1 = dsps+9;
	    dsss_m  = dsss;
	    dsss_m1 = dsss+6;
	    psps_m1 = psps+3;
	    // (xx,xx), (yy,yy), (zz,zz)
	    DSDS[0]+=QC[0]*dsps_m[0]+WQ[0]*dsps_m1[0]
		+eta2*(dsss_m[0]-re*dsss_m1[0])+ze22*psps_m1[0];
	    DSDS[1]+=QC[1]*dsps_m[4]+WQ[1]*dsps_m1[4]
		+eta2*(dsss_m[1]-re*dsss_m1[1])+ze22*psps_m1[1];
	    DSDS[2]+=QC[2]*dsps_m[8]+WQ[2]*dsps_m1[8]
		+eta2*(dsss_m[2]-re*dsss_m1[2])+ze22*psps_m1[2];
	    // (xy,xy), (yz,yz), (zx,zx)
	    DSDS[3]+=QC[0]*dsps_m[3]+WQ[0]*dsps_m1[3]+ze2*psps_m1[1];
	    DSDS[4]+=QC[1]*dsps_m[7]+WQ[1]*dsps_m1[7]+ze2*psps_m1[2];
	    DSDS[5]+=QC[2]*dsps_m[2]+WQ[2]*dsps_m1[2]+ze2*psps_m1[0];
	}	// klps
    }		// ijps
    // HRR
    // (PP,PS)
	// (PxP*|PxS)
    PPPS[0]=DSPS[0]-BA[0]*PSPS[0];
    PPPS[1]=DSPS[1]-BA[1]*PSPS[0];
    PPPS[2]=DSPS[2]-BA[2]*PSPS[0];
	// (PyP*|PyS)
    PPPS[3]=DSPS[3]-BA[0]*PSPS[1];
    PPPS[4]=DSPS[4]-BA[1]*PSPS[1];
    PPPS[5]=DSPS[5]-BA[2]*PSPS[1];
	// (PyP*|PyS)
    PPPS[6]=DSPS[6]-BA[0]*PSPS[2];
    PPPS[7]=DSPS[7]-BA[1]*PSPS[2];
    PPPS[8]=DSPS[8]-BA[2]*PSPS[2];
    // (PP,DS)
	// (PxP*|DS)
    PPDS[0]=DSDS[0]-BA[0]*PSDS[0];
    PPDS[1]=DSDS[3]-BA[1]*PSDS[1];
    PPDS[2]=DSDS[5]-BA[2]*PSDS[2];
	// (PyP*|DS)
    PPDS[3]=DSDS[3]-BA[0]*PSDS[3];
    PPDS[4]=DSDS[1]-BA[1]*PSDS[4];
    PPDS[5]=DSDS[4]-BA[2]*PSDS[5];
	// (PzP*|DS)
    PPDS[6]=DSDS[5]-BA[0]*PSDS[6];
    PPDS[7]=DSDS[4]-BA[1]*PSDS[7];
    PPDS[8]=DSDS[2]-BA[2]*PSDS[8];
    // (PP|PP)
    PPPP[0]=PPDS[0]-BA[0]*PPPS[0];
    PPPP[1]=PPDS[1]-BA[1]*PPPS[1];
    PPPP[2]=PPDS[2]-BA[2]*PPPS[2];
    PPPP[3]=PPDS[3]-BA[0]*PPPS[3];
    PPPP[4]=PPDS[4]-BA[1]*PPPS[4];
    PPPP[5]=PPDS[5]-BA[2]*PPPS[5];
    PPPP[6]=PPDS[6]-BA[0]*PPPS[6];
    PPPP[7]=PPDS[7]-BA[1]*PPPS[7];
    PPPP[8]=PPDS[8]-BA[2]*PPPS[8];

    for ( i=0; i<9; i++ ) PPPP[i]=fabs(PPPP[i]);
    maxval = PPPP[0];
    for ( i=1; i<9; i++ ) if ( PPPP[i]>maxval ) maxval = PPPP[i];

    return maxval;
}

/** (ds,ds)タイプのSchwarz積分計算関数
 * @ingroup core-cutoff
 * */
double schwarz_core_dsds_(
	const int npps, const double vzeta[], const double vdkps[],
	const double vxiza[], const double BA[3], const double AB2 ) {
    int i, m;
    int ijps, klps;
    double zeta, Kab, xiza, zeta2;
    double  eta, Kcd, xizc,   eta2;
    double Kabcd, sqrho, dxi, rho, T, cssss, rz, re, ze2, ze22, tmp1;
    double PA[3], QC[3], QPi, WP[3], WQ[3];
    double ssss[4+1], psss[(3+1)*3], dsss[(2+1)*6];
    double psps[3], dsps[(1+1)*6];
    double DSDS[6];
    double maxval;

    double *psss_m, *psss_m1, *dsss_m, *dsss_m1, *dsps_m, *dsps_m1;

    for (i=0; i<6; i++) DSDS[i] = ZERO;

    for ( ijps=0; ijps<npps; ijps++ ) {
	zeta  = vzeta[ijps];
	Kab   = vdkps[ijps];
	xiza  = vxiza[ijps];
	zeta2 = HALF*zeta;
	for ( i=0; i<3; i++ ) PA[i] = xiza*BA[i];
	for ( klps=0; klps<npps; klps++ ) {
	    eta  = vzeta[klps];
	    Kcd  = vdkps[klps];
	    xizc = vxiza[klps];
	    eta2 = HALF*eta;

	    Kabcd = Kab * Kcd;
	    sqrho = 1.e0 / sqrt( zeta + eta );
	    dxi   = xizc - xiza;
	    rho   = sqrho * sqrho;
	    T     = rho * dxi * dxi * AB2;
	    cssss = sqrho * Kabcd;
	    rz    = zeta * rho;
	    re    = eta  * rho;
	    ze2   = eta2 * rz;
	    ze22  = eta  * rz;
	    for ( i=0; i<3; i++ ) {
		QC[i] = xizc*BA[i];
		QPi   = QC[i] - PA[i];
		WP[i] = rz * QPi;
		WQ[i] = rz * QPi - QPi;
	    }
	    fmt(ssss, 4, T, cssss);
	    // (ps,ss) m=0-3
	    for ( m=0, psss_m=psss; m<=3; m++, psss_m+=3 ) {
		psss_m[0]=PA[0]*ssss[m]+WP[0]*ssss[m+1];
		psss_m[1]=PA[1]*ssss[m]+WP[1]*ssss[m+1];
		psss_m[2]=PA[2]*ssss[m]+WP[2]*ssss[m+1];
	    }
	    // (ds,ss) m=0-2
	    for ( m=0, dsss_m=dsss, psss_m=psss, psss_m1=psss+3;
		    m<=2;
		    m++, dsss_m+=6, psss_m+=3, psss_m1+=3) {
		tmp1 = ssss[m]-rz*ssss[m+1];
		dsss_m[0]=PA[0]*psss_m[0]+WP[0]*psss_m1[0]+zeta2*tmp1;
		dsss_m[1]=PA[1]*psss_m[1]+WP[1]*psss_m1[1]+zeta2*tmp1;
		dsss_m[2]=PA[2]*psss_m[2]+WP[2]*psss_m1[2]+zeta2*tmp1;
		dsss_m[3]=PA[0]*psss_m[1]+WP[0]*psss_m1[1];
		dsss_m[4]=PA[1]*psss_m[2]+WP[1]*psss_m1[2];
		dsss_m[5]=PA[2]*psss_m[0]+WP[2]*psss_m1[0];
	    }
	    // (ps,ps) m=1 ( only for diagonal term )
	    m = 1;
	    psps[0]=QC[0]*psss[3*m+0]+WQ[0]*psss[3*(m+1)+0]+ze2*ssss[m+1];
	    psps[1]=QC[1]*psss[3*m+1]+WQ[1]*psss[3*(m+1)+1]+ze2*ssss[m+1];
	    psps[2]=QC[2]*psss[3*m+2]+WQ[2]*psss[3*(m+1)+2]+ze2*ssss[m+1];
	    // (ds,ps) m=0,1 ( only for diagonal term )
	    for ( m=0, dsps_m=dsps, dsss_m=dsss, dsss_m1=dsss+6,
		    psss_m1=psss+3;
		    m<=1;
		    m++, dsps_m+=6, dsss_m+=6, dsss_m1+=6, psss_m1+=3 ) {
		// (xx|x), (yy|y), (zz|z)
		dsps_m[0]=QC[0]*dsss_m[0]+WQ[0]*dsss_m1[0]+ze22*psss_m1[0];
		dsps_m[1]=QC[1]*dsss_m[1]+WQ[1]*dsss_m1[1]+ze22*psss_m1[1];
		dsps_m[2]=QC[2]*dsss_m[2]+WQ[2]*dsss_m1[2]+ze22*psss_m1[2];
		// (xy,x), (yz,y), (zx,z)
		dsps_m[3]=QC[0]*dsss_m[3]+WQ[0]*dsss_m1[3]+ze22*psss_m1[1];
		dsps_m[4]=QC[1]*dsss_m[4]+WQ[1]*dsss_m1[4]+ze22*psss_m1[2];
		dsps_m[5]=QC[2]*dsss_m[5]+WQ[2]*dsss_m1[5]+ze22*psss_m1[0];
	    }
	    // (ds,ds) m=0
	    dsps_m  = dsps;
	    dsps_m1 = dsps+6;
	    dsss_m  = dsss;
	    dsss_m1 = dsss+6;
	    // (xx,xx), (yy,yy), (zz,zz)
	    DSDS[0]+=QC[0]*dsps_m[0]+WQ[0]*dsps_m1[0]
		+eta2*(dsss_m[0]-re*dsss_m1[0])+ze22*psps[0];
	    DSDS[1]+=QC[1]*dsps_m[1]+WQ[1]*dsps_m1[1]
		+eta2*(dsss_m[1]-re*dsss_m1[1])+ze22*psps[1];
	    DSDS[2]+=QC[2]*dsps_m[2]+WQ[2]*dsps_m1[2]
		+eta2*(dsss_m[2]-re*dsss_m1[2])+ze22*psps[2];
	    // (xy,xy), (yz,yz), (zx,zx)
	    DSDS[3]+=QC[1]*dsps_m[3]+WQ[1]*dsps_m1[3]+ze2*psps[1];
	    DSDS[4]+=QC[2]*dsps_m[4]+WQ[2]*dsps_m1[4]+ze2*psps[2];
	    DSDS[5]+=QC[0]*dsps_m[5]+WQ[0]*dsps_m1[5]+ze2*psps[0];
	}	// klps
    }		// ijps
    for ( i=3; i<6; i++ ) DSDS[i]*=3.e0;

    for ( i=0; i<6; i++ ) DSDS[i]=fabs(DSDS[i]);
    maxval = DSDS[0];
    for ( i=1; i<6; i++ ) if ( DSDS[i]>maxval ) maxval = DSDS[i];

    return maxval;
}

/** (dp,dp)タイプのSchwarz積分計算関数
 * @ingroup core-cutoff
 * */
double schwarz_core_dpdp_(
	const int npps, const double vzeta[], const double vdkps[],
	const double vxiza[], const double BA[3], const double AB2 ) {
    int i, m;
    int ijps, klps;
    double zeta, Kab, xiza, zeta2;
    double  eta, Kcd, xizc,   eta2;
    double Kabcd, sqrho, dxi, rho, T, cssss, rz, re, ze2, ze22, ze23, tmp1;
    double PA[3], QC[3], QPi, WP[3], WQ[3];
    double ssss[6+1], psss[(5+1)*3], dsss[(4+1)*6], fsss[(3+1)*10];
    double ssps[3], psps[2*3], dsps[(2+1)*9], fsps[(2+1)*12];
    double psds[9], dsds[(1+1)*6], fsds[(1+1)*18];
    double DSDS[6], FSDS[18], DSFS[18], FSFS[10];
    double maxval;

    double *psss_m, *psss_m1, *dsss_m, *dsss_m1, *fsss_m, *fsss_m1;
    double *psps_m, *psps_m1, *dsps_m, *dsps_m1, *fsps_m, *fsps_m1;
    double *dsds_m, *dsds_m1, *fsds_m, *fsds_m1;
    double DPDS[18], DPFS[18], DPDP[18];
    /*
    // debug
    double c_dsps, c_psds, c_psps;
    c_dsps = c_psds = c_psps = ZERO;
    */

    for (i=0; i<6; i++)  DSDS[i] = ZERO;
    for (i=0; i<18; i++) FSDS[i] = ZERO;
    for (i=0; i<18; i++) DSFS[i] = ZERO;
    for (i=0; i<10; i++) FSFS[i] = ZERO;

    for ( ijps=0; ijps<npps; ijps++ ) {
	zeta  = vzeta[ijps];
	Kab   = vdkps[ijps];
	xiza  = vxiza[ijps];
	zeta2 = HALF*zeta;
	for ( i=0; i<3; i++ ) PA[i] = xiza*BA[i];

	for ( klps=0; klps<npps; klps++ ) {
	    eta  = vzeta[klps];
	    Kcd  = vdkps[klps];
	    xizc = vxiza[klps];
	    eta2 = HALF*eta;

	    Kabcd = Kab * Kcd;
	    sqrho = 1.e0 / sqrt( zeta + eta );
	    dxi   = xizc - xiza;
	    rho   = sqrho * sqrho;
	    T     = rho * dxi * dxi * AB2;
	    cssss = sqrho * Kabcd;
	    rz    = zeta * rho;
	    re    = eta  * rho;
	    ze2   = eta2 * rz;
	    ze22  = eta  * rz;
	    ze23  = 3.e0 * ze2;
	    for ( i=0; i<3; i++ ) {
		QC[i] = xizc*BA[i];
		QPi   = QC[i] - PA[i];
		WP[i] = rz * QPi;
		WQ[i] = rz * QPi - QPi;
	    }
	    fmt(ssss, 6, T, cssss);
	    // (ps,ss) m=0-5
	    for ( m=0, psss_m=psss; m<=5; m++, psss_m+=3 ) {
		psss_m[0]=PA[0]*ssss[m]+WP[0]*ssss[m+1];
		psss_m[1]=PA[1]*ssss[m]+WP[1]*ssss[m+1];
		psss_m[2]=PA[2]*ssss[m]+WP[2]*ssss[m+1];
	    }
	    // (ds,ss) m=0-4
	    for ( m=0, dsss_m=dsss, psss_m=psss, psss_m1=psss+3;
		    m<=4;
		    m++, dsss_m+=6, psss_m+=3, psss_m1+=3) {
		tmp1 = ssss[m]-rz*ssss[m+1];
		dsss_m[0]=PA[0]*psss_m[0]+WP[0]*psss_m1[0]+zeta2*tmp1;
		dsss_m[1]=PA[1]*psss_m[1]+WP[1]*psss_m1[1]+zeta2*tmp1;
		dsss_m[2]=PA[2]*psss_m[2]+WP[2]*psss_m1[2]+zeta2*tmp1;
		dsss_m[3]=PA[0]*psss_m[1]+WP[0]*psss_m1[1];
		dsss_m[4]=PA[1]*psss_m[2]+WP[1]*psss_m1[2];
		dsss_m[5]=PA[2]*psss_m[0]+WP[2]*psss_m1[0];
	    }
	    // (fs,ss) m=0-3
	    for ( m=0, fsss_m=fsss, dsss_m=dsss, dsss_m1=dsss+6,
		    psss_m=psss, psss_m1=psss+3;
		    m<=3;
		    m++, fsss_m+=10, dsss_m+=6, dsss_m1+=6,
		    psss_m+=3, psss_m1+=3 ) {
		fsss_m[0]=PA[0]*dsss_m[0]+WP[0]*dsss_m1[0]
			 +zeta*(psss_m[0]-rz*psss_m1[0]);
		fsss_m[1]=PA[1]*dsss_m[1]+WP[1]*dsss_m1[1]
			 +zeta*(psss_m[1]-rz*psss_m1[1]);
		fsss_m[2]=PA[2]*dsss_m[2]+WP[2]*dsss_m1[2]
			 +zeta*(psss_m[2]-rz*psss_m1[2]);

		fsss_m[3]=PA[1]*dsss_m[0]+WP[1]*dsss_m1[0];
		fsss_m[4]=PA[2]*dsss_m[1]+WP[2]*dsss_m1[1];
		fsss_m[5]=PA[0]*dsss_m[2]+WP[0]*dsss_m1[2];

		fsss_m[6]=PA[0]*dsss_m[1]+WP[0]*dsss_m1[1];
		fsss_m[7]=PA[1]*dsss_m[2]+WP[1]*dsss_m1[2];
		fsss_m[8]=PA[2]*dsss_m[0]+WP[2]*dsss_m1[0];

		fsss_m[9]=PA[2]*dsss_m[3]+WP[2]*dsss_m1[3];
	    }
	    // (ss,ps) m=2
	    ssps[0]=QC[0]*ssss[2]+WQ[0]*ssss[3];
	    ssps[1]=QC[1]*ssss[2]+WQ[1]*ssss[3];
	    ssps[2]=QC[2]*ssss[2]+WQ[2]*ssss[3];

	    // (ps,ps) m=1,2 ( only for diagonal term )
	    for ( m=1, psps_m=psps, psss_m=psss+3, psss_m1=psss+6;
		    m<=2;
		    m++, psps_m+=3, psss_m+=3, psss_m1+=3 ) {
		psps_m[0]=QC[0]*psss_m[0]+WQ[0]*psss_m1[0]+ze2*ssss[m+1];
		psps_m[1]=QC[1]*psss_m[1]+WQ[1]*psss_m1[1]+ze2*ssss[m+1];
		psps_m[2]=QC[2]*psss_m[2]+WQ[2]*psss_m1[2]+ze2*ssss[m+1];
	    }
	    /*
	    // debug
	    for ( i=0; i<2*3; i++ ) c_psps += psps[i];
	    */
	    // (ds,ps) m=0,2
	    for ( m=0, dsps_m=dsps, dsss_m=dsss, dsss_m1=dsss+6,
		    psss_m1=psss+3;
		    m<=2;
		    m++, dsps_m+=9, dsss_m+=6, dsss_m1+=6, psss_m1+=3 ) {
		// (xx,x), (yy,y), (zz,z)
		dsps_m[0]=QC[0]*dsss_m[0]+WQ[0]*dsss_m1[0]+ze22*psss_m1[0];
		dsps_m[1]=QC[1]*dsss_m[1]+WQ[1]*dsss_m1[1]+ze22*psss_m1[1];
		dsps_m[2]=QC[2]*dsss_m[2]+WQ[2]*dsss_m1[2]+ze22*psss_m1[2];
		// (xy,x), (xy,y)
		dsps_m[3]=QC[0]*dsss_m[3]+WQ[0]*dsss_m1[3]+ze2*psss_m1[1];
		dsps_m[4]=QC[1]*dsss_m[3]+WQ[1]*dsss_m1[3]+ze2*psss_m1[0];
		// (yz,y), (yz,z)
		dsps_m[5]=QC[1]*dsss_m[4]+WQ[1]*dsss_m1[4]+ze2*psss_m1[2];
		dsps_m[6]=QC[2]*dsss_m[4]+WQ[2]*dsss_m1[4]+ze2*psss_m1[1];
		// (zx,z), (zx,x)
		dsps_m[7]=QC[2]*dsss_m[5]+WQ[2]*dsss_m1[5]+ze2*psss_m1[0];
		dsps_m[8]=QC[0]*dsss_m[5]+WQ[0]*dsss_m1[5]+ze2*psss_m1[2];
	    }
	    /*
	    // debug
	    for ( i=0; i<9; i++ ) c_dsps += dsps[i];
	    */

	    // (fs,ps) m=0,2
	    for ( m=0, fsps_m=fsps, fsss_m=fsss, fsss_m1=fsss+10,
		    dsss_m1=dsss+6;
		    m<=2;
		    fsps_m+=12, fsss_m+=10, fsss_m1+=10, dsss_m1+=6, m++ ) {
		// (xxx,x), (yyy,y), (zzz,z)
		fsps_m[ 0]=QC[0]*fsss_m[0]+WQ[0]*fsss_m1[0]+ze23*dsss_m1[0];
		fsps_m[ 1]=QC[1]*fsss_m[1]+WQ[1]*fsss_m1[1]+ze23*dsss_m1[1];
		fsps_m[ 2]=QC[2]*fsss_m[2]+WQ[2]*fsss_m1[2]+ze23*dsss_m1[2];
		// (xxy,x), (yyz,y), (zzx,z)
		fsps_m[ 3]=QC[0]*fsss_m[3]+WQ[0]*fsss_m1[3]+ze22*dsss_m1[3];
		fsps_m[ 4]=QC[1]*fsss_m[4]+WQ[1]*fsss_m1[4]+ze22*dsss_m1[4];
		fsps_m[ 5]=QC[2]*fsss_m[5]+WQ[2]*fsss_m1[5]+ze22*dsss_m1[5];
		// (xyy,y), (yzz,z), (zxx,x)
		fsps_m[ 6]=QC[1]*fsss_m[6]+WQ[1]*fsss_m1[6]+ze22*dsss_m1[3];
		fsps_m[ 7]=QC[2]*fsss_m[7]+WQ[2]*fsss_m1[7]+ze22*dsss_m1[4];
		fsps_m[ 8]=QC[0]*fsss_m[8]+WQ[0]*fsss_m1[8]+ze22*dsss_m1[5];
		// (xyz,x), (xyz,y), (xyz,z)
		fsps_m[ 9]=QC[0]*fsss_m[9]+WQ[0]*fsss_m1[9]+ze2*dsss_m1[4];
		fsps_m[10]=QC[1]*fsss_m[9]+WQ[1]*fsss_m1[9]+ze2*dsss_m1[5];
		fsps_m[11]=QC[2]*fsss_m[9]+WQ[2]*fsss_m1[9]+ze2*dsss_m1[3];
	    }
	    // (ps,ds) m=1
	    m=1;
	    psps_m  = psps + 0;	// psps is 1-offset
	    psps_m1 = psps + 3;
	    psss_m  = psss + 3;
	    psss_m1 = psss + 6;
	    // (x,xx), (y,yy), (z,zz)
	    psds[0]=QC[0]*psps_m[0]+WQ[0]*psps_m1[0]
		+eta2*(psss_m[0]-re*psss_m1[0]) + ze2*ssps[0];
	    psds[1]=QC[1]*psps_m[1]+WQ[1]*psps_m1[1]
		+eta2*(psss_m[1]-re*psss_m1[1]) + ze2*ssps[1];
	    psds[2]=QC[2]*psps_m[2]+WQ[2]*psps_m1[2]
		+eta2*(psss_m[2]-re*psss_m1[2]) + ze2*ssps[2];
	    // (y,xy), (x,xy)
	    psds[3]=QC[0]*psps_m[1]+WQ[0]*psps_m1[1];
	    psds[4]=QC[1]*psps_m[0]+WQ[1]*psps_m1[0];
	    // (z,yz), (y,yz)
	    psds[5]=QC[1]*psps_m[2]+WQ[1]*psps_m1[2];
	    psds[6]=QC[2]*psps_m[1]+WQ[2]*psps_m1[1];
	    // (x,zx), (z,zx)
	    psds[7]=QC[2]*psps_m[0]+WQ[2]*psps_m1[0];
	    psds[8]=QC[0]*psps_m[2]+WQ[0]*psps_m1[2];
	    /*
	    // debug
	    for ( i=0; i<9; i++ ) c_psds += psds[i];
	    */
	    // (ds,ds) m=0,1 (diagonal part)
	    for (m=0, dsds_m=dsds, dsps_m=dsps, dsps_m1=dsps+9,
		    dsss_m=dsss, dsss_m1=dsss+6, psps_m1=psps;
		    m<=1;
		    dsds_m+=6, dsps_m+=9, dsps_m1+=9, dsss_m+=6, dsss_m1+=6,
		    psps_m1+=3, m++ ) {
		dsds_m[0]=QC[0]*dsps_m[0]+WQ[0]*dsps_m1[0]
		    +eta2*(dsss_m[0]-re*dsss_m1[0])+ze22*psps_m1[0];
		dsds_m[1]=QC[1]*dsps_m[1]+WQ[1]*dsps_m1[1]
		    +eta2*(dsss_m[1]-re*dsss_m1[1])+ze22*psps_m1[1];
		dsds_m[2]=QC[2]*dsps_m[2]+WQ[2]*dsps_m1[2]
		    +eta2*(dsss_m[2]-re*dsss_m1[2])+ze22*psps_m1[2];
		dsds_m[3]=QC[0]*dsps_m[4]+WQ[0]*dsps_m1[4]+ze2*psps_m1[1];
		dsds_m[4]=QC[1]*dsps_m[6]+WQ[1]*dsps_m1[6]+ze2*psps_m1[2];
		dsds_m[5]=QC[2]*dsps_m[8]+WQ[2]*dsps_m1[8]+ze2*psps_m1[0];
	    }
	    for (i=0; i<6; i++) DSDS[i]+=dsds[i];
	    // (fs,ds) m=0,1
	    for ( m=0, fsds_m=fsds, fsps_m=fsps, fsps_m1=fsps+12,
		    fsss_m=fsss, fsss_m1=fsss+10, dsps_m1=dsps+9;
		    m<=1;
		    fsds_m+=18, fsps_m+=12, fsps_m1+=12,
		    fsss_m+=10, fsss_m1+=10, dsps_m1+=9, m++ ) {
		// (xxx,xx), (yyy,yy), (zzz,zz)
		fsds_m[ 0]=QC[0]*fsps_m[0]+WQ[0]*fsps_m1[0]+
		    eta2*(fsss_m[0]-re*fsss_m1[0])+ze23*dsps_m1[0];
		fsds_m[ 1]=QC[1]*fsps_m[1]+WQ[1]*fsps_m1[1]+
		    eta2*(fsss_m[1]-re*fsss_m1[1])+ze23*dsps_m1[1];
		fsds_m[ 2]=QC[2]*fsps_m[2]+WQ[2]*fsps_m1[2]+
		    eta2*(fsss_m[2]-re*fsss_m1[2])+ze23*dsps_m1[2];
		// (xxy,xx), (xxy,xy)
		fsds_m[ 3]=QC[0]*fsps_m[3]+WQ[0]*fsps_m1[3]+
		    eta2*(fsss_m[3]-re*fsss_m1[3])+ze22*dsps_m1[3];
		fsds_m[ 4]=QC[1]*fsps_m[3]+WQ[1]*fsps_m1[3]+ze2*dsps_m1[8];
		// (yyz,yy), (yyz,yz)
		fsds_m[ 5]=QC[1]*fsps_m[4]+WQ[1]*fsps_m1[4]+
		    eta2*(fsss_m[4]-re*fsss_m1[4])+ze22*dsps_m1[5];
		fsds_m[ 6]=QC[2]*fsps_m[4]+WQ[2]*fsps_m1[4]+ze2*dsps_m1[1];
		// (zzx,zz), (zzx,zx)
		fsds_m[ 7]=QC[2]*fsps_m[5]+WQ[2]*fsps_m1[5]+
		    eta2*(fsss_m[5]-re*fsss_m1[5])+ze22*dsps_m1[7];
		fsds_m[ 8]=QC[0]*fsps_m[5]+WQ[0]*fsps_m1[5]+ze2*dsps_m1[2];
		// (xyy,yy), (xyy,xy)
		fsds_m[ 9]=QC[1]*fsps_m[6]+WQ[1]*fsps_m1[6]+
		    eta2*(fsss_m[6]-re*fsss_m1[6])+ze22*dsps_m1[4];
		fsds_m[10]=QC[0]*fsps_m[6]+WQ[0]*fsps_m1[6]+ze2*dsps_m1[1];
		// (yzz,zz), (yzz,yz)
		fsds_m[11]=QC[2]*fsps_m[7]+WQ[2]*fsps_m1[7]+
		    eta2*(fsss_m[7]-re*fsss_m1[7])+ze22*dsps_m1[6];
		fsds_m[12]=QC[1]*fsps_m[7]+WQ[1]*fsps_m1[7]+ze2*dsps_m1[2];
		// (zxx,xx), (zxx,zx)
		fsds_m[13]=QC[0]*fsps_m[8]+WQ[0]*fsps_m1[8]+
		    eta2*(fsss_m[8]-re*fsss_m1[8])+ze22*dsps_m1[8];
		fsds_m[14]=QC[2]*fsps_m[8]+WQ[2]*fsps_m1[8]+ze2*dsps_m1[0];
		// (xyz,xy), (xyz,yz), (xyz,zx)
		fsds_m[15]=QC[1]*fsps_m[ 9]+WQ[1]*fsps_m1[ 9]+ze2*dsps_m1[8];
		fsds_m[16]=QC[2]*fsps_m[10]+WQ[2]*fsps_m1[10]+ze2*dsps_m1[4];
		fsds_m[17]=QC[0]*fsps_m[11]+WQ[0]*fsps_m1[11]+ze2*dsps_m1[6];
	    }
	    for (i=0; i<18; i++) FSDS[i]+=fsds[i];
	    // (ds,fs) m=0
	    dsds_m  = dsds;
	    dsds_m1 = dsds + 6;
	    dsps_m  = dsps;
	    dsps_m1 = dsps + 9;
	    // (xx,xxx), (xx,xxy), (xx,zxx)
	    DSFS[0]+=QC[0]*dsds_m[0]+WQ[0]*dsds_m1[0]
		+eta*(dsps_m[0]-re*dsps_m1[0])+ze22*psds[0];
	    DSFS[1]+=QC[1]*dsds_m[0]+WQ[1]*dsds_m1[0];
	    DSFS[2]+=QC[2]*dsds_m[0]+WQ[2]*dsds_m1[0];
	    // (yy,xyy), (yy,yyy), (yy,yyz)
	    DSFS[3]+=QC[0]*dsds_m[1]+WQ[0]*dsds_m1[1];
	    DSFS[4]+=QC[1]*dsds_m[1]+WQ[1]*dsds_m1[1]
		+eta*(dsps_m[1]-re*dsps_m1[1])+ze22*psds[1];
	    DSFS[5]+=QC[2]*dsds_m[1]+WQ[2]*dsds_m1[1];
	    // (zz,zzx), (zz,yzz), (zz,zzz)
	    DSFS[6]+=QC[0]*dsds_m[2]+WQ[0]*dsds_m1[2];
	    DSFS[7]+=QC[1]*dsds_m[2]+WQ[1]*dsds_m1[2];
	    DSFS[8]+=QC[2]*dsds_m[2]+WQ[2]*dsds_m1[2]
		+eta*(dsps_m[2]-re*dsps_m1[2])+ze22*psds[2];
	    // (xy,xxy), (xy,xyy), (xy,xyz)
	    DSFS[ 9]+=QC[0]*dsds_m[3]+WQ[0]*dsds_m1[3]
		+eta2*(dsps_m[4]-re*dsps_m1[4])+ze2*psds[3];
	    DSFS[10]+=QC[1]*dsds_m[3]+WQ[1]*dsds_m1[3]
		+eta2*(dsps_m[3]-re*dsps_m1[3])+ze2*psds[4];
	    DSFS[11]+=QC[2]*dsds_m[3]+WQ[2]*dsds_m1[3];
	    // (yz,zyz), (yz,yyz), (yz,yzz)
	    DSFS[12]+=QC[0]*dsds_m[4]+WQ[0]*dsds_m1[4];
	    DSFS[13]+=QC[1]*dsds_m[4]+WQ[1]*dsds_m1[4]
		+eta2*(dsps_m[6]-re*dsps_m1[6])+ze2*psds[5];
	    DSFS[14]+=QC[2]*dsds_m[4]+WQ[2]*dsds_m1[4]
		+eta2*(dsps_m[5]-re*dsps_m1[5])+ze2*psds[6];
	    // (zx,zxx), (zx,xyz), (zx,zzx)
	    DSFS[15]+=QC[0]*dsds_m[5]+WQ[0]*dsds_m1[5]
		+eta2*(dsps_m[7]-re*dsps_m1[7])+ze2*psds[8];
	    DSFS[16]+=QC[1]*dsds_m[5]+WQ[1]*dsds_m1[5];
	    DSFS[17]+=QC[2]*dsds_m[5]+WQ[2]*dsds_m1[5]
		+eta2*(dsps_m[8]-re*dsps_m1[8])+ze2*psds[7];
	    // (fs,fs) m=0
	    fsds_m  = fsds;
	    fsds_m1 = fsds+18;
	    fsps_m  = fsps;
	    fsps_m1 = fsps+12;
	    dsds_m1 = dsds+6;
	    // (xxx,xxx), (yyy,yyy), (zzz,zzz)
	    FSFS[0]+=QC[0]*fsds_m[0]+WQ[0]*fsds_m1[0]
		+eta*(fsps_m[0]-re*fsps_m1[0])+ze23*dsds_m1[0];
	    FSFS[1]+=QC[1]*fsds_m[1]+WQ[1]*fsds_m1[1]
		+eta*(fsps_m[1]-re*fsps_m1[1])+ze23*dsds_m1[1];
	    FSFS[2]+=QC[2]*fsds_m[2]+WQ[2]*fsds_m1[2]
		+eta*(fsps_m[2]-re*fsps_m1[2])+ze23*dsds_m1[2];
	    // (xxy,xxy), (yyz,yyz), (zzx,zzx)
	    FSFS[3]+=QC[1]*fsds_m[3]+WQ[1]*fsds_m1[3]+ze2*dsds_m1[0];
	    FSFS[4]+=QC[2]*fsds_m[5]+WQ[2]*fsds_m1[5]+ze2*dsds_m1[1];
	    FSFS[5]+=QC[0]*fsds_m[7]+WQ[0]*fsds_m1[7]+ze2*dsds_m1[2];
	    // (xyy,xyy), (yzz,yzz), (zxx,zxx)
	    FSFS[6]+=QC[0]*fsds_m[ 9]+WQ[0]*fsds_m1[ 9]+ze2*dsds_m1[1];
	    FSFS[7]+=QC[1]*fsds_m[11]+WQ[1]*fsds_m1[11]+ze2*dsds_m1[2];
	    FSFS[8]+=QC[2]*fsds_m[13]+WQ[2]*fsds_m1[13]+ze2*dsds_m1[0];
	    // (xyz,xyz)
	    FSFS[9]+=QC[2]*fsds_m[15]+WQ[2]*fsds_m1[15]+ze2*dsds_m1[3];
	}	// klps
    }		// ijps
    // HRR
    // (dp,ds) = (fs,ds)-BA(ds,ds)
    DPDS[ 0]=FSDS[ 0]-BA[0]*DSDS[0];
    DPDS[ 1]=FSDS[ 3]-BA[1]*DSDS[0];
    DPDS[ 2]=FSDS[13]-BA[2]*DSDS[0];

    DPDS[ 3]=FSDS[ 9]-BA[0]*DSDS[1];
    DPDS[ 4]=FSDS[ 1]-BA[1]*DSDS[1];
    DPDS[ 5]=FSDS[ 5]-BA[2]*DSDS[1];

    DPDS[ 6]=FSDS[ 7]-BA[0]*DSDS[2];
    DPDS[ 7]=FSDS[11]-BA[1]*DSDS[2];
    DPDS[ 8]=FSDS[ 2]-BA[2]*DSDS[2];

    DPDS[ 9]=FSDS[ 4]-BA[0]*DSDS[3];
    DPDS[10]=FSDS[10]-BA[1]*DSDS[3];
    DPDS[11]=FSDS[15]-BA[2]*DSDS[3];

    DPDS[12]=FSDS[16]-BA[0]*DSDS[4];
    DPDS[13]=FSDS[ 6]-BA[1]*DSDS[4];
    DPDS[14]=FSDS[12]-BA[2]*DSDS[4];

    DPDS[15]=FSDS[14]-BA[0]*DSDS[5];
    DPDS[16]=FSDS[17]-BA[1]*DSDS[5];
    DPDS[17]=FSDS[ 8]-BA[2]*DSDS[5];
    // ( dp,fs) = (fs,fs)-BA(ds,fs)
    DPFS[ 0]=FSFS[0]-BA[0]*DSFS[ 0];
    DPFS[ 1]=FSFS[3]-BA[1]*DSFS[ 1];
    DPFS[ 2]=FSFS[8]-BA[2]*DSFS[ 2];

    DPFS[ 3]=FSFS[6]-BA[0]*DSFS[ 3];
    DPFS[ 4]=FSFS[1]-BA[1]*DSFS[ 4];
    DPFS[ 5]=FSFS[4]-BA[2]*DSFS[ 5];

    DPFS[ 6]=FSFS[5]-BA[0]*DSFS[ 6];
    DPFS[ 7]=FSFS[7]-BA[1]*DSFS[ 7];
    DPFS[ 8]=FSFS[2]-BA[2]*DSFS[ 8];

    DPFS[ 9]=FSFS[3]-BA[0]*DSFS[ 9];
    DPFS[10]=FSFS[6]-BA[1]*DSFS[10];
    DPFS[11]=FSFS[9]-BA[2]*DSFS[11];

    DPFS[12]=FSFS[9]-BA[0]*DSFS[12];
    DPFS[13]=FSFS[4]-BA[1]*DSFS[13];
    DPFS[14]=FSFS[7]-BA[2]*DSFS[14];

    DPFS[15]=FSFS[8]-BA[0]*DSFS[15];
    DPFS[16]=FSFS[9]-BA[1]*DSFS[16];
    DPFS[17]=FSFS[5]-BA[2]*DSFS[17];
    // (dp,dp) = (dp,fs)-DC(dp,ds)
    for ( i=0; i<18; i+=3 ) {
	DPDP[i+0]=DPFS[i+0]-BA[0]*DPDS[i+0];
	DPDP[i+1]=DPFS[i+1]-BA[1]*DPDS[i+1];
	DPDP[i+2]=DPFS[i+2]-BA[2]*DPDS[i+2];
    }
    for ( i=9; i<18; i++ ) DPDP[i]*=3.e0;

    for ( i=0; i<18; i++ ) DPDP[i]=fabs(DPDP[i]);
    maxval = DPDP[0];
    for ( i=1; i<18; i++ ) if ( DPDP[i]>maxval ) maxval = DPDP[i];

    return maxval;
}
