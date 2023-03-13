/**
 * @file ofmo-cutoff-core-dddd.c
 * １つのSchwarz積分計算関数群（２）
 *
 * １つのSchwarz積分を計算する関数群のうち、(dd,dd)タイプの
 * 計算のための関数が定義してあるファイル。
 * 計算された複数のSchwarz積分の最大値を返す。
 * */
#include <stdio.h>
#include <math.h>

#define ZERO 0.e0
#define HALF .5e0

#ifndef false
#define false 0
#endif
#ifndef true
#define true 1
#endif

#define EPS_PS_PAIR 1.e-32

extern void fmt( double*, int, double, double );

/** (dd,dd)タイプのSchwarz積分計算関数
 * @ingroup core-cutoff
 * */
double schwarz_core_dddd_(
	const int npps, const double vzeta[], const double vdkps[],
	const double vxiza[], const double BA[3], const double AB2 ) {
    int i;
    int ijps, klps;
    double zeta, Kab, xiza, zeta2;
    double  eta, Kcd, xizc,   eta2;
    double Kabcd, sqrho, dxi, rho, T, cssss, rz, re, ze2, rze2[4+1];
    double PA[3], QC[3], QPi, WP[3], WQ[3];
    double maxval;
    // ========== intermediate integ. in VRR ==========
    double gsgs[15*15], ssss[8+1];
    double psss[7+1][3], dsss[6+1][6], fsss[5+1][10], gsss[4+1][15];
    double ssps[1+1][3], psps[2+1][3*3], dsps[3+1][6*3];
    double fsps[3+1][10*3], gsps[3+1][15*3], ssds[6];
    double psds[1+1][3*6], dsds[2+1][6*6], fsds[2+1][10*6];
    double gsds[2+1][15*6], psfs[3*10], dsfs[1+1][6*10];
    double fsfs[1+1][10*10], gsfs[1+1][15*10], dsgs[6*15];
    double fsgs[10*15];
    // ==========  work area used in VRR ==========
    union _temp_ {
        double entity[6];
        double ssss;
        double psss[3];
        double dsss[6];
        double fsss[10];
        double gsss[15];
        double psps[3*3];
        double dsps[6*3];
        double fsps[10*3];
        double gsps[15*3];
        double dsds[6*6];
        double fsds[10*6];
        double gsds[15*6];
    } tmp;
    // ========== contracted integ, calculated by VRR ==========
    double DSDS[6*6], FSDS[10*6], GSDS[15*6];
    double DSFS[6*10], FSFS[10*10], GSFS[15*10];
    double DSGS[6*15], FSGS[10*15], GSGS[15*15];
    // ========== intermediate integ. in HRR ==========
    double DPDS[6*3*6], DPFS[6*3*10], DPGS[6*3*15];
    double FPDS[10*3*6], FPFS[10*3*10], FPGS[10*3*15];
    double DDDS[6*6*6], DDFS[6*6*10], DDGS[6*6*15];
    double DDDP[6*6*6*3], DDFP[6*6*10*3], DDDD[6*6*6*6];
    // added in 2010.07.09
    double DDDD_diag[6*6];

    // ---- zero clear for contracted integrals ----
    for (int i=0; i<6*6; i++)   DSDS[i] = ZERO;
    for (int i=0; i<6*10; i++)  DSFS[i] = ZERO;
    for (int i=0; i<6*15; i++)  DSGS[i] = ZERO;
    for (int i=0; i<10*6; i++)  FSDS[i] = ZERO;
    for (int i=0; i<10*10; i++) FSFS[i] = ZERO;
    for (int i=0; i<10*15; i++) FSGS[i] = ZERO;
    for (int i=0; i<15*6; i++)  GSDS[i] = ZERO;
    for (int i=0; i<15*10; i++) GSFS[i] = ZERO;
    for (int i=0; i<15*15; i++) GSGS[i] = ZERO;

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
	    rze2[1] = ze2 = eta2 * rz;
	    rze2[2] = eta * rz;
	    rze2[3] = ze2*3.e0;
	    rze2[4] = ze2*4.e0;
	    for ( i=0; i<3; i++ ) {
		QC[i] = xizc*BA[i];
		QPi   = QC[i] - PA[i];
		WP[i] = rz * QPi;
		WQ[i] = rz * QPi - QPi;
	    }
	    fmt(ssss, 6, T, cssss);
            
	    // ========== VRR calculation (2-center) ==========
	    // (P,S|S,S)
	    for (int m=0; m<=7; m++) {
		int m1 = m+1;
		psss[m][0] = PA[0]*ssss[m] + WP[0]*ssss[m1];
		psss[m][1] = PA[1]*ssss[m] + WP[1]*ssss[m1];
		psss[m][2] = PA[2]*ssss[m] + WP[2]*ssss[m1];
	    }	// end m loop
	    // (D,S|S,S)
	    for (int m=0; m<=6; m++) {
		int m1 = m+1;
		tmp.ssss = zeta2*( ssss[m] - rz*ssss[m1] );

		dsss[m][0] = PA[0]*psss[m][0] + WP[0]*psss[m1][0]
			   +   tmp.ssss;
		dsss[m][1] = PA[1]*psss[m][1] + WP[1]*psss[m1][1]
			   +   tmp.ssss;
		dsss[m][2] = PA[2]*psss[m][2] + WP[2]*psss[m1][2]
			   +   tmp.ssss;
		dsss[m][3] = PA[0]*psss[m][1] + WP[0]*psss[m1][1];
		dsss[m][4] = PA[1]*psss[m][2] + WP[1]*psss[m1][2];
		dsss[m][5] = PA[2]*psss[m][0] + WP[2]*psss[m1][0];
	    }	// end m loop
	    // (F,S|S,S)
	    for (int m=0; m<=5; m++) {
		int m1 = m+1;
		for (int i=0; i<3; i++)
		    tmp.psss[i] = zeta2 * ( psss[m][i] - rz*psss[m1][i] );

		fsss[m][0] = PA[0]*dsss[m][0] + WP[0]*dsss[m1][0]
			   + 2*tmp.psss[0];
		fsss[m][1] = PA[1]*dsss[m][1] + WP[1]*dsss[m1][1]
			   + 2*tmp.psss[1];
		fsss[m][2] = PA[2]*dsss[m][2] + WP[2]*dsss[m1][2]
			   + 2*tmp.psss[2];
		fsss[m][3] = PA[0]*dsss[m][3] + WP[0]*dsss[m1][3]
			   +   tmp.psss[1];
		fsss[m][4] = PA[1]*dsss[m][4] + WP[1]*dsss[m1][4]
			   +   tmp.psss[2];
		fsss[m][5] = PA[2]*dsss[m][5] + WP[2]*dsss[m1][5]
			   +   tmp.psss[0];
		fsss[m][6] = PA[0]*dsss[m][1] + WP[0]*dsss[m1][1];
		fsss[m][7] = PA[1]*dsss[m][2] + WP[1]*dsss[m1][2];
		fsss[m][8] = PA[2]*dsss[m][0] + WP[2]*dsss[m1][0];
		fsss[m][9] = PA[0]*dsss[m][4] + WP[0]*dsss[m1][4];
	    }	// end m loop
	    // (G,S|S,S)
	    for (int m=0; m<=4; m++) {
		int m1 = m+1;
		for (int i=0; i<6; i++)
		    tmp.dsss[i] = zeta2 * ( dsss[m][i] - rz*dsss[m1][i] );

		gsss[m][ 0] = PA[0]*fsss[m][0] + WP[0]*fsss[m1][0]
			    + 3*tmp.dsss[0];
		gsss[m][ 1] = PA[1]*fsss[m][1] + WP[1]*fsss[m1][1]
			    + 3*tmp.dsss[1];
		gsss[m][ 2] = PA[2]*fsss[m][2] + WP[2]*fsss[m1][2]
			    + 3*tmp.dsss[2];
		gsss[m][ 3] = PA[0]*fsss[m][3] + WP[0]*fsss[m1][3]
			    + 2*tmp.dsss[3];
		gsss[m][ 4] = PA[1]*fsss[m][4] + WP[1]*fsss[m1][4]
			    + 2*tmp.dsss[4];
		gsss[m][ 5] = PA[2]*fsss[m][5] + WP[2]*fsss[m1][5]
			    + 2*tmp.dsss[5];
		gsss[m][ 6] = PA[0]*fsss[m][6] + WP[0]*fsss[m1][6]
			    +   tmp.dsss[1];
		gsss[m][ 7] = PA[1]*fsss[m][7] + WP[1]*fsss[m1][7]
			    +   tmp.dsss[2];
		gsss[m][ 8] = PA[2]*fsss[m][8] + WP[2]*fsss[m1][8]
			    +   tmp.dsss[0];
		gsss[m][ 9] = PA[0]*fsss[m][1] + WP[0]*fsss[m1][1];
		gsss[m][10] = PA[1]*fsss[m][2] + WP[1]*fsss[m1][2];
		gsss[m][11] = PA[2]*fsss[m][0] + WP[2]*fsss[m1][0];
		gsss[m][12] = PA[0]*fsss[m][9] + WP[0]*fsss[m1][9]
			    +   tmp.dsss[4];
		gsss[m][13] = PA[1]*fsss[m][9] + WP[1]*fsss[m1][9]
			    +   tmp.dsss[5];
		gsss[m][14] = PA[2]*fsss[m][9] + WP[2]*fsss[m1][9]
			    +   tmp.dsss[3];
	    }	// end m loop
	    // (S,S|P,S)
	    for (int m=2; m<=3; m++) {
		int m1 = m+1;
		ssps[m-2][0] = QC[0]*ssss[m] + WQ[0]*ssss[m1];
		ssps[m-2][1] = QC[1]*ssss[m] + WQ[1]*ssss[m1];
		ssps[m-2][2] = QC[2]*ssss[m] + WQ[2]*ssss[m1];
	    }	// end m loop
	    // (P,S|P,S)
	    for (int m=1; m<=3; m++) {
		int m1 = m+1;
		psps[m-1][0*3+0] = QC[0]*psss[m][0] + WQ[0]*psss[m1][0]
				 + rze2[1]*ssss[m1];
		psps[m-1][0*3+1] = QC[1]*psss[m][0] + WQ[1]*psss[m1][0];
		psps[m-1][0*3+2] = QC[2]*psss[m][0] + WQ[2]*psss[m1][0];
		psps[m-1][1*3+0] = QC[0]*psss[m][1] + WQ[0]*psss[m1][1];
		psps[m-1][1*3+1] = QC[1]*psss[m][1] + WQ[1]*psss[m1][1]
				 + rze2[1]*ssss[m1];
		psps[m-1][1*3+2] = QC[2]*psss[m][1] + WQ[2]*psss[m1][1];
		psps[m-1][2*3+0] = QC[0]*psss[m][2] + WQ[0]*psss[m1][2];
		psps[m-1][2*3+1] = QC[1]*psss[m][2] + WQ[1]*psss[m1][2];
		psps[m-1][2*3+2] = QC[2]*psss[m][2] + WQ[2]*psss[m1][2]
				 + rze2[1]*ssss[m1];
	    }	// end m loop
	    // (D,S|P,S)
	    for (int m=0; m<=3; m++) {
		int m1 = m+1;
		dsps[m][0*3+0] = QC[0]*dsss[m][0] + WQ[0]*dsss[m1][0]
			       + rze2[2]*psss[m1][0];
		dsps[m][0*3+1] = QC[1]*dsss[m][0] + WQ[1]*dsss[m1][0];
		dsps[m][0*3+2] = QC[2]*dsss[m][0] + WQ[2]*dsss[m1][0];
		dsps[m][1*3+0] = QC[0]*dsss[m][1] + WQ[0]*dsss[m1][1];
		dsps[m][1*3+1] = QC[1]*dsss[m][1] + WQ[1]*dsss[m1][1]
			       + rze2[2]*psss[m1][1];
		dsps[m][1*3+2] = QC[2]*dsss[m][1] + WQ[2]*dsss[m1][1];
		dsps[m][2*3+0] = QC[0]*dsss[m][2] + WQ[0]*dsss[m1][2];
		dsps[m][2*3+1] = QC[1]*dsss[m][2] + WQ[1]*dsss[m1][2];
		dsps[m][2*3+2] = QC[2]*dsss[m][2] + WQ[2]*dsss[m1][2]
			       + rze2[2]*psss[m1][2];
		dsps[m][3*3+0] = QC[0]*dsss[m][3] + WQ[0]*dsss[m1][3]
			       + rze2[1]*psss[m1][1];
		dsps[m][3*3+1] = QC[1]*dsss[m][3] + WQ[1]*dsss[m1][3]
			       + rze2[1]*psss[m1][0];
		dsps[m][3*3+2] = QC[2]*dsss[m][3] + WQ[2]*dsss[m1][3];
		dsps[m][4*3+0] = QC[0]*dsss[m][4] + WQ[0]*dsss[m1][4];
		dsps[m][4*3+1] = QC[1]*dsss[m][4] + WQ[1]*dsss[m1][4]
			       + rze2[1]*psss[m1][2];
		dsps[m][4*3+2] = QC[2]*dsss[m][4] + WQ[2]*dsss[m1][4]
			       + rze2[1]*psss[m1][1];
		dsps[m][5*3+0] = QC[0]*dsss[m][5] + WQ[0]*dsss[m1][5]
			       + rze2[1]*psss[m1][2];
		dsps[m][5*3+1] = QC[1]*dsss[m][5] + WQ[1]*dsss[m1][5];
		dsps[m][5*3+2] = QC[2]*dsss[m][5] + WQ[2]*dsss[m1][5]
			       + rze2[1]*psss[m1][0];
	    }	// end m loop
	    // (F,S|P,S)
	    for (int m=0; m<=3; m++) {
		int m1 = m+1;
		fsps[m][0*3+0] = QC[0]*fsss[m][0] + WQ[0]*fsss[m1][0]
			       + rze2[3]*dsss[m1][0];
		fsps[m][0*3+1] = QC[1]*fsss[m][0] + WQ[1]*fsss[m1][0];
		fsps[m][0*3+2] = QC[2]*fsss[m][0] + WQ[2]*fsss[m1][0];
		fsps[m][1*3+0] = QC[0]*fsss[m][1] + WQ[0]*fsss[m1][1];
		fsps[m][1*3+1] = QC[1]*fsss[m][1] + WQ[1]*fsss[m1][1]
			       + rze2[3]*dsss[m1][1];
		fsps[m][1*3+2] = QC[2]*fsss[m][1] + WQ[2]*fsss[m1][1];
		fsps[m][2*3+0] = QC[0]*fsss[m][2] + WQ[0]*fsss[m1][2];
		fsps[m][2*3+1] = QC[1]*fsss[m][2] + WQ[1]*fsss[m1][2];
		fsps[m][2*3+2] = QC[2]*fsss[m][2] + WQ[2]*fsss[m1][2]
			       + rze2[3]*dsss[m1][2];
		fsps[m][3*3+0] = QC[0]*fsss[m][3] + WQ[0]*fsss[m1][3]
			       + rze2[2]*dsss[m1][3];
		fsps[m][3*3+1] = QC[1]*fsss[m][3] + WQ[1]*fsss[m1][3]
			       + rze2[1]*dsss[m1][0];
		fsps[m][3*3+2] = QC[2]*fsss[m][3] + WQ[2]*fsss[m1][3];
		fsps[m][4*3+0] = QC[0]*fsss[m][4] + WQ[0]*fsss[m1][4];
		fsps[m][4*3+1] = QC[1]*fsss[m][4] + WQ[1]*fsss[m1][4]
			       + rze2[2]*dsss[m1][4];
		fsps[m][4*3+2] = QC[2]*fsss[m][4] + WQ[2]*fsss[m1][4]
			       + rze2[1]*dsss[m1][1];
		fsps[m][5*3+0] = QC[0]*fsss[m][5] + WQ[0]*fsss[m1][5]
			       + rze2[1]*dsss[m1][2];
		fsps[m][5*3+1] = QC[1]*fsss[m][5] + WQ[1]*fsss[m1][5];
		fsps[m][5*3+2] = QC[2]*fsss[m][5] + WQ[2]*fsss[m1][5]
			       + rze2[2]*dsss[m1][5];
		fsps[m][6*3+0] = QC[0]*fsss[m][6] + WQ[0]*fsss[m1][6]
			       + rze2[1]*dsss[m1][1];
		fsps[m][6*3+1] = QC[1]*fsss[m][6] + WQ[1]*fsss[m1][6]
			       + rze2[2]*dsss[m1][3];
		fsps[m][6*3+2] = QC[2]*fsss[m][6] + WQ[2]*fsss[m1][6];
		fsps[m][7*3+0] = QC[0]*fsss[m][7] + WQ[0]*fsss[m1][7];
		fsps[m][7*3+1] = QC[1]*fsss[m][7] + WQ[1]*fsss[m1][7]
			       + rze2[1]*dsss[m1][2];
		fsps[m][7*3+2] = QC[2]*fsss[m][7] + WQ[2]*fsss[m1][7]
			       + rze2[2]*dsss[m1][4];
		fsps[m][8*3+0] = QC[0]*fsss[m][8] + WQ[0]*fsss[m1][8]
			       + rze2[2]*dsss[m1][5];
		fsps[m][8*3+1] = QC[1]*fsss[m][8] + WQ[1]*fsss[m1][8];
		fsps[m][8*3+2] = QC[2]*fsss[m][8] + WQ[2]*fsss[m1][8]
			       + rze2[1]*dsss[m1][0];
		fsps[m][9*3+0] = QC[0]*fsss[m][9] + WQ[0]*fsss[m1][9]
			       + rze2[1]*dsss[m1][4];
		fsps[m][9*3+1] = QC[1]*fsss[m][9] + WQ[1]*fsss[m1][9]
			       + rze2[1]*dsss[m1][5];
		fsps[m][9*3+2] = QC[2]*fsss[m][9] + WQ[2]*fsss[m1][9]
			       + rze2[1]*dsss[m1][3];
	    }	// end m loop
	    // (G,S|P,S)
	    for (int m=0; m<=3; m++) {
		int m1 = m+1;
		gsps[m][ 0*3+0] = QC[0]*gsss[m][ 0] + WQ[0]*gsss[m1][ 0]
				+ rze2[4]*fsss[m1][0];
		gsps[m][ 0*3+1] = QC[1]*gsss[m][ 0] + WQ[1]*gsss[m1][ 0];
		gsps[m][ 0*3+2] = QC[2]*gsss[m][ 0] + WQ[2]*gsss[m1][ 0];
		gsps[m][ 1*3+0] = QC[0]*gsss[m][ 1] + WQ[0]*gsss[m1][ 1];
		gsps[m][ 1*3+1] = QC[1]*gsss[m][ 1] + WQ[1]*gsss[m1][ 1]
				+ rze2[4]*fsss[m1][1];
		gsps[m][ 1*3+2] = QC[2]*gsss[m][ 1] + WQ[2]*gsss[m1][ 1];
		gsps[m][ 2*3+0] = QC[0]*gsss[m][ 2] + WQ[0]*gsss[m1][ 2];
		gsps[m][ 2*3+1] = QC[1]*gsss[m][ 2] + WQ[1]*gsss[m1][ 2];
		gsps[m][ 2*3+2] = QC[2]*gsss[m][ 2] + WQ[2]*gsss[m1][ 2]
				+ rze2[4]*fsss[m1][2];
		gsps[m][ 3*3+0] = QC[0]*gsss[m][ 3] + WQ[0]*gsss[m1][ 3]
				+ rze2[3]*fsss[m1][3];
		gsps[m][ 3*3+1] = QC[1]*gsss[m][ 3] + WQ[1]*gsss[m1][ 3]
				+ rze2[1]*fsss[m1][0];
		gsps[m][ 3*3+2] = QC[2]*gsss[m][ 3] + WQ[2]*gsss[m1][ 3];
		gsps[m][ 4*3+0] = QC[0]*gsss[m][ 4] + WQ[0]*gsss[m1][ 4];
		gsps[m][ 4*3+1] = QC[1]*gsss[m][ 4] + WQ[1]*gsss[m1][ 4]
				+ rze2[3]*fsss[m1][4];
		gsps[m][ 4*3+2] = QC[2]*gsss[m][ 4] + WQ[2]*gsss[m1][ 4]
				+ rze2[1]*fsss[m1][1];
		gsps[m][ 5*3+0] = QC[0]*gsss[m][ 5] + WQ[0]*gsss[m1][ 5]
				+ rze2[1]*fsss[m1][2];
		gsps[m][ 5*3+1] = QC[1]*gsss[m][ 5] + WQ[1]*gsss[m1][ 5];
		gsps[m][ 5*3+2] = QC[2]*gsss[m][ 5] + WQ[2]*gsss[m1][ 5]
				+ rze2[3]*fsss[m1][5];
		gsps[m][ 6*3+0] = QC[0]*gsss[m][ 6] + WQ[0]*gsss[m1][ 6]
				+ rze2[2]*fsss[m1][6];
		gsps[m][ 6*3+1] = QC[1]*gsss[m][ 6] + WQ[1]*gsss[m1][ 6]
				+ rze2[2]*fsss[m1][3];
		gsps[m][ 6*3+2] = QC[2]*gsss[m][ 6] + WQ[2]*gsss[m1][ 6];
		gsps[m][ 7*3+0] = QC[0]*gsss[m][ 7] + WQ[0]*gsss[m1][ 7];
		gsps[m][ 7*3+1] = QC[1]*gsss[m][ 7] + WQ[1]*gsss[m1][ 7]
				+ rze2[2]*fsss[m1][7];
		gsps[m][ 7*3+2] = QC[2]*gsss[m][ 7] + WQ[2]*gsss[m1][ 7]
				+ rze2[2]*fsss[m1][4];
		gsps[m][ 8*3+0] = QC[0]*gsss[m][ 8] + WQ[0]*gsss[m1][ 8]
				+ rze2[2]*fsss[m1][5];
		gsps[m][ 8*3+1] = QC[1]*gsss[m][ 8] + WQ[1]*gsss[m1][ 8];
		gsps[m][ 8*3+2] = QC[2]*gsss[m][ 8] + WQ[2]*gsss[m1][ 8]
				+ rze2[2]*fsss[m1][8];
		gsps[m][ 9*3+0] = QC[0]*gsss[m][ 9] + WQ[0]*gsss[m1][ 9]
				+ rze2[1]*fsss[m1][1];
		gsps[m][ 9*3+1] = QC[1]*gsss[m][ 9] + WQ[1]*gsss[m1][ 9]
				+ rze2[3]*fsss[m1][6];
		gsps[m][ 9*3+2] = QC[2]*gsss[m][ 9] + WQ[2]*gsss[m1][ 9];
		gsps[m][10*3+0] = QC[0]*gsss[m][10] + WQ[0]*gsss[m1][10];
		gsps[m][10*3+1] = QC[1]*gsss[m][10] + WQ[1]*gsss[m1][10]
				+ rze2[1]*fsss[m1][2];
		gsps[m][10*3+2] = QC[2]*gsss[m][10] + WQ[2]*gsss[m1][10]
				+ rze2[3]*fsss[m1][7];
		gsps[m][11*3+0] = QC[0]*gsss[m][11] + WQ[0]*gsss[m1][11]
				+ rze2[3]*fsss[m1][8];
		gsps[m][11*3+1] = QC[1]*gsss[m][11] + WQ[1]*gsss[m1][11];
		gsps[m][11*3+2] = QC[2]*gsss[m][11] + WQ[2]*gsss[m1][11]
				+ rze2[1]*fsss[m1][0];
		gsps[m][12*3+0] = QC[0]*gsss[m][12] + WQ[0]*gsss[m1][12]
				+ rze2[2]*fsss[m1][9];
		gsps[m][12*3+1] = QC[1]*gsss[m][12] + WQ[1]*gsss[m1][12]
				+ rze2[1]*fsss[m1][8];
		gsps[m][12*3+2] = QC[2]*gsss[m][12] + WQ[2]*gsss[m1][12]
				+ rze2[1]*fsss[m1][3];
		gsps[m][13*3+0] = QC[0]*gsss[m][13] + WQ[0]*gsss[m1][13]
				+ rze2[1]*fsss[m1][4];
		gsps[m][13*3+1] = QC[1]*gsss[m][13] + WQ[1]*gsss[m1][13]
				+ rze2[2]*fsss[m1][9];
		gsps[m][13*3+2] = QC[2]*gsss[m][13] + WQ[2]*gsss[m1][13]
				+ rze2[1]*fsss[m1][6];
		gsps[m][14*3+0] = QC[0]*gsss[m][14] + WQ[0]*gsss[m1][14]
				+ rze2[1]*fsss[m1][7];
		gsps[m][14*3+1] = QC[1]*gsss[m][14] + WQ[1]*gsss[m1][14]
				+ rze2[1]*fsss[m1][5];
		gsps[m][14*3+2] = QC[2]*gsss[m][14] + WQ[2]*gsss[m1][14]
				+ rze2[2]*fsss[m1][9];
	    }	// end m loop
	    // (S,S|D,S)
	    tmp.ssss = eta2*( ssss[2] - re*ssss[3] );

	    ssds[0] = QC[0]*ssps[2-2][0] + WQ[0]*ssps[3-2][0]
		    +   tmp.ssss;
	    ssds[1] = QC[1]*ssps[2-2][1] + WQ[1]*ssps[3-2][1]
		    +   tmp.ssss;
	    ssds[2] = QC[2]*ssps[2-2][2] + WQ[2]*ssps[3-2][2]
		    +   tmp.ssss;
	    ssds[3] = QC[0]*ssps[2-2][1] + WQ[0]*ssps[3-2][1];
	    ssds[4] = QC[1]*ssps[2-2][2] + WQ[1]*ssps[3-2][2];
	    ssds[5] = QC[2]*ssps[2-2][0] + WQ[2]*ssps[3-2][0];
	    // (P,S|D,S)
	    for (int m=1; m<=2; m++) {
		int m1 = m+1;
		for (int i=0; i<3; i++)
		    tmp.psss[i] = eta2 * ( psss[m][i] - re*psss[m1][i] );

		psds[m-1][0*6+0] = QC[0]*psps[m-1][0*3+0] + WQ[0]*psps[m1-1][0*3+0]
				 +   tmp.psss[0] + rze2[1]*ssps[m1-2][0];
		psds[m-1][0*6+1] = QC[1]*psps[m-1][0*3+1] + WQ[1]*psps[m1-1][0*3+1]
				 +   tmp.psss[0];
		psds[m-1][0*6+2] = QC[2]*psps[m-1][0*3+2] + WQ[2]*psps[m1-1][0*3+2]
				 +   tmp.psss[0];
		psds[m-1][0*6+3] = QC[0]*psps[m-1][0*3+1] + WQ[0]*psps[m1-1][0*3+1]
				 + rze2[1]*ssps[m1-2][1];
		psds[m-1][0*6+4] = QC[1]*psps[m-1][0*3+2] + WQ[1]*psps[m1-1][0*3+2];
		psds[m-1][0*6+5] = QC[2]*psps[m-1][0*3+0] + WQ[2]*psps[m1-1][0*3+0];
		psds[m-1][1*6+0] = QC[0]*psps[m-1][1*3+0] + WQ[0]*psps[m1-1][1*3+0]
				 +   tmp.psss[1];
		psds[m-1][1*6+1] = QC[1]*psps[m-1][1*3+1] + WQ[1]*psps[m1-1][1*3+1]
				 +   tmp.psss[1] + rze2[1]*ssps[m1-2][1];
		psds[m-1][1*6+2] = QC[2]*psps[m-1][1*3+2] + WQ[2]*psps[m1-1][1*3+2]
				 +   tmp.psss[1];
		psds[m-1][1*6+3] = QC[0]*psps[m-1][1*3+1] + WQ[0]*psps[m1-1][1*3+1];
		psds[m-1][1*6+4] = QC[1]*psps[m-1][1*3+2] + WQ[1]*psps[m1-1][1*3+2]
				 + rze2[1]*ssps[m1-2][2];
		psds[m-1][1*6+5] = QC[2]*psps[m-1][1*3+0] + WQ[2]*psps[m1-1][1*3+0];
		psds[m-1][2*6+0] = QC[0]*psps[m-1][2*3+0] + WQ[0]*psps[m1-1][2*3+0]
				 +   tmp.psss[2];
		psds[m-1][2*6+1] = QC[1]*psps[m-1][2*3+1] + WQ[1]*psps[m1-1][2*3+1]
				 +   tmp.psss[2];
		psds[m-1][2*6+2] = QC[2]*psps[m-1][2*3+2] + WQ[2]*psps[m1-1][2*3+2]
				 +   tmp.psss[2] + rze2[1]*ssps[m1-2][2];
		psds[m-1][2*6+3] = QC[0]*psps[m-1][2*3+1] + WQ[0]*psps[m1-1][2*3+1];
		psds[m-1][2*6+4] = QC[1]*psps[m-1][2*3+2] + WQ[1]*psps[m1-1][2*3+2];
		psds[m-1][2*6+5] = QC[2]*psps[m-1][2*3+0] + WQ[2]*psps[m1-1][2*3+0]
				 + rze2[1]*ssps[m1-2][0];
	    }	// end m loop
	    // (D,S|D,S)
	    for (int m=0; m<=2; m++) {
		int m1 = m+1;
		for (int i=0; i<6; i++)
		    tmp.dsss[i] = eta2 * ( dsss[m][i] - re*dsss[m1][i] );

		dsds[m][0*6+0] = QC[0]*dsps[m][0*3+0] + WQ[0]*dsps[m1][0*3+0]
			       +   tmp.dsss[0] + rze2[2]*psps[m1-1][0*3+0];
		dsds[m][0*6+1] = QC[1]*dsps[m][0*3+1] + WQ[1]*dsps[m1][0*3+1]
			       +   tmp.dsss[0];
		dsds[m][0*6+2] = QC[2]*dsps[m][0*3+2] + WQ[2]*dsps[m1][0*3+2]
			       +   tmp.dsss[0];
		dsds[m][0*6+3] = QC[0]*dsps[m][0*3+1] + WQ[0]*dsps[m1][0*3+1]
			       + rze2[2]*psps[m1-1][0*3+1];
		dsds[m][0*6+4] = QC[1]*dsps[m][0*3+2] + WQ[1]*dsps[m1][0*3+2];
		dsds[m][0*6+5] = QC[2]*dsps[m][0*3+0] + WQ[2]*dsps[m1][0*3+0];
		dsds[m][1*6+0] = QC[0]*dsps[m][1*3+0] + WQ[0]*dsps[m1][1*3+0]
			       +   tmp.dsss[1];
		dsds[m][1*6+1] = QC[1]*dsps[m][1*3+1] + WQ[1]*dsps[m1][1*3+1]
			       +   tmp.dsss[1] + rze2[2]*psps[m1-1][1*3+1];
		dsds[m][1*6+2] = QC[2]*dsps[m][1*3+2] + WQ[2]*dsps[m1][1*3+2]
			       +   tmp.dsss[1];
		dsds[m][1*6+3] = QC[0]*dsps[m][1*3+1] + WQ[0]*dsps[m1][1*3+1];
		dsds[m][1*6+4] = QC[1]*dsps[m][1*3+2] + WQ[1]*dsps[m1][1*3+2]
			       + rze2[2]*psps[m1-1][1*3+2];
		dsds[m][1*6+5] = QC[2]*dsps[m][1*3+0] + WQ[2]*dsps[m1][1*3+0];
		dsds[m][2*6+0] = QC[0]*dsps[m][2*3+0] + WQ[0]*dsps[m1][2*3+0]
			       +   tmp.dsss[2];
		dsds[m][2*6+1] = QC[1]*dsps[m][2*3+1] + WQ[1]*dsps[m1][2*3+1]
			       +   tmp.dsss[2];
		dsds[m][2*6+2] = QC[2]*dsps[m][2*3+2] + WQ[2]*dsps[m1][2*3+2]
			       +   tmp.dsss[2] + rze2[2]*psps[m1-1][2*3+2];
		dsds[m][2*6+3] = QC[0]*dsps[m][2*3+1] + WQ[0]*dsps[m1][2*3+1];
		dsds[m][2*6+4] = QC[1]*dsps[m][2*3+2] + WQ[1]*dsps[m1][2*3+2];
		dsds[m][2*6+5] = QC[2]*dsps[m][2*3+0] + WQ[2]*dsps[m1][2*3+0]
			       + rze2[2]*psps[m1-1][2*3+0];
		dsds[m][3*6+0] = QC[0]*dsps[m][3*3+0] + WQ[0]*dsps[m1][3*3+0]
			       +   tmp.dsss[3] + rze2[1]*psps[m1-1][1*3+0];
		dsds[m][3*6+1] = QC[1]*dsps[m][3*3+1] + WQ[1]*dsps[m1][3*3+1]
			       +   tmp.dsss[3] + rze2[1]*psps[m1-1][0*3+1];
		dsds[m][3*6+2] = QC[2]*dsps[m][3*3+2] + WQ[2]*dsps[m1][3*3+2]
			       +   tmp.dsss[3];
		dsds[m][3*6+3] = QC[0]*dsps[m][3*3+1] + WQ[0]*dsps[m1][3*3+1]
			       + rze2[1]*psps[m1-1][1*3+1];
		dsds[m][3*6+4] = QC[1]*dsps[m][3*3+2] + WQ[1]*dsps[m1][3*3+2]
			       + rze2[1]*psps[m1-1][0*3+2];
		dsds[m][3*6+5] = QC[2]*dsps[m][3*3+0] + WQ[2]*dsps[m1][3*3+0];
		dsds[m][4*6+0] = QC[0]*dsps[m][4*3+0] + WQ[0]*dsps[m1][4*3+0]
			       +   tmp.dsss[4];
		dsds[m][4*6+1] = QC[1]*dsps[m][4*3+1] + WQ[1]*dsps[m1][4*3+1]
			       +   tmp.dsss[4] + rze2[1]*psps[m1-1][2*3+1];
		dsds[m][4*6+2] = QC[2]*dsps[m][4*3+2] + WQ[2]*dsps[m1][4*3+2]
			       +   tmp.dsss[4] + rze2[1]*psps[m1-1][1*3+2];
		dsds[m][4*6+3] = QC[0]*dsps[m][4*3+1] + WQ[0]*dsps[m1][4*3+1];
		dsds[m][4*6+4] = QC[1]*dsps[m][4*3+2] + WQ[1]*dsps[m1][4*3+2]
			       + rze2[1]*psps[m1-1][2*3+2];
		dsds[m][4*6+5] = QC[2]*dsps[m][4*3+0] + WQ[2]*dsps[m1][4*3+0]
			       + rze2[1]*psps[m1-1][1*3+0];
		dsds[m][5*6+0] = QC[0]*dsps[m][5*3+0] + WQ[0]*dsps[m1][5*3+0]
			       +   tmp.dsss[5] + rze2[1]*psps[m1-1][2*3+0];
		dsds[m][5*6+1] = QC[1]*dsps[m][5*3+1] + WQ[1]*dsps[m1][5*3+1]
			       +   tmp.dsss[5];
		dsds[m][5*6+2] = QC[2]*dsps[m][5*3+2] + WQ[2]*dsps[m1][5*3+2]
			       +   tmp.dsss[5] + rze2[1]*psps[m1-1][0*3+2];
		dsds[m][5*6+3] = QC[0]*dsps[m][5*3+1] + WQ[0]*dsps[m1][5*3+1]
			       + rze2[1]*psps[m1-1][2*3+1];
		dsds[m][5*6+4] = QC[1]*dsps[m][5*3+2] + WQ[1]*dsps[m1][5*3+2];
		dsds[m][5*6+5] = QC[2]*dsps[m][5*3+0] + WQ[2]*dsps[m1][5*3+0]
			       + rze2[1]*psps[m1-1][0*3+0];
	    }	// end m loop
	    // (F,S|D,S)
	    for (int m=0; m<=2; m++) {
		int m1 = m+1;
		for (int i=0; i<10; i++)
		    tmp.fsss[i] = eta2 * ( fsss[m][i] - re*fsss[m1][i] );

		fsds[m][0*6+0] = QC[0]*fsps[m][0*3+0] + WQ[0]*fsps[m1][0*3+0]
			       +   tmp.fsss[0] + rze2[3]*dsps[m1][0*3+0];
		fsds[m][0*6+1] = QC[1]*fsps[m][0*3+1] + WQ[1]*fsps[m1][0*3+1]
			       +   tmp.fsss[0];
		fsds[m][0*6+2] = QC[2]*fsps[m][0*3+2] + WQ[2]*fsps[m1][0*3+2]
			       +   tmp.fsss[0];
		fsds[m][0*6+3] = QC[0]*fsps[m][0*3+1] + WQ[0]*fsps[m1][0*3+1]
			       + rze2[3]*dsps[m1][0*3+1];
		fsds[m][0*6+4] = QC[1]*fsps[m][0*3+2] + WQ[1]*fsps[m1][0*3+2];
		fsds[m][0*6+5] = QC[2]*fsps[m][0*3+0] + WQ[2]*fsps[m1][0*3+0];
		fsds[m][1*6+0] = QC[0]*fsps[m][1*3+0] + WQ[0]*fsps[m1][1*3+0]
			       +   tmp.fsss[1];
		fsds[m][1*6+1] = QC[1]*fsps[m][1*3+1] + WQ[1]*fsps[m1][1*3+1]
			       +   tmp.fsss[1] + rze2[3]*dsps[m1][1*3+1];
		fsds[m][1*6+2] = QC[2]*fsps[m][1*3+2] + WQ[2]*fsps[m1][1*3+2]
			       +   tmp.fsss[1];
		fsds[m][1*6+3] = QC[0]*fsps[m][1*3+1] + WQ[0]*fsps[m1][1*3+1];
		fsds[m][1*6+4] = QC[1]*fsps[m][1*3+2] + WQ[1]*fsps[m1][1*3+2]
			       + rze2[3]*dsps[m1][1*3+2];
		fsds[m][1*6+5] = QC[2]*fsps[m][1*3+0] + WQ[2]*fsps[m1][1*3+0];
		fsds[m][2*6+0] = QC[0]*fsps[m][2*3+0] + WQ[0]*fsps[m1][2*3+0]
			       +   tmp.fsss[2];
		fsds[m][2*6+1] = QC[1]*fsps[m][2*3+1] + WQ[1]*fsps[m1][2*3+1]
			       +   tmp.fsss[2];
		fsds[m][2*6+2] = QC[2]*fsps[m][2*3+2] + WQ[2]*fsps[m1][2*3+2]
			       +   tmp.fsss[2] + rze2[3]*dsps[m1][2*3+2];
		fsds[m][2*6+3] = QC[0]*fsps[m][2*3+1] + WQ[0]*fsps[m1][2*3+1];
		fsds[m][2*6+4] = QC[1]*fsps[m][2*3+2] + WQ[1]*fsps[m1][2*3+2];
		fsds[m][2*6+5] = QC[2]*fsps[m][2*3+0] + WQ[2]*fsps[m1][2*3+0]
			       + rze2[3]*dsps[m1][2*3+0];
		fsds[m][3*6+0] = QC[0]*fsps[m][3*3+0] + WQ[0]*fsps[m1][3*3+0]
			       +   tmp.fsss[3] + rze2[2]*dsps[m1][3*3+0];
		fsds[m][3*6+1] = QC[1]*fsps[m][3*3+1] + WQ[1]*fsps[m1][3*3+1]
			       +   tmp.fsss[3] + rze2[1]*dsps[m1][0*3+1];
		fsds[m][3*6+2] = QC[2]*fsps[m][3*3+2] + WQ[2]*fsps[m1][3*3+2]
			       +   tmp.fsss[3];
		fsds[m][3*6+3] = QC[0]*fsps[m][3*3+1] + WQ[0]*fsps[m1][3*3+1]
			       + rze2[2]*dsps[m1][3*3+1];
		fsds[m][3*6+4] = QC[1]*fsps[m][3*3+2] + WQ[1]*fsps[m1][3*3+2]
			       + rze2[1]*dsps[m1][0*3+2];
		fsds[m][3*6+5] = QC[2]*fsps[m][3*3+0] + WQ[2]*fsps[m1][3*3+0];
		fsds[m][4*6+0] = QC[0]*fsps[m][4*3+0] + WQ[0]*fsps[m1][4*3+0]
			       +   tmp.fsss[4];
		fsds[m][4*6+1] = QC[1]*fsps[m][4*3+1] + WQ[1]*fsps[m1][4*3+1]
			       +   tmp.fsss[4] + rze2[2]*dsps[m1][4*3+1];
		fsds[m][4*6+2] = QC[2]*fsps[m][4*3+2] + WQ[2]*fsps[m1][4*3+2]
			       +   tmp.fsss[4] + rze2[1]*dsps[m1][1*3+2];
		fsds[m][4*6+3] = QC[0]*fsps[m][4*3+1] + WQ[0]*fsps[m1][4*3+1];
		fsds[m][4*6+4] = QC[1]*fsps[m][4*3+2] + WQ[1]*fsps[m1][4*3+2]
			       + rze2[2]*dsps[m1][4*3+2];
		fsds[m][4*6+5] = QC[2]*fsps[m][4*3+0] + WQ[2]*fsps[m1][4*3+0]
			       + rze2[1]*dsps[m1][1*3+0];
		fsds[m][5*6+0] = QC[0]*fsps[m][5*3+0] + WQ[0]*fsps[m1][5*3+0]
			       +   tmp.fsss[5] + rze2[1]*dsps[m1][2*3+0];
		fsds[m][5*6+1] = QC[1]*fsps[m][5*3+1] + WQ[1]*fsps[m1][5*3+1]
			       +   tmp.fsss[5];
		fsds[m][5*6+2] = QC[2]*fsps[m][5*3+2] + WQ[2]*fsps[m1][5*3+2]
			       +   tmp.fsss[5] + rze2[2]*dsps[m1][5*3+2];
		fsds[m][5*6+3] = QC[0]*fsps[m][5*3+1] + WQ[0]*fsps[m1][5*3+1]
			       + rze2[1]*dsps[m1][2*3+1];
		fsds[m][5*6+4] = QC[1]*fsps[m][5*3+2] + WQ[1]*fsps[m1][5*3+2];
		fsds[m][5*6+5] = QC[2]*fsps[m][5*3+0] + WQ[2]*fsps[m1][5*3+0]
			       + rze2[2]*dsps[m1][5*3+0];
		fsds[m][6*6+0] = QC[0]*fsps[m][6*3+0] + WQ[0]*fsps[m1][6*3+0]
			       +   tmp.fsss[6] + rze2[1]*dsps[m1][1*3+0];
		fsds[m][6*6+1] = QC[1]*fsps[m][6*3+1] + WQ[1]*fsps[m1][6*3+1]
			       +   tmp.fsss[6] + rze2[2]*dsps[m1][3*3+1];
		fsds[m][6*6+2] = QC[2]*fsps[m][6*3+2] + WQ[2]*fsps[m1][6*3+2]
			       +   tmp.fsss[6];
		fsds[m][6*6+3] = QC[0]*fsps[m][6*3+1] + WQ[0]*fsps[m1][6*3+1]
			       + rze2[1]*dsps[m1][1*3+1];
		fsds[m][6*6+4] = QC[1]*fsps[m][6*3+2] + WQ[1]*fsps[m1][6*3+2]
			       + rze2[2]*dsps[m1][3*3+2];
		fsds[m][6*6+5] = QC[2]*fsps[m][6*3+0] + WQ[2]*fsps[m1][6*3+0];
		fsds[m][7*6+0] = QC[0]*fsps[m][7*3+0] + WQ[0]*fsps[m1][7*3+0]
			       +   tmp.fsss[7];
		fsds[m][7*6+1] = QC[1]*fsps[m][7*3+1] + WQ[1]*fsps[m1][7*3+1]
			       +   tmp.fsss[7] + rze2[1]*dsps[m1][2*3+1];
		fsds[m][7*6+2] = QC[2]*fsps[m][7*3+2] + WQ[2]*fsps[m1][7*3+2]
			       +   tmp.fsss[7] + rze2[2]*dsps[m1][4*3+2];
		fsds[m][7*6+3] = QC[0]*fsps[m][7*3+1] + WQ[0]*fsps[m1][7*3+1];
		fsds[m][7*6+4] = QC[1]*fsps[m][7*3+2] + WQ[1]*fsps[m1][7*3+2]
			       + rze2[1]*dsps[m1][2*3+2];
		fsds[m][7*6+5] = QC[2]*fsps[m][7*3+0] + WQ[2]*fsps[m1][7*3+0]
			       + rze2[2]*dsps[m1][4*3+0];
		fsds[m][8*6+0] = QC[0]*fsps[m][8*3+0] + WQ[0]*fsps[m1][8*3+0]
			       +   tmp.fsss[8] + rze2[2]*dsps[m1][5*3+0];
		fsds[m][8*6+1] = QC[1]*fsps[m][8*3+1] + WQ[1]*fsps[m1][8*3+1]
			       +   tmp.fsss[8];
		fsds[m][8*6+2] = QC[2]*fsps[m][8*3+2] + WQ[2]*fsps[m1][8*3+2]
			       +   tmp.fsss[8] + rze2[1]*dsps[m1][0*3+2];
		fsds[m][8*6+3] = QC[0]*fsps[m][8*3+1] + WQ[0]*fsps[m1][8*3+1]
			       + rze2[2]*dsps[m1][5*3+1];
		fsds[m][8*6+4] = QC[1]*fsps[m][8*3+2] + WQ[1]*fsps[m1][8*3+2];
		fsds[m][8*6+5] = QC[2]*fsps[m][8*3+0] + WQ[2]*fsps[m1][8*3+0]
			       + rze2[1]*dsps[m1][0*3+0];
		fsds[m][9*6+0] = QC[0]*fsps[m][9*3+0] + WQ[0]*fsps[m1][9*3+0]
			       +   tmp.fsss[9] + rze2[1]*dsps[m1][4*3+0];
		fsds[m][9*6+1] = QC[1]*fsps[m][9*3+1] + WQ[1]*fsps[m1][9*3+1]
			       +   tmp.fsss[9] + rze2[1]*dsps[m1][5*3+1];
		fsds[m][9*6+2] = QC[2]*fsps[m][9*3+2] + WQ[2]*fsps[m1][9*3+2]
			       +   tmp.fsss[9] + rze2[1]*dsps[m1][3*3+2];
		fsds[m][9*6+3] = QC[0]*fsps[m][9*3+1] + WQ[0]*fsps[m1][9*3+1]
			       + rze2[1]*dsps[m1][4*3+1];
		fsds[m][9*6+4] = QC[1]*fsps[m][9*3+2] + WQ[1]*fsps[m1][9*3+2]
			       + rze2[1]*dsps[m1][5*3+2];
		fsds[m][9*6+5] = QC[2]*fsps[m][9*3+0] + WQ[2]*fsps[m1][9*3+0]
			       + rze2[1]*dsps[m1][3*3+0];
	    }	// end m loop
	    // (G,S|D,S)
	    for (int m=0; m<=2; m++) {
		int m1 = m+1;
		for (int i=0; i<15; i++)
		    tmp.gsss[i] = eta2 * ( gsss[m][i] - re*gsss[m1][i] );

		gsds[m][ 0*6+0] = QC[0]*gsps[m][ 0*3+0] + WQ[0]*gsps[m1][ 0*3+0]
				+   tmp.gsss[ 0] + rze2[4]*fsps[m1][0*3+0];
		gsds[m][ 0*6+1] = QC[1]*gsps[m][ 0*3+1] + WQ[1]*gsps[m1][ 0*3+1]
				+   tmp.gsss[ 0];
		gsds[m][ 0*6+2] = QC[2]*gsps[m][ 0*3+2] + WQ[2]*gsps[m1][ 0*3+2]
				+   tmp.gsss[ 0];
		gsds[m][ 0*6+3] = QC[0]*gsps[m][ 0*3+1] + WQ[0]*gsps[m1][ 0*3+1]
				+ rze2[4]*fsps[m1][0*3+1];
		gsds[m][ 0*6+4] = QC[1]*gsps[m][ 0*3+2] + WQ[1]*gsps[m1][ 0*3+2];
		gsds[m][ 0*6+5] = QC[2]*gsps[m][ 0*3+0] + WQ[2]*gsps[m1][ 0*3+0];
		gsds[m][ 1*6+0] = QC[0]*gsps[m][ 1*3+0] + WQ[0]*gsps[m1][ 1*3+0]
				+   tmp.gsss[ 1];
		gsds[m][ 1*6+1] = QC[1]*gsps[m][ 1*3+1] + WQ[1]*gsps[m1][ 1*3+1]
				+   tmp.gsss[ 1] + rze2[4]*fsps[m1][1*3+1];
		gsds[m][ 1*6+2] = QC[2]*gsps[m][ 1*3+2] + WQ[2]*gsps[m1][ 1*3+2]
				+   tmp.gsss[ 1];
		gsds[m][ 1*6+3] = QC[0]*gsps[m][ 1*3+1] + WQ[0]*gsps[m1][ 1*3+1];
		gsds[m][ 1*6+4] = QC[1]*gsps[m][ 1*3+2] + WQ[1]*gsps[m1][ 1*3+2]
				+ rze2[4]*fsps[m1][1*3+2];
		gsds[m][ 1*6+5] = QC[2]*gsps[m][ 1*3+0] + WQ[2]*gsps[m1][ 1*3+0];
		gsds[m][ 2*6+0] = QC[0]*gsps[m][ 2*3+0] + WQ[0]*gsps[m1][ 2*3+0]
				+   tmp.gsss[ 2];
		gsds[m][ 2*6+1] = QC[1]*gsps[m][ 2*3+1] + WQ[1]*gsps[m1][ 2*3+1]
				+   tmp.gsss[ 2];
		gsds[m][ 2*6+2] = QC[2]*gsps[m][ 2*3+2] + WQ[2]*gsps[m1][ 2*3+2]
				+   tmp.gsss[ 2] + rze2[4]*fsps[m1][2*3+2];
		gsds[m][ 2*6+3] = QC[0]*gsps[m][ 2*3+1] + WQ[0]*gsps[m1][ 2*3+1];
		gsds[m][ 2*6+4] = QC[1]*gsps[m][ 2*3+2] + WQ[1]*gsps[m1][ 2*3+2];
		gsds[m][ 2*6+5] = QC[2]*gsps[m][ 2*3+0] + WQ[2]*gsps[m1][ 2*3+0]
				+ rze2[4]*fsps[m1][2*3+0];
		gsds[m][ 3*6+0] = QC[0]*gsps[m][ 3*3+0] + WQ[0]*gsps[m1][ 3*3+0]
				+   tmp.gsss[ 3] + rze2[3]*fsps[m1][3*3+0];
		gsds[m][ 3*6+1] = QC[1]*gsps[m][ 3*3+1] + WQ[1]*gsps[m1][ 3*3+1]
				+   tmp.gsss[ 3] + rze2[1]*fsps[m1][0*3+1];
		gsds[m][ 3*6+2] = QC[2]*gsps[m][ 3*3+2] + WQ[2]*gsps[m1][ 3*3+2]
				+   tmp.gsss[ 3];
		gsds[m][ 3*6+3] = QC[0]*gsps[m][ 3*3+1] + WQ[0]*gsps[m1][ 3*3+1]
				+ rze2[3]*fsps[m1][3*3+1];
		gsds[m][ 3*6+4] = QC[1]*gsps[m][ 3*3+2] + WQ[1]*gsps[m1][ 3*3+2]
				+ rze2[1]*fsps[m1][0*3+2];
		gsds[m][ 3*6+5] = QC[2]*gsps[m][ 3*3+0] + WQ[2]*gsps[m1][ 3*3+0];
		gsds[m][ 4*6+0] = QC[0]*gsps[m][ 4*3+0] + WQ[0]*gsps[m1][ 4*3+0]
				+   tmp.gsss[ 4];
		gsds[m][ 4*6+1] = QC[1]*gsps[m][ 4*3+1] + WQ[1]*gsps[m1][ 4*3+1]
				+   tmp.gsss[ 4] + rze2[3]*fsps[m1][4*3+1];
		gsds[m][ 4*6+2] = QC[2]*gsps[m][ 4*3+2] + WQ[2]*gsps[m1][ 4*3+2]
				+   tmp.gsss[ 4] + rze2[1]*fsps[m1][1*3+2];
		gsds[m][ 4*6+3] = QC[0]*gsps[m][ 4*3+1] + WQ[0]*gsps[m1][ 4*3+1];
		gsds[m][ 4*6+4] = QC[1]*gsps[m][ 4*3+2] + WQ[1]*gsps[m1][ 4*3+2]
				+ rze2[3]*fsps[m1][4*3+2];
		gsds[m][ 4*6+5] = QC[2]*gsps[m][ 4*3+0] + WQ[2]*gsps[m1][ 4*3+0]
				+ rze2[1]*fsps[m1][1*3+0];
		gsds[m][ 5*6+0] = QC[0]*gsps[m][ 5*3+0] + WQ[0]*gsps[m1][ 5*3+0]
				+   tmp.gsss[ 5] + rze2[1]*fsps[m1][2*3+0];
		gsds[m][ 5*6+1] = QC[1]*gsps[m][ 5*3+1] + WQ[1]*gsps[m1][ 5*3+1]
				+   tmp.gsss[ 5];
		gsds[m][ 5*6+2] = QC[2]*gsps[m][ 5*3+2] + WQ[2]*gsps[m1][ 5*3+2]
				+   tmp.gsss[ 5] + rze2[3]*fsps[m1][5*3+2];
		gsds[m][ 5*6+3] = QC[0]*gsps[m][ 5*3+1] + WQ[0]*gsps[m1][ 5*3+1]
				+ rze2[1]*fsps[m1][2*3+1];
		gsds[m][ 5*6+4] = QC[1]*gsps[m][ 5*3+2] + WQ[1]*gsps[m1][ 5*3+2];
		gsds[m][ 5*6+5] = QC[2]*gsps[m][ 5*3+0] + WQ[2]*gsps[m1][ 5*3+0]
				+ rze2[3]*fsps[m1][5*3+0];
		gsds[m][ 6*6+0] = QC[0]*gsps[m][ 6*3+0] + WQ[0]*gsps[m1][ 6*3+0]
				+   tmp.gsss[ 6] + rze2[2]*fsps[m1][6*3+0];
		gsds[m][ 6*6+1] = QC[1]*gsps[m][ 6*3+1] + WQ[1]*gsps[m1][ 6*3+1]
				+   tmp.gsss[ 6] + rze2[2]*fsps[m1][3*3+1];
		gsds[m][ 6*6+2] = QC[2]*gsps[m][ 6*3+2] + WQ[2]*gsps[m1][ 6*3+2]
				+   tmp.gsss[ 6];
		gsds[m][ 6*6+3] = QC[0]*gsps[m][ 6*3+1] + WQ[0]*gsps[m1][ 6*3+1]
				+ rze2[2]*fsps[m1][6*3+1];
		gsds[m][ 6*6+4] = QC[1]*gsps[m][ 6*3+2] + WQ[1]*gsps[m1][ 6*3+2]
				+ rze2[2]*fsps[m1][3*3+2];
		gsds[m][ 6*6+5] = QC[2]*gsps[m][ 6*3+0] + WQ[2]*gsps[m1][ 6*3+0];
		gsds[m][ 7*6+0] = QC[0]*gsps[m][ 7*3+0] + WQ[0]*gsps[m1][ 7*3+0]
				+   tmp.gsss[ 7];
		gsds[m][ 7*6+1] = QC[1]*gsps[m][ 7*3+1] + WQ[1]*gsps[m1][ 7*3+1]
				+   tmp.gsss[ 7] + rze2[2]*fsps[m1][7*3+1];
		gsds[m][ 7*6+2] = QC[2]*gsps[m][ 7*3+2] + WQ[2]*gsps[m1][ 7*3+2]
				+   tmp.gsss[ 7] + rze2[2]*fsps[m1][4*3+2];
		gsds[m][ 7*6+3] = QC[0]*gsps[m][ 7*3+1] + WQ[0]*gsps[m1][ 7*3+1];
		gsds[m][ 7*6+4] = QC[1]*gsps[m][ 7*3+2] + WQ[1]*gsps[m1][ 7*3+2]
				+ rze2[2]*fsps[m1][7*3+2];
		gsds[m][ 7*6+5] = QC[2]*gsps[m][ 7*3+0] + WQ[2]*gsps[m1][ 7*3+0]
				+ rze2[2]*fsps[m1][4*3+0];
		gsds[m][ 8*6+0] = QC[0]*gsps[m][ 8*3+0] + WQ[0]*gsps[m1][ 8*3+0]
				+   tmp.gsss[ 8] + rze2[2]*fsps[m1][5*3+0];
		gsds[m][ 8*6+1] = QC[1]*gsps[m][ 8*3+1] + WQ[1]*gsps[m1][ 8*3+1]
				+   tmp.gsss[ 8];
		gsds[m][ 8*6+2] = QC[2]*gsps[m][ 8*3+2] + WQ[2]*gsps[m1][ 8*3+2]
				+   tmp.gsss[ 8] + rze2[2]*fsps[m1][8*3+2];
		gsds[m][ 8*6+3] = QC[0]*gsps[m][ 8*3+1] + WQ[0]*gsps[m1][ 8*3+1]
				+ rze2[2]*fsps[m1][5*3+1];
		gsds[m][ 8*6+4] = QC[1]*gsps[m][ 8*3+2] + WQ[1]*gsps[m1][ 8*3+2];
		gsds[m][ 8*6+5] = QC[2]*gsps[m][ 8*3+0] + WQ[2]*gsps[m1][ 8*3+0]
				+ rze2[2]*fsps[m1][8*3+0];
		gsds[m][ 9*6+0] = QC[0]*gsps[m][ 9*3+0] + WQ[0]*gsps[m1][ 9*3+0]
				+   tmp.gsss[ 9] + rze2[1]*fsps[m1][1*3+0];
		gsds[m][ 9*6+1] = QC[1]*gsps[m][ 9*3+1] + WQ[1]*gsps[m1][ 9*3+1]
				+   tmp.gsss[ 9] + rze2[3]*fsps[m1][6*3+1];
		gsds[m][ 9*6+2] = QC[2]*gsps[m][ 9*3+2] + WQ[2]*gsps[m1][ 9*3+2]
				+   tmp.gsss[ 9];
		gsds[m][ 9*6+3] = QC[0]*gsps[m][ 9*3+1] + WQ[0]*gsps[m1][ 9*3+1]
				+ rze2[1]*fsps[m1][1*3+1];
		gsds[m][ 9*6+4] = QC[1]*gsps[m][ 9*3+2] + WQ[1]*gsps[m1][ 9*3+2]
				+ rze2[3]*fsps[m1][6*3+2];
		gsds[m][ 9*6+5] = QC[2]*gsps[m][ 9*3+0] + WQ[2]*gsps[m1][ 9*3+0];
		gsds[m][10*6+0] = QC[0]*gsps[m][10*3+0] + WQ[0]*gsps[m1][10*3+0]
				+   tmp.gsss[10];
		gsds[m][10*6+1] = QC[1]*gsps[m][10*3+1] + WQ[1]*gsps[m1][10*3+1]
				+   tmp.gsss[10] + rze2[1]*fsps[m1][2*3+1];
		gsds[m][10*6+2] = QC[2]*gsps[m][10*3+2] + WQ[2]*gsps[m1][10*3+2]
				+   tmp.gsss[10] + rze2[3]*fsps[m1][7*3+2];
		gsds[m][10*6+3] = QC[0]*gsps[m][10*3+1] + WQ[0]*gsps[m1][10*3+1];
		gsds[m][10*6+4] = QC[1]*gsps[m][10*3+2] + WQ[1]*gsps[m1][10*3+2]
				+ rze2[1]*fsps[m1][2*3+2];
		gsds[m][10*6+5] = QC[2]*gsps[m][10*3+0] + WQ[2]*gsps[m1][10*3+0]
				+ rze2[3]*fsps[m1][7*3+0];
		gsds[m][11*6+0] = QC[0]*gsps[m][11*3+0] + WQ[0]*gsps[m1][11*3+0]
				+   tmp.gsss[11] + rze2[3]*fsps[m1][8*3+0];
		gsds[m][11*6+1] = QC[1]*gsps[m][11*3+1] + WQ[1]*gsps[m1][11*3+1]
				+   tmp.gsss[11];
		gsds[m][11*6+2] = QC[2]*gsps[m][11*3+2] + WQ[2]*gsps[m1][11*3+2]
				+   tmp.gsss[11] + rze2[1]*fsps[m1][0*3+2];
		gsds[m][11*6+3] = QC[0]*gsps[m][11*3+1] + WQ[0]*gsps[m1][11*3+1]
				+ rze2[3]*fsps[m1][8*3+1];
		gsds[m][11*6+4] = QC[1]*gsps[m][11*3+2] + WQ[1]*gsps[m1][11*3+2];
		gsds[m][11*6+5] = QC[2]*gsps[m][11*3+0] + WQ[2]*gsps[m1][11*3+0]
				+ rze2[1]*fsps[m1][0*3+0];
		gsds[m][12*6+0] = QC[0]*gsps[m][12*3+0] + WQ[0]*gsps[m1][12*3+0]
				+   tmp.gsss[12] + rze2[2]*fsps[m1][9*3+0];
		gsds[m][12*6+1] = QC[1]*gsps[m][12*3+1] + WQ[1]*gsps[m1][12*3+1]
				+   tmp.gsss[12] + rze2[1]*fsps[m1][8*3+1];
		gsds[m][12*6+2] = QC[2]*gsps[m][12*3+2] + WQ[2]*gsps[m1][12*3+2]
				+   tmp.gsss[12] + rze2[1]*fsps[m1][3*3+2];
		gsds[m][12*6+3] = QC[0]*gsps[m][12*3+1] + WQ[0]*gsps[m1][12*3+1]
				+ rze2[2]*fsps[m1][9*3+1];
		gsds[m][12*6+4] = QC[1]*gsps[m][12*3+2] + WQ[1]*gsps[m1][12*3+2]
				+ rze2[1]*fsps[m1][8*3+2];
		gsds[m][12*6+5] = QC[2]*gsps[m][12*3+0] + WQ[2]*gsps[m1][12*3+0]
				+ rze2[1]*fsps[m1][3*3+0];
		gsds[m][13*6+0] = QC[0]*gsps[m][13*3+0] + WQ[0]*gsps[m1][13*3+0]
				+   tmp.gsss[13] + rze2[1]*fsps[m1][4*3+0];
		gsds[m][13*6+1] = QC[1]*gsps[m][13*3+1] + WQ[1]*gsps[m1][13*3+1]
				+   tmp.gsss[13] + rze2[2]*fsps[m1][9*3+1];
		gsds[m][13*6+2] = QC[2]*gsps[m][13*3+2] + WQ[2]*gsps[m1][13*3+2]
				+   tmp.gsss[13] + rze2[1]*fsps[m1][6*3+2];
		gsds[m][13*6+3] = QC[0]*gsps[m][13*3+1] + WQ[0]*gsps[m1][13*3+1]
				+ rze2[1]*fsps[m1][4*3+1];
		gsds[m][13*6+4] = QC[1]*gsps[m][13*3+2] + WQ[1]*gsps[m1][13*3+2]
				+ rze2[2]*fsps[m1][9*3+2];
		gsds[m][13*6+5] = QC[2]*gsps[m][13*3+0] + WQ[2]*gsps[m1][13*3+0]
				+ rze2[1]*fsps[m1][6*3+0];
		gsds[m][14*6+0] = QC[0]*gsps[m][14*3+0] + WQ[0]*gsps[m1][14*3+0]
				+   tmp.gsss[14] + rze2[1]*fsps[m1][7*3+0];
		gsds[m][14*6+1] = QC[1]*gsps[m][14*3+1] + WQ[1]*gsps[m1][14*3+1]
				+   tmp.gsss[14] + rze2[1]*fsps[m1][5*3+1];
		gsds[m][14*6+2] = QC[2]*gsps[m][14*3+2] + WQ[2]*gsps[m1][14*3+2]
				+   tmp.gsss[14] + rze2[2]*fsps[m1][9*3+2];
		gsds[m][14*6+3] = QC[0]*gsps[m][14*3+1] + WQ[0]*gsps[m1][14*3+1]
				+ rze2[1]*fsps[m1][7*3+1];
		gsds[m][14*6+4] = QC[1]*gsps[m][14*3+2] + WQ[1]*gsps[m1][14*3+2]
				+ rze2[1]*fsps[m1][5*3+2];
		gsds[m][14*6+5] = QC[2]*gsps[m][14*3+0] + WQ[2]*gsps[m1][14*3+0]
				+ rze2[2]*fsps[m1][9*3+0];
	    }	// end m loop
	    // (P,S|F,S)
	    for (int i=0; i<3*3; i++)
		tmp.psps[i] = eta2 * ( psps[1-1][i] - re*psps[2-1][i] );

	    psfs[0*10+0] = QC[0]*psds[1-1][0*6+0] + WQ[0]*psds[2-1][0*6+0]
			 + 2*tmp.psps[0*3+0] + rze2[1]*ssds[0];
	    psfs[0*10+1] = QC[1]*psds[1-1][0*6+1] + WQ[1]*psds[2-1][0*6+1]
			 + 2*tmp.psps[0*3+1];
	    psfs[0*10+2] = QC[2]*psds[1-1][0*6+2] + WQ[2]*psds[2-1][0*6+2]
			 + 2*tmp.psps[0*3+2];
	    psfs[0*10+3] = QC[0]*psds[1-1][0*6+3] + WQ[0]*psds[2-1][0*6+3]
			 +   tmp.psps[0*3+1] + rze2[1]*ssds[3];
	    psfs[0*10+4] = QC[1]*psds[1-1][0*6+4] + WQ[1]*psds[2-1][0*6+4]
			 +   tmp.psps[0*3+2];
	    psfs[0*10+5] = QC[2]*psds[1-1][0*6+5] + WQ[2]*psds[2-1][0*6+5]
			 +   tmp.psps[0*3+0];
	    psfs[0*10+6] = QC[0]*psds[1-1][0*6+1] + WQ[0]*psds[2-1][0*6+1]
			 + rze2[1]*ssds[1];
	    psfs[0*10+7] = QC[1]*psds[1-1][0*6+2] + WQ[1]*psds[2-1][0*6+2];
	    psfs[0*10+8] = QC[2]*psds[1-1][0*6+0] + WQ[2]*psds[2-1][0*6+0];
	    psfs[0*10+9] = QC[0]*psds[1-1][0*6+4] + WQ[0]*psds[2-1][0*6+4]
			 + rze2[1]*ssds[4];
	    psfs[1*10+0] = QC[0]*psds[1-1][1*6+0] + WQ[0]*psds[2-1][1*6+0]
			 + 2*tmp.psps[1*3+0];
	    psfs[1*10+1] = QC[1]*psds[1-1][1*6+1] + WQ[1]*psds[2-1][1*6+1]
			 + 2*tmp.psps[1*3+1] + rze2[1]*ssds[1];
	    psfs[1*10+2] = QC[2]*psds[1-1][1*6+2] + WQ[2]*psds[2-1][1*6+2]
			 + 2*tmp.psps[1*3+2];
	    psfs[1*10+3] = QC[0]*psds[1-1][1*6+3] + WQ[0]*psds[2-1][1*6+3]
			 +   tmp.psps[1*3+1];
	    psfs[1*10+4] = QC[1]*psds[1-1][1*6+4] + WQ[1]*psds[2-1][1*6+4]
			 +   tmp.psps[1*3+2] + rze2[1]*ssds[4];
	    psfs[1*10+5] = QC[2]*psds[1-1][1*6+5] + WQ[2]*psds[2-1][1*6+5]
			 +   tmp.psps[1*3+0];
	    psfs[1*10+6] = QC[0]*psds[1-1][1*6+1] + WQ[0]*psds[2-1][1*6+1];
	    psfs[1*10+7] = QC[1]*psds[1-1][1*6+2] + WQ[1]*psds[2-1][1*6+2]
			 + rze2[1]*ssds[2];
	    psfs[1*10+8] = QC[2]*psds[1-1][1*6+0] + WQ[2]*psds[2-1][1*6+0];
	    psfs[1*10+9] = QC[0]*psds[1-1][1*6+4] + WQ[0]*psds[2-1][1*6+4];
	    psfs[2*10+0] = QC[0]*psds[1-1][2*6+0] + WQ[0]*psds[2-1][2*6+0]
			 + 2*tmp.psps[2*3+0];
	    psfs[2*10+1] = QC[1]*psds[1-1][2*6+1] + WQ[1]*psds[2-1][2*6+1]
			 + 2*tmp.psps[2*3+1];
	    psfs[2*10+2] = QC[2]*psds[1-1][2*6+2] + WQ[2]*psds[2-1][2*6+2]
			 + 2*tmp.psps[2*3+2] + rze2[1]*ssds[2];
	    psfs[2*10+3] = QC[0]*psds[1-1][2*6+3] + WQ[0]*psds[2-1][2*6+3]
			 +   tmp.psps[2*3+1];
	    psfs[2*10+4] = QC[1]*psds[1-1][2*6+4] + WQ[1]*psds[2-1][2*6+4]
			 +   tmp.psps[2*3+2];
	    psfs[2*10+5] = QC[2]*psds[1-1][2*6+5] + WQ[2]*psds[2-1][2*6+5]
			 +   tmp.psps[2*3+0] + rze2[1]*ssds[5];
	    psfs[2*10+6] = QC[0]*psds[1-1][2*6+1] + WQ[0]*psds[2-1][2*6+1];
	    psfs[2*10+7] = QC[1]*psds[1-1][2*6+2] + WQ[1]*psds[2-1][2*6+2];
	    psfs[2*10+8] = QC[2]*psds[1-1][2*6+0] + WQ[2]*psds[2-1][2*6+0]
			 + rze2[1]*ssds[0];
	    psfs[2*10+9] = QC[0]*psds[1-1][2*6+4] + WQ[0]*psds[2-1][2*6+4];
	    // (D,S|F,S)
	    for (int m=0; m<=1; m++) {
		int m1 = m+1;
		for (int i=0; i<6*3; i++)
		    tmp.dsps[i] = eta2 * ( dsps[m][i] - re*dsps[m1][i] );

		dsfs[m][0*10+0] = QC[0]*dsds[m][0*6+0] + WQ[0]*dsds[m1][0*6+0]
				+ 2*tmp.dsps[0*3+0] + rze2[2]*psds[m1-1][0*6+0];
		dsfs[m][0*10+1] = QC[1]*dsds[m][0*6+1] + WQ[1]*dsds[m1][0*6+1]
				+ 2*tmp.dsps[0*3+1];
		dsfs[m][0*10+2] = QC[2]*dsds[m][0*6+2] + WQ[2]*dsds[m1][0*6+2]
				+ 2*tmp.dsps[0*3+2];
		dsfs[m][0*10+3] = QC[0]*dsds[m][0*6+3] + WQ[0]*dsds[m1][0*6+3]
				+   tmp.dsps[0*3+1] + rze2[2]*psds[m1-1][0*6+3];
		dsfs[m][0*10+4] = QC[1]*dsds[m][0*6+4] + WQ[1]*dsds[m1][0*6+4]
				+   tmp.dsps[0*3+2];
		dsfs[m][0*10+5] = QC[2]*dsds[m][0*6+5] + WQ[2]*dsds[m1][0*6+5]
				+   tmp.dsps[0*3+0];
		dsfs[m][0*10+6] = QC[0]*dsds[m][0*6+1] + WQ[0]*dsds[m1][0*6+1]
				+ rze2[2]*psds[m1-1][0*6+1];
		dsfs[m][0*10+7] = QC[1]*dsds[m][0*6+2] + WQ[1]*dsds[m1][0*6+2];
		dsfs[m][0*10+8] = QC[2]*dsds[m][0*6+0] + WQ[2]*dsds[m1][0*6+0];
		dsfs[m][0*10+9] = QC[0]*dsds[m][0*6+4] + WQ[0]*dsds[m1][0*6+4]
				+ rze2[2]*psds[m1-1][0*6+4];
		dsfs[m][1*10+0] = QC[0]*dsds[m][1*6+0] + WQ[0]*dsds[m1][1*6+0]
				+ 2*tmp.dsps[1*3+0];
		dsfs[m][1*10+1] = QC[1]*dsds[m][1*6+1] + WQ[1]*dsds[m1][1*6+1]
				+ 2*tmp.dsps[1*3+1] + rze2[2]*psds[m1-1][1*6+1];
		dsfs[m][1*10+2] = QC[2]*dsds[m][1*6+2] + WQ[2]*dsds[m1][1*6+2]
				+ 2*tmp.dsps[1*3+2];
		dsfs[m][1*10+3] = QC[0]*dsds[m][1*6+3] + WQ[0]*dsds[m1][1*6+3]
				+   tmp.dsps[1*3+1];
		dsfs[m][1*10+4] = QC[1]*dsds[m][1*6+4] + WQ[1]*dsds[m1][1*6+4]
				+   tmp.dsps[1*3+2] + rze2[2]*psds[m1-1][1*6+4];
		dsfs[m][1*10+5] = QC[2]*dsds[m][1*6+5] + WQ[2]*dsds[m1][1*6+5]
				+   tmp.dsps[1*3+0];
		dsfs[m][1*10+6] = QC[0]*dsds[m][1*6+1] + WQ[0]*dsds[m1][1*6+1];
		dsfs[m][1*10+7] = QC[1]*dsds[m][1*6+2] + WQ[1]*dsds[m1][1*6+2]
				+ rze2[2]*psds[m1-1][1*6+2];
		dsfs[m][1*10+8] = QC[2]*dsds[m][1*6+0] + WQ[2]*dsds[m1][1*6+0];
		dsfs[m][1*10+9] = QC[0]*dsds[m][1*6+4] + WQ[0]*dsds[m1][1*6+4];
		dsfs[m][2*10+0] = QC[0]*dsds[m][2*6+0] + WQ[0]*dsds[m1][2*6+0]
				+ 2*tmp.dsps[2*3+0];
		dsfs[m][2*10+1] = QC[1]*dsds[m][2*6+1] + WQ[1]*dsds[m1][2*6+1]
				+ 2*tmp.dsps[2*3+1];
		dsfs[m][2*10+2] = QC[2]*dsds[m][2*6+2] + WQ[2]*dsds[m1][2*6+2]
				+ 2*tmp.dsps[2*3+2] + rze2[2]*psds[m1-1][2*6+2];
		dsfs[m][2*10+3] = QC[0]*dsds[m][2*6+3] + WQ[0]*dsds[m1][2*6+3]
				+   tmp.dsps[2*3+1];
		dsfs[m][2*10+4] = QC[1]*dsds[m][2*6+4] + WQ[1]*dsds[m1][2*6+4]
				+   tmp.dsps[2*3+2];
		dsfs[m][2*10+5] = QC[2]*dsds[m][2*6+5] + WQ[2]*dsds[m1][2*6+5]
				+   tmp.dsps[2*3+0] + rze2[2]*psds[m1-1][2*6+5];
		dsfs[m][2*10+6] = QC[0]*dsds[m][2*6+1] + WQ[0]*dsds[m1][2*6+1];
		dsfs[m][2*10+7] = QC[1]*dsds[m][2*6+2] + WQ[1]*dsds[m1][2*6+2];
		dsfs[m][2*10+8] = QC[2]*dsds[m][2*6+0] + WQ[2]*dsds[m1][2*6+0]
				+ rze2[2]*psds[m1-1][2*6+0];
		dsfs[m][2*10+9] = QC[0]*dsds[m][2*6+4] + WQ[0]*dsds[m1][2*6+4];
		dsfs[m][3*10+0] = QC[0]*dsds[m][3*6+0] + WQ[0]*dsds[m1][3*6+0]
				+ 2*tmp.dsps[3*3+0] + rze2[1]*psds[m1-1][1*6+0];
		dsfs[m][3*10+1] = QC[1]*dsds[m][3*6+1] + WQ[1]*dsds[m1][3*6+1]
				+ 2*tmp.dsps[3*3+1] + rze2[1]*psds[m1-1][0*6+1];
		dsfs[m][3*10+2] = QC[2]*dsds[m][3*6+2] + WQ[2]*dsds[m1][3*6+2]
				+ 2*tmp.dsps[3*3+2];
		dsfs[m][3*10+3] = QC[0]*dsds[m][3*6+3] + WQ[0]*dsds[m1][3*6+3]
				+   tmp.dsps[3*3+1] + rze2[1]*psds[m1-1][1*6+3];
		dsfs[m][3*10+4] = QC[1]*dsds[m][3*6+4] + WQ[1]*dsds[m1][3*6+4]
				+   tmp.dsps[3*3+2] + rze2[1]*psds[m1-1][0*6+4];
		dsfs[m][3*10+5] = QC[2]*dsds[m][3*6+5] + WQ[2]*dsds[m1][3*6+5]
				+   tmp.dsps[3*3+0];
		dsfs[m][3*10+6] = QC[0]*dsds[m][3*6+1] + WQ[0]*dsds[m1][3*6+1]
				+ rze2[1]*psds[m1-1][1*6+1];
		dsfs[m][3*10+7] = QC[1]*dsds[m][3*6+2] + WQ[1]*dsds[m1][3*6+2]
				+ rze2[1]*psds[m1-1][0*6+2];
		dsfs[m][3*10+8] = QC[2]*dsds[m][3*6+0] + WQ[2]*dsds[m1][3*6+0];
		dsfs[m][3*10+9] = QC[0]*dsds[m][3*6+4] + WQ[0]*dsds[m1][3*6+4]
				+ rze2[1]*psds[m1-1][1*6+4];
		dsfs[m][4*10+0] = QC[0]*dsds[m][4*6+0] + WQ[0]*dsds[m1][4*6+0]
				+ 2*tmp.dsps[4*3+0];
		dsfs[m][4*10+1] = QC[1]*dsds[m][4*6+1] + WQ[1]*dsds[m1][4*6+1]
				+ 2*tmp.dsps[4*3+1] + rze2[1]*psds[m1-1][2*6+1];
		dsfs[m][4*10+2] = QC[2]*dsds[m][4*6+2] + WQ[2]*dsds[m1][4*6+2]
				+ 2*tmp.dsps[4*3+2] + rze2[1]*psds[m1-1][1*6+2];
		dsfs[m][4*10+3] = QC[0]*dsds[m][4*6+3] + WQ[0]*dsds[m1][4*6+3]
				+   tmp.dsps[4*3+1];
		dsfs[m][4*10+4] = QC[1]*dsds[m][4*6+4] + WQ[1]*dsds[m1][4*6+4]
				+   tmp.dsps[4*3+2] + rze2[1]*psds[m1-1][2*6+4];
		dsfs[m][4*10+5] = QC[2]*dsds[m][4*6+5] + WQ[2]*dsds[m1][4*6+5]
				+   tmp.dsps[4*3+0] + rze2[1]*psds[m1-1][1*6+5];
		dsfs[m][4*10+6] = QC[0]*dsds[m][4*6+1] + WQ[0]*dsds[m1][4*6+1];
		dsfs[m][4*10+7] = QC[1]*dsds[m][4*6+2] + WQ[1]*dsds[m1][4*6+2]
				+ rze2[1]*psds[m1-1][2*6+2];
		dsfs[m][4*10+8] = QC[2]*dsds[m][4*6+0] + WQ[2]*dsds[m1][4*6+0]
				+ rze2[1]*psds[m1-1][1*6+0];
		dsfs[m][4*10+9] = QC[0]*dsds[m][4*6+4] + WQ[0]*dsds[m1][4*6+4];
		dsfs[m][5*10+0] = QC[0]*dsds[m][5*6+0] + WQ[0]*dsds[m1][5*6+0]
				+ 2*tmp.dsps[5*3+0] + rze2[1]*psds[m1-1][2*6+0];
		dsfs[m][5*10+1] = QC[1]*dsds[m][5*6+1] + WQ[1]*dsds[m1][5*6+1]
				+ 2*tmp.dsps[5*3+1];
		dsfs[m][5*10+2] = QC[2]*dsds[m][5*6+2] + WQ[2]*dsds[m1][5*6+2]
				+ 2*tmp.dsps[5*3+2] + rze2[1]*psds[m1-1][0*6+2];
		dsfs[m][5*10+3] = QC[0]*dsds[m][5*6+3] + WQ[0]*dsds[m1][5*6+3]
				+   tmp.dsps[5*3+1] + rze2[1]*psds[m1-1][2*6+3];
		dsfs[m][5*10+4] = QC[1]*dsds[m][5*6+4] + WQ[1]*dsds[m1][5*6+4]
				+   tmp.dsps[5*3+2];
		dsfs[m][5*10+5] = QC[2]*dsds[m][5*6+5] + WQ[2]*dsds[m1][5*6+5]
				+   tmp.dsps[5*3+0] + rze2[1]*psds[m1-1][0*6+5];
		dsfs[m][5*10+6] = QC[0]*dsds[m][5*6+1] + WQ[0]*dsds[m1][5*6+1]
				+ rze2[1]*psds[m1-1][2*6+1];
		dsfs[m][5*10+7] = QC[1]*dsds[m][5*6+2] + WQ[1]*dsds[m1][5*6+2];
		dsfs[m][5*10+8] = QC[2]*dsds[m][5*6+0] + WQ[2]*dsds[m1][5*6+0]
				+ rze2[1]*psds[m1-1][0*6+0];
		dsfs[m][5*10+9] = QC[0]*dsds[m][5*6+4] + WQ[0]*dsds[m1][5*6+4]
				+ rze2[1]*psds[m1-1][2*6+4];
	    }	// end m loop
	    // (F,S|F,S)
	    for (int m=0; m<=1; m++) {
		int m1 = m+1;
		for (int i=0; i<10*3; i++)
		    tmp.fsps[i] = eta2 * ( fsps[m][i] - re*fsps[m1][i] );

		fsfs[m][0*10+0] = QC[0]*fsds[m][0*6+0] + WQ[0]*fsds[m1][0*6+0]
				+ 2*tmp.fsps[0*3+0] + rze2[3]*dsds[m1][0*6+0];
		fsfs[m][0*10+1] = QC[1]*fsds[m][0*6+1] + WQ[1]*fsds[m1][0*6+1]
				+ 2*tmp.fsps[0*3+1];
		fsfs[m][0*10+2] = QC[2]*fsds[m][0*6+2] + WQ[2]*fsds[m1][0*6+2]
				+ 2*tmp.fsps[0*3+2];
		fsfs[m][0*10+3] = QC[0]*fsds[m][0*6+3] + WQ[0]*fsds[m1][0*6+3]
				+   tmp.fsps[0*3+1] + rze2[3]*dsds[m1][0*6+3];
		fsfs[m][0*10+4] = QC[1]*fsds[m][0*6+4] + WQ[1]*fsds[m1][0*6+4]
				+   tmp.fsps[0*3+2];
		fsfs[m][0*10+5] = QC[2]*fsds[m][0*6+5] + WQ[2]*fsds[m1][0*6+5]
				+   tmp.fsps[0*3+0];
		fsfs[m][0*10+6] = QC[0]*fsds[m][0*6+1] + WQ[0]*fsds[m1][0*6+1]
				+ rze2[3]*dsds[m1][0*6+1];
		fsfs[m][0*10+7] = QC[1]*fsds[m][0*6+2] + WQ[1]*fsds[m1][0*6+2];
		fsfs[m][0*10+8] = QC[2]*fsds[m][0*6+0] + WQ[2]*fsds[m1][0*6+0];
		fsfs[m][0*10+9] = QC[0]*fsds[m][0*6+4] + WQ[0]*fsds[m1][0*6+4]
				+ rze2[3]*dsds[m1][0*6+4];
		fsfs[m][1*10+0] = QC[0]*fsds[m][1*6+0] + WQ[0]*fsds[m1][1*6+0]
				+ 2*tmp.fsps[1*3+0];
		fsfs[m][1*10+1] = QC[1]*fsds[m][1*6+1] + WQ[1]*fsds[m1][1*6+1]
				+ 2*tmp.fsps[1*3+1] + rze2[3]*dsds[m1][1*6+1];
		fsfs[m][1*10+2] = QC[2]*fsds[m][1*6+2] + WQ[2]*fsds[m1][1*6+2]
				+ 2*tmp.fsps[1*3+2];
		fsfs[m][1*10+3] = QC[0]*fsds[m][1*6+3] + WQ[0]*fsds[m1][1*6+3]
				+   tmp.fsps[1*3+1];
		fsfs[m][1*10+4] = QC[1]*fsds[m][1*6+4] + WQ[1]*fsds[m1][1*6+4]
				+   tmp.fsps[1*3+2] + rze2[3]*dsds[m1][1*6+4];
		fsfs[m][1*10+5] = QC[2]*fsds[m][1*6+5] + WQ[2]*fsds[m1][1*6+5]
				+   tmp.fsps[1*3+0];
		fsfs[m][1*10+6] = QC[0]*fsds[m][1*6+1] + WQ[0]*fsds[m1][1*6+1];
		fsfs[m][1*10+7] = QC[1]*fsds[m][1*6+2] + WQ[1]*fsds[m1][1*6+2]
				+ rze2[3]*dsds[m1][1*6+2];
		fsfs[m][1*10+8] = QC[2]*fsds[m][1*6+0] + WQ[2]*fsds[m1][1*6+0];
		fsfs[m][1*10+9] = QC[0]*fsds[m][1*6+4] + WQ[0]*fsds[m1][1*6+4];
		fsfs[m][2*10+0] = QC[0]*fsds[m][2*6+0] + WQ[0]*fsds[m1][2*6+0]
				+ 2*tmp.fsps[2*3+0];
		fsfs[m][2*10+1] = QC[1]*fsds[m][2*6+1] + WQ[1]*fsds[m1][2*6+1]
				+ 2*tmp.fsps[2*3+1];
		fsfs[m][2*10+2] = QC[2]*fsds[m][2*6+2] + WQ[2]*fsds[m1][2*6+2]
				+ 2*tmp.fsps[2*3+2] + rze2[3]*dsds[m1][2*6+2];
		fsfs[m][2*10+3] = QC[0]*fsds[m][2*6+3] + WQ[0]*fsds[m1][2*6+3]
				+   tmp.fsps[2*3+1];
		fsfs[m][2*10+4] = QC[1]*fsds[m][2*6+4] + WQ[1]*fsds[m1][2*6+4]
				+   tmp.fsps[2*3+2];
		fsfs[m][2*10+5] = QC[2]*fsds[m][2*6+5] + WQ[2]*fsds[m1][2*6+5]
				+   tmp.fsps[2*3+0] + rze2[3]*dsds[m1][2*6+5];
		fsfs[m][2*10+6] = QC[0]*fsds[m][2*6+1] + WQ[0]*fsds[m1][2*6+1];
		fsfs[m][2*10+7] = QC[1]*fsds[m][2*6+2] + WQ[1]*fsds[m1][2*6+2];
		fsfs[m][2*10+8] = QC[2]*fsds[m][2*6+0] + WQ[2]*fsds[m1][2*6+0]
				+ rze2[3]*dsds[m1][2*6+0];
		fsfs[m][2*10+9] = QC[0]*fsds[m][2*6+4] + WQ[0]*fsds[m1][2*6+4];
		fsfs[m][3*10+0] = QC[0]*fsds[m][3*6+0] + WQ[0]*fsds[m1][3*6+0]
				+ 2*tmp.fsps[3*3+0] + rze2[2]*dsds[m1][3*6+0];
		fsfs[m][3*10+1] = QC[1]*fsds[m][3*6+1] + WQ[1]*fsds[m1][3*6+1]
				+ 2*tmp.fsps[3*3+1] + rze2[1]*dsds[m1][0*6+1];
		fsfs[m][3*10+2] = QC[2]*fsds[m][3*6+2] + WQ[2]*fsds[m1][3*6+2]
				+ 2*tmp.fsps[3*3+2];
		fsfs[m][3*10+3] = QC[0]*fsds[m][3*6+3] + WQ[0]*fsds[m1][3*6+3]
				+   tmp.fsps[3*3+1] + rze2[2]*dsds[m1][3*6+3];
		fsfs[m][3*10+4] = QC[1]*fsds[m][3*6+4] + WQ[1]*fsds[m1][3*6+4]
				+   tmp.fsps[3*3+2] + rze2[1]*dsds[m1][0*6+4];
		fsfs[m][3*10+5] = QC[2]*fsds[m][3*6+5] + WQ[2]*fsds[m1][3*6+5]
				+   tmp.fsps[3*3+0];
		fsfs[m][3*10+6] = QC[0]*fsds[m][3*6+1] + WQ[0]*fsds[m1][3*6+1]
				+ rze2[2]*dsds[m1][3*6+1];
		fsfs[m][3*10+7] = QC[1]*fsds[m][3*6+2] + WQ[1]*fsds[m1][3*6+2]
				+ rze2[1]*dsds[m1][0*6+2];
		fsfs[m][3*10+8] = QC[2]*fsds[m][3*6+0] + WQ[2]*fsds[m1][3*6+0];
		fsfs[m][3*10+9] = QC[0]*fsds[m][3*6+4] + WQ[0]*fsds[m1][3*6+4]
				+ rze2[2]*dsds[m1][3*6+4];
		fsfs[m][4*10+0] = QC[0]*fsds[m][4*6+0] + WQ[0]*fsds[m1][4*6+0]
				+ 2*tmp.fsps[4*3+0];
		fsfs[m][4*10+1] = QC[1]*fsds[m][4*6+1] + WQ[1]*fsds[m1][4*6+1]
				+ 2*tmp.fsps[4*3+1] + rze2[2]*dsds[m1][4*6+1];
		fsfs[m][4*10+2] = QC[2]*fsds[m][4*6+2] + WQ[2]*fsds[m1][4*6+2]
				+ 2*tmp.fsps[4*3+2] + rze2[1]*dsds[m1][1*6+2];
		fsfs[m][4*10+3] = QC[0]*fsds[m][4*6+3] + WQ[0]*fsds[m1][4*6+3]
				+   tmp.fsps[4*3+1];
		fsfs[m][4*10+4] = QC[1]*fsds[m][4*6+4] + WQ[1]*fsds[m1][4*6+4]
				+   tmp.fsps[4*3+2] + rze2[2]*dsds[m1][4*6+4];
		fsfs[m][4*10+5] = QC[2]*fsds[m][4*6+5] + WQ[2]*fsds[m1][4*6+5]
				+   tmp.fsps[4*3+0] + rze2[1]*dsds[m1][1*6+5];
		fsfs[m][4*10+6] = QC[0]*fsds[m][4*6+1] + WQ[0]*fsds[m1][4*6+1];
		fsfs[m][4*10+7] = QC[1]*fsds[m][4*6+2] + WQ[1]*fsds[m1][4*6+2]
				+ rze2[2]*dsds[m1][4*6+2];
		fsfs[m][4*10+8] = QC[2]*fsds[m][4*6+0] + WQ[2]*fsds[m1][4*6+0]
				+ rze2[1]*dsds[m1][1*6+0];
		fsfs[m][4*10+9] = QC[0]*fsds[m][4*6+4] + WQ[0]*fsds[m1][4*6+4];
		fsfs[m][5*10+0] = QC[0]*fsds[m][5*6+0] + WQ[0]*fsds[m1][5*6+0]
				+ 2*tmp.fsps[5*3+0] + rze2[1]*dsds[m1][2*6+0];
		fsfs[m][5*10+1] = QC[1]*fsds[m][5*6+1] + WQ[1]*fsds[m1][5*6+1]
				+ 2*tmp.fsps[5*3+1];
		fsfs[m][5*10+2] = QC[2]*fsds[m][5*6+2] + WQ[2]*fsds[m1][5*6+2]
				+ 2*tmp.fsps[5*3+2] + rze2[2]*dsds[m1][5*6+2];
		fsfs[m][5*10+3] = QC[0]*fsds[m][5*6+3] + WQ[0]*fsds[m1][5*6+3]
				+   tmp.fsps[5*3+1] + rze2[1]*dsds[m1][2*6+3];
		fsfs[m][5*10+4] = QC[1]*fsds[m][5*6+4] + WQ[1]*fsds[m1][5*6+4]
				+   tmp.fsps[5*3+2];
		fsfs[m][5*10+5] = QC[2]*fsds[m][5*6+5] + WQ[2]*fsds[m1][5*6+5]
				+   tmp.fsps[5*3+0] + rze2[2]*dsds[m1][5*6+5];
		fsfs[m][5*10+6] = QC[0]*fsds[m][5*6+1] + WQ[0]*fsds[m1][5*6+1]
				+ rze2[1]*dsds[m1][2*6+1];
		fsfs[m][5*10+7] = QC[1]*fsds[m][5*6+2] + WQ[1]*fsds[m1][5*6+2];
		fsfs[m][5*10+8] = QC[2]*fsds[m][5*6+0] + WQ[2]*fsds[m1][5*6+0]
				+ rze2[2]*dsds[m1][5*6+0];
		fsfs[m][5*10+9] = QC[0]*fsds[m][5*6+4] + WQ[0]*fsds[m1][5*6+4]
				+ rze2[1]*dsds[m1][2*6+4];
		fsfs[m][6*10+0] = QC[0]*fsds[m][6*6+0] + WQ[0]*fsds[m1][6*6+0]
				+ 2*tmp.fsps[6*3+0] + rze2[1]*dsds[m1][1*6+0];
		fsfs[m][6*10+1] = QC[1]*fsds[m][6*6+1] + WQ[1]*fsds[m1][6*6+1]
				+ 2*tmp.fsps[6*3+1] + rze2[2]*dsds[m1][3*6+1];
		fsfs[m][6*10+2] = QC[2]*fsds[m][6*6+2] + WQ[2]*fsds[m1][6*6+2]
				+ 2*tmp.fsps[6*3+2];
		fsfs[m][6*10+3] = QC[0]*fsds[m][6*6+3] + WQ[0]*fsds[m1][6*6+3]
				+   tmp.fsps[6*3+1] + rze2[1]*dsds[m1][1*6+3];
		fsfs[m][6*10+4] = QC[1]*fsds[m][6*6+4] + WQ[1]*fsds[m1][6*6+4]
				+   tmp.fsps[6*3+2] + rze2[2]*dsds[m1][3*6+4];
		fsfs[m][6*10+5] = QC[2]*fsds[m][6*6+5] + WQ[2]*fsds[m1][6*6+5]
				+   tmp.fsps[6*3+0];
		fsfs[m][6*10+6] = QC[0]*fsds[m][6*6+1] + WQ[0]*fsds[m1][6*6+1]
				+ rze2[1]*dsds[m1][1*6+1];
		fsfs[m][6*10+7] = QC[1]*fsds[m][6*6+2] + WQ[1]*fsds[m1][6*6+2]
				+ rze2[2]*dsds[m1][3*6+2];
		fsfs[m][6*10+8] = QC[2]*fsds[m][6*6+0] + WQ[2]*fsds[m1][6*6+0];
		fsfs[m][6*10+9] = QC[0]*fsds[m][6*6+4] + WQ[0]*fsds[m1][6*6+4]
				+ rze2[1]*dsds[m1][1*6+4];
		fsfs[m][7*10+0] = QC[0]*fsds[m][7*6+0] + WQ[0]*fsds[m1][7*6+0]
				+ 2*tmp.fsps[7*3+0];
		fsfs[m][7*10+1] = QC[1]*fsds[m][7*6+1] + WQ[1]*fsds[m1][7*6+1]
				+ 2*tmp.fsps[7*3+1] + rze2[1]*dsds[m1][2*6+1];
		fsfs[m][7*10+2] = QC[2]*fsds[m][7*6+2] + WQ[2]*fsds[m1][7*6+2]
				+ 2*tmp.fsps[7*3+2] + rze2[2]*dsds[m1][4*6+2];
		fsfs[m][7*10+3] = QC[0]*fsds[m][7*6+3] + WQ[0]*fsds[m1][7*6+3]
				+   tmp.fsps[7*3+1];
		fsfs[m][7*10+4] = QC[1]*fsds[m][7*6+4] + WQ[1]*fsds[m1][7*6+4]
				+   tmp.fsps[7*3+2] + rze2[1]*dsds[m1][2*6+4];
		fsfs[m][7*10+5] = QC[2]*fsds[m][7*6+5] + WQ[2]*fsds[m1][7*6+5]
				+   tmp.fsps[7*3+0] + rze2[2]*dsds[m1][4*6+5];
		fsfs[m][7*10+6] = QC[0]*fsds[m][7*6+1] + WQ[0]*fsds[m1][7*6+1];
		fsfs[m][7*10+7] = QC[1]*fsds[m][7*6+2] + WQ[1]*fsds[m1][7*6+2]
				+ rze2[1]*dsds[m1][2*6+2];
		fsfs[m][7*10+8] = QC[2]*fsds[m][7*6+0] + WQ[2]*fsds[m1][7*6+0]
				+ rze2[2]*dsds[m1][4*6+0];
		fsfs[m][7*10+9] = QC[0]*fsds[m][7*6+4] + WQ[0]*fsds[m1][7*6+4];
		fsfs[m][8*10+0] = QC[0]*fsds[m][8*6+0] + WQ[0]*fsds[m1][8*6+0]
				+ 2*tmp.fsps[8*3+0] + rze2[2]*dsds[m1][5*6+0];
		fsfs[m][8*10+1] = QC[1]*fsds[m][8*6+1] + WQ[1]*fsds[m1][8*6+1]
				+ 2*tmp.fsps[8*3+1];
		fsfs[m][8*10+2] = QC[2]*fsds[m][8*6+2] + WQ[2]*fsds[m1][8*6+2]
				+ 2*tmp.fsps[8*3+2] + rze2[1]*dsds[m1][0*6+2];
		fsfs[m][8*10+3] = QC[0]*fsds[m][8*6+3] + WQ[0]*fsds[m1][8*6+3]
				+   tmp.fsps[8*3+1] + rze2[2]*dsds[m1][5*6+3];
		fsfs[m][8*10+4] = QC[1]*fsds[m][8*6+4] + WQ[1]*fsds[m1][8*6+4]
				+   tmp.fsps[8*3+2];
		fsfs[m][8*10+5] = QC[2]*fsds[m][8*6+5] + WQ[2]*fsds[m1][8*6+5]
				+   tmp.fsps[8*3+0] + rze2[1]*dsds[m1][0*6+5];
		fsfs[m][8*10+6] = QC[0]*fsds[m][8*6+1] + WQ[0]*fsds[m1][8*6+1]
				+ rze2[2]*dsds[m1][5*6+1];
		fsfs[m][8*10+7] = QC[1]*fsds[m][8*6+2] + WQ[1]*fsds[m1][8*6+2];
		fsfs[m][8*10+8] = QC[2]*fsds[m][8*6+0] + WQ[2]*fsds[m1][8*6+0]
				+ rze2[1]*dsds[m1][0*6+0];
		fsfs[m][8*10+9] = QC[0]*fsds[m][8*6+4] + WQ[0]*fsds[m1][8*6+4]
				+ rze2[2]*dsds[m1][5*6+4];
		fsfs[m][9*10+0] = QC[0]*fsds[m][9*6+0] + WQ[0]*fsds[m1][9*6+0]
				+ 2*tmp.fsps[9*3+0] + rze2[1]*dsds[m1][4*6+0];
		fsfs[m][9*10+1] = QC[1]*fsds[m][9*6+1] + WQ[1]*fsds[m1][9*6+1]
				+ 2*tmp.fsps[9*3+1] + rze2[1]*dsds[m1][5*6+1];
		fsfs[m][9*10+2] = QC[2]*fsds[m][9*6+2] + WQ[2]*fsds[m1][9*6+2]
				+ 2*tmp.fsps[9*3+2] + rze2[1]*dsds[m1][3*6+2];
		fsfs[m][9*10+3] = QC[0]*fsds[m][9*6+3] + WQ[0]*fsds[m1][9*6+3]
				+   tmp.fsps[9*3+1] + rze2[1]*dsds[m1][4*6+3];
		fsfs[m][9*10+4] = QC[1]*fsds[m][9*6+4] + WQ[1]*fsds[m1][9*6+4]
				+   tmp.fsps[9*3+2] + rze2[1]*dsds[m1][5*6+4];
		fsfs[m][9*10+5] = QC[2]*fsds[m][9*6+5] + WQ[2]*fsds[m1][9*6+5]
				+   tmp.fsps[9*3+0] + rze2[1]*dsds[m1][3*6+5];
		fsfs[m][9*10+6] = QC[0]*fsds[m][9*6+1] + WQ[0]*fsds[m1][9*6+1]
				+ rze2[1]*dsds[m1][4*6+1];
		fsfs[m][9*10+7] = QC[1]*fsds[m][9*6+2] + WQ[1]*fsds[m1][9*6+2]
				+ rze2[1]*dsds[m1][5*6+2];
		fsfs[m][9*10+8] = QC[2]*fsds[m][9*6+0] + WQ[2]*fsds[m1][9*6+0]
				+ rze2[1]*dsds[m1][3*6+0];
		fsfs[m][9*10+9] = QC[0]*fsds[m][9*6+4] + WQ[0]*fsds[m1][9*6+4]
				+ rze2[1]*dsds[m1][4*6+4];
	    }	// end m loop
	    // (G,S|F,S)
	    for (int m=0; m<=1; m++) {
		int m1 = m+1;
		for (int i=0; i<15*3; i++)
		    tmp.gsps[i] = eta2 * ( gsps[m][i] - re*gsps[m1][i] );

		gsfs[m][ 0*10+0] = QC[0]*gsds[m][ 0*6+0] + WQ[0]*gsds[m1][ 0*6+0]
				 + 2*tmp.gsps[ 0*3+0] + rze2[4]*fsds[m1][0*6+0];
		gsfs[m][ 0*10+1] = QC[1]*gsds[m][ 0*6+1] + WQ[1]*gsds[m1][ 0*6+1]
				 + 2*tmp.gsps[ 0*3+1];
		gsfs[m][ 0*10+2] = QC[2]*gsds[m][ 0*6+2] + WQ[2]*gsds[m1][ 0*6+2]
				 + 2*tmp.gsps[ 0*3+2];
		gsfs[m][ 0*10+3] = QC[0]*gsds[m][ 0*6+3] + WQ[0]*gsds[m1][ 0*6+3]
				 +   tmp.gsps[ 0*3+1] + rze2[4]*fsds[m1][0*6+3];
		gsfs[m][ 0*10+4] = QC[1]*gsds[m][ 0*6+4] + WQ[1]*gsds[m1][ 0*6+4]
				 +   tmp.gsps[ 0*3+2];
		gsfs[m][ 0*10+5] = QC[2]*gsds[m][ 0*6+5] + WQ[2]*gsds[m1][ 0*6+5]
				 +   tmp.gsps[ 0*3+0];
		gsfs[m][ 0*10+6] = QC[0]*gsds[m][ 0*6+1] + WQ[0]*gsds[m1][ 0*6+1]
				 + rze2[4]*fsds[m1][0*6+1];
		gsfs[m][ 0*10+7] = QC[1]*gsds[m][ 0*6+2] + WQ[1]*gsds[m1][ 0*6+2];
		gsfs[m][ 0*10+8] = QC[2]*gsds[m][ 0*6+0] + WQ[2]*gsds[m1][ 0*6+0];
		gsfs[m][ 0*10+9] = QC[0]*gsds[m][ 0*6+4] + WQ[0]*gsds[m1][ 0*6+4]
				 + rze2[4]*fsds[m1][0*6+4];
		gsfs[m][ 1*10+0] = QC[0]*gsds[m][ 1*6+0] + WQ[0]*gsds[m1][ 1*6+0]
				 + 2*tmp.gsps[ 1*3+0];
		gsfs[m][ 1*10+1] = QC[1]*gsds[m][ 1*6+1] + WQ[1]*gsds[m1][ 1*6+1]
				 + 2*tmp.gsps[ 1*3+1] + rze2[4]*fsds[m1][1*6+1];
		gsfs[m][ 1*10+2] = QC[2]*gsds[m][ 1*6+2] + WQ[2]*gsds[m1][ 1*6+2]
				 + 2*tmp.gsps[ 1*3+2];
		gsfs[m][ 1*10+3] = QC[0]*gsds[m][ 1*6+3] + WQ[0]*gsds[m1][ 1*6+3]
				 +   tmp.gsps[ 1*3+1];
		gsfs[m][ 1*10+4] = QC[1]*gsds[m][ 1*6+4] + WQ[1]*gsds[m1][ 1*6+4]
				 +   tmp.gsps[ 1*3+2] + rze2[4]*fsds[m1][1*6+4];
		gsfs[m][ 1*10+5] = QC[2]*gsds[m][ 1*6+5] + WQ[2]*gsds[m1][ 1*6+5]
				 +   tmp.gsps[ 1*3+0];
		gsfs[m][ 1*10+6] = QC[0]*gsds[m][ 1*6+1] + WQ[0]*gsds[m1][ 1*6+1];
		gsfs[m][ 1*10+7] = QC[1]*gsds[m][ 1*6+2] + WQ[1]*gsds[m1][ 1*6+2]
				 + rze2[4]*fsds[m1][1*6+2];
		gsfs[m][ 1*10+8] = QC[2]*gsds[m][ 1*6+0] + WQ[2]*gsds[m1][ 1*6+0];
		gsfs[m][ 1*10+9] = QC[0]*gsds[m][ 1*6+4] + WQ[0]*gsds[m1][ 1*6+4];
		gsfs[m][ 2*10+0] = QC[0]*gsds[m][ 2*6+0] + WQ[0]*gsds[m1][ 2*6+0]
				 + 2*tmp.gsps[ 2*3+0];
		gsfs[m][ 2*10+1] = QC[1]*gsds[m][ 2*6+1] + WQ[1]*gsds[m1][ 2*6+1]
				 + 2*tmp.gsps[ 2*3+1];
		gsfs[m][ 2*10+2] = QC[2]*gsds[m][ 2*6+2] + WQ[2]*gsds[m1][ 2*6+2]
				 + 2*tmp.gsps[ 2*3+2] + rze2[4]*fsds[m1][2*6+2];
		gsfs[m][ 2*10+3] = QC[0]*gsds[m][ 2*6+3] + WQ[0]*gsds[m1][ 2*6+3]
				 +   tmp.gsps[ 2*3+1];
		gsfs[m][ 2*10+4] = QC[1]*gsds[m][ 2*6+4] + WQ[1]*gsds[m1][ 2*6+4]
				 +   tmp.gsps[ 2*3+2];
		gsfs[m][ 2*10+5] = QC[2]*gsds[m][ 2*6+5] + WQ[2]*gsds[m1][ 2*6+5]
				 +   tmp.gsps[ 2*3+0] + rze2[4]*fsds[m1][2*6+5];
		gsfs[m][ 2*10+6] = QC[0]*gsds[m][ 2*6+1] + WQ[0]*gsds[m1][ 2*6+1];
		gsfs[m][ 2*10+7] = QC[1]*gsds[m][ 2*6+2] + WQ[1]*gsds[m1][ 2*6+2];
		gsfs[m][ 2*10+8] = QC[2]*gsds[m][ 2*6+0] + WQ[2]*gsds[m1][ 2*6+0]
				 + rze2[4]*fsds[m1][2*6+0];
		gsfs[m][ 2*10+9] = QC[0]*gsds[m][ 2*6+4] + WQ[0]*gsds[m1][ 2*6+4];
		gsfs[m][ 3*10+0] = QC[0]*gsds[m][ 3*6+0] + WQ[0]*gsds[m1][ 3*6+0]
				 + 2*tmp.gsps[ 3*3+0] + rze2[3]*fsds[m1][3*6+0];
		gsfs[m][ 3*10+1] = QC[1]*gsds[m][ 3*6+1] + WQ[1]*gsds[m1][ 3*6+1]
				 + 2*tmp.gsps[ 3*3+1] + rze2[1]*fsds[m1][0*6+1];
		gsfs[m][ 3*10+2] = QC[2]*gsds[m][ 3*6+2] + WQ[2]*gsds[m1][ 3*6+2]
				 + 2*tmp.gsps[ 3*3+2];
		gsfs[m][ 3*10+3] = QC[0]*gsds[m][ 3*6+3] + WQ[0]*gsds[m1][ 3*6+3]
				 +   tmp.gsps[ 3*3+1] + rze2[3]*fsds[m1][3*6+3];
		gsfs[m][ 3*10+4] = QC[1]*gsds[m][ 3*6+4] + WQ[1]*gsds[m1][ 3*6+4]
				 +   tmp.gsps[ 3*3+2] + rze2[1]*fsds[m1][0*6+4];
		gsfs[m][ 3*10+5] = QC[2]*gsds[m][ 3*6+5] + WQ[2]*gsds[m1][ 3*6+5]
				 +   tmp.gsps[ 3*3+0];
		gsfs[m][ 3*10+6] = QC[0]*gsds[m][ 3*6+1] + WQ[0]*gsds[m1][ 3*6+1]
				 + rze2[3]*fsds[m1][3*6+1];
		gsfs[m][ 3*10+7] = QC[1]*gsds[m][ 3*6+2] + WQ[1]*gsds[m1][ 3*6+2]
				 + rze2[1]*fsds[m1][0*6+2];
		gsfs[m][ 3*10+8] = QC[2]*gsds[m][ 3*6+0] + WQ[2]*gsds[m1][ 3*6+0];
		gsfs[m][ 3*10+9] = QC[0]*gsds[m][ 3*6+4] + WQ[0]*gsds[m1][ 3*6+4]
				 + rze2[3]*fsds[m1][3*6+4];
		gsfs[m][ 4*10+0] = QC[0]*gsds[m][ 4*6+0] + WQ[0]*gsds[m1][ 4*6+0]
				 + 2*tmp.gsps[ 4*3+0];
		gsfs[m][ 4*10+1] = QC[1]*gsds[m][ 4*6+1] + WQ[1]*gsds[m1][ 4*6+1]
				 + 2*tmp.gsps[ 4*3+1] + rze2[3]*fsds[m1][4*6+1];
		gsfs[m][ 4*10+2] = QC[2]*gsds[m][ 4*6+2] + WQ[2]*gsds[m1][ 4*6+2]
				 + 2*tmp.gsps[ 4*3+2] + rze2[1]*fsds[m1][1*6+2];
		gsfs[m][ 4*10+3] = QC[0]*gsds[m][ 4*6+3] + WQ[0]*gsds[m1][ 4*6+3]
				 +   tmp.gsps[ 4*3+1];
		gsfs[m][ 4*10+4] = QC[1]*gsds[m][ 4*6+4] + WQ[1]*gsds[m1][ 4*6+4]
				 +   tmp.gsps[ 4*3+2] + rze2[3]*fsds[m1][4*6+4];
		gsfs[m][ 4*10+5] = QC[2]*gsds[m][ 4*6+5] + WQ[2]*gsds[m1][ 4*6+5]
				 +   tmp.gsps[ 4*3+0] + rze2[1]*fsds[m1][1*6+5];
		gsfs[m][ 4*10+6] = QC[0]*gsds[m][ 4*6+1] + WQ[0]*gsds[m1][ 4*6+1];
		gsfs[m][ 4*10+7] = QC[1]*gsds[m][ 4*6+2] + WQ[1]*gsds[m1][ 4*6+2]
				 + rze2[3]*fsds[m1][4*6+2];
		gsfs[m][ 4*10+8] = QC[2]*gsds[m][ 4*6+0] + WQ[2]*gsds[m1][ 4*6+0]
				 + rze2[1]*fsds[m1][1*6+0];
		gsfs[m][ 4*10+9] = QC[0]*gsds[m][ 4*6+4] + WQ[0]*gsds[m1][ 4*6+4];
		gsfs[m][ 5*10+0] = QC[0]*gsds[m][ 5*6+0] + WQ[0]*gsds[m1][ 5*6+0]
				 + 2*tmp.gsps[ 5*3+0] + rze2[1]*fsds[m1][2*6+0];
		gsfs[m][ 5*10+1] = QC[1]*gsds[m][ 5*6+1] + WQ[1]*gsds[m1][ 5*6+1]
				 + 2*tmp.gsps[ 5*3+1];
		gsfs[m][ 5*10+2] = QC[2]*gsds[m][ 5*6+2] + WQ[2]*gsds[m1][ 5*6+2]
				 + 2*tmp.gsps[ 5*3+2] + rze2[3]*fsds[m1][5*6+2];
		gsfs[m][ 5*10+3] = QC[0]*gsds[m][ 5*6+3] + WQ[0]*gsds[m1][ 5*6+3]
				 +   tmp.gsps[ 5*3+1] + rze2[1]*fsds[m1][2*6+3];
		gsfs[m][ 5*10+4] = QC[1]*gsds[m][ 5*6+4] + WQ[1]*gsds[m1][ 5*6+4]
				 +   tmp.gsps[ 5*3+2];
		gsfs[m][ 5*10+5] = QC[2]*gsds[m][ 5*6+5] + WQ[2]*gsds[m1][ 5*6+5]
				 +   tmp.gsps[ 5*3+0] + rze2[3]*fsds[m1][5*6+5];
		gsfs[m][ 5*10+6] = QC[0]*gsds[m][ 5*6+1] + WQ[0]*gsds[m1][ 5*6+1]
				 + rze2[1]*fsds[m1][2*6+1];
		gsfs[m][ 5*10+7] = QC[1]*gsds[m][ 5*6+2] + WQ[1]*gsds[m1][ 5*6+2];
		gsfs[m][ 5*10+8] = QC[2]*gsds[m][ 5*6+0] + WQ[2]*gsds[m1][ 5*6+0]
				 + rze2[3]*fsds[m1][5*6+0];
		gsfs[m][ 5*10+9] = QC[0]*gsds[m][ 5*6+4] + WQ[0]*gsds[m1][ 5*6+4]
				 + rze2[1]*fsds[m1][2*6+4];
		gsfs[m][ 6*10+0] = QC[0]*gsds[m][ 6*6+0] + WQ[0]*gsds[m1][ 6*6+0]
				 + 2*tmp.gsps[ 6*3+0] + rze2[2]*fsds[m1][6*6+0];
		gsfs[m][ 6*10+1] = QC[1]*gsds[m][ 6*6+1] + WQ[1]*gsds[m1][ 6*6+1]
				 + 2*tmp.gsps[ 6*3+1] + rze2[2]*fsds[m1][3*6+1];
		gsfs[m][ 6*10+2] = QC[2]*gsds[m][ 6*6+2] + WQ[2]*gsds[m1][ 6*6+2]
				 + 2*tmp.gsps[ 6*3+2];
		gsfs[m][ 6*10+3] = QC[0]*gsds[m][ 6*6+3] + WQ[0]*gsds[m1][ 6*6+3]
				 +   tmp.gsps[ 6*3+1] + rze2[2]*fsds[m1][6*6+3];
		gsfs[m][ 6*10+4] = QC[1]*gsds[m][ 6*6+4] + WQ[1]*gsds[m1][ 6*6+4]
				 +   tmp.gsps[ 6*3+2] + rze2[2]*fsds[m1][3*6+4];
		gsfs[m][ 6*10+5] = QC[2]*gsds[m][ 6*6+5] + WQ[2]*gsds[m1][ 6*6+5]
				 +   tmp.gsps[ 6*3+0];
		gsfs[m][ 6*10+6] = QC[0]*gsds[m][ 6*6+1] + WQ[0]*gsds[m1][ 6*6+1]
				 + rze2[2]*fsds[m1][6*6+1];
		gsfs[m][ 6*10+7] = QC[1]*gsds[m][ 6*6+2] + WQ[1]*gsds[m1][ 6*6+2]
				 + rze2[2]*fsds[m1][3*6+2];
		gsfs[m][ 6*10+8] = QC[2]*gsds[m][ 6*6+0] + WQ[2]*gsds[m1][ 6*6+0];
		gsfs[m][ 6*10+9] = QC[0]*gsds[m][ 6*6+4] + WQ[0]*gsds[m1][ 6*6+4]
				 + rze2[2]*fsds[m1][6*6+4];
		gsfs[m][ 7*10+0] = QC[0]*gsds[m][ 7*6+0] + WQ[0]*gsds[m1][ 7*6+0]
				 + 2*tmp.gsps[ 7*3+0];
		gsfs[m][ 7*10+1] = QC[1]*gsds[m][ 7*6+1] + WQ[1]*gsds[m1][ 7*6+1]
				 + 2*tmp.gsps[ 7*3+1] + rze2[2]*fsds[m1][7*6+1];
		gsfs[m][ 7*10+2] = QC[2]*gsds[m][ 7*6+2] + WQ[2]*gsds[m1][ 7*6+2]
				 + 2*tmp.gsps[ 7*3+2] + rze2[2]*fsds[m1][4*6+2];
		gsfs[m][ 7*10+3] = QC[0]*gsds[m][ 7*6+3] + WQ[0]*gsds[m1][ 7*6+3]
				 +   tmp.gsps[ 7*3+1];
		gsfs[m][ 7*10+4] = QC[1]*gsds[m][ 7*6+4] + WQ[1]*gsds[m1][ 7*6+4]
				 +   tmp.gsps[ 7*3+2] + rze2[2]*fsds[m1][7*6+4];
		gsfs[m][ 7*10+5] = QC[2]*gsds[m][ 7*6+5] + WQ[2]*gsds[m1][ 7*6+5]
				 +   tmp.gsps[ 7*3+0] + rze2[2]*fsds[m1][4*6+5];
		gsfs[m][ 7*10+6] = QC[0]*gsds[m][ 7*6+1] + WQ[0]*gsds[m1][ 7*6+1];
		gsfs[m][ 7*10+7] = QC[1]*gsds[m][ 7*6+2] + WQ[1]*gsds[m1][ 7*6+2]
				 + rze2[2]*fsds[m1][7*6+2];
		gsfs[m][ 7*10+8] = QC[2]*gsds[m][ 7*6+0] + WQ[2]*gsds[m1][ 7*6+0]
				 + rze2[2]*fsds[m1][4*6+0];
		gsfs[m][ 7*10+9] = QC[0]*gsds[m][ 7*6+4] + WQ[0]*gsds[m1][ 7*6+4];
		gsfs[m][ 8*10+0] = QC[0]*gsds[m][ 8*6+0] + WQ[0]*gsds[m1][ 8*6+0]
				 + 2*tmp.gsps[ 8*3+0] + rze2[2]*fsds[m1][5*6+0];
		gsfs[m][ 8*10+1] = QC[1]*gsds[m][ 8*6+1] + WQ[1]*gsds[m1][ 8*6+1]
				 + 2*tmp.gsps[ 8*3+1];
		gsfs[m][ 8*10+2] = QC[2]*gsds[m][ 8*6+2] + WQ[2]*gsds[m1][ 8*6+2]
				 + 2*tmp.gsps[ 8*3+2] + rze2[2]*fsds[m1][8*6+2];
		gsfs[m][ 8*10+3] = QC[0]*gsds[m][ 8*6+3] + WQ[0]*gsds[m1][ 8*6+3]
				 +   tmp.gsps[ 8*3+1] + rze2[2]*fsds[m1][5*6+3];
		gsfs[m][ 8*10+4] = QC[1]*gsds[m][ 8*6+4] + WQ[1]*gsds[m1][ 8*6+4]
				 +   tmp.gsps[ 8*3+2];
		gsfs[m][ 8*10+5] = QC[2]*gsds[m][ 8*6+5] + WQ[2]*gsds[m1][ 8*6+5]
				 +   tmp.gsps[ 8*3+0] + rze2[2]*fsds[m1][8*6+5];
		gsfs[m][ 8*10+6] = QC[0]*gsds[m][ 8*6+1] + WQ[0]*gsds[m1][ 8*6+1]
				 + rze2[2]*fsds[m1][5*6+1];
		gsfs[m][ 8*10+7] = QC[1]*gsds[m][ 8*6+2] + WQ[1]*gsds[m1][ 8*6+2];
		gsfs[m][ 8*10+8] = QC[2]*gsds[m][ 8*6+0] + WQ[2]*gsds[m1][ 8*6+0]
				 + rze2[2]*fsds[m1][8*6+0];
		gsfs[m][ 8*10+9] = QC[0]*gsds[m][ 8*6+4] + WQ[0]*gsds[m1][ 8*6+4]
				 + rze2[2]*fsds[m1][5*6+4];
		gsfs[m][ 9*10+0] = QC[0]*gsds[m][ 9*6+0] + WQ[0]*gsds[m1][ 9*6+0]
				 + 2*tmp.gsps[ 9*3+0] + rze2[1]*fsds[m1][1*6+0];
		gsfs[m][ 9*10+1] = QC[1]*gsds[m][ 9*6+1] + WQ[1]*gsds[m1][ 9*6+1]
				 + 2*tmp.gsps[ 9*3+1] + rze2[3]*fsds[m1][6*6+1];
		gsfs[m][ 9*10+2] = QC[2]*gsds[m][ 9*6+2] + WQ[2]*gsds[m1][ 9*6+2]
				 + 2*tmp.gsps[ 9*3+2];
		gsfs[m][ 9*10+3] = QC[0]*gsds[m][ 9*6+3] + WQ[0]*gsds[m1][ 9*6+3]
				 +   tmp.gsps[ 9*3+1] + rze2[1]*fsds[m1][1*6+3];
		gsfs[m][ 9*10+4] = QC[1]*gsds[m][ 9*6+4] + WQ[1]*gsds[m1][ 9*6+4]
				 +   tmp.gsps[ 9*3+2] + rze2[3]*fsds[m1][6*6+4];
		gsfs[m][ 9*10+5] = QC[2]*gsds[m][ 9*6+5] + WQ[2]*gsds[m1][ 9*6+5]
				 +   tmp.gsps[ 9*3+0];
		gsfs[m][ 9*10+6] = QC[0]*gsds[m][ 9*6+1] + WQ[0]*gsds[m1][ 9*6+1]
				 + rze2[1]*fsds[m1][1*6+1];
		gsfs[m][ 9*10+7] = QC[1]*gsds[m][ 9*6+2] + WQ[1]*gsds[m1][ 9*6+2]
				 + rze2[3]*fsds[m1][6*6+2];
		gsfs[m][ 9*10+8] = QC[2]*gsds[m][ 9*6+0] + WQ[2]*gsds[m1][ 9*6+0];
		gsfs[m][ 9*10+9] = QC[0]*gsds[m][ 9*6+4] + WQ[0]*gsds[m1][ 9*6+4]
				 + rze2[1]*fsds[m1][1*6+4];
		gsfs[m][10*10+0] = QC[0]*gsds[m][10*6+0] + WQ[0]*gsds[m1][10*6+0]
				 + 2*tmp.gsps[10*3+0];
		gsfs[m][10*10+1] = QC[1]*gsds[m][10*6+1] + WQ[1]*gsds[m1][10*6+1]
				 + 2*tmp.gsps[10*3+1] + rze2[1]*fsds[m1][2*6+1];
		gsfs[m][10*10+2] = QC[2]*gsds[m][10*6+2] + WQ[2]*gsds[m1][10*6+2]
				 + 2*tmp.gsps[10*3+2] + rze2[3]*fsds[m1][7*6+2];
		gsfs[m][10*10+3] = QC[0]*gsds[m][10*6+3] + WQ[0]*gsds[m1][10*6+3]
				 +   tmp.gsps[10*3+1];
		gsfs[m][10*10+4] = QC[1]*gsds[m][10*6+4] + WQ[1]*gsds[m1][10*6+4]
				 +   tmp.gsps[10*3+2] + rze2[1]*fsds[m1][2*6+4];
		gsfs[m][10*10+5] = QC[2]*gsds[m][10*6+5] + WQ[2]*gsds[m1][10*6+5]
				 +   tmp.gsps[10*3+0] + rze2[3]*fsds[m1][7*6+5];
		gsfs[m][10*10+6] = QC[0]*gsds[m][10*6+1] + WQ[0]*gsds[m1][10*6+1];
		gsfs[m][10*10+7] = QC[1]*gsds[m][10*6+2] + WQ[1]*gsds[m1][10*6+2]
				 + rze2[1]*fsds[m1][2*6+2];
		gsfs[m][10*10+8] = QC[2]*gsds[m][10*6+0] + WQ[2]*gsds[m1][10*6+0]
				 + rze2[3]*fsds[m1][7*6+0];
		gsfs[m][10*10+9] = QC[0]*gsds[m][10*6+4] + WQ[0]*gsds[m1][10*6+4];
		gsfs[m][11*10+0] = QC[0]*gsds[m][11*6+0] + WQ[0]*gsds[m1][11*6+0]
				 + 2*tmp.gsps[11*3+0] + rze2[3]*fsds[m1][8*6+0];
		gsfs[m][11*10+1] = QC[1]*gsds[m][11*6+1] + WQ[1]*gsds[m1][11*6+1]
				 + 2*tmp.gsps[11*3+1];
		gsfs[m][11*10+2] = QC[2]*gsds[m][11*6+2] + WQ[2]*gsds[m1][11*6+2]
				 + 2*tmp.gsps[11*3+2] + rze2[1]*fsds[m1][0*6+2];
		gsfs[m][11*10+3] = QC[0]*gsds[m][11*6+3] + WQ[0]*gsds[m1][11*6+3]
				 +   tmp.gsps[11*3+1] + rze2[3]*fsds[m1][8*6+3];
		gsfs[m][11*10+4] = QC[1]*gsds[m][11*6+4] + WQ[1]*gsds[m1][11*6+4]
				 +   tmp.gsps[11*3+2];
		gsfs[m][11*10+5] = QC[2]*gsds[m][11*6+5] + WQ[2]*gsds[m1][11*6+5]
				 +   tmp.gsps[11*3+0] + rze2[1]*fsds[m1][0*6+5];
		gsfs[m][11*10+6] = QC[0]*gsds[m][11*6+1] + WQ[0]*gsds[m1][11*6+1]
				 + rze2[3]*fsds[m1][8*6+1];
		gsfs[m][11*10+7] = QC[1]*gsds[m][11*6+2] + WQ[1]*gsds[m1][11*6+2];
		gsfs[m][11*10+8] = QC[2]*gsds[m][11*6+0] + WQ[2]*gsds[m1][11*6+0]
				 + rze2[1]*fsds[m1][0*6+0];
		gsfs[m][11*10+9] = QC[0]*gsds[m][11*6+4] + WQ[0]*gsds[m1][11*6+4]
				 + rze2[3]*fsds[m1][8*6+4];
		gsfs[m][12*10+0] = QC[0]*gsds[m][12*6+0] + WQ[0]*gsds[m1][12*6+0]
				 + 2*tmp.gsps[12*3+0] + rze2[2]*fsds[m1][9*6+0];
		gsfs[m][12*10+1] = QC[1]*gsds[m][12*6+1] + WQ[1]*gsds[m1][12*6+1]
				 + 2*tmp.gsps[12*3+1] + rze2[1]*fsds[m1][8*6+1];
		gsfs[m][12*10+2] = QC[2]*gsds[m][12*6+2] + WQ[2]*gsds[m1][12*6+2]
				 + 2*tmp.gsps[12*3+2] + rze2[1]*fsds[m1][3*6+2];
		gsfs[m][12*10+3] = QC[0]*gsds[m][12*6+3] + WQ[0]*gsds[m1][12*6+3]
				 +   tmp.gsps[12*3+1] + rze2[2]*fsds[m1][9*6+3];
		gsfs[m][12*10+4] = QC[1]*gsds[m][12*6+4] + WQ[1]*gsds[m1][12*6+4]
				 +   tmp.gsps[12*3+2] + rze2[1]*fsds[m1][8*6+4];
		gsfs[m][12*10+5] = QC[2]*gsds[m][12*6+5] + WQ[2]*gsds[m1][12*6+5]
				 +   tmp.gsps[12*3+0] + rze2[1]*fsds[m1][3*6+5];
		gsfs[m][12*10+6] = QC[0]*gsds[m][12*6+1] + WQ[0]*gsds[m1][12*6+1]
				 + rze2[2]*fsds[m1][9*6+1];
		gsfs[m][12*10+7] = QC[1]*gsds[m][12*6+2] + WQ[1]*gsds[m1][12*6+2]
				 + rze2[1]*fsds[m1][8*6+2];
		gsfs[m][12*10+8] = QC[2]*gsds[m][12*6+0] + WQ[2]*gsds[m1][12*6+0]
				 + rze2[1]*fsds[m1][3*6+0];
		gsfs[m][12*10+9] = QC[0]*gsds[m][12*6+4] + WQ[0]*gsds[m1][12*6+4]
				 + rze2[2]*fsds[m1][9*6+4];
		gsfs[m][13*10+0] = QC[0]*gsds[m][13*6+0] + WQ[0]*gsds[m1][13*6+0]
				 + 2*tmp.gsps[13*3+0] + rze2[1]*fsds[m1][4*6+0];
		gsfs[m][13*10+1] = QC[1]*gsds[m][13*6+1] + WQ[1]*gsds[m1][13*6+1]
				 + 2*tmp.gsps[13*3+1] + rze2[2]*fsds[m1][9*6+1];
		gsfs[m][13*10+2] = QC[2]*gsds[m][13*6+2] + WQ[2]*gsds[m1][13*6+2]
				 + 2*tmp.gsps[13*3+2] + rze2[1]*fsds[m1][6*6+2];
		gsfs[m][13*10+3] = QC[0]*gsds[m][13*6+3] + WQ[0]*gsds[m1][13*6+3]
				 +   tmp.gsps[13*3+1] + rze2[1]*fsds[m1][4*6+3];
		gsfs[m][13*10+4] = QC[1]*gsds[m][13*6+4] + WQ[1]*gsds[m1][13*6+4]
				 +   tmp.gsps[13*3+2] + rze2[2]*fsds[m1][9*6+4];
		gsfs[m][13*10+5] = QC[2]*gsds[m][13*6+5] + WQ[2]*gsds[m1][13*6+5]
				 +   tmp.gsps[13*3+0] + rze2[1]*fsds[m1][6*6+5];
		gsfs[m][13*10+6] = QC[0]*gsds[m][13*6+1] + WQ[0]*gsds[m1][13*6+1]
				 + rze2[1]*fsds[m1][4*6+1];
		gsfs[m][13*10+7] = QC[1]*gsds[m][13*6+2] + WQ[1]*gsds[m1][13*6+2]
				 + rze2[2]*fsds[m1][9*6+2];
		gsfs[m][13*10+8] = QC[2]*gsds[m][13*6+0] + WQ[2]*gsds[m1][13*6+0]
				 + rze2[1]*fsds[m1][6*6+0];
		gsfs[m][13*10+9] = QC[0]*gsds[m][13*6+4] + WQ[0]*gsds[m1][13*6+4]
				 + rze2[1]*fsds[m1][4*6+4];
		gsfs[m][14*10+0] = QC[0]*gsds[m][14*6+0] + WQ[0]*gsds[m1][14*6+0]
				 + 2*tmp.gsps[14*3+0] + rze2[1]*fsds[m1][7*6+0];
		gsfs[m][14*10+1] = QC[1]*gsds[m][14*6+1] + WQ[1]*gsds[m1][14*6+1]
				 + 2*tmp.gsps[14*3+1] + rze2[1]*fsds[m1][5*6+1];
		gsfs[m][14*10+2] = QC[2]*gsds[m][14*6+2] + WQ[2]*gsds[m1][14*6+2]
				 + 2*tmp.gsps[14*3+2] + rze2[2]*fsds[m1][9*6+2];
		gsfs[m][14*10+3] = QC[0]*gsds[m][14*6+3] + WQ[0]*gsds[m1][14*6+3]
				 +   tmp.gsps[14*3+1] + rze2[1]*fsds[m1][7*6+3];
		gsfs[m][14*10+4] = QC[1]*gsds[m][14*6+4] + WQ[1]*gsds[m1][14*6+4]
				 +   tmp.gsps[14*3+2] + rze2[1]*fsds[m1][5*6+4];
		gsfs[m][14*10+5] = QC[2]*gsds[m][14*6+5] + WQ[2]*gsds[m1][14*6+5]
				 +   tmp.gsps[14*3+0] + rze2[2]*fsds[m1][9*6+5];
		gsfs[m][14*10+6] = QC[0]*gsds[m][14*6+1] + WQ[0]*gsds[m1][14*6+1]
				 + rze2[1]*fsds[m1][7*6+1];
		gsfs[m][14*10+7] = QC[1]*gsds[m][14*6+2] + WQ[1]*gsds[m1][14*6+2]
				 + rze2[1]*fsds[m1][5*6+2];
		gsfs[m][14*10+8] = QC[2]*gsds[m][14*6+0] + WQ[2]*gsds[m1][14*6+0]
				 + rze2[2]*fsds[m1][9*6+0];
		gsfs[m][14*10+9] = QC[0]*gsds[m][14*6+4] + WQ[0]*gsds[m1][14*6+4]
				 + rze2[1]*fsds[m1][7*6+4];
	    }	// end m loop
	    // (D,S|G,S)
	    for (int i=0; i<6*6; i++)
		tmp.dsds[i] = eta2 * ( dsds[0][i] - re*dsds[1][i] );

	    dsgs[0*15+ 0] = QC[0]*dsfs[0][0*10+0] + WQ[0]*dsfs[1][0*10+0]
			  + 3*tmp.dsds[0*6+0] + rze2[2]*psfs[0*10+0];
	    dsgs[0*15+ 1] = QC[1]*dsfs[0][0*10+1] + WQ[1]*dsfs[1][0*10+1]
			  + 3*tmp.dsds[0*6+1];
	    dsgs[0*15+ 2] = QC[2]*dsfs[0][0*10+2] + WQ[2]*dsfs[1][0*10+2]
			  + 3*tmp.dsds[0*6+2];
	    dsgs[0*15+ 3] = QC[0]*dsfs[0][0*10+3] + WQ[0]*dsfs[1][0*10+3]
			  + 2*tmp.dsds[0*6+3] + rze2[2]*psfs[0*10+3];
	    dsgs[0*15+ 4] = QC[1]*dsfs[0][0*10+4] + WQ[1]*dsfs[1][0*10+4]
			  + 2*tmp.dsds[0*6+4];
	    dsgs[0*15+ 5] = QC[2]*dsfs[0][0*10+5] + WQ[2]*dsfs[1][0*10+5]
			  + 2*tmp.dsds[0*6+5];
	    dsgs[0*15+ 6] = QC[0]*dsfs[0][0*10+6] + WQ[0]*dsfs[1][0*10+6]
			  +   tmp.dsds[0*6+1] + rze2[2]*psfs[0*10+6];
	    dsgs[0*15+ 7] = QC[1]*dsfs[0][0*10+7] + WQ[1]*dsfs[1][0*10+7]
			  +   tmp.dsds[0*6+2];
	    dsgs[0*15+ 8] = QC[2]*dsfs[0][0*10+8] + WQ[2]*dsfs[1][0*10+8]
			  +   tmp.dsds[0*6+0];
	    dsgs[0*15+ 9] = QC[0]*dsfs[0][0*10+1] + WQ[0]*dsfs[1][0*10+1]
			  + rze2[2]*psfs[0*10+1];
	    dsgs[0*15+10] = QC[1]*dsfs[0][0*10+2] + WQ[1]*dsfs[1][0*10+2];
	    dsgs[0*15+11] = QC[2]*dsfs[0][0*10+0] + WQ[2]*dsfs[1][0*10+0];
	    dsgs[0*15+12] = QC[0]*dsfs[0][0*10+9] + WQ[0]*dsfs[1][0*10+9]
			  +   tmp.dsds[0*6+4] + rze2[2]*psfs[0*10+9];
	    dsgs[0*15+13] = QC[1]*dsfs[0][0*10+9] + WQ[1]*dsfs[1][0*10+9]
			  +   tmp.dsds[0*6+5];
	    dsgs[0*15+14] = QC[2]*dsfs[0][0*10+9] + WQ[2]*dsfs[1][0*10+9]
			  +   tmp.dsds[0*6+3];
	    dsgs[1*15+ 0] = QC[0]*dsfs[0][1*10+0] + WQ[0]*dsfs[1][1*10+0]
			  + 3*tmp.dsds[1*6+0];
	    dsgs[1*15+ 1] = QC[1]*dsfs[0][1*10+1] + WQ[1]*dsfs[1][1*10+1]
			  + 3*tmp.dsds[1*6+1] + rze2[2]*psfs[1*10+1];
	    dsgs[1*15+ 2] = QC[2]*dsfs[0][1*10+2] + WQ[2]*dsfs[1][1*10+2]
			  + 3*tmp.dsds[1*6+2];
	    dsgs[1*15+ 3] = QC[0]*dsfs[0][1*10+3] + WQ[0]*dsfs[1][1*10+3]
			  + 2*tmp.dsds[1*6+3];
	    dsgs[1*15+ 4] = QC[1]*dsfs[0][1*10+4] + WQ[1]*dsfs[1][1*10+4]
			  + 2*tmp.dsds[1*6+4] + rze2[2]*psfs[1*10+4];
	    dsgs[1*15+ 5] = QC[2]*dsfs[0][1*10+5] + WQ[2]*dsfs[1][1*10+5]
			  + 2*tmp.dsds[1*6+5];
	    dsgs[1*15+ 6] = QC[0]*dsfs[0][1*10+6] + WQ[0]*dsfs[1][1*10+6]
			  +   tmp.dsds[1*6+1];
	    dsgs[1*15+ 7] = QC[1]*dsfs[0][1*10+7] + WQ[1]*dsfs[1][1*10+7]
			  +   tmp.dsds[1*6+2] + rze2[2]*psfs[1*10+7];
	    dsgs[1*15+ 8] = QC[2]*dsfs[0][1*10+8] + WQ[2]*dsfs[1][1*10+8]
			  +   tmp.dsds[1*6+0];
	    dsgs[1*15+ 9] = QC[0]*dsfs[0][1*10+1] + WQ[0]*dsfs[1][1*10+1];
	    dsgs[1*15+10] = QC[1]*dsfs[0][1*10+2] + WQ[1]*dsfs[1][1*10+2]
			  + rze2[2]*psfs[1*10+2];
	    dsgs[1*15+11] = QC[2]*dsfs[0][1*10+0] + WQ[2]*dsfs[1][1*10+0];
	    dsgs[1*15+12] = QC[0]*dsfs[0][1*10+9] + WQ[0]*dsfs[1][1*10+9]
			  +   tmp.dsds[1*6+4];
	    dsgs[1*15+13] = QC[1]*dsfs[0][1*10+9] + WQ[1]*dsfs[1][1*10+9]
			  +   tmp.dsds[1*6+5] + rze2[2]*psfs[1*10+9];
	    dsgs[1*15+14] = QC[2]*dsfs[0][1*10+9] + WQ[2]*dsfs[1][1*10+9]
			  +   tmp.dsds[1*6+3];
	    dsgs[2*15+ 0] = QC[0]*dsfs[0][2*10+0] + WQ[0]*dsfs[1][2*10+0]
			  + 3*tmp.dsds[2*6+0];
	    dsgs[2*15+ 1] = QC[1]*dsfs[0][2*10+1] + WQ[1]*dsfs[1][2*10+1]
			  + 3*tmp.dsds[2*6+1];
	    dsgs[2*15+ 2] = QC[2]*dsfs[0][2*10+2] + WQ[2]*dsfs[1][2*10+2]
			  + 3*tmp.dsds[2*6+2] + rze2[2]*psfs[2*10+2];
	    dsgs[2*15+ 3] = QC[0]*dsfs[0][2*10+3] + WQ[0]*dsfs[1][2*10+3]
			  + 2*tmp.dsds[2*6+3];
	    dsgs[2*15+ 4] = QC[1]*dsfs[0][2*10+4] + WQ[1]*dsfs[1][2*10+4]
			  + 2*tmp.dsds[2*6+4];
	    dsgs[2*15+ 5] = QC[2]*dsfs[0][2*10+5] + WQ[2]*dsfs[1][2*10+5]
			  + 2*tmp.dsds[2*6+5] + rze2[2]*psfs[2*10+5];
	    dsgs[2*15+ 6] = QC[0]*dsfs[0][2*10+6] + WQ[0]*dsfs[1][2*10+6]
			  +   tmp.dsds[2*6+1];
	    dsgs[2*15+ 7] = QC[1]*dsfs[0][2*10+7] + WQ[1]*dsfs[1][2*10+7]
			  +   tmp.dsds[2*6+2];
	    dsgs[2*15+ 8] = QC[2]*dsfs[0][2*10+8] + WQ[2]*dsfs[1][2*10+8]
			  +   tmp.dsds[2*6+0] + rze2[2]*psfs[2*10+8];
	    dsgs[2*15+ 9] = QC[0]*dsfs[0][2*10+1] + WQ[0]*dsfs[1][2*10+1];
	    dsgs[2*15+10] = QC[1]*dsfs[0][2*10+2] + WQ[1]*dsfs[1][2*10+2];
	    dsgs[2*15+11] = QC[2]*dsfs[0][2*10+0] + WQ[2]*dsfs[1][2*10+0]
			  + rze2[2]*psfs[2*10+0];
	    dsgs[2*15+12] = QC[0]*dsfs[0][2*10+9] + WQ[0]*dsfs[1][2*10+9]
			  +   tmp.dsds[2*6+4];
	    dsgs[2*15+13] = QC[1]*dsfs[0][2*10+9] + WQ[1]*dsfs[1][2*10+9]
			  +   tmp.dsds[2*6+5];
	    dsgs[2*15+14] = QC[2]*dsfs[0][2*10+9] + WQ[2]*dsfs[1][2*10+9]
			  +   tmp.dsds[2*6+3] + rze2[2]*psfs[2*10+9];
	    dsgs[3*15+ 0] = QC[0]*dsfs[0][3*10+0] + WQ[0]*dsfs[1][3*10+0]
			  + 3*tmp.dsds[3*6+0] + rze2[1]*psfs[1*10+0];
	    dsgs[3*15+ 1] = QC[1]*dsfs[0][3*10+1] + WQ[1]*dsfs[1][3*10+1]
			  + 3*tmp.dsds[3*6+1] + rze2[1]*psfs[0*10+1];
	    dsgs[3*15+ 2] = QC[2]*dsfs[0][3*10+2] + WQ[2]*dsfs[1][3*10+2]
			  + 3*tmp.dsds[3*6+2];
	    dsgs[3*15+ 3] = QC[0]*dsfs[0][3*10+3] + WQ[0]*dsfs[1][3*10+3]
			  + 2*tmp.dsds[3*6+3] + rze2[1]*psfs[1*10+3];
	    dsgs[3*15+ 4] = QC[1]*dsfs[0][3*10+4] + WQ[1]*dsfs[1][3*10+4]
			  + 2*tmp.dsds[3*6+4] + rze2[1]*psfs[0*10+4];
	    dsgs[3*15+ 5] = QC[2]*dsfs[0][3*10+5] + WQ[2]*dsfs[1][3*10+5]
			  + 2*tmp.dsds[3*6+5];
	    dsgs[3*15+ 6] = QC[0]*dsfs[0][3*10+6] + WQ[0]*dsfs[1][3*10+6]
			  +   tmp.dsds[3*6+1] + rze2[1]*psfs[1*10+6];
	    dsgs[3*15+ 7] = QC[1]*dsfs[0][3*10+7] + WQ[1]*dsfs[1][3*10+7]
			  +   tmp.dsds[3*6+2] + rze2[1]*psfs[0*10+7];
	    dsgs[3*15+ 8] = QC[2]*dsfs[0][3*10+8] + WQ[2]*dsfs[1][3*10+8]
			  +   tmp.dsds[3*6+0];
	    dsgs[3*15+ 9] = QC[0]*dsfs[0][3*10+1] + WQ[0]*dsfs[1][3*10+1]
			  + rze2[1]*psfs[1*10+1];
	    dsgs[3*15+10] = QC[1]*dsfs[0][3*10+2] + WQ[1]*dsfs[1][3*10+2]
			  + rze2[1]*psfs[0*10+2];
	    dsgs[3*15+11] = QC[2]*dsfs[0][3*10+0] + WQ[2]*dsfs[1][3*10+0];
	    dsgs[3*15+12] = QC[0]*dsfs[0][3*10+9] + WQ[0]*dsfs[1][3*10+9]
			  +   tmp.dsds[3*6+4] + rze2[1]*psfs[1*10+9];
	    dsgs[3*15+13] = QC[1]*dsfs[0][3*10+9] + WQ[1]*dsfs[1][3*10+9]
			  +   tmp.dsds[3*6+5] + rze2[1]*psfs[0*10+9];
	    dsgs[3*15+14] = QC[2]*dsfs[0][3*10+9] + WQ[2]*dsfs[1][3*10+9]
			  +   tmp.dsds[3*6+3];
	    dsgs[4*15+ 0] = QC[0]*dsfs[0][4*10+0] + WQ[0]*dsfs[1][4*10+0]
			  + 3*tmp.dsds[4*6+0];
	    dsgs[4*15+ 1] = QC[1]*dsfs[0][4*10+1] + WQ[1]*dsfs[1][4*10+1]
			  + 3*tmp.dsds[4*6+1] + rze2[1]*psfs[2*10+1];
	    dsgs[4*15+ 2] = QC[2]*dsfs[0][4*10+2] + WQ[2]*dsfs[1][4*10+2]
			  + 3*tmp.dsds[4*6+2] + rze2[1]*psfs[1*10+2];
	    dsgs[4*15+ 3] = QC[0]*dsfs[0][4*10+3] + WQ[0]*dsfs[1][4*10+3]
			  + 2*tmp.dsds[4*6+3];
	    dsgs[4*15+ 4] = QC[1]*dsfs[0][4*10+4] + WQ[1]*dsfs[1][4*10+4]
			  + 2*tmp.dsds[4*6+4] + rze2[1]*psfs[2*10+4];
	    dsgs[4*15+ 5] = QC[2]*dsfs[0][4*10+5] + WQ[2]*dsfs[1][4*10+5]
			  + 2*tmp.dsds[4*6+5] + rze2[1]*psfs[1*10+5];
	    dsgs[4*15+ 6] = QC[0]*dsfs[0][4*10+6] + WQ[0]*dsfs[1][4*10+6]
			  +   tmp.dsds[4*6+1];
	    dsgs[4*15+ 7] = QC[1]*dsfs[0][4*10+7] + WQ[1]*dsfs[1][4*10+7]
			  +   tmp.dsds[4*6+2] + rze2[1]*psfs[2*10+7];
	    dsgs[4*15+ 8] = QC[2]*dsfs[0][4*10+8] + WQ[2]*dsfs[1][4*10+8]
			  +   tmp.dsds[4*6+0] + rze2[1]*psfs[1*10+8];
	    dsgs[4*15+ 9] = QC[0]*dsfs[0][4*10+1] + WQ[0]*dsfs[1][4*10+1];
	    dsgs[4*15+10] = QC[1]*dsfs[0][4*10+2] + WQ[1]*dsfs[1][4*10+2]
			  + rze2[1]*psfs[2*10+2];
	    dsgs[4*15+11] = QC[2]*dsfs[0][4*10+0] + WQ[2]*dsfs[1][4*10+0]
			  + rze2[1]*psfs[1*10+0];
	    dsgs[4*15+12] = QC[0]*dsfs[0][4*10+9] + WQ[0]*dsfs[1][4*10+9]
			  +   tmp.dsds[4*6+4];
	    dsgs[4*15+13] = QC[1]*dsfs[0][4*10+9] + WQ[1]*dsfs[1][4*10+9]
			  +   tmp.dsds[4*6+5] + rze2[1]*psfs[2*10+9];
	    dsgs[4*15+14] = QC[2]*dsfs[0][4*10+9] + WQ[2]*dsfs[1][4*10+9]
			  +   tmp.dsds[4*6+3] + rze2[1]*psfs[1*10+9];
	    dsgs[5*15+ 0] = QC[0]*dsfs[0][5*10+0] + WQ[0]*dsfs[1][5*10+0]
			  + 3*tmp.dsds[5*6+0] + rze2[1]*psfs[2*10+0];
	    dsgs[5*15+ 1] = QC[1]*dsfs[0][5*10+1] + WQ[1]*dsfs[1][5*10+1]
			  + 3*tmp.dsds[5*6+1];
	    dsgs[5*15+ 2] = QC[2]*dsfs[0][5*10+2] + WQ[2]*dsfs[1][5*10+2]
			  + 3*tmp.dsds[5*6+2] + rze2[1]*psfs[0*10+2];
	    dsgs[5*15+ 3] = QC[0]*dsfs[0][5*10+3] + WQ[0]*dsfs[1][5*10+3]
			  + 2*tmp.dsds[5*6+3] + rze2[1]*psfs[2*10+3];
	    dsgs[5*15+ 4] = QC[1]*dsfs[0][5*10+4] + WQ[1]*dsfs[1][5*10+4]
			  + 2*tmp.dsds[5*6+4];
	    dsgs[5*15+ 5] = QC[2]*dsfs[0][5*10+5] + WQ[2]*dsfs[1][5*10+5]
			  + 2*tmp.dsds[5*6+5] + rze2[1]*psfs[0*10+5];
	    dsgs[5*15+ 6] = QC[0]*dsfs[0][5*10+6] + WQ[0]*dsfs[1][5*10+6]
			  +   tmp.dsds[5*6+1] + rze2[1]*psfs[2*10+6];
	    dsgs[5*15+ 7] = QC[1]*dsfs[0][5*10+7] + WQ[1]*dsfs[1][5*10+7]
			  +   tmp.dsds[5*6+2];
	    dsgs[5*15+ 8] = QC[2]*dsfs[0][5*10+8] + WQ[2]*dsfs[1][5*10+8]
			  +   tmp.dsds[5*6+0] + rze2[1]*psfs[0*10+8];
	    dsgs[5*15+ 9] = QC[0]*dsfs[0][5*10+1] + WQ[0]*dsfs[1][5*10+1]
			  + rze2[1]*psfs[2*10+1];
	    dsgs[5*15+10] = QC[1]*dsfs[0][5*10+2] + WQ[1]*dsfs[1][5*10+2];
	    dsgs[5*15+11] = QC[2]*dsfs[0][5*10+0] + WQ[2]*dsfs[1][5*10+0]
			  + rze2[1]*psfs[0*10+0];
	    dsgs[5*15+12] = QC[0]*dsfs[0][5*10+9] + WQ[0]*dsfs[1][5*10+9]
			  +   tmp.dsds[5*6+4] + rze2[1]*psfs[2*10+9];
	    dsgs[5*15+13] = QC[1]*dsfs[0][5*10+9] + WQ[1]*dsfs[1][5*10+9]
			  +   tmp.dsds[5*6+5];
	    dsgs[5*15+14] = QC[2]*dsfs[0][5*10+9] + WQ[2]*dsfs[1][5*10+9]
			  +   tmp.dsds[5*6+3] + rze2[1]*psfs[0*10+9];
	    // (F,S|G,S)
	    for (int i=0; i<10*6; i++)
		tmp.fsds[i] = eta2 * ( fsds[0][i] - re*fsds[1][i] );

	    fsgs[0*15+ 0] = QC[0]*fsfs[0][0*10+0] + WQ[0]*fsfs[1][0*10+0]
			  + 3*tmp.fsds[0*6+0] + rze2[3]*dsfs[1][0*10+0];
	    fsgs[0*15+ 1] = QC[1]*fsfs[0][0*10+1] + WQ[1]*fsfs[1][0*10+1]
			  + 3*tmp.fsds[0*6+1];
	    fsgs[0*15+ 2] = QC[2]*fsfs[0][0*10+2] + WQ[2]*fsfs[1][0*10+2]
			  + 3*tmp.fsds[0*6+2];
	    fsgs[0*15+ 3] = QC[0]*fsfs[0][0*10+3] + WQ[0]*fsfs[1][0*10+3]
			  + 2*tmp.fsds[0*6+3] + rze2[3]*dsfs[1][0*10+3];
	    fsgs[0*15+ 4] = QC[1]*fsfs[0][0*10+4] + WQ[1]*fsfs[1][0*10+4]
			  + 2*tmp.fsds[0*6+4];
	    fsgs[0*15+ 5] = QC[2]*fsfs[0][0*10+5] + WQ[2]*fsfs[1][0*10+5]
			  + 2*tmp.fsds[0*6+5];
	    fsgs[0*15+ 6] = QC[0]*fsfs[0][0*10+6] + WQ[0]*fsfs[1][0*10+6]
			  +   tmp.fsds[0*6+1] + rze2[3]*dsfs[1][0*10+6];
	    fsgs[0*15+ 7] = QC[1]*fsfs[0][0*10+7] + WQ[1]*fsfs[1][0*10+7]
			  +   tmp.fsds[0*6+2];
	    fsgs[0*15+ 8] = QC[2]*fsfs[0][0*10+8] + WQ[2]*fsfs[1][0*10+8]
			  +   tmp.fsds[0*6+0];
	    fsgs[0*15+ 9] = QC[0]*fsfs[0][0*10+1] + WQ[0]*fsfs[1][0*10+1]
			  + rze2[3]*dsfs[1][0*10+1];
	    fsgs[0*15+10] = QC[1]*fsfs[0][0*10+2] + WQ[1]*fsfs[1][0*10+2];
	    fsgs[0*15+11] = QC[2]*fsfs[0][0*10+0] + WQ[2]*fsfs[1][0*10+0];
	    fsgs[0*15+12] = QC[0]*fsfs[0][0*10+9] + WQ[0]*fsfs[1][0*10+9]
			  +   tmp.fsds[0*6+4] + rze2[3]*dsfs[1][0*10+9];
	    fsgs[0*15+13] = QC[1]*fsfs[0][0*10+9] + WQ[1]*fsfs[1][0*10+9]
			  +   tmp.fsds[0*6+5];
	    fsgs[0*15+14] = QC[2]*fsfs[0][0*10+9] + WQ[2]*fsfs[1][0*10+9]
			  +   tmp.fsds[0*6+3];
	    fsgs[1*15+ 0] = QC[0]*fsfs[0][1*10+0] + WQ[0]*fsfs[1][1*10+0]
			  + 3*tmp.fsds[1*6+0];
	    fsgs[1*15+ 1] = QC[1]*fsfs[0][1*10+1] + WQ[1]*fsfs[1][1*10+1]
			  + 3*tmp.fsds[1*6+1] + rze2[3]*dsfs[1][1*10+1];
	    fsgs[1*15+ 2] = QC[2]*fsfs[0][1*10+2] + WQ[2]*fsfs[1][1*10+2]
			  + 3*tmp.fsds[1*6+2];
	    fsgs[1*15+ 3] = QC[0]*fsfs[0][1*10+3] + WQ[0]*fsfs[1][1*10+3]
			  + 2*tmp.fsds[1*6+3];
	    fsgs[1*15+ 4] = QC[1]*fsfs[0][1*10+4] + WQ[1]*fsfs[1][1*10+4]
			  + 2*tmp.fsds[1*6+4] + rze2[3]*dsfs[1][1*10+4];
	    fsgs[1*15+ 5] = QC[2]*fsfs[0][1*10+5] + WQ[2]*fsfs[1][1*10+5]
			  + 2*tmp.fsds[1*6+5];
	    fsgs[1*15+ 6] = QC[0]*fsfs[0][1*10+6] + WQ[0]*fsfs[1][1*10+6]
			  +   tmp.fsds[1*6+1];
	    fsgs[1*15+ 7] = QC[1]*fsfs[0][1*10+7] + WQ[1]*fsfs[1][1*10+7]
			  +   tmp.fsds[1*6+2] + rze2[3]*dsfs[1][1*10+7];
	    fsgs[1*15+ 8] = QC[2]*fsfs[0][1*10+8] + WQ[2]*fsfs[1][1*10+8]
			  +   tmp.fsds[1*6+0];
	    fsgs[1*15+ 9] = QC[0]*fsfs[0][1*10+1] + WQ[0]*fsfs[1][1*10+1];
	    fsgs[1*15+10] = QC[1]*fsfs[0][1*10+2] + WQ[1]*fsfs[1][1*10+2]
			  + rze2[3]*dsfs[1][1*10+2];
	    fsgs[1*15+11] = QC[2]*fsfs[0][1*10+0] + WQ[2]*fsfs[1][1*10+0];
	    fsgs[1*15+12] = QC[0]*fsfs[0][1*10+9] + WQ[0]*fsfs[1][1*10+9]
			  +   tmp.fsds[1*6+4];
	    fsgs[1*15+13] = QC[1]*fsfs[0][1*10+9] + WQ[1]*fsfs[1][1*10+9]
			  +   tmp.fsds[1*6+5] + rze2[3]*dsfs[1][1*10+9];
	    fsgs[1*15+14] = QC[2]*fsfs[0][1*10+9] + WQ[2]*fsfs[1][1*10+9]
			  +   tmp.fsds[1*6+3];
	    fsgs[2*15+ 0] = QC[0]*fsfs[0][2*10+0] + WQ[0]*fsfs[1][2*10+0]
			  + 3*tmp.fsds[2*6+0];
	    fsgs[2*15+ 1] = QC[1]*fsfs[0][2*10+1] + WQ[1]*fsfs[1][2*10+1]
			  + 3*tmp.fsds[2*6+1];
	    fsgs[2*15+ 2] = QC[2]*fsfs[0][2*10+2] + WQ[2]*fsfs[1][2*10+2]
			  + 3*tmp.fsds[2*6+2] + rze2[3]*dsfs[1][2*10+2];
	    fsgs[2*15+ 3] = QC[0]*fsfs[0][2*10+3] + WQ[0]*fsfs[1][2*10+3]
			  + 2*tmp.fsds[2*6+3];
	    fsgs[2*15+ 4] = QC[1]*fsfs[0][2*10+4] + WQ[1]*fsfs[1][2*10+4]
			  + 2*tmp.fsds[2*6+4];
	    fsgs[2*15+ 5] = QC[2]*fsfs[0][2*10+5] + WQ[2]*fsfs[1][2*10+5]
			  + 2*tmp.fsds[2*6+5] + rze2[3]*dsfs[1][2*10+5];
	    fsgs[2*15+ 6] = QC[0]*fsfs[0][2*10+6] + WQ[0]*fsfs[1][2*10+6]
			  +   tmp.fsds[2*6+1];
	    fsgs[2*15+ 7] = QC[1]*fsfs[0][2*10+7] + WQ[1]*fsfs[1][2*10+7]
			  +   tmp.fsds[2*6+2];
	    fsgs[2*15+ 8] = QC[2]*fsfs[0][2*10+8] + WQ[2]*fsfs[1][2*10+8]
			  +   tmp.fsds[2*6+0] + rze2[3]*dsfs[1][2*10+8];
	    fsgs[2*15+ 9] = QC[0]*fsfs[0][2*10+1] + WQ[0]*fsfs[1][2*10+1];
	    fsgs[2*15+10] = QC[1]*fsfs[0][2*10+2] + WQ[1]*fsfs[1][2*10+2];
	    fsgs[2*15+11] = QC[2]*fsfs[0][2*10+0] + WQ[2]*fsfs[1][2*10+0]
			  + rze2[3]*dsfs[1][2*10+0];
	    fsgs[2*15+12] = QC[0]*fsfs[0][2*10+9] + WQ[0]*fsfs[1][2*10+9]
			  +   tmp.fsds[2*6+4];
	    fsgs[2*15+13] = QC[1]*fsfs[0][2*10+9] + WQ[1]*fsfs[1][2*10+9]
			  +   tmp.fsds[2*6+5];
	    fsgs[2*15+14] = QC[2]*fsfs[0][2*10+9] + WQ[2]*fsfs[1][2*10+9]
			  +   tmp.fsds[2*6+3] + rze2[3]*dsfs[1][2*10+9];
	    fsgs[3*15+ 0] = QC[0]*fsfs[0][3*10+0] + WQ[0]*fsfs[1][3*10+0]
			  + 3*tmp.fsds[3*6+0] + rze2[2]*dsfs[1][3*10+0];
	    fsgs[3*15+ 1] = QC[1]*fsfs[0][3*10+1] + WQ[1]*fsfs[1][3*10+1]
			  + 3*tmp.fsds[3*6+1] + rze2[1]*dsfs[1][0*10+1];
	    fsgs[3*15+ 2] = QC[2]*fsfs[0][3*10+2] + WQ[2]*fsfs[1][3*10+2]
			  + 3*tmp.fsds[3*6+2];
	    fsgs[3*15+ 3] = QC[0]*fsfs[0][3*10+3] + WQ[0]*fsfs[1][3*10+3]
			  + 2*tmp.fsds[3*6+3] + rze2[2]*dsfs[1][3*10+3];
	    fsgs[3*15+ 4] = QC[1]*fsfs[0][3*10+4] + WQ[1]*fsfs[1][3*10+4]
			  + 2*tmp.fsds[3*6+4] + rze2[1]*dsfs[1][0*10+4];
	    fsgs[3*15+ 5] = QC[2]*fsfs[0][3*10+5] + WQ[2]*fsfs[1][3*10+5]
			  + 2*tmp.fsds[3*6+5];
	    fsgs[3*15+ 6] = QC[0]*fsfs[0][3*10+6] + WQ[0]*fsfs[1][3*10+6]
			  +   tmp.fsds[3*6+1] + rze2[2]*dsfs[1][3*10+6];
	    fsgs[3*15+ 7] = QC[1]*fsfs[0][3*10+7] + WQ[1]*fsfs[1][3*10+7]
			  +   tmp.fsds[3*6+2] + rze2[1]*dsfs[1][0*10+7];
	    fsgs[3*15+ 8] = QC[2]*fsfs[0][3*10+8] + WQ[2]*fsfs[1][3*10+8]
			  +   tmp.fsds[3*6+0];
	    fsgs[3*15+ 9] = QC[0]*fsfs[0][3*10+1] + WQ[0]*fsfs[1][3*10+1]
			  + rze2[2]*dsfs[1][3*10+1];
	    fsgs[3*15+10] = QC[1]*fsfs[0][3*10+2] + WQ[1]*fsfs[1][3*10+2]
			  + rze2[1]*dsfs[1][0*10+2];
	    fsgs[3*15+11] = QC[2]*fsfs[0][3*10+0] + WQ[2]*fsfs[1][3*10+0];
	    fsgs[3*15+12] = QC[0]*fsfs[0][3*10+9] + WQ[0]*fsfs[1][3*10+9]
			  +   tmp.fsds[3*6+4] + rze2[2]*dsfs[1][3*10+9];
	    fsgs[3*15+13] = QC[1]*fsfs[0][3*10+9] + WQ[1]*fsfs[1][3*10+9]
			  +   tmp.fsds[3*6+5] + rze2[1]*dsfs[1][0*10+9];
	    fsgs[3*15+14] = QC[2]*fsfs[0][3*10+9] + WQ[2]*fsfs[1][3*10+9]
			  +   tmp.fsds[3*6+3];
	    fsgs[4*15+ 0] = QC[0]*fsfs[0][4*10+0] + WQ[0]*fsfs[1][4*10+0]
			  + 3*tmp.fsds[4*6+0];
	    fsgs[4*15+ 1] = QC[1]*fsfs[0][4*10+1] + WQ[1]*fsfs[1][4*10+1]
			  + 3*tmp.fsds[4*6+1] + rze2[2]*dsfs[1][4*10+1];
	    fsgs[4*15+ 2] = QC[2]*fsfs[0][4*10+2] + WQ[2]*fsfs[1][4*10+2]
			  + 3*tmp.fsds[4*6+2] + rze2[1]*dsfs[1][1*10+2];
	    fsgs[4*15+ 3] = QC[0]*fsfs[0][4*10+3] + WQ[0]*fsfs[1][4*10+3]
			  + 2*tmp.fsds[4*6+3];
	    fsgs[4*15+ 4] = QC[1]*fsfs[0][4*10+4] + WQ[1]*fsfs[1][4*10+4]
			  + 2*tmp.fsds[4*6+4] + rze2[2]*dsfs[1][4*10+4];
	    fsgs[4*15+ 5] = QC[2]*fsfs[0][4*10+5] + WQ[2]*fsfs[1][4*10+5]
			  + 2*tmp.fsds[4*6+5] + rze2[1]*dsfs[1][1*10+5];
	    fsgs[4*15+ 6] = QC[0]*fsfs[0][4*10+6] + WQ[0]*fsfs[1][4*10+6]
			  +   tmp.fsds[4*6+1];
	    fsgs[4*15+ 7] = QC[1]*fsfs[0][4*10+7] + WQ[1]*fsfs[1][4*10+7]
			  +   tmp.fsds[4*6+2] + rze2[2]*dsfs[1][4*10+7];
	    fsgs[4*15+ 8] = QC[2]*fsfs[0][4*10+8] + WQ[2]*fsfs[1][4*10+8]
			  +   tmp.fsds[4*6+0] + rze2[1]*dsfs[1][1*10+8];
	    fsgs[4*15+ 9] = QC[0]*fsfs[0][4*10+1] + WQ[0]*fsfs[1][4*10+1];
	    fsgs[4*15+10] = QC[1]*fsfs[0][4*10+2] + WQ[1]*fsfs[1][4*10+2]
			  + rze2[2]*dsfs[1][4*10+2];
	    fsgs[4*15+11] = QC[2]*fsfs[0][4*10+0] + WQ[2]*fsfs[1][4*10+0]
			  + rze2[1]*dsfs[1][1*10+0];
	    fsgs[4*15+12] = QC[0]*fsfs[0][4*10+9] + WQ[0]*fsfs[1][4*10+9]
			  +   tmp.fsds[4*6+4];
	    fsgs[4*15+13] = QC[1]*fsfs[0][4*10+9] + WQ[1]*fsfs[1][4*10+9]
			  +   tmp.fsds[4*6+5] + rze2[2]*dsfs[1][4*10+9];
	    fsgs[4*15+14] = QC[2]*fsfs[0][4*10+9] + WQ[2]*fsfs[1][4*10+9]
			  +   tmp.fsds[4*6+3] + rze2[1]*dsfs[1][1*10+9];
	    fsgs[5*15+ 0] = QC[0]*fsfs[0][5*10+0] + WQ[0]*fsfs[1][5*10+0]
			  + 3*tmp.fsds[5*6+0] + rze2[1]*dsfs[1][2*10+0];
	    fsgs[5*15+ 1] = QC[1]*fsfs[0][5*10+1] + WQ[1]*fsfs[1][5*10+1]
			  + 3*tmp.fsds[5*6+1];
	    fsgs[5*15+ 2] = QC[2]*fsfs[0][5*10+2] + WQ[2]*fsfs[1][5*10+2]
			  + 3*tmp.fsds[5*6+2] + rze2[2]*dsfs[1][5*10+2];
	    fsgs[5*15+ 3] = QC[0]*fsfs[0][5*10+3] + WQ[0]*fsfs[1][5*10+3]
			  + 2*tmp.fsds[5*6+3] + rze2[1]*dsfs[1][2*10+3];
	    fsgs[5*15+ 4] = QC[1]*fsfs[0][5*10+4] + WQ[1]*fsfs[1][5*10+4]
			  + 2*tmp.fsds[5*6+4];
	    fsgs[5*15+ 5] = QC[2]*fsfs[0][5*10+5] + WQ[2]*fsfs[1][5*10+5]
			  + 2*tmp.fsds[5*6+5] + rze2[2]*dsfs[1][5*10+5];
	    fsgs[5*15+ 6] = QC[0]*fsfs[0][5*10+6] + WQ[0]*fsfs[1][5*10+6]
			  +   tmp.fsds[5*6+1] + rze2[1]*dsfs[1][2*10+6];
	    fsgs[5*15+ 7] = QC[1]*fsfs[0][5*10+7] + WQ[1]*fsfs[1][5*10+7]
			  +   tmp.fsds[5*6+2];
	    fsgs[5*15+ 8] = QC[2]*fsfs[0][5*10+8] + WQ[2]*fsfs[1][5*10+8]
			  +   tmp.fsds[5*6+0] + rze2[2]*dsfs[1][5*10+8];
	    fsgs[5*15+ 9] = QC[0]*fsfs[0][5*10+1] + WQ[0]*fsfs[1][5*10+1]
			  + rze2[1]*dsfs[1][2*10+1];
	    fsgs[5*15+10] = QC[1]*fsfs[0][5*10+2] + WQ[1]*fsfs[1][5*10+2];
	    fsgs[5*15+11] = QC[2]*fsfs[0][5*10+0] + WQ[2]*fsfs[1][5*10+0]
			  + rze2[2]*dsfs[1][5*10+0];
	    fsgs[5*15+12] = QC[0]*fsfs[0][5*10+9] + WQ[0]*fsfs[1][5*10+9]
			  +   tmp.fsds[5*6+4] + rze2[1]*dsfs[1][2*10+9];
	    fsgs[5*15+13] = QC[1]*fsfs[0][5*10+9] + WQ[1]*fsfs[1][5*10+9]
			  +   tmp.fsds[5*6+5];
	    fsgs[5*15+14] = QC[2]*fsfs[0][5*10+9] + WQ[2]*fsfs[1][5*10+9]
			  +   tmp.fsds[5*6+3] + rze2[2]*dsfs[1][5*10+9];
	    fsgs[6*15+ 0] = QC[0]*fsfs[0][6*10+0] + WQ[0]*fsfs[1][6*10+0]
			  + 3*tmp.fsds[6*6+0] + rze2[1]*dsfs[1][1*10+0];
	    fsgs[6*15+ 1] = QC[1]*fsfs[0][6*10+1] + WQ[1]*fsfs[1][6*10+1]
			  + 3*tmp.fsds[6*6+1] + rze2[2]*dsfs[1][3*10+1];
	    fsgs[6*15+ 2] = QC[2]*fsfs[0][6*10+2] + WQ[2]*fsfs[1][6*10+2]
			  + 3*tmp.fsds[6*6+2];
	    fsgs[6*15+ 3] = QC[0]*fsfs[0][6*10+3] + WQ[0]*fsfs[1][6*10+3]
			  + 2*tmp.fsds[6*6+3] + rze2[1]*dsfs[1][1*10+3];
	    fsgs[6*15+ 4] = QC[1]*fsfs[0][6*10+4] + WQ[1]*fsfs[1][6*10+4]
			  + 2*tmp.fsds[6*6+4] + rze2[2]*dsfs[1][3*10+4];
	    fsgs[6*15+ 5] = QC[2]*fsfs[0][6*10+5] + WQ[2]*fsfs[1][6*10+5]
			  + 2*tmp.fsds[6*6+5];
	    fsgs[6*15+ 6] = QC[0]*fsfs[0][6*10+6] + WQ[0]*fsfs[1][6*10+6]
			  +   tmp.fsds[6*6+1] + rze2[1]*dsfs[1][1*10+6];
	    fsgs[6*15+ 7] = QC[1]*fsfs[0][6*10+7] + WQ[1]*fsfs[1][6*10+7]
			  +   tmp.fsds[6*6+2] + rze2[2]*dsfs[1][3*10+7];
	    fsgs[6*15+ 8] = QC[2]*fsfs[0][6*10+8] + WQ[2]*fsfs[1][6*10+8]
			  +   tmp.fsds[6*6+0];
	    fsgs[6*15+ 9] = QC[0]*fsfs[0][6*10+1] + WQ[0]*fsfs[1][6*10+1]
			  + rze2[1]*dsfs[1][1*10+1];
	    fsgs[6*15+10] = QC[1]*fsfs[0][6*10+2] + WQ[1]*fsfs[1][6*10+2]
			  + rze2[2]*dsfs[1][3*10+2];
	    fsgs[6*15+11] = QC[2]*fsfs[0][6*10+0] + WQ[2]*fsfs[1][6*10+0];
	    fsgs[6*15+12] = QC[0]*fsfs[0][6*10+9] + WQ[0]*fsfs[1][6*10+9]
			  +   tmp.fsds[6*6+4] + rze2[1]*dsfs[1][1*10+9];
	    fsgs[6*15+13] = QC[1]*fsfs[0][6*10+9] + WQ[1]*fsfs[1][6*10+9]
			  +   tmp.fsds[6*6+5] + rze2[2]*dsfs[1][3*10+9];
	    fsgs[6*15+14] = QC[2]*fsfs[0][6*10+9] + WQ[2]*fsfs[1][6*10+9]
			  +   tmp.fsds[6*6+3];
	    fsgs[7*15+ 0] = QC[0]*fsfs[0][7*10+0] + WQ[0]*fsfs[1][7*10+0]
			  + 3*tmp.fsds[7*6+0];
	    fsgs[7*15+ 1] = QC[1]*fsfs[0][7*10+1] + WQ[1]*fsfs[1][7*10+1]
			  + 3*tmp.fsds[7*6+1] + rze2[1]*dsfs[1][2*10+1];
	    fsgs[7*15+ 2] = QC[2]*fsfs[0][7*10+2] + WQ[2]*fsfs[1][7*10+2]
			  + 3*tmp.fsds[7*6+2] + rze2[2]*dsfs[1][4*10+2];
	    fsgs[7*15+ 3] = QC[0]*fsfs[0][7*10+3] + WQ[0]*fsfs[1][7*10+3]
			  + 2*tmp.fsds[7*6+3];
	    fsgs[7*15+ 4] = QC[1]*fsfs[0][7*10+4] + WQ[1]*fsfs[1][7*10+4]
			  + 2*tmp.fsds[7*6+4] + rze2[1]*dsfs[1][2*10+4];
	    fsgs[7*15+ 5] = QC[2]*fsfs[0][7*10+5] + WQ[2]*fsfs[1][7*10+5]
			  + 2*tmp.fsds[7*6+5] + rze2[2]*dsfs[1][4*10+5];
	    fsgs[7*15+ 6] = QC[0]*fsfs[0][7*10+6] + WQ[0]*fsfs[1][7*10+6]
			  +   tmp.fsds[7*6+1];
	    fsgs[7*15+ 7] = QC[1]*fsfs[0][7*10+7] + WQ[1]*fsfs[1][7*10+7]
			  +   tmp.fsds[7*6+2] + rze2[1]*dsfs[1][2*10+7];
	    fsgs[7*15+ 8] = QC[2]*fsfs[0][7*10+8] + WQ[2]*fsfs[1][7*10+8]
			  +   tmp.fsds[7*6+0] + rze2[2]*dsfs[1][4*10+8];
	    fsgs[7*15+ 9] = QC[0]*fsfs[0][7*10+1] + WQ[0]*fsfs[1][7*10+1];
	    fsgs[7*15+10] = QC[1]*fsfs[0][7*10+2] + WQ[1]*fsfs[1][7*10+2]
			  + rze2[1]*dsfs[1][2*10+2];
	    fsgs[7*15+11] = QC[2]*fsfs[0][7*10+0] + WQ[2]*fsfs[1][7*10+0]
			  + rze2[2]*dsfs[1][4*10+0];
	    fsgs[7*15+12] = QC[0]*fsfs[0][7*10+9] + WQ[0]*fsfs[1][7*10+9]
			  +   tmp.fsds[7*6+4];
	    fsgs[7*15+13] = QC[1]*fsfs[0][7*10+9] + WQ[1]*fsfs[1][7*10+9]
			  +   tmp.fsds[7*6+5] + rze2[1]*dsfs[1][2*10+9];
	    fsgs[7*15+14] = QC[2]*fsfs[0][7*10+9] + WQ[2]*fsfs[1][7*10+9]
			  +   tmp.fsds[7*6+3] + rze2[2]*dsfs[1][4*10+9];
	    fsgs[8*15+ 0] = QC[0]*fsfs[0][8*10+0] + WQ[0]*fsfs[1][8*10+0]
			  + 3*tmp.fsds[8*6+0] + rze2[2]*dsfs[1][5*10+0];
	    fsgs[8*15+ 1] = QC[1]*fsfs[0][8*10+1] + WQ[1]*fsfs[1][8*10+1]
			  + 3*tmp.fsds[8*6+1];
	    fsgs[8*15+ 2] = QC[2]*fsfs[0][8*10+2] + WQ[2]*fsfs[1][8*10+2]
			  + 3*tmp.fsds[8*6+2] + rze2[1]*dsfs[1][0*10+2];
	    fsgs[8*15+ 3] = QC[0]*fsfs[0][8*10+3] + WQ[0]*fsfs[1][8*10+3]
			  + 2*tmp.fsds[8*6+3] + rze2[2]*dsfs[1][5*10+3];
	    fsgs[8*15+ 4] = QC[1]*fsfs[0][8*10+4] + WQ[1]*fsfs[1][8*10+4]
			  + 2*tmp.fsds[8*6+4];
	    fsgs[8*15+ 5] = QC[2]*fsfs[0][8*10+5] + WQ[2]*fsfs[1][8*10+5]
			  + 2*tmp.fsds[8*6+5] + rze2[1]*dsfs[1][0*10+5];
	    fsgs[8*15+ 6] = QC[0]*fsfs[0][8*10+6] + WQ[0]*fsfs[1][8*10+6]
			  +   tmp.fsds[8*6+1] + rze2[2]*dsfs[1][5*10+6];
	    fsgs[8*15+ 7] = QC[1]*fsfs[0][8*10+7] + WQ[1]*fsfs[1][8*10+7]
			  +   tmp.fsds[8*6+2];
	    fsgs[8*15+ 8] = QC[2]*fsfs[0][8*10+8] + WQ[2]*fsfs[1][8*10+8]
			  +   tmp.fsds[8*6+0] + rze2[1]*dsfs[1][0*10+8];
	    fsgs[8*15+ 9] = QC[0]*fsfs[0][8*10+1] + WQ[0]*fsfs[1][8*10+1]
			  + rze2[2]*dsfs[1][5*10+1];
	    fsgs[8*15+10] = QC[1]*fsfs[0][8*10+2] + WQ[1]*fsfs[1][8*10+2];
	    fsgs[8*15+11] = QC[2]*fsfs[0][8*10+0] + WQ[2]*fsfs[1][8*10+0]
			  + rze2[1]*dsfs[1][0*10+0];
	    fsgs[8*15+12] = QC[0]*fsfs[0][8*10+9] + WQ[0]*fsfs[1][8*10+9]
			  +   tmp.fsds[8*6+4] + rze2[2]*dsfs[1][5*10+9];
	    fsgs[8*15+13] = QC[1]*fsfs[0][8*10+9] + WQ[1]*fsfs[1][8*10+9]
			  +   tmp.fsds[8*6+5];
	    fsgs[8*15+14] = QC[2]*fsfs[0][8*10+9] + WQ[2]*fsfs[1][8*10+9]
			  +   tmp.fsds[8*6+3] + rze2[1]*dsfs[1][0*10+9];
	    fsgs[9*15+ 0] = QC[0]*fsfs[0][9*10+0] + WQ[0]*fsfs[1][9*10+0]
			  + 3*tmp.fsds[9*6+0] + rze2[1]*dsfs[1][4*10+0];
	    fsgs[9*15+ 1] = QC[1]*fsfs[0][9*10+1] + WQ[1]*fsfs[1][9*10+1]
			  + 3*tmp.fsds[9*6+1] + rze2[1]*dsfs[1][5*10+1];
	    fsgs[9*15+ 2] = QC[2]*fsfs[0][9*10+2] + WQ[2]*fsfs[1][9*10+2]
			  + 3*tmp.fsds[9*6+2] + rze2[1]*dsfs[1][3*10+2];
	    fsgs[9*15+ 3] = QC[0]*fsfs[0][9*10+3] + WQ[0]*fsfs[1][9*10+3]
			  + 2*tmp.fsds[9*6+3] + rze2[1]*dsfs[1][4*10+3];
	    fsgs[9*15+ 4] = QC[1]*fsfs[0][9*10+4] + WQ[1]*fsfs[1][9*10+4]
			  + 2*tmp.fsds[9*6+4] + rze2[1]*dsfs[1][5*10+4];
	    fsgs[9*15+ 5] = QC[2]*fsfs[0][9*10+5] + WQ[2]*fsfs[1][9*10+5]
			  + 2*tmp.fsds[9*6+5] + rze2[1]*dsfs[1][3*10+5];
	    fsgs[9*15+ 6] = QC[0]*fsfs[0][9*10+6] + WQ[0]*fsfs[1][9*10+6]
			  +   tmp.fsds[9*6+1] + rze2[1]*dsfs[1][4*10+6];
	    fsgs[9*15+ 7] = QC[1]*fsfs[0][9*10+7] + WQ[1]*fsfs[1][9*10+7]
			  +   tmp.fsds[9*6+2] + rze2[1]*dsfs[1][5*10+7];
	    fsgs[9*15+ 8] = QC[2]*fsfs[0][9*10+8] + WQ[2]*fsfs[1][9*10+8]
			  +   tmp.fsds[9*6+0] + rze2[1]*dsfs[1][3*10+8];
	    fsgs[9*15+ 9] = QC[0]*fsfs[0][9*10+1] + WQ[0]*fsfs[1][9*10+1]
			  + rze2[1]*dsfs[1][4*10+1];
	    fsgs[9*15+10] = QC[1]*fsfs[0][9*10+2] + WQ[1]*fsfs[1][9*10+2]
			  + rze2[1]*dsfs[1][5*10+2];
	    fsgs[9*15+11] = QC[2]*fsfs[0][9*10+0] + WQ[2]*fsfs[1][9*10+0]
			  + rze2[1]*dsfs[1][3*10+0];
	    fsgs[9*15+12] = QC[0]*fsfs[0][9*10+9] + WQ[0]*fsfs[1][9*10+9]
			  +   tmp.fsds[9*6+4] + rze2[1]*dsfs[1][4*10+9];
	    fsgs[9*15+13] = QC[1]*fsfs[0][9*10+9] + WQ[1]*fsfs[1][9*10+9]
			  +   tmp.fsds[9*6+5] + rze2[1]*dsfs[1][5*10+9];
	    fsgs[9*15+14] = QC[2]*fsfs[0][9*10+9] + WQ[2]*fsfs[1][9*10+9]
			  +   tmp.fsds[9*6+3] + rze2[1]*dsfs[1][3*10+9];
	    // (G,S|G,S)
	    for (int i=0; i<15*6; i++)
		tmp.gsds[i] = eta2 * ( gsds[0][i] - re*gsds[1][i] );

	    gsgs[ 0*15+ 0] = QC[0]*gsfs[0][ 0*10+0] + WQ[0]*gsfs[1][ 0*10+0]
			   + 3*tmp.gsds[ 0*6+0] + rze2[4]*fsfs[1][0*10+0];
	    gsgs[ 0*15+ 1] = QC[1]*gsfs[0][ 0*10+1] + WQ[1]*gsfs[1][ 0*10+1]
			   + 3*tmp.gsds[ 0*6+1];
	    gsgs[ 0*15+ 2] = QC[2]*gsfs[0][ 0*10+2] + WQ[2]*gsfs[1][ 0*10+2]
			   + 3*tmp.gsds[ 0*6+2];
	    gsgs[ 0*15+ 3] = QC[0]*gsfs[0][ 0*10+3] + WQ[0]*gsfs[1][ 0*10+3]
			   + 2*tmp.gsds[ 0*6+3] + rze2[4]*fsfs[1][0*10+3];
	    gsgs[ 0*15+ 4] = QC[1]*gsfs[0][ 0*10+4] + WQ[1]*gsfs[1][ 0*10+4]
			   + 2*tmp.gsds[ 0*6+4];
	    gsgs[ 0*15+ 5] = QC[2]*gsfs[0][ 0*10+5] + WQ[2]*gsfs[1][ 0*10+5]
			   + 2*tmp.gsds[ 0*6+5];
	    gsgs[ 0*15+ 6] = QC[0]*gsfs[0][ 0*10+6] + WQ[0]*gsfs[1][ 0*10+6]
			   +   tmp.gsds[ 0*6+1] + rze2[4]*fsfs[1][0*10+6];
	    gsgs[ 0*15+ 7] = QC[1]*gsfs[0][ 0*10+7] + WQ[1]*gsfs[1][ 0*10+7]
			   +   tmp.gsds[ 0*6+2];
	    gsgs[ 0*15+ 8] = QC[2]*gsfs[0][ 0*10+8] + WQ[2]*gsfs[1][ 0*10+8]
			   +   tmp.gsds[ 0*6+0];
	    gsgs[ 0*15+ 9] = QC[0]*gsfs[0][ 0*10+1] + WQ[0]*gsfs[1][ 0*10+1]
			   + rze2[4]*fsfs[1][0*10+1];
	    gsgs[ 0*15+10] = QC[1]*gsfs[0][ 0*10+2] + WQ[1]*gsfs[1][ 0*10+2];
	    gsgs[ 0*15+11] = QC[2]*gsfs[0][ 0*10+0] + WQ[2]*gsfs[1][ 0*10+0];
	    gsgs[ 0*15+12] = QC[0]*gsfs[0][ 0*10+9] + WQ[0]*gsfs[1][ 0*10+9]
			   +   tmp.gsds[ 0*6+4] + rze2[4]*fsfs[1][0*10+9];
	    gsgs[ 0*15+13] = QC[1]*gsfs[0][ 0*10+9] + WQ[1]*gsfs[1][ 0*10+9]
			   +   tmp.gsds[ 0*6+5];
	    gsgs[ 0*15+14] = QC[2]*gsfs[0][ 0*10+9] + WQ[2]*gsfs[1][ 0*10+9]
			   +   tmp.gsds[ 0*6+3];
	    gsgs[ 1*15+ 0] = QC[0]*gsfs[0][ 1*10+0] + WQ[0]*gsfs[1][ 1*10+0]
			   + 3*tmp.gsds[ 1*6+0];
	    gsgs[ 1*15+ 1] = QC[1]*gsfs[0][ 1*10+1] + WQ[1]*gsfs[1][ 1*10+1]
			   + 3*tmp.gsds[ 1*6+1] + rze2[4]*fsfs[1][1*10+1];
	    gsgs[ 1*15+ 2] = QC[2]*gsfs[0][ 1*10+2] + WQ[2]*gsfs[1][ 1*10+2]
			   + 3*tmp.gsds[ 1*6+2];
	    gsgs[ 1*15+ 3] = QC[0]*gsfs[0][ 1*10+3] + WQ[0]*gsfs[1][ 1*10+3]
			   + 2*tmp.gsds[ 1*6+3];
	    gsgs[ 1*15+ 4] = QC[1]*gsfs[0][ 1*10+4] + WQ[1]*gsfs[1][ 1*10+4]
			   + 2*tmp.gsds[ 1*6+4] + rze2[4]*fsfs[1][1*10+4];
	    gsgs[ 1*15+ 5] = QC[2]*gsfs[0][ 1*10+5] + WQ[2]*gsfs[1][ 1*10+5]
			   + 2*tmp.gsds[ 1*6+5];
	    gsgs[ 1*15+ 6] = QC[0]*gsfs[0][ 1*10+6] + WQ[0]*gsfs[1][ 1*10+6]
			   +   tmp.gsds[ 1*6+1];
	    gsgs[ 1*15+ 7] = QC[1]*gsfs[0][ 1*10+7] + WQ[1]*gsfs[1][ 1*10+7]
			   +   tmp.gsds[ 1*6+2] + rze2[4]*fsfs[1][1*10+7];
	    gsgs[ 1*15+ 8] = QC[2]*gsfs[0][ 1*10+8] + WQ[2]*gsfs[1][ 1*10+8]
			   +   tmp.gsds[ 1*6+0];
	    gsgs[ 1*15+ 9] = QC[0]*gsfs[0][ 1*10+1] + WQ[0]*gsfs[1][ 1*10+1];
	    gsgs[ 1*15+10] = QC[1]*gsfs[0][ 1*10+2] + WQ[1]*gsfs[1][ 1*10+2]
			   + rze2[4]*fsfs[1][1*10+2];
	    gsgs[ 1*15+11] = QC[2]*gsfs[0][ 1*10+0] + WQ[2]*gsfs[1][ 1*10+0];
	    gsgs[ 1*15+12] = QC[0]*gsfs[0][ 1*10+9] + WQ[0]*gsfs[1][ 1*10+9]
			   +   tmp.gsds[ 1*6+4];
	    gsgs[ 1*15+13] = QC[1]*gsfs[0][ 1*10+9] + WQ[1]*gsfs[1][ 1*10+9]
			   +   tmp.gsds[ 1*6+5] + rze2[4]*fsfs[1][1*10+9];
	    gsgs[ 1*15+14] = QC[2]*gsfs[0][ 1*10+9] + WQ[2]*gsfs[1][ 1*10+9]
			   +   tmp.gsds[ 1*6+3];
	    gsgs[ 2*15+ 0] = QC[0]*gsfs[0][ 2*10+0] + WQ[0]*gsfs[1][ 2*10+0]
			   + 3*tmp.gsds[ 2*6+0];
	    gsgs[ 2*15+ 1] = QC[1]*gsfs[0][ 2*10+1] + WQ[1]*gsfs[1][ 2*10+1]
			   + 3*tmp.gsds[ 2*6+1];
	    gsgs[ 2*15+ 2] = QC[2]*gsfs[0][ 2*10+2] + WQ[2]*gsfs[1][ 2*10+2]
			   + 3*tmp.gsds[ 2*6+2] + rze2[4]*fsfs[1][2*10+2];
	    gsgs[ 2*15+ 3] = QC[0]*gsfs[0][ 2*10+3] + WQ[0]*gsfs[1][ 2*10+3]
			   + 2*tmp.gsds[ 2*6+3];
	    gsgs[ 2*15+ 4] = QC[1]*gsfs[0][ 2*10+4] + WQ[1]*gsfs[1][ 2*10+4]
			   + 2*tmp.gsds[ 2*6+4];
	    gsgs[ 2*15+ 5] = QC[2]*gsfs[0][ 2*10+5] + WQ[2]*gsfs[1][ 2*10+5]
			   + 2*tmp.gsds[ 2*6+5] + rze2[4]*fsfs[1][2*10+5];
	    gsgs[ 2*15+ 6] = QC[0]*gsfs[0][ 2*10+6] + WQ[0]*gsfs[1][ 2*10+6]
			   +   tmp.gsds[ 2*6+1];
	    gsgs[ 2*15+ 7] = QC[1]*gsfs[0][ 2*10+7] + WQ[1]*gsfs[1][ 2*10+7]
			   +   tmp.gsds[ 2*6+2];
	    gsgs[ 2*15+ 8] = QC[2]*gsfs[0][ 2*10+8] + WQ[2]*gsfs[1][ 2*10+8]
			   +   tmp.gsds[ 2*6+0] + rze2[4]*fsfs[1][2*10+8];
	    gsgs[ 2*15+ 9] = QC[0]*gsfs[0][ 2*10+1] + WQ[0]*gsfs[1][ 2*10+1];
	    gsgs[ 2*15+10] = QC[1]*gsfs[0][ 2*10+2] + WQ[1]*gsfs[1][ 2*10+2];
	    gsgs[ 2*15+11] = QC[2]*gsfs[0][ 2*10+0] + WQ[2]*gsfs[1][ 2*10+0]
			   + rze2[4]*fsfs[1][2*10+0];
	    gsgs[ 2*15+12] = QC[0]*gsfs[0][ 2*10+9] + WQ[0]*gsfs[1][ 2*10+9]
			   +   tmp.gsds[ 2*6+4];
	    gsgs[ 2*15+13] = QC[1]*gsfs[0][ 2*10+9] + WQ[1]*gsfs[1][ 2*10+9]
			   +   tmp.gsds[ 2*6+5];
	    gsgs[ 2*15+14] = QC[2]*gsfs[0][ 2*10+9] + WQ[2]*gsfs[1][ 2*10+9]
			   +   tmp.gsds[ 2*6+3] + rze2[4]*fsfs[1][2*10+9];
	    gsgs[ 3*15+ 0] = QC[0]*gsfs[0][ 3*10+0] + WQ[0]*gsfs[1][ 3*10+0]
			   + 3*tmp.gsds[ 3*6+0] + rze2[3]*fsfs[1][3*10+0];
	    gsgs[ 3*15+ 1] = QC[1]*gsfs[0][ 3*10+1] + WQ[1]*gsfs[1][ 3*10+1]
			   + 3*tmp.gsds[ 3*6+1] + rze2[1]*fsfs[1][0*10+1];
	    gsgs[ 3*15+ 2] = QC[2]*gsfs[0][ 3*10+2] + WQ[2]*gsfs[1][ 3*10+2]
			   + 3*tmp.gsds[ 3*6+2];
	    gsgs[ 3*15+ 3] = QC[0]*gsfs[0][ 3*10+3] + WQ[0]*gsfs[1][ 3*10+3]
			   + 2*tmp.gsds[ 3*6+3] + rze2[3]*fsfs[1][3*10+3];
	    gsgs[ 3*15+ 4] = QC[1]*gsfs[0][ 3*10+4] + WQ[1]*gsfs[1][ 3*10+4]
			   + 2*tmp.gsds[ 3*6+4] + rze2[1]*fsfs[1][0*10+4];
	    gsgs[ 3*15+ 5] = QC[2]*gsfs[0][ 3*10+5] + WQ[2]*gsfs[1][ 3*10+5]
			   + 2*tmp.gsds[ 3*6+5];
	    gsgs[ 3*15+ 6] = QC[0]*gsfs[0][ 3*10+6] + WQ[0]*gsfs[1][ 3*10+6]
			   +   tmp.gsds[ 3*6+1] + rze2[3]*fsfs[1][3*10+6];
	    gsgs[ 3*15+ 7] = QC[1]*gsfs[0][ 3*10+7] + WQ[1]*gsfs[1][ 3*10+7]
			   +   tmp.gsds[ 3*6+2] + rze2[1]*fsfs[1][0*10+7];
	    gsgs[ 3*15+ 8] = QC[2]*gsfs[0][ 3*10+8] + WQ[2]*gsfs[1][ 3*10+8]
			   +   tmp.gsds[ 3*6+0];
	    gsgs[ 3*15+ 9] = QC[0]*gsfs[0][ 3*10+1] + WQ[0]*gsfs[1][ 3*10+1]
			   + rze2[3]*fsfs[1][3*10+1];
	    gsgs[ 3*15+10] = QC[1]*gsfs[0][ 3*10+2] + WQ[1]*gsfs[1][ 3*10+2]
			   + rze2[1]*fsfs[1][0*10+2];
	    gsgs[ 3*15+11] = QC[2]*gsfs[0][ 3*10+0] + WQ[2]*gsfs[1][ 3*10+0];
	    gsgs[ 3*15+12] = QC[0]*gsfs[0][ 3*10+9] + WQ[0]*gsfs[1][ 3*10+9]
			   +   tmp.gsds[ 3*6+4] + rze2[3]*fsfs[1][3*10+9];
	    gsgs[ 3*15+13] = QC[1]*gsfs[0][ 3*10+9] + WQ[1]*gsfs[1][ 3*10+9]
			   +   tmp.gsds[ 3*6+5] + rze2[1]*fsfs[1][0*10+9];
	    gsgs[ 3*15+14] = QC[2]*gsfs[0][ 3*10+9] + WQ[2]*gsfs[1][ 3*10+9]
			   +   tmp.gsds[ 3*6+3];
	    gsgs[ 4*15+ 0] = QC[0]*gsfs[0][ 4*10+0] + WQ[0]*gsfs[1][ 4*10+0]
			   + 3*tmp.gsds[ 4*6+0];
	    gsgs[ 4*15+ 1] = QC[1]*gsfs[0][ 4*10+1] + WQ[1]*gsfs[1][ 4*10+1]
			   + 3*tmp.gsds[ 4*6+1] + rze2[3]*fsfs[1][4*10+1];
	    gsgs[ 4*15+ 2] = QC[2]*gsfs[0][ 4*10+2] + WQ[2]*gsfs[1][ 4*10+2]
			   + 3*tmp.gsds[ 4*6+2] + rze2[1]*fsfs[1][1*10+2];
	    gsgs[ 4*15+ 3] = QC[0]*gsfs[0][ 4*10+3] + WQ[0]*gsfs[1][ 4*10+3]
			   + 2*tmp.gsds[ 4*6+3];
	    gsgs[ 4*15+ 4] = QC[1]*gsfs[0][ 4*10+4] + WQ[1]*gsfs[1][ 4*10+4]
			   + 2*tmp.gsds[ 4*6+4] + rze2[3]*fsfs[1][4*10+4];
	    gsgs[ 4*15+ 5] = QC[2]*gsfs[0][ 4*10+5] + WQ[2]*gsfs[1][ 4*10+5]
			   + 2*tmp.gsds[ 4*6+5] + rze2[1]*fsfs[1][1*10+5];
	    gsgs[ 4*15+ 6] = QC[0]*gsfs[0][ 4*10+6] + WQ[0]*gsfs[1][ 4*10+6]
			   +   tmp.gsds[ 4*6+1];
	    gsgs[ 4*15+ 7] = QC[1]*gsfs[0][ 4*10+7] + WQ[1]*gsfs[1][ 4*10+7]
			   +   tmp.gsds[ 4*6+2] + rze2[3]*fsfs[1][4*10+7];
	    gsgs[ 4*15+ 8] = QC[2]*gsfs[0][ 4*10+8] + WQ[2]*gsfs[1][ 4*10+8]
			   +   tmp.gsds[ 4*6+0] + rze2[1]*fsfs[1][1*10+8];
	    gsgs[ 4*15+ 9] = QC[0]*gsfs[0][ 4*10+1] + WQ[0]*gsfs[1][ 4*10+1];
	    gsgs[ 4*15+10] = QC[1]*gsfs[0][ 4*10+2] + WQ[1]*gsfs[1][ 4*10+2]
			   + rze2[3]*fsfs[1][4*10+2];
	    gsgs[ 4*15+11] = QC[2]*gsfs[0][ 4*10+0] + WQ[2]*gsfs[1][ 4*10+0]
			   + rze2[1]*fsfs[1][1*10+0];
	    gsgs[ 4*15+12] = QC[0]*gsfs[0][ 4*10+9] + WQ[0]*gsfs[1][ 4*10+9]
			   +   tmp.gsds[ 4*6+4];
	    gsgs[ 4*15+13] = QC[1]*gsfs[0][ 4*10+9] + WQ[1]*gsfs[1][ 4*10+9]
			   +   tmp.gsds[ 4*6+5] + rze2[3]*fsfs[1][4*10+9];
	    gsgs[ 4*15+14] = QC[2]*gsfs[0][ 4*10+9] + WQ[2]*gsfs[1][ 4*10+9]
			   +   tmp.gsds[ 4*6+3] + rze2[1]*fsfs[1][1*10+9];
	    gsgs[ 5*15+ 0] = QC[0]*gsfs[0][ 5*10+0] + WQ[0]*gsfs[1][ 5*10+0]
			   + 3*tmp.gsds[ 5*6+0] + rze2[1]*fsfs[1][2*10+0];
	    gsgs[ 5*15+ 1] = QC[1]*gsfs[0][ 5*10+1] + WQ[1]*gsfs[1][ 5*10+1]
			   + 3*tmp.gsds[ 5*6+1];
	    gsgs[ 5*15+ 2] = QC[2]*gsfs[0][ 5*10+2] + WQ[2]*gsfs[1][ 5*10+2]
			   + 3*tmp.gsds[ 5*6+2] + rze2[3]*fsfs[1][5*10+2];
	    gsgs[ 5*15+ 3] = QC[0]*gsfs[0][ 5*10+3] + WQ[0]*gsfs[1][ 5*10+3]
			   + 2*tmp.gsds[ 5*6+3] + rze2[1]*fsfs[1][2*10+3];
	    gsgs[ 5*15+ 4] = QC[1]*gsfs[0][ 5*10+4] + WQ[1]*gsfs[1][ 5*10+4]
			   + 2*tmp.gsds[ 5*6+4];
	    gsgs[ 5*15+ 5] = QC[2]*gsfs[0][ 5*10+5] + WQ[2]*gsfs[1][ 5*10+5]
			   + 2*tmp.gsds[ 5*6+5] + rze2[3]*fsfs[1][5*10+5];
	    gsgs[ 5*15+ 6] = QC[0]*gsfs[0][ 5*10+6] + WQ[0]*gsfs[1][ 5*10+6]
			   +   tmp.gsds[ 5*6+1] + rze2[1]*fsfs[1][2*10+6];
	    gsgs[ 5*15+ 7] = QC[1]*gsfs[0][ 5*10+7] + WQ[1]*gsfs[1][ 5*10+7]
			   +   tmp.gsds[ 5*6+2];
	    gsgs[ 5*15+ 8] = QC[2]*gsfs[0][ 5*10+8] + WQ[2]*gsfs[1][ 5*10+8]
			   +   tmp.gsds[ 5*6+0] + rze2[3]*fsfs[1][5*10+8];
	    gsgs[ 5*15+ 9] = QC[0]*gsfs[0][ 5*10+1] + WQ[0]*gsfs[1][ 5*10+1]
			   + rze2[1]*fsfs[1][2*10+1];
	    gsgs[ 5*15+10] = QC[1]*gsfs[0][ 5*10+2] + WQ[1]*gsfs[1][ 5*10+2];
	    gsgs[ 5*15+11] = QC[2]*gsfs[0][ 5*10+0] + WQ[2]*gsfs[1][ 5*10+0]
			   + rze2[3]*fsfs[1][5*10+0];
	    gsgs[ 5*15+12] = QC[0]*gsfs[0][ 5*10+9] + WQ[0]*gsfs[1][ 5*10+9]
			   +   tmp.gsds[ 5*6+4] + rze2[1]*fsfs[1][2*10+9];
	    gsgs[ 5*15+13] = QC[1]*gsfs[0][ 5*10+9] + WQ[1]*gsfs[1][ 5*10+9]
			   +   tmp.gsds[ 5*6+5];
	    gsgs[ 5*15+14] = QC[2]*gsfs[0][ 5*10+9] + WQ[2]*gsfs[1][ 5*10+9]
			   +   tmp.gsds[ 5*6+3] + rze2[3]*fsfs[1][5*10+9];
	    gsgs[ 6*15+ 0] = QC[0]*gsfs[0][ 6*10+0] + WQ[0]*gsfs[1][ 6*10+0]
			   + 3*tmp.gsds[ 6*6+0] + rze2[2]*fsfs[1][6*10+0];
	    gsgs[ 6*15+ 1] = QC[1]*gsfs[0][ 6*10+1] + WQ[1]*gsfs[1][ 6*10+1]
			   + 3*tmp.gsds[ 6*6+1] + rze2[2]*fsfs[1][3*10+1];
	    gsgs[ 6*15+ 2] = QC[2]*gsfs[0][ 6*10+2] + WQ[2]*gsfs[1][ 6*10+2]
			   + 3*tmp.gsds[ 6*6+2];
	    gsgs[ 6*15+ 3] = QC[0]*gsfs[0][ 6*10+3] + WQ[0]*gsfs[1][ 6*10+3]
			   + 2*tmp.gsds[ 6*6+3] + rze2[2]*fsfs[1][6*10+3];
	    gsgs[ 6*15+ 4] = QC[1]*gsfs[0][ 6*10+4] + WQ[1]*gsfs[1][ 6*10+4]
			   + 2*tmp.gsds[ 6*6+4] + rze2[2]*fsfs[1][3*10+4];
	    gsgs[ 6*15+ 5] = QC[2]*gsfs[0][ 6*10+5] + WQ[2]*gsfs[1][ 6*10+5]
			   + 2*tmp.gsds[ 6*6+5];
	    gsgs[ 6*15+ 6] = QC[0]*gsfs[0][ 6*10+6] + WQ[0]*gsfs[1][ 6*10+6]
			   +   tmp.gsds[ 6*6+1] + rze2[2]*fsfs[1][6*10+6];
	    gsgs[ 6*15+ 7] = QC[1]*gsfs[0][ 6*10+7] + WQ[1]*gsfs[1][ 6*10+7]
			   +   tmp.gsds[ 6*6+2] + rze2[2]*fsfs[1][3*10+7];
	    gsgs[ 6*15+ 8] = QC[2]*gsfs[0][ 6*10+8] + WQ[2]*gsfs[1][ 6*10+8]
			   +   tmp.gsds[ 6*6+0];
	    gsgs[ 6*15+ 9] = QC[0]*gsfs[0][ 6*10+1] + WQ[0]*gsfs[1][ 6*10+1]
			   + rze2[2]*fsfs[1][6*10+1];
	    gsgs[ 6*15+10] = QC[1]*gsfs[0][ 6*10+2] + WQ[1]*gsfs[1][ 6*10+2]
			   + rze2[2]*fsfs[1][3*10+2];
	    gsgs[ 6*15+11] = QC[2]*gsfs[0][ 6*10+0] + WQ[2]*gsfs[1][ 6*10+0];
	    gsgs[ 6*15+12] = QC[0]*gsfs[0][ 6*10+9] + WQ[0]*gsfs[1][ 6*10+9]
			   +   tmp.gsds[ 6*6+4] + rze2[2]*fsfs[1][6*10+9];
	    gsgs[ 6*15+13] = QC[1]*gsfs[0][ 6*10+9] + WQ[1]*gsfs[1][ 6*10+9]
			   +   tmp.gsds[ 6*6+5] + rze2[2]*fsfs[1][3*10+9];
	    gsgs[ 6*15+14] = QC[2]*gsfs[0][ 6*10+9] + WQ[2]*gsfs[1][ 6*10+9]
			   +   tmp.gsds[ 6*6+3];
	    gsgs[ 7*15+ 0] = QC[0]*gsfs[0][ 7*10+0] + WQ[0]*gsfs[1][ 7*10+0]
			   + 3*tmp.gsds[ 7*6+0];
	    gsgs[ 7*15+ 1] = QC[1]*gsfs[0][ 7*10+1] + WQ[1]*gsfs[1][ 7*10+1]
			   + 3*tmp.gsds[ 7*6+1] + rze2[2]*fsfs[1][7*10+1];
	    gsgs[ 7*15+ 2] = QC[2]*gsfs[0][ 7*10+2] + WQ[2]*gsfs[1][ 7*10+2]
			   + 3*tmp.gsds[ 7*6+2] + rze2[2]*fsfs[1][4*10+2];
	    gsgs[ 7*15+ 3] = QC[0]*gsfs[0][ 7*10+3] + WQ[0]*gsfs[1][ 7*10+3]
			   + 2*tmp.gsds[ 7*6+3];
	    gsgs[ 7*15+ 4] = QC[1]*gsfs[0][ 7*10+4] + WQ[1]*gsfs[1][ 7*10+4]
			   + 2*tmp.gsds[ 7*6+4] + rze2[2]*fsfs[1][7*10+4];
	    gsgs[ 7*15+ 5] = QC[2]*gsfs[0][ 7*10+5] + WQ[2]*gsfs[1][ 7*10+5]
			   + 2*tmp.gsds[ 7*6+5] + rze2[2]*fsfs[1][4*10+5];
	    gsgs[ 7*15+ 6] = QC[0]*gsfs[0][ 7*10+6] + WQ[0]*gsfs[1][ 7*10+6]
			   +   tmp.gsds[ 7*6+1];
	    gsgs[ 7*15+ 7] = QC[1]*gsfs[0][ 7*10+7] + WQ[1]*gsfs[1][ 7*10+7]
			   +   tmp.gsds[ 7*6+2] + rze2[2]*fsfs[1][7*10+7];
	    gsgs[ 7*15+ 8] = QC[2]*gsfs[0][ 7*10+8] + WQ[2]*gsfs[1][ 7*10+8]
			   +   tmp.gsds[ 7*6+0] + rze2[2]*fsfs[1][4*10+8];
	    gsgs[ 7*15+ 9] = QC[0]*gsfs[0][ 7*10+1] + WQ[0]*gsfs[1][ 7*10+1];
	    gsgs[ 7*15+10] = QC[1]*gsfs[0][ 7*10+2] + WQ[1]*gsfs[1][ 7*10+2]
			   + rze2[2]*fsfs[1][7*10+2];
	    gsgs[ 7*15+11] = QC[2]*gsfs[0][ 7*10+0] + WQ[2]*gsfs[1][ 7*10+0]
			   + rze2[2]*fsfs[1][4*10+0];
	    gsgs[ 7*15+12] = QC[0]*gsfs[0][ 7*10+9] + WQ[0]*gsfs[1][ 7*10+9]
			   +   tmp.gsds[ 7*6+4];
	    gsgs[ 7*15+13] = QC[1]*gsfs[0][ 7*10+9] + WQ[1]*gsfs[1][ 7*10+9]
			   +   tmp.gsds[ 7*6+5] + rze2[2]*fsfs[1][7*10+9];
	    gsgs[ 7*15+14] = QC[2]*gsfs[0][ 7*10+9] + WQ[2]*gsfs[1][ 7*10+9]
			   +   tmp.gsds[ 7*6+3] + rze2[2]*fsfs[1][4*10+9];
	    gsgs[ 8*15+ 0] = QC[0]*gsfs[0][ 8*10+0] + WQ[0]*gsfs[1][ 8*10+0]
			   + 3*tmp.gsds[ 8*6+0] + rze2[2]*fsfs[1][5*10+0];
	    gsgs[ 8*15+ 1] = QC[1]*gsfs[0][ 8*10+1] + WQ[1]*gsfs[1][ 8*10+1]
			   + 3*tmp.gsds[ 8*6+1];
	    gsgs[ 8*15+ 2] = QC[2]*gsfs[0][ 8*10+2] + WQ[2]*gsfs[1][ 8*10+2]
			   + 3*tmp.gsds[ 8*6+2] + rze2[2]*fsfs[1][8*10+2];
	    gsgs[ 8*15+ 3] = QC[0]*gsfs[0][ 8*10+3] + WQ[0]*gsfs[1][ 8*10+3]
			   + 2*tmp.gsds[ 8*6+3] + rze2[2]*fsfs[1][5*10+3];
	    gsgs[ 8*15+ 4] = QC[1]*gsfs[0][ 8*10+4] + WQ[1]*gsfs[1][ 8*10+4]
			   + 2*tmp.gsds[ 8*6+4];
	    gsgs[ 8*15+ 5] = QC[2]*gsfs[0][ 8*10+5] + WQ[2]*gsfs[1][ 8*10+5]
			   + 2*tmp.gsds[ 8*6+5] + rze2[2]*fsfs[1][8*10+5];
	    gsgs[ 8*15+ 6] = QC[0]*gsfs[0][ 8*10+6] + WQ[0]*gsfs[1][ 8*10+6]
			   +   tmp.gsds[ 8*6+1] + rze2[2]*fsfs[1][5*10+6];
	    gsgs[ 8*15+ 7] = QC[1]*gsfs[0][ 8*10+7] + WQ[1]*gsfs[1][ 8*10+7]
			   +   tmp.gsds[ 8*6+2];
	    gsgs[ 8*15+ 8] = QC[2]*gsfs[0][ 8*10+8] + WQ[2]*gsfs[1][ 8*10+8]
			   +   tmp.gsds[ 8*6+0] + rze2[2]*fsfs[1][8*10+8];
	    gsgs[ 8*15+ 9] = QC[0]*gsfs[0][ 8*10+1] + WQ[0]*gsfs[1][ 8*10+1]
			   + rze2[2]*fsfs[1][5*10+1];
	    gsgs[ 8*15+10] = QC[1]*gsfs[0][ 8*10+2] + WQ[1]*gsfs[1][ 8*10+2];
	    gsgs[ 8*15+11] = QC[2]*gsfs[0][ 8*10+0] + WQ[2]*gsfs[1][ 8*10+0]
			   + rze2[2]*fsfs[1][8*10+0];
	    gsgs[ 8*15+12] = QC[0]*gsfs[0][ 8*10+9] + WQ[0]*gsfs[1][ 8*10+9]
			   +   tmp.gsds[ 8*6+4] + rze2[2]*fsfs[1][5*10+9];
	    gsgs[ 8*15+13] = QC[1]*gsfs[0][ 8*10+9] + WQ[1]*gsfs[1][ 8*10+9]
			   +   tmp.gsds[ 8*6+5];
	    gsgs[ 8*15+14] = QC[2]*gsfs[0][ 8*10+9] + WQ[2]*gsfs[1][ 8*10+9]
			   +   tmp.gsds[ 8*6+3] + rze2[2]*fsfs[1][8*10+9];
	    gsgs[ 9*15+ 0] = QC[0]*gsfs[0][ 9*10+0] + WQ[0]*gsfs[1][ 9*10+0]
			   + 3*tmp.gsds[ 9*6+0] + rze2[1]*fsfs[1][1*10+0];
	    gsgs[ 9*15+ 1] = QC[1]*gsfs[0][ 9*10+1] + WQ[1]*gsfs[1][ 9*10+1]
			   + 3*tmp.gsds[ 9*6+1] + rze2[3]*fsfs[1][6*10+1];
	    gsgs[ 9*15+ 2] = QC[2]*gsfs[0][ 9*10+2] + WQ[2]*gsfs[1][ 9*10+2]
			   + 3*tmp.gsds[ 9*6+2];
	    gsgs[ 9*15+ 3] = QC[0]*gsfs[0][ 9*10+3] + WQ[0]*gsfs[1][ 9*10+3]
			   + 2*tmp.gsds[ 9*6+3] + rze2[1]*fsfs[1][1*10+3];
	    gsgs[ 9*15+ 4] = QC[1]*gsfs[0][ 9*10+4] + WQ[1]*gsfs[1][ 9*10+4]
			   + 2*tmp.gsds[ 9*6+4] + rze2[3]*fsfs[1][6*10+4];
	    gsgs[ 9*15+ 5] = QC[2]*gsfs[0][ 9*10+5] + WQ[2]*gsfs[1][ 9*10+5]
			   + 2*tmp.gsds[ 9*6+5];
	    gsgs[ 9*15+ 6] = QC[0]*gsfs[0][ 9*10+6] + WQ[0]*gsfs[1][ 9*10+6]
			   +   tmp.gsds[ 9*6+1] + rze2[1]*fsfs[1][1*10+6];
	    gsgs[ 9*15+ 7] = QC[1]*gsfs[0][ 9*10+7] + WQ[1]*gsfs[1][ 9*10+7]
			   +   tmp.gsds[ 9*6+2] + rze2[3]*fsfs[1][6*10+7];
	    gsgs[ 9*15+ 8] = QC[2]*gsfs[0][ 9*10+8] + WQ[2]*gsfs[1][ 9*10+8]
			   +   tmp.gsds[ 9*6+0];
	    gsgs[ 9*15+ 9] = QC[0]*gsfs[0][ 9*10+1] + WQ[0]*gsfs[1][ 9*10+1]
			   + rze2[1]*fsfs[1][1*10+1];
	    gsgs[ 9*15+10] = QC[1]*gsfs[0][ 9*10+2] + WQ[1]*gsfs[1][ 9*10+2]
			   + rze2[3]*fsfs[1][6*10+2];
	    gsgs[ 9*15+11] = QC[2]*gsfs[0][ 9*10+0] + WQ[2]*gsfs[1][ 9*10+0];
	    gsgs[ 9*15+12] = QC[0]*gsfs[0][ 9*10+9] + WQ[0]*gsfs[1][ 9*10+9]
			   +   tmp.gsds[ 9*6+4] + rze2[1]*fsfs[1][1*10+9];
	    gsgs[ 9*15+13] = QC[1]*gsfs[0][ 9*10+9] + WQ[1]*gsfs[1][ 9*10+9]
			   +   tmp.gsds[ 9*6+5] + rze2[3]*fsfs[1][6*10+9];
	    gsgs[ 9*15+14] = QC[2]*gsfs[0][ 9*10+9] + WQ[2]*gsfs[1][ 9*10+9]
			   +   tmp.gsds[ 9*6+3];
	    gsgs[10*15+ 0] = QC[0]*gsfs[0][10*10+0] + WQ[0]*gsfs[1][10*10+0]
			   + 3*tmp.gsds[10*6+0];
	    gsgs[10*15+ 1] = QC[1]*gsfs[0][10*10+1] + WQ[1]*gsfs[1][10*10+1]
			   + 3*tmp.gsds[10*6+1] + rze2[1]*fsfs[1][2*10+1];
	    gsgs[10*15+ 2] = QC[2]*gsfs[0][10*10+2] + WQ[2]*gsfs[1][10*10+2]
			   + 3*tmp.gsds[10*6+2] + rze2[3]*fsfs[1][7*10+2];
	    gsgs[10*15+ 3] = QC[0]*gsfs[0][10*10+3] + WQ[0]*gsfs[1][10*10+3]
			   + 2*tmp.gsds[10*6+3];
	    gsgs[10*15+ 4] = QC[1]*gsfs[0][10*10+4] + WQ[1]*gsfs[1][10*10+4]
			   + 2*tmp.gsds[10*6+4] + rze2[1]*fsfs[1][2*10+4];
	    gsgs[10*15+ 5] = QC[2]*gsfs[0][10*10+5] + WQ[2]*gsfs[1][10*10+5]
			   + 2*tmp.gsds[10*6+5] + rze2[3]*fsfs[1][7*10+5];
	    gsgs[10*15+ 6] = QC[0]*gsfs[0][10*10+6] + WQ[0]*gsfs[1][10*10+6]
			   +   tmp.gsds[10*6+1];
	    gsgs[10*15+ 7] = QC[1]*gsfs[0][10*10+7] + WQ[1]*gsfs[1][10*10+7]
			   +   tmp.gsds[10*6+2] + rze2[1]*fsfs[1][2*10+7];
	    gsgs[10*15+ 8] = QC[2]*gsfs[0][10*10+8] + WQ[2]*gsfs[1][10*10+8]
			   +   tmp.gsds[10*6+0] + rze2[3]*fsfs[1][7*10+8];
	    gsgs[10*15+ 9] = QC[0]*gsfs[0][10*10+1] + WQ[0]*gsfs[1][10*10+1];
	    gsgs[10*15+10] = QC[1]*gsfs[0][10*10+2] + WQ[1]*gsfs[1][10*10+2]
			   + rze2[1]*fsfs[1][2*10+2];
	    gsgs[10*15+11] = QC[2]*gsfs[0][10*10+0] + WQ[2]*gsfs[1][10*10+0]
			   + rze2[3]*fsfs[1][7*10+0];
	    gsgs[10*15+12] = QC[0]*gsfs[0][10*10+9] + WQ[0]*gsfs[1][10*10+9]
			   +   tmp.gsds[10*6+4];
	    gsgs[10*15+13] = QC[1]*gsfs[0][10*10+9] + WQ[1]*gsfs[1][10*10+9]
			   +   tmp.gsds[10*6+5] + rze2[1]*fsfs[1][2*10+9];
	    gsgs[10*15+14] = QC[2]*gsfs[0][10*10+9] + WQ[2]*gsfs[1][10*10+9]
			   +   tmp.gsds[10*6+3] + rze2[3]*fsfs[1][7*10+9];
	    gsgs[11*15+ 0] = QC[0]*gsfs[0][11*10+0] + WQ[0]*gsfs[1][11*10+0]
			   + 3*tmp.gsds[11*6+0] + rze2[3]*fsfs[1][8*10+0];
	    gsgs[11*15+ 1] = QC[1]*gsfs[0][11*10+1] + WQ[1]*gsfs[1][11*10+1]
			   + 3*tmp.gsds[11*6+1];
	    gsgs[11*15+ 2] = QC[2]*gsfs[0][11*10+2] + WQ[2]*gsfs[1][11*10+2]
			   + 3*tmp.gsds[11*6+2] + rze2[1]*fsfs[1][0*10+2];
	    gsgs[11*15+ 3] = QC[0]*gsfs[0][11*10+3] + WQ[0]*gsfs[1][11*10+3]
			   + 2*tmp.gsds[11*6+3] + rze2[3]*fsfs[1][8*10+3];
	    gsgs[11*15+ 4] = QC[1]*gsfs[0][11*10+4] + WQ[1]*gsfs[1][11*10+4]
			   + 2*tmp.gsds[11*6+4];
	    gsgs[11*15+ 5] = QC[2]*gsfs[0][11*10+5] + WQ[2]*gsfs[1][11*10+5]
			   + 2*tmp.gsds[11*6+5] + rze2[1]*fsfs[1][0*10+5];
	    gsgs[11*15+ 6] = QC[0]*gsfs[0][11*10+6] + WQ[0]*gsfs[1][11*10+6]
			   +   tmp.gsds[11*6+1] + rze2[3]*fsfs[1][8*10+6];
	    gsgs[11*15+ 7] = QC[1]*gsfs[0][11*10+7] + WQ[1]*gsfs[1][11*10+7]
			   +   tmp.gsds[11*6+2];
	    gsgs[11*15+ 8] = QC[2]*gsfs[0][11*10+8] + WQ[2]*gsfs[1][11*10+8]
			   +   tmp.gsds[11*6+0] + rze2[1]*fsfs[1][0*10+8];
	    gsgs[11*15+ 9] = QC[0]*gsfs[0][11*10+1] + WQ[0]*gsfs[1][11*10+1]
			   + rze2[3]*fsfs[1][8*10+1];
	    gsgs[11*15+10] = QC[1]*gsfs[0][11*10+2] + WQ[1]*gsfs[1][11*10+2];
	    gsgs[11*15+11] = QC[2]*gsfs[0][11*10+0] + WQ[2]*gsfs[1][11*10+0]
			   + rze2[1]*fsfs[1][0*10+0];
	    gsgs[11*15+12] = QC[0]*gsfs[0][11*10+9] + WQ[0]*gsfs[1][11*10+9]
			   +   tmp.gsds[11*6+4] + rze2[3]*fsfs[1][8*10+9];
	    gsgs[11*15+13] = QC[1]*gsfs[0][11*10+9] + WQ[1]*gsfs[1][11*10+9]
			   +   tmp.gsds[11*6+5];
	    gsgs[11*15+14] = QC[2]*gsfs[0][11*10+9] + WQ[2]*gsfs[1][11*10+9]
			   +   tmp.gsds[11*6+3] + rze2[1]*fsfs[1][0*10+9];
	    gsgs[12*15+ 0] = QC[0]*gsfs[0][12*10+0] + WQ[0]*gsfs[1][12*10+0]
			   + 3*tmp.gsds[12*6+0] + rze2[2]*fsfs[1][9*10+0];
	    gsgs[12*15+ 1] = QC[1]*gsfs[0][12*10+1] + WQ[1]*gsfs[1][12*10+1]
			   + 3*tmp.gsds[12*6+1] + rze2[1]*fsfs[1][8*10+1];
	    gsgs[12*15+ 2] = QC[2]*gsfs[0][12*10+2] + WQ[2]*gsfs[1][12*10+2]
			   + 3*tmp.gsds[12*6+2] + rze2[1]*fsfs[1][3*10+2];
	    gsgs[12*15+ 3] = QC[0]*gsfs[0][12*10+3] + WQ[0]*gsfs[1][12*10+3]
			   + 2*tmp.gsds[12*6+3] + rze2[2]*fsfs[1][9*10+3];
	    gsgs[12*15+ 4] = QC[1]*gsfs[0][12*10+4] + WQ[1]*gsfs[1][12*10+4]
			   + 2*tmp.gsds[12*6+4] + rze2[1]*fsfs[1][8*10+4];
	    gsgs[12*15+ 5] = QC[2]*gsfs[0][12*10+5] + WQ[2]*gsfs[1][12*10+5]
			   + 2*tmp.gsds[12*6+5] + rze2[1]*fsfs[1][3*10+5];
	    gsgs[12*15+ 6] = QC[0]*gsfs[0][12*10+6] + WQ[0]*gsfs[1][12*10+6]
			   +   tmp.gsds[12*6+1] + rze2[2]*fsfs[1][9*10+6];
	    gsgs[12*15+ 7] = QC[1]*gsfs[0][12*10+7] + WQ[1]*gsfs[1][12*10+7]
			   +   tmp.gsds[12*6+2] + rze2[1]*fsfs[1][8*10+7];
	    gsgs[12*15+ 8] = QC[2]*gsfs[0][12*10+8] + WQ[2]*gsfs[1][12*10+8]
			   +   tmp.gsds[12*6+0] + rze2[1]*fsfs[1][3*10+8];
	    gsgs[12*15+ 9] = QC[0]*gsfs[0][12*10+1] + WQ[0]*gsfs[1][12*10+1]
			   + rze2[2]*fsfs[1][9*10+1];
	    gsgs[12*15+10] = QC[1]*gsfs[0][12*10+2] + WQ[1]*gsfs[1][12*10+2]
			   + rze2[1]*fsfs[1][8*10+2];
	    gsgs[12*15+11] = QC[2]*gsfs[0][12*10+0] + WQ[2]*gsfs[1][12*10+0]
			   + rze2[1]*fsfs[1][3*10+0];
	    gsgs[12*15+12] = QC[0]*gsfs[0][12*10+9] + WQ[0]*gsfs[1][12*10+9]
			   +   tmp.gsds[12*6+4] + rze2[2]*fsfs[1][9*10+9];
	    gsgs[12*15+13] = QC[1]*gsfs[0][12*10+9] + WQ[1]*gsfs[1][12*10+9]
			   +   tmp.gsds[12*6+5] + rze2[1]*fsfs[1][8*10+9];
	    gsgs[12*15+14] = QC[2]*gsfs[0][12*10+9] + WQ[2]*gsfs[1][12*10+9]
			   +   tmp.gsds[12*6+3] + rze2[1]*fsfs[1][3*10+9];
	    gsgs[13*15+ 0] = QC[0]*gsfs[0][13*10+0] + WQ[0]*gsfs[1][13*10+0]
			   + 3*tmp.gsds[13*6+0] + rze2[1]*fsfs[1][4*10+0];
	    gsgs[13*15+ 1] = QC[1]*gsfs[0][13*10+1] + WQ[1]*gsfs[1][13*10+1]
			   + 3*tmp.gsds[13*6+1] + rze2[2]*fsfs[1][9*10+1];
	    gsgs[13*15+ 2] = QC[2]*gsfs[0][13*10+2] + WQ[2]*gsfs[1][13*10+2]
			   + 3*tmp.gsds[13*6+2] + rze2[1]*fsfs[1][6*10+2];
	    gsgs[13*15+ 3] = QC[0]*gsfs[0][13*10+3] + WQ[0]*gsfs[1][13*10+3]
			   + 2*tmp.gsds[13*6+3] + rze2[1]*fsfs[1][4*10+3];
	    gsgs[13*15+ 4] = QC[1]*gsfs[0][13*10+4] + WQ[1]*gsfs[1][13*10+4]
			   + 2*tmp.gsds[13*6+4] + rze2[2]*fsfs[1][9*10+4];
	    gsgs[13*15+ 5] = QC[2]*gsfs[0][13*10+5] + WQ[2]*gsfs[1][13*10+5]
			   + 2*tmp.gsds[13*6+5] + rze2[1]*fsfs[1][6*10+5];
	    gsgs[13*15+ 6] = QC[0]*gsfs[0][13*10+6] + WQ[0]*gsfs[1][13*10+6]
			   +   tmp.gsds[13*6+1] + rze2[1]*fsfs[1][4*10+6];
	    gsgs[13*15+ 7] = QC[1]*gsfs[0][13*10+7] + WQ[1]*gsfs[1][13*10+7]
			   +   tmp.gsds[13*6+2] + rze2[2]*fsfs[1][9*10+7];
	    gsgs[13*15+ 8] = QC[2]*gsfs[0][13*10+8] + WQ[2]*gsfs[1][13*10+8]
			   +   tmp.gsds[13*6+0] + rze2[1]*fsfs[1][6*10+8];
	    gsgs[13*15+ 9] = QC[0]*gsfs[0][13*10+1] + WQ[0]*gsfs[1][13*10+1]
			   + rze2[1]*fsfs[1][4*10+1];
	    gsgs[13*15+10] = QC[1]*gsfs[0][13*10+2] + WQ[1]*gsfs[1][13*10+2]
			   + rze2[2]*fsfs[1][9*10+2];
	    gsgs[13*15+11] = QC[2]*gsfs[0][13*10+0] + WQ[2]*gsfs[1][13*10+0]
			   + rze2[1]*fsfs[1][6*10+0];
	    gsgs[13*15+12] = QC[0]*gsfs[0][13*10+9] + WQ[0]*gsfs[1][13*10+9]
			   +   tmp.gsds[13*6+4] + rze2[1]*fsfs[1][4*10+9];
	    gsgs[13*15+13] = QC[1]*gsfs[0][13*10+9] + WQ[1]*gsfs[1][13*10+9]
			   +   tmp.gsds[13*6+5] + rze2[2]*fsfs[1][9*10+9];
	    gsgs[13*15+14] = QC[2]*gsfs[0][13*10+9] + WQ[2]*gsfs[1][13*10+9]
			   +   tmp.gsds[13*6+3] + rze2[1]*fsfs[1][6*10+9];
	    gsgs[14*15+ 0] = QC[0]*gsfs[0][14*10+0] + WQ[0]*gsfs[1][14*10+0]
			   + 3*tmp.gsds[14*6+0] + rze2[1]*fsfs[1][7*10+0];
	    gsgs[14*15+ 1] = QC[1]*gsfs[0][14*10+1] + WQ[1]*gsfs[1][14*10+1]
			   + 3*tmp.gsds[14*6+1] + rze2[1]*fsfs[1][5*10+1];
	    gsgs[14*15+ 2] = QC[2]*gsfs[0][14*10+2] + WQ[2]*gsfs[1][14*10+2]
			   + 3*tmp.gsds[14*6+2] + rze2[2]*fsfs[1][9*10+2];
	    gsgs[14*15+ 3] = QC[0]*gsfs[0][14*10+3] + WQ[0]*gsfs[1][14*10+3]
			   + 2*tmp.gsds[14*6+3] + rze2[1]*fsfs[1][7*10+3];
	    gsgs[14*15+ 4] = QC[1]*gsfs[0][14*10+4] + WQ[1]*gsfs[1][14*10+4]
			   + 2*tmp.gsds[14*6+4] + rze2[1]*fsfs[1][5*10+4];
	    gsgs[14*15+ 5] = QC[2]*gsfs[0][14*10+5] + WQ[2]*gsfs[1][14*10+5]
			   + 2*tmp.gsds[14*6+5] + rze2[2]*fsfs[1][9*10+5];
	    gsgs[14*15+ 6] = QC[0]*gsfs[0][14*10+6] + WQ[0]*gsfs[1][14*10+6]
			   +   tmp.gsds[14*6+1] + rze2[1]*fsfs[1][7*10+6];
	    gsgs[14*15+ 7] = QC[1]*gsfs[0][14*10+7] + WQ[1]*gsfs[1][14*10+7]
			   +   tmp.gsds[14*6+2] + rze2[1]*fsfs[1][5*10+7];
	    gsgs[14*15+ 8] = QC[2]*gsfs[0][14*10+8] + WQ[2]*gsfs[1][14*10+8]
			   +   tmp.gsds[14*6+0] + rze2[2]*fsfs[1][9*10+8];
	    gsgs[14*15+ 9] = QC[0]*gsfs[0][14*10+1] + WQ[0]*gsfs[1][14*10+1]
			   + rze2[1]*fsfs[1][7*10+1];
	    gsgs[14*15+10] = QC[1]*gsfs[0][14*10+2] + WQ[1]*gsfs[1][14*10+2]
			   + rze2[1]*fsfs[1][5*10+2];
	    gsgs[14*15+11] = QC[2]*gsfs[0][14*10+0] + WQ[2]*gsfs[1][14*10+0]
			   + rze2[2]*fsfs[1][9*10+0];
	    gsgs[14*15+12] = QC[0]*gsfs[0][14*10+9] + WQ[0]*gsfs[1][14*10+9]
			   +   tmp.gsds[14*6+4] + rze2[1]*fsfs[1][7*10+9];
	    gsgs[14*15+13] = QC[1]*gsfs[0][14*10+9] + WQ[1]*gsfs[1][14*10+9]
			   +   tmp.gsds[14*6+5] + rze2[1]*fsfs[1][5*10+9];
	    gsgs[14*15+14] = QC[2]*gsfs[0][14*10+9] + WQ[2]*gsfs[1][14*10+9]
			   +   tmp.gsds[14*6+3] + rze2[2]*fsfs[1][9*10+9];
			    
	    // ---- contraction for XSXS type integrals ----
	    for (int i=0; i<6*6; i++) DSDS[i] += dsds[0][i];
	    for (int i=0; i<6*10; i++) DSFS[i] += dsfs[0][i];
	    for (int i=0; i<6*15; i++) DSGS[i] += dsgs[i];
	    for (int i=0; i<10*6; i++) FSDS[i] += fsds[0][i];
	    for (int i=0; i<10*10; i++) FSFS[i] += fsfs[0][i];
	    for (int i=0; i<10*15; i++) FSGS[i] += fsgs[i];
	    for (int i=0; i<15*6; i++) GSDS[i] += gsds[0][i];
	    for (int i=0; i<15*10; i++) GSFS[i] += gsfs[0][i];
	    for (int i=0; i<15*15; i++) GSGS[i] += gsgs[i];
	}
    }
    // ========== HRR calculation ==========
    // (a,b+1|c,0) = (a+1,b|c,0) - BA(a,b|c,0)
    // (D,P|D,S)
    for (int c0=0; c0<6; c0++)
        DPDS[0*3*6+0*6+c0] = FSDS[0*6+c0] - BA[0]*DSDS[0*6+c0];
    for (int c0=0; c0<6; c0++)
        DPDS[0*3*6+1*6+c0] = FSDS[3*6+c0] - BA[1]*DSDS[0*6+c0];
    for (int c0=0; c0<6; c0++)
        DPDS[0*3*6+2*6+c0] = FSDS[8*6+c0] - BA[2]*DSDS[0*6+c0];
    for (int c0=0; c0<6; c0++)
        DPDS[1*3*6+0*6+c0] = FSDS[6*6+c0] - BA[0]*DSDS[1*6+c0];
    for (int c0=0; c0<6; c0++)
        DPDS[1*3*6+1*6+c0] = FSDS[1*6+c0] - BA[1]*DSDS[1*6+c0];
    for (int c0=0; c0<6; c0++)
        DPDS[1*3*6+2*6+c0] = FSDS[4*6+c0] - BA[2]*DSDS[1*6+c0];
    for (int c0=0; c0<6; c0++)
        DPDS[2*3*6+0*6+c0] = FSDS[5*6+c0] - BA[0]*DSDS[2*6+c0];
    for (int c0=0; c0<6; c0++)
        DPDS[2*3*6+1*6+c0] = FSDS[7*6+c0] - BA[1]*DSDS[2*6+c0];
    for (int c0=0; c0<6; c0++)
        DPDS[2*3*6+2*6+c0] = FSDS[2*6+c0] - BA[2]*DSDS[2*6+c0];
    for (int c0=0; c0<6; c0++)
        DPDS[3*3*6+0*6+c0] = FSDS[3*6+c0] - BA[0]*DSDS[3*6+c0];
    for (int c0=0; c0<6; c0++)
        DPDS[3*3*6+1*6+c0] = FSDS[6*6+c0] - BA[1]*DSDS[3*6+c0];
    for (int c0=0; c0<6; c0++)
        DPDS[3*3*6+2*6+c0] = FSDS[9*6+c0] - BA[2]*DSDS[3*6+c0];
    for (int c0=0; c0<6; c0++)
        DPDS[4*3*6+0*6+c0] = FSDS[9*6+c0] - BA[0]*DSDS[4*6+c0];
    for (int c0=0; c0<6; c0++)
        DPDS[4*3*6+1*6+c0] = FSDS[4*6+c0] - BA[1]*DSDS[4*6+c0];
    for (int c0=0; c0<6; c0++)
        DPDS[4*3*6+2*6+c0] = FSDS[7*6+c0] - BA[2]*DSDS[4*6+c0];
    for (int c0=0; c0<6; c0++)
        DPDS[5*3*6+0*6+c0] = FSDS[8*6+c0] - BA[0]*DSDS[5*6+c0];
    for (int c0=0; c0<6; c0++)
        DPDS[5*3*6+1*6+c0] = FSDS[9*6+c0] - BA[1]*DSDS[5*6+c0];
    for (int c0=0; c0<6; c0++)
        DPDS[5*3*6+2*6+c0] = FSDS[5*6+c0] - BA[2]*DSDS[5*6+c0];
    // (D,P|F,S)
    for (int c0=0; c0<10; c0++)
        DPFS[0*3*10+0*10+c0] = FSFS[0*10+c0] - BA[0]*DSFS[0*10+c0];
    for (int c0=0; c0<10; c0++)
        DPFS[0*3*10+1*10+c0] = FSFS[3*10+c0] - BA[1]*DSFS[0*10+c0];
    for (int c0=0; c0<10; c0++)
        DPFS[0*3*10+2*10+c0] = FSFS[8*10+c0] - BA[2]*DSFS[0*10+c0];
    for (int c0=0; c0<10; c0++)
        DPFS[1*3*10+0*10+c0] = FSFS[6*10+c0] - BA[0]*DSFS[1*10+c0];
    for (int c0=0; c0<10; c0++)
        DPFS[1*3*10+1*10+c0] = FSFS[1*10+c0] - BA[1]*DSFS[1*10+c0];
    for (int c0=0; c0<10; c0++)
        DPFS[1*3*10+2*10+c0] = FSFS[4*10+c0] - BA[2]*DSFS[1*10+c0];
    for (int c0=0; c0<10; c0++)
        DPFS[2*3*10+0*10+c0] = FSFS[5*10+c0] - BA[0]*DSFS[2*10+c0];
    for (int c0=0; c0<10; c0++)
        DPFS[2*3*10+1*10+c0] = FSFS[7*10+c0] - BA[1]*DSFS[2*10+c0];
    for (int c0=0; c0<10; c0++)
        DPFS[2*3*10+2*10+c0] = FSFS[2*10+c0] - BA[2]*DSFS[2*10+c0];
    for (int c0=0; c0<10; c0++)
        DPFS[3*3*10+0*10+c0] = FSFS[3*10+c0] - BA[0]*DSFS[3*10+c0];
    for (int c0=0; c0<10; c0++)
        DPFS[3*3*10+1*10+c0] = FSFS[6*10+c0] - BA[1]*DSFS[3*10+c0];
    for (int c0=0; c0<10; c0++)
        DPFS[3*3*10+2*10+c0] = FSFS[9*10+c0] - BA[2]*DSFS[3*10+c0];
    for (int c0=0; c0<10; c0++)
        DPFS[4*3*10+0*10+c0] = FSFS[9*10+c0] - BA[0]*DSFS[4*10+c0];
    for (int c0=0; c0<10; c0++)
        DPFS[4*3*10+1*10+c0] = FSFS[4*10+c0] - BA[1]*DSFS[4*10+c0];
    for (int c0=0; c0<10; c0++)
        DPFS[4*3*10+2*10+c0] = FSFS[7*10+c0] - BA[2]*DSFS[4*10+c0];
    for (int c0=0; c0<10; c0++)
        DPFS[5*3*10+0*10+c0] = FSFS[8*10+c0] - BA[0]*DSFS[5*10+c0];
    for (int c0=0; c0<10; c0++)
        DPFS[5*3*10+1*10+c0] = FSFS[9*10+c0] - BA[1]*DSFS[5*10+c0];
    for (int c0=0; c0<10; c0++)
        DPFS[5*3*10+2*10+c0] = FSFS[5*10+c0] - BA[2]*DSFS[5*10+c0];
    // (D,P|G,S)
    for (int c0=0; c0<15; c0++)
        DPGS[0*3*15+0*15+c0] = FSGS[0*15+c0] - BA[0]*DSGS[0*15+c0];
    for (int c0=0; c0<15; c0++)
        DPGS[0*3*15+1*15+c0] = FSGS[3*15+c0] - BA[1]*DSGS[0*15+c0];
    for (int c0=0; c0<15; c0++)
        DPGS[0*3*15+2*15+c0] = FSGS[8*15+c0] - BA[2]*DSGS[0*15+c0];
    for (int c0=0; c0<15; c0++)
        DPGS[1*3*15+0*15+c0] = FSGS[6*15+c0] - BA[0]*DSGS[1*15+c0];
    for (int c0=0; c0<15; c0++)
        DPGS[1*3*15+1*15+c0] = FSGS[1*15+c0] - BA[1]*DSGS[1*15+c0];
    for (int c0=0; c0<15; c0++)
        DPGS[1*3*15+2*15+c0] = FSGS[4*15+c0] - BA[2]*DSGS[1*15+c0];
    for (int c0=0; c0<15; c0++)
        DPGS[2*3*15+0*15+c0] = FSGS[5*15+c0] - BA[0]*DSGS[2*15+c0];
    for (int c0=0; c0<15; c0++)
        DPGS[2*3*15+1*15+c0] = FSGS[7*15+c0] - BA[1]*DSGS[2*15+c0];
    for (int c0=0; c0<15; c0++)
        DPGS[2*3*15+2*15+c0] = FSGS[2*15+c0] - BA[2]*DSGS[2*15+c0];
    for (int c0=0; c0<15; c0++)
        DPGS[3*3*15+0*15+c0] = FSGS[3*15+c0] - BA[0]*DSGS[3*15+c0];
    for (int c0=0; c0<15; c0++)
        DPGS[3*3*15+1*15+c0] = FSGS[6*15+c0] - BA[1]*DSGS[3*15+c0];
    for (int c0=0; c0<15; c0++)
        DPGS[3*3*15+2*15+c0] = FSGS[9*15+c0] - BA[2]*DSGS[3*15+c0];
    for (int c0=0; c0<15; c0++)
        DPGS[4*3*15+0*15+c0] = FSGS[9*15+c0] - BA[0]*DSGS[4*15+c0];
    for (int c0=0; c0<15; c0++)
        DPGS[4*3*15+1*15+c0] = FSGS[4*15+c0] - BA[1]*DSGS[4*15+c0];
    for (int c0=0; c0<15; c0++)
        DPGS[4*3*15+2*15+c0] = FSGS[7*15+c0] - BA[2]*DSGS[4*15+c0];
    for (int c0=0; c0<15; c0++)
        DPGS[5*3*15+0*15+c0] = FSGS[8*15+c0] - BA[0]*DSGS[5*15+c0];
    for (int c0=0; c0<15; c0++)
        DPGS[5*3*15+1*15+c0] = FSGS[9*15+c0] - BA[1]*DSGS[5*15+c0];
    for (int c0=0; c0<15; c0++)
        DPGS[5*3*15+2*15+c0] = FSGS[5*15+c0] - BA[2]*DSGS[5*15+c0];
    // (F,P|D,S)
    for (int c0=0; c0<6; c0++)
        FPDS[0*3*6+0*6+c0] = GSDS[ 0*6+c0] - BA[0]*FSDS[0*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[0*3*6+1*6+c0] = GSDS[ 3*6+c0] - BA[1]*FSDS[0*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[0*3*6+2*6+c0] = GSDS[11*6+c0] - BA[2]*FSDS[0*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[1*3*6+0*6+c0] = GSDS[ 9*6+c0] - BA[0]*FSDS[1*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[1*3*6+1*6+c0] = GSDS[ 1*6+c0] - BA[1]*FSDS[1*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[1*3*6+2*6+c0] = GSDS[ 4*6+c0] - BA[2]*FSDS[1*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[2*3*6+0*6+c0] = GSDS[ 5*6+c0] - BA[0]*FSDS[2*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[2*3*6+1*6+c0] = GSDS[10*6+c0] - BA[1]*FSDS[2*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[2*3*6+2*6+c0] = GSDS[ 2*6+c0] - BA[2]*FSDS[2*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[3*3*6+0*6+c0] = GSDS[ 3*6+c0] - BA[0]*FSDS[3*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[3*3*6+1*6+c0] = GSDS[ 6*6+c0] - BA[1]*FSDS[3*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[3*3*6+2*6+c0] = GSDS[12*6+c0] - BA[2]*FSDS[3*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[4*3*6+0*6+c0] = GSDS[13*6+c0] - BA[0]*FSDS[4*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[4*3*6+1*6+c0] = GSDS[ 4*6+c0] - BA[1]*FSDS[4*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[4*3*6+2*6+c0] = GSDS[ 7*6+c0] - BA[2]*FSDS[4*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[5*3*6+0*6+c0] = GSDS[ 8*6+c0] - BA[0]*FSDS[5*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[5*3*6+1*6+c0] = GSDS[14*6+c0] - BA[1]*FSDS[5*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[5*3*6+2*6+c0] = GSDS[ 5*6+c0] - BA[2]*FSDS[5*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[6*3*6+0*6+c0] = GSDS[ 6*6+c0] - BA[0]*FSDS[6*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[6*3*6+1*6+c0] = GSDS[ 9*6+c0] - BA[1]*FSDS[6*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[6*3*6+2*6+c0] = GSDS[13*6+c0] - BA[2]*FSDS[6*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[7*3*6+0*6+c0] = GSDS[14*6+c0] - BA[0]*FSDS[7*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[7*3*6+1*6+c0] = GSDS[ 7*6+c0] - BA[1]*FSDS[7*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[7*3*6+2*6+c0] = GSDS[10*6+c0] - BA[2]*FSDS[7*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[8*3*6+0*6+c0] = GSDS[11*6+c0] - BA[0]*FSDS[8*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[8*3*6+1*6+c0] = GSDS[12*6+c0] - BA[1]*FSDS[8*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[8*3*6+2*6+c0] = GSDS[ 8*6+c0] - BA[2]*FSDS[8*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[9*3*6+0*6+c0] = GSDS[12*6+c0] - BA[0]*FSDS[9*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[9*3*6+1*6+c0] = GSDS[13*6+c0] - BA[1]*FSDS[9*6+c0];
    for (int c0=0; c0<6; c0++)
        FPDS[9*3*6+2*6+c0] = GSDS[14*6+c0] - BA[2]*FSDS[9*6+c0];
    // (F,P|F,S)
    for (int c0=0; c0<10; c0++)
        FPFS[0*3*10+0*10+c0] = GSFS[ 0*10+c0] - BA[0]*FSFS[0*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[0*3*10+1*10+c0] = GSFS[ 3*10+c0] - BA[1]*FSFS[0*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[0*3*10+2*10+c0] = GSFS[11*10+c0] - BA[2]*FSFS[0*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[1*3*10+0*10+c0] = GSFS[ 9*10+c0] - BA[0]*FSFS[1*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[1*3*10+1*10+c0] = GSFS[ 1*10+c0] - BA[1]*FSFS[1*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[1*3*10+2*10+c0] = GSFS[ 4*10+c0] - BA[2]*FSFS[1*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[2*3*10+0*10+c0] = GSFS[ 5*10+c0] - BA[0]*FSFS[2*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[2*3*10+1*10+c0] = GSFS[10*10+c0] - BA[1]*FSFS[2*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[2*3*10+2*10+c0] = GSFS[ 2*10+c0] - BA[2]*FSFS[2*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[3*3*10+0*10+c0] = GSFS[ 3*10+c0] - BA[0]*FSFS[3*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[3*3*10+1*10+c0] = GSFS[ 6*10+c0] - BA[1]*FSFS[3*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[3*3*10+2*10+c0] = GSFS[12*10+c0] - BA[2]*FSFS[3*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[4*3*10+0*10+c0] = GSFS[13*10+c0] - BA[0]*FSFS[4*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[4*3*10+1*10+c0] = GSFS[ 4*10+c0] - BA[1]*FSFS[4*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[4*3*10+2*10+c0] = GSFS[ 7*10+c0] - BA[2]*FSFS[4*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[5*3*10+0*10+c0] = GSFS[ 8*10+c0] - BA[0]*FSFS[5*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[5*3*10+1*10+c0] = GSFS[14*10+c0] - BA[1]*FSFS[5*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[5*3*10+2*10+c0] = GSFS[ 5*10+c0] - BA[2]*FSFS[5*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[6*3*10+0*10+c0] = GSFS[ 6*10+c0] - BA[0]*FSFS[6*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[6*3*10+1*10+c0] = GSFS[ 9*10+c0] - BA[1]*FSFS[6*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[6*3*10+2*10+c0] = GSFS[13*10+c0] - BA[2]*FSFS[6*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[7*3*10+0*10+c0] = GSFS[14*10+c0] - BA[0]*FSFS[7*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[7*3*10+1*10+c0] = GSFS[ 7*10+c0] - BA[1]*FSFS[7*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[7*3*10+2*10+c0] = GSFS[10*10+c0] - BA[2]*FSFS[7*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[8*3*10+0*10+c0] = GSFS[11*10+c0] - BA[0]*FSFS[8*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[8*3*10+1*10+c0] = GSFS[12*10+c0] - BA[1]*FSFS[8*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[8*3*10+2*10+c0] = GSFS[ 8*10+c0] - BA[2]*FSFS[8*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[9*3*10+0*10+c0] = GSFS[12*10+c0] - BA[0]*FSFS[9*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[9*3*10+1*10+c0] = GSFS[13*10+c0] - BA[1]*FSFS[9*10+c0];
    for (int c0=0; c0<10; c0++)
        FPFS[9*3*10+2*10+c0] = GSFS[14*10+c0] - BA[2]*FSFS[9*10+c0];
    // (F,P|G,S)
    for (int c0=0; c0<15; c0++)
        FPGS[0*3*15+0*15+c0] = GSGS[ 0*15+c0] - BA[0]*FSGS[0*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[0*3*15+1*15+c0] = GSGS[ 3*15+c0] - BA[1]*FSGS[0*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[0*3*15+2*15+c0] = GSGS[11*15+c0] - BA[2]*FSGS[0*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[1*3*15+0*15+c0] = GSGS[ 9*15+c0] - BA[0]*FSGS[1*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[1*3*15+1*15+c0] = GSGS[ 1*15+c0] - BA[1]*FSGS[1*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[1*3*15+2*15+c0] = GSGS[ 4*15+c0] - BA[2]*FSGS[1*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[2*3*15+0*15+c0] = GSGS[ 5*15+c0] - BA[0]*FSGS[2*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[2*3*15+1*15+c0] = GSGS[10*15+c0] - BA[1]*FSGS[2*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[2*3*15+2*15+c0] = GSGS[ 2*15+c0] - BA[2]*FSGS[2*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[3*3*15+0*15+c0] = GSGS[ 3*15+c0] - BA[0]*FSGS[3*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[3*3*15+1*15+c0] = GSGS[ 6*15+c0] - BA[1]*FSGS[3*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[3*3*15+2*15+c0] = GSGS[12*15+c0] - BA[2]*FSGS[3*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[4*3*15+0*15+c0] = GSGS[13*15+c0] - BA[0]*FSGS[4*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[4*3*15+1*15+c0] = GSGS[ 4*15+c0] - BA[1]*FSGS[4*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[4*3*15+2*15+c0] = GSGS[ 7*15+c0] - BA[2]*FSGS[4*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[5*3*15+0*15+c0] = GSGS[ 8*15+c0] - BA[0]*FSGS[5*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[5*3*15+1*15+c0] = GSGS[14*15+c0] - BA[1]*FSGS[5*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[5*3*15+2*15+c0] = GSGS[ 5*15+c0] - BA[2]*FSGS[5*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[6*3*15+0*15+c0] = GSGS[ 6*15+c0] - BA[0]*FSGS[6*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[6*3*15+1*15+c0] = GSGS[ 9*15+c0] - BA[1]*FSGS[6*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[6*3*15+2*15+c0] = GSGS[13*15+c0] - BA[2]*FSGS[6*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[7*3*15+0*15+c0] = GSGS[14*15+c0] - BA[0]*FSGS[7*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[7*3*15+1*15+c0] = GSGS[ 7*15+c0] - BA[1]*FSGS[7*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[7*3*15+2*15+c0] = GSGS[10*15+c0] - BA[2]*FSGS[7*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[8*3*15+0*15+c0] = GSGS[11*15+c0] - BA[0]*FSGS[8*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[8*3*15+1*15+c0] = GSGS[12*15+c0] - BA[1]*FSGS[8*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[8*3*15+2*15+c0] = GSGS[ 8*15+c0] - BA[2]*FSGS[8*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[9*3*15+0*15+c0] = GSGS[12*15+c0] - BA[0]*FSGS[9*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[9*3*15+1*15+c0] = GSGS[13*15+c0] - BA[1]*FSGS[9*15+c0];
    for (int c0=0; c0<15; c0++)
        FPGS[9*3*15+2*15+c0] = GSGS[14*15+c0] - BA[2]*FSGS[9*15+c0];
    // (D,D|D,S)
    for (int c0=0; c0<6; c0++)
        DDDS[0*6*6+0*6+c0] = FPDS[0*3*6+0*6+c0] - BA[0]*DPDS[0*3*6+0*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[0*6*6+1*6+c0] = FPDS[3*3*6+1*6+c0] - BA[1]*DPDS[0*3*6+1*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[0*6*6+2*6+c0] = FPDS[8*3*6+2*6+c0] - BA[2]*DPDS[0*3*6+2*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[0*6*6+3*6+c0] = FPDS[0*3*6+1*6+c0] - BA[0]*DPDS[0*3*6+1*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[0*6*6+4*6+c0] = FPDS[3*3*6+2*6+c0] - BA[1]*DPDS[0*3*6+2*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[0*6*6+5*6+c0] = FPDS[8*3*6+0*6+c0] - BA[2]*DPDS[0*3*6+0*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[1*6*6+0*6+c0] = FPDS[6*3*6+0*6+c0] - BA[0]*DPDS[1*3*6+0*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[1*6*6+1*6+c0] = FPDS[1*3*6+1*6+c0] - BA[1]*DPDS[1*3*6+1*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[1*6*6+2*6+c0] = FPDS[4*3*6+2*6+c0] - BA[2]*DPDS[1*3*6+2*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[1*6*6+3*6+c0] = FPDS[6*3*6+1*6+c0] - BA[0]*DPDS[1*3*6+1*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[1*6*6+4*6+c0] = FPDS[1*3*6+2*6+c0] - BA[1]*DPDS[1*3*6+2*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[1*6*6+5*6+c0] = FPDS[4*3*6+0*6+c0] - BA[2]*DPDS[1*3*6+0*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[2*6*6+0*6+c0] = FPDS[5*3*6+0*6+c0] - BA[0]*DPDS[2*3*6+0*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[2*6*6+1*6+c0] = FPDS[7*3*6+1*6+c0] - BA[1]*DPDS[2*3*6+1*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[2*6*6+2*6+c0] = FPDS[2*3*6+2*6+c0] - BA[2]*DPDS[2*3*6+2*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[2*6*6+3*6+c0] = FPDS[5*3*6+1*6+c0] - BA[0]*DPDS[2*3*6+1*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[2*6*6+4*6+c0] = FPDS[7*3*6+2*6+c0] - BA[1]*DPDS[2*3*6+2*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[2*6*6+5*6+c0] = FPDS[2*3*6+0*6+c0] - BA[2]*DPDS[2*3*6+0*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[3*6*6+0*6+c0] = FPDS[3*3*6+0*6+c0] - BA[0]*DPDS[3*3*6+0*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[3*6*6+1*6+c0] = FPDS[6*3*6+1*6+c0] - BA[1]*DPDS[3*3*6+1*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[3*6*6+2*6+c0] = FPDS[9*3*6+2*6+c0] - BA[2]*DPDS[3*3*6+2*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[3*6*6+3*6+c0] = FPDS[3*3*6+1*6+c0] - BA[0]*DPDS[3*3*6+1*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[3*6*6+4*6+c0] = FPDS[6*3*6+2*6+c0] - BA[1]*DPDS[3*3*6+2*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[3*6*6+5*6+c0] = FPDS[9*3*6+0*6+c0] - BA[2]*DPDS[3*3*6+0*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[4*6*6+0*6+c0] = FPDS[9*3*6+0*6+c0] - BA[0]*DPDS[4*3*6+0*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[4*6*6+1*6+c0] = FPDS[4*3*6+1*6+c0] - BA[1]*DPDS[4*3*6+1*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[4*6*6+2*6+c0] = FPDS[7*3*6+2*6+c0] - BA[2]*DPDS[4*3*6+2*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[4*6*6+3*6+c0] = FPDS[9*3*6+1*6+c0] - BA[0]*DPDS[4*3*6+1*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[4*6*6+4*6+c0] = FPDS[4*3*6+2*6+c0] - BA[1]*DPDS[4*3*6+2*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[4*6*6+5*6+c0] = FPDS[7*3*6+0*6+c0] - BA[2]*DPDS[4*3*6+0*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[5*6*6+0*6+c0] = FPDS[8*3*6+0*6+c0] - BA[0]*DPDS[5*3*6+0*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[5*6*6+1*6+c0] = FPDS[9*3*6+1*6+c0] - BA[1]*DPDS[5*3*6+1*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[5*6*6+2*6+c0] = FPDS[5*3*6+2*6+c0] - BA[2]*DPDS[5*3*6+2*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[5*6*6+3*6+c0] = FPDS[8*3*6+1*6+c0] - BA[0]*DPDS[5*3*6+1*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[5*6*6+4*6+c0] = FPDS[9*3*6+2*6+c0] - BA[1]*DPDS[5*3*6+2*6+c0];
    for (int c0=0; c0<6; c0++)
        DDDS[5*6*6+5*6+c0] = FPDS[5*3*6+0*6+c0] - BA[2]*DPDS[5*3*6+0*6+c0];
    // (D,D|F,S)
    for (int c0=0; c0<10; c0++)
        DDFS[0*6*10+0*10+c0] = FPFS[0*3*10+0*10+c0] - BA[0]*DPFS[0*3*10+0*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[0*6*10+1*10+c0] = FPFS[3*3*10+1*10+c0] - BA[1]*DPFS[0*3*10+1*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[0*6*10+2*10+c0] = FPFS[8*3*10+2*10+c0] - BA[2]*DPFS[0*3*10+2*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[0*6*10+3*10+c0] = FPFS[0*3*10+1*10+c0] - BA[0]*DPFS[0*3*10+1*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[0*6*10+4*10+c0] = FPFS[3*3*10+2*10+c0] - BA[1]*DPFS[0*3*10+2*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[0*6*10+5*10+c0] = FPFS[8*3*10+0*10+c0] - BA[2]*DPFS[0*3*10+0*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[1*6*10+0*10+c0] = FPFS[6*3*10+0*10+c0] - BA[0]*DPFS[1*3*10+0*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[1*6*10+1*10+c0] = FPFS[1*3*10+1*10+c0] - BA[1]*DPFS[1*3*10+1*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[1*6*10+2*10+c0] = FPFS[4*3*10+2*10+c0] - BA[2]*DPFS[1*3*10+2*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[1*6*10+3*10+c0] = FPFS[6*3*10+1*10+c0] - BA[0]*DPFS[1*3*10+1*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[1*6*10+4*10+c0] = FPFS[1*3*10+2*10+c0] - BA[1]*DPFS[1*3*10+2*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[1*6*10+5*10+c0] = FPFS[4*3*10+0*10+c0] - BA[2]*DPFS[1*3*10+0*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[2*6*10+0*10+c0] = FPFS[5*3*10+0*10+c0] - BA[0]*DPFS[2*3*10+0*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[2*6*10+1*10+c0] = FPFS[7*3*10+1*10+c0] - BA[1]*DPFS[2*3*10+1*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[2*6*10+2*10+c0] = FPFS[2*3*10+2*10+c0] - BA[2]*DPFS[2*3*10+2*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[2*6*10+3*10+c0] = FPFS[5*3*10+1*10+c0] - BA[0]*DPFS[2*3*10+1*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[2*6*10+4*10+c0] = FPFS[7*3*10+2*10+c0] - BA[1]*DPFS[2*3*10+2*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[2*6*10+5*10+c0] = FPFS[2*3*10+0*10+c0] - BA[2]*DPFS[2*3*10+0*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[3*6*10+0*10+c0] = FPFS[3*3*10+0*10+c0] - BA[0]*DPFS[3*3*10+0*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[3*6*10+1*10+c0] = FPFS[6*3*10+1*10+c0] - BA[1]*DPFS[3*3*10+1*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[3*6*10+2*10+c0] = FPFS[9*3*10+2*10+c0] - BA[2]*DPFS[3*3*10+2*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[3*6*10+3*10+c0] = FPFS[3*3*10+1*10+c0] - BA[0]*DPFS[3*3*10+1*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[3*6*10+4*10+c0] = FPFS[6*3*10+2*10+c0] - BA[1]*DPFS[3*3*10+2*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[3*6*10+5*10+c0] = FPFS[9*3*10+0*10+c0] - BA[2]*DPFS[3*3*10+0*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[4*6*10+0*10+c0] = FPFS[9*3*10+0*10+c0] - BA[0]*DPFS[4*3*10+0*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[4*6*10+1*10+c0] = FPFS[4*3*10+1*10+c0] - BA[1]*DPFS[4*3*10+1*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[4*6*10+2*10+c0] = FPFS[7*3*10+2*10+c0] - BA[2]*DPFS[4*3*10+2*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[4*6*10+3*10+c0] = FPFS[9*3*10+1*10+c0] - BA[0]*DPFS[4*3*10+1*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[4*6*10+4*10+c0] = FPFS[4*3*10+2*10+c0] - BA[1]*DPFS[4*3*10+2*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[4*6*10+5*10+c0] = FPFS[7*3*10+0*10+c0] - BA[2]*DPFS[4*3*10+0*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[5*6*10+0*10+c0] = FPFS[8*3*10+0*10+c0] - BA[0]*DPFS[5*3*10+0*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[5*6*10+1*10+c0] = FPFS[9*3*10+1*10+c0] - BA[1]*DPFS[5*3*10+1*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[5*6*10+2*10+c0] = FPFS[5*3*10+2*10+c0] - BA[2]*DPFS[5*3*10+2*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[5*6*10+3*10+c0] = FPFS[8*3*10+1*10+c0] - BA[0]*DPFS[5*3*10+1*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[5*6*10+4*10+c0] = FPFS[9*3*10+2*10+c0] - BA[1]*DPFS[5*3*10+2*10+c0];
    for (int c0=0; c0<10; c0++)
        DDFS[5*6*10+5*10+c0] = FPFS[5*3*10+0*10+c0] - BA[2]*DPFS[5*3*10+0*10+c0];
    // (D,D|G,S)
    for (int c0=0; c0<15; c0++)
        DDGS[0*6*15+0*15+c0] = FPGS[0*3*15+0*15+c0] - BA[0]*DPGS[0*3*15+0*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[0*6*15+1*15+c0] = FPGS[3*3*15+1*15+c0] - BA[1]*DPGS[0*3*15+1*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[0*6*15+2*15+c0] = FPGS[8*3*15+2*15+c0] - BA[2]*DPGS[0*3*15+2*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[0*6*15+3*15+c0] = FPGS[0*3*15+1*15+c0] - BA[0]*DPGS[0*3*15+1*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[0*6*15+4*15+c0] = FPGS[3*3*15+2*15+c0] - BA[1]*DPGS[0*3*15+2*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[0*6*15+5*15+c0] = FPGS[8*3*15+0*15+c0] - BA[2]*DPGS[0*3*15+0*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[1*6*15+0*15+c0] = FPGS[6*3*15+0*15+c0] - BA[0]*DPGS[1*3*15+0*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[1*6*15+1*15+c0] = FPGS[1*3*15+1*15+c0] - BA[1]*DPGS[1*3*15+1*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[1*6*15+2*15+c0] = FPGS[4*3*15+2*15+c0] - BA[2]*DPGS[1*3*15+2*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[1*6*15+3*15+c0] = FPGS[6*3*15+1*15+c0] - BA[0]*DPGS[1*3*15+1*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[1*6*15+4*15+c0] = FPGS[1*3*15+2*15+c0] - BA[1]*DPGS[1*3*15+2*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[1*6*15+5*15+c0] = FPGS[4*3*15+0*15+c0] - BA[2]*DPGS[1*3*15+0*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[2*6*15+0*15+c0] = FPGS[5*3*15+0*15+c0] - BA[0]*DPGS[2*3*15+0*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[2*6*15+1*15+c0] = FPGS[7*3*15+1*15+c0] - BA[1]*DPGS[2*3*15+1*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[2*6*15+2*15+c0] = FPGS[2*3*15+2*15+c0] - BA[2]*DPGS[2*3*15+2*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[2*6*15+3*15+c0] = FPGS[5*3*15+1*15+c0] - BA[0]*DPGS[2*3*15+1*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[2*6*15+4*15+c0] = FPGS[7*3*15+2*15+c0] - BA[1]*DPGS[2*3*15+2*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[2*6*15+5*15+c0] = FPGS[2*3*15+0*15+c0] - BA[2]*DPGS[2*3*15+0*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[3*6*15+0*15+c0] = FPGS[3*3*15+0*15+c0] - BA[0]*DPGS[3*3*15+0*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[3*6*15+1*15+c0] = FPGS[6*3*15+1*15+c0] - BA[1]*DPGS[3*3*15+1*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[3*6*15+2*15+c0] = FPGS[9*3*15+2*15+c0] - BA[2]*DPGS[3*3*15+2*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[3*6*15+3*15+c0] = FPGS[3*3*15+1*15+c0] - BA[0]*DPGS[3*3*15+1*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[3*6*15+4*15+c0] = FPGS[6*3*15+2*15+c0] - BA[1]*DPGS[3*3*15+2*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[3*6*15+5*15+c0] = FPGS[9*3*15+0*15+c0] - BA[2]*DPGS[3*3*15+0*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[4*6*15+0*15+c0] = FPGS[9*3*15+0*15+c0] - BA[0]*DPGS[4*3*15+0*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[4*6*15+1*15+c0] = FPGS[4*3*15+1*15+c0] - BA[1]*DPGS[4*3*15+1*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[4*6*15+2*15+c0] = FPGS[7*3*15+2*15+c0] - BA[2]*DPGS[4*3*15+2*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[4*6*15+3*15+c0] = FPGS[9*3*15+1*15+c0] - BA[0]*DPGS[4*3*15+1*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[4*6*15+4*15+c0] = FPGS[4*3*15+2*15+c0] - BA[1]*DPGS[4*3*15+2*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[4*6*15+5*15+c0] = FPGS[7*3*15+0*15+c0] - BA[2]*DPGS[4*3*15+0*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[5*6*15+0*15+c0] = FPGS[8*3*15+0*15+c0] - BA[0]*DPGS[5*3*15+0*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[5*6*15+1*15+c0] = FPGS[9*3*15+1*15+c0] - BA[1]*DPGS[5*3*15+1*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[5*6*15+2*15+c0] = FPGS[5*3*15+2*15+c0] - BA[2]*DPGS[5*3*15+2*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[5*6*15+3*15+c0] = FPGS[8*3*15+1*15+c0] - BA[0]*DPGS[5*3*15+1*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[5*6*15+4*15+c0] = FPGS[9*3*15+2*15+c0] - BA[1]*DPGS[5*3*15+2*15+c0];
    for (int c0=0; c0<15; c0++)
        DDGS[5*6*15+5*15+c0] = FPGS[5*3*15+0*15+c0] - BA[2]*DPGS[5*3*15+0*15+c0];
    // (a,b|c,d+1) = (a,b|c+1,d) - DC(a,b|c,d)
    // (D,D|D,P)
    for (int ab=0, ab01=0, ab10=0, ab00=0;
            ab<6*6; ab++, ab01 += 6*3, ab10 += 10, ab00 += 6) {
        DDDP[ab01+0*3+0] = DDFS[ab10+0] - BA[0]*DDDS[ab00+0];
        DDDP[ab01+0*3+1] = DDFS[ab10+3] - BA[1]*DDDS[ab00+0];
        DDDP[ab01+0*3+2] = DDFS[ab10+8] - BA[2]*DDDS[ab00+0];
        DDDP[ab01+1*3+0] = DDFS[ab10+6] - BA[0]*DDDS[ab00+1];
        DDDP[ab01+1*3+1] = DDFS[ab10+1] - BA[1]*DDDS[ab00+1];
        DDDP[ab01+1*3+2] = DDFS[ab10+4] - BA[2]*DDDS[ab00+1];
        DDDP[ab01+2*3+0] = DDFS[ab10+5] - BA[0]*DDDS[ab00+2];
        DDDP[ab01+2*3+1] = DDFS[ab10+7] - BA[1]*DDDS[ab00+2];
        DDDP[ab01+2*3+2] = DDFS[ab10+2] - BA[2]*DDDS[ab00+2];
        DDDP[ab01+3*3+0] = DDFS[ab10+3] - BA[0]*DDDS[ab00+3];
        DDDP[ab01+3*3+1] = DDFS[ab10+6] - BA[1]*DDDS[ab00+3];
        DDDP[ab01+3*3+2] = DDFS[ab10+9] - BA[2]*DDDS[ab00+3];
        DDDP[ab01+4*3+0] = DDFS[ab10+9] - BA[0]*DDDS[ab00+4];
        DDDP[ab01+4*3+1] = DDFS[ab10+4] - BA[1]*DDDS[ab00+4];
        DDDP[ab01+4*3+2] = DDFS[ab10+7] - BA[2]*DDDS[ab00+4];
        DDDP[ab01+5*3+0] = DDFS[ab10+8] - BA[0]*DDDS[ab00+5];
        DDDP[ab01+5*3+1] = DDFS[ab10+9] - BA[1]*DDDS[ab00+5];
        DDDP[ab01+5*3+2] = DDFS[ab10+5] - BA[2]*DDDS[ab00+5];
    }
    // (D,D|F,P)
    for (int ab=0, ab01=0, ab10=0, ab00=0;
            ab<6*6; ab++, ab01 += 10*3, ab10 += 15, ab00 += 10) {
        DDFP[ab01+0*3+0] = DDGS[ab10+ 0] - BA[0]*DDFS[ab00+0];
        DDFP[ab01+0*3+1] = DDGS[ab10+ 3] - BA[1]*DDFS[ab00+0];
        DDFP[ab01+0*3+2] = DDGS[ab10+11] - BA[2]*DDFS[ab00+0];
        DDFP[ab01+1*3+0] = DDGS[ab10+ 9] - BA[0]*DDFS[ab00+1];
        DDFP[ab01+1*3+1] = DDGS[ab10+ 1] - BA[1]*DDFS[ab00+1];
        DDFP[ab01+1*3+2] = DDGS[ab10+ 4] - BA[2]*DDFS[ab00+1];
        DDFP[ab01+2*3+0] = DDGS[ab10+ 5] - BA[0]*DDFS[ab00+2];
        DDFP[ab01+2*3+1] = DDGS[ab10+10] - BA[1]*DDFS[ab00+2];
        DDFP[ab01+2*3+2] = DDGS[ab10+ 2] - BA[2]*DDFS[ab00+2];
        DDFP[ab01+3*3+0] = DDGS[ab10+ 3] - BA[0]*DDFS[ab00+3];
        DDFP[ab01+3*3+1] = DDGS[ab10+ 6] - BA[1]*DDFS[ab00+3];
        DDFP[ab01+3*3+2] = DDGS[ab10+12] - BA[2]*DDFS[ab00+3];
        DDFP[ab01+4*3+0] = DDGS[ab10+13] - BA[0]*DDFS[ab00+4];
        DDFP[ab01+4*3+1] = DDGS[ab10+ 4] - BA[1]*DDFS[ab00+4];
        DDFP[ab01+4*3+2] = DDGS[ab10+ 7] - BA[2]*DDFS[ab00+4];
        DDFP[ab01+5*3+0] = DDGS[ab10+ 8] - BA[0]*DDFS[ab00+5];
        DDFP[ab01+5*3+1] = DDGS[ab10+14] - BA[1]*DDFS[ab00+5];
        DDFP[ab01+5*3+2] = DDGS[ab10+ 5] - BA[2]*DDFS[ab00+5];
        DDFP[ab01+6*3+0] = DDGS[ab10+ 6] - BA[0]*DDFS[ab00+6];
        DDFP[ab01+6*3+1] = DDGS[ab10+ 9] - BA[1]*DDFS[ab00+6];
        DDFP[ab01+6*3+2] = DDGS[ab10+13] - BA[2]*DDFS[ab00+6];
        DDFP[ab01+7*3+0] = DDGS[ab10+14] - BA[0]*DDFS[ab00+7];
        DDFP[ab01+7*3+1] = DDGS[ab10+ 7] - BA[1]*DDFS[ab00+7];
        DDFP[ab01+7*3+2] = DDGS[ab10+10] - BA[2]*DDFS[ab00+7];
        DDFP[ab01+8*3+0] = DDGS[ab10+11] - BA[0]*DDFS[ab00+8];
        DDFP[ab01+8*3+1] = DDGS[ab10+12] - BA[1]*DDFS[ab00+8];
        DDFP[ab01+8*3+2] = DDGS[ab10+ 8] - BA[2]*DDFS[ab00+8];
        DDFP[ab01+9*3+0] = DDGS[ab10+12] - BA[0]*DDFS[ab00+9];
        DDFP[ab01+9*3+1] = DDGS[ab10+13] - BA[1]*DDFS[ab00+9];
        DDFP[ab01+9*3+2] = DDGS[ab10+14] - BA[2]*DDFS[ab00+9];
    }
    // (D,D|D,D)
    for (int ab=0, ab01=0, ab10=0, ab00=0;
            ab<6*6; ab++, ab01 += 6*6, ab10 += 10*3, ab00 += 6*3) {
        DDDD[ab01+0*6+0] = DDFP[ab10+0*3+0] - BA[0]*DDDP[ab00+0*3+0];
        DDDD[ab01+0*6+1] = DDFP[ab10+3*3+1] - BA[1]*DDDP[ab00+0*3+1];
        DDDD[ab01+0*6+2] = DDFP[ab10+8*3+2] - BA[2]*DDDP[ab00+0*3+2];
        DDDD[ab01+0*6+3] = DDFP[ab10+0*3+1] - BA[0]*DDDP[ab00+0*3+1];
        DDDD[ab01+0*6+4] = DDFP[ab10+3*3+2] - BA[1]*DDDP[ab00+0*3+2];
        DDDD[ab01+0*6+5] = DDFP[ab10+8*3+0] - BA[2]*DDDP[ab00+0*3+0];
        DDDD[ab01+1*6+0] = DDFP[ab10+6*3+0] - BA[0]*DDDP[ab00+1*3+0];
        DDDD[ab01+1*6+1] = DDFP[ab10+1*3+1] - BA[1]*DDDP[ab00+1*3+1];
        DDDD[ab01+1*6+2] = DDFP[ab10+4*3+2] - BA[2]*DDDP[ab00+1*3+2];
        DDDD[ab01+1*6+3] = DDFP[ab10+6*3+1] - BA[0]*DDDP[ab00+1*3+1];
        DDDD[ab01+1*6+4] = DDFP[ab10+1*3+2] - BA[1]*DDDP[ab00+1*3+2];
        DDDD[ab01+1*6+5] = DDFP[ab10+4*3+0] - BA[2]*DDDP[ab00+1*3+0];
        DDDD[ab01+2*6+0] = DDFP[ab10+5*3+0] - BA[0]*DDDP[ab00+2*3+0];
        DDDD[ab01+2*6+1] = DDFP[ab10+7*3+1] - BA[1]*DDDP[ab00+2*3+1];
        DDDD[ab01+2*6+2] = DDFP[ab10+2*3+2] - BA[2]*DDDP[ab00+2*3+2];
        DDDD[ab01+2*6+3] = DDFP[ab10+5*3+1] - BA[0]*DDDP[ab00+2*3+1];
        DDDD[ab01+2*6+4] = DDFP[ab10+7*3+2] - BA[1]*DDDP[ab00+2*3+2];
        DDDD[ab01+2*6+5] = DDFP[ab10+2*3+0] - BA[2]*DDDP[ab00+2*3+0];
        DDDD[ab01+3*6+0] = DDFP[ab10+3*3+0] - BA[0]*DDDP[ab00+3*3+0];
        DDDD[ab01+3*6+1] = DDFP[ab10+6*3+1] - BA[1]*DDDP[ab00+3*3+1];
        DDDD[ab01+3*6+2] = DDFP[ab10+9*3+2] - BA[2]*DDDP[ab00+3*3+2];
        DDDD[ab01+3*6+3] = DDFP[ab10+3*3+1] - BA[0]*DDDP[ab00+3*3+1];
        DDDD[ab01+3*6+4] = DDFP[ab10+6*3+2] - BA[1]*DDDP[ab00+3*3+2];
        DDDD[ab01+3*6+5] = DDFP[ab10+9*3+0] - BA[2]*DDDP[ab00+3*3+0];
        DDDD[ab01+4*6+0] = DDFP[ab10+9*3+0] - BA[0]*DDDP[ab00+4*3+0];
        DDDD[ab01+4*6+1] = DDFP[ab10+4*3+1] - BA[1]*DDDP[ab00+4*3+1];
        DDDD[ab01+4*6+2] = DDFP[ab10+7*3+2] - BA[2]*DDDP[ab00+4*3+2];
        DDDD[ab01+4*6+3] = DDFP[ab10+9*3+1] - BA[0]*DDDP[ab00+4*3+1];
        DDDD[ab01+4*6+4] = DDFP[ab10+4*3+2] - BA[1]*DDDP[ab00+4*3+2];
        DDDD[ab01+4*6+5] = DDFP[ab10+7*3+0] - BA[2]*DDDP[ab00+4*3+0];
        DDDD[ab01+5*6+0] = DDFP[ab10+8*3+0] - BA[0]*DDDP[ab00+5*3+0];
        DDDD[ab01+5*6+1] = DDFP[ab10+9*3+1] - BA[1]*DDDP[ab00+5*3+1];
        DDDD[ab01+5*6+2] = DDFP[ab10+5*3+2] - BA[2]*DDDP[ab00+5*3+2];
        DDDD[ab01+5*6+3] = DDFP[ab10+8*3+1] - BA[0]*DDDP[ab00+5*3+1];
        DDDD[ab01+5*6+4] = DDFP[ab10+9*3+2] - BA[1]*DDDP[ab00+5*3+2];
        DDDD[ab01+5*6+5] = DDFP[ab10+5*3+0] - BA[2]*DDDP[ab00+5*3+0];
    }
    // select integrals
    {
	int i,j, ijij, ij, n;
	double coe_i, coe_j;
	n=0;
	for (i=0; i<6; i++) {
	    coe_i = (i<3? 1.e0 : 3.e0 );
	    for ( j=0; j<6; j++) {
		coe_j = (j<3? 1.e0 : 3.e0 );
		ij   = i*6+j;
		ijij = ij*36+ij;
		DDDD_diag[n] = coe_i * coe_j * DDDD[ijij];
		n++;
	    }
	}
    }

    for ( i=0; i<36; i++ ) DDDD_diag[i]=fabs(DDDD_diag[i]);
    maxval = DDDD_diag[0];
    for ( i=1; i<36; i++ )
	if ( DDDD_diag[i]>maxval ) maxval = DDDD_diag[i];

    return maxval;
}
