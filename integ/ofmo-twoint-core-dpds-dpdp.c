/**
 * @file ofmo-twoint-core-dpds-dpdp.c
 * １つのCS４重対に対する２電子積分を計算する関数群。
 * 2011/06/16現在、(ss,ss)～(dd,dd)までの２１種類の２電子積分
 * 計算をサポートしている。
 * このファイルには、(dp,ds)～(dp,dp)タイプの２電子積分計算を
 * 行う関数が記述されている。
 *
 * 引数の関しては、@see ofmo-twoint-core-ssss-dspp.c に記述してある
 *
 * */
#include <stdio.h>
#include <math.h>
#include "fmt.h"

#define ZERO 0.e0
#define ONE  1.e0
#define HALF .5e0

extern double _twoint_inv2_;
extern double _twoint_inv3_;
extern double _twoint_spi2_;

extern double *FMT_fmt_table0;
extern double *FMT_fmt_table1;
extern double *FMT_fmt_table2;
extern double *FMT_fmt_table3;
extern double *FMT_fmt_table4;
extern double *FMT_fmt_table5;
extern double *FMT_fmt_table6;

extern double FMT_fmt_step_size;
extern double FMT_fmt_inv_step_size;
extern double FMT_pi_div2;

/** １つのCS４重対に対して(dp,ds)タイプの２電子積分を計算する関数
 * @ingroup core-twoint
 * */
void twoint_core_dpds__(
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double DPDS[6*3*6] ) {
    int ijps, klps, i, j, k, ix, m, m1;
    int c0;
    double cssss, zeta, dkab, xiza, eta, xizc, dk;
    double rz, re, ze2, ze22, ze23, zeta2, eta2;
    double tmp, tmp3[3], tmp6[6], tmp10[10];
    double sqrho, rho, PQ2, T;
    double PC[3], PA[3], WP[3], QP[3], QC[3], WQ[3];
    double ssss[5+1], psss[4+1][3], dsss[3+1][6], fsss[2+1][10];
    double psps[3*3], dsps[1+1][6*3], fsps[1+1][10*3];
    double DSDS[6*6], FSDS[10*6];
    double sqr3, coe, coe0;
    sqr3 = sqrt(3.e0);
    for ( i=0; i< 6*6; i++ ) DSDS[i] = ZERO;
    for ( i=0; i<10*6; i++ ) FSDS[i] = ZERO;
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
	    ze23  = 3.e0 * ze2;
	    for ( i=0; i<3; i++ ) {
		WP[i] = rz*QP[i];
		WQ[i] = rz*QP[i] - QP[i];
	    }
	    T     = rho * PQ2;
	    cssss = sqrho * dk;
	    {
		int it0, pos;
		double dT, t_inv, st_inv, dT2, dT3;
		if ( T < 46e0 ) {
		    it0 = (int)(0.5e0 + T * FMT_fmt_inv_step_size);
		    dT  = it0 * FMT_fmt_step_size - T;
		    dT2 = dT * _twoint_inv2_;
		    dT3 = dT * _twoint_inv3_;
		    pos = it0 * (5+4);
		    ssss[0] = cssss*(((FMT_fmt_table5[pos+3] * dT3
				    + FMT_fmt_table5[pos+2] ) * dT2
				    + FMT_fmt_table5[pos+1] ) * dT
				    + FMT_fmt_table5[pos+0] );
		    ssss[1] = cssss*(((FMT_fmt_table5[pos+4] * dT3
				    + FMT_fmt_table5[pos+3] ) * dT2
				    + FMT_fmt_table5[pos+2] ) * dT
				    + FMT_fmt_table5[pos+1] );
		    ssss[2] = cssss*(((FMT_fmt_table5[pos+5] * dT3
				    + FMT_fmt_table5[pos+4] ) * dT2
				    + FMT_fmt_table5[pos+3] ) * dT
				    + FMT_fmt_table5[pos+2] );
		    ssss[3] = cssss*(((FMT_fmt_table5[pos+6] * dT3
				    + FMT_fmt_table5[pos+5] ) * dT2
				    + FMT_fmt_table5[pos+4] ) * dT
				    + FMT_fmt_table5[pos+3] );
		    ssss[4] = cssss*(((FMT_fmt_table5[pos+7] * dT3
				    + FMT_fmt_table5[pos+6] ) * dT2
				    + FMT_fmt_table5[pos+5] ) * dT
				    + FMT_fmt_table5[pos+4] );
		    ssss[5] = cssss*(((FMT_fmt_table5[pos+8] * dT3
				    + FMT_fmt_table5[pos+7] ) * dT2
				    + FMT_fmt_table5[pos+6] ) * dT
				    + FMT_fmt_table5[pos+5] );
		} else {
		    st_inv = sqrt( 1.e0 / (T+T) );
		    t_inv  = st_inv * st_inv;
		    ssss[0] = cssss * _twoint_spi2_ * st_inv;
		    ssss[1] = t_inv * ssss[0];
		    ssss[2] = 3.e0 * t_inv * ssss[1];
		    ssss[3] = 5.e0 * t_inv * ssss[2];
		    ssss[4] = 7.e0 * t_inv * ssss[3];
		    ssss[5] = 9.e0 * t_inv * ssss[4];
		}
	    }
	    //fmt( ssss, 5, T, cssss );
	    // psss (m=0,4)
	    for ( m=0; m<=4; m++) {
		m1 = m+1;
		psss[m][0] = PA[0]*ssss[m] + WP[0]*ssss[m1];
		psss[m][1] = PA[1]*ssss[m] + WP[1]*ssss[m1];
		psss[m][2] = PA[2]*ssss[m] + WP[2]*ssss[m1];
	    }
	    // dsss (m=0,3)
	    for ( m=0; m<=3; m++) {
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
	    // fsss (m=0,2)
	    for ( m=0; m<=2; m++) {
		m1 = m+1;
		for ( i=0; i<3; i++) tmp3[i] = psss[m][i] - rz*psss[m1][i];

		fsss[m][0] = PA[0]*dsss[m][0] + WP[0]*dsss[m1][0]
			   + zeta*tmp3[0];
		fsss[m][1] = PA[1]*dsss[m][1] + WP[1]*dsss[m1][1]
			   + zeta*tmp3[1];
		fsss[m][2] = PA[2]*dsss[m][2] + WP[2]*dsss[m1][2]
			   + zeta*tmp3[2];
		fsss[m][3] = PA[0]*dsss[m][3] + WP[0]*dsss[m1][3]
			   + zeta2*tmp3[1];
		fsss[m][4] = PA[1]*dsss[m][4] + WP[1]*dsss[m1][4]
			   + zeta2*tmp3[2];
		fsss[m][5] = PA[2]*dsss[m][5] + WP[2]*dsss[m1][5]
			   + zeta2*tmp3[0];
		fsss[m][6] = PA[0]*dsss[m][1] + WP[0]*dsss[m1][1];
		fsss[m][7] = PA[1]*dsss[m][2] + WP[1]*dsss[m1][2];
		fsss[m][8] = PA[2]*dsss[m][0] + WP[2]*dsss[m1][0];
		fsss[m][9] = PA[0]*dsss[m][4] + WP[0]*dsss[m1][4];
	    }
	    // psps (m=1)
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
	    // dsps (m=0,1)
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
	    }
	    // fsps (m=0,1)
	    for ( m=0; m<=1; m++) {
		m1 = m+1;
		fsps[m][0*3+0] = QC[0]*fsss[m][0] + WQ[0]*fsss[m1][0]
			       + ze23*dsss[m1][0];
		fsps[m][0*3+1] = QC[1]*fsss[m][0] + WQ[1]*fsss[m1][0];
		fsps[m][0*3+2] = QC[2]*fsss[m][0] + WQ[2]*fsss[m1][0];
		fsps[m][1*3+0] = QC[0]*fsss[m][1] + WQ[0]*fsss[m1][1];
		fsps[m][1*3+1] = QC[1]*fsss[m][1] + WQ[1]*fsss[m1][1]
			       + ze23*dsss[m1][1];
		fsps[m][1*3+2] = QC[2]*fsss[m][1] + WQ[2]*fsss[m1][1];
		fsps[m][2*3+0] = QC[0]*fsss[m][2] + WQ[0]*fsss[m1][2];
		fsps[m][2*3+1] = QC[1]*fsss[m][2] + WQ[1]*fsss[m1][2];
		fsps[m][2*3+2] = QC[2]*fsss[m][2] + WQ[2]*fsss[m1][2]
			       + ze23*dsss[m1][2];
		fsps[m][3*3+0] = QC[0]*fsss[m][3] + WQ[0]*fsss[m1][3]
			       + ze22*dsss[m1][3];
		fsps[m][3*3+1] = QC[1]*fsss[m][3] + WQ[1]*fsss[m1][3]
			       + ze2*dsss[m1][0];
		fsps[m][3*3+2] = QC[2]*fsss[m][3] + WQ[2]*fsss[m1][3];
		fsps[m][4*3+0] = QC[0]*fsss[m][4] + WQ[0]*fsss[m1][4];
		fsps[m][4*3+1] = QC[1]*fsss[m][4] + WQ[1]*fsss[m1][4]
			       + ze22*dsss[m1][4];
		fsps[m][4*3+2] = QC[2]*fsss[m][4] + WQ[2]*fsss[m1][4]
			       + ze2*dsss[m1][1];
		fsps[m][5*3+0] = QC[0]*fsss[m][5] + WQ[0]*fsss[m1][5]
			       + ze2*dsss[m1][2];
		fsps[m][5*3+1] = QC[1]*fsss[m][5] + WQ[1]*fsss[m1][5];
		fsps[m][5*3+2] = QC[2]*fsss[m][5] + WQ[2]*fsss[m1][5]
			       + ze22*dsss[m1][5];
		fsps[m][6*3+0] = QC[0]*fsss[m][6] + WQ[0]*fsss[m1][6]
			       + ze2*dsss[m1][1];
		fsps[m][6*3+1] = QC[1]*fsss[m][6] + WQ[1]*fsss[m1][6]
			       + ze22*dsss[m1][3];
		fsps[m][6*3+2] = QC[2]*fsss[m][6] + WQ[2]*fsss[m1][6];
		fsps[m][7*3+0] = QC[0]*fsss[m][7] + WQ[0]*fsss[m1][7];
		fsps[m][7*3+1] = QC[1]*fsss[m][7] + WQ[1]*fsss[m1][7]
			       + ze2*dsss[m1][2];
		fsps[m][7*3+2] = QC[2]*fsss[m][7] + WQ[2]*fsss[m1][7]
			       + ze22*dsss[m1][4];
		fsps[m][8*3+0] = QC[0]*fsss[m][8] + WQ[0]*fsss[m1][8]
			       + ze22*dsss[m1][5];
		fsps[m][8*3+1] = QC[1]*fsss[m][8] + WQ[1]*fsss[m1][8];
		fsps[m][8*3+2] = QC[2]*fsss[m][8] + WQ[2]*fsss[m1][8]
			       + ze2*dsss[m1][0];
		fsps[m][9*3+0] = QC[0]*fsss[m][9] + WQ[0]*fsss[m1][9]
			       + ze2*dsss[m1][4];
		fsps[m][9*3+1] = QC[1]*fsss[m][9] + WQ[1]*fsss[m1][9]
			       + ze2*dsss[m1][5];
		fsps[m][9*3+2] = QC[2]*fsss[m][9] + WQ[2]*fsss[m1][9]
			       + ze2*dsss[m1][3];
	    }	// end m loop
	    // dsds
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
	    // fsds
	    for ( i=0; i<10; i++) tmp10[i] = fsss[0][i] - re*fsss[1][i];
	    FSDS[0*6+0] += QC[0]*fsps[0][0*3+0] + WQ[0]*fsps[1][0*3+0]
			+ eta2*tmp10[0] + ze23*dsps[1][0*3+0];
	    FSDS[0*6+1] += QC[1]*fsps[0][0*3+1] + WQ[1]*fsps[1][0*3+1]
			+ eta2*tmp10[0];
	    FSDS[0*6+2] += QC[2]*fsps[0][0*3+2] + WQ[2]*fsps[1][0*3+2]
			+ eta2*tmp10[0];
	    FSDS[0*6+3] += QC[0]*fsps[0][0*3+1] + WQ[0]*fsps[1][0*3+1]
			+ ze23*dsps[1][0*3+1];
	    FSDS[0*6+4] += QC[1]*fsps[0][0*3+2] + WQ[1]*fsps[1][0*3+2];
	    FSDS[0*6+5] += QC[2]*fsps[0][0*3+0] + WQ[2]*fsps[1][0*3+0];
	    FSDS[1*6+0] += QC[0]*fsps[0][1*3+0] + WQ[0]*fsps[1][1*3+0]
			+ eta2*tmp10[1];
	    FSDS[1*6+1] += QC[1]*fsps[0][1*3+1] + WQ[1]*fsps[1][1*3+1]
			+ eta2*tmp10[1] + ze23*dsps[1][1*3+1];
	    FSDS[1*6+2] += QC[2]*fsps[0][1*3+2] + WQ[2]*fsps[1][1*3+2]
			+ eta2*tmp10[1];
	    FSDS[1*6+3] += QC[0]*fsps[0][1*3+1] + WQ[0]*fsps[1][1*3+1];
	    FSDS[1*6+4] += QC[1]*fsps[0][1*3+2] + WQ[1]*fsps[1][1*3+2]
			+ ze23*dsps[1][1*3+2];
	    FSDS[1*6+5] += QC[2]*fsps[0][1*3+0] + WQ[2]*fsps[1][1*3+0];
	    FSDS[2*6+0] += QC[0]*fsps[0][2*3+0] + WQ[0]*fsps[1][2*3+0]
			+ eta2*tmp10[2];
	    FSDS[2*6+1] += QC[1]*fsps[0][2*3+1] + WQ[1]*fsps[1][2*3+1]
			+ eta2*tmp10[2];
	    FSDS[2*6+2] += QC[2]*fsps[0][2*3+2] + WQ[2]*fsps[1][2*3+2]
			+ eta2*tmp10[2] + ze23*dsps[1][2*3+2];
	    FSDS[2*6+3] += QC[0]*fsps[0][2*3+1] + WQ[0]*fsps[1][2*3+1];
	    FSDS[2*6+4] += QC[1]*fsps[0][2*3+2] + WQ[1]*fsps[1][2*3+2];
	    FSDS[2*6+5] += QC[2]*fsps[0][2*3+0] + WQ[2]*fsps[1][2*3+0]
			+ ze23*dsps[1][2*3+0];
	    FSDS[3*6+0] += QC[0]*fsps[0][3*3+0] + WQ[0]*fsps[1][3*3+0]
			+ eta2*tmp10[3] + ze22*dsps[1][3*3+0];
	    FSDS[3*6+1] += QC[1]*fsps[0][3*3+1] + WQ[1]*fsps[1][3*3+1]
			+ eta2*tmp10[3] + ze2*dsps[1][0*3+1];
	    FSDS[3*6+2] += QC[2]*fsps[0][3*3+2] + WQ[2]*fsps[1][3*3+2]
			+ eta2*tmp10[3];
	    FSDS[3*6+3] += QC[0]*fsps[0][3*3+1] + WQ[0]*fsps[1][3*3+1]
			+ ze22*dsps[1][3*3+1];
	    FSDS[3*6+4] += QC[1]*fsps[0][3*3+2] + WQ[1]*fsps[1][3*3+2]
			+ ze2*dsps[1][0*3+2];
	    FSDS[3*6+5] += QC[2]*fsps[0][3*3+0] + WQ[2]*fsps[1][3*3+0];
	    FSDS[4*6+0] += QC[0]*fsps[0][4*3+0] + WQ[0]*fsps[1][4*3+0]
			+ eta2*tmp10[4];
	    FSDS[4*6+1] += QC[1]*fsps[0][4*3+1] + WQ[1]*fsps[1][4*3+1]
			+ eta2*tmp10[4] + ze22*dsps[1][4*3+1];
	    FSDS[4*6+2] += QC[2]*fsps[0][4*3+2] + WQ[2]*fsps[1][4*3+2]
			+ eta2*tmp10[4] + ze2*dsps[1][1*3+2];
	    FSDS[4*6+3] += QC[0]*fsps[0][4*3+1] + WQ[0]*fsps[1][4*3+1];
	    FSDS[4*6+4] += QC[1]*fsps[0][4*3+2] + WQ[1]*fsps[1][4*3+2]
			+ ze22*dsps[1][4*3+2];
	    FSDS[4*6+5] += QC[2]*fsps[0][4*3+0] + WQ[2]*fsps[1][4*3+0]
			+ ze2*dsps[1][1*3+0];
	    FSDS[5*6+0] += QC[0]*fsps[0][5*3+0] + WQ[0]*fsps[1][5*3+0]
			+ eta2*tmp10[5] + ze2*dsps[1][2*3+0];
	    FSDS[5*6+1] += QC[1]*fsps[0][5*3+1] + WQ[1]*fsps[1][5*3+1]
			+ eta2*tmp10[5];
	    FSDS[5*6+2] += QC[2]*fsps[0][5*3+2] + WQ[2]*fsps[1][5*3+2]
			+ eta2*tmp10[5] + ze22*dsps[1][5*3+2];
	    FSDS[5*6+3] += QC[0]*fsps[0][5*3+1] + WQ[0]*fsps[1][5*3+1]
			+ ze2*dsps[1][2*3+1];
	    FSDS[5*6+4] += QC[1]*fsps[0][5*3+2] + WQ[1]*fsps[1][5*3+2];
	    FSDS[5*6+5] += QC[2]*fsps[0][5*3+0] + WQ[2]*fsps[1][5*3+0]
			+ ze22*dsps[1][5*3+0];
	    FSDS[6*6+0] += QC[0]*fsps[0][6*3+0] + WQ[0]*fsps[1][6*3+0]
			+ eta2*tmp10[6] + ze2*dsps[1][1*3+0];
	    FSDS[6*6+1] += QC[1]*fsps[0][6*3+1] + WQ[1]*fsps[1][6*3+1]
			+ eta2*tmp10[6] + ze22*dsps[1][3*3+1];
	    FSDS[6*6+2] += QC[2]*fsps[0][6*3+2] + WQ[2]*fsps[1][6*3+2]
			+ eta2*tmp10[6];
	    FSDS[6*6+3] += QC[0]*fsps[0][6*3+1] + WQ[0]*fsps[1][6*3+1]
			+ ze2*dsps[1][1*3+1];
	    FSDS[6*6+4] += QC[1]*fsps[0][6*3+2] + WQ[1]*fsps[1][6*3+2]
			+ ze22*dsps[1][3*3+2];
	    FSDS[6*6+5] += QC[2]*fsps[0][6*3+0] + WQ[2]*fsps[1][6*3+0];
	    FSDS[7*6+0] += QC[0]*fsps[0][7*3+0] + WQ[0]*fsps[1][7*3+0]
			+ eta2*tmp10[7];
	    FSDS[7*6+1] += QC[1]*fsps[0][7*3+1] + WQ[1]*fsps[1][7*3+1]
			+ eta2*tmp10[7] + ze2*dsps[1][2*3+1];
	    FSDS[7*6+2] += QC[2]*fsps[0][7*3+2] + WQ[2]*fsps[1][7*3+2]
			+ eta2*tmp10[7] + ze22*dsps[1][4*3+2];
	    FSDS[7*6+3] += QC[0]*fsps[0][7*3+1] + WQ[0]*fsps[1][7*3+1];
	    FSDS[7*6+4] += QC[1]*fsps[0][7*3+2] + WQ[1]*fsps[1][7*3+2]
			+ ze2*dsps[1][2*3+2];
	    FSDS[7*6+5] += QC[2]*fsps[0][7*3+0] + WQ[2]*fsps[1][7*3+0]
			+ ze22*dsps[1][4*3+0];
	    FSDS[8*6+0] += QC[0]*fsps[0][8*3+0] + WQ[0]*fsps[1][8*3+0]
			+ eta2*tmp10[8] + ze22*dsps[1][5*3+0];
	    FSDS[8*6+1] += QC[1]*fsps[0][8*3+1] + WQ[1]*fsps[1][8*3+1]
			+ eta2*tmp10[8];
	    FSDS[8*6+2] += QC[2]*fsps[0][8*3+2] + WQ[2]*fsps[1][8*3+2]
			+ eta2*tmp10[8] + ze2*dsps[1][0*3+2];
	    FSDS[8*6+3] += QC[0]*fsps[0][8*3+1] + WQ[0]*fsps[1][8*3+1]
			+ ze22*dsps[1][5*3+1];
	    FSDS[8*6+4] += QC[1]*fsps[0][8*3+2] + WQ[1]*fsps[1][8*3+2];
	    FSDS[8*6+5] += QC[2]*fsps[0][8*3+0] + WQ[2]*fsps[1][8*3+0]
			+ ze2*dsps[1][0*3+0];
	    FSDS[9*6+0] += QC[0]*fsps[0][9*3+0] + WQ[0]*fsps[1][9*3+0]
			+ eta2*tmp10[9] + ze2*dsps[1][4*3+0];
	    FSDS[9*6+1] += QC[1]*fsps[0][9*3+1] + WQ[1]*fsps[1][9*3+1]
			+ eta2*tmp10[9] + ze2*dsps[1][5*3+1];
	    FSDS[9*6+2] += QC[2]*fsps[0][9*3+2] + WQ[2]*fsps[1][9*3+2]
			+ eta2*tmp10[9] + ze2*dsps[1][3*3+2];
	    FSDS[9*6+3] += QC[0]*fsps[0][9*3+1] + WQ[0]*fsps[1][9*3+1]
			+ ze2*dsps[1][4*3+1];
	    FSDS[9*6+4] += QC[1]*fsps[0][9*3+2] + WQ[1]*fsps[1][9*3+2]
			+ ze2*dsps[1][5*3+2];
	    FSDS[9*6+5] += QC[2]*fsps[0][9*3+0] + WQ[2]*fsps[1][9*3+0]
			+ ze2*dsps[1][3*3+0];
	}	// klps
    }	// ijps
    // (D,P|D,S)
    for ( c0=0; c0<6; c0++) {
        DPDS[0*3*6+0*6+c0] = FSDS[0*6+c0] - BA[0]*DSDS[0*6+c0];
        DPDS[0*3*6+1*6+c0] = FSDS[3*6+c0] - BA[1]*DSDS[0*6+c0];
        DPDS[0*3*6+2*6+c0] = FSDS[8*6+c0] - BA[2]*DSDS[0*6+c0];
        DPDS[1*3*6+0*6+c0] = FSDS[6*6+c0] - BA[0]*DSDS[1*6+c0];
        DPDS[1*3*6+1*6+c0] = FSDS[1*6+c0] - BA[1]*DSDS[1*6+c0];
        DPDS[1*3*6+2*6+c0] = FSDS[4*6+c0] - BA[2]*DSDS[1*6+c0];
        DPDS[2*3*6+0*6+c0] = FSDS[5*6+c0] - BA[0]*DSDS[2*6+c0];
        DPDS[2*3*6+1*6+c0] = FSDS[7*6+c0] - BA[1]*DSDS[2*6+c0];
        DPDS[2*3*6+2*6+c0] = FSDS[2*6+c0] - BA[2]*DSDS[2*6+c0];
        DPDS[3*3*6+0*6+c0] = FSDS[3*6+c0] - BA[0]*DSDS[3*6+c0];
        DPDS[3*3*6+1*6+c0] = FSDS[6*6+c0] - BA[1]*DSDS[3*6+c0];
        DPDS[3*3*6+2*6+c0] = FSDS[9*6+c0] - BA[2]*DSDS[3*6+c0];
        DPDS[4*3*6+0*6+c0] = FSDS[9*6+c0] - BA[0]*DSDS[4*6+c0];
        DPDS[4*3*6+1*6+c0] = FSDS[4*6+c0] - BA[1]*DSDS[4*6+c0];
        DPDS[4*3*6+2*6+c0] = FSDS[7*6+c0] - BA[2]*DSDS[4*6+c0];
        DPDS[5*3*6+0*6+c0] = FSDS[8*6+c0] - BA[0]*DSDS[5*6+c0];
        DPDS[5*3*6+1*6+c0] = FSDS[9*6+c0] - BA[1]*DSDS[5*6+c0];
        DPDS[5*3*6+2*6+c0] = FSDS[5*6+c0] - BA[2]*DSDS[5*6+c0];
    }
    for ( i=0, ix=0; i<6; i++ ) {
	coe0 = (i<3? ONE : sqr3 );
	for ( j=0; j<3; j++ ) {
	    for ( k=0; k<6; k++, ix++ ) {
		coe = coe0 * ( k<3? ONE : sqr3 );
		DPDS[ix] *= coe;
	    }
	}
    }
}

/** １つのCS４重対に対して(dp,dp)タイプの２電子積分を計算する関数
 * @ingroup core-twoint
 * */
void twoint_core_dpdp__(
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double DPDP[6*3*6*3] ) {
    int ijps, klps, i, j, k, l, ix, m, m1;
    int ab, ab01, ab10, ab00, c0;
    double cssss, zeta, dkab, xiza, eta, xizc, dk;
    double rz, re, ze2, ze22, ze23, zeta2, eta2;
    union _temp_ {
        double entity[30];
        double ssss;
        double psss[3];
        double dsss[6];
        double fsss[10];
        double dsps[6*3];
        double fsps[10*3];
    } tmp;
    double sqrho, rho, PQ2, T;
    double PC[3], PA[3], WP[3], QP[3], QC[3], WQ[3];
    double ssss[6+1], psss[5+1][3], dsss[4+1][6], fsss[3+1][10];
    double ssps[3], psps[2+0][3*3], dsps[2+1][6*3], fsps[2+1][10*3];
    double psds[3*6], dsds[1+1][6*6], fsds[1+1][10*6];
    double DSFS[6*10], FSFS[10*10], DSDS[6*6], FSDS[10*6];
    double DPFS[6*3*10], DPDS[6*3*6];
    double sqr3, coe0, coe;
    sqr3 = sqrt(3.e0);
    for ( i=0; i< 6* 6; i++ ) DSDS[i] = ZERO;
    for ( i=0; i<10* 6; i++ ) FSDS[i] = ZERO;
    for ( i=0; i< 6*10; i++ ) DSFS[i] = ZERO;
    for ( i=0; i<10*10; i++ ) FSFS[i] = ZERO;
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
	    ze23  = 3.e0 * ze2;
	    for ( i=0; i<3; i++ ) {
		WP[i] = rz*QP[i];
		WQ[i] = rz*QP[i] - QP[i];
	    }
	    T     = rho * PQ2;
	    cssss = sqrho * dk;
	    {
		int it0, pos;
		double dT, t_inv, st_inv, dT2, dT3;
		if ( T < 48.e0 ) {
		    it0 = (int)(0.5e0 + T * FMT_fmt_inv_step_size);
		    dT  = it0 * FMT_fmt_step_size - T;
		    dT2 = dT * _twoint_inv2_;
		    dT3 = dT * _twoint_inv3_;
		    pos = it0 * (6+4);
		    ssss[0] = cssss*(((FMT_fmt_table6[pos+3]   * dT3
				     + FMT_fmt_table6[pos+2] ) * dT2
				     + FMT_fmt_table6[pos+1] ) * dT
				     + FMT_fmt_table6[pos+0] );
		    ssss[1] = cssss*(((FMT_fmt_table6[pos+4]   * dT3
				     + FMT_fmt_table6[pos+3] ) * dT2
				     + FMT_fmt_table6[pos+2] ) * dT
				     + FMT_fmt_table6[pos+1] );
		    ssss[2] = cssss*(((FMT_fmt_table6[pos+5]   * dT3
				     + FMT_fmt_table6[pos+4] ) * dT2
				     + FMT_fmt_table6[pos+3] ) * dT
				     + FMT_fmt_table6[pos+2] );
		    ssss[3] = cssss*(((FMT_fmt_table6[pos+6]   * dT3
				     + FMT_fmt_table6[pos+5] ) * dT2
				     + FMT_fmt_table6[pos+4] ) * dT
				     + FMT_fmt_table6[pos+3] );
		    ssss[4] = cssss*(((FMT_fmt_table6[pos+7]   * dT3
				     + FMT_fmt_table6[pos+6] ) * dT2
				     + FMT_fmt_table6[pos+5] ) * dT
				     + FMT_fmt_table6[pos+4] );
		    ssss[5] = cssss*(((FMT_fmt_table6[pos+8]   * dT3
				     + FMT_fmt_table6[pos+7] ) * dT2
				     + FMT_fmt_table6[pos+6] ) * dT
				     + FMT_fmt_table6[pos+5] );
		    ssss[6] = cssss*(((FMT_fmt_table6[pos+9]   * dT3
				     + FMT_fmt_table6[pos+8] ) * dT2
				     + FMT_fmt_table6[pos+7] ) * dT
				     + FMT_fmt_table6[pos+6] );
		} else {
		    st_inv = sqrt( 1.e0 / (T+T) );
		    t_inv  = st_inv * st_inv;
		    ssss[0] = cssss * _twoint_spi2_ * st_inv;
		    ssss[1] =         t_inv * ssss[0];
		    ssss[2] =  3.e0 * t_inv * ssss[1];
		    ssss[3] =  5.e0 * t_inv * ssss[2];
		    ssss[4] =  7.e0 * t_inv * ssss[3];
		    ssss[5] =  9.e0 * t_inv * ssss[4];
		    ssss[6] = 11.e0 * t_inv * ssss[5];
		}
	    }
	    //fmt( ssss, 6, T, cssss );
	    // psss (m=0,5)
	    for (m=0, m1=1; m<=5; m++, m1++) {
		psss[m][0] = PA[0]*ssss[m] + WP[0]*ssss[m1];
		psss[m][1] = PA[1]*ssss[m] + WP[1]*ssss[m1];
		psss[m][2] = PA[2]*ssss[m] + WP[2]*ssss[m1];
	    }
	    // dsss (m=0,4)
	    for (m=0, m1=1; m<=4; m++, m1++) {
		tmp.ssss = ssss[m] - rz*ssss[m1];

		dsss[m][0] = PA[0]*psss[m][0] + WP[0]*psss[m1][0]
			   + zeta2*tmp.ssss;
		dsss[m][1] = PA[1]*psss[m][1] + WP[1]*psss[m1][1]
			   + zeta2*tmp.ssss;
		dsss[m][2] = PA[2]*psss[m][2] + WP[2]*psss[m1][2]
			   + zeta2*tmp.ssss;
		dsss[m][3] = PA[0]*psss[m][1] + WP[0]*psss[m1][1];
		dsss[m][4] = PA[1]*psss[m][2] + WP[1]*psss[m1][2];
		dsss[m][5] = PA[2]*psss[m][0] + WP[2]*psss[m1][0];
	    }	// end m loop
	    // fsss (m=0,3)
	    for (m=0, m1=1; m<=3; m++, m1++) {
		for (i=0; i<3; i++) tmp.psss[i]=psss[m][i]-rz*psss[m1][i];
		fsss[m][0] = PA[0]*dsss[m][0] + WP[0]*dsss[m1][0]
			   + zeta*tmp.psss[0];
		fsss[m][1] = PA[1]*dsss[m][1] + WP[1]*dsss[m1][1]
			   + zeta*tmp.psss[1];
		fsss[m][2] = PA[2]*dsss[m][2] + WP[2]*dsss[m1][2]
			   + zeta*tmp.psss[2];
		fsss[m][3] = PA[0]*dsss[m][3] + WP[0]*dsss[m1][3]
			   + zeta2*tmp.psss[1];
		fsss[m][4] = PA[1]*dsss[m][4] + WP[1]*dsss[m1][4]
			   + zeta2*tmp.psss[2];
		fsss[m][5] = PA[2]*dsss[m][5] + WP[2]*dsss[m1][5]
			   + zeta2*tmp.psss[0];
		fsss[m][6] = PA[0]*dsss[m][1] + WP[0]*dsss[m1][1];
		fsss[m][7] = PA[1]*dsss[m][2] + WP[1]*dsss[m1][2];
		fsss[m][8] = PA[2]*dsss[m][0] + WP[2]*dsss[m1][0];
		fsss[m][9] = PA[0]*dsss[m][4] + WP[0]*dsss[m1][4];
	    }	// end m loop
	    // ssps (m=2)
	    ssps[0] = QC[0]*ssss[2] + WQ[0]*ssss[3];
	    ssps[1] = QC[1]*ssss[2] + WQ[1]*ssss[3];
	    ssps[2] = QC[2]*ssss[2] + WQ[2]*ssss[3];
	    // psps (m=1,2)
	    for (m=1, m1=2; m<=2; m++, m1++) {
		psps[m-1][0*3+0] = QC[0]*psss[m][0] + WQ[0]*psss[m1][0]
				 + ze2*ssss[m1];
		psps[m-1][0*3+1] = QC[1]*psss[m][0] + WQ[1]*psss[m1][0];
		psps[m-1][0*3+2] = QC[2]*psss[m][0] + WQ[2]*psss[m1][0];
		psps[m-1][1*3+0] = QC[0]*psss[m][1] + WQ[0]*psss[m1][1];
		psps[m-1][1*3+1] = QC[1]*psss[m][1] + WQ[1]*psss[m1][1]
				 + ze2*ssss[m1];
		psps[m-1][1*3+2] = QC[2]*psss[m][1] + WQ[2]*psss[m1][1];
		psps[m-1][2*3+0] = QC[0]*psss[m][2] + WQ[0]*psss[m1][2];
		psps[m-1][2*3+1] = QC[1]*psss[m][2] + WQ[1]*psss[m1][2];
		psps[m-1][2*3+2] = QC[2]*psss[m][2] + WQ[2]*psss[m1][2]
				 + ze2*ssss[m1];
	    }
	    // dsps (m=0,2)
	    for (m=0, m1=1; m<=2; m++, m1++) {
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
	    }
	    // fsps (m=0,2)
	    for (m=0, m1=1; m<=2; m++, m1++) {
		fsps[m][0*3+0] = QC[0]*fsss[m][0] + WQ[0]*fsss[m1][0]
			       + ze23*dsss[m1][0];
		fsps[m][0*3+1] = QC[1]*fsss[m][0] + WQ[1]*fsss[m1][0];
		fsps[m][0*3+2] = QC[2]*fsss[m][0] + WQ[2]*fsss[m1][0];
		fsps[m][1*3+0] = QC[0]*fsss[m][1] + WQ[0]*fsss[m1][1];
		fsps[m][1*3+1] = QC[1]*fsss[m][1] + WQ[1]*fsss[m1][1]
			       + ze23*dsss[m1][1];
		fsps[m][1*3+2] = QC[2]*fsss[m][1] + WQ[2]*fsss[m1][1];
		fsps[m][2*3+0] = QC[0]*fsss[m][2] + WQ[0]*fsss[m1][2];
		fsps[m][2*3+1] = QC[1]*fsss[m][2] + WQ[1]*fsss[m1][2];
		fsps[m][2*3+2] = QC[2]*fsss[m][2] + WQ[2]*fsss[m1][2]
			       + ze23*dsss[m1][2];
		fsps[m][3*3+0] = QC[0]*fsss[m][3] + WQ[0]*fsss[m1][3]
			       + ze22*dsss[m1][3];
		fsps[m][3*3+1] = QC[1]*fsss[m][3] + WQ[1]*fsss[m1][3]
			       + ze2*dsss[m1][0];
		fsps[m][3*3+2] = QC[2]*fsss[m][3] + WQ[2]*fsss[m1][3];
		fsps[m][4*3+0] = QC[0]*fsss[m][4] + WQ[0]*fsss[m1][4];
		fsps[m][4*3+1] = QC[1]*fsss[m][4] + WQ[1]*fsss[m1][4]
			       + ze22*dsss[m1][4];
		fsps[m][4*3+2] = QC[2]*fsss[m][4] + WQ[2]*fsss[m1][4]
			       + ze2*dsss[m1][1];
		fsps[m][5*3+0] = QC[0]*fsss[m][5] + WQ[0]*fsss[m1][5]
			       + ze2*dsss[m1][2];
		fsps[m][5*3+1] = QC[1]*fsss[m][5] + WQ[1]*fsss[m1][5];
		fsps[m][5*3+2] = QC[2]*fsss[m][5] + WQ[2]*fsss[m1][5]
			       + ze22*dsss[m1][5];
		fsps[m][6*3+0] = QC[0]*fsss[m][6] + WQ[0]*fsss[m1][6]
			       + ze2*dsss[m1][1];
		fsps[m][6*3+1] = QC[1]*fsss[m][6] + WQ[1]*fsss[m1][6]
			       + ze22*dsss[m1][3];
		fsps[m][6*3+2] = QC[2]*fsss[m][6] + WQ[2]*fsss[m1][6];
		fsps[m][7*3+0] = QC[0]*fsss[m][7] + WQ[0]*fsss[m1][7];
		fsps[m][7*3+1] = QC[1]*fsss[m][7] + WQ[1]*fsss[m1][7]
			       + ze2*dsss[m1][2];
		fsps[m][7*3+2] = QC[2]*fsss[m][7] + WQ[2]*fsss[m1][7]
			       + ze22*dsss[m1][4];
		fsps[m][8*3+0] = QC[0]*fsss[m][8] + WQ[0]*fsss[m1][8]
			       + ze22*dsss[m1][5];
		fsps[m][8*3+1] = QC[1]*fsss[m][8] + WQ[1]*fsss[m1][8];
		fsps[m][8*3+2] = QC[2]*fsss[m][8] + WQ[2]*fsss[m1][8]
			       + ze2*dsss[m1][0];
		fsps[m][9*3+0] = QC[0]*fsss[m][9] + WQ[0]*fsss[m1][9]
			       + ze2*dsss[m1][4];
		fsps[m][9*3+1] = QC[1]*fsss[m][9] + WQ[1]*fsss[m1][9]
			       + ze2*dsss[m1][5];
		fsps[m][9*3+2] = QC[2]*fsss[m][9] + WQ[2]*fsss[m1][9]
			       + ze2*dsss[m1][3];
	    }
	    // psds (m=1)
	    for (i=0; i<3; i++) tmp.psss[i]=psss[1][i]-re*psss[2][i];
	    psds[0*6+0] = QC[0]*psps[1-1][0*3+0] + WQ[0]*psps[2-1][0*3+0]
			+ eta2*tmp.psss[0] + ze2*ssps[0];
	    psds[0*6+1] = QC[1]*psps[1-1][0*3+1] + WQ[1]*psps[2-1][0*3+1]
			+ eta2*tmp.psss[0];
	    psds[0*6+2] = QC[2]*psps[1-1][0*3+2] + WQ[2]*psps[2-1][0*3+2]
			+ eta2*tmp.psss[0];
	    psds[0*6+3] = QC[0]*psps[1-1][0*3+1] + WQ[0]*psps[2-1][0*3+1]
			+ ze2*ssps[1];
	    psds[0*6+4] = QC[1]*psps[1-1][0*3+2] + WQ[1]*psps[2-1][0*3+2];
	    psds[0*6+5] = QC[2]*psps[1-1][0*3+0] + WQ[2]*psps[2-1][0*3+0];
	    psds[1*6+0] = QC[0]*psps[1-1][1*3+0] + WQ[0]*psps[2-1][1*3+0]
			+ eta2*tmp.psss[1];
	    psds[1*6+1] = QC[1]*psps[1-1][1*3+1] + WQ[1]*psps[2-1][1*3+1]
			+ eta2*tmp.psss[1] + ze2*ssps[1];
	    psds[1*6+2] = QC[2]*psps[1-1][1*3+2] + WQ[2]*psps[2-1][1*3+2]
			+ eta2*tmp.psss[1];
	    psds[1*6+3] = QC[0]*psps[1-1][1*3+1] + WQ[0]*psps[2-1][1*3+1];
	    psds[1*6+4] = QC[1]*psps[1-1][1*3+2] + WQ[1]*psps[2-1][1*3+2]
			+ ze2*ssps[2];
	    psds[1*6+5] = QC[2]*psps[1-1][1*3+0] + WQ[2]*psps[2-1][1*3+0];
	    psds[2*6+0] = QC[0]*psps[1-1][2*3+0] + WQ[0]*psps[2-1][2*3+0]
			+ eta2*tmp.psss[2];
	    psds[2*6+1] = QC[1]*psps[1-1][2*3+1] + WQ[1]*psps[2-1][2*3+1]
			+ eta2*tmp.psss[2];
	    psds[2*6+2] = QC[2]*psps[1-1][2*3+2] + WQ[2]*psps[2-1][2*3+2]
			+ eta2*tmp.psss[2] + ze2*ssps[2];
	    psds[2*6+3] = QC[0]*psps[1-1][2*3+1] + WQ[0]*psps[2-1][2*3+1];
	    psds[2*6+4] = QC[1]*psps[1-1][2*3+2] + WQ[1]*psps[2-1][2*3+2];
	    psds[2*6+5] = QC[2]*psps[1-1][2*3+0] + WQ[2]*psps[2-1][2*3+0]
			+ ze2*ssps[0];
	    // dsds (m=0,1)
	    for (m=0, m1=1; m<=1; m++, m1++) {
		for (i=0; i<6; i++) tmp.dsss[i]=dsss[m][i]-re*dsss[m1][i];
		dsds[m][0*6+0]=QC[0]*dsps[m][0*3+0]+WQ[0]*dsps[m1][0*3+0]
			       + eta2*tmp.dsss[0] + ze22*psps[m1-1][0*3+0];
		dsds[m][0*6+1]=QC[1]*dsps[m][0*3+1]+WQ[1]*dsps[m1][0*3+1]
			       + eta2*tmp.dsss[0];
		dsds[m][0*6+2]=QC[2]*dsps[m][0*3+2]+WQ[2]*dsps[m1][0*3+2]
			       + eta2*tmp.dsss[0];
		dsds[m][0*6+3]=QC[0]*dsps[m][0*3+1]+WQ[0]*dsps[m1][0*3+1]
			       + ze22*psps[m1-1][0*3+1];
		dsds[m][0*6+4]=QC[1]*dsps[m][0*3+2]+WQ[1]*dsps[m1][0*3+2];
		dsds[m][0*6+5]=QC[2]*dsps[m][0*3+0]+WQ[2]*dsps[m1][0*3+0];
		dsds[m][1*6+0]=QC[0]*dsps[m][1*3+0]+WQ[0]*dsps[m1][1*3+0]
			       + eta2*tmp.dsss[1];
		dsds[m][1*6+1]=QC[1]*dsps[m][1*3+1]+WQ[1]*dsps[m1][1*3+1]
			       + eta2*tmp.dsss[1] + ze22*psps[m1-1][1*3+1];
		dsds[m][1*6+2]=QC[2]*dsps[m][1*3+2]+WQ[2]*dsps[m1][1*3+2]
			       + eta2*tmp.dsss[1];
		dsds[m][1*6+3]=QC[0]*dsps[m][1*3+1]+WQ[0]*dsps[m1][1*3+1];
		dsds[m][1*6+4]=QC[1]*dsps[m][1*3+2]+WQ[1]*dsps[m1][1*3+2]
			       + ze22*psps[m1-1][1*3+2];
		dsds[m][1*6+5]=QC[2]*dsps[m][1*3+0]+WQ[2]*dsps[m1][1*3+0];
		dsds[m][2*6+0]=QC[0]*dsps[m][2*3+0]+WQ[0]*dsps[m1][2*3+0]
			       + eta2*tmp.dsss[2];
		dsds[m][2*6+1]=QC[1]*dsps[m][2*3+1]+WQ[1]*dsps[m1][2*3+1]
			       + eta2*tmp.dsss[2];
		dsds[m][2*6+2]=QC[2]*dsps[m][2*3+2]+WQ[2]*dsps[m1][2*3+2]
			       + eta2*tmp.dsss[2] + ze22*psps[m1-1][2*3+2];
		dsds[m][2*6+3]=QC[0]*dsps[m][2*3+1]+WQ[0]*dsps[m1][2*3+1];
		dsds[m][2*6+4]=QC[1]*dsps[m][2*3+2]+WQ[1]*dsps[m1][2*3+2];
		dsds[m][2*6+5]=QC[2]*dsps[m][2*3+0]+WQ[2]*dsps[m1][2*3+0]
			       + ze22*psps[m1-1][2*3+0];
		dsds[m][3*6+0]=QC[0]*dsps[m][3*3+0]+WQ[0]*dsps[m1][3*3+0]
			       + eta2*tmp.dsss[3] + ze2*psps[m1-1][1*3+0];
		dsds[m][3*6+1]=QC[1]*dsps[m][3*3+1]+WQ[1]*dsps[m1][3*3+1]
			       + eta2*tmp.dsss[3] + ze2*psps[m1-1][0*3+1];
		dsds[m][3*6+2]=QC[2]*dsps[m][3*3+2]+WQ[2]*dsps[m1][3*3+2]
			       + eta2*tmp.dsss[3];
		dsds[m][3*6+3]=QC[0]*dsps[m][3*3+1]+WQ[0]*dsps[m1][3*3+1]
			       + ze2*psps[m1-1][1*3+1];
		dsds[m][3*6+4]=QC[1]*dsps[m][3*3+2]+WQ[1]*dsps[m1][3*3+2]
			       + ze2*psps[m1-1][0*3+2];
		dsds[m][3*6+5]=QC[2]*dsps[m][3*3+0]+WQ[2]*dsps[m1][3*3+0];
		dsds[m][4*6+0]=QC[0]*dsps[m][4*3+0]+WQ[0]*dsps[m1][4*3+0]
			       + eta2*tmp.dsss[4];
		dsds[m][4*6+1]=QC[1]*dsps[m][4*3+1]+WQ[1]*dsps[m1][4*3+1]
			       + eta2*tmp.dsss[4] + ze2*psps[m1-1][2*3+1];
		dsds[m][4*6+2]=QC[2]*dsps[m][4*3+2]+WQ[2]*dsps[m1][4*3+2]
			       + eta2*tmp.dsss[4] + ze2*psps[m1-1][1*3+2];
		dsds[m][4*6+3]=QC[0]*dsps[m][4*3+1]+WQ[0]*dsps[m1][4*3+1];
		dsds[m][4*6+4]=QC[1]*dsps[m][4*3+2]+WQ[1]*dsps[m1][4*3+2]
			       + ze2*psps[m1-1][2*3+2];
		dsds[m][4*6+5]=QC[2]*dsps[m][4*3+0]+WQ[2]*dsps[m1][4*3+0]
			       + ze2*psps[m1-1][1*3+0];
		dsds[m][5*6+0]=QC[0]*dsps[m][5*3+0]+WQ[0]*dsps[m1][5*3+0]
			       + eta2*tmp.dsss[5] + ze2*psps[m1-1][2*3+0];
		dsds[m][5*6+1]=QC[1]*dsps[m][5*3+1]+WQ[1]*dsps[m1][5*3+1]
			       + eta2*tmp.dsss[5];
		dsds[m][5*6+2]=QC[2]*dsps[m][5*3+2]+WQ[2]*dsps[m1][5*3+2]
			       + eta2*tmp.dsss[5] + ze2*psps[m1-1][0*3+2];
		dsds[m][5*6+3]=QC[0]*dsps[m][5*3+1]+WQ[0]*dsps[m1][5*3+1]
			       + ze2*psps[m1-1][2*3+1];
		dsds[m][5*6+4]=QC[1]*dsps[m][5*3+2]+WQ[1]*dsps[m1][5*3+2];
		dsds[m][5*6+5]=QC[2]*dsps[m][5*3+0]+WQ[2]*dsps[m1][5*3+0]
			       + ze2*psps[m1-1][0*3+0];
	    }
	    for (i=0; i<6*6; i++) DSDS[i] += dsds[0][i];
	    // fsds (m=0,1)
	    for (m=0, m1=1; m<=1; m++, m1++) {
		for (i=0; i<10; i++) tmp.fsss[i]=fsss[m][i]-re*fsss[m1][i];
		fsds[m][0*6+0]=QC[0]*fsps[m][0*3+0]+WQ[0]*fsps[m1][0*3+0]
			       + eta2*tmp.fsss[0] + ze23*dsps[m1][0*3+0];
		fsds[m][0*6+1]=QC[1]*fsps[m][0*3+1]+WQ[1]*fsps[m1][0*3+1]
			       + eta2*tmp.fsss[0];
		fsds[m][0*6+2]=QC[2]*fsps[m][0*3+2]+WQ[2]*fsps[m1][0*3+2]
			       + eta2*tmp.fsss[0];
		fsds[m][0*6+3]=QC[0]*fsps[m][0*3+1]+WQ[0]*fsps[m1][0*3+1]
			       + ze23*dsps[m1][0*3+1];
		fsds[m][0*6+4]=QC[1]*fsps[m][0*3+2]+WQ[1]*fsps[m1][0*3+2];
		fsds[m][0*6+5]=QC[2]*fsps[m][0*3+0]+WQ[2]*fsps[m1][0*3+0];
		fsds[m][1*6+0]=QC[0]*fsps[m][1*3+0]+WQ[0]*fsps[m1][1*3+0]
			       + eta2*tmp.fsss[1];
		fsds[m][1*6+1]=QC[1]*fsps[m][1*3+1]+WQ[1]*fsps[m1][1*3+1]
			       + eta2*tmp.fsss[1] + ze23*dsps[m1][1*3+1];
		fsds[m][1*6+2]=QC[2]*fsps[m][1*3+2]+WQ[2]*fsps[m1][1*3+2]
			       + eta2*tmp.fsss[1];
		fsds[m][1*6+3]=QC[0]*fsps[m][1*3+1]+WQ[0]*fsps[m1][1*3+1];
		fsds[m][1*6+4]=QC[1]*fsps[m][1*3+2]+WQ[1]*fsps[m1][1*3+2]
			       + ze23*dsps[m1][1*3+2];
		fsds[m][1*6+5]=QC[2]*fsps[m][1*3+0]+WQ[2]*fsps[m1][1*3+0];
		fsds[m][2*6+0]=QC[0]*fsps[m][2*3+0]+WQ[0]*fsps[m1][2*3+0]
			       + eta2*tmp.fsss[2];
		fsds[m][2*6+1]=QC[1]*fsps[m][2*3+1]+WQ[1]*fsps[m1][2*3+1]
			       + eta2*tmp.fsss[2];
		fsds[m][2*6+2]=QC[2]*fsps[m][2*3+2]+WQ[2]*fsps[m1][2*3+2]
			       + eta2*tmp.fsss[2] + ze23*dsps[m1][2*3+2];
		fsds[m][2*6+3]=QC[0]*fsps[m][2*3+1]+WQ[0]*fsps[m1][2*3+1];
		fsds[m][2*6+4]=QC[1]*fsps[m][2*3+2]+WQ[1]*fsps[m1][2*3+2];
		fsds[m][2*6+5]=QC[2]*fsps[m][2*3+0]+WQ[2]*fsps[m1][2*3+0]
			       + ze23*dsps[m1][2*3+0];
		fsds[m][3*6+0]=QC[0]*fsps[m][3*3+0]+WQ[0]*fsps[m1][3*3+0]
			       + eta2*tmp.fsss[3] + ze22*dsps[m1][3*3+0];
		fsds[m][3*6+1]=QC[1]*fsps[m][3*3+1]+WQ[1]*fsps[m1][3*3+1]
			       + eta2*tmp.fsss[3] + ze2*dsps[m1][0*3+1];
		fsds[m][3*6+2]=QC[2]*fsps[m][3*3+2]+WQ[2]*fsps[m1][3*3+2]
			       + eta2*tmp.fsss[3];
		fsds[m][3*6+3]=QC[0]*fsps[m][3*3+1]+WQ[0]*fsps[m1][3*3+1]
			       + ze22*dsps[m1][3*3+1];
		fsds[m][3*6+4]=QC[1]*fsps[m][3*3+2]+WQ[1]*fsps[m1][3*3+2]
			       + ze2*dsps[m1][0*3+2];
		fsds[m][3*6+5]=QC[2]*fsps[m][3*3+0]+WQ[2]*fsps[m1][3*3+0];
		fsds[m][4*6+0]=QC[0]*fsps[m][4*3+0]+WQ[0]*fsps[m1][4*3+0]
			       + eta2*tmp.fsss[4];
		fsds[m][4*6+1]=QC[1]*fsps[m][4*3+1]+WQ[1]*fsps[m1][4*3+1]
			       + eta2*tmp.fsss[4] + ze22*dsps[m1][4*3+1];
		fsds[m][4*6+2]=QC[2]*fsps[m][4*3+2]+WQ[2]*fsps[m1][4*3+2]
			       + eta2*tmp.fsss[4] + ze2*dsps[m1][1*3+2];
		fsds[m][4*6+3]=QC[0]*fsps[m][4*3+1]+WQ[0]*fsps[m1][4*3+1];
		fsds[m][4*6+4]=QC[1]*fsps[m][4*3+2]+WQ[1]*fsps[m1][4*3+2]
			       + ze22*dsps[m1][4*3+2];
		fsds[m][4*6+5]=QC[2]*fsps[m][4*3+0]+WQ[2]*fsps[m1][4*3+0]
			       + ze2*dsps[m1][1*3+0];
		fsds[m][5*6+0]=QC[0]*fsps[m][5*3+0]+WQ[0]*fsps[m1][5*3+0]
			       + eta2*tmp.fsss[5] + ze2*dsps[m1][2*3+0];
		fsds[m][5*6+1]=QC[1]*fsps[m][5*3+1]+WQ[1]*fsps[m1][5*3+1]
			       + eta2*tmp.fsss[5];
		fsds[m][5*6+2]=QC[2]*fsps[m][5*3+2]+WQ[2]*fsps[m1][5*3+2]
			       + eta2*tmp.fsss[5] + ze22*dsps[m1][5*3+2];
		fsds[m][5*6+3]=QC[0]*fsps[m][5*3+1]+WQ[0]*fsps[m1][5*3+1]
			       + ze2*dsps[m1][2*3+1];
		fsds[m][5*6+4]=QC[1]*fsps[m][5*3+2]+WQ[1]*fsps[m1][5*3+2];
		fsds[m][5*6+5]=QC[2]*fsps[m][5*3+0]+WQ[2]*fsps[m1][5*3+0]
			       + ze22*dsps[m1][5*3+0];
		fsds[m][6*6+0]=QC[0]*fsps[m][6*3+0]+WQ[0]*fsps[m1][6*3+0]
			       + eta2*tmp.fsss[6] + ze2*dsps[m1][1*3+0];
		fsds[m][6*6+1]=QC[1]*fsps[m][6*3+1]+WQ[1]*fsps[m1][6*3+1]
			       + eta2*tmp.fsss[6] + ze22*dsps[m1][3*3+1];
		fsds[m][6*6+2]=QC[2]*fsps[m][6*3+2]+WQ[2]*fsps[m1][6*3+2]
			       + eta2*tmp.fsss[6];
		fsds[m][6*6+3]=QC[0]*fsps[m][6*3+1]+WQ[0]*fsps[m1][6*3+1]
			       + ze2*dsps[m1][1*3+1];
		fsds[m][6*6+4]=QC[1]*fsps[m][6*3+2]+WQ[1]*fsps[m1][6*3+2]
			       + ze22*dsps[m1][3*3+2];
		fsds[m][6*6+5]=QC[2]*fsps[m][6*3+0]+WQ[2]*fsps[m1][6*3+0];
		fsds[m][7*6+0]=QC[0]*fsps[m][7*3+0]+WQ[0]*fsps[m1][7*3+0]
			       + eta2*tmp.fsss[7];
		fsds[m][7*6+1]=QC[1]*fsps[m][7*3+1]+WQ[1]*fsps[m1][7*3+1]
			       + eta2*tmp.fsss[7] + ze2*dsps[m1][2*3+1];
		fsds[m][7*6+2]=QC[2]*fsps[m][7*3+2]+WQ[2]*fsps[m1][7*3+2]
			       + eta2*tmp.fsss[7] + ze22*dsps[m1][4*3+2];
		fsds[m][7*6+3]=QC[0]*fsps[m][7*3+1]+WQ[0]*fsps[m1][7*3+1];
		fsds[m][7*6+4]=QC[1]*fsps[m][7*3+2]+WQ[1]*fsps[m1][7*3+2]
			       + ze2*dsps[m1][2*3+2];
		fsds[m][7*6+5]=QC[2]*fsps[m][7*3+0]+WQ[2]*fsps[m1][7*3+0]
			       + ze22*dsps[m1][4*3+0];
		fsds[m][8*6+0]=QC[0]*fsps[m][8*3+0]+WQ[0]*fsps[m1][8*3+0]
			       + eta2*tmp.fsss[8] + ze22*dsps[m1][5*3+0];
		fsds[m][8*6+1]=QC[1]*fsps[m][8*3+1]+WQ[1]*fsps[m1][8*3+1]
			       + eta2*tmp.fsss[8];
		fsds[m][8*6+2]=QC[2]*fsps[m][8*3+2]+WQ[2]*fsps[m1][8*3+2]
			       + eta2*tmp.fsss[8] + ze2*dsps[m1][0*3+2];
		fsds[m][8*6+3]=QC[0]*fsps[m][8*3+1]+WQ[0]*fsps[m1][8*3+1]
			       + ze22*dsps[m1][5*3+1];
		fsds[m][8*6+4]=QC[1]*fsps[m][8*3+2]+WQ[1]*fsps[m1][8*3+2];
		fsds[m][8*6+5]=QC[2]*fsps[m][8*3+0]+WQ[2]*fsps[m1][8*3+0]
			       + ze2*dsps[m1][0*3+0];
		fsds[m][9*6+0]=QC[0]*fsps[m][9*3+0]+WQ[0]*fsps[m1][9*3+0]
			       + eta2*tmp.fsss[9] + ze2*dsps[m1][4*3+0];
		fsds[m][9*6+1]=QC[1]*fsps[m][9*3+1]+WQ[1]*fsps[m1][9*3+1]
			       + eta2*tmp.fsss[9] + ze2*dsps[m1][5*3+1];
		fsds[m][9*6+2]=QC[2]*fsps[m][9*3+2]+WQ[2]*fsps[m1][9*3+2]
			       + eta2*tmp.fsss[9] + ze2*dsps[m1][3*3+2];
		fsds[m][9*6+3]=QC[0]*fsps[m][9*3+1]+WQ[0]*fsps[m1][9*3+1]
			       + ze2*dsps[m1][4*3+1];
		fsds[m][9*6+4]=QC[1]*fsps[m][9*3+2]+WQ[1]*fsps[m1][9*3+2]
			       + ze2*dsps[m1][5*3+2];
		fsds[m][9*6+5]=QC[2]*fsps[m][9*3+0]+WQ[2]*fsps[m1][9*3+0]
			       + ze2*dsps[m1][3*3+0];
	    }
	    for (i=0; i<10*6; i++) FSDS[i] += fsds[0][i];
	    // dsfs (m=0)
	    for (i=0; i<6*3; i++) tmp.dsps[i]=dsps[0][i]-re*dsps[1][i];
	    DSFS[0*10+0]+=QC[0]*dsds[0][0*6+0] + WQ[0]*dsds[1][0*6+0]
			 + eta*tmp.dsps[0*3+0] + ze22*psds[0*6+0];
	    DSFS[0*10+1]+=QC[1]*dsds[0][0*6+1] + WQ[1]*dsds[1][0*6+1]
			 + eta*tmp.dsps[0*3+1];
	    DSFS[0*10+2]+=QC[2]*dsds[0][0*6+2] + WQ[2]*dsds[1][0*6+2]
			 + eta*tmp.dsps[0*3+2];
	    DSFS[0*10+3]+=QC[0]*dsds[0][0*6+3] + WQ[0]*dsds[1][0*6+3]
			 + eta2*tmp.dsps[0*3+1] + ze22*psds[0*6+3];
	    DSFS[0*10+4]+=QC[1]*dsds[0][0*6+4] + WQ[1]*dsds[1][0*6+4]
			 + eta2*tmp.dsps[0*3+2];
	    DSFS[0*10+5]+=QC[2]*dsds[0][0*6+5] + WQ[2]*dsds[1][0*6+5]
			 + eta2*tmp.dsps[0*3+0];
	    DSFS[0*10+6]+=QC[0]*dsds[0][0*6+1] + WQ[0]*dsds[1][0*6+1]
			 + ze22*psds[0*6+1];
	    DSFS[0*10+7]+=QC[1]*dsds[0][0*6+2] + WQ[1]*dsds[1][0*6+2];
	    DSFS[0*10+8]+=QC[2]*dsds[0][0*6+0] + WQ[2]*dsds[1][0*6+0];
	    DSFS[0*10+9]+=QC[0]*dsds[0][0*6+4] + WQ[0]*dsds[1][0*6+4]
			 + ze22*psds[0*6+4];
	    DSFS[1*10+0]+=QC[0]*dsds[0][1*6+0] + WQ[0]*dsds[1][1*6+0]
			 + eta*tmp.dsps[1*3+0];
	    DSFS[1*10+1]+=QC[1]*dsds[0][1*6+1] + WQ[1]*dsds[1][1*6+1]
			 + eta*tmp.dsps[1*3+1] + ze22*psds[1*6+1];
	    DSFS[1*10+2]+=QC[2]*dsds[0][1*6+2] + WQ[2]*dsds[1][1*6+2]
			 + eta*tmp.dsps[1*3+2];
	    DSFS[1*10+3]+=QC[0]*dsds[0][1*6+3] + WQ[0]*dsds[1][1*6+3]
			 + eta2*tmp.dsps[1*3+1];
	    DSFS[1*10+4]+=QC[1]*dsds[0][1*6+4] + WQ[1]*dsds[1][1*6+4]
			 + eta2*tmp.dsps[1*3+2] + ze22*psds[1*6+4];
	    DSFS[1*10+5]+=QC[2]*dsds[0][1*6+5] + WQ[2]*dsds[1][1*6+5]
			 + eta2*tmp.dsps[1*3+0];
	    DSFS[1*10+6]+=QC[0]*dsds[0][1*6+1] + WQ[0]*dsds[1][1*6+1];
	    DSFS[1*10+7]+=QC[1]*dsds[0][1*6+2] + WQ[1]*dsds[1][1*6+2]
			 + ze22*psds[1*6+2];
	    DSFS[1*10+8]+=QC[2]*dsds[0][1*6+0] + WQ[2]*dsds[1][1*6+0];
	    DSFS[1*10+9]+=QC[0]*dsds[0][1*6+4] + WQ[0]*dsds[1][1*6+4];
	    DSFS[2*10+0]+=QC[0]*dsds[0][2*6+0] + WQ[0]*dsds[1][2*6+0]
			 + eta*tmp.dsps[2*3+0];
	    DSFS[2*10+1]+=QC[1]*dsds[0][2*6+1] + WQ[1]*dsds[1][2*6+1]
			 + eta*tmp.dsps[2*3+1];
	    DSFS[2*10+2]+=QC[2]*dsds[0][2*6+2] + WQ[2]*dsds[1][2*6+2]
			 + eta*tmp.dsps[2*3+2] + ze22*psds[2*6+2];
	    DSFS[2*10+3]+=QC[0]*dsds[0][2*6+3] + WQ[0]*dsds[1][2*6+3]
			 + eta2*tmp.dsps[2*3+1];
	    DSFS[2*10+4]+=QC[1]*dsds[0][2*6+4] + WQ[1]*dsds[1][2*6+4]
			 + eta2*tmp.dsps[2*3+2];
	    DSFS[2*10+5]+=QC[2]*dsds[0][2*6+5] + WQ[2]*dsds[1][2*6+5]
			 + eta2*tmp.dsps[2*3+0] + ze22*psds[2*6+5];
	    DSFS[2*10+6]+=QC[0]*dsds[0][2*6+1] + WQ[0]*dsds[1][2*6+1];
	    DSFS[2*10+7]+=QC[1]*dsds[0][2*6+2] + WQ[1]*dsds[1][2*6+2];
	    DSFS[2*10+8]+=QC[2]*dsds[0][2*6+0] + WQ[2]*dsds[1][2*6+0]
			 + ze22*psds[2*6+0];
	    DSFS[2*10+9]+=QC[0]*dsds[0][2*6+4] + WQ[0]*dsds[1][2*6+4];
	    DSFS[3*10+0]+=QC[0]*dsds[0][3*6+0] + WQ[0]*dsds[1][3*6+0]
			 + eta*tmp.dsps[3*3+0] + ze2*psds[1*6+0];
	    DSFS[3*10+1]+=QC[1]*dsds[0][3*6+1] + WQ[1]*dsds[1][3*6+1]
			 + eta*tmp.dsps[3*3+1] + ze2*psds[0*6+1];
	    DSFS[3*10+2]+=QC[2]*dsds[0][3*6+2] + WQ[2]*dsds[1][3*6+2]
			 + eta*tmp.dsps[3*3+2];
	    DSFS[3*10+3]+=QC[0]*dsds[0][3*6+3] + WQ[0]*dsds[1][3*6+3]
			 + eta2*tmp.dsps[3*3+1] + ze2*psds[1*6+3];
	    DSFS[3*10+4]+=QC[1]*dsds[0][3*6+4] + WQ[1]*dsds[1][3*6+4]
			 + eta2*tmp.dsps[3*3+2] + ze2*psds[0*6+4];
	    DSFS[3*10+5]+=QC[2]*dsds[0][3*6+5] + WQ[2]*dsds[1][3*6+5]
			 + eta2*tmp.dsps[3*3+0];
	    DSFS[3*10+6]+=QC[0]*dsds[0][3*6+1] + WQ[0]*dsds[1][3*6+1]
			 + ze2*psds[1*6+1];
	    DSFS[3*10+7]+=QC[1]*dsds[0][3*6+2] + WQ[1]*dsds[1][3*6+2]
			 + ze2*psds[0*6+2];
	    DSFS[3*10+8]+=QC[2]*dsds[0][3*6+0] + WQ[2]*dsds[1][3*6+0];
	    DSFS[3*10+9]+=QC[0]*dsds[0][3*6+4] + WQ[0]*dsds[1][3*6+4]
			 + ze2*psds[1*6+4];
	    DSFS[4*10+0]+=QC[0]*dsds[0][4*6+0] + WQ[0]*dsds[1][4*6+0]
			 + eta*tmp.dsps[4*3+0];
	    DSFS[4*10+1]+=QC[1]*dsds[0][4*6+1] + WQ[1]*dsds[1][4*6+1]
			 + eta*tmp.dsps[4*3+1] + ze2*psds[2*6+1];
	    DSFS[4*10+2]+=QC[2]*dsds[0][4*6+2] + WQ[2]*dsds[1][4*6+2]
			 + eta*tmp.dsps[4*3+2] + ze2*psds[1*6+2];
	    DSFS[4*10+3]+=QC[0]*dsds[0][4*6+3] + WQ[0]*dsds[1][4*6+3]
			 + eta2*tmp.dsps[4*3+1];
	    DSFS[4*10+4]+=QC[1]*dsds[0][4*6+4] + WQ[1]*dsds[1][4*6+4]
			 + eta2*tmp.dsps[4*3+2] + ze2*psds[2*6+4];
	    DSFS[4*10+5]+=QC[2]*dsds[0][4*6+5] + WQ[2]*dsds[1][4*6+5]
			 + eta2*tmp.dsps[4*3+0] + ze2*psds[1*6+5];
	    DSFS[4*10+6]+=QC[0]*dsds[0][4*6+1] + WQ[0]*dsds[1][4*6+1];
	    DSFS[4*10+7]+=QC[1]*dsds[0][4*6+2] + WQ[1]*dsds[1][4*6+2]
			 + ze2*psds[2*6+2];
	    DSFS[4*10+8]+=QC[2]*dsds[0][4*6+0] + WQ[2]*dsds[1][4*6+0]
			 + ze2*psds[1*6+0];
	    DSFS[4*10+9]+=QC[0]*dsds[0][4*6+4] + WQ[0]*dsds[1][4*6+4];
	    DSFS[5*10+0]+=QC[0]*dsds[0][5*6+0] + WQ[0]*dsds[1][5*6+0]
			 + eta*tmp.dsps[5*3+0] + ze2*psds[2*6+0];
	    DSFS[5*10+1]+=QC[1]*dsds[0][5*6+1] + WQ[1]*dsds[1][5*6+1]
			 + eta*tmp.dsps[5*3+1];
	    DSFS[5*10+2]+=QC[2]*dsds[0][5*6+2] + WQ[2]*dsds[1][5*6+2]
			 + eta*tmp.dsps[5*3+2] + ze2*psds[0*6+2];
	    DSFS[5*10+3]+=QC[0]*dsds[0][5*6+3] + WQ[0]*dsds[1][5*6+3]
			 + eta2*tmp.dsps[5*3+1] + ze2*psds[2*6+3];
	    DSFS[5*10+4]+=QC[1]*dsds[0][5*6+4] + WQ[1]*dsds[1][5*6+4]
			 + eta2*tmp.dsps[5*3+2];
	    DSFS[5*10+5]+=QC[2]*dsds[0][5*6+5] + WQ[2]*dsds[1][5*6+5]
			 + eta2*tmp.dsps[5*3+0] + ze2*psds[0*6+5];
	    DSFS[5*10+6]+=QC[0]*dsds[0][5*6+1] + WQ[0]*dsds[1][5*6+1]
			 + ze2*psds[2*6+1];
	    DSFS[5*10+7]+=QC[1]*dsds[0][5*6+2] + WQ[1]*dsds[1][5*6+2];
	    DSFS[5*10+8]+=QC[2]*dsds[0][5*6+0] + WQ[2]*dsds[1][5*6+0]
			 + ze2*psds[0*6+0];
	    DSFS[5*10+9]+=QC[0]*dsds[0][5*6+4] + WQ[0]*dsds[1][5*6+4]
			 + ze2*psds[2*6+4];
	    // fsfs (m=0)
	    for (i=0; i<10*3; i++) tmp.fsps[i]=fsps[0][i]-re*fsps[1][i];
	    FSFS[0*10+0]+=QC[0]*fsds[0][0*6+0] + WQ[0]*fsds[1][0*6+0]
			 + eta*tmp.fsps[0*3+0] + ze23*dsds[1][0*6+0];
	    FSFS[0*10+1]+=QC[1]*fsds[0][0*6+1] + WQ[1]*fsds[1][0*6+1]
			 + eta*tmp.fsps[0*3+1];
	    FSFS[0*10+2]+=QC[2]*fsds[0][0*6+2] + WQ[2]*fsds[1][0*6+2]
			 + eta*tmp.fsps[0*3+2];
	    FSFS[0*10+3]+=QC[0]*fsds[0][0*6+3] + WQ[0]*fsds[1][0*6+3]
			 + eta2*tmp.fsps[0*3+1] + ze23*dsds[1][0*6+3];
	    FSFS[0*10+4]+=QC[1]*fsds[0][0*6+4] + WQ[1]*fsds[1][0*6+4]
			 + eta2*tmp.fsps[0*3+2];
	    FSFS[0*10+5]+=QC[2]*fsds[0][0*6+5] + WQ[2]*fsds[1][0*6+5]
			 + eta2*tmp.fsps[0*3+0];
	    FSFS[0*10+6]+=QC[0]*fsds[0][0*6+1] + WQ[0]*fsds[1][0*6+1]
			 + ze23*dsds[1][0*6+1];
	    FSFS[0*10+7]+=QC[1]*fsds[0][0*6+2] + WQ[1]*fsds[1][0*6+2];
	    FSFS[0*10+8]+=QC[2]*fsds[0][0*6+0] + WQ[2]*fsds[1][0*6+0];
	    FSFS[0*10+9]+=QC[0]*fsds[0][0*6+4] + WQ[0]*fsds[1][0*6+4]
			 + ze23*dsds[1][0*6+4];
	    FSFS[1*10+0]+=QC[0]*fsds[0][1*6+0] + WQ[0]*fsds[1][1*6+0]
			 + eta*tmp.fsps[1*3+0];
	    FSFS[1*10+1]+=QC[1]*fsds[0][1*6+1] + WQ[1]*fsds[1][1*6+1]
			 + eta*tmp.fsps[1*3+1] + ze23*dsds[1][1*6+1];
	    FSFS[1*10+2]+=QC[2]*fsds[0][1*6+2] + WQ[2]*fsds[1][1*6+2]
			 + eta*tmp.fsps[1*3+2];
	    FSFS[1*10+3]+=QC[0]*fsds[0][1*6+3] + WQ[0]*fsds[1][1*6+3]
			 + eta2*tmp.fsps[1*3+1];
	    FSFS[1*10+4]+=QC[1]*fsds[0][1*6+4] + WQ[1]*fsds[1][1*6+4]
			 + eta2*tmp.fsps[1*3+2] + ze23*dsds[1][1*6+4];
	    FSFS[1*10+5]+=QC[2]*fsds[0][1*6+5] + WQ[2]*fsds[1][1*6+5]
			 + eta2*tmp.fsps[1*3+0];
	    FSFS[1*10+6]+=QC[0]*fsds[0][1*6+1] + WQ[0]*fsds[1][1*6+1];
	    FSFS[1*10+7]+=QC[1]*fsds[0][1*6+2] + WQ[1]*fsds[1][1*6+2]
			 + ze23*dsds[1][1*6+2];
	    FSFS[1*10+8]+=QC[2]*fsds[0][1*6+0] + WQ[2]*fsds[1][1*6+0];
	    FSFS[1*10+9]+=QC[0]*fsds[0][1*6+4] + WQ[0]*fsds[1][1*6+4];
	    FSFS[2*10+0]+=QC[0]*fsds[0][2*6+0] + WQ[0]*fsds[1][2*6+0]
			 + eta*tmp.fsps[2*3+0];
	    FSFS[2*10+1]+=QC[1]*fsds[0][2*6+1] + WQ[1]*fsds[1][2*6+1]
			 + eta*tmp.fsps[2*3+1];
	    FSFS[2*10+2]+=QC[2]*fsds[0][2*6+2] + WQ[2]*fsds[1][2*6+2]
			 + eta*tmp.fsps[2*3+2] + ze23*dsds[1][2*6+2];
	    FSFS[2*10+3]+=QC[0]*fsds[0][2*6+3] + WQ[0]*fsds[1][2*6+3]
			 + eta2*tmp.fsps[2*3+1];
	    FSFS[2*10+4]+=QC[1]*fsds[0][2*6+4] + WQ[1]*fsds[1][2*6+4]
			 + eta2*tmp.fsps[2*3+2];
	    FSFS[2*10+5]+=QC[2]*fsds[0][2*6+5] + WQ[2]*fsds[1][2*6+5]
			 + eta2*tmp.fsps[2*3+0] + ze23*dsds[1][2*6+5];
	    FSFS[2*10+6]+=QC[0]*fsds[0][2*6+1] + WQ[0]*fsds[1][2*6+1];
	    FSFS[2*10+7]+=QC[1]*fsds[0][2*6+2] + WQ[1]*fsds[1][2*6+2];
	    FSFS[2*10+8]+=QC[2]*fsds[0][2*6+0] + WQ[2]*fsds[1][2*6+0]
			 + ze23*dsds[1][2*6+0];
	    FSFS[2*10+9]+=QC[0]*fsds[0][2*6+4] + WQ[0]*fsds[1][2*6+4];
	    FSFS[3*10+0]+=QC[0]*fsds[0][3*6+0] + WQ[0]*fsds[1][3*6+0]
			 + eta*tmp.fsps[3*3+0] + ze22*dsds[1][3*6+0];
	    FSFS[3*10+1]+=QC[1]*fsds[0][3*6+1] + WQ[1]*fsds[1][3*6+1]
			 + eta*tmp.fsps[3*3+1] + ze2*dsds[1][0*6+1];
	    FSFS[3*10+2]+=QC[2]*fsds[0][3*6+2] + WQ[2]*fsds[1][3*6+2]
			 + eta*tmp.fsps[3*3+2];
	    FSFS[3*10+3]+=QC[0]*fsds[0][3*6+3] + WQ[0]*fsds[1][3*6+3]
			 + eta2*tmp.fsps[3*3+1] + ze22*dsds[1][3*6+3];
	    FSFS[3*10+4]+=QC[1]*fsds[0][3*6+4] + WQ[1]*fsds[1][3*6+4]
			 + eta2*tmp.fsps[3*3+2] + ze2*dsds[1][0*6+4];
	    FSFS[3*10+5]+=QC[2]*fsds[0][3*6+5] + WQ[2]*fsds[1][3*6+5]
			 + eta2*tmp.fsps[3*3+0];
	    FSFS[3*10+6]+=QC[0]*fsds[0][3*6+1] + WQ[0]*fsds[1][3*6+1]
			 + ze22*dsds[1][3*6+1];
	    FSFS[3*10+7]+=QC[1]*fsds[0][3*6+2] + WQ[1]*fsds[1][3*6+2]
			 + ze2*dsds[1][0*6+2];
	    FSFS[3*10+8]+=QC[2]*fsds[0][3*6+0] + WQ[2]*fsds[1][3*6+0];
	    FSFS[3*10+9]+=QC[0]*fsds[0][3*6+4] + WQ[0]*fsds[1][3*6+4]
			 + ze22*dsds[1][3*6+4];
	    FSFS[4*10+0]+=QC[0]*fsds[0][4*6+0] + WQ[0]*fsds[1][4*6+0]
			 + eta*tmp.fsps[4*3+0];
	    FSFS[4*10+1]+=QC[1]*fsds[0][4*6+1] + WQ[1]*fsds[1][4*6+1]
			 + eta*tmp.fsps[4*3+1] + ze22*dsds[1][4*6+1];
	    FSFS[4*10+2]+=QC[2]*fsds[0][4*6+2] + WQ[2]*fsds[1][4*6+2]
			 + eta*tmp.fsps[4*3+2] + ze2*dsds[1][1*6+2];
	    FSFS[4*10+3]+=QC[0]*fsds[0][4*6+3] + WQ[0]*fsds[1][4*6+3]
			 + eta2*tmp.fsps[4*3+1];
	    FSFS[4*10+4]+=QC[1]*fsds[0][4*6+4] + WQ[1]*fsds[1][4*6+4]
			 + eta2*tmp.fsps[4*3+2] + ze22*dsds[1][4*6+4];
	    FSFS[4*10+5]+=QC[2]*fsds[0][4*6+5] + WQ[2]*fsds[1][4*6+5]
			 + eta2*tmp.fsps[4*3+0] + ze2*dsds[1][1*6+5];
	    FSFS[4*10+6]+=QC[0]*fsds[0][4*6+1] + WQ[0]*fsds[1][4*6+1];
	    FSFS[4*10+7]+=QC[1]*fsds[0][4*6+2] + WQ[1]*fsds[1][4*6+2]
			 + ze22*dsds[1][4*6+2];
	    FSFS[4*10+8]+=QC[2]*fsds[0][4*6+0] + WQ[2]*fsds[1][4*6+0]
			 + ze2*dsds[1][1*6+0];
	    FSFS[4*10+9]+=QC[0]*fsds[0][4*6+4] + WQ[0]*fsds[1][4*6+4];
	    FSFS[5*10+0]+=QC[0]*fsds[0][5*6+0] + WQ[0]*fsds[1][5*6+0]
			 + eta*tmp.fsps[5*3+0] + ze2*dsds[1][2*6+0];
	    FSFS[5*10+1]+=QC[1]*fsds[0][5*6+1] + WQ[1]*fsds[1][5*6+1]
			 + eta*tmp.fsps[5*3+1];
	    FSFS[5*10+2]+=QC[2]*fsds[0][5*6+2] + WQ[2]*fsds[1][5*6+2]
			 + eta*tmp.fsps[5*3+2] + ze22*dsds[1][5*6+2];
	    FSFS[5*10+3]+=QC[0]*fsds[0][5*6+3] + WQ[0]*fsds[1][5*6+3]
			 + eta2*tmp.fsps[5*3+1] + ze2*dsds[1][2*6+3];
	    FSFS[5*10+4]+=QC[1]*fsds[0][5*6+4] + WQ[1]*fsds[1][5*6+4]
			 + eta2*tmp.fsps[5*3+2];
	    FSFS[5*10+5]+=QC[2]*fsds[0][5*6+5] + WQ[2]*fsds[1][5*6+5]
			 + eta2*tmp.fsps[5*3+0] + ze22*dsds[1][5*6+5];
	    FSFS[5*10+6]+=QC[0]*fsds[0][5*6+1] + WQ[0]*fsds[1][5*6+1]
			 + ze2*dsds[1][2*6+1];
	    FSFS[5*10+7]+=QC[1]*fsds[0][5*6+2] + WQ[1]*fsds[1][5*6+2];
	    FSFS[5*10+8]+=QC[2]*fsds[0][5*6+0] + WQ[2]*fsds[1][5*6+0]
			 + ze22*dsds[1][5*6+0];
	    FSFS[5*10+9]+=QC[0]*fsds[0][5*6+4] + WQ[0]*fsds[1][5*6+4]
			 + ze2*dsds[1][2*6+4];
	    FSFS[6*10+0]+=QC[0]*fsds[0][6*6+0] + WQ[0]*fsds[1][6*6+0]
			 + eta*tmp.fsps[6*3+0] + ze2*dsds[1][1*6+0];
	    FSFS[6*10+1]+=QC[1]*fsds[0][6*6+1] + WQ[1]*fsds[1][6*6+1]
			 + eta*tmp.fsps[6*3+1] + ze22*dsds[1][3*6+1];
	    FSFS[6*10+2]+=QC[2]*fsds[0][6*6+2] + WQ[2]*fsds[1][6*6+2]
			 + eta*tmp.fsps[6*3+2];
	    FSFS[6*10+3]+=QC[0]*fsds[0][6*6+3] + WQ[0]*fsds[1][6*6+3]
			 + eta2*tmp.fsps[6*3+1] + ze2*dsds[1][1*6+3];
	    FSFS[6*10+4]+=QC[1]*fsds[0][6*6+4] + WQ[1]*fsds[1][6*6+4]
			 + eta2*tmp.fsps[6*3+2] + ze22*dsds[1][3*6+4];
	    FSFS[6*10+5]+=QC[2]*fsds[0][6*6+5] + WQ[2]*fsds[1][6*6+5]
			 + eta2*tmp.fsps[6*3+0];
	    FSFS[6*10+6]+=QC[0]*fsds[0][6*6+1] + WQ[0]*fsds[1][6*6+1]
			 + ze2*dsds[1][1*6+1];
	    FSFS[6*10+7]+=QC[1]*fsds[0][6*6+2] + WQ[1]*fsds[1][6*6+2]
			 + ze22*dsds[1][3*6+2];
	    FSFS[6*10+8]+=QC[2]*fsds[0][6*6+0] + WQ[2]*fsds[1][6*6+0];
	    FSFS[6*10+9]+=QC[0]*fsds[0][6*6+4] + WQ[0]*fsds[1][6*6+4]
			 + ze2*dsds[1][1*6+4];
	    FSFS[7*10+0]+=QC[0]*fsds[0][7*6+0] + WQ[0]*fsds[1][7*6+0]
			 + eta*tmp.fsps[7*3+0];
	    FSFS[7*10+1]+=QC[1]*fsds[0][7*6+1] + WQ[1]*fsds[1][7*6+1]
			 + eta*tmp.fsps[7*3+1] + ze2*dsds[1][2*6+1];
	    FSFS[7*10+2]+=QC[2]*fsds[0][7*6+2] + WQ[2]*fsds[1][7*6+2]
			 + eta*tmp.fsps[7*3+2] + ze22*dsds[1][4*6+2];
	    FSFS[7*10+3]+=QC[0]*fsds[0][7*6+3] + WQ[0]*fsds[1][7*6+3]
			 + eta2*tmp.fsps[7*3+1];
	    FSFS[7*10+4]+=QC[1]*fsds[0][7*6+4] + WQ[1]*fsds[1][7*6+4]
			 + eta2*tmp.fsps[7*3+2] + ze2*dsds[1][2*6+4];
	    FSFS[7*10+5]+=QC[2]*fsds[0][7*6+5] + WQ[2]*fsds[1][7*6+5]
			 + eta2*tmp.fsps[7*3+0] + ze22*dsds[1][4*6+5];
	    FSFS[7*10+6]+=QC[0]*fsds[0][7*6+1] + WQ[0]*fsds[1][7*6+1];
	    FSFS[7*10+7]+=QC[1]*fsds[0][7*6+2] + WQ[1]*fsds[1][7*6+2]
			 + ze2*dsds[1][2*6+2];
	    FSFS[7*10+8]+=QC[2]*fsds[0][7*6+0] + WQ[2]*fsds[1][7*6+0]
			 + ze22*dsds[1][4*6+0];
	    FSFS[7*10+9]+=QC[0]*fsds[0][7*6+4] + WQ[0]*fsds[1][7*6+4];
	    FSFS[8*10+0]+=QC[0]*fsds[0][8*6+0] + WQ[0]*fsds[1][8*6+0]
			 + eta*tmp.fsps[8*3+0] + ze22*dsds[1][5*6+0];
	    FSFS[8*10+1]+=QC[1]*fsds[0][8*6+1] + WQ[1]*fsds[1][8*6+1]
			 + eta*tmp.fsps[8*3+1];
	    FSFS[8*10+2]+=QC[2]*fsds[0][8*6+2] + WQ[2]*fsds[1][8*6+2]
			 + eta*tmp.fsps[8*3+2] + ze2*dsds[1][0*6+2];
	    FSFS[8*10+3]+=QC[0]*fsds[0][8*6+3] + WQ[0]*fsds[1][8*6+3]
			 + eta2*tmp.fsps[8*3+1] + ze22*dsds[1][5*6+3];
	    FSFS[8*10+4]+=QC[1]*fsds[0][8*6+4] + WQ[1]*fsds[1][8*6+4]
			 + eta2*tmp.fsps[8*3+2];
	    FSFS[8*10+5]+=QC[2]*fsds[0][8*6+5] + WQ[2]*fsds[1][8*6+5]
			 + eta2*tmp.fsps[8*3+0] + ze2*dsds[1][0*6+5];
	    FSFS[8*10+6]+=QC[0]*fsds[0][8*6+1] + WQ[0]*fsds[1][8*6+1]
			 + ze22*dsds[1][5*6+1];
	    FSFS[8*10+7]+=QC[1]*fsds[0][8*6+2] + WQ[1]*fsds[1][8*6+2];
	    FSFS[8*10+8]+=QC[2]*fsds[0][8*6+0] + WQ[2]*fsds[1][8*6+0]
			 + ze2*dsds[1][0*6+0];
	    FSFS[8*10+9]+=QC[0]*fsds[0][8*6+4] + WQ[0]*fsds[1][8*6+4]
			 + ze22*dsds[1][5*6+4];
	    FSFS[9*10+0]+=QC[0]*fsds[0][9*6+0] + WQ[0]*fsds[1][9*6+0]
			 + eta*tmp.fsps[9*3+0] + ze2*dsds[1][4*6+0];
	    FSFS[9*10+1]+=QC[1]*fsds[0][9*6+1] + WQ[1]*fsds[1][9*6+1]
			 + eta*tmp.fsps[9*3+1] + ze2*dsds[1][5*6+1];
	    FSFS[9*10+2]+=QC[2]*fsds[0][9*6+2] + WQ[2]*fsds[1][9*6+2]
			 + eta*tmp.fsps[9*3+2] + ze2*dsds[1][3*6+2];
	    FSFS[9*10+3]+=QC[0]*fsds[0][9*6+3] + WQ[0]*fsds[1][9*6+3]
			 + eta2*tmp.fsps[9*3+1] + ze2*dsds[1][4*6+3];
	    FSFS[9*10+4]+=QC[1]*fsds[0][9*6+4] + WQ[1]*fsds[1][9*6+4]
			 + eta2*tmp.fsps[9*3+2] + ze2*dsds[1][5*6+4];
	    FSFS[9*10+5]+=QC[2]*fsds[0][9*6+5] + WQ[2]*fsds[1][9*6+5]
			 + eta2*tmp.fsps[9*3+0] + ze2*dsds[1][3*6+5];
	    FSFS[9*10+6]+=QC[0]*fsds[0][9*6+1] + WQ[0]*fsds[1][9*6+1]
			 + ze2*dsds[1][4*6+1];
	    FSFS[9*10+7]+=QC[1]*fsds[0][9*6+2] + WQ[1]*fsds[1][9*6+2]
			 + ze2*dsds[1][5*6+2];
	    FSFS[9*10+8]+=QC[2]*fsds[0][9*6+0] + WQ[2]*fsds[1][9*6+0]
			 + ze2*dsds[1][3*6+0];
	    FSFS[9*10+9]+=QC[0]*fsds[0][9*6+4] + WQ[0]*fsds[1][9*6+4]
			 + ze2*dsds[1][4*6+4];
	}	// klps
    }	// ijps
    //
    // dpds
    for (c0=0; c0<6; c0++) {
        DPDS[0*3*6+0*6+c0] = FSDS[0*6+c0] - BA[0]*DSDS[0*6+c0];
        DPDS[0*3*6+1*6+c0] = FSDS[3*6+c0] - BA[1]*DSDS[0*6+c0];
        DPDS[0*3*6+2*6+c0] = FSDS[8*6+c0] - BA[2]*DSDS[0*6+c0];
        DPDS[1*3*6+0*6+c0] = FSDS[6*6+c0] - BA[0]*DSDS[1*6+c0];
        DPDS[1*3*6+1*6+c0] = FSDS[1*6+c0] - BA[1]*DSDS[1*6+c0];
        DPDS[1*3*6+2*6+c0] = FSDS[4*6+c0] - BA[2]*DSDS[1*6+c0];
        DPDS[2*3*6+0*6+c0] = FSDS[5*6+c0] - BA[0]*DSDS[2*6+c0];
        DPDS[2*3*6+1*6+c0] = FSDS[7*6+c0] - BA[1]*DSDS[2*6+c0];
        DPDS[2*3*6+2*6+c0] = FSDS[2*6+c0] - BA[2]*DSDS[2*6+c0];
        DPDS[3*3*6+0*6+c0] = FSDS[3*6+c0] - BA[0]*DSDS[3*6+c0];
        DPDS[3*3*6+1*6+c0] = FSDS[6*6+c0] - BA[1]*DSDS[3*6+c0];
        DPDS[3*3*6+2*6+c0] = FSDS[9*6+c0] - BA[2]*DSDS[3*6+c0];
        DPDS[4*3*6+0*6+c0] = FSDS[9*6+c0] - BA[0]*DSDS[4*6+c0];
        DPDS[4*3*6+1*6+c0] = FSDS[4*6+c0] - BA[1]*DSDS[4*6+c0];
        DPDS[4*3*6+2*6+c0] = FSDS[7*6+c0] - BA[2]*DSDS[4*6+c0];
        DPDS[5*3*6+0*6+c0] = FSDS[8*6+c0] - BA[0]*DSDS[5*6+c0];
        DPDS[5*3*6+1*6+c0] = FSDS[9*6+c0] - BA[1]*DSDS[5*6+c0];
        DPDS[5*3*6+2*6+c0] = FSDS[5*6+c0] - BA[2]*DSDS[5*6+c0];
    }
    // dpfs
    for (c0=0; c0<10; c0++) {
        DPFS[0*3*10+0*10+c0] = FSFS[0*10+c0] - BA[0]*DSFS[0*10+c0];
        DPFS[0*3*10+1*10+c0] = FSFS[3*10+c0] - BA[1]*DSFS[0*10+c0];
        DPFS[0*3*10+2*10+c0] = FSFS[8*10+c0] - BA[2]*DSFS[0*10+c0];
        DPFS[1*3*10+0*10+c0] = FSFS[6*10+c0] - BA[0]*DSFS[1*10+c0];
        DPFS[1*3*10+1*10+c0] = FSFS[1*10+c0] - BA[1]*DSFS[1*10+c0];
        DPFS[1*3*10+2*10+c0] = FSFS[4*10+c0] - BA[2]*DSFS[1*10+c0];
        DPFS[2*3*10+0*10+c0] = FSFS[5*10+c0] - BA[0]*DSFS[2*10+c0];
        DPFS[2*3*10+1*10+c0] = FSFS[7*10+c0] - BA[1]*DSFS[2*10+c0];
        DPFS[2*3*10+2*10+c0] = FSFS[2*10+c0] - BA[2]*DSFS[2*10+c0];
        DPFS[3*3*10+0*10+c0] = FSFS[3*10+c0] - BA[0]*DSFS[3*10+c0];
        DPFS[3*3*10+1*10+c0] = FSFS[6*10+c0] - BA[1]*DSFS[3*10+c0];
        DPFS[3*3*10+2*10+c0] = FSFS[9*10+c0] - BA[2]*DSFS[3*10+c0];
        DPFS[4*3*10+0*10+c0] = FSFS[9*10+c0] - BA[0]*DSFS[4*10+c0];
        DPFS[4*3*10+1*10+c0] = FSFS[4*10+c0] - BA[1]*DSFS[4*10+c0];
        DPFS[4*3*10+2*10+c0] = FSFS[7*10+c0] - BA[2]*DSFS[4*10+c0];
        DPFS[5*3*10+0*10+c0] = FSFS[8*10+c0] - BA[0]*DSFS[5*10+c0];
        DPFS[5*3*10+1*10+c0] = FSFS[9*10+c0] - BA[1]*DSFS[5*10+c0];
        DPFS[5*3*10+2*10+c0] = FSFS[5*10+c0] - BA[2]*DSFS[5*10+c0];
    }
    // dpdp
    for (ab=0, ab01=0, ab10=0, ab00=0;
            ab<6*3; ab++, ab01 += 6*3, ab10 += 10, ab00 += 6) {
        DPDP[ab01+0*3+0] = DPFS[ab10+0] - DC[0]*DPDS[ab00+0];
        DPDP[ab01+0*3+1] = DPFS[ab10+3] - DC[1]*DPDS[ab00+0];
        DPDP[ab01+0*3+2] = DPFS[ab10+8] - DC[2]*DPDS[ab00+0];
        DPDP[ab01+1*3+0] = DPFS[ab10+6] - DC[0]*DPDS[ab00+1];
        DPDP[ab01+1*3+1] = DPFS[ab10+1] - DC[1]*DPDS[ab00+1];
        DPDP[ab01+1*3+2] = DPFS[ab10+4] - DC[2]*DPDS[ab00+1];
        DPDP[ab01+2*3+0] = DPFS[ab10+5] - DC[0]*DPDS[ab00+2];
        DPDP[ab01+2*3+1] = DPFS[ab10+7] - DC[1]*DPDS[ab00+2];
        DPDP[ab01+2*3+2] = DPFS[ab10+2] - DC[2]*DPDS[ab00+2];
        DPDP[ab01+3*3+0] = DPFS[ab10+3] - DC[0]*DPDS[ab00+3];
        DPDP[ab01+3*3+1] = DPFS[ab10+6] - DC[1]*DPDS[ab00+3];
        DPDP[ab01+3*3+2] = DPFS[ab10+9] - DC[2]*DPDS[ab00+3];
        DPDP[ab01+4*3+0] = DPFS[ab10+9] - DC[0]*DPDS[ab00+4];
        DPDP[ab01+4*3+1] = DPFS[ab10+4] - DC[1]*DPDS[ab00+4];
        DPDP[ab01+4*3+2] = DPFS[ab10+7] - DC[2]*DPDS[ab00+4];
        DPDP[ab01+5*3+0] = DPFS[ab10+8] - DC[0]*DPDS[ab00+5];
        DPDP[ab01+5*3+1] = DPFS[ab10+9] - DC[1]*DPDS[ab00+5];
        DPDP[ab01+5*3+2] = DPFS[ab10+5] - DC[2]*DPDS[ab00+5];
    }
    //
    for ( i=0, ix=0; i<6; i++ ) {
	coe0 = (i<3? ONE : sqr3 );
	for ( j=0; j<3; j++ ) {
	    for ( k=0; k<6; k++ ) {
		coe = coe0 * (k<3? ONE : sqr3 );
		for ( l=0; l<3; l++, ix++ ) DPDP[ix] *= coe;
	    }
	}
    }
}
