/**
 * @file ofmo-twoint-core-dddd.c
 * １つのCS４重対に対する２電子積分を計算する関数群。
 * 2011/06/16現在、(ss,ss)～(dd,dd)までの２１種類の２電子積分
 * 計算をサポートしている。
 * このファイルには、(dd,dd)タイプの２電子積分計算を
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

extern double *FMT_fmt_table8;

extern double FMT_fmt_step_size;
extern double FMT_fmt_inv_step_size;
extern double FMT_pi_div2;

/** １つのCS４重対に対して(dd,dd)タイプの２電子積分を計算する関数
 * @ingroup core-twoint
 * */
void twoint_core_dddd__(
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double DDDD[6*6*6*6] ) {
    int ijps, klps, i, j, k, l, ix, m, m1;
    int ab, ab01, ab10, ab00, c0;
    double cssss, zeta, dkab, xiza, eta, xizc, dk;
    double rz, re, ze2, ze22, ze23, ze24, zeta2, eta2, zeta23, eta23;
    union _temp_ {
	double entity[15*6];
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
    double sqrho, rho, PQ2, T;
    double PC[3], PA[3], WP[3], QP[3], QC[3], WQ[3];
    double ssss[8+1], psss[7+1][3], dsss[6+1][6], fsss[5+1][10];
    double gsss[4+1][15];
    double ssps[3-1][3], psps[3+0][3*3], dsps[3+1][6*3], fsps[3+1][10*3];
    double gsps[3+1][15*3];
    double ssds[6], psds[2+0][3*6], dsds[2+1][6*6], fsds[2+1][10*6];
    double gsds[2+1][15*6];
    double psfs[3*10], dsfs[1+1][6*10], fsfs[1+1][10*10], gsfs[1+1][15*10];
    double DSDS[6* 6], FSDS[10* 6], GSDS[15* 6];
    double DSFS[6*10], FSFS[10*10], GSFS[15*10];
    double DSGS[6*15], FSGS[10*15], GSGS[15*15];
    double DDDS[ 6*6*6], DDFS[ 6*6*10], DDGS[ 6*6*15];
    double FPDS[10*3*6], FPFS[10*3*10], FPGS[10*3*15];
    double DPDS[15*3*6], DPFS[15*3*10], DPGS[15*3*15];
    double DDDP[6*6*6*3], DDFP[6*6*10*3];
    double sqr3, coe0, coe1, coe2, coe;
    sqr3 = sqrt(3.e0);
    for ( i=0; i< 6* 6; i++ ) DSDS[i] = ZERO;
    for ( i=0; i<10* 6; i++ ) FSDS[i] = ZERO;
    for ( i=0; i<15* 6; i++ ) GSDS[i] = ZERO;
    for ( i=0; i< 6*10; i++ ) DSFS[i] = ZERO;
    for ( i=0; i<10*10; i++ ) FSFS[i] = ZERO;
    for ( i=0; i<15*10; i++ ) GSFS[i] = ZERO;
    for ( i=0; i< 6*15; i++ ) DSGS[i] = ZERO;
    for ( i=0; i<10*15; i++ ) FSGS[i] = ZERO;
    for ( i=0; i<15*15; i++ ) GSGS[i] = ZERO;
    for ( ijps=0; ijps<(*nijps); ijps++ ) {
	zeta = vzeta[ijps];
	dkab = vdkab[ijps];
	xiza = vxiza[ijps];
	zeta2 = HALF * zeta;
	zeta23 = zeta + zeta2;
	for ( i=0; i<3; i++ ) {
	    PC[i] = AC[i] + xiza*BA[i];
	    PA[i] = xiza * BA[i];
	}
	for ( klps=0; klps<(*nklps); klps++ ) {
	    eta  = veta[klps];
	    dk   = dkab * vdkcd[klps];
	    xizc = vxizc[klps];
	    eta2 = HALF * eta;
	    eta23 = 1.5e0 * eta;
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
	    ze23  = 3.e0 * ze2;
	    ze24  = ze22 * 2.e0;
	    for ( i=0; i<3; i++ ) {
		WP[i] = rz*QP[i];
		WQ[i] = rz*QP[i] - QP[i];
	    }
	    T     = rho * PQ2;
	    cssss = sqrho * dk;
	    {
		int it0, pos;
		double dT, t_inv, st_inv, dT2, dT3;
		if ( T < 52.e0 ) {
		    it0 = (int)(0.5e0 + T * FMT_fmt_inv_step_size);
		    dT  = it0 * FMT_fmt_step_size - T;
		    dT2 = dT * _twoint_inv2_;
		    dT3 = dT * _twoint_inv3_;
		    pos = it0 * (8+4);
		    ssss[0] = cssss*(((FMT_fmt_table8[pos+3]   * dT3
				     + FMT_fmt_table8[pos+2] ) * dT2
				     + FMT_fmt_table8[pos+1] ) * dT
				     + FMT_fmt_table8[pos+0] );
		    ssss[1] = cssss*(((FMT_fmt_table8[pos+4]   * dT3
				     + FMT_fmt_table8[pos+3] ) * dT2
				     + FMT_fmt_table8[pos+2] ) * dT
				     + FMT_fmt_table8[pos+1] );
		    ssss[2] = cssss*(((FMT_fmt_table8[pos+5]   * dT3
				     + FMT_fmt_table8[pos+4] ) * dT2
				     + FMT_fmt_table8[pos+3] ) * dT
				     + FMT_fmt_table8[pos+2] );
		    ssss[3] = cssss*(((FMT_fmt_table8[pos+6]   * dT3
				     + FMT_fmt_table8[pos+5] ) * dT2
				     + FMT_fmt_table8[pos+4] ) * dT
				     + FMT_fmt_table8[pos+3] );
		    ssss[4] = cssss*(((FMT_fmt_table8[pos+7]   * dT3
				     + FMT_fmt_table8[pos+6] ) * dT2
				     + FMT_fmt_table8[pos+5] ) * dT
				     + FMT_fmt_table8[pos+4] );
		    ssss[5] = cssss*(((FMT_fmt_table8[pos+8]   * dT3
				     + FMT_fmt_table8[pos+7] ) * dT2
				     + FMT_fmt_table8[pos+6] ) * dT
				     + FMT_fmt_table8[pos+5] );
		    ssss[6] = cssss*(((FMT_fmt_table8[pos+9]   * dT3
				     + FMT_fmt_table8[pos+8] ) * dT2
				     + FMT_fmt_table8[pos+7] ) * dT
				     + FMT_fmt_table8[pos+6] );
		    ssss[7] = cssss*(((FMT_fmt_table8[pos+10]   * dT3
				     + FMT_fmt_table8[pos+ 9] ) * dT2
				     + FMT_fmt_table8[pos+ 8] ) * dT
				     + FMT_fmt_table8[pos+ 7] );
		    ssss[8] = cssss*(((FMT_fmt_table8[pos+11]   * dT3
				     + FMT_fmt_table8[pos+10] ) * dT2
				     + FMT_fmt_table8[pos+ 9] ) * dT
				     + FMT_fmt_table8[pos+ 8] );
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
		    ssss[7] = 13.e0 * t_inv * ssss[6];
		    ssss[8] = 15.e0 * t_inv * ssss[7];
		}
	    }
	    //fmt( ssss, 8, T, cssss );
	    // psss (m=0,7)
	    for (m=0, m1=1; m<=7; m++, m1++) {
		psss[m][0] = PA[0]*ssss[m] + WP[0]*ssss[m1];
		psss[m][1] = PA[1]*ssss[m] + WP[1]*ssss[m1];
		psss[m][2] = PA[2]*ssss[m] + WP[2]*ssss[m1];
	    }
	    // dsss (m=0,6)
	    for (m=0, m1=1; m<=6; m++, m1++) {
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
	    }
	    // fsss (m=0,5)
	    for (m=0, m1=1; m<=5; m++, m1++) {
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
	    }
	    // gsss (m=0,4)
	    for (m=0, m1=1; m<=4; m++, m1++) {
		for (i=0; i<6; i++) tmp.dsss[i]=dsss[m][i]-rz*dsss[m1][i];
		gsss[m][ 0] = PA[0]*fsss[m][0] + WP[0]*fsss[m1][0]
			    + zeta23*tmp.dsss[0];
		gsss[m][ 1] = PA[1]*fsss[m][1] + WP[1]*fsss[m1][1]
			    + zeta23*tmp.dsss[1];
		gsss[m][ 2] = PA[2]*fsss[m][2] + WP[2]*fsss[m1][2]
			    + zeta23*tmp.dsss[2];
		gsss[m][ 3] = PA[0]*fsss[m][3] + WP[0]*fsss[m1][3]
			    + zeta*tmp.dsss[3];
		gsss[m][ 4] = PA[1]*fsss[m][4] + WP[1]*fsss[m1][4]
			    + zeta*tmp.dsss[4];
		gsss[m][ 5] = PA[2]*fsss[m][5] + WP[2]*fsss[m1][5]
			    + zeta*tmp.dsss[5];
		gsss[m][ 6] = PA[0]*fsss[m][6] + WP[0]*fsss[m1][6]
			    + zeta2*tmp.dsss[1];
		gsss[m][ 7] = PA[1]*fsss[m][7] + WP[1]*fsss[m1][7]
			    + zeta2*tmp.dsss[2];
		gsss[m][ 8] = PA[2]*fsss[m][8] + WP[2]*fsss[m1][8]
			    + zeta2*tmp.dsss[0];
		gsss[m][ 9] = PA[0]*fsss[m][1] + WP[0]*fsss[m1][1];
		gsss[m][10] = PA[1]*fsss[m][2] + WP[1]*fsss[m1][2];
		gsss[m][11] = PA[2]*fsss[m][0] + WP[2]*fsss[m1][0];
		gsss[m][12] = PA[0]*fsss[m][9] + WP[0]*fsss[m1][9]
			    + zeta2*tmp.dsss[4];
		gsss[m][13] = PA[1]*fsss[m][9] + WP[1]*fsss[m1][9]
			    + zeta2*tmp.dsss[5];
		gsss[m][14] = PA[2]*fsss[m][9] + WP[2]*fsss[m1][9]
			    + zeta2*tmp.dsss[3];
	    }
	    // ssps (m=2,3)
	    for (m=2, m1=3; m<=3; m++, m1++) {
		ssps[m-2][0] = QC[0]*ssss[m] + WQ[0]*ssss[m1];
		ssps[m-2][1] = QC[1]*ssss[m] + WQ[1]*ssss[m1];
		ssps[m-2][2] = QC[2]*ssss[m] + WQ[2]*ssss[m1];
	    }	// end ssps loop
	    // psps (m=1,3)
	    for (m=1, m1=2; m<=3; m++, m1++) {
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
	    // dsps (m=0,3)
	    for (m=0, m1=1; m<=3; m++, m1++) {
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
	    // fsps (m=0,3)
	    for (m=0, m1=1; m<=3; m++, m1++) {
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
	    // gsps (m=0,3)
	    for (m=0, m1=1; m<=3; m++, m1++) {
		gsps[m][ 0*3+0] = QC[0]*gsss[m][ 0] + WQ[0]*gsss[m1][ 0]
				+ ze24*fsss[m1][0];
		gsps[m][ 0*3+1] = QC[1]*gsss[m][ 0] + WQ[1]*gsss[m1][ 0];
		gsps[m][ 0*3+2] = QC[2]*gsss[m][ 0] + WQ[2]*gsss[m1][ 0];
		gsps[m][ 1*3+0] = QC[0]*gsss[m][ 1] + WQ[0]*gsss[m1][ 1];
		gsps[m][ 1*3+1] = QC[1]*gsss[m][ 1] + WQ[1]*gsss[m1][ 1]
				+ ze24*fsss[m1][1];
		gsps[m][ 1*3+2] = QC[2]*gsss[m][ 1] + WQ[2]*gsss[m1][ 1];
		gsps[m][ 2*3+0] = QC[0]*gsss[m][ 2] + WQ[0]*gsss[m1][ 2];
		gsps[m][ 2*3+1] = QC[1]*gsss[m][ 2] + WQ[1]*gsss[m1][ 2];
		gsps[m][ 2*3+2] = QC[2]*gsss[m][ 2] + WQ[2]*gsss[m1][ 2]
				+ ze24*fsss[m1][2];
		gsps[m][ 3*3+0] = QC[0]*gsss[m][ 3] + WQ[0]*gsss[m1][ 3]
				+ ze23*fsss[m1][3];
		gsps[m][ 3*3+1] = QC[1]*gsss[m][ 3] + WQ[1]*gsss[m1][ 3]
				+ ze2*fsss[m1][0];
		gsps[m][ 3*3+2] = QC[2]*gsss[m][ 3] + WQ[2]*gsss[m1][ 3];
		gsps[m][ 4*3+0] = QC[0]*gsss[m][ 4] + WQ[0]*gsss[m1][ 4];
		gsps[m][ 4*3+1] = QC[1]*gsss[m][ 4] + WQ[1]*gsss[m1][ 4]
				+ ze23*fsss[m1][4];
		gsps[m][ 4*3+2] = QC[2]*gsss[m][ 4] + WQ[2]*gsss[m1][ 4]
				+ ze2*fsss[m1][1];
		gsps[m][ 5*3+0] = QC[0]*gsss[m][ 5] + WQ[0]*gsss[m1][ 5]
				+ ze2*fsss[m1][2];
		gsps[m][ 5*3+1] = QC[1]*gsss[m][ 5] + WQ[1]*gsss[m1][ 5];
		gsps[m][ 5*3+2] = QC[2]*gsss[m][ 5] + WQ[2]*gsss[m1][ 5]
				+ ze23*fsss[m1][5];
		gsps[m][ 6*3+0] = QC[0]*gsss[m][ 6] + WQ[0]*gsss[m1][ 6]
				+ ze22*fsss[m1][6];
		gsps[m][ 6*3+1] = QC[1]*gsss[m][ 6] + WQ[1]*gsss[m1][ 6]
				+ ze22*fsss[m1][3];
		gsps[m][ 6*3+2] = QC[2]*gsss[m][ 6] + WQ[2]*gsss[m1][ 6];
		gsps[m][ 7*3+0] = QC[0]*gsss[m][ 7] + WQ[0]*gsss[m1][ 7];
		gsps[m][ 7*3+1] = QC[1]*gsss[m][ 7] + WQ[1]*gsss[m1][ 7]
				+ ze22*fsss[m1][7];
		gsps[m][ 7*3+2] = QC[2]*gsss[m][ 7] + WQ[2]*gsss[m1][ 7]
				+ ze22*fsss[m1][4];
		gsps[m][ 8*3+0] = QC[0]*gsss[m][ 8] + WQ[0]*gsss[m1][ 8]
				+ ze22*fsss[m1][5];
		gsps[m][ 8*3+1] = QC[1]*gsss[m][ 8] + WQ[1]*gsss[m1][ 8];
		gsps[m][ 8*3+2] = QC[2]*gsss[m][ 8] + WQ[2]*gsss[m1][ 8]
				+ ze22*fsss[m1][8];
		gsps[m][ 9*3+0] = QC[0]*gsss[m][ 9] + WQ[0]*gsss[m1][ 9]
				+ ze2*fsss[m1][1];
		gsps[m][ 9*3+1] = QC[1]*gsss[m][ 9] + WQ[1]*gsss[m1][ 9]
				+ ze23*fsss[m1][6];
		gsps[m][ 9*3+2] = QC[2]*gsss[m][ 9] + WQ[2]*gsss[m1][ 9];
		gsps[m][10*3+0] = QC[0]*gsss[m][10] + WQ[0]*gsss[m1][10];
		gsps[m][10*3+1] = QC[1]*gsss[m][10] + WQ[1]*gsss[m1][10]
				+ ze2*fsss[m1][2];
		gsps[m][10*3+2] = QC[2]*gsss[m][10] + WQ[2]*gsss[m1][10]
				+ ze23*fsss[m1][7];
		gsps[m][11*3+0] = QC[0]*gsss[m][11] + WQ[0]*gsss[m1][11]
				+ ze23*fsss[m1][8];
		gsps[m][11*3+1] = QC[1]*gsss[m][11] + WQ[1]*gsss[m1][11];
		gsps[m][11*3+2] = QC[2]*gsss[m][11] + WQ[2]*gsss[m1][11]
				+ ze2*fsss[m1][0];
		gsps[m][12*3+0] = QC[0]*gsss[m][12] + WQ[0]*gsss[m1][12]
				+ ze22*fsss[m1][9];
		gsps[m][12*3+1] = QC[1]*gsss[m][12] + WQ[1]*gsss[m1][12]
				+ ze2*fsss[m1][8];
		gsps[m][12*3+2] = QC[2]*gsss[m][12] + WQ[2]*gsss[m1][12]
				+ ze2*fsss[m1][3];
		gsps[m][13*3+0] = QC[0]*gsss[m][13] + WQ[0]*gsss[m1][13]
				+ ze2*fsss[m1][4];
		gsps[m][13*3+1] = QC[1]*gsss[m][13] + WQ[1]*gsss[m1][13]
				+ ze22*fsss[m1][9];
		gsps[m][13*3+2] = QC[2]*gsss[m][13] + WQ[2]*gsss[m1][13]
				+ ze2*fsss[m1][6];
		gsps[m][14*3+0] = QC[0]*gsss[m][14] + WQ[0]*gsss[m1][14]
				+ ze2*fsss[m1][7];
		gsps[m][14*3+1] = QC[1]*gsss[m][14] + WQ[1]*gsss[m1][14]
				+ ze2*fsss[m1][5];
		gsps[m][14*3+2] = QC[2]*gsss[m][14] + WQ[2]*gsss[m1][14]
				+ ze22*fsss[m1][9];
	    }	// end do gsps
	    // ssds (m=2)
	    tmp.ssss = ssss[2] - re*ssss[3];
	    ssds[0] = QC[0]*ssps[2-2][0] + WQ[0]*ssps[3-2][0]
		    + eta2*tmp.ssss;
	    ssds[1] = QC[1]*ssps[2-2][1] + WQ[1]*ssps[3-2][1]
		    + eta2*tmp.ssss;
	    ssds[2] = QC[2]*ssps[2-2][2] + WQ[2]*ssps[3-2][2]
		    + eta2*tmp.ssss;
	    ssds[3] = QC[0]*ssps[2-2][1] + WQ[0]*ssps[3-2][1];
	    ssds[4] = QC[1]*ssps[2-2][2] + WQ[1]*ssps[3-2][2];
	    ssds[5] = QC[2]*ssps[2-2][0] + WQ[2]*ssps[3-2][0];
	    // psds (m=1,2)
	    for (m=1, m1=2; m<=2; m++, m1++) {
		for (i=0; i<3; i++) tmp.psss[i]=psss[m][i]-re*psss[m1][i];
		psds[m-1][0*6+0] = QC[0]*psps[m-1][0*3+0]
		    + WQ[0]*psps[m1-1][0*3+0] + eta2*tmp.psss[0]
		    + ze2*ssps[m1-2][0];
		psds[m-1][0*6+1] = QC[1]*psps[m-1][0*3+1]
		    + WQ[1]*psps[m1-1][0*3+1] + eta2*tmp.psss[0];
		psds[m-1][0*6+2] = QC[2]*psps[m-1][0*3+2]
		    + WQ[2]*psps[m1-1][0*3+2] + eta2*tmp.psss[0];
		psds[m-1][0*6+3] = QC[0]*psps[m-1][0*3+1]
		    + WQ[0]*psps[m1-1][0*3+1] + ze2*ssps[m1-2][1];
		psds[m-1][0*6+4] = QC[1]*psps[m-1][0*3+2]
		    + WQ[1]*psps[m1-1][0*3+2];
		psds[m-1][0*6+5] = QC[2]*psps[m-1][0*3+0]
		    + WQ[2]*psps[m1-1][0*3+0];
		psds[m-1][1*6+0] = QC[0]*psps[m-1][1*3+0]
		    + WQ[0]*psps[m1-1][1*3+0] + eta2*tmp.psss[1];
		psds[m-1][1*6+1] = QC[1]*psps[m-1][1*3+1]
		    + WQ[1]*psps[m1-1][1*3+1] + eta2*tmp.psss[1]
		    + ze2*ssps[m1-2][1];
		psds[m-1][1*6+2] = QC[2]*psps[m-1][1*3+2]
		    + WQ[2]*psps[m1-1][1*3+2] + eta2*tmp.psss[1];
		psds[m-1][1*6+3] = QC[0]*psps[m-1][1*3+1]
		    + WQ[0]*psps[m1-1][1*3+1];
		psds[m-1][1*6+4] = QC[1]*psps[m-1][1*3+2]
		    + WQ[1]*psps[m1-1][1*3+2] + ze2*ssps[m1-2][2];
		psds[m-1][1*6+5] = QC[2]*psps[m-1][1*3+0]
		    + WQ[2]*psps[m1-1][1*3+0];
		psds[m-1][2*6+0] = QC[0]*psps[m-1][2*3+0]
		    + WQ[0]*psps[m1-1][2*3+0] + eta2*tmp.psss[2];
		psds[m-1][2*6+1] = QC[1]*psps[m-1][2*3+1]
		    + WQ[1]*psps[m1-1][2*3+1] + eta2*tmp.psss[2];
		psds[m-1][2*6+2] = QC[2]*psps[m-1][2*3+2]
		    + WQ[2]*psps[m1-1][2*3+2] + eta2*tmp.psss[2]
		    + ze2*ssps[m1-2][2];
		psds[m-1][2*6+3] = QC[0]*psps[m-1][2*3+1]
		    + WQ[0]*psps[m1-1][2*3+1];
		psds[m-1][2*6+4] = QC[1]*psps[m-1][2*3+2]
		    + WQ[1]*psps[m1-1][2*3+2];
		psds[m-1][2*6+5] = QC[2]*psps[m-1][2*3+0]
		    + WQ[2]*psps[m1-1][2*3+0] + ze2*ssps[m1-2][0];
	    }	// end psds loop
	    // dsds (m=0,2)
	    for ( m=0, m1=1; m<=2; m++, m1++) {
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
	    }	// end dsds loop
	    for (i=0; i<6*6; i++) DSDS[i] += dsds[0][i];
	    // fsds (m=0,2)
	    for (m=0, m1=1; m<=2; m++, m1++) {
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
	    }	// end fsds loop
	    for (i=0; i<10*6; i++) FSDS[i] += fsds[0][i];
	    // gsds (m=0,2)
	    for (m=0, m1=1; m<=2; m++, m1++) {
		for (i=0; i<15; i++) tmp.gsss[i]=gsss[m][i]-re*gsss[m1][i];
		gsds[m][0*6+0]=QC[0]*gsps[m][0*3+0]+WQ[0]*gsps[m1][0*3+0]
				+ eta2*tmp.gsss[0] + ze24*fsps[m1][0*3+0];
		gsds[m][0*6+1]=QC[1]*gsps[m][0*3+1]+WQ[1]*gsps[m1][0*3+1]
				+ eta2*tmp.gsss[0];
		gsds[m][0*6+2]=QC[2]*gsps[m][0*3+2]+WQ[2]*gsps[m1][0*3+2]
				+ eta2*tmp.gsss[0];
		gsds[m][0*6+3]=QC[0]*gsps[m][0*3+1]+WQ[0]*gsps[m1][0*3+1]
				+ ze24*fsps[m1][0*3+1];
		gsds[m][0*6+4]=QC[1]*gsps[m][0*3+2]+WQ[1]*gsps[m1][0*3+2];
		gsds[m][0*6+5]=QC[2]*gsps[m][0*3+0]+WQ[2]*gsps[m1][0*3+0];
		gsds[m][1*6+0]=QC[0]*gsps[m][1*3+0]+WQ[0]*gsps[m1][1*3+0]
				+ eta2*tmp.gsss[1];
		gsds[m][1*6+1]=QC[1]*gsps[m][1*3+1]+WQ[1]*gsps[m1][1*3+1]
				+ eta2*tmp.gsss[1] + ze24*fsps[m1][1*3+1];
		gsds[m][1*6+2]=QC[2]*gsps[m][1*3+2]+WQ[2]*gsps[m1][1*3+2]
				+ eta2*tmp.gsss[1];
		gsds[m][1*6+3]=QC[0]*gsps[m][1*3+1]+WQ[0]*gsps[m1][1*3+1];
		gsds[m][1*6+4]=QC[1]*gsps[m][1*3+2]+WQ[1]*gsps[m1][1*3+2]
				+ ze24*fsps[m1][1*3+2];
		gsds[m][1*6+5]=QC[2]*gsps[m][1*3+0]+WQ[2]*gsps[m1][1*3+0];
		gsds[m][2*6+0]=QC[0]*gsps[m][2*3+0]+WQ[0]*gsps[m1][2*3+0]
				+ eta2*tmp.gsss[2];
		gsds[m][2*6+1]=QC[1]*gsps[m][2*3+1]+WQ[1]*gsps[m1][2*3+1]
				+ eta2*tmp.gsss[2];
		gsds[m][2*6+2]=QC[2]*gsps[m][2*3+2]+WQ[2]*gsps[m1][2*3+2]
				+ eta2*tmp.gsss[2] + ze24*fsps[m1][2*3+2];
		gsds[m][2*6+3]=QC[0]*gsps[m][2*3+1]+WQ[0]*gsps[m1][2*3+1];
		gsds[m][2*6+4]=QC[1]*gsps[m][2*3+2]+WQ[1]*gsps[m1][2*3+2];
		gsds[m][2*6+5]=QC[2]*gsps[m][2*3+0]+WQ[2]*gsps[m1][2*3+0]
				+ ze24*fsps[m1][2*3+0];
		gsds[m][3*6+0]=QC[0]*gsps[m][3*3+0]+WQ[0]*gsps[m1][3*3+0]
				+ eta2*tmp.gsss[3] + ze23*fsps[m1][3*3+0];
		gsds[m][3*6+1]=QC[1]*gsps[m][3*3+1]+WQ[1]*gsps[m1][3*3+1]
				+ eta2*tmp.gsss[3] + ze2*fsps[m1][0*3+1];
		gsds[m][3*6+2]=QC[2]*gsps[m][3*3+2]+WQ[2]*gsps[m1][3*3+2]
				+ eta2*tmp.gsss[3];
		gsds[m][3*6+3]=QC[0]*gsps[m][3*3+1]+WQ[0]*gsps[m1][3*3+1]
				+ ze23*fsps[m1][3*3+1];
		gsds[m][3*6+4]=QC[1]*gsps[m][3*3+2]+WQ[1]*gsps[m1][3*3+2]
				+ ze2*fsps[m1][0*3+2];
		gsds[m][3*6+5]=QC[2]*gsps[m][3*3+0]+WQ[2]*gsps[m1][3*3+0];
		gsds[m][4*6+0]=QC[0]*gsps[m][4*3+0]+WQ[0]*gsps[m1][4*3+0]
				+ eta2*tmp.gsss[4];
		gsds[m][4*6+1]=QC[1]*gsps[m][4*3+1]+WQ[1]*gsps[m1][4*3+1]
				+ eta2*tmp.gsss[4] + ze23*fsps[m1][4*3+1];
		gsds[m][4*6+2]=QC[2]*gsps[m][4*3+2]+WQ[2]*gsps[m1][4*3+2]
				+ eta2*tmp.gsss[4] + ze2*fsps[m1][1*3+2];
		gsds[m][4*6+3]=QC[0]*gsps[m][4*3+1]+WQ[0]*gsps[m1][4*3+1];
		gsds[m][4*6+4]=QC[1]*gsps[m][4*3+2]+WQ[1]*gsps[m1][4*3+2]
				+ ze23*fsps[m1][4*3+2];
		gsds[m][4*6+5]=QC[2]*gsps[m][4*3+0]+WQ[2]*gsps[m1][4*3+0]
				+ ze2*fsps[m1][1*3+0];
		gsds[m][5*6+0]=QC[0]*gsps[m][5*3+0]+WQ[0]*gsps[m1][5*3+0]
				+ eta2*tmp.gsss[5] + ze2*fsps[m1][2*3+0];
		gsds[m][5*6+1]=QC[1]*gsps[m][5*3+1]+WQ[1]*gsps[m1][5*3+1]
				+ eta2*tmp.gsss[5];
		gsds[m][5*6+2]=QC[2]*gsps[m][5*3+2]+WQ[2]*gsps[m1][5*3+2]
				+ eta2*tmp.gsss[5] + ze23*fsps[m1][5*3+2];
		gsds[m][5*6+3]=QC[0]*gsps[m][5*3+1]+WQ[0]*gsps[m1][5*3+1]
				+ ze2*fsps[m1][2*3+1];
		gsds[m][5*6+4]=QC[1]*gsps[m][5*3+2]+WQ[1]*gsps[m1][5*3+2];
		gsds[m][5*6+5]=QC[2]*gsps[m][5*3+0]+WQ[2]*gsps[m1][5*3+0]
				+ ze23*fsps[m1][5*3+0];
		gsds[m][6*6+0]=QC[0]*gsps[m][6*3+0]+WQ[0]*gsps[m1][6*3+0]
				+ eta2*tmp.gsss[6] + ze22*fsps[m1][6*3+0];
		gsds[m][6*6+1]=QC[1]*gsps[m][6*3+1]+WQ[1]*gsps[m1][6*3+1]
				+ eta2*tmp.gsss[6] + ze22*fsps[m1][3*3+1];
		gsds[m][6*6+2]=QC[2]*gsps[m][6*3+2]+WQ[2]*gsps[m1][6*3+2]
				+ eta2*tmp.gsss[6];
		gsds[m][6*6+3]=QC[0]*gsps[m][6*3+1]+WQ[0]*gsps[m1][6*3+1]
				+ ze22*fsps[m1][6*3+1];
		gsds[m][6*6+4]=QC[1]*gsps[m][6*3+2]+WQ[1]*gsps[m1][6*3+2]
				+ ze22*fsps[m1][3*3+2];
		gsds[m][6*6+5]=QC[2]*gsps[m][6*3+0]+WQ[2]*gsps[m1][6*3+0];
		gsds[m][7*6+0]=QC[0]*gsps[m][7*3+0]+WQ[0]*gsps[m1][7*3+0]
				+ eta2*tmp.gsss[7];
		gsds[m][7*6+1]=QC[1]*gsps[m][7*3+1]+WQ[1]*gsps[m1][7*3+1]
				+ eta2*tmp.gsss[7] + ze22*fsps[m1][7*3+1];
		gsds[m][7*6+2]=QC[2]*gsps[m][7*3+2]+WQ[2]*gsps[m1][7*3+2]
				+ eta2*tmp.gsss[7] + ze22*fsps[m1][4*3+2];
		gsds[m][7*6+3]=QC[0]*gsps[m][7*3+1]+WQ[0]*gsps[m1][7*3+1];
		gsds[m][7*6+4]=QC[1]*gsps[m][7*3+2]+WQ[1]*gsps[m1][7*3+2]
				+ ze22*fsps[m1][7*3+2];
		gsds[m][7*6+5]=QC[2]*gsps[m][7*3+0]+WQ[2]*gsps[m1][7*3+0]
				+ ze22*fsps[m1][4*3+0];
		gsds[m][8*6+0]=QC[0]*gsps[m][8*3+0]+WQ[0]*gsps[m1][8*3+0]
				+ eta2*tmp.gsss[8] + ze22*fsps[m1][5*3+0];
		gsds[m][8*6+1]=QC[1]*gsps[m][8*3+1]+WQ[1]*gsps[m1][8*3+1]
				+ eta2*tmp.gsss[8];
		gsds[m][8*6+2]=QC[2]*gsps[m][8*3+2]+WQ[2]*gsps[m1][8*3+2]
				+ eta2*tmp.gsss[8] + ze22*fsps[m1][8*3+2];
		gsds[m][8*6+3]=QC[0]*gsps[m][8*3+1]+WQ[0]*gsps[m1][8*3+1]
				+ ze22*fsps[m1][5*3+1];
		gsds[m][8*6+4]=QC[1]*gsps[m][8*3+2]+WQ[1]*gsps[m1][8*3+2];
		gsds[m][8*6+5]=QC[2]*gsps[m][8*3+0]+WQ[2]*gsps[m1][8*3+0]
				+ ze22*fsps[m1][8*3+0];
		gsds[m][9*6+0]=QC[0]*gsps[m][9*3+0]+WQ[0]*gsps[m1][9*3+0]
				+ eta2*tmp.gsss[9] + ze2*fsps[m1][1*3+0];
		gsds[m][9*6+1]=QC[1]*gsps[m][9*3+1]+WQ[1]*gsps[m1][9*3+1]
				+ eta2*tmp.gsss[9] + ze23*fsps[m1][6*3+1];
		gsds[m][9*6+2]=QC[2]*gsps[m][9*3+2]+WQ[2]*gsps[m1][9*3+2]
				+ eta2*tmp.gsss[9];
		gsds[m][9*6+3]=QC[0]*gsps[m][9*3+1]+WQ[0]*gsps[m1][9*3+1]
				+ ze2*fsps[m1][1*3+1];
		gsds[m][9*6+4]=QC[1]*gsps[m][9*3+2]+WQ[1]*gsps[m1][9*3+2]
				+ ze23*fsps[m1][6*3+2];
		gsds[m][9*6+5]=QC[2]*gsps[m][9*3+0]+WQ[2]*gsps[m1][9*3+0];
		gsds[m][10*6+0]=QC[0]*gsps[m][10*3+0]+WQ[0]*gsps[m1][10*3+0]
				+ eta2*tmp.gsss[10];
		gsds[m][10*6+1]=QC[1]*gsps[m][10*3+1]+WQ[1]*gsps[m1][10*3+1]
				+ eta2*tmp.gsss[10] + ze2*fsps[m1][2*3+1];
		gsds[m][10*6+2]=QC[2]*gsps[m][10*3+2]+WQ[2]*gsps[m1][10*3+2]
				+ eta2*tmp.gsss[10] + ze23*fsps[m1][7*3+2];
		gsds[m][10*6+3]=QC[0]*gsps[m][10*3+1]+WQ[0]*gsps[m1][10*3+1];
		gsds[m][10*6+4]=QC[1]*gsps[m][10*3+2]+WQ[1]*gsps[m1][10*3+2]
				+ ze2*fsps[m1][2*3+2];
		gsds[m][10*6+5]=QC[2]*gsps[m][10*3+0]+WQ[2]*gsps[m1][10*3+0]
				+ ze23*fsps[m1][7*3+0];
		gsds[m][11*6+0]=QC[0]*gsps[m][11*3+0]+WQ[0]*gsps[m1][11*3+0]
				+ eta2*tmp.gsss[11] + ze23*fsps[m1][8*3+0];
		gsds[m][11*6+1]=QC[1]*gsps[m][11*3+1]+WQ[1]*gsps[m1][11*3+1]
				+ eta2*tmp.gsss[11];
		gsds[m][11*6+2]=QC[2]*gsps[m][11*3+2]+WQ[2]*gsps[m1][11*3+2]
				+ eta2*tmp.gsss[11] + ze2*fsps[m1][0*3+2];
		gsds[m][11*6+3]=QC[0]*gsps[m][11*3+1]+WQ[0]*gsps[m1][11*3+1]
				+ ze23*fsps[m1][8*3+1];
		gsds[m][11*6+4]=QC[1]*gsps[m][11*3+2]+WQ[1]*gsps[m1][11*3+2];
		gsds[m][11*6+5]=QC[2]*gsps[m][11*3+0]+WQ[2]*gsps[m1][11*3+0]
				+ ze2*fsps[m1][0*3+0];
		gsds[m][12*6+0]=QC[0]*gsps[m][12*3+0]+WQ[0]*gsps[m1][12*3+0]
				+ eta2*tmp.gsss[12] + ze22*fsps[m1][9*3+0];
		gsds[m][12*6+1]=QC[1]*gsps[m][12*3+1]+WQ[1]*gsps[m1][12*3+1]
				+ eta2*tmp.gsss[12] + ze2*fsps[m1][8*3+1];
		gsds[m][12*6+2]=QC[2]*gsps[m][12*3+2]+WQ[2]*gsps[m1][12*3+2]
				+ eta2*tmp.gsss[12] + ze2*fsps[m1][3*3+2];
		gsds[m][12*6+3]=QC[0]*gsps[m][12*3+1]+WQ[0]*gsps[m1][12*3+1]
				+ ze22*fsps[m1][9*3+1];
		gsds[m][12*6+4]=QC[1]*gsps[m][12*3+2]+WQ[1]*gsps[m1][12*3+2]
				+ ze2*fsps[m1][8*3+2];
		gsds[m][12*6+5]=QC[2]*gsps[m][12*3+0]+WQ[2]*gsps[m1][12*3+0]
				+ ze2*fsps[m1][3*3+0];
		gsds[m][13*6+0]=QC[0]*gsps[m][13*3+0]+WQ[0]*gsps[m1][13*3+0]
				+ eta2*tmp.gsss[13] + ze2*fsps[m1][4*3+0];
		gsds[m][13*6+1]=QC[1]*gsps[m][13*3+1]+WQ[1]*gsps[m1][13*3+1]
				+ eta2*tmp.gsss[13] + ze22*fsps[m1][9*3+1];
		gsds[m][13*6+2]=QC[2]*gsps[m][13*3+2]+WQ[2]*gsps[m1][13*3+2]
				+ eta2*tmp.gsss[13] + ze2*fsps[m1][6*3+2];
		gsds[m][13*6+3]=QC[0]*gsps[m][13*3+1]+WQ[0]*gsps[m1][13*3+1]
				+ ze2*fsps[m1][4*3+1];
		gsds[m][13*6+4]=QC[1]*gsps[m][13*3+2]+WQ[1]*gsps[m1][13*3+2]
				+ ze22*fsps[m1][9*3+2];
		gsds[m][13*6+5]=QC[2]*gsps[m][13*3+0]+WQ[2]*gsps[m1][13*3+0]
				+ ze2*fsps[m1][6*3+0];
		gsds[m][14*6+0]=QC[0]*gsps[m][14*3+0]+WQ[0]*gsps[m1][14*3+0]
				+ eta2*tmp.gsss[14] + ze2*fsps[m1][7*3+0];
		gsds[m][14*6+1]=QC[1]*gsps[m][14*3+1]+WQ[1]*gsps[m1][14*3+1]
				+ eta2*tmp.gsss[14] + ze2*fsps[m1][5*3+1];
		gsds[m][14*6+2]=QC[2]*gsps[m][14*3+2]+WQ[2]*gsps[m1][14*3+2]
				+ eta2*tmp.gsss[14] + ze22*fsps[m1][9*3+2];
		gsds[m][14*6+3]=QC[0]*gsps[m][14*3+1]+WQ[0]*gsps[m1][14*3+1]
				+ ze2*fsps[m1][7*3+1];
		gsds[m][14*6+4]=QC[1]*gsps[m][14*3+2]+WQ[1]*gsps[m1][14*3+2]
				+ ze2*fsps[m1][5*3+2];
		gsds[m][14*6+5]=QC[2]*gsps[m][14*3+0]+WQ[2]*gsps[m1][14*3+0]
				+ ze22*fsps[m1][9*3+0];
	    }	// end gsds loop
	    for (i=0; i<15*6; i++) GSDS[i] += gsds[0][i];
	    // psfs (m=1)
	    for (i=0; i<3*3; i++) tmp.psps[i]=psps[1-1][i]-re*psps[2-1][i];
	    psfs[0*10+0] = QC[0]*psds[1-1][0*6+0] + WQ[0]*psds[2-1][0*6+0]
			 + eta*tmp.psps[0*3+0] + ze2*ssds[0];
	    psfs[0*10+1] = QC[1]*psds[1-1][0*6+1] + WQ[1]*psds[2-1][0*6+1]
			 + eta*tmp.psps[0*3+1];
	    psfs[0*10+2] = QC[2]*psds[1-1][0*6+2] + WQ[2]*psds[2-1][0*6+2]
			 + eta*tmp.psps[0*3+2];
	    psfs[0*10+3] = QC[0]*psds[1-1][0*6+3] + WQ[0]*psds[2-1][0*6+3]
			 + eta2*tmp.psps[0*3+1] + ze2*ssds[3];
	    psfs[0*10+4] = QC[1]*psds[1-1][0*6+4] + WQ[1]*psds[2-1][0*6+4]
			 + eta2*tmp.psps[0*3+2];
	    psfs[0*10+5] = QC[2]*psds[1-1][0*6+5] + WQ[2]*psds[2-1][0*6+5]
			 + eta2*tmp.psps[0*3+0];
	    psfs[0*10+6] = QC[0]*psds[1-1][0*6+1] + WQ[0]*psds[2-1][0*6+1]
			 + ze2*ssds[1];
	    psfs[0*10+7] = QC[1]*psds[1-1][0*6+2] + WQ[1]*psds[2-1][0*6+2];
	    psfs[0*10+8] = QC[2]*psds[1-1][0*6+0] + WQ[2]*psds[2-1][0*6+0];
	    psfs[0*10+9] = QC[0]*psds[1-1][0*6+4] + WQ[0]*psds[2-1][0*6+4]
			 + ze2*ssds[4];
	    psfs[1*10+0] = QC[0]*psds[1-1][1*6+0] + WQ[0]*psds[2-1][1*6+0]
			 + eta*tmp.psps[1*3+0];
	    psfs[1*10+1] = QC[1]*psds[1-1][1*6+1] + WQ[1]*psds[2-1][1*6+1]
			 + eta*tmp.psps[1*3+1] + ze2*ssds[1];
	    psfs[1*10+2] = QC[2]*psds[1-1][1*6+2] + WQ[2]*psds[2-1][1*6+2]
			 + eta*tmp.psps[1*3+2];
	    psfs[1*10+3] = QC[0]*psds[1-1][1*6+3] + WQ[0]*psds[2-1][1*6+3]
			 + eta2*tmp.psps[1*3+1];
	    psfs[1*10+4] = QC[1]*psds[1-1][1*6+4] + WQ[1]*psds[2-1][1*6+4]
			 + eta2*tmp.psps[1*3+2] + ze2*ssds[4];
	    psfs[1*10+5] = QC[2]*psds[1-1][1*6+5] + WQ[2]*psds[2-1][1*6+5]
			 + eta2*tmp.psps[1*3+0];
	    psfs[1*10+6] = QC[0]*psds[1-1][1*6+1] + WQ[0]*psds[2-1][1*6+1];
	    psfs[1*10+7] = QC[1]*psds[1-1][1*6+2] + WQ[1]*psds[2-1][1*6+2]
			 + ze2*ssds[2];
	    psfs[1*10+8] = QC[2]*psds[1-1][1*6+0] + WQ[2]*psds[2-1][1*6+0];
	    psfs[1*10+9] = QC[0]*psds[1-1][1*6+4] + WQ[0]*psds[2-1][1*6+4];
	    psfs[2*10+0] = QC[0]*psds[1-1][2*6+0] + WQ[0]*psds[2-1][2*6+0]
			 + eta*tmp.psps[2*3+0];
	    psfs[2*10+1] = QC[1]*psds[1-1][2*6+1] + WQ[1]*psds[2-1][2*6+1]
			 + eta*tmp.psps[2*3+1];
	    psfs[2*10+2] = QC[2]*psds[1-1][2*6+2] + WQ[2]*psds[2-1][2*6+2]
			 + eta*tmp.psps[2*3+2] + ze2*ssds[2];
	    psfs[2*10+3] = QC[0]*psds[1-1][2*6+3] + WQ[0]*psds[2-1][2*6+3]
			 + eta2*tmp.psps[2*3+1];
	    psfs[2*10+4] = QC[1]*psds[1-1][2*6+4] + WQ[1]*psds[2-1][2*6+4]
			 + eta2*tmp.psps[2*3+2];
	    psfs[2*10+5] = QC[2]*psds[1-1][2*6+5] + WQ[2]*psds[2-1][2*6+5]
			 + eta2*tmp.psps[2*3+0] + ze2*ssds[5];
	    psfs[2*10+6] = QC[0]*psds[1-1][2*6+1] + WQ[0]*psds[2-1][2*6+1];
	    psfs[2*10+7] = QC[1]*psds[1-1][2*6+2] + WQ[1]*psds[2-1][2*6+2];
	    psfs[2*10+8] = QC[2]*psds[1-1][2*6+0] + WQ[2]*psds[2-1][2*6+0]
			 + ze2*ssds[0];
	    psfs[2*10+9] = QC[0]*psds[1-1][2*6+4] + WQ[0]*psds[2-1][2*6+4];
	    // dsfs (m=0,1)
	    for (m=0, m1=1; m<=1; m++, m1++) {
		for (i=0; i<6*3; i++)
		    tmp.dsps[i]=dsps[m][i]-re*dsps[m1][i];

		dsfs[m][0*10+0]=QC[0]*dsds[m][0*6+0]+WQ[0]*dsds[m1][0*6+0]
		    +eta*tmp.dsps[0*3+0]+ze22*psds[m1-1][0*6+0];
		dsfs[m][0*10+1]=QC[1]*dsds[m][0*6+1]+WQ[1]*dsds[m1][0*6+1]
		    +eta*tmp.dsps[0*3+1];
		dsfs[m][0*10+2]=QC[2]*dsds[m][0*6+2]+WQ[2]*dsds[m1][0*6+2]
		    +eta*tmp.dsps[0*3+2];
		dsfs[m][0*10+3]=QC[0]*dsds[m][0*6+3]+WQ[0]*dsds[m1][0*6+3]
		    +eta2*tmp.dsps[0*3+1]+ze22*psds[m1-1][0*6+3];
		dsfs[m][0*10+4]=QC[1]*dsds[m][0*6+4]+WQ[1]*dsds[m1][0*6+4]
		    +eta2*tmp.dsps[0*3+2];
		dsfs[m][0*10+5]=QC[2]*dsds[m][0*6+5]+WQ[2]*dsds[m1][0*6+5]
		    +eta2*tmp.dsps[0*3+0];
		dsfs[m][0*10+6]=QC[0]*dsds[m][0*6+1]+WQ[0]*dsds[m1][0*6+1]
		    +ze22*psds[m1-1][0*6+1];
		dsfs[m][0*10+7]=QC[1]*dsds[m][0*6+2]+WQ[1]*dsds[m1][0*6+2];
		dsfs[m][0*10+8]=QC[2]*dsds[m][0*6+0]+WQ[2]*dsds[m1][0*6+0];
		dsfs[m][0*10+9]=QC[0]*dsds[m][0*6+4]+WQ[0]*dsds[m1][0*6+4]
		    + ze22*psds[m1-1][0*6+4];
		dsfs[m][1*10+0]=QC[0]*dsds[m][1*6+0]+WQ[0]*dsds[m1][1*6+0]
		    + eta*tmp.dsps[1*3+0];
		dsfs[m][1*10+1]=QC[1]*dsds[m][1*6+1]+WQ[1]*dsds[m1][1*6+1]
		    +eta*tmp.dsps[1*3+1]+ze22*psds[m1-1][1*6+1];
		dsfs[m][1*10+2]=QC[2]*dsds[m][1*6+2]+WQ[2]*dsds[m1][1*6+2]
		    + eta*tmp.dsps[1*3+2];
		dsfs[m][1*10+3]=QC[0]*dsds[m][1*6+3]+WQ[0]*dsds[m1][1*6+3]
		    + eta2*tmp.dsps[1*3+1];
		dsfs[m][1*10+4]=QC[1]*dsds[m][1*6+4]+WQ[1]*dsds[m1][1*6+4]
		    + eta2*tmp.dsps[1*3+2] + ze22*psds[m1-1][1*6+4];
		dsfs[m][1*10+5]=QC[2]*dsds[m][1*6+5]+WQ[2]*dsds[m1][1*6+5]
		    + eta2*tmp.dsps[1*3+0];
		dsfs[m][1*10+6]=QC[0]*dsds[m][1*6+1]+WQ[0]*dsds[m1][1*6+1];
		dsfs[m][1*10+7]=QC[1]*dsds[m][1*6+2]+WQ[1]*dsds[m1][1*6+2]
		    + ze22*psds[m1-1][1*6+2];
		dsfs[m][1*10+8]=QC[2]*dsds[m][1*6+0]+WQ[2]*dsds[m1][1*6+0];
		dsfs[m][1*10+9]=QC[0]*dsds[m][1*6+4]+WQ[0]*dsds[m1][1*6+4];
		dsfs[m][2*10+0]=QC[0]*dsds[m][2*6+0]+WQ[0]*dsds[m1][2*6+0]
		    + eta*tmp.dsps[2*3+0];
		dsfs[m][2*10+1]=QC[1]*dsds[m][2*6+1]+WQ[1]*dsds[m1][2*6+1]
		    + eta*tmp.dsps[2*3+1];
		dsfs[m][2*10+2]=QC[2]*dsds[m][2*6+2]+WQ[2]*dsds[m1][2*6+2]
		    + eta*tmp.dsps[2*3+2] + ze22*psds[m1-1][2*6+2];
		dsfs[m][2*10+3]=QC[0]*dsds[m][2*6+3]+WQ[0]*dsds[m1][2*6+3]
		    + eta2*tmp.dsps[2*3+1];
		dsfs[m][2*10+4]=QC[1]*dsds[m][2*6+4]+WQ[1]*dsds[m1][2*6+4]
		    + eta2*tmp.dsps[2*3+2];
		dsfs[m][2*10+5]=QC[2]*dsds[m][2*6+5]+WQ[2]*dsds[m1][2*6+5]
		    + eta2*tmp.dsps[2*3+0] + ze22*psds[m1-1][2*6+5];
		dsfs[m][2*10+6]=QC[0]*dsds[m][2*6+1]+WQ[0]*dsds[m1][2*6+1];
		dsfs[m][2*10+7]=QC[1]*dsds[m][2*6+2]+WQ[1]*dsds[m1][2*6+2];
		dsfs[m][2*10+8]=QC[2]*dsds[m][2*6+0]+WQ[2]*dsds[m1][2*6+0]
		    + ze22*psds[m1-1][2*6+0];
		dsfs[m][2*10+9]=QC[0]*dsds[m][2*6+4]+WQ[0]*dsds[m1][2*6+4];
		dsfs[m][3*10+0]=QC[0]*dsds[m][3*6+0]+WQ[0]*dsds[m1][3*6+0]
		    + eta*tmp.dsps[3*3+0] + ze2*psds[m1-1][1*6+0];
		dsfs[m][3*10+1]=QC[1]*dsds[m][3*6+1]+WQ[1]*dsds[m1][3*6+1]
		    + eta*tmp.dsps[3*3+1] + ze2*psds[m1-1][0*6+1];
		dsfs[m][3*10+2]=QC[2]*dsds[m][3*6+2]+WQ[2]*dsds[m1][3*6+2]
		    + eta*tmp.dsps[3*3+2];
		dsfs[m][3*10+3]=QC[0]*dsds[m][3*6+3]+WQ[0]*dsds[m1][3*6+3]
		    + eta2*tmp.dsps[3*3+1] + ze2*psds[m1-1][1*6+3];
		dsfs[m][3*10+4]=QC[1]*dsds[m][3*6+4]+WQ[1]*dsds[m1][3*6+4]
		    + eta2*tmp.dsps[3*3+2] + ze2*psds[m1-1][0*6+4];
		dsfs[m][3*10+5]=QC[2]*dsds[m][3*6+5]+WQ[2]*dsds[m1][3*6+5]
		    + eta2*tmp.dsps[3*3+0];
		dsfs[m][3*10+6]=QC[0]*dsds[m][3*6+1]+WQ[0]*dsds[m1][3*6+1]
		    + ze2*psds[m1-1][1*6+1];
		dsfs[m][3*10+7]=QC[1]*dsds[m][3*6+2]+WQ[1]*dsds[m1][3*6+2]
		    + ze2*psds[m1-1][0*6+2];
		dsfs[m][3*10+8]=QC[2]*dsds[m][3*6+0]+WQ[2]*dsds[m1][3*6+0];
		dsfs[m][3*10+9]=QC[0]*dsds[m][3*6+4]+WQ[0]*dsds[m1][3*6+4]
		    + ze2*psds[m1-1][1*6+4];
		dsfs[m][4*10+0]=QC[0]*dsds[m][4*6+0]+WQ[0]*dsds[m1][4*6+0]
		    + eta*tmp.dsps[4*3+0];
		dsfs[m][4*10+1]=QC[1]*dsds[m][4*6+1]+WQ[1]*dsds[m1][4*6+1]
		    + eta*tmp.dsps[4*3+1] + ze2*psds[m1-1][2*6+1];
		dsfs[m][4*10+2]=QC[2]*dsds[m][4*6+2]+WQ[2]*dsds[m1][4*6+2]
		    + eta*tmp.dsps[4*3+2] + ze2*psds[m1-1][1*6+2];
		dsfs[m][4*10+3]=QC[0]*dsds[m][4*6+3]+WQ[0]*dsds[m1][4*6+3]
		    + eta2*tmp.dsps[4*3+1];
		dsfs[m][4*10+4]=QC[1]*dsds[m][4*6+4]+WQ[1]*dsds[m1][4*6+4]
		    + eta2*tmp.dsps[4*3+2] + ze2*psds[m1-1][2*6+4];
		dsfs[m][4*10+5]=QC[2]*dsds[m][4*6+5]+WQ[2]*dsds[m1][4*6+5]
		    + eta2*tmp.dsps[4*3+0] + ze2*psds[m1-1][1*6+5];
		dsfs[m][4*10+6]=QC[0]*dsds[m][4*6+1]+WQ[0]*dsds[m1][4*6+1];
		dsfs[m][4*10+7]=QC[1]*dsds[m][4*6+2]+WQ[1]*dsds[m1][4*6+2]
		    + ze2*psds[m1-1][2*6+2];
		dsfs[m][4*10+8]=QC[2]*dsds[m][4*6+0]+WQ[2]*dsds[m1][4*6+0]
		    + ze2*psds[m1-1][1*6+0];
		dsfs[m][4*10+9]=QC[0]*dsds[m][4*6+4]+WQ[0]*dsds[m1][4*6+4];
		dsfs[m][5*10+0]=QC[0]*dsds[m][5*6+0]+WQ[0]*dsds[m1][5*6+0]
		    + eta*tmp.dsps[5*3+0] + ze2*psds[m1-1][2*6+0];
		dsfs[m][5*10+1]=QC[1]*dsds[m][5*6+1]+WQ[1]*dsds[m1][5*6+1]
		    + eta*tmp.dsps[5*3+1];
		dsfs[m][5*10+2]=QC[2]*dsds[m][5*6+2]+WQ[2]*dsds[m1][5*6+2]
		    + eta*tmp.dsps[5*3+2] + ze2*psds[m1-1][0*6+2];
		dsfs[m][5*10+3]=QC[0]*dsds[m][5*6+3]+WQ[0]*dsds[m1][5*6+3]
		    + eta2*tmp.dsps[5*3+1] + ze2*psds[m1-1][2*6+3];
		dsfs[m][5*10+4]=QC[1]*dsds[m][5*6+4]+WQ[1]*dsds[m1][5*6+4]
		    + eta2*tmp.dsps[5*3+2];
		dsfs[m][5*10+5]=QC[2]*dsds[m][5*6+5]+WQ[2]*dsds[m1][5*6+5]
		    + eta2*tmp.dsps[5*3+0] + ze2*psds[m1-1][0*6+5];
		dsfs[m][5*10+6]=QC[0]*dsds[m][5*6+1]+WQ[0]*dsds[m1][5*6+1]
		    + ze2*psds[m1-1][2*6+1];
		dsfs[m][5*10+7]=QC[1]*dsds[m][5*6+2]+WQ[1]*dsds[m1][5*6+2];
		dsfs[m][5*10+8]=QC[2]*dsds[m][5*6+0]+WQ[2]*dsds[m1][5*6+0]
		    + ze2*psds[m1-1][0*6+0];
		dsfs[m][5*10+9]=QC[0]*dsds[m][5*6+4]+WQ[0]*dsds[m1][5*6+4]
		    + ze2*psds[m1-1][2*6+4];
	    }	// end dsfs loop
	    for (i=0; i<6*10; i++) DSFS[i] += dsfs[0][i];
	    // fsfs (m=0,1)
	    for (m=0, m1=1; m<=1; m++, m1++) {
		for (i=0; i<10*3; i++)
		    tmp.fsps[i]=fsps[m][i]-re*fsps[m1][i];

		fsfs[m][0*10+0]=QC[0]*fsds[m][0*6+0]+WQ[0]*fsds[m1][0*6+0]
				+eta*tmp.fsps[0*3+0]+ze23*dsds[m1][0*6+0];
		fsfs[m][0*10+1]=QC[1]*fsds[m][0*6+1]+WQ[1]*fsds[m1][0*6+1]
				+eta*tmp.fsps[0*3+1];
		fsfs[m][0*10+2]=QC[2]*fsds[m][0*6+2]+WQ[2]*fsds[m1][0*6+2]
				+eta*tmp.fsps[0*3+2];
		fsfs[m][0*10+3]=QC[0]*fsds[m][0*6+3]+WQ[0]*fsds[m1][0*6+3]
				+eta2*tmp.fsps[0*3+1]+ze23*dsds[m1][0*6+3];
		fsfs[m][0*10+4]=QC[1]*fsds[m][0*6+4]+WQ[1]*fsds[m1][0*6+4]
				+eta2*tmp.fsps[0*3+2];
		fsfs[m][0*10+5]=QC[2]*fsds[m][0*6+5]+WQ[2]*fsds[m1][0*6+5]
				+eta2*tmp.fsps[0*3+0];
		fsfs[m][0*10+6]=QC[0]*fsds[m][0*6+1]+WQ[0]*fsds[m1][0*6+1]
				+ze23*dsds[m1][0*6+1];
		fsfs[m][0*10+7]=QC[1]*fsds[m][0*6+2]+WQ[1]*fsds[m1][0*6+2];
		fsfs[m][0*10+8]=QC[2]*fsds[m][0*6+0]+WQ[2]*fsds[m1][0*6+0];
		fsfs[m][0*10+9]=QC[0]*fsds[m][0*6+4]+WQ[0]*fsds[m1][0*6+4]
				+ze23*dsds[m1][0*6+4];
		fsfs[m][1*10+0]=QC[0]*fsds[m][1*6+0]+WQ[0]*fsds[m1][1*6+0]
				+eta*tmp.fsps[1*3+0];
		fsfs[m][1*10+1]=QC[1]*fsds[m][1*6+1]+WQ[1]*fsds[m1][1*6+1]
				+eta*tmp.fsps[1*3+1]+ze23*dsds[m1][1*6+1];
		fsfs[m][1*10+2]=QC[2]*fsds[m][1*6+2]+WQ[2]*fsds[m1][1*6+2]
				+eta*tmp.fsps[1*3+2];
		fsfs[m][1*10+3]=QC[0]*fsds[m][1*6+3]+WQ[0]*fsds[m1][1*6+3]
				+eta2*tmp.fsps[1*3+1];
		fsfs[m][1*10+4]=QC[1]*fsds[m][1*6+4]+WQ[1]*fsds[m1][1*6+4]
				+eta2*tmp.fsps[1*3+2]+ze23*dsds[m1][1*6+4];
		fsfs[m][1*10+5]=QC[2]*fsds[m][1*6+5]+WQ[2]*fsds[m1][1*6+5]
				+eta2*tmp.fsps[1*3+0];
		fsfs[m][1*10+6]=QC[0]*fsds[m][1*6+1]+WQ[0]*fsds[m1][1*6+1];
		fsfs[m][1*10+7]=QC[1]*fsds[m][1*6+2]+WQ[1]*fsds[m1][1*6+2]
				+ze23*dsds[m1][1*6+2];
		fsfs[m][1*10+8]=QC[2]*fsds[m][1*6+0]+WQ[2]*fsds[m1][1*6+0];
		fsfs[m][1*10+9]=QC[0]*fsds[m][1*6+4]+WQ[0]*fsds[m1][1*6+4];
		fsfs[m][2*10+0]=QC[0]*fsds[m][2*6+0]+WQ[0]*fsds[m1][2*6+0]
				+eta*tmp.fsps[2*3+0];
		fsfs[m][2*10+1]=QC[1]*fsds[m][2*6+1]+WQ[1]*fsds[m1][2*6+1]
				+eta*tmp.fsps[2*3+1];
		fsfs[m][2*10+2]=QC[2]*fsds[m][2*6+2]+WQ[2]*fsds[m1][2*6+2]
				+eta*tmp.fsps[2*3+2]+ze23*dsds[m1][2*6+2];
		fsfs[m][2*10+3]=QC[0]*fsds[m][2*6+3]+WQ[0]*fsds[m1][2*6+3]
				+eta2*tmp.fsps[2*3+1];
		fsfs[m][2*10+4]=QC[1]*fsds[m][2*6+4]+WQ[1]*fsds[m1][2*6+4]
				+eta2*tmp.fsps[2*3+2];
		fsfs[m][2*10+5]=QC[2]*fsds[m][2*6+5]+WQ[2]*fsds[m1][2*6+5]
				+eta2*tmp.fsps[2*3+0]+ze23*dsds[m1][2*6+5];
		fsfs[m][2*10+6]=QC[0]*fsds[m][2*6+1]+WQ[0]*fsds[m1][2*6+1];
		fsfs[m][2*10+7]=QC[1]*fsds[m][2*6+2]+WQ[1]*fsds[m1][2*6+2];
		fsfs[m][2*10+8]=QC[2]*fsds[m][2*6+0]+WQ[2]*fsds[m1][2*6+0]
				+ze23*dsds[m1][2*6+0];
		fsfs[m][2*10+9]=QC[0]*fsds[m][2*6+4]+WQ[0]*fsds[m1][2*6+4];
		fsfs[m][3*10+0]=QC[0]*fsds[m][3*6+0]+WQ[0]*fsds[m1][3*6+0]
				+eta*tmp.fsps[3*3+0]+ze22*dsds[m1][3*6+0];
		fsfs[m][3*10+1]=QC[1]*fsds[m][3*6+1]+WQ[1]*fsds[m1][3*6+1]
				+eta*tmp.fsps[3*3+1]+ze2*dsds[m1][0*6+1];
		fsfs[m][3*10+2]=QC[2]*fsds[m][3*6+2]+WQ[2]*fsds[m1][3*6+2]
				+eta*tmp.fsps[3*3+2];
		fsfs[m][3*10+3]=QC[0]*fsds[m][3*6+3]+WQ[0]*fsds[m1][3*6+3]
				+eta2*tmp.fsps[3*3+1]+ze22*dsds[m1][3*6+3];
		fsfs[m][3*10+4]=QC[1]*fsds[m][3*6+4]+WQ[1]*fsds[m1][3*6+4]
				+eta2*tmp.fsps[3*3+2]+ze2*dsds[m1][0*6+4];
		fsfs[m][3*10+5]=QC[2]*fsds[m][3*6+5]+WQ[2]*fsds[m1][3*6+5]
				+eta2*tmp.fsps[3*3+0];
		fsfs[m][3*10+6]=QC[0]*fsds[m][3*6+1]+WQ[0]*fsds[m1][3*6+1]
				+ze22*dsds[m1][3*6+1];
		fsfs[m][3*10+7]=QC[1]*fsds[m][3*6+2]+WQ[1]*fsds[m1][3*6+2]
				+ze2*dsds[m1][0*6+2];
		fsfs[m][3*10+8]=QC[2]*fsds[m][3*6+0]+WQ[2]*fsds[m1][3*6+0];
		fsfs[m][3*10+9]=QC[0]*fsds[m][3*6+4]+WQ[0]*fsds[m1][3*6+4]
				+ze22*dsds[m1][3*6+4];
		fsfs[m][4*10+0]=QC[0]*fsds[m][4*6+0]+WQ[0]*fsds[m1][4*6+0]
				+eta*tmp.fsps[4*3+0];
		fsfs[m][4*10+1]=QC[1]*fsds[m][4*6+1]+WQ[1]*fsds[m1][4*6+1]
				+eta*tmp.fsps[4*3+1]+ze22*dsds[m1][4*6+1];
		fsfs[m][4*10+2]=QC[2]*fsds[m][4*6+2]+WQ[2]*fsds[m1][4*6+2]
				+eta*tmp.fsps[4*3+2]+ze2*dsds[m1][1*6+2];
		fsfs[m][4*10+3]=QC[0]*fsds[m][4*6+3]+WQ[0]*fsds[m1][4*6+3]
				+eta2*tmp.fsps[4*3+1];
		fsfs[m][4*10+4]=QC[1]*fsds[m][4*6+4]+WQ[1]*fsds[m1][4*6+4]
				+eta2*tmp.fsps[4*3+2]+ze22*dsds[m1][4*6+4];
		fsfs[m][4*10+5]=QC[2]*fsds[m][4*6+5]+WQ[2]*fsds[m1][4*6+5]
				+eta2*tmp.fsps[4*3+0]+ze2*dsds[m1][1*6+5];
		fsfs[m][4*10+6]=QC[0]*fsds[m][4*6+1]+WQ[0]*fsds[m1][4*6+1];
		fsfs[m][4*10+7]=QC[1]*fsds[m][4*6+2]+WQ[1]*fsds[m1][4*6+2]
				+ze22*dsds[m1][4*6+2];
		fsfs[m][4*10+8]=QC[2]*fsds[m][4*6+0]+WQ[2]*fsds[m1][4*6+0]
				+ze2*dsds[m1][1*6+0];
		fsfs[m][4*10+9]=QC[0]*fsds[m][4*6+4]+WQ[0]*fsds[m1][4*6+4];
		fsfs[m][5*10+0]=QC[0]*fsds[m][5*6+0]+WQ[0]*fsds[m1][5*6+0]
				+eta*tmp.fsps[5*3+0]+ze2*dsds[m1][2*6+0];
		fsfs[m][5*10+1]=QC[1]*fsds[m][5*6+1]+WQ[1]*fsds[m1][5*6+1]
				+eta*tmp.fsps[5*3+1];
		fsfs[m][5*10+2]=QC[2]*fsds[m][5*6+2]+WQ[2]*fsds[m1][5*6+2]
				+eta*tmp.fsps[5*3+2]+ze22*dsds[m1][5*6+2];
		fsfs[m][5*10+3]=QC[0]*fsds[m][5*6+3]+WQ[0]*fsds[m1][5*6+3]
				+eta2*tmp.fsps[5*3+1]+ze2*dsds[m1][2*6+3];
		fsfs[m][5*10+4]=QC[1]*fsds[m][5*6+4]+WQ[1]*fsds[m1][5*6+4]
				+eta2*tmp.fsps[5*3+2];
		fsfs[m][5*10+5]=QC[2]*fsds[m][5*6+5]+WQ[2]*fsds[m1][5*6+5]
				+eta2*tmp.fsps[5*3+0]+ze22*dsds[m1][5*6+5];
		fsfs[m][5*10+6]=QC[0]*fsds[m][5*6+1]+WQ[0]*fsds[m1][5*6+1]
				+ze2*dsds[m1][2*6+1];
		fsfs[m][5*10+7]=QC[1]*fsds[m][5*6+2]+WQ[1]*fsds[m1][5*6+2];
		fsfs[m][5*10+8]=QC[2]*fsds[m][5*6+0]+WQ[2]*fsds[m1][5*6+0]
				+ze22*dsds[m1][5*6+0];
		fsfs[m][5*10+9]=QC[0]*fsds[m][5*6+4]+WQ[0]*fsds[m1][5*6+4]
				+ze2*dsds[m1][2*6+4];
		fsfs[m][6*10+0]=QC[0]*fsds[m][6*6+0]+WQ[0]*fsds[m1][6*6+0]
				+eta*tmp.fsps[6*3+0]+ze2*dsds[m1][1*6+0];
		fsfs[m][6*10+1]=QC[1]*fsds[m][6*6+1]+WQ[1]*fsds[m1][6*6+1]
				+eta*tmp.fsps[6*3+1]+ze22*dsds[m1][3*6+1];
		fsfs[m][6*10+2]=QC[2]*fsds[m][6*6+2]+WQ[2]*fsds[m1][6*6+2]
				+eta*tmp.fsps[6*3+2];
		fsfs[m][6*10+3]=QC[0]*fsds[m][6*6+3]+WQ[0]*fsds[m1][6*6+3]
				+eta2*tmp.fsps[6*3+1]+ze2*dsds[m1][1*6+3];
		fsfs[m][6*10+4]=QC[1]*fsds[m][6*6+4]+WQ[1]*fsds[m1][6*6+4]
				+eta2*tmp.fsps[6*3+2]+ze22*dsds[m1][3*6+4];
		fsfs[m][6*10+5]=QC[2]*fsds[m][6*6+5]+WQ[2]*fsds[m1][6*6+5]
				+eta2*tmp.fsps[6*3+0];
		fsfs[m][6*10+6]=QC[0]*fsds[m][6*6+1]+WQ[0]*fsds[m1][6*6+1]
				+ze2*dsds[m1][1*6+1];
		fsfs[m][6*10+7]=QC[1]*fsds[m][6*6+2]+WQ[1]*fsds[m1][6*6+2]
				+ze22*dsds[m1][3*6+2];
		fsfs[m][6*10+8]=QC[2]*fsds[m][6*6+0]+WQ[2]*fsds[m1][6*6+0];
		fsfs[m][6*10+9]=QC[0]*fsds[m][6*6+4]+WQ[0]*fsds[m1][6*6+4]
				+ze2*dsds[m1][1*6+4];
		fsfs[m][7*10+0]=QC[0]*fsds[m][7*6+0]+WQ[0]*fsds[m1][7*6+0]
				+eta*tmp.fsps[7*3+0];
		fsfs[m][7*10+1]=QC[1]*fsds[m][7*6+1]+WQ[1]*fsds[m1][7*6+1]
				+eta*tmp.fsps[7*3+1]+ze2*dsds[m1][2*6+1];
		fsfs[m][7*10+2]=QC[2]*fsds[m][7*6+2]+WQ[2]*fsds[m1][7*6+2]
				+eta*tmp.fsps[7*3+2]+ze22*dsds[m1][4*6+2];
		fsfs[m][7*10+3]=QC[0]*fsds[m][7*6+3]+WQ[0]*fsds[m1][7*6+3]
				+eta2*tmp.fsps[7*3+1];
		fsfs[m][7*10+4]=QC[1]*fsds[m][7*6+4]+WQ[1]*fsds[m1][7*6+4]
				+eta2*tmp.fsps[7*3+2]+ze2*dsds[m1][2*6+4];
		fsfs[m][7*10+5]=QC[2]*fsds[m][7*6+5]+WQ[2]*fsds[m1][7*6+5]
				+eta2*tmp.fsps[7*3+0]+ze22*dsds[m1][4*6+5];
		fsfs[m][7*10+6]=QC[0]*fsds[m][7*6+1]+WQ[0]*fsds[m1][7*6+1];
		fsfs[m][7*10+7]=QC[1]*fsds[m][7*6+2]+WQ[1]*fsds[m1][7*6+2]
				+ze2*dsds[m1][2*6+2];
		fsfs[m][7*10+8]=QC[2]*fsds[m][7*6+0]+WQ[2]*fsds[m1][7*6+0]
				+ze22*dsds[m1][4*6+0];
		fsfs[m][7*10+9]=QC[0]*fsds[m][7*6+4]+WQ[0]*fsds[m1][7*6+4];
		fsfs[m][8*10+0]=QC[0]*fsds[m][8*6+0]+WQ[0]*fsds[m1][8*6+0]
				+eta*tmp.fsps[8*3+0]+ze22*dsds[m1][5*6+0];
		fsfs[m][8*10+1]=QC[1]*fsds[m][8*6+1]+WQ[1]*fsds[m1][8*6+1]
				+eta*tmp.fsps[8*3+1];
		fsfs[m][8*10+2]=QC[2]*fsds[m][8*6+2]+WQ[2]*fsds[m1][8*6+2]
				+eta*tmp.fsps[8*3+2]+ze2*dsds[m1][0*6+2];
		fsfs[m][8*10+3]=QC[0]*fsds[m][8*6+3]+WQ[0]*fsds[m1][8*6+3]
				+eta2*tmp.fsps[8*3+1]+ze22*dsds[m1][5*6+3];
		fsfs[m][8*10+4]=QC[1]*fsds[m][8*6+4]+WQ[1]*fsds[m1][8*6+4]
				+eta2*tmp.fsps[8*3+2];
		fsfs[m][8*10+5]=QC[2]*fsds[m][8*6+5]+WQ[2]*fsds[m1][8*6+5]
				+eta2*tmp.fsps[8*3+0]+ze2*dsds[m1][0*6+5];
		fsfs[m][8*10+6]=QC[0]*fsds[m][8*6+1]+WQ[0]*fsds[m1][8*6+1]
				+ze22*dsds[m1][5*6+1];
		fsfs[m][8*10+7]=QC[1]*fsds[m][8*6+2]+WQ[1]*fsds[m1][8*6+2];
		fsfs[m][8*10+8]=QC[2]*fsds[m][8*6+0]+WQ[2]*fsds[m1][8*6+0]
				+ze2*dsds[m1][0*6+0];
		fsfs[m][8*10+9]=QC[0]*fsds[m][8*6+4]+WQ[0]*fsds[m1][8*6+4]
				+ze22*dsds[m1][5*6+4];
		fsfs[m][9*10+0]=QC[0]*fsds[m][9*6+0]+WQ[0]*fsds[m1][9*6+0]
				+eta*tmp.fsps[9*3+0]+ze2*dsds[m1][4*6+0];
		fsfs[m][9*10+1]=QC[1]*fsds[m][9*6+1]+WQ[1]*fsds[m1][9*6+1]
				+eta*tmp.fsps[9*3+1]+ze2*dsds[m1][5*6+1];
		fsfs[m][9*10+2]=QC[2]*fsds[m][9*6+2]+WQ[2]*fsds[m1][9*6+2]
				+eta*tmp.fsps[9*3+2]+ze2*dsds[m1][3*6+2];
		fsfs[m][9*10+3]=QC[0]*fsds[m][9*6+3]+WQ[0]*fsds[m1][9*6+3]
				+eta2*tmp.fsps[9*3+1]+ze2*dsds[m1][4*6+3];
		fsfs[m][9*10+4]=QC[1]*fsds[m][9*6+4]+WQ[1]*fsds[m1][9*6+4]
				+eta2*tmp.fsps[9*3+2]+ze2*dsds[m1][5*6+4];
		fsfs[m][9*10+5]=QC[2]*fsds[m][9*6+5]+WQ[2]*fsds[m1][9*6+5]
				+eta2*tmp.fsps[9*3+0]+ze2*dsds[m1][3*6+5];
		fsfs[m][9*10+6]=QC[0]*fsds[m][9*6+1]+WQ[0]*fsds[m1][9*6+1]
				+ze2*dsds[m1][4*6+1];
		fsfs[m][9*10+7]=QC[1]*fsds[m][9*6+2]+WQ[1]*fsds[m1][9*6+2]
				+ze2*dsds[m1][5*6+2];
		fsfs[m][9*10+8]=QC[2]*fsds[m][9*6+0]+WQ[2]*fsds[m1][9*6+0]
				+ze2*dsds[m1][3*6+0];
		fsfs[m][9*10+9]=QC[0]*fsds[m][9*6+4]+WQ[0]*fsds[m1][9*6+4]
				+ze2*dsds[m1][4*6+4];
	    }	// end fsfs loop
	    for (i=0; i<10*10; i++) FSFS[i] += fsfs[0][i];
	    // gsfs (m=0,1)
	    for (m=0, m1=1; m<=1; m++, m1++) {
		for (i=0; i<15*3; i++)
		  tmp.gsps[i]=gsps[m][i]-re*gsps[m1][i];

		gsfs[m][ 0*10+0] = QC[0]*gsds[m][ 0*6+0] + WQ[0]*gsds[m1][ 0*6+0]
				 + eta*tmp.gsps[ 0*3+0] + ze24*fsds[m1][0*6+0];
		gsfs[m][ 0*10+1] = QC[1]*gsds[m][ 0*6+1] + WQ[1]*gsds[m1][ 0*6+1]
				 + eta*tmp.gsps[ 0*3+1];
		gsfs[m][ 0*10+2] = QC[2]*gsds[m][ 0*6+2] + WQ[2]*gsds[m1][ 0*6+2]
				 + eta*tmp.gsps[ 0*3+2];
		gsfs[m][ 0*10+3] = QC[0]*gsds[m][ 0*6+3] + WQ[0]*gsds[m1][ 0*6+3]
				 + eta2*tmp.gsps[ 0*3+1] + ze24*fsds[m1][0*6+3];
		gsfs[m][ 0*10+4] = QC[1]*gsds[m][ 0*6+4] + WQ[1]*gsds[m1][ 0*6+4]
				 + eta2*tmp.gsps[ 0*3+2];
		gsfs[m][ 0*10+5] = QC[2]*gsds[m][ 0*6+5] + WQ[2]*gsds[m1][ 0*6+5]
				 + eta2*tmp.gsps[ 0*3+0];
		gsfs[m][ 0*10+6] = QC[0]*gsds[m][ 0*6+1] + WQ[0]*gsds[m1][ 0*6+1]
				 + ze24*fsds[m1][0*6+1];
		gsfs[m][ 0*10+7] = QC[1]*gsds[m][ 0*6+2] + WQ[1]*gsds[m1][ 0*6+2];
		gsfs[m][ 0*10+8] = QC[2]*gsds[m][ 0*6+0] + WQ[2]*gsds[m1][ 0*6+0];
		gsfs[m][ 0*10+9] = QC[0]*gsds[m][ 0*6+4] + WQ[0]*gsds[m1][ 0*6+4]
				 + ze24*fsds[m1][0*6+4];
		gsfs[m][ 1*10+0] = QC[0]*gsds[m][ 1*6+0] + WQ[0]*gsds[m1][ 1*6+0]
				 + eta*tmp.gsps[ 1*3+0];
		gsfs[m][ 1*10+1] = QC[1]*gsds[m][ 1*6+1] + WQ[1]*gsds[m1][ 1*6+1]
				 + eta*tmp.gsps[ 1*3+1] + ze24*fsds[m1][1*6+1];
		gsfs[m][ 1*10+2] = QC[2]*gsds[m][ 1*6+2] + WQ[2]*gsds[m1][ 1*6+2]
				 + eta*tmp.gsps[ 1*3+2];
		gsfs[m][ 1*10+3] = QC[0]*gsds[m][ 1*6+3] + WQ[0]*gsds[m1][ 1*6+3]
				 + eta2*tmp.gsps[ 1*3+1];
		gsfs[m][ 1*10+4] = QC[1]*gsds[m][ 1*6+4] + WQ[1]*gsds[m1][ 1*6+4]
				 + eta2*tmp.gsps[ 1*3+2] + ze24*fsds[m1][1*6+4];
		gsfs[m][ 1*10+5] = QC[2]*gsds[m][ 1*6+5] + WQ[2]*gsds[m1][ 1*6+5]
				 + eta2*tmp.gsps[ 1*3+0];
		gsfs[m][ 1*10+6] = QC[0]*gsds[m][ 1*6+1] + WQ[0]*gsds[m1][ 1*6+1];
		gsfs[m][ 1*10+7] = QC[1]*gsds[m][ 1*6+2] + WQ[1]*gsds[m1][ 1*6+2]
				 + ze24*fsds[m1][1*6+2];
		gsfs[m][ 1*10+8] = QC[2]*gsds[m][ 1*6+0] + WQ[2]*gsds[m1][ 1*6+0];
		gsfs[m][ 1*10+9] = QC[0]*gsds[m][ 1*6+4] + WQ[0]*gsds[m1][ 1*6+4];
		gsfs[m][ 2*10+0] = QC[0]*gsds[m][ 2*6+0] + WQ[0]*gsds[m1][ 2*6+0]
				 + eta*tmp.gsps[ 2*3+0];
		gsfs[m][ 2*10+1] = QC[1]*gsds[m][ 2*6+1] + WQ[1]*gsds[m1][ 2*6+1]
				 + eta*tmp.gsps[ 2*3+1];
		gsfs[m][ 2*10+2] = QC[2]*gsds[m][ 2*6+2] + WQ[2]*gsds[m1][ 2*6+2]
				 + eta*tmp.gsps[ 2*3+2] + ze24*fsds[m1][2*6+2];
		gsfs[m][ 2*10+3] = QC[0]*gsds[m][ 2*6+3] + WQ[0]*gsds[m1][ 2*6+3]
				 + eta2*tmp.gsps[ 2*3+1];
		gsfs[m][ 2*10+4] = QC[1]*gsds[m][ 2*6+4] + WQ[1]*gsds[m1][ 2*6+4]
				 + eta2*tmp.gsps[ 2*3+2];
		gsfs[m][ 2*10+5] = QC[2]*gsds[m][ 2*6+5] + WQ[2]*gsds[m1][ 2*6+5]
				 + eta2*tmp.gsps[ 2*3+0] + ze24*fsds[m1][2*6+5];
		gsfs[m][ 2*10+6] = QC[0]*gsds[m][ 2*6+1] + WQ[0]*gsds[m1][ 2*6+1];
		gsfs[m][ 2*10+7] = QC[1]*gsds[m][ 2*6+2] + WQ[1]*gsds[m1][ 2*6+2];
		gsfs[m][ 2*10+8] = QC[2]*gsds[m][ 2*6+0] + WQ[2]*gsds[m1][ 2*6+0]
				 + ze24*fsds[m1][2*6+0];
		gsfs[m][ 2*10+9] = QC[0]*gsds[m][ 2*6+4] + WQ[0]*gsds[m1][ 2*6+4];
		gsfs[m][ 3*10+0] = QC[0]*gsds[m][ 3*6+0] + WQ[0]*gsds[m1][ 3*6+0]
				 + eta*tmp.gsps[ 3*3+0] + ze23*fsds[m1][3*6+0];
		gsfs[m][ 3*10+1] = QC[1]*gsds[m][ 3*6+1] + WQ[1]*gsds[m1][ 3*6+1]
				 + eta*tmp.gsps[ 3*3+1] + ze2*fsds[m1][0*6+1];
		gsfs[m][ 3*10+2] = QC[2]*gsds[m][ 3*6+2] + WQ[2]*gsds[m1][ 3*6+2]
				 + eta*tmp.gsps[ 3*3+2];
		gsfs[m][ 3*10+3] = QC[0]*gsds[m][ 3*6+3] + WQ[0]*gsds[m1][ 3*6+3]
				 + eta2*tmp.gsps[ 3*3+1] + ze23*fsds[m1][3*6+3];
		gsfs[m][ 3*10+4] = QC[1]*gsds[m][ 3*6+4] + WQ[1]*gsds[m1][ 3*6+4]
				 + eta2*tmp.gsps[ 3*3+2] + ze2*fsds[m1][0*6+4];
		gsfs[m][ 3*10+5] = QC[2]*gsds[m][ 3*6+5] + WQ[2]*gsds[m1][ 3*6+5]
				 + eta2*tmp.gsps[ 3*3+0];
		gsfs[m][ 3*10+6] = QC[0]*gsds[m][ 3*6+1] + WQ[0]*gsds[m1][ 3*6+1]
				 + ze23*fsds[m1][3*6+1];
		gsfs[m][ 3*10+7] = QC[1]*gsds[m][ 3*6+2] + WQ[1]*gsds[m1][ 3*6+2]
				 + ze2*fsds[m1][0*6+2];
		gsfs[m][ 3*10+8] = QC[2]*gsds[m][ 3*6+0] + WQ[2]*gsds[m1][ 3*6+0];
		gsfs[m][ 3*10+9] = QC[0]*gsds[m][ 3*6+4] + WQ[0]*gsds[m1][ 3*6+4]
				 + ze23*fsds[m1][3*6+4];
		gsfs[m][ 4*10+0] = QC[0]*gsds[m][ 4*6+0] + WQ[0]*gsds[m1][ 4*6+0]
				 + eta*tmp.gsps[ 4*3+0];
		gsfs[m][ 4*10+1] = QC[1]*gsds[m][ 4*6+1] + WQ[1]*gsds[m1][ 4*6+1]
				 + eta*tmp.gsps[ 4*3+1] + ze23*fsds[m1][4*6+1];
		gsfs[m][ 4*10+2] = QC[2]*gsds[m][ 4*6+2] + WQ[2]*gsds[m1][ 4*6+2]
				 + eta*tmp.gsps[ 4*3+2] + ze2*fsds[m1][1*6+2];
		gsfs[m][ 4*10+3] = QC[0]*gsds[m][ 4*6+3] + WQ[0]*gsds[m1][ 4*6+3]
				 + eta2*tmp.gsps[ 4*3+1];
		gsfs[m][ 4*10+4] = QC[1]*gsds[m][ 4*6+4] + WQ[1]*gsds[m1][ 4*6+4]
				 + eta2*tmp.gsps[ 4*3+2] + ze23*fsds[m1][4*6+4];
		gsfs[m][ 4*10+5] = QC[2]*gsds[m][ 4*6+5] + WQ[2]*gsds[m1][ 4*6+5]
				 + eta2*tmp.gsps[ 4*3+0] + ze2*fsds[m1][1*6+5];
		gsfs[m][ 4*10+6] = QC[0]*gsds[m][ 4*6+1] + WQ[0]*gsds[m1][ 4*6+1];
		gsfs[m][ 4*10+7] = QC[1]*gsds[m][ 4*6+2] + WQ[1]*gsds[m1][ 4*6+2]
				 + ze23*fsds[m1][4*6+2];
		gsfs[m][ 4*10+8] = QC[2]*gsds[m][ 4*6+0] + WQ[2]*gsds[m1][ 4*6+0]
				 + ze2*fsds[m1][1*6+0];
		gsfs[m][ 4*10+9] = QC[0]*gsds[m][ 4*6+4] + WQ[0]*gsds[m1][ 4*6+4];
		gsfs[m][ 5*10+0] = QC[0]*gsds[m][ 5*6+0] + WQ[0]*gsds[m1][ 5*6+0]
				 + eta*tmp.gsps[ 5*3+0] + ze2*fsds[m1][2*6+0];
		gsfs[m][ 5*10+1] = QC[1]*gsds[m][ 5*6+1] + WQ[1]*gsds[m1][ 5*6+1]
				 + eta*tmp.gsps[ 5*3+1];
		gsfs[m][ 5*10+2] = QC[2]*gsds[m][ 5*6+2] + WQ[2]*gsds[m1][ 5*6+2]
				 + eta*tmp.gsps[ 5*3+2] + ze23*fsds[m1][5*6+2];
		gsfs[m][ 5*10+3] = QC[0]*gsds[m][ 5*6+3] + WQ[0]*gsds[m1][ 5*6+3]
				 + eta2*tmp.gsps[ 5*3+1] + ze2*fsds[m1][2*6+3];
		gsfs[m][ 5*10+4] = QC[1]*gsds[m][ 5*6+4] + WQ[1]*gsds[m1][ 5*6+4]
				 + eta2*tmp.gsps[ 5*3+2];
		gsfs[m][ 5*10+5] = QC[2]*gsds[m][ 5*6+5] + WQ[2]*gsds[m1][ 5*6+5]
				 + eta2*tmp.gsps[ 5*3+0] + ze23*fsds[m1][5*6+5];
		gsfs[m][ 5*10+6] = QC[0]*gsds[m][ 5*6+1] + WQ[0]*gsds[m1][ 5*6+1]
				 + ze2*fsds[m1][2*6+1];
		gsfs[m][ 5*10+7] = QC[1]*gsds[m][ 5*6+2] + WQ[1]*gsds[m1][ 5*6+2];
		gsfs[m][ 5*10+8] = QC[2]*gsds[m][ 5*6+0] + WQ[2]*gsds[m1][ 5*6+0]
				 + ze23*fsds[m1][5*6+0];
		gsfs[m][ 5*10+9] = QC[0]*gsds[m][ 5*6+4] + WQ[0]*gsds[m1][ 5*6+4]
				 + ze2*fsds[m1][2*6+4];
		gsfs[m][ 6*10+0] = QC[0]*gsds[m][ 6*6+0] + WQ[0]*gsds[m1][ 6*6+0]
				 + eta*tmp.gsps[ 6*3+0] + ze22*fsds[m1][6*6+0];
		gsfs[m][ 6*10+1] = QC[1]*gsds[m][ 6*6+1] + WQ[1]*gsds[m1][ 6*6+1]
				 + eta*tmp.gsps[ 6*3+1] + ze22*fsds[m1][3*6+1];
		gsfs[m][ 6*10+2] = QC[2]*gsds[m][ 6*6+2] + WQ[2]*gsds[m1][ 6*6+2]
				 + eta*tmp.gsps[ 6*3+2];
		gsfs[m][ 6*10+3] = QC[0]*gsds[m][ 6*6+3] + WQ[0]*gsds[m1][ 6*6+3]
				 + eta2*tmp.gsps[ 6*3+1] + ze22*fsds[m1][6*6+3];
		gsfs[m][ 6*10+4] = QC[1]*gsds[m][ 6*6+4] + WQ[1]*gsds[m1][ 6*6+4]
				 + eta2*tmp.gsps[ 6*3+2] + ze22*fsds[m1][3*6+4];
		gsfs[m][ 6*10+5] = QC[2]*gsds[m][ 6*6+5] + WQ[2]*gsds[m1][ 6*6+5]
				 + eta2*tmp.gsps[ 6*3+0];
		gsfs[m][ 6*10+6] = QC[0]*gsds[m][ 6*6+1] + WQ[0]*gsds[m1][ 6*6+1]
				 + ze22*fsds[m1][6*6+1];
		gsfs[m][ 6*10+7] = QC[1]*gsds[m][ 6*6+2] + WQ[1]*gsds[m1][ 6*6+2]
				 + ze22*fsds[m1][3*6+2];
		gsfs[m][ 6*10+8] = QC[2]*gsds[m][ 6*6+0] + WQ[2]*gsds[m1][ 6*6+0];
		gsfs[m][ 6*10+9] = QC[0]*gsds[m][ 6*6+4] + WQ[0]*gsds[m1][ 6*6+4]
				 + ze22*fsds[m1][6*6+4];
		gsfs[m][ 7*10+0] = QC[0]*gsds[m][ 7*6+0] + WQ[0]*gsds[m1][ 7*6+0]
				 + eta*tmp.gsps[ 7*3+0];
		gsfs[m][ 7*10+1] = QC[1]*gsds[m][ 7*6+1] + WQ[1]*gsds[m1][ 7*6+1]
				 + eta*tmp.gsps[ 7*3+1] + ze22*fsds[m1][7*6+1];
		gsfs[m][ 7*10+2] = QC[2]*gsds[m][ 7*6+2] + WQ[2]*gsds[m1][ 7*6+2]
				 + eta*tmp.gsps[ 7*3+2] + ze22*fsds[m1][4*6+2];
		gsfs[m][ 7*10+3] = QC[0]*gsds[m][ 7*6+3] + WQ[0]*gsds[m1][ 7*6+3]
				 + eta2*tmp.gsps[ 7*3+1];
		gsfs[m][ 7*10+4] = QC[1]*gsds[m][ 7*6+4] + WQ[1]*gsds[m1][ 7*6+4]
				 + eta2*tmp.gsps[ 7*3+2] + ze22*fsds[m1][7*6+4];
		gsfs[m][ 7*10+5] = QC[2]*gsds[m][ 7*6+5] + WQ[2]*gsds[m1][ 7*6+5]
				 + eta2*tmp.gsps[ 7*3+0] + ze22*fsds[m1][4*6+5];
		gsfs[m][ 7*10+6] = QC[0]*gsds[m][ 7*6+1] + WQ[0]*gsds[m1][ 7*6+1];
		gsfs[m][ 7*10+7] = QC[1]*gsds[m][ 7*6+2] + WQ[1]*gsds[m1][ 7*6+2]
				 + ze22*fsds[m1][7*6+2];
		gsfs[m][ 7*10+8] = QC[2]*gsds[m][ 7*6+0] + WQ[2]*gsds[m1][ 7*6+0]
				 + ze22*fsds[m1][4*6+0];
		gsfs[m][ 7*10+9] = QC[0]*gsds[m][ 7*6+4] + WQ[0]*gsds[m1][ 7*6+4];
		gsfs[m][ 8*10+0] = QC[0]*gsds[m][ 8*6+0] + WQ[0]*gsds[m1][ 8*6+0]
				 + eta*tmp.gsps[ 8*3+0] + ze22*fsds[m1][5*6+0];
		gsfs[m][ 8*10+1] = QC[1]*gsds[m][ 8*6+1] + WQ[1]*gsds[m1][ 8*6+1]
				 + eta*tmp.gsps[ 8*3+1];
		gsfs[m][ 8*10+2] = QC[2]*gsds[m][ 8*6+2] + WQ[2]*gsds[m1][ 8*6+2]
				 + eta*tmp.gsps[ 8*3+2] + ze22*fsds[m1][8*6+2];
		gsfs[m][ 8*10+3] = QC[0]*gsds[m][ 8*6+3] + WQ[0]*gsds[m1][ 8*6+3]
				 + eta2*tmp.gsps[ 8*3+1] + ze22*fsds[m1][5*6+3];
		gsfs[m][ 8*10+4] = QC[1]*gsds[m][ 8*6+4] + WQ[1]*gsds[m1][ 8*6+4]
				 + eta2*tmp.gsps[ 8*3+2];
		gsfs[m][ 8*10+5] = QC[2]*gsds[m][ 8*6+5] + WQ[2]*gsds[m1][ 8*6+5]
				 + eta2*tmp.gsps[ 8*3+0] + ze22*fsds[m1][8*6+5];
		gsfs[m][ 8*10+6] = QC[0]*gsds[m][ 8*6+1] + WQ[0]*gsds[m1][ 8*6+1]
				 + ze22*fsds[m1][5*6+1];
		gsfs[m][ 8*10+7] = QC[1]*gsds[m][ 8*6+2] + WQ[1]*gsds[m1][ 8*6+2];
		gsfs[m][ 8*10+8] = QC[2]*gsds[m][ 8*6+0] + WQ[2]*gsds[m1][ 8*6+0]
				 + ze22*fsds[m1][8*6+0];
		gsfs[m][ 8*10+9] = QC[0]*gsds[m][ 8*6+4] + WQ[0]*gsds[m1][ 8*6+4]
				 + ze22*fsds[m1][5*6+4];
		gsfs[m][ 9*10+0] = QC[0]*gsds[m][ 9*6+0] + WQ[0]*gsds[m1][ 9*6+0]
				 + eta*tmp.gsps[ 9*3+0] + ze2*fsds[m1][1*6+0];
		gsfs[m][ 9*10+1] = QC[1]*gsds[m][ 9*6+1] + WQ[1]*gsds[m1][ 9*6+1]
				 + eta*tmp.gsps[ 9*3+1] + ze23*fsds[m1][6*6+1];
		gsfs[m][ 9*10+2] = QC[2]*gsds[m][ 9*6+2] + WQ[2]*gsds[m1][ 9*6+2]
				 + eta*tmp.gsps[ 9*3+2];
		gsfs[m][ 9*10+3] = QC[0]*gsds[m][ 9*6+3] + WQ[0]*gsds[m1][ 9*6+3]
				 + eta2*tmp.gsps[ 9*3+1] + ze2*fsds[m1][1*6+3];
		gsfs[m][ 9*10+4] = QC[1]*gsds[m][ 9*6+4] + WQ[1]*gsds[m1][ 9*6+4]
				 + eta2*tmp.gsps[ 9*3+2] + ze23*fsds[m1][6*6+4];
		gsfs[m][ 9*10+5] = QC[2]*gsds[m][ 9*6+5] + WQ[2]*gsds[m1][ 9*6+5]
				 + eta2*tmp.gsps[ 9*3+0];
		gsfs[m][ 9*10+6] = QC[0]*gsds[m][ 9*6+1] + WQ[0]*gsds[m1][ 9*6+1]
				 + ze2*fsds[m1][1*6+1];
		gsfs[m][ 9*10+7] = QC[1]*gsds[m][ 9*6+2] + WQ[1]*gsds[m1][ 9*6+2]
				 + ze23*fsds[m1][6*6+2];
		gsfs[m][ 9*10+8] = QC[2]*gsds[m][ 9*6+0] + WQ[2]*gsds[m1][ 9*6+0];
		gsfs[m][ 9*10+9] = QC[0]*gsds[m][ 9*6+4] + WQ[0]*gsds[m1][ 9*6+4]
				 + ze2*fsds[m1][1*6+4];
		gsfs[m][10*10+0] = QC[0]*gsds[m][10*6+0] + WQ[0]*gsds[m1][10*6+0]
				 + eta*tmp.gsps[10*3+0];
		gsfs[m][10*10+1] = QC[1]*gsds[m][10*6+1] + WQ[1]*gsds[m1][10*6+1]
				 + eta*tmp.gsps[10*3+1] + ze2*fsds[m1][2*6+1];
		gsfs[m][10*10+2] = QC[2]*gsds[m][10*6+2] + WQ[2]*gsds[m1][10*6+2]
				 + eta*tmp.gsps[10*3+2] + ze23*fsds[m1][7*6+2];
		gsfs[m][10*10+3] = QC[0]*gsds[m][10*6+3] + WQ[0]*gsds[m1][10*6+3]
				 + eta2*tmp.gsps[10*3+1];
		gsfs[m][10*10+4] = QC[1]*gsds[m][10*6+4] + WQ[1]*gsds[m1][10*6+4]
				 + eta2*tmp.gsps[10*3+2] + ze2*fsds[m1][2*6+4];
		gsfs[m][10*10+5] = QC[2]*gsds[m][10*6+5] + WQ[2]*gsds[m1][10*6+5]
				 + eta2*tmp.gsps[10*3+0] + ze23*fsds[m1][7*6+5];
		gsfs[m][10*10+6] = QC[0]*gsds[m][10*6+1] + WQ[0]*gsds[m1][10*6+1];
		gsfs[m][10*10+7] = QC[1]*gsds[m][10*6+2] + WQ[1]*gsds[m1][10*6+2]
				 + ze2*fsds[m1][2*6+2];
		gsfs[m][10*10+8] = QC[2]*gsds[m][10*6+0] + WQ[2]*gsds[m1][10*6+0]
				 + ze23*fsds[m1][7*6+0];
		gsfs[m][10*10+9] = QC[0]*gsds[m][10*6+4] + WQ[0]*gsds[m1][10*6+4];
		gsfs[m][11*10+0] = QC[0]*gsds[m][11*6+0] + WQ[0]*gsds[m1][11*6+0]
				 + eta*tmp.gsps[11*3+0] + ze23*fsds[m1][8*6+0];
		gsfs[m][11*10+1] = QC[1]*gsds[m][11*6+1] + WQ[1]*gsds[m1][11*6+1]
				 + eta*tmp.gsps[11*3+1];
		gsfs[m][11*10+2] = QC[2]*gsds[m][11*6+2] + WQ[2]*gsds[m1][11*6+2]
				 + eta*tmp.gsps[11*3+2] + ze2*fsds[m1][0*6+2];
		gsfs[m][11*10+3] = QC[0]*gsds[m][11*6+3] + WQ[0]*gsds[m1][11*6+3]
				 + eta2*tmp.gsps[11*3+1] + ze23*fsds[m1][8*6+3];
		gsfs[m][11*10+4] = QC[1]*gsds[m][11*6+4] + WQ[1]*gsds[m1][11*6+4]
				 + eta2*tmp.gsps[11*3+2];
		gsfs[m][11*10+5] = QC[2]*gsds[m][11*6+5] + WQ[2]*gsds[m1][11*6+5]
				 + eta2*tmp.gsps[11*3+0] + ze2*fsds[m1][0*6+5];
		gsfs[m][11*10+6] = QC[0]*gsds[m][11*6+1] + WQ[0]*gsds[m1][11*6+1]
				 + ze23*fsds[m1][8*6+1];
		gsfs[m][11*10+7] = QC[1]*gsds[m][11*6+2] + WQ[1]*gsds[m1][11*6+2];
		gsfs[m][11*10+8] = QC[2]*gsds[m][11*6+0] + WQ[2]*gsds[m1][11*6+0]
				 + ze2*fsds[m1][0*6+0];
		gsfs[m][11*10+9] = QC[0]*gsds[m][11*6+4] + WQ[0]*gsds[m1][11*6+4]
				 + ze23*fsds[m1][8*6+4];
		gsfs[m][12*10+0] = QC[0]*gsds[m][12*6+0] + WQ[0]*gsds[m1][12*6+0]
				 + eta*tmp.gsps[12*3+0] + ze22*fsds[m1][9*6+0];
		gsfs[m][12*10+1] = QC[1]*gsds[m][12*6+1] + WQ[1]*gsds[m1][12*6+1]
				 + eta*tmp.gsps[12*3+1] + ze2*fsds[m1][8*6+1];
		gsfs[m][12*10+2] = QC[2]*gsds[m][12*6+2] + WQ[2]*gsds[m1][12*6+2]
				 + eta*tmp.gsps[12*3+2] + ze2*fsds[m1][3*6+2];
		gsfs[m][12*10+3] = QC[0]*gsds[m][12*6+3] + WQ[0]*gsds[m1][12*6+3]
				 + eta2*tmp.gsps[12*3+1] + ze22*fsds[m1][9*6+3];
		gsfs[m][12*10+4] = QC[1]*gsds[m][12*6+4] + WQ[1]*gsds[m1][12*6+4]
				 + eta2*tmp.gsps[12*3+2] + ze2*fsds[m1][8*6+4];
		gsfs[m][12*10+5] = QC[2]*gsds[m][12*6+5] + WQ[2]*gsds[m1][12*6+5]
				 + eta2*tmp.gsps[12*3+0] + ze2*fsds[m1][3*6+5];
		gsfs[m][12*10+6] = QC[0]*gsds[m][12*6+1] + WQ[0]*gsds[m1][12*6+1]
				 + ze22*fsds[m1][9*6+1];
		gsfs[m][12*10+7] = QC[1]*gsds[m][12*6+2] + WQ[1]*gsds[m1][12*6+2]
				 + ze2*fsds[m1][8*6+2];
		gsfs[m][12*10+8] = QC[2]*gsds[m][12*6+0] + WQ[2]*gsds[m1][12*6+0]
				 + ze2*fsds[m1][3*6+0];
		gsfs[m][12*10+9] = QC[0]*gsds[m][12*6+4] + WQ[0]*gsds[m1][12*6+4]
				 + ze22*fsds[m1][9*6+4];
		gsfs[m][13*10+0] = QC[0]*gsds[m][13*6+0] + WQ[0]*gsds[m1][13*6+0]
				 + eta*tmp.gsps[13*3+0] + ze2*fsds[m1][4*6+0];
		gsfs[m][13*10+1] = QC[1]*gsds[m][13*6+1] + WQ[1]*gsds[m1][13*6+1]
				 + eta*tmp.gsps[13*3+1] + ze22*fsds[m1][9*6+1];
		gsfs[m][13*10+2] = QC[2]*gsds[m][13*6+2] + WQ[2]*gsds[m1][13*6+2]
				 + eta*tmp.gsps[13*3+2] + ze2*fsds[m1][6*6+2];
		gsfs[m][13*10+3] = QC[0]*gsds[m][13*6+3] + WQ[0]*gsds[m1][13*6+3]
				 + eta2*tmp.gsps[13*3+1] + ze2*fsds[m1][4*6+3];
		gsfs[m][13*10+4] = QC[1]*gsds[m][13*6+4] + WQ[1]*gsds[m1][13*6+4]
				 + eta2*tmp.gsps[13*3+2] + ze22*fsds[m1][9*6+4];
		gsfs[m][13*10+5] = QC[2]*gsds[m][13*6+5] + WQ[2]*gsds[m1][13*6+5]
				 + eta2*tmp.gsps[13*3+0] + ze2*fsds[m1][6*6+5];
		gsfs[m][13*10+6] = QC[0]*gsds[m][13*6+1] + WQ[0]*gsds[m1][13*6+1]
				 + ze2*fsds[m1][4*6+1];
		gsfs[m][13*10+7] = QC[1]*gsds[m][13*6+2] + WQ[1]*gsds[m1][13*6+2]
				 + ze22*fsds[m1][9*6+2];
		gsfs[m][13*10+8] = QC[2]*gsds[m][13*6+0] + WQ[2]*gsds[m1][13*6+0]
				 + ze2*fsds[m1][6*6+0];
		gsfs[m][13*10+9] = QC[0]*gsds[m][13*6+4] + WQ[0]*gsds[m1][13*6+4]
				 + ze2*fsds[m1][4*6+4];
		gsfs[m][14*10+0] = QC[0]*gsds[m][14*6+0] + WQ[0]*gsds[m1][14*6+0]
				 + eta*tmp.gsps[14*3+0] + ze2*fsds[m1][7*6+0];
		gsfs[m][14*10+1] = QC[1]*gsds[m][14*6+1] + WQ[1]*gsds[m1][14*6+1]
				 + eta*tmp.gsps[14*3+1] + ze2*fsds[m1][5*6+1];
		gsfs[m][14*10+2] = QC[2]*gsds[m][14*6+2] + WQ[2]*gsds[m1][14*6+2]
				 + eta*tmp.gsps[14*3+2] + ze22*fsds[m1][9*6+2];
		gsfs[m][14*10+3] = QC[0]*gsds[m][14*6+3] + WQ[0]*gsds[m1][14*6+3]
				 + eta2*tmp.gsps[14*3+1] + ze2*fsds[m1][7*6+3];
		gsfs[m][14*10+4] = QC[1]*gsds[m][14*6+4] + WQ[1]*gsds[m1][14*6+4]
				 + eta2*tmp.gsps[14*3+2] + ze2*fsds[m1][5*6+4];
		gsfs[m][14*10+5] = QC[2]*gsds[m][14*6+5] + WQ[2]*gsds[m1][14*6+5]
				 + eta2*tmp.gsps[14*3+0] + ze22*fsds[m1][9*6+5];
		gsfs[m][14*10+6] = QC[0]*gsds[m][14*6+1] + WQ[0]*gsds[m1][14*6+1]
				 + ze2*fsds[m1][7*6+1];
		gsfs[m][14*10+7] = QC[1]*gsds[m][14*6+2] + WQ[1]*gsds[m1][14*6+2]
				 + ze2*fsds[m1][5*6+2];
		gsfs[m][14*10+8] = QC[2]*gsds[m][14*6+0] + WQ[2]*gsds[m1][14*6+0]
				 + ze22*fsds[m1][9*6+0];
		gsfs[m][14*10+9] = QC[0]*gsds[m][14*6+4] + WQ[0]*gsds[m1][14*6+4]
				 + ze2*fsds[m1][7*6+4];
	    }	// end gsfs loop
	    for (i=0; i<15*10; i++) GSFS[i] += gsfs[0][i];
	    // dsgs (m=0)
	    for (i=0; i<6*6; i++) tmp.dsds[i]=dsds[0][i]-re*dsds[1][i];
	    DSGS[0*15+ 0]+=QC[0]*dsfs[0][0*10+0] + WQ[0]*dsfs[1][0*10+0]
			  + eta23*tmp.dsds[0*6+0] + ze22*psfs[0*10+0];
	    DSGS[0*15+ 1]+=QC[1]*dsfs[0][0*10+1] + WQ[1]*dsfs[1][0*10+1]
			  + eta23*tmp.dsds[0*6+1];
	    DSGS[0*15+ 2]+=QC[2]*dsfs[0][0*10+2] + WQ[2]*dsfs[1][0*10+2]
			  + eta23*tmp.dsds[0*6+2];
	    DSGS[0*15+ 3]+=QC[0]*dsfs[0][0*10+3] + WQ[0]*dsfs[1][0*10+3]
			  + eta*tmp.dsds[0*6+3] + ze22*psfs[0*10+3];
	    DSGS[0*15+ 4]+=QC[1]*dsfs[0][0*10+4] + WQ[1]*dsfs[1][0*10+4]
			  + eta*tmp.dsds[0*6+4];
	    DSGS[0*15+ 5]+=QC[2]*dsfs[0][0*10+5] + WQ[2]*dsfs[1][0*10+5]
			  + eta*tmp.dsds[0*6+5];
	    DSGS[0*15+ 6]+=QC[0]*dsfs[0][0*10+6] + WQ[0]*dsfs[1][0*10+6]
			  + eta2*tmp.dsds[0*6+1] + ze22*psfs[0*10+6];
	    DSGS[0*15+ 7]+=QC[1]*dsfs[0][0*10+7] + WQ[1]*dsfs[1][0*10+7]
			  + eta2*tmp.dsds[0*6+2];
	    DSGS[0*15+ 8]+=QC[2]*dsfs[0][0*10+8] + WQ[2]*dsfs[1][0*10+8]
			  + eta2*tmp.dsds[0*6+0];
	    DSGS[0*15+ 9]+=QC[0]*dsfs[0][0*10+1] + WQ[0]*dsfs[1][0*10+1]
			  + ze22*psfs[0*10+1];
	    DSGS[0*15+10]+=QC[1]*dsfs[0][0*10+2] + WQ[1]*dsfs[1][0*10+2];
	    DSGS[0*15+11]+=QC[2]*dsfs[0][0*10+0] + WQ[2]*dsfs[1][0*10+0];
	    DSGS[0*15+12]+=QC[0]*dsfs[0][0*10+9] + WQ[0]*dsfs[1][0*10+9]
			  + eta2*tmp.dsds[0*6+4] + ze22*psfs[0*10+9];
	    DSGS[0*15+13]+=QC[1]*dsfs[0][0*10+9] + WQ[1]*dsfs[1][0*10+9]
			  + eta2*tmp.dsds[0*6+5];
	    DSGS[0*15+14]+=QC[2]*dsfs[0][0*10+9] + WQ[2]*dsfs[1][0*10+9]
			  + eta2*tmp.dsds[0*6+3];
	    DSGS[1*15+ 0]+=QC[0]*dsfs[0][1*10+0] + WQ[0]*dsfs[1][1*10+0]
			  + eta23*tmp.dsds[1*6+0];
	    DSGS[1*15+ 1]+=QC[1]*dsfs[0][1*10+1] + WQ[1]*dsfs[1][1*10+1]
			  + eta23*tmp.dsds[1*6+1] + ze22*psfs[1*10+1];
	    DSGS[1*15+ 2]+=QC[2]*dsfs[0][1*10+2] + WQ[2]*dsfs[1][1*10+2]
			  + eta23*tmp.dsds[1*6+2];
	    DSGS[1*15+ 3]+=QC[0]*dsfs[0][1*10+3] + WQ[0]*dsfs[1][1*10+3]
			  + eta*tmp.dsds[1*6+3];
	    DSGS[1*15+ 4]+=QC[1]*dsfs[0][1*10+4] + WQ[1]*dsfs[1][1*10+4]
			  + eta*tmp.dsds[1*6+4] + ze22*psfs[1*10+4];
	    DSGS[1*15+ 5]+=QC[2]*dsfs[0][1*10+5] + WQ[2]*dsfs[1][1*10+5]
			  + eta*tmp.dsds[1*6+5];
	    DSGS[1*15+ 6]+=QC[0]*dsfs[0][1*10+6] + WQ[0]*dsfs[1][1*10+6]
			  + eta2*tmp.dsds[1*6+1];
	    DSGS[1*15+ 7]+=QC[1]*dsfs[0][1*10+7] + WQ[1]*dsfs[1][1*10+7]
			  + eta2*tmp.dsds[1*6+2] + ze22*psfs[1*10+7];
	    DSGS[1*15+ 8]+=QC[2]*dsfs[0][1*10+8] + WQ[2]*dsfs[1][1*10+8]
			  + eta2*tmp.dsds[1*6+0];
	    DSGS[1*15+ 9]+=QC[0]*dsfs[0][1*10+1] + WQ[0]*dsfs[1][1*10+1];
	    DSGS[1*15+10]+=QC[1]*dsfs[0][1*10+2] + WQ[1]*dsfs[1][1*10+2]
			  + ze22*psfs[1*10+2];
	    DSGS[1*15+11]+=QC[2]*dsfs[0][1*10+0] + WQ[2]*dsfs[1][1*10+0];
	    DSGS[1*15+12]+=QC[0]*dsfs[0][1*10+9] + WQ[0]*dsfs[1][1*10+9]
			  + eta2*tmp.dsds[1*6+4];
	    DSGS[1*15+13]+=QC[1]*dsfs[0][1*10+9] + WQ[1]*dsfs[1][1*10+9]
			  + eta2*tmp.dsds[1*6+5] + ze22*psfs[1*10+9];
	    DSGS[1*15+14]+=QC[2]*dsfs[0][1*10+9] + WQ[2]*dsfs[1][1*10+9]
			  + eta2*tmp.dsds[1*6+3];
	    DSGS[2*15+ 0]+=QC[0]*dsfs[0][2*10+0] + WQ[0]*dsfs[1][2*10+0]
			  + eta23*tmp.dsds[2*6+0];
	    DSGS[2*15+ 1]+=QC[1]*dsfs[0][2*10+1] + WQ[1]*dsfs[1][2*10+1]
			  + eta23*tmp.dsds[2*6+1];
	    DSGS[2*15+ 2]+=QC[2]*dsfs[0][2*10+2] + WQ[2]*dsfs[1][2*10+2]
			  + eta23*tmp.dsds[2*6+2] + ze22*psfs[2*10+2];
	    DSGS[2*15+ 3]+=QC[0]*dsfs[0][2*10+3] + WQ[0]*dsfs[1][2*10+3]
			  + eta*tmp.dsds[2*6+3];
	    DSGS[2*15+ 4]+=QC[1]*dsfs[0][2*10+4] + WQ[1]*dsfs[1][2*10+4]
			  + eta*tmp.dsds[2*6+4];
	    DSGS[2*15+ 5]+=QC[2]*dsfs[0][2*10+5] + WQ[2]*dsfs[1][2*10+5]
			  + eta*tmp.dsds[2*6+5] + ze22*psfs[2*10+5];
	    DSGS[2*15+ 6]+=QC[0]*dsfs[0][2*10+6] + WQ[0]*dsfs[1][2*10+6]
			  + eta2*tmp.dsds[2*6+1];
	    DSGS[2*15+ 7]+=QC[1]*dsfs[0][2*10+7] + WQ[1]*dsfs[1][2*10+7]
			  + eta2*tmp.dsds[2*6+2];
	    DSGS[2*15+ 8]+=QC[2]*dsfs[0][2*10+8] + WQ[2]*dsfs[1][2*10+8]
			  + eta2*tmp.dsds[2*6+0] + ze22*psfs[2*10+8];
	    DSGS[2*15+ 9]+=QC[0]*dsfs[0][2*10+1] + WQ[0]*dsfs[1][2*10+1];
	    DSGS[2*15+10]+=QC[1]*dsfs[0][2*10+2] + WQ[1]*dsfs[1][2*10+2];
	    DSGS[2*15+11]+=QC[2]*dsfs[0][2*10+0] + WQ[2]*dsfs[1][2*10+0]
			  + ze22*psfs[2*10+0];
	    DSGS[2*15+12]+=QC[0]*dsfs[0][2*10+9] + WQ[0]*dsfs[1][2*10+9]
			  + eta2*tmp.dsds[2*6+4];
	    DSGS[2*15+13]+=QC[1]*dsfs[0][2*10+9] + WQ[1]*dsfs[1][2*10+9]
			  + eta2*tmp.dsds[2*6+5];
	    DSGS[2*15+14]+=QC[2]*dsfs[0][2*10+9] + WQ[2]*dsfs[1][2*10+9]
			  + eta2*tmp.dsds[2*6+3] + ze22*psfs[2*10+9];
	    DSGS[3*15+ 0]+=QC[0]*dsfs[0][3*10+0] + WQ[0]*dsfs[1][3*10+0]
			  + eta23*tmp.dsds[3*6+0] + ze2*psfs[1*10+0];
	    DSGS[3*15+ 1]+=QC[1]*dsfs[0][3*10+1] + WQ[1]*dsfs[1][3*10+1]
			  + eta23*tmp.dsds[3*6+1] + ze2*psfs[0*10+1];
	    DSGS[3*15+ 2]+=QC[2]*dsfs[0][3*10+2] + WQ[2]*dsfs[1][3*10+2]
			  + eta23*tmp.dsds[3*6+2];
	    DSGS[3*15+ 3]+=QC[0]*dsfs[0][3*10+3] + WQ[0]*dsfs[1][3*10+3]
			  + eta*tmp.dsds[3*6+3] + ze2*psfs[1*10+3];
	    DSGS[3*15+ 4]+=QC[1]*dsfs[0][3*10+4] + WQ[1]*dsfs[1][3*10+4]
			  + eta*tmp.dsds[3*6+4] + ze2*psfs[0*10+4];
	    DSGS[3*15+ 5]+=QC[2]*dsfs[0][3*10+5] + WQ[2]*dsfs[1][3*10+5]
			  + eta*tmp.dsds[3*6+5];
	    DSGS[3*15+ 6]+=QC[0]*dsfs[0][3*10+6] + WQ[0]*dsfs[1][3*10+6]
			  + eta2*tmp.dsds[3*6+1] + ze2*psfs[1*10+6];
	    DSGS[3*15+ 7]+=QC[1]*dsfs[0][3*10+7] + WQ[1]*dsfs[1][3*10+7]
			  + eta2*tmp.dsds[3*6+2] + ze2*psfs[0*10+7];
	    DSGS[3*15+ 8]+=QC[2]*dsfs[0][3*10+8] + WQ[2]*dsfs[1][3*10+8]
			  + eta2*tmp.dsds[3*6+0];
	    DSGS[3*15+ 9]+=QC[0]*dsfs[0][3*10+1] + WQ[0]*dsfs[1][3*10+1]
			  + ze2*psfs[1*10+1];
	    DSGS[3*15+10]+=QC[1]*dsfs[0][3*10+2] + WQ[1]*dsfs[1][3*10+2]
			  + ze2*psfs[0*10+2];
	    DSGS[3*15+11]+=QC[2]*dsfs[0][3*10+0] + WQ[2]*dsfs[1][3*10+0];
	    DSGS[3*15+12]+=QC[0]*dsfs[0][3*10+9] + WQ[0]*dsfs[1][3*10+9]
			  + eta2*tmp.dsds[3*6+4] + ze2*psfs[1*10+9];
	    DSGS[3*15+13]+=QC[1]*dsfs[0][3*10+9] + WQ[1]*dsfs[1][3*10+9]
			  + eta2*tmp.dsds[3*6+5] + ze2*psfs[0*10+9];
	    DSGS[3*15+14]+=QC[2]*dsfs[0][3*10+9] + WQ[2]*dsfs[1][3*10+9]
			  + eta2*tmp.dsds[3*6+3];
	    DSGS[4*15+ 0]+=QC[0]*dsfs[0][4*10+0] + WQ[0]*dsfs[1][4*10+0]
			  + eta23*tmp.dsds[4*6+0];
	    DSGS[4*15+ 1]+=QC[1]*dsfs[0][4*10+1] + WQ[1]*dsfs[1][4*10+1]
			  + eta23*tmp.dsds[4*6+1] + ze2*psfs[2*10+1];
	    DSGS[4*15+ 2]+=QC[2]*dsfs[0][4*10+2] + WQ[2]*dsfs[1][4*10+2]
			  + eta23*tmp.dsds[4*6+2] + ze2*psfs[1*10+2];
	    DSGS[4*15+ 3]+=QC[0]*dsfs[0][4*10+3] + WQ[0]*dsfs[1][4*10+3]
			  + eta*tmp.dsds[4*6+3];
	    DSGS[4*15+ 4]+=QC[1]*dsfs[0][4*10+4] + WQ[1]*dsfs[1][4*10+4]
			  + eta*tmp.dsds[4*6+4] + ze2*psfs[2*10+4];
	    DSGS[4*15+ 5]+=QC[2]*dsfs[0][4*10+5] + WQ[2]*dsfs[1][4*10+5]
			  + eta*tmp.dsds[4*6+5] + ze2*psfs[1*10+5];
	    DSGS[4*15+ 6]+=QC[0]*dsfs[0][4*10+6] + WQ[0]*dsfs[1][4*10+6]
			  + eta2*tmp.dsds[4*6+1];
	    DSGS[4*15+ 7]+=QC[1]*dsfs[0][4*10+7] + WQ[1]*dsfs[1][4*10+7]
			  + eta2*tmp.dsds[4*6+2] + ze2*psfs[2*10+7];
	    DSGS[4*15+ 8]+=QC[2]*dsfs[0][4*10+8] + WQ[2]*dsfs[1][4*10+8]
			  + eta2*tmp.dsds[4*6+0] + ze2*psfs[1*10+8];
	    DSGS[4*15+ 9]+=QC[0]*dsfs[0][4*10+1] + WQ[0]*dsfs[1][4*10+1];
	    DSGS[4*15+10]+=QC[1]*dsfs[0][4*10+2] + WQ[1]*dsfs[1][4*10+2]
			  + ze2*psfs[2*10+2];
	    DSGS[4*15+11]+=QC[2]*dsfs[0][4*10+0] + WQ[2]*dsfs[1][4*10+0]
			  + ze2*psfs[1*10+0];
	    DSGS[4*15+12]+=QC[0]*dsfs[0][4*10+9] + WQ[0]*dsfs[1][4*10+9]
			  + eta2*tmp.dsds[4*6+4];
	    DSGS[4*15+13]+=QC[1]*dsfs[0][4*10+9] + WQ[1]*dsfs[1][4*10+9]
			  + eta2*tmp.dsds[4*6+5] + ze2*psfs[2*10+9];
	    DSGS[4*15+14]+=QC[2]*dsfs[0][4*10+9] + WQ[2]*dsfs[1][4*10+9]
			  + eta2*tmp.dsds[4*6+3] + ze2*psfs[1*10+9];
	    DSGS[5*15+ 0]+=QC[0]*dsfs[0][5*10+0] + WQ[0]*dsfs[1][5*10+0]
			  + eta23*tmp.dsds[5*6+0] + ze2*psfs[2*10+0];
	    DSGS[5*15+ 1]+=QC[1]*dsfs[0][5*10+1] + WQ[1]*dsfs[1][5*10+1]
			  + eta23*tmp.dsds[5*6+1];
	    DSGS[5*15+ 2]+=QC[2]*dsfs[0][5*10+2] + WQ[2]*dsfs[1][5*10+2]
			  + eta23*tmp.dsds[5*6+2] + ze2*psfs[0*10+2];
	    DSGS[5*15+ 3]+=QC[0]*dsfs[0][5*10+3] + WQ[0]*dsfs[1][5*10+3]
			  + eta*tmp.dsds[5*6+3] + ze2*psfs[2*10+3];
	    DSGS[5*15+ 4]+=QC[1]*dsfs[0][5*10+4] + WQ[1]*dsfs[1][5*10+4]
			  + eta*tmp.dsds[5*6+4];
	    DSGS[5*15+ 5]+=QC[2]*dsfs[0][5*10+5] + WQ[2]*dsfs[1][5*10+5]
			  + eta*tmp.dsds[5*6+5] + ze2*psfs[0*10+5];
	    DSGS[5*15+ 6]+=QC[0]*dsfs[0][5*10+6] + WQ[0]*dsfs[1][5*10+6]
			  + eta2*tmp.dsds[5*6+1] + ze2*psfs[2*10+6];
	    DSGS[5*15+ 7]+=QC[1]*dsfs[0][5*10+7] + WQ[1]*dsfs[1][5*10+7]
			  + eta2*tmp.dsds[5*6+2];
	    DSGS[5*15+ 8]+=QC[2]*dsfs[0][5*10+8] + WQ[2]*dsfs[1][5*10+8]
			  + eta2*tmp.dsds[5*6+0] + ze2*psfs[0*10+8];
	    DSGS[5*15+ 9]+=QC[0]*dsfs[0][5*10+1] + WQ[0]*dsfs[1][5*10+1]
			  + ze2*psfs[2*10+1];
	    DSGS[5*15+10]+=QC[1]*dsfs[0][5*10+2] + WQ[1]*dsfs[1][5*10+2];
	    DSGS[5*15+11]+=QC[2]*dsfs[0][5*10+0] + WQ[2]*dsfs[1][5*10+0]
			  + ze2*psfs[0*10+0];
	    DSGS[5*15+12]+=QC[0]*dsfs[0][5*10+9] + WQ[0]*dsfs[1][5*10+9]
			  + eta2*tmp.dsds[5*6+4] + ze2*psfs[2*10+9];
	    DSGS[5*15+13]+=QC[1]*dsfs[0][5*10+9] + WQ[1]*dsfs[1][5*10+9]
			  + eta2*tmp.dsds[5*6+5];
	    DSGS[5*15+14]+=QC[2]*dsfs[0][5*10+9] + WQ[2]*dsfs[1][5*10+9]
			  + eta2*tmp.dsds[5*6+3] + ze2*psfs[0*10+9];
	    // fsgs (m=0)
	    for (i=0; i<10*6; i++) tmp.fsds[i]=fsds[0][i]-re*fsds[1][i];
	    FSGS[0*15+ 0]+=QC[0]*fsfs[0][0*10+0] + WQ[0]*fsfs[1][0*10+0]
			  + eta23*tmp.fsds[0*6+0] + ze23*dsfs[1][0*10+0];
	    FSGS[0*15+ 1]+=QC[1]*fsfs[0][0*10+1] + WQ[1]*fsfs[1][0*10+1]
			  + eta23*tmp.fsds[0*6+1];
	    FSGS[0*15+ 2]+=QC[2]*fsfs[0][0*10+2] + WQ[2]*fsfs[1][0*10+2]
			  + eta23*tmp.fsds[0*6+2];
	    FSGS[0*15+ 3]+=QC[0]*fsfs[0][0*10+3] + WQ[0]*fsfs[1][0*10+3]
			  + eta*tmp.fsds[0*6+3] + ze23*dsfs[1][0*10+3];
	    FSGS[0*15+ 4]+=QC[1]*fsfs[0][0*10+4] + WQ[1]*fsfs[1][0*10+4]
			  + eta*tmp.fsds[0*6+4];
	    FSGS[0*15+ 5]+=QC[2]*fsfs[0][0*10+5] + WQ[2]*fsfs[1][0*10+5]
			  + eta*tmp.fsds[0*6+5];
	    FSGS[0*15+ 6]+=QC[0]*fsfs[0][0*10+6] + WQ[0]*fsfs[1][0*10+6]
			  + eta2*tmp.fsds[0*6+1] + ze23*dsfs[1][0*10+6];
	    FSGS[0*15+ 7]+=QC[1]*fsfs[0][0*10+7] + WQ[1]*fsfs[1][0*10+7]
			  + eta2*tmp.fsds[0*6+2];
	    FSGS[0*15+ 8]+=QC[2]*fsfs[0][0*10+8] + WQ[2]*fsfs[1][0*10+8]
			  + eta2*tmp.fsds[0*6+0];
	    FSGS[0*15+ 9]+=QC[0]*fsfs[0][0*10+1] + WQ[0]*fsfs[1][0*10+1]
			  + ze23*dsfs[1][0*10+1];
	    FSGS[0*15+10]+=QC[1]*fsfs[0][0*10+2] + WQ[1]*fsfs[1][0*10+2];
	    FSGS[0*15+11]+=QC[2]*fsfs[0][0*10+0] + WQ[2]*fsfs[1][0*10+0];
	    FSGS[0*15+12]+=QC[0]*fsfs[0][0*10+9] + WQ[0]*fsfs[1][0*10+9]
			  + eta2*tmp.fsds[0*6+4] + ze23*dsfs[1][0*10+9];
	    FSGS[0*15+13]+=QC[1]*fsfs[0][0*10+9] + WQ[1]*fsfs[1][0*10+9]
			  + eta2*tmp.fsds[0*6+5];
	    FSGS[0*15+14]+=QC[2]*fsfs[0][0*10+9] + WQ[2]*fsfs[1][0*10+9]
			  + eta2*tmp.fsds[0*6+3];
	    FSGS[1*15+ 0]+=QC[0]*fsfs[0][1*10+0] + WQ[0]*fsfs[1][1*10+0]
			  + eta23*tmp.fsds[1*6+0];
	    FSGS[1*15+ 1]+=QC[1]*fsfs[0][1*10+1] + WQ[1]*fsfs[1][1*10+1]
			  + eta23*tmp.fsds[1*6+1] + ze23*dsfs[1][1*10+1];
	    FSGS[1*15+ 2]+=QC[2]*fsfs[0][1*10+2] + WQ[2]*fsfs[1][1*10+2]
			  + eta23*tmp.fsds[1*6+2];
	    FSGS[1*15+ 3]+=QC[0]*fsfs[0][1*10+3] + WQ[0]*fsfs[1][1*10+3]
			  + eta*tmp.fsds[1*6+3];
	    FSGS[1*15+ 4]+=QC[1]*fsfs[0][1*10+4] + WQ[1]*fsfs[1][1*10+4]
			  + eta*tmp.fsds[1*6+4] + ze23*dsfs[1][1*10+4];
	    FSGS[1*15+ 5]+=QC[2]*fsfs[0][1*10+5] + WQ[2]*fsfs[1][1*10+5]
			  + eta*tmp.fsds[1*6+5];
	    FSGS[1*15+ 6]+=QC[0]*fsfs[0][1*10+6] + WQ[0]*fsfs[1][1*10+6]
			  + eta2*tmp.fsds[1*6+1];
	    FSGS[1*15+ 7]+=QC[1]*fsfs[0][1*10+7] + WQ[1]*fsfs[1][1*10+7]
			  + eta2*tmp.fsds[1*6+2] + ze23*dsfs[1][1*10+7];
	    FSGS[1*15+ 8]+=QC[2]*fsfs[0][1*10+8] + WQ[2]*fsfs[1][1*10+8]
			  + eta2*tmp.fsds[1*6+0];
	    FSGS[1*15+ 9]+=QC[0]*fsfs[0][1*10+1] + WQ[0]*fsfs[1][1*10+1];
	    FSGS[1*15+10]+=QC[1]*fsfs[0][1*10+2] + WQ[1]*fsfs[1][1*10+2]
			  + ze23*dsfs[1][1*10+2];
	    FSGS[1*15+11]+=QC[2]*fsfs[0][1*10+0] + WQ[2]*fsfs[1][1*10+0];
	    FSGS[1*15+12]+=QC[0]*fsfs[0][1*10+9] + WQ[0]*fsfs[1][1*10+9]
			  + eta2*tmp.fsds[1*6+4];
	    FSGS[1*15+13]+=QC[1]*fsfs[0][1*10+9] + WQ[1]*fsfs[1][1*10+9]
			  + eta2*tmp.fsds[1*6+5] + ze23*dsfs[1][1*10+9];
	    FSGS[1*15+14]+=QC[2]*fsfs[0][1*10+9] + WQ[2]*fsfs[1][1*10+9]
			  + eta2*tmp.fsds[1*6+3];
	    FSGS[2*15+ 0]+=QC[0]*fsfs[0][2*10+0] + WQ[0]*fsfs[1][2*10+0]
			  + eta23*tmp.fsds[2*6+0];
	    FSGS[2*15+ 1]+=QC[1]*fsfs[0][2*10+1] + WQ[1]*fsfs[1][2*10+1]
			  + eta23*tmp.fsds[2*6+1];
	    FSGS[2*15+ 2]+=QC[2]*fsfs[0][2*10+2] + WQ[2]*fsfs[1][2*10+2]
			  + eta23*tmp.fsds[2*6+2] + ze23*dsfs[1][2*10+2];
	    FSGS[2*15+ 3]+=QC[0]*fsfs[0][2*10+3] + WQ[0]*fsfs[1][2*10+3]
			  + eta*tmp.fsds[2*6+3];
	    FSGS[2*15+ 4]+=QC[1]*fsfs[0][2*10+4] + WQ[1]*fsfs[1][2*10+4]
			  + eta*tmp.fsds[2*6+4];
	    FSGS[2*15+ 5]+=QC[2]*fsfs[0][2*10+5] + WQ[2]*fsfs[1][2*10+5]
			  + eta*tmp.fsds[2*6+5] + ze23*dsfs[1][2*10+5];
	    FSGS[2*15+ 6]+=QC[0]*fsfs[0][2*10+6] + WQ[0]*fsfs[1][2*10+6]
			  + eta2*tmp.fsds[2*6+1];
	    FSGS[2*15+ 7]+=QC[1]*fsfs[0][2*10+7] + WQ[1]*fsfs[1][2*10+7]
			  + eta2*tmp.fsds[2*6+2];
	    FSGS[2*15+ 8]+=QC[2]*fsfs[0][2*10+8] + WQ[2]*fsfs[1][2*10+8]
			  + eta2*tmp.fsds[2*6+0] + ze23*dsfs[1][2*10+8];
	    FSGS[2*15+ 9]+=QC[0]*fsfs[0][2*10+1] + WQ[0]*fsfs[1][2*10+1];
	    FSGS[2*15+10]+=QC[1]*fsfs[0][2*10+2] + WQ[1]*fsfs[1][2*10+2];
	    FSGS[2*15+11]+=QC[2]*fsfs[0][2*10+0] + WQ[2]*fsfs[1][2*10+0]
			  + ze23*dsfs[1][2*10+0];
	    FSGS[2*15+12]+=QC[0]*fsfs[0][2*10+9] + WQ[0]*fsfs[1][2*10+9]
			  + eta2*tmp.fsds[2*6+4];
	    FSGS[2*15+13]+=QC[1]*fsfs[0][2*10+9] + WQ[1]*fsfs[1][2*10+9]
			  + eta2*tmp.fsds[2*6+5];
	    FSGS[2*15+14]+=QC[2]*fsfs[0][2*10+9] + WQ[2]*fsfs[1][2*10+9]
			  + eta2*tmp.fsds[2*6+3] + ze23*dsfs[1][2*10+9];
	    FSGS[3*15+ 0]+=QC[0]*fsfs[0][3*10+0] + WQ[0]*fsfs[1][3*10+0]
			  + eta23*tmp.fsds[3*6+0] + ze22*dsfs[1][3*10+0];
	    FSGS[3*15+ 1]+=QC[1]*fsfs[0][3*10+1] + WQ[1]*fsfs[1][3*10+1]
			  + eta23*tmp.fsds[3*6+1] + ze2*dsfs[1][0*10+1];
	    FSGS[3*15+ 2]+=QC[2]*fsfs[0][3*10+2] + WQ[2]*fsfs[1][3*10+2]
			  + eta23*tmp.fsds[3*6+2];
	    FSGS[3*15+ 3]+=QC[0]*fsfs[0][3*10+3] + WQ[0]*fsfs[1][3*10+3]
			  + eta*tmp.fsds[3*6+3] + ze22*dsfs[1][3*10+3];
	    FSGS[3*15+ 4]+=QC[1]*fsfs[0][3*10+4] + WQ[1]*fsfs[1][3*10+4]
			  + eta*tmp.fsds[3*6+4] + ze2*dsfs[1][0*10+4];
	    FSGS[3*15+ 5]+=QC[2]*fsfs[0][3*10+5] + WQ[2]*fsfs[1][3*10+5]
			  + eta*tmp.fsds[3*6+5];
	    FSGS[3*15+ 6]+=QC[0]*fsfs[0][3*10+6] + WQ[0]*fsfs[1][3*10+6]
			  + eta2*tmp.fsds[3*6+1] + ze22*dsfs[1][3*10+6];
	    FSGS[3*15+ 7]+=QC[1]*fsfs[0][3*10+7] + WQ[1]*fsfs[1][3*10+7]
			  + eta2*tmp.fsds[3*6+2] + ze2*dsfs[1][0*10+7];
	    FSGS[3*15+ 8]+=QC[2]*fsfs[0][3*10+8] + WQ[2]*fsfs[1][3*10+8]
			  + eta2*tmp.fsds[3*6+0];
	    FSGS[3*15+ 9]+=QC[0]*fsfs[0][3*10+1] + WQ[0]*fsfs[1][3*10+1]
			  + ze22*dsfs[1][3*10+1];
	    FSGS[3*15+10]+=QC[1]*fsfs[0][3*10+2] + WQ[1]*fsfs[1][3*10+2]
			  + ze2*dsfs[1][0*10+2];
	    FSGS[3*15+11]+=QC[2]*fsfs[0][3*10+0] + WQ[2]*fsfs[1][3*10+0];
	    FSGS[3*15+12]+=QC[0]*fsfs[0][3*10+9] + WQ[0]*fsfs[1][3*10+9]
			  + eta2*tmp.fsds[3*6+4] + ze22*dsfs[1][3*10+9];
	    FSGS[3*15+13]+=QC[1]*fsfs[0][3*10+9] + WQ[1]*fsfs[1][3*10+9]
			  + eta2*tmp.fsds[3*6+5] + ze2*dsfs[1][0*10+9];
	    FSGS[3*15+14]+=QC[2]*fsfs[0][3*10+9] + WQ[2]*fsfs[1][3*10+9]
			  + eta2*tmp.fsds[3*6+3];
	    FSGS[4*15+ 0]+=QC[0]*fsfs[0][4*10+0] + WQ[0]*fsfs[1][4*10+0]
			  + eta23*tmp.fsds[4*6+0];
	    FSGS[4*15+ 1]+=QC[1]*fsfs[0][4*10+1] + WQ[1]*fsfs[1][4*10+1]
			  + eta23*tmp.fsds[4*6+1] + ze22*dsfs[1][4*10+1];
	    FSGS[4*15+ 2]+=QC[2]*fsfs[0][4*10+2] + WQ[2]*fsfs[1][4*10+2]
			  + eta23*tmp.fsds[4*6+2] + ze2*dsfs[1][1*10+2];
	    FSGS[4*15+ 3]+=QC[0]*fsfs[0][4*10+3] + WQ[0]*fsfs[1][4*10+3]
			  + eta*tmp.fsds[4*6+3];
	    FSGS[4*15+ 4]+=QC[1]*fsfs[0][4*10+4] + WQ[1]*fsfs[1][4*10+4]
			  + eta*tmp.fsds[4*6+4] + ze22*dsfs[1][4*10+4];
	    FSGS[4*15+ 5]+=QC[2]*fsfs[0][4*10+5] + WQ[2]*fsfs[1][4*10+5]
			  + eta*tmp.fsds[4*6+5] + ze2*dsfs[1][1*10+5];
	    FSGS[4*15+ 6]+=QC[0]*fsfs[0][4*10+6] + WQ[0]*fsfs[1][4*10+6]
			  + eta2*tmp.fsds[4*6+1];
	    FSGS[4*15+ 7]+=QC[1]*fsfs[0][4*10+7] + WQ[1]*fsfs[1][4*10+7]
			  + eta2*tmp.fsds[4*6+2] + ze22*dsfs[1][4*10+7];
	    FSGS[4*15+ 8]+=QC[2]*fsfs[0][4*10+8] + WQ[2]*fsfs[1][4*10+8]
			  + eta2*tmp.fsds[4*6+0] + ze2*dsfs[1][1*10+8];
	    FSGS[4*15+ 9]+=QC[0]*fsfs[0][4*10+1] + WQ[0]*fsfs[1][4*10+1];
	    FSGS[4*15+10]+=QC[1]*fsfs[0][4*10+2] + WQ[1]*fsfs[1][4*10+2]
			  + ze22*dsfs[1][4*10+2];
	    FSGS[4*15+11]+=QC[2]*fsfs[0][4*10+0] + WQ[2]*fsfs[1][4*10+0]
			  + ze2*dsfs[1][1*10+0];
	    FSGS[4*15+12]+=QC[0]*fsfs[0][4*10+9] + WQ[0]*fsfs[1][4*10+9]
			  + eta2*tmp.fsds[4*6+4];
	    FSGS[4*15+13]+=QC[1]*fsfs[0][4*10+9] + WQ[1]*fsfs[1][4*10+9]
			  + eta2*tmp.fsds[4*6+5] + ze22*dsfs[1][4*10+9];
	    FSGS[4*15+14]+=QC[2]*fsfs[0][4*10+9] + WQ[2]*fsfs[1][4*10+9]
			  + eta2*tmp.fsds[4*6+3] + ze2*dsfs[1][1*10+9];
	    FSGS[5*15+ 0]+=QC[0]*fsfs[0][5*10+0] + WQ[0]*fsfs[1][5*10+0]
			  + eta23*tmp.fsds[5*6+0] + ze2*dsfs[1][2*10+0];
	    FSGS[5*15+ 1]+=QC[1]*fsfs[0][5*10+1] + WQ[1]*fsfs[1][5*10+1]
			  + eta23*tmp.fsds[5*6+1];
	    FSGS[5*15+ 2]+=QC[2]*fsfs[0][5*10+2] + WQ[2]*fsfs[1][5*10+2]
			  + eta23*tmp.fsds[5*6+2] + ze22*dsfs[1][5*10+2];
	    FSGS[5*15+ 3]+=QC[0]*fsfs[0][5*10+3] + WQ[0]*fsfs[1][5*10+3]
			  + eta*tmp.fsds[5*6+3] + ze2*dsfs[1][2*10+3];
	    FSGS[5*15+ 4]+=QC[1]*fsfs[0][5*10+4] + WQ[1]*fsfs[1][5*10+4]
			  + eta*tmp.fsds[5*6+4];
	    FSGS[5*15+ 5]+=QC[2]*fsfs[0][5*10+5] + WQ[2]*fsfs[1][5*10+5]
			  + eta*tmp.fsds[5*6+5] + ze22*dsfs[1][5*10+5];
	    FSGS[5*15+ 6]+=QC[0]*fsfs[0][5*10+6] + WQ[0]*fsfs[1][5*10+6]
			  + eta2*tmp.fsds[5*6+1] + ze2*dsfs[1][2*10+6];
	    FSGS[5*15+ 7]+=QC[1]*fsfs[0][5*10+7] + WQ[1]*fsfs[1][5*10+7]
			  + eta2*tmp.fsds[5*6+2];
	    FSGS[5*15+ 8]+=QC[2]*fsfs[0][5*10+8] + WQ[2]*fsfs[1][5*10+8]
			  + eta2*tmp.fsds[5*6+0] + ze22*dsfs[1][5*10+8];
	    FSGS[5*15+ 9]+=QC[0]*fsfs[0][5*10+1] + WQ[0]*fsfs[1][5*10+1]
			  + ze2*dsfs[1][2*10+1];
	    FSGS[5*15+10]+=QC[1]*fsfs[0][5*10+2] + WQ[1]*fsfs[1][5*10+2];
	    FSGS[5*15+11]+=QC[2]*fsfs[0][5*10+0] + WQ[2]*fsfs[1][5*10+0]
			  + ze22*dsfs[1][5*10+0];
	    FSGS[5*15+12]+=QC[0]*fsfs[0][5*10+9] + WQ[0]*fsfs[1][5*10+9]
			  + eta2*tmp.fsds[5*6+4] + ze2*dsfs[1][2*10+9];
	    FSGS[5*15+13]+=QC[1]*fsfs[0][5*10+9] + WQ[1]*fsfs[1][5*10+9]
			  + eta2*tmp.fsds[5*6+5];
	    FSGS[5*15+14]+=QC[2]*fsfs[0][5*10+9] + WQ[2]*fsfs[1][5*10+9]
			  + eta2*tmp.fsds[5*6+3] + ze22*dsfs[1][5*10+9];
	    FSGS[6*15+ 0]+=QC[0]*fsfs[0][6*10+0] + WQ[0]*fsfs[1][6*10+0]
			  + eta23*tmp.fsds[6*6+0] + ze2*dsfs[1][1*10+0];
	    FSGS[6*15+ 1]+=QC[1]*fsfs[0][6*10+1] + WQ[1]*fsfs[1][6*10+1]
			  + eta23*tmp.fsds[6*6+1] + ze22*dsfs[1][3*10+1];
	    FSGS[6*15+ 2]+=QC[2]*fsfs[0][6*10+2] + WQ[2]*fsfs[1][6*10+2]
			  + eta23*tmp.fsds[6*6+2];
	    FSGS[6*15+ 3]+=QC[0]*fsfs[0][6*10+3] + WQ[0]*fsfs[1][6*10+3]
			  + eta*tmp.fsds[6*6+3] + ze2*dsfs[1][1*10+3];
	    FSGS[6*15+ 4]+=QC[1]*fsfs[0][6*10+4] + WQ[1]*fsfs[1][6*10+4]
			  + eta*tmp.fsds[6*6+4] + ze22*dsfs[1][3*10+4];
	    FSGS[6*15+ 5]+=QC[2]*fsfs[0][6*10+5] + WQ[2]*fsfs[1][6*10+5]
			  + eta*tmp.fsds[6*6+5];
	    FSGS[6*15+ 6]+=QC[0]*fsfs[0][6*10+6] + WQ[0]*fsfs[1][6*10+6]
			  + eta2*tmp.fsds[6*6+1] + ze2*dsfs[1][1*10+6];
	    FSGS[6*15+ 7]+=QC[1]*fsfs[0][6*10+7] + WQ[1]*fsfs[1][6*10+7]
			  + eta2*tmp.fsds[6*6+2] + ze22*dsfs[1][3*10+7];
	    FSGS[6*15+ 8]+=QC[2]*fsfs[0][6*10+8] + WQ[2]*fsfs[1][6*10+8]
			  + eta2*tmp.fsds[6*6+0];
	    FSGS[6*15+ 9]+=QC[0]*fsfs[0][6*10+1] + WQ[0]*fsfs[1][6*10+1]
			  + ze2*dsfs[1][1*10+1];
	    FSGS[6*15+10]+=QC[1]*fsfs[0][6*10+2] + WQ[1]*fsfs[1][6*10+2]
			  + ze22*dsfs[1][3*10+2];
	    FSGS[6*15+11]+=QC[2]*fsfs[0][6*10+0] + WQ[2]*fsfs[1][6*10+0];
	    FSGS[6*15+12]+=QC[0]*fsfs[0][6*10+9] + WQ[0]*fsfs[1][6*10+9]
			  + eta2*tmp.fsds[6*6+4] + ze2*dsfs[1][1*10+9];
	    FSGS[6*15+13]+=QC[1]*fsfs[0][6*10+9] + WQ[1]*fsfs[1][6*10+9]
			  + eta2*tmp.fsds[6*6+5] + ze22*dsfs[1][3*10+9];
	    FSGS[6*15+14]+=QC[2]*fsfs[0][6*10+9] + WQ[2]*fsfs[1][6*10+9]
			  + eta2*tmp.fsds[6*6+3];
	    FSGS[7*15+ 0]+=QC[0]*fsfs[0][7*10+0] + WQ[0]*fsfs[1][7*10+0]
			  + eta23*tmp.fsds[7*6+0];
	    FSGS[7*15+ 1]+=QC[1]*fsfs[0][7*10+1] + WQ[1]*fsfs[1][7*10+1]
			  + eta23*tmp.fsds[7*6+1] + ze2*dsfs[1][2*10+1];
	    FSGS[7*15+ 2]+=QC[2]*fsfs[0][7*10+2] + WQ[2]*fsfs[1][7*10+2]
			  + eta23*tmp.fsds[7*6+2] + ze22*dsfs[1][4*10+2];
	    FSGS[7*15+ 3]+=QC[0]*fsfs[0][7*10+3] + WQ[0]*fsfs[1][7*10+3]
			  + eta*tmp.fsds[7*6+3];
	    FSGS[7*15+ 4]+=QC[1]*fsfs[0][7*10+4] + WQ[1]*fsfs[1][7*10+4]
			  + eta*tmp.fsds[7*6+4] + ze2*dsfs[1][2*10+4];
	    FSGS[7*15+ 5]+=QC[2]*fsfs[0][7*10+5] + WQ[2]*fsfs[1][7*10+5]
			  + eta*tmp.fsds[7*6+5] + ze22*dsfs[1][4*10+5];
	    FSGS[7*15+ 6]+=QC[0]*fsfs[0][7*10+6] + WQ[0]*fsfs[1][7*10+6]
			  + eta2*tmp.fsds[7*6+1];
	    FSGS[7*15+ 7]+=QC[1]*fsfs[0][7*10+7] + WQ[1]*fsfs[1][7*10+7]
			  + eta2*tmp.fsds[7*6+2] + ze2*dsfs[1][2*10+7];
	    FSGS[7*15+ 8]+=QC[2]*fsfs[0][7*10+8] + WQ[2]*fsfs[1][7*10+8]
			  + eta2*tmp.fsds[7*6+0] + ze22*dsfs[1][4*10+8];
	    FSGS[7*15+ 9]+=QC[0]*fsfs[0][7*10+1] + WQ[0]*fsfs[1][7*10+1];
	    FSGS[7*15+10]+=QC[1]*fsfs[0][7*10+2] + WQ[1]*fsfs[1][7*10+2]
			  + ze2*dsfs[1][2*10+2];
	    FSGS[7*15+11]+=QC[2]*fsfs[0][7*10+0] + WQ[2]*fsfs[1][7*10+0]
			  + ze22*dsfs[1][4*10+0];
	    FSGS[7*15+12]+=QC[0]*fsfs[0][7*10+9] + WQ[0]*fsfs[1][7*10+9]
			  + eta2*tmp.fsds[7*6+4];
	    FSGS[7*15+13]+=QC[1]*fsfs[0][7*10+9] + WQ[1]*fsfs[1][7*10+9]
			  + eta2*tmp.fsds[7*6+5] + ze2*dsfs[1][2*10+9];
	    FSGS[7*15+14]+=QC[2]*fsfs[0][7*10+9] + WQ[2]*fsfs[1][7*10+9]
			  + eta2*tmp.fsds[7*6+3] + ze22*dsfs[1][4*10+9];
	    FSGS[8*15+ 0]+=QC[0]*fsfs[0][8*10+0] + WQ[0]*fsfs[1][8*10+0]
			  + eta23*tmp.fsds[8*6+0] + ze22*dsfs[1][5*10+0];
	    FSGS[8*15+ 1]+=QC[1]*fsfs[0][8*10+1] + WQ[1]*fsfs[1][8*10+1]
			  + eta23*tmp.fsds[8*6+1];
	    FSGS[8*15+ 2]+=QC[2]*fsfs[0][8*10+2] + WQ[2]*fsfs[1][8*10+2]
			  + eta23*tmp.fsds[8*6+2] + ze2*dsfs[1][0*10+2];
	    FSGS[8*15+ 3]+=QC[0]*fsfs[0][8*10+3] + WQ[0]*fsfs[1][8*10+3]
			  + eta*tmp.fsds[8*6+3] + ze22*dsfs[1][5*10+3];
	    FSGS[8*15+ 4]+=QC[1]*fsfs[0][8*10+4] + WQ[1]*fsfs[1][8*10+4]
			  + eta*tmp.fsds[8*6+4];
	    FSGS[8*15+ 5]+=QC[2]*fsfs[0][8*10+5] + WQ[2]*fsfs[1][8*10+5]
			  + eta*tmp.fsds[8*6+5] + ze2*dsfs[1][0*10+5];
	    FSGS[8*15+ 6]+=QC[0]*fsfs[0][8*10+6] + WQ[0]*fsfs[1][8*10+6]
			  + eta2*tmp.fsds[8*6+1] + ze22*dsfs[1][5*10+6];
	    FSGS[8*15+ 7]+=QC[1]*fsfs[0][8*10+7] + WQ[1]*fsfs[1][8*10+7]
			  + eta2*tmp.fsds[8*6+2];
	    FSGS[8*15+ 8]+=QC[2]*fsfs[0][8*10+8] + WQ[2]*fsfs[1][8*10+8]
			  + eta2*tmp.fsds[8*6+0] + ze2*dsfs[1][0*10+8];
	    FSGS[8*15+ 9]+=QC[0]*fsfs[0][8*10+1] + WQ[0]*fsfs[1][8*10+1]
			  + ze22*dsfs[1][5*10+1];
	    FSGS[8*15+10]+=QC[1]*fsfs[0][8*10+2] + WQ[1]*fsfs[1][8*10+2];
	    FSGS[8*15+11]+=QC[2]*fsfs[0][8*10+0] + WQ[2]*fsfs[1][8*10+0]
			  + ze2*dsfs[1][0*10+0];
	    FSGS[8*15+12]+=QC[0]*fsfs[0][8*10+9] + WQ[0]*fsfs[1][8*10+9]
			  + eta2*tmp.fsds[8*6+4] + ze22*dsfs[1][5*10+9];
	    FSGS[8*15+13]+=QC[1]*fsfs[0][8*10+9] + WQ[1]*fsfs[1][8*10+9]
			  + eta2*tmp.fsds[8*6+5];
	    FSGS[8*15+14]+=QC[2]*fsfs[0][8*10+9] + WQ[2]*fsfs[1][8*10+9]
			  + eta2*tmp.fsds[8*6+3] + ze2*dsfs[1][0*10+9];
	    FSGS[9*15+ 0]+=QC[0]*fsfs[0][9*10+0] + WQ[0]*fsfs[1][9*10+0]
			  + eta23*tmp.fsds[9*6+0] + ze2*dsfs[1][4*10+0];
	    FSGS[9*15+ 1]+=QC[1]*fsfs[0][9*10+1] + WQ[1]*fsfs[1][9*10+1]
			  + eta23*tmp.fsds[9*6+1] + ze2*dsfs[1][5*10+1];
	    FSGS[9*15+ 2]+=QC[2]*fsfs[0][9*10+2] + WQ[2]*fsfs[1][9*10+2]
			  + eta23*tmp.fsds[9*6+2] + ze2*dsfs[1][3*10+2];
	    FSGS[9*15+ 3]+=QC[0]*fsfs[0][9*10+3] + WQ[0]*fsfs[1][9*10+3]
			  + eta*tmp.fsds[9*6+3] + ze2*dsfs[1][4*10+3];
	    FSGS[9*15+ 4]+=QC[1]*fsfs[0][9*10+4] + WQ[1]*fsfs[1][9*10+4]
			  + eta*tmp.fsds[9*6+4] + ze2*dsfs[1][5*10+4];
	    FSGS[9*15+ 5]+=QC[2]*fsfs[0][9*10+5] + WQ[2]*fsfs[1][9*10+5]
			  + eta*tmp.fsds[9*6+5] + ze2*dsfs[1][3*10+5];
	    FSGS[9*15+ 6]+=QC[0]*fsfs[0][9*10+6] + WQ[0]*fsfs[1][9*10+6]
			  + eta2*tmp.fsds[9*6+1] + ze2*dsfs[1][4*10+6];
	    FSGS[9*15+ 7]+=QC[1]*fsfs[0][9*10+7] + WQ[1]*fsfs[1][9*10+7]
			  + eta2*tmp.fsds[9*6+2] + ze2*dsfs[1][5*10+7];
	    FSGS[9*15+ 8]+=QC[2]*fsfs[0][9*10+8] + WQ[2]*fsfs[1][9*10+8]
			  + eta2*tmp.fsds[9*6+0] + ze2*dsfs[1][3*10+8];
	    FSGS[9*15+ 9]+=QC[0]*fsfs[0][9*10+1] + WQ[0]*fsfs[1][9*10+1]
			  + ze2*dsfs[1][4*10+1];
	    FSGS[9*15+10]+=QC[1]*fsfs[0][9*10+2] + WQ[1]*fsfs[1][9*10+2]
			  + ze2*dsfs[1][5*10+2];
	    FSGS[9*15+11]+=QC[2]*fsfs[0][9*10+0] + WQ[2]*fsfs[1][9*10+0]
			  + ze2*dsfs[1][3*10+0];
	    FSGS[9*15+12]+=QC[0]*fsfs[0][9*10+9] + WQ[0]*fsfs[1][9*10+9]
			  + eta2*tmp.fsds[9*6+4] + ze2*dsfs[1][4*10+9];
	    FSGS[9*15+13]+=QC[1]*fsfs[0][9*10+9] + WQ[1]*fsfs[1][9*10+9]
			  + eta2*tmp.fsds[9*6+5] + ze2*dsfs[1][5*10+9];
	    FSGS[9*15+14]+=QC[2]*fsfs[0][9*10+9] + WQ[2]*fsfs[1][9*10+9]
			  + eta2*tmp.fsds[9*6+3] + ze2*dsfs[1][3*10+9];
	    // (G,S|G,S)
	    for (i=0; i<15*6; i++) tmp.gsds[i]=gsds[0][i]-re*gsds[1][i];
	    GSGS[ 0*15+ 0]+=QC[0]*gsfs[0][ 0*10+0]+WQ[0]*gsfs[1][ 0*10+0]
			   + eta23*tmp.gsds[ 0*6+0] + ze24*fsfs[1][0*10+0];
	    GSGS[ 0*15+ 1]+=QC[1]*gsfs[0][ 0*10+1]+WQ[1]*gsfs[1][ 0*10+1]
			   + eta23*tmp.gsds[ 0*6+1];
	    GSGS[ 0*15+ 2]+=QC[2]*gsfs[0][ 0*10+2]+WQ[2]*gsfs[1][ 0*10+2]
			   + eta23*tmp.gsds[ 0*6+2];
	    GSGS[ 0*15+ 3]+=QC[0]*gsfs[0][ 0*10+3]+WQ[0]*gsfs[1][ 0*10+3]
			   + eta*tmp.gsds[ 0*6+3] + ze24*fsfs[1][0*10+3];
	    GSGS[ 0*15+ 4]+=QC[1]*gsfs[0][ 0*10+4]+WQ[1]*gsfs[1][ 0*10+4]
			   + eta*tmp.gsds[ 0*6+4];
	    GSGS[ 0*15+ 5]+=QC[2]*gsfs[0][ 0*10+5]+WQ[2]*gsfs[1][ 0*10+5]
			   + eta*tmp.gsds[ 0*6+5];
	    GSGS[ 0*15+ 6]+=QC[0]*gsfs[0][ 0*10+6]+WQ[0]*gsfs[1][ 0*10+6]
			   + eta2*tmp.gsds[ 0*6+1] + ze24*fsfs[1][0*10+6];
	    GSGS[ 0*15+ 7]+=QC[1]*gsfs[0][ 0*10+7]+WQ[1]*gsfs[1][ 0*10+7]
			   + eta2*tmp.gsds[ 0*6+2];
	    GSGS[ 0*15+ 8]+=QC[2]*gsfs[0][ 0*10+8]+WQ[2]*gsfs[1][ 0*10+8]
			   + eta2*tmp.gsds[ 0*6+0];
	    GSGS[ 0*15+ 9]+=QC[0]*gsfs[0][ 0*10+1]+WQ[0]*gsfs[1][ 0*10+1]
			   + ze24*fsfs[1][0*10+1];
	    GSGS[ 0*15+10]+=QC[1]*gsfs[0][ 0*10+2]+WQ[1]*gsfs[1][ 0*10+2];
	    GSGS[ 0*15+11]+=QC[2]*gsfs[0][ 0*10+0]+WQ[2]*gsfs[1][ 0*10+0];
	    GSGS[ 0*15+12]+=QC[0]*gsfs[0][ 0*10+9]+WQ[0]*gsfs[1][ 0*10+9]
			   + eta2*tmp.gsds[ 0*6+4] + ze24*fsfs[1][0*10+9];
	    GSGS[ 0*15+13]+=QC[1]*gsfs[0][ 0*10+9]+WQ[1]*gsfs[1][ 0*10+9]
			   + eta2*tmp.gsds[ 0*6+5];
	    GSGS[ 0*15+14]+=QC[2]*gsfs[0][ 0*10+9]+WQ[2]*gsfs[1][ 0*10+9]
			   + eta2*tmp.gsds[ 0*6+3];
	    GSGS[ 1*15+ 0]+=QC[0]*gsfs[0][ 1*10+0]+WQ[0]*gsfs[1][ 1*10+0]
			   + eta23*tmp.gsds[ 1*6+0];
	    GSGS[ 1*15+ 1]+=QC[1]*gsfs[0][ 1*10+1]+WQ[1]*gsfs[1][ 1*10+1]
			   + eta23*tmp.gsds[ 1*6+1] + ze24*fsfs[1][1*10+1];
	    GSGS[ 1*15+ 2]+=QC[2]*gsfs[0][ 1*10+2]+WQ[2]*gsfs[1][ 1*10+2]
			   + eta23*tmp.gsds[ 1*6+2];
	    GSGS[ 1*15+ 3]+=QC[0]*gsfs[0][ 1*10+3]+WQ[0]*gsfs[1][ 1*10+3]
			   + eta*tmp.gsds[ 1*6+3];
	    GSGS[ 1*15+ 4]+=QC[1]*gsfs[0][ 1*10+4]+WQ[1]*gsfs[1][ 1*10+4]
			   + eta*tmp.gsds[ 1*6+4] + ze24*fsfs[1][1*10+4];
	    GSGS[ 1*15+ 5]+=QC[2]*gsfs[0][ 1*10+5]+WQ[2]*gsfs[1][ 1*10+5]
			   + eta*tmp.gsds[ 1*6+5];
	    GSGS[ 1*15+ 6]+=QC[0]*gsfs[0][ 1*10+6]+WQ[0]*gsfs[1][ 1*10+6]
			   + eta2*tmp.gsds[ 1*6+1];
	    GSGS[ 1*15+ 7]+=QC[1]*gsfs[0][ 1*10+7]+WQ[1]*gsfs[1][ 1*10+7]
			   + eta2*tmp.gsds[ 1*6+2] + ze24*fsfs[1][1*10+7];
	    GSGS[ 1*15+ 8]+=QC[2]*gsfs[0][ 1*10+8]+WQ[2]*gsfs[1][ 1*10+8]
			   + eta2*tmp.gsds[ 1*6+0];
	    GSGS[ 1*15+ 9]+=QC[0]*gsfs[0][ 1*10+1]+WQ[0]*gsfs[1][ 1*10+1];
	    GSGS[ 1*15+10]+=QC[1]*gsfs[0][ 1*10+2]+WQ[1]*gsfs[1][ 1*10+2]
			   + ze24*fsfs[1][1*10+2];
	    GSGS[ 1*15+11]+=QC[2]*gsfs[0][ 1*10+0]+WQ[2]*gsfs[1][ 1*10+0];
	    GSGS[ 1*15+12]+=QC[0]*gsfs[0][ 1*10+9]+WQ[0]*gsfs[1][ 1*10+9]
			   + eta2*tmp.gsds[ 1*6+4];
	    GSGS[ 1*15+13]+=QC[1]*gsfs[0][ 1*10+9]+WQ[1]*gsfs[1][ 1*10+9]
			   + eta2*tmp.gsds[ 1*6+5] + ze24*fsfs[1][1*10+9];
	    GSGS[ 1*15+14]+=QC[2]*gsfs[0][ 1*10+9]+WQ[2]*gsfs[1][ 1*10+9]
			   + eta2*tmp.gsds[ 1*6+3];
	    GSGS[ 2*15+ 0]+=QC[0]*gsfs[0][ 2*10+0]+WQ[0]*gsfs[1][ 2*10+0]
			   + eta23*tmp.gsds[ 2*6+0];
	    GSGS[ 2*15+ 1]+=QC[1]*gsfs[0][ 2*10+1]+WQ[1]*gsfs[1][ 2*10+1]
			   + eta23*tmp.gsds[ 2*6+1];
	    GSGS[ 2*15+ 2]+=QC[2]*gsfs[0][ 2*10+2]+WQ[2]*gsfs[1][ 2*10+2]
			   + eta23*tmp.gsds[ 2*6+2] + ze24*fsfs[1][2*10+2];
	    GSGS[ 2*15+ 3]+=QC[0]*gsfs[0][ 2*10+3]+WQ[0]*gsfs[1][ 2*10+3]
			   + eta*tmp.gsds[ 2*6+3];
	    GSGS[ 2*15+ 4]+=QC[1]*gsfs[0][ 2*10+4]+WQ[1]*gsfs[1][ 2*10+4]
			   + eta*tmp.gsds[ 2*6+4];
	    GSGS[ 2*15+ 5]+=QC[2]*gsfs[0][ 2*10+5]+WQ[2]*gsfs[1][ 2*10+5]
			   + eta*tmp.gsds[ 2*6+5] + ze24*fsfs[1][2*10+5];
	    GSGS[ 2*15+ 6]+=QC[0]*gsfs[0][ 2*10+6]+WQ[0]*gsfs[1][ 2*10+6]
			   + eta2*tmp.gsds[ 2*6+1];
	    GSGS[ 2*15+ 7]+=QC[1]*gsfs[0][ 2*10+7]+WQ[1]*gsfs[1][ 2*10+7]
			   + eta2*tmp.gsds[ 2*6+2];
	    GSGS[ 2*15+ 8]+=QC[2]*gsfs[0][ 2*10+8]+WQ[2]*gsfs[1][ 2*10+8]
			   + eta2*tmp.gsds[ 2*6+0] + ze24*fsfs[1][2*10+8];
	    GSGS[ 2*15+ 9]+=QC[0]*gsfs[0][ 2*10+1]+WQ[0]*gsfs[1][ 2*10+1];
	    GSGS[ 2*15+10]+=QC[1]*gsfs[0][ 2*10+2]+WQ[1]*gsfs[1][ 2*10+2];
	    GSGS[ 2*15+11]+=QC[2]*gsfs[0][ 2*10+0]+WQ[2]*gsfs[1][ 2*10+0]
			   + ze24*fsfs[1][2*10+0];
	    GSGS[ 2*15+12]+=QC[0]*gsfs[0][ 2*10+9]+WQ[0]*gsfs[1][ 2*10+9]
			   + eta2*tmp.gsds[ 2*6+4];
	    GSGS[ 2*15+13]+=QC[1]*gsfs[0][ 2*10+9]+WQ[1]*gsfs[1][ 2*10+9]
			   + eta2*tmp.gsds[ 2*6+5];
	    GSGS[ 2*15+14]+=QC[2]*gsfs[0][ 2*10+9]+WQ[2]*gsfs[1][ 2*10+9]
			   + eta2*tmp.gsds[ 2*6+3] + ze24*fsfs[1][2*10+9];
	    GSGS[ 3*15+ 0]+=QC[0]*gsfs[0][ 3*10+0]+WQ[0]*gsfs[1][ 3*10+0]
			   + eta23*tmp.gsds[ 3*6+0] + ze23*fsfs[1][3*10+0];
	    GSGS[ 3*15+ 1]+=QC[1]*gsfs[0][ 3*10+1]+WQ[1]*gsfs[1][ 3*10+1]
			   + eta23*tmp.gsds[ 3*6+1] + ze2*fsfs[1][0*10+1];
	    GSGS[ 3*15+ 2]+=QC[2]*gsfs[0][ 3*10+2]+WQ[2]*gsfs[1][ 3*10+2]
			   + eta23*tmp.gsds[ 3*6+2];
	    GSGS[ 3*15+ 3]+=QC[0]*gsfs[0][ 3*10+3]+WQ[0]*gsfs[1][ 3*10+3]
			   + eta*tmp.gsds[ 3*6+3] + ze23*fsfs[1][3*10+3];
	    GSGS[ 3*15+ 4]+=QC[1]*gsfs[0][ 3*10+4]+WQ[1]*gsfs[1][ 3*10+4]
			   + eta*tmp.gsds[ 3*6+4] + ze2*fsfs[1][0*10+4];
	    GSGS[ 3*15+ 5]+=QC[2]*gsfs[0][ 3*10+5]+WQ[2]*gsfs[1][ 3*10+5]
			   + eta*tmp.gsds[ 3*6+5];
	    GSGS[ 3*15+ 6]+=QC[0]*gsfs[0][ 3*10+6]+WQ[0]*gsfs[1][ 3*10+6]
			   + eta2*tmp.gsds[ 3*6+1] + ze23*fsfs[1][3*10+6];
	    GSGS[ 3*15+ 7]+=QC[1]*gsfs[0][ 3*10+7]+WQ[1]*gsfs[1][ 3*10+7]
			   + eta2*tmp.gsds[ 3*6+2] + ze2*fsfs[1][0*10+7];
	    GSGS[ 3*15+ 8]+=QC[2]*gsfs[0][ 3*10+8]+WQ[2]*gsfs[1][ 3*10+8]
			   + eta2*tmp.gsds[ 3*6+0];
	    GSGS[ 3*15+ 9]+=QC[0]*gsfs[0][ 3*10+1]+WQ[0]*gsfs[1][ 3*10+1]
			   + ze23*fsfs[1][3*10+1];
	    GSGS[ 3*15+10]+=QC[1]*gsfs[0][ 3*10+2]+WQ[1]*gsfs[1][ 3*10+2]
			   + ze2*fsfs[1][0*10+2];
	    GSGS[ 3*15+11]+=QC[2]*gsfs[0][ 3*10+0]+WQ[2]*gsfs[1][ 3*10+0];
	    GSGS[ 3*15+12]+=QC[0]*gsfs[0][ 3*10+9]+WQ[0]*gsfs[1][ 3*10+9]
			   + eta2*tmp.gsds[ 3*6+4] + ze23*fsfs[1][3*10+9];
	    GSGS[ 3*15+13]+=QC[1]*gsfs[0][ 3*10+9]+WQ[1]*gsfs[1][ 3*10+9]
			   + eta2*tmp.gsds[ 3*6+5] + ze2*fsfs[1][0*10+9];
	    GSGS[ 3*15+14]+=QC[2]*gsfs[0][ 3*10+9]+WQ[2]*gsfs[1][ 3*10+9]
			   + eta2*tmp.gsds[ 3*6+3];
	    GSGS[ 4*15+ 0]+=QC[0]*gsfs[0][ 4*10+0]+WQ[0]*gsfs[1][ 4*10+0]
			   + eta23*tmp.gsds[ 4*6+0];
	    GSGS[ 4*15+ 1]+=QC[1]*gsfs[0][ 4*10+1]+WQ[1]*gsfs[1][ 4*10+1]
			   + eta23*tmp.gsds[ 4*6+1] + ze23*fsfs[1][4*10+1];
	    GSGS[ 4*15+ 2]+=QC[2]*gsfs[0][ 4*10+2]+WQ[2]*gsfs[1][ 4*10+2]
			   + eta23*tmp.gsds[ 4*6+2] + ze2*fsfs[1][1*10+2];
	    GSGS[ 4*15+ 3]+=QC[0]*gsfs[0][ 4*10+3]+WQ[0]*gsfs[1][ 4*10+3]
			   + eta*tmp.gsds[ 4*6+3];
	    GSGS[ 4*15+ 4]+=QC[1]*gsfs[0][ 4*10+4]+WQ[1]*gsfs[1][ 4*10+4]
			   + eta*tmp.gsds[ 4*6+4] + ze23*fsfs[1][4*10+4];
	    GSGS[ 4*15+ 5]+=QC[2]*gsfs[0][ 4*10+5]+WQ[2]*gsfs[1][ 4*10+5]
			   + eta*tmp.gsds[ 4*6+5] + ze2*fsfs[1][1*10+5];
	    GSGS[ 4*15+ 6]+=QC[0]*gsfs[0][ 4*10+6]+WQ[0]*gsfs[1][ 4*10+6]
			   + eta2*tmp.gsds[ 4*6+1];
	    GSGS[ 4*15+ 7]+=QC[1]*gsfs[0][ 4*10+7]+WQ[1]*gsfs[1][ 4*10+7]
			   + eta2*tmp.gsds[ 4*6+2] + ze23*fsfs[1][4*10+7];
	    GSGS[ 4*15+ 8]+=QC[2]*gsfs[0][ 4*10+8]+WQ[2]*gsfs[1][ 4*10+8]
			   + eta2*tmp.gsds[ 4*6+0] + ze2*fsfs[1][1*10+8];
	    GSGS[ 4*15+ 9]+=QC[0]*gsfs[0][ 4*10+1]+WQ[0]*gsfs[1][ 4*10+1];
	    GSGS[ 4*15+10]+=QC[1]*gsfs[0][ 4*10+2]+WQ[1]*gsfs[1][ 4*10+2]
			   + ze23*fsfs[1][4*10+2];
	    GSGS[ 4*15+11]+=QC[2]*gsfs[0][ 4*10+0]+WQ[2]*gsfs[1][ 4*10+0]
			   + ze2*fsfs[1][1*10+0];
	    GSGS[ 4*15+12]+=QC[0]*gsfs[0][ 4*10+9]+WQ[0]*gsfs[1][ 4*10+9]
			   + eta2*tmp.gsds[ 4*6+4];
	    GSGS[ 4*15+13]+=QC[1]*gsfs[0][ 4*10+9]+WQ[1]*gsfs[1][ 4*10+9]
			   + eta2*tmp.gsds[ 4*6+5] + ze23*fsfs[1][4*10+9];
	    GSGS[ 4*15+14]+=QC[2]*gsfs[0][ 4*10+9]+WQ[2]*gsfs[1][ 4*10+9]
			   + eta2*tmp.gsds[ 4*6+3] + ze2*fsfs[1][1*10+9];
	    GSGS[ 5*15+ 0]+=QC[0]*gsfs[0][ 5*10+0]+WQ[0]*gsfs[1][ 5*10+0]
			   + eta23*tmp.gsds[ 5*6+0] + ze2*fsfs[1][2*10+0];
	    GSGS[ 5*15+ 1]+=QC[1]*gsfs[0][ 5*10+1]+WQ[1]*gsfs[1][ 5*10+1]
			   + eta23*tmp.gsds[ 5*6+1];
	    GSGS[ 5*15+ 2]+=QC[2]*gsfs[0][ 5*10+2]+WQ[2]*gsfs[1][ 5*10+2]
			   + eta23*tmp.gsds[ 5*6+2] + ze23*fsfs[1][5*10+2];
	    GSGS[ 5*15+ 3]+=QC[0]*gsfs[0][ 5*10+3]+WQ[0]*gsfs[1][ 5*10+3]
			   + eta*tmp.gsds[ 5*6+3] + ze2*fsfs[1][2*10+3];
	    GSGS[ 5*15+ 4]+=QC[1]*gsfs[0][ 5*10+4]+WQ[1]*gsfs[1][ 5*10+4]
			   + eta*tmp.gsds[ 5*6+4];
	    GSGS[ 5*15+ 5]+=QC[2]*gsfs[0][ 5*10+5]+WQ[2]*gsfs[1][ 5*10+5]
			   + eta*tmp.gsds[ 5*6+5] + ze23*fsfs[1][5*10+5];
	    GSGS[ 5*15+ 6]+=QC[0]*gsfs[0][ 5*10+6]+WQ[0]*gsfs[1][ 5*10+6]
			   + eta2*tmp.gsds[ 5*6+1] + ze2*fsfs[1][2*10+6];
	    GSGS[ 5*15+ 7]+=QC[1]*gsfs[0][ 5*10+7]+WQ[1]*gsfs[1][ 5*10+7]
			   + eta2*tmp.gsds[ 5*6+2];
	    GSGS[ 5*15+ 8]+=QC[2]*gsfs[0][ 5*10+8]+WQ[2]*gsfs[1][ 5*10+8]
			   + eta2*tmp.gsds[ 5*6+0] + ze23*fsfs[1][5*10+8];
	    GSGS[ 5*15+ 9]+=QC[0]*gsfs[0][ 5*10+1]+WQ[0]*gsfs[1][ 5*10+1]
			   + ze2*fsfs[1][2*10+1];
	    GSGS[ 5*15+10]+=QC[1]*gsfs[0][ 5*10+2]+WQ[1]*gsfs[1][ 5*10+2];
	    GSGS[ 5*15+11]+=QC[2]*gsfs[0][ 5*10+0]+WQ[2]*gsfs[1][ 5*10+0]
			   + ze23*fsfs[1][5*10+0];
	    GSGS[ 5*15+12]+=QC[0]*gsfs[0][ 5*10+9]+WQ[0]*gsfs[1][ 5*10+9]
			   + eta2*tmp.gsds[ 5*6+4] + ze2*fsfs[1][2*10+9];
	    GSGS[ 5*15+13]+=QC[1]*gsfs[0][ 5*10+9]+WQ[1]*gsfs[1][ 5*10+9]
			   + eta2*tmp.gsds[ 5*6+5];
	    GSGS[ 5*15+14]+=QC[2]*gsfs[0][ 5*10+9]+WQ[2]*gsfs[1][ 5*10+9]
			   + eta2*tmp.gsds[ 5*6+3] + ze23*fsfs[1][5*10+9];
	    GSGS[ 6*15+ 0]+=QC[0]*gsfs[0][ 6*10+0]+WQ[0]*gsfs[1][ 6*10+0]
			   + eta23*tmp.gsds[ 6*6+0] + ze22*fsfs[1][6*10+0];
	    GSGS[ 6*15+ 1]+=QC[1]*gsfs[0][ 6*10+1]+WQ[1]*gsfs[1][ 6*10+1]
			   + eta23*tmp.gsds[ 6*6+1] + ze22*fsfs[1][3*10+1];
	    GSGS[ 6*15+ 2]+=QC[2]*gsfs[0][ 6*10+2]+WQ[2]*gsfs[1][ 6*10+2]
			   + eta23*tmp.gsds[ 6*6+2];
	    GSGS[ 6*15+ 3]+=QC[0]*gsfs[0][ 6*10+3]+WQ[0]*gsfs[1][ 6*10+3]
			   + eta*tmp.gsds[ 6*6+3] + ze22*fsfs[1][6*10+3];
	    GSGS[ 6*15+ 4]+=QC[1]*gsfs[0][ 6*10+4]+WQ[1]*gsfs[1][ 6*10+4]
			   + eta*tmp.gsds[ 6*6+4] + ze22*fsfs[1][3*10+4];
	    GSGS[ 6*15+ 5]+=QC[2]*gsfs[0][ 6*10+5]+WQ[2]*gsfs[1][ 6*10+5]
			   + eta*tmp.gsds[ 6*6+5];
	    GSGS[ 6*15+ 6]+=QC[0]*gsfs[0][ 6*10+6]+WQ[0]*gsfs[1][ 6*10+6]
			   + eta2*tmp.gsds[ 6*6+1] + ze22*fsfs[1][6*10+6];
	    GSGS[ 6*15+ 7]+=QC[1]*gsfs[0][ 6*10+7]+WQ[1]*gsfs[1][ 6*10+7]
			   + eta2*tmp.gsds[ 6*6+2] + ze22*fsfs[1][3*10+7];
	    GSGS[ 6*15+ 8]+=QC[2]*gsfs[0][ 6*10+8]+WQ[2]*gsfs[1][ 6*10+8]
			   + eta2*tmp.gsds[ 6*6+0];
	    GSGS[ 6*15+ 9]+=QC[0]*gsfs[0][ 6*10+1]+WQ[0]*gsfs[1][ 6*10+1]
			   + ze22*fsfs[1][6*10+1];
	    GSGS[ 6*15+10]+=QC[1]*gsfs[0][ 6*10+2]+WQ[1]*gsfs[1][ 6*10+2]
			   + ze22*fsfs[1][3*10+2];
	    GSGS[ 6*15+11]+=QC[2]*gsfs[0][ 6*10+0]+WQ[2]*gsfs[1][ 6*10+0];
	    GSGS[ 6*15+12]+=QC[0]*gsfs[0][ 6*10+9]+WQ[0]*gsfs[1][ 6*10+9]
			   + eta2*tmp.gsds[ 6*6+4] + ze22*fsfs[1][6*10+9];
	    GSGS[ 6*15+13]+=QC[1]*gsfs[0][ 6*10+9]+WQ[1]*gsfs[1][ 6*10+9]
			   + eta2*tmp.gsds[ 6*6+5] + ze22*fsfs[1][3*10+9];
	    GSGS[ 6*15+14]+=QC[2]*gsfs[0][ 6*10+9]+WQ[2]*gsfs[1][ 6*10+9]
			   + eta2*tmp.gsds[ 6*6+3];
	    GSGS[ 7*15+ 0]+=QC[0]*gsfs[0][ 7*10+0]+WQ[0]*gsfs[1][ 7*10+0]
			   + eta23*tmp.gsds[ 7*6+0];
	    GSGS[ 7*15+ 1]+=QC[1]*gsfs[0][ 7*10+1]+WQ[1]*gsfs[1][ 7*10+1]
			   + eta23*tmp.gsds[ 7*6+1] + ze22*fsfs[1][7*10+1];
	    GSGS[ 7*15+ 2]+=QC[2]*gsfs[0][ 7*10+2]+WQ[2]*gsfs[1][ 7*10+2]
			   + eta23*tmp.gsds[ 7*6+2] + ze22*fsfs[1][4*10+2];
	    GSGS[ 7*15+ 3]+=QC[0]*gsfs[0][ 7*10+3]+WQ[0]*gsfs[1][ 7*10+3]
			   + eta*tmp.gsds[ 7*6+3];
	    GSGS[ 7*15+ 4]+=QC[1]*gsfs[0][ 7*10+4]+WQ[1]*gsfs[1][ 7*10+4]
			   + eta*tmp.gsds[ 7*6+4] + ze22*fsfs[1][7*10+4];
	    GSGS[ 7*15+ 5]+=QC[2]*gsfs[0][ 7*10+5]+WQ[2]*gsfs[1][ 7*10+5]
			   + eta*tmp.gsds[ 7*6+5] + ze22*fsfs[1][4*10+5];
	    GSGS[ 7*15+ 6]+=QC[0]*gsfs[0][ 7*10+6]+WQ[0]*gsfs[1][ 7*10+6]
			   + eta2*tmp.gsds[ 7*6+1];
	    GSGS[ 7*15+ 7]+=QC[1]*gsfs[0][ 7*10+7]+WQ[1]*gsfs[1][ 7*10+7]
			   + eta2*tmp.gsds[ 7*6+2] + ze22*fsfs[1][7*10+7];
	    GSGS[ 7*15+ 8]+=QC[2]*gsfs[0][ 7*10+8]+WQ[2]*gsfs[1][ 7*10+8]
			   + eta2*tmp.gsds[ 7*6+0] + ze22*fsfs[1][4*10+8];
	    GSGS[ 7*15+ 9]+=QC[0]*gsfs[0][ 7*10+1]+WQ[0]*gsfs[1][ 7*10+1];
	    GSGS[ 7*15+10]+=QC[1]*gsfs[0][ 7*10+2]+WQ[1]*gsfs[1][ 7*10+2]
			   + ze22*fsfs[1][7*10+2];
	    GSGS[ 7*15+11]+=QC[2]*gsfs[0][ 7*10+0]+WQ[2]*gsfs[1][ 7*10+0]
			   + ze22*fsfs[1][4*10+0];
	    GSGS[ 7*15+12]+=QC[0]*gsfs[0][ 7*10+9]+WQ[0]*gsfs[1][ 7*10+9]
			   + eta2*tmp.gsds[ 7*6+4];
	    GSGS[ 7*15+13]+=QC[1]*gsfs[0][ 7*10+9]+WQ[1]*gsfs[1][ 7*10+9]
			   + eta2*tmp.gsds[ 7*6+5] + ze22*fsfs[1][7*10+9];
	    GSGS[ 7*15+14]+=QC[2]*gsfs[0][ 7*10+9]+WQ[2]*gsfs[1][ 7*10+9]
			   + eta2*tmp.gsds[ 7*6+3] + ze22*fsfs[1][4*10+9];
	    GSGS[ 8*15+ 0]+=QC[0]*gsfs[0][ 8*10+0]+WQ[0]*gsfs[1][ 8*10+0]
			   + eta23*tmp.gsds[ 8*6+0] + ze22*fsfs[1][5*10+0];
	    GSGS[ 8*15+ 1]+=QC[1]*gsfs[0][ 8*10+1]+WQ[1]*gsfs[1][ 8*10+1]
			   + eta23*tmp.gsds[ 8*6+1];
	    GSGS[ 8*15+ 2]+=QC[2]*gsfs[0][ 8*10+2]+WQ[2]*gsfs[1][ 8*10+2]
			   + eta23*tmp.gsds[ 8*6+2] + ze22*fsfs[1][8*10+2];
	    GSGS[ 8*15+ 3]+=QC[0]*gsfs[0][ 8*10+3]+WQ[0]*gsfs[1][ 8*10+3]
			   + eta*tmp.gsds[ 8*6+3] + ze22*fsfs[1][5*10+3];
	    GSGS[ 8*15+ 4]+=QC[1]*gsfs[0][ 8*10+4]+WQ[1]*gsfs[1][ 8*10+4]
			   + eta*tmp.gsds[ 8*6+4];
	    GSGS[ 8*15+ 5]+=QC[2]*gsfs[0][ 8*10+5]+WQ[2]*gsfs[1][ 8*10+5]
			   + eta*tmp.gsds[ 8*6+5] + ze22*fsfs[1][8*10+5];
	    GSGS[ 8*15+ 6]+=QC[0]*gsfs[0][ 8*10+6]+WQ[0]*gsfs[1][ 8*10+6]
			   + eta2*tmp.gsds[ 8*6+1] + ze22*fsfs[1][5*10+6];
	    GSGS[ 8*15+ 7]+=QC[1]*gsfs[0][ 8*10+7]+WQ[1]*gsfs[1][ 8*10+7]
			   + eta2*tmp.gsds[ 8*6+2];
	    GSGS[ 8*15+ 8]+=QC[2]*gsfs[0][ 8*10+8]+WQ[2]*gsfs[1][ 8*10+8]
			   + eta2*tmp.gsds[ 8*6+0] + ze22*fsfs[1][8*10+8];
	    GSGS[ 8*15+ 9]+=QC[0]*gsfs[0][ 8*10+1]+WQ[0]*gsfs[1][ 8*10+1]
			   + ze22*fsfs[1][5*10+1];
	    GSGS[ 8*15+10]+=QC[1]*gsfs[0][ 8*10+2]+WQ[1]*gsfs[1][ 8*10+2];
	    GSGS[ 8*15+11]+=QC[2]*gsfs[0][ 8*10+0]+WQ[2]*gsfs[1][ 8*10+0]
			   + ze22*fsfs[1][8*10+0];
	    GSGS[ 8*15+12]+=QC[0]*gsfs[0][ 8*10+9]+WQ[0]*gsfs[1][ 8*10+9]
			   + eta2*tmp.gsds[ 8*6+4] + ze22*fsfs[1][5*10+9];
	    GSGS[ 8*15+13]+=QC[1]*gsfs[0][ 8*10+9]+WQ[1]*gsfs[1][ 8*10+9]
			   + eta2*tmp.gsds[ 8*6+5];
	    GSGS[ 8*15+14]+=QC[2]*gsfs[0][ 8*10+9]+WQ[2]*gsfs[1][ 8*10+9]
			   + eta2*tmp.gsds[ 8*6+3] + ze22*fsfs[1][8*10+9];
	    GSGS[ 9*15+ 0]+=QC[0]*gsfs[0][ 9*10+0]+WQ[0]*gsfs[1][ 9*10+0]
			   + eta23*tmp.gsds[ 9*6+0] + ze2*fsfs[1][1*10+0];
	    GSGS[ 9*15+ 1]+=QC[1]*gsfs[0][ 9*10+1]+WQ[1]*gsfs[1][ 9*10+1]
			   + eta23*tmp.gsds[ 9*6+1] + ze23*fsfs[1][6*10+1];
	    GSGS[ 9*15+ 2]+=QC[2]*gsfs[0][ 9*10+2]+WQ[2]*gsfs[1][ 9*10+2]
			   + eta23*tmp.gsds[ 9*6+2];
	    GSGS[ 9*15+ 3]+=QC[0]*gsfs[0][ 9*10+3]+WQ[0]*gsfs[1][ 9*10+3]
			   + eta*tmp.gsds[ 9*6+3] + ze2*fsfs[1][1*10+3];
	    GSGS[ 9*15+ 4]+=QC[1]*gsfs[0][ 9*10+4]+WQ[1]*gsfs[1][ 9*10+4]
			   + eta*tmp.gsds[ 9*6+4] + ze23*fsfs[1][6*10+4];
	    GSGS[ 9*15+ 5]+=QC[2]*gsfs[0][ 9*10+5]+WQ[2]*gsfs[1][ 9*10+5]
			   + eta*tmp.gsds[ 9*6+5];
	    GSGS[ 9*15+ 6]+=QC[0]*gsfs[0][ 9*10+6]+WQ[0]*gsfs[1][ 9*10+6]
			   + eta2*tmp.gsds[ 9*6+1] + ze2*fsfs[1][1*10+6];
	    GSGS[ 9*15+ 7]+=QC[1]*gsfs[0][ 9*10+7]+WQ[1]*gsfs[1][ 9*10+7]
			   + eta2*tmp.gsds[ 9*6+2] + ze23*fsfs[1][6*10+7];
	    GSGS[ 9*15+ 8]+=QC[2]*gsfs[0][ 9*10+8]+WQ[2]*gsfs[1][ 9*10+8]
			   + eta2*tmp.gsds[ 9*6+0];
	    GSGS[ 9*15+ 9]+=QC[0]*gsfs[0][ 9*10+1]+WQ[0]*gsfs[1][ 9*10+1]
			   + ze2*fsfs[1][1*10+1];
	    GSGS[ 9*15+10]+=QC[1]*gsfs[0][ 9*10+2]+WQ[1]*gsfs[1][ 9*10+2]
			   + ze23*fsfs[1][6*10+2];
	    GSGS[ 9*15+11]+=QC[2]*gsfs[0][ 9*10+0]+WQ[2]*gsfs[1][ 9*10+0];
	    GSGS[ 9*15+12]+=QC[0]*gsfs[0][ 9*10+9]+WQ[0]*gsfs[1][ 9*10+9]
			   + eta2*tmp.gsds[ 9*6+4] + ze2*fsfs[1][1*10+9];
	    GSGS[ 9*15+13]+=QC[1]*gsfs[0][ 9*10+9]+WQ[1]*gsfs[1][ 9*10+9]
			   + eta2*tmp.gsds[ 9*6+5] + ze23*fsfs[1][6*10+9];
	    GSGS[ 9*15+14]+=QC[2]*gsfs[0][ 9*10+9]+WQ[2]*gsfs[1][ 9*10+9]
			   + eta2*tmp.gsds[ 9*6+3];
	    GSGS[10*15+ 0]+=QC[0]*gsfs[0][10*10+0]+WQ[0]*gsfs[1][10*10+0]
			   + eta23*tmp.gsds[10*6+0];
	    GSGS[10*15+ 1]+=QC[1]*gsfs[0][10*10+1]+WQ[1]*gsfs[1][10*10+1]
			   + eta23*tmp.gsds[10*6+1] + ze2*fsfs[1][2*10+1];
	    GSGS[10*15+ 2]+=QC[2]*gsfs[0][10*10+2]+WQ[2]*gsfs[1][10*10+2]
			   + eta23*tmp.gsds[10*6+2] + ze23*fsfs[1][7*10+2];
	    GSGS[10*15+ 3]+=QC[0]*gsfs[0][10*10+3]+WQ[0]*gsfs[1][10*10+3]
			   + eta*tmp.gsds[10*6+3];
	    GSGS[10*15+ 4]+=QC[1]*gsfs[0][10*10+4]+WQ[1]*gsfs[1][10*10+4]
			   + eta*tmp.gsds[10*6+4] + ze2*fsfs[1][2*10+4];
	    GSGS[10*15+ 5]+=QC[2]*gsfs[0][10*10+5]+WQ[2]*gsfs[1][10*10+5]
			   + eta*tmp.gsds[10*6+5] + ze23*fsfs[1][7*10+5];
	    GSGS[10*15+ 6]+=QC[0]*gsfs[0][10*10+6]+WQ[0]*gsfs[1][10*10+6]
			   + eta2*tmp.gsds[10*6+1];
	    GSGS[10*15+ 7]+=QC[1]*gsfs[0][10*10+7]+WQ[1]*gsfs[1][10*10+7]
			   + eta2*tmp.gsds[10*6+2] + ze2*fsfs[1][2*10+7];
	    GSGS[10*15+ 8]+=QC[2]*gsfs[0][10*10+8]+WQ[2]*gsfs[1][10*10+8]
			   + eta2*tmp.gsds[10*6+0] + ze23*fsfs[1][7*10+8];
	    GSGS[10*15+ 9]+=QC[0]*gsfs[0][10*10+1]+WQ[0]*gsfs[1][10*10+1];
	    GSGS[10*15+10]+=QC[1]*gsfs[0][10*10+2]+WQ[1]*gsfs[1][10*10+2]
			   + ze2*fsfs[1][2*10+2];
	    GSGS[10*15+11]+=QC[2]*gsfs[0][10*10+0]+WQ[2]*gsfs[1][10*10+0]
			   + ze23*fsfs[1][7*10+0];
	    GSGS[10*15+12]+=QC[0]*gsfs[0][10*10+9]+WQ[0]*gsfs[1][10*10+9]
			   + eta2*tmp.gsds[10*6+4];
	    GSGS[10*15+13]+=QC[1]*gsfs[0][10*10+9]+WQ[1]*gsfs[1][10*10+9]
			   + eta2*tmp.gsds[10*6+5] + ze2*fsfs[1][2*10+9];
	    GSGS[10*15+14]+=QC[2]*gsfs[0][10*10+9]+WQ[2]*gsfs[1][10*10+9]
			   + eta2*tmp.gsds[10*6+3] + ze23*fsfs[1][7*10+9];
	    GSGS[11*15+ 0]+=QC[0]*gsfs[0][11*10+0]+WQ[0]*gsfs[1][11*10+0]
			   + eta23*tmp.gsds[11*6+0] + ze23*fsfs[1][8*10+0];
	    GSGS[11*15+ 1]+=QC[1]*gsfs[0][11*10+1]+WQ[1]*gsfs[1][11*10+1]
			   + eta23*tmp.gsds[11*6+1];
	    GSGS[11*15+ 2]+=QC[2]*gsfs[0][11*10+2]+WQ[2]*gsfs[1][11*10+2]
			   + eta23*tmp.gsds[11*6+2] + ze2*fsfs[1][0*10+2];
	    GSGS[11*15+ 3]+=QC[0]*gsfs[0][11*10+3]+WQ[0]*gsfs[1][11*10+3]
			   + eta*tmp.gsds[11*6+3] + ze23*fsfs[1][8*10+3];
	    GSGS[11*15+ 4]+=QC[1]*gsfs[0][11*10+4]+WQ[1]*gsfs[1][11*10+4]
			   + eta*tmp.gsds[11*6+4];
	    GSGS[11*15+ 5]+=QC[2]*gsfs[0][11*10+5]+WQ[2]*gsfs[1][11*10+5]
			   + eta*tmp.gsds[11*6+5] + ze2*fsfs[1][0*10+5];
	    GSGS[11*15+ 6]+=QC[0]*gsfs[0][11*10+6]+WQ[0]*gsfs[1][11*10+6]
			   + eta2*tmp.gsds[11*6+1] + ze23*fsfs[1][8*10+6];
	    GSGS[11*15+ 7]+=QC[1]*gsfs[0][11*10+7]+WQ[1]*gsfs[1][11*10+7]
			   + eta2*tmp.gsds[11*6+2];
	    GSGS[11*15+ 8]+=QC[2]*gsfs[0][11*10+8]+WQ[2]*gsfs[1][11*10+8]
			   + eta2*tmp.gsds[11*6+0] + ze2*fsfs[1][0*10+8];
	    GSGS[11*15+ 9]+=QC[0]*gsfs[0][11*10+1]+WQ[0]*gsfs[1][11*10+1]
			   + ze23*fsfs[1][8*10+1];
	    GSGS[11*15+10]+=QC[1]*gsfs[0][11*10+2]+WQ[1]*gsfs[1][11*10+2];
	    GSGS[11*15+11]+=QC[2]*gsfs[0][11*10+0]+WQ[2]*gsfs[1][11*10+0]
			   + ze2*fsfs[1][0*10+0];
	    GSGS[11*15+12]+=QC[0]*gsfs[0][11*10+9]+WQ[0]*gsfs[1][11*10+9]
			   + eta2*tmp.gsds[11*6+4] + ze23*fsfs[1][8*10+9];
	    GSGS[11*15+13]+=QC[1]*gsfs[0][11*10+9]+WQ[1]*gsfs[1][11*10+9]
			   + eta2*tmp.gsds[11*6+5];
	    GSGS[11*15+14]+=QC[2]*gsfs[0][11*10+9]+WQ[2]*gsfs[1][11*10+9]
			   + eta2*tmp.gsds[11*6+3] + ze2*fsfs[1][0*10+9];
	    GSGS[12*15+ 0]+=QC[0]*gsfs[0][12*10+0]+WQ[0]*gsfs[1][12*10+0]
			   + eta23*tmp.gsds[12*6+0] + ze22*fsfs[1][9*10+0];
	    GSGS[12*15+ 1]+=QC[1]*gsfs[0][12*10+1]+WQ[1]*gsfs[1][12*10+1]
			   + eta23*tmp.gsds[12*6+1] + ze2*fsfs[1][8*10+1];
	    GSGS[12*15+ 2]+=QC[2]*gsfs[0][12*10+2]+WQ[2]*gsfs[1][12*10+2]
			   + eta23*tmp.gsds[12*6+2] + ze2*fsfs[1][3*10+2];
	    GSGS[12*15+ 3]+=QC[0]*gsfs[0][12*10+3]+WQ[0]*gsfs[1][12*10+3]
			   + eta*tmp.gsds[12*6+3] + ze22*fsfs[1][9*10+3];
	    GSGS[12*15+ 4]+=QC[1]*gsfs[0][12*10+4]+WQ[1]*gsfs[1][12*10+4]
			   + eta*tmp.gsds[12*6+4] + ze2*fsfs[1][8*10+4];
	    GSGS[12*15+ 5]+=QC[2]*gsfs[0][12*10+5]+WQ[2]*gsfs[1][12*10+5]
			   + eta*tmp.gsds[12*6+5] + ze2*fsfs[1][3*10+5];
	    GSGS[12*15+ 6]+=QC[0]*gsfs[0][12*10+6]+WQ[0]*gsfs[1][12*10+6]
			   + eta2*tmp.gsds[12*6+1] + ze22*fsfs[1][9*10+6];
	    GSGS[12*15+ 7]+=QC[1]*gsfs[0][12*10+7]+WQ[1]*gsfs[1][12*10+7]
			   + eta2*tmp.gsds[12*6+2] + ze2*fsfs[1][8*10+7];
	    GSGS[12*15+ 8]+=QC[2]*gsfs[0][12*10+8]+WQ[2]*gsfs[1][12*10+8]
			   + eta2*tmp.gsds[12*6+0] + ze2*fsfs[1][3*10+8];
	    GSGS[12*15+ 9]+=QC[0]*gsfs[0][12*10+1]+WQ[0]*gsfs[1][12*10+1]
			   + ze22*fsfs[1][9*10+1];
	    GSGS[12*15+10]+=QC[1]*gsfs[0][12*10+2]+WQ[1]*gsfs[1][12*10+2]
			   + ze2*fsfs[1][8*10+2];
	    GSGS[12*15+11]+=QC[2]*gsfs[0][12*10+0]+WQ[2]*gsfs[1][12*10+0]
			   + ze2*fsfs[1][3*10+0];
	    GSGS[12*15+12]+=QC[0]*gsfs[0][12*10+9]+WQ[0]*gsfs[1][12*10+9]
			   + eta2*tmp.gsds[12*6+4] + ze22*fsfs[1][9*10+9];
	    GSGS[12*15+13]+=QC[1]*gsfs[0][12*10+9]+WQ[1]*gsfs[1][12*10+9]
			   + eta2*tmp.gsds[12*6+5] + ze2*fsfs[1][8*10+9];
	    GSGS[12*15+14]+=QC[2]*gsfs[0][12*10+9]+WQ[2]*gsfs[1][12*10+9]
			   + eta2*tmp.gsds[12*6+3] + ze2*fsfs[1][3*10+9];
	    GSGS[13*15+ 0]+=QC[0]*gsfs[0][13*10+0]+WQ[0]*gsfs[1][13*10+0]
			   + eta23*tmp.gsds[13*6+0] + ze2*fsfs[1][4*10+0];
	    GSGS[13*15+ 1]+=QC[1]*gsfs[0][13*10+1]+WQ[1]*gsfs[1][13*10+1]
			   + eta23*tmp.gsds[13*6+1] + ze22*fsfs[1][9*10+1];
	    GSGS[13*15+ 2]+=QC[2]*gsfs[0][13*10+2]+WQ[2]*gsfs[1][13*10+2]
			   + eta23*tmp.gsds[13*6+2] + ze2*fsfs[1][6*10+2];
	    GSGS[13*15+ 3]+=QC[0]*gsfs[0][13*10+3]+WQ[0]*gsfs[1][13*10+3]
			   + eta*tmp.gsds[13*6+3] + ze2*fsfs[1][4*10+3];
	    GSGS[13*15+ 4]+=QC[1]*gsfs[0][13*10+4]+WQ[1]*gsfs[1][13*10+4]
			   + eta*tmp.gsds[13*6+4] + ze22*fsfs[1][9*10+4];
	    GSGS[13*15+ 5]+=QC[2]*gsfs[0][13*10+5]+WQ[2]*gsfs[1][13*10+5]
			   + eta*tmp.gsds[13*6+5] + ze2*fsfs[1][6*10+5];
	    GSGS[13*15+ 6]+=QC[0]*gsfs[0][13*10+6]+WQ[0]*gsfs[1][13*10+6]
			   + eta2*tmp.gsds[13*6+1] + ze2*fsfs[1][4*10+6];
	    GSGS[13*15+ 7]+=QC[1]*gsfs[0][13*10+7]+WQ[1]*gsfs[1][13*10+7]
			   + eta2*tmp.gsds[13*6+2] + ze22*fsfs[1][9*10+7];
	    GSGS[13*15+ 8]+=QC[2]*gsfs[0][13*10+8]+WQ[2]*gsfs[1][13*10+8]
			   + eta2*tmp.gsds[13*6+0] + ze2*fsfs[1][6*10+8];
	    GSGS[13*15+ 9]+=QC[0]*gsfs[0][13*10+1]+WQ[0]*gsfs[1][13*10+1]
			   + ze2*fsfs[1][4*10+1];
	    GSGS[13*15+10]+=QC[1]*gsfs[0][13*10+2]+WQ[1]*gsfs[1][13*10+2]
			   + ze22*fsfs[1][9*10+2];
	    GSGS[13*15+11]+=QC[2]*gsfs[0][13*10+0]+WQ[2]*gsfs[1][13*10+0]
			   + ze2*fsfs[1][6*10+0];
	    GSGS[13*15+12]+=QC[0]*gsfs[0][13*10+9]+WQ[0]*gsfs[1][13*10+9]
			   + eta2*tmp.gsds[13*6+4] + ze2*fsfs[1][4*10+9];
	    GSGS[13*15+13]+=QC[1]*gsfs[0][13*10+9]+WQ[1]*gsfs[1][13*10+9]
			   + eta2*tmp.gsds[13*6+5] + ze22*fsfs[1][9*10+9];
	    GSGS[13*15+14]+=QC[2]*gsfs[0][13*10+9]+WQ[2]*gsfs[1][13*10+9]
			   + eta2*tmp.gsds[13*6+3] + ze2*fsfs[1][6*10+9];
	    GSGS[14*15+ 0]+=QC[0]*gsfs[0][14*10+0]+WQ[0]*gsfs[1][14*10+0]
			   + eta23*tmp.gsds[14*6+0] + ze2*fsfs[1][7*10+0];
	    GSGS[14*15+ 1]+=QC[1]*gsfs[0][14*10+1]+WQ[1]*gsfs[1][14*10+1]
			   + eta23*tmp.gsds[14*6+1] + ze2*fsfs[1][5*10+1];
	    GSGS[14*15+ 2]+=QC[2]*gsfs[0][14*10+2]+WQ[2]*gsfs[1][14*10+2]
			   + eta23*tmp.gsds[14*6+2] + ze22*fsfs[1][9*10+2];
	    GSGS[14*15+ 3]+=QC[0]*gsfs[0][14*10+3]+WQ[0]*gsfs[1][14*10+3]
			   + eta*tmp.gsds[14*6+3] + ze2*fsfs[1][7*10+3];
	    GSGS[14*15+ 4]+=QC[1]*gsfs[0][14*10+4]+WQ[1]*gsfs[1][14*10+4]
			   + eta*tmp.gsds[14*6+4] + ze2*fsfs[1][5*10+4];
	    GSGS[14*15+ 5]+=QC[2]*gsfs[0][14*10+5]+WQ[2]*gsfs[1][14*10+5]
			   + eta*tmp.gsds[14*6+5] + ze22*fsfs[1][9*10+5];
	    GSGS[14*15+ 6]+=QC[0]*gsfs[0][14*10+6]+WQ[0]*gsfs[1][14*10+6]
			   + eta2*tmp.gsds[14*6+1] + ze2*fsfs[1][7*10+6];
	    GSGS[14*15+ 7]+=QC[1]*gsfs[0][14*10+7]+WQ[1]*gsfs[1][14*10+7]
			   + eta2*tmp.gsds[14*6+2] + ze2*fsfs[1][5*10+7];
	    GSGS[14*15+ 8]+=QC[2]*gsfs[0][14*10+8]+WQ[2]*gsfs[1][14*10+8]
			   + eta2*tmp.gsds[14*6+0] + ze22*fsfs[1][9*10+8];
	    GSGS[14*15+ 9]+=QC[0]*gsfs[0][14*10+1]+WQ[0]*gsfs[1][14*10+1]
			   + ze2*fsfs[1][7*10+1];
	    GSGS[14*15+10]+=QC[1]*gsfs[0][14*10+2]+WQ[1]*gsfs[1][14*10+2]
			   + ze2*fsfs[1][5*10+2];
	    GSGS[14*15+11]+=QC[2]*gsfs[0][14*10+0]+WQ[2]*gsfs[1][14*10+0]
			   + ze22*fsfs[1][9*10+0];
	    GSGS[14*15+12]+=QC[0]*gsfs[0][14*10+9]+WQ[0]*gsfs[1][14*10+9]
			   + eta2*tmp.gsds[14*6+4] + ze2*fsfs[1][7*10+9];
	    GSGS[14*15+13]+=QC[1]*gsfs[0][14*10+9]+WQ[1]*gsfs[1][14*10+9]
			   + eta2*tmp.gsds[14*6+5] + ze2*fsfs[1][5*10+9];
	    GSGS[14*15+14]+=QC[2]*gsfs[0][14*10+9]+WQ[2]*gsfs[1][14*10+9]
			   + eta2*tmp.gsds[14*6+3] + ze22*fsfs[1][9*10+9];
	}	// klps
    }	// ijps
    // (D,P|D,S)
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
    // (D,P|F,S)
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
    // (D,P|G,S)
    for (c0=0; c0<15; c0++) {
        DPGS[0*3*15+0*15+c0] = FSGS[0*15+c0] - BA[0]*DSGS[0*15+c0];
        DPGS[0*3*15+1*15+c0] = FSGS[3*15+c0] - BA[1]*DSGS[0*15+c0];
        DPGS[0*3*15+2*15+c0] = FSGS[8*15+c0] - BA[2]*DSGS[0*15+c0];
        DPGS[1*3*15+0*15+c0] = FSGS[6*15+c0] - BA[0]*DSGS[1*15+c0];
        DPGS[1*3*15+1*15+c0] = FSGS[1*15+c0] - BA[1]*DSGS[1*15+c0];
        DPGS[1*3*15+2*15+c0] = FSGS[4*15+c0] - BA[2]*DSGS[1*15+c0];
        DPGS[2*3*15+0*15+c0] = FSGS[5*15+c0] - BA[0]*DSGS[2*15+c0];
        DPGS[2*3*15+1*15+c0] = FSGS[7*15+c0] - BA[1]*DSGS[2*15+c0];
        DPGS[2*3*15+2*15+c0] = FSGS[2*15+c0] - BA[2]*DSGS[2*15+c0];
        DPGS[3*3*15+0*15+c0] = FSGS[3*15+c0] - BA[0]*DSGS[3*15+c0];
        DPGS[3*3*15+1*15+c0] = FSGS[6*15+c0] - BA[1]*DSGS[3*15+c0];
        DPGS[3*3*15+2*15+c0] = FSGS[9*15+c0] - BA[2]*DSGS[3*15+c0];
        DPGS[4*3*15+0*15+c0] = FSGS[9*15+c0] - BA[0]*DSGS[4*15+c0];
        DPGS[4*3*15+1*15+c0] = FSGS[4*15+c0] - BA[1]*DSGS[4*15+c0];
        DPGS[4*3*15+2*15+c0] = FSGS[7*15+c0] - BA[2]*DSGS[4*15+c0];
        DPGS[5*3*15+0*15+c0] = FSGS[8*15+c0] - BA[0]*DSGS[5*15+c0];
        DPGS[5*3*15+1*15+c0] = FSGS[9*15+c0] - BA[1]*DSGS[5*15+c0];
        DPGS[5*3*15+2*15+c0] = FSGS[5*15+c0] - BA[2]*DSGS[5*15+c0];
    }
    // (F,P|D,S)
    for (c0=0; c0<6; c0++) {
        FPDS[0*3*6+0*6+c0] = GSDS[ 0*6+c0] - BA[0]*FSDS[0*6+c0];
        FPDS[0*3*6+1*6+c0] = GSDS[ 3*6+c0] - BA[1]*FSDS[0*6+c0];
        FPDS[0*3*6+2*6+c0] = GSDS[11*6+c0] - BA[2]*FSDS[0*6+c0];
        FPDS[1*3*6+0*6+c0] = GSDS[ 9*6+c0] - BA[0]*FSDS[1*6+c0];
        FPDS[1*3*6+1*6+c0] = GSDS[ 1*6+c0] - BA[1]*FSDS[1*6+c0];
        FPDS[1*3*6+2*6+c0] = GSDS[ 4*6+c0] - BA[2]*FSDS[1*6+c0];
        FPDS[2*3*6+0*6+c0] = GSDS[ 5*6+c0] - BA[0]*FSDS[2*6+c0];
        FPDS[2*3*6+1*6+c0] = GSDS[10*6+c0] - BA[1]*FSDS[2*6+c0];
        FPDS[2*3*6+2*6+c0] = GSDS[ 2*6+c0] - BA[2]*FSDS[2*6+c0];
        FPDS[3*3*6+0*6+c0] = GSDS[ 3*6+c0] - BA[0]*FSDS[3*6+c0];
        FPDS[3*3*6+1*6+c0] = GSDS[ 6*6+c0] - BA[1]*FSDS[3*6+c0];
        FPDS[3*3*6+2*6+c0] = GSDS[12*6+c0] - BA[2]*FSDS[3*6+c0];
        FPDS[4*3*6+0*6+c0] = GSDS[13*6+c0] - BA[0]*FSDS[4*6+c0];
        FPDS[4*3*6+1*6+c0] = GSDS[ 4*6+c0] - BA[1]*FSDS[4*6+c0];
        FPDS[4*3*6+2*6+c0] = GSDS[ 7*6+c0] - BA[2]*FSDS[4*6+c0];
        FPDS[5*3*6+0*6+c0] = GSDS[ 8*6+c0] - BA[0]*FSDS[5*6+c0];
        FPDS[5*3*6+1*6+c0] = GSDS[14*6+c0] - BA[1]*FSDS[5*6+c0];
        FPDS[5*3*6+2*6+c0] = GSDS[ 5*6+c0] - BA[2]*FSDS[5*6+c0];
        FPDS[6*3*6+0*6+c0] = GSDS[ 6*6+c0] - BA[0]*FSDS[6*6+c0];
        FPDS[6*3*6+1*6+c0] = GSDS[ 9*6+c0] - BA[1]*FSDS[6*6+c0];
        FPDS[6*3*6+2*6+c0] = GSDS[13*6+c0] - BA[2]*FSDS[6*6+c0];
        FPDS[7*3*6+0*6+c0] = GSDS[14*6+c0] - BA[0]*FSDS[7*6+c0];
        FPDS[7*3*6+1*6+c0] = GSDS[ 7*6+c0] - BA[1]*FSDS[7*6+c0];
        FPDS[7*3*6+2*6+c0] = GSDS[10*6+c0] - BA[2]*FSDS[7*6+c0];
        FPDS[8*3*6+0*6+c0] = GSDS[11*6+c0] - BA[0]*FSDS[8*6+c0];
        FPDS[8*3*6+1*6+c0] = GSDS[12*6+c0] - BA[1]*FSDS[8*6+c0];
        FPDS[8*3*6+2*6+c0] = GSDS[ 8*6+c0] - BA[2]*FSDS[8*6+c0];
        FPDS[9*3*6+0*6+c0] = GSDS[12*6+c0] - BA[0]*FSDS[9*6+c0];
        FPDS[9*3*6+1*6+c0] = GSDS[13*6+c0] - BA[1]*FSDS[9*6+c0];
        FPDS[9*3*6+2*6+c0] = GSDS[14*6+c0] - BA[2]*FSDS[9*6+c0];
    }
    // (F,P|F,S)
    for (c0=0; c0<10; c0++) {
        FPFS[0*3*10+0*10+c0] = GSFS[ 0*10+c0] - BA[0]*FSFS[0*10+c0];
        FPFS[0*3*10+1*10+c0] = GSFS[ 3*10+c0] - BA[1]*FSFS[0*10+c0];
        FPFS[0*3*10+2*10+c0] = GSFS[11*10+c0] - BA[2]*FSFS[0*10+c0];
        FPFS[1*3*10+0*10+c0] = GSFS[ 9*10+c0] - BA[0]*FSFS[1*10+c0];
        FPFS[1*3*10+1*10+c0] = GSFS[ 1*10+c0] - BA[1]*FSFS[1*10+c0];
        FPFS[1*3*10+2*10+c0] = GSFS[ 4*10+c0] - BA[2]*FSFS[1*10+c0];
        FPFS[2*3*10+0*10+c0] = GSFS[ 5*10+c0] - BA[0]*FSFS[2*10+c0];
        FPFS[2*3*10+1*10+c0] = GSFS[10*10+c0] - BA[1]*FSFS[2*10+c0];
        FPFS[2*3*10+2*10+c0] = GSFS[ 2*10+c0] - BA[2]*FSFS[2*10+c0];
        FPFS[3*3*10+0*10+c0] = GSFS[ 3*10+c0] - BA[0]*FSFS[3*10+c0];
        FPFS[3*3*10+1*10+c0] = GSFS[ 6*10+c0] - BA[1]*FSFS[3*10+c0];
        FPFS[3*3*10+2*10+c0] = GSFS[12*10+c0] - BA[2]*FSFS[3*10+c0];
        FPFS[4*3*10+0*10+c0] = GSFS[13*10+c0] - BA[0]*FSFS[4*10+c0];
        FPFS[4*3*10+1*10+c0] = GSFS[ 4*10+c0] - BA[1]*FSFS[4*10+c0];
        FPFS[4*3*10+2*10+c0] = GSFS[ 7*10+c0] - BA[2]*FSFS[4*10+c0];
        FPFS[5*3*10+0*10+c0] = GSFS[ 8*10+c0] - BA[0]*FSFS[5*10+c0];
        FPFS[5*3*10+1*10+c0] = GSFS[14*10+c0] - BA[1]*FSFS[5*10+c0];
        FPFS[5*3*10+2*10+c0] = GSFS[ 5*10+c0] - BA[2]*FSFS[5*10+c0];
        FPFS[6*3*10+0*10+c0] = GSFS[ 6*10+c0] - BA[0]*FSFS[6*10+c0];
        FPFS[6*3*10+1*10+c0] = GSFS[ 9*10+c0] - BA[1]*FSFS[6*10+c0];
        FPFS[6*3*10+2*10+c0] = GSFS[13*10+c0] - BA[2]*FSFS[6*10+c0];
        FPFS[7*3*10+0*10+c0] = GSFS[14*10+c0] - BA[0]*FSFS[7*10+c0];
        FPFS[7*3*10+1*10+c0] = GSFS[ 7*10+c0] - BA[1]*FSFS[7*10+c0];
        FPFS[7*3*10+2*10+c0] = GSFS[10*10+c0] - BA[2]*FSFS[7*10+c0];
        FPFS[8*3*10+0*10+c0] = GSFS[11*10+c0] - BA[0]*FSFS[8*10+c0];
        FPFS[8*3*10+1*10+c0] = GSFS[12*10+c0] - BA[1]*FSFS[8*10+c0];
        FPFS[8*3*10+2*10+c0] = GSFS[ 8*10+c0] - BA[2]*FSFS[8*10+c0];
        FPFS[9*3*10+0*10+c0] = GSFS[12*10+c0] - BA[0]*FSFS[9*10+c0];
        FPFS[9*3*10+1*10+c0] = GSFS[13*10+c0] - BA[1]*FSFS[9*10+c0];
        FPFS[9*3*10+2*10+c0] = GSFS[14*10+c0] - BA[2]*FSFS[9*10+c0];
    }
    // (F,P|G,S)
    for (c0=0; c0<15; c0++) {
        FPGS[0*3*15+0*15+c0] = GSGS[ 0*15+c0] - BA[0]*FSGS[0*15+c0];
        FPGS[0*3*15+1*15+c0] = GSGS[ 3*15+c0] - BA[1]*FSGS[0*15+c0];
        FPGS[0*3*15+2*15+c0] = GSGS[11*15+c0] - BA[2]*FSGS[0*15+c0];
        FPGS[1*3*15+0*15+c0] = GSGS[ 9*15+c0] - BA[0]*FSGS[1*15+c0];
        FPGS[1*3*15+1*15+c0] = GSGS[ 1*15+c0] - BA[1]*FSGS[1*15+c0];
        FPGS[1*3*15+2*15+c0] = GSGS[ 4*15+c0] - BA[2]*FSGS[1*15+c0];
        FPGS[2*3*15+0*15+c0] = GSGS[ 5*15+c0] - BA[0]*FSGS[2*15+c0];
        FPGS[2*3*15+1*15+c0] = GSGS[10*15+c0] - BA[1]*FSGS[2*15+c0];
        FPGS[2*3*15+2*15+c0] = GSGS[ 2*15+c0] - BA[2]*FSGS[2*15+c0];
        FPGS[3*3*15+0*15+c0] = GSGS[ 3*15+c0] - BA[0]*FSGS[3*15+c0];
        FPGS[3*3*15+1*15+c0] = GSGS[ 6*15+c0] - BA[1]*FSGS[3*15+c0];
        FPGS[3*3*15+2*15+c0] = GSGS[12*15+c0] - BA[2]*FSGS[3*15+c0];
        FPGS[4*3*15+0*15+c0] = GSGS[13*15+c0] - BA[0]*FSGS[4*15+c0];
        FPGS[4*3*15+1*15+c0] = GSGS[ 4*15+c0] - BA[1]*FSGS[4*15+c0];
        FPGS[4*3*15+2*15+c0] = GSGS[ 7*15+c0] - BA[2]*FSGS[4*15+c0];
        FPGS[5*3*15+0*15+c0] = GSGS[ 8*15+c0] - BA[0]*FSGS[5*15+c0];
        FPGS[5*3*15+1*15+c0] = GSGS[14*15+c0] - BA[1]*FSGS[5*15+c0];
        FPGS[5*3*15+2*15+c0] = GSGS[ 5*15+c0] - BA[2]*FSGS[5*15+c0];
        FPGS[6*3*15+0*15+c0] = GSGS[ 6*15+c0] - BA[0]*FSGS[6*15+c0];
        FPGS[6*3*15+1*15+c0] = GSGS[ 9*15+c0] - BA[1]*FSGS[6*15+c0];
        FPGS[6*3*15+2*15+c0] = GSGS[13*15+c0] - BA[2]*FSGS[6*15+c0];
        FPGS[7*3*15+0*15+c0] = GSGS[14*15+c0] - BA[0]*FSGS[7*15+c0];
        FPGS[7*3*15+1*15+c0] = GSGS[ 7*15+c0] - BA[1]*FSGS[7*15+c0];
        FPGS[7*3*15+2*15+c0] = GSGS[10*15+c0] - BA[2]*FSGS[7*15+c0];
        FPGS[8*3*15+0*15+c0] = GSGS[11*15+c0] - BA[0]*FSGS[8*15+c0];
        FPGS[8*3*15+1*15+c0] = GSGS[12*15+c0] - BA[1]*FSGS[8*15+c0];
        FPGS[8*3*15+2*15+c0] = GSGS[ 8*15+c0] - BA[2]*FSGS[8*15+c0];
        FPGS[9*3*15+0*15+c0] = GSGS[12*15+c0] - BA[0]*FSGS[9*15+c0];
        FPGS[9*3*15+1*15+c0] = GSGS[13*15+c0] - BA[1]*FSGS[9*15+c0];
        FPGS[9*3*15+2*15+c0] = GSGS[14*15+c0] - BA[2]*FSGS[9*15+c0];
    }
    // (D,D|D,S)
    for (c0=0; c0<6; c0++) {
        DDDS[0*6*6+0*6+c0] = FPDS[0*3*6+0*6+c0] - BA[0]*DPDS[0*3*6+0*6+c0];
        DDDS[0*6*6+1*6+c0] = FPDS[3*3*6+1*6+c0] - BA[1]*DPDS[0*3*6+1*6+c0];
        DDDS[0*6*6+2*6+c0] = FPDS[8*3*6+2*6+c0] - BA[2]*DPDS[0*3*6+2*6+c0];
        DDDS[0*6*6+3*6+c0] = FPDS[0*3*6+1*6+c0] - BA[0]*DPDS[0*3*6+1*6+c0];
        DDDS[0*6*6+4*6+c0] = FPDS[3*3*6+2*6+c0] - BA[1]*DPDS[0*3*6+2*6+c0];
        DDDS[0*6*6+5*6+c0] = FPDS[8*3*6+0*6+c0] - BA[2]*DPDS[0*3*6+0*6+c0];
        DDDS[1*6*6+0*6+c0] = FPDS[6*3*6+0*6+c0] - BA[0]*DPDS[1*3*6+0*6+c0];
        DDDS[1*6*6+1*6+c0] = FPDS[1*3*6+1*6+c0] - BA[1]*DPDS[1*3*6+1*6+c0];
        DDDS[1*6*6+2*6+c0] = FPDS[4*3*6+2*6+c0] - BA[2]*DPDS[1*3*6+2*6+c0];
        DDDS[1*6*6+3*6+c0] = FPDS[6*3*6+1*6+c0] - BA[0]*DPDS[1*3*6+1*6+c0];
        DDDS[1*6*6+4*6+c0] = FPDS[1*3*6+2*6+c0] - BA[1]*DPDS[1*3*6+2*6+c0];
        DDDS[1*6*6+5*6+c0] = FPDS[4*3*6+0*6+c0] - BA[2]*DPDS[1*3*6+0*6+c0];
        DDDS[2*6*6+0*6+c0] = FPDS[5*3*6+0*6+c0] - BA[0]*DPDS[2*3*6+0*6+c0];
        DDDS[2*6*6+1*6+c0] = FPDS[7*3*6+1*6+c0] - BA[1]*DPDS[2*3*6+1*6+c0];
        DDDS[2*6*6+2*6+c0] = FPDS[2*3*6+2*6+c0] - BA[2]*DPDS[2*3*6+2*6+c0];
        DDDS[2*6*6+3*6+c0] = FPDS[5*3*6+1*6+c0] - BA[0]*DPDS[2*3*6+1*6+c0];
        DDDS[2*6*6+4*6+c0] = FPDS[7*3*6+2*6+c0] - BA[1]*DPDS[2*3*6+2*6+c0];
        DDDS[2*6*6+5*6+c0] = FPDS[2*3*6+0*6+c0] - BA[2]*DPDS[2*3*6+0*6+c0];
        DDDS[3*6*6+0*6+c0] = FPDS[3*3*6+0*6+c0] - BA[0]*DPDS[3*3*6+0*6+c0];
        DDDS[3*6*6+1*6+c0] = FPDS[6*3*6+1*6+c0] - BA[1]*DPDS[3*3*6+1*6+c0];
        DDDS[3*6*6+2*6+c0] = FPDS[9*3*6+2*6+c0] - BA[2]*DPDS[3*3*6+2*6+c0];
        DDDS[3*6*6+3*6+c0] = FPDS[3*3*6+1*6+c0] - BA[0]*DPDS[3*3*6+1*6+c0];
        DDDS[3*6*6+4*6+c0] = FPDS[6*3*6+2*6+c0] - BA[1]*DPDS[3*3*6+2*6+c0];
        DDDS[3*6*6+5*6+c0] = FPDS[9*3*6+0*6+c0] - BA[2]*DPDS[3*3*6+0*6+c0];
        DDDS[4*6*6+0*6+c0] = FPDS[9*3*6+0*6+c0] - BA[0]*DPDS[4*3*6+0*6+c0];
        DDDS[4*6*6+1*6+c0] = FPDS[4*3*6+1*6+c0] - BA[1]*DPDS[4*3*6+1*6+c0];
        DDDS[4*6*6+2*6+c0] = FPDS[7*3*6+2*6+c0] - BA[2]*DPDS[4*3*6+2*6+c0];
        DDDS[4*6*6+3*6+c0] = FPDS[9*3*6+1*6+c0] - BA[0]*DPDS[4*3*6+1*6+c0];
        DDDS[4*6*6+4*6+c0] = FPDS[4*3*6+2*6+c0] - BA[1]*DPDS[4*3*6+2*6+c0];
        DDDS[4*6*6+5*6+c0] = FPDS[7*3*6+0*6+c0] - BA[2]*DPDS[4*3*6+0*6+c0];
        DDDS[5*6*6+0*6+c0] = FPDS[8*3*6+0*6+c0] - BA[0]*DPDS[5*3*6+0*6+c0];
        DDDS[5*6*6+1*6+c0] = FPDS[9*3*6+1*6+c0] - BA[1]*DPDS[5*3*6+1*6+c0];
        DDDS[5*6*6+2*6+c0] = FPDS[5*3*6+2*6+c0] - BA[2]*DPDS[5*3*6+2*6+c0];
        DDDS[5*6*6+3*6+c0] = FPDS[8*3*6+1*6+c0] - BA[0]*DPDS[5*3*6+1*6+c0];
        DDDS[5*6*6+4*6+c0] = FPDS[9*3*6+2*6+c0] - BA[1]*DPDS[5*3*6+2*6+c0];
        DDDS[5*6*6+5*6+c0] = FPDS[5*3*6+0*6+c0] - BA[2]*DPDS[5*3*6+0*6+c0];
    }
    // (D,D|F,S)
    for (c0=0; c0<10; c0++) {
        DDFS[0*6*10+0*10+c0]=FPFS[0*3*10+0*10+c0]-BA[0]*DPFS[0*3*10+0*10+c0];
        DDFS[0*6*10+1*10+c0]=FPFS[3*3*10+1*10+c0]-BA[1]*DPFS[0*3*10+1*10+c0];
        DDFS[0*6*10+2*10+c0]=FPFS[8*3*10+2*10+c0]-BA[2]*DPFS[0*3*10+2*10+c0];
        DDFS[0*6*10+3*10+c0]=FPFS[0*3*10+1*10+c0]-BA[0]*DPFS[0*3*10+1*10+c0];
        DDFS[0*6*10+4*10+c0]=FPFS[3*3*10+2*10+c0]-BA[1]*DPFS[0*3*10+2*10+c0];
        DDFS[0*6*10+5*10+c0]=FPFS[8*3*10+0*10+c0]-BA[2]*DPFS[0*3*10+0*10+c0];
        DDFS[1*6*10+0*10+c0]=FPFS[6*3*10+0*10+c0]-BA[0]*DPFS[1*3*10+0*10+c0];
        DDFS[1*6*10+1*10+c0]=FPFS[1*3*10+1*10+c0]-BA[1]*DPFS[1*3*10+1*10+c0];
        DDFS[1*6*10+2*10+c0]=FPFS[4*3*10+2*10+c0]-BA[2]*DPFS[1*3*10+2*10+c0];
        DDFS[1*6*10+3*10+c0]=FPFS[6*3*10+1*10+c0]-BA[0]*DPFS[1*3*10+1*10+c0];
        DDFS[1*6*10+4*10+c0]=FPFS[1*3*10+2*10+c0]-BA[1]*DPFS[1*3*10+2*10+c0];
        DDFS[1*6*10+5*10+c0]=FPFS[4*3*10+0*10+c0]-BA[2]*DPFS[1*3*10+0*10+c0];
        DDFS[2*6*10+0*10+c0]=FPFS[5*3*10+0*10+c0]-BA[0]*DPFS[2*3*10+0*10+c0];
        DDFS[2*6*10+1*10+c0]=FPFS[7*3*10+1*10+c0]-BA[1]*DPFS[2*3*10+1*10+c0];
        DDFS[2*6*10+2*10+c0]=FPFS[2*3*10+2*10+c0]-BA[2]*DPFS[2*3*10+2*10+c0];
        DDFS[2*6*10+3*10+c0]=FPFS[5*3*10+1*10+c0]-BA[0]*DPFS[2*3*10+1*10+c0];
        DDFS[2*6*10+4*10+c0]=FPFS[7*3*10+2*10+c0]-BA[1]*DPFS[2*3*10+2*10+c0];
        DDFS[2*6*10+5*10+c0]=FPFS[2*3*10+0*10+c0]-BA[2]*DPFS[2*3*10+0*10+c0];
        DDFS[3*6*10+0*10+c0]=FPFS[3*3*10+0*10+c0]-BA[0]*DPFS[3*3*10+0*10+c0];
        DDFS[3*6*10+1*10+c0]=FPFS[6*3*10+1*10+c0]-BA[1]*DPFS[3*3*10+1*10+c0];
        DDFS[3*6*10+2*10+c0]=FPFS[9*3*10+2*10+c0]-BA[2]*DPFS[3*3*10+2*10+c0];
        DDFS[3*6*10+3*10+c0]=FPFS[3*3*10+1*10+c0]-BA[0]*DPFS[3*3*10+1*10+c0];
        DDFS[3*6*10+4*10+c0]=FPFS[6*3*10+2*10+c0]-BA[1]*DPFS[3*3*10+2*10+c0];
        DDFS[3*6*10+5*10+c0]=FPFS[9*3*10+0*10+c0]-BA[2]*DPFS[3*3*10+0*10+c0];
        DDFS[4*6*10+0*10+c0]=FPFS[9*3*10+0*10+c0]-BA[0]*DPFS[4*3*10+0*10+c0];
        DDFS[4*6*10+1*10+c0]=FPFS[4*3*10+1*10+c0]-BA[1]*DPFS[4*3*10+1*10+c0];
        DDFS[4*6*10+2*10+c0]=FPFS[7*3*10+2*10+c0]-BA[2]*DPFS[4*3*10+2*10+c0];
        DDFS[4*6*10+3*10+c0]=FPFS[9*3*10+1*10+c0]-BA[0]*DPFS[4*3*10+1*10+c0];
        DDFS[4*6*10+4*10+c0]=FPFS[4*3*10+2*10+c0]-BA[1]*DPFS[4*3*10+2*10+c0];
        DDFS[4*6*10+5*10+c0]=FPFS[7*3*10+0*10+c0]-BA[2]*DPFS[4*3*10+0*10+c0];
        DDFS[5*6*10+0*10+c0]=FPFS[8*3*10+0*10+c0]-BA[0]*DPFS[5*3*10+0*10+c0];
        DDFS[5*6*10+1*10+c0]=FPFS[9*3*10+1*10+c0]-BA[1]*DPFS[5*3*10+1*10+c0];
        DDFS[5*6*10+2*10+c0]=FPFS[5*3*10+2*10+c0]-BA[2]*DPFS[5*3*10+2*10+c0];
        DDFS[5*6*10+3*10+c0]=FPFS[8*3*10+1*10+c0]-BA[0]*DPFS[5*3*10+1*10+c0];
        DDFS[5*6*10+4*10+c0]=FPFS[9*3*10+2*10+c0]-BA[1]*DPFS[5*3*10+2*10+c0];
        DDFS[5*6*10+5*10+c0]=FPFS[5*3*10+0*10+c0]-BA[2]*DPFS[5*3*10+0*10+c0];
    }
    // (D,D|G,S)
    for (c0=0; c0<15; c0++) {
        DDGS[0*6*15+0*15+c0]=FPGS[0*3*15+0*15+c0]-BA[0]*DPGS[0*3*15+0*15+c0];
        DDGS[0*6*15+1*15+c0]=FPGS[3*3*15+1*15+c0]-BA[1]*DPGS[0*3*15+1*15+c0];
        DDGS[0*6*15+2*15+c0]=FPGS[8*3*15+2*15+c0]-BA[2]*DPGS[0*3*15+2*15+c0];
        DDGS[0*6*15+3*15+c0]=FPGS[0*3*15+1*15+c0]-BA[0]*DPGS[0*3*15+1*15+c0];
        DDGS[0*6*15+4*15+c0]=FPGS[3*3*15+2*15+c0]-BA[1]*DPGS[0*3*15+2*15+c0];
        DDGS[0*6*15+5*15+c0]=FPGS[8*3*15+0*15+c0]-BA[2]*DPGS[0*3*15+0*15+c0];
        DDGS[1*6*15+0*15+c0]=FPGS[6*3*15+0*15+c0]-BA[0]*DPGS[1*3*15+0*15+c0];
        DDGS[1*6*15+1*15+c0]=FPGS[1*3*15+1*15+c0]-BA[1]*DPGS[1*3*15+1*15+c0];
        DDGS[1*6*15+2*15+c0]=FPGS[4*3*15+2*15+c0]-BA[2]*DPGS[1*3*15+2*15+c0];
        DDGS[1*6*15+3*15+c0]=FPGS[6*3*15+1*15+c0]-BA[0]*DPGS[1*3*15+1*15+c0];
        DDGS[1*6*15+4*15+c0]=FPGS[1*3*15+2*15+c0]-BA[1]*DPGS[1*3*15+2*15+c0];
        DDGS[1*6*15+5*15+c0]=FPGS[4*3*15+0*15+c0]-BA[2]*DPGS[1*3*15+0*15+c0];
        DDGS[2*6*15+0*15+c0]=FPGS[5*3*15+0*15+c0]-BA[0]*DPGS[2*3*15+0*15+c0];
        DDGS[2*6*15+1*15+c0]=FPGS[7*3*15+1*15+c0]-BA[1]*DPGS[2*3*15+1*15+c0];
        DDGS[2*6*15+2*15+c0]=FPGS[2*3*15+2*15+c0]-BA[2]*DPGS[2*3*15+2*15+c0];
        DDGS[2*6*15+3*15+c0]=FPGS[5*3*15+1*15+c0]-BA[0]*DPGS[2*3*15+1*15+c0];
        DDGS[2*6*15+4*15+c0]=FPGS[7*3*15+2*15+c0]-BA[1]*DPGS[2*3*15+2*15+c0];
        DDGS[2*6*15+5*15+c0]=FPGS[2*3*15+0*15+c0]-BA[2]*DPGS[2*3*15+0*15+c0];
        DDGS[3*6*15+0*15+c0]=FPGS[3*3*15+0*15+c0]-BA[0]*DPGS[3*3*15+0*15+c0];
        DDGS[3*6*15+1*15+c0]=FPGS[6*3*15+1*15+c0]-BA[1]*DPGS[3*3*15+1*15+c0];
        DDGS[3*6*15+2*15+c0]=FPGS[9*3*15+2*15+c0]-BA[2]*DPGS[3*3*15+2*15+c0];
        DDGS[3*6*15+3*15+c0]=FPGS[3*3*15+1*15+c0]-BA[0]*DPGS[3*3*15+1*15+c0];
        DDGS[3*6*15+4*15+c0]=FPGS[6*3*15+2*15+c0]-BA[1]*DPGS[3*3*15+2*15+c0];
        DDGS[3*6*15+5*15+c0]=FPGS[9*3*15+0*15+c0]-BA[2]*DPGS[3*3*15+0*15+c0];
        DDGS[4*6*15+0*15+c0]=FPGS[9*3*15+0*15+c0]-BA[0]*DPGS[4*3*15+0*15+c0];
        DDGS[4*6*15+1*15+c0]=FPGS[4*3*15+1*15+c0]-BA[1]*DPGS[4*3*15+1*15+c0];
        DDGS[4*6*15+2*15+c0]=FPGS[7*3*15+2*15+c0]-BA[2]*DPGS[4*3*15+2*15+c0];
        DDGS[4*6*15+3*15+c0]=FPGS[9*3*15+1*15+c0]-BA[0]*DPGS[4*3*15+1*15+c0];
        DDGS[4*6*15+4*15+c0]=FPGS[4*3*15+2*15+c0]-BA[1]*DPGS[4*3*15+2*15+c0];
        DDGS[4*6*15+5*15+c0]=FPGS[7*3*15+0*15+c0]-BA[2]*DPGS[4*3*15+0*15+c0];
        DDGS[5*6*15+0*15+c0]=FPGS[8*3*15+0*15+c0]-BA[0]*DPGS[5*3*15+0*15+c0];
        DDGS[5*6*15+1*15+c0]=FPGS[9*3*15+1*15+c0]-BA[1]*DPGS[5*3*15+1*15+c0];
        DDGS[5*6*15+2*15+c0]=FPGS[5*3*15+2*15+c0]-BA[2]*DPGS[5*3*15+2*15+c0];
        DDGS[5*6*15+3*15+c0]=FPGS[8*3*15+1*15+c0]-BA[0]*DPGS[5*3*15+1*15+c0];
        DDGS[5*6*15+4*15+c0]=FPGS[9*3*15+2*15+c0]-BA[1]*DPGS[5*3*15+2*15+c0];
        DDGS[5*6*15+5*15+c0]=FPGS[5*3*15+0*15+c0]-BA[2]*DPGS[5*3*15+0*15+c0];
    }
    // (a,b|c,d+1)=(a,b|c+1,d)-DC(a,b|c,d)
    // (D,D|D,P)
    for (ab=0, ab01=0, ab10=0, ab00=0;
            ab<6*6; ab++, ab01 += 6*3, ab10 += 10, ab00 += 6) {
        DDDP[ab01+0*3+0]=DDFS[ab10+0]-DC[0]*DDDS[ab00+0];
        DDDP[ab01+0*3+1]=DDFS[ab10+3]-DC[1]*DDDS[ab00+0];
        DDDP[ab01+0*3+2]=DDFS[ab10+8]-DC[2]*DDDS[ab00+0];
        DDDP[ab01+1*3+0]=DDFS[ab10+6]-DC[0]*DDDS[ab00+1];
        DDDP[ab01+1*3+1]=DDFS[ab10+1]-DC[1]*DDDS[ab00+1];
        DDDP[ab01+1*3+2]=DDFS[ab10+4]-DC[2]*DDDS[ab00+1];
        DDDP[ab01+2*3+0]=DDFS[ab10+5]-DC[0]*DDDS[ab00+2];
        DDDP[ab01+2*3+1]=DDFS[ab10+7]-DC[1]*DDDS[ab00+2];
        DDDP[ab01+2*3+2]=DDFS[ab10+2]-DC[2]*DDDS[ab00+2];
        DDDP[ab01+3*3+0]=DDFS[ab10+3]-DC[0]*DDDS[ab00+3];
        DDDP[ab01+3*3+1]=DDFS[ab10+6]-DC[1]*DDDS[ab00+3];
        DDDP[ab01+3*3+2]=DDFS[ab10+9]-DC[2]*DDDS[ab00+3];
        DDDP[ab01+4*3+0]=DDFS[ab10+9]-DC[0]*DDDS[ab00+4];
        DDDP[ab01+4*3+1]=DDFS[ab10+4]-DC[1]*DDDS[ab00+4];
        DDDP[ab01+4*3+2]=DDFS[ab10+7]-DC[2]*DDDS[ab00+4];
        DDDP[ab01+5*3+0]=DDFS[ab10+8]-DC[0]*DDDS[ab00+5];
        DDDP[ab01+5*3+1]=DDFS[ab10+9]-DC[1]*DDDS[ab00+5];
        DDDP[ab01+5*3+2]=DDFS[ab10+5]-DC[2]*DDDS[ab00+5];
    }
    // (D,D|F,P)
    for (ab=0, ab01=0, ab10=0, ab00=0;
            ab<6*6; ab++, ab01 += 10*3, ab10 += 15, ab00 += 10) {
        DDFP[ab01+0*3+0]=DDGS[ab10+ 0]-DC[0]*DDFS[ab00+0];
        DDFP[ab01+0*3+1]=DDGS[ab10+ 3]-DC[1]*DDFS[ab00+0];
        DDFP[ab01+0*3+2]=DDGS[ab10+11]-DC[2]*DDFS[ab00+0];
        DDFP[ab01+1*3+0]=DDGS[ab10+ 9]-DC[0]*DDFS[ab00+1];
        DDFP[ab01+1*3+1]=DDGS[ab10+ 1]-DC[1]*DDFS[ab00+1];
        DDFP[ab01+1*3+2]=DDGS[ab10+ 4]-DC[2]*DDFS[ab00+1];
        DDFP[ab01+2*3+0]=DDGS[ab10+ 5]-DC[0]*DDFS[ab00+2];
        DDFP[ab01+2*3+1]=DDGS[ab10+10]-DC[1]*DDFS[ab00+2];
        DDFP[ab01+2*3+2]=DDGS[ab10+ 2]-DC[2]*DDFS[ab00+2];
        DDFP[ab01+3*3+0]=DDGS[ab10+ 3]-DC[0]*DDFS[ab00+3];
        DDFP[ab01+3*3+1]=DDGS[ab10+ 6]-DC[1]*DDFS[ab00+3];
        DDFP[ab01+3*3+2]=DDGS[ab10+12]-DC[2]*DDFS[ab00+3];
        DDFP[ab01+4*3+0]=DDGS[ab10+13]-DC[0]*DDFS[ab00+4];
        DDFP[ab01+4*3+1]=DDGS[ab10+ 4]-DC[1]*DDFS[ab00+4];
        DDFP[ab01+4*3+2]=DDGS[ab10+ 7]-DC[2]*DDFS[ab00+4];
        DDFP[ab01+5*3+0]=DDGS[ab10+ 8]-DC[0]*DDFS[ab00+5];
        DDFP[ab01+5*3+1]=DDGS[ab10+14]-DC[1]*DDFS[ab00+5];
        DDFP[ab01+5*3+2]=DDGS[ab10+ 5]-DC[2]*DDFS[ab00+5];
        DDFP[ab01+6*3+0]=DDGS[ab10+ 6]-DC[0]*DDFS[ab00+6];
        DDFP[ab01+6*3+1]=DDGS[ab10+ 9]-DC[1]*DDFS[ab00+6];
        DDFP[ab01+6*3+2]=DDGS[ab10+13]-DC[2]*DDFS[ab00+6];
        DDFP[ab01+7*3+0]=DDGS[ab10+14]-DC[0]*DDFS[ab00+7];
        DDFP[ab01+7*3+1]=DDGS[ab10+ 7]-DC[1]*DDFS[ab00+7];
        DDFP[ab01+7*3+2]=DDGS[ab10+10]-DC[2]*DDFS[ab00+7];
        DDFP[ab01+8*3+0]=DDGS[ab10+11]-DC[0]*DDFS[ab00+8];
        DDFP[ab01+8*3+1]=DDGS[ab10+12]-DC[1]*DDFS[ab00+8];
        DDFP[ab01+8*3+2]=DDGS[ab10+ 8]-DC[2]*DDFS[ab00+8];
        DDFP[ab01+9*3+0]=DDGS[ab10+12]-DC[0]*DDFS[ab00+9];
        DDFP[ab01+9*3+1]=DDGS[ab10+13]-DC[1]*DDFS[ab00+9];
        DDFP[ab01+9*3+2]=DDGS[ab10+14]-DC[2]*DDFS[ab00+9];
    }
    // (D,D|D,D)
    for (ab=0, ab01=0, ab10=0, ab00=0;
            ab<6*6; ab++, ab01 += 6*6, ab10 += 10*3, ab00 += 6*3) {
        DDDD[ab01+0*6+0]=DDFP[ab10+0*3+0]-DC[0]*DDDP[ab00+0*3+0];
        DDDD[ab01+0*6+1]=DDFP[ab10+3*3+1]-DC[1]*DDDP[ab00+0*3+1];
        DDDD[ab01+0*6+2]=DDFP[ab10+8*3+2]-DC[2]*DDDP[ab00+0*3+2];
        DDDD[ab01+0*6+3]=DDFP[ab10+0*3+1]-DC[0]*DDDP[ab00+0*3+1];
        DDDD[ab01+0*6+4]=DDFP[ab10+3*3+2]-DC[1]*DDDP[ab00+0*3+2];
        DDDD[ab01+0*6+5]=DDFP[ab10+8*3+0]-DC[2]*DDDP[ab00+0*3+0];
        DDDD[ab01+1*6+0]=DDFP[ab10+6*3+0]-DC[0]*DDDP[ab00+1*3+0];
        DDDD[ab01+1*6+1]=DDFP[ab10+1*3+1]-DC[1]*DDDP[ab00+1*3+1];
        DDDD[ab01+1*6+2]=DDFP[ab10+4*3+2]-DC[2]*DDDP[ab00+1*3+2];
        DDDD[ab01+1*6+3]=DDFP[ab10+6*3+1]-DC[0]*DDDP[ab00+1*3+1];
        DDDD[ab01+1*6+4]=DDFP[ab10+1*3+2]-DC[1]*DDDP[ab00+1*3+2];
        DDDD[ab01+1*6+5]=DDFP[ab10+4*3+0]-DC[2]*DDDP[ab00+1*3+0];
        DDDD[ab01+2*6+0]=DDFP[ab10+5*3+0]-DC[0]*DDDP[ab00+2*3+0];
        DDDD[ab01+2*6+1]=DDFP[ab10+7*3+1]-DC[1]*DDDP[ab00+2*3+1];
        DDDD[ab01+2*6+2]=DDFP[ab10+2*3+2]-DC[2]*DDDP[ab00+2*3+2];
        DDDD[ab01+2*6+3]=DDFP[ab10+5*3+1]-DC[0]*DDDP[ab00+2*3+1];
        DDDD[ab01+2*6+4]=DDFP[ab10+7*3+2]-DC[1]*DDDP[ab00+2*3+2];
        DDDD[ab01+2*6+5]=DDFP[ab10+2*3+0]-DC[2]*DDDP[ab00+2*3+0];
        DDDD[ab01+3*6+0]=DDFP[ab10+3*3+0]-DC[0]*DDDP[ab00+3*3+0];
        DDDD[ab01+3*6+1]=DDFP[ab10+6*3+1]-DC[1]*DDDP[ab00+3*3+1];
        DDDD[ab01+3*6+2]=DDFP[ab10+9*3+2]-DC[2]*DDDP[ab00+3*3+2];
        DDDD[ab01+3*6+3]=DDFP[ab10+3*3+1]-DC[0]*DDDP[ab00+3*3+1];
        DDDD[ab01+3*6+4]=DDFP[ab10+6*3+2]-DC[1]*DDDP[ab00+3*3+2];
        DDDD[ab01+3*6+5]=DDFP[ab10+9*3+0]-DC[2]*DDDP[ab00+3*3+0];
        DDDD[ab01+4*6+0]=DDFP[ab10+9*3+0]-DC[0]*DDDP[ab00+4*3+0];
        DDDD[ab01+4*6+1]=DDFP[ab10+4*3+1]-DC[1]*DDDP[ab00+4*3+1];
        DDDD[ab01+4*6+2]=DDFP[ab10+7*3+2]-DC[2]*DDDP[ab00+4*3+2];
        DDDD[ab01+4*6+3]=DDFP[ab10+9*3+1]-DC[0]*DDDP[ab00+4*3+1];
        DDDD[ab01+4*6+4]=DDFP[ab10+4*3+2]-DC[1]*DDDP[ab00+4*3+2];
        DDDD[ab01+4*6+5]=DDFP[ab10+7*3+0]-DC[2]*DDDP[ab00+4*3+0];
        DDDD[ab01+5*6+0]=DDFP[ab10+8*3+0]-DC[0]*DDDP[ab00+5*3+0];
        DDDD[ab01+5*6+1]=DDFP[ab10+9*3+1]-DC[1]*DDDP[ab00+5*3+1];
        DDDD[ab01+5*6+2]=DDFP[ab10+5*3+2]-DC[2]*DDDP[ab00+5*3+2];
        DDDD[ab01+5*6+3]=DDFP[ab10+8*3+1]-DC[0]*DDDP[ab00+5*3+1];
        DDDD[ab01+5*6+4]=DDFP[ab10+9*3+2]-DC[1]*DDDP[ab00+5*3+2];
        DDDD[ab01+5*6+5]=DDFP[ab10+5*3+0]-DC[2]*DDDP[ab00+5*3+0];
    }
    //
    for ( i=0, ix=0; i<6; i++ ) {
	coe0 = (i<3? ONE : sqr3 );
	for ( j=0; j<6; j++ ) {
	    coe1 = coe0 * (j<3? ONE : sqr3 );
	    for ( k=0; k<6; k++ ) {
		coe2 = coe1 * (k<3? ONE : sqr3 );
		for ( l=0; l<6; l++, ix++ ) {
		    coe = coe2 * (l<3? ONE : sqr3 );
		    DDDD[ix] *= coe;
		}
	    }
	}
    }
}
