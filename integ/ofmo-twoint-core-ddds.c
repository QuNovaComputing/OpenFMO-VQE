/**
 * @file ofmo-twoint-core-ddds.c
 * １つのCS４重対に対する２電子積分を計算する関数群。
 * 2011/06/16現在、(ss,ss)～(dd,dd)までの２１種類の２電子積分
 * 計算をサポートしている。
 * このファイルには、(dd,ds)タイプの２電子積分計算を
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

extern double *FMT_fmt_table6;

extern double FMT_fmt_step_size;
extern double FMT_fmt_inv_step_size;
extern double FMT_pi_div2;

/** １つのCS４重対に対して(dd,ds)タイプの２電子積分を計算する関数
 * @ingroup core-twoint
 * */
void twoint_core_ddds__(
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double DDDS[6*6*6] ) {
    int ijps, klps, i, j, k, ix, m, m1, c0;
    double cssss, zeta, dkab, xiza, eta, xizc, dk;
    double rz, re, ze2, ze22, ze23, ze24, zeta2, eta2, zeta23;
    union _temp_ {
	double entity[15];
	double ssss;
	double psss[3];
	double dsss[6];
	double fsss[10];
	double gsss[15];
    } tmp;
    double sqrho, rho, PQ2, T;
    double PC[3], PA[3], WP[3], QP[3], QC[3], WQ[3];
    double ssss[6+1], psss[5+1][3], dsss[4+1][6], fsss[3+1][10];
    double gsss[2+1][15];
    double psps[3*3], dsps[1+1][6*3], fsps[1+1][10*3], gsps[1+1][15*3];
    double DSDS[6*6], FSDS[10*6], GSDS[15*6];
    double DPDS[6*3*6], FPDS[10*3*6];
    double sqr3, coe0, coe1, coe;
    sqr3 = sqrt(3.e0);
    for ( i=0; i< 6*6; i++ ) DSDS[i] = ZERO;
    for ( i=0; i<10*6; i++ ) FSDS[i] = ZERO;
    for ( i=0; i<15*6; i++ ) GSDS[i] = ZERO;
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
	    }
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
	    }
	    // gsss (m=0,2)
	    for (m=0, m1=1; m<=2; m++, m1++) {
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
	    // psps (m=1)
	    psps[0*3+0]=QC[0]*psss[1][0]+WQ[0]*psss[2][0]+ze2*ssss[2];
	    psps[0*3+1]=QC[1]*psss[1][0]+WQ[1]*psss[2][0];
	    psps[0*3+2]=QC[2]*psss[1][0]+WQ[2]*psss[2][0];
	    psps[1*3+0]=QC[0]*psss[1][1]+WQ[0]*psss[2][1];
	    psps[1*3+1]=QC[1]*psss[1][1]+WQ[1]*psss[2][1]+ze2*ssss[2];
	    psps[1*3+2]=QC[2]*psss[1][1]+WQ[2]*psss[2][1];
	    psps[2*3+0]=QC[0]*psss[1][2]+WQ[0]*psss[2][2];
	    psps[2*3+1]=QC[1]*psss[1][2]+WQ[1]*psss[2][2];
	    psps[2*3+2]=QC[2]*psss[1][2]+WQ[2]*psss[2][2]+ze2*ssss[2];
	    // dsps (m=0,1)
	    for (m=0, m1=1; m<=1; m++, m1++) {
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
	    for (m=0, m1=1; m<=1; m++, m1++) {
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
	    // gsps (m=0,1)
	    for (m=0, m1=1; m<=1; m++, m1++) {
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
	    }
	    // dsds (m=0)
	    for (i=0; i<6; i++) tmp.dsss[i]=dsss[0][i]-re*dsss[1][i];
	    DSDS[0*6+0]+=QC[0]*dsps[0][0*3+0] + WQ[0]*dsps[1][0*3+0]
			+ eta2*tmp.dsss[0] + ze22*psps[0*3+0];
	    DSDS[0*6+1]+=QC[1]*dsps[0][0*3+1] + WQ[1]*dsps[1][0*3+1]
			+ eta2*tmp.dsss[0];
	    DSDS[0*6+2]+=QC[2]*dsps[0][0*3+2] + WQ[2]*dsps[1][0*3+2]
			+ eta2*tmp.dsss[0];
	    DSDS[0*6+3]+=QC[0]*dsps[0][0*3+1] + WQ[0]*dsps[1][0*3+1]
			+ ze22*psps[0*3+1];
	    DSDS[0*6+4]+=QC[1]*dsps[0][0*3+2] + WQ[1]*dsps[1][0*3+2];
	    DSDS[0*6+5]+=QC[2]*dsps[0][0*3+0] + WQ[2]*dsps[1][0*3+0];
	    DSDS[1*6+0]+=QC[0]*dsps[0][1*3+0] + WQ[0]*dsps[1][1*3+0]
			+ eta2*tmp.dsss[1];
	    DSDS[1*6+1]+=QC[1]*dsps[0][1*3+1] + WQ[1]*dsps[1][1*3+1]
			+ eta2*tmp.dsss[1] + ze22*psps[1*3+1];
	    DSDS[1*6+2]+=QC[2]*dsps[0][1*3+2] + WQ[2]*dsps[1][1*3+2]
			+ eta2*tmp.dsss[1];
	    DSDS[1*6+3]+=QC[0]*dsps[0][1*3+1] + WQ[0]*dsps[1][1*3+1];
	    DSDS[1*6+4]+=QC[1]*dsps[0][1*3+2] + WQ[1]*dsps[1][1*3+2]
			+ ze22*psps[1*3+2];
	    DSDS[1*6+5]+=QC[2]*dsps[0][1*3+0] + WQ[2]*dsps[1][1*3+0];
	    DSDS[2*6+0]+=QC[0]*dsps[0][2*3+0] + WQ[0]*dsps[1][2*3+0]
			+ eta2*tmp.dsss[2];
	    DSDS[2*6+1]+=QC[1]*dsps[0][2*3+1] + WQ[1]*dsps[1][2*3+1]
			+ eta2*tmp.dsss[2];
	    DSDS[2*6+2]+=QC[2]*dsps[0][2*3+2] + WQ[2]*dsps[1][2*3+2]
			+ eta2*tmp.dsss[2] + ze22*psps[2*3+2];
	    DSDS[2*6+3]+=QC[0]*dsps[0][2*3+1] + WQ[0]*dsps[1][2*3+1];
	    DSDS[2*6+4]+=QC[1]*dsps[0][2*3+2] + WQ[1]*dsps[1][2*3+2];
	    DSDS[2*6+5]+=QC[2]*dsps[0][2*3+0] + WQ[2]*dsps[1][2*3+0]
			+ ze22*psps[2*3+0];
	    DSDS[3*6+0]+=QC[0]*dsps[0][3*3+0] + WQ[0]*dsps[1][3*3+0]
			+ eta2*tmp.dsss[3] + ze2*psps[1*3+0];
	    DSDS[3*6+1]+=QC[1]*dsps[0][3*3+1] + WQ[1]*dsps[1][3*3+1]
			+ eta2*tmp.dsss[3] + ze2*psps[0*3+1];
	    DSDS[3*6+2]+=QC[2]*dsps[0][3*3+2] + WQ[2]*dsps[1][3*3+2]
			+ eta2*tmp.dsss[3];
	    DSDS[3*6+3]+=QC[0]*dsps[0][3*3+1] + WQ[0]*dsps[1][3*3+1]
			+ ze2*psps[1*3+1];
	    DSDS[3*6+4]+=QC[1]*dsps[0][3*3+2] + WQ[1]*dsps[1][3*3+2]
			+ ze2*psps[0*3+2];
	    DSDS[3*6+5]+=QC[2]*dsps[0][3*3+0] + WQ[2]*dsps[1][3*3+0];
	    DSDS[4*6+0]+=QC[0]*dsps[0][4*3+0] + WQ[0]*dsps[1][4*3+0]
			+ eta2*tmp.dsss[4];
	    DSDS[4*6+1]+=QC[1]*dsps[0][4*3+1] + WQ[1]*dsps[1][4*3+1]
			+ eta2*tmp.dsss[4] + ze2*psps[2*3+1];
	    DSDS[4*6+2]+=QC[2]*dsps[0][4*3+2] + WQ[2]*dsps[1][4*3+2]
			+ eta2*tmp.dsss[4] + ze2*psps[1*3+2];
	    DSDS[4*6+3]+=QC[0]*dsps[0][4*3+1] + WQ[0]*dsps[1][4*3+1];
	    DSDS[4*6+4]+=QC[1]*dsps[0][4*3+2] + WQ[1]*dsps[1][4*3+2]
			+ ze2*psps[2*3+2];
	    DSDS[4*6+5]+=QC[2]*dsps[0][4*3+0] + WQ[2]*dsps[1][4*3+0]
			+ ze2*psps[1*3+0];
	    DSDS[5*6+0]+=QC[0]*dsps[0][5*3+0] + WQ[0]*dsps[1][5*3+0]
			+ eta2*tmp.dsss[5] + ze2*psps[2*3+0];
	    DSDS[5*6+1]+=QC[1]*dsps[0][5*3+1] + WQ[1]*dsps[1][5*3+1]
			+ eta2*tmp.dsss[5];
	    DSDS[5*6+2]+=QC[2]*dsps[0][5*3+2] + WQ[2]*dsps[1][5*3+2]
			+ eta2*tmp.dsss[5] + ze2*psps[0*3+2];
	    DSDS[5*6+3]+=QC[0]*dsps[0][5*3+1] + WQ[0]*dsps[1][5*3+1]
			+ ze2*psps[2*3+1];
	    DSDS[5*6+4]+=QC[1]*dsps[0][5*3+2] + WQ[1]*dsps[1][5*3+2];
	    DSDS[5*6+5]+=QC[2]*dsps[0][5*3+0] + WQ[2]*dsps[1][5*3+0]
			+ ze2*psps[0*3+0];
	    // fsds (m=0)
	    for (i=0; i<10; i++) tmp.fsss[i]=fsss[0][i]-re*fsss[1][i];
	    FSDS[0*6+0]+=QC[0]*fsps[0][0*3+0] + WQ[0]*fsps[1][0*3+0]
			+ eta2*tmp.fsss[0] + ze23*dsps[1][0*3+0];
	    FSDS[0*6+1]+=QC[1]*fsps[0][0*3+1] + WQ[1]*fsps[1][0*3+1]
			+ eta2*tmp.fsss[0];
	    FSDS[0*6+2]+=QC[2]*fsps[0][0*3+2] + WQ[2]*fsps[1][0*3+2]
			+ eta2*tmp.fsss[0];
	    FSDS[0*6+3]+=QC[0]*fsps[0][0*3+1] + WQ[0]*fsps[1][0*3+1]
			+ ze23*dsps[1][0*3+1];
	    FSDS[0*6+4]+=QC[1]*fsps[0][0*3+2] + WQ[1]*fsps[1][0*3+2];
	    FSDS[0*6+5]+=QC[2]*fsps[0][0*3+0] + WQ[2]*fsps[1][0*3+0];
	    FSDS[1*6+0]+=QC[0]*fsps[0][1*3+0] + WQ[0]*fsps[1][1*3+0]
			+ eta2*tmp.fsss[1];
	    FSDS[1*6+1]+=QC[1]*fsps[0][1*3+1] + WQ[1]*fsps[1][1*3+1]
			+ eta2*tmp.fsss[1] + ze23*dsps[1][1*3+1];
	    FSDS[1*6+2]+=QC[2]*fsps[0][1*3+2] + WQ[2]*fsps[1][1*3+2]
			+ eta2*tmp.fsss[1];
	    FSDS[1*6+3]+=QC[0]*fsps[0][1*3+1] + WQ[0]*fsps[1][1*3+1];
	    FSDS[1*6+4]+=QC[1]*fsps[0][1*3+2] + WQ[1]*fsps[1][1*3+2]
			+ ze23*dsps[1][1*3+2];
	    FSDS[1*6+5]+=QC[2]*fsps[0][1*3+0] + WQ[2]*fsps[1][1*3+0];
	    FSDS[2*6+0]+=QC[0]*fsps[0][2*3+0] + WQ[0]*fsps[1][2*3+0]
			+ eta2*tmp.fsss[2];
	    FSDS[2*6+1]+=QC[1]*fsps[0][2*3+1] + WQ[1]*fsps[1][2*3+1]
			+ eta2*tmp.fsss[2];
	    FSDS[2*6+2]+=QC[2]*fsps[0][2*3+2] + WQ[2]*fsps[1][2*3+2]
			+ eta2*tmp.fsss[2] + ze23*dsps[1][2*3+2];
	    FSDS[2*6+3]+=QC[0]*fsps[0][2*3+1] + WQ[0]*fsps[1][2*3+1];
	    FSDS[2*6+4]+=QC[1]*fsps[0][2*3+2] + WQ[1]*fsps[1][2*3+2];
	    FSDS[2*6+5]+=QC[2]*fsps[0][2*3+0] + WQ[2]*fsps[1][2*3+0]
			+ ze23*dsps[1][2*3+0];
	    FSDS[3*6+0]+=QC[0]*fsps[0][3*3+0] + WQ[0]*fsps[1][3*3+0]
			+ eta2*tmp.fsss[3] + ze22*dsps[1][3*3+0];
	    FSDS[3*6+1]+=QC[1]*fsps[0][3*3+1] + WQ[1]*fsps[1][3*3+1]
			+ eta2*tmp.fsss[3] + ze2*dsps[1][0*3+1];
	    FSDS[3*6+2]+=QC[2]*fsps[0][3*3+2] + WQ[2]*fsps[1][3*3+2]
			+ eta2*tmp.fsss[3];
	    FSDS[3*6+3]+=QC[0]*fsps[0][3*3+1] + WQ[0]*fsps[1][3*3+1]
			+ ze22*dsps[1][3*3+1];
	    FSDS[3*6+4]+=QC[1]*fsps[0][3*3+2] + WQ[1]*fsps[1][3*3+2]
			+ ze2*dsps[1][0*3+2];
	    FSDS[3*6+5]+=QC[2]*fsps[0][3*3+0] + WQ[2]*fsps[1][3*3+0];
	    FSDS[4*6+0]+=QC[0]*fsps[0][4*3+0] + WQ[0]*fsps[1][4*3+0]
			+ eta2*tmp.fsss[4];
	    FSDS[4*6+1]+=QC[1]*fsps[0][4*3+1] + WQ[1]*fsps[1][4*3+1]
			+ eta2*tmp.fsss[4] + ze22*dsps[1][4*3+1];
	    FSDS[4*6+2]+=QC[2]*fsps[0][4*3+2] + WQ[2]*fsps[1][4*3+2]
			+ eta2*tmp.fsss[4] + ze2*dsps[1][1*3+2];
	    FSDS[4*6+3]+=QC[0]*fsps[0][4*3+1] + WQ[0]*fsps[1][4*3+1];
	    FSDS[4*6+4]+=QC[1]*fsps[0][4*3+2] + WQ[1]*fsps[1][4*3+2]
			+ ze22*dsps[1][4*3+2];
	    FSDS[4*6+5]+=QC[2]*fsps[0][4*3+0] + WQ[2]*fsps[1][4*3+0]
			+ ze2*dsps[1][1*3+0];
	    FSDS[5*6+0]+=QC[0]*fsps[0][5*3+0] + WQ[0]*fsps[1][5*3+0]
			+ eta2*tmp.fsss[5] + ze2*dsps[1][2*3+0];
	    FSDS[5*6+1]+=QC[1]*fsps[0][5*3+1] + WQ[1]*fsps[1][5*3+1]
			+ eta2*tmp.fsss[5];
	    FSDS[5*6+2]+=QC[2]*fsps[0][5*3+2] + WQ[2]*fsps[1][5*3+2]
			+ eta2*tmp.fsss[5] + ze22*dsps[1][5*3+2];
	    FSDS[5*6+3]+=QC[0]*fsps[0][5*3+1] + WQ[0]*fsps[1][5*3+1]
			+ ze2*dsps[1][2*3+1];
	    FSDS[5*6+4]+=QC[1]*fsps[0][5*3+2] + WQ[1]*fsps[1][5*3+2];
	    FSDS[5*6+5]+=QC[2]*fsps[0][5*3+0] + WQ[2]*fsps[1][5*3+0]
			+ ze22*dsps[1][5*3+0];
	    FSDS[6*6+0]+=QC[0]*fsps[0][6*3+0] + WQ[0]*fsps[1][6*3+0]
			+ eta2*tmp.fsss[6] + ze2*dsps[1][1*3+0];
	    FSDS[6*6+1]+=QC[1]*fsps[0][6*3+1] + WQ[1]*fsps[1][6*3+1]
			+ eta2*tmp.fsss[6] + ze22*dsps[1][3*3+1];
	    FSDS[6*6+2]+=QC[2]*fsps[0][6*3+2] + WQ[2]*fsps[1][6*3+2]
			+ eta2*tmp.fsss[6];
	    FSDS[6*6+3]+=QC[0]*fsps[0][6*3+1] + WQ[0]*fsps[1][6*3+1]
			+ ze2*dsps[1][1*3+1];
	    FSDS[6*6+4]+=QC[1]*fsps[0][6*3+2] + WQ[1]*fsps[1][6*3+2]
			+ ze22*dsps[1][3*3+2];
	    FSDS[6*6+5]+=QC[2]*fsps[0][6*3+0] + WQ[2]*fsps[1][6*3+0];
	    FSDS[7*6+0]+=QC[0]*fsps[0][7*3+0] + WQ[0]*fsps[1][7*3+0]
			+ eta2*tmp.fsss[7];
	    FSDS[7*6+1]+=QC[1]*fsps[0][7*3+1] + WQ[1]*fsps[1][7*3+1]
			+ eta2*tmp.fsss[7] + ze2*dsps[1][2*3+1];
	    FSDS[7*6+2]+=QC[2]*fsps[0][7*3+2] + WQ[2]*fsps[1][7*3+2]
			+ eta2*tmp.fsss[7] + ze22*dsps[1][4*3+2];
	    FSDS[7*6+3]+=QC[0]*fsps[0][7*3+1] + WQ[0]*fsps[1][7*3+1];
	    FSDS[7*6+4]+=QC[1]*fsps[0][7*3+2] + WQ[1]*fsps[1][7*3+2]
			+ ze2*dsps[1][2*3+2];
	    FSDS[7*6+5]+=QC[2]*fsps[0][7*3+0] + WQ[2]*fsps[1][7*3+0]
			+ ze22*dsps[1][4*3+0];
	    FSDS[8*6+0]+=QC[0]*fsps[0][8*3+0] + WQ[0]*fsps[1][8*3+0]
			+ eta2*tmp.fsss[8] + ze22*dsps[1][5*3+0];
	    FSDS[8*6+1]+=QC[1]*fsps[0][8*3+1] + WQ[1]*fsps[1][8*3+1]
			+ eta2*tmp.fsss[8];
	    FSDS[8*6+2]+=QC[2]*fsps[0][8*3+2] + WQ[2]*fsps[1][8*3+2]
			+ eta2*tmp.fsss[8] + ze2*dsps[1][0*3+2];
	    FSDS[8*6+3]+=QC[0]*fsps[0][8*3+1] + WQ[0]*fsps[1][8*3+1]
			+ ze22*dsps[1][5*3+1];
	    FSDS[8*6+4]+=QC[1]*fsps[0][8*3+2] + WQ[1]*fsps[1][8*3+2];
	    FSDS[8*6+5]+=QC[2]*fsps[0][8*3+0] + WQ[2]*fsps[1][8*3+0]
			+ ze2*dsps[1][0*3+0];
	    FSDS[9*6+0]+=QC[0]*fsps[0][9*3+0] + WQ[0]*fsps[1][9*3+0]
			+ eta2*tmp.fsss[9] + ze2*dsps[1][4*3+0];
	    FSDS[9*6+1]+=QC[1]*fsps[0][9*3+1] + WQ[1]*fsps[1][9*3+1]
			+ eta2*tmp.fsss[9] + ze2*dsps[1][5*3+1];
	    FSDS[9*6+2]+=QC[2]*fsps[0][9*3+2] + WQ[2]*fsps[1][9*3+2]
			+ eta2*tmp.fsss[9] + ze2*dsps[1][3*3+2];
	    FSDS[9*6+3]+=QC[0]*fsps[0][9*3+1] + WQ[0]*fsps[1][9*3+1]
			+ ze2*dsps[1][4*3+1];
	    FSDS[9*6+4]+=QC[1]*fsps[0][9*3+2] + WQ[1]*fsps[1][9*3+2]
			+ ze2*dsps[1][5*3+2];
	    FSDS[9*6+5]+=QC[2]*fsps[0][9*3+0] + WQ[2]*fsps[1][9*3+0]
			+ ze2*dsps[1][3*3+0];
	    // gsds (m=0)
	    for (i=0; i<15; i++) tmp.gsss[i]=gsss[0][i]-re*gsss[1][i];
	    GSDS[ 0*6+0]+=QC[0]*gsps[0][ 0*3+0] + WQ[0]*gsps[1][ 0*3+0]
			 + eta2*tmp.gsss[ 0] + ze24*fsps[1][0*3+0];
	    GSDS[ 0*6+1]+=QC[1]*gsps[0][ 0*3+1] + WQ[1]*gsps[1][ 0*3+1]
			 + eta2*tmp.gsss[ 0];
	    GSDS[ 0*6+2]+=QC[2]*gsps[0][ 0*3+2] + WQ[2]*gsps[1][ 0*3+2]
			 + eta2*tmp.gsss[ 0];
	    GSDS[ 0*6+3]+=QC[0]*gsps[0][ 0*3+1] + WQ[0]*gsps[1][ 0*3+1]
			 + ze24*fsps[1][0*3+1];
	    GSDS[ 0*6+4]+=QC[1]*gsps[0][ 0*3+2] + WQ[1]*gsps[1][ 0*3+2];
	    GSDS[ 0*6+5]+=QC[2]*gsps[0][ 0*3+0] + WQ[2]*gsps[1][ 0*3+0];
	    GSDS[ 1*6+0]+=QC[0]*gsps[0][ 1*3+0] + WQ[0]*gsps[1][ 1*3+0]
			 + eta2*tmp.gsss[ 1];
	    GSDS[ 1*6+1]+=QC[1]*gsps[0][ 1*3+1] + WQ[1]*gsps[1][ 1*3+1]
			 + eta2*tmp.gsss[ 1] + ze24*fsps[1][1*3+1];
	    GSDS[ 1*6+2]+=QC[2]*gsps[0][ 1*3+2] + WQ[2]*gsps[1][ 1*3+2]
			 + eta2*tmp.gsss[ 1];
	    GSDS[ 1*6+3]+=QC[0]*gsps[0][ 1*3+1] + WQ[0]*gsps[1][ 1*3+1];
	    GSDS[ 1*6+4]+=QC[1]*gsps[0][ 1*3+2] + WQ[1]*gsps[1][ 1*3+2]
			 + ze24*fsps[1][1*3+2];
	    GSDS[ 1*6+5]+=QC[2]*gsps[0][ 1*3+0] + WQ[2]*gsps[1][ 1*3+0];
	    GSDS[ 2*6+0]+=QC[0]*gsps[0][ 2*3+0] + WQ[0]*gsps[1][ 2*3+0]
			 + eta2*tmp.gsss[ 2];
	    GSDS[ 2*6+1]+=QC[1]*gsps[0][ 2*3+1] + WQ[1]*gsps[1][ 2*3+1]
			 + eta2*tmp.gsss[ 2];
	    GSDS[ 2*6+2]+=QC[2]*gsps[0][ 2*3+2] + WQ[2]*gsps[1][ 2*3+2]
			 + eta2*tmp.gsss[ 2] + ze24*fsps[1][2*3+2];
	    GSDS[ 2*6+3]+=QC[0]*gsps[0][ 2*3+1] + WQ[0]*gsps[1][ 2*3+1];
	    GSDS[ 2*6+4]+=QC[1]*gsps[0][ 2*3+2] + WQ[1]*gsps[1][ 2*3+2];
	    GSDS[ 2*6+5]+=QC[2]*gsps[0][ 2*3+0] + WQ[2]*gsps[1][ 2*3+0]
			 + ze24*fsps[1][2*3+0];
	    GSDS[ 3*6+0]+=QC[0]*gsps[0][ 3*3+0] + WQ[0]*gsps[1][ 3*3+0]
			 + eta2*tmp.gsss[ 3] + ze23*fsps[1][3*3+0];
	    GSDS[ 3*6+1]+=QC[1]*gsps[0][ 3*3+1] + WQ[1]*gsps[1][ 3*3+1]
			 + eta2*tmp.gsss[ 3] + ze2*fsps[1][0*3+1];
	    GSDS[ 3*6+2]+=QC[2]*gsps[0][ 3*3+2] + WQ[2]*gsps[1][ 3*3+2]
			 + eta2*tmp.gsss[ 3];
	    GSDS[ 3*6+3]+=QC[0]*gsps[0][ 3*3+1] + WQ[0]*gsps[1][ 3*3+1]
			 + ze23*fsps[1][3*3+1];
	    GSDS[ 3*6+4]+=QC[1]*gsps[0][ 3*3+2] + WQ[1]*gsps[1][ 3*3+2]
			 + ze2*fsps[1][0*3+2];
	    GSDS[ 3*6+5]+=QC[2]*gsps[0][ 3*3+0] + WQ[2]*gsps[1][ 3*3+0];
	    GSDS[ 4*6+0]+=QC[0]*gsps[0][ 4*3+0] + WQ[0]*gsps[1][ 4*3+0]
			 + eta2*tmp.gsss[ 4];
	    GSDS[ 4*6+1]+=QC[1]*gsps[0][ 4*3+1] + WQ[1]*gsps[1][ 4*3+1]
			 + eta2*tmp.gsss[ 4] + ze23*fsps[1][4*3+1];
	    GSDS[ 4*6+2]+=QC[2]*gsps[0][ 4*3+2] + WQ[2]*gsps[1][ 4*3+2]
			 + eta2*tmp.gsss[ 4] + ze2*fsps[1][1*3+2];
	    GSDS[ 4*6+3]+=QC[0]*gsps[0][ 4*3+1] + WQ[0]*gsps[1][ 4*3+1];
	    GSDS[ 4*6+4]+=QC[1]*gsps[0][ 4*3+2] + WQ[1]*gsps[1][ 4*3+2]
			 + ze23*fsps[1][4*3+2];
	    GSDS[ 4*6+5]+=QC[2]*gsps[0][ 4*3+0] + WQ[2]*gsps[1][ 4*3+0]
			 + ze2*fsps[1][1*3+0];
	    GSDS[ 5*6+0]+=QC[0]*gsps[0][ 5*3+0] + WQ[0]*gsps[1][ 5*3+0]
			 + eta2*tmp.gsss[ 5] + ze2*fsps[1][2*3+0];
	    GSDS[ 5*6+1]+=QC[1]*gsps[0][ 5*3+1] + WQ[1]*gsps[1][ 5*3+1]
			 + eta2*tmp.gsss[ 5];
	    GSDS[ 5*6+2]+=QC[2]*gsps[0][ 5*3+2] + WQ[2]*gsps[1][ 5*3+2]
			 + eta2*tmp.gsss[ 5] + ze23*fsps[1][5*3+2];
	    GSDS[ 5*6+3]+=QC[0]*gsps[0][ 5*3+1] + WQ[0]*gsps[1][ 5*3+1]
			 + ze2*fsps[1][2*3+1];
	    GSDS[ 5*6+4]+=QC[1]*gsps[0][ 5*3+2] + WQ[1]*gsps[1][ 5*3+2];
	    GSDS[ 5*6+5]+=QC[2]*gsps[0][ 5*3+0] + WQ[2]*gsps[1][ 5*3+0]
			 + ze23*fsps[1][5*3+0];
	    GSDS[ 6*6+0]+=QC[0]*gsps[0][ 6*3+0] + WQ[0]*gsps[1][ 6*3+0]
			 + eta2*tmp.gsss[ 6] + ze22*fsps[1][6*3+0];
	    GSDS[ 6*6+1]+=QC[1]*gsps[0][ 6*3+1] + WQ[1]*gsps[1][ 6*3+1]
			 + eta2*tmp.gsss[ 6] + ze22*fsps[1][3*3+1];
	    GSDS[ 6*6+2]+=QC[2]*gsps[0][ 6*3+2] + WQ[2]*gsps[1][ 6*3+2]
			 + eta2*tmp.gsss[ 6];
	    GSDS[ 6*6+3]+=QC[0]*gsps[0][ 6*3+1] + WQ[0]*gsps[1][ 6*3+1]
			 + ze22*fsps[1][6*3+1];
	    GSDS[ 6*6+4]+=QC[1]*gsps[0][ 6*3+2] + WQ[1]*gsps[1][ 6*3+2]
			 + ze22*fsps[1][3*3+2];
	    GSDS[ 6*6+5]+=QC[2]*gsps[0][ 6*3+0] + WQ[2]*gsps[1][ 6*3+0];
	    GSDS[ 7*6+0]+=QC[0]*gsps[0][ 7*3+0] + WQ[0]*gsps[1][ 7*3+0]
			 + eta2*tmp.gsss[ 7];
	    GSDS[ 7*6+1]+=QC[1]*gsps[0][ 7*3+1] + WQ[1]*gsps[1][ 7*3+1]
			 + eta2*tmp.gsss[ 7] + ze22*fsps[1][7*3+1];
	    GSDS[ 7*6+2]+=QC[2]*gsps[0][ 7*3+2] + WQ[2]*gsps[1][ 7*3+2]
			 + eta2*tmp.gsss[ 7] + ze22*fsps[1][4*3+2];
	    GSDS[ 7*6+3]+=QC[0]*gsps[0][ 7*3+1] + WQ[0]*gsps[1][ 7*3+1];
	    GSDS[ 7*6+4]+=QC[1]*gsps[0][ 7*3+2] + WQ[1]*gsps[1][ 7*3+2]
			 + ze22*fsps[1][7*3+2];
	    GSDS[ 7*6+5]+=QC[2]*gsps[0][ 7*3+0] + WQ[2]*gsps[1][ 7*3+0]
			 + ze22*fsps[1][4*3+0];
	    GSDS[ 8*6+0]+=QC[0]*gsps[0][ 8*3+0] + WQ[0]*gsps[1][ 8*3+0]
			 + eta2*tmp.gsss[ 8] + ze22*fsps[1][5*3+0];
	    GSDS[ 8*6+1]+=QC[1]*gsps[0][ 8*3+1] + WQ[1]*gsps[1][ 8*3+1]
			 + eta2*tmp.gsss[ 8];
	    GSDS[ 8*6+2]+=QC[2]*gsps[0][ 8*3+2] + WQ[2]*gsps[1][ 8*3+2]
			 + eta2*tmp.gsss[ 8] + ze22*fsps[1][8*3+2];
	    GSDS[ 8*6+3]+=QC[0]*gsps[0][ 8*3+1] + WQ[0]*gsps[1][ 8*3+1]
			 + ze22*fsps[1][5*3+1];
	    GSDS[ 8*6+4]+=QC[1]*gsps[0][ 8*3+2] + WQ[1]*gsps[1][ 8*3+2];
	    GSDS[ 8*6+5]+=QC[2]*gsps[0][ 8*3+0] + WQ[2]*gsps[1][ 8*3+0]
			 + ze22*fsps[1][8*3+0];
	    GSDS[ 9*6+0]+=QC[0]*gsps[0][ 9*3+0] + WQ[0]*gsps[1][ 9*3+0]
			 + eta2*tmp.gsss[ 9] + ze2*fsps[1][1*3+0];
	    GSDS[ 9*6+1]+=QC[1]*gsps[0][ 9*3+1] + WQ[1]*gsps[1][ 9*3+1]
			 + eta2*tmp.gsss[ 9] + ze23*fsps[1][6*3+1];
	    GSDS[ 9*6+2]+=QC[2]*gsps[0][ 9*3+2] + WQ[2]*gsps[1][ 9*3+2]
			 + eta2*tmp.gsss[ 9];
	    GSDS[ 9*6+3]+=QC[0]*gsps[0][ 9*3+1] + WQ[0]*gsps[1][ 9*3+1]
			 + ze2*fsps[1][1*3+1];
	    GSDS[ 9*6+4]+=QC[1]*gsps[0][ 9*3+2] + WQ[1]*gsps[1][ 9*3+2]
			 + ze23*fsps[1][6*3+2];
	    GSDS[ 9*6+5]+=QC[2]*gsps[0][ 9*3+0] + WQ[2]*gsps[1][ 9*3+0];
	    GSDS[10*6+0]+=QC[0]*gsps[0][10*3+0] + WQ[0]*gsps[1][10*3+0]
			 + eta2*tmp.gsss[10];
	    GSDS[10*6+1]+=QC[1]*gsps[0][10*3+1] + WQ[1]*gsps[1][10*3+1]
			 + eta2*tmp.gsss[10] + ze2*fsps[1][2*3+1];
	    GSDS[10*6+2]+=QC[2]*gsps[0][10*3+2] + WQ[2]*gsps[1][10*3+2]
			 + eta2*tmp.gsss[10] + ze23*fsps[1][7*3+2];
	    GSDS[10*6+3]+=QC[0]*gsps[0][10*3+1] + WQ[0]*gsps[1][10*3+1];
	    GSDS[10*6+4]+=QC[1]*gsps[0][10*3+2] + WQ[1]*gsps[1][10*3+2]
			 + ze2*fsps[1][2*3+2];
	    GSDS[10*6+5]+=QC[2]*gsps[0][10*3+0] + WQ[2]*gsps[1][10*3+0]
			 + ze23*fsps[1][7*3+0];
	    GSDS[11*6+0]+=QC[0]*gsps[0][11*3+0] + WQ[0]*gsps[1][11*3+0]
			 + eta2*tmp.gsss[11] + ze23*fsps[1][8*3+0];
	    GSDS[11*6+1]+=QC[1]*gsps[0][11*3+1] + WQ[1]*gsps[1][11*3+1]
			 + eta2*tmp.gsss[11];
	    GSDS[11*6+2]+=QC[2]*gsps[0][11*3+2] + WQ[2]*gsps[1][11*3+2]
			 + eta2*tmp.gsss[11] + ze2*fsps[1][0*3+2];
	    GSDS[11*6+3]+=QC[0]*gsps[0][11*3+1] + WQ[0]*gsps[1][11*3+1]
			 + ze23*fsps[1][8*3+1];
	    GSDS[11*6+4]+=QC[1]*gsps[0][11*3+2] + WQ[1]*gsps[1][11*3+2];
	    GSDS[11*6+5]+=QC[2]*gsps[0][11*3+0] + WQ[2]*gsps[1][11*3+0]
			 + ze2*fsps[1][0*3+0];
	    GSDS[12*6+0]+=QC[0]*gsps[0][12*3+0] + WQ[0]*gsps[1][12*3+0]
			 + eta2*tmp.gsss[12] + ze22*fsps[1][9*3+0];
	    GSDS[12*6+1]+=QC[1]*gsps[0][12*3+1] + WQ[1]*gsps[1][12*3+1]
			 + eta2*tmp.gsss[12] + ze2*fsps[1][8*3+1];
	    GSDS[12*6+2]+=QC[2]*gsps[0][12*3+2] + WQ[2]*gsps[1][12*3+2]
			 + eta2*tmp.gsss[12] + ze2*fsps[1][3*3+2];
	    GSDS[12*6+3]+=QC[0]*gsps[0][12*3+1] + WQ[0]*gsps[1][12*3+1]
			 + ze22*fsps[1][9*3+1];
	    GSDS[12*6+4]+=QC[1]*gsps[0][12*3+2] + WQ[1]*gsps[1][12*3+2]
			 + ze2*fsps[1][8*3+2];
	    GSDS[12*6+5]+=QC[2]*gsps[0][12*3+0] + WQ[2]*gsps[1][12*3+0]
			 + ze2*fsps[1][3*3+0];
	    GSDS[13*6+0]+=QC[0]*gsps[0][13*3+0] + WQ[0]*gsps[1][13*3+0]
			 + eta2*tmp.gsss[13] + ze2*fsps[1][4*3+0];
	    GSDS[13*6+1]+=QC[1]*gsps[0][13*3+1] + WQ[1]*gsps[1][13*3+1]
			 + eta2*tmp.gsss[13] + ze22*fsps[1][9*3+1];
	    GSDS[13*6+2]+=QC[2]*gsps[0][13*3+2] + WQ[2]*gsps[1][13*3+2]
			 + eta2*tmp.gsss[13] + ze2*fsps[1][6*3+2];
	    GSDS[13*6+3]+=QC[0]*gsps[0][13*3+1] + WQ[0]*gsps[1][13*3+1]
			 + ze2*fsps[1][4*3+1];
	    GSDS[13*6+4]+=QC[1]*gsps[0][13*3+2] + WQ[1]*gsps[1][13*3+2]
			 + ze22*fsps[1][9*3+2];
	    GSDS[13*6+5]+=QC[2]*gsps[0][13*3+0] + WQ[2]*gsps[1][13*3+0]
			 + ze2*fsps[1][6*3+0];
	    GSDS[14*6+0]+=QC[0]*gsps[0][14*3+0] + WQ[0]*gsps[1][14*3+0]
			 + eta2*tmp.gsss[14] + ze2*fsps[1][7*3+0];
	    GSDS[14*6+1]+=QC[1]*gsps[0][14*3+1] + WQ[1]*gsps[1][14*3+1]
			 + eta2*tmp.gsss[14] + ze2*fsps[1][5*3+1];
	    GSDS[14*6+2]+=QC[2]*gsps[0][14*3+2] + WQ[2]*gsps[1][14*3+2]
			 + eta2*tmp.gsss[14] + ze22*fsps[1][9*3+2];
	    GSDS[14*6+3]+=QC[0]*gsps[0][14*3+1] + WQ[0]*gsps[1][14*3+1]
			 + ze2*fsps[1][7*3+1];
	    GSDS[14*6+4]+=QC[1]*gsps[0][14*3+2] + WQ[1]*gsps[1][14*3+2]
			 + ze2*fsps[1][5*3+2];
	    GSDS[14*6+5]+=QC[2]*gsps[0][14*3+0] + WQ[2]*gsps[1][14*3+0]
			 + ze22*fsps[1][9*3+0];
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
    //
    for ( i=0, ix=0; i<6; i++ ) {
	coe0 = (i<3? ONE : sqr3 );
	for ( j=0; j<6; j++ ) {
	    coe1 = coe0 * (j<3? ONE : sqr3 );
	    for ( k=0; k<6; k++, ix++ ) {
		coe = coe1 * (k<3? ONE : sqr3 );
		DDDS[ix] *= coe;
	    }
	}
    }
}
