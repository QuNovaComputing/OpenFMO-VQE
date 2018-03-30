/**
 * @file ofmo-twoint-core-ddss-ddpp.c
 * １つのCS４重対に対する２電子積分を計算する関数群。
 * 2011/06/16現在、(ss,ss)～(dd,dd)までの２１種類の２電子積分
 * 計算をサポートしている。
 * このファイルには、(dd,ss)～(dd,pp)タイプの２電子積分計算を
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
/** １つのCS４重対に対して(dd,ss)タイプの２電子積分を計算する関数
 * @ingroup core-twoint
 * */
void twoint_core_ddss__(
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double DDSS[6*6] ) {
    int ijps, klps, i, j, ix, m, m1;
    double cssss, zeta, dkab, xiza, eta, xizc, dk;
    double rz, zeta2, zeta23;
    double sqrho, rho, PQ2, T;
    double PC[3], PA[3], WP[3], QP[3];
    double ssss[4+1], psss[3+1][3], dsss[2+1][6], fsss[1+1][10];
    double DSSS[6], FSSS[10], GSSS[15];
    double FPSS[10*3], DPSS[6*3];
    union _temp_ {
        double entity[6];
        double ssss;
        double psss[3];
        double dsss[6];
    } tmp;
    double sqr3, coe0, coe;
    sqr3 = sqrt(3.e0);
    for ( i=0; i<6; i++ )  DSSS[i] = ZERO;
    for ( i=0; i<10; i++ ) FSSS[i] = ZERO;
    for ( i=0; i<15; i++ ) GSSS[i] = ZERO;
    for ( ijps=0; ijps<(*nijps); ijps++ ) {
	zeta   = vzeta[ijps];
	dkab   = vdkab[ijps];
	xiza   = vxiza[ijps];
	zeta2  = HALF * zeta;
	zeta23 = zeta + zeta2;
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
	    // psss (m=0,3)
	    for ( m=0; m<=3; m++) {
		m1 = m+1;
		psss[m][0] = PA[0]*ssss[m] + WP[0]*ssss[m1];
		psss[m][1] = PA[1]*ssss[m] + WP[1]*ssss[m1];
		psss[m][2] = PA[2]*ssss[m] + WP[2]*ssss[m1];
	    }
	    // dsss (m=0,2)
	    for ( m=0; m<=2; m++) {
		m1 = m+1;
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
	    for ( i=0; i<6; i++ ) DSSS[i] += dsss[0][i];
	    // fsss (m=0,1)
	    for (m=0; m<=1; m++) {
		m1 = m+1;
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
	    for (i=0; i<10; i++) FSSS[i] += fsss[0][i];
	    // gsss (m=0)
	    for (i=0; i<6; i++) tmp.dsss[i]=dsss[0][i]-rz*dsss[1][i];
	    GSSS[ 0]+=PA[0]*fsss[0][0] + WP[0]*fsss[1][0]
		     + zeta23*tmp.dsss[0];
	    GSSS[ 1]+=PA[1]*fsss[0][1] + WP[1]*fsss[1][1]
		     + zeta23*tmp.dsss[1];
	    GSSS[ 2]+=PA[2]*fsss[0][2] + WP[2]*fsss[1][2]
		     + zeta23*tmp.dsss[2];
	    GSSS[ 3]+=PA[0]*fsss[0][3] + WP[0]*fsss[1][3]
		     + zeta*tmp.dsss[3];
	    GSSS[ 4]+=PA[1]*fsss[0][4] + WP[1]*fsss[1][4]
		     + zeta*tmp.dsss[4];
	    GSSS[ 5]+=PA[2]*fsss[0][5] + WP[2]*fsss[1][5]
		     + zeta*tmp.dsss[5];
	    GSSS[ 6]+=PA[0]*fsss[0][6] + WP[0]*fsss[1][6]
		     + zeta2*tmp.dsss[1];
	    GSSS[ 7]+=PA[1]*fsss[0][7] + WP[1]*fsss[1][7]
		     + zeta2*tmp.dsss[2];
	    GSSS[ 8]+=PA[2]*fsss[0][8] + WP[2]*fsss[1][8]
		     + zeta2*tmp.dsss[0];
	    GSSS[ 9]+=PA[0]*fsss[0][1] + WP[0]*fsss[1][1];
	    GSSS[10]+=PA[1]*fsss[0][2] + WP[1]*fsss[1][2];
	    GSSS[11]+=PA[2]*fsss[0][0] + WP[2]*fsss[1][0];
	    GSSS[12]+=PA[0]*fsss[0][9] + WP[0]*fsss[1][9]
		     + zeta2*tmp.dsss[4];
	    GSSS[13]+=PA[1]*fsss[0][9] + WP[1]*fsss[1][9]
		     + zeta2*tmp.dsss[5];
	    GSSS[14]+=PA[2]*fsss[0][9] + WP[2]*fsss[1][9]
		     + zeta2*tmp.dsss[3];
	}	// klps
    }	// ijps
    // (D,P|S,S)
    DPSS[0*3+0] = FSSS[0] - BA[0]*DSSS[0];
    DPSS[0*3+1] = FSSS[3] - BA[1]*DSSS[0];
    DPSS[0*3+2] = FSSS[8] - BA[2]*DSSS[0];
    DPSS[1*3+0] = FSSS[6] - BA[0]*DSSS[1];
    DPSS[1*3+1] = FSSS[1] - BA[1]*DSSS[1];
    DPSS[1*3+2] = FSSS[4] - BA[2]*DSSS[1];
    DPSS[2*3+0] = FSSS[5] - BA[0]*DSSS[2];
    DPSS[2*3+1] = FSSS[7] - BA[1]*DSSS[2];
    DPSS[2*3+2] = FSSS[2] - BA[2]*DSSS[2];
    DPSS[3*3+0] = FSSS[3] - BA[0]*DSSS[3];
    DPSS[3*3+1] = FSSS[6] - BA[1]*DSSS[3];
    DPSS[3*3+2] = FSSS[9] - BA[2]*DSSS[3];
    DPSS[4*3+0] = FSSS[9] - BA[0]*DSSS[4];
    DPSS[4*3+1] = FSSS[4] - BA[1]*DSSS[4];
    DPSS[4*3+2] = FSSS[7] - BA[2]*DSSS[4];
    DPSS[5*3+0] = FSSS[8] - BA[0]*DSSS[5];
    DPSS[5*3+1] = FSSS[9] - BA[1]*DSSS[5];
    DPSS[5*3+2] = FSSS[5] - BA[2]*DSSS[5];
    // (F,P|S,S)
    FPSS[0*3+0] = GSSS[ 0] - BA[0]*FSSS[0];
    FPSS[0*3+1] = GSSS[ 3] - BA[1]*FSSS[0];
    FPSS[0*3+2] = GSSS[11] - BA[2]*FSSS[0];
    FPSS[1*3+0] = GSSS[ 9] - BA[0]*FSSS[1];
    FPSS[1*3+1] = GSSS[ 1] - BA[1]*FSSS[1];
    FPSS[1*3+2] = GSSS[ 4] - BA[2]*FSSS[1];
    FPSS[2*3+0] = GSSS[ 5] - BA[0]*FSSS[2];
    FPSS[2*3+1] = GSSS[10] - BA[1]*FSSS[2];
    FPSS[2*3+2] = GSSS[ 2] - BA[2]*FSSS[2];
    FPSS[3*3+0] = GSSS[ 3] - BA[0]*FSSS[3];
    FPSS[3*3+1] = GSSS[ 6] - BA[1]*FSSS[3];
    FPSS[3*3+2] = GSSS[12] - BA[2]*FSSS[3];
    FPSS[4*3+0] = GSSS[13] - BA[0]*FSSS[4];
    FPSS[4*3+1] = GSSS[ 4] - BA[1]*FSSS[4];
    FPSS[4*3+2] = GSSS[ 7] - BA[2]*FSSS[4];
    FPSS[5*3+0] = GSSS[ 8] - BA[0]*FSSS[5];
    FPSS[5*3+1] = GSSS[14] - BA[1]*FSSS[5];
    FPSS[5*3+2] = GSSS[ 5] - BA[2]*FSSS[5];
    FPSS[6*3+0] = GSSS[ 6] - BA[0]*FSSS[6];
    FPSS[6*3+1] = GSSS[ 9] - BA[1]*FSSS[6];
    FPSS[6*3+2] = GSSS[13] - BA[2]*FSSS[6];
    FPSS[7*3+0] = GSSS[14] - BA[0]*FSSS[7];
    FPSS[7*3+1] = GSSS[ 7] - BA[1]*FSSS[7];
    FPSS[7*3+2] = GSSS[10] - BA[2]*FSSS[7];
    FPSS[8*3+0] = GSSS[11] - BA[0]*FSSS[8];
    FPSS[8*3+1] = GSSS[12] - BA[1]*FSSS[8];
    FPSS[8*3+2] = GSSS[ 8] - BA[2]*FSSS[8];
    FPSS[9*3+0] = GSSS[12] - BA[0]*FSSS[9];
    FPSS[9*3+1] = GSSS[13] - BA[1]*FSSS[9];
    FPSS[9*3+2] = GSSS[14] - BA[2]*FSSS[9];
    // (D,D|S,S)
    DDSS[0*6+0] = FPSS[0*3+0] - BA[0]*DPSS[0*3+0];
    DDSS[0*6+1] = FPSS[3*3+1] - BA[1]*DPSS[0*3+1];
    DDSS[0*6+2] = FPSS[8*3+2] - BA[2]*DPSS[0*3+2];
    DDSS[0*6+3] = FPSS[0*3+1] - BA[0]*DPSS[0*3+1];
    DDSS[0*6+4] = FPSS[3*3+2] - BA[1]*DPSS[0*3+2];
    DDSS[0*6+5] = FPSS[8*3+0] - BA[2]*DPSS[0*3+0];
    DDSS[1*6+0] = FPSS[6*3+0] - BA[0]*DPSS[1*3+0];
    DDSS[1*6+1] = FPSS[1*3+1] - BA[1]*DPSS[1*3+1];
    DDSS[1*6+2] = FPSS[4*3+2] - BA[2]*DPSS[1*3+2];
    DDSS[1*6+3] = FPSS[6*3+1] - BA[0]*DPSS[1*3+1];
    DDSS[1*6+4] = FPSS[1*3+2] - BA[1]*DPSS[1*3+2];
    DDSS[1*6+5] = FPSS[4*3+0] - BA[2]*DPSS[1*3+0];
    DDSS[2*6+0] = FPSS[5*3+0] - BA[0]*DPSS[2*3+0];
    DDSS[2*6+1] = FPSS[7*3+1] - BA[1]*DPSS[2*3+1];
    DDSS[2*6+2] = FPSS[2*3+2] - BA[2]*DPSS[2*3+2];
    DDSS[2*6+3] = FPSS[5*3+1] - BA[0]*DPSS[2*3+1];
    DDSS[2*6+4] = FPSS[7*3+2] - BA[1]*DPSS[2*3+2];
    DDSS[2*6+5] = FPSS[2*3+0] - BA[2]*DPSS[2*3+0];
    DDSS[3*6+0] = FPSS[3*3+0] - BA[0]*DPSS[3*3+0];
    DDSS[3*6+1] = FPSS[6*3+1] - BA[1]*DPSS[3*3+1];
    DDSS[3*6+2] = FPSS[9*3+2] - BA[2]*DPSS[3*3+2];
    DDSS[3*6+3] = FPSS[3*3+1] - BA[0]*DPSS[3*3+1];
    DDSS[3*6+4] = FPSS[6*3+2] - BA[1]*DPSS[3*3+2];
    DDSS[3*6+5] = FPSS[9*3+0] - BA[2]*DPSS[3*3+0];
    DDSS[4*6+0] = FPSS[9*3+0] - BA[0]*DPSS[4*3+0];
    DDSS[4*6+1] = FPSS[4*3+1] - BA[1]*DPSS[4*3+1];
    DDSS[4*6+2] = FPSS[7*3+2] - BA[2]*DPSS[4*3+2];
    DDSS[4*6+3] = FPSS[9*3+1] - BA[0]*DPSS[4*3+1];
    DDSS[4*6+4] = FPSS[4*3+2] - BA[1]*DPSS[4*3+2];
    DDSS[4*6+5] = FPSS[7*3+0] - BA[2]*DPSS[4*3+0];
    DDSS[5*6+0] = FPSS[8*3+0] - BA[0]*DPSS[5*3+0];
    DDSS[5*6+1] = FPSS[9*3+1] - BA[1]*DPSS[5*3+1];
    DDSS[5*6+2] = FPSS[5*3+2] - BA[2]*DPSS[5*3+2];
    DDSS[5*6+3] = FPSS[8*3+1] - BA[0]*DPSS[5*3+1];
    DDSS[5*6+4] = FPSS[9*3+2] - BA[1]*DPSS[5*3+2];
    DDSS[5*6+5] = FPSS[5*3+0] - BA[2]*DPSS[5*3+0];
    //
    for ( i=0, ix=0; i<6; i++ ) {
	coe0 = (i<3? ONE : sqr3 );
	for ( j=0; j<6; j++, ix++ ) {
	    coe = coe0 * (j<3? ONE : sqr3 );
	    DDSS[ix] *= coe;
	}
    }
}

/** １つのCS４重対に対して(dd,ps)タイプの２電子積分を計算する関数
 * @ingroup core-twoint
 * */
void twoint_core_ddps__(
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double DDPS[6*6*3] ) {
    int ijps, klps, i, j, k, ix, m, m1, c0;
    double cssss, zeta, dkab, xiza, eta, xizc, dk;
    double rz, re, ze2, ze22, ze23, ze24, zeta2, zeta23;
    union _temp_ {
	double entity[6];
	double ssss;
	double psss[3];
	double dsss[6];
    } tmp;
    double sqrho, rho, PQ2, T;
    double PC[3], PA[3], WP[3], QP[3], QC[3], WQ[3];
    double ssss[5+1], psss[4+1][3], dsss[3+1][6], fsss[2+1][10];
    double gsss[1+1][15];
    double DSPS[6*3], FSPS[10*3], GSPS[15*3];
    double DPPS[6*3*3], FPPS[10*3*3];
    double sqr3, coe0, coe;
    sqr3 = sqrt(3.e0);
    for ( i=0; i< 6*3; i++ ) DSPS[i] = ZERO;
    for ( i=0; i<10*3; i++ ) FSPS[i] = ZERO;
    for ( i=0; i<15*3; i++ ) GSPS[i] = ZERO;
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
	    for (m=0; m<=4; m++) {
		m1 = m+1;
		psss[m][0] = PA[0]*ssss[m] + WP[0]*ssss[m1];
		psss[m][1] = PA[1]*ssss[m] + WP[1]*ssss[m1];
		psss[m][2] = PA[2]*ssss[m] + WP[2]*ssss[m1];
	    }
	    // dsss (m=0,3)
	    for (m=0; m<=3; m++) {
		m1 = m+1;
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
	    // fsss (m=0,2)
	    for (m=0; m<=2; m++) {
		m1 = m+1;
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
	    // gsss (m=0,1)
	    for (m=0; m<=1; m++) {
		m1 = m+1;
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
	    // dsps (m=0)
	    DSPS[0*3+0]+=QC[0]*dsss[0][0]+WQ[0]*dsss[1][0]+ze22*psss[1][0];
	    DSPS[0*3+1]+=QC[1]*dsss[0][0]+WQ[1]*dsss[1][0];
	    DSPS[0*3+2]+=QC[2]*dsss[0][0]+WQ[2]*dsss[1][0];
	    DSPS[1*3+0]+=QC[0]*dsss[0][1]+WQ[0]*dsss[1][1];
	    DSPS[1*3+1]+=QC[1]*dsss[0][1]+WQ[1]*dsss[1][1]+ze22*psss[1][1];
	    DSPS[1*3+2]+=QC[2]*dsss[0][1]+WQ[2]*dsss[1][1];
	    DSPS[2*3+0]+=QC[0]*dsss[0][2]+WQ[0]*dsss[1][2];
	    DSPS[2*3+1]+=QC[1]*dsss[0][2]+WQ[1]*dsss[1][2];
	    DSPS[2*3+2]+=QC[2]*dsss[0][2]+WQ[2]*dsss[1][2]+ze22*psss[1][2];
	    DSPS[3*3+0]+=QC[0]*dsss[0][3]+WQ[0]*dsss[1][3]+ze2*psss[1][1];
	    DSPS[3*3+1]+=QC[1]*dsss[0][3]+WQ[1]*dsss[1][3]+ze2*psss[1][0];
	    DSPS[3*3+2]+=QC[2]*dsss[0][3]+WQ[2]*dsss[1][3];
	    DSPS[4*3+0]+=QC[0]*dsss[0][4]+WQ[0]*dsss[1][4];
	    DSPS[4*3+1]+=QC[1]*dsss[0][4]+WQ[1]*dsss[1][4]+ze2*psss[1][2];
	    DSPS[4*3+2]+=QC[2]*dsss[0][4]+WQ[2]*dsss[1][4]+ze2*psss[1][1];
	    DSPS[5*3+0]+=QC[0]*dsss[0][5]+WQ[0]*dsss[1][5]+ze2*psss[1][2];
	    DSPS[5*3+1]+=QC[1]*dsss[0][5]+WQ[1]*dsss[1][5];
	    DSPS[5*3+2]+=QC[2]*dsss[0][5]+WQ[2]*dsss[1][5]+ze2*psss[1][0];
	    // fsps (m=0)
	    FSPS[0*3+0]+=QC[0]*fsss[0][0]+WQ[0]*fsss[1][0]+ze23*dsss[1][0];
	    FSPS[0*3+1]+=QC[1]*fsss[0][0]+WQ[1]*fsss[1][0];
	    FSPS[0*3+2]+=QC[2]*fsss[0][0]+WQ[2]*fsss[1][0];
	    FSPS[1*3+0]+=QC[0]*fsss[0][1]+WQ[0]*fsss[1][1];
	    FSPS[1*3+1]+=QC[1]*fsss[0][1]+WQ[1]*fsss[1][1]+ze23*dsss[1][1];
	    FSPS[1*3+2]+=QC[2]*fsss[0][1]+WQ[2]*fsss[1][1];
	    FSPS[2*3+0]+=QC[0]*fsss[0][2]+WQ[0]*fsss[1][2];
	    FSPS[2*3+1]+=QC[1]*fsss[0][2]+WQ[1]*fsss[1][2];
	    FSPS[2*3+2]+=QC[2]*fsss[0][2]+WQ[2]*fsss[1][2]+ze23*dsss[1][2];
	    FSPS[3*3+0]+=QC[0]*fsss[0][3]+WQ[0]*fsss[1][3]+ze22*dsss[1][3];
	    FSPS[3*3+1]+=QC[1]*fsss[0][3]+WQ[1]*fsss[1][3]+ze2*dsss[1][0];
	    FSPS[3*3+2]+=QC[2]*fsss[0][3]+WQ[2]*fsss[1][3];
	    FSPS[4*3+0]+=QC[0]*fsss[0][4]+WQ[0]*fsss[1][4];
	    FSPS[4*3+1]+=QC[1]*fsss[0][4]+WQ[1]*fsss[1][4]+ze22*dsss[1][4];
	    FSPS[4*3+2]+=QC[2]*fsss[0][4]+WQ[2]*fsss[1][4]+ze2*dsss[1][1];
	    FSPS[5*3+0]+=QC[0]*fsss[0][5]+WQ[0]*fsss[1][5]+ze2*dsss[1][2];
	    FSPS[5*3+1]+=QC[1]*fsss[0][5]+WQ[1]*fsss[1][5];
	    FSPS[5*3+2]+=QC[2]*fsss[0][5]+WQ[2]*fsss[1][5]+ze22*dsss[1][5];
	    FSPS[6*3+0]+=QC[0]*fsss[0][6]+WQ[0]*fsss[1][6]+ze2*dsss[1][1];
	    FSPS[6*3+1]+=QC[1]*fsss[0][6]+WQ[1]*fsss[1][6]+ze22*dsss[1][3];
	    FSPS[6*3+2]+=QC[2]*fsss[0][6]+WQ[2]*fsss[1][6];
	    FSPS[7*3+0]+=QC[0]*fsss[0][7]+WQ[0]*fsss[1][7];
	    FSPS[7*3+1]+=QC[1]*fsss[0][7]+WQ[1]*fsss[1][7]+ze2*dsss[1][2];
	    FSPS[7*3+2]+=QC[2]*fsss[0][7]+WQ[2]*fsss[1][7]+ze22*dsss[1][4];
	    FSPS[8*3+0]+=QC[0]*fsss[0][8]+WQ[0]*fsss[1][8]+ze22*dsss[1][5];
	    FSPS[8*3+1]+=QC[1]*fsss[0][8]+WQ[1]*fsss[1][8];
	    FSPS[8*3+2]+=QC[2]*fsss[0][8]+WQ[2]*fsss[1][8]+ze2*dsss[1][0];
	    FSPS[9*3+0]+=QC[0]*fsss[0][9]+WQ[0]*fsss[1][9]+ze2*dsss[1][4];
	    FSPS[9*3+1]+=QC[1]*fsss[0][9]+WQ[1]*fsss[1][9]+ze2*dsss[1][5];
	    FSPS[9*3+2]+=QC[2]*fsss[0][9]+WQ[2]*fsss[1][9]+ze2*dsss[1][3];
	    // gsps (m=0)
	    GSPS[ 0*3+0]+=QC[0]*gsss[0][ 0]+WQ[0]*gsss[1][ 0]
			 + ze24*fsss[1][0];
	    GSPS[ 0*3+1]+=QC[1]*gsss[0][ 0]+WQ[1]*gsss[1][ 0];
	    GSPS[ 0*3+2]+=QC[2]*gsss[0][ 0]+WQ[2]*gsss[1][ 0];
	    GSPS[ 1*3+0]+=QC[0]*gsss[0][ 1]+WQ[0]*gsss[1][ 1];
	    GSPS[ 1*3+1]+=QC[1]*gsss[0][ 1]+WQ[1]*gsss[1][ 1]
			 + ze24*fsss[1][1];
	    GSPS[ 1*3+2]+=QC[2]*gsss[0][ 1]+WQ[2]*gsss[1][ 1];
	    GSPS[ 2*3+0]+=QC[0]*gsss[0][ 2]+WQ[0]*gsss[1][ 2];
	    GSPS[ 2*3+1]+=QC[1]*gsss[0][ 2]+WQ[1]*gsss[1][ 2];
	    GSPS[ 2*3+2]+=QC[2]*gsss[0][ 2]+WQ[2]*gsss[1][ 2]
			 + ze24*fsss[1][2];
	    GSPS[ 3*3+0]+=QC[0]*gsss[0][ 3]+WQ[0]*gsss[1][ 3]
			 + ze23*fsss[1][3];
	    GSPS[ 3*3+1]+=QC[1]*gsss[0][ 3]+WQ[1]*gsss[1][ 3]
			 + ze2*fsss[1][0];
	    GSPS[ 3*3+2]+=QC[2]*gsss[0][ 3]+WQ[2]*gsss[1][ 3];
	    GSPS[ 4*3+0]+=QC[0]*gsss[0][ 4]+WQ[0]*gsss[1][ 4];
	    GSPS[ 4*3+1]+=QC[1]*gsss[0][ 4]+WQ[1]*gsss[1][ 4]
			 + ze23*fsss[1][4];
	    GSPS[ 4*3+2]+=QC[2]*gsss[0][ 4]+WQ[2]*gsss[1][ 4]
			 + ze2*fsss[1][1];
	    GSPS[ 5*3+0]+=QC[0]*gsss[0][ 5]+WQ[0]*gsss[1][ 5]
			 + ze2*fsss[1][2];
	    GSPS[ 5*3+1]+=QC[1]*gsss[0][ 5]+WQ[1]*gsss[1][ 5];
	    GSPS[ 5*3+2]+=QC[2]*gsss[0][ 5]+WQ[2]*gsss[1][ 5]
			 + ze23*fsss[1][5];
	    GSPS[ 6*3+0]+=QC[0]*gsss[0][ 6]+WQ[0]*gsss[1][ 6]
			 + ze22*fsss[1][6];
	    GSPS[ 6*3+1]+=QC[1]*gsss[0][ 6]+WQ[1]*gsss[1][ 6]
			 + ze22*fsss[1][3];
	    GSPS[ 6*3+2]+=QC[2]*gsss[0][ 6]+WQ[2]*gsss[1][ 6];
	    GSPS[ 7*3+0]+=QC[0]*gsss[0][ 7]+WQ[0]*gsss[1][ 7];
	    GSPS[ 7*3+1]+=QC[1]*gsss[0][ 7]+WQ[1]*gsss[1][ 7]
			 + ze22*fsss[1][7];
	    GSPS[ 7*3+2]+=QC[2]*gsss[0][ 7]+WQ[2]*gsss[1][ 7]
			 + ze22*fsss[1][4];
	    GSPS[ 8*3+0]+=QC[0]*gsss[0][ 8]+WQ[0]*gsss[1][ 8]
			 + ze22*fsss[1][5];
	    GSPS[ 8*3+1]+=QC[1]*gsss[0][ 8]+WQ[1]*gsss[1][ 8];
	    GSPS[ 8*3+2]+=QC[2]*gsss[0][ 8]+WQ[2]*gsss[1][ 8]
			 + ze22*fsss[1][8];
	    GSPS[ 9*3+0]+=QC[0]*gsss[0][ 9]+WQ[0]*gsss[1][ 9]
			 + ze2*fsss[1][1];
	    GSPS[ 9*3+1]+=QC[1]*gsss[0][ 9]+WQ[1]*gsss[1][ 9]
			 + ze23*fsss[1][6];
	    GSPS[ 9*3+2]+=QC[2]*gsss[0][ 9]+WQ[2]*gsss[1][ 9];
	    GSPS[10*3+0]+=QC[0]*gsss[0][10]+WQ[0]*gsss[1][10];
	    GSPS[10*3+1]+=QC[1]*gsss[0][10]+WQ[1]*gsss[1][10]
			 + ze2*fsss[1][2];
	    GSPS[10*3+2]+=QC[2]*gsss[0][10]+WQ[2]*gsss[1][10]
			 + ze23*fsss[1][7];
	    GSPS[11*3+0]+=QC[0]*gsss[0][11]+WQ[0]*gsss[1][11]
			 + ze23*fsss[1][8];
	    GSPS[11*3+1]+=QC[1]*gsss[0][11]+WQ[1]*gsss[1][11];
	    GSPS[11*3+2]+=QC[2]*gsss[0][11]+WQ[2]*gsss[1][11]
			 + ze2*fsss[1][0];
	    GSPS[12*3+0]+=QC[0]*gsss[0][12]+WQ[0]*gsss[1][12]
			 + ze22*fsss[1][9];
	    GSPS[12*3+1]+=QC[1]*gsss[0][12]+WQ[1]*gsss[1][12]
			 + ze2*fsss[1][8];
	    GSPS[12*3+2]+=QC[2]*gsss[0][12]+WQ[2]*gsss[1][12]
			 + ze2*fsss[1][3];
	    GSPS[13*3+0]+=QC[0]*gsss[0][13]+WQ[0]*gsss[1][13]
			 + ze2*fsss[1][4];
	    GSPS[13*3+1]+=QC[1]*gsss[0][13]+WQ[1]*gsss[1][13]
			 + ze22*fsss[1][9];
	    GSPS[13*3+2]+=QC[2]*gsss[0][13]+WQ[2]*gsss[1][13]
			 + ze2*fsss[1][6];
	    GSPS[14*3+0]+=QC[0]*gsss[0][14]+WQ[0]*gsss[1][14]
			 + ze2*fsss[1][7];
	    GSPS[14*3+1]+=QC[1]*gsss[0][14]+WQ[1]*gsss[1][14]
			 + ze2*fsss[1][5];
	    GSPS[14*3+2]+=QC[2]*gsss[0][14]+WQ[2]*gsss[1][14]
			 + ze22*fsss[1][9];
	}	// klps
    }	// ijps
    // (D,P|P,S)
    for (c0=0; c0<3; c0++) {
        DPPS[0*3*3+0*3+c0]=FSPS[0*3+c0]-BA[0]*DSPS[0*3+c0];
        DPPS[0*3*3+1*3+c0]=FSPS[3*3+c0]-BA[1]*DSPS[0*3+c0];
        DPPS[0*3*3+2*3+c0]=FSPS[8*3+c0]-BA[2]*DSPS[0*3+c0];
        DPPS[1*3*3+0*3+c0]=FSPS[6*3+c0]-BA[0]*DSPS[1*3+c0];
        DPPS[1*3*3+1*3+c0]=FSPS[1*3+c0]-BA[1]*DSPS[1*3+c0];
        DPPS[1*3*3+2*3+c0]=FSPS[4*3+c0]-BA[2]*DSPS[1*3+c0];
        DPPS[2*3*3+0*3+c0]=FSPS[5*3+c0]-BA[0]*DSPS[2*3+c0];
        DPPS[2*3*3+1*3+c0]=FSPS[7*3+c0]-BA[1]*DSPS[2*3+c0];
        DPPS[2*3*3+2*3+c0]=FSPS[2*3+c0]-BA[2]*DSPS[2*3+c0];
        DPPS[3*3*3+0*3+c0]=FSPS[3*3+c0]-BA[0]*DSPS[3*3+c0];
        DPPS[3*3*3+1*3+c0]=FSPS[6*3+c0]-BA[1]*DSPS[3*3+c0];
        DPPS[3*3*3+2*3+c0]=FSPS[9*3+c0]-BA[2]*DSPS[3*3+c0];
        DPPS[4*3*3+0*3+c0]=FSPS[9*3+c0]-BA[0]*DSPS[4*3+c0];
        DPPS[4*3*3+1*3+c0]=FSPS[4*3+c0]-BA[1]*DSPS[4*3+c0];
        DPPS[4*3*3+2*3+c0]=FSPS[7*3+c0]-BA[2]*DSPS[4*3+c0];
        DPPS[5*3*3+0*3+c0]=FSPS[8*3+c0]-BA[0]*DSPS[5*3+c0];
        DPPS[5*3*3+1*3+c0]=FSPS[9*3+c0]-BA[1]*DSPS[5*3+c0];
        DPPS[5*3*3+2*3+c0]=FSPS[5*3+c0]-BA[2]*DSPS[5*3+c0];
    }
    // (F,P|P,S)
    for (c0=0; c0<3; c0++) {
        FPPS[0*3*3+0*3+c0]=GSPS[ 0*3+c0]-BA[0]*FSPS[0*3+c0];
        FPPS[0*3*3+1*3+c0]=GSPS[ 3*3+c0]-BA[1]*FSPS[0*3+c0];
        FPPS[0*3*3+2*3+c0]=GSPS[11*3+c0]-BA[2]*FSPS[0*3+c0];
        FPPS[1*3*3+0*3+c0]=GSPS[ 9*3+c0]-BA[0]*FSPS[1*3+c0];
        FPPS[1*3*3+1*3+c0]=GSPS[ 1*3+c0]-BA[1]*FSPS[1*3+c0];
        FPPS[1*3*3+2*3+c0]=GSPS[ 4*3+c0]-BA[2]*FSPS[1*3+c0];
        FPPS[2*3*3+0*3+c0]=GSPS[ 5*3+c0]-BA[0]*FSPS[2*3+c0];
        FPPS[2*3*3+1*3+c0]=GSPS[10*3+c0]-BA[1]*FSPS[2*3+c0];
        FPPS[2*3*3+2*3+c0]=GSPS[ 2*3+c0]-BA[2]*FSPS[2*3+c0];
        FPPS[3*3*3+0*3+c0]=GSPS[ 3*3+c0]-BA[0]*FSPS[3*3+c0];
        FPPS[3*3*3+1*3+c0]=GSPS[ 6*3+c0]-BA[1]*FSPS[3*3+c0];
        FPPS[3*3*3+2*3+c0]=GSPS[12*3+c0]-BA[2]*FSPS[3*3+c0];
        FPPS[4*3*3+0*3+c0]=GSPS[13*3+c0]-BA[0]*FSPS[4*3+c0];
        FPPS[4*3*3+1*3+c0]=GSPS[ 4*3+c0]-BA[1]*FSPS[4*3+c0];
        FPPS[4*3*3+2*3+c0]=GSPS[ 7*3+c0]-BA[2]*FSPS[4*3+c0];
        FPPS[5*3*3+0*3+c0]=GSPS[ 8*3+c0]-BA[0]*FSPS[5*3+c0];
        FPPS[5*3*3+1*3+c0]=GSPS[14*3+c0]-BA[1]*FSPS[5*3+c0];
        FPPS[5*3*3+2*3+c0]=GSPS[ 5*3+c0]-BA[2]*FSPS[5*3+c0];
        FPPS[6*3*3+0*3+c0]=GSPS[ 6*3+c0]-BA[0]*FSPS[6*3+c0];
        FPPS[6*3*3+1*3+c0]=GSPS[ 9*3+c0]-BA[1]*FSPS[6*3+c0];
        FPPS[6*3*3+2*3+c0]=GSPS[13*3+c0]-BA[2]*FSPS[6*3+c0];
        FPPS[7*3*3+0*3+c0]=GSPS[14*3+c0]-BA[0]*FSPS[7*3+c0];
        FPPS[7*3*3+1*3+c0]=GSPS[ 7*3+c0]-BA[1]*FSPS[7*3+c0];
        FPPS[7*3*3+2*3+c0]=GSPS[10*3+c0]-BA[2]*FSPS[7*3+c0];
        FPPS[8*3*3+0*3+c0]=GSPS[11*3+c0]-BA[0]*FSPS[8*3+c0];
        FPPS[8*3*3+1*3+c0]=GSPS[12*3+c0]-BA[1]*FSPS[8*3+c0];
        FPPS[8*3*3+2*3+c0]=GSPS[ 8*3+c0]-BA[2]*FSPS[8*3+c0];
        FPPS[9*3*3+0*3+c0]=GSPS[12*3+c0]-BA[0]*FSPS[9*3+c0];
        FPPS[9*3*3+1*3+c0]=GSPS[13*3+c0]-BA[1]*FSPS[9*3+c0];
        FPPS[9*3*3+2*3+c0]=GSPS[14*3+c0]-BA[2]*FSPS[9*3+c0];
    }
    // (D,D|P,S)
    for (c0=0; c0<3; c0++) {
        DDPS[0*6*3+0*3+c0]=FPPS[0*3*3+0*3+c0]-BA[0]*DPPS[0*3*3+0*3+c0];
        DDPS[0*6*3+1*3+c0]=FPPS[3*3*3+1*3+c0]-BA[1]*DPPS[0*3*3+1*3+c0];
        DDPS[0*6*3+2*3+c0]=FPPS[8*3*3+2*3+c0]-BA[2]*DPPS[0*3*3+2*3+c0];
        DDPS[0*6*3+3*3+c0]=FPPS[0*3*3+1*3+c0]-BA[0]*DPPS[0*3*3+1*3+c0];
        DDPS[0*6*3+4*3+c0]=FPPS[3*3*3+2*3+c0]-BA[1]*DPPS[0*3*3+2*3+c0];
        DDPS[0*6*3+5*3+c0]=FPPS[8*3*3+0*3+c0]-BA[2]*DPPS[0*3*3+0*3+c0];
        DDPS[1*6*3+0*3+c0]=FPPS[6*3*3+0*3+c0]-BA[0]*DPPS[1*3*3+0*3+c0];
        DDPS[1*6*3+1*3+c0]=FPPS[1*3*3+1*3+c0]-BA[1]*DPPS[1*3*3+1*3+c0];
        DDPS[1*6*3+2*3+c0]=FPPS[4*3*3+2*3+c0]-BA[2]*DPPS[1*3*3+2*3+c0];
        DDPS[1*6*3+3*3+c0]=FPPS[6*3*3+1*3+c0]-BA[0]*DPPS[1*3*3+1*3+c0];
        DDPS[1*6*3+4*3+c0]=FPPS[1*3*3+2*3+c0]-BA[1]*DPPS[1*3*3+2*3+c0];
        DDPS[1*6*3+5*3+c0]=FPPS[4*3*3+0*3+c0]-BA[2]*DPPS[1*3*3+0*3+c0];
        DDPS[2*6*3+0*3+c0]=FPPS[5*3*3+0*3+c0]-BA[0]*DPPS[2*3*3+0*3+c0];
        DDPS[2*6*3+1*3+c0]=FPPS[7*3*3+1*3+c0]-BA[1]*DPPS[2*3*3+1*3+c0];
        DDPS[2*6*3+2*3+c0]=FPPS[2*3*3+2*3+c0]-BA[2]*DPPS[2*3*3+2*3+c0];
        DDPS[2*6*3+3*3+c0]=FPPS[5*3*3+1*3+c0]-BA[0]*DPPS[2*3*3+1*3+c0];
        DDPS[2*6*3+4*3+c0]=FPPS[7*3*3+2*3+c0]-BA[1]*DPPS[2*3*3+2*3+c0];
        DDPS[2*6*3+5*3+c0]=FPPS[2*3*3+0*3+c0]-BA[2]*DPPS[2*3*3+0*3+c0];
        DDPS[3*6*3+0*3+c0]=FPPS[3*3*3+0*3+c0]-BA[0]*DPPS[3*3*3+0*3+c0];
        DDPS[3*6*3+1*3+c0]=FPPS[6*3*3+1*3+c0]-BA[1]*DPPS[3*3*3+1*3+c0];
        DDPS[3*6*3+2*3+c0]=FPPS[9*3*3+2*3+c0]-BA[2]*DPPS[3*3*3+2*3+c0];
        DDPS[3*6*3+3*3+c0]=FPPS[3*3*3+1*3+c0]-BA[0]*DPPS[3*3*3+1*3+c0];
        DDPS[3*6*3+4*3+c0]=FPPS[6*3*3+2*3+c0]-BA[1]*DPPS[3*3*3+2*3+c0];
        DDPS[3*6*3+5*3+c0]=FPPS[9*3*3+0*3+c0]-BA[2]*DPPS[3*3*3+0*3+c0];
        DDPS[4*6*3+0*3+c0]=FPPS[9*3*3+0*3+c0]-BA[0]*DPPS[4*3*3+0*3+c0];
        DDPS[4*6*3+1*3+c0]=FPPS[4*3*3+1*3+c0]-BA[1]*DPPS[4*3*3+1*3+c0];
        DDPS[4*6*3+2*3+c0]=FPPS[7*3*3+2*3+c0]-BA[2]*DPPS[4*3*3+2*3+c0];
        DDPS[4*6*3+3*3+c0]=FPPS[9*3*3+1*3+c0]-BA[0]*DPPS[4*3*3+1*3+c0];
        DDPS[4*6*3+4*3+c0]=FPPS[4*3*3+2*3+c0]-BA[1]*DPPS[4*3*3+2*3+c0];
        DDPS[4*6*3+5*3+c0]=FPPS[7*3*3+0*3+c0]-BA[2]*DPPS[4*3*3+0*3+c0];
        DDPS[5*6*3+0*3+c0]=FPPS[8*3*3+0*3+c0]-BA[0]*DPPS[5*3*3+0*3+c0];
        DDPS[5*6*3+1*3+c0]=FPPS[9*3*3+1*3+c0]-BA[1]*DPPS[5*3*3+1*3+c0];
        DDPS[5*6*3+2*3+c0]=FPPS[5*3*3+2*3+c0]-BA[2]*DPPS[5*3*3+2*3+c0];
        DDPS[5*6*3+3*3+c0]=FPPS[8*3*3+1*3+c0]-BA[0]*DPPS[5*3*3+1*3+c0];
        DDPS[5*6*3+4*3+c0]=FPPS[9*3*3+2*3+c0]-BA[1]*DPPS[5*3*3+2*3+c0];
        DDPS[5*6*3+5*3+c0]=FPPS[5*3*3+0*3+c0]-BA[2]*DPPS[5*3*3+0*3+c0];
    }
    for ( i=0, ix=0; i<6; i++ ) {
	coe0 = (i<3 ? ONE : sqr3 );
	for ( j=0; j<6; j++ ) {
	    coe = coe0 * (j<3? ONE : sqr3 );
	    for ( k=0; k<3; k++, ix++ ) DDPS[ix] *= coe;
	}
    }
}

/** １つのCS４重対に対して(dd,pp)タイプの２電子積分を計算する関数
 * @ingroup core-twoint
 * */
void twoint_core_ddpp__(
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double DDPP[6*6*3*3] ) {
    int ijps, klps, i, j, kl, ix, m, m1;
    int ab, ab01, ab10, ab00, c0;
    double cssss, zeta, dkab, xiza, eta, xizc, dk;
    double rz, re, ze2, ze22, ze23, ze24, zeta2, zeta23, eta2;
    double sqrho, rho, PQ2, T;
    double PC[3], PA[3], WP[3], QP[3], QC[3], WQ[3];
    double ssss[6+1], psss[5+1][3], dsss[4+1][6], fsss[3+1][10];
    double gsss[2+1][15];
    double psps[3*3], dsps[1+1][6*3], fsps[1+1][10*3], gsps[1+1][15*3];
    double DSDS[6*6], FSDS[10*6], GSDS[15*6];
    double DSPS[6*3], FSPS[10*3], GSPS[15*3];
    double DPPS[6*3*3], FPPS[10*3*3], DDPS[6*6*3];
    double DPDS[6*3*6], FPDS[10*3*6], DDDS[6*6*6];
    union _temp_ {
        double entity[6];
        double ssss;
        double psss[3];
        double dsss[6];
        double fsss[10];
        double gsss[15];
    } tmp;
    double sqr3, coe0, coe;
    sqr3 = sqrt(3.e0);
    for ( i=0; i< 6*3; i++ ) DSPS[i] = ZERO;
    for ( i=0; i<10*3; i++ ) FSPS[i] = ZERO;
    for ( i=0; i<15*3; i++ ) GSPS[i] = ZERO;
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
	    ze2   = rz  * eta2;
	    ze22  = rz  * eta;
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
	    for (m=0; m<=5; m++) {
		m1 = m+1;
		psss[m][0] = PA[0]*ssss[m] + WP[0]*ssss[m1];
		psss[m][1] = PA[1]*ssss[m] + WP[1]*ssss[m1];
		psss[m][2] = PA[2]*ssss[m] + WP[2]*ssss[m1];
	    }
	    // dsss (m=0,4)
	    for (m=0; m<=4; m++) {
		m1 = m+1;
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
	    for (m=0; m<=3; m++) {
		m1 = m+1;
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
	    for (m=0; m<=2; m++) {
		m1 = m+1;
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
	    for (m=0; m<=1; m++) {
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
	    for ( i=0; i<6*3; i++ ) DSPS[i] += dsps[0][i];
	    // fsps (m=0,1)
	    for (m=0; m<=1; m++) {
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
	    }
	    for ( i=0; i<10*3; i++ ) FSPS[i] += fsps[0][i];
	    // gsps (m=0,1)
	    for (m=0; m<=1; m++) {
		m1 = m+1;
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
	    for ( i=0; i<15*3; i++ ) GSPS[i] += gsps[0][i];
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
	    // (F,S|D,S)
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
    // (D,P|P,S)
    for (c0=0; c0<3; c0++) {
        DPPS[0*3*3+0*3+c0] = FSPS[0*3+c0] - BA[0]*DSPS[0*3+c0];
        DPPS[0*3*3+1*3+c0] = FSPS[3*3+c0] - BA[1]*DSPS[0*3+c0];
        DPPS[0*3*3+2*3+c0] = FSPS[8*3+c0] - BA[2]*DSPS[0*3+c0];
        DPPS[1*3*3+0*3+c0] = FSPS[6*3+c0] - BA[0]*DSPS[1*3+c0];
        DPPS[1*3*3+1*3+c0] = FSPS[1*3+c0] - BA[1]*DSPS[1*3+c0];
        DPPS[1*3*3+2*3+c0] = FSPS[4*3+c0] - BA[2]*DSPS[1*3+c0];
        DPPS[2*3*3+0*3+c0] = FSPS[5*3+c0] - BA[0]*DSPS[2*3+c0];
        DPPS[2*3*3+1*3+c0] = FSPS[7*3+c0] - BA[1]*DSPS[2*3+c0];
        DPPS[2*3*3+2*3+c0] = FSPS[2*3+c0] - BA[2]*DSPS[2*3+c0];
        DPPS[3*3*3+0*3+c0] = FSPS[3*3+c0] - BA[0]*DSPS[3*3+c0];
        DPPS[3*3*3+1*3+c0] = FSPS[6*3+c0] - BA[1]*DSPS[3*3+c0];
        DPPS[3*3*3+2*3+c0] = FSPS[9*3+c0] - BA[2]*DSPS[3*3+c0];
        DPPS[4*3*3+0*3+c0] = FSPS[9*3+c0] - BA[0]*DSPS[4*3+c0];
        DPPS[4*3*3+1*3+c0] = FSPS[4*3+c0] - BA[1]*DSPS[4*3+c0];
        DPPS[4*3*3+2*3+c0] = FSPS[7*3+c0] - BA[2]*DSPS[4*3+c0];
        DPPS[5*3*3+0*3+c0] = FSPS[8*3+c0] - BA[0]*DSPS[5*3+c0];
        DPPS[5*3*3+1*3+c0] = FSPS[9*3+c0] - BA[1]*DSPS[5*3+c0];
        DPPS[5*3*3+2*3+c0] = FSPS[5*3+c0] - BA[2]*DSPS[5*3+c0];
    }
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
    // (F,P|P,S)
    for (c0=0; c0<3; c0++) {
        FPPS[0*3*3+0*3+c0] = GSPS[ 0*3+c0] - BA[0]*FSPS[0*3+c0];
        FPPS[0*3*3+1*3+c0] = GSPS[ 3*3+c0] - BA[1]*FSPS[0*3+c0];
        FPPS[0*3*3+2*3+c0] = GSPS[11*3+c0] - BA[2]*FSPS[0*3+c0];
        FPPS[1*3*3+0*3+c0] = GSPS[ 9*3+c0] - BA[0]*FSPS[1*3+c0];
        FPPS[1*3*3+1*3+c0] = GSPS[ 1*3+c0] - BA[1]*FSPS[1*3+c0];
        FPPS[1*3*3+2*3+c0] = GSPS[ 4*3+c0] - BA[2]*FSPS[1*3+c0];
        FPPS[2*3*3+0*3+c0] = GSPS[ 5*3+c0] - BA[0]*FSPS[2*3+c0];
        FPPS[2*3*3+1*3+c0] = GSPS[10*3+c0] - BA[1]*FSPS[2*3+c0];
        FPPS[2*3*3+2*3+c0] = GSPS[ 2*3+c0] - BA[2]*FSPS[2*3+c0];
        FPPS[3*3*3+0*3+c0] = GSPS[ 3*3+c0] - BA[0]*FSPS[3*3+c0];
        FPPS[3*3*3+1*3+c0] = GSPS[ 6*3+c0] - BA[1]*FSPS[3*3+c0];
        FPPS[3*3*3+2*3+c0] = GSPS[12*3+c0] - BA[2]*FSPS[3*3+c0];
        FPPS[4*3*3+0*3+c0] = GSPS[13*3+c0] - BA[0]*FSPS[4*3+c0];
        FPPS[4*3*3+1*3+c0] = GSPS[ 4*3+c0] - BA[1]*FSPS[4*3+c0];
        FPPS[4*3*3+2*3+c0] = GSPS[ 7*3+c0] - BA[2]*FSPS[4*3+c0];
        FPPS[5*3*3+0*3+c0] = GSPS[ 8*3+c0] - BA[0]*FSPS[5*3+c0];
        FPPS[5*3*3+1*3+c0] = GSPS[14*3+c0] - BA[1]*FSPS[5*3+c0];
        FPPS[5*3*3+2*3+c0] = GSPS[ 5*3+c0] - BA[2]*FSPS[5*3+c0];
        FPPS[6*3*3+0*3+c0] = GSPS[ 6*3+c0] - BA[0]*FSPS[6*3+c0];
        FPPS[6*3*3+1*3+c0] = GSPS[ 9*3+c0] - BA[1]*FSPS[6*3+c0];
        FPPS[6*3*3+2*3+c0] = GSPS[13*3+c0] - BA[2]*FSPS[6*3+c0];
        FPPS[7*3*3+0*3+c0] = GSPS[14*3+c0] - BA[0]*FSPS[7*3+c0];
        FPPS[7*3*3+1*3+c0] = GSPS[ 7*3+c0] - BA[1]*FSPS[7*3+c0];
        FPPS[7*3*3+2*3+c0] = GSPS[10*3+c0] - BA[2]*FSPS[7*3+c0];
        FPPS[8*3*3+0*3+c0] = GSPS[11*3+c0] - BA[0]*FSPS[8*3+c0];
        FPPS[8*3*3+1*3+c0] = GSPS[12*3+c0] - BA[1]*FSPS[8*3+c0];
        FPPS[8*3*3+2*3+c0] = GSPS[ 8*3+c0] - BA[2]*FSPS[8*3+c0];
        FPPS[9*3*3+0*3+c0] = GSPS[12*3+c0] - BA[0]*FSPS[9*3+c0];
        FPPS[9*3*3+1*3+c0] = GSPS[13*3+c0] - BA[1]*FSPS[9*3+c0];
        FPPS[9*3*3+2*3+c0] = GSPS[14*3+c0] - BA[2]*FSPS[9*3+c0];
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
    // (D,D|P,S)
    for (c0=0; c0<3; c0++) {
        DDPS[0*6*3+0*3+c0]=FPPS[0*3*3+0*3+c0]-BA[0]*DPPS[0*3*3+0*3+c0];
        DDPS[0*6*3+1*3+c0]=FPPS[3*3*3+1*3+c0]-BA[1]*DPPS[0*3*3+1*3+c0];
        DDPS[0*6*3+2*3+c0]=FPPS[8*3*3+2*3+c0]-BA[2]*DPPS[0*3*3+2*3+c0];
        DDPS[0*6*3+3*3+c0]=FPPS[0*3*3+1*3+c0]-BA[0]*DPPS[0*3*3+1*3+c0];
        DDPS[0*6*3+4*3+c0]=FPPS[3*3*3+2*3+c0]-BA[1]*DPPS[0*3*3+2*3+c0];
        DDPS[0*6*3+5*3+c0]=FPPS[8*3*3+0*3+c0]-BA[2]*DPPS[0*3*3+0*3+c0];
        DDPS[1*6*3+0*3+c0]=FPPS[6*3*3+0*3+c0]-BA[0]*DPPS[1*3*3+0*3+c0];
        DDPS[1*6*3+1*3+c0]=FPPS[1*3*3+1*3+c0]-BA[1]*DPPS[1*3*3+1*3+c0];
        DDPS[1*6*3+2*3+c0]=FPPS[4*3*3+2*3+c0]-BA[2]*DPPS[1*3*3+2*3+c0];
        DDPS[1*6*3+3*3+c0]=FPPS[6*3*3+1*3+c0]-BA[0]*DPPS[1*3*3+1*3+c0];
        DDPS[1*6*3+4*3+c0]=FPPS[1*3*3+2*3+c0]-BA[1]*DPPS[1*3*3+2*3+c0];
        DDPS[1*6*3+5*3+c0]=FPPS[4*3*3+0*3+c0]-BA[2]*DPPS[1*3*3+0*3+c0];
        DDPS[2*6*3+0*3+c0]=FPPS[5*3*3+0*3+c0]-BA[0]*DPPS[2*3*3+0*3+c0];
        DDPS[2*6*3+1*3+c0]=FPPS[7*3*3+1*3+c0]-BA[1]*DPPS[2*3*3+1*3+c0];
        DDPS[2*6*3+2*3+c0]=FPPS[2*3*3+2*3+c0]-BA[2]*DPPS[2*3*3+2*3+c0];
        DDPS[2*6*3+3*3+c0]=FPPS[5*3*3+1*3+c0]-BA[0]*DPPS[2*3*3+1*3+c0];
        DDPS[2*6*3+4*3+c0]=FPPS[7*3*3+2*3+c0]-BA[1]*DPPS[2*3*3+2*3+c0];
        DDPS[2*6*3+5*3+c0]=FPPS[2*3*3+0*3+c0]-BA[2]*DPPS[2*3*3+0*3+c0];
        DDPS[3*6*3+0*3+c0]=FPPS[3*3*3+0*3+c0]-BA[0]*DPPS[3*3*3+0*3+c0];
        DDPS[3*6*3+1*3+c0]=FPPS[6*3*3+1*3+c0]-BA[1]*DPPS[3*3*3+1*3+c0];
        DDPS[3*6*3+2*3+c0]=FPPS[9*3*3+2*3+c0]-BA[2]*DPPS[3*3*3+2*3+c0];
        DDPS[3*6*3+3*3+c0]=FPPS[3*3*3+1*3+c0]-BA[0]*DPPS[3*3*3+1*3+c0];
        DDPS[3*6*3+4*3+c0]=FPPS[6*3*3+2*3+c0]-BA[1]*DPPS[3*3*3+2*3+c0];
        DDPS[3*6*3+5*3+c0]=FPPS[9*3*3+0*3+c0]-BA[2]*DPPS[3*3*3+0*3+c0];
        DDPS[4*6*3+0*3+c0]=FPPS[9*3*3+0*3+c0]-BA[0]*DPPS[4*3*3+0*3+c0];
        DDPS[4*6*3+1*3+c0]=FPPS[4*3*3+1*3+c0]-BA[1]*DPPS[4*3*3+1*3+c0];
        DDPS[4*6*3+2*3+c0]=FPPS[7*3*3+2*3+c0]-BA[2]*DPPS[4*3*3+2*3+c0];
        DDPS[4*6*3+3*3+c0]=FPPS[9*3*3+1*3+c0]-BA[0]*DPPS[4*3*3+1*3+c0];
        DDPS[4*6*3+4*3+c0]=FPPS[4*3*3+2*3+c0]-BA[1]*DPPS[4*3*3+2*3+c0];
        DDPS[4*6*3+5*3+c0]=FPPS[7*3*3+0*3+c0]-BA[2]*DPPS[4*3*3+0*3+c0];
        DDPS[5*6*3+0*3+c0]=FPPS[8*3*3+0*3+c0]-BA[0]*DPPS[5*3*3+0*3+c0];
        DDPS[5*6*3+1*3+c0]=FPPS[9*3*3+1*3+c0]-BA[1]*DPPS[5*3*3+1*3+c0];
        DDPS[5*6*3+2*3+c0]=FPPS[5*3*3+2*3+c0]-BA[2]*DPPS[5*3*3+2*3+c0];
        DDPS[5*6*3+3*3+c0]=FPPS[8*3*3+1*3+c0]-BA[0]*DPPS[5*3*3+1*3+c0];
        DDPS[5*6*3+4*3+c0]=FPPS[9*3*3+2*3+c0]-BA[1]*DPPS[5*3*3+2*3+c0];
        DDPS[5*6*3+5*3+c0]=FPPS[5*3*3+0*3+c0]-BA[2]*DPPS[5*3*3+0*3+c0];
    }
    // (D,D|D,S)
    for (c0=0; c0<6; c0++) {
        DDDS[0*6*6+0*6+c0]=FPDS[0*3*6+0*6+c0]-BA[0]*DPDS[0*3*6+0*6+c0];
        DDDS[0*6*6+1*6+c0]=FPDS[3*3*6+1*6+c0]-BA[1]*DPDS[0*3*6+1*6+c0];
        DDDS[0*6*6+2*6+c0]=FPDS[8*3*6+2*6+c0]-BA[2]*DPDS[0*3*6+2*6+c0];
        DDDS[0*6*6+3*6+c0]=FPDS[0*3*6+1*6+c0]-BA[0]*DPDS[0*3*6+1*6+c0];
        DDDS[0*6*6+4*6+c0]=FPDS[3*3*6+2*6+c0]-BA[1]*DPDS[0*3*6+2*6+c0];
        DDDS[0*6*6+5*6+c0]=FPDS[8*3*6+0*6+c0]-BA[2]*DPDS[0*3*6+0*6+c0];
        DDDS[1*6*6+0*6+c0]=FPDS[6*3*6+0*6+c0]-BA[0]*DPDS[1*3*6+0*6+c0];
        DDDS[1*6*6+1*6+c0]=FPDS[1*3*6+1*6+c0]-BA[1]*DPDS[1*3*6+1*6+c0];
        DDDS[1*6*6+2*6+c0]=FPDS[4*3*6+2*6+c0]-BA[2]*DPDS[1*3*6+2*6+c0];
        DDDS[1*6*6+3*6+c0]=FPDS[6*3*6+1*6+c0]-BA[0]*DPDS[1*3*6+1*6+c0];
        DDDS[1*6*6+4*6+c0]=FPDS[1*3*6+2*6+c0]-BA[1]*DPDS[1*3*6+2*6+c0];
        DDDS[1*6*6+5*6+c0]=FPDS[4*3*6+0*6+c0]-BA[2]*DPDS[1*3*6+0*6+c0];
        DDDS[2*6*6+0*6+c0]=FPDS[5*3*6+0*6+c0]-BA[0]*DPDS[2*3*6+0*6+c0];
        DDDS[2*6*6+1*6+c0]=FPDS[7*3*6+1*6+c0]-BA[1]*DPDS[2*3*6+1*6+c0];
        DDDS[2*6*6+2*6+c0]=FPDS[2*3*6+2*6+c0]-BA[2]*DPDS[2*3*6+2*6+c0];
        DDDS[2*6*6+3*6+c0]=FPDS[5*3*6+1*6+c0]-BA[0]*DPDS[2*3*6+1*6+c0];
        DDDS[2*6*6+4*6+c0]=FPDS[7*3*6+2*6+c0]-BA[1]*DPDS[2*3*6+2*6+c0];
        DDDS[2*6*6+5*6+c0]=FPDS[2*3*6+0*6+c0]-BA[2]*DPDS[2*3*6+0*6+c0];
        DDDS[3*6*6+0*6+c0]=FPDS[3*3*6+0*6+c0]-BA[0]*DPDS[3*3*6+0*6+c0];
        DDDS[3*6*6+1*6+c0]=FPDS[6*3*6+1*6+c0]-BA[1]*DPDS[3*3*6+1*6+c0];
        DDDS[3*6*6+2*6+c0]=FPDS[9*3*6+2*6+c0]-BA[2]*DPDS[3*3*6+2*6+c0];
        DDDS[3*6*6+3*6+c0]=FPDS[3*3*6+1*6+c0]-BA[0]*DPDS[3*3*6+1*6+c0];
        DDDS[3*6*6+4*6+c0]=FPDS[6*3*6+2*6+c0]-BA[1]*DPDS[3*3*6+2*6+c0];
        DDDS[3*6*6+5*6+c0]=FPDS[9*3*6+0*6+c0]-BA[2]*DPDS[3*3*6+0*6+c0];
        DDDS[4*6*6+0*6+c0]=FPDS[9*3*6+0*6+c0]-BA[0]*DPDS[4*3*6+0*6+c0];
        DDDS[4*6*6+1*6+c0]=FPDS[4*3*6+1*6+c0]-BA[1]*DPDS[4*3*6+1*6+c0];
        DDDS[4*6*6+2*6+c0]=FPDS[7*3*6+2*6+c0]-BA[2]*DPDS[4*3*6+2*6+c0];
        DDDS[4*6*6+3*6+c0]=FPDS[9*3*6+1*6+c0]-BA[0]*DPDS[4*3*6+1*6+c0];
        DDDS[4*6*6+4*6+c0]=FPDS[4*3*6+2*6+c0]-BA[1]*DPDS[4*3*6+2*6+c0];
        DDDS[4*6*6+5*6+c0]=FPDS[7*3*6+0*6+c0]-BA[2]*DPDS[4*3*6+0*6+c0];
        DDDS[5*6*6+0*6+c0]=FPDS[8*3*6+0*6+c0]-BA[0]*DPDS[5*3*6+0*6+c0];
        DDDS[5*6*6+1*6+c0]=FPDS[9*3*6+1*6+c0]-BA[1]*DPDS[5*3*6+1*6+c0];
        DDDS[5*6*6+2*6+c0]=FPDS[5*3*6+2*6+c0]-BA[2]*DPDS[5*3*6+2*6+c0];
        DDDS[5*6*6+3*6+c0]=FPDS[8*3*6+1*6+c0]-BA[0]*DPDS[5*3*6+1*6+c0];
        DDDS[5*6*6+4*6+c0]=FPDS[9*3*6+2*6+c0]-BA[1]*DPDS[5*3*6+2*6+c0];
        DDDS[5*6*6+5*6+c0]=FPDS[5*3*6+0*6+c0]-BA[2]*DPDS[5*3*6+0*6+c0];
    }
    // (a,b|c,d+1) = (a,b|c+1,d) - DC(a,b|c,d)
    // (D,D|P,P)
    for (ab=0, ab01=0, ab10=0, ab00=0;
            ab<6*6; ab++, ab01 += 3*3, ab10 += 6, ab00 += 3) {
        DDPP[ab01+0*3+0] = DDDS[ab10+0] - DC[0]*DDPS[ab00+0];
        DDPP[ab01+0*3+1] = DDDS[ab10+3] - DC[1]*DDPS[ab00+0];
        DDPP[ab01+0*3+2] = DDDS[ab10+5] - DC[2]*DDPS[ab00+0];
        DDPP[ab01+1*3+0] = DDDS[ab10+3] - DC[0]*DDPS[ab00+1];
        DDPP[ab01+1*3+1] = DDDS[ab10+1] - DC[1]*DDPS[ab00+1];
        DDPP[ab01+1*3+2] = DDDS[ab10+4] - DC[2]*DDPS[ab00+1];
        DDPP[ab01+2*3+0] = DDDS[ab10+5] - DC[0]*DDPS[ab00+2];
        DDPP[ab01+2*3+1] = DDDS[ab10+4] - DC[1]*DDPS[ab00+2];
        DDPP[ab01+2*3+2] = DDDS[ab10+2] - DC[2]*DDPS[ab00+2];
    }
    //
    for ( i=0, ix=0; i<6; i++ ) {
	coe0 = (i<3 ? ONE : sqr3 );
	for ( j=0; j<6; j++ ) {
	    coe = coe0 * (j<3? ONE : sqr3 );
	    for ( kl=0; kl<3*3; kl++, ix++ ) DDPP[ix] *= coe;
	}
    }
}
