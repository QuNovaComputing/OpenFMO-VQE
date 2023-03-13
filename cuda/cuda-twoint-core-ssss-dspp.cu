#include <stdio.h>
#include <math.h>

#include "cuda-twoint-core.h"
#include "cuda-fmt.h"

#ifndef ZERO
#define ZERO 0.e0
#endif
#ifndef HALF
#define HALF .5e0
#endif
    
/** １つのCS４重対に対して(ss,ss)タイプの２電子積分を計算する関数
 * @ingroup core-twoint
 * */
__device__ void gpu_twoint_core_ssss_(
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double SSSS[1] ) {
//    int ijps, klps, i;
//    double ssss, cssss, zeta, dkab, xiza, eta, xizc, dk;
//    double sqrho, rho, PQ2, T;
//    double PC[3], PQ[3];
    SSSS[0] = ZERO;
    for (int ijps=0; ijps<(*nijps); ijps++ ) {
      double zeta, dkab, xiza;
      double PC[3];
	zeta = LDG(vzeta[ijps]);
	dkab = LDG(vdkab[ijps]);
	xiza = LDG(vxiza[ijps]);
#pragma unroll
	for (int i=0; i<3; i++ ) PC[i] = AC[i] + xiza*BA[i];

	for (int klps=0; klps<(*nklps); klps++ ) {
          double eta, dk, xizc;
          double ssss, cssss;
          double sqrho, rho, PQ2, T;
	    eta  = veta[klps];
	    dk   = dkab * LDG(vdkcd[klps]);
	    sqrho = sqrt(1.e0/(zeta+eta));
	    cssss = sqrho * dk;
	    rho   = sqrho*sqrho;
	    xizc = LDG(vxizc[klps]);
//#pragma unroll
//	    for ( i=0; i<3; i++ ) PQ[i] = PC[i] - xizc*DC[i];
//	    PQ2   = PQ[0]*PQ[0] + PQ[1]*PQ[1] + PQ[2]*PQ[2];
//	    T     = rho * PQ2;
	    T = ZERO;
#pragma unroll
	    for (int i=0; i<3; i++ ) {
              double PQ = PC[i] - xizc*DC[i];
	      T   += PQ*PQ;
            }
	    T     *= rho;
	    gpu_fmt( &ssss, 0, T, cssss );
#if 0
#if   CUDA_FMT_M == 3
	    gpu_fmt0_method3( T, cssss, &ssss );
#elif CUDA_FMT_M == 2
	    gpu_fmt0_method2( T, cssss, &ssss );
#elif CUDA_FMT_M == 1
	    gpu_fmt0_method1( T, cssss, &ssss );
#else
	    gpu_fmt0( &ssss, T, cssss );
#endif
#endif
	    SSSS[0] += ssss;
	}
    }
}

/** １つのCS４重対に対して(ps,ss)タイプの２電子積分を計算する関数
 * @ingroup core-twoint
 * */
__device__ void gpu_twoint_core_psss_(
	 const int *nijps,  const double vzeta[],  const double vdkab[],
	 const double vxiza[],  const double BA[3],
	 const int *nklps,  const double veta[],  const double vdkcd[],
	 const double vxizc[],  const double DC[3],  const double AC[3],
	 double PSSS[3] ) {
    int ijps, klps;
//    double ssss[1+1], cssss, zeta, dkab, xiza, eta, xizc, dk;
    double zeta, dkab, xiza;
//    double rz;
//    double sqrho, rho, PQ2, T;
//    double PC[3], PA[3], WP[3], QP[3];
    double PC[3], PA[3];

#pragma unroll
    for (int i=0; i<3; i++ ) PSSS[i] = ZERO;
    for (int ijps=0; ijps<(*nijps); ijps++ ) {
	zeta = LDG(vzeta[ijps]);
	dkab = LDG(vdkab[ijps]);
	xiza = LDG(vxiza[ijps]);
#pragma unroll
	for (int i=0; i<3; i++ ) {
	    PA[i] = xiza * BA[i];
	    PC[i] = AC[i] + PA[i];
	}
	for (int klps=0; klps<(*nklps); klps++ ) {
          double WP[3];
//          double ssss[1+1], cssss, eta, xizc, dk;
          double cssss, eta, xizc, dk;
          double sqrho, rho, T;
          double rz;
	    eta  = LDG(veta[klps]);
	    dk   = dkab * LDG(vdkcd[klps]);
	    xizc = LDG(vxizc[klps]);
	    sqrho = sqrt(1.e0/(zeta+eta));
	    cssss = sqrho * dk;
	    rho   = sqrho*sqrho;
	    rz    = rho * zeta;
            T=ZERO;
#pragma unroll
	    for (int i=0; i<3; i++ ) {
		double qp = xizc*DC[i] - PC[i];
		WP[i] = rz * qp;
		T  += qp*qp;
	    }
	    T     *= rho;
            {
            double ssss[1+1];
	    //gpu_fmt( ssss, 1, T, cssss );
#if defined(CUDA_FMT_M_K1) && CUDA_ARCH >= 350 && CUDA_FMT_M_NEXP == 6
                        gpu_fmt1_method3( T, cssss, ssss );
#else
#if   CUDA_FMT_M == 3
	    gpu_fmt1_method3( T, cssss, ssss );
#elif CUDA_FMT_M == 2
	    gpu_fmt1_method2( T, cssss, ssss );
#elif CUDA_FMT_M == 1
	    gpu_fmt1_method1( T, cssss, ssss );
#else
	    gpu_fmt1( ssss, T, cssss );
#endif
#endif
#pragma unroll
	    for (int i=0; i<3; i++ ) PSSS[i] += PA[i]*ssss[0]+WP[i]*ssss[1];
            }
	}
    }
}

/** １つのCS４重対に対して(ps,ps)タイプの２電子積分を計算する関数
 * @ingroup core-twoint
 * */
__device__ void gpu_twoint_core_psps_(
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double PSPS[3*3]) {
//    int ijps, klps, i, m;
//    double cssss, zeta, dkab, xiza, eta, xizc, dk;
//    double rz, ze2;
//    double sqrho, rho, PQ2, T;
//    double PC[3], PA[3], WP[3], QP[3], QC[3], WQ[3];
    double PC[3], PA[3];
//    double ssss[1+2], psss[1+1][3];
#pragma unroll
    for (int i=0; i<3*3; i++ ) PSPS[i] = ZERO;
    for (int ijps=0; ijps<(*nijps); ijps++ ) {
	double zeta = LDG(vzeta[ijps]);
	double dkab = LDG(vdkab[ijps]);
	double xiza = LDG(vxiza[ijps]);
#pragma unroll
	for (int i=0; i<3; i++ ) {
	    PA[i] = xiza * BA[i];
	    PC[i] = AC[i] + PA[i];
	}
	for (int klps=0; klps<(*nklps); klps++ ) {
          double eta, dk, xizc;
          double PQ2;
          double sqrho, rho, rz, ze2;
          double T, cssss;
          double WP[3], QP[3], QC[3], WQ[3];
//          double ssss[1+2], psss[1+1][3];
          double ssss[1+2];
	    eta  = LDG(veta[klps]);
	    dk   = dkab * LDG(vdkcd[klps]);
	    xizc = LDG(vxizc[klps]);
	    sqrho = sqrt(1.e0/(zeta+eta));
	    rho   = sqrho*sqrho;
	    cssss = sqrho * dk;
	    rz    = rho * zeta;
	    ze2   = HALF * rz * eta;
	    T  = ZERO;
#pragma unroll
	    for (int i=0; i<3; i++ ) {
		QC[i] = xizc*DC[i];
		QP[i] = QC[i] - PC[i];
		T  += QP[i]*QP[i];
	    }
	    T     *= rho;
	    //gpu_fmt( ssss, 2, T, cssss );
#if   CUDA_FMT_M == 3
	    gpu_fmt2_method3( T, cssss, ssss );
#elif CUDA_FMT_M == 2
	    gpu_fmt2_method2( T, cssss, ssss );
#elif CUDA_FMT_M == 1
	    gpu_fmt2_method1( T, cssss, ssss );
#else
	    gpu_fmt2( ssss, T, cssss );
#endif
#pragma unroll
	    for (int i=0; i<3; i++ ) {
		WP[i]=rz*QP[i];
		WQ[i] = WP[i] - QP[i];
	    }
	    // psss
            /*
#pragma unroll
	    for ( m=0; m<=1; m++ ) {
#pragma unroll
		for ( i=0; i<3; i++ )
		    psss[m][i] = PA[i]*ssss[m]+WP[i]*ssss[m+1];
	    }
	    // psps
	    PSPS[0*3+0] += QC[0]*psss[0][0]+WQ[0]*psss[1][0]+ze2*ssss[1];
	    PSPS[0*3+1] += QC[1]*psss[0][0]+WQ[1]*psss[1][0];
	    PSPS[0*3+2] += QC[2]*psss[0][0]+WQ[2]*psss[1][0];
	    PSPS[1*3+0] += QC[0]*psss[0][1]+WQ[0]*psss[1][1];
	    PSPS[1*3+1] += QC[1]*psss[0][1]+WQ[1]*psss[1][1]+ze2*ssss[1];;
	    PSPS[1*3+2] += QC[2]*psss[0][1]+WQ[2]*psss[1][1];
	    PSPS[2*3+0] += QC[0]*psss[0][2]+WQ[0]*psss[1][2];
	    PSPS[2*3+1] += QC[1]*psss[0][2]+WQ[1]*psss[1][2];
	    PSPS[2*3+2] += QC[2]*psss[0][2]+WQ[2]*psss[1][2]+ze2*ssss[1];;;
            */
            /*
            {
            double d2;
            double d0, d1;
	    d0 = PA[0]*ssss[0]+WP[0]*ssss[1];
	    d1 = PA[0]*ssss[1]+WP[0]*ssss[2];
            d2 = ze2*ssss[1];
	    PSPS[0*3+0] += QC[0]*d0+WQ[0]*d1+d2;
	    PSPS[0*3+1] += QC[1]*d0+WQ[1]*d1;
	    PSPS[0*3+2] += QC[2]*d0+WQ[2]*d1;
	    d0 = PA[1]*ssss[0]+WP[1]*ssss[1];
	    d1 = PA[1]*ssss[1]+WP[1]*ssss[2];
	    PSPS[1*3+0] += QC[0]*d0+WQ[0]*d1;
	    PSPS[1*3+1] += QC[1]*d0+WQ[1]*d1+d2;
	    PSPS[1*3+2] += QC[2]*d0+WQ[2]*d1;
	    d0 = PA[2]*ssss[0]+WP[2]*ssss[1];
	    d1 = PA[2]*ssss[1]+WP[2]*ssss[2];
	    PSPS[2*3+0] += QC[0]*d0+WQ[0]*d1;
	    PSPS[2*3+1] += QC[1]*d0+WQ[1]*d1;
	    PSPS[2*3+2] += QC[2]*d0+WQ[2]*d1+d2;
            }
            */
            {
            double d2 = ze2*ssss[1];
#pragma unroll
            for (int i=0; i<3; i++) {
	      double d0 = PA[i]*ssss[0]+WP[i]*ssss[1];
	      double d1 = PA[i]*ssss[1]+WP[i]*ssss[2];
	      PSPS[i*3+0] += QC[0]*d0+WQ[0]*d1;
	      PSPS[i*3+1] += QC[1]*d0+WQ[1]*d1;
	      PSPS[i*3+2] += QC[2]*d0+WQ[2]*d1;
	      PSPS[i*3+i] += d2;
            }
            }
	}
    }
}

/** １つのCS４重対に対して(pp,ss)タイプの２電子積分を計算する関数
 * @ingroup core-twoint
 * */
__device__ void gpu_twoint_core_ppss_(
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double PPSS[3*3] ) {
//    int ijps, klps, i, m;
//    double cssss, zeta, dkab, xiza, eta, xizc, dk;
//    double rz, zeta2, tmp;
//    double sqrho, rho, PQ2, T;
//    double PC[3], PA[3], WP[3], QP[3];
//    double ssss[1+2], psss[1+1][3];
    double PSSS[3], DSSS[6];
#pragma unroll
    for (int i=0; i<3; i++ ) PSSS[i] = ZERO;
#pragma unroll
    for (int i=0; i<6; i++ ) DSSS[i] = ZERO;
    for (int ijps=0; ijps<(*nijps); ijps++ ) {
        double PC[3], PA[3];
	double zeta = LDG(vzeta[ijps]);
	double dkab = LDG(vdkab[ijps]);
	double xiza = LDG(vxiza[ijps]);
	double zeta2 = HALF * zeta;
#pragma unroll
	for (int i=0; i<3; i++ ) {
//	    PC[i] = AC[i] + xiza*BA[i];
//	    PA[i] = xiza * BA[i];
	    PC[i] = AC[i] + (PA[i] = xiza * BA[i]);
	}
	for (int klps=0; klps<(*nklps); klps++ ) {
            double WP[3], QP[3];
            double ssss[1+2];
            double cssss, T;
            double sqrho, rho, rz;
            double PQ2 = ZERO;
	    double eta  = LDG(veta[klps]);
	    double xizc = LDG(vxizc[klps]);
	    double dk   = dkab * LDG(vdkcd[klps]);
#pragma unroll
	    for (int i=0; i<3; i++ ) {
		QP[i] = xizc*DC[i] - PC[i];
		PQ2  += QP[i]*QP[i];
	    }
	    sqrho = sqrt(1.e0/(zeta+eta));
	    rho   = sqrho*sqrho;
	    rz    = rho * zeta;
#pragma unroll
	    for (int i=0; i<3; i++ ) WP[i]=rz*QP[i];
	    T = rho * PQ2;
	    cssss = sqrho * dk;
	    //gpu_fmt( ssss, 2, T, cssss );
#if   CUDA_FMT_M == 3
	    gpu_fmt2_method3( T, cssss, ssss );
#elif CUDA_FMT_M == 2
	    gpu_fmt2_method2( T, cssss, ssss );
#elif CUDA_FMT_M == 1
	    gpu_fmt2_method1( T, cssss, ssss );
#else
	    gpu_fmt2( ssss, T, cssss );
#endif
            /*
	    // psss
#pragma unroll
	    for ( m=0; m<=1; m++ ) {
#pragma unroll
		for ( i=0; i<3; i++ )
		    psss[m][i] = PA[i]*ssss[m]+WP[i]*ssss[m+1];
	    }
#pragma unroll
	    for (i=0; i<3; i++ ) PSSS[i] += psss[0][i];
	    // dsss
	    tmp = ssss[0] - rz*ssss[1];
	    DSSS[0] += PA[0]*psss[0][0] + WP[0]*psss[1][0] + zeta2*tmp;
	    DSSS[1] += PA[1]*psss[0][1] + WP[1]*psss[1][1] + zeta2*tmp;
	    DSSS[2] += PA[2]*psss[0][2] + WP[2]*psss[1][2] + zeta2*tmp;
	    DSSS[3] += PA[0]*psss[0][1] + WP[0]*psss[1][1];
	    DSSS[4] += PA[1]*psss[0][2] + WP[1]*psss[1][2];
	    DSSS[5] += PA[2]*psss[0][0] + WP[2]*psss[1][0];
            */
            /*
            {
	    double tmp = (ssss[0] - rz*ssss[1])*zeta2;
            double d0, d1;
            d0 = PA[0]*ssss[0]+WP[0]*ssss[0+1];
            d1 = PA[0]*ssss[1]+WP[0]*ssss[1+1];
            PSSS[0] += d0;
	    DSSS[0] += PA[0]*d0 + WP[0]*d1 + tmp;
	    DSSS[5] += PA[2]*d0 + WP[2]*d1;
            d0 = PA[1]*ssss[0]+WP[1]*ssss[0+1];
            d1 = PA[1]*ssss[1]+WP[1]*ssss[1+1];
            PSSS[1] += d0;
	    DSSS[1] += PA[1]*d0 + WP[1]*d1 + tmp;
	    DSSS[3] += PA[0]*d0 + WP[0]*d1;
            d0 = PA[2]*ssss[0]+WP[2]*ssss[0+1];
            d1 = PA[2]*ssss[1]+WP[2]*ssss[1+1];
            PSSS[2] += d0;
	    DSSS[2] += PA[2]*d0 + WP[2]*d1 + tmp;
	    DSSS[4] += PA[1]*d0 + WP[1]*d1;
            }
            */
            {
	    double tmp = (ssss[0] - rz*ssss[1])*zeta2;
#pragma unroll
            for (int i=0; i<3; i++) {
            int j;
            double d0, d1;
            d0 = PA[i]*ssss[0]+WP[i]*ssss[0+1];
            d1 = PA[i]*ssss[1]+WP[i]*ssss[1+1];
            PSSS[i] += d0;
	    DSSS[i] += PA[i]*d0 + WP[i]*d1 + tmp;
            j=(i+2)%3;
	    DSSS[j+3] += PA[j]*d0 + WP[j]*d1;
            }
            }
	}	// klps
    }	// ijps
    // (P,P|S,S)
    /*
    PPSS[0*3+0] = DSSS[0] - BA[0]*PSSS[0];
    PPSS[0*3+1] = DSSS[3] - BA[1]*PSSS[0];
    PPSS[0*3+2] = DSSS[5] - BA[2]*PSSS[0];
    PPSS[1*3+0] = DSSS[3] - BA[0]*PSSS[1];
    PPSS[1*3+1] = DSSS[1] - BA[1]*PSSS[1];
    PPSS[1*3+2] = DSSS[4] - BA[2]*PSSS[1];
    PPSS[2*3+0] = DSSS[5] - BA[0]*PSSS[2];
    PPSS[2*3+1] = DSSS[4] - BA[1]*PSSS[2];
    PPSS[2*3+2] = DSSS[2] - BA[2]*PSSS[2];
    */
    // 035
    // 314
    // 542
    {
    int ix3[9]={0,3,5,3,1,4,5,4,2};
#pragma unroll
    for (int i=0; i<3; i++) {
#pragma unroll
      for (int j=0; j<3; j++) {
        int ix=i*3+j;
        PPSS[ix] = DSSS[ix3[ix]] - BA[j]*PSSS[i];
      }
    }
    }
}

/** １つのCS４重対に対して(pp,ps)タイプの２電子積分を計算する関数
 * @ingroup core-twoint
 * */
__device__ void gpu_twoint_core_ppps_(
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double PPPS[3*3*3] ) {
//    int ijps, klps, i, m, m1, c0;
//    int ijps, klps, i;
//    double cssss, zeta, dkab, xiza, eta, xizc, dk;
//    double rz, re, ze2, ze22, zeta2, tmp;
//    double sqrho, rho, PQ2, T;
//    double PC[3], PA[3], WP[3], QP[3], QC[3], WQ[3];
//    double ssss[3+1], psss[2+1][3], dsss[1+1][6];
    double PSPS[3*3], DSPS[6*3];
#pragma unroll
    for (int i=0; i<3*3; i++ ) PSPS[i] = ZERO;
#pragma unroll
    for (int i=0; i<6*3; i++ ) DSPS[i] = ZERO;
    for (int ijps=0; ijps<(*nijps); ijps++ ) {
        double PC[3], PA[3];
//	double zeta = vzeta[ijps];
	double dkab = LDG(vdkab[ijps]);
//	double xiza = vxiza[ijps];
	double xiza = LDG(vxiza[ijps]);
//	double zeta2 = HALF * zeta;
	double zeta2 = HALF * LDG(vzeta[ijps]);
#pragma unroll
	for (int i=0; i<3; i++ ) {
//	    PC[i] = AC[i] + xiza*BA[i];
//	    PA[i] = xiza * BA[i];
	    PC[i] = AC[i] + (PA[i] = xiza * BA[i]);
	}
	for (int klps=0; klps<(*nklps); klps++ ) {
            double ze2, ze22;
            double WP[3], QC[3], WQ[3];
            double ssss[3+1];
            //double cssss, T;
            double zs0, zs1;
            {
            double rz;
            double QP[3];
            //double re;
	    double eta  = LDG(veta[klps]);
	    //double dk   = dkab * vdkcd[klps];
	    double xizc = LDG(vxizc[klps]);
            //double sqrho, rho;
            double rho;
	    //double PQ2  = ZERO;
            double dtmp;
            dtmp = ZERO;
#pragma unroll
	    for (int i=0; i<3; i++ ) {
//		QC[i] = xizc*DC[i];
//		QP[i] = xizc*DC[i] - PC[i];
		QP[i] = (QC[i] = xizc*DC[i]) - PC[i];
		//PQ2  += QP[i]*QP[i];
		dtmp  += QP[i]*QP[i];
	    }
            zs0 = dtmp;
            dtmp = zeta2 * 2;
	    //sqrho = sqrt(1.e0/(zeta+eta));
	    //sqrho = sqrt(1.e0/(dtmp+eta));
	    zs1 = sqrt(1.e0/(dtmp+eta));
	    rho   = zs1*zs1;
	    //zs1 = sqrho * dkab * vdkcd[klps];
	    zs1 *= dkab * LDG(vdkcd[klps]);
	    //rho   = sqrho*sqrho;
            zs0 *= rho;
	    //zs0 = rho * PQ2;
//	    rz    = rho * zeta;
//	    re    = rho * eta;
//	    ze2   = re  * zeta2;
//	    ze22  = re  * zeta;
	    dtmp *= rho;
	    ze22  = dtmp * eta;
	    ze2   = ze22 * HALF;
#pragma unroll
	    for (int i=0; i<3; i++ ) {
//		WP[i] = rz*QP[i];
//		WQ[i] = rz*QP[i] - QP[i];
//		WQ[i] = (WP[i] = rz*QP[i]) - QP[i];
		WQ[i] = (WP[i] = dtmp*QP[i]) - QP[i];
	    }
	    //T     = rho * PQ2;
	    //cssss = sqrho * dk;
	    //cssss = sqrho * dkab * vdkcd[klps];
	    //gpu_fmt( ssss, 3, T, cssss );
	    //zs0 = rho * PQ2;
	    //zs1 = sqrho * dkab * vdkcd[klps];
	    //gpu_fmt( ssss, 3, zs0, zs1 );
#if   CUDA_FMT_M == 3
	    gpu_fmt3_method3( zs0, zs1, ssss );
#elif CUDA_FMT_M == 2
	    gpu_fmt3_method2( zs0, zs1, ssss );
#elif CUDA_FMT_M == 1
	    gpu_fmt3_method1( zs0, zs1, ssss );
#else
	    gpu_fmt3( ssss, zs0, zs1 );
#endif
            //zs0 = zeta2*(ssss[0] - rz*ssss[1]);
            //zs1 = zeta2*(ssss[1] - rz*ssss[2]);
            zs0 = zeta2*(ssss[0] - dtmp*ssss[1]);
            zs1 = zeta2*(ssss[1] - dtmp*ssss[2]);
            }
            /*
	    // psss
#pragma unroll
	    for ( m=0; m<=2; m++ ) {
		m1 = m+1;
#pragma unroll
		for ( i=0; i<3; i++ )
		    psss[m][i] = PA[i]*ssss[m]+WP[i]*ssss[m1];
	    }
            */
            /*
	    // dsss
#pragma unroll
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
            */

            /*
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
            */
            /*
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
            */
//#pragma unroll
            for (int j=0; j<3; j++)
            {
              double d0, d1;
              double p[3];
              int j2;
              //int jid[]={0,1,2,0,1,2};

#pragma unroll
              for (int i=0; i<3; i++)
                p[i] = PA[j]*ssss[i]+WP[j]*ssss[i+1];
#pragma unroll
              for (int i=0; i<3; i++)
                PSPS[j*3+i] += QC[i]*p[0] + WQ[i]*p[1];
              PSPS[j*3+j] += ze2*ssss[1];

              d0 = PA[j]*p[0] + WP[j]*p[1] + zs0;
              d1 = PA[j]*p[1] + WP[j]*p[2] + zs1;
#pragma unroll
              for (int i=0; i<3; i++)
                DSPS[j*3+i] += QC[i]*d0 + WQ[i]*d1;
	      DSPS[j*3+j] += ze22*p[1];

              j2 = (j+2)%3;
              //j2 = jid[j+2];
              d0 = PA[j2]*p[0] + WP[j2]*p[1];
              d1 = PA[j2]*p[1] + WP[j2]*p[2];
#pragma unroll
              for (int i=0; i<3; i++)
                DSPS[(j2+3)*3+i] += QC[i]*d0 + WQ[i]*d1;
	      DSPS[(j2+3)*3+j2] += ze2*p[1];

              j2 = (j+1)%3;
              //j2 = jid[j+1];
	      DSPS[(j+3)*3+j2] += ze2*p[1];
            }
	}	// klps
    }	// ijps
    // ppps
    /*
#pragma unroll
    for (int i=0; i<3; i++)
        PPPS[0*3*3+0*3+i] = DSPS[0*3+i] - BA[0]*PSPS[0*3+i];
#pragma unroll
    for (int i=0; i<3; i++)
        PPPS[0*3*3+1*3+i] = DSPS[3*3+i] - BA[1]*PSPS[0*3+i];
#pragma unroll
    for (int i=0; i<3; i++)
        PPPS[0*3*3+2*3+i] = DSPS[5*3+i] - BA[2]*PSPS[0*3+i];
#pragma unroll
    for (int i=0; i<3; i++)
        PPPS[1*3*3+0*3+i] = DSPS[3*3+i] - BA[0]*PSPS[1*3+i];
#pragma unroll
    for (int i=0; i<3; i++)
        PPPS[1*3*3+1*3+i] = DSPS[1*3+i] - BA[1]*PSPS[1*3+i];
#pragma unroll
    for (int i=0; i<3; i++)
        PPPS[1*3*3+2*3+i] = DSPS[4*3+i] - BA[2]*PSPS[1*3+i];
#pragma unroll
    for (int i=0; i<3; i++)
        PPPS[2*3*3+0*3+i] = DSPS[5*3+i] - BA[0]*PSPS[2*3+i];
#pragma unroll
    for (int i=0; i<3; i++)
        PPPS[2*3*3+1*3+i] = DSPS[4*3+i] - BA[1]*PSPS[2*3+i];
#pragma unroll
    for (int i=0; i<3; i++)
        PPPS[2*3*3+2*3+i] = DSPS[2*3+i] - BA[2]*PSPS[2*3+i];
    */
    {
    int ix3[9]={0,3,5,3,1,4,5,4,2};
#pragma unroll
    for (int k=0; k<3; k++)
#pragma unroll
    for (int j=0; j<3; j++)
    {
      int ix=ix3[k*3+j];
#pragma unroll
      for (int i=0; i<3; i++)
        PPPS[(k*3+j)*3+i] = DSPS[ix*3+i] - BA[j]*PSPS[k*3+i];
    }
    }
}

/** １つのCS４重対に対して(pp,pp)タイプの２電子積分を計算する関数
 * @ingroup core-twoint
 * */
__device__ void gpu_twoint_core_pppp_(
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double PPPP[3*3*3*3] ) {
//    int ijps, klps, i, m, m1, c0;
//    int ab, ab01, ab10, ab00;
//    double cssss, zeta, dkab, xiza, eta, xizc, dk;
//    double rz, re, ze2, ze22, zeta2, eta2;
//    double tmp, tmp3[3], tmp6[6];
//    double sqrho, rho, PQ2, T;
//    double PC[3], PA[3], WP[3], QP[3], QC[3], WQ[3];
//    double ssss[4+1], psss[3+1][3], dsss[2+1][6];
//    double ssps[3], psps[1+1][3*3], dsps[1+1][6*3];
    double PSPS[3*3], DSPS[6*3], PSDS[3*6], DSDS[6*6];
    double PPPS[3*3*3], PPDS[3*3*6];
#pragma unroll
    for (int i=0; i<3*3; i++ ) PSPS[i] = ZERO;
#pragma unroll
    for (int i=0; i<6*3; i++ ) DSPS[i] = ZERO;
#pragma unroll
    for (int i=0; i<3*6; i++ ) PSDS[i] = ZERO;
#pragma unroll
    for (int i=0; i<6*6; i++ ) DSDS[i] = ZERO;
    for (int ijps=0; ijps<(*nijps); ijps++ ) {
        double PC[3], PA[3];
        double zeta, dkab, xiza, zeta2;
	zeta = LDG(vzeta[ijps]);
	dkab = LDG(vdkab[ijps]);
	xiza = LDG(vxiza[ijps]);
	zeta2 = HALF * zeta;
#pragma unroll
	for (int i=0; i<3; i++ ) {
	    PC[i] = AC[i] + xiza*BA[i];
	    PA[i] = xiza * BA[i];
	}
	for (int klps=0; klps<(*nklps); klps++ ) {
            double ssss[4+1], psss[3+1][3];
            double psps[1+1][3*3];
            double WP[3], QP[3], QC[3], WQ[3];
            double cssss, eta, xizc, dk;
            double sqrho, rho, PQ2, T;
            double rz, re, ze2, ze22, eta2;
	    eta  = LDG(veta[klps]);
	    dk   = dkab * LDG(vdkcd[klps]);
	    xizc = LDG(vxizc[klps]);
	    eta2 = eta * HALF;
	    PQ2  = ZERO;
#pragma unroll
	    for (int i=0; i<3; i++ ) {
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
#pragma unroll
	    for (int i=0; i<3; i++ ) {
		WP[i] = rz*QP[i];
		WQ[i] = rz*QP[i] - QP[i];
	    }
	    T     = rho * PQ2;
	    cssss = sqrho * dk;
#if   CUDA_FMT_M == 3
	    gpu_fmt4_method3( T, cssss, ssss );
#elif CUDA_FMT_M == 2
	    gpu_fmt4_method2( T, cssss, ssss );
#elif CUDA_FMT_M == 1
	    gpu_fmt4_method1( T, cssss, ssss );
#else
	    gpu_fmt( ssss, 4, T, cssss );
#endif
	    // (P,S|S,S)
#pragma unroll
	    for (int m=0, m1=1; m<=3; m++, m1++) {
		psss[m][0] = PA[0]*ssss[m] + WP[0]*ssss[m1];
		psss[m][1] = PA[1]*ssss[m] + WP[1]*ssss[m1];
		psss[m][2] = PA[2]*ssss[m] + WP[2]*ssss[m1];
	    }	// end m loop
	    // (D,S|S,S)
            /*
#pragma unroll
	    for (int m=0, m1=1; m<=2; m++, m1++) {
		double tmp = zeta2*(ssss[m] - rz*ssss[m1]);
		dsss[m][0] = PA[0]*psss[m][0] + WP[0]*psss[m1][0]
			   + tmp;
		dsss[m][1] = PA[1]*psss[m][1] + WP[1]*psss[m1][1]
			   + tmp;
		dsss[m][2] = PA[2]*psss[m][2] + WP[2]*psss[m1][2]
			   + tmp;
		dsss[m][3] = PA[0]*psss[m][1] + WP[0]*psss[m1][1];
		dsss[m][4] = PA[1]*psss[m][2] + WP[1]*psss[m1][2];
		dsss[m][5] = PA[2]*psss[m][0] + WP[2]*psss[m1][0];
	    }	// end m loop
            */
/*
	    // (S,S|P,S)
	    ssps[0] = QC[0]*ssss[1] + WQ[0]*ssss[2];
	    ssps[1] = QC[1]*ssss[1] + WQ[1]*ssss[2];
	    ssps[2] = QC[2]*ssss[1] + WQ[2]*ssss[2];
*/
	    // (P,S|P,S)
            /*
#pragma unroll
	    for (int m=0, m1=1; m<=1; m++, m1++) {
              double tmp = ze2*ssss[m1];
		psps[m][0*3+0] = QC[0]*psss[m][0] + WQ[0]*psss[m1][0]
			       + tmp;
//			       + ze2*ssss[m1];
		psps[m][0*3+1] = QC[1]*psss[m][0] + WQ[1]*psss[m1][0];
		psps[m][0*3+2] = QC[2]*psss[m][0] + WQ[2]*psss[m1][0];
		psps[m][1*3+0] = QC[0]*psss[m][1] + WQ[0]*psss[m1][1];
		psps[m][1*3+1] = QC[1]*psss[m][1] + WQ[1]*psss[m1][1]
			       + tmp;
//			       + ze2*ssss[m1];
		psps[m][1*3+2] = QC[2]*psss[m][1] + WQ[2]*psss[m1][1];
		psps[m][2*3+0] = QC[0]*psss[m][2] + WQ[0]*psss[m1][2];
		psps[m][2*3+1] = QC[1]*psss[m][2] + WQ[1]*psss[m1][2];
		psps[m][2*3+2] = QC[2]*psss[m][2] + WQ[2]*psss[m1][2]
			       + tmp;
//			       + ze2*ssss[m1];
	    }	// end m loop
#pragma unroll
	    for (int i=0; i<3*3; i++) PSPS[i] += psps[0][i];
            */
#pragma unroll
            for (int i=0; i<3; i++) {
#pragma unroll
              for (int j=0; j<3; j++) {
                psps[0][i*3+j] = QC[j]*psss[0][i] + WQ[j]*psss[1][i];
                psps[1][i*3+j] = QC[j]*psss[1][i] + WQ[j]*psss[2][i];
              }
              psps[0][i*3+i] += ze2*ssss[1];
              psps[1][i*3+i] += ze2*ssss[2];
#pragma unroll
              for (int j=0; j<3; j++)
                PSPS[i*3+j] += psps[0][i*3+j];
            }
/*
	    // (D,S|P,S)
#pragma unroll
	    for (int m=0; m<=1; m++) {
		int m1 = m+1;
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
#pragma unroll
	    for (int i=0; i<6*3; i++) DSPS[i] += dsps[0][i];
*/
	    // (P,S|D,S)
/*
#pragma unroll
	    for (int i=0; i<3; i++) tmp3[i] = psss[0][i] - re*psss[1][i];

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
*/
#pragma unroll
            for (int i=0; i<3; i++) {
	      double tmp;
#pragma unroll
              for (int j=0; j<6; j++) {
                int j0 = j%3;
                int j1 = (j<3)? j: (j+1)%3;
                PSDS[i*6+j] += QC[j0]*psps[0][i*3+j1] + WQ[j0]*psps[1][i*3+j1];
              }
              tmp = eta2*(psss[0][i] - re*psss[1][i]);
#pragma unroll
              for (int j=0; j<3; j++)
                PSDS[i*6+j] += tmp;
              //PSDS[i*6+i  ] += ze2*ssps[i];
              //PSDS[i*6+i+3] += ze2*ssps[(i+1)%3];
            }
            {
              double ssps[3];
              ssps[0] = QC[0]*ssss[1] + WQ[0]*ssss[2];
              ssps[1] = QC[1]*ssss[1] + WQ[1]*ssss[2];
              ssps[2] = QC[2]*ssss[1] + WQ[2]*ssss[2];
#pragma unroll
              for (int i=0; i<3; i++) {
              PSDS[i*6+i  ] += ze2*ssps[i];
              PSDS[i*6+i+3] += ze2*ssps[(i+1)%3];
              }
            }
	    // (D,S|D,S)
/*
#pragma unroll
	    for (int i=0; i<6; i++) tmp6[i] = dsss[0][i] - re*dsss[1][i];

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
*/
	    // (D,S|P,S)
            /*
#pragma unroll
	    for (int m=0; m<=1; m++) {
		int m1 = m+1;
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
#pragma unroll
	    for (int i=0; i<6*3; i++) DSPS[i] += dsps[0][i];
            */
#pragma unroll
            for (int i=0; i<6; i++) {
              double dsps[3][1+1];
              double dsss[2+1];
              double zt;
              int i0 = i%3;
              int i1 = (i<3)? i: (i+1)%3;
              int i2 = (i+2)%3;

              zt = (i<3)? zeta2: 0;
              dsss[0] = PA[i0]*psss[0][i1] + WP[i0]*psss[1][i1] + zt*(ssss[0] - rz*ssss[1]);
              dsss[1] = PA[i0]*psss[1][i1] + WP[i0]*psss[2][i1] + zt*(ssss[1] - rz*ssss[2]);
              dsss[2] = PA[i0]*psss[2][i1] + WP[i0]*psss[3][i1] + zt*(ssss[2] - rz*ssss[3]);

#pragma unroll
              for (int j=0; j<3; j++) {
                dsps[j][0] = QC[j]*dsss[0] + WQ[j]*dsss[1];
                dsps[j][1] = QC[j]*dsss[1] + WQ[j]*dsss[2];
              }
              zt = (i<3)? ze22: ze2;
		dsps[i0][0] += zt*psss[1][i1];
		dsps[i0][1] += zt*psss[2][i1];
              if (i>=3) {
		dsps[i1][0] += zt*psss[1][i0];
		dsps[i1][1] += zt*psss[2][i0];
              }
#pragma unroll
              for (int j=0; j<3; j++) 
                DSPS[i*3+j] += dsps[j][0];
#pragma unroll
              for (int j=0; j<6; j++) {
                int j0 = j%3;
                int j1 = (j<3)? j: (j+1)%3;
                DSDS[i*6+j] += QC[j0]*dsps[j1][0] + WQ[j0]*dsps[j1][1];
              }
              {
                double tmp;
                tmp = eta2*(dsss[0] - re*dsss[1]);
#pragma unroll
                for (int j=0; j<3; j++) DSDS[i*6+j] += tmp;
              }
                DSDS[i*6+i0  ] += zt*psps[1][i1*3+i0];
              if (i<3) {
                //DSDS[i*6+i  ] += ze22*psps[1][i*3+i];
                DSDS[i*6+i+3] += ze22*psps[1][i*3+(i+1)%3];
              } else {
                //DSDS[i*6+i0  ] += ze2*psps[1][i1*3+i0];
                DSDS[i*6+i1  ] += ze2*psps[1][i0*3+i1];
                DSDS[i*6+i   ] += ze2*psps[1][i1*3+i1];
                DSDS[i*6+i1+3] += ze2*psps[1][i0*3+i2];
              }
            }
//#pragma unroll
//	    for (int i=0; i<3*3; i++) PSPS[i] += psps[0][i];
//#pragma unroll
//	    for (int i=0; i<6*3; i++) DSPS[i] += dsps[0][i];
	}	// klps
    }	// ijps
    // (P,P|P,S)
#pragma unroll
    for (int c0=0; c0<3; c0++)
        PPPS[0*3*3+0*3+c0] = DSPS[0*3+c0] - BA[0]*PSPS[0*3+c0];
#pragma unroll
    for (int c0=0; c0<3; c0++)
        PPPS[0*3*3+1*3+c0] = DSPS[3*3+c0] - BA[1]*PSPS[0*3+c0];
#pragma unroll
    for (int c0=0; c0<3; c0++)
        PPPS[0*3*3+2*3+c0] = DSPS[5*3+c0] - BA[2]*PSPS[0*3+c0];
#pragma unroll
    for (int c0=0; c0<3; c0++)
        PPPS[1*3*3+0*3+c0] = DSPS[3*3+c0] - BA[0]*PSPS[1*3+c0];
#pragma unroll
    for (int c0=0; c0<3; c0++)
        PPPS[1*3*3+1*3+c0] = DSPS[1*3+c0] - BA[1]*PSPS[1*3+c0];
#pragma unroll
    for (int c0=0; c0<3; c0++)
        PPPS[1*3*3+2*3+c0] = DSPS[4*3+c0] - BA[2]*PSPS[1*3+c0];
#pragma unroll
    for (int c0=0; c0<3; c0++)
        PPPS[2*3*3+0*3+c0] = DSPS[5*3+c0] - BA[0]*PSPS[2*3+c0];
#pragma unroll
    for (int c0=0; c0<3; c0++)
        PPPS[2*3*3+1*3+c0] = DSPS[4*3+c0] - BA[1]*PSPS[2*3+c0];
#pragma unroll
    for (int c0=0; c0<3; c0++)
        PPPS[2*3*3+2*3+c0] = DSPS[2*3+c0] - BA[2]*PSPS[2*3+c0];
    // (P,P|D,S)
#pragma unroll
    for (int c0=0; c0<6; c0++)
        PPDS[0*3*6+0*6+c0] = DSDS[0*6+c0] - BA[0]*PSDS[0*6+c0];
#pragma unroll
    for (int c0=0; c0<6; c0++)
        PPDS[0*3*6+1*6+c0] = DSDS[3*6+c0] - BA[1]*PSDS[0*6+c0];
#pragma unroll
    for (int c0=0; c0<6; c0++)
        PPDS[0*3*6+2*6+c0] = DSDS[5*6+c0] - BA[2]*PSDS[0*6+c0];
#pragma unroll
    for (int c0=0; c0<6; c0++)
        PPDS[1*3*6+0*6+c0] = DSDS[3*6+c0] - BA[0]*PSDS[1*6+c0];
#pragma unroll
    for (int c0=0; c0<6; c0++)
        PPDS[1*3*6+1*6+c0] = DSDS[1*6+c0] - BA[1]*PSDS[1*6+c0];
#pragma unroll
    for (int c0=0; c0<6; c0++)
        PPDS[1*3*6+2*6+c0] = DSDS[4*6+c0] - BA[2]*PSDS[1*6+c0];
#pragma unroll
    for (int c0=0; c0<6; c0++)
        PPDS[2*3*6+0*6+c0] = DSDS[5*6+c0] - BA[0]*PSDS[2*6+c0];
#pragma unroll
    for (int c0=0; c0<6; c0++)
        PPDS[2*3*6+1*6+c0] = DSDS[4*6+c0] - BA[1]*PSDS[2*6+c0];
#pragma unroll
    for (int c0=0; c0<6; c0++)
        PPDS[2*3*6+2*6+c0] = DSDS[2*6+c0] - BA[2]*PSDS[2*6+c0];
    // (a,b|c,d+1) = (a,b|c+1,d) - DC(a,b|c,d)
    // (P,P|P,P)
#pragma unroll
    for (int ab=0, ab01=0, ab10=0, ab00=0;
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
__device__ void gpu_twoint_core_dsss_(
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double DSSS[6] ) {
//    int ijps, klps, i, m;
//    double cssss, zeta, dkab, xiza, eta, xizc, dk;
//    double zeta, dkab, xiza;
    double zeta, dkab, xiza;
//    double rz, zeta2, tmp;
    double zeta2;
//    double PC[3], PA[3], WP[3], QP[3];
    double PC[3], PA[3];
//    double ssss[1+2], psss[1+1][3];
    double dsss[6];

#pragma unroll
    for (int i=0; i<6; i++ ) dsss[i] = ZERO;
//    for (int i=0; i<6; i++ ) DSSS[i] = ZERO;
    for (int ijps=0; ijps<(*nijps); ijps++ ) {
	zeta = LDG(vzeta[ijps]);
	dkab = LDG(vdkab[ijps]);
	xiza = LDG(vxiza[ijps]);
	zeta2 = HALF * zeta;
#pragma unroll
	for (int i=0; i<3; i++ ) {
//	    PC[i] = AC[i] + xiza*BA[i];
//	    PA[i] = xiza * BA[i];
	    PC[i] = AC[i] + (PA[i] = xiza * BA[i]);
	}
	for (int klps=0; klps<(*nklps); klps++ ) {
//          double WP[3],QP[3];
          double WP[3];
//          double ssss[1+2], psss[1+1][3];
          double ssss[1+2];
          double rz;
          {
          double cssss, eta, xizc, dk;
          double sqrho, rho, PQ2, T;
	    eta  = LDG(veta[klps]);
	    dk   = dkab * LDG(vdkcd[klps]);
	    xizc = LDG(vxizc[klps]);
	    PQ2  = ZERO;
#pragma unroll
	    for (int i=0; i<3; i++ ) {
		//QP[i] = xizc*DC[i] - PC[i];
		//PQ2  += QP[i]*QP[i];
		WP[i] = xizc*DC[i] - PC[i];
		PQ2  += WP[i]*WP[i];
	    }
	    sqrho = sqrt(1.e0/(zeta+eta));
	    rho   = sqrho*sqrho;
	    rz    = rho * zeta;
#pragma unroll
	    for (int i=0; i<3; i++ ) WP[i] *= rz;
//	    for (i=0; i<3; i++ ) WP[i]=rz*QP[i];
	    T     = rho * PQ2;
	    cssss = sqrho * dk;
	    //gpu_fmt( ssss, 2, T, cssss );
#if   CUDA_FMT_M == 3
	    gpu_fmt2_method3( T, cssss, ssss );
#elif CUDA_FMT_M == 2
	    gpu_fmt2_method2( T, cssss, ssss );
#elif CUDA_FMT_M == 1
	    gpu_fmt2_method1( T, cssss, ssss );
#else
	    gpu_fmt2( ssss, T, cssss );
#endif
          }
          /*
	    // psss
#pragma unroll
	    for (int m=0; m<=1; m++ ) {
#pragma unroll
		for (int i=0; i<3; i++ )
		    psss[m][i] = PA[i]*ssss[m]+WP[i]*ssss[m+1];
	    }
	    // dsss
            {
            double tmp;
	    tmp = ssss[0] - rz*ssss[1];
	    DSSS[0] += PA[0]*psss[0][0] + WP[0]*psss[1][0] + zeta2*tmp;
	    DSSS[1] += PA[1]*psss[0][1] + WP[1]*psss[1][1] + zeta2*tmp;
	    DSSS[2] += PA[2]*psss[0][2] + WP[2]*psss[1][2] + zeta2*tmp;
	    DSSS[3] += PA[0]*psss[0][1] + WP[0]*psss[1][1];
	    DSSS[5] += PA[1]*psss[0][2] + WP[1]*psss[1][2];
	    DSSS[4] += PA[2]*psss[0][0] + WP[2]*psss[1][0];
            }
          */
          {
	    double tmp = zeta2 * (ssss[0] - rz*ssss[1]);
            /*
#pragma unroll
          for (int i=0; i<3; i++)
          {
            double psss0,psss1;
            int i2 = (i+2)%3;

	    psss0 = PA[i]*ssss[0]+WP[i]*ssss[1];
            psss1 = PA[i]*ssss[1]+WP[i]*ssss[2];
//	    DSSS[i   ] += PA[i ]*psss0 + WP[i ]*psss1 + tmp;
//	    DSSS[i2+3] += PA[i2]*psss0 + WP[i2]*psss1;
	    dsss[i   ] += PA[i ]*psss0 + WP[i ]*psss1 + tmp;
	    dsss[i2+3] += PA[i2]*psss0 + WP[i2]*psss1;

          }
          */
          {
            double psss0,psss1;
	    psss0 = PA[0]*ssss[0]+WP[0]*ssss[1];
            psss1 = PA[0]*ssss[1]+WP[0]*ssss[2];
	    dsss[0] += PA[0]*psss0 + WP[0]*psss1 + tmp;
	    dsss[4] += PA[2]*psss0 + WP[2]*psss1;
	    //DSSS[0] += PA[0]*psss0 + WP[0]*psss1 + tmp;
	    //DSSS[4] += PA[2]*psss0 + WP[2]*psss1;

	    psss0 = PA[1]*ssss[0]+WP[1]*ssss[1];
            psss1 = PA[1]*ssss[1]+WP[1]*ssss[2];
	    dsss[1] += PA[1]*psss0 + WP[1]*psss1 + tmp;
	    dsss[3] += PA[0]*psss0 + WP[0]*psss1;
	    //DSSS[1] += PA[1]*psss0 + WP[1]*psss1 + tmp;
	    //DSSS[3] += PA[0]*psss0 + WP[0]*psss1;

	    psss0 = PA[2]*ssss[0]+WP[2]*ssss[1];
            psss1 = PA[2]*ssss[1]+WP[2]*ssss[2];
	    dsss[2] += PA[2]*psss0 + WP[2]*psss1 + tmp;
	    dsss[5] += PA[1]*psss0 + WP[1]*psss1;
	    //DSSS[2] += PA[2]*psss0 + WP[2]*psss1 + tmp;
	    //DSSS[5] += PA[1]*psss0 + WP[1]*psss1;
          }
          }
	}	// klps
    }	// ijps
    {
    double sqr3 = sqrt(3.e0);
#pragma unroll
//    for (int i=3; i<6; i++ ) DSSS[i] *= sqr3;
    for (int i=0; i<3; i++ ) DSSS[i] = dsss[i];
#pragma unroll
    for (int i=3; i<6; i++ ) DSSS[i] = dsss[i] * sqr3;
    }
}

#if 0
/** １つのCS４重対に対して(ds,ps)タイプの２電子積分を計算する関数
 * @ingroup core-twoint
 * */
__device__ void gpu_twoint_core_dsps_(
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double DSPS[6*3] ) {
//    int ijps, klps, i, m, m1;
    double cssss, zeta, dkab, xiza, eta, xizc, dk;
    double rz, re, ze2, ze22, zeta2, tmp;
    double sqrho, rho, PQ2, T;
    double PC[3], PA[3], WP[3], QP[3], QC[3], WQ[3];
    double ssss[3+1], psss[2+1][3], dsss[1+1][6];
    double sqr3;
    sqr3 = sqrt(3.e0);
#pragma unroll
    for (int i=0; i<6*3; i++ ) DSPS[i] = ZERO;
    for (int ijps=0; ijps<(*nijps); ijps++ ) {
	zeta = vzeta[ijps];
	dkab = vdkab[ijps];
	xiza = vxiza[ijps];
	zeta2 = HALF * zeta;
#pragma unroll
	for (int i=0; i<3; i++ ) {
	    PC[i] = AC[i] + xiza*BA[i];
	    PA[i] = xiza * BA[i];
	}
	for (int klps=0; klps<(*nklps); klps++ ) {
	    eta  = veta[klps];
	    dk   = dkab * vdkcd[klps];
	    xizc = vxizc[klps];
	    PQ2  = ZERO;
#pragma unroll
	    for (int i=0; i<3; i++ ) {
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
#pragma unroll
	    for (int i=0; i<3; i++ ) {
		WP[i] = rz*QP[i];
		WQ[i] = rz*QP[i] - QP[i];
	    }
	    T     = rho * PQ2;
	    cssss = sqrho * dk;
	    //gpu_fmt( ssss, 3, T, cssss );
#if   CUDA_FMT_M == 3
	    gpu_fmt3_method3( T, cssss, ssss );
#elif CUDA_FMT_M == 2
	    gpu_fmt3_method2( T, cssss, ssss );
#elif CUDA_FMT_M == 1
	    gpu_fmt3_method1( T, cssss, ssss );
#else
	    gpu_fmt3( ssss, T, cssss );
#endif
/*
	    // psss
#pragma unroll
	    for (int m=0; m<=2; m++ ) {
		m1 = m+1;
#pragma unroll
		for ( i=0; i<3; i++ )
		    psss[m][i] = PA[i]*ssss[m]+WP[i]*ssss[m1];
	    }
	    // dsss
#pragma unroll
	    for (int m=0; m<=1; m++) {
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
*/
            {
              double z0, z1;
              rz = -rz;
              z0 = zeta2*(ssss[0] + rz*ssss[1]);
              z1 = zeta2*(ssss[1] + rz*ssss[2]);
#pragma unroll
            for (int j=0; j<3; j++)
            {
              double dsss0, dsss1;
              double psss[3];
              int j2,j3;
	    // dsps
#pragma unroll
            for (int i=0; i<3; i++) 
              psss[i] = PA[j]*ssss[i]+WP[j]*ssss[i+1];
            dsss0 = PA[j]*psss[0] + WP[j]*psss[1] + z0;
            dsss1 = PA[j]*psss[1] + WP[j]*psss[2] + z1;
#pragma unroll
            for (int i=0; i<3; i++) 
              DSPS[j*3+i] += QC[i]*dsss0 + WQ[i]*dsss1;
	    DSPS[j*3+j] += ze22*psss[1];

            j2 = (j+2)%3;
            j3 = (j2+3)*3;
            dsss0 = PA[j2]*psss[0] + WP[j2]*psss[1];
            dsss1 = PA[j2]*psss[1] + WP[j2]*psss[2];
#pragma unroll
            for (int i=0; i<3; i++) 
              DSPS[j3+i] += QC[i]*dsss0 + WQ[i]*dsss1;
	    DSPS[j3+j2] += ze2*psss[1];

            j2 = (j+1)%3;
	    DSPS[(j+3)*3+j2] += ze2*psss[1];
            }
            }
	}	// klps
    }	// ijps
#pragma unroll
    for (int  i=3*3; i<6*3; i++) DSPS[i] *= sqr3;
}

/** １つのCS４重対に対して(ds,pp)タイプの２電子積分を計算する関数
 * @ingroup core-twoint
 * */
__device__ void gpu_twoint_core_dspp_(
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double DSPP[6*3*3] ) {
//    int ijps, klps, i, m, m1;
//    int ab, ab01, ab10, ab00;
//    double cssss, zeta, dkab, xiza, eta, xizc, dk;
//    double rz, re, ze2, ze22, zeta2, eta2;
//    double tmp, tmp6[6];
//    double sqrho, rho, PQ2, T;
//    double PC[3], PA[3], WP[3], QP[3], QC[3], WQ[3];
//    double ssss[4+1], psss[3+1][3], dsss[2+1][6];
//    double psps[3*3], dsps[1+1][6*3];
    double DSPS[6*3], DSDS[6*6];
#pragma unroll
    for (int i=0; i<6*3; i++ ) DSPS[i] = ZERO;
#pragma unroll
    for (int i=0; i<6*6; i++ ) DSDS[i] = ZERO;
    for (int ijps=0; ijps<(*nijps); ijps++ ) {
      double zeta, dkab, xiza, zeta2;
      double PC[3], PA[3];
	zeta = vzeta[ijps];
	dkab = vdkab[ijps];
	xiza = vxiza[ijps];
	zeta2 = HALF * zeta;
#pragma unroll
	for (int i=0; i<3; i++ ) {
	    PC[i] = AC[i] + xiza*BA[i];
	    PA[i] = xiza * BA[i];
	}
	for (int klps=0; klps<(*nklps); klps++ ) {
          double cssss, eta, xizc, dk;
          double rz, re, ze2, ze22, eta2;
          double sqrho, rho, PQ2, T;
          double ssss[4+1];
          double WP[3], QP[3], QC[3], WQ[3];
	    eta  = veta[klps];
	    dk   = dkab * vdkcd[klps];
	    xizc = vxizc[klps];
	    eta2 = HALF * eta;
	    PQ2  = ZERO;
#pragma unroll
	    for (int i=0; i<3; i++ ) {
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
#pragma unroll
	    for (int i=0; i<3; i++ ) {
		WP[i] = rz*QP[i];
		WQ[i] = rz*QP[i] - QP[i];
	    }
	    T     = rho * PQ2;
	    cssss = sqrho * dk;
	    //gpu_fmt( ssss, 4, T, cssss );
#if   CUDA_FMT_M == 3
	    gpu_fmt4_method3( T, cssss, ssss );
#elif CUDA_FMT_M == 2
	    gpu_fmt4_method2( T, cssss, ssss );
#elif CUDA_FMT_M == 1
	    gpu_fmt4_method1( T, cssss, ssss );
#else
	    gpu_fmt3( ssss, T, cssss );
#endif
            {
            //  double psss[3+1][3], dsss[2+1][6];
              double psss[3+1][3];
            //  double psps[3*3], dsps[1+1][6*3];
	    // (P,S|S,S)
#pragma unroll
	    for (int m=0; m<=3; m++) {
		int m1 = m+1;
		psss[m][0] = PA[0]*ssss[m] + WP[0]*ssss[m1];
		psss[m][1] = PA[1]*ssss[m] + WP[1]*ssss[m1];
		psss[m][2] = PA[2]*ssss[m] + WP[2]*ssss[m1];
	    }	// end m loop
	    // (D,S|S,S)
            /*
#pragma unroll
	    for (int m=0; m<=2; m++) {
              double tmp;
		int m1 = m+1;
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
#pragma unroll
	    for (int i=0; i<6; i++) {
              double tmp = dsss[0][i] - re*dsss[1][i];
	      DSDS[i*6+0] += eta2*tmp;
	      DSDS[i*6+1] += eta2*tmp;
	      DSDS[i*6+2] += eta2*tmp;
            }
            */
	    // (P,S|P,S)
            /*
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
            */
	    // (D,S|P,S)
            /*
#pragma unroll
	    for (int m=0; m<=1; m++) {
		int m1 = m+1;
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
            */

            /*
#pragma unroll
	    for (int i=0; i<6*3; i++ ) DSPS[i] += dsps[0][i];
            */
	    // (D,S|D,S)
            /*
#pragma unroll
//	    for (int i=0; i<6; i++) tmp6[i] = dsss[0][i] - re*dsss[1][i];
	    for (int i=0; i<6; i++) {
              double tmp = dsss[0][i] - re*dsss[1][i];
	      DSDS[i*6+0] += eta2*tmp;
	      DSDS[i*6+1] += eta2*tmp;
	      DSDS[i*6+2] += eta2*tmp;
            }
            */

            /*
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
            */
            /*
            {
            double dsps[6*3][1+1];
            int i3x[]={0,1,2,1,2,0};
            for (int j=0; j<6; j++) {
#pragma unroll
              for (int i=0; i<3; i++) {
                dsps[j*3+i][0] = QC[i]*dsss[0][j] + WQ[i]*dsss[1][j];
                dsps[j*3+i][1] = QC[i]*dsss[1][j] + WQ[i]*dsss[2][j];
              }
              {
                double zz=(j<3)? ze22: ze2;
                int j0 = (j<3)? j: j%3;
                int j1 = (j<3)? j: (j+1)%3;
		dsps[j*3+j0][0] += zz*psss[1][j1];
		dsps[j*3+j0][1] += zz*psss[2][j1];
                if (j>=3) {
		dsps[j*3+j1][0] += zz*psss[1][j0];
		dsps[j*3+j1][1] += zz*psss[2][j0];
                }
              }
              for (int i=0; i<3; i++) 
	        DSPS[j*3+i] += dsps[j*3+i][0];
#pragma unroll
              for (int i=0; i<6; i++) { 
                //DSDS[j*6+i] += QC[i%3]*dsps[0][j*3+i3x[i]] + WQ[i%3]*dsps[1][j*3+i3x[i]];
                DSDS[j*6+i] += QC[i%3]*dsps[j*3+i3x[i]][0] + WQ[i%3]*dsps[j*3+i3x[i]][1];
              }
            } 
            } //dsps
            */

            {
              double dsss[2+1];
              double dsps[3][1+1];
#pragma unroll
            for (int j=0; j<6; j++) {
              int j0 = j%3;
              int j1 = (j<3)? j: (j+1)%3;
              double zt = (j<3)? zeta2: 0;
#pragma unroll
	    for (int m=0; m<=2; m++) {
              double tmp;
		int m1 = m+1;
		tmp = ssss[m] - rz*ssss[m1];

		dsss[m] = PA[j0]*psss[m][j1] + WP[j0]*psss[m1][j1] + zt*tmp;
	    }	// end m loop
            {
              double tmp = dsss[0] - re*dsss[1];
	      DSDS[j*6+0] += eta2*tmp;
	      DSDS[j*6+1] += eta2*tmp;
	      DSDS[j*6+2] += eta2*tmp;
            }

#pragma unroll
              for (int i=0; i<3; i++) {
                dsps[i][0] = QC[i]*dsss[0] + WQ[i]*dsss[1];
                dsps[i][1] = QC[i]*dsss[1] + WQ[i]*dsss[2];
              }
                zt=(j<3)? ze22: ze2;
		dsps[j0][0] += zt*psss[1][j1];
		dsps[j0][1] += zt*psss[2][j1];
                if (j>=3) {
		dsps[j1][0] += zt*psss[1][j0];
		dsps[j1][1] += zt*psss[2][j0];
                }
#pragma unroll
              for (int i=0; i<3; i++) 
	        DSPS[j*3+i] += dsps[i][0];
#pragma unroll
              for (int i=0; i<6; i++) { 
                //DSDS[j*6+i] += QC[i%3]*dsps[i3x[i]][0] + WQ[i%3]*dsps[i3x[i]][1];
                int i1 = (i<3)? i: (i+1)%3;
                DSDS[j*6+i] += QC[i%3]*dsps[i1][0] + WQ[i%3]*dsps[i1][1];
              }
            }
            } //dsss, dsps

            // Should be parallelize with 6 threads
#pragma unroll
            for (int j=0; j<3; j++) {
              double psps[3];
              int j1=(j+1)%3;
              int j2=(j+2)%3;
#pragma unroll
              for (int i=0; i<3; i++) {
	        psps[i] = QC[i]*psss[1][j] + WQ[i]*psss[2][j];
              }
	      psps[j] += ze2*ssss[2];

              DSDS[j*6+j]   += ze22*psps[j];
              DSDS[j*6+j+3] += ze22*psps[j1];
	      DSDS[(j2+3)*6+j2] += ze2*psps[j2];
	      DSDS[(j2+3)*6+j2+3] += ze2*psps[j];
	      DSDS[(j+3)*6+j1] += ze2*psps[j1];
	      DSDS[(j+3)*6+j1+3] += ze2*psps[j2];
            }
            /*
	    DSDS[0*6+0] += ze22*psps[0*3+0];
	    DSDS[0*6+3] += ze22*psps[0*3+1];
	    DSDS[1*6+1] += ze22*psps[1*3+1];
	    DSDS[1*6+4] += ze22*psps[1*3+2];
	    DSDS[2*6+2] += ze22*psps[2*3+2];
	    DSDS[2*6+5] += ze22*psps[2*3+0];

	    DSDS[3*6+0] += ze2*psps[1*3+0];
	    DSDS[3*6+3] += ze2*psps[1*3+1];
	    DSDS[4*6+1] += ze2*psps[2*3+1];
	    DSDS[4*6+4] += ze2*psps[2*3+2];
	    DSDS[5*6+0] += ze2*psps[2*3+0];
	    DSDS[5*6+3] += ze2*psps[2*3+1];

	    DSDS[3*6+1] += ze2*psps[0*3+1];
	    DSDS[3*6+4] += ze2*psps[0*3+2];
	    DSDS[4*6+2] += ze2*psps[1*3+2];
	    DSDS[4*6+5] += ze2*psps[1*3+0];
	    DSDS[5*6+2] += ze2*psps[0*3+2];
	    DSDS[5*6+5] += ze2*psps[0*3+0];
            */

            }
	}	// klps
    }	// ijps
    // (D,S|P,P)
#pragma unroll
    for (int ab=0, ab01=0, ab10=0, ab00=0;
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
    {
    double sqr3;
    sqr3 = sqrt(3.e0);
#pragma unroll
    for (int i=3*3*3; i<6*3*3; i++ ) DSPP[i] *= sqr3;
    }
}

#endif /* 0 */
