// dspp
#include "cuda-twoint-core-os.h"

__device__ void gpu_hrr_clear_os_dspp( double *eh ) {
    int i;
    // (DS|PS)
#pragma unroll
    for ( i=0; i<(0+18); i++ ) eh[i] = 0.e0;
    // (DS|DS)
#pragma unroll
    for ( i=18; i<(18+36); i++ ) eh[i] = 0.e0;
}

__device__ void gpu_hrr_calc_os_dspp(
//        const double BA[3], const double DC[3], double *eh ) {
        const double BA[3], const double DC[3], const double *ev, double *eh ) {
    // HRR for (XX|XX)-type integral (center CD)
    // (DS,PP)
#pragma unroll
    for (int i=0,j=0,ix=18,iy=0;i<6;i++) {
      eh[iy++] = ev[ix+0] - DC[0]*ev[j];
      eh[iy++] = ev[ix+3] - DC[1]*ev[j];
      eh[iy++] = ev[ix+4] - DC[2]*ev[j++];
      eh[iy++] = ev[ix+3] - DC[0]*ev[j];
      eh[iy++] = ev[ix+1] - DC[1]*ev[j];
      eh[iy++] = ev[ix+5] - DC[2]*ev[j++];
      eh[iy++] = ev[ix+4] - DC[0]*ev[j];
      eh[iy++] = ev[ix+5] - DC[1]*ev[j];
      eh[iy++] = ev[ix+2] - DC[2]*ev[j++];
      ix+=6;
    }
/*
#pragma unroll
    for (int i=0,iy=0;i<18;i++) {
#pragma unroll
      for (int j=0;j<3;j++) {
        eh[iy++] = -DC[j]*eh[i];
      }
    }

#pragma unroll
    for (int i=0,ix=18,iy=0;i<6;i++) {
      eh[iy++] += eh[ix+0];
      eh[iy++] += eh[ix+3];
      eh[iy++] += eh[ix+4];
      eh[iy++] += eh[ix+3];
      eh[iy++] += eh[ix+1];
      eh[iy++] += eh[ix+5];
      eh[iy++] += eh[ix+4];
      eh[iy++] += eh[ix+5];
      eh[iy++] += eh[ix+2];
      ix+=6;
    }
    */
    /*
    eh[  54] = eh[  18] - DC[0]*eh[   0];
    eh[  55] = eh[  21] - DC[1]*eh[   0];
    eh[  56] = eh[  22] - DC[2]*eh[   0];
    eh[  57] = eh[  21] - DC[0]*eh[   1];
    eh[  58] = eh[  19] - DC[1]*eh[   1];
    eh[  59] = eh[  23] - DC[2]*eh[   1];
    eh[  60] = eh[  22] - DC[0]*eh[   2];
    eh[  61] = eh[  23] - DC[1]*eh[   2];
    eh[  62] = eh[  20] - DC[2]*eh[   2];
    eh[  63] = eh[  24] - DC[0]*eh[   3];
    eh[  64] = eh[  27] - DC[1]*eh[   3];
    eh[  65] = eh[  28] - DC[2]*eh[   3];
    eh[  66] = eh[  27] - DC[0]*eh[   4];
    eh[  67] = eh[  25] - DC[1]*eh[   4];
    eh[  68] = eh[  29] - DC[2]*eh[   4];
    eh[  69] = eh[  28] - DC[0]*eh[   5];
    eh[  70] = eh[  29] - DC[1]*eh[   5];
    eh[  71] = eh[  26] - DC[2]*eh[   5];
    eh[  72] = eh[  30] - DC[0]*eh[   6];
    eh[  73] = eh[  33] - DC[1]*eh[   6];
    eh[  74] = eh[  34] - DC[2]*eh[   6];
    eh[  75] = eh[  33] - DC[0]*eh[   7];
    eh[  76] = eh[  31] - DC[1]*eh[   7];
    eh[  77] = eh[  35] - DC[2]*eh[   7];
    eh[  78] = eh[  34] - DC[0]*eh[   8];
    eh[  79] = eh[  35] - DC[1]*eh[   8];
    eh[  80] = eh[  32] - DC[2]*eh[   8];
    eh[  81] = eh[  36] - DC[0]*eh[   9];
    eh[  82] = eh[  39] - DC[1]*eh[   9];
    eh[  83] = eh[  40] - DC[2]*eh[   9];
    eh[  84] = eh[  39] - DC[0]*eh[  10];
    eh[  85] = eh[  37] - DC[1]*eh[  10];
    eh[  86] = eh[  41] - DC[2]*eh[  10];
    eh[  87] = eh[  40] - DC[0]*eh[  11];
    eh[  88] = eh[  41] - DC[1]*eh[  11];
    eh[  89] = eh[  38] - DC[2]*eh[  11];
    eh[  90] = eh[  42] - DC[0]*eh[  12];
    eh[  91] = eh[  45] - DC[1]*eh[  12];
    eh[  92] = eh[  46] - DC[2]*eh[  12];
    eh[  93] = eh[  45] - DC[0]*eh[  13];
    eh[  94] = eh[  43] - DC[1]*eh[  13];
    eh[  95] = eh[  47] - DC[2]*eh[  13];
    eh[  96] = eh[  46] - DC[0]*eh[  14];
    eh[  97] = eh[  47] - DC[1]*eh[  14];
    eh[  98] = eh[  44] - DC[2]*eh[  14];
    eh[  99] = eh[  48] - DC[0]*eh[  15];
    eh[ 100] = eh[  51] - DC[1]*eh[  15];
    eh[ 101] = eh[  52] - DC[2]*eh[  15];
    eh[ 102] = eh[  51] - DC[0]*eh[  16];
    eh[ 103] = eh[  49] - DC[1]*eh[  16];
    eh[ 104] = eh[  53] - DC[2]*eh[  16];
    eh[ 105] = eh[  52] - DC[0]*eh[  17];
    eh[ 106] = eh[  53] - DC[1]*eh[  17];
    eh[ 107] = eh[  50] - DC[2]*eh[  17];
    */
}

__device__ void gpu_hrr_coef_os_dspp(
        //const double *eh, double *DINT ) {
        const double *DC, const double *eh, double *DINT ) {
    int i, j, k, l, iao, jao, kao, lao, ix, iy;
    double coef_a, coef_ab, coef_abc;
    ix = 54;
    iy = 0;

    /*
#pragma unroll
    for ( i=0, iao=4; i<6; i++, iao++ ) {
        coef_a = LDG(DFACT[iao]);
        coef_ab = coef_a;
#pragma unroll
        for ( k=0, kao=1; k<3; k++, kao++ ) {
            coef_abc = coef_ab * LDG(DFACT[kao]);
#pragma unroll
            for ( l=0, lao=1; l<3; l++, lao++ ) {
                //DINT[iy] = coef_abc * DFACT[lao] * eh[ix];
                DINT[iy] *= coef_abc * LDG(DFACT[lao]);
                iy++;
                ix++;
            }
        }
    }
    */
    const int ix9[]={0,3,4,3,1,5,4,5,2};
    int iz;
#pragma unroll
    for ( i=0, iao=4, ix=18, iz=0; i<6; i++, iao++ ) {
        coef_a = LDG(DFACT[iao]);
        coef_ab = coef_a;
#pragma unroll
        for ( k=0, kao=1, j=0; k<3; k++, kao++, iz++ ) {
            coef_abc = coef_ab * LDG(DFACT[kao]);
#pragma unroll
            for ( l=0, lao=1; l<3; l++, lao++ ) {
                DINT[iy++] = coef_abc * LDG(DFACT[lao]) *
                  (eh[ix+ix9[j++]] - DC[l]*eh[iz]);
            }
        }
        ix+=6;
    }
}

__device__ void gpu_vrr_calc_os_dspp(
        const double T, const double cssss,
        const double zeta2, const double eta2, const double ze2,
        const double rz, const double re,
        const double PA[3], const double WP[3],
        const double QC[3], const double WQ[3],
        double *eh ) {
    double ev[27], ew[36], d[6];
    // (ss|ss) m=0,4
    //gpu_fmt( &ev[0], 4, T, cssss );
#if   CUDA_FMT_M == 3
    gpu_fmt4_method3( T, cssss, &ev[0] );
#elif CUDA_FMT_M == 2
    gpu_fmt4_method2( T, cssss, &ev[0] );
#elif CUDA_FMT_M == 1
    gpu_fmt4_method1( T, cssss, &ev[0] );
#else
    gpu_fmt4( &ev[0], T, cssss );
#endif
    // (ps|ss) m=0,3
    ew[ 0]=PA[0]*ev[0]+WP[0]*ev[1];
    ew[ 1]=PA[1]*ev[0]+WP[1]*ev[1];
    ew[ 2]=PA[2]*ev[0]+WP[2]*ev[1];
    ew[ 3]=PA[0]*ev[1]+WP[0]*ev[2];
    ew[ 4]=PA[1]*ev[1]+WP[1]*ev[2];
    ew[ 5]=PA[2]*ev[1]+WP[2]*ev[2];
    ew[ 6]=PA[0]*ev[2]+WP[0]*ev[3];
    ew[ 7]=PA[1]*ev[2]+WP[1]*ev[3];
    ew[ 8]=PA[2]*ev[2]+WP[2]*ev[3];
    ew[ 9]=PA[0]*ev[3]+WP[0]*ev[4];
    ew[10]=PA[1]*ev[3]+WP[1]*ev[4];
    ew[11]=PA[2]*ev[3]+WP[2]*ev[4];

    d[0] = zeta2*(ev[0]-rz*ev[1]);
    d[1] = zeta2*(ev[1]-rz*ev[2]);
    d[2] = zeta2*(ev[2]-rz*ev[3]);
    d[3] = ze2*ev[2];
    // (ds|ss) m=0,2
    ev[ 0]=PA[0]*ew[ 0]+WP[0]*ew[ 3]+d[0];
    ev[ 1]=PA[1]*ew[ 1]+WP[1]*ew[ 4]+d[0];
    ev[ 2]=PA[2]*ew[ 2]+WP[2]*ew[ 5]+d[0];
    ev[ 3]=PA[0]*ew[ 1]+WP[0]*ew[ 4];
    ev[ 4]=PA[0]*ew[ 2]+WP[0]*ew[ 5];
    ev[ 5]=PA[1]*ew[ 2]+WP[1]*ew[ 5];
    ev[ 6]=PA[0]*ew[ 3]+WP[0]*ew[ 6]+d[1];
    ev[ 7]=PA[1]*ew[ 4]+WP[1]*ew[ 7]+d[1];
    ev[ 8]=PA[2]*ew[ 5]+WP[2]*ew[ 8]+d[1];
    ev[ 9]=PA[0]*ew[ 4]+WP[0]*ew[ 7];
    ev[10]=PA[0]*ew[ 5]+WP[0]*ew[ 8];
    ev[11]=PA[1]*ew[ 5]+WP[1]*ew[ 8];
    ev[12]=PA[0]*ew[ 6]+WP[0]*ew[ 9]+d[2];
    ev[13]=PA[1]*ew[ 7]+WP[1]*ew[10]+d[2];
    ev[14]=PA[2]*ew[ 8]+WP[2]*ew[11]+d[2];
    ev[15]=PA[0]*ew[ 7]+WP[0]*ew[10];
    ev[16]=PA[0]*ew[ 8]+WP[0]*ew[11];
    ev[17]=PA[1]*ew[ 8]+WP[1]*ew[11];
    // (ps|ps) m=[1,1]
    ev[18]=QC[0]*ew[ 3]+WQ[0]*ew[ 6]+d[3];
    ev[19]=QC[1]*ew[ 3]+WQ[1]*ew[ 6];
    ev[20]=QC[2]*ew[ 3]+WQ[2]*ew[ 6];
    ev[21]=QC[0]*ew[ 4]+WQ[0]*ew[ 7];
    ev[22]=QC[1]*ew[ 4]+WQ[1]*ew[ 7]+d[3];
    ev[23]=QC[2]*ew[ 4]+WQ[2]*ew[ 7];
    ev[24]=QC[0]*ew[ 5]+WQ[0]*ew[ 8];
    ev[25]=QC[1]*ew[ 5]+WQ[1]*ew[ 8];
    ev[26]=QC[2]*ew[ 5]+WQ[2]*ew[ 8]+d[3];
#pragma unroll
    for (int i=0,j=3; i<6; i++) d[i] = ze2*ew[j++];
    // (ds|ps) m=[0,1]
    ew[ 0]=QC[0]*ev[ 0]+WQ[0]*ev[ 6]+2.e0*d[ 0];
    ew[ 1]=QC[1]*ev[ 0]+WQ[1]*ev[ 6];
    ew[ 2]=QC[2]*ev[ 0]+WQ[2]*ev[ 6];
    ew[ 3]=QC[0]*ev[ 1]+WQ[0]*ev[ 7];
    ew[ 4]=QC[1]*ev[ 1]+WQ[1]*ev[ 7]+2.e0*d[ 1];
    ew[ 5]=QC[2]*ev[ 1]+WQ[2]*ev[ 7];
    ew[ 6]=QC[0]*ev[ 2]+WQ[0]*ev[ 8];
    ew[ 7]=QC[1]*ev[ 2]+WQ[1]*ev[ 8];
    ew[ 8]=QC[2]*ev[ 2]+WQ[2]*ev[ 8]+2.e0*d[ 2];
    ew[ 9]=QC[0]*ev[ 3]+WQ[0]*ev[ 9]+     d[ 1];
    ew[10]=QC[1]*ev[ 3]+WQ[1]*ev[ 9]+     d[ 0];
    ew[11]=QC[2]*ev[ 3]+WQ[2]*ev[ 9];
    ew[12]=QC[0]*ev[ 4]+WQ[0]*ev[10]+     d[ 2];
    ew[13]=QC[1]*ev[ 4]+WQ[1]*ev[10];
    ew[14]=QC[2]*ev[ 4]+WQ[2]*ev[10]+     d[ 0];
    ew[15]=QC[0]*ev[ 5]+WQ[0]*ev[11];
    ew[16]=QC[1]*ev[ 5]+WQ[1]*ev[11]+     d[ 2];
    ew[17]=QC[2]*ev[ 5]+WQ[2]*ev[11]+     d[ 1];
    ew[18]=QC[0]*ev[ 6]+WQ[0]*ev[12]+2.e0*d[ 3];
    ew[19]=QC[1]*ev[ 6]+WQ[1]*ev[12];
    ew[20]=QC[2]*ev[ 6]+WQ[2]*ev[12];
    ew[21]=QC[0]*ev[ 7]+WQ[0]*ev[13];
    ew[22]=QC[1]*ev[ 7]+WQ[1]*ev[13]+2.e0*d[ 4];
    ew[23]=QC[2]*ev[ 7]+WQ[2]*ev[13];
    ew[24]=QC[0]*ev[ 8]+WQ[0]*ev[14];
    ew[25]=QC[1]*ev[ 8]+WQ[1]*ev[14];
    ew[26]=QC[2]*ev[ 8]+WQ[2]*ev[14]+2.e0*d[ 5];
    ew[27]=QC[0]*ev[ 9]+WQ[0]*ev[15]+     d[ 4];
    ew[28]=QC[1]*ev[ 9]+WQ[1]*ev[15]+     d[ 3];
    ew[29]=QC[2]*ev[ 9]+WQ[2]*ev[15];
    ew[30]=QC[0]*ev[10]+WQ[0]*ev[16]+     d[ 5];
    ew[31]=QC[1]*ev[10]+WQ[1]*ev[16];
    ew[32]=QC[2]*ev[10]+WQ[2]*ev[16]+     d[ 3];
    ew[33]=QC[0]*ev[11]+WQ[0]*ev[17];
    ew[34]=QC[1]*ev[11]+WQ[1]*ev[17]+     d[ 5];
    ew[35]=QC[2]*ev[11]+WQ[2]*ev[17]+     d[ 4];
#pragma unroll
    for (int i=0; i<6; i++) d[i] = eta2*(ev[ i]-re*ev[i+6]);
#pragma unroll
    for (int i=0; i<18; i++ ) eh[i]+=ew[i];
    // (DS|DS)
    // (ds|ds) m=[0,0]
    eh[ 18]+=QC[0]*ew[ 0]+WQ[0]*ew[18]+d[0] +2.e0*ze2*ev[18];
    eh[ 19]+=QC[1]*ew[ 1]+WQ[1]*ew[19]+d[0];
    eh[ 20]+=QC[2]*ew[ 2]+WQ[2]*ew[20]+d[0];
    eh[ 21]+=QC[0]*ew[ 1]+WQ[0]*ew[19]+2.e0*ze2*ev[19];
    eh[ 22]+=QC[0]*ew[ 2]+WQ[0]*ew[20]+2.e0*ze2*ev[20];
    eh[ 23]+=QC[1]*ew[ 2]+WQ[1]*ew[20];
    eh[ 24]+=QC[0]*ew[ 3]+WQ[0]*ew[21]+d[1];
    eh[ 25]+=QC[1]*ew[ 4]+WQ[1]*ew[22]+d[1] +2.e0*ze2*ev[22];
    eh[ 26]+=QC[2]*ew[ 5]+WQ[2]*ew[23]+d[1];
    eh[ 27]+=QC[0]*ew[ 4]+WQ[0]*ew[22];
    eh[ 28]+=QC[0]*ew[ 5]+WQ[0]*ew[23];
    eh[ 29]+=QC[1]*ew[ 5]+WQ[1]*ew[23]+2.e0*ze2*ev[23];
    eh[ 30]+=QC[0]*ew[ 6]+WQ[0]*ew[24]+d[2];
    eh[ 31]+=QC[1]*ew[ 7]+WQ[1]*ew[25]+d[2];
    eh[ 32]+=QC[2]*ew[ 8]+WQ[2]*ew[26]+d[2] +2.e0*ze2*ev[26];
    eh[ 33]+=QC[0]*ew[ 7]+WQ[0]*ew[25];
    eh[ 34]+=QC[0]*ew[ 8]+WQ[0]*ew[26];
    eh[ 35]+=QC[1]*ew[ 8]+WQ[1]*ew[26];
    eh[ 36]+=QC[0]*ew[ 9]+WQ[0]*ew[27]+d[3] +     ze2*ev[21];
    eh[ 37]+=QC[1]*ew[10]+WQ[1]*ew[28]+d[3] +     ze2*ev[19];
    eh[ 38]+=QC[2]*ew[11]+WQ[2]*ew[29]+d[3];
    eh[ 39]+=QC[0]*ew[10]+WQ[0]*ew[28]+     ze2*ev[22];
    eh[ 40]+=QC[0]*ew[11]+WQ[0]*ew[29]+     ze2*ev[23];
    eh[ 41]+=QC[1]*ew[11]+WQ[1]*ew[29]+     ze2*ev[20];
    eh[ 42]+=QC[0]*ew[12]+WQ[0]*ew[30]+d[4] +     ze2*ev[24];
    eh[ 43]+=QC[1]*ew[13]+WQ[1]*ew[31]+d[4];
    eh[ 44]+=QC[2]*ew[14]+WQ[2]*ew[32]+d[4] +     ze2*ev[20];
    eh[ 45]+=QC[0]*ew[13]+WQ[0]*ew[31]+     ze2*ev[25];
    eh[ 46]+=QC[0]*ew[14]+WQ[0]*ew[32]+     ze2*ev[26];
    eh[ 47]+=QC[1]*ew[14]+WQ[1]*ew[32];
    eh[ 48]+=QC[0]*ew[15]+WQ[0]*ew[33]+d[5];
    eh[ 49]+=QC[1]*ew[16]+WQ[1]*ew[34]+d[5] +     ze2*ev[25];
    eh[ 50]+=QC[2]*ew[17]+WQ[2]*ew[35]+d[5] +     ze2*ev[23];
    eh[ 51]+=QC[0]*ew[16]+WQ[0]*ew[34];
    eh[ 52]+=QC[0]*ew[17]+WQ[0]*ew[35];
    eh[ 53]+=QC[1]*ew[17]+WQ[1]*ew[35]+     ze2*ev[26];
}

__device__ void gpu_vrr_cint_os_dspp( const double *ev, double *eh ) {
    int La=2, Lb=0, Lc=1, Ld=1;
    int i, ih, iv;
    // (DS|PS)
#pragma unroll
    for ( i=0, iv=44, ih=0; i<18; i++, iv++, ih++ ) eh[ih]+=ev[iv];
    // (DS|DS)
#pragma unroll
    for ( i=0, iv=80, ih=18; i<36; i++, iv++, ih++ ) eh[ih]+=ev[iv];
}

__device__ void gpu_twoint_core_os_dspp(
//        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *nijps, const double vzeta[], const double vdkab[],
        const double vxiza[], const double BA[3],
        const int *nklps, const double veta[], const double vdkcd[],
        const double vxizc[], const double DC[3], const double AC[3],
        double *DINT ) {
    int ijps, klps, i;
    double cssss, zeta, dkab, xiza, eta, xizc, dk, T;
    double zeta2, eta2, ze2, rz, re, PA[3], WP[3], QC[3], WQ[3];
    double PQ2, sqrho, rho, PC[3], QP[3];
//    double ev[116], eh[108];
    double eh[54];
//    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;

    gpu_hrr_clear_os_dspp( eh );
    for ( ijps=0; ijps<(*nijps); ijps++ ) {
        zeta  = LDG(vzeta[ijps]);
        dkab  = LDG(vdkab[ijps]);
        xiza  = LDG(vxiza[ijps]);
        zeta2 = HALF * zeta;
#pragma unroll
        for ( i=0; i<3; i++ ) {
            PC[i] = AC[i] + xiza*BA[i];
            PA[i] = xiza * BA[i];
        }
        for ( klps=0; klps<(*nklps); klps++ ) {
            eta  = LDG(veta[klps]);
            dk   = dkab * LDG(vdkcd[klps]);
            xizc = LDG(vxizc[klps]);
            eta2 = HALF * eta;
            PQ2  = ZERO;
#pragma unroll
            for ( i=0; i<3; i++ ) {
                QC[i] = xizc*DC[i];
                QP[i] = xizc*DC[i] - PC[i];
                PQ2  += QP[i]*QP[i];
            }
            sqrho = sqrt(1.e0/(zeta+eta));
            rho   = sqrho*sqrho;
            rz    = rho * zeta;
            re    = rho * eta;
            ze2 = rz * eta2;
#pragma unroll
            for ( i=0; i<3; i++ ) {
                WP[i] = rz*QP[i];
                WQ[i] = rz*QP[i] - QP[i];
            }
            T     = rho * PQ2;
            cssss = sqrho * dk;
            gpu_vrr_calc_os_dspp(
                    T, cssss, zeta2, eta2, ze2, rz, re, PA, WP, QC, WQ,
                    //ev );
                    eh );
            //gpu_vrr_cint_os_dspp( ev, eh );
        }	// for (klps)
    }	// for (ijps)
    //gpu_hrr_calc_os_dspp( BA, DC, eh );
    //gpu_hrr_calc_os_dspp( BA, DC, eh, DINT );
    //gpu_hrr_coef_os_dspp( eh, DINT );
    gpu_hrr_coef_os_dspp( DC, eh, DINT );
}

#if 0
int gpu_twoint_os_dspp(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int shel_atm[], const int shel_ini[],
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for partially direct SCF
        const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
        double ebuf_val[], short int ebuf_ind4[],
        int *last_ijcs, int *last_klcs ) {
    int Lab, Lcd, i, j, k, l, ipat, ix;
    int I2, IJ, K2, KL;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1, max_klcs;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe, coe0;
    double DINTEG[6*1*3*3];
    long nzeri, max_nzeri, nzeri4;
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    long ebuf_max_nzeri = *pebuf_max_nzeri;
    int mythread;

    mythread = omp_get_thread_num();
    if ( DFACT == NULL ) DFACT = gpu_getadd_dfact();

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs0 = leading_cs_pair[Lab];
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];
    nzeri     = *ebuf_non_zero_eri;
    max_nzeri = ebuf_max_nzeri - 6*1*3*3;
    nzeri4    = nzeri*4;
    if ( nzeri >= max_nzeri ) {
        *last_ijcs = ijcs0+workerid;
        *last_klcs = klcs0 - 1;
        *ebuf_non_zero_eri = nzeri;
        return OFMO_EBUF_FULL;
    }

    for ( ijcs=ijcs0+workerid; ijcs<ijcs1; ijcs+=nworkers ) {
        val_ab = csp_schwarz[ijcs];
        ics    = csp_ics[ijcs];
        jcs    = csp_jcs[ijcs];
        ijps0  = csp_leading_ps_pair[ijcs];
        nijps  = csp_leading_ps_pair[ijcs+1]-ijps0;
        iat    = shel_atm[ics];
        jat    = shel_atm[jcs];
        iao0   = shel_ini[ics];
        jao0   = shel_ini[jcs];
        A[0]=atom_x[iat]; A[1]=atom_y[iat]; A[2]=atom_z[iat];
        B[0]=atom_x[jat]; B[1]=atom_y[jat]; B[2]=atom_z[jat];
        for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];
        max_klcs = ( Lab == Lcd ? ijcs+1 : klcs1 );
        for ( klcs=klcs0; klcs<max_klcs; klcs++ ) {
            val_cd = csp_schwarz[klcs];
            if ( val_ab*val_cd < EPS_PS4 ) continue;
            kcs    = csp_ics[klcs];
            lcs    = csp_jcs[klcs];
            klps0  = csp_leading_ps_pair[klcs];
            nklps  = csp_leading_ps_pair[klcs+1]-klps0;
            kat    = shel_atm[kcs];
            lat    = shel_atm[lcs];
            kao0   = shel_ini[kcs];
            lao0   = shel_ini[lcs];
            C[0]=atom_x[kat]; C[1]=atom_y[kat]; C[2]=atom_z[kat];
            D[0]=atom_x[lat]; D[1]=atom_y[lat]; D[2]=atom_z[lat];
            for ( i=0; i<3; i++ ) {
                AC[i] = A[i] - C[i];
                DC[i] = D[i] - C[i];
            }
            gpu_twoint_core_os_dspp(
                    &La, &Lb, &Lc, &Ld,
                    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
                    &psp_xiza[ijps0], BA,
                    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
                    &psp_xiza[klps0], DC,   AC,      DINTEG );
            ipat=((Lab != Lcd)||(ics==kcs && jcs>lcs) ? true : false );
            for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
                I2 = (iao*iao+iao)>>1;
                for ( j=0, jao=jao0; j<1; j++, jao++ ) {
                    if ( jao>iao ) { ix+=3*3; continue; }
                    IJ = I2 + jao;
                    coe0 = ( iao==jao ? HALF : ONE );
                    for ( k=0, kao=kao0; k<3; k++, kao++ ) {
                        K2 = (kao*kao+kao)>>1;
                        for ( l=0, lao=lao0; l<3; l++, lao++, ix++ ) {
                            if ( lao>kao ) continue;
                            if ( fabs(DINTEG[ix]) > EPS_ERI ) {
                                KL = K2 + lao;
                                if ( IJ >= KL ) {
                                    coe = coe0;
                                    if ( kao==lao ) coe *= HALF;
                                    if ( KL == IJ ) coe *= HALF;
                                    ebuf_val[nzeri]     = coe*DINTEG[ix];
                                    ebuf_ind4[nzeri4+0] = (short int)iao;
                                    ebuf_ind4[nzeri4+1] = (short int)jao;
                                    ebuf_ind4[nzeri4+2] = (short int)kao;
                                    ebuf_ind4[nzeri4+3] = (short int)lao;
                                    nzeri++;
                                    nzeri4+=4;
                                } else if ( ipat ) {
                                    coe = coe0;
                                    if ( kao==lao ) coe*=HALF;
                                    ebuf_val[nzeri]     = coe*DINTEG[ix];
                                    ebuf_ind4[nzeri4+0] = (short int)kao;
                                    ebuf_ind4[nzeri4+1] = (short int)lao;
                                    ebuf_ind4[nzeri4+2] = (short int)iao;
                                    ebuf_ind4[nzeri4+3] = (short int)jao;
                                    nzeri++;
                                    nzeri4+=4;
                                }
                            }
                        }	// l
                    }	// k
                }	// j
            }	// i
            if ( nzeri >= max_nzeri ) {
                *last_ijcs = ijcs;
                *last_klcs = klcs;
                *ebuf_non_zero_eri = nzeri;
                return OFMO_EBUF_FULL;
            }
        }	// for (klcs)
    }	// for (ijcs)
    *ebuf_non_zero_eri = nzeri;
    return OFMO_EBUF_NOFULL;
}

int gpu_twoint_direct_os_dspp(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int shel_atm[], const int shel_ini[],
        const double atom_x[], const double atom_y[],
        const double atom_z[], const int leading_cs_pair[],
        const double csp_schwarz[],
        const int csp_ics[], const int csp_jcs[],
        const int csp_leading_ps_pair[],
        const double psp_zeta[], const double psp_dkps[],
        const double psp_xiza[],
        // for direct SCF
        const long *petmp_max_nzeri, long *petmp_non_zero_eri,
        double etmp_val[], short int etmp_ind4[],
        const int *plast_ijcs, const int *plast_klcs,
        // density matrix & G-matrix data
        const int *pnao, const double Ds[], double G[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int last_ijcs=*plast_ijcs, last_klcs=*plast_klcs, nao=*pnao;
    long max_nzeri=*petmp_max_nzeri;
    long nzeri4, nzeri=*petmp_non_zero_eri;
    //
    int Lab, Lcd, i, j, k, l, ipat, ix;
    int I2, IJ, K2, KL;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1, max_klcs;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe, coe0;
    double DINTEG[6*1*3*3];
    int mythread;

    mythread = omp_get_thread_num();
    if ( DFACT == NULL ) DFACT = gpu_getadd_dfact();

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];
    if ( last_ijcs != -1 ) {
        ijcs = last_ijcs;
        klcs = last_klcs+1;
    } else {
        ijcs = leading_cs_pair[Lab] + workerid;
        klcs = klcs0;
    }

    max_nzeri -= 6*1*3*3;
    nzeri4    = nzeri*4;
    if ( nzeri >= max_nzeri ) {
        gpu_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4, Ds, G );
        nzeri = nzeri4 = 0;
    }

    for ( ; ijcs<ijcs1; ijcs+=nworkers ) {
        val_ab = csp_schwarz[ijcs];
        ics    = csp_ics[ijcs];
        jcs    = csp_jcs[ijcs];
        ijps0  = csp_leading_ps_pair[ijcs];
        nijps  = csp_leading_ps_pair[ijcs+1]-ijps0;
        iat    = shel_atm[ics];
        jat    = shel_atm[jcs];
        iao0   = shel_ini[ics];
        jao0   = shel_ini[jcs];
        A[0]=atom_x[iat]; A[1]=atom_y[iat]; A[2]=atom_z[iat];
        B[0]=atom_x[jat]; B[1]=atom_y[jat]; B[2]=atom_z[jat];
        for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];
        max_klcs = ( Lab == Lcd ? ijcs+1 : klcs1 );
        for ( ; klcs<max_klcs; klcs++ ) {
            val_cd = csp_schwarz[klcs];
            if ( val_ab*val_cd < EPS_PS4 ) continue;
            kcs    = csp_ics[klcs];
            lcs    = csp_jcs[klcs];
            klps0  = csp_leading_ps_pair[klcs];
            nklps  = csp_leading_ps_pair[klcs+1]-klps0;
            kat    = shel_atm[kcs];
            lat    = shel_atm[lcs];
            kao0   = shel_ini[kcs];
            lao0   = shel_ini[lcs];
            C[0]=atom_x[kat]; C[1]=atom_y[kat]; C[2]=atom_z[kat];
            D[0]=atom_x[lat]; D[1]=atom_y[lat]; D[2]=atom_z[lat];
            for ( i=0; i<3; i++ ) {
                AC[i] = A[i] - C[i];
                DC[i] = D[i] - C[i];
            }
            gpu_twoint_core_os_dspp(
                    &La, &Lb, &Lc, &Ld,
                    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
                    &psp_xiza[ijps0], BA,
                    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
                    &psp_xiza[klps0], DC,   AC,      DINTEG );
            ipat=((Lab != Lcd)||(ics==kcs && jcs>lcs) ? true : false );
            for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
                I2 = (iao*iao+iao)>>1;
                for ( j=0, jao=jao0; j<1; j++, jao++ ) {
                    if ( jao>iao ) { ix+=3*3; continue; }
                    IJ = I2 + jao;
                    coe0 = ( iao==jao ? HALF : ONE );
                    for ( k=0, kao=kao0; k<3; k++, kao++ ) {
                        K2 = (kao*kao+kao)>>1;
                        for ( l=0, lao=lao0; l<3; l++, lao++, ix++ ) {
                            if ( lao>kao ) continue;
                            if ( fabs(DINTEG[ix]) > EPS_ERI ) {
                                KL = K2 + lao;
                                if ( IJ >= KL ) {
                                    coe = coe0;
                                    if ( kao==lao ) coe *= HALF;
                                    if ( KL == IJ ) coe *= HALF;
                                    etmp_val[nzeri]     = coe*DINTEG[ix];
                                    etmp_ind4[nzeri4+0] = (short int)iao;
                                    etmp_ind4[nzeri4+1] = (short int)jao;
                                    etmp_ind4[nzeri4+2] = (short int)kao;
                                    etmp_ind4[nzeri4+3] = (short int)lao;
                                    nzeri++;
                                    nzeri4+=4;
                                } else if ( ipat ) {
                                    coe = coe0;
                                    if ( kao==lao ) coe*=HALF;
                                    etmp_val[nzeri]     = coe*DINTEG[ix];
                                    etmp_ind4[nzeri4+0] = (short int)kao;
                                    etmp_ind4[nzeri4+1] = (short int)lao;
                                    etmp_ind4[nzeri4+2] = (short int)iao;
                                    etmp_ind4[nzeri4+3] = (short int)jao;
                                    nzeri++;
                                    nzeri4+=4;
                                }
                            }
                        }	// l
                    }	// k
                }	// j
            }	// i
            if ( nzeri >= max_nzeri ) {
                gpu_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4,
                        Ds, G );
                nzeri = nzeri4= 0;
            }
        }	// for (klcs)
        klcs = klcs0;
    }	// for (ijcs)
    *petmp_non_zero_eri = nzeri;
    return 0;
}
#endif // 0
