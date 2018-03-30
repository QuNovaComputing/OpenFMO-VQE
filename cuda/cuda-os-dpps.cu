// dpps
#include "cuda-twoint-core-os.h"

__device__ void gpu_hrr_clear_dpps( double *eh, double *ef ) {
    int i;
    // (DS|PS)
#pragma unroll
    for ( i=0; i<(0+18); i++ ) eh[i] = 0.e0;
    // (FS|PS)
#pragma unroll
    //for ( i=18; i<(18+30); i++ ) eh[i] = 0.e0;
    for ( i=0; i<(30); i++ ) ef[i] = 0.e0;
}

__device__ void gpu_hrr_calc_dpps(
        //const double BA[3], const double DC[3], double *eh ) {
        const double BA[3], const double DC[3],
        const double eh[18], const double ef[30], double *ew ) {
    // (DP,PS)
    ew[  0] = ef[   0] - BA[0]*eh[   0];
    ew[  1] = ef[   1] - BA[0]*eh[   1];
    ew[  2] = ef[   2] - BA[0]*eh[   2];
    ew[  3] = ef[   9] - BA[1]*eh[   0];
    ew[  4] = ef[  10] - BA[1]*eh[   1];
    ew[  5] = ef[  11] - BA[1]*eh[   2];
    ew[  6] = ef[  12] - BA[2]*eh[   0];
    ew[  7] = ef[  13] - BA[2]*eh[   1];
    ew[  8] = ef[  14] - BA[2]*eh[   2];
    ew[  9] = ef[  15] - BA[0]*eh[   3];
    ew[ 10] = ef[  16] - BA[0]*eh[   4];
    ew[ 11] = ef[  17] - BA[0]*eh[   5];
    ew[ 12] = ef[   3] - BA[1]*eh[   3];
    ew[ 13] = ef[   4] - BA[1]*eh[   4];
    ew[ 14] = ef[   5] - BA[1]*eh[   5];
    ew[ 15] = ef[  24] - BA[2]*eh[   3];
    ew[ 16] = ef[  25] - BA[2]*eh[   4];
    ew[ 17] = ef[  26] - BA[2]*eh[   5];
    ew[ 18] = ef[  18] - BA[0]*eh[   6];
    ew[ 19] = ef[  19] - BA[0]*eh[   7];
    ew[ 20] = ef[  20] - BA[0]*eh[   8];
    ew[ 21] = ef[  27] - BA[1]*eh[   6];
    ew[ 22] = ef[  28] - BA[1]*eh[   7];
    ew[ 23] = ef[  29] - BA[1]*eh[   8];
    ew[ 24] = ef[   6] - BA[2]*eh[   6];
    ew[ 25] = ef[   7] - BA[2]*eh[   7];
    ew[ 26] = ef[   8] - BA[2]*eh[   8];
    ew[ 27] = ef[   9] - BA[0]*eh[   9];
    ew[ 28] = ef[  10] - BA[0]*eh[  10];
    ew[ 29] = ef[  11] - BA[0]*eh[  11];
    ew[ 30] = ef[  15] - BA[1]*eh[   9];
    ew[ 31] = ef[  16] - BA[1]*eh[  10];
    ew[ 32] = ef[  17] - BA[1]*eh[  11];
    ew[ 33] = ef[  21] - BA[2]*eh[   9];
    ew[ 34] = ef[  22] - BA[2]*eh[  10];
    ew[ 35] = ef[  23] - BA[2]*eh[  11];
    ew[ 36] = ef[  12] - BA[0]*eh[  12];
    ew[ 37] = ef[  13] - BA[0]*eh[  13];
    ew[ 38] = ef[  14] - BA[0]*eh[  14];
    ew[ 39] = ef[  21] - BA[1]*eh[  12];
    ew[ 40] = ef[  22] - BA[1]*eh[  13];
    ew[ 41] = ef[  23] - BA[1]*eh[  14];
    ew[ 42] = ef[  18] - BA[2]*eh[  12];
    ew[ 43] = ef[  19] - BA[2]*eh[  13];
    ew[ 44] = ef[  20] - BA[2]*eh[  14];
    ew[ 45] = ef[  21] - BA[0]*eh[  15];
    ew[ 46] = ef[  22] - BA[0]*eh[  16];
    ew[ 47] = ef[  23] - BA[0]*eh[  17];
    ew[ 48] = ef[  24] - BA[1]*eh[  15];
    ew[ 49] = ef[  25] - BA[1]*eh[  16];
    ew[ 50] = ef[  26] - BA[1]*eh[  17];
    ew[ 51] = ef[  27] - BA[2]*eh[  15];
    ew[ 52] = ef[  28] - BA[2]*eh[  16];
    ew[ 53] = ef[  29] - BA[2]*eh[  17];
    // HRR for (XX|XX)-type integral (center CD)
}

__device__ void gpu_hrr_coef_dpps(
        //const double *eh, double *DINT ) {
        double *DINT ) {
    int i, j, k, l, iao, jao, kao, lao, ix, iy;
    double coef_a, coef_ab, coef_abc;
    ix = 48;
    iy = 0;

#pragma unroll
    for ( i=0, iao=4; i<6; i++, iao++ ) {
        coef_a = LDG(DFACT[iao]);
#pragma unroll
        for ( j=0, jao=1; j<3; j++, jao++ ) {
            coef_ab = coef_a * LDG(DFACT[jao]);
#pragma unroll
            for ( k=0, kao=1; k<3; k++, kao++ ) {
                coef_abc = coef_ab * LDG(DFACT[kao]);
                //DINT[iy] = coef_abc * eh[ix];
                DINT[iy] *= coef_abc;
                iy++;
                ix++;
            }
        }
    }
}

__device__ void gpu_vrr_calc_dpps(
        const double T, const double cssss,
        const double zeta2, const double eta2, const double ze2,
        const double rz, const double re,
        const double PA[3], const double WP[3],
        const double QC[3], const double WQ[3],
        double *eh, double *ef ) {
    double d[6], ew[20], ev[16];
    // (ss|ss) m=0,4
    //gpu_fmt( &ev[0], 4, T, cssss );
#if   CUDA_FMT_M == 3
    gpu_fmt4_method3( T, cssss, &d[0] );
#elif CUDA_FMT_M == 2
    gpu_fmt4_method2( T, cssss, &d[0] );
#elif CUDA_FMT_M == 1
    gpu_fmt4_method1( T, cssss, &d[0] );
#else
    gpu_fmt4( &d[0], T, cssss );
#endif
    // (ps|ss) m=0,3
#pragma unroll
    for (int i=0; i<3; i++) {
#pragma unroll
      for (int j=0; j<4; j++) {
        ew[i*4+j] = PA[i]*d[j]+WP[i]*d[j+1];
      }
    }
    d[0]=zeta2*(d[0]-rz*d[1]);
    d[1]=zeta2*(d[1]-rz*d[2]);
    d[2]=zeta2*(d[2]-rz*d[3]);
    // (ds|ss) m=0,2
    ev[ 0]=PA[0]*ew[ 0]+WP[0]*ew[ 1]+d[0];
    ev[ 1]=PA[0]*ew[ 1]+WP[0]*ew[ 2]+d[1];
    ev[ 2]=PA[0]*ew[ 2]+WP[0]*ew[ 3]+d[2];
    ev[ 3]=PA[1]*ew[ 4]+WP[1]*ew[ 5]+d[0];
    ev[ 4]=PA[1]*ew[ 5]+WP[1]*ew[ 6]+d[1];
    ev[ 5]=PA[1]*ew[ 6]+WP[1]*ew[ 7]+d[2];
    ev[ 6]=PA[2]*ew[ 8]+WP[2]*ew[ 9]+d[0];
    ev[ 7]=PA[2]*ew[ 9]+WP[2]*ew[10]+d[1];
    ev[ 8]=PA[2]*ew[10]+WP[2]*ew[11]+d[2];
    ev[ 9]=PA[1]*ew[ 8]+WP[1]*ew[ 9];
    ev[10]=PA[1]*ew[ 9]+WP[1]*ew[10];
    ev[11]=PA[1]*ew[10]+WP[1]*ew[11];
    ev[12]=PA[0]*ew[ 4]+WP[0]*ew[ 5];
    ev[13]=PA[0]*ew[ 5]+WP[0]*ew[ 6];
    ev[14]=PA[0]*ew[ 8]+WP[0]*ew[ 9];
    ev[15]=PA[0]*ew[ 9]+WP[0]*ew[10];
    //ev[32]=PA[0]*ew[ 7]+WP[0]*ew[10];
    //ev[33]=PA[0]*ew[ 8]+WP[0]*ew[11];
    d[0] = ze2*ew[1];
    d[1] = ze2*ew[5];
    d[2] = ze2*ew[9];
    // (ds|ps) m=[0,0]
    eh[ 0]+=QC[0]*ev[ 0]+WQ[0]*ev[ 1]+2.e0*d[0];
    eh[ 1]+=QC[1]*ev[ 0]+WQ[1]*ev[ 1];
    eh[ 2]+=QC[2]*ev[ 0]+WQ[2]*ev[ 1];
    eh[ 3]+=QC[0]*ev[ 3]+WQ[0]*ev[ 4];
    eh[ 4]+=QC[1]*ev[ 3]+WQ[1]*ev[ 4]+2.e0*d[1];
    eh[ 5]+=QC[2]*ev[ 3]+WQ[2]*ev[ 4];
    eh[ 6]+=QC[0]*ev[ 6]+WQ[0]*ev[ 7];
    eh[ 7]+=QC[1]*ev[ 6]+WQ[1]*ev[ 7];
    eh[ 8]+=QC[2]*ev[ 6]+WQ[2]*ev[ 7]+2.e0*d[2];
    eh[ 9]+=QC[0]*ev[12]+WQ[0]*ev[13]+     d[1];
    eh[10]+=QC[1]*ev[12]+WQ[1]*ev[13]+     d[0];
    eh[11]+=QC[2]*ev[12]+WQ[2]*ev[13];
    eh[12]+=QC[0]*ev[14]+WQ[0]*ev[15]+     d[2];
    eh[13]+=QC[1]*ev[14]+WQ[1]*ev[15];
    eh[14]+=QC[2]*ev[14]+WQ[2]*ev[15]+     d[0];
    eh[15]+=QC[0]*ev[ 9]+WQ[0]*ev[10];
    eh[16]+=QC[1]*ev[ 9]+WQ[1]*ev[10]+     d[2];
    eh[17]+=QC[2]*ev[ 9]+WQ[2]*ev[10]+     d[1];
    // (fs|ss) m=0,1
    d[0] = 2.e0*zeta2*(ew[ 0]-rz*ew[ 1]);
    d[1] = 2.e0*zeta2*(ew[ 1]-rz*ew[ 2]);
    d[2] = 2.e0*zeta2*(ew[ 4]-rz*ew[ 5]);
    d[3] = 2.e0*zeta2*(ew[ 5]-rz*ew[ 6]);
    d[4] = 2.e0*zeta2*(ew[ 8]-rz*ew[ 9]);
    d[5] = 2.e0*zeta2*(ew[ 9]-rz*ew[10]);
    ew[ 0]=PA[0]*ev[ 0]+WP[0]*ev[ 1]+d[0];
    ew[ 1]=PA[0]*ev[ 1]+WP[0]*ev[ 2]+d[1];
    ew[ 2]=PA[1]*ev[ 3]+WP[1]*ev[ 4]+d[2];
    ew[ 3]=PA[1]*ev[ 4]+WP[1]*ev[ 5]+d[3];
    ew[ 4]=PA[2]*ev[ 6]+WP[2]*ev[ 7]+d[4];
    ew[ 5]=PA[2]*ev[ 7]+WP[2]*ev[ 8]+d[5];
    ew[ 6]=PA[1]*ev[ 0]+WP[1]*ev[ 1];
    ew[ 7]=PA[1]*ev[ 1]+WP[1]*ev[ 2];
    ew[ 8]=PA[2]*ev[ 0]+WP[2]*ev[ 1];
    ew[ 9]=PA[2]*ev[ 1]+WP[2]*ev[ 2];
    ew[10]=PA[0]*ev[ 3]+WP[0]*ev[ 4];
    ew[11]=PA[0]*ev[ 4]+WP[0]*ev[ 5];
    ew[12]=PA[0]*ev[ 6]+WP[0]*ev[ 7];
    ew[13]=PA[0]*ev[ 7]+WP[0]*ev[ 8];
    ew[14]=PA[0]*ev[ 9]+WP[0]*ev[10];
    ew[15]=PA[0]*ev[10]+WP[0]*ev[11];
    ew[16]=PA[2]*ev[ 3]+WP[2]*ev[ 4];
    ew[17]=PA[2]*ev[ 4]+WP[2]*ev[ 5];
    ew[18]=PA[1]*ev[ 6]+WP[1]*ev[ 7];
    ew[19]=PA[1]*ev[ 7]+WP[1]*ev[ 8];
    // (fs|ps) m+=[0,0]
    ef[ 0]+=QC[0]*ew[ 0]+WQ[0]*ew[ 1]+3.e0*ze2*ev[ 1];
    ef[ 1]+=QC[1]*ew[ 0]+WQ[1]*ew[ 1];
    ef[ 2]+=QC[2]*ew[ 0]+WQ[2]*ew[ 1];
    ef[ 3]+=QC[0]*ew[ 2]+WQ[0]*ew[ 3];
    ef[ 4]+=QC[1]*ew[ 2]+WQ[1]*ew[ 3]+3.e0*ze2*ev[ 4];
    ef[ 5]+=QC[2]*ew[ 2]+WQ[2]*ew[ 3];
    ef[ 6]+=QC[0]*ew[ 4]+WQ[0]*ew[ 5];
    ef[ 7]+=QC[1]*ew[ 4]+WQ[1]*ew[ 5];
    ef[ 8]+=QC[2]*ew[ 4]+WQ[2]*ew[ 5]+3.e0*ze2*ev[ 7];
    ef[ 9]+=QC[0]*ew[ 6]+WQ[0]*ew[ 7]+2.e0*ze2*ev[13];
    ef[10]+=QC[1]*ew[ 6]+WQ[1]*ew[ 7]+     ze2*ev[ 1];
    ef[11]+=QC[2]*ew[ 6]+WQ[2]*ew[ 7];
    ef[12]+=QC[0]*ew[ 8]+WQ[0]*ew[ 9]+2.e0*ze2*ev[15];
    ef[13]+=QC[1]*ew[ 8]+WQ[1]*ew[ 9];
    ef[14]+=QC[2]*ew[ 8]+WQ[2]*ew[ 9]+     ze2*ev[ 1];
    ef[15]+=QC[0]*ew[10]+WQ[0]*ew[11]+     ze2*ev[ 4];
    ef[16]+=QC[1]*ew[10]+WQ[1]*ew[11]+2.e0*ze2*ev[13];
    ef[17]+=QC[2]*ew[10]+WQ[2]*ew[11];
    ef[18]+=QC[0]*ew[12]+WQ[0]*ew[13]+     ze2*ev[ 7];
    ef[19]+=QC[1]*ew[12]+WQ[1]*ew[13];
    ef[20]+=QC[2]*ew[12]+WQ[2]*ew[13]+2.e0*ze2*ev[15];
    ef[21]+=QC[0]*ew[14]+WQ[0]*ew[15]+     ze2*ev[10];
    ef[22]+=QC[1]*ew[14]+WQ[1]*ew[15]+     ze2*ev[15];
    ef[23]+=QC[2]*ew[14]+WQ[2]*ew[15]+     ze2*ev[13];
    ef[24]+=QC[0]*ew[16]+WQ[0]*ew[17];
    ef[25]+=QC[1]*ew[16]+WQ[1]*ew[17]+2.e0*ze2*ev[10];
    ef[26]+=QC[2]*ew[16]+WQ[2]*ew[17]+     ze2*ev[ 4];
    ef[27]+=QC[0]*ew[18]+WQ[0]*ew[19];
    ef[28]+=QC[1]*ew[18]+WQ[1]*ew[19]+     ze2*ev[ 7];
    ef[29]+=QC[2]*ew[18]+WQ[2]*ew[19]+2.e0*ze2*ev[10];
}

__device__ void gpu_vrr_cint_dpps( const double *ev, double *eh ) {
    int La=2, Lb=1, Lc=1, Ld=0;
    int i, ih, iv;
    // (DS|PS)
#pragma unroll
    for ( i=0, iv=55, ih=0; i<18; i++, iv++, ih++ ) eh[ih]+=ev[iv];
    // (FS|PS)
#pragma unroll
    for ( i=0, iv=73, ih=18; i<30; i++, iv++, ih++ ) eh[ih]+=ev[iv];
}

__device__ void gpu_twoint_core_os_dpps(
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
//    double ev[103], eh[102];
    //double eh[48];
    double eh[18], ef[30];
//    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;

    gpu_hrr_clear_dpps( eh, ef );
    for ( ijps=0; ijps<(*nijps); ijps++ ) {
        zeta  = LDG(vzeta[ijps]);
        dkab  = LDG(vdkab[ijps]);
        xiza  = LDG(vxiza[ijps]);
        zeta2 = HALF * zeta;
        for ( i=0; i<3; i++ ) {
            PC[i] = AC[i] + xiza*BA[i];
            PA[i] = xiza * BA[i];
        }
        for ( klps=0; klps<(*nklps); klps++ ) {
            eta  = LDG(veta[klps]);
            dk   = dkab * LDG(vdkcd[klps]);
            xizc = LDG(vxizc[klps]);
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
            ze2 = rz * eta * HALF;
#pragma unroll
            for ( i=0; i<3; i++ ) {
                WP[i] = rz*QP[i];
                WQ[i] = rz*QP[i] - QP[i];
            }
            T     = rho * PQ2;
            cssss = sqrho * dk;
            gpu_vrr_calc_dpps(
                    T, cssss, zeta2, eta2, ze2, rz, re, PA, WP, QC, WQ,
                    eh, ef );
//                    ev );
//            gpu_vrr_cint_dpps( ev, eh );
        }	// for (klps)
    }	// for (ijps)
    //gpu_hrr_calc_dpps( BA, DC, eh );
    //gpu_hrr_coef_dpps( eh, DINT );
    gpu_hrr_calc_dpps( BA, DC, eh, ef, DINT );
    gpu_hrr_coef_dpps( DINT );
}

#if 0
int gpu_twoint_os_dpps(
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
    double DINTEG[6*3*3*1];
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
    max_nzeri = ebuf_max_nzeri - 6*3*3*1;
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
            gpu_twoint_core_os_dpps(
                    &La, &Lb, &Lc, &Ld,
                    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
                    &psp_xiza[ijps0], BA,
                    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
                    &psp_xiza[klps0], DC,   AC,      DINTEG );
            ipat=((Lab != Lcd)||(ics==kcs && jcs>lcs) ? true : false );
            for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
                I2 = (iao*iao+iao)>>1;
                for ( j=0, jao=jao0; j<3; j++, jao++ ) {
                    if ( jao>iao ) { ix+=3*1; continue; }
                    IJ = I2 + jao;
                    coe0 = ( iao==jao ? HALF : ONE );
                    for ( k=0, kao=kao0; k<3; k++, kao++ ) {
                        K2 = (kao*kao+kao)>>1;
                        for ( l=0, lao=lao0; l<1; l++, lao++, ix++ ) {
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

int gpu_twoint_direct_os_dpps(
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
    double DINTEG[6*3*3*1];
    int mythread;
    float eps_eri = gpu_twoint_eps_eri(0);
    float eps_ps4 = gpu_twoint_eps_ps4(0);
    float eps_sch = gpu_twoint_eps_sch(0);

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
    
    max_nzeri -= 6*3*3*1;
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
            if ( val_ab*val_cd < eps_ps4 ) continue;
            kcs    = csp_ics[klcs];
            lcs    = csp_jcs[klcs];
            if ( val_ab*val_cd*gpu_twoint_dmax6(ics,jcs,kcs,lcs) < eps_sch ) continue;
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
            gpu_twoint_core_os_dpps(
                    &La, &Lb, &Lc, &Ld,
                    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
                    &psp_xiza[ijps0], BA,
                    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
                    &psp_xiza[klps0], DC,   AC,      DINTEG );
            ipat=((Lab != Lcd)||(ics==kcs && jcs>lcs) ? true : false );
            for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
                I2 = (iao*iao+iao)>>1;
                for ( j=0, jao=jao0; j<3; j++, jao++ ) {
                    if ( jao>iao ) { ix+=3*1; continue; }
                    IJ = I2 + jao;
                    coe0 = ( iao==jao ? HALF : ONE );
                    for ( k=0, kao=kao0; k<3; k++, kao++ ) {
                        K2 = (kao*kao+kao)>>1;
                        for ( l=0, lao=lao0; l<1; l++, lao++, ix++ ) {
                            if ( lao>kao ) continue;
                            if ( fabs(DINTEG[ix]) > eps_eri ) {
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
