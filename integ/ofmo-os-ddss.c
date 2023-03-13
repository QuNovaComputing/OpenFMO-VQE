#include <stdio.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "ofmo-twoint.h"
#include "fmt-m.h"

#ifndef false
#define false 0
#endif
#ifndef true
#define true 1
#endif

#define HALF    0.5e0
#define ONE     1.e0
#define ZERO    0.e0


#define OFMO_EBUF_FULL          1
#define OFMO_EBUF_NOFULL        0

#define EPS_PS_PAIR     1.e-32
#define EPS_CS_PAIR2    1.e-30
#define MAXNPSPAIR 100

extern double* ofmo_os_getadd_vrr( const int mythread );
extern double* ofmo_os_getadd_hrr( const int mythread );
extern double* ofmo_os_getadd_eri( const int mythread );

extern void fmt( double F[],
        const int m, const double T, const double cssss );
extern int ofmo_integ_add_fock( const int nao, const size_t nstored_eri,
        const double eri_val[], const short int eri_ind4[],
        const double D[], double G[] );

extern double* ofmo_getadd_dfact();
static double *DFACT = NULL;

static void ofmo_hrr_clear_ddss( double *eh ) {
    int i;
    // (DS|SS)
    for ( i=0; i<(0+6); i++ ) eh[i] = 0.e0;
    // (FS|SS)
    for ( i=6; i<(6+10); i++ ) eh[i] = 0.e0;
    // (GS|SS)
    for ( i=16; i<(16+15); i++ ) eh[i] = 0.e0;
}

static void ofmo_hrr_calc_ddss(
        const double BA[3], const double DC[3], double *eh ) {
    // (DP,SS)
    eh[  31] = eh[   6] - BA[0]*eh[   0];
    eh[  32] = eh[   9] - BA[1]*eh[   0];
    eh[  33] = eh[  10] - BA[2]*eh[   0];
    eh[  34] = eh[  11] - BA[0]*eh[   1];
    eh[  35] = eh[   7] - BA[1]*eh[   1];
    eh[  36] = eh[  14] - BA[2]*eh[   1];
    eh[  37] = eh[  12] - BA[0]*eh[   2];
    eh[  38] = eh[  15] - BA[1]*eh[   2];
    eh[  39] = eh[   8] - BA[2]*eh[   2];
    eh[  40] = eh[   9] - BA[0]*eh[   3];
    eh[  41] = eh[  11] - BA[1]*eh[   3];
    eh[  42] = eh[  13] - BA[2]*eh[   3];
    eh[  43] = eh[  10] - BA[0]*eh[   4];
    eh[  44] = eh[  13] - BA[1]*eh[   4];
    eh[  45] = eh[  12] - BA[2]*eh[   4];
    eh[  46] = eh[  13] - BA[0]*eh[   5];
    eh[  47] = eh[  14] - BA[1]*eh[   5];
    eh[  48] = eh[  15] - BA[2]*eh[   5];
    // (FP,SS)
    eh[  49] = eh[  16] - BA[0]*eh[   6];
    eh[  50] = eh[  19] - BA[1]*eh[   6];
    eh[  51] = eh[  20] - BA[2]*eh[   6];
    eh[  52] = eh[  24] - BA[0]*eh[   7];
    eh[  53] = eh[  17] - BA[1]*eh[   7];
    eh[  54] = eh[  28] - BA[2]*eh[   7];
    eh[  55] = eh[  25] - BA[0]*eh[   8];
    eh[  56] = eh[  29] - BA[1]*eh[   8];
    eh[  57] = eh[  18] - BA[2]*eh[   8];
    eh[  58] = eh[  19] - BA[0]*eh[   9];
    eh[  59] = eh[  21] - BA[1]*eh[   9];
    eh[  60] = eh[  23] - BA[2]*eh[   9];
    eh[  61] = eh[  20] - BA[0]*eh[  10];
    eh[  62] = eh[  23] - BA[1]*eh[  10];
    eh[  63] = eh[  22] - BA[2]*eh[  10];
    eh[  64] = eh[  21] - BA[0]*eh[  11];
    eh[  65] = eh[  24] - BA[1]*eh[  11];
    eh[  66] = eh[  26] - BA[2]*eh[  11];
    eh[  67] = eh[  22] - BA[0]*eh[  12];
    eh[  68] = eh[  27] - BA[1]*eh[  12];
    eh[  69] = eh[  25] - BA[2]*eh[  12];
    eh[  70] = eh[  23] - BA[0]*eh[  13];
    eh[  71] = eh[  26] - BA[1]*eh[  13];
    eh[  72] = eh[  27] - BA[2]*eh[  13];
    eh[  73] = eh[  26] - BA[0]*eh[  14];
    eh[  74] = eh[  28] - BA[1]*eh[  14];
    eh[  75] = eh[  30] - BA[2]*eh[  14];
    eh[  76] = eh[  27] - BA[0]*eh[  15];
    eh[  77] = eh[  30] - BA[1]*eh[  15];
    eh[  78] = eh[  29] - BA[2]*eh[  15];
    // (DD,SS)
    eh[  79] = eh[  49] - BA[0]*eh[  31];
    eh[  80] = eh[  59] - BA[1]*eh[  32];
    eh[  81] = eh[  63] - BA[2]*eh[  33];
    eh[  82] = eh[  50] - BA[0]*eh[  32];
    eh[  83] = eh[  51] - BA[0]*eh[  33];
    eh[  84] = eh[  60] - BA[1]*eh[  33];
    eh[  85] = eh[  64] - BA[0]*eh[  34];
    eh[  86] = eh[  53] - BA[1]*eh[  35];
    eh[  87] = eh[  75] - BA[2]*eh[  36];
    eh[  88] = eh[  65] - BA[0]*eh[  35];
    eh[  89] = eh[  66] - BA[0]*eh[  36];
    eh[  90] = eh[  54] - BA[1]*eh[  36];
    eh[  91] = eh[  67] - BA[0]*eh[  37];
    eh[  92] = eh[  77] - BA[1]*eh[  38];
    eh[  93] = eh[  57] - BA[2]*eh[  39];
    eh[  94] = eh[  68] - BA[0]*eh[  38];
    eh[  95] = eh[  69] - BA[0]*eh[  39];
    eh[  96] = eh[  78] - BA[1]*eh[  39];
    eh[  97] = eh[  58] - BA[0]*eh[  40];
    eh[  98] = eh[  65] - BA[1]*eh[  41];
    eh[  99] = eh[  72] - BA[2]*eh[  42];
    eh[ 100] = eh[  59] - BA[0]*eh[  41];
    eh[ 101] = eh[  60] - BA[0]*eh[  42];
    eh[ 102] = eh[  66] - BA[1]*eh[  42];
    eh[ 103] = eh[  61] - BA[0]*eh[  43];
    eh[ 104] = eh[  71] - BA[1]*eh[  44];
    eh[ 105] = eh[  69] - BA[2]*eh[  45];
    eh[ 106] = eh[  62] - BA[0]*eh[  44];
    eh[ 107] = eh[  63] - BA[0]*eh[  45];
    eh[ 108] = eh[  72] - BA[1]*eh[  45];
    eh[ 109] = eh[  70] - BA[0]*eh[  46];
    eh[ 110] = eh[  74] - BA[1]*eh[  47];
    eh[ 111] = eh[  78] - BA[2]*eh[  48];
    eh[ 112] = eh[  71] - BA[0]*eh[  47];
    eh[ 113] = eh[  72] - BA[0]*eh[  48];
    eh[ 114] = eh[  75] - BA[1]*eh[  48];
    // HRR for (XX|XX)-type integral (center CD)
}

static void ofmo_hrr_coef_ddss(
        const double *eh, double *DINT ) {
    int i, j, k, l, iao, jao, kao, lao, ix, iy;
    double coef_a, coef_ab, coef_abc;
    ix = 79;
    iy = 0;

    for ( i=0, iao=4; i<6; i++, iao++ ) {
        coef_a = DFACT[iao];
        for ( j=0, jao=4; j<6; j++, jao++ ) {
            coef_ab = coef_a * DFACT[jao];
            coef_abc = coef_ab;
            DINT[iy] = coef_abc * eh[ix];
            iy++;
            ix++;
        }
    }
}

static void ofmo_vrr_calc_ddss(
        const double T, const double cssss,
        const double zeta2, const double eta2, const double ze2,
        const double rz, const double re,
        const double PA[3], const double WP[3],
        const double QC[3], const double WQ[3],
        double *ev ) {
    // (ss|ss) m=0,4
    //fmt( &ev[0], 4, T, cssss );
    OFMO_FMT( &ev[0], 4, T, cssss );
    // (ps|ss) m=0,3
    ev[ 5]=PA[0]*ev[0]+WP[0]*ev[1];
    ev[ 6]=PA[1]*ev[0]+WP[1]*ev[1];
    ev[ 7]=PA[2]*ev[0]+WP[2]*ev[1];
    ev[ 8]=PA[0]*ev[1]+WP[0]*ev[2];
    ev[ 9]=PA[1]*ev[1]+WP[1]*ev[2];
    ev[10]=PA[2]*ev[1]+WP[2]*ev[2];
    ev[11]=PA[0]*ev[2]+WP[0]*ev[3];
    ev[12]=PA[1]*ev[2]+WP[1]*ev[3];
    ev[13]=PA[2]*ev[2]+WP[2]*ev[3];
    ev[14]=PA[0]*ev[3]+WP[0]*ev[4];
    ev[15]=PA[1]*ev[3]+WP[1]*ev[4];
    ev[16]=PA[2]*ev[3]+WP[2]*ev[4];
    // (ds|ss) m=0,2
    ev[17]=PA[0]*ev[ 5]+WP[0]*ev[ 8]+zeta2*(ev[0]-rz*ev[1]);
    ev[18]=PA[1]*ev[ 6]+WP[1]*ev[ 9]+zeta2*(ev[0]-rz*ev[1]);
    ev[19]=PA[2]*ev[ 7]+WP[2]*ev[10]+zeta2*(ev[0]-rz*ev[1]);
    ev[20]=PA[0]*ev[ 6]+WP[0]*ev[ 9];
    ev[21]=PA[0]*ev[ 7]+WP[0]*ev[10];
    ev[22]=PA[1]*ev[ 7]+WP[1]*ev[10];
    ev[23]=PA[0]*ev[ 8]+WP[0]*ev[11]+zeta2*(ev[1]-rz*ev[2]);
    ev[24]=PA[1]*ev[ 9]+WP[1]*ev[12]+zeta2*(ev[1]-rz*ev[2]);
    ev[25]=PA[2]*ev[10]+WP[2]*ev[13]+zeta2*(ev[1]-rz*ev[2]);
    ev[26]=PA[0]*ev[ 9]+WP[0]*ev[12];
    ev[27]=PA[0]*ev[10]+WP[0]*ev[13];
    ev[28]=PA[1]*ev[10]+WP[1]*ev[13];
    ev[29]=PA[0]*ev[11]+WP[0]*ev[14]+zeta2*(ev[2]-rz*ev[3]);
    ev[30]=PA[1]*ev[12]+WP[1]*ev[15]+zeta2*(ev[2]-rz*ev[3]);
    ev[31]=PA[2]*ev[13]+WP[2]*ev[16]+zeta2*(ev[2]-rz*ev[3]);
    ev[32]=PA[0]*ev[12]+WP[0]*ev[15];
    ev[33]=PA[0]*ev[13]+WP[0]*ev[16];
    ev[34]=PA[1]*ev[13]+WP[1]*ev[16];
    // (fs|ss) m=0,1
    ev[35]=PA[0]*ev[17]+WP[0]*ev[23]+2.e0*zeta2*(ev[ 5]-rz*ev[ 8]);
    ev[36]=PA[1]*ev[18]+WP[1]*ev[24]+2.e0*zeta2*(ev[ 6]-rz*ev[ 9]);
    ev[37]=PA[2]*ev[19]+WP[2]*ev[25]+2.e0*zeta2*(ev[ 7]-rz*ev[10]);
    ev[38]=PA[1]*ev[17]+WP[1]*ev[23];
    ev[39]=PA[2]*ev[17]+WP[2]*ev[23];
    ev[40]=PA[0]*ev[18]+WP[0]*ev[24];
    ev[41]=PA[0]*ev[19]+WP[0]*ev[25];
    ev[42]=PA[0]*ev[22]+WP[0]*ev[28];
    ev[43]=PA[2]*ev[18]+WP[2]*ev[24];
    ev[44]=PA[1]*ev[19]+WP[1]*ev[25];
    ev[45]=PA[0]*ev[23]+WP[0]*ev[29]+2.e0*zeta2*(ev[ 8]-rz*ev[11]);
    ev[46]=PA[1]*ev[24]+WP[1]*ev[30]+2.e0*zeta2*(ev[ 9]-rz*ev[12]);
    ev[47]=PA[2]*ev[25]+WP[2]*ev[31]+2.e0*zeta2*(ev[10]-rz*ev[13]);
    ev[48]=PA[1]*ev[23]+WP[1]*ev[29];
    ev[49]=PA[2]*ev[23]+WP[2]*ev[29];
    ev[50]=PA[0]*ev[24]+WP[0]*ev[30];
    ev[51]=PA[0]*ev[25]+WP[0]*ev[31];
    ev[52]=PA[0]*ev[28]+WP[0]*ev[34];
    ev[53]=PA[2]*ev[24]+WP[2]*ev[30];
    ev[54]=PA[1]*ev[25]+WP[1]*ev[31];
    // (gs|ss) m=0,0
    ev[55]=PA[0]*ev[35]+WP[0]*ev[45]+3.e0*zeta2*(ev[17]-rz*ev[23]);
    ev[56]=PA[1]*ev[36]+WP[1]*ev[46]+3.e0*zeta2*(ev[18]-rz*ev[24]);
    ev[57]=PA[2]*ev[37]+WP[2]*ev[47]+3.e0*zeta2*(ev[19]-rz*ev[25]);
    ev[58]=PA[1]*ev[35]+WP[1]*ev[45];
    ev[59]=PA[2]*ev[35]+WP[2]*ev[45];
    ev[60]=PA[0]*ev[40]+WP[0]*ev[50]+     zeta2*(ev[18]-rz*ev[24]);
    ev[61]=PA[0]*ev[41]+WP[0]*ev[51]+     zeta2*(ev[19]-rz*ev[25]);
    ev[62]=PA[1]*ev[39]+WP[1]*ev[49];
    ev[63]=PA[0]*ev[36]+WP[0]*ev[46];
    ev[64]=PA[0]*ev[37]+WP[0]*ev[47];
    ev[65]=PA[0]*ev[43]+WP[0]*ev[53];
    ev[66]=PA[0]*ev[44]+WP[0]*ev[54];
    ev[67]=PA[2]*ev[36]+WP[2]*ev[46];
    ev[68]=PA[1]*ev[37]+WP[1]*ev[47];
    ev[69]=PA[1]*ev[44]+WP[1]*ev[54]+     zeta2*(ev[19]-rz*ev[25]);
}

static void ofmo_vrr_cint_ddss( const double *ev, double *eh ) {
    int La=2, Lb=2, Lc=0, Ld=0;
    int i, ih, iv;
    // (DS|SS)
    for ( i=0, iv=17, ih=0; i<6; i++, iv++, ih++ ) eh[ih]+=ev[iv];
    // (FS|SS)
    for ( i=0, iv=35, ih=6; i<10; i++, iv++, ih++ ) eh[ih]+=ev[iv];
    // (GS|SS)
    for ( i=0, iv=55, ih=16; i<15; i++, iv++, ih++ ) eh[ih]+=ev[iv];
}

void ofmo_twoint_core_os_ddss(
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *nijps, const double vzeta[], const double vdkab[],
        const double vxiza[], const double BA[3],
        const int *nklps, const double veta[], const double vdkcd[],
        const double vxizc[], const double DC[3], const double AC[3],
        double *DINT ) {
    int ijps, klps, i;
    double cssss, zeta, dkab, xiza, eta, xizc, dk, T;
    double zeta2, eta2, ze2, rz, re, PA[3], WP[3], QC[3], WQ[3];
    double PQ2, sqrho, rho, PC[3], QP[3];
    double ev[70], eh[115];
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;

    DFACT = ofmo_getadd_dfact();
    ofmo_hrr_clear_ddss( eh );
    for ( ijps=0; ijps<(*nijps); ijps++ ) {
        zeta  = vzeta[ijps];
        dkab  = vdkab[ijps];
        xiza  = vxiza[ijps];
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
            ze2 = rz * eta * HALF;
            for ( i=0; i<3; i++ ) {
                WP[i] = rz*QP[i];
            }
            T     = rho * PQ2;
            cssss = sqrho * dk;
            ofmo_vrr_calc_ddss(
                    T, cssss, zeta2, eta2, ze2, rz, re, PA, WP, QC, WQ,
                    ev );
            ofmo_vrr_cint_ddss( ev, eh );
        }	// for (klps)
    }	// for (ijps)
    ofmo_hrr_calc_ddss( BA, DC, eh );
    ofmo_hrr_coef_ddss( eh, DINT );
}

int ofmo_twoint_os_ddss(
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
    double DINTEG[6*6*1*1];
    long nzeri, max_nzeri, nzeri4;
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    long ebuf_max_nzeri = *pebuf_max_nzeri;
    int mythread;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);

    
    mythread = omp_get_thread_num();
    if ( DFACT == NULL ) DFACT = ofmo_getadd_dfact();
    
    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs0 = leading_cs_pair[Lab];
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];
    nzeri     = *ebuf_non_zero_eri;
    max_nzeri = ebuf_max_nzeri - 6*6*1*1;
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
            if ( val_ab*val_cd < eps_ps4 ) continue;
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
            ofmo_twoint_core_os_ddss(
                    &La, &Lb, &Lc, &Ld,
                    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
                    &psp_xiza[ijps0], BA,
                    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
                    &psp_xiza[klps0], DC,   AC,      DINTEG );
            ipat=((Lab != Lcd)||(ics==kcs && jcs>lcs) ? true : false );
            for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
                I2 = (iao*iao+iao)>>1;
                for ( j=0, jao=jao0; j<6; j++, jao++ ) {
                    if ( jao>iao ) { ix+=1*1; continue; }
                    IJ = I2 + jao;
                    coe0 = ( iao==jao ? HALF : ONE );
                    for ( k=0, kao=kao0; k<1; k++, kao++ ) {
                        K2 = (kao*kao+kao)>>1;
                        for ( l=0, lao=lao0; l<1; l++, lao++, ix++ ) {
                            if ( lao>kao ) continue;
                            if ( fabs(DINTEG[ix]) > eps_eri ) {
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

int ofmo_twoint_direct_os_ddss(
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
    double DINTEG[6*6*1*1];
    int mythread;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);
    float eps_sch = ofmo_twoint_eps_sch(0);

    mythread = omp_get_thread_num();
    if ( DFACT == NULL ) DFACT = ofmo_getadd_dfact();
    
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
    
    max_nzeri -= 6*6*1*1;
    nzeri4    = nzeri*4;
    if ( nzeri >= max_nzeri ) {
        ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4, Ds, G );
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
            if ( val_ab*val_cd*ofmo_twoint_dmax6(ics,jcs,kcs,lcs) < eps_sch ) continue;
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
            ofmo_twoint_core_os_ddss(
                    &La, &Lb, &Lc, &Ld,
                    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
                    &psp_xiza[ijps0], BA,
                    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
                    &psp_xiza[klps0], DC,   AC,      DINTEG );
            ipat=((Lab != Lcd)||(ics==kcs && jcs>lcs) ? true : false );
            for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
                I2 = (iao*iao+iao)>>1;
                for ( j=0, jao=jao0; j<6; j++, jao++ ) {
                    if ( jao>iao ) { ix+=1*1; continue; }
                    IJ = I2 + jao;
                    coe0 = ( iao==jao ? HALF : ONE );
                    for ( k=0, kao=kao0; k<1; k++, kao++ ) {
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
                ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4,
                        Ds, G );
                nzeri = nzeri4= 0;
            }
        }	// for (klcs)
        klcs = klcs0;
    }	// for (ijcs)
    *petmp_non_zero_eri = nzeri;
    return 0;
}
