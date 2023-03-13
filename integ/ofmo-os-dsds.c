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

static void ofmo_hrr_clear_dsds( double *eh ) {
    int i;
    // (DS|DS)
    for ( i=0; i<(0+36); i++ ) eh[i] = 0.e0;
}

static void ofmo_hrr_coef_dsds(
        const double *eh, double *DINT ) {
    int i, j, k, l, iao, jao, kao, lao, ix, iy;
    double coef_a, coef_ab, coef_abc;
    ix = 0;
    for ( i=0, iao=4; i<6; i++, iao++ ) {
        coef_a = DFACT[iao];
        coef_ab = coef_a;
        for ( k=0, kao=4; k<6; k++, kao++ ) {
            coef_abc = coef_ab * DFACT[kao];
            DINT[ix] = coef_abc * eh[ix];
            ix++;
        }
    }
}

static void ofmo_vrr_calc_dsds(
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
    // (ps|ps) m=[1,1]
    ev[35]=QC[0]*ev[ 8]+WQ[0]*ev[11]+ze2*ev[2];
    ev[36]=QC[1]*ev[ 8]+WQ[1]*ev[11];
    ev[37]=QC[2]*ev[ 8]+WQ[2]*ev[11];
    ev[38]=QC[0]*ev[ 9]+WQ[0]*ev[12];
    ev[39]=QC[1]*ev[ 9]+WQ[1]*ev[12]+ze2*ev[2];
    ev[40]=QC[2]*ev[ 9]+WQ[2]*ev[12];
    ev[41]=QC[0]*ev[10]+WQ[0]*ev[13];
    ev[42]=QC[1]*ev[10]+WQ[1]*ev[13];
    ev[43]=QC[2]*ev[10]+WQ[2]*ev[13]+ze2*ev[2];
    // (ds|ps) m=[0,1]
    ev[44]=QC[0]*ev[17]+WQ[0]*ev[23]+2.e0*ze2*ev[ 8];
    ev[45]=QC[1]*ev[17]+WQ[1]*ev[23];
    ev[46]=QC[2]*ev[17]+WQ[2]*ev[23];
    ev[47]=QC[0]*ev[18]+WQ[0]*ev[24];
    ev[48]=QC[1]*ev[18]+WQ[1]*ev[24]+2.e0*ze2*ev[ 9];
    ev[49]=QC[2]*ev[18]+WQ[2]*ev[24];
    ev[50]=QC[0]*ev[19]+WQ[0]*ev[25];
    ev[51]=QC[1]*ev[19]+WQ[1]*ev[25];
    ev[52]=QC[2]*ev[19]+WQ[2]*ev[25]+2.e0*ze2*ev[10];
    ev[53]=QC[0]*ev[20]+WQ[0]*ev[26]+     ze2*ev[ 9];
    ev[54]=QC[1]*ev[20]+WQ[1]*ev[26]+     ze2*ev[ 8];
    ev[55]=QC[2]*ev[20]+WQ[2]*ev[26];
    ev[56]=QC[0]*ev[21]+WQ[0]*ev[27]+     ze2*ev[10];
    ev[57]=QC[1]*ev[21]+WQ[1]*ev[27];
    ev[58]=QC[2]*ev[21]+WQ[2]*ev[27]+     ze2*ev[ 8];
    ev[59]=QC[0]*ev[22]+WQ[0]*ev[28];
    ev[60]=QC[1]*ev[22]+WQ[1]*ev[28]+     ze2*ev[10];
    ev[61]=QC[2]*ev[22]+WQ[2]*ev[28]+     ze2*ev[ 9];
    ev[62]=QC[0]*ev[23]+WQ[0]*ev[29]+2.e0*ze2*ev[11];
    ev[63]=QC[1]*ev[23]+WQ[1]*ev[29];
    ev[64]=QC[2]*ev[23]+WQ[2]*ev[29];
    ev[65]=QC[0]*ev[24]+WQ[0]*ev[30];
    ev[66]=QC[1]*ev[24]+WQ[1]*ev[30]+2.e0*ze2*ev[12];
    ev[67]=QC[2]*ev[24]+WQ[2]*ev[30];
    ev[68]=QC[0]*ev[25]+WQ[0]*ev[31];
    ev[69]=QC[1]*ev[25]+WQ[1]*ev[31];
    ev[70]=QC[2]*ev[25]+WQ[2]*ev[31]+2.e0*ze2*ev[13];
    ev[71]=QC[0]*ev[26]+WQ[0]*ev[32]+     ze2*ev[12];
    ev[72]=QC[1]*ev[26]+WQ[1]*ev[32]+     ze2*ev[11];
    ev[73]=QC[2]*ev[26]+WQ[2]*ev[32];
    ev[74]=QC[0]*ev[27]+WQ[0]*ev[33]+     ze2*ev[13];
    ev[75]=QC[1]*ev[27]+WQ[1]*ev[33];
    ev[76]=QC[2]*ev[27]+WQ[2]*ev[33]+     ze2*ev[11];
    ev[77]=QC[0]*ev[28]+WQ[0]*ev[34];
    ev[78]=QC[1]*ev[28]+WQ[1]*ev[34]+     ze2*ev[13];
    ev[79]=QC[2]*ev[28]+WQ[2]*ev[34]+     ze2*ev[12];
    // (ds|ds) m=[0,0]
    ev[ 80]=QC[0]*ev[44]+WQ[0]*ev[62]+eta2*(ev[17]-re*ev[23])
            +2.e0*ze2*ev[35];
    ev[ 81]=QC[1]*ev[45]+WQ[1]*ev[63]+eta2*(ev[17]-re*ev[23]);
    ev[ 82]=QC[2]*ev[46]+WQ[2]*ev[64]+eta2*(ev[17]-re*ev[23]);
    ev[ 83]=QC[0]*ev[45]+WQ[0]*ev[63]+2.e0*ze2*ev[36];
    ev[ 84]=QC[0]*ev[46]+WQ[0]*ev[64]+2.e0*ze2*ev[37];
    ev[ 85]=QC[1]*ev[46]+WQ[1]*ev[64];
    ev[ 86]=QC[0]*ev[47]+WQ[0]*ev[65]+eta2*(ev[18]-re*ev[24]);
    ev[ 87]=QC[1]*ev[48]+WQ[1]*ev[66]+eta2*(ev[18]-re*ev[24])
            +2.e0*ze2*ev[39];
    ev[ 88]=QC[2]*ev[49]+WQ[2]*ev[67]+eta2*(ev[18]-re*ev[24]);
    ev[ 89]=QC[0]*ev[48]+WQ[0]*ev[66];
    ev[ 90]=QC[0]*ev[49]+WQ[0]*ev[67];
    ev[ 91]=QC[1]*ev[49]+WQ[1]*ev[67]+2.e0*ze2*ev[40];
    ev[ 92]=QC[0]*ev[50]+WQ[0]*ev[68]+eta2*(ev[19]-re*ev[25]);
    ev[ 93]=QC[1]*ev[51]+WQ[1]*ev[69]+eta2*(ev[19]-re*ev[25]);
    ev[ 94]=QC[2]*ev[52]+WQ[2]*ev[70]+eta2*(ev[19]-re*ev[25])
            +2.e0*ze2*ev[43];
    ev[ 95]=QC[0]*ev[51]+WQ[0]*ev[69];
    ev[ 96]=QC[0]*ev[52]+WQ[0]*ev[70];
    ev[ 97]=QC[1]*ev[52]+WQ[1]*ev[70];
    ev[ 98]=QC[0]*ev[53]+WQ[0]*ev[71]+eta2*(ev[20]-re*ev[26])
            +     ze2*ev[38];
    ev[ 99]=QC[1]*ev[54]+WQ[1]*ev[72]+eta2*(ev[20]-re*ev[26])
            +     ze2*ev[36];
    ev[100]=QC[2]*ev[55]+WQ[2]*ev[73]+eta2*(ev[20]-re*ev[26]);
    ev[101]=QC[0]*ev[54]+WQ[0]*ev[72]+     ze2*ev[39];
    ev[102]=QC[0]*ev[55]+WQ[0]*ev[73]+     ze2*ev[40];
    ev[103]=QC[1]*ev[55]+WQ[1]*ev[73]+     ze2*ev[37];
    ev[104]=QC[0]*ev[56]+WQ[0]*ev[74]+eta2*(ev[21]-re*ev[27])
            +     ze2*ev[41];
    ev[105]=QC[1]*ev[57]+WQ[1]*ev[75]+eta2*(ev[21]-re*ev[27]);
    ev[106]=QC[2]*ev[58]+WQ[2]*ev[76]+eta2*(ev[21]-re*ev[27])
            +     ze2*ev[37];
    ev[107]=QC[0]*ev[57]+WQ[0]*ev[75]+     ze2*ev[42];
    ev[108]=QC[0]*ev[58]+WQ[0]*ev[76]+     ze2*ev[43];
    ev[109]=QC[1]*ev[58]+WQ[1]*ev[76];
    ev[110]=QC[0]*ev[59]+WQ[0]*ev[77]+eta2*(ev[22]-re*ev[28]);
    ev[111]=QC[1]*ev[60]+WQ[1]*ev[78]+eta2*(ev[22]-re*ev[28])
            +     ze2*ev[42];
    ev[112]=QC[2]*ev[61]+WQ[2]*ev[79]+eta2*(ev[22]-re*ev[28])
            +     ze2*ev[40];
    ev[113]=QC[0]*ev[60]+WQ[0]*ev[78];
    ev[114]=QC[0]*ev[61]+WQ[0]*ev[79];
    ev[115]=QC[1]*ev[61]+WQ[1]*ev[79]+     ze2*ev[43];
}

static void ofmo_vrr_cint_dsds( const double *ev, double *eh ) {
    int La=2, Lb=0, Lc=2, Ld=0;
    int i, ih, iv;
    // (DS|DS)
    for ( i=0, iv=80, ih=0; i<36; i++, iv++, ih++ ) eh[ih]+=ev[iv];
}

void ofmo_twoint_core_os_dsds(
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
    double ev[116];
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;

    DFACT = ofmo_getadd_dfact();
    ofmo_hrr_clear_dsds( DINT );
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
            ze2 = rz * eta2;
            for ( i=0; i<3; i++ ) {
                WP[i] = rz*QP[i];
                WQ[i] = rz*QP[i] - QP[i];
            }
            T     = rho * PQ2;
            cssss = sqrho * dk;
            ofmo_vrr_calc_dsds(
                    T, cssss, zeta2, eta2, ze2, rz, re, PA, WP, QC, WQ,
                    ev );
            ofmo_vrr_cint_dsds( ev, DINT );
        }	// for (klps)
    }	// for (ijps)
    ofmo_hrr_coef_dsds( DINT, DINT );
}

int ofmo_twoint_os_dsds(
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
    double DINTEG[6*1*6*1];
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
    max_nzeri = ebuf_max_nzeri - 6*1*6*1;
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
            ofmo_twoint_core_os_dsds(
                    &La, &Lb, &Lc, &Ld,
                    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
                    &psp_xiza[ijps0], BA,
                    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
                    &psp_xiza[klps0], DC,   AC,      DINTEG );
            ipat=((Lab != Lcd)||(ics==kcs && jcs>lcs) ? true : false );
#ifdef SORT_CSP
            int ijgekl = (ics>kcs);
            if (ics==kcs) ijgekl = (jcs>=lcs);
            if (!ijgekl) ipat = ( (ics==kcs && jcs<lcs) ? true : false);
#endif
            for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
                I2 = (iao*iao+iao)>>1;
                for ( j=0, jao=jao0; j<1; j++, jao++ ) {
                    if ( jao>iao ) { ix+=6*1; continue; }
                    IJ = I2 + jao;
                    coe0 = ( iao==jao ? HALF : ONE );
                    for ( k=0, kao=kao0; k<6; k++, kao++ ) {
                        K2 = (kao*kao+kao)>>1;
                        for ( l=0, lao=lao0; l<1; l++, lao++, ix++ ) {
                            if ( lao>kao ) continue;
                            if ( fabs(DINTEG[ix]) > eps_eri ) {
                                KL = K2 + lao;
#ifndef SORT_CSP
                                if ( IJ >= KL ) {
#else
                                if ((ijgekl&&IJ>=KL) || (!ijgekl&&KL>=IJ)) {
#endif
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

int ofmo_twoint_direct_os_dsds(
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
    double DINTEG[6*1*6*1];
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
    
    max_nzeri -= 6*1*6*1;
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
            ofmo_twoint_core_os_dsds(
                    &La, &Lb, &Lc, &Ld,
                    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
                    &psp_xiza[ijps0], BA,
                    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
                    &psp_xiza[klps0], DC,   AC,      DINTEG );
            ipat=((Lab != Lcd)||(ics==kcs && jcs>lcs) ? true : false );
#ifdef SORT_CSP
            int ijgekl = (ics>kcs);
            if (ics==kcs) ijgekl = (jcs>=lcs);
            if (!ijgekl) ipat = ( (ics==kcs && jcs<lcs) ? true : false);
#endif
            for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
                I2 = (iao*iao+iao)>>1;
                for ( j=0, jao=jao0; j<1; j++, jao++ ) {
                    if ( jao>iao ) { ix+=6*1; continue; }
                    IJ = I2 + jao;
                    coe0 = ( iao==jao ? HALF : ONE );
                    for ( k=0, kao=kao0; k<6; k++, kao++ ) {
                        K2 = (kao*kao+kao)>>1;
                        for ( l=0, lao=lao0; l<1; l++, lao++, ix++ ) {
                            if ( lao>kao ) continue;
                            if ( fabs(DINTEG[ix]) > eps_eri ) {
                                KL = K2 + lao;
#ifndef SORT_CSP
                                if ( IJ >= KL ) {
#else
                                if ((ijgekl&&IJ>=KL) || (!ijgekl&&KL>=IJ)) {
#endif
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
