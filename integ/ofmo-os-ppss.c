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

static void ofmo_hrr_clear_ppss( double *eh ) {
    int i;
    // (PS|SS)
    for ( i=0; i<(0+3); i++ ) eh[i] = 0.e0;
    // (DS|SS)
    for ( i=3; i<(3+6); i++ ) eh[i] = 0.e0;
}

static void ofmo_hrr_calc_ppss(
        const double BA[3], const double DC[3], double *eh ) {
    // (PP,SS)
    eh[   9] = eh[   3] - BA[0]*eh[   0];
    eh[  10] = eh[   6] - BA[1]*eh[   0];
    eh[  11] = eh[   7] - BA[2]*eh[   0];
    eh[  12] = eh[   6] - BA[0]*eh[   1];
    eh[  13] = eh[   4] - BA[1]*eh[   1];
    eh[  14] = eh[   8] - BA[2]*eh[   1];
    eh[  15] = eh[   7] - BA[0]*eh[   2];
    eh[  16] = eh[   8] - BA[1]*eh[   2];
    eh[  17] = eh[   5] - BA[2]*eh[   2];
    // HRR for (XX|XX)-type integral (center CD)
}

static void ofmo_hrr_coef_ppss(
        const double *eh, double *DINT ) {
    memcpy( DINT, &eh[9], sizeof(double)*3*3*1*1 );
}

static void ofmo_vrr_calc_ppss(
        const double T, const double cssss,
        const double zeta2, const double eta2, const double ze2,
        const double rz, const double re,
        const double PA[3], const double WP[3],
        const double QC[3], const double WQ[3],
        double *ev ) {
    // (ss|ss) m=0,2
    //fmt( &ev[0], 2, T, cssss );
    OFMO_FMT( &ev[0], 2, T, cssss );
    // (ps|ss) m=0,1
    ev[3]=PA[0]*ev[0]+WP[0]*ev[1];
    ev[4]=PA[1]*ev[0]+WP[1]*ev[1];
    ev[5]=PA[2]*ev[0]+WP[2]*ev[1];
    ev[6]=PA[0]*ev[1]+WP[0]*ev[2];
    ev[7]=PA[1]*ev[1]+WP[1]*ev[2];
    ev[8]=PA[2]*ev[1]+WP[2]*ev[2];
    // (ds|ss) m=0,0
    ev[ 9]=PA[0]*ev[3]+WP[0]*ev[6]+zeta2*(ev[0]-rz*ev[1]);
    ev[10]=PA[1]*ev[4]+WP[1]*ev[7]+zeta2*(ev[0]-rz*ev[1]);
    ev[11]=PA[2]*ev[5]+WP[2]*ev[8]+zeta2*(ev[0]-rz*ev[1]);
    ev[12]=PA[0]*ev[4]+WP[0]*ev[7];
    ev[13]=PA[0]*ev[5]+WP[0]*ev[8];
    ev[14]=PA[1]*ev[5]+WP[1]*ev[8];
}

static void ofmo_vrr_cint_ppss( const double *ev, double *eh ) {
    int La=1, Lb=1, Lc=0, Ld=0;
    int i, ih, iv;
    // (PS|SS)
    for ( i=0, iv=3, ih=0; i<3; i++, iv++, ih++ ) eh[ih]+=ev[iv];
    // (DS|SS)
    for ( i=0, iv=9, ih=3; i<6; i++, iv++, ih++ ) eh[ih]+=ev[iv];
}

void ofmo_twoint_core_os_ppss(
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
    double ev[15], eh[18];
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;

    ofmo_hrr_clear_ppss( eh );
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
            ofmo_vrr_calc_ppss(
                    T, cssss, zeta2, eta2, ze2, rz, re, PA, WP, QC, WQ,
                    ev );
            ofmo_vrr_cint_ppss( ev, eh );
        }	// for (klps)
    }	// for (ijps)
    ofmo_hrr_calc_ppss( BA, DC, eh );
    ofmo_hrr_coef_ppss( eh, DINT );
}

int ofmo_twoint_os_ppss(
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
    double DINTEG[3*3*1*1];
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
    max_nzeri = ebuf_max_nzeri - 3*3*1*1;
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
            ofmo_twoint_core_os_ppss(
                    &La, &Lb, &Lc, &Ld,
                    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
                    &psp_xiza[ijps0], BA,
                    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
                    &psp_xiza[klps0], DC,   AC,      DINTEG );
            ipat=((Lab != Lcd)||(ics==kcs && jcs>lcs) ? true : false );
            for ( i=0, iao=iao0, ix=0; i<3; i++, iao++ ) {
                I2 = (iao*iao+iao)>>1;
                for ( j=0, jao=jao0; j<3; j++, jao++ ) {
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

int ofmo_twoint_direct_os_ppss(
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
    double DINTEG[3*3*1*1];
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
    
    max_nzeri -= 3*3*1*1;
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
            ofmo_twoint_core_os_ppss(
                    &La, &Lb, &Lc, &Ld,
                    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
                    &psp_xiza[ijps0], BA,
                    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
                    &psp_xiza[klps0], DC,   AC,      DINTEG );
            ipat=((Lab != Lcd)||(ics==kcs && jcs>lcs) ? true : false );
            for ( i=0, iao=iao0, ix=0; i<3; i++, iao++ ) {
                I2 = (iao*iao+iao)>>1;
                for ( j=0, jao=jao0; j<3; j++, jao++ ) {
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
