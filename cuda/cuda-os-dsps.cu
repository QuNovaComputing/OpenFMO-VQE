// dsps
#include "cuda-twoint-core-os.h"

__device__ void gpu_hrr_clear_os_dsps( double *eh ) {
    int i;
    // (DS|PS)
#pragma unroll
    for ( i=0; i<(0+18); i++ ) eh[i] = 0.e0;
}

__device__ void gpu_hrr_coef_os_dsps(
        const double *eh, double *DINT ) {
    int i, j, k, l, iao, jao, kao, lao, ix, iy;
    double coef_a, coef_ab, coef_abc;
    ix = 0;
#pragma unroll
    for ( i=0, iao=4; i<6; i++, iao++ ) {
        coef_a = LDG(DFACT[iao]);
        coef_ab = coef_a;
#pragma unroll
        for ( k=0, kao=1; k<3; k++, kao++ ) {
            coef_abc = coef_ab * LDG(DFACT[kao]);
            //DINT[ix] = coef_abc * eh[ix];
            //ix++;
            DINT[ix++] *= coef_abc;
        }
    }
}

__device__ void gpu_vrr_calc_os_dsps(
        const double T, const double cssss,
        const double zeta2, const double eta2, const double ze2,
        const double rz, const double re,
        const double PA[3], const double WP[3],
        const double QC[3], const double WQ[3],
        double *eh ) {
    double tmp0,tmp1;
    double ev[12];
    double ep[9];
    // (ss|ss) m=0,3
    //gpu_fmt( &ev[0], 3, T, cssss );
#if   CUDA_FMT_M == 3
    gpu_fmt3_method3( T, cssss, ev );
#elif CUDA_FMT_M == 2
    gpu_fmt3_method2( T, cssss, ev );
#elif CUDA_FMT_M == 1
    gpu_fmt3_method1( T, cssss, ev );
#else
    gpu_fmt3( ev, T, cssss );
#endif
    tmp0 = zeta2*(ev[0]-rz*ev[1]);
    tmp1 = zeta2*(ev[1]-rz*ev[2]);
    // (ps|ss) m=0,2
    ep[ 0]=PA[0]*ev[0]+WP[0]*ev[1];
    ep[ 1]=PA[1]*ev[0]+WP[1]*ev[1];
    ep[ 2]=PA[2]*ev[0]+WP[2]*ev[1];
    ep[ 3]=PA[0]*ev[1]+WP[0]*ev[2];
    ep[ 4]=PA[1]*ev[1]+WP[1]*ev[2];
    ep[ 5]=PA[2]*ev[1]+WP[2]*ev[2];
    ep[ 6]=PA[0]*ev[2]+WP[0]*ev[3];
    ep[ 7]=PA[1]*ev[2]+WP[1]*ev[3];
    ep[ 8]=PA[2]*ev[2]+WP[2]*ev[3];
    // (ds|ss) m=0,1
    ev[ 0]=PA[0]*ep[0]+WP[0]*ep[ 3]+tmp0;
    ev[ 1]=PA[0]*ep[3]+WP[0]*ep[ 6]+tmp1;
    ev[ 2]=PA[1]*ep[1]+WP[1]*ep[ 4]+tmp0;
    ev[ 3]=PA[1]*ep[4]+WP[1]*ep[ 7]+tmp1;
    ev[ 4]=PA[2]*ep[2]+WP[2]*ep[ 5]+tmp0;
    ev[ 5]=PA[2]*ep[5]+WP[2]*ep[ 8]+tmp1;
    ev[ 6]=PA[0]*ep[1]+WP[0]*ep[ 4];
    ev[ 7]=PA[0]*ep[4]+WP[0]*ep[ 7];
    ev[ 8]=PA[0]*ep[2]+WP[0]*ep[ 5];
    ev[ 9]=PA[0]*ep[5]+WP[0]*ep[ 8];
    ev[10]=PA[1]*ep[2]+WP[1]*ep[ 5];
    ev[11]=PA[1]*ep[5]+WP[1]*ep[ 8];
    // (ds|ps) m=[0,0]
    eh[ 0]+=QC[0]*ev[ 0]+WQ[0]*ev[ 1]+2.e0*ze2*ep[3];
    eh[ 1]+=QC[1]*ev[ 0]+WQ[1]*ev[ 1];
    eh[ 2]+=QC[2]*ev[ 0]+WQ[2]*ev[ 1];
    eh[ 3]+=QC[0]*ev[ 2]+WQ[0]*ev[ 3];
    eh[ 4]+=QC[1]*ev[ 2]+WQ[1]*ev[ 3]+2.e0*ze2*ep[4];
    eh[ 5]+=QC[2]*ev[ 2]+WQ[2]*ev[ 3];
    eh[ 6]+=QC[0]*ev[ 4]+WQ[0]*ev[ 5];
    eh[ 7]+=QC[1]*ev[ 4]+WQ[1]*ev[ 5];
    eh[ 8]+=QC[2]*ev[ 4]+WQ[2]*ev[ 5]+2.e0*ze2*ep[5];
    eh[ 9]+=QC[0]*ev[ 6]+WQ[0]*ev[ 7]+     ze2*ep[4];
    eh[10]+=QC[1]*ev[ 6]+WQ[1]*ev[ 7]+     ze2*ep[3];
    eh[11]+=QC[2]*ev[ 6]+WQ[2]*ev[ 7];
    eh[12]+=QC[0]*ev[ 8]+WQ[0]*ev[ 9]+     ze2*ep[5];
    eh[13]+=QC[1]*ev[ 8]+WQ[1]*ev[ 9];
    eh[14]+=QC[2]*ev[ 8]+WQ[2]*ev[ 9]+     ze2*ep[3];
    eh[15]+=QC[0]*ev[10]+WQ[0]*ev[11];
    eh[16]+=QC[1]*ev[10]+WQ[1]*ev[11]+     ze2*ep[5];
    eh[17]+=QC[2]*ev[10]+WQ[2]*ev[11]+     ze2*ep[4];
}

__device__ void gpu_vrr_cint_os_dsps( const double *ev, double *eh ) {
    int La=2, Lb=0, Lc=1, Ld=0;
    int i, ih, iv;
    // (DS|PS)
#pragma unroll
    //for ( i=0, iv=25, ih=0; i<18; i++, iv++, ih++ ) eh[ih]+=ev[iv];
    for ( i=0; i<18; i++) eh[i]+=ev[i];
}

__device__ void gpu_twoint_core_os_dsps(
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
//    double ev[18];
//    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;

    gpu_hrr_clear_os_dsps( DINT );
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
            gpu_vrr_calc_os_dsps(
                    T, cssss, zeta2, eta2, ze2, rz, re, PA, WP, QC, WQ,
                    //ev );
                    DINT );
            //gpu_vrr_cint_os_dsps( ev, DINT );
        }	// for (klps)
    }	// for (ijps)
    gpu_hrr_coef_os_dsps( DINT, DINT );
}

#if 0
int gpu_twoint_os_dsps(
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
    double DINTEG[6*1*3*1];
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
    max_nzeri = ebuf_max_nzeri - 6*1*3*1;
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
            gpu_twoint_core_os_dsps(
                    &La, &Lb, &Lc, &Ld,
                    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
                    &psp_xiza[ijps0], BA,
                    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
                    &psp_xiza[klps0], DC,   AC,      DINTEG );
            ipat=((Lab != Lcd)||(ics==kcs && jcs>lcs) ? true : false );
            for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
                I2 = (iao*iao+iao)>>1;
                for ( j=0, jao=jao0; j<1; j++, jao++ ) {
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

int gpu_twoint_direct_os_dsps(
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
    double DINTEG[6*1*3*1];
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
    
    max_nzeri -= 6*1*3*1;
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
            gpu_twoint_core_os_dsps(
                    &La, &Lb, &Lc, &Ld,
                    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
                    &psp_xiza[ijps0], BA,
                    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
                    &psp_xiza[klps0], DC,   AC,      DINTEG );
            ipat=((Lab != Lcd)||(ics==kcs && jcs>lcs) ? true : false );
            for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
                I2 = (iao*iao+iao)>>1;
                for ( j=0, jao=jao0; j<1; j++, jao++ ) {
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
