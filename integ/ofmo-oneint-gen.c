#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ofmo-index.h"

#ifdef _OPENMP
#include <omp.h>
#else
#include "omp-dummy.h"
#endif

extern void fmt( double F[],
	const int m, const double T, const double cssss );

#ifndef false
#define false 0
#endif

#ifndef true
#define true 1
#endif

#define HALF	0.5e0
#define ZERO	0.e0

#define EPS_PS_PAIR 1.e-32

// 2次元整数配列の確保
static int** ofmo_alloc_imatrix( const int na, const int nb ) {
    int **ip, i;
    ip = (int**)malloc( sizeof(int*) * na );
    ip[0] = (int*)malloc( sizeof(int) * na * nb );
    for (i=1; i<na; i++ ) ip[i] = ip[i-1] + nb;
    return ip;
}

// 二次元整数配列の解放
static void ofmo_free_imatrix( int** ip ) {
    if ( ip ) {
	if ( ip[0] ) free( ip[0] );
	free( ip );
    }
}

/* ofmo-index.cで生成される変数 */
static int *NNAO;
static int *LAOT;
static int *INDX;
static int **ANGM;
static int **NAM;
static double *DFACT;

/* 1電子積分のVRRで必要な変数 */
static int ***V_OVIADD = NULL;
static int ***V_NAIADD = NULL;
static int ***V_MINDEX = NULL;
static double **V_OVI  = NULL;
static double **V_KEI  = NULL;
static double **V_NAI  = NULL;

/* 縮約積分の保存領域 */
static double **OVI_MASTER = NULL;
static double **HCORE_MASTER = NULL;

/* 定数など */
static double _PI_32_;
static double _2PI_;
static double EPS_PS_PAIR_NAI;

/* 各種アドレス計算 */
static int ofmo_oneint_gen_make_add(
	const int mythread, const int La, const int Lb,
	int *OVI_MEM, int *NAI_MEM ) {
    int ovi_mem, nai_mem;
    int Lab, ma, mb, na, nb, nab, lmin, mmax;
    int **OVIADD, **NAIADD, **MINDEX;
    Lab = La + Lb;
    OVIADD = V_OVIADD[mythread];
    NAIADD = V_NAIADD[mythread];
    MINDEX = V_MINDEX[mythread];
    ovi_mem = nai_mem = 0;
    // (a|s)タイプ
    for ( ma=0; ma<=La; ma++ ) {
	na            = NNAO[ma];
	OVIADD[ma][0] = ovi_mem;
	NAIADD[ma][0] = nai_mem;
	MINDEX[ma][0] = na;
	ovi_mem += na;
	nai_mem += (Lab-ma+1) * na;
    }
    // (a|b)タイプ
    for ( mb=1; mb<=Lb; mb++ ) {
	lmin = La - Lb + mb;
	if ( lmin < 0 ) lmin = 0;
	mmax = Lb - mb;
	nb   = NNAO[mb];
	for ( ma=lmin; ma<=La; ma++ ) {
	    nab = nb * NNAO[ma];
	    OVIADD[ma][mb] = ovi_mem;
	    NAIADD[ma][mb] = nai_mem;
	    MINDEX[ma][mb] = nab;
	    ovi_mem += nab;
	    nai_mem += (mmax+1)*nab;
	}
    }
    *OVI_MEM = ovi_mem;
    *NAI_MEM = nai_mem;
    return 0;
}

static void ofmo_oneint_gen_finalize() {
    int i, nthreads;
    nthreads = omp_get_max_threads();
    for ( i=0; i<nthreads; i++ ) {
	if ( V_OVIADD[i] ) ofmo_free_imatrix( V_OVIADD[i] );
	if ( V_NAIADD[i] ) ofmo_free_imatrix( V_NAIADD[i] );
	if ( V_MINDEX[i] ) ofmo_free_imatrix( V_MINDEX[i] );
	if ( V_OVI[i] ) free( V_OVI[i] );
	if ( V_KEI[i] ) free( V_KEI[i] );
	if ( V_NAI[i] ) free( V_NAI[i] );
	if ( OVI_MASTER[i] )   free( OVI_MASTER[i] );
	if ( HCORE_MASTER[i] ) free( HCORE_MASTER[i] );
    }
    free( V_OVIADD );
    free( V_NAIADD );
    free( V_MINDEX );
    free( V_OVI );
    free( V_KEI );
    free( V_NAI );
    free( OVI_MASTER );
    free( HCORE_MASTER );
    V_OVIADD = NULL;
    V_NAIADD = NULL;
    V_MINDEX = NULL;
    V_OVI  = NULL;
    V_KEI  = NULL;
    V_NAI  = NULL;
    OVI_MASTER = NULL;
    HCORE_MASTER = NULL;
}

int ofmo_oneint_gen_init( const int maxlqn ) {
    int nthreads;
    double pi;
    static int called = false;
    if ( called ) return 0;
    pi = 4.e0 * atan( 1.e0 );
    _PI_32_ = pi * sqrt(pi);
    _2PI_   = 2.e0 * pi;
    EPS_PS_PAIR_NAI = EPS_PS_PAIR;
    ofmo_index_init( 2*maxlqn );
    NNAO = ofmo_getadd_nnao();
    LAOT = ofmo_getadd_laot();
    ANGM = ofmo_getadd_angm();
    INDX = ofmo_getadd_indx();
    NAM  = ofmo_getadd_nam();
    DFACT = ofmo_getadd_dfact();

    nthreads = omp_get_max_threads();
    V_OVIADD = (int***)malloc( sizeof(int**) * nthreads );
    V_NAIADD = (int***)malloc( sizeof(int**) * nthreads );
    V_MINDEX = (int***)malloc( sizeof(int**) * nthreads );
    V_OVI    = (double**)malloc( sizeof(double*) * nthreads );
    V_KEI    = (double**)malloc( sizeof(double*) * nthreads );
    V_NAI    = (double**)malloc( sizeof(double*) * nthreads );
    OVI_MASTER   = (double**)malloc( sizeof(double*) * nthreads );
    HCORE_MASTER = (double**)malloc( sizeof(double*) * nthreads );
#pragma omp parallel
    {
	int mythread;
	int ovi_mem, nai_mem;
	int n, n2;
	n = NNAO[maxlqn];
	n2 = n*n;
	mythread = omp_get_thread_num();
	V_OVIADD[mythread] = ofmo_alloc_imatrix( maxlqn+1, maxlqn+1 );
	V_NAIADD[mythread] = ofmo_alloc_imatrix( maxlqn+1, maxlqn+1 );
	V_MINDEX[mythread] = ofmo_alloc_imatrix( maxlqn+1, maxlqn+1 );
	ofmo_oneint_gen_make_add( mythread, maxlqn, maxlqn,
		&ovi_mem, &nai_mem );
	V_OVI[mythread] = (double*)malloc( sizeof(double) * ovi_mem );
	V_KEI[mythread] = (double*)malloc( sizeof(double) * ovi_mem );
	V_NAI[mythread] = (double*)malloc( sizeof(double) * nai_mem );
	OVI_MASTER[mythread]   = (double*)malloc( sizeof(double) * n2 );
	HCORE_MASTER[mythread] = (double*)malloc( sizeof(double) * n2 );
    }
    atexit( ofmo_oneint_gen_finalize );
    called = true;
    return 0;
}

static int oneint_core_xx(
	const int mythread, const int *pLa, const int *pLb,
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const int atomic_number[],
	double OVI[], double HCORE[] ) {
    int i, ij, m, ips, jps, kat, ips1, jps1;
    double zeta_a, zeta_b, zeta, coef_a, coef_b, coef;
    double sqrzi, zi, xiza, xizb, xi, xiab2, xi2, exp_ab, css, zeta2;
    double oviss, keiss;
    double BA[3], P[3], PA[3], PC[3], PB[3], AB2;
    double PC2, dq, U;
    double _pi32_ = _PI_32_, _2pi_ = _2PI_;
    double coe_a, coe;
    int La=*pLa, Lb=*pLb, Lab, nab, na, nb, ma, mb;
    int **OVIADD, **NAIADD, **MINDEX;
    double *ovi, *kei, *nai, *p, *ov, *ke;
    int ip0, i0p, i00, i10, i01, ix, nia, nib, lmin, mmax;
    int ip00, i0p0, i000, i001, i100, i101, i010, i011;
    int iaop, iaop0, iaop1, iao, iao0, iao1, iaom, iaom0;
    int jaop, jaop0, jaop1, jao, jao0, jao1, jaom, jaom0;

    ovi = V_OVI[mythread];
    kei = V_KEI[mythread];
    nai = V_NAI[mythread];
    OVIADD = V_OVIADD[mythread];
    NAIADD = V_NAIADD[mythread];
    MINDEX = V_MINDEX[mythread];
    na = NNAO[La];
    nb = NNAO[Lb];
    nab = na * nb;
    Lab = La + Lb;
    for ( i=0; i<nab; i++ ) OVI[i] = HCORE[i] = ZERO;
    AB2 = ZERO;
    for ( i=0; i<3; i++ ) {
	BA[i] = B[i]-A[i];
	AB2 += BA[i]*BA[i];
    }
    ips1 = *ips0 + (*nps_i);
    jps1 = *jps0 + (*nps_j);

    for ( ips=(*ips0); ips<ips1; ips++ ) {
	zeta_a = prim_exp[ips];
	coef_a = prim_coe[ips];
	for ( jps=(*jps0); jps<jps1; jps++ ) {
	    zeta_b = prim_exp[jps];
	    coef_b = prim_coe[jps];

	    zeta   = zeta_a + zeta_b;
	    sqrzi  = sqrt( 1.e0 / zeta );
	    zi     = sqrzi  * sqrzi;
	    coef   = coef_a * coef_b;

	    xiza  = zeta_b * zi;
	    xizb  = zeta_a * zi;
	    xi    = xiza * zeta_a;
	    xiab2 = xi * AB2;
	    xi2   = 2.e0 * xi;
	    zeta2 = HALF * zi;

	    exp_ab = coef * exp( -xiab2 );
	    oviss = _pi32_ * zi * sqrzi * exp_ab;
	    keiss = xi * ( 3.e0 - 2.e0 * xiab2 ) * oviss;
	    css   = _2pi_ * zi * exp_ab;

	    if ( fabs(css) < EPS_PS_PAIR_NAI ) continue;
	    for ( i=0; i<3; i++ ) {
		P[i]  = xiza*BA[i] + A[i];
		PA[i] = xiza*BA[i];
		PB[i] = xiza*BA[i] - BA[i];
	    }
	    // OVI and KEI
	    ovi[0] = oviss;
	    kei[0] = keiss;
	    // (x||s) and (x|T|s)
	    for ( ma=1; ma<=La; ma++ ) {
		iaop0 = LAOT[ma];
		iaop1 = iaop0 + NNAO[ma];
		iao0  = LAOT[ma-1];
		if ( ma>=2 ) iaom0 = LAOT[ma-2];
		for ( iaop=iaop0; iaop<iaop1; iaop++ ) {
		    ix  = INDX[iaop];
		    iao = NAM[iaop][ix];
		    nia = ANGM[iao][ix];
		    ip0 = OVIADD[ma  ][0] + (iaop - iaop0);
		    i00 = OVIADD[ma-1][0] + (iao  - iao0 );
		    ovi[ip0] = PA[ix]*ovi[i00];
		    kei[ip0] = PA[ix]*kei[i00];
		    if ( nia > 0 ) {
			iaom = NAM[iao][ix];
			i10  = OVIADD[ma-2][0] + (iaom - iaom0);
			ovi[ip0] += nia * zeta2 * ovi[i10];
			kei[ip0] += nia * (zeta2*kei[i10]-xiza*ovi[i10]);
		    }
		    kei[ip0] += xi2 * ovi[ip0];
		}
	    }	// for ( ma )
	    // (x||y) and (x|T|y)
	    for ( mb=1; mb<=Lb; mb++ ) {
		lmin = La - Lb + mb;
		jaop0 = LAOT[mb];
		jao0  = LAOT[mb-1];
		jaop1 = jaop0 + NNAO[mb];
		if ( lmin < 0 ) lmin = 0;
		if ( mb>=2 ) jaom0 = LAOT[mb-2];
		for ( ma=lmin; ma<=La; ma++ ) {
		    iao0 = LAOT[ma];
		    iao1 = iao0 + NNAO[ma];
		    if ( ma>=1 ) iaom0 = LAOT[ma-1];
		    for ( iao=iao0; iao<iao1; iao++ ) {
			for ( jaop=jaop0; jaop<jaop1; jaop++ ) {
			    ix  = INDX[jaop];
			    jao = NAM[jaop][ix];
			    nib = ANGM[jao][ix];
			    nia = ANGM[iao][ix];
			    i0p = OVIADD[ma][mb]
				+ (iao-iao0)*NNAO[mb]+(jaop-jaop0);
			    i00 = OVIADD[ma][mb-1]
				+ (iao-iao0)*NNAO[mb-1]+(jao-jao0);
			    ovi[i0p] = PB[ix]*ovi[i00];
			    kei[i0p] = PB[ix]*kei[i00];
			    if ( nib > 0 ) {
				jaom = NAM[jao][ix];
				i01  = OVIADD[ma][mb-2]
				    + (iao-iao0)*NNAO[mb-2] + (jaom-jaom0);
				ovi[i0p] += nib * zeta2 * ovi[i01];
				kei[i0p] +=
				    nib * (zeta2*kei[i01]-xizb*ovi[i01]);
			    }
			    if ( nia > 0 ) {
				iaom = NAM[iao][ix];
				i10  = OVIADD[ma-1][mb-1]
				    + (iaom-iaom0)*NNAO[mb-1]+(jao-jao0);
				ovi[i0p] += nia * zeta2 * ovi[i10];
				kei[i0p] += nia * zeta2 * kei[i10];
			    }
			    kei[i0p] += xi2 * ovi[i0p];
			}	// for (jao)
		    }	// for ( iao )
		}	// for ( ma )
	    }	// for ( mb );
	    // NAI
	    for ( kat=0; kat<(*nat); kat++ ) {
		dq = -(double)atomic_number[kat];
		PC[0] = P[0]-atom_x[kat];
		PC[1] = P[1]-atom_y[kat];
		PC[2] = P[2]-atom_z[kat];
		PC2 = PC[0]*PC[0] + PC[1]*PC[1] + PC[2]*PC[2];
		U = zeta * PC2;
		fmt( &nai[0], Lab, U, dq*css );
		// (x|A|s) type
		for ( ma=1; ma<=La; ma++ ) {
		    mmax  = Lab - ma;
		    iaop0 = LAOT[ma];
		    iaop1 = iaop0 + NNAO[ma];
		    iao0  = LAOT[ma-1];
		    if ( ma>=2 ) iaom0 = LAOT[ma-2];
		    for ( m=0; m<=mmax; m++ ) {
			for ( iaop=iaop0; iaop<iaop1; iaop++ ) {
			    ix   = INDX[iaop];
			    iao  = NAM[iaop][ix];
			    nia  = ANGM[iao][ix];
			    ip00 = NAIADD[ma][0]
				+ m*NNAO[ma]   + (iaop-iaop0);
			    i000 = NAIADD[ma-1][0]
				+ m*NNAO[ma-1] + (iao - iao0);
			    i001 = i000 + NNAO[ma-1];
			    nai[ip00] = PA[ix]*nai[i000]-PC[ix]*nai[i001];
			    if ( nia > 0 ) {
				iaom = NAM[iao][ix];
				i100 = NAIADD[ma-2][0]
				    + m*NNAO[ma-2] + (iaom-iaom0);
				i101 = i100 + NNAO[ma-2];
				nai[ip00] +=
				    nia * zeta2 * (nai[i100]-nai[i101]);
			    }
			}	// for (iaop)
		    }	// for (m)
		}	// for (ma)
		// (x|A|y) type
		for ( mb=1; mb<=Lb; mb++ ) {
		    lmin = La - Lb + mb;
		    if ( lmin < 0 ) lmin = 0;
		    mmax = Lb - mb;
		    jaop0 = LAOT[mb];
		    jaop1 = jaop0 + NNAO[mb];
		    jao0  = LAOT[mb-1];
		    if ( mb >= 2 ) jaom0 = LAOT[mb-2];
		    for ( ma=lmin; ma<=La; ma++ ) {
			iao0 = LAOT[ma];
			iao1 = iao0 + NNAO[ma];
			if ( ma > 0 ) iaom0 = LAOT[ma-1];
			for ( m=0; m<=mmax; m++ ) {
			    for ( iao=iao0; iao<iao1; iao++ ) {
				for ( jaop=jaop0; jaop<jaop1; jaop++ ) {
				    ix  = INDX[jaop];
				    jao = NAM[jaop][ix];
				    nib = ANGM[jao][ix];
				    nia = ANGM[iao][ix];
				    i0p0 =
					NAIADD[ma][mb]+m*MINDEX[ma][mb] +
					(iao-iao0)*NNAO[mb] + (jaop-jaop0);
				    i000 =
					NAIADD[ma][mb-1]+m*MINDEX[ma][mb-1]+
					(iao-iao0)*NNAO[mb-1] + (jao-jao0);
				    i001 = i000 + MINDEX[ma][mb-1];
				    nai[i0p0] =
					PB[ix]*nai[i000] - PC[ix]*nai[i001];
				    if ( nib > 0 ) {
					jaom = NAM[jao][ix];
					i010 = NAIADD[ma][mb-2]
					    + m*MINDEX[ma][mb-2]
					    + (iao-iao0)*NNAO[mb-2]
					    + (jaom-jaom0);
					i011 = i010 + MINDEX[ma][mb-2];
					nai[i0p0] +=
					    nib * zeta2 * ( nai[i010] - nai[i011] );
				    }
				    if ( nia > 0 ) {
					iaom = NAM[iao][ix];
					i100 = NAIADD[ma-1][mb-1]
					    + m * MINDEX[ma-1][mb-1]
					    + (iaom-iaom0)*NNAO[mb-1]
					    + (jao-jao0);
					i101 = i100 + MINDEX[ma-1][mb-1];
					nai[i0p0] +=
					    nia * zeta2 * ( nai[i100] - nai[i101] );
				    }
				}	// for (jaop)
			    }	// for (iao)
			}	// for (m)
		    }	// for (ma)
		}	// for (mb)
		// contraction of NAI
		p = &nai[NAIADD[La][Lb]];
		for ( i=0; i<nab; i++ ) HCORE[i] += p[i];
	    }	// for (kat)
	    // contraction of OVI and KEI
	    ov = &ovi[OVIADD[La][Lb]];
	    ke = &kei[OVIADD[La][Lb]];
	    for ( i=0; i<nab; i++ ) {
		OVI[i]   += ov[i];
		HCORE[i] += ke[i];
	    }
	}	// for ( jps )
    }	// for ( ips )
    // multiply coefficients
    iao0 = LAOT[La];
    iao1 = iao0 + na;
    jao0 = LAOT[Lb];
    jao1 = jao0 + nb;
    ij   = 0;
    //for ( iao=iao0, ij=0; iao<iao1; iao++ ) {
    for ( iao=iao0; iao<iao1; iao++ ) {
	coe_a = DFACT[iao];
	for ( jao=jao0; jao<jao1; jao++ ) {
	    coe = coe_a * DFACT[jao];
	    // debug
	    //if ( coe != 1.e0 ) printf("coe= %10.5f\n", coe );

	    OVI[ij]   *= coe;
	    HCORE[ij] *= coe;
	    ij++;
	}
    }
    return 0;
}

int ofmo_oneint_xx(
	const int *pnworkers, const int *pworkerid,
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[],
	const int shel_add[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	const int *pnat, const int atomic_number[],
	double S[], double H[] ) {
    int ij, iao2, ijao;
    int ips0, ics, ics0, ics1, iat, iao, iao0, iao1;
    int jps0, jcs, jcs0, jcs1, jat, jao, jao0, jao1, jcs_max;
    int nps_i, nps_j;
    double A[3], B[3];
    double *ovixx, *hcorexx;
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, nat=*pnat, dum1, dum2;
    //
    int na, nb, mythread;
    mythread = omp_get_thread_num();

    ovixx   = OVI_MASTER[mythread];
    hcorexx = HCORE_MASTER[mythread];

    ofmo_oneint_gen_make_add( mythread, La, Lb, &dum1, &dum2 );
    na = NNAO[La];
    nb = NNAO[Lb];
    ics0 = leading_cs[La];
    jcs0 = leading_cs[Lb];
    ics1 = leading_cs[La+1];
    jcs1 = leading_cs[Lb+1];
    for ( ics=ics0+workerid; ics<ics1; ics+=nworkers ) {
	jcs_max = ( La==Lb ? ics+1 : jcs1 );
	ips0  = shel_add[ics];
	iat   = shel_atm[ics];
	iao0  = shel_ini[ics];
	nps_i = shel_tem[ics];
	A[0]=atom_x[ iat ]; A[1]=atom_y[ iat ]; A[2]=atom_z[ iat ];
	for ( jcs=jcs0; jcs<jcs_max; jcs++ ) {
	    jps0  = shel_add[jcs];
	    jat   = shel_atm[jcs];
	    jao0  = shel_ini[jcs];
	    nps_j = shel_tem[jcs];
	    B[0]=atom_x[ jat ]; B[1]=atom_y[ jat ]; B[2]=atom_z[ jat ];
	    oneint_core_xx(
		    mythread, &La, &Lb,
		    &ips0, &nps_i, A,
		    &jps0, &nps_j, B,
		    prim_exp, prim_coe,
		    &nat, atom_x, atom_y, atom_z, atomic_number,
		    ovixx, hcorexx );
	    iao1 = iao0 + na;
	    jao1 = jao0 + nb;
	    ij = -1;
	    for ( iao=iao0; iao<iao1; iao++ ) {
		iao2 = (iao*iao+iao)>>1;
		for ( jao=jao0; jao<jao1; jao++ ) {
		    ij++;
		    if ( jao>iao ) continue;
		    ijao    = iao2 + jao;
		    S[ijao] = ovixx[ij];
		    H[ijao] = hcorexx[ij];
		}
	    }
	}	// for (jcs)
    }	// for (ics)
    return 0;
}
