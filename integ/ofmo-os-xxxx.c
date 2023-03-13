#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ofmo-index.h"
#include "ofmo-twoint.h"

#ifdef _OPENMP
#include <omp.h>
#else
#include "omp-dummy.h"
#endif

extern void fmt( double F[],
	const int m, const double T, const double cssss );

extern int ofmo_integ_add_fock( const int nao, const size_t nstored_eri,
	const double eri_val[], const short int eri_ind4[],
	const double D[], double G[] );

extern int ofmo_twoint_core_rys_xxxx(
	const int mythread,
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double DINT[] );

#ifndef false
#define false 0
#endif

#ifndef true
#define true 1
#endif

#define HALF	0.5e0
#define ZERO	0.e0
#define ONE	1.e0


#define OFMO_EBUF_FULL		1
#define OFMO_EBUF_NOFULL	0

#define MAXNPSPAIR 100
#define EPS_PS_PAIR	1.e-32
#define EPS_CS_PAIR2	1.e-30


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

// 4次元整数配列の確保
static int**** ofmo_alloc_i4d( int na, int nb, int nc, int nd ) {
    int ****ip, i, j, k;
    ip = (int****)malloc( sizeof(int***) * na );
    ip[0] = (int***)malloc( sizeof(int**) * na * nb );
    ip[0][0] = (int**)malloc( sizeof(int*) * na * nb * nc );
    ip[0][0][0] = (int*   )malloc( sizeof(int   ) * na * nb * nc * nd );

    for ( i=1; i<na; i++ ) ip[i]       = ip[i-1]       + nb;

    for ( j=1; j<nb; j++ ) ip[0][j]    = ip[0][j-1]    + nc;
    for ( i=1; i<na; i++ ) {
	for ( j=0; j<nb; j++ ) ip[i][j] = ip[i-1][j] + nb * nc;
    }

    for ( k=1; k<nc; k++ ) ip[0][0][k] = ip[0][0][k-1] + nd;
    for ( j=1; j<nb; j++ ) {
	for ( k=0; k<nc; k++ ) ip[0][j][k] = ip[0][j-1][k] + nc * nd;
    }
    for ( i=1; i<na; i++ ) {
	for ( j=0; j<nb; j++ ) {
	    for ( k=0; k<nc; k++ )
		ip[i][j][k] = ip[i-1][j][k] + nb * nc * nd;
	}
    }

    return ip;
}

// 4次元整数配列の解放
static void ofmo_free_i4d( int**** ip ) {
    if ( ip ) {
	if ( ip[0] ) {
	    if ( ip[0][0] ) {
		if ( ip[0][0][0] ) free( ip[0][0][0] );
		free ( ip[0][0] );
	    }
	    free( ip[0] );
	}
	free( ip );
    }
}

/* ofmo-index.cで生成される変数 */
static int *NNAO;
static int *LAOT;
static int *INDX;
static int **ANGM;
static int **NAM;
static int **NAP;
static double *DFACT;

/* VRR関数群で使用する変数 */
static int ***V_VADD = NULL;
static int ***V_Mindex = NULL;
static int ***V_KLindex = NULL;
static int ***V_IJindex = NULL;
static int ***V_Mmin = NULL;
static double **V_ev = NULL;

/* HRR計算で使用する変数 */
static int *****V_HADD = NULL;
static double **V_eh = NULL;

/* 縮約分子積分の格納場所 */
static double **DINTEG_MASTER = NULL;

/* カットオフテーブル関連 */
static double _CK_;

// VRR計算で必要となるアドレス情報などを設定する
static int ofmo_vrr_make_add(
	const int mythread,
	const int La, const int Lb, const int Lc, const int Ld,
	const int mmajor
	) {
    int vrr_mem, Lab, Lcd, Labcd, mab, mcd;
    int imin, mmin, mmax, nm;
    int **Mindex, **KLindex, **IJindex, **Mmin, **VADD;
    Mindex = V_Mindex[mythread];
    KLindex = V_KLindex[mythread];
    IJindex = V_IJindex[mythread];
    Mmin    = V_Mmin[mythread];
    VADD    = V_VADD[mythread];
    Lab = La + Lb;
    Lcd = Lc + Ld;
    Labcd = Lab + Lcd;
    vrr_mem = 0;
    // (xs,ss) type
    // m*Mmindex[mab][mcd] + (iao-iao0)*IJindex[mab][mcd]
    mcd = 0;
    for ( mab=0; mab<=Lab; mab++ ) {
	if ( mmajor ) {
	    Mindex[mab][mcd]  = NNAO[mab];
	    KLindex[mab][mcd] = 0;
	    IJindex[mab][mcd] = 1;
	} else {
	    Mindex[mab][mcd]  = 1;
	    KLindex[mab][mcd] = 0;
	    IJindex[mab][mcd] = (Labcd - mab + 1);
	}
	Mmin[mab][mcd] = 0;
	VADD[mab][mcd] = vrr_mem;
	vrr_mem += (Labcd - mab + 1) * NNAO[mab];
    }
    // (xs,ys) type
    // m*Mindex[mab][mcd] + (kao-kao0)*KLindex[mab][mcd]
    //     + (iao-iao0)*IJindex[mab][mcd]
    for ( mcd=1; mcd<=Lcd; mcd++ ) {
	mmax = Lcd - mcd;
	imin = La - Lcd + mcd;
	if ( imin < 0 ) imin = 0;
	for ( mab=imin; mab<=Lab; mab++ ) {
	    mmin = La - mab;
	    if ( mmin < 0 ) mmin = 0;
	    Mmin[mab][mcd] = mmin;
	    nm = mmax - mmin + 1;
	    if ( mmajor ) {
		Mindex[mab][mcd]  = NNAO[mab] * NNAO[mcd];
		KLindex[mab][mcd] = 1;
		IJindex[mab][mcd] = NNAO[mcd];
	    } else {
		Mindex[mab][mcd]  = 1;
		KLindex[mab][mcd] = nm;
		IJindex[mab][mcd] = nm * NNAO[mcd];
	    }
	    VADD[mab][mcd] = vrr_mem;
	    vrr_mem += nm * NNAO[mab] * NNAO[mcd];
	}
    }
    return vrr_mem;
}

// HRR計算で必要となるアドレス情報などを設定する
static int ofmo_hrr_make_add(
	const int mythread,
	const int La, const int Lb, const int Lc, const int Ld ) {
    int hrr_mem, ma, mb, mc, md;
    int na, nb, nab, nabd;
    int Lab, Lcd;
    int ****HADD;
    HADD = V_HADD[mythread];
    Lab = La + Lb;
    Lcd = Lc + Ld;
    hrr_mem = 0;
    // VRRで生成される縮約積分のアドレス
    mb = md = 0;
    for ( ma=La; ma<=Lab; ma++ ) {
	na = NNAO[ma];
	for ( mc=Lc; mc<=Lcd; mc++ ) {
	    HADD[ma][mb][mc][md] = hrr_mem;
	    hrr_mem += ( na*NNAO[mc] );
	}
    }
    // ABに対するHRRのアドレス
    md = 0;
    for ( mb=1; mb<=Lb; mb++ ) {
	nb = NNAO[mb];
	for ( ma=La; ma<=(Lab-mb); ma++ ) {
	    nab = nb * NNAO[ma];
	    for ( mc=Lc; mc<=Lcd; mc++ ) {
		HADD[ma][mb][mc][md] = hrr_mem;
		hrr_mem += ( nab * NNAO[mc]);
	    }
	}
    }
    // CDに対するHRRのアドレス
    ma = La;
    mb = Lb;
    nab = NNAO[La]*NNAO[Lb];
    for ( md=1; md<=Ld; md++ ) {
	nabd = nab*NNAO[md];
	for ( mc=Lc; mc<=(Lcd-md); mc++ ) {
	    HADD[ma][mb][mc][md] = hrr_mem;
	    hrr_mem += nabd * NNAO[mc];
	}
    }
    return hrr_mem;
}

static void ofmo_vrr_finalize() {
    int i, nthreads;
    nthreads = omp_get_max_threads();
    for ( i=0; i<nthreads; i++ ) {
	if ( V_VADD[i]    ) ofmo_free_imatrix( V_VADD[i] );
	if ( V_Mmin[i]    ) ofmo_free_imatrix( V_Mmin[i] );
	if ( V_Mindex[i]  ) ofmo_free_imatrix( V_Mindex[i] );
	if ( V_IJindex[i] ) ofmo_free_imatrix( V_IJindex[i] );
	if ( V_KLindex[i] ) ofmo_free_imatrix( V_KLindex[i] );
	if ( V_ev[i]      ) free( V_ev[i] );
	if ( DINTEG_MASTER[i] ) free( DINTEG_MASTER[i] );
    }
    free( V_VADD );
    free( V_Mmin );
    free( V_Mindex );
    free( V_IJindex );
    free( V_KLindex );
    free( V_ev );
    free( DINTEG_MASTER );
    V_VADD    = NULL;
    V_Mmin    = NULL;
    V_Mindex  = NULL;
    V_IJindex = NULL;
    V_KLindex = NULL;
    V_ev      = NULL;
    DINTEG_MASTER = NULL;
}

/* VRR初期化関数
 * １回だけ呼び出す必要あり */
static int ofmo_vrr_init( const int maxlqn ) {
    int Lab;
    int nthreads;
    static int called = false;
    if ( called ) return 0;
    ofmo_index_init( 2*maxlqn );
    NNAO = ofmo_getadd_nnao();
    LAOT = ofmo_getadd_laot();
    ANGM = ofmo_getadd_angm();
    INDX = ofmo_getadd_indx();
    NAM  = ofmo_getadd_nam();

    Lab     = maxlqn + maxlqn;
    nthreads = omp_get_max_threads();
    V_VADD    = (int***)malloc( sizeof(int**) * nthreads );
    V_Mindex  = (int***)malloc( sizeof(int**) * nthreads );
    V_KLindex = (int***)malloc( sizeof(int**) * nthreads );
    V_IJindex = (int***)malloc( sizeof(int**) * nthreads );
    V_Mmin    = (int***)malloc( sizeof(int**) * nthreads );
    V_ev      = (double**)malloc( sizeof(double*) * nthreads );

    DINTEG_MASTER = (double**)malloc( sizeof(double*) * nthreads );
#pragma omp parallel
    {
	int mythread, vrr_mem, n, n4;
	n = NNAO[maxlqn];
	n4 = n*n*n*n;
	mythread = omp_get_thread_num();
	V_VADD[mythread]    = ofmo_alloc_imatrix( Lab+1, Lab+1 );
	V_Mmin[mythread]    = ofmo_alloc_imatrix( Lab+1, Lab+1 );
	V_Mindex[mythread]  = ofmo_alloc_imatrix( Lab+1, Lab+1 );
	V_KLindex[mythread] = ofmo_alloc_imatrix( Lab+1, Lab+1 );
	V_IJindex[mythread] = ofmo_alloc_imatrix( Lab+1, Lab+1 );
	vrr_mem = ofmo_vrr_make_add( mythread,
		maxlqn, maxlqn, maxlqn, maxlqn, false );
	V_ev[mythread] = (double*)malloc( sizeof(double) * vrr_mem );
	DINTEG_MASTER[mythread] = (double*)malloc( sizeof(double) * n4 );
    }
    atexit( ofmo_vrr_finalize );
    called = true;
    return 0;
}

static void ofmo_hrr_finalize() {
    int nthreads, i;
    nthreads = omp_get_max_threads();
    for ( i=0; i<nthreads; i++ ) {
	if ( V_eh[i]   ) free( V_eh[i] );
	if ( V_HADD[i] ) ofmo_free_i4d( V_HADD[i] );
    }
    free( V_eh );
    free( V_HADD );
    V_eh = NULL;
    V_HADD = NULL;
}

// HRRの初期化関数（１回だけ呼び出せばよい）
static int ofmo_hrr_init( const int maxlqn ) {
    int nthreads;
    static int called = false;
    if ( called ) return 0;
    nthreads = omp_get_max_threads();
    V_HADD = (int*****)malloc( sizeof(int****) * nthreads );
    V_eh   = (double**)malloc( sizeof(double*) * nthreads );
    NAP   = ofmo_getadd_nap();
    DFACT = ofmo_getadd_dfact();
    atexit( ofmo_hrr_finalize );
#pragma omp parallel
    {
	int mythread, hrr_mem;
	mythread = omp_get_thread_num();
	V_HADD[mythread] =
	    ofmo_alloc_i4d( 2*maxlqn+1, maxlqn+1, 2*maxlqn+1, maxlqn+1);
	hrr_mem = ofmo_hrr_make_add( mythread,
		maxlqn, maxlqn, maxlqn, maxlqn );
	V_eh[mythread] = (double*)malloc( sizeof(double) * hrr_mem );
    }
    called = true;
    return 0;
}

// 確保した配列のアドレスを返す関数
double* ofmo_os_getadd_eri( const int mythread ) {
    return DINTEG_MASTER[mythread];
}

double* ofmo_os_getadd_vrr( const int mythread ) {
    return V_ev[mythread];
}

double* ofmo_os_getadd_hrr( const int mythread ) {
    return V_eh[mythread];
}

// 1つの原始積分に対するVRR計算を行う（縮約はしていない）
static int ofmo_vrr_calc(
	const int mythread,
	const int La, const int Lb, const int Lc, const int Ld,
	const double T, const double cssss,
	const double zeta2, const double eta2, const double ze2,
	const double rz, const double re,
	const double PA[3], const double WP[3],
	const double QC[3], const double WQ[3] ) {
    int mab, mcd, mmax, m, ix;
    int Lab, Lcd, Labcd, Lmin;
    int iao, iao0, iao1, iaom, iaom0, iaomm, iaomm0;
    int jao, jao0, jao1, jaom, jaom0, jaomm, jaomm0;
    int Ip00, I000, I001, I100, I101, I011;
    int I0p0,             I010;
    int nia, nic;
    int **Mindex, **KLindex, **IJindex, **Mmin, **VADD;
    double *ev;
    Mindex = V_Mindex[mythread];
    KLindex = V_KLindex[mythread];
    IJindex = V_IJindex[mythread];
    Mmin    = V_Mmin[mythread];
    VADD    = V_VADD[mythread];
    ev      = V_ev[mythread];

    Lab = La + Lb;
    Lcd = Lc + Ld;
    Labcd = Lab + Lcd;
    // (ss,ss)
    fmt( &ev[VADD[0][0]], Labcd, T, cssss );
    // (xs,ss) (x>=p)
    for ( mab=1; mab<=Lab; mab++ ) {
	mmax   = Labcd - mab;
	iao0   = LAOT[mab];
	iao1   = iao0 + NNAO[mab];
	iaom0  = LAOT[mab-1];
	if ( mab > 1 ) iaomm0 = LAOT[mab-2];
	for ( m=0; m<=mmax; m++ ) {
	    for ( iao=iao0; iao<iao1; iao++ ) {
		ix = INDX[iao];
		iaom = NAM[iao][ix];
		nia  = ANGM[iaom][ix];

		Ip00 = VADD[mab][0]
		    + m*Mindex[mab][0] + (iao-iao0)*IJindex[mab][0];
		I000 = VADD[mab-1][0]
		    + m*Mindex[mab-1][0] + (iaom-iaom0)*IJindex[mab-1][0];
		I001 = I000 + Mindex[mab-1][0];
		ev[Ip00] = PA[ix]*ev[I000] + WP[ix]*ev[I001];
		if ( nia > 0 ) {
		    iaomm = NAM[iaom][ix];
		    I100 = VADD[mab-2][0]
			+ m*Mindex[mab-2][0]
			+ (iaomm-iaomm0)*IJindex[mab-2][0];
		    I101 = I100 + Mindex[mab-2][0];
		    ev[Ip00] += (double)nia*zeta2*(ev[I100] - rz*ev[I101] );
		}
	    }
	}
    }
    // (xs,ys) (y>=p)
    for ( mcd=1; mcd<=Lcd; mcd++ ) {
	mmax = Lcd - mcd;
	jao0 = LAOT[mcd];
	jao1 = jao0 + NNAO[mcd];
	jaom0 = LAOT[mcd-1];
	if ( mcd > 1 ) jaomm0 = LAOT[mcd-2];
	Lmin = La - (Lcd-mcd);
	if ( Lmin < 0 ) Lmin = 0;
	for ( mab=Lmin; mab<=Lab; mab++ ) {
	    iao0 = LAOT[mab];
	    iao1 = iao0 + NNAO[mab];
	    if ( mab>0 ) iaom0 = LAOT[mab-1];
	    for ( m=Mmin[mab][mcd]; m<=mmax; m++ ) {
		for ( iao=iao0; iao<iao1; iao++ ) {
		    for ( jao=jao0; jao<jao1; jao++ ) {
			ix = INDX[jao];
			jaom = NAM[jao][ix];
			I0p0 = VADD[mab][mcd]
			    + (m-Mmin[mab][mcd])*Mindex[mab][mcd]
			    + (jao-jao0)*KLindex[mab][mcd]
			    + (iao-iao0)*IJindex[mab][mcd];
			I000 = VADD[mab][mcd-1]
			    + (m-Mmin[mab][mcd-1])*Mindex[mab][mcd-1]
			    + (jaom-jaom0)*KLindex[mab][mcd-1]
			    + (iao-iao0)*IJindex[mab][mcd-1];
			I001 = I000 + Mindex[mab][mcd-1];
			ev[I0p0] = QC[ix]*ev[I000] + WQ[ix]*ev[I001];
			nic = ANGM[jaom][ix];
			if ( nic > 0 ) {
			    jaomm = NAM[jaom][ix];
			    I010 = VADD[mab][mcd-2]
				+ (m-Mmin[mab][mcd-2])*Mindex[mab][mcd-2]
				+ (jaomm - jaomm0) * KLindex[mab][mcd-2]
				+ (iao-iao0) * IJindex[mab][mcd-2];
			    I011 = I010 + Mindex[mab][mcd-2];
			    ev[I0p0] +=
				(double)nic*eta2*(ev[I010] - re*ev[I011]);
			}
			nia = ANGM[iao][ix];
			if ( nia > 0 ) {
			    iaom = NAM[iao][ix];
			    I101 = VADD[mab-1][mcd-1]
				+ (m+1-Mmin[mab-1][mcd-1])*Mindex[mab-1][mcd-1]
				+ (jaom-jaom0) * KLindex[mab-1][mcd-1]
				+ (iaom-iaom0) * IJindex[mab-1][mcd-1];
			    ev[I0p0] += (double)nia*ze2*ev[I101];
			}
		    }
		}
	    }
	}
    }
    return 0;
}


/** HRRを行う関数
 * 外部変数
 * LAOT[CS type] = CSに含まれる先頭AO番号
 * NNAO[CS type] = CSに含まれるAO数
 * INDX[AO type] = 添字
 * HADD[ma][mb][mc][md] = 各軌道量子数4重対の先頭アドレス
 * eh[] = HRRに関連する縮約積分保存に使用する配列
 * */
static int ofmo_hrr_calc(
	const int mythread,
	const int La, const int Lb, const int Lc, const int Ld,
	const double BA[3], const double DC[3]
	) {
    int ma, mb, mc, md;
    int Lab, Lcd;
    int ix;
    int iao, iao0, iao1, iaop, iaop0;
    int jao, jao0, jao1, jaom, jaom0;
    int kao, kao0, kao1, kaop, kaop0, k;
    int lao, lao0, lao1, laom, laom0;
    int add01, add10, add00;
    int I01, I10, I00, IJ01, IJ10, IJ00, IJK01, IJK00;
    int IJKL01, IJKL10, IJKL00;
    double *d01, *d10, *d00;
    int ****HADD;
    double *eh;
    HADD = V_HADD[mythread];
    eh   = V_eh[mythread];

    Lab = La + Lb;
    Lcd = Lc + Ld;
    // ABに対するHRR
    for ( mb=1; mb<=Lb; mb++ ) {
	jao0  = LAOT[mb];
	jao1  = jao0 + NNAO[mb];
	jaom0 = LAOT[mb-1];
	for ( ma=La; ma<=(Lab-mb); ma++ ) {
	    iao0  = LAOT[ma];
	    iao1  = iao0 + NNAO[ma];
	    iaop0 = LAOT[ma+1];
	    for ( mc=Lc; mc<=Lcd; mc++ ) {
		kao0 = LAOT[mc];
		kao1 = kao0 + NNAO[mc];
		add01 = HADD[ma  ][mb  ][mc][0];
		add10 = HADD[ma+1][mb-1][mc][0];
		add00 = HADD[ma  ][mb-1][mc][0];
		for ( iao=iao0; iao<iao1; iao++ ) {
		    I01 = add01 + (iao-iao0)*NNAO[mb  ]*NNAO[mc];
		    I00 = add00 + (iao-iao0)*NNAO[mb-1]*NNAO[mc];
		    for ( jao=jao0; jao<jao1; jao++ ) {
			ix   = INDX[jao];
			jaom = NAM[jao][ix];
			iaop = NAP[iao][ix];
			IJ01 = I01 + (jao-jao0)*NNAO[mc];
			IJ10 = add10 + (iaop-iaop0)*NNAO[mb-1]*NNAO[mc]
			    + (jaom-jaom0)*NNAO[mc];
			IJ00 = I00 + (jaom-jaom0)*NNAO[mc];
			d01 = &eh[IJ01];
			d10 = &eh[IJ10];
			d00 = &eh[IJ00];
			for ( kao=kao0, k=0; kao<kao1; kao++, k++ )
			    d01[k] = d10[k] - BA[ix]*d00[k];
		    }
		}
	    }	// for (mc)
	}	// for (ma)
    }	// for (mb);
    // CDに対するHRR
    ma = La;
    mb = Lb;
    iao0 = LAOT[ma];
    iao1 = iao0 + NNAO[ma];
    jao0 = LAOT[mb];
    jao1 = jao0 + NNAO[mb];
    for ( md=1; md<=Ld; md++ ) {
	lao0  = LAOT[md];
	lao1  = lao0 + NNAO[md];
	laom0 = LAOT[md-1];
	for ( mc=Lc; mc<=(Lcd-md); mc++ ) {
	    kao0  = LAOT[mc];
	    kao1  = kao0 + NNAO[mc];
	    kaop0 = LAOT[mc+1];
	    add01 = HADD[ma][mb][mc  ][md  ];
	    add10 = HADD[ma][mb][mc+1][md-1];
	    add00 = HADD[ma][mb][mc  ][md-1];
	    for ( iao=iao0; iao<iao1; iao++ ) {
		I01 = add01 + (iao-iao0)*NNAO[mb]*NNAO[mc  ]*NNAO[md  ];
		I10 = add10 + (iao-iao0)*NNAO[mb]*NNAO[mc+1]*NNAO[md-1];
		I00 = add00 + (iao-iao0)*NNAO[mb]*NNAO[mc  ]*NNAO[md-1];
		for ( jao=jao0; jao<jao1; jao++ ) {
		    IJ01 = I01 + (jao-jao0)*NNAO[mc  ]*NNAO[md  ];
		    IJ10 = I10 + (jao-jao0)*NNAO[mc+1]*NNAO[md-1];
		    IJ00 = I00 + (jao-jao0)*NNAO[mc  ]*NNAO[md-1];
		    for ( kao=kao0; kao<kao1; kao++ ) {
			IJK01 = IJ01 + (kao-kao0)*NNAO[md  ];
			IJK00 = IJ00 + (kao-kao0)*NNAO[md-1];
			for ( lao=lao0; lao<lao1; lao++ ) {
			    ix = INDX[lao];
			    laom = NAM[lao][ix];
			    kaop = NAP[kao][ix];
			    IJKL01 = IJK01 + (lao-lao0);
			    IJKL10 = IJ10
				+ (kaop-kaop0)*NNAO[md-1] + (laom-laom0);
			    IJKL00 = IJK00 + (laom-laom0);
			    eh[IJKL01] = eh[IJKL10] - DC[ix]*eh[IJKL00];
			}	// for (lao)
		    }	// for (kao)
		}	// for (jao)
	    }	// for (iao)
	}	// for (mc)
    }	// for (md)
    return 0;
}

static int ofmo_hrr_clear(
	const int mythread,
	const int La, const int Lb, const int Lc, const int Ld ) {
    int mab, mcd, Lab, Lcd;
    int nab, nabcd, i;
    int ****HADD;
    double *th, *eh;
    eh   = V_eh[mythread];
    HADD = V_HADD[mythread];
    Lab = La + Lb;
    Lcd = Lc + Ld;
    for ( mab=La; mab<=Lab; mab++ ) {
	nab = NNAO[mab];
	for ( mcd=Lc; mcd<=Lcd; mcd++ ) {
	    nabcd = nab * NNAO[mcd];
	    th = &eh[ HADD[mab][0][mcd][0] ];
	    for ( i=0; i<nabcd; i++ ) th[i] = 0.e0;
	}
    }
    return 0;
}

static int ofmo_vrr_cint(
	const int mythread,
	const int La, const int Lb, const int Lc, const int Ld ) {
    int mab, mcd, Lab, Lcd;
    int nab, nabcd, i;
    double *tv, *th;
    double *ev, *eh;
    int **VADD, ****HADD;
    ev = V_ev[mythread];
    eh = V_eh[mythread];
    VADD = V_VADD[mythread];
    HADD = V_HADD[mythread];

    Lab = La + Lb;
    Lcd = Lc + Ld;
    for ( mab=La; mab<=Lab; mab++ ) {
	nab = NNAO[mab];
	for ( mcd=Lc; mcd<=Lcd; mcd++ ) {
	    nabcd = nab * NNAO[mcd];
	    tv = &ev[ VADD[mab][mcd] ];
	    th = &eh[ HADD[mab][0][mcd][0] ];
	    for ( i=0; i<nabcd; i++ ) th[i] += tv[i];
	}
    }
    return 0;
}

static int ofmo_hrr_coef(
	const int mythread,
	const int La, const int Lb, const int Lc, const int Ld,
	double DINT[] ) {
    int i, j, k, l, iao, jao, kao, lao, ix;
    double *th, coef_a, coef_ab, coef_abc;
    //
    double *eh;
    int ****HADD;
    eh = V_eh[mythread];
    HADD = V_HADD[mythread];
    th = &eh[ HADD[La][Lb][Lc][Ld] ];
    ix = 0;
    for ( i=0, iao=LAOT[La]; i<NNAO[La]; i++, iao++ ) {
	coef_a = DFACT[iao];
	for ( j=0, jao=LAOT[Lb]; j<NNAO[Lb]; j++, jao++ ) {
	    coef_ab = coef_a * DFACT[jao];
	    for ( k=0, kao=LAOT[Lc]; k<NNAO[Lc]; k++, kao++ ) {
		coef_abc = coef_ab * DFACT[kao];
		for ( l=0, lao=LAOT[Ld]; l<NNAO[Ld]; l++, lao++ ) {
		    DINT[ix] = coef_abc * DFACT[lao] * th[ix];
		    ix++;
		}
	    }
	}
    }
    return 0;
}

static int ofmo_twoint_core_xxxx(
	const int mythread,
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double DINT[] ) {
    int ijps, klps, i;
    double cssss, zeta, dkab, xiza, eta, xizc, dk, T;
    double zeta2, eta2, ze2, rz, re, PA[3], WP[3], QC[3], WQ[3];
    double PQ2, sqrho, rho, PC[3], QP[3];
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;

    ofmo_hrr_clear( mythread, La, Lb, Lc, Ld );
    for ( ijps=0; ijps<(*nijps); ijps++ ) {
	zeta = vzeta[ijps];
	dkab = vdkab[ijps];
	xiza = vxiza[ijps];
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
	    ze2   = re  * zeta2;
	    for ( i=0; i<3; i++ ) {
		WP[i] = rz*QP[i];
		WQ[i] = rz*QP[i] - QP[i];
	    }
	    T     = rho * PQ2;
	    cssss = sqrho * dk;
	    ofmo_vrr_calc( mythread, La, Lb, Lc, Ld,
		    T, cssss, zeta2, eta2, ze2, rz, re, PA, WP, QC, WQ );
	    ofmo_vrr_cint( mythread, La, Lb, Lc, Ld );
	}
    }
    ofmo_hrr_calc( mythread, La, Lb, Lc, Ld, BA, DC );
    ofmo_hrr_coef( mythread, La, Lb, Lc, Ld, DINT );
    return 0;
}

static int schwarz_calc_ps_pair_params(
	const double prim_exp[], const double prim_coe[],
	const int ips0, const int ips1,
	const int jps0, const int jps1,
	const double AB2,
	double vzeta[], double vdkps[], double vxiza[] ) {
    double zeta_a, coef_a, zeta_b, coef_b;
    double sqrz, rz, xiza, Kab, zeta, coef;
    int npps, ips, jps;
    npps = 0;
    for ( ips=ips0; ips<ips1; ips++ ) {
	zeta_a = prim_exp[ips];
	coef_a = prim_coe[ips];
	for ( jps=jps0; jps<jps1; jps++ ) {
	    zeta_b = prim_exp[jps];
	    coef_b = prim_coe[jps];

	    zeta = zeta_a + zeta_b;
	    coef = coef_a * coef_b;

	    sqrz = 1.e0 / sqrt( zeta );
	    rz   = sqrz * sqrz;
	    xiza = zeta_b * rz;
	    Kab  = _CK_ * coef * sqrz * rz * exp( -xiza * zeta_a * AB2 );
	    if ( fabs(Kab) > EPS_PS_PAIR ) {
		vzeta[npps] = rz;
		vdkps[npps] = Kab;
		vxiza[npps] = xiza;
		npps++;
	    }
	}
    }
    return npps;
}

static double ofmo_schwarz_core_xxxx(
	const int mythread,
	const int *pLa, const int *pLb,
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3], const double AB2 ) {
    int ijps, klps, i;
    double cssss, zeta, dkab, xiza, eta, xizc, dk, T, dxi;
    double zeta2, eta2, ze2, rz, re, PA[3], WP[3], QC[3], WQ[3];
    double PQ2, sqrho, rho, QP[3];
    int La=*pLa, Lb=*pLb, Lc=*pLa, Ld=*pLb;
    double *DINTEG;
    int na, nb, nab1, j, ij, i0;
    double dmaxint;

    na   = NNAO[La];
    nb   = NNAO[Lb];
    nab1 = na*nb+1;

    DINTEG = DINTEG_MASTER[mythread];

    ofmo_hrr_clear( mythread, La, Lb, Lc, Ld );

    for ( ijps=0; ijps<(*nijps); ijps++ ) {
	zeta = vzeta[ijps];
	dkab = vdkab[ijps];
	xiza = vxiza[ijps];
	zeta2 = HALF * zeta;
	for ( i=0; i<3; i++ ) PA[i] = xiza * BA[i];
	for ( klps=0; klps<(*nijps); klps++ ) {
	    eta  = vzeta[klps];
	    dk   = dkab * vdkab[klps];
	    xizc = vxiza[klps];
	    dxi  = xizc - xiza;
	    eta2 = HALF * eta;
	    PQ2  = dxi*dxi*AB2;
	    for ( i=0; i<3; i++ ) {
		QC[i] = xizc* BA[i];
		QP[i] = dxi * BA[i];
	    }
	    sqrho = sqrt(1.e0/(zeta+eta));
	    rho   = sqrho*sqrho;
	    rz    = rho * zeta;
	    re    = rho * eta;
	    ze2   = re  * zeta2;
	    for ( i=0; i<3; i++ ) {
		WP[i] = rz*QP[i];
		WQ[i] = rz*QP[i] - QP[i];
	    }
	    T     = rho * PQ2;
	    cssss = sqrho * dk;
	    ofmo_vrr_calc( mythread, La, Lb, Lc, Ld,
		    T, cssss, zeta2, eta2, ze2, rz, re, PA, WP, QC, WQ );
	    ofmo_vrr_cint( mythread, La, Lb, Lc, Ld );
	}
    }

    ofmo_hrr_calc( mythread, La, Lb, Lc, Ld, BA, BA );
    ofmo_hrr_coef( mythread, La, Lb, Lc, Ld, DINTEG );
    //
    dmaxint=0.e0;
    for ( i=0; i<na; i++ ) {
	i0 = i*nb*nab1;
	for ( j=0; j<nb; j++ ) {
	    ij = i0 + j*nab1;
	    if ( fabs(DINTEG[ij]) > dmaxint ) dmaxint = fabs(DINTEG[ij]);
	}
    }
    return dmaxint;
}

int ofmo_cutoff_xx(
	// input arguments
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[], const int shel_add[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	// output arguments
	int leading_cs_pair[],
	double csp_schwarz[], int csp_ics[], int csp_jcs[],
	int csp_leading_ps_pair[],
	double psp_zeta[], double psp_dkps[], double psp_xiza[] ) {
    int ics, ics0, ics1, ips0, ips1, iat;
    int jcs, jcs0, jcs1, jps0, jps1, jat;
    int npps, ipps, ncs_pair, nps_pair;
    int i, Lab, jcsmax;
    double BA[3], A[3], B[3], AB2;
    double max_eri;
    double vzeta[MAXNPSPAIR], vdkps[MAXNPSPAIR], vxiza[MAXNPSPAIR];
    int La=*pLa, Lb=*pLb, mythread;

    mythread  = omp_get_thread_num();
    ofmo_vrr_make_add( mythread, La, Lb, La, Lb, true );
    ofmo_hrr_make_add( mythread, La, Lb, La, Lb );

    Lab       = La*(La+1)/2 + Lb;
    ncs_pair  = leading_cs_pair[Lab];
    nps_pair  = csp_leading_ps_pair[ncs_pair];

    ics0 = leading_cs[La];	// La=2
    ics1 = leading_cs[La+1];
    jcs0 = leading_cs[Lb];	// Lb=2
    jcs1 = leading_cs[Lb+1];
    for ( ics=ics0; ics<ics1; ics++ ) {
	ips0 = shel_add[ics];
	ips1 = ips0 + shel_tem[ics];
	iat  = shel_atm[ics];
	A[0]=atom_x[ iat ]; A[1]=atom_y[ iat ]; A[2]=atom_z[ iat ];
	jcsmax = ( Lb==La ? ics+1 : jcs1 );
	for ( jcs=jcs0; jcs<jcsmax; jcs++ ) {
	    jps0 = shel_add[jcs];
	    jps1 = jps0 + shel_tem[jcs];
	    jat  = shel_atm[jcs];
	    B[0]=atom_x[ jat ]; B[1]=atom_y[ jat ]; B[2]=atom_z[ jat ];

	    AB2 = 0.e0;
	    for ( i=0; i<3; i++ ) {
		BA[i] = B[i] - A[i];
		AB2 += BA[i]*BA[i];
	    }
	    
	    npps = schwarz_calc_ps_pair_params(
		    prim_exp, prim_coe, ips0, ips1, jps0, jps1, AB2,
		    vzeta, vdkps, vxiza );
	    if ( npps == 0 ) continue;
	    
	    max_eri = ofmo_schwarz_core_xxxx(
		    mythread, pLa, pLb,
		    &npps, vzeta, vdkps, vxiza, BA, AB2 );

	    if ( max_eri > EPS_CS_PAIR2 ) {
		csp_schwarz[ncs_pair] = sqrt( max_eri );
		csp_ics[ncs_pair] = ics;
		csp_jcs[ncs_pair] = jcs;
		for ( ipps=0; ipps<npps; ipps++ ) {
		    psp_zeta[nps_pair] = vzeta[ipps];
		    psp_dkps[nps_pair] = vdkps[ipps];
		    psp_xiza[nps_pair] = vxiza[ipps];
		    nps_pair++;
		}
		csp_leading_ps_pair[ncs_pair+1] = nps_pair;
		ncs_pair++;
	    }
	}	// jcs
    }		// ics
    leading_cs_pair[Lab+1] = ncs_pair;
    return 0;
}

// 初期化関数（１回だけ呼び出せばよい）
int ofmo_OS_integ_init( const int maxlqn ) {
    static int called = false;
    if ( called ) return 0;
    ofmo_vrr_init( maxlqn );
    ofmo_hrr_init( maxlqn );
    //
    double pi, t;
    pi = 4.e0 * atan( 1.e0 );
    t = 2.e0 * pi * pi * sqrt( pi );
    _CK_ = sqrt( t );
    called = true;
    return 0;
}

// 縮約積分計算関数
int ofmo_twoint_xxxx(
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
    int klcs, klcs0, klcs1, max_klcs ;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe, coe0;
    double *DINTEG;
    long nzeri, max_nzeri, nzeri4;
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    long ebuf_max_nzeri = *pebuf_max_nzeri;
    //
    int na, nb, nc, nd;
    //
    int Labcd, lambda;
    //
    int mythread;
    mythread = omp_get_thread_num();
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);

    ofmo_vrr_make_add( mythread, La, Lb, Lc, Ld, true );
    ofmo_hrr_make_add( mythread, La, Lb, Lc, Ld );
    DINTEG = DINTEG_MASTER[mythread];

    na   = NNAO[La];
    nb   = NNAO[Lb];
    nc   = NNAO[Lc];
    nd   = NNAO[Ld];
    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    //
    Labcd = Lab*(Lab+1)/2 + Lcd;
    lambda = La + Lb + Lc + Ld;

    ijcs0 = leading_cs_pair[Lab];
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];

    nzeri     = *ebuf_non_zero_eri;
    max_nzeri = ebuf_max_nzeri - na*nb*nc*nd;
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

	    ofmo_twoint_core_xxxx( mythread,
		    &La, &Lb, &Lc, &Ld,
		    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
		    &psp_xiza[ijps0], BA,
		    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
		    &psp_xiza[klps0], DC,   AC,      DINTEG );

	    /*// debug
#pragma omp master
	    {
		if ( Lab==1 && Lcd==0 && fabs(DINTEG[0]) > 1.e-7 ) {
		    printf("ijcs, klcs= %4d, %4d E[0]= %10.7f\n",
			    ijcs, klcs, DINTEG[0] );
		    fflush(stdout);
		}
	    }*/

	    ipat = ((Lab != Lcd) || (ics==kcs && jcs>lcs) ? true : false);
#ifdef SORT_CSP
            int ijgekl = (ics>kcs);
            if (ics==kcs) ijgekl = (jcs>=lcs);
            if (!ijgekl) ipat = ( (ics==kcs && jcs<lcs) ? true : false);
#endif
	    for ( i=0, iao=iao0, ix=0; i<na; i++, iao++ ) {
		I2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<nb; j++, jao++ ) {
		    if ( jao>iao ) { ix+=nc*nd; continue; }
		    IJ = I2 + jao;
		    coe0 = ( iao==jao ? HALF : ONE );
		    for ( k=0, kao=kao0; k<nc; k++, kao++ ) {
			K2 = (kao*kao+kao)>>1;
			for ( l=0, lao=lao0; l<nd; l++, lao++, ix++ ) {
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
		    }		// k
		}		// j
	    }			// i
	    if ( nzeri >= max_nzeri ) {
		*last_ijcs = ijcs;
		*last_klcs = klcs;
		*ebuf_non_zero_eri = nzeri;
		return OFMO_EBUF_FULL;
	    }
	}	// for ( klcs );
    }		// for ( ijcs );
    *ebuf_non_zero_eri = nzeri;
    return OFMO_EBUF_NOFULL;
}
//
// 縮約積分計算関数
int ofmo_twoint_direct_xxxx(
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
    int klcs, klcs0, klcs1, max_klcs ;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe, coe0;
    double *DINTEG;
    //
    int na, nb, nc, nd;
    //
    int mythread;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);
    float eps_sch = ofmo_twoint_eps_sch(0);

    mythread = omp_get_thread_num();
    ofmo_vrr_make_add( mythread, La, Lb, Lc, Ld, true );
    ofmo_hrr_make_add( mythread, La, Lb, Lc, Ld );
    DINTEG = DINTEG_MASTER[mythread];

    na   = NNAO[La];
    nb   = NNAO[Lb];
    nc   = NNAO[Lc];
    nd   = NNAO[Ld];

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

    max_nzeri -= na*nb*nc*nd;
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

	    ofmo_twoint_core_xxxx( mythread,
		    &La, &Lb, &Lc, &Ld,
		    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
		    &psp_xiza[ijps0], BA,
		    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
		    &psp_xiza[klps0], DC,   AC,      DINTEG );

	    ipat = ((Lab != Lcd) || (ics==kcs && jcs>lcs) ? true : false);
#ifdef SORT_CSP
            int ijgekl = (ics>kcs);
            if (ics==kcs) ijgekl = (jcs>=lcs);
            if (!ijgekl) ipat = ( (ics==kcs && jcs<lcs) ? true : false);
#endif
	    for ( i=0, iao=iao0, ix=0; i<na; i++, iao++ ) {
		I2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<nb; j++, jao++ ) {
		    if ( jao>iao ) { ix+=nc*nd; continue; }
		    IJ = I2 + jao;
		    coe0 = ( iao==jao ? HALF : ONE );
		    for ( k=0, kao=kao0; k<nc; k++, kao++ ) {
			K2 = (kao*kao+kao)>>1;
			for ( l=0, lao=lao0; l<nd; l++, lao++, ix++ ) {
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
		    }		// k
		}		// j
	    }			// i
	    if ( nzeri >= max_nzeri ) {
		ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4,
			Ds, G );
		nzeri = nzeri4= 0;
	    }
	}	// for ( klcs );
	klcs = klcs0;
    }		// for ( ijcs );
    *petmp_non_zero_eri = nzeri;
    return 0;
}
