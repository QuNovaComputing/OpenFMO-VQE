#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#else
#include "omp-dummy.h"
#endif

#include "ofmo-index.h"
#include "ofmo-twoint.h"

#ifndef false
#define false 0
#endif
#ifndef true
#define true 1
#endif

#define HALF	0.5e0
#define ONE	1.e0
#define ZERO	0.e0

#define EPS_PS4 1.e-30
#define EPS_ERI 1.e-15

#define OFMO_EBUF_FULL		1
#define OFMO_EBUF_NOFULL	0

#define MAXNPSPAIR 100
#define EPS_PS_PAIR	1.e-32
#define EPS_CS_PAIR2	1.e-30

extern void calc_root( const int nroot, const double T,
	double *U, double *W );

extern int ofmo_integ_add_fock( const int nao, const size_t nstored_eri,
	const double eri_val[], const short int eri_ind4[],
	const double D[], double G[] );


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

static int *NNAO;
static int **ANGM;
static int *LAOT;
static int *INDX;
static int **NAM;
static int **NAP;
static double *DFACT;

/* Rys積分で用いる変数 */
static double **V_XINT;
static double **V_YINT;
static double **V_ZINT;

/* HRR計算で使用する変数 */
static int *****V_HADD = NULL;
static double **V_eh = NULL;

/* indx */
static int *NROOTS;
static int **INS;

/* 縮約分子積分の格納場所 */
static double **DINTEG_MASTER = NULL;

/* カットオフテーブル関連 */
static double _CK_;

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

static void ofmo_rys_finalize() {
    int nthreads, i;
    nthreads = omp_get_max_threads();
    for ( i=0; i<nthreads; i++ ) {
	free( V_XINT[i] );
	free( V_YINT[i] );
	free( V_ZINT[i] );
	free( INS[i] );
	free( DINTEG_MASTER[i] );
    }
    free( V_XINT );
    free( V_YINT );
    free( V_ZINT );
    free( INS );
    free( NROOTS );
    free( DINTEG_MASTER );
}

static int ofmo_rys_init( const int maxlqn ) {
    int nthreads, nroot, maxlqn2;
    maxlqn2 = 2*maxlqn;
    ofmo_index_init( maxlqn2 );
    NNAO  = ofmo_getadd_nnao();
    LAOT  = ofmo_getadd_laot();
    ANGM  = ofmo_getadd_angm();
    INDX  = ofmo_getadd_indx();
    NAM   = ofmo_getadd_nam();
    nroot = ((4*maxlqn)>>1) + 1;
    nthreads = omp_get_max_threads();
    V_XINT = (double**)malloc( sizeof(double*) * nthreads );
    V_YINT = (double**)malloc( sizeof(double*) * nthreads );
    V_ZINT = (double**)malloc( sizeof(double*) * nthreads );
    NROOTS = (int*)malloc( sizeof(int) * nthreads );
    INS    = (int**)malloc( sizeof(int*) * nthreads );
    DINTEG_MASTER = (double**)malloc( sizeof(double*) * nthreads );
#pragma omp parallel
    {
	int mythread, nint, n, n4;
	n  = NNAO[maxlqn];
	n4 = n*n*n*n;
	mythread = omp_get_thread_num();
	nint     = (maxlqn2+1)*(maxlqn2+1)*nroot;
	V_XINT[mythread] = (double*)malloc( sizeof(double) * nint );
	V_YINT[mythread] = (double*)malloc( sizeof(double) * nint );
	V_ZINT[mythread] = (double*)malloc( sizeof(double) * nint );
	INS[mythread]    = (int*)malloc( sizeof(int) * 3 );
	DINTEG_MASTER[mythread] = (double*)malloc( sizeof(double) * n4 );
    }
    atexit( ofmo_rys_finalize );
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
    return 0;
}


// 初期化関数（１回だけ呼び出せばよい）
int ofmo_Rys_integ_init( const int maxlqn ) {
    static int called = false;
    double pi, t;
    if ( called ) return 0;
    ofmo_rys_init( maxlqn );
    ofmo_hrr_init( maxlqn );
    pi = 4.e0 * atan( 1.e0 );
    t = 2.e0 * pi * pi * sqrt( pi );
    _CK_ = sqrt( t );
    called = true;
    return 0;
}

static int ofmo_hrr_clear(
	const int La, const int Lb, const int Lc, const int Ld,
	int ****HADD, double *eh ) {
    int mab, mcd, Lab, Lcd;
    int nab, nabcd, i;
    double *th;
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

static int ofmo_hrr_coef(
	const int La, const int Lb, const int Lc, const int Ld,
	double DINT[], 
	int ****HADD, double *eh ) {
    int i, j, k, l, iao, jao, kao, lao, ix;
    double *th, coef_a, coef_ab, coef_abc;
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

/* 確保した並列のアドレス取得に関する関数群 */
double* ofmo_integ_getadd_xint( const int mythread ) {
    return V_XINT[mythread];
}

double* ofmo_integ_getadd_yint( const int mythread ) {
    return V_YINT[mythread];
}

double* ofmo_integ_getadd_zint( const int mythread ) {
    return V_ZINT[mythread];
}

double* ofmo_integ_getadd_eh( const int mythread ) {
    return V_eh[mythread];
}

double* ofmo_integ_getadd_eri( const int mythread ) {
    return DINTEG_MASTER[mythread];
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
	const int La, const int Lb, const int Lc, const int Ld,
	const double BA[3], const double DC[3],
	int ****HADD, double *eh ) {
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

static void set_indx( const int mythread,
	const int La, const int Lb, const int Lc, const int Ld ) {
    int Lab, Lcd;
    Lab = La + Lb;
    Lcd = Lc + Ld;
    NROOTS[mythread] = ( (Lab+Lcd)>>1 ) + 1;
    INS[mythread][0] = 1;
    INS[mythread][1] = NROOTS[mythread] * (Lcd+1);
    INS[mythread][2] = NROOTS[mythread];
}

static int indx( const int mythread,
	const int m, const int N, const int M ) {
    return m*INS[mythread][0]+N*INS[mythread][1]+M*INS[mythread][2];
}

static void ofmo_form( const int mythread,
	const int La, const int Lb, const int Lc, const int Ld,
	const int nroot, const double *xint, const double *yint,
	const double *zint, int ****HADD, double *eh ) {
    int Lab, Lcd;
    int mab, mcd, iao, iao0, iao1, kao, kao0, kao1, m;
    int ix, iy, iz, kx, ky, kz;
    int IJKL, ncd;
    double *XSXS;
    Lab  = La + Lb;
    Lcd  = Lc + Ld;
    for ( mab=La; mab<=Lab; mab++ ) {
	iao0 = LAOT[mab];
	iao1 = iao0 + NNAO[mab];
	for ( mcd=Lc; mcd<=Lcd; mcd++ ) {
	    kao0 = LAOT[mcd];
	    ncd  = NNAO[mcd];
	    kao1 = kao0 + ncd;
	    XSXS = &eh[ HADD[mab][0][mcd][0] ];
	    for ( iao=iao0; iao<iao1; iao++ ) {
		ix = ANGM[iao][0];
		iy = ANGM[iao][1];
		iz = ANGM[iao][2];
		for ( kao=kao0; kao<kao1; kao++ ) {
		    kx = ANGM[kao][0];
		    ky = ANGM[kao][1];
		    kz = ANGM[kao][2];
		    IJKL = (iao-iao0)*ncd + (kao-kao0);
		    for ( m=0; m<nroot; m++ )
			XSXS[IJKL] += xint[ indx(mythread,m,ix,kx) ] *
			    yint[ indx(mythread,m,iy,ky) ] *
			    zint[ indx(mythread,m,iz,kz) ];
		}
	    }
	}
    }
}

static void ofmo_xyzint_v( const int mythread,
	const int La, const int Lb, const int Lc, const int Ld,
	const double *F00, const double *B00, const double *B10,
	const double *B01, const double *C00, const double *CP00,
	const int nroot, double *xint, double *yint, double *zint ) {
    int Lab, Lcd;
    int m, m3, N, M, ix3, ix2, ix1, ix0;
    double C10[13], CP10[13], CP01[13], C01[13];

    Lab   = La + Lb;
    Lcd   = Lc + Ld;

    // (0,0)
    for ( m=0; m<nroot; m++ ) {
	ix0 = indx( mythread, m, 0, 0 );
	xint[ix0] = 1.e0;
	yint[ix0] = 1.e0;
	zint[ix0] = F00[m];
    }
    // (1,0)
    if ( Lab > 0 ) {
	for ( m=m3=0; m<nroot; m++, m3+=3 ) {
	    ix0 = indx( mythread, m, 1, 0 );
	    xint[ix0] = C00[m3+0];
	    yint[ix0] = C00[m3+1];
	    zint[ix0] = F00[m]*C00[m3+2];
	}
    }
    // (0,1)
    if ( Lcd > 0 ) {
	for ( m=m3=0; m<nroot; m++, m3+=3 ) {
	    ix0 = indx( mythread, m, 0, 1 );
	    xint[ix0] = CP00[m3+0];
	    yint[ix0] = CP00[m3+1];
	    zint[ix0] = F00[m] * CP00[m3+2];
	}
    }
    // (1,1) = C'00*(1,0)+B00*(0,0)
    if ( Lab>0 && Lcd>0 ) {
	for ( m=m3=0; m<nroot; m++, m3+=3 ) {
	    ix1 = indx( mythread, m, 1, 1 );
	    ix0 = indx( mythread, m, 1, 0 );
	    xint[ix1] = CP00[m3+0]*xint[ix0] + B00[m];
	    yint[ix1] = CP00[m3+1]*yint[ix0] + B00[m];
	    zint[ix1] = CP00[m3+2]*zint[ix0] + B00[m]*F00[m];
	}
    }
    // (N,0) and (N,1) (N>=2)
    // (N,0) = C00 *(N-1,0) + (N-1)*B10*(N-2,0)
    // (N,1) = C'00*(N,0)   +     N*B00*(N-1,0)
    if ( Lab > 1 ) {
	for ( m=0; m<nroot; m++ ) {
	    C10[m]  = 0.e0;
	    CP10[m] = B00[m];
	}
	for ( N=2; N<=Lab; N++ ) {
	    for ( m=m3=0; m<nroot; m++, m3+=3 ) {
		C10[m] += B10[m];
		ix2 = indx( mythread, m, N  , 0 );
		ix1 = indx( mythread, m, N-1, 0 );
		ix0 = indx( mythread, m, N-2, 0 );
		xint[ix2]=C00[m3+0]*xint[ix1] + C10[m]*xint[ix0];
		yint[ix2]=C00[m3+1]*yint[ix1] + C10[m]*yint[ix0];
		zint[ix2]=C00[m3+2]*zint[ix1] + C10[m]*zint[ix0];
	    }
	    if ( Lcd>0 ) {
		for ( m=m3=0; m<nroot; m++, m3+=3 ) {
		    CP10[m] += B00[m];
		    ix2 = indx( mythread, m, N  , 1 );
		    ix1 = indx( mythread, m, N  , 0 );
		    ix0 = indx( mythread, m, N-1, 0 );
		    xint[ix2] = CP00[m3+0]*xint[ix1] + CP10[m]*xint[ix0];
		    yint[ix2] = CP00[m3+1]*yint[ix1] + CP10[m]*yint[ix0];
		    zint[ix2] = CP00[m3+2]*zint[ix1] + CP10[m]*zint[ix0];
		}
	    }
	}
    }	// if (Lab>1)
    // (0,M) and (1,M) (M>=2)
    // (0,M) = C'00*(0,M-1) + (M-1)*B'01*(0,M-2)
    // (1,M) = C00 *(0,M)   +     M* B00*(0,M-1)
    if ( Lcd > 1 ) {
	for ( m=0; m<nroot; m++ ) {
	    CP01[m] = 0.e0;
	    C01[m]  = B00[m];
	}
	for ( M=2; M<=Lcd; M++ ) {
	    for ( m=m3=0; m<nroot; m++, m3+=3 ) {
		CP01[m] += B01[m];
		ix2 = indx( mythread, m, 0, M );
		ix1 = indx( mythread, m, 0, M-1 );
		ix0 = indx( mythread, m, 0, M-2 );
		xint[ix2] = CP00[m3+0]*xint[ix1] + CP01[m]*xint[ix0];
		yint[ix2] = CP00[m3+1]*yint[ix1] + CP01[m]*yint[ix0];
		zint[ix2] = CP00[m3+2]*zint[ix1] + CP01[m]*zint[ix0];
	    }
	    if ( Lab>0 ) {
		for ( m=m3=0; m<nroot; m++, m3+=3 ) {
		    C01[m] += B00[m];
		    ix2 = indx( mythread, m, 1, M );
		    ix1 = indx( mythread, m, 0, M );
		    ix0 = indx( mythread, m, 0, M-1 );
		    xint[ix2] = C00[m3+0]*xint[ix1] + C01[m]*xint[ix0];
		    yint[ix2] = C00[m3+1]*yint[ix1] + C01[m]*yint[ix0];
		    zint[ix2] = C00[m3+2]*zint[ix1] + C01[m]*zint[ix0];
		}
	    }
	}
    }
    // (N,M) (N>=2 and M>=2)
    // (N,M) = C00*(N-1,M) + (N-1)*B10*(N-2,M) + M*B00*(N-1,M-1)
    if ( Lab>1 && Lcd>1 ) {
	for ( m=0; m<nroot; m++ ) C01[m] = B00[m];
	for ( M=2; M<=Lcd; M++ ) {
	    for ( m=0; m<nroot; m++ ) {
		C01[m] += B00[m];
		C10[m]  = B10[m];
	    }
	    for ( N=2; N<=Lab; N++ ) {
		for ( m=m3=0; m<nroot; m++, m3+=3 ) {
		    ix3 = indx( mythread, m, N  , M   );
		    ix2 = indx( mythread, m, N-1, M   );
		    ix1 = indx( mythread, m, N-2, M   );
		    ix0 = indx( mythread, m, N-1, M-1 );
		    xint[ix3] = C00[m3+0]*xint[ix2]+C10[m]*xint[ix1]
			+C01[m]*xint[ix0];
		    yint[ix3] = C00[m3+1]*yint[ix2]+C10[m]*yint[ix1]
			+C01[m]*yint[ix0];
		    zint[ix3] = C00[m3+2]*zint[ix2]+C10[m]*zint[ix1]
			+C01[m]*zint[ix0];
		    C10[m] += B10[m];
		}
	    }
	}
    }
}

static int ofmo_twoint_core_rys_xxxx(
	const int mythread,
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	const int *nijps, const double vzeta[], const double vdkab[],
	const double vxiza[], const double BA[3],
	const int *nklps, const double veta[], const double vdkcd[],
	const double vxizc[], const double DC[3], const double AC[3],
	double DINT[] ) {
    int ijps, klps, i;
    double cssss, zeta, dkab, xiza, eta, xizc, dk, T;
    double zeta2, eta2, rz, PA[3], QC[3];
    double PQ2, sqrho, rho, PC[3], QP[3];
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    //
    double C00[13*3], CP00[13*3], B00[13], B10[13], B01[13], F00[13];
    double rrho, rze, W[13], U[13];
    double u2, duminv, dm2inv, dum;
    int m, m3;
    //
    int ****HADD, nroot;
    double *xint, *yint, *zint, *eh;

    HADD = V_HADD[mythread];
    xint = V_XINT[mythread];
    yint = V_YINT[mythread];
    zint = V_ZINT[mythread];
    eh   = V_eh[mythread];
    nroot = NROOTS[mythread];

    ofmo_hrr_clear( La, Lb, Lc, Ld, HADD, eh );

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
	    rrho  = zeta + eta;
	    rze   = zeta * eta;
	    sqrho = sqrt(1.e0/rrho);
	    rho   = sqrho * sqrho;
	    rz    = rho * zeta;
	    T     = rho * PQ2;
	    cssss = sqrho * dk;
	    calc_root( nroot, T, U, W );

	    for ( m=m3=0; m<nroot; m++, m3+=3 ) {
		u2     = rho * U[m];
		F00[m]    = cssss * W[m];
		duminv = 1.e0 / ( 1.e0 + rrho * u2 );
		dm2inv = 0.5e0 * duminv;
		B00[m]    = dm2inv * rze * u2;
		B10[m]    = dm2inv * ( zeta + rze*u2 );
		B01[m]    = dm2inv * ( eta  + rze*u2 );
		dum    = zeta * u2 * duminv;
		for ( i=0; i<3; i++ ) C00[m3+i]  = PA[i] + dum * QP[i];
		dum    = eta * u2 * duminv;
		for ( i=0; i<3; i++ ) CP00[m3+i] = QC[i] - dum * QP[i];
	    }
	    ofmo_xyzint_v( mythread, La, Lb, Lc, Ld,
		    F00, B00, B10, B01, C00, CP00,
		    nroot, xint, yint, zint );
	    ofmo_form( mythread, La, Lb, Lc, Ld,
		    nroot, xint, yint, zint, HADD, eh );
	}
    }

    ofmo_hrr_calc( La, Lb, Lc, Ld, BA, DC, HADD, eh );
    ofmo_hrr_coef( La, Lb, Lc, Ld, DINT,   HADD, eh );
    return 0;
}

// 縮約積分計算関数
int ofmo_twoint_rys_xxxx(
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
    int mythread;
    mythread = omp_get_thread_num();

    ofmo_hrr_make_add( mythread, La, Lb, Lc, Ld );
    set_indx( mythread, La, Lb, Lc, Ld );
    DINTEG = DINTEG_MASTER[mythread];

    na   = NNAO[La];
    nb   = NNAO[Lb];
    nc   = NNAO[Lc];
    nd   = NNAO[Ld];
    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;

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

	    ofmo_twoint_core_rys_xxxx( mythread,
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
			    if ( fabs(DINTEG[ix]) > EPS_ERI ) {
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
int ofmo_twoint_direct_rys_xxxx(
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
    ofmo_hrr_make_add( mythread, La, Lb, Lc, Ld );
    set_indx( mythread, La, Lb, Lc, Ld );
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

	    ofmo_twoint_core_rys_xxxx( mythread,
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
