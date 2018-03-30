#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cudalib.h"
#include "cuda-fmt-m.h"
#include "cuda-fmt-drv.h"

#ifndef CUDA_FMT_M
int cuda_fmt_m_finalize(void) {return 0;};
int cuda_fmt_m_init(void) {return 0;};
#else

#ifndef false
#define false 0
#endif
#ifndef true 
#define true 1
#endif
#define T_TOL 90.e0

static double gau_slow( const int m, const double t ) {
    int k;
    double gau, coef, t2;
    double PI = 4.e0 * atan( 1.e0 );
    gau = 0.5e0 * sqrt( PI / t );
    if ( m > 0 ) {
        t2 = 0.5e0 / t;
        coef = 1.e0;
        for ( k=1; k<=m; k++ ) {
            gau *= (coef*t2);
            coef += 2.e0;
        }
    }
    return gau;
}

static double fmt_slow( const int m, const double t ) {
    double t2, fac, term, fmt, expt, thr_zero=1.e-17, eps;
    if ( t >= T_TOL ) return gau_slow( m, t );
    expt = exp( -t );
    t2  = 2.e0 * t;
    eps = (expt/t2) * thr_zero;
    fac = (double)(2*m+1);
    term = 1.e0 / fac;
    fmt  = term;
    while (1) {
        fac += 2.e0;
        term *= (t2/fac);
        fmt += term;
        if ( term < eps ) break;
    }
    return fmt*expt; 
}

static double* fmt_make_table_method1( const int m,
        const int tmax, const int nexp, const int ndiv ) {
    double *table_p, *table_i;
    int step, i, k;
    double d, d2, t0;
    size_t table_size;
    d      = 1.e0 / (double)ndiv;
    d2     = 0.5e0 * d;
    step       = m+ nexp;
    table_size = step*tmax*ndiv*sizeof(double);
    table_p = (double*)malloc( table_size );
    for ( i=0, table_i=table_p; i<(tmax*ndiv); i++, table_i+=step ) {
        t0 = d*(double)i + d2;
        for ( k=0; k<step; k++ ) table_i[k] = fmt_slow( k, t0 );
    }
    return table_p;
}

static double* fmt_make_table_method2( const int m, const int tmax,
        const int nexp, const int ndiv ) {
    double *table_p, *table_i;
    int step, i, k;
    double d, d2, t0;
    size_t table_size;
    d      = 1.e0 / (double)ndiv;
    d2     = 0.5e0 * d;
    step       = nexp + 1;
    table_size = step*tmax*ndiv*sizeof(double);
    table_p = (double*)malloc( table_size );
    for ( i=0, table_i=table_p; i<(tmax*ndiv); i++, table_i+=step ) {
        t0 = d*(double)i + d2;
        for ( k=0; k<nexp; k++ ) table_i[k] = fmt_slow( m+k, t0 );
        table_i[nexp] = exp( -t0 );
    }
    return table_p;
}

static double* fmt_make_table_method3( const int m, const int tmax,
        const int nexp, const int ndiv ) {
    double *table_p, *table_i;
    int step, i;
    double d, d2, t0;
    size_t table_size;
    d      = 1.e0 / (double)ndiv;
    d2     = 0.5e0 * d;
    step       = 1 + 1;
    table_size = step*tmax*ndiv*sizeof(double);
    table_p = (double*)malloc( table_size );
    for ( i=0, table_i=table_p; i<(tmax*ndiv); i++, table_i+=step ) {
        t0 = d*(double)i + d2;
        table_i[0] = fmt_slow( (m+nexp-1), t0 );
        table_i[1] = exp( -t0 );
    }
    return table_p;
}

int cuda_fmt_m_finalize(void) {
  int ret=0;
#pragma omp master
  ret = cuda_FMT_m_Finalize();
  return ret;
}

int cuda_fmt_m_init(void) {
    int k, i, ip;
    int nexp;
    int m, ndiv, tmax;
    double fac[10];
    /*
    table_sizes = (size_t*)malloc(sizeof(size_t) * (8+1) );
    NDIVS = (int*)malloc( sizeof(int) * (8+1) );
    TMAXS = (int*)malloc( sizeof(int) * (8+1) );
    */
    size_t table_sizes[8+1];
    // method[23]

    // NEPS=12
#if CUDA_FMT_M_NEXP == 6
    int NEXPS[8+1]={6,6,6,6, 6,6,6,6, 6};
    int NDIVS[8+1]={16,16,16,16, 16,16,16,16, 16};
#elif CUDA_FMT_M_NEXP == 8
    int NEXPS[8+1]={8,8,8,8, 8,8,8,8, 8};
    int NDIVS[8+1]={4,4,4,4, 4,4,4,4, 4};
#elif CUDA_FMT_M_NEXP == 10
    int NEXPS[8+1]={10,10,10,10, 10,10,10,10, 10};
    int NDIVS[8+1]={2,2,2,2, 2,2,2,2, 2};
#else // CUDA_FMT_M == 0
    int NEXPS[8+1]={10,10,10,10, 10,10,10,10, 10};
    int NDIVS[8+1]={2,2,2,2, 2,2,2,2, 2};
#endif
    int TMAXS[8+1]={26, 30, 33, 36, 39, 41, 43, 45, 48};

#if CUDA_FMT_M == 3
    int METHD[8+1]={3,3,3,3, 3,3,3,3, 3};
#elif CUDA_FMT_M == 2
    int METHD[8+1]={2,2,2,2, 2,2,2,2, 2};
#elif CUDA_FMT_M == 1
    int METHD[8+1]={1,1,1,1, 1,1,1,1, 1};
#else  // CUDA_FMT_M == 0 dummy data
    int METHD[8+1]={0,0,0,0, 0,0,0,0, 0};
#endif
#if defined(CUDA_FMT_M_K1) && CUDA_ARCH >= 350 && CUDA_FMT_M_NEXP == 6
    METHD[1] = 3;
#endif

    double *fmt_m_tbl[8+1];

    int NDEV = cuda_get_numDevice();
    if (NDEV<=0) return 2;

    // table assignment
    fac[0] = fac[1] = 1.e0;
    for ( k=2; k<10; k++ ) fac[k] = fac[k-1] / (double)k;

    for (m=0; m<=8; m++) {
      double *tbl;
      nexp=NEXPS[m];
      ndiv=NDIVS[m];
      tmax=TMAXS[m];
      if (METHD[m]==3) {
        tbl = fmt_make_table_method3( m, tmax, nexp, ndiv );
        table_sizes[m] = (1+1) * tmax * ndiv;
      } else if (METHD[m]==2) {
        tbl = fmt_make_table_method2( m, tmax, nexp, ndiv );
        table_sizes[m] = (nexp+1) * tmax * ndiv;
        for ( i=0; i<(tmax*ndiv); i++ ) {
          ip = i * (1+nexp);
          for ( k=2; k<nexp; k++ ) tbl[ip+k] *= fac[k];
        }
      } else if (METHD[m]==1) {
        tbl = fmt_make_table_method1( m, tmax, nexp, ndiv );
        table_sizes[m] = (m+nexp) * tmax * ndiv;
        if (m==0) {
          for ( i=0; i<(tmax*ndiv); i++ ) {
            ip = i * (0+nexp);
            for ( k=2; k<nexp; k++ ) tbl[ip+k] *= fac[k];
          }
        }
      } else {
        table_sizes[m] = 0;
        tbl = NULL;
      }
      fmt_m_tbl[m]=tbl;
    }

    // CUDA
    int ret;
    int mmax=8;
    ret = cuda_FMT_m_Init(fmt_m_tbl, table_sizes, NDIVS, mmax);
    if (ret<0) return -1;

    for (m=0; m<=8; m++) if ( fmt_m_tbl[m] ) free( fmt_m_tbl[m] );
    return 0;
}
#endif /* CUDA_FMT_M */
