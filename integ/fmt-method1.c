/*
 * fmt-method1 for Intel PC
 * eps=12
 * mmax=8
 * nexp=6 for m=2,3,5,6
 * nexp=8 for m=0,1,4
 * nexp=10 for m=7,8
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fmt-method1.h"

#ifndef false
#define false 0
#endif
#ifndef true 
#define true 1
#endif
#define T_TOL 90.e0

static size_t *table_sizes = NULL;

static int *NDIVS = NULL;

static int *TMAXS = NULL;

static double* fmt0_table = NULL;
static double* fmt1_table = NULL;
static double* fmt2_table = NULL;
static double* fmt3_table = NULL;
static double* fmt4_table = NULL;
static double* fmt5_table = NULL;
static double* fmt6_table = NULL;
static double* fmt7_table = NULL;
static double* fmt8_table = NULL;
static double delta0, dhalf0;
static double delta1, dhalf1;
static double delta2, dhalf2;
static double delta3, dhalf3;
static double delta4, dhalf4;
static double delta5, dhalf5;
static double delta6, dhalf6;
static double delta7, dhalf7;
static double delta8, dhalf8;
static double dinv[10];
static double sqrt_pi_2;
static double sqrt_pi2;
static double PI;

double fmt_method1_get_eps() { return 1.0e-12; }

int fmt_method1_get_mmax() { return 8; }

int fmt_method1_get_nexp() { return 10; }

int fmt_method1_get_ndiv( const int m ) { return NDIVS[m]; }

int fmt_method1_get_tmax( const int m ) { return TMAXS[m]; }

size_t fmt_method1_get_size( const int m ) { return table_sizes[m]; }

static int init() {
    static int called = false;
    if ( called ) return 0;
    PI = 4.e0 * atan( 1.e0 );
    called = true;
    return 0;
}

static double gau_slow( const int m, const double t ) {
    int k;
    double gau, coef, t2;
    init();
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
    init();
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

static void fmt_method1_finalize() {
    if ( fmt0_table ) free( fmt0_table );
    fmt0_table = NULL;
    if ( fmt1_table ) free( fmt1_table );
    fmt1_table = NULL;
    if ( fmt2_table ) free( fmt2_table );
    fmt2_table = NULL;
    if ( fmt3_table ) free( fmt3_table );
    fmt3_table = NULL;
    if ( fmt4_table ) free( fmt4_table );
    fmt4_table = NULL;
    if ( fmt5_table ) free( fmt5_table );
    fmt5_table = NULL;
    if ( fmt6_table ) free( fmt6_table );
    fmt6_table = NULL;
    if ( fmt7_table ) free( fmt7_table );
    fmt7_table = NULL;
    if ( fmt8_table ) free( fmt8_table );
    fmt8_table = NULL;
    if ( table_sizes ) free ( table_sizes );
    table_sizes = NULL;
    if ( NDIVS ) free ( NDIVS );
    NDIVS = NULL;
    if ( TMAXS ) free ( TMAXS );
    TMAXS = NULL;
}

int fmt_method1_init() {
    int k, i, ip;
    int nexp;
    int m, ndiv, tmax;
    double fac[10];
    static int called = false;
    if ( called ) return 0; 
    table_sizes = (size_t*)malloc(sizeof(size_t) * (8+1) );
    NDIVS = (int*)malloc( sizeof(int) * (8+1) );
    TMAXS = (int*)malloc( sizeof(int) * (8+1) );
    //nexps[8+1]={8,8,6,6, 8,6,6,10, 10};
    int m2idx[8+1]={1,1,0,0, 1,0,0,2, 2};
    int idx2nexp[3]={6,8,10};
    int idx2ndiv[3]={16,4,2};
    int tmax1[8+1]={26, 30, 33, 36, 39, 41, 43, 45, 48};
    double *tbl[8+1];
    int mmax=8;

    // table assignment
    for (m=0; m<=mmax; m++) {
      nexp = idx2nexp[m2idx[m]];
      ndiv = idx2ndiv[m2idx[m]];
      tmax = tmax1[m];
      tbl[m] = fmt_make_table_method1( m, tmax, nexp, ndiv );
      table_sizes[m] = sizeof(double) * (m+nexp) * tmax * ndiv;
      NDIVS[m] = ndiv;
      TMAXS[m] = tmax;
      if (m==0) {
        // special procedure for m=0
        fac[0] = fac[1] = 1.e0;
        for ( k=2; k<nexp; k++ ) fac[k] = fac[k-1] / (double)k;
        for ( i=0; i<(tmax*ndiv); i++ ) {
          ip = i * (0+nexp);
          for ( k=2; k<nexp; k++ ) tbl[m][ip+k] *= fac[k];
        }
      }
    }
    fmt0_table = tbl[0];
    fmt1_table = tbl[1];
    fmt2_table = tbl[2];
    fmt3_table = tbl[3];
    fmt4_table = tbl[4];
    fmt5_table = tbl[5];
    fmt6_table = tbl[6];
    fmt7_table = tbl[7];
    fmt8_table = tbl[8];
    delta0 = 1.e0 / (double)NDIVS[0];
    dhalf0 = 0.5e0 * delta0;
    delta1 = 1.e0 / (double)NDIVS[1];
    dhalf1 = 0.5e0 * delta1;
    delta2 = 1.e0 / (double)NDIVS[2];
    dhalf2 = 0.5e0 * delta2;
    delta3 = 1.e0 / (double)NDIVS[3];
    dhalf3 = 0.5e0 * delta3;
    delta4 = 1.e0 / (double)NDIVS[4];
    dhalf4 = 0.5e0 * delta4;
    delta5 = 1.e0 / (double)NDIVS[5];
    dhalf5 = 0.5e0 * delta5;
    delta6 = 1.e0 / (double)NDIVS[6];
    dhalf6 = 0.5e0 * delta6;
    delta7 = 1.e0 / (double)NDIVS[7];
    dhalf7 = 0.5e0 * delta7;
    delta8 = 1.e0 / (double)NDIVS[8];
    dhalf8 = 0.5e0 * delta8;

    // assign constants
    sqrt_pi_2 = sqrt( atan( 1.e0 ) );
    sqrt_pi2  = sqrt( 2.e0*atan( 1.e0 ) );
    dinv[0] = dinv[1] = 1.e0;
    for ( k=2; k<10; k++ ) dinv[k] = 1.e0 / (double)k;
    atexit( fmt_method1_finalize );
    called = true;
    return 0;
}

// automatically generated by gen_fmt_method1 function
// m = 0, nexp = 8, eps = 1.0e-12
void fmt0_method1( const double t, const double coef, double fmt[] ) {
    int it;
    double t0, dt, *f;
    if ( t >= 26 ) {
        fmt[0] = coef * sqrt_pi_2 * sqrt(1.e0/t);
    } else {
        it = (int)(t*4);
        t0 = delta0 * (double)it + dhalf0;
        dt = t0 - t;
        f  = &fmt0_table[it*(0+8)];
        fmt[0] = f[0] + dt *( f[1] + dt *( f[2] + dt *( f[3] + dt *
                ( f[4] + dt *( f[5] + dt *( f[6] + dt * f[7]))))));
        fmt[0] *= coef;
    }
}

// automatically generated by gen_fmt_method1 function
// m = 1, nexp = 8, eps = 1.0e-12
void fmt1_method1( const double t, const double coef, double fmt[] ) {
    int it, k;
    double t0, dt, *f;
    double t2, sqrt2, dtk[8];
    if ( t >= 30 ) {
        sqrt2  = sqrt( 0.5e0 / t );
        t2     = sqrt2*sqrt2;
        fmt[0] = coef * sqrt_pi2 * sqrt2;
        fmt[1] =        t2 * fmt[0];
    } else {
        it = (int)(t*4);
        t0 = delta1 * (double)it + dhalf1;
        dt = t0 - t;
        f  = &fmt1_table[it*(1+8)];
        dtk[1] = dt;
        for ( k=2; k<8; k++ ) dtk[k] = dtk[k-1]*(dinv[k]*dt);
        fmt[0] = f[0] + dtk[1]*f[1] + dtk[2]*f[2] + dtk[3]*f[3]
                 + dtk[4]*f[4] + dtk[5]*f[5] + dtk[6]*f[6] + dtk[7]*f[7];
        fmt[0] *= coef;
        fmt[1] = f[1] + dtk[1]*f[2] + dtk[2]*f[3] + dtk[3]*f[4]
                 + dtk[4]*f[5] + dtk[5]*f[6] + dtk[6]*f[7] + dtk[7]*f[8];
        fmt[1] *= coef;
    }
}

// automatically generated by gen_fmt_method1 function
// m = 2, nexp = 6, eps = 1.0e-12
void fmt2_method1( const double t, const double coef, double fmt[] ) {
    int it, k;
    double t0, dt, *f;
    double t2, sqrt2, dtk[6];
    if ( t >= 33 ) {
        sqrt2  = sqrt( 0.5e0 / t );
        t2     = sqrt2*sqrt2;
        fmt[0] = coef * sqrt_pi2 * sqrt2;
        fmt[1] =        t2 * fmt[0];
        fmt[2] = 3.e0 * t2 * fmt[1];
    } else {
        it = (int)(t*16);
        t0 = delta2 * (double)it + dhalf2;
        dt = t0 - t;
        f  = &fmt2_table[it*(2+6)];
        dtk[1] = dt;
        for ( k=2; k<6; k++ ) dtk[k] = dtk[k-1]*(dinv[k]*dt);
        fmt[0] = f[0] + dtk[1]*f[1] + dtk[2]*f[2] + dtk[3]*f[3]
                 + dtk[4]*f[4] + dtk[5]*f[5];
        fmt[0] *= coef;
        fmt[1] = f[1] + dtk[1]*f[2] + dtk[2]*f[3] + dtk[3]*f[4]
                 + dtk[4]*f[5] + dtk[5]*f[6];
        fmt[1] *= coef;
        fmt[2] = f[2] + dtk[1]*f[3] + dtk[2]*f[4] + dtk[3]*f[5]
                 + dtk[4]*f[6] + dtk[5]*f[7];
        fmt[2] *= coef;
    }
}

// automatically generated by gen_fmt_method1 function
// m = 3, nexp = 6, eps = 1.0e-12
void fmt3_method1( const double t, const double coef, double fmt[] ) {
    int it, k;
    double t0, dt, *f;
    double t2, sqrt2, dtk[6];
    if ( t >= 36 ) {
        sqrt2  = sqrt( 0.5e0 / t );
        t2     = sqrt2*sqrt2;
        fmt[0] = coef * sqrt_pi2 * sqrt2;
        fmt[1] =        t2 * fmt[0];
        fmt[2] = 3.e0 * t2 * fmt[1];
        fmt[3] = 5.e0 * t2 * fmt[2];
    } else {
        it = (int)(t*16);
        t0 = delta3 * (double)it + dhalf3;
        dt = t0 - t;
        f  = &fmt3_table[it*(3+6)];
        dtk[1] = dt;
        for ( k=2; k<6; k++ ) dtk[k] = dtk[k-1]*(dinv[k]*dt);
        fmt[0] = f[0] + dtk[1]*f[1] + dtk[2]*f[2] + dtk[3]*f[3]
                 + dtk[4]*f[4] + dtk[5]*f[5];
        fmt[0] *= coef;
        fmt[1] = f[1] + dtk[1]*f[2] + dtk[2]*f[3] + dtk[3]*f[4]
                 + dtk[4]*f[5] + dtk[5]*f[6];
        fmt[1] *= coef;
        fmt[2] = f[2] + dtk[1]*f[3] + dtk[2]*f[4] + dtk[3]*f[5]
                 + dtk[4]*f[6] + dtk[5]*f[7];
        fmt[2] *= coef;
        fmt[3] = f[3] + dtk[1]*f[4] + dtk[2]*f[5] + dtk[3]*f[6]
                 + dtk[4]*f[7] + dtk[5]*f[8];
        fmt[3] *= coef;
    }
}

// automatically generated by gen_fmt_method1 function
// m = 4, nexp = 8, eps = 1.0e-12
void fmt4_method1( const double t, const double coef, double fmt[] ) {
    int it, k;
    double t0, dt, *f;
    double t2, sqrt2, dtk[8];
    if ( t >= 39 ) {
        sqrt2  = sqrt( 0.5e0 / t );
        t2     = sqrt2*sqrt2;
        fmt[0] = coef * sqrt_pi2 * sqrt2;
        fmt[1] =        t2 * fmt[0];
        fmt[2] = 3.e0 * t2 * fmt[1];
        fmt[3] = 5.e0 * t2 * fmt[2];
        fmt[4] = 7.e0 * t2 * fmt[3];
    } else {
        it = (int)(t*4);
        t0 = delta4 * (double)it + dhalf4;
        dt = t0 - t;
        f  = &fmt4_table[it*(4+8)];
        dtk[1] = dt;
        for ( k=2; k<8; k++ ) dtk[k] = dtk[k-1]*(dinv[k]*dt);
        fmt[0] = f[0] + dtk[1]*f[1] + dtk[2]*f[2] + dtk[3]*f[3]
                 + dtk[4]*f[4] + dtk[5]*f[5] + dtk[6]*f[6] + dtk[7]*f[7];
        fmt[0] *= coef;
        fmt[1] = f[1] + dtk[1]*f[2] + dtk[2]*f[3] + dtk[3]*f[4]
                 + dtk[4]*f[5] + dtk[5]*f[6] + dtk[6]*f[7] + dtk[7]*f[8];
        fmt[1] *= coef;
        fmt[2] = f[2] + dtk[1]*f[3] + dtk[2]*f[4] + dtk[3]*f[5]
                 + dtk[4]*f[6] + dtk[5]*f[7] + dtk[6]*f[8] + dtk[7]*f[9];
        fmt[2] *= coef;
        fmt[3] = f[3] + dtk[1]*f[4] + dtk[2]*f[5] + dtk[3]*f[6]
                 + dtk[4]*f[7] + dtk[5]*f[8] + dtk[6]*f[9] + dtk[7]*f[10];
        fmt[3] *= coef;
        fmt[4] = f[4] + dtk[1]*f[5] + dtk[2]*f[6] + dtk[3]*f[7]
                 + dtk[4]*f[8] + dtk[5]*f[9] + dtk[6]*f[10] + dtk[7]*f[11]
                ;
        fmt[4] *= coef;
    }
}

// automatically generated by gen_fmt_method1 function
// m = 5, nexp = 6, eps = 1.0e-12
void fmt5_method1( const double t, const double coef, double fmt[] ) {
    int it, k;
    double t0, dt, *f;
    double t2, sqrt2, dtk[6];
    if ( t >= 41 ) {
        sqrt2  = sqrt( 0.5e0 / t );
        t2     = sqrt2*sqrt2;
        fmt[0] = coef * sqrt_pi2 * sqrt2;
        fmt[1] =         t2 * fmt[0];
        fmt[2] = 3.e0 * t2 * fmt[1];
        fmt[3] = 5.e0 * t2 * fmt[2];
        fmt[4] = 7.e0 * t2 * fmt[3];
        fmt[5] = 9.e0 * t2 * fmt[4];
    } else {
        it = (int)(t*16);
        t0 = delta5 * (double)it + dhalf5;
        dt = t0 - t;
        f  = &fmt5_table[it*(5+6)];
        dtk[1] = dt;
        for ( k=2; k<6; k++ ) dtk[k] = dtk[k-1]*(dinv[k]*dt);
        fmt[0] = f[0] + dtk[1]*f[1] + dtk[2]*f[2] + dtk[3]*f[3]
                 + dtk[4]*f[4] + dtk[5]*f[5];
        fmt[0] *= coef;
        fmt[1] = f[1] + dtk[1]*f[2] + dtk[2]*f[3] + dtk[3]*f[4]
                 + dtk[4]*f[5] + dtk[5]*f[6];
        fmt[1] *= coef;
        fmt[2] = f[2] + dtk[1]*f[3] + dtk[2]*f[4] + dtk[3]*f[5]
                 + dtk[4]*f[6] + dtk[5]*f[7];
        fmt[2] *= coef;
        fmt[3] = f[3] + dtk[1]*f[4] + dtk[2]*f[5] + dtk[3]*f[6]
                 + dtk[4]*f[7] + dtk[5]*f[8];
        fmt[3] *= coef;
        fmt[4] = f[4] + dtk[1]*f[5] + dtk[2]*f[6] + dtk[3]*f[7]
                 + dtk[4]*f[8] + dtk[5]*f[9];
        fmt[4] *= coef;
        fmt[5] = f[5] + dtk[1]*f[6] + dtk[2]*f[7] + dtk[3]*f[8]
                 + dtk[4]*f[9] + dtk[5]*f[10];
        fmt[5] *= coef;
    }
}

// automatically generated by gen_fmt_method1 function
// m = 6, nexp = 6, eps = 1.0e-12
void fmt6_method1( const double t, const double coef, double fmt[] ) {
    int it, k;
    double t0, dt, *f;
    double t2, sqrt2, dtk[6];
    if ( t >= 43 ) {
        sqrt2  = sqrt( 0.5e0 / t );
        t2     = sqrt2*sqrt2;
        fmt[0] = coef * sqrt_pi2 * sqrt2;
        fmt[1] =         t2 * fmt[0];
        fmt[2] = 3.e0 * t2 * fmt[1];
        fmt[3] = 5.e0 * t2 * fmt[2];
        fmt[4] = 7.e0 * t2 * fmt[3];
        fmt[5] = 9.e0 * t2 * fmt[4];
        fmt[6] = 11.e0 * t2 * fmt[5];
    } else {
        it = (int)(t*16);
        t0 = delta6 * (double)it + dhalf6;
        dt = t0 - t;
        f  = &fmt6_table[it*(6+6)];
        dtk[1] = dt;
        for ( k=2; k<6; k++ ) dtk[k] = dtk[k-1]*(dinv[k]*dt);
        fmt[0] = f[0] + dtk[1]*f[1] + dtk[2]*f[2] + dtk[3]*f[3]
                 + dtk[4]*f[4] + dtk[5]*f[5];
        fmt[0] *= coef;
        fmt[1] = f[1] + dtk[1]*f[2] + dtk[2]*f[3] + dtk[3]*f[4]
                 + dtk[4]*f[5] + dtk[5]*f[6];
        fmt[1] *= coef;
        fmt[2] = f[2] + dtk[1]*f[3] + dtk[2]*f[4] + dtk[3]*f[5]
                 + dtk[4]*f[6] + dtk[5]*f[7];
        fmt[2] *= coef;
        fmt[3] = f[3] + dtk[1]*f[4] + dtk[2]*f[5] + dtk[3]*f[6]
                 + dtk[4]*f[7] + dtk[5]*f[8];
        fmt[3] *= coef;
        fmt[4] = f[4] + dtk[1]*f[5] + dtk[2]*f[6] + dtk[3]*f[7]
                 + dtk[4]*f[8] + dtk[5]*f[9];
        fmt[4] *= coef;
        fmt[5] = f[5] + dtk[1]*f[6] + dtk[2]*f[7] + dtk[3]*f[8]
                 + dtk[4]*f[9] + dtk[5]*f[10];
        fmt[5] *= coef;
        fmt[6] = f[6] + dtk[1]*f[7] + dtk[2]*f[8] + dtk[3]*f[9]
                 + dtk[4]*f[10] + dtk[5]*f[11];
        fmt[6] *= coef;
    }
}

// automatically generated by gen_fmt_method1 function
// m = 7, nexp = 10, eps = 1.0e-12
void fmt7_method1( const double t, const double coef, double fmt[] ) {
    int it, k;
    double t0, dt, *f;
    double t2, sqrt2, dtk[10];
    if ( t >= 45 ) {
        sqrt2  = sqrt( 0.5e0 / t );
        t2     = sqrt2*sqrt2;
        fmt[0] = coef * sqrt_pi2 * sqrt2;
        fmt[1] =         t2 * fmt[0];
        fmt[2] = 3.e0 * t2 * fmt[1];
        fmt[3] = 5.e0 * t2 * fmt[2];
        fmt[4] = 7.e0 * t2 * fmt[3];
        fmt[5] = 9.e0 * t2 * fmt[4];
        fmt[6] = 11.e0 * t2 * fmt[5];
        fmt[7] = 13.e0 * t2 * fmt[6];
    } else {
        it = (int)(t*2);
        t0 = delta7 * (double)it + dhalf7;
        dt = t0 - t;
        f  = &fmt7_table[it*(7+10)];
        dtk[1] = dt;
        for ( k=2; k<10; k++ ) dtk[k] = dtk[k-1]*(dinv[k]*dt);
        fmt[0] = f[0] + dtk[1]*f[1] + dtk[2]*f[2] + dtk[3]*f[3]
                 + dtk[4]*f[4] + dtk[5]*f[5] + dtk[6]*f[6] + dtk[7]*f[7]
                 + dtk[8]*f[8] + dtk[9]*f[9];
        fmt[0] *= coef;
        fmt[1] = f[1] + dtk[1]*f[2] + dtk[2]*f[3] + dtk[3]*f[4]
                 + dtk[4]*f[5] + dtk[5]*f[6] + dtk[6]*f[7] + dtk[7]*f[8]
                 + dtk[8]*f[9] + dtk[9]*f[10];
        fmt[1] *= coef;
        fmt[2] = f[2] + dtk[1]*f[3] + dtk[2]*f[4] + dtk[3]*f[5]
                 + dtk[4]*f[6] + dtk[5]*f[7] + dtk[6]*f[8] + dtk[7]*f[9]
                 + dtk[8]*f[10] + dtk[9]*f[11];
        fmt[2] *= coef;
        fmt[3] = f[3] + dtk[1]*f[4] + dtk[2]*f[5] + dtk[3]*f[6]
                 + dtk[4]*f[7] + dtk[5]*f[8] + dtk[6]*f[9] + dtk[7]*f[10]
                 + dtk[8]*f[11] + dtk[9]*f[12];
        fmt[3] *= coef;
        fmt[4] = f[4] + dtk[1]*f[5] + dtk[2]*f[6] + dtk[3]*f[7]
                 + dtk[4]*f[8] + dtk[5]*f[9] + dtk[6]*f[10] + dtk[7]*f[11]
                 + dtk[8]*f[12] + dtk[9]*f[13];
        fmt[4] *= coef;
        fmt[5] = f[5] + dtk[1]*f[6] + dtk[2]*f[7] + dtk[3]*f[8]
                 + dtk[4]*f[9] + dtk[5]*f[10] + dtk[6]*f[11]
                 + dtk[7]*f[12] + dtk[8]*f[13] + dtk[9]*f[14];
        fmt[5] *= coef;
        fmt[6] = f[6] + dtk[1]*f[7] + dtk[2]*f[8] + dtk[3]*f[9]
                 + dtk[4]*f[10] + dtk[5]*f[11] + dtk[6]*f[12]
                 + dtk[7]*f[13] + dtk[8]*f[14] + dtk[9]*f[15];
        fmt[6] *= coef;
        fmt[7] = f[7] + dtk[1]*f[8] + dtk[2]*f[9] + dtk[3]*f[10]
                 + dtk[4]*f[11] + dtk[5]*f[12] + dtk[6]*f[13]
                 + dtk[7]*f[14] + dtk[8]*f[15] + dtk[9]*f[16];
        fmt[7] *= coef;
    }
}

// automatically generated by gen_fmt_method1 function
// m = 8, nexp = 10, eps = 1.0e-12
void fmt8_method1( const double t, const double coef, double fmt[] ) {
    int it, k;
    double t0, dt, *f;
    double t2, sqrt2, dtk[10];
    if ( t >= 48 ) {
        sqrt2  = sqrt( 0.5e0 / t );
        t2     = sqrt2*sqrt2;
        fmt[0] = coef * sqrt_pi2 * sqrt2;
        fmt[1] =         t2 * fmt[0];
        fmt[2] = 3.e0 * t2 * fmt[1];
        fmt[3] = 5.e0 * t2 * fmt[2];
        fmt[4] = 7.e0 * t2 * fmt[3];
        fmt[5] = 9.e0 * t2 * fmt[4];
        fmt[6] = 11.e0 * t2 * fmt[5];
        fmt[7] = 13.e0 * t2 * fmt[6];
        fmt[8] = 15.e0 * t2 * fmt[7];
    } else {
        it = (int)(t*2);
        t0 = delta8 * (double)it + dhalf8;
        dt = t0 - t;
        f  = &fmt8_table[it*(8+10)];
        dtk[1] = dt;
        for ( k=2; k<10; k++ ) dtk[k] = dtk[k-1]*(dinv[k]*dt);
        fmt[0] = f[0] + dtk[1]*f[1] + dtk[2]*f[2] + dtk[3]*f[3]
                 + dtk[4]*f[4] + dtk[5]*f[5] + dtk[6]*f[6] + dtk[7]*f[7]
                 + dtk[8]*f[8] + dtk[9]*f[9];
        fmt[0] *= coef;
        fmt[1] = f[1] + dtk[1]*f[2] + dtk[2]*f[3] + dtk[3]*f[4]
                 + dtk[4]*f[5] + dtk[5]*f[6] + dtk[6]*f[7] + dtk[7]*f[8]
                 + dtk[8]*f[9] + dtk[9]*f[10];
        fmt[1] *= coef;
        fmt[2] = f[2] + dtk[1]*f[3] + dtk[2]*f[4] + dtk[3]*f[5]
                 + dtk[4]*f[6] + dtk[5]*f[7] + dtk[6]*f[8] + dtk[7]*f[9]
                 + dtk[8]*f[10] + dtk[9]*f[11];
        fmt[2] *= coef;
        fmt[3] = f[3] + dtk[1]*f[4] + dtk[2]*f[5] + dtk[3]*f[6]
                 + dtk[4]*f[7] + dtk[5]*f[8] + dtk[6]*f[9] + dtk[7]*f[10]
                 + dtk[8]*f[11] + dtk[9]*f[12];
        fmt[3] *= coef;
        fmt[4] = f[4] + dtk[1]*f[5] + dtk[2]*f[6] + dtk[3]*f[7]
                 + dtk[4]*f[8] + dtk[5]*f[9] + dtk[6]*f[10] + dtk[7]*f[11]
                 + dtk[8]*f[12] + dtk[9]*f[13];
        fmt[4] *= coef;
        fmt[5] = f[5] + dtk[1]*f[6] + dtk[2]*f[7] + dtk[3]*f[8]
                 + dtk[4]*f[9] + dtk[5]*f[10] + dtk[6]*f[11]
                 + dtk[7]*f[12] + dtk[8]*f[13] + dtk[9]*f[14];
        fmt[5] *= coef;
        fmt[6] = f[6] + dtk[1]*f[7] + dtk[2]*f[8] + dtk[3]*f[9]
                 + dtk[4]*f[10] + dtk[5]*f[11] + dtk[6]*f[12]
                 + dtk[7]*f[13] + dtk[8]*f[14] + dtk[9]*f[15];
        fmt[6] *= coef;
        fmt[7] = f[7] + dtk[1]*f[8] + dtk[2]*f[9] + dtk[3]*f[10]
                 + dtk[4]*f[11] + dtk[5]*f[12] + dtk[6]*f[13]
                 + dtk[7]*f[14] + dtk[8]*f[15] + dtk[9]*f[16];
        fmt[7] *= coef;
        fmt[8] = f[8] + dtk[1]*f[9] + dtk[2]*f[10] + dtk[3]*f[11]
                 + dtk[4]*f[12] + dtk[5]*f[13] + dtk[6]*f[14]
                 + dtk[7]*f[15] + dtk[8]*f[16] + dtk[9]*f[17];
        fmt[8] *= coef;
    }
}
