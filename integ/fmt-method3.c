#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef false
#define false 0
#endif
#ifndef true 
#define true 1
#endif
#define T_TOL 90.e0

static size_t* table_sizes = NULL;

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
static double dinv[8];
static double dinv2[8+6];
static double sqrt_pi2, sqrt_pi_2;
static double PI;

double fmt_method3_get_eps() { return 1.0e-12; }

int fmt_method3_get_mmax() { return 8; }

int fmt_method3_get_nexp() { return 6; }

int fmt_method3_get_ndiv( const int m ) { return NDIVS[m]; }

int fmt_method3_get_tmax( const int m ) { return TMAXS[m]; }

size_t fmt_method3_get_size( const int m ) { return table_sizes[m]; }

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

static void fmt_method3_finalize() {
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
}

int fmt_method3_init() {
    int k;
    int nexp = 6;
    int m, ndiv, tmax;
    double fac[6];
    static int called = false;
    if ( called ) return 0; 
    table_sizes = (size_t*)malloc(sizeof(size_t) * (8+1) );
    NDIVS = (int*)malloc( sizeof(int) * (8+1) );
    TMAXS = (int*)malloc( sizeof(int) * (8+1) );
    // table assignment
    fac[0] = fac[1] = 1.e0;
    for ( k=2; k<6; k++ ) fac[k] = fac[k-1] / (double)k;
    // m=0
    m=0;
    nexp=6;
    ndiv=16;
    tmax=26;
    fmt0_table = fmt_make_table_method3( m, tmax, nexp, ndiv );
    table_sizes[m] = sizeof(double) * (1+1) * tmax * ndiv;
    NDIVS[m] = ndiv;
    TMAXS[m] = tmax;
    delta0 = 1.e0 / (double)ndiv;
    dhalf0 = 0.5e0 * delta0;
    // m=1
    m=1;
    nexp=6;
    ndiv=16;
    tmax=30;
    fmt1_table = fmt_make_table_method3( m, tmax, nexp, ndiv );
    table_sizes[m] = sizeof(double) * (1+1) * tmax * ndiv;
    NDIVS[m] = ndiv;
    TMAXS[m] = tmax;
    delta1 = 1.e0 / (double)ndiv;
    dhalf1 = 0.5e0 * delta1;
    // m=2
    m=2;
    nexp=6;
    ndiv=16;
    tmax=33;
    fmt2_table = fmt_make_table_method3( m, tmax, nexp, ndiv );
    table_sizes[m] = sizeof(double) * (1+1) * tmax * ndiv;
    NDIVS[m] = ndiv;
    TMAXS[m] = tmax;
    delta2 = 1.e0 / (double)ndiv;
    dhalf2 = 0.5e0 * delta2;
    // m=3
    m=3;
    nexp=6;
    ndiv=16;
    tmax=36;
    fmt3_table = fmt_make_table_method3( m, tmax, nexp, ndiv );
    table_sizes[m] = sizeof(double) * (1+1) * tmax * ndiv;
    NDIVS[m] = ndiv;
    TMAXS[m] = tmax;
    delta3 = 1.e0 / (double)ndiv;
    dhalf3 = 0.5e0 * delta3;
    // m=4
    m=4;
    nexp=6;
    ndiv=16;
    tmax=39;
    fmt4_table = fmt_make_table_method3( m, tmax, nexp, ndiv );
    table_sizes[m] = sizeof(double) * (1+1) * tmax * ndiv;
    NDIVS[m] = ndiv;
    TMAXS[m] = tmax;
    delta4 = 1.e0 / (double)ndiv;
    dhalf4 = 0.5e0 * delta4;
    // m=5
    m=5;
    nexp=6;
    ndiv=16;
    tmax=41;
    fmt5_table = fmt_make_table_method3( m, tmax, nexp, ndiv );
    table_sizes[m] = sizeof(double) * (1+1) * tmax * ndiv;
    NDIVS[m] = ndiv;
    TMAXS[m] = tmax;
    delta5 = 1.e0 / (double)ndiv;
    dhalf5 = 0.5e0 * delta5;
    // m=6
    m=6;
    nexp=6;
    ndiv=16;
    tmax=43;
    fmt6_table = fmt_make_table_method3( m, tmax, nexp, ndiv );
    table_sizes[m] = sizeof(double) * (1+1) * tmax * ndiv;
    NDIVS[m] = ndiv;
    TMAXS[m] = tmax;
    delta6 = 1.e0 / (double)ndiv;
    dhalf6 = 0.5e0 * delta6;
    // m=7
    m=7;
    nexp=6;
    ndiv=16;
    tmax=45;
    fmt7_table = fmt_make_table_method3( m, tmax, nexp, ndiv );
    table_sizes[m] = sizeof(double) * (1+1) * tmax * ndiv;
    NDIVS[m] = ndiv;
    TMAXS[m] = tmax;
    delta7 = 1.e0 / (double)ndiv;
    dhalf7 = 0.5e0 * delta7;
    // m=8
    m=8;
    nexp=6;
    ndiv=16;
    tmax=48;
    fmt8_table = fmt_make_table_method3( m, tmax, nexp, ndiv );
    table_sizes[m] = sizeof(double) * (1+1) * tmax * ndiv;
    NDIVS[m] = ndiv;
    TMAXS[m] = tmax;
    delta8 = 1.e0 / (double)ndiv;
    dhalf8 = 0.5e0 * delta8;
    // assign constants
    sqrt_pi_2 = sqrt( atan( 1.e0 ) );
    sqrt_pi2  = sqrt( 2.e0*atan( 1.e0 ) );
    dinv[0] = dinv[1] = 1.e0;
    for ( k=2; k<8; k++ ) dinv[k] = 1.e0 / (double)k;
    dinv2[0] = 1.e0;
    for ( k=1; k<(8+6); k++ ) dinv2[k] = 1.e0/ (double)(2*k+1);
    atexit( fmt_method3_finalize );
    called = true;
    return 0;
}


// automatically generated by gen_fmt_method3 function
// m = 0, nexp = 6, eps = 1.0e-12
void fmt0_method3( const double t, const double coef, double fmt[] ) {
    int it;
    double t0, dt, *f, expt0, t2, ff[6];
    if ( t >= 26 ) {
        fmt[0] = coef * sqrt_pi_2 * sqrt(1.e0/t);
    } else {
        it = (int)(t*16);
        t0 = delta0 * (double)it + dhalf0;
        t2 = t0 + t0;
        dt = t0 - t;
        f  = &fmt0_table[it*2];
        ff[6-1] = f[0];
        expt0       = f[1];
        // F(m+k)(T0) (k=0, n-2)
        ff[4] = dinv2[0+4] * ( expt0 + t2*ff[5] );
        ff[3] = dinv2[0+3] * ( expt0 + t2*ff[4] );
        ff[2] = dinv2[0+2] * ( expt0 + t2*ff[3] );
        ff[1] = dinv2[0+1] * ( expt0 + t2*ff[2] );
        ff[0] = dinv2[0+0] * ( expt0 + t2*ff[1] );
        // F[0]
        fmt[0] = ff[0] + dt*( ff[1] + dinv[2]*dt*( ff[2] + dinv[3]*dt*
                ( ff[3] + dinv[4]*dt*( ff[4] + dinv[5]*dt*ff[5]))));
        fmt[0] *= coef;
    }
}

// automatically generated by gen_fmt_method3 function
// m = 1, nexp = 6, eps = 1.0e-12
void fmt1_method3( const double t, const double coef, double fmt[] ) {
    int it, k;
    double t0, dt, *f, expt0, t2, ff[6], sqrt2, dtk[8], expt;
    if ( t >= 30 ) {
        sqrt2  = sqrt( 0.5e0 / t );
        t2     = sqrt2*sqrt2;
        fmt[0] = coef * sqrt_pi2 * sqrt2;
        fmt[1] =        t2 * fmt[0];
    } else {
        it = (int)(t*16);
        t0 = delta1 * (double)it + dhalf1;
        t2 = t0 + t0;
        dt = t0 - t;
        f  = &fmt1_table[it*2];
        ff[6-1] = f[0];
        expt0       = f[1];
        // (-dt)^k/k!
        dtk[1] = dt;
        for ( k=2; k<8; k++ ) dtk[k] = dt*dinv[k]*dtk[k-1];
        // exp(-T)
        expt = 1.e0 + dtk[1] + dtk[2] + dtk[3] + dtk[4] + dtk[5] + dtk[6]
                 + dtk[7];
        expt *= expt0;
        // F(m+k)(T0) (k=0, n-2)
        ff[4] = dinv2[1+4] * ( expt0 + t2*ff[5] );
        ff[3] = dinv2[1+3] * ( expt0 + t2*ff[4] );
        ff[2] = dinv2[1+2] * ( expt0 + t2*ff[3] );
        ff[1] = dinv2[1+1] * ( expt0 + t2*ff[2] );
        ff[0] = dinv2[1+0] * ( expt0 + t2*ff[1] );
        // F[1]
        fmt[1] = ff[0] + dtk[1]*ff[1] + dtk[2]*ff[2] + dtk[3]*ff[3]
                 + dtk[4]*ff[4] + dtk[5]*ff[5];
        fmt[1] *= coef;
        // F[0]-F[0]
        t2 = t + t;
        fmt[0] = coef*expt + t2*fmt[1];
    }
}

// automatically generated by gen_fmt_method3 function
// m = 2, nexp = 6, eps = 1.0e-12
void fmt2_method3( const double t, const double coef, double fmt[] ) {
    int it, k;
    double t0, dt, *f, expt0, t2, ff[6], sqrt2, dtk[8], expt;
    if ( t >= 33 ) {
        sqrt2  = sqrt( 0.5e0 / t );
        t2     = sqrt2*sqrt2;
        fmt[0] = coef * sqrt_pi2 * sqrt2;
        fmt[1] =        t2 * fmt[0];
        fmt[2] = 3.e0 * t2 * fmt[1];
    } else {
        it = (int)(t*16);
        t0 = delta2 * (double)it + dhalf2;
        t2 = t0 + t0;
        dt = t0 - t;
        f  = &fmt2_table[it*2];
        ff[6-1] = f[0];
        expt0       = f[1];
        // (-dt)^k/k!
        dtk[1] = dt;
        for ( k=2; k<8; k++ ) dtk[k] = dt*dinv[k]*dtk[k-1];
        // exp(-T)
        expt = 1.e0 + dtk[1] + dtk[2] + dtk[3] + dtk[4] + dtk[5] + dtk[6]
                 + dtk[7];
        expt *= expt0;
        // F(m+k)(T0) (k=0, n-2)
        ff[4] = dinv2[2+4] * ( expt0 + t2*ff[5] );
        ff[3] = dinv2[2+3] * ( expt0 + t2*ff[4] );
        ff[2] = dinv2[2+2] * ( expt0 + t2*ff[3] );
        ff[1] = dinv2[2+1] * ( expt0 + t2*ff[2] );
        ff[0] = dinv2[2+0] * ( expt0 + t2*ff[1] );
        // F[2]
        fmt[2] = ff[0] + dtk[1]*ff[1] + dtk[2]*ff[2] + dtk[3]*ff[3]
                 + dtk[4]*ff[4] + dtk[5]*ff[5];
        fmt[2] *= coef;
        // F[1]-F[0]
        t2 = t + t;
        expt      *= coef;
        fmt[1] = dinv2[1] * ( expt + t2 * fmt[2] );
        fmt[0] =              expt + t2 * fmt[1];
    }
}

// automatically generated by gen_fmt_method3 function
// m = 3, nexp = 6, eps = 1.0e-12
void fmt3_method3( const double t, const double coef, double fmt[] ) {
    int it, k;
    double t0, dt, *f, expt0, t2, ff[6], sqrt2, dtk[8], expt;
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
        t2 = t0 + t0;
        dt = t0 - t;
        f  = &fmt3_table[it*2];
        ff[6-1] = f[0];
        expt0       = f[1];
        // (-dt)^k/k!
        dtk[1] = dt;
        for ( k=2; k<8; k++ ) dtk[k] = dt*dinv[k]*dtk[k-1];
        // exp(-T)
        expt = 1.e0 + dtk[1] + dtk[2] + dtk[3] + dtk[4] + dtk[5] + dtk[6]
                 + dtk[7];
        expt *= expt0;
        // F(m+k)(T0) (k=0, n-2)
        ff[4] = dinv2[3+4] * ( expt0 + t2*ff[5] );
        ff[3] = dinv2[3+3] * ( expt0 + t2*ff[4] );
        ff[2] = dinv2[3+2] * ( expt0 + t2*ff[3] );
        ff[1] = dinv2[3+1] * ( expt0 + t2*ff[2] );
        ff[0] = dinv2[3+0] * ( expt0 + t2*ff[1] );
        // F[3]
        fmt[3] = ff[0] + dtk[1]*ff[1] + dtk[2]*ff[2] + dtk[3]*ff[3]
                 + dtk[4]*ff[4] + dtk[5]*ff[5];
        fmt[3] *= coef;
        // F[2]-F[0]
        t2 = t + t;
        expt      *= coef;
        fmt[2] = dinv2[2] * ( expt + t2 * fmt[3] );
        fmt[1] = dinv2[1] * ( expt + t2 * fmt[2] );
        fmt[0] =              expt + t2 * fmt[1];
    }
}

// automatically generated by gen_fmt_method3 function
// m = 4, nexp = 6, eps = 1.0e-12
void fmt4_method3( const double t, const double coef, double fmt[] ) {
    int it, k;
    double t0, dt, *f, expt0, t2, ff[6], sqrt2, dtk[8], expt;
    if ( t >= 39 ) {
        sqrt2  = sqrt( 0.5e0 / t );
        t2     = sqrt2*sqrt2;
        fmt[0] = coef * sqrt_pi2 * sqrt2;
        fmt[1] =        t2 * fmt[0];
        fmt[2] = 3.e0 * t2 * fmt[1];
        fmt[3] = 5.e0 * t2 * fmt[2];
        fmt[4] = 7.e0 * t2 * fmt[3];
    } else {
        it = (int)(t*16);
        t0 = delta4 * (double)it + dhalf4;
        t2 = t0 + t0;
        dt = t0 - t;
        f  = &fmt4_table[it*2];
        ff[6-1] = f[0];
        expt0       = f[1];
        // (-dt)^k/k!
        dtk[1] = dt;
        for ( k=2; k<8; k++ ) dtk[k] = dt*dinv[k]*dtk[k-1];
        // exp(-T)
        expt = 1.e0 + dtk[1] + dtk[2] + dtk[3] + dtk[4] + dtk[5] + dtk[6]
                 + dtk[7];
        expt *= expt0;
        // F(m+k)(T0) (k=0, n-2)
        ff[4] = dinv2[4+4] * ( expt0 + t2*ff[5] );
        ff[3] = dinv2[4+3] * ( expt0 + t2*ff[4] );
        ff[2] = dinv2[4+2] * ( expt0 + t2*ff[3] );
        ff[1] = dinv2[4+1] * ( expt0 + t2*ff[2] );
        ff[0] = dinv2[4+0] * ( expt0 + t2*ff[1] );
        // F[4]
        fmt[4] = ff[0] + dtk[1]*ff[1] + dtk[2]*ff[2] + dtk[3]*ff[3]
                 + dtk[4]*ff[4] + dtk[5]*ff[5];
        fmt[4] *= coef;
        // F[3]-F[0]
        t2 = t + t;
        expt      *= coef;
        fmt[3] = dinv2[3] * ( expt + t2 * fmt[4] );
        fmt[2] = dinv2[2] * ( expt + t2 * fmt[3] );
        fmt[1] = dinv2[1] * ( expt + t2 * fmt[2] );
        fmt[0] =              expt + t2 * fmt[1];
    }
}

// automatically generated by gen_fmt_method3 function
// m = 5, nexp = 6, eps = 1.0e-12
void fmt5_method3( const double t, const double coef, double fmt[] ) {
    int it, k;
    double t0, dt, *f, expt0, t2, ff[6], sqrt2, dtk[8], expt;
    if ( t >= 41 ) {
        sqrt2  = sqrt( 0.5e0 / t );
        t2     = sqrt2*sqrt2;
        fmt[0] = coef * sqrt_pi2 * sqrt2;
        fmt[1] =         t2 * fmt[0];
        fmt[2] =  3.e0 * t2 * fmt[1];
        fmt[3] =  5.e0 * t2 * fmt[2];
        fmt[4] =  7.e0 * t2 * fmt[3];
        fmt[5] =  9.e0 * t2 * fmt[4];
    } else {
        it = (int)(t*16);
        t0 = delta5 * (double)it + dhalf5;
        t2 = t0 + t0;
        dt = t0 - t;
        f  = &fmt5_table[it*2];
        ff[6-1] = f[0];
        expt0       = f[1];
        // (-dt)^k/k!
        dtk[1] = dt;
        for ( k=2; k<8; k++ ) dtk[k] = dt*dinv[k]*dtk[k-1];
        // exp(-T)
        expt = 1.e0 + dtk[1] + dtk[2] + dtk[3] + dtk[4] + dtk[5] + dtk[6]
                 + dtk[7];
        expt *= expt0;
        // F(m+k)(T0) (k=0, n-2)
        ff[4] = dinv2[5+4] * ( expt0 + t2*ff[5] );
        ff[3] = dinv2[5+3] * ( expt0 + t2*ff[4] );
        ff[2] = dinv2[5+2] * ( expt0 + t2*ff[3] );
        ff[1] = dinv2[5+1] * ( expt0 + t2*ff[2] );
        ff[0] = dinv2[5+0] * ( expt0 + t2*ff[1] );
        // F[5]
        fmt[5] = ff[0] + dtk[1]*ff[1] + dtk[2]*ff[2] + dtk[3]*ff[3]
                 + dtk[4]*ff[4] + dtk[5]*ff[5];
        fmt[5] *= coef;
        // F[4]-F[0]
        t2 = t + t;
        expt      *= coef;
        fmt[4] = dinv2[4] * ( expt + t2 * fmt[5] );
        fmt[3] = dinv2[3] * ( expt + t2 * fmt[4] );
        fmt[2] = dinv2[2] * ( expt + t2 * fmt[3] );
        fmt[1] = dinv2[1] * ( expt + t2 * fmt[2] );
        fmt[0] =              expt + t2 * fmt[1];
    }
}

// automatically generated by gen_fmt_method3 function
// m = 6, nexp = 6, eps = 1.0e-12
void fmt6_method3( const double t, const double coef, double fmt[] ) {
    int it, k;
    double t0, dt, *f, expt0, t2, ff[6], sqrt2, dtk[8], expt;
    if ( t >= 43 ) {
        sqrt2  = sqrt( 0.5e0 / t );
        t2     = sqrt2*sqrt2;
        fmt[0] = coef * sqrt_pi2 * sqrt2;
        fmt[1] =         t2 * fmt[0];
        fmt[2] =  3.e0 * t2 * fmt[1];
        fmt[3] =  5.e0 * t2 * fmt[2];
        fmt[4] =  7.e0 * t2 * fmt[3];
        fmt[5] =  9.e0 * t2 * fmt[4];
        fmt[6] = 11.e0 * t2 * fmt[5];
    } else {
        it = (int)(t*16);
        t0 = delta6 * (double)it + dhalf6;
        t2 = t0 + t0;
        dt = t0 - t;
        f  = &fmt6_table[it*2];
        ff[6-1] = f[0];
        expt0       = f[1];
        // (-dt)^k/k!
        dtk[1] = dt;
        for ( k=2; k<8; k++ ) dtk[k] = dt*dinv[k]*dtk[k-1];
        // exp(-T)
        expt = 1.e0 + dtk[1] + dtk[2] + dtk[3] + dtk[4] + dtk[5] + dtk[6]
                 + dtk[7];
        expt *= expt0;
        // F(m+k)(T0) (k=0, n-2)
        ff[4] = dinv2[6+4] * ( expt0 + t2*ff[5] );
        ff[3] = dinv2[6+3] * ( expt0 + t2*ff[4] );
        ff[2] = dinv2[6+2] * ( expt0 + t2*ff[3] );
        ff[1] = dinv2[6+1] * ( expt0 + t2*ff[2] );
        ff[0] = dinv2[6+0] * ( expt0 + t2*ff[1] );
        // F[6]
        fmt[6] = ff[0] + dtk[1]*ff[1] + dtk[2]*ff[2] + dtk[3]*ff[3]
                 + dtk[4]*ff[4] + dtk[5]*ff[5];
        fmt[6] *= coef;
        // F[5]-F[0]
        t2 = t + t;
        expt      *= coef;
        fmt[5] = dinv2[5] * ( expt + t2 * fmt[6] );
        fmt[4] = dinv2[4] * ( expt + t2 * fmt[5] );
        fmt[3] = dinv2[3] * ( expt + t2 * fmt[4] );
        fmt[2] = dinv2[2] * ( expt + t2 * fmt[3] );
        fmt[1] = dinv2[1] * ( expt + t2 * fmt[2] );
        fmt[0] =              expt + t2 * fmt[1];
    }
}

// automatically generated by gen_fmt_method3 function
// m = 7, nexp = 6, eps = 1.0e-12
void fmt7_method3( const double t, const double coef, double fmt[] ) {
    int it, k;
    double t0, dt, *f, expt0, t2, ff[6], sqrt2, dtk[8], expt;
    if ( t >= 45 ) {
        sqrt2  = sqrt( 0.5e0 / t );
        t2     = sqrt2*sqrt2;
        fmt[0] = coef * sqrt_pi2 * sqrt2;
        fmt[1] =         t2 * fmt[0];
        fmt[2] =  3.e0 * t2 * fmt[1];
        fmt[3] =  5.e0 * t2 * fmt[2];
        fmt[4] =  7.e0 * t2 * fmt[3];
        fmt[5] =  9.e0 * t2 * fmt[4];
        fmt[6] = 11.e0 * t2 * fmt[5];
        fmt[7] = 13.e0 * t2 * fmt[6];
    } else {
        it = (int)(t*16);
        t0 = delta7 * (double)it + dhalf7;
        t2 = t0 + t0;
        dt = t0 - t;
        f  = &fmt7_table[it*2];
        ff[6-1] = f[0];
        expt0       = f[1];
        // (-dt)^k/k!
        dtk[1] = dt;
        for ( k=2; k<8; k++ ) dtk[k] = dt*dinv[k]*dtk[k-1];
        // exp(-T)
        expt = 1.e0 + dtk[1] + dtk[2] + dtk[3] + dtk[4] + dtk[5] + dtk[6]
                 + dtk[7];
        expt *= expt0;
        // F(m+k)(T0) (k=0, n-2)
        ff[4] = dinv2[7+4] * ( expt0 + t2*ff[5] );
        ff[3] = dinv2[7+3] * ( expt0 + t2*ff[4] );
        ff[2] = dinv2[7+2] * ( expt0 + t2*ff[3] );
        ff[1] = dinv2[7+1] * ( expt0 + t2*ff[2] );
        ff[0] = dinv2[7+0] * ( expt0 + t2*ff[1] );
        // F[7]
        fmt[7] = ff[0] + dtk[1]*ff[1] + dtk[2]*ff[2] + dtk[3]*ff[3]
                 + dtk[4]*ff[4] + dtk[5]*ff[5];
        fmt[7] *= coef;
        // F[6]-F[0]
        t2 = t + t;
        expt      *= coef;
        fmt[6] = dinv2[6] * ( expt + t2 * fmt[7] );
        fmt[5] = dinv2[5] * ( expt + t2 * fmt[6] );
        fmt[4] = dinv2[4] * ( expt + t2 * fmt[5] );
        fmt[3] = dinv2[3] * ( expt + t2 * fmt[4] );
        fmt[2] = dinv2[2] * ( expt + t2 * fmt[3] );
        fmt[1] = dinv2[1] * ( expt + t2 * fmt[2] );
        fmt[0] =              expt + t2 * fmt[1];
    }
}

// automatically generated by gen_fmt_method3 function
// m = 8, nexp = 6, eps = 1.0e-12
void fmt8_method3( const double t, const double coef, double fmt[] ) {
    int it, k;
    double t0, dt, *f, expt0, t2, ff[6], sqrt2, dtk[8], expt;
    if ( t >= 48 ) {
        sqrt2  = sqrt( 0.5e0 / t );
        t2     = sqrt2*sqrt2;
        fmt[0] = coef * sqrt_pi2 * sqrt2;
        fmt[1] =         t2 * fmt[0];
        fmt[2] =  3.e0 * t2 * fmt[1];
        fmt[3] =  5.e0 * t2 * fmt[2];
        fmt[4] =  7.e0 * t2 * fmt[3];
        fmt[5] =  9.e0 * t2 * fmt[4];
        fmt[6] = 11.e0 * t2 * fmt[5];
        fmt[7] = 13.e0 * t2 * fmt[6];
        fmt[8] = 15.e0 * t2 * fmt[7];
    } else {
        it = (int)(t*16);
        t0 = delta8 * (double)it + dhalf8;
        t2 = t0 + t0;
        dt = t0 - t;
        f  = &fmt8_table[it*2];
        ff[6-1] = f[0];
        expt0       = f[1];
        // (-dt)^k/k!
        dtk[1] = dt;
        for ( k=2; k<8; k++ ) dtk[k] = dt*dinv[k]*dtk[k-1];
        // exp(-T)
        expt = 1.e0 + dtk[1] + dtk[2] + dtk[3] + dtk[4] + dtk[5] + dtk[6]
                 + dtk[7];
        expt *= expt0;
        // F(m+k)(T0) (k=0, n-2)
        ff[4] = dinv2[8+4] * ( expt0 + t2*ff[5] );
        ff[3] = dinv2[8+3] * ( expt0 + t2*ff[4] );
        ff[2] = dinv2[8+2] * ( expt0 + t2*ff[3] );
        ff[1] = dinv2[8+1] * ( expt0 + t2*ff[2] );
        ff[0] = dinv2[8+0] * ( expt0 + t2*ff[1] );
        // F[8]
        fmt[8] = ff[0] + dtk[1]*ff[1] + dtk[2]*ff[2] + dtk[3]*ff[3]
                 + dtk[4]*ff[4] + dtk[5]*ff[5];
        fmt[8] *= coef;
        // F[7]-F[0]
        t2 = t + t;
        expt      *= coef;
        fmt[7] = dinv2[7] * ( expt + t2 * fmt[8] );
        fmt[6] = dinv2[6] * ( expt + t2 * fmt[7] );
        fmt[5] = dinv2[5] * ( expt + t2 * fmt[6] );
        fmt[4] = dinv2[4] * ( expt + t2 * fmt[5] );
        fmt[3] = dinv2[3] * ( expt + t2 * fmt[4] );
        fmt[2] = dinv2[2] * ( expt + t2 * fmt[3] );
        fmt[1] = dinv2[1] * ( expt + t2 * fmt[2] );
        fmt[0] =              expt + t2 * fmt[1];
    }
}
