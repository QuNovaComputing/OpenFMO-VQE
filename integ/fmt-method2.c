/*
 * fmt-method2 for Intel PC
 * eps=12
 * mmax=8
 * nexp=6 for m=0-4
 * nexp=8 for m=5
 * nexp=10 for m=6-8
*/

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
static double dinv[11];
static double dinv2[8];
static double sqrt_pi2, sqrt_pi_2;
static double PI;

double fmt_method2_get_eps() { return 1.0e-12; }

int fmt_method2_get_mmax() { return 8; }

int fmt_method2_get_nexp() { return 10; }

int fmt_method2_get_ndiv( const int m ) { return NDIVS[m]; }

int fmt_method2_get_tmax( const int m ) { return TMAXS[m]; }

size_t fmt_method2_get_size( const int m ) { return table_sizes[m]; }

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

static void fmt_method2_finalize() {
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

int fmt_method2_init() {
    int k, i, ip;
    int nexp;
    int m, ndiv, tmax;
    double fac[10];
    static int called = false;
    if ( called ) return 0; 
    table_sizes = (size_t*)malloc(sizeof(size_t) * (8+1) );
    NDIVS = (int*)malloc( sizeof(int) * (8+1) );
    TMAXS = (int*)malloc( sizeof(int) * (8+1) );
    //nexps[8+1]={6,6,6,6, 6,8,10,10, 10};
    int m2idx[8+1]={0,0,0,0, 0,1,2,2, 2};
    int idx2nexp[3]={6,8,10};
    int idx2ndiv[3]={16,4,2};
    int tmax1[8+1]={26, 30, 33, 36, 39, 41, 43, 45, 48};
    double *tbl[8+1];
    int mmax=8;

    // table assignment
    fac[0] = fac[1] = 1.e0;
    for ( k=2; k<10; k++ ) fac[k] = fac[k-1] / (double)k;

    for (m=0; m<=mmax; m++) {
      nexp=idx2nexp[m2idx[m]];
      ndiv=idx2ndiv[m2idx[m]];
      tmax=tmax1[m];
      tbl[m] = fmt_make_table_method2( m, tmax, nexp, ndiv );
      table_sizes[m] = sizeof(double) * (nexp+1) * tmax * ndiv;
      NDIVS[m] = ndiv;
      TMAXS[m] = tmax;
      for ( i=0; i<(tmax*ndiv); i++ ) {
        ip = i * (1+nexp);
        for ( k=2; k<nexp; k++ ) tbl[m][ip+k] *= fac[k];
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
    for ( k=2; k<11; k++ ) dinv[k] = 1.e0 / (double)k;
    dinv2[0] = 1.e0;
    for ( k=1; k<8; k++ ) dinv2[k] = 1.e0/ (double)(2*k+1);
    atexit( fmt_method2_finalize );
    called = true;
    return 0;
}

// automatically generated by gen_fmt_method2 function
// m = 0, nexp = 6, eps = 1.0e-12
void fmt0_method2( const double t, const double coef, double fmt[] ) {
    int it;
    double t0, dt, *f;
    if ( t >= 26 ) {
        fmt[0] = coef * sqrt_pi_2 * sqrt(1.e0/t);
    } else {
        it = (int)(t*16);
        t0 = delta0 * (double)it + dhalf0;
        dt = t0 - t;
        f  = &fmt0_table[it*(6+1)];
        // F0(T)
        fmt[0] = f[0] + dt *( f[1] + dt *( f[2] + dt *( f[3] + dt *
                ( f[4] + dt * f[5]))));
        fmt[0] *= coef;
    }
}

// automatically generated by gen_fmt_method2 function
// m = 1, nexp = 6, eps = 1.0e-12
void fmt1_method2( const double t, const double coef, double fmt[] ) {
    int it;
    double t0, dt, *f, t2, sqrt2, expt;
    if ( t >= 30 ) {
        sqrt2  = sqrt( 0.5e0 / t );
        t2     = sqrt2*sqrt2;
        fmt[0] = coef * sqrt_pi2 * sqrt2;
        fmt[1] =        t2 * fmt[0];
    } else {
        it = (int)(t*16);
        t0 = delta1 * (double)it + dhalf1;
        t2 = t + t;
        dt = t0 - t;
        f  = &fmt1_table[it*(6+1)];
        // F1(T)
        fmt[1] = f[0] + dt *( f[1] + dt *( f[2] + dt *( f[3] + dt *
                ( f[4] + dt * f[5]))));
        // exp(-T)
        expt = 1.e0 + dt*(1.e0 + dinv[2]*dt*(1.e0 + dinv[3]*dt
                *(1.e0 + dinv[4]*dt*(1.e0 + dinv[5]*dt*(1.e0 + dinv[6]*dt
                *(1.e0 + dinv[7]*dt))))));
        expt *= f[7-1];
        // F[0]-F[0]
        fmt[1] *= coef;
        fmt[0] = coef*expt + t2*fmt[1];
    }
}

// automatically generated by gen_fmt_method2 function
// m = 2, nexp = 6, eps = 1.0e-12
void fmt2_method2( const double t, const double coef, double fmt[] ) {
    int it;
    double t0, dt, *f, t2, sqrt2, expt;
    if ( t >= 33 ) {
        sqrt2  = sqrt( 0.5e0 / t );
        t2     = sqrt2*sqrt2;
        fmt[0] = coef * sqrt_pi2 * sqrt2;
        fmt[1] =        t2 * fmt[0];
        fmt[2] = 3.e0 * t2 * fmt[1];
    } else {
        it = (int)(t*16);
        t0 = delta2 * (double)it + dhalf2;
        t2 = t + t;
        dt = t0 - t;
        f  = &fmt2_table[it*(6+1)];
        // F2(T)
        fmt[2] = f[0] + dt *( f[1] + dt *( f[2] + dt *( f[3] + dt *
                ( f[4] + dt * f[5]))));
        // exp(-T)
        expt = 1.e0 + dt*(1.e0 + dinv[2]*dt*(1.e0 + dinv[3]*dt
                *(1.e0 + dinv[4]*dt*(1.e0 + dinv[5]*dt*(1.e0 + dinv[6]*dt
                *(1.e0 + dinv[7]*dt))))));
        expt *= f[7-1];
        // F[1]-F[0]
        fmt[2] *= coef;
        expt   *= coef;
        fmt[1] = dinv2[1]*( expt + t2*fmt[2] );
        fmt[0] =            expt + t2*fmt[1];
    }
}

// automatically generated by gen_fmt_method2 function
// m = 3, nexp = 6, eps = 1.0e-12
void fmt3_method2( const double t, const double coef, double fmt[] ) {
    int it;
    double t0, dt, *f, t2, sqrt2, expt;
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
        t2 = t + t;
        dt = t0 - t;
        f  = &fmt3_table[it*(6+1)];
        // F3(T)
        fmt[3] = f[0] + dt *( f[1] + dt *( f[2] + dt *( f[3] + dt *
                ( f[4] + dt * f[5]))));
        // exp(-T)
        expt = 1.e0 + dt*(1.e0 + dinv[2]*dt*(1.e0 + dinv[3]*dt
                *(1.e0 + dinv[4]*dt*(1.e0 + dinv[5]*dt*(1.e0 + dinv[6]*dt
                *(1.e0 + dinv[7]*dt))))));
        expt *= f[7-1];
        // F[2]-F[0]
        fmt[3] *= coef;
        expt   *= coef;
        fmt[2] = dinv2[2]*( expt + t2*fmt[3] );
        fmt[1] = dinv2[1]*( expt + t2*fmt[2] );
        fmt[0] =            expt + t2*fmt[1];
    }
}

// automatically generated by gen_fmt_method2 function
// m = 4, nexp = 6, eps = 1.0e-12
void fmt4_method2( const double t, const double coef, double fmt[] ) {
    int it;
    double t0, dt, *f, t2, sqrt2, expt;
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
        t2 = t + t;
        dt = t0 - t;
        f  = &fmt4_table[it*(6+1)];
        // F4(T)
        fmt[4] = f[0] + dt *( f[1] + dt *( f[2] + dt *( f[3] + dt *
                ( f[4] + dt * f[5]))));
        // exp(-T)
        expt = 1.e0 + dt*(1.e0 + dinv[2]*dt*(1.e0 + dinv[3]*dt
                *(1.e0 + dinv[4]*dt*(1.e0 + dinv[5]*dt*(1.e0 + dinv[6]*dt
                *(1.e0 + dinv[7]*dt))))));
        expt *= f[7-1];
        // F[3]-F[0]
        fmt[4] *= coef;
        expt   *= coef;
        fmt[3] = dinv2[3]*( expt + t2*fmt[4] );
        fmt[2] = dinv2[2]*( expt + t2*fmt[3] );
        fmt[1] = dinv2[1]*( expt + t2*fmt[2] );
        fmt[0] =            expt + t2*fmt[1];
    }
}

// automatically generated by gen_fmt_method2 function
// m = 5, nexp = 8, eps = 1.0e-12
void fmt5_method2( const double t, const double coef, double fmt[] ) {
    int it;
    double t0, dt, *f, t2, sqrt2, expt;
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
        it = (int)(t*4);
        t0 = delta5 * (double)it + dhalf5;
        t2 = t + t;
        dt = t0 - t;
        f  = &fmt5_table[it*(8+1)];
        // F5(T)
        fmt[5] = f[0] + dt *( f[1] + dt *( f[2] + dt *( f[3] + dt *
                ( f[4] + dt *( f[5] + dt *( f[6] + dt * f[7]))))));
        // exp(-T)
        expt = 1.e0 + dt*(1.e0 + dinv[2]*dt*(1.e0 + dinv[3]*dt
                *(1.e0 + dinv[4]*dt*(1.e0 + dinv[5]*dt*(1.e0 + dinv[6]*dt
                *(1.e0 + dinv[7]*dt*(1.e0 + dinv[8]*dt*(1.e0 + dinv[9]*dt)
                )))))));
        expt *= f[9-1];
        // F[4]-F[0]
        fmt[5] *= coef;
        expt   *= coef;
        fmt[4] = dinv2[4]*( expt + t2*fmt[5] );
        fmt[3] = dinv2[3]*( expt + t2*fmt[4] );
        fmt[2] = dinv2[2]*( expt + t2*fmt[3] );
        fmt[1] = dinv2[1]*( expt + t2*fmt[2] );
        fmt[0] =            expt + t2*fmt[1];
    }
}

// automatically generated by gen_fmt_method2 function
// m = 6, nexp = 10, eps = 1.0e-12
void fmt6_method2( const double t, const double coef, double fmt[] ) {
    int it;
    double t0, dt, *f, t2, sqrt2, expt;
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
        it = (int)(t*2);
        t0 = delta6 * (double)it + dhalf6;
        t2 = t + t;
        dt = t0 - t;
        f  = &fmt6_table[it*(10+1)];
        // F6(T)
        fmt[6] = f[0] + dt *( f[1] + dt *( f[2] + dt *( f[3] + dt *
                ( f[4] + dt *( f[5] + dt *( f[6] + dt *( f[7] + dt *
                ( f[8] + dt * f[9]))))))));
        // exp(-T)
        expt = 1.e0 + dt*(1.e0 + dinv[2]*dt*(1.e0 + dinv[3]*dt
                *(1.e0 + dinv[4]*dt*(1.e0 + dinv[5]*dt*(1.e0 + dinv[6]*dt
                *(1.e0 + dinv[7]*dt*(1.e0 + dinv[8]*dt*(1.e0 + dinv[9]*dt
                *(1.e0 + dinv[10]*dt)))))))));
        expt *= f[11-1];
        // F[5]-F[0]
        fmt[6] *= coef;
        expt   *= coef;
        fmt[5] = dinv2[5]*( expt + t2*fmt[6] );
        fmt[4] = dinv2[4]*( expt + t2*fmt[5] );
        fmt[3] = dinv2[3]*( expt + t2*fmt[4] );
        fmt[2] = dinv2[2]*( expt + t2*fmt[3] );
        fmt[1] = dinv2[1]*( expt + t2*fmt[2] );
        fmt[0] =            expt + t2*fmt[1];
    }
}

// automatically generated by gen_fmt_method2 function
// m = 7, nexp = 10, eps = 1.0e-12
void fmt7_method2( const double t, const double coef, double fmt[] ) {
    int it;
    double t0, dt, *f, t2, sqrt2, expt;
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
        it = (int)(t*2);
        t0 = delta7 * (double)it + dhalf7;
        t2 = t + t;
        dt = t0 - t;
        f  = &fmt7_table[it*(10+1)];
        // F7(T)
        fmt[7] = f[0] + dt *( f[1] + dt *( f[2] + dt *( f[3] + dt *
                ( f[4] + dt *( f[5] + dt *( f[6] + dt *( f[7] + dt *
                ( f[8] + dt * f[9]))))))));
        // exp(-T)
        expt = 1.e0 + dt*(1.e0 + dinv[2]*dt*(1.e0 + dinv[3]*dt
                *(1.e0 + dinv[4]*dt*(1.e0 + dinv[5]*dt*(1.e0 + dinv[6]*dt
                *(1.e0 + dinv[7]*dt*(1.e0 + dinv[8]*dt*(1.e0 + dinv[9]*dt
                *(1.e0 + dinv[10]*dt)))))))));
        expt *= f[11-1];
        // F[6]-F[0]
        fmt[7] *= coef;
        expt   *= coef;
        fmt[6] = dinv2[6]*( expt + t2*fmt[7] );
        fmt[5] = dinv2[5]*( expt + t2*fmt[6] );
        fmt[4] = dinv2[4]*( expt + t2*fmt[5] );
        fmt[3] = dinv2[3]*( expt + t2*fmt[4] );
        fmt[2] = dinv2[2]*( expt + t2*fmt[3] );
        fmt[1] = dinv2[1]*( expt + t2*fmt[2] );
        fmt[0] =            expt + t2*fmt[1];
    }
}

// automatically generated by gen_fmt_method2 function
// m = 8, nexp = 10, eps = 1.0e-12
void fmt8_method2( const double t, const double coef, double fmt[] ) {
    int it;
    double t0, dt, *f, t2, sqrt2, expt;
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
        it = (int)(t*2);
        t0 = delta8 * (double)it + dhalf8;
        t2 = t + t;
        dt = t0 - t;
        f  = &fmt8_table[it*(10+1)];
        // F8(T)
        fmt[8] = f[0] + dt *( f[1] + dt *( f[2] + dt *( f[3] + dt *
                ( f[4] + dt *( f[5] + dt *( f[6] + dt *( f[7] + dt *
                ( f[8] + dt * f[9]))))))));
        // exp(-T)
        expt = 1.e0 + dt*(1.e0 + dinv[2]*dt*(1.e0 + dinv[3]*dt
                *(1.e0 + dinv[4]*dt*(1.e0 + dinv[5]*dt*(1.e0 + dinv[6]*dt
                *(1.e0 + dinv[7]*dt*(1.e0 + dinv[8]*dt*(1.e0 + dinv[9]*dt
                *(1.e0 + dinv[10]*dt)))))))));
        expt *= f[11-1];
        // F[7]-F[0]
        fmt[8] *= coef;
        expt   *= coef;
        fmt[7] = dinv2[7]*( expt + t2*fmt[8] );
        fmt[6] = dinv2[6]*( expt + t2*fmt[7] );
        fmt[5] = dinv2[5]*( expt + t2*fmt[6] );
        fmt[4] = dinv2[4]*( expt + t2*fmt[5] );
        fmt[3] = dinv2[3]*( expt + t2*fmt[4] );
        fmt[2] = dinv2[2]*( expt + t2*fmt[3] );
        fmt[1] = dinv2[1]*( expt + t2*fmt[2] );
        fmt[0] =            expt + t2*fmt[1];
    }
}

