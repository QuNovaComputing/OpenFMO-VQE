#ifndef _FMT_METHOD2_H_
#define _FMT_METHOD2_H_
#include <stdio.h>
extern double fmt_method2_get_eps();

extern int fmt_method2_get_mmax();

extern int fmt_method2_get_nexp();

extern int fmt_method2_init();

extern int fmt_method2_get_ndiv( const int m );

extern int fmt_method2_get_tmax( const int m );

extern size_t fmt_method2_get_size( const int m );

extern void fmt0_method2( const double t, const double coef, double fmt[] );
extern void fmt1_method2( const double t, const double coef, double fmt[] );
extern void fmt2_method2( const double t, const double coef, double fmt[] );
extern void fmt3_method2( const double t, const double coef, double fmt[] );
extern void fmt4_method2( const double t, const double coef, double fmt[] );
extern void fmt5_method2( const double t, const double coef, double fmt[] );
extern void fmt6_method2( const double t, const double coef, double fmt[] );
extern void fmt7_method2( const double t, const double coef, double fmt[] );
extern void fmt8_method2( const double t, const double coef, double fmt[] );

static void (*fmt_method2[])
        ( const double t, const double coef, double fmt[] ) = {
    fmt0_method2,
    fmt1_method2,
    fmt2_method2,
    fmt3_method2,
    fmt4_method2,
    fmt5_method2,
    fmt6_method2,
    fmt7_method2,
    fmt8_method2,
};

#endif
