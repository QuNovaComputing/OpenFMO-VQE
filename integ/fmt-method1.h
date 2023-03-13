#ifndef _FMT_METHOD1_H_
#define _FMT_METHOD1_H_
#include <stdio.h>

extern double fmt_method1_get_eps();

extern int fmt_method1_get_mmax();

extern int fmt_method1_get_nexp();

extern int fmt_method1_get_ndiv( const int m );

extern int fmt_method1_get_tmax( const int m );

extern size_t fmt_method1_get_size( const int m );

extern int fmt_method1_init();

extern void fmt0_method1(const double t, const double coef, double fmt[]);
extern void fmt1_method1(const double t, const double coef, double fmt[]);
extern void fmt2_method1(const double t, const double coef, double fmt[]);
extern void fmt3_method1(const double t, const double coef, double fmt[]);
extern void fmt4_method1(const double t, const double coef, double fmt[]);
extern void fmt5_method1(const double t, const double coef, double fmt[]);
extern void fmt6_method1(const double t, const double coef, double fmt[]);
extern void fmt7_method1(const double t, const double coef, double fmt[]);
extern void fmt8_method1(const double t, const double coef, double fmt[]);

static void (*fmt_method1[])
        ( const double t, const double coef, double fmt[] ) = {
    fmt0_method1,
    fmt1_method1,
    fmt2_method1,
    fmt3_method1,
    fmt4_method1,
    fmt5_method1,
    fmt6_method1,
    fmt7_method1,
    fmt8_method1,
};

#endif
