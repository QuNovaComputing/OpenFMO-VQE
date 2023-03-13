#ifndef _FMT_METHOD3_H_
#define _FMT_METHOD3_H_
#include <stdio.h>

extern double fmt_method3_get_eps();

extern int fmt_method3_get_mmax();

extern int fmt_method3_get_nexp();

extern int fmt_method3_get_ndiv( const int m );

extern int fmt_method3_get_tmax( const int m );

extern size_t fmt_method3_get_size( const int m );

extern int fmt_method3_init();

extern void fmt0_method3(const double t, const double coef, double fmt[]);
extern void fmt1_method3(const double t, const double coef, double fmt[]);
extern void fmt2_method3(const double t, const double coef, double fmt[]);
extern void fmt3_method3(const double t, const double coef, double fmt[]);
extern void fmt4_method3(const double t, const double coef, double fmt[]);
extern void fmt5_method3(const double t, const double coef, double fmt[]);
extern void fmt6_method3(const double t, const double coef, double fmt[]);
extern void fmt7_method3(const double t, const double coef, double fmt[]);
extern void fmt8_method3(const double t, const double coef, double fmt[]);

static void (*fmt_method3[])
        ( const double t, const double coef, double fmt[] ) = {
    fmt0_method3,
    fmt1_method3,
    fmt2_method3,
    fmt3_method3,
    fmt4_method3,
    fmt5_method3,
    fmt6_method3,
    fmt7_method3,
    fmt8_method3,
};

#endif
