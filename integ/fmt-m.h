// driver routines for fmt funcs created by gen-fmt

#ifndef _FMT_M_H_
#define _FMT_M_H_

#include "fmt-method1.h"
#include "fmt-method2.h"
#include "fmt-method3.h"

extern int fmt_m_init(void);

//void (*fmt_m[]) ( const double t, const double coef, double f[] );
/*
static void (*fmt_m[]) ( const double t, const double coef, double f[] ) = {
  fmt0_method1,
  fmt1_method1,
  fmt2_method1,
  fmt3_method1,
  fmt4_method1,
  fmt5_method1,
  fmt6_method3,
  fmt7_method3,
  fmt8_method3,
};
*/

void fmt_mm( const int m, const double t, const double coef, double f[] );

//#define OFMO_FMT(a,b,c,d) fmt((a), (b), (c), (d))
//#define OFMO_FMT(a,b,c,d) fmt_m[(b)]((c), (d), (a))
//#define OFMO_FMT(a,b,c,d) fmt_mm((b), (c), (d), (a))
//#define OFMO_FMT(a,b,c,d) {if ((b)==0) fmt((a), 0, (c), (d)); else fmt_m[(b)]((c), (d), (a));}
//#define OFMO_FMT(a,b,c,d) {if ((b)==0) fmt((a), 0, (c), (d)); else fmt_mm((b), (c), (d), (a));}
#define OFMO_FMT(a,b,c,d) {if ((b)<=5) fmt##b##_method1((c), (d), (a)); else fmt##b##_method3((c), (d), (a));}

#endif /* _FMTM_H_ */
