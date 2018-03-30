// driver routines for fmt funcs created by gen-fmt

#include <stdlib.h>
#include "fmt.h"
#include "fmt-m.h"
#ifdef USE_CUDA
#include "cuda/cuda-fmt-drv.h"
#endif

int fmt_m_init(void) {
  int ret=-1;
  ret = fmt_method1_init();
  if (ret<0) return -1;
  ret = fmt_method2_init();
  if (ret<0) return -2;
  ret = fmt_method3_init();
  if (ret<0) return -3;

#ifdef USE_CUDA
#pragma omp master
  ret = cuda_fmt_m_init();
  if (ret<0) return -4;
//  atexit( cuda_fmt_m_finalize() );
#endif

  return 0;
};

void (*fmt_m[]) ( const double t, const double coef, double f[] ) = {
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

void fmt_mm(const int m, const double t, const double coef, double f[] )
{
  fmt_m[m](t, coef, f);
#if 0
  if (m==0) fmt(f, 0, t, coef);
  else fmt_m[m](t, coef, f);
#endif
};

