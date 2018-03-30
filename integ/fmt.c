/**
 * @file fmt.c
 * 各種分子積分計算で用いる誤差関数の計算に関する関数群
 *
 * */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ofmo-def.h"
#ifdef USE_CUDA
//#include "cuda/cudalib.h"
#include "cuda/cuda-fmt.h"
#endif

//#define Free(a) if ( a != NULL ) free( a ); a = NULL

static double FMT_inv2 = 1.0 / 2.0;
static double FMT_inv6 = 1.0 / 6.0;
static int FMT_fmt_m;
static int FMT_fmt_inv_d = 512;	// 2^9
static int FMT_fmt_max_m;
static double FMT_fmt_t;
static int FMT_fmt_n_step;
static double *FMT_fmt_table;
// 変更
double FMT_fmt_step_size;
double FMT_fmt_inv_step_size;
double FMT_pi_div2;

static int FMT_fmt_max_m1;

static int FMT_MAXLQN = -1;

// 各m専用のテーブル
//
double *FMT_fmt_table0;
double *FMT_fmt_table1;
double *FMT_fmt_table2;
double *FMT_fmt_table3;
double *FMT_fmt_table4;
double *FMT_fmt_table5;
double *FMT_fmt_table6;
double *FMT_fmt_table7;
double *FMT_fmt_table8;

/* ---------------------------------------------------------------
    不完全Γ関数（誤差関数）計算の後処理ルーチン
------------------------------------------------------------------ */
static void fmt_finalize() {
    int ret;
    if (FMT_MAXLQN < 0) return;
    if (FMT_fmt_table  != NULL) Free( FMT_fmt_table );
    if (FMT_fmt_table0 != NULL) {
      Free( FMT_fmt_table0 );
      Free( FMT_fmt_table1 );
      Free( FMT_fmt_table2 );
      Free( FMT_fmt_table3 );
      Free( FMT_fmt_table4 );
      Free( FMT_fmt_table5 );
      Free( FMT_fmt_table6 );
      Free( FMT_fmt_table7 );
      Free( FMT_fmt_table8 );
    }
    FMT_fmt_table  = NULL;
    FMT_fmt_table0 = NULL;
    FMT_fmt_table1 = NULL;
    FMT_fmt_table2 = NULL;
    FMT_fmt_table3 = NULL;
    FMT_fmt_table4 = NULL;
    FMT_fmt_table5 = NULL;
    FMT_fmt_table6 = NULL;
    FMT_fmt_table7 = NULL;
    FMT_fmt_table8 = NULL;
    FMT_MAXLQN = -1;
}

/** 誤差関数計算関数の初期化ルーチン
 *
 * 誤差関数\f$ F_m(T)=\int_0^1 t^{2m}\exp[-Tt^2]\f$の計算で用いる
 * 誤差関数テーブルを作成する
 *
 * @param[in] max_lqn 最大軌道量子数
 *
 * @return なし
 *
 * @ingroup integ-fmt
 * */
void fmt_initialize(int max_lqn) {
    double thr_zero = 1.0e-17;
    int i,j,m,nu;
    double eps,t,expt,t2,term,func;
    int it0;

#ifdef USE_CUDA
    if (max_lqn<2) max_lqn = 2;
#endif
    if (max_lqn <= FMT_MAXLQN) return;
    if (FMT_MAXLQN < 0) atexit(fmt_finalize);	// 最初に呼び出されたとき

    if ( FMT_fmt_table != NULL ) free( FMT_fmt_table );
    //fmt_finalize();

    // 基本的な変数の定義
    FMT_fmt_m = 4 * max_lqn + 2;
    FMT_fmt_max_m = FMT_fmt_m + 3;
    FMT_fmt_t = 2 * FMT_fmt_m + 36;
    FMT_fmt_n_step = FMT_fmt_t * FMT_fmt_inv_d;
    FMT_fmt_step_size = FMT_fmt_t / FMT_fmt_n_step;
    FMT_fmt_inv_step_size= 1.0 / FMT_fmt_step_size;
    FMT_fmt_max_m1 = FMT_fmt_max_m + 1;
    FMT_pi_div2 = 2.0 * atan(1.0);

    // 誤差関数計算のテーブルのためのメモリ確保
    FMT_fmt_table = (double*)malloc(sizeof(double)*(FMT_fmt_max_m+1)*(FMT_fmt_n_step+1));
    // added
    if (FMT_fmt_table0==NULL) {
    FMT_fmt_table0 = (double*)malloc(sizeof(double)*(0+4)*(36*FMT_fmt_inv_d+1) );
    FMT_fmt_table1 = (double*)malloc(sizeof(double)*(1+4)*(38*FMT_fmt_inv_d+1) );
    FMT_fmt_table2 = (double*)malloc(sizeof(double)*(2+4)*(40*FMT_fmt_inv_d+1) );
    FMT_fmt_table3 = (double*)malloc(sizeof(double)*(3+4)*(42*FMT_fmt_inv_d+1) );
    FMT_fmt_table4 = (double*)malloc(sizeof(double)*(4+4)*(44*FMT_fmt_inv_d+1) );
    FMT_fmt_table5 = (double*)malloc(sizeof(double)*(5+4)*(46*FMT_fmt_inv_d+1) );
    FMT_fmt_table6 = (double*)malloc(sizeof(double)*(6+4)*(48*FMT_fmt_inv_d+1) );
    FMT_fmt_table7 = (double*)malloc(sizeof(double)*(7+4)*(50*FMT_fmt_inv_d+1) );
    FMT_fmt_table8 = (double*)malloc(sizeof(double)*(8+4)*(52*FMT_fmt_inv_d+1) );
    }

    // 計算開始
	// T=0の場合の計算 Fm(0) = 1/(2m+1)
    for (m=0; m<=FMT_fmt_max_m; m++) FMT_fmt_table[m] = 1.0 / (2*m+1);

	// T>0の場合の計算 Fm(T) (T>0)
    m=FMT_fmt_max_m;
    for (j=1, it0=FMT_fmt_max_m1; j<=FMT_fmt_n_step; j++, it0 += FMT_fmt_max_m1){
	t = FMT_fmt_step_size * j;
	expt = exp(-t);
	nu = 2*m+1;
	t2 = 2.0 * t;
	eps = (expt/t2) * thr_zero;
	term = 1.0 / nu;
	func = term;
	i = nu;
	while (1) {
	    i += 2;
	    term *= t2 / i;
	    func += term;
	    if (term <= eps) break;
	}
	FMT_fmt_table[it0 + m] = expt * func;
	for (i=m-1; i>=0; i--){
	    nu -= 2;
	    FMT_fmt_table[it0 + i] = (expt + t2*FMT_fmt_table[it0 + (i+1)]) / nu;
	}
    }
    // added
    double *dtmp[8+1+1];
    dtmp[0] = FMT_fmt_table0;
    dtmp[1] = FMT_fmt_table1;
    dtmp[2] = FMT_fmt_table2;
    dtmp[3] = FMT_fmt_table3;
    dtmp[4] = FMT_fmt_table4;
    dtmp[5] = FMT_fmt_table5;
    dtmp[6] = FMT_fmt_table6;
    dtmp[7] = FMT_fmt_table7;
    dtmp[8] = FMT_fmt_table8;
    for ( m=0; m<=max_lqn*4; m++ ) {
	for ( j=0; j<=(36+2*m)*FMT_fmt_inv_d; j++ ) {
	    for ( i=0; i<(m+4); i++ ) {
		dtmp[m][j*(m+4)+i] = FMT_fmt_table[j*FMT_fmt_max_m1 + i];
	    }
	}
    }
    // special procedure at m=0
    double c[4];
    c[0] = 1.e0;
    c[1] = 1.e0;
    c[2] = 1.e0 / 2.e0;
    c[3] = 1.e0 / (2.e0*3.e0);
    for ( j=0; j<=36*FMT_fmt_inv_d; j++ ) {
	for ( i=0; i<4; i++ ) FMT_fmt_table0[j*4+i] *= c[i];
    }
    FMT_MAXLQN = max_lqn;
}

#ifdef USE_CUDA
int cuda_fmt_initialize(void)
{
  int ret = 0;
  int max_lqn = FMT_MAXLQN;

#pragma omp master
  {
    double *dtmp[8+1+1];
    size_t mtmp[8+1+1];
    dtmp[0] = FMT_fmt_table0;
    dtmp[1] = FMT_fmt_table1;
    dtmp[2] = FMT_fmt_table2;
    dtmp[3] = FMT_fmt_table3;
    dtmp[4] = FMT_fmt_table4;
    dtmp[5] = FMT_fmt_table5;
    dtmp[6] = FMT_fmt_table6;
    dtmp[7] = FMT_fmt_table7;
    dtmp[8] = FMT_fmt_table8;
    dtmp[9] = FMT_fmt_table;
    for (int m=0; m<=max_lqn*4; m++) mtmp[m] = (m+4)*((36+2*m)*FMT_fmt_inv_d+1);
    mtmp[9] = (FMT_fmt_max_m+1)*(FMT_fmt_n_step+1);
    ret = cuda_FMT_Init(dtmp, mtmp, FMT_fmt_step_size, FMT_fmt_max_m1);
  }
  return ret;
}
#endif

/** 誤差関数を計算する関数
 *
 * \c fmt_initialize 関数で作成された誤差関数テーブルを用いるなどして
 * 誤差関数を計算する
 *
 * @attention
 * @li 事前に初期化ルーチンを \c fmt_initialize を呼び出しておく必要がある
 *
 * @param[out] f[] 計算した誤差関数を代入する配列
 *     （\f$ \tt{f[0]}\sim \tt{f[m]} \f$
 * @param[in] m 誤差関数の次数
 * @param[in] t 誤差関数の引数
 * @param[in] coef 誤差関数に掛ける定数
 *
 * @return なし
 *
 * @ingroup integ-fmt
 * */

void fmt(double f[], const int m, const double t, const double coef){
    int i,ts;
    double d,t_inv,nu;
    int ts0;
    // main
    if (t <= (2*m+36)) {
	ts = 0.5 + t * FMT_fmt_inv_step_size;
	d = ts * FMT_fmt_step_size - t;
	ts0 = ts * FMT_fmt_max_m1;
	for (i=0; i<=m; i++){
	    f[i] = coef * (((FMT_fmt_table[ts0 + i + 3] * FMT_inv6   * d
			   + FMT_fmt_table[ts0 + i + 2] * FMT_inv2 ) * d
			   + FMT_fmt_table[ts0 + i + 1]            ) * d
			   + FMT_fmt_table[ts0 + i + 0]            );
	}
    } else {
	t_inv = FMT_inv2 / t;
	f[0] = coef * sqrt(FMT_pi_div2*t_inv);
	nu = 1.0;
	for (i=1; i<=m; i++){
	    f[i] = t_inv * nu * f[i-1];
	    nu += 2.0;
	}
    }
}


void fmt_(double f[], int *pm, double *pt, double *pcoef) {
    int i,ts, m=*pm;
    double d,t_inv,nu, t=*pt, coef=*pcoef;
    int ts0;
    // main
    if (t <= (2*m+36)) {
	ts = 0.5 + t * FMT_fmt_inv_step_size;
	d = ts * FMT_fmt_step_size - t;
	ts0 = ts * FMT_fmt_max_m1;
	for (i=0; i<=m; i++){
	    f[i] = coef * (((FMT_fmt_table[ts0 + i + 3] * FMT_inv6   * d
			   + FMT_fmt_table[ts0 + i + 2] * FMT_inv2 ) * d
			   + FMT_fmt_table[ts0 + i + 1]            ) * d
			   + FMT_fmt_table[ts0 + i + 0]            );
	}
    } else {
	t_inv = FMT_inv2 / t;
	f[0] = coef * sqrt(FMT_pi_div2*t_inv);
	nu = 1.0;
	for (i=1; i<=m; i++){
	    f[i] = t_inv * nu * f[i-1];
	    nu += 2.0;
	}
    }
}
