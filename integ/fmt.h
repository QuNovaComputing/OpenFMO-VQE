/**
 * @file fmt.h
 * 誤差関数計算に関する関数のプロトタイプ宣言を行っている
 * ヘッダファイル
 *
 * */
#ifndef _FMT_H_
#define _FMT_H_

extern void fmt_initialize( int maxlqn );
extern void fmt(double FmT[], const int m, const double T, const double COEF);

#endif
