/**
 * @file ofmo-diis.h
 *
 * DIISに関する関数群のプロトタイプ宣言を行っている
 * ヘッダファイル
 * */
#ifndef _OFMO_DIIS_H_
#define _OFMO_DIIS_H_

extern int ofmo_diis_alloc( const int maxdiis, const int maxnao );
extern int ofmo_diis_init( const int nao, const double Sp[] );
extern int ofmo_diis_update( const int nao, double Dp[], const double Fp[],
	double Fdiis[], const int dodiis );
extern double ofmo_diis_profiling(
	const int nao, const double Dp[], const double Fp[]);

#endif
