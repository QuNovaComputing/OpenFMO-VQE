/**
 * @file ofmo-cutoff-core.h
 * １つのSchwarz積分を計算するための関数のプロトタイプ宣言を行っている
 * ファイル。
 * */
#ifndef _OFMO_CUTOFF_CORE_H_
#define _OFMO_CUTOFF_CORE_H_

extern double schwarz_core_ssss_(
	const int npps, const double vzeta[], const double vdkps[],
	const double vxiza[], const double BA[3], const double AB2 );
extern double schwarz_core_psps_(
	const int npps, const double vzeta[], const double vdkps[],
	const double vxiza[], const double BA[3], const double AB2 );
extern double schwarz_core_pppp_(
	const int npps, const double vzeta[], const double vdkps[],
	const double vxiza[], const double BA[3], const double AB2 );
extern double schwarz_core_dsds_(
	const int npps, const double vzeta[], const double vdkps[],
	const double vxiza[], const double BA[3], const double AB2 );
extern double schwarz_core_dpdp_(
	const int npps, const double vzeta[], const double vdkps[],
	const double vxiza[], const double BA[3], const double AB2 );
extern double schwarz_core_dddd_(
	const int npps, const double vzeta[], const double vdkps[],
	const double vxiza[], const double BA[3], const double AB2 );

#endif
