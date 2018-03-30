/**
 * @file ofmo-oneint.h
 *
 * 通常のHartree-Fock SCF計算に現れる１電子積分（重なり行列、
 * １電子ハミルトン行列）を計算するための関数群の
 * プロトタイプ宣言を行うヘッダファイル
 *
 * */
#ifndef _OFMO_ONEINT_H_
#define _OFMO_ONEINT_H_

extern int ofmo_oneint_init();

extern int ofmo_oneint_set_sum_atomic_numbers(const int sum_atomic_number);

extern int ofmo_oneint_ss_(
	const int *pnworkers, const int *pworkerid,
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[],
	const int shel_add[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	const int *pnat, const int atomic_number[],
	double S[], double H[] );

extern int ofmo_oneint_ps_(
	const int *pnworkers, const int *pworkerid,
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[],
	const int shel_add[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	const int *pnat, const int atomic_number[],
	double S[], double H[] );

extern int ofmo_oneint_pp_(
	const int *pnworkers, const int *pworkerid,
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[],
	const int shel_add[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	const int *pnat, const int atomic_number[],
	double S[], double H[] );

extern int ofmo_oneint_ds_(
	const int *pnworkers, const int *pworkerid,
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[],
	const int shel_add[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	const int *pnat, const int atomic_number[],
	double S[], double H[] );

extern int ofmo_oneint_dp_(
	const int *pnworkers, const int *pworkerid,
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[],
	const int shel_add[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	const int *pnat, const int atomic_number[],
	double S[], double H[] );

extern int ofmo_oneint_dd_(
	const int *pnworkers, const int *pworkerid,
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[],
	const int shel_add[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	const int *pnat, const int atomic_number[],
	double S[], double H[] );

extern int ofmo_oneint_ss__(
	const int *pnworkers, const int *pworkerid,
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[],
	const int shel_add[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	const int *pnat, const int atomic_number[],
	double S[], double H[] );

extern int ofmo_oneint_ps__(
	const int *pnworkers, const int *pworkerid,
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[],
	const int shel_add[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	const int *pnat, const int atomic_number[],
	double S[], double H[] );

extern int ofmo_oneint_pp__(
	const int *pnworkers, const int *pworkerid,
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[],
	const int shel_add[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	const int *pnat, const int atomic_number[],
	double S[], double H[] );

extern int ofmo_oneint_ds__(
	const int *pnworkers, const int *pworkerid,
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[],
	const int shel_add[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	const int *pnat, const int atomic_number[],
	double S[], double H[] );

extern int ofmo_oneint_dp__(
	const int *pnworkers, const int *pworkerid,
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[],
	const int shel_add[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	const int *pnat, const int atomic_number[],
	double S[], double H[] );

extern int ofmo_oneint_dd__(
	const int *pnworkers, const int *pworkerid,
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[],
	const int shel_add[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	const int *pnat, const int atomic_number[],
	double S[], double H[] );
#endif
