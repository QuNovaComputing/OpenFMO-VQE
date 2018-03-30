/**
 * @file ofmo-oneint-core.h
 * １つのCSペアに対する１電子積分を計算する関数群の
 * プロトタイプ宣言を行っているヘッダファイル
 * */
#ifndef _OFMO_ONEINT_CORE_H_
#define _OFMO_ONEINT_CORE_H_

extern void oneint_core_ss_(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const int atomic_number[],
	double OVI[], double HCORE[] );
extern void oneint_core_ps_(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const int atomic_number[],
	double OVI[], double HCORE[] );
extern void oneint_core_pp_(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const int atomic_number[],
	double OVI[], double HCORE[] );
extern void oneint_core_ds_(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const int atomic_number[],
	double OVI[], double HCORE[] );
extern void oneint_core_dp_(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const int atomic_number[],
	double OVI[], double HCORE[] );
extern void oneint_core_dd_(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const int atomic_number[],
	double OVI[], double HCORE[] );

extern void oneint_core_ss__(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const int atomic_number[],
	double OVI[], double HCORE[] );
extern void oneint_core_ps__(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const int atomic_number[],
	double OVI[], double HCORE[] );
extern void oneint_core_pp__(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const int atomic_number[],
	double OVI[], double HCORE[] );
extern void oneint_core_ds__(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const int atomic_number[],
	double OVI[], double HCORE[] );
extern void oneint_core_dp__(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const int atomic_number[],
	double OVI[], double HCORE[] );
extern void oneint_core_dd__(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const int atomic_number[],
	double OVI[], double HCORE[] );

#endif
