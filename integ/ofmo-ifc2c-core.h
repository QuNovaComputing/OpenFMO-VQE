/**
 * @file ofmo-ifc2c-core.h
 * １つのCSペアに対する２中心クーロン積分を計算する関数群の
 * プロトタイプ宣言を行っているヘッダファイル
 * */
#ifndef _OFMO_IFC2C_CORE_H_
#define _OFMO_IFC2C_CORE_H_ 

extern int ifc2c_core_ss__(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const double atm_pop[],
	double SS[1] ); 

extern int ifc2c_core_ps__(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const double atm_pop[],
	double PS[3] ); 

extern int ifc2c_core_pp__(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const double atm_pop[],
	double PP[3*3] ); 

extern int ifc2c_core_ds__(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const double atm_pop[],
	double DS[6] ); 

extern int ifc2c_core_dp__(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const double atm_pop[],
	double DP[6*3] ); 

extern int ifc2c_core_dd__(
	const int *ips0, const int *nps_i, const double A[3],
	const int *jps0, const int *nps_j, const double B[3],
	const double prim_exp[], const double prim_coe[],
	const int *nat, const double atom_x[], const double atom_y[],
	const double atom_z[], const double atm_pop[],
	double DD[6*6] ); 

#endif 
