/**
 * @file ofmo-ifc2c.h
 * ２中心クーロン相互作用項を計算する中位の関数群の
 * プロトタイプ宣言を行っているヘッダファイル
 *
 * */
#ifndef _OFMO_IFC2C_H_
#define _OFMO_IFC2C_H_

extern int ofmo_ifc2c_init();

extern int ofmo_ifc2c_ss__(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb,
	// basis set data for fragment
	const int leading_cs_frg[],
	const int shel_tem_frg[], const int shel_atm_frg[],
	const int shel_add_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[],
	const double prim_exp_frg[], const double prim_coe_frg[],
	// atomic charge and atomic coordinate data for monomer
	const int *pnat_mon, const double atom_x_mon[],
	const double atom_y_mon[], const double atom_z_mon[],
	const double atm_pop_mon[],
	// output
	double V_frg[] );

extern int ofmo_ifc2c_ps__(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb,
	// basis set data for fragment
	 const int leading_cs_frg[],
	const int shel_tem_frg[], const int shel_atm_frg[],
	const int shel_add_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[],
	const double prim_exp_frg[], const double prim_coe_frg[],
	// atomic charge and atomic coordinate data for monomer
	const int *pnat_mon, const double atom_x_mon[],
	const double atom_y_mon[], const double atom_z_mon[],
	const double atm_pop_mon[],
	// output
	double V_frg[] );

extern int ofmo_ifc2c_pp__(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb,
	// basis set data for fragment
	const int leading_cs_frg[],
	const int shel_tem_frg[], const int shel_atm_frg[],
	const int shel_add_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[],
	const double prim_exp_frg[], const double prim_coe_frg[],
	// atomic charge and atomic coordinate data for monomer
	const int *pnat_mon, const double atom_x_mon[],
	const double atom_y_mon[], const double atom_z_mon[],
	const double atm_pop_mon[],
	// output
	double V_frg[] );

extern int ofmo_ifc2c_ds__(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb,
	// basis set data for fragment
	const int leading_cs_frg[],
	const int shel_tem_frg[], const int shel_atm_frg[],
	const int shel_add_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[],
	const double prim_exp_frg[], const double prim_coe_frg[],
	// atomic charge and atomic coordinate data for monomer
	const int *pnat_mon, const double atom_x_mon[],
	const double atom_y_mon[], const double atom_z_mon[],
	const double atm_pop_mon[],
	// output
	double V_frg[] );

extern int ofmo_ifc2c_dp__(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb,
	// basis set data for fragment
	const int leading_cs_frg[],
	const int shel_tem_frg[], const int shel_atm_frg[],
	const int shel_add_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[],
	const double prim_exp_frg[], const double prim_coe_frg[],
	// atomic charge and atomic coordinate data for monomer
	const int *pnat_mon, const double atom_x_mon[],
	const double atom_y_mon[], const double atom_z_mon[],
	const double atm_pop_mon[],
	// output
	double V_frg[] );

extern int ofmo_ifc2c_dd__(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb,
	// basis set data for fragment
	const int leading_cs_frg[],
	const int shel_tem_frg[], const int shel_atm_frg[],
	const int shel_add_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[],
	const double prim_exp_frg[], const double prim_coe_frg[],
	// atomic charge and atomic coordinate data for monomer
	const int *pnat_mon, const double atom_x_mon[],
	const double atom_y_mon[], const double atom_z_mon[],
	const double atm_pop_mon[],
	// output
	double V_frg[] );

#endif
