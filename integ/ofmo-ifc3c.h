#ifndef _OFMO_IFC3C_H_
#define _OFMO_IFC3C_H_

extern  int ofmo_ifc3c_os_init();
extern  int ofmo_ifc3c_rys_init();

extern int ofmo_ifc3c_os_ssss(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_os_sspp(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_os_ssdd(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_os_psss(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_os_pspp(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_os_psdd(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_os_ppss(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_os_pppp(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_os_ppdd(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_os_dsss(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_os_dspp(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_os_dsdd(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_os_dpss(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_os_dppp(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_os_dpdd(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_os_ddss(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_os_ddpp(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_os_dddd(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );
	//

extern int ofmo_ifc3c_rys_ssss(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_rys_sspp(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_rys_ssdd(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_rys_psss(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_rys_pspp(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_rys_psdd(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_rys_ppss(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_rys_pppp(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_rys_ppdd(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_rys_dsss(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_rys_dspp(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_rys_dsdd(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_rys_dpss(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_rys_dppp(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_rys_dpdd(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_rys_ddss(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_rys_ddpp(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );

extern int ofmo_ifc3c_rys_dddd(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] );
#endif
