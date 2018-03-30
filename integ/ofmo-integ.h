/**
 * @file ofmo-integ.h
 * 分子積分計算の最上位APIのプロトタイプ宣言を行っている
 * ヘッダファイル
 *
 * */
#ifndef _OFMO_INTEG_H_
#define _OFMO_INTEG_H_

extern int ofmo_integ_init( const int maxlqn );

/*  **** cutoff table ****
 * This function is thread-safe if output arrayes
 * (such as leading_cs_pair[]) are not shared by multiple threads */
extern int ofmo_cutoff_make_table(
	// input arguments
	const int maxlqn, const int leading_cs[],
	const int shel_tem[], const int shel_atm[], const int shel_add[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	// output arguments
	int leading_cs_pair[],
	double csp_schwarz[], int csp_ics[], int csp_jcs[],
	int csp_leading_ps_pair[],
	double psp_zeta[], double psp_dkps[], double psp_xiza[] );

/*  **** one-electron integral ****
 *  This function is thread-parallelized code.
 *  Overlap matrix (S), 1-electron Hamiltonian matrix (Hcore) are to be
 *  calculated in this routine.
 *  keyword: static load balancing
 * */
extern int ofmo_integ_oneint_sorted(
	// parallelization
	const int nworkers, const int workerid,
	// basis set data
	const int maxlqn, const int leading_cs[],
	const int shel_tem[], const int shel_atm[], const int shel_add[],
	const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	// atom data
	const int nat, const int atomic_number[],
	// output
	double S[], double H[]);

/*  **** two-electron integral ****
 *  (1) ofmo_integ_twoint_first
 *  This function is hybrid parallelized code and to be called from
 *  thread-parallel region.
 *  This function calculate electron repulsion integrals
 *  and stored in buffer, but don't calculate Fock matrix.
 * */
extern int ofmo_integ_twoint_first(
	// parallelization
	const int nworkers, const int workerid,
	// buffer size
	const size_t ebuf_buffer_size_mb,
	// basis set & cutoff table data
	const int maxlqn, const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const int leading_cs_pair[], const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[] );

/*  **** two-electron integral ****
 *  (2) ofmo_integ_gen_gmat
 *  This function is hybrid parallelized code and to be called from
 *  non-thread-parallel region.
 *  This function calculate ERIs not calculated in ofmo_integ_twoint_first
 *  function, and generate Fock matrix using stored and calculated ERIs.
 *  */
extern int ofmo_integ_gen_gmat(
	// parallelization
	const int nworkers, const int workerid,
	// basis set & cutoff table data
	const int maxlqn, const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const int leading_cs_pair[],
	const int leading_cs[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// density matrix data & G-matrix (output)
	const int nao, const double D[], double G[] );

/*  **** inter-fragment 4-center Clulomb integral (ifc4c) ****
 *  This function is hybrid parallelized code and to be called from
 *  thread-parallel region.
 * */
extern int ofmo_integ_ifc4c_sorted_partial(
	// parallelization
	const int nworkers, const int workerid,
	// basis and cutoff table data for fragment
	const int maxlqn,
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
        //
	const int leading_cs_mon[],
	// density matrix of monomer
	const int nao_mon, const double D_mon[],
	// (output) Coulomb potential
	double V[] );

/*  **** inter-fragment 3-center Clulomb integral (ifc3c) ****
 *  This function is hybrid parallelized code and to be called from
 *  thread-parallel region.
 * */
extern int ofmo_integ_ifc3c_sorted_partial(
	// parallelization
	const int nworkers, const int workerid,
	// basis and cutoff table data for fragment
	const int maxlqn,
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const int leading_cs_pair[],
	//const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
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
	double V[] );

/*  **** inter-fragment 2-center Clulomb integral (ifc2c) ****
 *  This function is hybrid parallelized code and to be called from
 *  thread-parallel region.
 * */
extern int ofmo_integ_ifc2c_sorted_partial(
	// parallelization
	const int nworkers, const int workerid,
	// input data of fragment
	const int maxlqn, const int leading_cs[],
	const int shel_tem[], const int shel_atm[],
	const int shel_add[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	// input data of counter monomer
	const int nat_mon, const double atom_x_mon[],
	const double atom_y_mon[], const double atom_z_mon[],
	const double atm_pop_mon[],
	// output data
	double V[] );

/*extern int ofmo_integ_ifc2c_sorted_partial_g(
	// parallelization
	const int nworkers, const int workerid,
	// for precise load-balancing
	const int type,
	// input data of fragment
	const int maxlqn, const int leading_cs_frg[],
	const int shel_tem_frg[], const int shel_atm_frg[],
	const int shel_add_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[],
	const double prim_exp_frg[], const double prim_coe_frg[],
	// input data of counter monomer
	const int nat_mon, const double atom_x_mon[],
	const double atom_y_mon[], const double atom_z_mon[],
	const double atm_pop_mon[],
	// output data
	double V_frg[] );*/
/*
   functions to control load-balancing
   */
extern size_t ofmo_integ_get_loop_offset( const int mythread );
extern void ofmo_integ_set_loop_offset( const int mythread,
	const size_t offset );
extern void ofmo_integ_set_target_type( const int mythread,
	const int ttype );

/*  **** initialize function ****
 *  This function initialize all integral codes to be provided by
 *  this library.
 *  This function must be called from non-thread-parallel region
 *  before any other integral code is called.
 * */
//extern int ofmo_integ_init();
/*extern int ofmo_integ_twoint_first_g(
	// parallelization
	const int nworkers, const int workerid,
	// added
	int *OFFSET,
	// buffer size
	const size_t ebuf_buffer_size_mb,
	// basis set & cutoff table data
	const int maxlqn, const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const int leading_cs_pair[], const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[] );

extern int ofmo_integ_ifc4c_sorted_partial_g(
	// parallelization
	const int nworkers, const int workerid,
	// added
	int *OFFSET,
	// basis and cutoff table data for fragment
	const int maxlqn,
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const int nao_mon, const double D_mon[],
	// (output) Coulomb potential
	double V[] );

extern int ofmo_integ_ifc3c_sorted_partial_g(
	// parallelization
	const int nworkers, const int workerid,
	// added
	int *OFFSET,
	// basis and cutoff table data for fragment
	const int maxlqn,
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const int leading_cs_pair[],
	//const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
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
	double V[] );*/
#endif
