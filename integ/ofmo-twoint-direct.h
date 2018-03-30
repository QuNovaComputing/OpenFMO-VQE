/**
 * @file ofmo-twoint-direct.h
 * ２電子積分計算、および、２電子ハミルトン行列計算を
 * 行う関数群のプロトタイプ定義をしているヘッダファイル
 * */
#ifndef _OFMO_TWOINT_DIRECT_H_
#define _OFMO_TWOINT_DIRECT_H_

extern int ofmo_twoint_direct_ssss__(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_psss__(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_psps__(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_ppss__(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_ppps__(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_pppp__(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_dsss__(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_dsps__(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_dspp__(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_dsds__(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_dpss__(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_dpps__(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_dppp__(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_dpds__(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_dpdp__(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_ddss__(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_ddps__(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_ddpp__(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_ddds__(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_dddp__(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_dddd__(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

// Fortran番
extern int ofmo_twoint_direct_ssss_(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_psss_(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_psps_(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_ppss_(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_ppps_(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_pppp_(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_dsss_(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_dsps_(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_dspp_(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_dsds_(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_dpss_(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_dpps_(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_dppp_(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_dpds_(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_dpdp_(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_ddss_(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_ddps_(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_ddpp_(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_ddds_(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_dddp_(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_direct_dddd_(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );
#define OFMO_EBUF_FULL   1
#define OFMO_EBUF_NOFULL 0

#endif
