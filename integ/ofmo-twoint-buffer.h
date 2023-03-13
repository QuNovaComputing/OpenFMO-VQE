/**
 * @file ofmo-twoint-buffer.h
 * 各タイプの２電子積分を計算して、バッファに保存する処理を行う
 * 関数群のプロトタイプ宣言を行っているヘッダファイル。
 * */
#ifndef _OFMO_TWOINT_BUFFER_H_
#define _OFMO_TWOINT_BUFFER_H_

/* one underscore version */
extern int ofmo_twoint_buffer_ssss_(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_psss_(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_psps_(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_ppss_(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_ppps_(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_pppp_(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_dsss_(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_dsps_(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_dspp_(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_dsds_(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_dpss_(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_dpps_(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_dppp_(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_dpds_(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_dpdp_(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_ddss_(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_ddps_(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_ddpp_(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_ddds_(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_dddp_(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_dddd_(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

/* two-underscore version */
extern int ofmo_twoint_buffer_ssss__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_psss__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_psps__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_ppss__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_ppps__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_pppp__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_dsss__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_dsps__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_dspp__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_dsds__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_dpss__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_dpps__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_dppp__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_dpds__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_dpdp__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_ddss__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_ddps__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_ddpp__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_ddds__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_dddp__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

extern int ofmo_twoint_buffer_dddd__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs);

#endif
