/**
 * @file ofmo-cutoff.h
 * 各タイプのカットオフテーブルを計算するための関数群の
 * プロトタイプ宣言を行っているファイル。
 *
 * */
#ifndef _OFMO_CUTOFF_H_
#define _OFMO_CUTOFF_H_

extern int ofmo_cutoff_init();

extern int ofmo_cutoff_ss_(
	const int *pLa, const int *pLb, const int leading_cs[], 
	const int shel_tem[], const int shel_atm[], const int shel_add[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	int leading_cs_pair[],
	double csp_schwarz[], int csp_ics[], int csp_jcs[],
	int csp_leading_ps_pair[],
	double psp_zeta[], double psp_dkps[], double psp_xiza[] );

extern int ofmo_cutoff_ps_(
	const int *pLa, const int *pLb, const int leading_cs[], 
	const int shel_tem[], const int shel_atm[], const int shel_add[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	int leading_cs_pair[],
	double csp_schwarz[], int csp_ics[], int csp_jcs[],
	int csp_leading_ps_pair[],
	double psp_zeta[], double psp_dkps[], double psp_xiza[] );

extern int ofmo_cutoff_pp_(
	const int *pLa, const int *pLb, const int leading_cs[], 
	const int shel_tem[], const int shel_atm[], const int shel_add[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	int leading_cs_pair[],
	double csp_schwarz[], int csp_ics[], int csp_jcs[],
	int csp_leading_ps_pair[],
	double psp_zeta[], double psp_dkps[], double psp_xiza[] );

extern int ofmo_cutoff_ds_(
	const int *pLa, const int *pLb, const int leading_cs[], 
	const int shel_tem[], const int shel_atm[], const int shel_add[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	int leading_cs_pair[],
	double csp_schwarz[], int csp_ics[], int csp_jcs[],
	int csp_leading_ps_pair[],
	double psp_zeta[], double psp_dkps[], double psp_xiza[] );

extern int ofmo_cutoff_dp_(
	const int *pLa, const int *pLb, const int leading_cs[], 
	const int shel_tem[], const int shel_atm[], const int shel_add[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	int leading_cs_pair[],
	double csp_schwarz[], int csp_ics[], int csp_jcs[],
	int csp_leading_ps_pair[],
	double psp_zeta[], double psp_dkps[], double psp_xiza[] );

extern int ofmo_cutoff_dd_(
	const int *pLa, const int *pLb, const int leading_cs[], 
	const int shel_tem[], const int shel_atm[], const int shel_add[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	int leading_cs_pair[],
	double csp_schwarz[], int csp_ics[], int csp_jcs[],
	int csp_leading_ps_pair[],
	double psp_zeta[], double psp_dkps[], double psp_xiza[] );

/* ------------------------------------- */
extern int ofmo_cutoff_show_table( const int maxlqn,
        const int leading_cs[], const int shel_tem[],
        const int leading_cs_pair[],
        const int csp_leading_ps_pair[] );
/* ------------------------------------- */
/* sort csp */
#if defined(SORT_CSP_SCHWARZ)||defined(SORT_CSP_REV)
#ifndef SORT_CSP
#define SORT_CSP
#endif
#endif

#ifdef SORT_CSP

struct ofmo_cutoff_sort_t {
  int ic;
  int iv;
  double vab;
};

extern int ofmo_cutoff_sort_(
    const int La, const int Lb, const int leading_cs[],
    const int shel_tem[], const int shel_atm[], const int shel_add[],
    const double atom_x[], const double atom_y[],
    const double atom_z[],
    const double prim_exp[], const double prim_coe[],
    int leading_cs_pair[],
    double csp_schwarz[], int csp_ics[], int csp_jcs[],
    int csp_leading_ps_pair[],
    double psp_zeta[], double psp_dkps[], double psp_xiza[] );

#endif /* SORT_CSP */

#endif
