#ifndef _CUDA_IFC4C_CALC_H_
#define _CUDA_IFC4C_CALC_H_

#ifdef __cplusplus
extern "C" {
#endif
extern int dim_ifc4c[][2];

/* ------------------------------------- */
void cuda_print_wifc4c(void);

int cuda_ifc4c_calc(const int idev,
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
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
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] );
/* ------------------------------------- */


#ifdef __cplusplus
}
#endif

#endif /* _CUDA_IFC4C_CALC_H_ */
