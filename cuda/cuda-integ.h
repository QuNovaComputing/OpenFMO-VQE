#ifndef _CUDA_INTEG_H_
#define _CUDA_INTEG_H_

#ifdef __cplusplus
extern "C" {
#endif
extern int dim2e[][2];
extern char *cuda_s2e[];

int cuda_twoint_direct_ssss_(
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

int cuda_twoint_direct_psss_(
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

/* ------------------------------------- */
typedef enum {
  INTTYPE_QUERY    =-1,
  INTTYPE_OS_XXXX  = 0,
  INTTYPE_OS       = 1,
  INTTYPE_RYS_XXXX = 2,
  INTTYPE_RYS      = 3
} inttype;

inttype cuda_twoint_inttype(inttype type);

/* ------------------------------------- */
int cuda_get_num_types(void);
void cuda_print_w2e(void);

int cuda_calc_twoint_direct
(
        const int Labcd,
        // paralleization
        const int nworkers, const int workerid,
        // integral type data
        const int La, const int Lb, const int Lc, const int Ld,
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
        const long max_nzeri, long *p_nzeri,
        double etmp_val[], short int etmp_ind4[],
        const int last_ijcs, const int last_klcs,
        // density matrix & G-matrix data
        const int nao, const double Ds[], double G[]);

#ifdef __cplusplus
}
#endif

#endif /* _CUDA_INTEG_H_ */
