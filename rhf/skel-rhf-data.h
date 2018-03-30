#ifndef _SKEL_RHF_DATA_H_
#define _SKEL_RHF_DATA_H_

#include "ofmo-def.h"

struct skel_rhf_data_st {
  int natom;
  int charge;
  int* atomic_number;
  double* atom_x;
  double* atom_y;
  double* atom_z;
  char basis_name[MAXSTRLEN];
  int nelec;
  int nocc;

  int maxlqn;
  int ncs;
  int nao;
  int nps;
  int* leading_cs;
  int* shel_tem;
  int* shel_atm;
  int* shel_add;
  int* shel_ini;
  double* prim_exp;
  double* prim_coe;
  int* sao2uao;
  int npspair;
  int* leading_cs_pair;
  int* csp_leading_ps_pair;
  int* csp_ics;
  int* csp_jcs;
  double* csp_schwarz;
  double* psp_zeta;
  double* psp_dkps;
  double* psp_xiza;
};
typedef struct skel_rhf_data_st skel_rhf_data_t;

int read_packed_matrix_binary(const int nao, const char *filename,
            double *ap);
int skel_rhf_init_data( skel_rhf_config_t* config, skel_rhf_data_t* data );
void skel_rhf_fin_data( skel_rhf_data_t* data );
int skel_rhf_cutoff_make_table(
        skel_rhf_config_t* config, skel_rhf_data_t* data );

#endif /* _SKEL_RHF_DATA_H_ */

