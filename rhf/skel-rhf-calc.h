#ifndef _SKEL_RHF_CALC_H_
#define _SKEL_RHF_CALC_H_

#include "ofmo-def.h"

struct skel_rhf_array_st {
  double* D;
  double* S;
  double* H;
  double* C;
  double* moe;
  double* TEMP;
  double* aopop;
  double* atpop;
};
typedef struct skel_rhf_array_st skel_rhf_array_t;

int skel_rhf_init_array(
    skel_rhf_config_t* config, skel_rhf_data_t* data,
    skel_rhf_array_t* array );
void skel_rhf_fin_array( skel_rhf_array_t* array );
int skel_rhf_oneint(
    skel_rhf_config_t* config, skel_rhf_data_t* data,
    skel_rhf_array_t* array );
int skel_rhf_init_density(
    skel_rhf_config_t* config, skel_rhf_data_t* data,
    skel_rhf_array_t* array );
int skel_rhf_twoint_first(
    skel_rhf_config_t* config, skel_rhf_data_t* data,
    skel_rhf_array_t* array );
double skel_rhf_scf(
    skel_rhf_config_t* config, skel_rhf_data_t* data,
    skel_rhf_array_t* array );

#endif /* _SKEL_RHF_CALC_H_ */

