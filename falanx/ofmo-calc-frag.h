
#ifndef _OFMO_CALC_FRAGMENT_H_
#define _OFMO_CALC_FRAGMENT_H_

#include <mpi.h>

extern int ofmo_monomer_init_density(const int* imsg, MPI_Comm comm);
extern void ofmo_set_new_scc_flag();
extern int ofmo_calc_fragment_electronic_state(
    MPI_Comm comm, int nmonomer, int monomer_list[], int level,
    double tolscf,
    double *energy, double *energy0, double *ddv,
    int *fnao, double daopop[], int sao2tao[],
    int *fnat, double datpop[], int fatom2tatom[]);

extern int ofmo_calc_es_dimer(MPI_Comm comm, int njob, int joblist[], double *energy0);

#endif

