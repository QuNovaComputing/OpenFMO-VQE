/**
 * @file ofmo-vqe.h
 *
 * Header file that declares prototypes of functions
 * related to VQE call
 * 
 * */
#ifndef _OFMO_VQE_H_
#define _OFMO_VQE_H_

#include <stdlib.h>

extern int ofmo_amp_alloc( const int nfrag );

extern int ofmo_amp_dealloc();

extern int ofmo_vqe_call( const int mythread, const int nmonomer, const int monomer_list[], const int nao, const double H[],
    const double mo_tei[], const double S[], const double C[], const int nelec, const double Enuc,
    double *energy, const int iscc, const double ev[]);

extern int ofmo_get_amps( const int ifrag, int *namp, double **alpha, int **fock_vec);

extern int ofmo_get_oldamps( const int ifrag, int *namp, double **old_alpha, int **old_fock_vec);

extern int ofmo_update_amps();


#endif
