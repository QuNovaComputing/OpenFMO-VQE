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

extern int ofmo_vqe_call( const int mythread, const int nmonomer, const int monomer_list[], const int nao, const double H[],
    const double U[], const double mo_tei[], const double S[], const double C[], const int nelec, const double Enuc,
    const double energy, const int iscc, const double ev[], char * desc, const int homo, const int lumo, const int ent);

extern int ofmo_vqe_get_amplitudes( const int ifrag, const int iscc, const int nso, int * namps, double ** amp, char *** fock, char * desc);

extern int ofmo_vqe_get_energy( const int nmonomer, const int monomer_list[], const int iscc, double * energy, double * dv, char * desc );

extern int ofmo_vqe_posthf_density( const int na, const double * A, char ** fock_vec,
    const double C[], const int nao, double D[]);
#endif
