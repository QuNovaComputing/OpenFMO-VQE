#ifndef _OFMO_MONOMER_DATA_H_
#define _OFMO_MONOMER_DATA_H_

#include <stdio.h>
extern void ofmo_reset_monomer_data_workerp();

extern int ofmo_get_monomer_density( MPI_Comm comm, const int ifrag,
	double D[] );
//extern int ofmo_get_monomer_density( int ifrag, double D[] );
extern int ofmo_put_monomer_density( int ifrag, double D[] );

extern int ofmo_get_monomer_aopop( int ifrag, double aopop[] );
extern int ofmo_put_monomer_aopop( int ifrag, double aopop[] );

extern int ofmo_get_monomer_atpop( int ifrag, double atpop[] );
extern int ofmo_put_monomer_atpop( int ifrag, double atpop[] );

extern int ofmo_update_monomer_data();

extern int ofmo_get_interfrag_distance( int ifrag, double dist[] );
extern int ofmo_put_interfrag_distance( int ifrag, double dist[] );

extern int ofmo_put_monomer_energy( int ifrag, double src[] );
extern int ofmo_get_monomer_energy( int ifrag, double dist[] );

extern int ofmo_put_monomer_energy0( int ifrag, double src[] );
extern int ofmo_get_monomer_energy0( int ifrag, double dist[] );

extern int ofmo_acc_dimer_aopop( int nao, double *aopop, int *iconv );
extern int ofmo_acc_dimer_atpop( int nat, double *atpop, int *iconv );

extern void ofmo_show_cache_prof();
#endif
