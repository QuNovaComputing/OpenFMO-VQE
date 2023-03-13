#ifndef _OFMO_PUT_GET_MASTER_H_
#define _OFMO_PUT_GET_MASTER_H_

#include <stdio.h>
#include <mpi.h>

extern int ofmo_master_get( MPI_Comm comm_mserv, int data_id,
	int ifrag, double D[] );
extern int ofmo_master_put( MPI_Comm comm_mserv, int data_id,
	int ifrag, double *src );
extern int ofmo_master_put_all( int nmservs, MPI_Comm comm_mservs[],
	int data_id, int ifrag, double *src );
 
extern int ofmo_get_monomer_density_master( MPI_Comm comm_mserv,
	int ifrag, double D[] );
extern int ofmo_get_monomer_aopop_master( MPI_Comm comm_mserv,
	int ifrag, double aopop[] );
extern int ofmo_get_monomer_atpop_master( MPI_Comm comm_mserv,
	int ifrag, double atpop[] );
extern int ofmo_update_monomer_data_master();
extern void ofmo_reset_monomer_data_master();

extern int ofmo_read_and_put_all( int nmservs, MPI_Comm comm_mservs[],
	char *fileheader );
extern int ofmo_read_and_put( MPI_Comm comm_mserv, char *fileheader );
extern int ofmo_get_and_write( MPI_Comm comm_mserv, char *fileheader );

extern int ofmo_copy_mserv_data( MPI_Comm comm_mserv_src,
	MPI_Comm comm_mserv_dst );

#endif
