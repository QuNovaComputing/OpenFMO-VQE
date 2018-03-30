#ifndef _OFMO_MSERV_CONT_H_
#define _OFMO_MSERV_CONT_H_

#include <stdio.h>
#include <mpi.h>

extern int ofmo_mserv_container_init(
	int def_mserv_id, int ndata, int nfrag,
	int** data_target, int nmservs, char* service_name_header );

extern int ofmo_connect_mserv( int mserv_id, MPI_Comm comm_worker );
extern int ofmo_disconnect_mserv( int mserv_id );
extern int ofmo_worker_put( int data_id, int ifrag, double *src );
extern int ofmo_worker_get( int data_id, int ifrag, double *dest );
extern int ofmo_worker_acc( int data_id, int nelem, double *data, int *iconv );


extern char* ofmo_getadd_service_name( int mserv_id );
extern char* ofmo_getadd_port_name( int mserv_id );
extern int ofmo_set_comm_mserv( int mserv_id, MPI_Comm comm_mserv );
extern int ofmo_set_default_mserv( int mserv_id );
extern int ofmo_worker_barrierp();
#endif
