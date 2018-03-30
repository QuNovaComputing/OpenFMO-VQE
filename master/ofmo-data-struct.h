#ifndef _OFMO_DATA_STRUCT_H_
#define _OFMO_DATA_STRUCT_H_
#include <stdio.h>

extern int* ofmo_getadd_data_target( int data_id );
extern int* ofmo_getadd_data_nelems( int data_id );
extern int* ofmo_getadd_data_offset( int data_id );

extern int ofmo_get_data_target( int data_id, int ifrag );
extern int ofmo_get_data_nelems( int data_id, int ifrag );
extern int ofmo_get_data_offset( int data_id, int ifrag );
extern int ofmo_get_nfrag();
extern int ofmo_get_ndata();

extern int ofmo_make_data_struct( int ndata, int nfrag, int mserv_size );

#endif
