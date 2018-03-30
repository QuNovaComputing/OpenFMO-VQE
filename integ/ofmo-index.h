#ifndef _OFMO_INDEX_2_H
#define _OFMO_INDEX_2_H

extern int ofmo_index_init( const int maxlqn );
extern int** ofmo_getadd_angm();
extern int* ofmo_getadd_laot();
extern int* ofmo_getadd_nnao();
extern int** ofmo_getadd_nam();
extern int** ofmo_getadd_nap();
extern int* ofmo_getadd_indx();
extern int** ofmo_getadd_nam2();
extern int** ofmo_getadd_nap2();
extern double* ofmo_getadd_dfact();

#endif
