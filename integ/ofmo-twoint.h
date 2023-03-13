#ifndef _OFMO_TWOINT_H_
#define _OFMO_TWOINT_H_

#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// to be called from non-thread-parallel region
extern int ofmo_twoint_init();

// to be called from thread-parallel region
extern size_t ofmo_twoint_set_buffer_size( const int mythread,
	const size_t ebuf_buffer_size_mb );

extern double* ofmo_twoint_alloc_local_gmat( const int mythread,
	const int maxnao );

extern double* ofmo_twoint_alloc_square_density( const int maxnao );

extern void ofmo_twoint_set_last_ijcs( const int mythread,
	const int last_ijcs );

extern void ofmo_twoint_set_last_klcs( const int mythread,
	const int last_klcs );

extern void ofmo_twoint_set_last_eri_type( const int mythread,
	const int last_eri_type );

extern void ofmo_twoint_set_stored_nzeri( const int mythread,
	const size_t ebuf_non_zero_eri );

extern int ofmo_twoint_get_last_eri_type( const int mythread );

extern int ofmo_twoint_get_last_ijcs( const int mythread );

extern int ofmo_twoint_get_last_klcs( const int mythread );

extern double* ofmo_twoint_get_ebuf_eri( const int mythread );

extern short int* ofmo_twoint_get_ebuf_ind4( const int mythread );

extern size_t ofmo_twoint_get_stored_nzeri( const int mythread );

// for debug
extern size_t ofmo_twoint_get_max_nzeri( const int mythread );

// added to use small buffer
extern size_t ofmo_twoint_get_max_stored_integ( const int mythread );
extern void ofmo_twoint_set_stored_integ( const int mythread, const size_t ninteg );
extern size_t ofmo_twoint_get_stored_integ( const int mythread );
extern double* ofmo_twoint_getadd_integ_val( const int mythread );
extern short* ofmo_twoint_getadd_integ_ind4( const int mythread );

//
extern void ofmo_twoint_set_global_last_eri_type( const int last_eri_type );
extern int ofmo_twoint_get_global_last_eri_type( void );

// Dcs
extern float* ofmo_twoint_get_Dcs( void );
extern float* ofmo_twoint_alloc_Dcs( const int maxncs );
extern void ofmo_twoint_free_Dcs( void );
extern float* ofmo_twoint_gen_Dcs( const int maxlqn, const int nao,
    const int leading_cs[], const double D[]);
extern float ofmo_twoint_dmax6(const int i, const int j,
    const int k, const int l);
extern float ofmo_twoint_dmax2(const int k, const int l);

// EPS
extern float ofmo_twoint_eps_eri(float eps);
extern float ofmo_twoint_eps_ps4(float eps);
extern float ofmo_twoint_eps_sch(float eps);

#endif
