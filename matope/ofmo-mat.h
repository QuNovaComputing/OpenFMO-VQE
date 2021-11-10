/**
 * @file ofmo-mat.h
 *
 * Header file that declares various matrix operations and function 
 * prototypes for matrix operations used in quantum chemistry calculations
 *
 * */
#ifndef _OFMO_MAT_H_
#define _OFMO_MAT_H_

#include <stdio.h>
#include "ofmo-def.h"

extern int ofmo_mat_init( const int n );

/* BLAS, LAPACKともに使用しない関数群 */
extern void ofmo_unpack_matrix(const int n, const double U[], double S[]);
extern void ofmo_pack_matrix(const int n, const double S[], double U[]);
extern void ofmo_fold_matrix( const int n, const double S[], double U[] );
extern void ofmo_expand_U2L( const int n, double A[] );
extern void ofmo_transpose_matrix(const int n, const double S[],
	double St[]);
extern void ofmo_scale_diag(const int n, const double alpha, double X[]);
extern double ofmo_max_diff( const int n, const double v1[],
	const double v2[] );
extern void ofmo_dcopy( const int n, const double src[], double dest[] );
extern void ofmo_dscale( const int n, const double alpha, double x[] );
extern int ofmo_isum( const int n, const int iv[] );

/* BLASの関数を利用する関数群 */
extern double ofmo_dot_product( const int n, const double x[],
	const double y[] );
extern void ofmo_daxpy( const int n, const double alpha,
	const double x[], double y[] );
extern void ofmo_dgemm( const int n, const char *transa, const char *transb,
	const double alpha, const double A[], const double B[],
	const double beta, double C[] );
extern int ofmo_dsygst( const int itype, const int n, double A[],
	const double U[] );

/* LAPACKの関数を利用する関数群 */
extern int ofmo_chodec( const int n, double S[] );
extern int ofmo_solv_GSEP( const int n, const double U[],
	double A[], double e[] );
extern int ofmo_solv_SEP( const int n, double A[], double e[] );

extern int ofmo_solv_leq( const int n, double A[], double b[] );

#endif
