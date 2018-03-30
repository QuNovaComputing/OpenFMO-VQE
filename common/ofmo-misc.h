#ifndef _OFMO_MISC_H_
#define _OFMO_MISC_H_

extern int ofmo_isum2( const int n, const int ia[] );
extern int ofmo_imax( const int n, const int ia[] );

extern int** ofmo_alloc_imatrixv( const int nrow, const int ncols[] );
extern long** ofmo_alloc_lmatrix( const int nrow, const int ncol );
extern double** ofmo_alloc_dmatrixv( const int nrow, const int ncols[] );
extern int** ofmo_alloc_imatrix( const int nrow, const int ncol );
extern double** ofmo_alloc_dmatrix( const int nrow, const int ncol );
extern char** ofmo_alloc_cmatrix( const int nrow, const int ncol );

extern void ofmo_free_imatrix( int **t );
extern void ofmo_free_lmatrix( long **t );
extern void ofmo_free_dmatrix( double **t );
extern void ofmo_free_cmatrix( char **t );

/* ------------------------------------- */
double ofmo_assignRes(int nres, int ntask, double w[]);
int ofmo_getMyColorFromAssignedRes(int me, int res[]);
int ofmo_getMyColorAndRankFromAssignedRes(int me, int res[], int *newrank);
int ofmo_aggregateTask(int nres, int ntask, double *weight, int *atask);

#endif
