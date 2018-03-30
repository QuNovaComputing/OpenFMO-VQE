#include <stdio.h>
#include <stdlib.h>

/** 整数並列要素の最大値を返す
 */
int ofmo_imax( const int n, const int ia[] ) {
    int i, imax;
    imax = ia[0];
    for ( i=1; i<n; i++ ) if ( ia[i] > imax ) imax = ia[i];
    return imax;
}

/** 配列の要素の和を求める
 * */
int ofmo_isum2( const int n, const int ia[] ) {
    int i, sum = 0;
    for ( i=0; i<n; i++ ) sum += ia[i];
    return sum;
}

/* 列サイズの異なる二次元整数配列を確保する */
int** ofmo_alloc_imatrixv( const int nrow, const int ncols[] ) {
    int total, i, **t;
    for ( i=0, total=0; i<nrow; i++ ) total += ncols[i];
    if ( total < 1 ) return NULL;
    t = (int**)malloc( sizeof(int*) * nrow );
    t[0] = (int*)malloc( sizeof(int) * total );
    for ( i=1; i<nrow; i++ ) t[i] = t[i-1] + ncols[i-1];
    return t;
}

/* 列サイズの異なる二次元実数配列を確保する */
double** ofmo_alloc_dmatrixv( const int nrow, const int ncols[] ) {
    int total, i;
    double **t;
    for ( i=0, total=0; i<nrow; i++ ) total += ncols[i];
    if ( total < 1 ) return NULL;
    t = (double**)malloc( sizeof(double*) * nrow );
    t[0] = (double*)malloc( sizeof(double) * total );
    for ( i=1; i<nrow; i++ ) t[i] = t[i-1] + ncols[i-1];
    return t;
}

/* 矩形の二次元整数配列(nrow x ncol)を確保する */
int** ofmo_alloc_imatrix( const int nrow, const int ncol ) {
    int i, **t;
    t    = (int**)malloc( sizeof(int*) * nrow );
    t[0] = (int*)malloc( sizeof(int) * nrow * ncol );
    for ( i=1; i<nrow; i++ ) t[i] = t[i-1] + ncol;
    return t;
}

/* 矩形の二次長整数配列(nrow x ncol)を確保する */
long** ofmo_alloc_lmatrix( const int nrow, const int ncol ) {
    int i;
    long **t;
    t    = (long**)malloc( sizeof(long*) * nrow );
    t[0] = (long*)malloc( sizeof(long) * nrow * ncol );
    for ( i=1; i<nrow; i++ ) t[i] = t[i-1] + ncol;
    return t;
}

/** 矩形の二次元実数配列(nrow x ncol)を確保する */
double** ofmo_alloc_dmatrix( const int nrow, const int ncol ) {
    int i;
    double **t;
    t    = (double**)malloc( sizeof(double*) * nrow );
    t[0] = (double*)malloc( sizeof(double) * nrow * ncol );
    for ( i=1; i<nrow; i++ ) t[i] = t[i-1] + ncol;
    return t;
}

/** 矩形の二次元文字配列(nrow x ncol)を確保する */
char** ofmo_alloc_cmatrix( const int nrow, const int ncol ) {
    int i;
    char **t;
    t  = (char**)malloc( sizeof(char*) * nrow );
    t[0] = (char*)malloc( sizeof(char) * nrow * ncol );
    for ( i=1; i<nrow; i++ ) t[i] = t[i-1] + ncol;
    return t;
}

/** 二次元整数配列を開放する */
void ofmo_free_imatrix( int **t ) {
    if ( t != NULL ) {
	if ( t[0] != NULL ) free( t[0] );
	free( t );
    }
}

/** 二次元実数配列を開放する */
void ofmo_free_dmatrix( double **t ) {
    if ( t != NULL ) {
	if ( t[0] != NULL ) free( t[0] );
	free( t );
    }
}

/** 二次元文字配列を開放する */
void ofmo_free_cmatrix( char **t ) {
    if ( t != NULL ) {
	if ( t[0] != NULL ) free( t[0] );
	free( t );
    }
}

/** 二次元長整数配列を開放する */
void ofmo_free_lmatrix( long **t ) {
    if ( t != NULL ) {
	if ( t[0] != NULL ) free( t[0] );
	free( t );
    }
}

/* ------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#ifndef MAX2
#define MAX2(a, b) (((a)>(b))? (a): (b))
#endif

// Assign resources (such as procs) to tasks by their weight
// nres: # resources (IN)
// ntask: # tasks (IN)
// w[]: task weight (IN), # assigned resources (OUT)
// (!) Must be nres>=ntask
double ofmo_assignRes(int nres, int ntask, double w[])
{
  double r[ntask];
  int nr[ntask];
  double rr=1.0/nres;
  int i;
  assert(nres>=ntask);

  double tw=0;
  for (i=0; i<ntask; i++ ) tw += w[i];
  for (i=0; i<ntask; i++ ) r[i]=w[i]/tw;

  int nnr=0;
  //for ( i=0; i<ntask; i++ ) nnr+=(nr[i]=1);
  for ( i=0; i<ntask; i++ ) nnr+=(nr[i]=(int)round(MAX2(1,(nres-ntask)*r[i])));
  //for ( i=0; i<ntask; i++ ) nnr+=(nr[i]=(int)(MAX2(1,(nres-ntask)*r[i])));
//  printf("%d for %d\n",nnr, nres);

  while (nnr!=nres) {
    if (nnr>nres) {
      int m=-1;
      double rm=1.0;
      for ( i=0; i<ntask; i++ ) {
        double rn=r[i]/nr[i];
        if (rn<rm && nr[i]>1) { m=i; rm=rn;}
      }
      nr[m]--;
      nnr--;
    } else {
      int m=-1;
      double rm=0.0;
      for ( i=0; i<ntask; i++ ) {
        double rn=r[i]/nr[i];
        if (rn>rm) { m=i; rm=rn;}
      }
      nr[m]++;
      nnr++;
    }
  }
//  for ( i=0; i<ntask; i++ ) printf("%4d [%.3f]: %2d (%.3f) %.3f\n",(int)w[i],r[i],nr[i],r[i]/nr[i], rr/(r[i]/nr[i]));

  double rm=0.0;
  for ( i=0; i<ntask; i++ ) {
    rm=MAX2(rm, r[i]/nr[i]);
    w[i] = nr[i];
  }

//  printf("%.3f / %.3f = %.3f\n",rr, rm, rr/rm);
  return rr/rm;
}

// Get my task number from assigned resouces
// me: my rank (IN)
// res[]: # assigned resources for task (IN)
// (!) res[] is int array, not as w[] assginRes() output
int ofmo_getMyColorFromAssignedRes(int me, int res[])
{
  int color=0, n=0;
  while (n+res[color]<=me) n+=res[color++];

  return color;
}

// Get my task number & new rank from assigned resouces
// me: my rank (IN)
// res[]: # assigned resources for task (IN)
// newrank: new rank (OUT)
// (!) res[] is int array, not as w[] assginRes() output
int ofmo_getMyColorAndRankFromAssignedRes(int me, int res[], int *newrank)
{
  int color=0, n=0;
  while (n+res[color]<=me) n+=res[color++];
  *newrank=me-n;

  return color;
}

static double *wcmpA;
static int wcmp(const void *a, const void *b)
{
  return wcmpA[*(int*)a]<wcmpA[*(int*)b];
}

// Aggregate small tasks
// nres: # resources (IN)
// ntask: # tasks (IN)
// weight[]: task weight (IN), task weight after aggregation (OUT)
// atask[]: conv. table from old task to new task (OUT)
// return value: # tasks after aggregation
int ofmo_aggregateTask(int nres, int ntask, double *weight, int *atask)
{
  int i,j;

  double *w;
  int *gidx,*idx,*ridx;
  w=(double*)malloc(ntask*sizeof(double)+ntask*sizeof(int));
  gidx=(int*)(w+ntask);
  idx=atask;
  ridx=(int*)w;

  int tw=0;
  for (i=0;i<ntask;i++) tw+=weight[i];
  double aw=tw/nres;

  double th;
  int ith;
  int iw=-1;
  int imw=-1;
  int mw=tw;
  int nt=ntask;
//  for (ith=9;ith<=12;ith++) {
  for (ith=20;ith<=22;ith++) {
    th=ith/10.0;

    for (i=0;i<ntask;i++) w[i]=weight[i];

    for (i=0;i<ntask;i++) gidx[i]=idx[i]=i;

    wcmpA=w;
    qsort(idx,ntask,sizeof(int),wcmp);

    iw=-1;
    imw=-1;
    mw=tw;
    nt=ntask;
    for (j=0;j<ntask;j++) {
      int i=idx[j];
      if(w[i]<aw*th) {
        if (iw<0) iw=i;
        else {
          w[iw]+=w[i];
          gidx[i]=iw;
          w[i]=0;
          nt--;
          if(w[iw]>aw*th) {
            if(w[iw]<mw) {imw=iw; mw=w[iw];}
            iw=-1;
          }
        }
      }
      if(w[i]>aw*th && w[i]<mw) {imw=i; mw=w[i];}
    }
    //if (iw>=0 && w[iw]<aw*0.9) {
    if (iw>=0 && w[iw]<aw*0.8) {
      w[imw]+=w[iw];
      w[iw]=0;
      nt--;
      gidx[iw]=imw;
      for (j=0;j<ntask;j++)
        if (gidx[j]==iw) gidx[j]=imw;
    }
    if (nt<=nres) break;
  }

  for (i=0;i<ntask;i++) idx[i]=i;
  wcmpA=w;
  qsort(idx,ntask,sizeof(int),wcmp);
  for (i=0;i<ntask;i++) weight[i]=w[idx[i]];
  for (i=0;i<nt;i++) ridx[idx[i]]=i;
  for (i=0;i<ntask;i++) atask[i]=ridx[gidx[i]];
  free(w);
  return nt;
}

#if 0
int main(int argc, char *argv[])
{
  int nproc=8;
  int nfcs[]={1000,440, 3200,70,1200,400,320, 100, 150};
//  int nfcs[]={1000,400,500,120,400,300,70};
  int nifc4c;
  int i,j;

  if (argc>1) nproc = atoi(argv[1]);


  nifc4c=sizeof(nfcs)/sizeof(int);
  double *w=(double*)malloc(nifc4c*sizeof(double));
  for (i=0;i<nifc4c;i++) w[i]=nfcs[i];

  int *atask=(int*)malloc(nifc4c*sizeof(int));
  int natask;
  natask=ofmo_aggregateTask(nproc, nifc4c, w, atask);
//  for (i=0;i<nifc4c;i++) printf("(%2d) %4d : %2d\n",i, nfcs[i], atask[i]);


  ofmo_assignRes(nproc, natask, w);

  int *res=(int*)malloc(natask*sizeof(int));
  for (i=0;i<natask;i++) res[i]=(int)w[i];

  for (i=0;i<natask;i++) {
//    printf("(%2d) ", i);
//    for (j=0;j<nifc4c;j++) if (i==atask[j]) printf("%2d ", j);
    printf("(%2d) %2d: ", i, res[i]);
    for (j=0;j<nifc4c;j++) if (i==atask[j]) printf("%4d ", nfcs[j]);
    printf("\n");
  }

  for (i=0; i<nproc; i++) {
    int nrank;
    int color=ofmo_getMyColorFromAssignedRes(i,natask,res);
    nrank=res[color];
//    printf("(%2d) %2d: %1d/%1d\n", i, color,nrank);
    printf("%2d ", color);
  }
  printf("\n");

  return 0;
}
#endif
/* ------------------------------------- */
