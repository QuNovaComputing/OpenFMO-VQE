/**
 * @file ofmo-mat.c
 * BLAS、LAPACKの関数を呼び出して、行列演算を行う関数群。
 *
 * BLAS、LAPACKの関数呼び出しを、このファイルで定義された関数だけで
 * 行いたいために、用意している。
 * ほぼ、すべてのOpenFMOのクラスから参照される。
 *
 * @attention
 * - 対称行列に対する操作を行う際には、ほぼすべての場合で、U形式での
 *   操作を行っている。つまり、BLAS、LAPACKの関数呼び出し時に、
 *   \c UPLO を指定する必要がある場合、この関数内部では、ほぼ常に、
 *   "U"を渡して演算を行っている。
 *
 * */
/**
 * @defgroup ofmo-mat 行列演算関連の関数群
 * 
 * @ingroup ofmo
 * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ofmo-def.h"

#define FORT_US1

#ifdef FORT_US1
#define dcopy dcopy_
#define ddot ddot_
#define daxpy daxpy_
#define dscal dscal_
#define dgemm dgemm_
#define dsymm dsymm_
#define dtrsv dtrsv_
#define dtrsm dtrsm_
#define dgesvd dgesvd_

#define dpotrf dpotrf_
#define ilaenv ilaenv_
#define dsyev dsyev_
#define dsygst dsygst_
#define dgetrs dgetrs_
#define dgetrf dgetrf_

#define dsysv dsysv_
#define dgesv dgesv_
#endif

#ifdef FORT_DOTUPPER
#define dcopy .DCOPY
#define ddot .DDOT
#define daxpy .DAXPY
#define dscal .DSCAL
#define dgemm .DGEMM
#define dsymm .DSYMM
#define dtrsv .DTRSV
#define dtrsm .DTRSM

#define dpotrf .DPOTRF
#define ilaenv .ILAENV
#define dsyev .DSYEV
#define dsygst .DSYGST
#define dgetrs .DGETRS
#define dgetrf .DGETRF

#define dsysv .DSYSV
#define dgesv .DGESV
#endif

extern double ddot(int* N, const double x[], int* incx,
	const double y[], int* incy);
extern void daxpy(const int* N, const double* A, const double x[], const int* incx,
	double y[], const int* incy);
extern void dscal(const int* N, const double* A, double x[], const int* incx);
extern void dgemm ( const char *transa, const char *transb,
	const int *m, const int *n, const int *k, const double *alpha,
	const double A[], const int *lda, const double B[],
	const int *ldb, const double *beta, double C[], const int *ldc );

extern void dsygst(const int* ITYPE, const char* UPLO, const int* N,
	double* A, const int* LDA, const double* B, const int* LDB,
	int* INFO);
extern void dpotrf(char* UPLO, int* N, double* A, int* LDA, int* INFO);
extern void dsyev(char* JOBZ, char* UPLO, int* N, double* A,
	int* LDA, double* W, double* WORK, int* LWORK, int* INFO);
extern void dtrsm(char* SIDE, char* UPLO, char* TRANSA, char* DIAG,
	int* M, int* N, double* alpha, const double A[], int* LDA,
	double B[], int* LDB );
extern void dgesv( const int* N, const int* HRHS, double A[],
	const int* LDA, int ipiv[],
	double B[], const int* LDB, int* INFO );

/** Substitute the compressed matrix (U format) into a square matrix
 *
 * Substitute (expand) a symmetric matrix saved in compressed U format
 * into a square matrix as a perfect symmetric matrix
 *
 * @param[in] n The size of the symmetric matrix to be manipulated
 * @param[in] U[n*(n+1)/2] An array containing the elements of a
 *    compressed U-style symmetric matrix
 * @param[out] S[n*n] An array to which the elements of the expanded square matrix are assigned
 *
 * */
void ofmo_unpack_matrix(const int n, const double U[], double S[]) {
    int i, j, ij, I, J;
    ij = 0;
    for (i=0, I=0; i<n; i++, I+=n) {
	for (j=0, J=0; j<=i; j++, J+=n) {
	    S[I+j] = S[J+i] = U[ij];
	    ij++;
	}
    }
}

/**  Substitute the U part of a square matrix into a compressed matrix
 *
 * Substitute (compress) the upper right triangular matrix (U) of the square
 * matrix into the compressed U-format symmetric matrix.
 *
 * @param[in] n The size of the symmetric matrix to be manipulated
 * @param[in] S[n*n] An array to which the elements of a square matrix symmetric
 *     matrix are assigned. Only the upper right triangular part including the diagonal element is referenced.
 * @param[out] U[n*(n+1)/2] An array to which the compressed symmetric matrix elements are assigned
 *
 * @ingroup ofmo-mat
 * * */
void ofmo_pack_matrix(const int n, const double S[], double U[]) {
    int i, j, ij, J, IJ;
    ij = 0;
    for (j=0, J=0; j<n; j++, J+=n) {
	for (i=0, IJ=J; i<=j; i++, IJ++) {
	    U[ij] = S[IJ];
	    ij++;
	}
    }
}

/** Fold a square matrix into a compressed matrix (U format)
 *
 * 正方行列を圧縮行列に折りたたむ。対角成分は2倍になる。
 *
 * つまり、元になる正方行列を\f$ S \f$、圧縮後の行列を\f$ U \f$とした
 * 場合、\f$ U_{ij} = S_{ij} + S_{ji},\ \ \ (0\le i\le j < \tt{n}) \f$を
 * 計算する
 *
 * @param[in] n 正方行列のサイズ
 * @param[in] S[n*n] 正方行列の要素が代入されている配列
 * @param[out] U[n*(n+1)/2] 折りたたみ後の圧縮U形式の対称行列要素が
 *     代入されている配列
 *
 * @ingroup ofmo-mat
 * */
void ofmo_fold_matrix( const int n, const double S[], double U[] ) {
    int i, j, IJ;
    IJ = 0;
    for ( i=0; i<n; i++ ) {
	for ( j=0; j<=i; j++ ) U[IJ] = S[i*n+j] + S[j*n+i];
    }
}

/** Copy the upper right triangular matrix (U) part of the square matrix to the lower left triangular part
 *
 * 正方行列の右上三角行列(U)部分を左下三角部分にコピーする。
 * LAPACKの行列演算関数の中には、正方行列に対する対称変換
 * などの演算を、右上（あるいは、左下）三角行列部分にのみ
 * 適用するものがある（dsygstなど）。そのような関数を呼び出した後に、
 * 完全な正方行列の形にするために、この関数を準備した。
 *
 * @param[in] n 正方行列のサイズ
 * @param[in] A[n*n] （入力時）対角成分を含む右上三角行列部分にのみ
 *     正しい値が代入されている正方行列。\n
 *     （出力時）対角要素の除く左下三角部分にも、対応する右上参加九部分の
 *     値が代入された正方行列形式の完全な対称行列
 *
 * @ingroup ofmo-mat
 * */
void ofmo_expand_U2L( const int n, double A[] ) {
    int i,j, jn, ij, ji;
    for (j=1, jn=n; j<n; j++, jn+=n) {
	for (i=0, ij=jn, ji=j; i<j; i++, ij++, ji+=n) {
	    //ij = i + j*n; ji = j + i*n;
	    A[ji] = A[ij];
	}
    }
}

/** A function that returns the maximum absolute value of the difference between two real vectors
 *
 * ２つの実数ベクトルの差の最大絶対値を返す関数。
 * 
 * @param[in] n ベクトルのサイズ
 * @param[in] v1[n] １つ目のベクトル \f$ \bold{v}_1 \f$
 * @param[in] v2[n] ２つ目のベクトル \f$ \bold{v}_2 \f$
 *
 * @return \f$ \max( | v_1(i)-v_2(i) |),\ \ \ (0\le i<\tt{n})\f$
 *
 * @ingroup ofmo-mat
 * */
double ofmo_max_diff( const int n, const double v1[], const double v2[] ) {
    int i;
    double d, max_diff = 0.e0;
    for ( i=0; i<n; i++ ) {
	d = fabs( v1[i] - v2[i] );
	if ( d > max_diff ) max_diff = d;
    }
    return max_diff;
}

/** Transpose a square matrix
 *
 * 正方行列の転置を行う。
 *
 * @param[in] n 正方行列のサイズ
 * @param[in] S[n*n] 転置前の正方行列
 * @param[out] St[n*n] 転置後の正方行列
 *
 * @ingroup ofmo-mat
 * */
void ofmo_transpose_matrix(const int n, const double S[], double St[]) {
    int i, j, in, jn, ij;
    for (i=0, in=0; i<n; i++, in+=n) {
	for (j=0, jn=0, ij=in; j<n; j++, jn+=n, ij++) {
	    St[j*n + i] = S[i*n + j];
	}
    }
}

/** Multiply the diagonal component of a compressed U-format matrix by a constant
 *
 * 圧縮U形式で１次元配列に格納された倍精度浮動小数対称行列の対角要素を
 * 定数倍する関数
 *
 * @param[in] n 正方行列のサイズ
 * @param[in] alpha 対角要素のかける定数
 * @param[in,out] X[] 操作対象となる対称行列（圧縮U形式）の行列要素が
 *     代入されている１次元配列
 *
 * @ingroup ofmo-mat
 * */
void ofmo_scale_diag(const int n, const double alpha, double X[]) {
    int i,j;
    i = -1;
    for (j=1; j<=n; j++) {
        i += j;
        X[i] *= alpha;
    }
}

/** Copy of real vector
 *
 * 実数ベクトルのコピー
 *
 * @param[in] n ベクトルのサイズ
 * @param[in] src[n] コピー元のベクトル要素が代入された配列
 * @param[out] dest[n] コピー先のベクトル要素が代入された配列
 *
 * @ingroup ofmo-mat
 * */
void ofmo_dcopy( const int n, const double src[], double dest[] ) {
    memcpy( dest, src, sizeof(double)*n );
}

/** Multiply each element of a real vector by a constant
 *
 * 実数ベクトル各要素を定数倍する
 * \f$ \bold{x} \leftarrow \alpha\bold{x}\f$
 *
 * @param[in] n ベクトルサイズ
 * @param[in] alpha 各要素にかける定数 \f$ \alpha \f$
 * @param[in,out] x[n] 操作対象となるベクトルの要素が代入された配列
 *
 * @ingroup ofmo-mat
 * */
void ofmo_dscale( const int n, const double alpha, double x[] ) {
    for ( int i=0; i<n; i++ ) x[i] *= alpha;
  //int incx = 1;
  //dscal( &n, &alpha, x, &incx);
}

/** daxpy calculation
 *
 * vec(y) <- alpha * vec(x) + vec(y)
 * \f$ \bold{y} \leftarrow \alpha\bold{x} + \bold{y} \f$の計算を行う
 *
 * @param[in] n ベクトルサイズ
 * @param[in] alpha 定数\f$ \alpha \f$
 * @param[in] x[n] ベクトル\f$ \bold{x}\f$の要素が代入
 *     された配列
 * @param[in,out] y[n] ベクトル\f$\bold{y}\f$の要素が代入
 *     された配列。出力時には各要素に\f$\alpha\bold{x}\f$
 *     が加算されている。
 *
 * @ingroup ofmo-mat
 * */
void ofmo_daxpy( const int n, const double alpha, const double x[],
	double y[] ) {
  //  for ( int i=0; i<n; i++ ) y[i] += alpha*x[i];
  int incx = 1;
  int incy = 1;
  daxpy( &n, &alpha, x, &incx, y, &incy);
}

/** Take the sum of integer arrays
 *
 * 整数配列の和をとる関数
 *
 * @ingroup ofmo-mat
 * */
int ofmo_isum( const int n, const int iv[] ) {
    int sum = 0, i;
    for ( i=0; i<n; i++ ) sum += iv[i];
    return sum;
}

/** Matrix product
 *
 * Perform DGEMM calculations on square matrices of the same size
 * C <- alpha * A * B + C
 * \f$ \bold{C} \leftarrow \alpha\bold{AB}+\beta\bold{C} \f$
 *
 * @ingroup ofmo-mat
 * */
void ofmo_dgemm( const int n, const char *transa, const char *transb,
	const double alpha, const double A[], const double B[],
	const double beta, double C[] ) {
    dgemm( transa, transb, &n, &n, &n, &alpha, A, &n, B, &n, &beta, C, &n );
}

/** Matrix transformation from the generalized eigenvalue problem to the standard eigenvalue problem
 *
 * Only for compressed U format
 * @ingroup ofmo-mat
 * */
int ofmo_dsygst( const int itype, const int n, double A[],
	const double U[] ) {
    int info;
    dsygst( &itype, "U", &n, A, &n, U, &n, &info );
    return info;
}

static int RESERVED = 0;
static int LWORK    = 0;
//static int ILWORK   = 0;
static int *IPIV    = NULL;
static double *WORK = NULL;

static void dealloc_work_area() {
    if ( IPIV != NULL ) free( IPIV );
    if ( WORK != NULL ) free( WORK );
    IPIV = 0;
    WORK = 0;
    LWORK = 0;
    //ILWORK = 0;
    RESERVED = 0;
}

static int alloc_work_area( const int n ) {
    int nb = 64;
    dealloc_work_area();
    IPIV = (int*)malloc( sizeof(int) * n );
    WORK = (double*)malloc( sizeof(double) * (nb+2) * n );
    if ( IPIV == NULL || WORK == NULL ) {
	dealloc_work_area();
	return -1;
    }
    LWORK = (nb+2) * n;
    //ILWORK = n;
    RESERVED = n;
    return 0;
}

/** Initialize this function group
 *
 * Secure work areas used for eigenvalue calculation, etc.
 *
 * @param[in] n Number of dimensions of the matrix to handle
 *
 * @note If you call this function with the maximum dimension
 * of the matrix to be calculated as an argument at the first call,
 * you can avoid calling unnecessary \c free and \c malloc .
 *
 * @ingroup ofmo-mat
 * */
int ofmo_mat_init( const int n ) {
    static int called = false;
    if ( RESERVED >= n ) return 0;
    if ( alloc_work_area( n ) != 0 ) return -1;
    if ( called == false ) {
	atexit( dealloc_work_area );
	called = true;
    }
    return 0;
}

/** Take the inner product of two real vectors
 *
 * ２つの倍精度浮動小数ベクトルの内積をとる関数。
 *
 * @param[in] n ベクトルの次元
 * @param[in] x[] １つ目のベクトル要素が代入された配列
 * @param[in] y[] ２つ目のベクトル要素が代入された配列
 *
 * @return 内積の値
 *
 * @ingroup ofmo-mat
 * */
double ofmo_dot_product( const int n, const double x[], const double y[] ) {
    int incx=1, N=n;
    return ddot( &N, x, &incx, y, &incx );
}

/** Functions that perform Cholesky decomposition
 *
 * Perform Cholesky decomposition on a positive-definite matrix (square matrix format).
 *
 * @attention
 * @li 内部でLapackの関数 \c DPOTRF を呼んでいる
 *
 * @param[in] n Matrix size
 * @param[in,out] S[] (At the time of input) Elements of positive-definite symmetric matrix
 *     (At the time of output) Cholesky Matrix element of the decomposed result.
 *     The array size must be at least n^2 (\ f $ \ tt {n} ^ 2 \ f $).
 *
 * @retval 0 正常終了（ちゃんとCholesky分解できた）
 * @retval 0以外 異常終了（途中で固有値が負になって計算が出来なくなった）
 *
 * @ingroup ofmo-mat
 * */
int ofmo_chodec( const int n, double S[] ) {
    int info, N=n;
    dpotrf( "U", &N, S, &N, &info );
    return info;
}

/** Solve the generalized symmetric eigenvalue problem
 * Ax = lambda Bx
 * @f$ \bold{Ax}=\lambda\bold{Bx}@f$
 * Solve the type generalization eigenvalue problem.
 * The coefficient matrix @f $ \ bold {A} @f $ is a symmetric matrix,
 * and the matrix @f $ \ bold {B} @f $ is a definite matrix.
 * Find all eigenvalues ​​and eigenvectors.
 *
 * @param[in] n Number of dimensions of the matrix to solve
 * @param[in] U[] An array to which the result of Cholesky decomposition of
 *    the positive-definite symmetric matrix @f $ \ bold {B} @f $ is assigned
 *    (square matrix format).
 *    This data was obtained by calling \c ofmo_chodec .
 * @param[in,out] A[] (At the time of input) An array to which the coefficient matrix
 *    (@f $ \ bold {A} @f $, symmetric matrix) is assigned.
 *    (At the time of output) An array to which the eigenvector is assigned. Square matrix format.
 * @param[out] e[] Array to which eigenvalues ​​are assigned
 *
 * @ingroup ofmo-mat
 * */
int ofmo_solv_GSEP( const int n, const double U[],
	double A[], double e[] ) {
    int info, itype=1, N=n, lwork;
    double alpha=1.e0;
    ofmo_mat_init( n );
    lwork = LWORK;
    dsygst( &itype, "U", &N, A, &N, U, &N, &info );
    if ( info != 0 ) {
	printf("error in dsygst\n");
	return -1;
    }
    dsyev( "V", "U", &N, A, &N, e, WORK, &lwork, &info );
    if ( info != 0 ) {
	printf("error in dsyev\n");
	return -1;
    }
    dtrsm("L", "U", "N", "N", &N, &N, &alpha, U, &N, A, &N );
    return 0;
}

/** Solve the standard eigenvalue problem
 *
 * @f$\bold{Ax}=\lambda\bold{x}@f$タイプの標準固有値問題を解く。
 * ただし、係数行列@f$\bold{A}@f$は対称行列である。
 *
 * @param[in] n 行列のサイズ
 * @param[in,out] A[] （入力時）係数行列。（出力時）固有ベクトル。
 * 正方行列形式。
 * @param[out] e[] 固有ベクトルが代入される配列
 *
 * @ingroup ofmo-mat
 * */
int ofmo_solv_SEP( const int n, double A[], double e[] ) {
    int info, lwork, N=n;
    ofmo_mat_init( n );
    lwork = LWORK;
    dsyev( "V", "U", &N, A, &N, e, WORK, &lwork, &info );
    if ( info != 0 ) {
	printf("error in dsyev\n");
	return -1;
    }
    return 0;
}

/** Solve simultaneous linear equations
 *
 * 連立一次方程式 \f$\bold{Ax}=\bold{b}\f$ を解く
 *
 * @param[in] n 行列のサイズ
 * @param[in,out] A[] 係数行列（正方行列形式）
 * @param[out] b[] （入力時）ベクトル@f$\bold{b}@f$。
 * （出力時）解ベクトル@f$\bold{x}@f$
 *
 * @ingroup ofmo-mat
 * */
int ofmo_solv_leq( const int n, double A[], double b[] ) {
    int info, m=1;
    ofmo_mat_init( n );
    dgesv( &n, &m, A, &n, IPIV, b, &n, &info );
    return info;
}

int ofmo_svd( const int n, double A[], double S[], double U[], double VT[]){
    int info, lwork;
    int ret_n=n, ret_m=n, ret_lda=n, ret_ldu=n, ldvt=n;
    ofmo_mat_init( n );
    lwork = LWORK;
    dgesvd("All", "All", &ret_n, &ret_m, A, &ret_lda, S, U, &ret_ldu, VT,
           &ldvt, WORK, &lwork, &info);
    return info;
}
