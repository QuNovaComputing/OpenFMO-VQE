/* 積分計算に必要なindexの作成
 * 2次元配列を使用するバージョン
 * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>

#ifndef false
#define false 0
#endif
#ifndef true
#define true 1
#endif

static int** alloc_imatrix( const int m, const int n ) {
    int **im, i;
    im = (int**)malloc( sizeof(int*) * m );
    im[0] = (int*)malloc( sizeof(int) * m * n );
    for ( i=1; i<m; i++ ) im[i] = im[i-1] + n;
    return im;
}

static void free_imatrix( int **im ) {
    if ( im ) {
	if ( im[0] ) free( im[0] );
	free( im );
    }
}

// もっとも基本的な配列変数
static int **angm = NULL;
static int *laot = NULL;
static int *nnao = NULL;
static int MAXLQN = -1;

// 派生データ(1)
static int **nam = NULL;	// a-1i
static int **nap = NULL;	// a+1i
static int *indx = NULL;	// i

// 派生データ(2)
static int **nam2 = NULL;
static int **nap2 = NULL;

// 派生データ(3)
static double *dfact = NULL;	//

static void ofmo_index_finalize() {
    if ( laot ) free( laot );
    if ( nnao ) free( nnao );
    free_imatrix( angm );
    free_imatrix( nam );
    free_imatrix( nap );
    free_imatrix( nam2 );
    free_imatrix( nap2 );
    if ( indx ) free( indx );
    if ( dfact ) free( dfact );
    angm = NULL;
    laot = NULL;
    nnao = NULL;
    MAXLQN = -1;
    nam = NULL;
    nap = NULL;
    indx = NULL;
    nam2 = NULL;
    nap2 = NULL;
    dfact = NULL;
}

// 以下の値を代入する関数
// nnao[cs_type] = 各CSに含まれるAO数
// laot[cs_type] = 各CSタイプの先頭AOタイプ番号
// angm[ao_type][i] = 各AOタイプの軌道量子指数の各成分
// 
// ------------------------------------------------
// CSタイプ     S P D F G H I
// CSタイプ番号 0 1 2 3 4 5 6
//
// AOタイプ     s px py pz dxx dyy dzz dxy dxz dyz
// AOタイプ番号 0 1  2  3   4   5   6   7   8   9
//
// AOタイプ     fxxx fyyy fzzz fxxy fxxz fxyy fxzz fxyz fyyz fzzy
// AOタイプ番号  10   11   12   13   14   15   16   17   18   19
//
// 配列は、この関数内部で領域を確保する（後で解放すること）
static int assign_angular_momentum_index( const int maxlqn ) {
    int naot, lqn, iaot, ix, iy, iz;
    nnao = (int*)malloc( sizeof(int) * (maxlqn+1) );
    laot = (int*)malloc( sizeof(int) * (maxlqn+1+1) );
    naot = 0;
    for ( lqn=0; lqn<=maxlqn; lqn++ ) {
	laot[lqn] = naot;
	nnao[lqn] = (lqn+1)*(lqn+2)/2;
	naot += nnao[lqn];
    }
    laot[maxlqn+1] = naot;
    // laot, angmの代入
    angm = alloc_imatrix( naot, 3 );
    iaot = 0;
    // sタイプ
    laot[0] = iaot;
    angm[iaot][0] = angm[iaot][1] = angm[iaot][2] = 0;
    iaot++;
    // pタイプ以降
    for ( lqn=1; lqn<=maxlqn; lqn++ ) {
	laot[lqn] = iaot;
	for ( ix=lqn; ix>=0; ix-- ) {
	    for ( iy=lqn-ix; iy>=0; iy-- ) {
		iz = lqn - (ix + iy);
		if ( ix == lqn ) {		// ix=lqn, iy=iz=0
		    angm[iaot][0] = ix;
		    angm[iaot][1] = 0;
		    angm[iaot][2] = 0;
		    iaot++;
		    angm[iaot][0] = 0;
		    angm[iaot][1] = ix;
		    angm[iaot][2] = 0;
		    iaot++;
		    angm[iaot][0] = 0;
		    angm[iaot][1] = 0;
		    angm[iaot][2] = ix;
		    iaot++;
		} else if ( iy == lqn ) {
		    continue;
		} else if ( iz > iy ) {
		    continue;
		} else if ( iy == iz ) {	// iy=iz
		    angm[iaot][0] = ix;
		    angm[iaot][1] = iy;
		    angm[iaot][2] = iz;
		    iaot++;
		} else {			// lqn > iy > iz
		    angm[iaot][0] = ix;
		    angm[iaot][1] = iy;
		    angm[iaot][2] = iz;
		    iaot++;
		    angm[iaot][0] = ix;
		    angm[iaot][1] = iz;
		    angm[iaot][2] = iy;
		    iaot++;
		}
	    }	// for ( iy )
	}	// for ( ix )
    }	// for ( lqn )
    return 0;
}



// 与えられた軌道量子指数を持つAOタイプを返す
static int get_aotype( const int nx, const int ny, const int nz ) {
    int istart, iend, i, lambda;
    lambda = nx+ny+nz;
    istart = laot[lambda];
    iend   = istart + nnao[lambda];
    for ( i=istart; i<iend; i++ ) {
	if ( angm[i][0] == nx && angm[i][1] == ny && angm[i][2] == nz )
	    return i;
    }
    return -1;
}

// 三要素のうち、非ゼロ要素で最も小さな値を持つ位置を返す
static int min_nonzero_loc( const int iv[] ) {
    int min=1000000, loc=-1;
    if ( iv[0] > 0 ) {
	min = iv[0];
	loc = 0;
    }
    if ( iv[1] > 0 ) {
	if ( iv[1] < min ) {
	    min = iv[1];
	    loc = 1;
	}
    }
    if ( iv[2] > 0 ) {
	if ( iv[2] < min ) loc = 2;
    }
    return loc;
}

/** a-1i, a+1i, 変化させる成分(Indx)のリストを作成
 * メモリも内部で確保する
 * */
static int make_nam() {
    int lqn, maxlqn = MAXLQN, iao, iao0, iao1, i;
    int a[3], naotype;
    // AOタイプの数を計算
    naotype = 0;
    for ( lqn=0; lqn<=maxlqn; lqn++ ) naotype += nnao[lqn];
    // メモリ確保
    indx = (int*)malloc( sizeof(int) * naotype );
    nam  = alloc_imatrix( naotype, 3 );
    nap  = alloc_imatrix( naotype, 3 );
    // まず、indxの作成
    for ( lqn=0; lqn<=maxlqn; lqn++ ) {
	iao0 = laot[lqn];
	iao1 = iao0 + nnao[lqn];
	for ( iao=iao0; iao<iao1; iao++ )
	    indx[iao] = min_nonzero_loc( angm[iao] );
    }

    // a-1iの作成
    for ( i=0; i<3; i++ ) nam[0][i] = -1;
    for ( lqn=1; lqn<=maxlqn; lqn++ ) {
	iao0 = laot[lqn];
	iao1 = iao0 + nnao[lqn];
	for ( iao=iao0; iao<iao1; iao++ ) {
	    memcpy( a, angm[iao], sizeof(int)*3 );
	    for ( i=0; i<3; i++ ) {
		a[i]--;
		nam[iao][i] = get_aotype( a[0], a[1], a[2] );
		a[i]++;
	    }
	}
    }
    // a+1iの作成
    for ( lqn=0; lqn<maxlqn; lqn++ ) {
	iao0 = laot[lqn];
	iao1 = iao0 + nnao[lqn];
	for ( iao=iao0; iao<iao1; iao++ ) {
	    memcpy( a, angm[iao], sizeof(int)*3 );
	    for ( i=0; i<3; i++ ) {
		a[i]++;
		nap[iao][i] = get_aotype( a[0], a[1], a[2] );
		a[i]--;
	    }
	}
    }
    iao0 = laot[maxlqn];
    iao1 = iao0 + nnao[maxlqn];
    for ( iao=iao0; iao<iao1; iao++ ) {
	for ( i=0; i<3; i++ ) nap[iao][i] = -1;
    }
    return 0;
}

/** nam, napが示す値から、先頭AOタイプ番号を引いたものにする
 * */
static int make_nam2() {
    int lqn, maxlqn = MAXLQN, iao, iao0, iao1, i;
    int naotype, lao;;
    // AOタイプの数を計算
    naotype = 0;
    for ( lqn=0; lqn<=maxlqn; lqn++ ) naotype += nnao[lqn];
    // メモリ確保
    nam2  = alloc_imatrix( naotype, 3 );
    nap2  = alloc_imatrix( naotype, 3 );
    // nam2の作成
    for ( lqn=0; lqn<=maxlqn; lqn++ ) {
	iao0 = laot[lqn];
	iao1 = iao0 + nnao[lqn];
	lao  = ( lqn==0 ? 0 : laot[lqn-1] );
	for ( iao=iao0; iao<iao1; iao++ ) {
	    for ( i=0; i<3; i++ ) {
		nam2[iao][i] = ( nam[iao][i] >= 0 ? nam[iao][i] - lao : -1 );
	    }
	}
    }
    // nap2の作成
    for ( lqn=0; lqn<=maxlqn; lqn++ ) {
	iao0 = laot[lqn];
	iao1 = iao0 + nnao[lqn];
	lao  = ( lqn==maxlqn ? 0 : laot[lqn+1] );
	for ( iao=iao0; iao<iao1; iao++ ) {
	    for ( i=0; i<3; i++ ) {
		nap2[iao][i] = ( nap[iao][i] >= 0 ? nap[iao][i] - lao : -1 );
	    }
	}
    }
    return 0;
}
/** (2n-1)!!を計算する
 * */
static int fact2( const int n ) {
    int istart, i, ifac;
    if ( n<2 ) return 1;
    istart = 2*n-1;
    ifac = 1;
    for ( i=istart; i>1; i-=2 ) ifac *= i;
    return ifac;
}

/** (2nx-1)!!(2ny-1)!!(2nz-1)!!を計算する
 * */
static int fact3( const int iv[] ) {
    return fact2(iv[0])*fact2(iv[1])*fact2(iv[2]);
}
/** 積分にかける定数の計算
 *
 * */
static int make_prefactor() {
    int lqn, maxlqn = MAXLQN, iao, iao0, iao1, i;
    int naotype;
    double base, dfac;
    // AOタイプの数を計算
    naotype = 0;
    for ( lqn=0; lqn<=maxlqn; lqn++ ) naotype += nnao[lqn];
    // メモリ確保
    dfact = (double*)malloc( sizeof(double) * naotype );
    // 計算開始
    // sについて
    dfact[0] = 1.e0;
    for ( lqn=1; lqn<=maxlqn; lqn++ ) {
	iao0 = laot[lqn];
	iao1 = iao0 + nnao[lqn];
	base = (double)fact3( angm[iao0] );
	for ( iao=iao0; iao<iao1; iao++ ) {
	    dfac = (double)fact3( angm[iao] );
	    dfact[iao] = sqrt( base / dfac );
	}
    }
    return 0;
}

/* 初期化関数 */
int ofmo_index_init( const int maxlqn ) {
    static int called = false;
    if ( called ) return 0;
    MAXLQN = 2*maxlqn;
    assign_angular_momentum_index( MAXLQN );
    make_nam();
    make_nam2();
    make_prefactor();
    atexit( ofmo_index_finalize );
    called = true;
    return 0;
}

/* --------------------------------------------------------
 * このクラス内部で定義されている変数を参照するための関数
 * -------------------------------------------------------- */
int** ofmo_getadd_angm() { return angm; }
int* ofmo_getadd_laot() { return laot; }
int* ofmo_getadd_nnao() { return nnao; }
int** ofmo_getadd_nam() { return nam; }
int** ofmo_getadd_nap() { return nap; }
int* ofmo_getadd_indx() { return indx; }
int** ofmo_getadd_nam2() { return nam2; }
int** ofmo_getadd_nap2() { return nap2; }
double* ofmo_getadd_dfact() { return dfact; }

// ofmo_index_initのテスト
/*int main( int argc, char* argv[] ) {
    int maxlqn;
    int nao, naot;
    int lqn, iao0, iao1, iao;
    char *ctype = "SPDFGHIJKLMNO";
    if ( argc == 1 ) {
	maxlqn = 2;
    } else {
	maxlqn = atoi( argv[1] );
	if ( maxlqn < 1 ) maxlqn = 2;
    }
    printf("maxlqn = %d\n", maxlqn );
    ofmo_index_init( maxlqn );
    // output
    for ( lqn=0; lqn<=MAXLQN; lqn++ ) {
	iao0 = laot[lqn];
	iao1 = iao0 + nnao[lqn];
	printf("%c-type (lqn= %d, nao= %2d, laot= %2d)\n",
		ctype[lqn], lqn, nnao[lqn], laot[lqn] );
	for ( iao=iao0; iao<iao1; iao++ ) {
	    printf(" type= %3d  angm= %2d %2d %2d",
		    iao, angm[iao][0], angm[iao][1], angm[iao][2] );
	    printf(" indx=%2d  fact=%6.3f nam= %2d %2d %2d  nap= %2d %2d %2d",
		    indx[iao], dfact[iao],
		    nam[iao][0], nam[iao][1], nam[iao][2],
		    nap[iao][0], nap[iao][1], nap[iao][2] );
	    printf(" nam2= %2d %2d %2d  nap2= %2d %2d %2d\n",
		    nam2[iao][0], nam2[iao][1], nam2[iao][2],
		    nap2[iao][0], nap2[iao][1], nap2[iao][2] );
	}
    }

    free( nnao );
    free( laot );
    free_imatrix( angm );
    return 0;
}*/
