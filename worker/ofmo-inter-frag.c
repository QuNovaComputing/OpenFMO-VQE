/**
 * @file ofmo-inter-frag.c
 * フラグメント電子状態計算のクーロン相互作用の近似レベルなどを
 * 決めるモノマー間距離の計算などに関する関数群を定義している。
 *
 * TODO:
 * モノマー間距離はdouble型ではなく、float型で十分では？
 *
 * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ofmo-def.h"
#include "ofmo-data.h"
#include "ofmo-mserv-cont.h"

/* Wan der Waals半径（１オフセット）*/
static double Van_Der_Waals_Radius[] = {0.0,
    1.20,                              1.20,
    1.37,1.45,1.45,1.50,1.50,1.40,1.35,1.30, 
    1.57,1.36,1.24,1.17,1.80,1.75,1.70,2.50, 
    2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50, 
    2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50, 
    2.30,2.50, 
    2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50, 
    2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50, 
    2.50,2.50, 
    2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50, 
    2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50, 
    2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50, 
    2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50, 
    2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50, 
    2.50,2.50,2.50,2.50,2.50,2.50 };

static void init_vdw() {
    double inv_bohr;
    int n, i;
    inv_bohr = 1.e0 / BOHR_RADIUS;
    n = sizeof(Van_Der_Waals_Radius) / sizeof(double);
    for ( i=0; i<n; i++ ) Van_Der_Waals_Radius[i] *= inv_bohr;
}

/** 特定のモノマーが関係するモノマー間距離を求める関数
 *
 *
 * */
static int ofmo_calc_inter_fragment_distance_v0(
	const int ifrag, double dist[] ) {
    static int called = false;
    static double *atom_x, *atom_y, *atom_z;
    static int *atomic_number, **ifatom, *nfatom;
    static int nfrag;
    if ( !called ) {
	int ierr;
	ierr = ofmo_data_get_vals( "nfrag atn atx aty atz ifatom nfatom",
		&nfrag, &atomic_number, &atom_x, &atom_y, &atom_z,
		&ifatom, &nfatom );
	init_vdw();
	if ( ierr != 0 ) return -1;
	called = true;
    }
    int jfrag;
#pragma omp parallel for
    for ( jfrag=0; jfrag<nfrag; jfrag++ ) {
	int iat, iatm, iatn, jat, jatm, jatn;
	double rix, riy, riz, rijx, rijy, rijz, rij2, rij;
	double vdw_min, vdwi, vdwj, vdw_dist;
	if ( jfrag == ifrag ) {
	    dist[jfrag] = 0.e0;
	} else {
	    vdw_min = HUGE_VAL;
	    for ( iat=0; iat<nfatom[ifrag]; iat++ ) {
		iatm = ifatom[ifrag][iat];
		iatn = atomic_number[iatm];
		rix  = atom_x[iatm];
		riy  = atom_y[iatm];
		riz  = atom_z[iatm];
		vdwi = Van_Der_Waals_Radius[iatn];
		for ( jat=0; jat<nfatom[jfrag]; jat++ ) {
		    jatm = ifatom[jfrag][jat];
		    jatn = atomic_number[jatm];
		    vdwj = Van_Der_Waals_Radius[jatn];
		    rijx = rix - atom_x[jatm];
		    rijy = riy - atom_y[jatm];
		    rijz = riz - atom_z[jatm];
		    rij2 = rijx*rijx + rijy*rijy + rijz*rijz;
		    rij  = sqrt( rij2 );
		    vdw_dist = rij / (vdwi+vdwj);
		    if ( vdw_dist < vdw_min ) vdw_min = vdw_dist;
		}
	    }
	    dist[jfrag] = vdw_min;
	}
    }
    return 0;
}

static int    *itmp1    = NULL;
static double *dtmp1	= NULL;
static double *dtmp2    = NULL;

static void dealloc() {
    Free( itmp1 );
    Free( dtmp1 );
    Free( dtmp2 );
}

static int alloc( int nfrag ) {
    static int called = false;
    if ( !called ) {
	dealloc();
	itmp1 = (int*)malloc( sizeof(int) * nfrag );
	dtmp1 = (double*)malloc( sizeof(double) * nfrag );
	dtmp2 = (double*)malloc( sizeof(double) * nfrag );
	atexit( dealloc );
	called = true;
    }
    return 0;
}

/** モノマー間距離を計算する
  ワーカー（コミュニケータ）内のプロセスを用いて、モノマー間距離の
  リストを作成する。
  一回だけ呼び出せばよい。また、１つのワーカーから呼び出すだけでよい。
 */
int ofmo_make_inter_frag_distance_list( MPI_Comm comm ) {
    int nfrag, ierr, ifrag;
    int myrank, nprocs;
    double *inter_fragment_distance = NULL;
    if ( ofmo_data_get_vals("nfrag", &nfrag ) != 0 ) return -1;
    alloc( nfrag );
    inter_fragment_distance = dtmp1;
    MPI_Comm_rank( comm, &myrank );
    MPI_Comm_size( comm, &nprocs );
    /* モノマー間距離の計算 */
    for ( ifrag=myrank; ifrag<nfrag; ifrag+=nprocs ) {
	ierr = ofmo_calc_inter_fragment_distance_v0( ifrag,
		inter_fragment_distance );
	if ( ierr != 0 ) return -1;
	ofmo_worker_put( OFMO_DISTA, ifrag, inter_fragment_distance );
    }
    MPI_Barrier( comm );
    return 0;
}

/* Is monomer ifrag in fragment (true) or not (false) */
static int is_in_fragment( const int ifrag,
	const int nmonomer, const int monomer_list[] ) {
    for ( int i=0; i<nmonomer; i++ ) {
	if ( monomer_list[i] == ifrag ) return true;
    }
    return false;
}

/** 環境ポテンシャル計算における近似レベルを返す関数
 *
 * FMO計算におけるフラグメント（モノマー、ダイマーなど(\f$x\f$)）
 * 電子状態計算を行う場合には、周辺モノマー(\f$I\f$)からの静電相互作用項
 * @f${}^x\mbox{\boldmath$V$}_I@f$
 * を計算する必要がある。
 * この環境ポテンシャル項は、計算量の削減を行うために、電子状態計算を
 * 行うフラグメントと環境ポテンシャル計算を行う対象モノマー間との距離に
 * 応じて、近似計算を行う。
 * @f{eqnarray*}{
 * \left( {}^x\mbox{\boldmath$V$}_I \right)_{\mu\nu} &=&
 * \sum_{\sigma\lambda} (\mbox{\boldmath$D$}_{I})_{\sigma\lambda}
 * (\mu\nu|\sigma\lambda) \\
 *  &\approx& \sum_{\sigma} \left( \mbox{\boldmath$D$}_I
 *  \mbox{\boldmath$S$}_I \right)_{\sigma\sigma} (\mu\nu|\sigma\sigma) \\
 *  &\approx& \sum_{A\in I} \int dr \phi_\mu(r) \frac{Q_A}{|r-A|}
 *  \phi_\nu(r)
 * @f}
 * 
 * 近似なしの場合には通常の２電子積分と同じオーダーの計算量の
 * ４中心クーロン積分を計算する必要があるが、1つ目の近似（pop近似）を
 * 行うと計算量が少ない３中心クーロン積分の計算だけで済む。さらに、
 * 二つ目の近似（点電荷近似）を用いると、１電子積分の1種である２中心
 * クーロン積分だけの計算で済む。この切り替えは、フラグメントと
 * モノマー間の距離を元にして行う。
 * この関数内部で、\c nmonomer と\c monomer_list[] で指定された
 * フラグメントと周辺モノマー間の環境ポテンシャルの近似レベルリスト
 * が生成される。このリストは\c ofmo_get_approx_level 関数を用いて
 * 参照することが出来る。
 * また、この関数の引数のポインタ変数を経由して、４中心クーロン積分が
 * 必要な相手モノマー数と相手モノマー番号のリスト（\c *nifc4c,
 * \c joblist_ifc4c[] ）と、３中心クーロン積分が必要な相手モノマー数と
 * 相手モノマー番号のリスト（\c *nifc3c, \c joblist_ifc3c[] ）とを、
 * 呼び出し元に返すことが出来る。
 *
 * @param[in] nmonomer フラグメントを構成するモノマー数
 * @param[in] monomer_list[] フラグメントを構成するモノマー番号のリスト
 * @param[out] *nifc4c ４中心積分が必要な近似なしの相手モノマー数
 * @param[out] joblist_ifc4c[] ４中心積分が必要な近似なしの相手モノマー
 * 番号のリスト
 * @param[out] *nifc3c ３中心積分が必要なpop近似を行う相手モノマー数
 * @param[out] joblist_ifc3c[] ３中心積分が必要なpop近似を行う相手モノマー
 * 番号のリスト
 *
 * @ingroup ofmo-calc
 *
 * */
int ofmo_make_approx_level( const int nmonomer, const int monomer_list[],
	int *nifc4c, int joblist_ifc4c[],
	int *nifc3c, int joblist_ifc3c[],
	MPI_Comm comm ) {
    static int nfrag, called = false;
    static double laop, lptc;
    if ( !called ) {
	ofmo_data_get_vals("nfrag laop lptc", &nfrag, &laop, &lptc);
	alloc( nfrag );
	called = true;
    }
    double *rmin = dtmp1, *dist = dtmp2;
    int *approx_level = itmp1;
    int i, ifrag, jfrag, idum[2], root=0;
    /* make rmin list for fragment */
    int myrank;
    MPI_Comm_rank( comm, &myrank );
    if ( myrank == root ) {
	ifrag = monomer_list[0];
	ofmo_worker_get( OFMO_DISTA, ifrag, rmin );
	for ( i=1; i<nmonomer; i++ ) {
	    ifrag = monomer_list[i];
	    ofmo_worker_get( OFMO_DISTA, ifrag, dist );
	    for ( jfrag=0; jfrag<nfrag; jfrag++ ) {
		if ( rmin[jfrag] > dist[jfrag] ) rmin[jfrag] = dist[jfrag];
	    }
	}
	/* */
	int n4, n3;
	n4 = n3 = 0;
	for ( ifrag=0; ifrag<nfrag; ifrag++ ) {
	    if ( is_in_fragment( ifrag, nmonomer, monomer_list ) ) {
		approx_level[ifrag] = OFMO_IFC0C;
	    } else if ( rmin[ifrag] < laop ) {
		joblist_ifc4c[n4] = ifrag;
		approx_level[ifrag] = OFMO_IFC4C;
		n4++;
	    } else if ( rmin[ifrag] < lptc ) {
		joblist_ifc3c[n3] = ifrag;
		approx_level[ifrag] = OFMO_IFC3C;
		n3++;
	    } else {
		approx_level[ifrag] = OFMO_IFC2C;
	    }
	}
	*nifc4c = n4;
	*nifc3c = n3;
	idum[0] = n4;
	idum[1] = n3;
    }
    MPI_Bcast( idum, 2, MPI_INT, root, comm );
    if ( myrank != root ) {
	*nifc4c = idum[0];
	*nifc3c = idum[1];
    }
    MPI_Bcast( joblist_ifc4c, *nifc4c, MPI_INT, root, comm );
    MPI_Bcast( joblist_ifc3c, *nifc3c, MPI_INT, root, comm );
    MPI_Bcast( approx_level, nfrag, MPI_INT, root, comm );
    return 0;
}

/** 現在計算中のフラグメントと指定されたモノマー間の環境ポテンシャルの
 * 近似レベルを返す
 *
 * 直近で呼ばれた\c ofmo_make_approx_level 関数で与えられたフラグメント
 * と、この関数の引数\c ifrag で指定されたモノマーとの間の環境ポテンシャル
 * の近似レベルを返す。
 *
 * @param[in] ifrag 相手モノマーの番号
 *
 * @return 環境ポテンシャルの近似レベル
 *
 * @ingroup ofmo-calc
 *
 * */
int ofmo_get_approx_level( const int ifrag ) {
    return itmp1[ifrag];
}
