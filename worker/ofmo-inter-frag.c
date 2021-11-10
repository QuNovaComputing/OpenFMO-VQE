/**
 * @file ofmo-inter-frag.c
 * It defines a group of functions related to the calculation of intermonomer distance,
 * which determines the approximate level of Coulomb interaction in fragment electronic state calculation.
 *
 * TODO:
 * Is the float type sufficient for the inter-monomer distance instead of the double type?
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

/** A function that returns the approximation level in the environmental potential calculation
 *
 * Fragment (monomer, dimer, etc. (\ f $ x \ f $)) in FMO calculation
 * When performing electronic state calculation, the electrostatic interaction term
 * @f $ {from the peripheral monomer (\ f $ I \ f $) } ^ x \ mbox {\ boldmath $ V $} _I @ f $
 * needs to be calculated.
 * In order to reduce the amount of calculation, this environmental potential term performs
 * an approximate calculation according to the distance between the fragment for which
 * the electronic state calculation is performed and the target monomer for which
 * the environmental potential calculation is performed.
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
 * If there is no approximation, it is necessary to calculate the 4-center Coulomb integral
 * with the same order of calculation as the normal 2-electron integral, but if the first
 * approximation (pop approximation) is performed, the amount of calculation is small.
 * All you have to do is calculate. Furthermore, if the second approximation
 * (point charge approximation) is used, only the two-center Coulomb integral,
 * which is a kind of one-electron integral, can be calculated.
 * This switching is based on the distance between the fragment and the monomer.
 * Inside this function, an approximate level list of the environmental potential between
 * the fragment specified by \ c nmonomer and \ c monomer_list [] and the surrounding
 * monomers is generated. This list can be referenced using the \ c ofmo_get_approx_level function.
 * Also, via the pointer variable of the argument of this function,
 * a list of the number of partner monomers and partner monomer numbers
 * that require 4-center Coulomb integration (\ c * nifc4c, \ c joblist_ifc4c [])
 * and 3-center Coulomb integration are required. The number of partner monomers
 * and the list of partner monomer numbers (\ c * nifc3c, \ c joblist_ifc3c [])
 * can be returned to the caller.
 *
 * @param[in] nmonomer Number of monomers that make up the fragment
 * @param[in] monomer_list[] List of monomer numbers that make up the fragment
 * @param[out] *nifc4c Number of non-approximate partner monomers requiring 4 central integration
 * @param[out] joblist_ifc4c[] List of non-approximate monomer numbers that require 4-center integration
 * @param[out] *nifc3c Number of partner monomers that perform pop approximation that requires 3-center integration
 * @param[out] joblist_ifc3c[] List of partner monomer numbers for pop approximation requiring 3-center integration
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
