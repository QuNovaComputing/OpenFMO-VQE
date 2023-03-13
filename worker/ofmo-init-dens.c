#include <stdio.h>
#include <stdlib.h>

#include <time.h>

#include <mpi.h>

#ifdef _OPENMP
#include <omp.h>
#else
#include "omp-dummy.h"
#endif

#include "ofmo-def.h"
#include "ofmo-data.h"
#include "ofmo-integ.h"
#include "ofmo-scf.h"
#include "ofmo-mat.h"

#ifdef DEBUG_MODE
extern FILE *fp_debug;
#endif

static double *_S_ = NULL;
static double *_H_ = NULL;
static int _maxnfao_ = 0;

static void dealloc() {
    Free( _S_ );
    Free( _H_ );
    _maxnfao_ = 0;
}

static int alloc( int nao ) {
    static int called = false;
    if ( nao > _maxnfao_ ) {
	int nao2;
	dealloc();
	nao2 = (nao*nao+nao)>>1;
	_S_ = (double*)malloc( sizeof(double) * nao2 );
	_H_ = (double*)malloc( sizeof(double) * nao2 );
	_maxnfao_ = nao;
    }
    if ( !called ) {
	atexit( dealloc );
	called = true;
    }
    return 0;
}

/** 指定されたモノマーの初期密度行列を計算する
 *
 * @arg[in] myfrag モノマー番号
 * @arg[out] D[] 初期密度行列（圧縮形式）
 * @arg[out] aop[] AO population
 * @arg[out] atp[] atomic population
 *
 * */
int ofmo_monomer_initial_density( int myfrag, double D[],
	double aop[], double atp[] ) {
    static int nfrag, *nfao, maxnfao, maxnfatom;
    static int maxlqn, **mleading_cs;
    static int **mshel_tem, **mshel_atm, **mshel_add, **mshel_ini;
    static double **matom_x, **matom_y, **matom_z;
    static double **mprim_exp, **mprim_coe;
    static int *nfatom, **matomic_number, *nfcs, *icharge;
    static int **msao2tuao, natom, nao, **ifatom;
    //double *SP;
    static int called = false;

    if ( !called ) {
	int ierr;
	ierr = ofmo_data_get_vals(
		"nfrag nfao maxnfao maxlqn maxnfatom "
		"mlcs mshel_tem mshel_atm mshel_add mshel_ini "
		"mprim_exp mprim_coe icharg "
		"matx maty matz matn nfatom "
		"msao2tuao nao natom ifatom nfcs",
		&nfrag, &nfao, &maxnfao, &maxlqn, &maxnfatom,
		&mleading_cs, &mshel_tem, &mshel_atm, &mshel_add,
		&mshel_ini,
		&mprim_exp, &mprim_coe, &icharge,
		&matom_x, &matom_y, &matom_z, &matomic_number, &nfatom,
		&msao2tuao, &nao, &natom, &ifatom, &nfcs );
	if ( ierr != 0 ) return -1;
	alloc( maxnfao );
	called = true;
    }
    if ( myfrag < 0 || myfrag >= nfrag ) { return -1; }
    double *S = _S_, *H=_H_;
    int nocc;
    nocc = ofmo_isum( nfatom[myfrag], matomic_number[myfrag] );
    nocc -= icharge[myfrag];
    nocc >>= 1;

#pragma omp parallel
    {
	int nthreads, mythread;
	nthreads = omp_get_num_threads();
	mythread = omp_get_thread_num();
	ofmo_integ_oneint_sorted( nthreads, mythread, maxlqn,
		mleading_cs[myfrag], mshel_tem[myfrag], mshel_atm[myfrag],
		mshel_add[myfrag], mshel_ini[myfrag],
		matom_x[myfrag], matom_y[myfrag], matom_z[myfrag],
		mprim_exp[myfrag], mprim_coe[myfrag], nfatom[myfrag],
		matomic_number[myfrag], S, H );
    }
    ofmo_scf_init_density_ehuckel( nfatom[myfrag], nfcs[myfrag],
	    nfao[myfrag], maxlqn, nocc,
	    matomic_number[myfrag], mleading_cs[myfrag],
	    mshel_atm[myfrag], mshel_ini[myfrag],
	    S, D, aop, atp );

    return 0;
}
