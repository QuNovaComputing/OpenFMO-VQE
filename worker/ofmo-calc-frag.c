/**
 * @file ofmo-calc-frag.c
 * @brief A file that describes functions related to fragment electronic state calculation
 * 
 * A file that describes functions related to fragment electronic state calculation
 * The calculation of the 2-center Coulomb integral is now divided for each integral type.
 * The 4-center Coulomb interaction term is now calculated even with worker = 0, except for the last one.
 * 
 * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <limits.h>

#ifdef USE_MPI
#include <mpi.h>
#else
#include "mpi-dummy.h"
#endif

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

#include "ofmo-monomer-data.h"
#include "ofmo-prof.h"
#include "ofmo-misc.h"

#include "ofmo-twoint.h"

#include "ofmo-vqe.h"

#ifdef USE_CUDA
#include "cuda/cuda-drv.h"
#include "cuda/cudalib.h"
#include "cuda/cuda-ifc4c.h"
void cuda_print_wifc4c(void);
#endif

//#define EPS_PS4 1.e-30
//#define EPS_ERI 1.e-15
//#define EPS_PS4 1.e-20
//#define EPS_ERI 1.e-12
//#define EPS_SCH 1.e-12
#define EPS_FAC_IFC4C 0.5e0

#define NPARTIAL 1000

#ifdef DEBUG_MODE
extern FILE* fp_debug;
#endif

// global counter
extern int ofmo_gc_init( const int id,
	MPI_Comm comm, const int init_val,
	const int njobs );
extern int ofmo_gc_nxtval( const int id );
//    debug
extern void ofmo_gc_set_debug_mode();
extern void ofmo_gc_reset_debug_mode();

extern int ofmo_projection_operator(
	const int nmonomer, const int monomer_list[],
	const int nao, const int sao2uao[], const double Ss[],
	double Ps[]);

extern int ofmo_fragment_init( int nmonomer, int monomer_list[] );

extern int ofmo_make_approx_level(
	const int nmonomer, const int monomer_list[],
	int *nifc4c, int joblist_ifc4c[],
	int *nifc3c, int joblist_ifc3c[],  MPI_Comm comm );
extern int ofmo_get_approx_level( const int ifrag );
extern int ofmo_monomer_initial_density(int frag, double D[], double aop[],
	double atp[] );


/* calc. and return the sum of integer array elements */
static int isum( const int n, const int ix[] ) {
    //int i, sum=0;
    //for ( i=0; i<n; i++ ) sum += ix[i];
    //return sum;
    return ofmo_isum(n, ix);
}

/* B += A */
static void acc_array( const int n, const double A[], double B[] ) {
    //for ( int i=0; i<n; i++ ) B[i] += A[i];
    ofmo_daxpy(n, 1.0e0, A, B);
}

static int is_in_fragment( const int ifrag, const int nmonomer,
	const int monomer_list[] ) {
    for ( int i=0; i<nmonomer; i++ )
	if ( monomer_list[i] == ifrag ) return true;
    return false;
}

/* --------------------------------------------------------------
 * ハイブリッド並列時のワーカー数、ワーカーIDの取得に関わる関数群
 * Function group related to acquisition of worker number and worker ID at the time of hybrid parallel
 * -------------------------------------------------------------- */
static int _nworkers_ = 1;
static int *_workerid_ = NULL;
static int *_nprocs_  = NULL;
static int *_itemp_   = NULL;

static void ofmo_mt_finalize() {
    Free( _nprocs_ );
    Free( _workerid_ );
    Free( _itemp_ );
}

static int ofmo_mt_init( MPI_Comm comm ) {
    static int called = false, maxthreads;
    if ( !called ) {
	int nprocs;
	MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
	_nprocs_ = (int*)malloc( sizeof(int) * nprocs );
	_itemp_  = (int*)malloc( sizeof(int) * nprocs );
	maxthreads = omp_get_max_threads();
	_workerid_ = (int*)malloc( sizeof(int) * maxthreads );
	atexit( ofmo_mt_finalize );
	called = true;
    }
    int myrank, nprocs;
    MPI_Comm_size( comm, &nprocs );
    MPI_Comm_rank( comm, &myrank );
    /* get # of threads in each process (_nprocs_[rank]) */
    for ( int i=0; i<nprocs; i++ ) _itemp_[i] = 0;
#pragma omp parallel
    {
	int nthreads;
	nthreads = omp_get_num_threads();
#pragma omp master
	_itemp_[myrank] = nthreads;
    }
    MPI_Allreduce( _itemp_, _nprocs_, nprocs, MPI_INT, MPI_SUM, comm );

    /* determine nworkers and workerid */
    int lwkid;
    _nworkers_ = isum( nprocs, _nprocs_ );
    if ( myrank==0 ) lwkid = 0;
    else             lwkid = isum( myrank, _nprocs_ );
    for ( int i=0; i<maxthreads; i++ ) _workerid_[i] = lwkid + i;
    return 0;
}

static int ofmo_mt_get_nworkers() { return _nworkers_; }

static int ofmo_mt_get_workerid( const int mythread ) {
    return _workerid_[mythread];
}
/* -----------------------------------------------------
 * モノマーをAO数の大きい順に並び替えたリストの作成
 * Creating a list of monomers sorted in descending order of AO number
 * ----------------------------------------------------- */
/** ２つの整数を比較する関数
 *  A function that compares two integers
 * */
static int comp2( const void *p1, const void *p2 ) {
    return ( (*(int*)p2) - (*(int*)p1) );	// good
}

static int *frag_order = NULL;

static void dealloc_frag_order() {
    if ( frag_order != NULL ) free( frag_order );
    frag_order = NULL;
}

static int init_frag_order() {
    static int called = false;
    if ( called ) return 0;
    int nfrag, i, i2, *nfao;
    if ( ofmo_data_get_vals("nfrag nfao", &nfrag, &nfao) != 0 ) {
	dbg("error\n");
	return -1;
    }
    frag_order = (int*)malloc( sizeof(int) * nfrag *2 );
    if ( frag_order == NULL ) return -1;
    for ( i=0, i2=0; i<nfrag; i++, i2+=2 ) {
	frag_order[i2+0] = nfao[i];
	frag_order[i2+1] = i;
    }
    qsort( frag_order, nfrag, sizeof(int)*2, comp2 );
    for ( i=0, i2=0; i<nfrag; i++, i2+=2 )
	frag_order[i] = frag_order[i2+1];
    atexit( dealloc_frag_order );
    called = true;
    return 0;
}

static double *_Sfrag_ = NULL;	/* overlap matrix of fragment */
static double *_Hfrag_ = NULL;	/* 1e-Hamilton matrix of fragment */
static double *_Ufrag_ = NULL;	/* external potential matrix of fragment */
static double *_Pfrag_ = NULL;	/* proj. operator matrix of fragment */
static double *_D0_    = NULL;	/* initial density matrix of fragment */
static double *_D_     = NULL;	/* density matrix of fragment */
static double *_ev_    = NULL;
static double *_WORK_  = NULL;
static double *_C_     = NULL;	/* MO coefficient matrix */

static double **_Dmon_ = NULL;	/* monomer density matrices for IFC4C */

/* cutoff table data of fragment */
static int *lcs_pair       = NULL;
static int *csp_ics        = NULL;
static int *csp_jcs        = NULL;
static int *csp_lps_pair   = NULL;
static double *csp_schwarz = NULL;
//
static double *psp_zeta    = NULL;
static double *psp_dkps    = NULL;
static double *psp_xiza    = NULL;
/* cutoff table data of monomers of IFC4C */
static int **lcs_pair_mon       = NULL;
static int **csp_ics_mon        = NULL;
static int **csp_jcs_mon        = NULL;
static int **csp_lps_pair_mon   = NULL;
static double **csp_schwarz_mon = NULL;
//
static double **psp_zeta_mon    = NULL;
static double **psp_dkps_mon    = NULL;
static double **psp_xiza_mon    = NULL;

/* Array that needs to be allocated for each thread */
//static int _maxnifc4c_;
static double **dUmaster = NULL;
static double **dUmaster2 = NULL;
static double **_atpop_local_ = NULL;

// AO population and atomic population of all monomers
static double **aopop_mon = NULL;
static double **atpop_mon = NULL;

// Fragment population data
// This data frees space with ofmo-data
static double *_daopop_ = NULL;
static double *_datpop_ = NULL;

/* vector of integer used to make initial density matrix */
static int *_aoconv_ = NULL;

/* job list of IFC4C and IFC3C */
static int *_joblist_ifc4c_ = NULL;
static int *_joblist_ifc3c_ = NULL;

/* monomer energy0 */
static double *_menergy0_   = NULL;

static void dealloc() {
    Free( _Sfrag_ );
    Free( _Hfrag_ );
    Free( _Ufrag_ );
    Free( _Pfrag_ );
    Free( _WORK_ );
    Free( _D0_ );
    Free( _D_ );
    Free( _ev_ );
    Free( _C_ );
    ofmo_free_dmatrix( _Dmon_ ); _Dmon_ = NULL;

    Free( lcs_pair );
    Free( csp_ics );
    Free( csp_jcs );
    Free( csp_lps_pair );
    Free( csp_schwarz );
    Free( psp_zeta );
    Free( psp_dkps );
    Free( psp_xiza );
    ofmo_free_imatrix( lcs_pair_mon ); lcs_pair_mon = NULL;
    ofmo_free_imatrix( csp_ics_mon );  csp_ics_mon = NULL;
    ofmo_free_imatrix( csp_jcs_mon );  csp_jcs_mon = NULL;
    ofmo_free_imatrix( csp_lps_pair_mon );  csp_lps_pair_mon = NULL;
    ofmo_free_dmatrix( csp_schwarz_mon );  csp_schwarz_mon = NULL;
    ofmo_free_dmatrix( psp_zeta_mon ); psp_zeta_mon = NULL;
    ofmo_free_dmatrix( psp_dkps_mon ); psp_dkps_mon = NULL;
    ofmo_free_dmatrix( psp_xiza_mon ); psp_xiza_mon = NULL;

    ofmo_free_dmatrix( dUmaster );
    ofmo_free_dmatrix( dUmaster2 );
    ofmo_free_dmatrix( aopop_mon );
    ofmo_free_dmatrix( atpop_mon );
    ofmo_free_dmatrix( _atpop_local_ );

    Free( _aoconv_ );

    Free( _joblist_ifc4c_ );
    Free( _joblist_ifc3c_ );

    Free( _menergy0_ );

    Free( _daopop_ );
    Free( _datpop_ );
}

static int alloc() {
    static int called = false;
    if ( called ) return 0;
    int ierr, maxlqn, maxnfatom, maxnfao, maxnfcs, maxnfps, maxnpspair;
    int nbody, nfrag, nao_total;
    int *nfao, *nfatom, total_nfao, total_nfatom;
    size_t total = 0, t;
    ierr = ofmo_data_get_vals(
	    "maxlqn maxnfatom maxnfcs maxnfao maxnfps "
	    "maxnpspair nbody nao nfrag nfao nfatom",
	    &maxlqn, &maxnfatom, &maxnfcs, &maxnfao, &maxnfps,
	    &maxnpspair, &nbody, &nao_total, &nfrag,
	    &nfao, &nfatom );
    if ( ierr != 0 ) {
	if ( fp_prof ) {
	    fdbg( fp_prof, "error\n");
	    fflush( fp_prof );
	}
	return -1;
    }
    int fnao2, nao2, maxlqn2, ncs2, fncs2, fnpspair;

    total_nfao   = ofmo_isum( nfrag, nfao );
    total_nfatom = ofmo_isum( nfrag, nfatom );

    nao2     = maxnfao * (maxnfao+1) / 2;
    fnao2    = nbody*maxnfao * (nbody*maxnfao + 1) / 2;
    ncs2     = maxnfcs * (maxnfcs+1) / 2;
    fncs2    = nbody*maxnfcs * (nbody*maxnfcs+1) / 2;
    maxlqn2  = (maxlqn+1) * (maxlqn+2) / 2;
    fnpspair = nbody * nbody * maxnpspair;
    /* memory allocation */
    // fragment
    t = sizeof(double) * fnao2;
    _Sfrag_ = (double*)malloc( t ); total += t;
    _Hfrag_ = (double*)malloc( t ); total += t;
    _Ufrag_ = (double*)malloc( t ); total += t;
    _Pfrag_ = (double*)malloc( t ); total += t;
    _WORK_  = (double*)malloc( t ); total += t;
    _D0_    = (double*)malloc( t ); total += t;
    _D_     = (double*)malloc( t ); total += t;
    t = sizeof(double) * nbody * maxnfao;
    _ev_    = (double*)malloc( t ); total += t;
    t = sizeof(double) * nbody * maxnfatom;
    t = sizeof(double) * nbody * nbody * maxnfao * maxnfao;
    _C_     = (double*)malloc( t ); total += t;
    // IFC4C monomers
    t = sizeof(double) * nao2;
    _Dmon_ = ofmo_alloc_dmatrix( MAXNIFC4C, nao2 );
    total += (t * MAXNIFC4C);
    /* cutoff table of fragment */
    t             = sizeof(int) * (maxlqn2+1+1);
    lcs_pair      = (int*)malloc( t ); total += t;
    t             = sizeof(int) * fncs2;
    csp_ics       = (int*)malloc( t ); total += t;
    csp_jcs       = (int*)malloc( t ); total += t;
    csp_lps_pair  = (int*)malloc( t ); total += t;
    t             = sizeof(double) * fncs2;
    csp_schwarz   = (double*)malloc( t ); total += t;
    t             = sizeof(double) * fnpspair;
    psp_zeta      = (double*)malloc( t ); total += t;
    psp_dkps      = (double*)malloc( t ); total += t;
    psp_xiza      = (double*)malloc( t ); total += t;
    /* cutoff table of IFC4C monomers */
    t = sizeof(int) * MAXNIFC4C * (maxlqn2+1+1);
    lcs_pair_mon     = ofmo_alloc_imatrix( MAXNIFC4C, (maxlqn2+1+1) );
    total           += t;
    t                = sizeof(int) * MAXNIFC4C * ncs2;
    csp_ics_mon      = ofmo_alloc_imatrix( MAXNIFC4C, ncs2 ); total += t;
    csp_jcs_mon      = ofmo_alloc_imatrix( MAXNIFC4C, ncs2 ); total += t;
    csp_lps_pair_mon = ofmo_alloc_imatrix( MAXNIFC4C, ncs2 ); total += t;
    t                = sizeof(double) * MAXNIFC4C * ncs2;
    csp_schwarz_mon  = ofmo_alloc_dmatrix( MAXNIFC4C, ncs2 ); total += t;
    t                = sizeof(double) * MAXNIFC4C * maxnpspair;
    psp_zeta_mon     = ofmo_alloc_dmatrix( MAXNIFC4C, maxnpspair);
    total += t;
    psp_dkps_mon     = ofmo_alloc_dmatrix( MAXNIFC4C, maxnpspair);
    total += t;
    psp_xiza_mon     = ofmo_alloc_dmatrix( MAXNIFC4C, maxnpspair);
    total += t;
    /* Array that needs to be allocated for each thread */
    int maxthreads;
    maxthreads = omp_get_max_threads();
    dUmaster = ofmo_alloc_dmatrix( maxthreads, fnao2 );
    dUmaster2 = ofmo_alloc_dmatrix( maxthreads, fnao2 );
    total += (maxthreads * fnao2)*sizeof(double);
    _atpop_local_ = ofmo_alloc_dmatrix( maxthreads, maxnfatom );
    total += (maxthreads * maxnfatom) * sizeof(double);

    /* Population of all monomers */
    aopop_mon = ofmo_alloc_dmatrixv( nfrag, nfao );
    total += total_nfao * sizeof(double);
    atpop_mon = ofmo_alloc_dmatrixv( nfrag, nfatom );
    total += total_nfatom * sizeof(double);
    memset(aopop_mon[0], '\0', sizeof(double)*total_nfao );
    memset(atpop_mon[0], '\0', sizeof(double)*total_nfatom );

    /* For AO order conversion */
    t = sizeof(int) * nao_total;
    _aoconv_ = (int*)malloc( t ); total += t;

    /* Job list in fragment electronic state calculation */
    t = sizeof(int) * nfrag;
    _joblist_ifc4c_ = (int*)malloc( t ); total += t;
    _joblist_ifc3c_ = (int*)malloc( t ); total += t;

    /* Monomer energy (excluding environmental potential term) */
    t = sizeof(double) * nfrag;
    _menergy0_ = (double*)malloc( t ); total += t;

    /* Fragment population data */
    _daopop_ = (double*)malloc( sizeof(double) * nbody * maxnfao );
    _datpop_ = (double*)malloc( sizeof(double) * nbody * maxnfatom );
    total += ( sizeof(double) * nbody * (maxnfao + maxnfatom) );


    // information for memory allocation
    if ( fp_prof ) {
	double dsize;
	dsize = (double)total / (double)(1024*1024);
	fprintf( fp_prof,
		"== allocd memory size in ofmo-calc-frag.c = %10.3f MB\n",
		dsize );
    }
    atexit( dealloc );
    called = true;
    return 0;
}

int ofmo_frag_init() {
    int maxnfao, nbody;
    ofmo_data_get_vals("maxnfao nbody", &maxnfao, &nbody );
    ofmo_scf_init( maxnfao*nbody );
    alloc();
    init_frag_order();
    return 0;
}

int ofmo_monomer_init_density( const int *imsg, MPI_Comm comm ) {
    int nprocs, myrank;
    int myfrag;
    static int called=false, nfrag, *nfao, *nfatom, maxnfao, maxnfatom;
    double *D, *aop, *atp;
    if ( !called ) {
	ofmo_data_get_vals("nfrag nfao nfatom maxnfao maxnfatom",
		&nfrag, &nfao, &nfatom, &maxnfao, &maxnfatom );
	called = true;
    }
    MPI_Comm_rank( comm, &myrank );
    MPI_Comm_size( comm, &nprocs );
    if ( imsg[6] >= nfrag ) return 0;
    D   = _D_;
    aop = _daopop_;
    atp = _datpop_;

    myfrag = imsg[6] + myrank;
    {
	if ( fp_prof ) {
	    int start, end;
	    start = imsg[6];
	    end   = start + nprocs;
	    if ( end > nfrag ) end = nfrag;
	    fprintf(fp_prof, "frag= ");
	    for ( int ifrag=start; ifrag<end; ifrag++ ) {
		fprintf( fp_prof, " %3d", ifrag );
	    }
	    fprintf( fp_prof, "\n");
	    fflush( fp_prof );
	}
    }
    if ( myfrag < nfrag ) {
	ofmo_monomer_initial_density( myfrag, D, aop, atp );
	ofmo_put_monomer_density( myfrag, D );
	ofmo_put_monomer_aopop( myfrag, aop );
	ofmo_put_monomer_atpop( myfrag, atp );
    }
    return 0;
}

static int ofmo_construct_init_density( MPI_Comm comm,
	const int nmonomer, const int monomer_list[],
	const int nao, const int fsao2tuao[],
	const int nao_total, int **msao2tuao, const int nfao[],
	double D[] ) {
    int myrank;
    MPI_Comm_rank( comm, &myrank );
    //is_root = ( myrank == root );
    if ( nmonomer == 1 ) {
	int ifrag;
	//nao2 = (nao*nao+nao)>>1;
	ifrag = monomer_list[0];
	ofmo_get_monomer_density( comm, ifrag, D );
	return 0;
    }
    /* initial condition */
    int *tuao2fsao=_aoconv_, nao2;
    for ( int iao=0; iao<nao_total; iao++ ) tuao2fsao[iao] = -1;
    for ( int iao=0; iao<nao; iao++ ) tuao2fsao[ fsao2tuao[iao] ] = iao;
    nao2 = (nao*nao+nao)>>1;
    memset( D, '\0', sizeof(double)*nao2 );
    /* construct fragment (dimer, trimer, ...) initial density */
    int i, ifrag, iao, jao, ijao, I, J, IJ, I2;
    double *Dmon = _Sfrag_; /* _Sfrag_ is used as temporary array */
    for ( i=0; i<nmonomer; i++ ) {
	ifrag = monomer_list[i];
	ofmo_get_monomer_density( comm, ifrag, Dmon );
	for ( iao=0, ijao=0; iao<nfao[ifrag]; iao++ ) {
	    if ( (I=tuao2fsao[ msao2tuao[ifrag][iao] ]) < 0 ) {
		dbg("error(ifrag=%d, iao=%d)\n", ifrag, iao );
		fflush(stdout);
		return -1;
	    }
	    I2 = (I*I+I)>>1; /* = I*(I+1)/2 */
	    for ( jao=0; jao<=iao; jao++, ijao++ ) {
		if ( (J = tuao2fsao[ msao2tuao[ifrag][jao] ]) < 0 ) {
		    dbg("error(ifrag=%d, iao=%d)\n", ifrag, jao );
		    fflush(stdout);
		    return -1;
		}
		IJ = (I>=J ? I2+J : ( ((J*J+J)>>1)+I ) );
		D[IJ] += Dmon[ijao];
	    }	/* for ( jao ) */
	}	/* for ( iao ) */
    }		/* for ( i<nmonomer ) */
    return 0;
}

extern size_t ofmo_twoint_get_max_nzeri( const int mythread );
extern size_t ofmo_twoint_get_stored_nzeri( const int mythread );

static int NEW_SCC_STEP = false;
void ofmo_set_new_scc_flag() { NEW_SCC_STEP = true; }
/** Function that performs RHF calculation of fragments
 *
 * In the FMO calculation, the RHF electronic state calculation is performed
 * in consideration of the environmental potential term from the peripheral
 * monomer and the projection operator term.
 * Inside this function, the following processing is mainly performed.
 *
 * - Cut-off table calculation for 2-electron integral and 4-center Coulomb integral
 * - Acquisition of required monomer density matrix data
 * - First two-electron integral calculation for the buffered SCF method
 *   (two-electron integral calculation is performed only in the memory of the fixed capacity)
 * - Calculation of 4-center Coulomb potential term
 * - Calculation of 3-center Coulomb potential term
 * - Calculation of 2-center Coulomb potential term
 * - (Normal) 1-electron integral calculation
 * - Calculation of projection operator terms
 * - SCF function call
 *
 * Of these, dynamic load balancing is applied to the 2-center Coulomb potential portion, 
 * and static load balancing is used for the other portions.
 *
 * @param[in] comm Communicator of worker group that performs electronic state calculation
 * @param[in] nmonomer Number of monomers that make up the fragment
 * @param[in] monomer_list[] List of monomer numbers that make up the fragment
 * @param[out] D[] Obtained density matrix (compressed U format)
 * @param[out] C[] Obtained MO coefficient matrix (square matrix)
 * @param[out] e[] Obtained MO energy (vector)
 * @param[out] aopop[] AO population by Mulliken population analysis
 * @param[out] atpop[] Atomic population by Mulliken population analysis
 * @param[out] *energy Energy of the fragment containing the environmental potential term
 * @param[out] *energy0 Fragment energy excluding the environmental potential term
 * @param[out] *ddv Amount of change in environmental potential in dimer SCF calculation
 * @param[out] daopop[] Amount of change in AO population in dimer SCF calculation
 * @param[out] datpop[] Amount of change in atomic population in dimer SCF calculation
 * \f[
 *   \tt{*ddv} = \rm{Tr}\,
 *   \left( {\boldmath D \unboldmath}_{IJ} -
 *   {\boldmath D \unboldmath}_I - {\boldmath D \unboldmath}_J \right)
 *   {\boldmath U \unboldmath}_{IJ}
 *
 * \f]
 * @param[in] iscc Num SCC iteration
 *
 * @ingroup ofmo-calc
 * */
int ofmo_calc_fragment_electronic_state(
	MPI_Comm comm, int nmonomer, int monomer_list[], int level,
	double tolscf,
	double *energy, double *energy0, double *ddv,
	int *fnao, double daopop[], int sao2tao[],
	int *fnat, double datpop[], int fatom2tatom[], const int iscc ) {
    static int nfrag, maxnfatom, maxnfcs, maxnfao, maxnfps, maxscf;
    static int total_nfao, total_nfatom;
    static int maxlqn, **mlcs;
    static int *nfatom, *nfcs, *nfao, *nfps;
    static int **mshel_tem, **mshel_atm, **mshel_add, **mshel_ini;
    static double **matom_x, **matom_y, **matom_z;
    static double **mprim_exp, **mprim_coe;
    static int **matomic_number, *icharge, **ifatom;
    static int maxnpspair, nbody, **msao2muao, **msao2tuao;
    static int nao_total;
    static int maxlqn2;
    static int itypes[6*3];
    static int type4c[6*6];
    static int called = false;
    static long nintic = 0;
    static int itol=30, icut=15;

    if(fp_prof) {fprintf(fp_prof,"--- calc_fragment^%d (tolscf=%8.2e) ---\n",nmonomer,tolscf);fflush(fp_prof);}
    if ( !called ) {
	int ierr;
	ierr = ofmo_data_get_vals(
                "nintic itol icut "
		"nfrag maxlqn maxnfatom maxnfcs maxnfao maxnfps "
		"nfatom nfcs nfao nfps "
		"mlcs mshel_tem mshel_atm mshel_add mshel_ini "
		"mprim_exp mprim_coe icharg "
		"matx maty matz matn "
		"maxnpspair nbody ifatom msao2muao msao2tuao maxscf nao",
                &nintic, &itol, &icut,
		&nfrag, &maxlqn, &maxnfatom, &maxnfcs, &maxnfao, &maxnfps,
		&nfatom, &nfcs, &nfao, &nfps,
		&mlcs, &mshel_tem, &mshel_atm, &mshel_add, &mshel_ini,
		&mprim_exp, &mprim_coe, &icharge,
		&matom_x, &matom_y, &matom_z, &matomic_number,
		&maxnpspair, &nbody, &ifatom, &msao2muao, &msao2tuao,
		&maxscf, &nao_total);
	if ( ierr != 0 ) return -1;
    // FIXME : need modification for dimer here.
    if ( level == OFMO_RHF_VQE || level == OFMO_VQE_RHF ) level = OFMO_VQE;
	total_nfao   = ofmo_isum( nfrag, nfao );
	total_nfatom = ofmo_isum( nfrag, nfatom );
	maxlqn2 = ((maxlqn+1)*(maxlqn+2))>>1;
	if ( maxlqn == 0 ) {
	    itypes[0] = type4c[0] = 0;
	} else if ( maxlqn <= 2 ) {
	    int Lab, Lcd, Lc, ix;
	    ix=0;
	    for ( Lab=0; Lab<maxlqn2; Lab++ ) {
		for ( Lc=0; Lc<=maxlqn; Lc++ ) {
		    itypes[ix] = Lab*3 + Lc;
		    ix++;
		}
	    }
	    ix=0;
	    for ( Lab=0; Lab<maxlqn2; Lab++ ) {
		for ( Lcd=0; Lcd<maxlqn2; Lcd++ ) {
		    type4c[ix] = Lab*6 + Lcd;
		    ix++;
		}
	    }
	} else return -1;
#pragma omp parallel
	{
	    int mythread;
	    mythread = omp_get_thread_num();
	    (void)ofmo_twoint_alloc_local_gmat( mythread, maxnfao*nbody );
	}
	called = true;
    }
    /* Acquisition of MPI information */
    int myrank, nprocs;
    int root = 0, is_root;
    MPI_Comm_rank( comm, &myrank );
    MPI_Comm_size( comm, &nprocs );
    is_root = ( myrank == root );
    /* Get profile ID */
    static int cid_cutoff, cid_eri, cid_4c, cid_3c, cid_2c, cid_buf;
    static int tid_init, tid_cutoff, tid_integ, tid_comm;	// 詳細
    static int tid_Init, tid_Integ, tid_SCF, tid_Total;
    cid_cutoff = ofmo_create_thread_timer( "CUTOFF", 0 );
    cid_eri    = ofmo_create_thread_timer( "ERI", 0 );
    cid_4c     = ofmo_create_thread_timer( "IFC4C", 0 );
    cid_3c     = ofmo_create_thread_timer( "IFC3C", 0 );
    cid_2c     = ofmo_create_thread_timer( "IFC2C", 0 );
    cid_buf    = ofmo_create_thread_timer( "BUF", 1 );

    tid_init   = ofmo_create_proc_timer( "init", 0 );
    tid_cutoff = ofmo_create_proc_timer( "cutoff", 0 );
    tid_integ  = ofmo_create_proc_timer( "integ", 0 );
    tid_comm   = ofmo_create_proc_timer( "comm", 0 );
    
    tid_Init   = ofmo_create_proc_timer( "INIT", 1 );
    tid_Integ  = ofmo_create_proc_timer( "INTEG", 1 );
    tid_SCF    = ofmo_create_proc_timer( "SCF", 1 );
    tid_Total  = ofmo_create_proc_timer( "TOTAL", 1 );


    ofmo_start_proc_timer( tid_init );
    ofmo_start_proc_timer( tid_Init );
    ofmo_start_proc_timer( tid_Total );

    /* Acquisition of atomic data and basis function data of fragments */
    int nat, ncs, nao, nps, *atomic_number, *fat2tat, ierr;
    int *flcs, *shel_tem, *shel_atm, *shel_add, *shel_ini, *fsao2tuao;
    int *fsao2fuao;
    double *atom_x, *atom_y, *atom_z;
    double *prim_exp, *prim_coe;
    int charge;
    ofmo_fragment_init( nmonomer, monomer_list );
    ierr = ofmo_data_get_vals("fnatom fncs fnao fnps "
	    "fatn fat2tat fatx faty fatz "
	    "flcs fshel_tem fshel_atm fshel_add fshel_ini "
	    "fprim_exp fprim_coe fsao2tuao fsao2fuao",
	    &nat, &ncs, &nao, &nps,
	    &atomic_number, &fat2tat, &atom_x, &atom_y, &atom_z,
	    &flcs, &shel_tem, &shel_atm, &shel_add, &shel_ini,
	    &prim_exp, &prim_coe, &fsao2tuao, &fsao2fuao );
    if ( ierr != 0 ) {
	if ( fp_prof ) fdbg(fp_prof, "error\n");
	return -1;
    }
    charge = 0;
    for ( int i=0; i<nmonomer; i++ )
	charge += icharge[ monomer_list[i] ];

    /* number of electrons */
    int nelec, nocc, nao2;
    nao2 = nao*(nao+1)/2;
    nelec = isum( nat, atomic_number ) - charge;
    nocc = nelec / 2;
    if ( (nelec%2) != 0) {
        if ( fp_prof ) {
	    fprintf( fp_prof, "ERROR: fragment(");
	    for ( int i=0; i<nmonomer; i++ )
		fprintf( fp_prof, "%d%s", monomer_list[i],
			( (i==(nmonomer-1)? ") " : "," ) ) );
	    fprintf( fp_prof, ": odd number of electron (%d)\n", nelec );
        }
        //return -1;
    }

    double *S = _Sfrag_, *H = _Hfrag_, *P = _Pfrag_, *U = _Ufrag_;
    double *D = _D_, *C = _C_;
    double *ev = _ev_;
    memset( U, '\0', sizeof(double)*nao2 );
    // temporary
    //int eribfsz = 4096;
    //int eribfsz = 1200;
    //int eribfsz = 1024;
    //int eribfsz = 480;
    //int eribfsz = 2;
    size_t eribfsz = 0;
    //if (nintic>0) eribfsz = nintic*(sizeof(double)+4*sizeof(short))/1024/1024;
    if (nintic>0) eribfsz = nintic;
    else eribfsz = -nintic*8/1024/1024;
    if ( fp_prof ) {fprintf(fp_prof,"buffer size  = %ld\n", eribfsz); fflush(fp_prof);}
    eribfsz /= omp_get_max_threads();

    /* make joblist */
    int *joblist_ifc4c = _joblist_ifc4c_;
    int *joblist_ifc3c = _joblist_ifc3c_;
    int nifc4c, nifc3c;
    int njob_ifc3c, njob_ifc2c, uifc3c;
    int njob_ifc4c, uifc4c;
    ofmo_make_approx_level( nmonomer, monomer_list,
	    &nifc4c, joblist_ifc4c, &nifc3c, joblist_ifc3c, comm );
    /*// global counterを用いた動的負荷分散の準備
    uifc4c     = maxlqn2 * maxlqn2;
    njob_ifc4c = ((nifc4c-1) * uifc4c)<<7; // << 7 means *128
    if ( nifc4c == 0 ) njob_ifc4c = 0;
    uifc3c     = (maxlqn+1) * maxlqn2;
    njob_ifc3c = (nifc3c * uifc3c)<<4; // <<4 means *16
    njob_ifc2c = nfrag * maxlqn2;
    ofmo_gc_init( 2, comm, 0, njob_ifc4c );
    ofmo_gc_init( 0, comm, 0, njob_ifc3c );
    ofmo_gc_init( 1, comm, 0, njob_ifc2c );*/
    ofmo_mt_init( comm );

    ofmo_acc_proc_timer( tid_init );
    /* thread-parallelized calculation of cutoff tables */
    ofmo_start_proc_timer( tid_cutoff );
    double **Dmons=_Dmon_;
    double *D0 = _D0_;
    double Enuc;
    int mc = -2;	/* global counter in local process */
#ifndef PARA_SUB
#pragma omp parallel
#else
    int nnewjob;
    int *newjoblist;
    int amyrank, anprocs;
    int color;
    double *wjob;
#pragma omp parallel shared(nnewjob, newjoblist, amyrank, anprocs, color, wjob)
#endif
    {
#pragma omp critical // FIXME after finding the bug
    {
	int i, ifrag, flag=false;
	//int n;
	// for profile
	int mythread, nthreads;
	mythread = omp_get_thread_num();
	nthreads = omp_get_num_threads();
	ofmo_start_thread_timer( cid_cutoff, mythread );

    /*
    // not known bug
    int loop_idx = 0;
    if(monomer_list[0] == 1){// not known bug
        printf("%d/%d/%d loop start\n", monomer_list[0], iscc, mythread);
        fflush(stdout);
    }
    */
	while (1) {
#pragma omp master
	    {
		if ( NEW_SCC_STEP ) {
		    if ( is_root ) {
			ofmo_get_monomer_aopop( -1, aopop_mon[0] );
			ofmo_get_monomer_atpop( -1, atpop_mon[0] );
		    }
		    MPI_Bcast( aopop_mon[0], total_nfao, MPI_DOUBLE,
			    root, comm );
		    MPI_Bcast( atpop_mon[0], total_nfatom, MPI_DOUBLE,
			    root, comm );
		    NEW_SCC_STEP = false;
		}
        /*
        if(monomer_list[0] == 1){// not known bug
            printf("%d/%d init_density_loop0_%d\n", monomer_list[0], iscc, loop_idx);
            fflush(stdout);
        }
        */
		ofmo_construct_init_density( comm, nmonomer, monomer_list,
			nao, fsao2tuao,
			nao_total, msao2tuao, nfao,
			D );
        /*
        if(monomer_list[0] == 1){// not known bug
            printf("%d/%d init_density_loop1_%d\n", monomer_list[0], iscc, loop_idx);
            fflush(stdout);
        }
        */
		if ( nmonomer > 1 ) memcpy( D0, D, sizeof(double) * nao2 );
		/* read density matrices to be used ifc4c calculations */
		for ( i=0; i<nifc4c; i++ ) {
		    ifrag = joblist_ifc4c[i];
		    ofmo_get_monomer_density( comm, ifrag, Dmons[i] );
		}
        /*
        if(monomer_list[0] == 1){// not known bug
            printf("%d/%d init_density_loop2_%d\n", monomer_list[0], iscc, loop_idx);
            fflush(stdout);
        }
        */
		flag = ( nthreads == 1 ? false : true );
	    }	// pragma omp master
	    if ( flag == true ) break;
//#pragma omp critical
	    { i = mc; mc++; }
        /* Cut-off table calculation for 2-electron integral and 4-center Coulomb integral */
	    if ( i == -2 ) {
		/* nuclear repulsion */
        /*
        if(monomer_list[0] == 1){// not known bug
            printf("%d/%d init_density_loop3_%d\n", monomer_list[0], iscc, loop_idx);
            loop_idx += 1;
            fflush(stdout);
        }
        */
		Enuc = ofmo_calc_nuclear_repulsion( nat, atomic_number,
			atom_x, atom_y, atom_z );
	    } else if ( i == -1 ) {
		ofmo_cutoff_make_table( maxlqn, flcs, shel_tem,
			shel_atm, shel_add, atom_x, atom_y, atom_z,
			prim_exp, prim_coe,
			lcs_pair, csp_schwarz, csp_ics, csp_jcs,
			csp_lps_pair, psp_zeta, psp_dkps, psp_xiza );
	    } else {
            if ( i >= nifc4c ) break;
            ifrag = joblist_ifc4c[i];
            ofmo_cutoff_make_table( maxlqn, mlcs[ifrag],
                mshel_tem[ifrag], mshel_atm[ifrag],
                mshel_add[ifrag],
                matom_x[ifrag], matom_y[ifrag], matom_z[ifrag],
                mprim_exp[ifrag], mprim_coe[ifrag],

                lcs_pair_mon[i], csp_schwarz_mon[i],
                csp_ics_mon[i], csp_jcs_mon[i],
                csp_lps_pair_mon[i], psp_zeta_mon[i],
                psp_dkps_mon[i], psp_xiza_mon[i] );
	    }
    //
	}	// while
    /*
    if(monomer_list[0] == 1){// not known bug
        printf("%d/%d/%d loop break\n", monomer_list[0], iscc, mythread);
        fflush(stdout);
    }
    */
	ofmo_acc_thread_timer( cid_cutoff, mythread );
    /*
    if(monomer_list[0] == 1){// not known bug
        printf("%d/%d/%d thread timer\n", monomer_list[0], iscc, mythread);
        fflush(stdout);
    }
    */
    }   // pragma omp critical (FIXME after find the bug)
    }	// pragma omp parallel
    /*
    if(monomer_list[0] == 1){// not known bug
        printf("%d/%d timer0\n", monomer_list[0], iscc);
        fflush(stdout);
    }
    */
    ofmo_acc_proc_timer( tid_cutoff );
    /*
    if(monomer_list[0] == 1){// not known bug
        printf("%d/%d timer1\n", monomer_list[0], iscc);
        fflush(stdout);
    }
    */
    ofmo_acc_proc_timer( tid_Init );
    /*
    if(monomer_list[0] == 1){// not known bug
        printf("%d/%d timer2\n", monomer_list[0], iscc);
        fflush(stdout);
    }
    */
    ofmo_start_proc_timer( tid_integ );
    /*
    if(monomer_list[0] == 1){// not known bug
        printf("%d/%d timer3\n", monomer_list[0], iscc);
        fflush(stdout);
    }
    */
    ofmo_start_proc_timer( tid_Integ );
    /*
    if(monomer_list[0] == 1){// not known bug
        printf("%d/%d timer done\n", monomer_list[0], iscc);
        fflush(stdout);
    }
    */
    double scfd=tolscf, scfe;
    if      ( scfd <= 1.e-4 ) scfe = scfd * 1.e-2;
    else if ( scfd <= 1.e-3 ) scfe = scfd * 1.e-1;
    else                      scfe = scfd * 1.e-0;
    if (scfe<1e-10) scfe=1e-10;

#if 0
    float eps_ps4 = EPS_PS4;
    float eps_eri = EPS_ERI;
    float eps_sch = EPS_SCH;
#else
    float eps_ps4 = pow(10.0,-itol);
    float eps_sch = pow(10.0,-icut);
    float eps_eri = eps_sch;
#endif
//    float eps_fac = (scfe>1e-8)? scfe/1e-8: 1.0;
    float eps_fac = (scfe>1e-7)? scfe/1e-8: 0.1F; // Increase accuracy for direct & fdiff run
    eps_ps4 *= eps_fac;
    eps_eri *= eps_fac;
    eps_sch *= eps_fac;

    if(fp_prof) {fprintf(fp_prof, "scfd scfe eps_ps4 eps_sch eps_eri: %8.2e %8.2e %8.2e %8.2e %8.2e\n", scfd, scfe, eps_ps4, eps_sch, eps_eri); fflush(fp_prof);}
    ofmo_twoint_eps_ps4(eps_ps4*EPS_FAC_IFC4C);
    ofmo_twoint_eps_sch(eps_sch*EPS_FAC_IFC4C);
    ofmo_twoint_eps_eri(eps_eri*EPS_FAC_IFC4C);
    //printf("%d/%d init_density\n", monomer_list[0], iscc);
    //fflush(stdout);

#pragma omp parallel
    {
	int nworkers, workerid;
	int j, jfrag;
        int jn;
	int iat, k;
	double *dU, *atpop_mt;
	int local_id;
	int mythread, nthreads;
	// for control load-balancing
	int offset = 0;
	//
	int mytype, tfrag, mypos, base;
	mythread = omp_get_thread_num();
	nthreads = omp_get_num_threads();
	dU       = dUmaster[mythread];
	double *dUtmp    = dUmaster2[mythread];
	atpop_mt = _atpop_local_[mythread];
	memset( dU, '\0', sizeof(double)*nao2 );

	nworkers = ofmo_mt_get_nworkers();
	workerid = ofmo_mt_get_workerid( mythread );

	ofmo_start_thread_timer( cid_eri, mythread );
	// for control load-balancing
	ofmo_integ_set_loop_offset( mythread, offset );
	/* ERI calculation ( 1st time ) */
    /* - First two-electron integral calculation for the buffered SCF method */

    //printf("%d/%d eri_first\n", monomer_list[0], iscc);
    //fflush(stdout);

	ofmo_integ_twoint_first(
		nworkers, workerid, eribfsz,
		maxlqn, shel_atm, shel_ini, atom_x, atom_y, atom_z,
		lcs_pair, csp_schwarz, csp_ics, csp_jcs, csp_lps_pair,
		psp_zeta, psp_dkps, psp_xiza );
	ofmo_acc_thread_timer( cid_eri, mythread );
	{
	    int last_eri_type;
	    double dnmax, dnzeri, rate;
	    dnmax = (double)ofmo_twoint_get_max_nzeri( mythread );
	    dnzeri = (double)ofmo_twoint_get_stored_nzeri( mythread );
	    last_eri_type = ofmo_twoint_get_last_eri_type( mythread );
        //if (fp_prof) {fprintf(fp_prof,"(%d) %d %ld/%ld\n",mythread,last_eri_type,(long)dnzeri,(long)dnmax);fflush(fp_prof);};
	    if ( last_eri_type < 21 ) {	// d関数までを仮定
		rate = (double)last_eri_type;
	    } else if ( dnmax > 0.e0 ) rate = dnzeri / dnmax * 100.e0;
	    else                rate = 100.1e0;
	    ofmo_set_thread_timer( cid_buf, mythread, rate );
	}

    //printf("%d/%d pot_env_start\n", monomer_list[0], iscc);
    //fflush(stdout);

	/* environmental potential */
	/* 4-centered inter-fragment Coulomb term */
	ofmo_start_thread_timer( cid_4c, mythread );
	if ( myrank == 0 ) {
	    if ( mythread == 0 ) {
		memset( S, '\0', sizeof(double)*nao2 );
		memset( H, '\0', sizeof(double)*nao2 );
	    }
#pragma omp barrier
	    ofmo_integ_oneint_sorted( nthreads, mythread, maxlqn, flcs,
		    shel_tem, shel_atm, shel_add, shel_ini,
		    atom_x, atom_y, atom_z, prim_exp, prim_coe,
		    nat, atomic_number, S, H );
#pragma omp barrier
	    if ( mythread == 0 ) {
		ofmo_projection_operator( nmonomer, monomer_list,
			nao, fsao2tuao, S, P );
		acc_array(nao2, P, H);
		ofmo_scale_diag( nao, 0.5e0, H );
	    }
#pragma omp barrier
	}
#pragma omp master
    if(fp_prof) {fprintf(fp_prof,"(%d) calc_fragment ifc4c (%d)\n",myrank,nifc4c);}
#ifdef PARA_SUB
    int njob = nifc4c;
#pragma omp master
    {
      newjoblist=(int*)malloc(njob*sizeof(int));
      wjob=(double *)malloc(njob*sizeof(double)+njob*2*sizeof(int));
    }
#if 0
#pragma omp master
      for ( j=0; j<njob; j++ ) {
        int Lab = maxlqn*(maxlqn+1)/2 + maxlqn;
        wjob[j] = (double)(lcs_pair_mon[j])[Lab+1];
      }
#else
#pragma omp barrier
#pragma omp for
      for ( j=0; j<njob; j++ ) {
        int Lab = maxlqn*(maxlqn+1)/2 + maxlqn;
        int ncspair_f = lcs_pair[Lab+1];
        int ncspair_m = (lcs_pair_mon[j])[Lab+1];
        wjob[j]=(double)ncspair_f*ncspair_m;
        //double eps_ps4 = ofmo_twoint_eps_ps4(0);
        double eps_ps4_ifc4c = ofmo_twoint_eps_ps4(0);
        size_t nps4=0;
        for (int ii=0; ii<ncspair_f; ii++) {
          for (int jj=0; jj<ncspair_m; jj++)
            if (csp_schwarz[ii]*(csp_schwarz_mon[j])[jj]>=eps_ps4_ifc4c) nps4++;
        }
        wjob[j]=(double)nps4;
      }
#endif
#pragma omp master
    {
      int *atask=(int *)(wjob+njob);
      int *res=atask+njob;
      int natask=ofmo_aggregateTask(nprocs, njob, wjob, atask);
      //if (myrank==0) printf("%d %d -> %d %d\n",nprocs, njob, nprocs, natask);
      ofmo_assignRes(nprocs, natask, wjob);
      for ( j=0; j<natask; j++ ) res[j]=(int)wjob[j];
      color=ofmo_getMyColorAndRankFromAssignedRes(myrank,res,&amyrank);
      anprocs = res[color];
      nnewjob=0;
      for (j=0; j<njob; j++) if (atask[j]==color) newjoblist[nnewjob++]=j;
      free(wjob);
    }
#ifdef USE_CUDA
      cuda_Reconfig(amyrank, anprocs, comm);
#endif
#pragma omp barrier
    ofmo_integ_set_loop_offset( mythread, 0 );

#endif /* PARA_SUB */
#ifdef USE_CUDA
        float *csp_schwarz_f;
        {
          int nat_f=nat, nao_f=nao, ncs_f=ncs;
          int nat_m=0, nao_m=0, ncs_m=0;
          int Lab = maxlqn*(maxlqn+1)/2 + maxlqn;
          int ncspair_f = lcs_pair[Lab+1];
          int npspair_f = csp_lps_pair[ncspair_f];
          int ncspair_m = 0;
          int npspair_m = 0;
          int max_num_klcs = 0;
#ifndef PARA_SUB
          for ( j=0; j<nifc4c; j++ ) {
#else
          for (int aj=0; aj<nnewjob; aj++ ) {
            j = newjoblist[aj];
#endif
            jfrag = joblist_ifc4c[j];
            nat_m = MAX2(nat_m, nfatom[jfrag]);
            ncs_m = MAX2(ncs_m, nfcs[jfrag]);
            nao_m = MAX2(nao_m, nfao[jfrag]);
            int ncspair_m0 = (lcs_pair_mon[j])[Lab+1];
            ncspair_m = MAX2(ncspair_m, ncspair_m0);
            npspair_m = MAX2(npspair_m, (csp_lps_pair_mon[j])[ncspair_m0]);
            max_num_klcs = MAX2(max_num_klcs, cuda_max_num_klcs(maxlqn, lcs_pair_mon[j]));
          }
          csp_schwarz_f = (float*)malloc(sizeof(float)*MAX2(ncspair_f,ncspair_m));

          cuda_ifc4c_Init(maxlqn, max_num_klcs,
              nat_f, ncs_f, nao_f, ncspair_f, npspair_f,
              nat_m, ncs_m, nao_m, ncspair_m, npspair_m );

          for (int ii=0; ii<ncspair_f; ii++)
            csp_schwarz_f[ii]=(float)csp_schwarz[ii];
          cuda_ifc4c_SetData(0,
              maxlqn, nat_f, ncs_f, nao_f, ncspair_f, npspair_f,
              shel_atm, shel_ini, atom_x, atom_y, atom_z,
              lcs_pair, csp_lps_pair, csp_ics, csp_jcs,
              psp_zeta, psp_dkps, psp_xiza,
              csp_schwarz_f, NULL);
        }
#endif

	memset( dUtmp, '\0', sizeof(double)*nao2 );
#ifndef PARA_SUB
	for ( j=0; j<nifc4c; j++ ) {
#else
        for (int aj=0; aj<nnewjob; aj++ ) {
            j = newjoblist[aj];
#endif
	    jfrag = joblist_ifc4c[j];
#ifdef USE_CUDA
            int nat_mon = nfatom[jfrag];
            int ncs_mon = nfcs[jfrag];
            int nao_mon = nfao[jfrag];
            int Lab = maxlqn*(maxlqn+1)/2 + maxlqn;
            int ncspair_mon = (lcs_pair_mon[j])[Lab+1];
            int npspair_mon = (csp_lps_pair_mon[j])[ncspair_mon];
            for (int ii=0; ii<ncspair_mon; ii++)
              csp_schwarz_f[ii]=(float)(csp_schwarz_mon[j])[ii];
            cuda_ifc4c_SetData(1,
                maxlqn, nat_mon, ncs_mon, nao_mon,
                ncspair_mon, npspair_mon,
                mshel_atm[jfrag], mshel_ini[jfrag],
                matom_x[jfrag], matom_y[jfrag], matom_z[jfrag],
                lcs_pair_mon[j], csp_lps_pair_mon[j],
                csp_ics_mon[j], csp_jcs_mon[j],
                psp_zeta_mon[j], psp_dkps_mon[j], psp_xiza_mon[j],
                csp_schwarz_f, Dmons[j]);
#endif
#ifndef PARA_SUB
	    ofmo_integ_ifc4c_sorted_partial( nworkers, workerid,
#else
            int anworkers = anprocs * omp_get_num_threads();
            int aworkerid = amyrank * omp_get_num_threads() + omp_get_thread_num();
	    ofmo_integ_ifc4c_sorted_partial( anworkers, aworkerid,
#endif
		    maxlqn, shel_atm, shel_ini, atom_x, atom_y, atom_z,
		    lcs_pair, csp_schwarz, csp_ics, csp_jcs,
		    csp_lps_pair, psp_zeta, psp_dkps, psp_xiza,
		    mshel_atm[jfrag], mshel_ini[jfrag],
		    matom_x[jfrag], matom_y[jfrag], matom_z[jfrag],

		    lcs_pair_mon[j],
		    csp_schwarz_mon[j], csp_ics_mon[j], csp_jcs_mon[j],
		    csp_lps_pair_mon[j],
		    psp_zeta_mon[j], psp_dkps_mon[j], psp_xiza_mon[j],
                    mlcs[jfrag],
		    nfao[jfrag], Dmons[j], dUtmp );
	}
#ifdef USE_CUDA
        Free(csp_schwarz_f);
        cuda_ifc4c_GetVfrg(nao, dUtmp);
        cuda_ifc4c_Finalize();
#endif
#ifdef PARA_SUB
#pragma omp barrier
#pragma omp master
        free(newjoblist);
#ifdef USE_CUDA
        cuda_Reconfig(myrank, nprocs, comm);
#endif
        ofmo_integ_set_loop_offset( mythread, 0 );
#endif /* PARA_SUB */
        acc_array( nao2, dUtmp, dU );
	//for ( int ii=0; ii<nao2; ii++ ) dU[ii] *= 2.e0;
        ofmo_dscale( nao2, 2.e0, dU );

	ofmo_acc_thread_timer( cid_4c, mythread );

    //printf("%d/%d pot_env3_start\n", monomer_list[0], iscc);
    //fflush(stdout);

	/* 3-centered inter-fragment Coulomb term */
	ofmo_start_thread_timer( cid_3c, mythread );
    //if(fp_prof) {fprintf(fp_prof,"calc_fragment ifc3c\n");fflush(fp_prof);}
	memset( dUtmp, '\0', sizeof(double)*nao2 );
	for ( j=0,jn=0; j<nifc3c; j++ ) {
	    jfrag = joblist_ifc3c[j];
	    ofmo_integ_ifc3c_sorted_partial( nworkers, workerid,
		    maxlqn, shel_atm, shel_ini, atom_x, atom_y, atom_z,
		    lcs_pair, csp_ics, csp_jcs, csp_lps_pair,
		    psp_zeta, psp_dkps, psp_xiza,

		    mlcs[jfrag],
		    mshel_tem[jfrag], mshel_atm[jfrag],
		    mshel_add[jfrag], mshel_ini[jfrag],
		    matom_x[jfrag], matom_y[jfrag], matom_z[jfrag],
		    mprim_exp[jfrag], mprim_coe[jfrag],
		    aopop_mon[jfrag], dUtmp );
	}
	acc_array( nao2, dUtmp, dU );
	ofmo_acc_thread_timer( cid_3c, mythread );
	// debug
	//for ( int ii=0; ii<nao2; ii++ ) dU[ii] = 0.e0;

    //printf("%d/%d pot_env2_start\n", monomer_list[0], iscc);
    //fflush(stdout);
	/* 2-centered inter-fragment Coulomb term */
	ofmo_start_thread_timer( cid_2c, mythread );
    //if(fp_prof) {fprintf(fp_prof,"calc_fragment ifc2c\n");fflush(fp_prof);}
	memset( dUtmp, '\0', sizeof(double)*nao2 );
	for ( j=myrank, jn=0; j<nfrag; j+=nprocs ) {
            if (jn++>NPARTIAL) {
              acc_array( nao2, dUtmp, dU );
              memset( dUtmp, '\0', sizeof(double)*nao2 );
              jn=0;
            }
	    jfrag = frag_order[j];
	    if ( is_in_fragment( jfrag, nmonomer, monomer_list) ) continue;
	    //nat_mon = nfatom[jfrag];
	    switch( ofmo_get_approx_level(jfrag) ) {
		case OFMO_IFC4C:
		case OFMO_IFC3C:
		    for ( iat=0; iat<nfatom[jfrag]; iat++ )
			atpop_mt[iat] = (double)matomic_number[jfrag][iat];
		    break;
		case OFMO_IFC2C:
		    for ( iat=0; iat<nfatom[jfrag]; iat++ )
			atpop_mt[iat] = (double)matomic_number[jfrag][iat]
			    - atpop_mon[jfrag][iat];
	    }
	    ofmo_integ_ifc2c_sorted_partial( nthreads, mythread, maxlqn,
		    flcs, shel_tem, shel_atm, shel_add, shel_ini,
		    atom_x, atom_y, atom_z, prim_exp, prim_coe,
		    nfatom[jfrag], matom_x[jfrag], matom_y[jfrag],
		    matom_z[jfrag], atpop_mt, dUtmp );
		    //matom_z[jfrag], atpop, _WORK_ );
	}
	acc_array( nao2, dUtmp, dU );
	// accumulate env.pot in local process
#pragma omp critical
	acc_array( nao2, dU, U );

	ofmo_acc_thread_timer( cid_2c, mythread );
    }	// #pragma omp parallel
    ofmo_scale_diag( nao, 0.5e0, U );
    ofmo_acc_proc_timer( tid_integ );
    
    ofmo_start_proc_timer( tid_integ );

    //printf("%d/%d pot_env_done\n", monomer_list[0], iscc);
    //fflush(stdout);

    //if(fp_prof) {fprintf(fp_prof,"calc_fragment reduction\n");fflush(fp_prof);}
    /* reduction */
    ofmo_start_proc_timer( tid_comm );
    double *WORK = _WORK_;
    MPI_Allreduce( U, WORK, nao2, MPI_DOUBLE, MPI_SUM, comm );
    memcpy( U, WORK, sizeof(double)*nao2 );
    MPI_Bcast( S, nao2, MPI_DOUBLE, 0, comm );
    MPI_Bcast( H, nao2, MPI_DOUBLE, 0, comm );

    acc_array(nao2, U, H);
    ofmo_acc_proc_timer( tid_comm );
    ofmo_acc_proc_timer( tid_Integ );

    // profile
    if ( fp_prof ) {
	if ( nmonomer == 1 )
	    fprintf(fp_prof,
		    "**M frag= %4d nat= %3d ncs= %3d nao= %3d"
		    " nifc4c= %3d nifc3c= %3d\n",
		    monomer_list[0], nat, ncs, nao, nifc4c, nifc3c );
	else if ( nmonomer == 2 ) {
	    fprintf(fp_prof,
		    "**D frag= %4d %4d nat= %3d ncs= %3d nao= %3d"
		    " nifc4c= %3d nifc3c= %3d\n",
		    monomer_list[0], monomer_list[1], nat, ncs, nao,
		    nifc4c, nifc3c );
	}
        fflush(fp_prof);
    }

    /* SCF calculation */
    ofmo_start_proc_timer( tid_SCF );
    //if(fp_prof) {fprintf(fp_prof,"calc_fragment scf_rhf\n");fflush(fp_prof);}
    /*
    double scfd=tolscf, scfe;
    if      ( scfd <= 1.e-4 ) scfe = scfd * 1.e-3;
    else if ( scfd <= 1.e-3 ) scfe = scfd * 1.e-2;
    else                      scfe = scfd * 1.e-1;
    */
    ofmo_twoint_eps_ps4(eps_ps4);
    ofmo_twoint_eps_eri(eps_eri);
    ofmo_twoint_eps_sch(eps_sch);

    double * mo_tei = NULL;

    if ( level == OFMO_VQE ){
        const int nao_4 = nao * nao * nao * nao;
        mo_tei = (double *)malloc(sizeof(double) * nao_4);
    }

    if ( level == OFMO_RHF || level == OFMO_VQE ) {
        ofmo_scf_set_convType((nmonomer==1)? scc: scf); // Should be in args.
	ofmo_scf_rhf(
		comm, maxlqn, Enuc, ncs, nao,
		flcs, shel_atm, shel_ini, atom_x, atom_y, atom_z,
		lcs_pair, csp_schwarz, csp_ics, csp_jcs,
		csp_lps_pair, psp_zeta, psp_dkps, psp_xiza,
		nat, nocc, S, H, maxscf, scfe, scfd,
		D, C, mo_tei, ev, energy );
    } else {
	if ( fp_prof )
	    fdbg( fp_prof, "ERROR: method %d is not supported\n", level );
	return -1;
    }

    double _dv;
    if(level == OFMO_VQE){
        //double * X = (double *)malloc(sizeof(double)*nao*nao); //Orth mat
        double * H_MO = (double *)malloc(sizeof(double)*nao2); // MO H
        double * U_MO = (double *)malloc(sizeof(double)*nao2); // MO basis U (environmental potential)
        //int orth_ret = ofmo_symm_orth(nao, S, C, X);
        char * desc;
        int monhomo, monlumo, dimhomo, dimlumo, monent, diment;
        ofmo_data_get_vals("desc monhomo monlumo dimhomo dimlumo monent diment",
            &desc, &monhomo, &monlumo, &dimhomo, &dimlumo, &monent, &diment);

        /*if(orth_ret == 0){
            ofmo_orth_C(nao, X, C);
        }else if(orth_ret < 0){
            printf("ERROR orth.\n");
            return -1;
        }*/
        ofmo_ao2mo_H(nao, H, C, H_MO); // Diagonal components mult by 2.
        ofmo_ao2mo_H(nao, U, C, U_MO); // Diagonal components mult by 2.
        int mythread = omp_get_thread_num();
        if(nmonomer == 1)
            ierr = ofmo_vqe_call(myrank, nmonomer, monomer_list, nao, H_MO, U_MO, mo_tei, S, C, nelec, Enuc, *energy, iscc, ev, desc,
                monhomo, monlumo, monent );
        else if(nmonomer == 2)
            ierr = ofmo_vqe_call(myrank, nmonomer, monomer_list, nao, H_MO, U_MO, mo_tei, S, C, nelec, Enuc, *energy, iscc, ev, desc,
                dimhomo, dimlumo, diment );
    	if ( ierr != 0 ){
            if(nmonomer == 1)
                printf("%d/%d Err in VQE call.\n", monomer_list[0], iscc);
            else
                printf("%d-%d/%d Err in VQE call.\n", monomer_list[0], monomer_list[1], iscc);
            fflush(stdout);
            return -1;
        }
        ierr = ofmo_vqe_get_energy(nmonomer, monomer_list, iscc, energy, &_dv, desc);
    	if ( ierr != 0 ){
            if(nmonomer == 1)
                printf("%d/%d Err in Energy acq.\n", monomer_list[0], iscc);
            else
                printf("%d-%d/%d Err Energy acq.\n", monomer_list[0], monomer_list[1], iscc);
            fflush(stdout);
            return -1;
        }

        free(H_MO);
        free(U_MO);

        if(nmonomer == 1){
            double * amps;
            char ** fock;
            int namps, i, j, ij;
            double * corr_mat;

            //double * hf_D = (double *) malloc (sizeof(double) * nao2);
            //memcpy(hf_D, D, sizeof(double) * nao2);

            ierr = ofmo_vqe_get_amplitudes(monomer_list[0], iscc, nao, &namps, &amps, &fock, &corr_mat, desc);
            if (ierr != 0) return -1;
            ofmo_vqe_posthf_density(namps, amps, fock, corr_mat, C, nao, D);
            free(amps);
            for(i=0; i<namps; i++){
                free(fock[i]);
            }
            free(fock);
            free(corr_mat);
/*
            printf("===== D _ HF =====\n");
            for(i=0, ij=0; i<nao; i++){
                for(j=0; j<=i; j++, ij++){
                    printf("%f\t", hf_D[ij]);
                }
                printf("\n");
            }
            free(hf_D);

            printf("===== D =====\n");
            for(i=0, ij=0; i<nao; i++){
                for(j=0; j<=i; j++, ij++){
                    printf("%f\t", D[ij]);
                }
                printf("\n");
            }
*/

        }
    }

    if(mo_tei) free(mo_tei);

    /* energies */
    double dv;
    if(level == OFMO_RHF) dv = 4.0e0*ofmo_dot_product( nao2, D, U );
    else if(level == OFMO_VQE) dv = _dv;
    else{
        printf("ERROR : invalid method\n");
        fflush(stdout);
    }
    *energy0 = *energy - dv;
    if(nmonomer == 1){
        printf("it=%d\tmon=[%d]\tenergy=%f\tenergy_0=%f\n", iscc, monomer_list[0], *energy, *energy0);
    }else if(nmonomer == 2){
        printf("it=%d\tmon=[%d, %d]\tenergy=%f\tenergy_0=%f\n", iscc, monomer_list[0], monomer_list[1], *energy, *energy0);
    }


    /* Mulliken population */
    double *aopop = _daopop_, *atpop = _datpop_;
    if ( nmonomer == 1 ) {
	int ifrag = monomer_list[0];
	ofmo_scf_mulliken_population(
		nat, nao, maxlqn, flcs, shel_atm, shel_ini, S, D,
		aopop, atpop );
	if ( is_root ) {
	    ofmo_put_monomer_density( ifrag, D );
	    ofmo_put_monomer_aopop( ifrag, aopop );
	    ofmo_put_monomer_atpop( ifrag, atpop );
	    ofmo_put_monomer_energy( ifrag, energy );
	    ofmo_put_monomer_energy0( ifrag, energy0 );
	}
	if ( fp_prof ) {
	    fprintf(fp_prof, "# E= %16.10f E-DV= %16.10f DV= %16.10f\n",
		    *energy, *energy0, dv);
	    fflush( fp_prof );
	}
    } else if (nmonomer == 2) {
	int i;
	for ( i=0; i<nao2; i++ ) D0[i]-=D[i];
	ofmo_scf_mulliken_population(
		nat, nao, maxlqn, flcs, shel_atm, shel_ini, S, D0,
		aopop, atpop );
	for ( i=0; i<nat; i++ ) atpop[i] *= -1.e0;
	for ( i=0; i<nao; i++ ) aopop[i] *= -1.e0;
	*ddv = -4.0e0 * ofmo_dot_product( nao2, D0, U );
    printf("# E= %16.10f E-DV= %16.10f DV= %16.10f "
        "DDV= %16.10f\n", *energy, *energy0, dv, *ddv);
	if ( fp_prof ) {
	    fprintf( fp_prof,
		    "# E= %16.10f E-DV= %16.10f DV= %16.10f "
		    "DDV= %16.10f\n", *energy, *energy0, dv, *ddv);
	    fflush( fp_prof );
	}
	// ofmo-worker-mainに渡す変数のコピー（rootランクだけでもよい）
	*fnao = nao;
	*fnat = nat;
	memcpy( daopop, aopop, sizeof(double) * nao );
	memcpy( sao2tao, fsao2tuao, sizeof(int) * nao );
	memcpy( datpop, atpop, sizeof(double) * nat );
	memcpy( fatom2tatom, fat2tat, sizeof(int) * nat );
    }

    ofmo_acc_proc_timer( tid_SCF );
    ofmo_acc_proc_timer( tid_Total );
    // timer
    ofmo_show_thread_timer_all();
    ofmo_show_proc_timer_all();
    return 0;
}

static double ofmo_calc_enucij(
	const int nati, const int atni[],
	const double xi[], const double yi[], const double zi[],
	const int natj, const int atnj[],
	const double xj[], const double yj[], const double zj[] ) {
    int iat, jat;
    double A[3], AB[3], AB2;
    double qi, qj, enuc = 0.e0;
    for ( iat=0; iat<nati; iat++ ) {
	qi = (double)atni[iat];
	A[0] = xi[iat];
	A[1] = yi[iat];
	A[2] = zi[iat];
	for ( jat=0; jat<natj; jat++ ) {
	    qj = (double)atnj[jat];
	    AB[0] = A[0] - xj[jat];
	    AB[1] = A[1] - yj[jat];
	    AB[2] = A[2] - zj[jat];
	    AB2 = AB[0]*AB[0] + AB[1]*AB[1] + AB[2]*AB[2];
	    if ( AB2 < 1.e-8 ) continue;
	    enuc += qi*qj * sqrt(1.e0/AB2);
	}
    }
    return enuc;
}

/** A function that performs all approximate dimer calculations
 *
 * いくつかのダイマーをまとめて
 * 計算するようにしている。
 *
 * @param[out] *total_dimer_es_energy ES dimer energy sum
 *
 * @ingroup ofmo-calc
 * */
int ofmo_calc_es_dimer( MPI_Comm comm, int njob, int joblist[],
	double *energy0 ) {
    static int nfrag, *nfao, maxnfao, maxnfatom, maxnfcs;
    static int maxlqn, **mlcs;
    static int **mshel_tem, **mshel_atm, **mshel_add, **mshel_ini;
    static double **matom_x, **matom_y, **matom_z;
    static double **mprim_exp, **mprim_coe;
    static int *nfatom, **matomic_number, *nfcs, *icharge;
    static int **msao2tuao, natom, nao, **ifatom, dold;
    static int maxnpspair, nbody;
    static double ldimer;
    static int type4c[6*6];
    static int called = false;
#if 1
#ifdef USE_CUDA
    static int first = true;
    if (first) {
      cuda_print_wifc4c();
    }
#endif
#endif
    //MPI_Status status;
    //int tag = 101;
    if ( !called ) {
	int ierr, maxlqn2;
	ierr = ofmo_data_get_vals(
		"nfrag nfao maxnfao maxnfcs maxlqn maxnfatom "
		"mlcs mshel_tem mshel_atm mshel_add mshel_ini "
		"mprim_exp mprim_coe icharg "
		"matx maty matz matn nfatom "
		"msao2tuao dold nao natom ifatom "
		"nfcs maxnpspair nbody "
		"ldim",
		&nfrag, &nfao, &maxnfao, &maxnfcs, &maxlqn, &maxnfatom,
		&mlcs, &mshel_tem, &mshel_atm, &mshel_add,
		&mshel_ini,
		&mprim_exp, &mprim_coe, &icharge,
		&matom_x, &matom_y, &matom_z, &matomic_number, &nfatom,
		&msao2tuao, &dold, &nao, &natom, &ifatom,
		&nfcs, &maxnpspair, &nbody,
		&ldimer);
	if ( ierr != 0 ) return -1;
	alloc();
	maxlqn2 = ((maxlqn+1)*(maxlqn+2))>>1;
	if ( maxlqn == 0 ) {
	    type4c[0] = 0;
	} else if ( maxlqn <= 2 ) {
	    int Lab, Lcd, ix;
	    ix=0;
	    for ( Lab=0; Lab<maxlqn2; Lab++ ) {
		for ( Lcd=0; Lcd<maxlqn2; Lcd++ ) {
		    type4c[ix] = Lab*6 + Lcd;
		    ix++;
		}
	    }
	} else return -1;
	called = true;
    }
    // for profile
    int cid_init, cid_integ;
    int tid_init, tid_integ, tid_intred, tid_all;
    cid_init   = ofmo_create_thread_timer( "INIT", 0 );
    cid_integ  = ofmo_create_thread_timer( "INTEG", 0 );

    tid_init   = ofmo_create_proc_timer( "init", 0 );
    tid_integ  = ofmo_create_proc_timer( "integ", 0 );
    tid_intred = ofmo_create_proc_timer( "integ+reduce", 0 );
    tid_all    = ofmo_create_proc_timer( "all", 0 );

    ofmo_start_proc_timer( tid_init );
    ofmo_start_proc_timer( tid_all );

    //int count;
    double t0[MAXNJOB*3], t[MAXNJOB*3];
    double *dE, *UjDi, *UiDj, Etot[MAXNJOB];
    double enucij[MAXNJOB];
    double **Di, **Dj;
    double energy_es_dimer = 0.e0;
    dE   = &t0[MAXNJOB*0];
    UjDi = &t0[MAXNJOB*1];
    UiDj = &t0[MAXNJOB*2];
    /* for cutoff table */
    /* ifrags */
    int **lcs_pair_i;
    int **csp_ics_i, **csp_jcs_i, **csp_lps_pair_i;
    double **csp_schwarz_i;
    double **psp_zeta_i, **psp_dkps_i, **psp_xiza_i;
    /* jfrags */
    int **lcs_pair_j;
    int **csp_ics_j, **csp_jcs_j, **csp_lps_pair_j;
    double **csp_schwarz_j;
    double **psp_zeta_j, **psp_dkps_j, **psp_xiza_j;
    lcs_pair_i     = &lcs_pair_mon[0];
    csp_ics_i      = &csp_ics_mon[0];
    csp_jcs_i      = &csp_jcs_mon[0];
    csp_lps_pair_i = &csp_lps_pair_mon[0];
    csp_schwarz_i  = &csp_schwarz_mon[0];
    psp_zeta_i     = &psp_zeta_mon[0];
    psp_dkps_i     = &psp_dkps_mon[0];
    psp_xiza_i     = &psp_xiza_mon[0];

    lcs_pair_j     = &lcs_pair_mon[MAXNJOB];
    csp_ics_j      = &csp_ics_mon[MAXNJOB];
    csp_jcs_j      = &csp_jcs_mon[MAXNJOB];
    csp_lps_pair_j = &csp_lps_pair_mon[MAXNJOB];
    csp_schwarz_j  = &csp_schwarz_mon[MAXNJOB];
    psp_zeta_j     = &psp_zeta_mon[MAXNJOB];
    psp_dkps_j     = &psp_dkps_mon[MAXNJOB];
    psp_xiza_j     = &psp_xiza_mon[MAXNJOB];
    // monomer density matrices
    Di = &_Dmon_[0];
    Dj = &_Dmon_[MAXNJOB];

    /* communicator */
    int nprocs, myrank, root=0, is_root;
    MPI_Comm_size( comm, &nprocs );
    MPI_Comm_rank( comm, &myrank );
    is_root = ( myrank == root );

    /* */
    if ( njob > MAXNJOB ) {
	if ( fp_prof )
	    fdbg( fp_prof, "Illegal number of jobs (%d)\n", njob );
	MPI_Abort( MPI_COMM_WORLD, 1 );
    } else if ( fp_prof ) {
	fprintf( fp_prof, "#ES njob= %2d\n", njob );
	fflush( fp_prof );
    }

    // Initialization
    for ( int ii=0; ii<njob; ii++ )
	dE[ii] = UiDj[ii] = UjDi[ii] = enucij[ii] = 0.e0;
    *energy0 = 0.e0;

    // 動的負荷分散
    // 各積分タイプの計算を128分割する
    /*int njob_ifc4c, uifc4c, maxlqn2;
    maxlqn2 = ((maxlqn+1)*(maxlqn+2))>>1;
    uifc4c = maxlqn2 * maxlqn2;
    njob_ifc4c = (njob*uifc4c)<<7;	// << 7 means "*128"
    ofmo_mt_init( comm );
    ofmo_gc_init( 2, comm, 0, njob_ifc4c );*/
    //ofmo_gc_set_debug_mode();
    // カットオフテーブル計算、および、密度行列データ取得
    int icut, ncut;
    ncut = njob * 2;	// カットオフテーブル計算回数
    icut = 0;
    ofmo_acc_proc_timer( tid_init );
    ofmo_start_proc_timer( tid_integ );
    ofmo_start_proc_timer( tid_intred );
#ifndef PARA_SUB
#pragma omp parallel
#else
    int nnewjob;
    int *newjoblist;
    int amyrank, anprocs;
    int color;
    double *wjob;
#pragma omp parallel shared(nnewjob, newjoblist, amyrank, anprocs, color, wjob)
#endif
    {
	int mythread, nthreads;
	int ifrag, jfrag, ni, nj, ni2, nj2;
	int ijob, mycut, is_odd;
	double *Utmp, *charge;
	double dE_tmp;
	// for dynamic load-balancing
	//int offset;
	int nworkers, workerid, vnworkers, vworkerid;
	int k, mytype, mypos, base;

	mythread = omp_get_thread_num();
	nthreads = omp_get_num_threads();

	ofmo_start_thread_timer( cid_init, mythread );
	Utmp     = dUmaster[mythread];
	charge   = _atpop_local_[mythread];
	//
	nworkers = ofmo_mt_get_nworkers();
	workerid = ofmo_mt_get_workerid( mythread );
	vnworkers = nworkers - 1;
	vworkerid = workerid - 1;
	// Acquisition of monomer density matrix data in master thread
	if ( mythread == 0 ) {
	    for ( ijob=0; ijob<njob; ijob++ ) {
		ifrag = joblist[ijob*2+0];
		jfrag = joblist[ijob*2+1];
		ni    = nfao[ifrag];
		nj    = nfao[jfrag];
		ni2   = (ni*ni+ni)>>1;
		nj2   = (nj*nj+nj)>>1;
		ofmo_get_monomer_density( comm, ifrag, Di[ijob] );
		ofmo_get_monomer_density( comm, jfrag, Dj[ijob] );
	    }
	}
#pragma omp critical
	{
	    mycut = icut;
	    icut++;
	}
	while ( mycut < ncut ) {
	    ifrag  = joblist[mycut];
	    ijob   = mycut>>1;
	    is_odd = ( mycut & 0x01 );
	    if ( is_odd ) {
		/* make cutoff table j */
		ofmo_cutoff_make_table( maxlqn, mlcs[ifrag], 
			mshel_tem[ifrag], mshel_atm[ifrag],
			mshel_add[ifrag],
			matom_x[ifrag], matom_y[ifrag],
			matom_z[ifrag],
			mprim_exp[ifrag], mprim_coe[ifrag],
			lcs_pair_j[ijob],
			csp_schwarz_j[ijob], csp_ics_j[ijob],
			csp_jcs_j[ijob], csp_lps_pair_j[ijob],
			psp_zeta_j[ijob], psp_dkps_j[ijob],
			psp_xiza_j[ijob] );
	    } else {
		/* make cutoff table i */
		ofmo_cutoff_make_table( maxlqn, mlcs[ifrag], 
			mshel_tem[ifrag], mshel_atm[ifrag],
			mshel_add[ifrag],
			matom_x[ifrag], matom_y[ifrag],
			matom_z[ifrag],
			mprim_exp[ifrag], mprim_coe[ifrag],
			lcs_pair_i[ijob],
			csp_schwarz_i[ijob], csp_ics_i[ijob],
			csp_jcs_i[ijob], csp_lps_pair_i[ijob],
			psp_zeta_i[ijob], psp_dkps_i[ijob],
			psp_xiza_i[ijob] );
	    }
#pragma omp critical
	    {
		mycut = icut;
		icut++;
	    }
	}	// while ( mycut < ncut );
	// ---- ここまでで、カットオフテーブル計算、密度行列データ
	// ---- 取得が終わっている
	ofmo_acc_thread_timer( cid_init, mythread );
#pragma omp barrier
	ofmo_start_thread_timer( cid_integ, mythread );
	// カウンタマスタースレッド以外のすべてのスレッドで計算を行う
	if ( workerid != 0 ) {
	    // まず、２中心クーロン相互作用計算を静的負荷分散で行う
	    for ( ijob=0; ijob<njob; ijob++ ) {
		ifrag = joblist[ijob*2+0];
		jfrag = joblist[ijob*2+1];
		ni    = nfao[ifrag];
		nj    = nfao[jfrag];
		ni2   = (ni*ni+ni)>>1;
		nj2   = (nj*nj+nj)>>1;
		/* UjDi */
		for ( k=0; k<nfatom[jfrag]; k++ )
		    charge[k] = (double)matomic_number[jfrag][k];
		memset( Utmp, '\0', sizeof(double)*ni2 );
		ofmo_integ_ifc2c_sorted_partial( vnworkers, vworkerid,
			maxlqn, mlcs[ifrag], mshel_tem[ifrag],
			mshel_atm[ifrag], mshel_add[ifrag],
			mshel_ini[ifrag],
			matom_x[ifrag], matom_y[ifrag], matom_z[ifrag],
			mprim_exp[ifrag], mprim_coe[ifrag],
			nfatom[jfrag],
			matom_x[jfrag], matom_y[jfrag], matom_z[jfrag],
			charge, Utmp );
		ofmo_scale_diag( ni, 0.5e0, Utmp );
		dE_tmp = 4.e0 * ofmo_dot_product( ni2, Di[ijob], Utmp );
#ifdef FJ_COMP
#pragma omp critical
		UjDi[ijob] += dE_tmp;
#else
#pragma omp atomic
		UjDi[ijob] += dE_tmp;
#endif
		/* UiDj */
		for ( k=0; k<nfatom[ifrag]; k++ )
		    charge[k] = (double)matomic_number[ifrag][k];
		memset( Utmp, '\0', sizeof(double)*nj2 );
		ofmo_integ_ifc2c_sorted_partial( vnworkers, vworkerid,
			maxlqn, mlcs[jfrag], mshel_tem[jfrag],
			mshel_atm[jfrag], mshel_add[jfrag],
			mshel_ini[jfrag],
			matom_x[jfrag], matom_y[jfrag], matom_z[jfrag],
			mprim_exp[jfrag], mprim_coe[jfrag],
			nfatom[ifrag],
			matom_x[ifrag], matom_y[ifrag], matom_z[ifrag],
			charge, Utmp );
		ofmo_scale_diag( nj, 0.5e0, Utmp );
		dE_tmp = 4.e0 * ofmo_dot_product( nj2, Dj[ijob], Utmp );
#ifdef FJ_COMP
#pragma omp critical
		UiDj[ijob] += dE_tmp;
#else
#pragma omp atomic
		UiDj[ijob] += dE_tmp;
#endif
	    }	// for ( ijob=0; ;  )
	} else {	// if ( workerid != 0 )
	    // カウンタマスタースレッドでは核間反発エネルギーを計算
	    for ( ijob=0; ijob<njob; ijob++ ) {
		ifrag = joblist[ijob*2+0];
		jfrag = joblist[ijob*2+1];
		enucij[ijob] = ofmo_calc_enucij(
			nfatom[ifrag], matomic_number[ifrag],
			matom_x[ifrag], matom_y[ifrag],
			matom_z[ifrag],
			nfatom[jfrag], matomic_number[jfrag],
			matom_x[jfrag], matom_y[jfrag],
			matom_z[jfrag] );
	    }
	}	// if ( workerid != 0 )

	// 次に、4中心クーロン相互作用項を静的負荷分散で計算
	// for control load-balancing
	ofmo_integ_set_loop_offset( mythread, 0 );
#ifdef PARA_SUB
#pragma omp master
    {
      newjoblist=(int*)malloc(njob*sizeof(int));
      wjob=(double *)malloc(njob*sizeof(double)+njob*2*sizeof(int));
    }
#pragma omp barrier
#pragma omp for
      for ( ijob=0; ijob<njob; ijob++ ) {
        int Lab = maxlqn*(maxlqn+1)/2 + maxlqn;
        int ncspair_i = (lcs_pair_i[ijob])[Lab+1];
        int ncspair_j = (lcs_pair_j[ijob])[Lab+1];
        wjob[ijob]=(double)ncspair_i*ncspair_j;
#if 1
        double eps_ps4 = ofmo_twoint_eps_ps4(0);
        size_t nps4=0;
        for (int ii=0;ii<ncspair_i;ii++) {
          for (int jj=0;jj<ncspair_j;jj++)
            if ((csp_schwarz_i[ijob])[ii]*(csp_schwarz_j[ijob])[jj]>=eps_ps4) nps4++;
        }
        wjob[ijob]=(double)nps4;
#endif
      }
#pragma omp master
    {
      int *atask=(int *)(wjob+njob);
      int *res=atask+njob;
      int natask=ofmo_aggregateTask(nprocs, njob, wjob, atask);
      //if (myrank==0) printf("%d %d -> %d %d\n",nprocs, njob, nprocs, natask);
      ofmo_assignRes(nprocs, natask, wjob);
      for (int j=0; j<natask; j++ ) res[j]=(int)wjob[j];
      color=ofmo_getMyColorAndRankFromAssignedRes(myrank,res,&amyrank);
      anprocs = res[color];
      nnewjob=0;
      for (int j=0; j<njob; j++) if (atask[j]==color) newjoblist[nnewjob++]=j;
      free(wjob);
    }
#ifdef USE_CUDA
      cuda_Reconfig(amyrank, anprocs, comm);
#endif
#pragma omp barrier
    ofmo_integ_set_loop_offset( mythread, 0 );

#endif /* PARA_SUB */
#ifdef USE_CUDA
        float *csp_schwarz_f;
        {
          int nat_i=0, ncs_i=0, nao_i=0, ncspair_i=0, npspair_i=0;
          int nat_j=0, ncs_j=0, nao_j=0, ncspair_j=0, npspair_j=0;
          int Lab = maxlqn*(maxlqn+1)/2 + maxlqn;
          int max_num_klcs = 0;
#ifndef PARA_SUB
          for ( ijob=0; ijob<njob; ijob++ ) {
#else
          for (int aj=0; aj<nnewjob; aj++ ) {
            ijob = newjoblist[aj];
#endif
            ifrag = joblist[ijob*2+0];
            jfrag = joblist[ijob*2+1];
            nat_i = MAX2(nat_i, nfatom[ifrag]);
            nat_j = MAX2(nat_j, nfatom[jfrag]);
            ncs_i = MAX2(ncs_i, nfcs[ifrag]);
            ncs_j = MAX2(ncs_j, nfcs[jfrag]);
            nao_i = MAX2(nao_i, nfao[ifrag]);
            nao_j = MAX2(nao_j, nfao[jfrag]);
            int ncspair_t;
            ncspair_t = (lcs_pair_i[ijob])[Lab+1];
            ncspair_i = MAX2(ncspair_i, ncspair_t);
            npspair_i = MAX2(npspair_i, (csp_lps_pair_i[ijob])[ncspair_t]);
            ncspair_t = (lcs_pair_j[ijob])[Lab+1];
            ncspair_j = MAX2(ncspair_j, ncspair_t);
            npspair_j = MAX2(npspair_j, (csp_lps_pair_j[ijob])[ncspair_t]);
            int max_num_klcs_t = cuda_max_num_klcs(maxlqn, lcs_pair_j[ijob]);
            max_num_klcs = MAX2(max_num_klcs, max_num_klcs_t);
          }

          cuda_ifc4c_Init(maxlqn, max_num_klcs,
              nat_i, ncs_i, nao_i, ncspair_i, npspair_i,
              nat_j, ncs_j, nao_j, ncspair_j, npspair_j);

          csp_schwarz_f = (float*)malloc(sizeof(float)*MAX2(ncspair_i,ncspair_j));
        }
#endif
	//while ( (k=ofmo_gc_nxtval(2)) < njob_ifc4c ) {
#ifndef PARA_SUB
	for ( ijob=0; ijob<njob; ijob++ ) {
#else
        for (int aj=0; aj<nnewjob; aj++ ) {
            ijob = newjoblist[aj];
#endif
	    ifrag = joblist[ijob*2+0];
	    jfrag = joblist[ijob*2+1];
	    ni    = nfao[ifrag];
	    nj    = nfao[jfrag];
	    ni2   = ni*(ni+1)/2;
	    nj2   = nj*(nj+1)/2;
	    memset( Utmp, '\0', sizeof(double)*ni2 );

#ifdef USE_CUDA
            int Lab = maxlqn*(maxlqn+1)/2 + maxlqn;
            int nat_i = nfatom[ifrag];
            int ncs_i = nfcs[ifrag];
            int nao_i = nfao[ifrag];
            int ncspair_i = (lcs_pair_i[ijob])[Lab+1];
            int npspair_i = (csp_lps_pair_i[ijob])[ncspair_i];
            for (int ii=0; ii<ncspair_i; ii++)
              csp_schwarz_f[ii]=(float)(csp_schwarz_i[ijob])[ii];
            cuda_ifc4c_SetData(0,
                maxlqn, nat_i, ncs_i, nao_i,
                ncspair_i, npspair_i,
                mshel_atm[ifrag], mshel_ini[ifrag],
                matom_x[ifrag], matom_y[ifrag], matom_z[ifrag],
                lcs_pair_i[ijob], csp_lps_pair_i[ijob],
                csp_ics_i[ijob], csp_jcs_i[ijob],
                psp_zeta_i[ijob], psp_dkps_i[ijob], psp_xiza_i[ijob],
                csp_schwarz_f, NULL);
            int nat_j = nfatom[jfrag];
            int ncs_j = nfcs[jfrag];
            int nao_j = nfao[jfrag];
            int ncspair_j = (lcs_pair_j[ijob])[Lab+1];
            int npspair_j = (csp_lps_pair_j[ijob])[ncspair_j];
            for (int ii=0; ii<ncspair_j; ii++)
              csp_schwarz_f[ii]=(float)(csp_schwarz_j[ijob])[ii];
            cuda_ifc4c_SetData(1,
                maxlqn, nat_j, ncs_j, nao_j,
                ncspair_j, npspair_j,
                mshel_atm[jfrag], mshel_ini[jfrag],
                matom_x[jfrag], matom_y[jfrag], matom_z[jfrag],
                lcs_pair_j[ijob], csp_lps_pair_j[ijob],
                csp_ics_j[ijob], csp_jcs_j[ijob],
                psp_zeta_j[ijob], psp_dkps_j[ijob], psp_xiza_j[ijob],
                csp_schwarz_f, Dj[ijob]);
#endif

#ifndef PARA_SUB
	    ofmo_integ_ifc4c_sorted_partial( nworkers, workerid,
#else
            int anworkers = anprocs * omp_get_num_threads();
            int aworkerid = amyrank * omp_get_num_threads() + omp_get_thread_num();
	    ofmo_integ_ifc4c_sorted_partial( anworkers, aworkerid,
#endif
		    maxlqn,
		    mshel_atm[ifrag], mshel_ini[ifrag],
		    matom_x[ifrag], matom_y[ifrag], matom_z[ifrag],
		    lcs_pair_i[ijob],
		    csp_schwarz_i[ijob], csp_ics_i[ijob],
		    csp_jcs_i[ijob], csp_lps_pair_i[ijob],
		    psp_zeta_i[ijob], psp_dkps_i[ijob],
		    psp_xiza_i[ijob],

		    mshel_atm[jfrag], mshel_ini[jfrag],
		    matom_x[jfrag], matom_y[jfrag], matom_z[jfrag],
		    lcs_pair_j[ijob],
		    csp_schwarz_j[ijob], csp_ics_j[ijob],
		    csp_jcs_j[ijob], csp_lps_pair_j[ijob],
		    psp_zeta_j[ijob], psp_dkps_j[ijob],
		    psp_xiza_j[ijob],

                    mlcs[jfrag],
		    nj, Dj[ijob], Utmp );
#ifdef USE_CUDA
            cuda_ifc4c_GetVfrg(ni, Utmp);
#endif

	    ofmo_scale_diag( ni, 0.5e0, Utmp );
	    dE_tmp = 8.e0 * ofmo_dot_product( ni2, Di[ijob], Utmp );
#ifdef FJ_COMP
#pragma omp critical
	    dE[ijob] += dE_tmp;
#else
#pragma omp atomic
	    dE[ijob] += dE_tmp;
#endif
	}	// for ( ijob=0 );
#ifdef USE_CUDA
        Free(csp_schwarz_f);
        cuda_ifc4c_Finalize();
#endif
#ifdef PARA_SUB
#pragma omp barrier
#pragma omp master
        free(newjoblist);
#ifdef USE_CUDA
        cuda_Reconfig(myrank, nprocs, comm);
#endif
        ofmo_integ_set_loop_offset( mythread, 0 );
#endif /* PARA_SUB */
	ofmo_acc_thread_timer( cid_integ, mythread );
    }	/* pragma omp parallel */
    ofmo_acc_proc_timer( tid_integ );
    MPI_Reduce( t0, t, MAXNJOB*3, MPI_DOUBLE, MPI_SUM, root, comm ); 
    memcpy( t0, t, sizeof(double)*MAXNJOB*3 );
    ofmo_acc_proc_timer( tid_intred );

    if ( is_root ) {
	for ( int i=0; i<njob; i++ ) {
	    fprintf(fp_prof, "#ES %3d - %3d dE, UjDi, UiDj, enucij ="
		    " %16.10f %16.10f %16.10f %16.10f\n", joblist[i*2+0],
		    joblist[i*2+1], dE[i], UjDi[i], UiDj[i], enucij[i] );
	}
	fflush( fp_prof );
        int jn=0;
        double pEtot=0.0e0;
	for ( int i=0; i<njob; i++ ) {
	    Etot[i] = dE[i]+ enucij[i] + UjDi[i] + UiDj[i];
	    //energy_es_dimer += Etot[i];
	    //
	    //*energy0 += Etot[i];
            pEtot += Etot[i];
            if (++jn>NPARTIAL) {
              energy_es_dimer += pEtot;
	      *energy0 += pEtot;
              pEtot=0.0e0;
              jn=0;
            }
	}
        energy_es_dimer += pEtot;
        *energy0 += pEtot;
    }
    ofmo_acc_proc_timer( tid_all );
    ofmo_show_thread_timer_all();
    ofmo_show_proc_timer_all();
    return 0;
}
