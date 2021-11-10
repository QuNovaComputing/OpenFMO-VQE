/** モノマー密度行列のキャッシュ機構を導入した
  */
#include <stdio.h>
#include <time.h>
#include <string.h>

#include "ofmo-def.h"
#include "ofmo-data.h"
#include "ofmo-mserv-cont.h"
#include "ofmo-prof.h"

#define MAXDENCACHE	16

static time_t	*LastUse = NULL;
static int	*CachedFrag = NULL;
static double	**CachedDens = NULL;
static int      _maxdencache_ = 0;

/* profile information */
static int NCALLED = 0;	/* # of calling of ofmo_get_monomer_density */
static int NHITS   = 0; /* # of cache hits in ofmo_get_monomer_density */

/* deallocaion of cache for monomer density matrices */
static void dealloc_cache() {
    Free( CachedDens[0] );
    Free( CachedDens );
    Free( CachedFrag );
    Free( LastUse );
}

void ofmo_show_cache_prof() {
    if ( fp_prof ) {
	double rate;
	rate = (double)NHITS/(double)NCALLED * 100.e0;
	fprintf( fp_prof,
		"CALLED=%6d, HITS=%6d, rate=%7.2f(%%)\n",
		NCALLED, NHITS, rate);
	fflush( fp_prof );
    }
}

/* memory allocation of cache for monomer density matrices */
static int alloc_cache() {
    static int called = false;
    if ( called ) return 0;
    int maxdencache;
    char *str;
    if ( (str=getenv("OFMO_MAXDENCACHE")) != NULL ) {
	maxdencache = atoi( str );
	if ( maxdencache < 0 ) maxdencache = 0;
    } else {
	maxdencache = MAXDENCACHE;
    }
    if ( maxdencache < 1 ) {
	_maxdencache_ = 0;
	if ( fp_prof )
	    fprintf(fp_prof, "== density matrices is not cached\n");
	called = true;
	return 0;
    } else {
	if ( fp_prof )
	    fprintf( fp_prof, "== density cache = %d ==\n", maxdencache );
    }
    int maxnfao, n2;
    if ( ofmo_data_get_vals("maxnfao", &maxnfao) != 0 )
	return -1;
    n2 = maxnfao*(maxnfao+1)/2;
    LastUse    = (time_t*)malloc( sizeof(time_t) * maxdencache );
    CachedFrag = (int*)malloc( sizeof(int) * maxdencache );
    CachedDens = (double**)malloc( sizeof(double*) * maxdencache );
    CachedDens[0] = (double*)malloc( sizeof(double) * n2 * maxdencache );
    if ( LastUse == NULL || CachedFrag == NULL || CachedDens == NULL
	   || CachedDens[0] == NULL ) {
	dbg("error: failure in memory allocation\n");
	return -1;
    }
    for ( int i=1; i<maxdencache; i++ )
	CachedDens[i] = CachedDens[i-1]+n2;
    for ( int i=0; i<maxdencache; i++ ) {
	LastUse[i] = 0;
	CachedFrag[i] = -1;
    }
    // information of memory allocation
    if ( fp_prof ) {
	double dsize;
	dsize = (double)(maxdencache * n2 * sizeof(double)) /
	    (double)(1024*1024);
	fprintf(fp_prof,
		"== allocd memory size in ofmo-monomer-data.c = "
		"%10.3f MB\n", dsize );
    }
    atexit( dealloc_cache );
    _maxdencache_ = maxdencache;
    called = true;
    return 0;
}

/* clear cache */
static void clear_cache() {
    if ( _maxdencache_ < 1 ) return;
    for ( int i=0; i<_maxdencache_; i++ ) CachedFrag[i] = -1;
}

/* density matrix of ifrag is cached (>=0) or not (-1) */
static int is_cached( const int ifrag ) {
    if ( _maxdencache_ < 1 ) return -1;
    for ( int i=0; i<_maxdencache_; i++ ) {
	if ( CachedFrag[i] == ifrag ) return i;
    }
    return -1;
}

/* get cache location of oldest reffered one */
static int is_oldest() {
    int    oldest_loc;
    time_t oldest_time;
    if ( _maxdencache_ < 1 ) return -1;
    oldest_loc  = 0;
    oldest_time = LastUse[0];
    for ( int i=1; i<_maxdencache_; i++ ) {
	if ( LastUse[i] < oldest_time ) {
	    oldest_time = LastUse[i];
	    oldest_loc  = i;
	}
    }
    return oldest_loc;
}

static int dold = OFMO_DENS1;
static int dnew = OFMO_DENS2;
static int aopold = OFMO_AOPOP1;
static int aopnew = OFMO_AOPOP2;
static int atpold = OFMO_ATPOP1;
static int atpnew = OFMO_ATPOP2;

void ofmo_reset_monomer_data_workerp() {
    dold = OFMO_DENS1;
    dnew = OFMO_DENS2;
    aopold = OFMO_AOPOP1;
    aopnew = OFMO_AOPOP2;
    atpold = OFMO_ATPOP1;
    atpnew = OFMO_ATPOP2;
}

/** Function to get the monomer density matrix (group call)
 *
 * A function to get the monomer density matrix.
 * At the end, all processes belonging to the communicator given in
 * the argument get the density matrix data. Since some recently
 * referenced monomer density matrix data is cached inside this function,
 * the cached data is copied without communication at the time of cache
 * hit. If a cache miss occurs, the specified monomer density matrix
 * data is acquired by communication and cached in place of the oldest
 * referenced data.
 *
 * @li \c MAXDENCACHE Maximum number of density matrix data to cache
 *
 * @ingroup ofmo-calc
 *
 * */
int ofmo_get_monomer_density( MPI_Comm comm, const int ifrag,
	double D[] ) {
    static int called = false;
    static int *nfao;
    if ( !called ) {
	int ierr;
	ierr = ofmo_data_get_vals("nfao", &nfao );
	if ( ierr != 0 ) {
	    dbg(" ");
	    return -1;
	}
	alloc_cache();
	called = true;
    }
    int n, n2, id;
    int root=0, myrank;
    MPI_Comm_rank( comm, &myrank );
    n = nfao[ifrag];
    n2 = ((n*n+n)>>1);
    id = is_cached(ifrag);
    if ( id < 0 ) {
	if ( myrank == root ) {
	    id = is_oldest();
	    ofmo_worker_get( dold, ifrag, D );
	}
	MPI_Bcast( &id, 1, MPI_INT, root, comm );
	MPI_Bcast( D, n2, MPI_DOUBLE, root, comm );
	if ( id >= 0 ) {
	    memcpy( CachedDens[id], D, sizeof(double)*n2 );
	    LastUse[id]    = time(NULL);
	    CachedFrag[id] = ifrag;
	}
    } else {
	memcpy( D, CachedDens[id], sizeof(double)*n2 );
	LastUse[id] = time(NULL);
	NHITS++;
    }
    /* for profiling */
    NCALLED++;
    return 0;
}


/*// 呼び出したプロセスのみ、モノマー密度行列データを取得
int ofmo_get_monomer_density( int ifrag, double D[] ) {
    ofmo_worker_get( dold, ifrag, D );
    return 0;
}*/

// 呼び出したプロセスのみ、モノマー密度行列データを更新
int ofmo_put_monomer_density( int ifrag, double D[] ) {
    ofmo_worker_put( dnew, ifrag, D );
    return 0;
}

// 呼び出したプロセスのみ、モノマーAOPOPを取得
int ofmo_get_monomer_aopop( int ifrag, double aopop[] ) {
    ofmo_worker_get( aopold, ifrag, aopop );
    return 0;
}

// 呼び出したプロセスのみ、モノマーAOPOPを更新
int ofmo_put_monomer_aopop( int ifrag, double aopop[] ) {
    ofmo_worker_put( aopnew, ifrag, aopop );
    return 0;
}

// 呼び出したプロセスのみ、モノマーのatomic populationを取得
int ofmo_get_monomer_atpop( int ifrag, double atpop[] ) {
    ofmo_worker_get( atpold, ifrag, atpop );
    return 0;
}

// 呼び出したプロセスのみ、モノマーのatomic populationを更新
int ofmo_put_monomer_atpop( int ifrag, double atpop[] ) {
    ofmo_worker_put( atpnew, ifrag, atpop );
    return 0;
}


// データの更新（データIDの交換）を行う
int ofmo_update_monomer_data() {
    int itmp;
    itmp = dold;     dold = dnew;     dnew = itmp;
    itmp = aopold; aopold = aopnew; aopnew = itmp;
    itmp = atpold; atpold = atpnew; atpnew = itmp;
    clear_cache();
    return 0;
}

int ofmo_get_interfrag_distance( int ifrag, double dist[] ) {
    ofmo_worker_get( OFMO_DISTA, ifrag, dist );
    return 0;
}

int ofmo_put_interfrag_distance( int ifrag, double dist[] ) {
    ofmo_worker_put( OFMO_DISTA, ifrag, dist );
    return 0;
}

int ofmo_put_monomer_energy( int ifrag, double src[] ) {
    ofmo_worker_put( OFMO_ENERGY, ifrag, src );
    return 0;
}

int ofmo_get_monomer_energy( int ifrag, double dist[] ) {
    ofmo_worker_get( OFMO_ENERGY, ifrag, dist );
    return 0;
}

int ofmo_put_monomer_energy0( int ifrag, double src[] ) {
    ofmo_worker_put( OFMO_ENERGY0, ifrag, src );
    return 0;
}

int ofmo_get_monomer_energy0( int ifrag, double dist[] ) {
    ofmo_worker_get( OFMO_ENERGY0, ifrag, dist );
    return 0;
}

// dimer SCF時のpopulationに対するaccumulation関連の関数
int ofmo_acc_dimer_aopop( int nao, double *aopop, int *iconv ) {
    int data_id = OFMO_TOTAL_AOPOP;
    ofmo_worker_acc( data_id, nao, aopop, iconv );
    return 0;
}

int ofmo_acc_dimer_atpop( int nat, double *atpop, int *iconv ) {
    int data_id = OFMO_TOTAL_ATPOP;
    ofmo_worker_acc( data_id, nat, atpop, iconv );
    return 0;
}
