/**
 * モノマー密度行列のキャッシュ機構を導入した
 *
 */

#include <assert.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#include "ofmo-def.h"
#include "ofmo-data.h"
#include "ofmo-datastore.h"
#include "ofmo-prof.h"
#include "ofmo-storage.h"

#define MAXDENCACHE 16

static time_t*  LastUse       = NULL;
static int*     CachedFrag    = NULL;
static double** CachedDens    = NULL;
static int      _maxdencache_ = 0;

/* profile information */
static int NCALLED = 0; /* # of calling of ofmo_get_monomer_density */
static int NHITS   = 0; /* # of cache hits in ofmo_get_monomer_density */

/* deallocaion of cache for monomer density matrices */
static void
dealloc_cache()
{
    Free(CachedDens[0]);
    Free(CachedDens);
    Free(CachedFrag);
    Free(LastUse);
}

void
ofmo_show_cache_prof()
{
    if (fp_prof) {
        double rate;
        rate = (double) NHITS / (double) NCALLED * 100.e0;
        fprintf(fp_prof, "CALLED=%6d, HITS=%6d, rate=%7.2f(%%)\n", NCALLED, NHITS, rate);
        fflush(fp_prof);
    }
}

/* memory allocation of cache for monomer density matrices */
static int
alloc_cache()
{
    static int called = false;
    if (called) { return 0; }
    int   maxdencache;
    char* str;
    if ((str = getenv("OFMO_MAXDENCACHE")) != NULL) {
        maxdencache = atoi(str);
        if (maxdencache < 0) { maxdencache = 0; }
    } else {
        maxdencache = MAXDENCACHE;
    }
    if (maxdencache < 1) {
        _maxdencache_ = 0;
        if (fp_prof) {
            fprintf(fp_prof, "== density matrices is not cached\n");
        }
        called = true;
        return 0;
    } else {
        if (fp_prof) {
            fprintf(fp_prof, "== density cache = %d ==\n", maxdencache);
        }
    }
    int maxnfao, n2;
    if (ofmo_data_get_vals("maxnfao", &maxnfao) != 0) {
        return -1;
    }
    n2 = maxnfao * (maxnfao + 1) / 2;
    LastUse       = (time_t*) malloc(sizeof(time_t) * maxdencache);
    CachedFrag    = (int*) malloc(sizeof(int) * maxdencache);
    CachedDens    = (double**) malloc(sizeof(double*) * maxdencache);
    CachedDens[0] = (double*) malloc(sizeof(double) * n2 * maxdencache);
    if (LastUse == NULL || CachedFrag == NULL || CachedDens == NULL ||
        CachedDens[0] == NULL) {
        dbg("error: failure in memory allocation\n");
        return -1;
    }
    for (int i = 1; i < maxdencache; i++) {
        CachedDens[i] = CachedDens[i - 1] + n2;
    }
    for (int i = 0; i < maxdencache; i++) {
        LastUse[i]    = 0;
        CachedFrag[i] = -1;
    }
    // information of memory allocation
    if (fp_prof) {
        double dsize;
        dsize = (double) (maxdencache * n2 * sizeof(double)) / (double) (1024 * 1024);
        fprintf(fp_prof,
                "== allocd memory size in ofmo-monomer-data.c = "
                "%10.3f MB\n", dsize);
    }
    atexit(dealloc_cache);
    _maxdencache_ = maxdencache;
    called = true;
    return 0;
}

/* clear cache */
static void
clear_cache()
{
    int i;
    if (_maxdencache_ < 1) { return; }
    for (i = 0; i < _maxdencache_; i++) {
        CachedFrag[i] = -1;
    }
    if (fp_prof) {
        fprintf(fp_prof, "clear cache\n");
    }
}

/* density matrix of ifrag is cached (>=0) or not (-1) */
static int
is_cached(const int ifrag)
{
    if (_maxdencache_ < 1) { return -1; }
    for (int i = 0; i < _maxdencache_; i++) {
        if (CachedFrag[i] == ifrag) { return i; }
    }
    return -1;
}

/* get cache location of oldest reffered one */
static int
is_oldest()
{
    int    oldest_loc;
    time_t oldest_time;
    if (_maxdencache_ < 1) { return -1; }
    oldest_loc  = 0;
    oldest_time = LastUse[0];
    for (int i = 1; i < _maxdencache_; i++) {
        if (LastUse[i] < oldest_time) {
            oldest_time = LastUse[i];
            oldest_loc  = i;
        }
    }
    return oldest_loc;
}

void
ofmo_reset_monomer_data_workerp()
{
}

/* データの更新（データIDの交換）を行う */
int
ofmo_update_monomer_data()
{
    clear_cache();
    return 0;
}

/** モノマー密度行列を取得する関数（集団呼び出し）
 *
 * モノマー密度行列を取得する関数。終了時には、引数で与えられた
 * コミュニケータに属するすべてのプロセスが、密度行列データを取得する。
 * この関数内部で、最近参照したいくつかのモノマー密度行列データが
 * キャッシュされているので、キャッシュヒット時には、通信せずに
 * キャッシュされたデータをコピーしてもちいる。キャッシュミスした
 * 場合には、指定されたモノマー密度行列データを通信で
 * 取得して、もっとも古く参照したデータの代わりにキャッシュする。
 *
 * @li \c MAXDENCACHE キャッシュする密度行列データの最大数
 *
 * @ingroup ofmo-calc
 *
 * @parm comm worker groupのコミュニケータ
 **/
int
ofmo_get_monomer_density(MPI_Comm comm, const int ifrag, double D[])
{
    static int  called = false;
    static int* nfao;
    if (!called) {
        int ierr;
        ierr = ofmo_data_get_vals("nfao", &nfao);
        if (ierr != 0) {
            return -1;
        }
        alloc_cache();
        called = true;
    }
    int n, n2, id;
    int root = 0, myrank;

    int data_ids[6];
    falanx_context_t* context = NULL;
    data_store_t* ds = NULL;

    ofmo_storage_get_data_ids(data_ids);
    assert((data_ids[0] == OFMO_DENS1)||(data_ids[0] == OFMO_DENS2));
    context = ofmo_storage_get_context();
    assert(context != NULL);
    ds = falanx_context_get_data_store(context);
    assert(ds != NULL);

    MPI_Comm_rank(comm, &myrank);

    n  = nfao[ifrag];
    n2 = ((n * n + n) >> 1);
    id = is_cached(ifrag);
    if (id < 0) {
        if (myrank == root) {
            int _rc;
            id = is_oldest();
            /* dold */
            _rc = ofmo_master_get(ds, data_ids[0], ifrag, D);
            assert(_rc == 0);
        }
        MPI_Bcast(&id, 1, MPI_INT, root, comm);
        MPI_Bcast(D, n2, MPI_DOUBLE, root, comm);
        if (id >= 0) {
            memcpy(CachedDens[id], D, sizeof(double) * n2);
            LastUse[id]    = time(NULL);
            CachedFrag[id] = ifrag;
        }
    } else {
        memcpy(D, CachedDens[id], sizeof(double) * n2);
        LastUse[id] = time(NULL);
        NHITS++;
    }
    /* for profiling */
    NCALLED++;
    return 0;
}

/* 呼び出したプロセスのみ、モノマー密度行列データを更新 */
int
ofmo_put_monomer_density(int ifrag, double D[])
{
    int rc;
    int data_ids[6];
    falanx_context_t* context = NULL;
    data_store_t* ds = NULL;

    ofmo_storage_get_data_ids(data_ids);
    assert((data_ids[1] == OFMO_DENS1)||(data_ids[1] == OFMO_DENS2));
    context = ofmo_storage_get_context();
    assert(context != NULL);
    ds = falanx_context_get_data_store(context);
    assert(ds != NULL);

    /* dnew */
    rc = ofmo_master_put(ds, data_ids[1], ifrag, D);
    return rc;
}

static int
ofmo_get_monomer_aopop_all(data_store_t* ds, int aopold, double aopop[])
{
    int i;
    int rc;
    int offset = 0;
    static int called = false;
    static int nfrag, *nfao;

    if (!called) {
        ofmo_data_get_vals("nfrag nfao", &nfrag, &nfao);
        called = true;
    }

    for (i = 0; i < nfrag; i++) {
        rc = ofmo_master_get_val(ds, aopold, i,
                                (char*)&aopop[offset],
                                (sizeof(double) * nfao[i]));
        if (rc != 0) {
            break;
        }
        offset += nfao[i];
    }
    return rc;
}

/* 呼び出したプロセスのみ、モノマーAOPOPを取得 */
int
ofmo_get_monomer_aopop(int ifrag, double aopop[])
{
    int rc;
    int data_ids[6];
    falanx_context_t* context = NULL;
    data_store_t* ds = NULL;

    ofmo_storage_get_data_ids(data_ids);
    assert((data_ids[2] == OFMO_AOPOP1)||(data_ids[2] == OFMO_AOPOP2));
    context = ofmo_storage_get_context();
    assert(context != NULL);
    ds = falanx_context_get_data_store(context);
    assert(ds != NULL);

    if (ifrag < 0) {
        /* 全データを取り出す */
        rc = ofmo_get_monomer_aopop_all(ds, data_ids[2], aopop);
        return rc;
    }

    /* aopold */
    rc = ofmo_master_get(ds, data_ids[2], ifrag, aopop);
    return rc;
}

/* 呼び出したプロセスのみ、モノマーAOPOPを更新 */
int
ofmo_put_monomer_aopop(int ifrag, double aopop[])
{
    int rc;
    int data_ids[6];
    falanx_context_t* context = NULL;
    data_store_t* ds = NULL;

    ofmo_storage_get_data_ids(data_ids);
    assert((data_ids[3] == OFMO_AOPOP1)||(data_ids[3] == OFMO_AOPOP2));
    context = ofmo_storage_get_context();
    assert(context != NULL);
    ds = falanx_context_get_data_store(context);
    assert(ds != NULL);

    /* aopnew */
    rc = ofmo_master_put(ds, data_ids[3], ifrag, aopop);
    return rc;
}


static int
ofmo_get_monomer_atpop_all(data_store_t* ds, int atpold, double atpop[])
{
    int i;
    int rc;
    int offset = 0;
    static int called = false;
    static int nfrag, * nfatom;

    if (!called) {
        ofmo_data_get_vals("nfrag nfatom", &nfrag, &nfatom);
        called = true;
    }

    for (i = 0; i < nfrag; i++) {
        rc = ofmo_master_get_val(ds, atpold, i,
                                (char*)&atpop[offset],
                                (sizeof(double) * nfatom[i]));
        if (rc != 0) {
            break;
        }
        offset += nfatom[i];
    }
    return rc;
}

/* 呼び出したプロセスのみ、モノマーのatomic populationを取得 */
int
ofmo_get_monomer_atpop(int ifrag, double atpop[])
{
    int rc;
    int data_ids[6];
    falanx_context_t* context = NULL;
    data_store_t* ds = NULL;

    ofmo_storage_get_data_ids(data_ids);
    assert((data_ids[4] == OFMO_ATPOP1)||(data_ids[4] == OFMO_ATPOP2));
    context = ofmo_storage_get_context();
    assert(context != NULL);
    ds = falanx_context_get_data_store(context);
    assert(ds != NULL);

    if (ifrag < 0) {
        rc = ofmo_get_monomer_atpop_all(ds, data_ids[4], atpop);
        return rc;
    }

    /* atpold */
    rc = ofmo_master_get(ds, data_ids[4], ifrag, atpop);
    return rc;
}

/* 呼び出したプロセスのみ、モノマーのatomic populationを更新 */
int
ofmo_put_monomer_atpop(int ifrag, double atpop[])
{
    int rc;
    int data_ids[6];
    falanx_context_t* context = NULL;
    data_store_t* ds = NULL;

    ofmo_storage_get_data_ids(data_ids);
    assert((data_ids[5] == OFMO_ATPOP1)||(data_ids[5] == OFMO_ATPOP2));
    context = ofmo_storage_get_context();
    assert(context != NULL);
    ds = falanx_context_get_data_store(context);
    assert(ds != NULL);

    /* atpnew */
    rc = ofmo_master_put(ds, data_ids[5], ifrag, atpop);
    return rc;
}

int
ofmo_put_monomer_energy(int ifrag, double src[])
{
    int rc;
    falanx_context_t* context = NULL;
    data_store_t* ds = NULL;

    context = ofmo_storage_get_context();
    assert(context != NULL);
    ds = falanx_context_get_data_store(context);
    assert(ds != NULL);

    rc = ofmo_master_put(ds, OFMO_ENERGY, ifrag, src);
    return rc;
}

int
ofmo_get_monomer_energy(int ifrag, double dist[])
{
    int rc;
    falanx_context_t* context = NULL;
    data_store_t* ds = NULL;

    context = ofmo_storage_get_context();
    assert(context != NULL);
    ds = falanx_context_get_data_store(context);
    assert(ds != NULL);


    rc = ofmo_master_get(ds, OFMO_ENERGY, ifrag, dist);
    return rc;
}

int
ofmo_put_monomer_energy0(int ifrag, double src[])
{
    int rc;
    falanx_context_t* context = NULL;
    data_store_t* ds = NULL;

    context = ofmo_storage_get_context();
    assert(context != NULL);
    ds = falanx_context_get_data_store(context);
    assert(ds != NULL);

    rc = ofmo_master_put(ds, OFMO_ENERGY0, ifrag, src);
    return rc;
}

int
ofmo_get_monomer_energy0(int ifrag, double dist[])
{
    int rc;
    falanx_context_t* context = NULL;
    data_store_t* ds = NULL;

    context = ofmo_storage_get_context();
    assert(context != NULL);
    ds = falanx_context_get_data_store(context);
    assert(ds != NULL);

    rc = ofmo_master_get(ds, OFMO_ENERGY0, ifrag, dist);
    return rc;
}






