/*
 * ワーカーからデータストアへのアクセスを 行うための関数群
 *
 * ofmo_worker_put/get() をfalanx用に変更した.
 *
 */

#include <assert.h>

#include <falanx.h>
#include <datastore.h>

#include "ofmo-datastore.h"
#include "ofmo-storage.h"


/*
 * データストアにデータを書き込む関数
 */
int
ofmo_worker_put(int data_id, int ifrag, double *src )
{
    int rc;
    falanx_context_t* context = NULL;
    data_store_t* ds = NULL;

    context = ofmo_storage_get_context();
    assert(context != NULL);
    ds = falanx_context_get_data_store(context);
    assert(ds != NULL);

    rc = ofmo_master_put(ds, data_id, ifrag, src);
    return rc;
}


/*
 * データストアからデータを読み込む関数
 */
int
ofmo_worker_get( int data_id, int ifrag, double *dest )
{
    int rc;
    falanx_context_t* context = NULL;
    data_store_t* ds = NULL;

    context = ofmo_storage_get_context();
    assert(context != NULL);
    ds = falanx_context_get_data_store(context);
    assert(ds != NULL);

    rc = ofmo_master_get(ds, data_id, ifrag, dest);
    return rc;
}

