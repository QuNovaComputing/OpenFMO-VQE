
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <assert.h>

#include <mpi.h>

#include "ofmo-def.h"
#include "ofmo-data-struct.h"
#include "ofmo-data.h"

#include "ofmo-misc.h"
#include "ofmo-prof.h"

#include <falanx.h>
#include <datastore.h>

#include "ofmo-task-util.h"
#include "ofmo-datastore.h"

#define FILE_NAME_FORMAT "%s-%s.dat"

int ofmo_master_get_energy_all(data_store_t* ds, const int data_id, double* outBuff, int nfrag);

typedef struct ofmo_monomer_file_st
{
    int     type;
    int     nfrag;
    int*    nelems;
    double* data;
    FILE*   fp;
    char    filename[MAXSTRLEN];
} ofmo_monomer_file_t;

const static int DATA_IDS[6] = {
    OFMO_DENS1,
    OFMO_AOPOP1, OFMO_ATPOP1,
    OFMO_ENERGY, OFMO_ENERGY0,
    OFMO_DISTA
};

static const char* SUFFIX[] = {
    "dens", "aop", "atp", "ene", "ene0", "dist"
};

static int dold   = OFMO_DENS1;
static int dnew   = OFMO_DENS2;
static int aopold = OFMO_AOPOP1;
static int aopnew = OFMO_AOPOP2;
static int atpold = OFMO_ATPOP1;
static int atpnew = OFMO_ATPOP2;

static void
ofmo_reset_monomer_data_master(void)
{
    dold   = OFMO_DENS1;
    dnew   = OFMO_DENS2;
    aopold = OFMO_AOPOP1;
    aopnew = OFMO_AOPOP2;
    atpold = OFMO_ATPOP1;
    atpnew = OFMO_ATPOP2;
}

#define swap_i(OLD, NEW) do { int t = (OLD); (OLD)=(NEW); (NEW)=t; } while (0)
int
ofmo_update_monomer_data_master(void)
{
    swap_i(dold,   dnew);
    swap_i(aopold, aopnew);
    swap_i(atpold, atpnew);
    return 0;
}

void
ofmo_get_monomer_data_master(int data_ids[6])
{
    data_ids[0] = dold;
    data_ids[1] = dnew;
    data_ids[2] = aopold;
    data_ids[3] = aopnew;
    data_ids[4] = atpold;
    data_ids[5] = atpnew;
}


/*
 * データストア用のキーを生成する.
 *  - キーの書式は "%d:%d", data_id, ifrag とした.
 */
static int
ofmo_ds_make_key(char key[MAX_KEY_LENGTH], int data_id, int ifrag)
{
    int n = snprintf(key, MAX_KEY_LENGTH, "%d:%d", data_id, ifrag);
    if (n >= MAX_KEY_LENGTH) {
        return -1;
    }
    return 0;
}

/*
 * データストアからデータを取得する関数
 *
 * @param ds データストアクライアントハンドル
 * @param data_id マスタからの依頼のID (see ofmo-def.h)
 * @param ifrag フラグメント番号
 * @param outBuff 出力バッファ
 */
int
ofmo_master_get(
    data_store_t* ds,
    int      data_id,
    int      ifrag,
    double   outBuff[])
{
    int rc;
    char key[MAX_KEY_LENGTH];
    char* value = NULL;
    size_t valueSize = 0;

    assert(ifrag >= 0);

    ofmo_ds_make_key(key, data_id, ifrag); 
    rc = data_store_get(ds, key, strlen(key), &value, &valueSize);
    if (rc == FLX_SUCCESS) {
        memcpy(outBuff, value, valueSize);
        free(value);
    } else {
        rc = -1;
    }
    return rc;
}

/* データストアからデータを取得する関数
 */
int
ofmo_master_get_val(
    data_store_t* ds,
    int      data_id,
    int      ifrag,
    char*    outBuff,
    size_t   buffSize)
{
    int rc;
    char key[MAX_KEY_LENGTH];

    ofmo_ds_make_key(key, data_id, ifrag);
    rc = data_store_get_value(ds, key, strlen(key), outBuff, buffSize);
    if (rc != FLX_SUCCESS) {
        rc = -1;
    }
    return rc;
}


static int
ofmo_master_get_num_of_elements(
    const int  data_id,
    const int  ifrag,
    const int* nfao,
    const int  nfao_total,
    const int* nfatom,
    const int  nfatom_total,
    const int  nfrag,
    const int  nao,
    const int  natom)
{
    int nelem = 0;
    switch (data_id) {
    case OFMO_DENS1:
    case OFMO_DENS2:
        if (ifrag >= 0) {
            int n  = nfao[ifrag];
            nelem  = (n * n + n) >> 1;
        }
        break;
    case OFMO_AOPOP1:
    case OFMO_AOPOP2:
        nelem  = (ifrag < 0 ? nfao_total : nfao[ifrag]);
        break;
    case OFMO_ATPOP1:
    case OFMO_ATPOP2:
        nelem  = (ifrag < 0 ? nfatom_total : nfatom[ifrag]);
        break;
    case OFMO_DISTA:
        nelem  = nfrag;
        break;
    case OFMO_ENERGY:
    case OFMO_ENERGY0:
        nelem  = (ifrag < 0 ? nfrag : 1);
        break;
    case OFMO_TOTAL_AOPOP:
        nelem  = nao;
        break;
    case OFMO_TOTAL_ATPOP:
        nelem  = natom;
        break;
    }
    return nelem;
}

/*
 * データストアにデータを書き込む.
 */
int
ofmo_master_put(
    data_store_t* ds,
    int           data_id,
    int           ifrag,
    double*       src)
{
    int rc;
    int nelem;
    char key[MAX_KEY_LENGTH];

    static int nao, natom;
    static int nfrag, * nfao, * nfatom, nfao_total, nfatom_total;
    static int called = false;
    if (!called) {
        ofmo_data_get_vals("nfrag nfao nfatom nao natom",
                           &nfrag, &nfao, &nfatom, &nao, &natom);
        nfao_total   = ofmo_isum2(nfrag, nfao);
        nfatom_total = ofmo_isum2(nfrag, nfatom);
        called       = true;
    }

    assert(ifrag >= 0);

    if (ifrag >= nfrag) {
        if (fp_prof) {
            fdbg(fp_prof, "ERROR: Illegal monomer number (%d)\n", ifrag);
            fflush(fp_prof);
        }
        return -1;
    }

    nelem = ofmo_master_get_num_of_elements(data_id, ifrag,
                                            nfao, nfao_total,
                                            nfatom, nfatom_total, nfrag, nao, natom);
    if (nelem <= 0) {
        if (fp_prof) {
            fdbg(fp_prof, "ERROR: Illegal parameter (ifrag=%d, data ID=%d)\n",
                 ifrag, data_id);
            fflush(fp_prof);
        }
        return -1;
    }

    ofmo_ds_make_key(key, data_id, ifrag);
    rc = data_store_set(ds, key, strlen(key), (char*)src, sizeof(double)*nelem);
    assert(rc == FLX_SUCCESS);

    return 0;
}

static int
ofmo_master_put_all(
    data_store_t* ds,
    int           data_id,
    double*       src,
    int           nfrag)
{
    int i;
    int rc;
    for (i = 0; i < nfrag; i++) {
        rc = ofmo_master_put(ds, data_id, i, src);
        if (rc < 0) {
            break;
        }
    }
    return rc;
}

int
ofmo_get_monomer_aopop_master(data_store_t* ds, int ifrag, double aopop[])
{
    return ofmo_master_get(ds, aopold, ifrag, aopop);
}

int
ofmo_get_monomer_atpop_master(data_store_t* ds, int ifrag, double atpop[])
{
    return ofmo_master_get(ds, atpold, ifrag, atpop);
}

static int
file_exists(const char* path)
{
    struct stat st;

    if (stat(path, &st) < 0) {
        if (errno != ENOENT) {
            perror("stat");
        }
        return 0;
    }
    return 1;
}

static void
ofmo_monomer_file_clear(ofmo_monomer_file_t* file)
{
    file->type   = 0;
    file->nfrag  = 0;
    file->nelems = NULL;
    file->data   = NULL;
    file->fp     = NULL;
    file->filename[0] = '\0';
}

static int
ofmo_monomer_file_get_max_data_size(ofmo_monomer_file_t* file)
{
    int i, size = 0;
    switch (file->type) {
    case 0:
        size = file->nfrag;
        break;
    case 1:
    case 2:
        for (i = 0; i < file->nfrag; i++) {
            if (size < file->nelems[i]) {
                size = file->nelems[i];
            }
        }
        break;
    default:
        MPI_Abort(MPI_COMM_WORLD, 1);
        break;
    }
    return size;
}

/*
 * ファイルからデータを読みこんでデータストアに書き込む.
 */
static int
ofmo_monomer_file_read(
    ofmo_monomer_file_t* file,
    data_store_t*        ds,
    const int            nfrag,
    const int            data_id)
{
    int rc = 0; /* success */
    int i, n, nr, size, ifrag;

    /* Read first fields */
    fread(&file->type, sizeof(int), 2, file->fp);
    if (file->nfrag != nfrag) {
        dbg("ERROR:%s: Illgel value of nfrag (given=%d, file=%d)\n",
            file->filename, nfrag, file->nfrag);
        return -1;
    }

    /* Read num of elems.
     * It needs for the type 1 or 2 file.
     */
    switch (file->type) {
    case 1:
        file->nelems = (int*)malloc(sizeof(int) * file->nfrag);
        assert(file->nelems != NULL);
        fread(&size, sizeof(int), 1, file->fp);
        for (i = 0; i < file->nfrag; i++) {
            file->nelems[i] = size;
        }
        break;
    case 2:
        file->nelems = (int*)malloc(sizeof(int) * file->nfrag);
        assert(file->nelems != NULL);
        nr = fread(file->nelems, sizeof(int), file->nfrag, file->fp);
        if (nr != nfrag) {
            dbg("ERROR:%s: Unexected EOF in reading nelems" " nfrag=%d, nread=%d\n",
                file->filename, file->nfrag, nr);
            free(file->nelems);
            return -1;
        }
        break;
    }

    /* Allocate data buffer */
    n = ofmo_monomer_file_get_max_data_size(file);
    file->data = (double*)malloc(sizeof(double) * n);
    assert(file->data != NULL);

    /* Read the data */
    switch (file->type) {
    case 0:
        nr = fread(file->data, sizeof(double), file->nfrag, file->fp);
        if (nr != nfrag) {
            dbg("ERROR:%s: Bad number of elemnts (%d vs %d )\n", file->filename, nr, nfrag);
            return -1;
        }
        ofmo_master_put_all(ds, data_id, file->data, file->nfrag);
        break;
    case 1:
    case 2:
        for (ifrag = 0; ifrag < file->nfrag; ifrag++) {
            nr = fread(file->data, sizeof(double), file->nelems[ifrag], file->fp);
            if (nr != file->nelems[ifrag]) {
                dbg("ERROR:%s: Bad number of elemnts (%d vs %d )\n",
                    file->filename, nr, file->nelems[ifrag]);
                rc = -1;
                break;
            }
            ofmo_master_put(ds, data_id, ifrag, file->data);
        }
        break;
    default:
        rc = -1;
        break;
    }

    free(file->nelems);
    free(file->data);

    return rc;
}

int
ofmo_read_and_put_all(
    ofmo_master_config_t* config,
    falanx_context_t*     context,
    char*                 fileheader)
{
    int i, rc;
    ofmo_monomer_file_t mono_file;
    data_store_t* ds = NULL;

    /* 必要なデータのファイルがあるかチェック */
    for (i = 0; i < 6; i++) {
        sprintf(mono_file.filename, FILE_NAME_FORMAT, fileheader, SUFFIX[i]);
        if (!file_exists(mono_file.filename)) {
            dbg("file[%s] is not found\n", mono_file.filename);
            return -2;
        }
    }

    ds = falanx_context_get_data_store(context);
    assert(ds != NULL);
    for (i = 0; i < 6; i++) {
        ofmo_monomer_file_clear(&mono_file);
        sprintf(mono_file.filename, FILE_NAME_FORMAT, fileheader, SUFFIX[i]);
        mono_file.fp = fopen(mono_file.filename, "r");
        assert(mono_file.fp != NULL);
        rc = ofmo_monomer_file_read(&mono_file, ds, config->nfrag, DATA_IDS[i]);
        fclose(mono_file.fp);
        if (rc < 0) {
            return -1;
        }
    }

    ofmo_reset_monomer_data_master();

    return 0;
}


/*
 * データをファイルに書き出す.
 * ファイルフォーマットは以下の通り.
 *
 * type 0 の場合
 *   1: type  (int) タイプ
 *   2: nfrag (int) フラグメント数
 *   3: data  (double * nfrag) データ
 * type 1 の場合
 *   1: type  (int) タイプ
 *   2: nfrag (int) フラグメント数
 *   3: size  (int) データ要素数(?)
 *   4: data  ( [(double * nelems[x]) for x in range(nfrag)] ) データ
 * type 2 の場合
 *   1: type  (int) タイプ
 *   2: nfrag (int) フラグメント数
 *   3: size  (int * nfrag) 各フラグメント毎の要素数(?)
 *   4: data  ( [(double * nelems[x]) for x in range(nfrag)] ) データ
 */
static int
ofmo_monomer_file_write(
    ofmo_monomer_file_t* file,
    data_store_t* ds,
    const int     data_id)
{
    int i, n = 0;
    double *buff = NULL;

    n = ofmo_monomer_file_get_max_data_size(file);
    buff = (double*)malloc(sizeof(double) * n);
    assert(buff != NULL);

    /* 先頭2要素を書き込み */
    fwrite(&file->type, sizeof(int), 2, file->fp);

    switch (file->type) {
    case 0:
        assert(data_id == OFMO_ENERGY || data_id == OFMO_ENERGY0); 
        ofmo_master_get_energy_all(ds, data_id, buff, file->nfrag);
        fwrite(buff, sizeof(double), file->nfrag, file->fp);
        break;
    case 1:
        fwrite(&file->nfrag, sizeof(int), 1, file->fp);
        for (i = 0; i < file->nfrag; i++) {
            ofmo_master_get(ds, data_id, i, buff);
            fwrite(buff, sizeof(double), file->nelems[i], file->fp);
        }
        break;
    case 2:
        fwrite(file->nelems, sizeof(int), file->nfrag, file->fp);
        for (i = 0; i < file->nfrag; i++) {
            ofmo_master_get(ds, data_id, i, buff);
            fwrite(buff, sizeof(double), file->nelems[i], file->fp);
        }
        break;
    default:
        return -1;
        break;
    }
    free(buff);
    return 0;
}

/** 計算で得られた（FMO計算続行に必要な）モノマーに間するデータを
 *  ファイルに出力する関数
 */
int
ofmo_get_and_write(falanx_context_t* context, char* fileheader, int nfrag)
{
    int  i;
    int rc;
    int  types[6] = { 2, 2, 2, 0, 0, 1 };
    int  data_ids[6];
    data_store_t* ds = NULL;
    ofmo_monomer_file_t mono_file;

    ds = falanx_context_get_data_store(context);
    assert(ds != NULL);

    data_ids[0] = dold;
    data_ids[1] = aopold;
    data_ids[2] = atpold;
    data_ids[3] = OFMO_ENERGY;
    data_ids[4] = OFMO_ENERGY0;
    data_ids[5] = OFMO_DISTA;

    for (i = 0; i < 6; i++) {
        ofmo_monomer_file_clear(&mono_file);

        mono_file.type   = types[i];
        mono_file.nfrag  = nfrag;
        mono_file.nelems = ofmo_getadd_data_nelems(data_ids[i]);

        sprintf(mono_file.filename, FILE_NAME_FORMAT, fileheader, SUFFIX[i]);
        if ((mono_file.fp = fopen(mono_file.filename, "w")) == NULL) {
            dbg("ERROR: Failure in open file %s\n", mono_file.filename);
            return -1;
        }
        rc = ofmo_monomer_file_write(&mono_file, ds, data_ids[i]);
        fclose(mono_file.fp);
        if (rc < 0) {
            return -1;
        }
    }
    return 0;
}

/*
 * ENERGY/ENERGY0データ全てを取り出す.
 * データストアに部分更新がないため、1つずつ取り出している.
 */
int
ofmo_master_get_energy_all(
    data_store_t* ds,
    const int     data_id,
    double*       outBuff,
    int           nfrag)
{
    int i;
    int rc;

    if (!(data_id == OFMO_ENERGY || data_id == OFMO_ENERGY0)) {
        return -1;
    }

    for (i = 0; i < nfrag; i++) {
        rc = ofmo_master_get_val(ds, data_id, i, (char*)&outBuff[i], sizeof(double));
        if (rc != 0) {
            break;
        }
    }
    return rc;
}

int
ofmo_master_get_scf_ao_pop(
    data_store_t* ds,
    char*         key,
    size_t        keysize,
    ofmo_scf_output_ao_t* out)
{
    int rc;
    char* data = NULL;
    size_t size;

    rc = data_store_get(ds, key, keysize, &data, &size);
    if (rc != FLX_SUCCESS) {
        return -1;
    }
    rc = ofmo_scf_output_ao_unpack(data, size, out);
    return rc;
}

int
ofmo_master_get_scf_at_pop(
    data_store_t* ds,
    char*         key,
    size_t        keysize,
    ofmo_scf_output_at_t* out)
{
    int rc;
    char* data = NULL;
    size_t size;

    rc = data_store_get(ds, key, keysize, &data, &size);
    if (rc != FLX_SUCCESS) {
        return -1;
    }
    rc = ofmo_scf_output_at_unpack(data, size, out);
    return rc;
}


