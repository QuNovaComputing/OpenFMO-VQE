
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "ofmo-def.h"
#include "ofmo-data.h"
#include "ofmo-monomer-data.h"
#include "ofmo-task-util.h"

/*/
 * 初期密度行列計算の入力キーフォーマット
 *   "<#fragment>:<dold><dnew><aopold><aopnew><atpold><atpnew>"
 */
#define INIT_DENS_INPUT_FORMAT "%d:%1d%1d%1d%1d%1d%1d"

/*
 * SCFタスクの入力キーフォーマット
 *   "scfin:<method>:<nmon>:<scc>:<scfconf>:<mon1>:<mon2>:<mon3>:<data_ids>"
 */
#define SCF_INPUT_FORMAT  "scfin:%d:%d:%d:%d"   \
                          ":%d:%d:%d"           \
                          ":%1d%1d%1d%1d%1d%1d"

/* SCFタスクに必要なパラメータの数.  */
#define NUM_OF_SCF_PARAMS (13)

/*
 * SCFタスクの出力キーフォーマット
 *   "scfout:<method>:<nmon>:<scc>:<mon1>:<mon2>:<mon3>"
 */
#define SCF_OUTPUT_FORMAT "scfout:%d:%d:%d:%d:%d:%d"

/*
 * Approxタスクの入力キーフォーマット
 */
#define APPROX_INPUT_FORMAT "aprxin:%1d:%d:%d"

/*
 * Approxタスクの出力キーフォーマット
 */
#define APPROX_OUTPUT_FORMAT "aprxout:%1d:%d:%d"

/*
 * 初期密度行列計算用入力キーを生成する
 */
int
ofmo_init_dens_inkey_make(char* key, size_t size, int ifrag, int data_ids[6])
{
    int n;
    n = snprintf(key, size, INIT_DENS_INPUT_FORMAT,
                 ifrag,
                 data_ids[0], data_ids[1], data_ids[2], data_ids[3], data_ids[4], data_ids[5]);
    if (n < 0 || n >= size) {
        return -1;
    }
    return 0;
}

void
ofmo_scf_output_en_zero(ofmo_scf_output_en_t* out)
{
    out->energy = out->energy0 = out->ddv = 0.e0;
}

/*
 * 文字列形式の入力から ofmo_scf_input_t構造体のフィールドを埋める
 */
int
ofmo_scf_input_read(ofmo_scf_input_t* in, char* key, size_t len)
{
    int n;

    UNUSED(len);

    n = sscanf(key, SCF_INPUT_FORMAT,
               &in->method, &in->nmonomer, &in->iscc, &in->itol,
               &in->monomer_list[0], &in->monomer_list[1], &in->monomer_list[2],
               &in->data_ids[0], &in->data_ids[1], &in->data_ids[2],
               &in->data_ids[3], &in->data_ids[4], &in->data_ids[5]);
    if (n != NUM_OF_SCF_PARAMS) {
        return -1;
    }
    return 0;
}

void
ofmo_scf_input_show(ofmo_scf_input_t* scfin)
{
    printf("SCF input\nmethod: %d\nnmonomer: %d\niscc: %d\nitol: %d\n"
           "mon1: %d\nmon2: %d\nmon3: %d\n"
           "data_id[0]: %d\ndata_id[1]: %d\n"
           "data_id[2]: %d\ndata_id[3]: %d\n"
           "data_id[4]: %d\ndata_id[5]: %d\n",
           scfin->method, scfin->nmonomer, scfin->iscc, scfin->itol,
           scfin->monomer_list[0], scfin->monomer_list[1], scfin->monomer_list[2],
           scfin->data_ids[0], scfin->data_ids[1], scfin->data_ids[2],
           scfin->data_ids[3], scfin->data_ids[4], scfin->data_ids[5]);
}

/*
 * OFMO_SCF の入力キーを生成する関数.
 */
int
ofmo_scf_inkey_make(
    char* inputkey, size_t size,
    int method,
    int nmon,
    int scc,
    int conv,
    int mon1, int mon2, int mon3,
    int data_ids[6])
{
    int n;

    assert(is_valid_data_id(data_ids));

    n = snprintf(inputkey, size,
                 SCF_INPUT_FORMAT,
                 method, nmon, scc, conv,
                 mon1, mon2, mon3,
                 data_ids[0], data_ids[1], data_ids[2],
                 data_ids[3], data_ids[4], data_ids[5]);
    if (n != NUM_OF_SCF_PARAMS) {
        return -1;
    }
    return 0;
}

/*
 * OFMO_SCF の出力キーを生成する関数
 */
int
ofmo_scf_outkey_make(char* outkey, size_t size, ofmo_scf_input_t* scfin)
{
    int n;
    n = snprintf(outkey, size,
                SCF_OUTPUT_FORMAT,
                scfin->method, scfin->nmonomer, scfin->iscc,
                scfin->monomer_list[0], scfin->monomer_list[1], scfin->monomer_list[2]);
    if (n < 0 || n >= size) {
        return -1;
    }
    return 0;
}

/*
 * OFMO_SCF の出力キーを解析する
 */
int
ofmo_scf_outkey_parse(
    char* outkey, size_t size,
    int* outMethod, int* outNmon, int* outScc, int* outMonomerList)
{
    int n;
    int method, nmon, scc, mon1, mon2, mon3;

    UNUSED(size);

    n = sscanf(outkey, SCF_OUTPUT_FORMAT, &method, &nmon, &scc, &mon1, &mon2, &mon3);
    if (n != 6) {
        return -1;
    }

    if (outMethod) { *outMethod = method; }
    if (outNmon) { *outNmon = nmon; }
    if (outScc) { *outScc = scc;  }
    if (outMonomerList) {
        outMonomerList[0] = mon1;
        outMonomerList[1] = mon2;
        outMonomerList[2] = mon3;
    }
    return 0;
}

/* 構造体の生成用関数 */
ofmo_scf_output_ao_t*
ofmo_scf_output_ao_new(int nbody, int maxnfao)
{
    ofmo_scf_output_ao_t* p = NULL;
    p = (ofmo_scf_output_ao_t*)malloc(sizeof(ofmo_scf_output_ao_t));
    if (p) {
        p->fnao = 0;
        p->daopop = (double*)malloc(sizeof(double)*nbody*maxnfao);
        p->fsao2tuao = (int*)malloc(sizeof(int)*nbody*maxnfao);
        if (p->daopop == NULL) { goto error; }
        if (p->fsao2tuao == NULL) { goto error; }
    }
    return p;

error:
    ofmo_scf_output_ao_free(p);
    return NULL;
}

ofmo_scf_output_at_t*
ofmo_scf_output_at_new(int nbody, int maxnfatom)
{
    ofmo_scf_output_at_t* p = NULL;
    p = (ofmo_scf_output_at_t*)malloc(sizeof(ofmo_scf_output_at_t));
    if (p) {
        p->fnat = 0;
        p->datpop = (double*)malloc(sizeof(double)*nbody*maxnfatom);
        p->fatom2tatom = (int*)malloc(sizeof(int)*nbody*maxnfatom);
        if (p->datpop == NULL) { goto error; }
        if (p->fatom2tatom == NULL) { goto error; }
    }
    return p;

error:
    ofmo_scf_output_at_free(p);
    return NULL;
}

/* 構造体の解放用関数 */
void
ofmo_scf_output_ao_free(ofmo_scf_output_ao_t* scfout)
{
    if (scfout) {
        free(scfout->daopop);
        free(scfout->fsao2tuao);
        free(scfout);
    }
}

void
ofmo_scf_output_at_free(ofmo_scf_output_at_t* scfout)
{
    if (scfout) {
        free(scfout->datpop);
        free(scfout->fatom2tatom);
        free(scfout);
    }
}

#if defined(MPICH) || defined(MPICH2)
#define Pack(IN, INCNT, TYPE, OUT, OUTCNT, POS)     \
    MPI_Pack(IN, INCNT, TYPE, OUT, OUTCNT, POS, MPI_COMM_WORLD)
#define Unpack(IN, INCNT, POS, OUT, OUTCNT, TYPE)   \
    MPI_Unpack(IN, INCNT, POS, OUT, OUTCNT, TYPE, MPI_COMM_WORLD)
#define position_t int
#else
#define Pack(IN, INCNT, TYPE, OUT, OUTSZ, POS)      \
    MPI_Pack_external("external32", IN, INCNT, TYPE, OUT, OUTSZ, POS)
#define Unpack(IN, INCNT, POS, OUT, OUTCNT, TYPE)   \
    MPI_Unpack_external("external32", IN, INCNT, POS, OUT, OUTCNT, TYPE)
#define position_t MPI_Aint
#endif

/* SCF output 構造体のpack処理  */
int
ofmo_scf_output_ao_pack(ofmo_scf_output_ao_t* scfout, char** outBuf, size_t* outSize)
{
    char* p = NULL;
    size_t sz;
    int rc;
    position_t position = 0;

    sz = sizeof(int);
    sz += sizeof(double) * scfout->fnao;
    sz += sizeof(int) * scfout->fnao;

    p = (char*)malloc(sz);
    if (p) {
        rc = Pack(&scfout->fnao, 1, MPI_INT, p, sz, &position);
        if (rc != MPI_SUCCESS) { goto error; }
        rc = Pack(scfout->daopop, scfout->fnao, MPI_DOUBLE, p, sz, &position);
        if (rc != MPI_SUCCESS) { goto error; }
        rc = Pack(scfout->fsao2tuao, scfout->fnao, MPI_INT, p, sz, &position);
        if (rc != MPI_SUCCESS) { goto error; }
        assert(sz == position);

        *outBuf  = p;
        *outSize = sz;

        return 0;
    }

error:
    free(p);
    return -1;
}

int
ofmo_scf_output_at_pack(ofmo_scf_output_at_t* scfout, char** outBuf, size_t* outSize)
{
    char* p = NULL;
    size_t sz;
    int rc;
    position_t position = 0;

    sz = sizeof(int);
    sz += sizeof(double) * scfout->fnat;
    sz += sizeof(int) * scfout->fnat;

    p = (char*)malloc(sz);
    if (p) {
        rc = Pack(&scfout->fnat, 1, MPI_INT, p, sz, &position);
        if (rc != MPI_SUCCESS) { goto error; }
        rc = Pack(scfout->datpop, scfout->fnat, MPI_DOUBLE, p, sz, &position);
        if (rc != MPI_SUCCESS) { goto error; }
        rc = Pack(scfout->fatom2tatom, scfout->fnat, MPI_INT, p, sz, &position);
        if (rc != MPI_SUCCESS) { goto error; }

        assert(sz == position);

        *outBuf  = p;
        *outSize = sz;

        return 0;
    }

error:
    free(p);
    return -1;
}

int
ofmo_scf_output_ao_unpack(char* data, size_t size, ofmo_scf_output_ao_t* scfout)
{
    int rc;
    position_t position = 0;

    rc = Unpack(data, size, &position, &scfout->fnao, 1, MPI_INT);
    if (rc != MPI_SUCCESS) { return -1; }
    rc = Unpack(data, size, &position, scfout->daopop, scfout->fnao, MPI_DOUBLE);
    if (rc != MPI_SUCCESS) { return -1; }
    rc = Unpack(data, size, &position, scfout->fsao2tuao, scfout->fnao, MPI_INT);
    if (rc != MPI_SUCCESS) { return -1; }

    return 0;
}

int
ofmo_scf_output_at_unpack(char* data, size_t size, ofmo_scf_output_at_t* scfout)
{
    int rc;
    position_t position = 0;

    rc = Unpack(data, size, &position, &scfout->fnat, 1, MPI_INT);
    if (rc != MPI_SUCCESS) { return -1; }
    rc = Unpack(data, size, &position, scfout->datpop , scfout->fnat, MPI_DOUBLE);
    if (rc != MPI_SUCCESS) { return -1; }
    rc = Unpack(data, size, &position, scfout->fatom2tatom, scfout->fnat, MPI_INT);
    if (rc != MPI_SUCCESS) { return -1; }
    
    return 0;
}

int
ofmo_scf_get_ao_pop_key(const char* key, char* outKey, size_t size)
{
    int n;
    n = snprintf(outKey, size, "%s:ao", key);
    if (n < 0 || n >= size) {
        return -1;
    }
    return 0;
}

int
ofmo_scf_get_at_pop_key(const char* key, char* outKey, size_t size)
{
    int n;
    n = snprintf(outKey, size, "%s:at", key);
    if (n < 0 || n >= size) {
        return -1;
    }
    return 0;
}

int
ofmo_approx_inkey_make(char* outkey, size_t size, int data_id, int njob, int seq)
{
    int n;
    
    n = snprintf(outkey, size, APPROX_INPUT_FORMAT, data_id, njob, seq);
    if (n < 0 || n >= size) {
        return -1;
    }
    return 0;
}

int
ofmo_approx_outkey_make(char* outkey, size_t size, ofmo_approx_input_t* aprxin)
{
    int n;

    n = snprintf(outkey, size, APPROX_OUTPUT_FORMAT,
                 aprxin->data_id, aprxin->njob, aprxin->seq);
    if (n < 0 || n >= size) {
        return -1;
    }
    return 0;
}

int
ofmo_approx_inkey_parse(
    char* inputkey,
    size_t size,
    int* outDataId,
    int* outNjob,
    int* outSeq)
{
    int n;
    UNUSED(size);
    n = sscanf(inputkey, APPROX_INPUT_FORMAT, outDataId, outNjob, outSeq);
    if (n != 3) {
        return -1;
    }
    return 0;
}

int
ofmo_is_dimer_scf_input_key(char* key, size_t size)
{
/*
 * Dimerタスクの入力キーは文字列の最初が "scfin:1:2:" となるためそれを確認する.
 */
    return ((size >= 9) && (strncmp(key, "scfin:1:2", 9) == 0));
}

int
ofmo_is_dimer_scf_output_key(char* key, size_t size)
{
/*
 * Dimerタスクの出力キーは文字列の最初が "scfout:1:2:" となるためそれを確認する.
 */
    return ((size >= 10) && (strncmp(key, "scfout:1:2", 10) == 0));
}

int
ofmo_is_approx_input_key(char* key, size_t size)
{
    return ((size >= 7) && (strncmp(key, "aprxin:", 7) == 0));
}
  
int
ofmo_is_approx_output_key(char* key, size_t size)
{
    return ((size >= 8) && (strncmp(key, "aprxout:", 8) == 0));
}



