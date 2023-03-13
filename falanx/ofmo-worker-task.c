
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#include <mpi.h>
#include <falanx_worker.h>
#include <datastore.h>

#include "ofmo-def.h"
#include "ofmo-prof.h"
#include "ofmo-data.h"
#include "ofmo-monomer-data.h"
#include "ofmo-task-util.h"
#include "ofmo-calc-frag.h"
#include "ofmo-inter-frag.h"
#include "ofmo-storage.h"
#include "ofmo-tlog.h"

#define ofmo_clear_cache() ofmo_update_monomer_data()
static int iscc_prev = -1;

/*
 * 初期密度行列計算 タスク
 *
 * パラメータ
 *  "<I_MON1(int)>"
 *
 * 出力
 *  ""
 *
 * */
void
ofmo_init_dens_task(falanx_context_t* context, MPI_Comm worker)
{
    int rc;
    int rank, size;
    char key[MAX_KEY_LENGTH];
    size_t len;
    int n;
    /* ofmo_monomer_init_density() が monomer_list[6] を利用するため */
    int monomer_list[7];
    int data_ids[6];
    data_store_t* ds = NULL;

    MPI_Comm_rank(worker, &rank);
    MPI_Comm_size(worker, &size);
    ds = falanx_context_get_data_store(context);
    assert(ds != NULL);

    if (rank == 0) {
        rc = falanx_get_request_key_value(context, key, MAX_KEY_LENGTH, &len);
        assert(rc == FLX_SUCCESS);
        //dbg("%d:%d:arg:[%s]\n", rank, size, key);
    }
    MPI_Bcast(key, MAX_KEY_LENGTH, MPI_CHAR, 0, worker);

    n = sscanf(key, "%d:%1d%1d%1d%1d%1d%1d",
                &monomer_list[6],
                &data_ids[0],
                &data_ids[1],
                &data_ids[2],
                &data_ids[3],
                &data_ids[4],
                &data_ids[5]);
    assert(n == 7);

    ofmo_storage_set_context(context);
    ofmo_storage_set_data_ids(data_ids);

    ofmo_monomer_init_density(monomer_list, worker); /* ofmo-calc-frag.c */

    MPI_Barrier(worker);

    if (rank == 0) {
        len = snprintf(key, MAX_KEY_LENGTH, "%d", monomer_list[6]);
        falanx_set_response_key(context, key, len);
        //dbg("%d:%d:key:%s\n", rank, size, key);
    }

    ofmo_storage_set_context(NULL);
    ofmo_storage_set_data_ids(OFMO_STORAGE_DATA_IDS_NULL);
}

void
ofmo_distance_task(falanx_context_t* context, MPI_Comm worker)
{
    int rank;

    MPI_Comm_rank(worker, &rank);

    ofmo_storage_set_context(context);

    ofmo_make_inter_frag_distance_list(worker);

    if (rank == 0) {
        falanx_set_response_key(context, "", 1);
    }

    ofmo_storage_set_context(NULL);
}

/*
 * SCF計算タスク
 *
 * input:
 *
 * output:
 *
 * DSへの書き込み:
 *
 */
void
ofmo_scf_task(falanx_context_t* context, MPI_Comm worker)
{
    int rc;
    int rank, size;
    char key[MAX_KEY_LENGTH];
    size_t len;
    data_store_t* ds = NULL;

    /* iterationに変更があるか確認するためisccパラメータをstatic領域に保存する. */
//    static int iscc_prev = -1;
    static int nbody;
    static int maxnfao;
    static int maxnfatom;
    static bool called = false;

    double               tolscf;
    ofmo_scf_input_t     scfin;
    ofmo_scf_output_en_t scfouten;
    ofmo_scf_output_ao_t* scfoutao = NULL;
    ofmo_scf_output_at_t* scfoutat = NULL;

    if (!called) {
        /* common block */
        ofmo_data_get_vals("nbody maxnfao maxnfatom", &nbody, &maxnfao, &maxnfatom);
        called = true;
    }

    scfoutao = ofmo_scf_output_ao_new(nbody, maxnfao);
    scfoutat = ofmo_scf_output_at_new(nbody, maxnfatom);
    assert(scfoutao != NULL);
    assert(scfoutat != NULL);

    MPI_Comm_rank(worker, &rank);
    MPI_Comm_size(worker, &size);

    /* experiment */
    ofmo_storage_set_context(context);

    if (rank == 0) {
        rc = falanx_get_request_key_value(context, key, MAX_KEY_LENGTH, &len);
        assert(rc == FLX_SUCCESS);
        //dbg("start:key:[%s]:%zd\n", key, len);
    }
    MPI_Bcast(key, MAX_KEY_LENGTH, MPI_CHAR, 0, worker);

    /* ---- begin -----  */
    rc = ofmo_scf_input_read(&scfin, key, len);
    assert(rc == 0);
    /* ofmo_scf_input_show(&scfin); */ /* for debug */
    ofmo_storage_set_data_ids(scfin.data_ids);
    TLOG_LOG_IN(2+scfin.nmonomer);

    if (scfin.nmonomer > nbody || scfin.nmonomer == 0){
        ofmo_scf_output_en_zero(&scfouten); 
    } else {
        if (scfin.iscc != iscc_prev) {
            /* 新しいステップ */
            ofmo_set_new_scc_flag();
            ofmo_clear_cache();
    TLOG_LOG_EVENT(5);
            iscc_prev = scfin.iscc;
            if (scfin.nmonomer == 1 && fp_prof) {
                fprintf(fp_prof, "==== SCC itera= %d ====\n", scfin.iscc);
                fflush(fp_prof);
            }
        }
        tolscf = pow(0.1e0, (double)scfin.itol);
        if (scfin.nmonomer == 1) {
            /* 都度取りにいく */
            if (rank == 0) {
                ofmo_get_monomer_energy(scfin.monomer_list[0], &scfouten.energy);
                //dbg("before ofmo_calc_fragment_electronic_state:energy:%d:%lf\n",
                //    scfin.monomer_list[0], scfouten.energy);
            }
            MPI_Bcast(&scfouten.energy, 1, MPI_DOUBLE, 0, worker);
        }
        /* 計算 */
        //if (rank == 0) {
        //    dbg("ofmo_calc_fragment_electronic_state(con, comm, data_ids, %d, %d, %d, %lf, %lf, %lf, %lf, *)\n",
        //        scfin.nmonomer, scfin.monomer_list[0], scfin.method, tolscf,
        //        scfouten.energy,scfouten.energy0, scfouten.ddv);
        //}
        ofmo_calc_fragment_electronic_state(worker,
            scfin.nmonomer, scfin.monomer_list, scfin.method,
            tolscf,
            &scfouten.energy, &scfouten.energy0, &scfouten.ddv,
            &scfoutao->fnao, scfoutao->daopop, scfoutao->fsao2tuao,
            &scfoutat->fnat, scfoutat->datpop, scfoutat->fatom2tatom);
    }
    /* -----  end ----- */

    if (rank == 0) {
        char outkey[MAX_KEY_LENGTH];

        rc = ofmo_scf_outkey_make(outkey, MAX_KEY_LENGTH, &scfin);
        assert(rc == 0);

        ds = falanx_context_get_data_store(context);
        assert(ds != NULL);
        data_store_set(ds, outkey, strlen(outkey)+1, (char*)&scfouten.energy, sizeof(scfouten));

        if (!(scfin.nmonomer > nbody || scfin.nmonomer == 0)) {
            char outaokey[MAX_KEY_LENGTH];
            char outatkey[MAX_KEY_LENGTH];
            char* aobuf = NULL;
            char* atbuf = NULL;
            size_t aobuf_sz;
            size_t atbuf_sz;

            ofmo_scf_get_ao_pop_key(outkey, outaokey, MAX_KEY_LENGTH);
            ofmo_scf_get_at_pop_key(outkey, outatkey, MAX_KEY_LENGTH);

            rc = ofmo_scf_output_ao_pack(scfoutao, &aobuf, &aobuf_sz);
            assert(rc == 0);
            rc = ofmo_scf_output_at_pack(scfoutat, &atbuf, &atbuf_sz);
            assert(rc == 0);

            data_store_set(ds, outaokey, strlen(outaokey)+1, aobuf, aobuf_sz);
            data_store_set(ds, outatkey, strlen(outatkey)+1, atbuf, atbuf_sz);


            free(aobuf);
            free(atbuf);
        }

        falanx_set_response_key(context, outkey, strlen(outkey)+1);
        //dbg("end:%d:energy:%lf:energy0:%lf:ddv:%lf\n",
        //    scfin.monomer_list[0], scfouten.energy, scfouten.energy0, scfouten.ddv);
    }

    ofmo_scf_output_ao_free(scfoutao);
    ofmo_scf_output_at_free(scfoutat);
    ofmo_storage_set_context(NULL);
    ofmo_storage_set_data_ids(OFMO_STORAGE_DATA_IDS_NULL);
    TLOG_LOG_OUT(2+scfin.nmonomer);

}


/*
 * ESdimerタスク
 */
void
ofmo_approx_task(falanx_context_t* context, MPI_Comm worker)
{
    int rc;
    int rank;
    char inputkey[MAX_KEY_LENGTH];
    size_t keysize;
    data_store_t* ds = NULL;
    ofmo_approx_input_t aprin;
    double eesdimer0;
    int data_ids[6];

    static int nbody;
    static bool called = false;

    {
      int maxscc;
      ofmo_data_get_vals("maxscc", &maxscc);
      if (iscc_prev<=maxscc) {
        ofmo_clear_cache();
        iscc_prev = maxscc + 100;
      }
    }
    TLOG_LOG_IN(5);
    if (!called) {
        /* common block */
        ofmo_data_get_vals("nbody", &nbody);
        called = true;
    }

    MPI_Comm_rank(worker, &rank);

    ofmo_storage_set_context(context);
    ds = falanx_context_get_data_store(context);
    assert(ds != NULL);

    if (rank == 0) {
        rc = falanx_get_request_key_value(context, inputkey, MAX_KEY_LENGTH, &keysize);
        assert(rc == FLX_SUCCESS);
        ofmo_approx_inkey_parse(inputkey, keysize,
                                &aprin.data_id, &aprin.njob, &aprin.seq);
        data_store_get_value(ds, inputkey, strlen(inputkey)+1,
                            (char*)aprin.joblist, (sizeof(int) * MAXNJOB*2));
    }
    MPI_Bcast(&aprin, sizeof(ofmo_approx_input_t)/sizeof(int), MPI_INT, 0, worker);


    if (nbody < 2) {
        eesdimer0 = 0.e0;
    } else {
        if (aprin.njob == 0) {
            eesdimer0 = 0.e0;
        } else {
            data_ids[0] = aprin.data_id;
            ofmo_storage_set_data_ids(data_ids);
            ofmo_calc_es_dimer(worker, aprin.njob, aprin.joblist, &eesdimer0);
        }
    }

    if (rank == 0) {
        char outkey[MAX_KEY_LENGTH];
        ofmo_approx_outkey_make(outkey, MAX_KEY_LENGTH, &aprin);
        data_store_set(ds, outkey, strlen(outkey)+1, (char*)&eesdimer0, sizeof(double));
        falanx_set_response_key(context, outkey, strlen(outkey)+1);
    }

    ofmo_storage_set_context(NULL);
    ofmo_storage_set_data_ids(OFMO_STORAGE_DATA_IDS_NULL);
    TLOG_LOG_OUT(5);
}

