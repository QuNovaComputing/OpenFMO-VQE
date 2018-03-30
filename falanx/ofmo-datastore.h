#ifndef _OFMO_DS_PUT_GET_MASTER_H_
#define _OFMO_DS_PUT_GET_MASTER_H_

#include <stdio.h>
#include <mpi.h>

#include <falanx.h>
#include <datastore.h>

#include "ofmo-falanx-main.h"
#include "ofmo-task-util.h"

extern int ofmo_read_and_put_all(ofmo_master_config_t* config,
                                 falanx_context_t*     context,
                                 char*                 fileheader);

extern int ofmo_get_and_write(falanx_context_t* context, char *fileheader, int nfrag );
extern int ofmo_get_monomer_aopop_master(data_store_t* ds , int ifrag, double aopop[] );
extern int ofmo_get_monomer_atpop_master(data_store_t* ds , int ifrag, double atpop[] );

extern int ofmo_update_monomer_data_master(void);
extern int ofmo_master_get(data_store_t* ds, int data_id, int ifrag, double outBuff[]);
extern int ofmo_master_put(data_store_t* ds, int data_id, int ifrag, double* src);


extern void ofmo_get_monomer_data_master(int data_ids[6]);

extern int ofmo_master_get_val(data_store_t* ds, int data_id, int ifrag, char* outBuff, size_t buffSize);

/*
 * OFMO_ENERGY/ENERGY0 データをすべて取り出す.
 *
 * @param ds データストアクライアントハンドル
 * @param data_id データタイプ(OFMO_ENERGYもしくはOFMO_ENERGY0)
 * @param buff 出力バッファへのポインタ
 * @param nfrag 要素数(フラグメント数)
 */
extern int ofmo_master_get_energy_all(data_store_t* ds, const int data_id, double* outBuff, int nfrag);

extern int ofmo_master_get_scf_ao_pop(data_store_t* ds, char* key, size_t keysize, ofmo_scf_output_ao_t* out);
extern int ofmo_master_get_scf_at_pop(data_store_t* ds, char* key, size_t keysize, ofmo_scf_output_at_t* out);


#endif

