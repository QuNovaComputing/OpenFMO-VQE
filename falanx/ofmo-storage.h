
#ifndef _OFMO_STORAGE_H_
#define _OFMO_STORAGE_H_

/*
 * オリジナルOpenFMOの関数に必要なデータをこのオブジェクトファイルを経由して渡す.
 */

#include <falanx.h>

#define OFMO_STORAGE_DATA_IDS_NULL g_ofmo_storage_data_ids_null
extern const int const g_ofmo_storage_data_ids_null[6];

extern void ofmo_storage_set_context(falanx_context_t* context);
extern falanx_context_t* ofmo_storage_get_context(void);

extern void ofmo_storage_set_data_ids(const int data_ids[6]);
extern void ofmo_storage_get_data_ids(int data_ids[6]);

#endif

