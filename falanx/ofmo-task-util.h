
#ifndef _OFMO_SCF_UTIL_H_
#define _OFMO_SCF_UTIL_H_

#include <stdlib.h>

#include "ofmo-def.h"

#ifndef UNUSED
#define UNUSED(X) ((void)X)
#endif

#define MAX_KEY_LENGTH 256

/*
 * SCFタスクの入力
 */
struct ofmo_scf_input_st {
    int method;
    int nmonomer;
    int iscc;
    int itol;
    int data_ids[6];
    int monomer_list[3];
};
typedef struct ofmo_scf_input_st ofmo_scf_input_t;

/*
 * SCFタスクの出力用構造体
 */
struct ofmo_scf_output_en_st {
    double energy;
    double energy0;
    double ddv;
};
typedef struct ofmo_scf_output_en_st ofmo_scf_output_en_t;

/*
 * SCFタスクの出力用構造体
 */
struct ofmo_scf_output_ao_st {
    int     fnao;
    double* daopop;
    int*    fsao2tuao;
};
typedef struct ofmo_scf_output_ao_st ofmo_scf_output_ao_t;

/*
 * SCFタスクの出力用構造体
 */
struct ofmo_scf_output_at_st {
    int     fnat;
    double* datpop;
    int*    fatom2tatom;
};
typedef struct ofmo_scf_output_at_st ofmo_scf_output_at_t;

/*
 * APPROXタスク用の入力
 */
struct ofmo_approx_input_st {
    int seq;
    int data_id;            /* データタイプ (dold) */
    int njob;               /* ジョブ数 */
    int joblist[MAXNJOB*2]; /* ジョブリスト(固定長) */
};
typedef struct ofmo_approx_input_st ofmo_approx_input_t;

#define is_valid_data_id(IDS) \
    ( (((IDS[0]) == OFMO_DENS1)  || ((IDS[0]) == OFMO_DENS2))  && \
      (((IDS[1]) == OFMO_DENS1)  || ((IDS[1]) == OFMO_DENS2))  && \
      (((IDS[2]) == OFMO_AOPOP1) || ((IDS[2]) == OFMO_AOPOP2)) && \
      (((IDS[3]) == OFMO_AOPOP1) || ((IDS[3]) == OFMO_AOPOP2)) && \
      (((IDS[4]) == OFMO_ATPOP1) || ((IDS[4]) == OFMO_ATPOP2)) && \
      (((IDS[5]) == OFMO_ATPOP1) || ((IDS[5]) == OFMO_ATPOP2)) )

/*
 * 初期密度行列計算用入力キーを生成する.
 */
int ofmo_init_dens_inkey_make(char* key, size_t size, int ifrag, int data_ids[6]);


void ofmo_scf_output_en_zero(ofmo_scf_output_en_t* out);

/*
 * 文字列形式の入力から ofmo_scf_input_t構造体のフィールドを埋める
 */
int ofmo_scf_input_read(ofmo_scf_input_t* in, char* key, size_t len);

void ofmo_scf_input_show(ofmo_scf_input_t* scfin);

/*
 * OFMO_SCF の入力キーを生成する関数.
 */
int ofmo_scf_inkey_make(char* inputkey, size_t size,
                        int method, int nmon, int scc, int conv,
                        int mon1, int mon2, int mon3,
                        int data_ids[6]);

/*
 * OFMO_SCF の出力キーを生成する関数.
 */
int ofmo_scf_outkey_make(char* outkey, size_t size, ofmo_scf_input_t* scfin);

/*
 * OFMO_SCF の出力キーを解析する
 *
 * 出力変数は全部optional.
 * outMonomerListを指定する場合は、int[3]の配列へのポインタを渡す.
 */
int ofmo_scf_outkey_parse(char* outkey, size_t size,
                          int* outMethod, int* outNmon, int* outScc, int* outMonomerList);

ofmo_scf_output_ao_t* ofmo_scf_output_ao_new(int nbody, int maxnfao);
ofmo_scf_output_at_t* ofmo_scf_output_at_new(int nbody, int maxnfatom);
void ofmo_scf_output_ao_free(ofmo_scf_output_ao_t* scfout);
void ofmo_scf_output_at_free(ofmo_scf_output_at_t* scfout);

/*
 * Pack関数
 *
 * @param [in] scfout 入力用構造体
 * @param [out] outBuf Pack先データ
 * @param [out] outSize outBufサイズ
 */
int ofmo_scf_output_ao_pack(ofmo_scf_output_ao_t* scfout, char** outBuf, size_t* outSize);
int ofmo_scf_output_at_pack(ofmo_scf_output_at_t* scfout, char** outBuf, size_t* outSize);

/*
 * Unpack関数
 *
 * @param [in] data Falan Data StoreSから取り出したデータ
 * @param [in] size dataサイズ
 * @param [out] scfout 出力先構造体. 適切な領域が割り当てられていること.
 */
int ofmo_scf_output_ao_unpack(char* data, size_t size, ofmo_scf_output_ao_t* scfout);
int ofmo_scf_output_at_unpack(char* data, size_t size, ofmo_scf_output_at_t* scfout);

/*
 * SCFタスクの出力キーからAO populationデータ用のキーを用意する.
 *
 * @param [in] key SCFタスクの出力キー
 * @param [out] outKey 出力先.
 * @param [in] size 出力先のサイズ.
 */
int ofmo_scf_get_ao_pop_key(const char* key, char* outKey, size_t size);
int ofmo_scf_get_at_pop_key(const char* key, char* outKey, size_t size);

/*
 * Approxタスクの入力キーを生成する
 */
int ofmo_approx_inkey_make(char* outkey, size_t size, int data_id, int njob, int seq);

/*
 * Approxタスクの出力キーを生成する
 */
int ofmo_approx_outkey_make(char* outkey, size_t size, ofmo_approx_input_t* aprin);

/*
 * Approxタスクの入力キーを解析する.
 */
int ofmo_approx_inkey_parse(char* inputkey, size_t size, int* outDataId, int* outNjob, int* outSeq);

/*
 * Dimer SCFタスクの入力キーかテストする.
 * @return dimer SCFタスクの入力キーなら True、そうでなければFalseを返す.
 */
int ofmo_is_dimer_scf_input_key(char* key, size_t size);

/*
 * Dimer SCFタスクの出力キーかテストする.
 * @return dimer SCFタスクの出力キーなら True、そうでなければFalseを返す.
 */
int ofmo_is_dimer_scf_output_key(char* key, size_t size);

/*
 * Approxタスクの入力キーかテストする.
 * @return Approxタスクの入力キーなら True、そうでなければFalseを返す.
 */
int ofmo_is_approx_input_key(char* key, size_t size);

/*
 * Approxタスクの出力キーかテストする.
 * @return Approxタスクの出力キーなら True、そうでなければFalseを返す.
 */
int ofmo_is_approx_output_key(char* key, size_t size);

#endif

