/*
   メモリサーバーに登録するデータ構造を決める関数群
   */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#include <unistd.h>

#include <mpi.h>

#include "ofmo-def.h"
#include "ofmo-data.h"
#include "ofmo-misc.h"

#define MAXTOKLEN MAXSTRLEN

/** ２つの整数を比較する関数 */
static int comp2( const void *p1, const void *p2 ) {
    return ( ((int*)p2)[0] - ((int*)p1)[0] );
}

/** 整数配列の最小値を持つ添字を返す */
static int min_loc( const int nprocs, const int sz[] ) {
    int loc, irank;
    int min;
    min = sz[0];
    loc = 0;
    for ( irank=1; irank<nprocs; irank++ )
	if ( sz[irank] < min ) loc = irank;
    return loc;
}

/* nfrag個のデータをnprocsプロセスにできるだけ公平に分割する */
/* 引数
 * nfrag (in) : データ数
 * nelem[ifrag] (in) : データ番号ifragの要素数
 * nprocs (in) : プロセス数
 *
 * target[ifrag] (out) : ifrag番目のデータの格納先（ターゲットプロセス番号）
 * offset[ifrag] (out) : ifrag番目のデータのoffset
 * total_nelem[irank] (out) : irankプロセスに割り当てられた要素数
 *
 * *iwork (in) : データ操作のワーク領域（nfrag*2個の要素が必要）
 * */
static int divide_n_data_by_p_proc(
	const int nfrag, const int nelem[], const int nprocs, 
	int target[], int offset[], int total_nelem[],
	int *iwork ) {
    int ifrag2, ifrag, id, id2, elems, irank;
    for ( ifrag=0, ifrag2=0; ifrag<nfrag; ifrag++, ifrag2+=2 ) { 
	iwork[ifrag2+0] = nelem[ifrag];
	iwork[ifrag2+1] = ifrag;
    }
    qsort( iwork, nfrag, 2*sizeof(int), comp2 );
    for ( irank=0; irank<nprocs; irank++ ) total_nelem[irank] = 0;
    for ( id=0, id2=0; id<nfrag; id++, id2+=2 ) {
	elems = iwork[id2+0];
	ifrag = iwork[id2+1];
	irank = min_loc( nprocs, total_nelem );
	offset[ifrag]       = total_nelem[irank];
	target[ifrag]       = irank;
	total_nelem[irank] += elems;
    }
    return 0;
}

static int NDATA = 0;		// データ塊の数
static int NFRAG = 0;
static int **DATA_TARGET = NULL;	// 各データのターゲットプロセス
static int **DATA_NELEMS = NULL;
static int **DATA_OFFSET = NULL;

static void dealloc() {
    ofmo_free_imatrix( DATA_TARGET );
    ofmo_free_imatrix( DATA_NELEMS );
    ofmo_free_imatrix( DATA_OFFSET );
    NDATA = 0;
    NFRAG = 0;
}

static int alloc( int ndata, int nfrag ) {
    static int called = false;
    if ( called ) return 0;
    DATA_TARGET = ofmo_alloc_imatrix( ndata, nfrag );
    DATA_NELEMS = ofmo_alloc_imatrix( ndata, nfrag );
    DATA_OFFSET = ofmo_alloc_imatrix( ndata, nfrag );
    NDATA = ndata;
    NFRAG = nfrag;
    atexit( dealloc );
    called = true;
    return 0;
}

static int invalid_data_id( int data_id ) {
    if ( data_id < 0 || data_id >= NDATA ) return true;
    return false;
}

static int invalid_frag( int data_id, int ifrag ) {
    if ( ifrag < 0 || ifrag >= NFRAG ) return true;
    return false;
}

int* ofmo_getadd_data_target( int data_id ) {
    if ( invalid_data_id( data_id ) ) return NULL;
    return DATA_TARGET[data_id];
}

int* ofmo_getadd_data_nelems( int data_id ) {
    if ( invalid_data_id( data_id ) ) return NULL;
    return DATA_NELEMS[data_id];
}

int* ofmo_getadd_data_offset( int data_id  ) {
    if ( invalid_data_id( data_id ) ) return NULL;
    return DATA_OFFSET[data_id];
}

int ofmo_get_data_target( int data_id, int ifrag ) {
    if ( invalid_data_id( data_id ) ) return -1;
    if ( invalid_frag( data_id, ifrag ) ) return -1;
    return DATA_TARGET[data_id][ifrag];
}


int ofmo_get_data_nelems( int data_id, int ifrag ) {
    if ( invalid_data_id( data_id ) ) return -1;
    if ( invalid_frag( data_id, ifrag ) ) return -1;
    return DATA_NELEMS[data_id][ifrag];
}

int ofmo_get_data_offset( int data_id, int ifrag ) {
    if ( invalid_data_id( data_id ) ) return -1;
    if ( invalid_frag( data_id, ifrag ) ) return -1;
    return DATA_OFFSET[data_id][ifrag];
}

int ofmo_get_nfrag() { return NFRAG; }
int ofmo_get_ndata() { return NDATA; }

/* メモリサーバーへのデータの格納方法を決める
    2013/02/13現在、以下の５つのデータのメモリサーバーへの保存を予定
    0. OFMO_DENS1	モノマー密度行列データ１
    1. OFMO_DENS2	モノマー密度行列データ２
    2. OFMO_AOPOP1	モノマーのAO population
    3. OFMO_AOPOP2	モノマーのAO population
    4. OFMO_ATPOP1	モノマーのatomic population
    5. OFMO_ATPOP2	モノマーのatomic population
    6. OFMO_DISTA	モノマー間距離
    7. OFMO_ENERGY	モノマーエネルギー
    8. OFMO_ENERGY0	環境ポテンシャル抜きのモノマーエネルギー
    9. OFMO_TOTAL_AOPOP	分子全体のAO population
   10. OFMO_TOTAL_ATPOP	分子全体のatomic population

    注意点：
    1. すべてのデータは倍精度実数型である
	（他のデータ型には対応していない）
    2. AO populationとatomic populationのデータは１つのプロセスに
       保存されているため、全データを一度に取り出すことが容易
    3. すべて、要素数はnfragである

    */
int ofmo_make_data_struct( int ndata, int nfrag, int mserv_size ) {
    int *nfao, *nfatom;
    int *iwork, *subtot_nelem, *total_nelems;
    int data_id, *nelems, *target, *offset;
    //
    int nao, irank, ifrag, total, tot_nao, tot_natom;
    ofmo_data_get_vals("nfao nfatom natom nao",
	    &nfao, &nfatom, &tot_natom, &tot_nao );
    alloc( ndata, nfrag );
    // 一時配列の確保
    iwork  = (int*)malloc( sizeof(int) * nfrag * 2 );
    subtot_nelem = (int*)malloc( sizeof(int) * mserv_size );
    total_nelems = (int*)malloc( sizeof(int) * mserv_size );
    for ( irank=0; irank<mserv_size; irank++ ) total_nelems[irank] = 0;

    data_id=0;
    // モノマー密度行列でーた１
    nelems = ofmo_getadd_data_nelems( data_id );
    target = ofmo_getadd_data_target( data_id );
    offset = ofmo_getadd_data_offset( data_id );
    for ( ifrag=0; ifrag<nfrag; ifrag++ ) {
	nao = nfao[ifrag];
	nelems[ifrag] = (nao*nao+nao)>>1;
    }
    divide_n_data_by_p_proc( nfrag, nelems, mserv_size, target, offset,
	    subtot_nelem, iwork );
    for ( irank=0; irank<mserv_size; irank++ )
	total_nelems[irank] += subtot_nelem[irank];
    data_id++;
    // モノマー密度行列データ２
    nelems = ofmo_getadd_data_nelems( data_id );
    target = ofmo_getadd_data_target( data_id );
    offset = ofmo_getadd_data_offset( data_id );
    for ( ifrag=0; ifrag<nfrag; ifrag++ ) {
	nao = nfao[ifrag];
	nelems[ifrag] = (nao*nao+nao)>>1;
    }
    divide_n_data_by_p_proc( nfrag, nelems, mserv_size, target, offset,
	    subtot_nelem, iwork );
    for ( irank=0; irank<mserv_size; irank++ )
	total_nelems[irank] += subtot_nelem[irank];
    data_id++;
    // AO population１
    nelems = ofmo_getadd_data_nelems( data_id );
    target = ofmo_getadd_data_target( data_id );
    offset = ofmo_getadd_data_offset( data_id );
    for ( ifrag=0; ifrag<nfrag; ifrag++ ) nelems[ifrag] = nfao[ifrag];
    irank = min_loc( mserv_size, total_nelems );
    total = 0;
    for ( ifrag=0; ifrag<nfrag; ifrag++ ) {
	offset[ifrag] = total;
	target[ifrag] = irank;
	total        += nelems[ifrag];
    }
    total_nelems[irank] += total;
    data_id++;
    // AO population２
    nelems = ofmo_getadd_data_nelems( data_id );
    target = ofmo_getadd_data_target( data_id );
    offset = ofmo_getadd_data_offset( data_id );
    for ( ifrag=0; ifrag<nfrag; ifrag++ ) nelems[ifrag] = nfao[ifrag];
    irank = min_loc( mserv_size, total_nelems );
    total = 0;
    for ( ifrag=0; ifrag<nfrag; ifrag++ ) {
	offset[ifrag] = total;
	target[ifrag] = irank;
	total        += nelems[ifrag];
    }
    total_nelems[irank] += total;
    data_id++;
    // Atomic population１
    nelems = ofmo_getadd_data_nelems( data_id );
    target = ofmo_getadd_data_target( data_id );
    offset = ofmo_getadd_data_offset( data_id );
    for ( ifrag=0; ifrag<nfrag; ifrag++ ) nelems[ifrag] = nfatom[ifrag];
    irank = min_loc( mserv_size, total_nelems );
    total = 0;
    for ( ifrag=0; ifrag<nfrag; ifrag++ ) {
	offset[ifrag] = total;
	target[ifrag] = irank;
	total        += nelems[ifrag];
    }
    total_nelems[irank] += total;
    data_id++;
    // Atomic population２
    nelems = ofmo_getadd_data_nelems( data_id );
    target = ofmo_getadd_data_target( data_id );
    offset = ofmo_getadd_data_offset( data_id );
    for ( ifrag=0; ifrag<nfrag; ifrag++ ) nelems[ifrag] = nfatom[ifrag];
    irank = min_loc( mserv_size, total_nelems );
    total = 0;
    for ( ifrag=0; ifrag<nfrag; ifrag++ ) {
	offset[ifrag] = total;
	target[ifrag] = irank;
	total        += nelems[ifrag];
    }
    total_nelems[irank] += total;
    data_id++;
    // フラグメント間距離
    nelems = ofmo_getadd_data_nelems( data_id );
    target = ofmo_getadd_data_target( data_id );
    offset = ofmo_getadd_data_offset( data_id );
    for ( ifrag=0; ifrag<nfrag; ifrag++ ) nelems[ifrag] = nfrag;
    divide_n_data_by_p_proc( nfrag, nelems, mserv_size, target, offset,
	    subtot_nelem, iwork );
    for ( irank=0; irank<mserv_size; irank++ )
	total_nelems[irank] += subtot_nelem[irank];
    data_id++;
    // モノマーのエネルギー（環境ポテンシャル項あり）
    nelems = ofmo_getadd_data_nelems( data_id );
    target = ofmo_getadd_data_target( data_id );
    offset = ofmo_getadd_data_offset( data_id );
    for ( ifrag=0; ifrag<nfrag; ifrag++ ) nelems[ifrag] = 1;
    irank = min_loc( mserv_size, total_nelems );
    total = 0;
    for ( ifrag=0; ifrag<nfrag; ifrag++ ) {
	offset[ifrag] = total;
	target[ifrag] = irank;
	total        += nelems[ifrag];
    }
    total_nelems[irank] += total;
    data_id++;
    // モノマーのエネルギー（環境ポテンシャル項なし）
    nelems = ofmo_getadd_data_nelems( data_id );
    target = ofmo_getadd_data_target( data_id );
    offset = ofmo_getadd_data_offset( data_id );
    for ( ifrag=0; ifrag<nfrag; ifrag++ ) nelems[ifrag] = 1;
    irank = min_loc( mserv_size, total_nelems );
    total = 0;
    for ( ifrag=0; ifrag<nfrag; ifrag++ ) {
	offset[ifrag] = total;
	target[ifrag] = irank;
	total        += nelems[ifrag];
    }
    total_nelems[irank] += total;
    data_id++;
    // 全体のAO population
    nelems = ofmo_getadd_data_nelems( data_id );
    target = ofmo_getadd_data_target( data_id );
    offset = ofmo_getadd_data_offset( data_id );
    for ( ifrag=1; ifrag<nfrag; ifrag++ ) nelems[ifrag] = 0;
    nelems[0] = tot_nao;
    irank = min_loc( mserv_size, total_nelems );
    total = 0;
    for ( ifrag=0; ifrag<nfrag; ifrag++ ) {
	offset[ifrag] = total;
	target[ifrag] = irank;
	total        += nelems[ifrag];
    }
    total_nelems[irank] += total;
    data_id++;
    // 全体のAtomic population
    nelems = ofmo_getadd_data_nelems( data_id );
    target = ofmo_getadd_data_target( data_id );
    offset = ofmo_getadd_data_offset( data_id );
    for ( ifrag=1; ifrag<nfrag; ifrag++ ) nelems[ifrag] = 0;
    nelems[0] = tot_natom;
    irank = min_loc( mserv_size, total_nelems );
    total = 0;
    for ( ifrag=0; ifrag<nfrag; ifrag++ ) {
	offset[ifrag] = total;
	target[ifrag] = irank;
	total        += nelems[ifrag];
    }
    total_nelems[irank] += total;
    data_id++;

    Free( iwork );
    Free( subtot_nelem );
    Free( total_nelems );
    return 0;
}
