/**
 * @file ofmo-data.c
 * FMOで扱う各種データを保存して、参照、更新するための関数群
 *
 * FMO計算で必要となるほとんどの情報は、このファイルで定義された
 * 関数群を通してアクセスできる。
 *
 * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "ofmo-def.h"

#define POINTER 128

#define INT	1
#define CHAR	2
#define DBLE	3
#define LONG	4

#define P_INT		INT + POINTER
#define P_CHAR		CHAR + POINTER
#define P_DBLE		DBLE + POINTER
#define P_LONG		LONG + POINTER

#define PP_INT		P_INT + POINTER
#define PP_CHAR		P_CHAR + POINTER
#define PP_DBLE		P_DBLE + POINTER
#define PP_LONG		P_LONG + POINTER

static char *typename[] = {
    "undefined", "int", "char", "double", "long",
};

struct _val_carrier_ {
    char*	identifier;
    int 	datatype;
    void*	val;
    char*	explanation;
};

static struct _val_carrier_ vc[] = {
    /* ------------ data from input file ------------- */
    // in $contrl section
    {"maxscf", INT, NULL, "maximum SCF cycles"},
    {"nprint", INT, NULL, "output level"},
    {"itol", INT, NULL, "primitive cutoff factor ( 10^{-n} )"},
    {"icut", INT, NULL, "integrals less than 10^{-n} are ignored"},
    // in $basis section
    {"ngauss", INT, NULL, "# of Gaussians (N)"},
    {"diffsp", INT, NULL, "flag of SP diffuse function for heavy atom"},
    {"diffs", INT, NULL, "flag of S diffuse function for H atom"},
    {"ndfunc", INT, NULL, "# of d-functions for heavy atom"},
    {"npfunc", INT, NULL, "# of p-functions for H atom"},
    {"nffunc", INT, NULL, "# of f-functions for heavy atom"},
    {"gbasis", P_CHAR, NULL, "name of basis set categoly"},
    {"atombs", P_INT, NULL, "basis set ID for each atom"},
    // in $intgrl section
    {"nintmx", LONG, NULL, "maximum # of integrals in a record block"},
    {"nintic", LONG, NULL, "size of buffer memory in bytes"},
    // in $gddi section
    {"ngroup", INT, NULL, "# of worker groups"},
    {"niogroup", INT, NULL, "# of I/O groups"},
    {"nioprocs", INT, NULL, "# of I/O processes"},
    {"master", INT, NULL, "rank of master process"},
    {"iopos", INT, NULL, "position of I/O processes in I/O group"},
    // in $scf section
    {"diis", INT, NULL, "flag whether DIIS is to be used or not"},
    {"npunch", INT, NULL, "output level after SCF"},
    {"scfconv", INT, NULL, "convergence criterion of SCF(10^{-n})"},
    // in $fmoxyz section
    {"natom", INT, NULL, "# of atoms in molecule"},
    {"atn", P_INT, NULL, "atomic number of each atom in molecule"},
    {"atx", P_DBLE, NULL, "x-coodinates of each atom in molecule"},
    {"aty", P_DBLE, NULL, "y-coodinates of each atom in molecule"},
    {"atz", P_DBLE, NULL, "z-coodinates of each atom in molecule"},
    // in $fmo section
    {"nfrag", INT, NULL, "# of fragments in molecule"},
    {"icharg", P_INT, NULL, "net charge of each monomer"},
    {"at2frg", P_INT, NULL, "fragment number in which each atom is"},
    {"nbody", INT, NULL, "nbody approximation is used"},
    {"lptc", DBLE, NULL, "threshold for Mulliken point charge approx."},
    {"laop", DBLE, NULL, "threshold for AO population approx."},
    {"ldim", DBLE, NULL, "threshold for separated dimer approx."},
    // in $fmolmo section
    {"nlmo", INT, NULL, "# of hybrid orbitals in bond detached atom"},
    {"nfunc", INT, NULL, "# of basis functions in hybrid orbital"},
    {"clmo", P_DBLE, NULL, "coefficients of hybrid orbitals"},
    {"ban", P_INT, NULL, "bond assignment number"},
    {"lmobs", P_CHAR, NULL, "basis name of LMO"},
    // in $fmobnd section
    {"nbond",  INT, NULL, "# of detached bonds"},
    {"welec", P_INT, NULL, "SN of detached atom with electron"},
    {"woelec", P_INT, NULL, "SN of detached atom w/o electron"},
    // in $fmoprp
    {"sccconv",  INT, NULL, "SCC convergence criteria"},
    {"maxscc",  INT, NULL, "maximum number of SCC cycle"},
    /* ================ data derived from input data =============== */
    {"nfatom", P_INT, NULL, "# of atoms in each monomer"},
    {"nfbond", P_INT, NULL, "# of BDA in each monomer"},
    {"maxnfatom", INT, NULL, "maximum # of atoms in monomer"},
    {"maxnfbond", INT, NULL, "maximum # of BDA in monomer"},
    {"ifatom", PP_INT, NULL, "SN of atom for i-th monomer, j-th atom"},
    {"matn", PP_INT, NULL, "atomic number of atom in each monomer"},
    {"fbsn1", PP_INT, NULL, "total SN of BDA in SN of BDA in monomer"},
    {"fbda", PP_INT, NULL, "monomer SN of BDA"},
    {"matx", PP_DBLE, NULL, "x-coordinate of atoms in each monomer"},
    {"maty", PP_DBLE, NULL, "y-coordinate of atoms in each monomer"},
    {"matz", PP_DBLE, NULL, "z-coordinate of atoms in each monomer"},
    /* 分子全体の基底関数データ */
    {"nsbs", INT, NULL, "# of kind of basis set used in calculation"},
    {"bslst", PP_CHAR, NULL, "list of basis set names"},
    {"maxlqn", INT, NULL, "maximum orbital quantum number"},
    {"ncs", INT, NULL, "# of CS of entire molecule"},
    {"nao", INT, NULL, "# of AO of entire molecule"},
    {"nps", INT, NULL, "# of PS of entire molecule"},
    {"ushel_lqn", P_INT, NULL, "orbital quantum number of CS"},
    {"ushel_ini", P_INT, NULL, "leading AO number of CS"},
    {"atm_lcs", P_INT, NULL, "leading CS number of atom"},
    /* モノマー基底関数に関するデータ */
    {"maxnfcs", INT, NULL, "maximum number of CS in each monomer"},
    {"maxnfao", INT, NULL, "maximum number of AO in each monomer"},
    {"maxnfps", INT, NULL, "maximum number of PS in each monomer"},
    {"maxnpspair", INT, NULL, "maximum number of PS pair"},
    {"nfcs", P_INT, NULL, "# of CS in each monomer"},
    {"nfao", P_INT, NULL, "# of AO in each monomer"},
    {"nfps", P_INT, NULL, "# of PS in each monomer"},
    {"mlcs", PP_INT, NULL, "leading CS of each CS type"},
    {"mshel_tem", PP_INT, NULL, "orbital quantum number for each monomer"},
    {"mshel_atm", PP_INT, NULL, "SN of atom which target CS is belonging"},
    {"mshel_add", PP_INT, NULL, "offset of leading PS"},
    {"mshel_ini", PP_INT, NULL, "leading AO number of each monomer"},
    {"mprim_exp", PP_DBLE, NULL, "orbital exponent"},
    {"mprim_coe", PP_DBLE, NULL, "contraction coef. with norm. coef."},
    /* モノマーAO変換テーブル */
    {"msao2muao", PP_INT, NULL, "monomer sorted AO to monomer unsorted AO"},
    {"msao2tuao", PP_INT, NULL, "monomer sorted AO to total unsorted AO"},
    {"muao2tuao", PP_INT, NULL, "monomer unsorted AO to total unsorted AO"},
    /* モノマーpopulationデータ */
    {"maopop", PP_DBLE, NULL, "AO population data for monomers"},
    {"matpop", PP_DBLE, NULL, "atomic population data for monomers"},
    /* フラグメント（モノマー以外）のデータ */
    {"fnatom", INT, NULL, "# of atom in fragment"},
    {"fncs", INT, NULL, "# of CS in fragment"},
    {"fnao", INT, NULL, "# of AO in fragment"},
    {"fnps", INT, NULL, "# of PS in fragment"},
    {"fatn", P_INT, NULL, "atomic number of atoms in fragment"},
    {"fat2tat", P_INT, NULL, "fragment atom to entire molecule"},
    {"fatx", P_DBLE, NULL, "x-coordinate of atoms in fragment"},
    {"faty", P_DBLE, NULL, "y-coordinate of atoms in fragment"},
    {"fatz", P_DBLE, NULL, "z-coordinate of atoms in fragment"},
    {"flcs", P_INT, NULL, "leading CS of each CS type for fragment"},
    {"fshel_tem", P_INT, NULL, "orbital quantum number for fragment"},
    {"fshel_atm", P_INT, NULL, "SN of atom which target CS is belonging"},
    {"fshel_add", P_INT, NULL, "offset of leading PS"},
    {"fshel_ini", P_INT, NULL, "leading AO number of fragment"},
    {"fprim_exp", P_DBLE, NULL, "orbital exponent"},
    {"fprim_coe", P_DBLE, NULL, "contraction coef. with norm. coef."},
    /* フラグメントAO変換テーブル */
    {"fsao2fuao", P_INT, NULL, "frag. sorted AO to frag. unsorted AO"},
    {"fsao2tuao", P_INT, NULL, "frag. sorted AO to total unsorted AO"},
    {"fuao2tuao", P_INT, NULL, "frag. unsorted AO to total unsorted AO"},
    {"dold", INT, NULL, "storage ID of old monomer density matrices"},
    {"dnew", INT, NULL, "storage ID of old monomer density matrices"},
    /* フラグメントのAO population, Atomic population（の変位）*/
    {"daopop", P_DBLE, NULL, "delta AO pop. of fragment"},
    {"datpop", P_DBLE, NULL, "delta atomic pop. of fragment"},
};

static int NCOUNT = 0;	// # of variables dealt with this class

static void finalize_data() {
    int id, dtype;
    for ( id=0; id<NCOUNT; id++ ) {
	dtype = vc[id].datatype;
	if        ( dtype < 2*POINTER ) {
	    if ( vc[id].val != NULL ) free( vc[id].val );
	} else if ( dtype < 3*POINTER ) {
	    if ( vc[id].val != NULL ) {
		if        ( dtype == PP_INT ) {
		    if ( ((int**)vc[id].val)[0] != NULL )
			free( ((int**)vc[id].val)[0] );
		} else if ( dtype == PP_CHAR ) {
		    if ( ((char**)vc[id].val)[0] != NULL )
			free( ((char**)vc[id].val)[0] );
		} else if ( dtype == PP_DBLE ) {
		    if ( ((double**)vc[id].val)[0] != NULL )
			free( ((double**)vc[id].val)[0] );
		} else if ( dtype == PP_LONG ) {
		    if ( ((long**)vc[id].val)[0] != NULL )
			free( ((long**)vc[id].val)[0] );
		}
		free( vc[id].val );
	    }
	}
	vc[id].val = NULL;
    }
}

static int init_data() {
    static int called = false;
    int id, dtype;
    if ( called ) return 0;
    NCOUNT = sizeof(vc) / sizeof(struct _val_carrier_);
    for ( id=0; id<NCOUNT; id++ ) {
	dtype = vc[id].datatype;
	if ( dtype < POINTER ) {
	    if ( dtype == INT ) {
		vc[id].val = (int*)malloc( sizeof(int) );
		*((int*)vc[id].val) = -1;
	    } else if ( dtype == CHAR ) {
		vc[id].val = (char*)malloc( sizeof(char) );
		*((char*)vc[id].val) = '\0';
	    } else if ( dtype == DBLE ) {
		vc[id].val = (double*)malloc( sizeof(double) );
		*((double*)vc[id].val) = 0.e0;
	    } else if ( dtype == LONG ) {
		vc[id].val = (long*)malloc( sizeof(long) );
		*((long*)vc[id].val) = -1;
	    } else {
		printf("Illegal datatype (%d) in %d-th element\n",
			dtype, (id+1) );
	    }
	} else {
	    vc[id].val = NULL;
	}
    }
    atexit( finalize_data );
    called = true;
    return 0;
}

// 与えられた識別子と一致する変数のIDを返す
// 正常終了時＝非負整数
// 異常終了時＝-1
static int get_id( const char *name ) {
    int id;
    init_data();
    for ( id=0; id<NCOUNT; id++ ) {
	if ( strcmp( name, vc[id].identifier ) == 0 ) return id;
    }
    return -1;
}

/** 変数の一覧を出力する
 *
 * このファイルで扱う各種変数の情報を表示する。スカラ変数は現在の値を、
 * また、配列変数についてはポインタの値を表示する。
 *
 * @note デバッグ用の関数として準備している
 *
 * @ingroup ofmo-input
 *
 * */
int ofmo_data_show_all() {
    init_data();
    int id, dtype;
    char stype[MAXSTRLEN];
    printf("======================= all data =========================\n");
    printf(" # of variables = %d\n", NCOUNT );
    printf("----------------------------------------------------------\n");
    printf(" %10s %11s %12s : %s\n",
	    "name", "data type", "val", "explanation" );
    printf("----------------------------------------------------------\n");
    for ( id=0; id<NCOUNT; id++ ) {
	dtype = vc[id].datatype;
	// data typeを表す文字列
	if ( dtype > 2*POINTER )
	    sprintf( stype, "**%s", typename[dtype-2*POINTER] );
	else if ( dtype > POINTER )
	    sprintf( stype, "*%s", typename[dtype-POINTER] );
	else if ( dtype < POINTER && dtype > 0 )
	    sprintf( stype, "%s", typename[dtype] );
	else
	    sprintf( stype, "%s", "undefined" );
	// output
	printf(" %10s %11s ", vc[id].identifier, stype );
	if ( dtype > POINTER ) printf("%12p", vc[id].val );
	else if ( dtype == INT ) printf("%12d", *((int*)vc[id].val) );
	else if ( dtype == CHAR )
	    printf("%11s%c", " ", *((char*)vc[id].val) );
	else if ( dtype == DBLE ) printf("%12f", *((double*)vc[id].val) );
	else if ( dtype == LONG ) printf("%12ld", *((long*)vc[id].val) );
	printf(" : %s\n", vc[id].explanation );
    }
    printf("----------------------------------------------------------\n");
    return 0;
}

/** この関数群で保持している値を参照するための関数
 *
 * @param[in] format 参照したい変数名を空白区切りで列挙した文字列
 * @param[out] ... 参照する値が代入される変数のポインタたち
 *
 * @note \c format に書かれた文字列の順に、変数（のポインタ）を並べる
 * 必要がある
 *
 * @ingroup ofmo-input
 * */
int ofmo_data_get_vals( const char* format, ...) {
    int narg, len, id, dtype;
    char *p, *sep = " \t\n\r", strtmp[MAXSTRLEN];
    va_list parg;
    init_data();
    va_start(parg, format);
    len = strlen(format);
    if (len >= (MAXSTRLEN-1)) {
        dbg("Too long length of \'format\' string!\n");
        va_end(parg);
        return -1;
    }
    strncpy(strtmp, format, len);
    strtmp[len] = '\0';
    narg = 0;
    p = strtok(strtmp, sep);
    do {
        if (p == NULL) break;
        if ( (id=get_id(p)) < 0) {
	    dbg("ERROR\n  Illegal identifier (%d-th element)"
		    " in format string! (%s)\n", (narg+1), p );
            va_end(parg);
            return -1;
        }
        dtype = vc[id].datatype;
        if ( dtype < POINTER ) {
            if ( dtype == INT ) 
		*((int*)va_arg(parg, int*)) = *((int*)vc[id].val);
	    else if ( dtype == CHAR )
		*((char*)va_arg(parg, char*)) = *((char*)vc[id].val);
	    else if ( dtype == DBLE )
		*((double*)va_arg(parg, double*)) = *((double*)vc[id].val);
	    else if ( dtype == LONG )
		*((long*)va_arg(parg, long*)) = *((long*)vc[id].val);
        } else {
            *((void**)va_arg(parg, void**)) = vc[id].val;
        }
        narg++;
    } while ( (p=strtok(NULL, sep)) != NULL );
    va_end(parg);
    return 0;
}

/** この関数群で保持している値を更新するための関数
 *
 * @param[in] format 更新したい変数名を空白区切りで列挙した文字列
 * @param[out] ... 更新する値が代入される変数
 *
 * @note \c format に書かれた文字列の順に、変数を並べる必要がある
 *
 * @ingroup ofmo-input
 * */
int ofmo_data_put_vals( const char* format, ... ) {
    int narg, len, id, dtype, i;
    char *p, *sep = " \t\n\r", strtmp[MAXSTRLEN];
    va_list parg;
    init_data();
    va_start(parg, format);
    len = strlen(format);
    if (len >= (MAXSTRLEN-1)) {
        dbg("Too long length of \'format\' string!\n");
        va_end(parg);
        return -1;
    }
    strncpy(strtmp, format, len);
    strtmp[len] = '\0';
    narg = 0;
    p = strtok(strtmp, sep);
    do {
        if (p == NULL) break;
        if ( ( id=get_id(p) ) < 0) {
	    dbg("ERROR\n  Illegal identifier (%d-th element)"
		    " in format string! (%s)\n", (narg+1), p );
            va_end(parg);
            return -1;
        }
        dtype = vc[id].datatype;
        if ( dtype < POINTER ) {
            if ( dtype == INT ) 
		*((int*)vc[id].val) = (int)va_arg(parg, int);
	    else if ( dtype == CHAR ) {
		i = (int)va_arg(parg, int);
		*((char*)vc[id].val) = (char)(i%256);
	    } else if ( dtype == DBLE )
		*((double*)vc[id].val) = (double)va_arg(parg, double);
	    else if ( dtype == LONG )
		*((long*)vc[id].val) = (long)va_arg(parg, long);
        } else {
	    vc[id].val = (void*)va_arg(parg, void*);
        }
        narg++;
    } while ( (p=strtok(NULL, sep)) != NULL );
    va_end(parg);
    return 0;
}

