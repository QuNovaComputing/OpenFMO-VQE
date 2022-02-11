/**
 * @file ofmo-input.c
 * 入力データを読み込むための関数群を定義したファイル。
 *
 * OpenFMOでは、GAMESSの入力形式に近い入力データを読み込む。
 * このファイルに記述された関数群は、各セクションの入力を担当する
 * 関数である。
 *
 * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <math.h>

#include "ofmo-def.h"
#include "ofmo-string.h"
#include "ofmo-data.h"

static char **local;
static char **elems;
static char line[MAXSTRLEN];

static int NATOM = 0;
static int NFRAG = 1;

/* ------------------------------------------------------
 * $CONTRL部分の読み込み
 * 以下の変数を読み込む
 *
 * maxit = Maximum number of repetitions of SCF
 * nprint = 出力レベル
 * itol   = primitiveのカットオフファクタ
 * icut   = 2電子積分のカットオフファクタ
 * method_list = ...
 * method = calculation method (RHF by default, VQE)
 * vqescr = vqe python script file name
 * ------------------------------------------------------- */
static int input_contrl( FILE *fp, char **token, const int ntok ) {
    //int maxit=30, nprint=0, itol=20, icut=15;
    int maxit=30, nprint=0, itol=20, icut=12;
	int method=OFMO_RHF;
	int * method_list = NULL;
	int * method_dim_idx = NULL;
	int * method_dim = NULL;
	int method_dim_n = 0;
	char *vqescr=(char*)malloc( sizeof(char) * MAXSTRLEN );
	char *desc=(char*)malloc( sizeof(char) * MAXSTRLEN );
	vqescr[0] = '\0';
	desc[0] = '\0';
    // local variables
    int i, nline = 0, n, ne, errpos = -1, flag_end = false;

    if ( fp == NULL ) {
	ofmo_data_put_vals("maxscf nprint itol icut method method_list method_dim_idx method_dim method_dim_n vqescr desc",
		                maxit, nprint,itol,icut,method,method_list,method_dim_idx,method_dim,method_dim_n,vqescr,desc);
	return 0;
    }

    ofmo_strcpy_2d( token, local, ntok, MAXSTRLEN ); // delimeters for token are " \t\n\r".
    n = ntok;
    while ( errpos < 0 && !flag_end ) {
	nline++;
	errpos = -1;
	for ( i=0; i<n; i++ ) {
	    if ( strcmp( local[i], "$contrl" ) == 0 ) continue;
	    if ( strcmp( local[i], "$end" ) == 0 ) {	// 終了
		flag_end = true;
		break;
	    }
	    ne = ofmo_split( local[i], "=", elems, MAXTOKEN, MAXSTRLEN );
	    if ( ne<2 ) {
		errpos = i;
		break;
	    }
	    if ( strcmp( elems[0], "maxit" ) == 0 )
		maxit = atoi( elems[1] );
	    else if ( strcmp( elems[0], "nprint" ) == 0 )
		nprint = atoi( elems[1] );
	    else if ( strcmp( elems[0], "itol" ) == 0 )
		itol = atoi( elems[1] );
	    else if ( strcmp( elems[0], "icut" ) == 0 )
		icut = atoi( elems[1] );
		else if ( strcmp( elems[0], "method" ) == 0){
			if (strcmp(elems[1], "rhf") == 0){
				method = OFMO_RHF;
			}
			else if (strcmp(elems[1], "vqe") == 0){
				method = OFMO_VQE;
			}
			else if (ofmo_find_char(elems[1], ',') > -1){
				int j;
				char** tmp_tok = (char **) malloc(NFRAG * sizeof(char *));
				for(j=0; j<NFRAG; j++) tmp_tok[j] = (char *)malloc(MAXSTRLEN);
				int num_tok = ofmo_split(elems[1], ",", tmp_tok, NFRAG, MAXSTRLEN);
				if(num_tok != NFRAG){
					for(j=0; j<NFRAG; j++) free(tmp_tok[j]);
					free(tmp_tok);
					dbg("line=%d, elem=%d : %s different from nfrag.\n", nline, ++i, elems[1] );
					break;
				}
				method = OFMO_UNDEF;
				method_list = (int *)malloc(NFRAG * sizeof(int));
				for(j=0; j<NFRAG; j++){
					if (strcmp(tmp_tok[j], "rhf") == 0) method_list[j] = OFMO_RHF;
					else if (strcmp(tmp_tok[j], "vqe") == 0) method_list[j] = OFMO_VQE;
					else{
						for(j=0; j<NFRAG; j++) free(tmp_tok[j]);
						free(tmp_tok);
						free(method_list);
						method_list = NULL;
						dbg("line=%d, elem=%d : %s different from nfrag.\n", nline, ++i, elems[1] );
						break;
					}
				}
			}
			else{
				dbg("line=%d, elem=%d : %s not valid\n", nline, ++i, elems[1] );
	    		break;
			}
		}
		else if ( strcmp( elems[0], "method_dim" ) == 0){
			int j;
			int end_dim = false;
			const int max_n_dim = NFRAG * (NFRAG - 1) / 2;
			char * tmp_elem1 = (char*) malloc( MAXSTRLEN );
			char ** contents = (char**) malloc( sizeof(char*) * max_n_dim);
			int dbg_flag = false;
			for(j=0; j<max_n_dim; j++) contents[j] = (char*) malloc( MAXSTRLEN );
			
			strcpy(tmp_elem1, elems[1]);
			int dim_count = ofmo_split(tmp_elem1, ";", contents, max_n_dim, MAXSTRLEN);
			method_dim_idx = (int *) malloc(sizeof(int) * 2 * dim_count);
			method_dim = (int *) malloc(sizeof(int) * dim_count);
			method_dim_n = dim_count;

			for(j=0; j<dim_count; j++){
				int dim1, dim2;
				char str_method[10];
				sscanf(contents[j], "(%d,%d):%s", &dim1, &dim2, str_method);
				method_dim_idx[2*j] = dim1-1;
				method_dim_idx[2*j+1] = dim2-1;
				if      (strcmp( str_method, "rhf" ) == 0) method_dim[j] = OFMO_RHF;
				else if (strcmp( str_method, "vqe" ) == 0) method_dim[j] = OFMO_VQE;
				else{
					dbg_flag = true;
					dbg("line=%d, elem=%d : %s not valid\n", nline, ++i, elems[1] );
					break;
				}
			}
			free(tmp_elem1);
			for(j=0; j<max_n_dim; j++) free(contents[j]);
			free(contents);
			if(dbg_flag){
				break;
			}
		}
		else if ( strcmp( elems[0], "vqescr" ) == 0)
			strcpy(vqescr, elems[1]);
		else if ( strcmp( elems[0], "desc") == 0 )
			strcpy(desc, elems[1]);
	}
	if ( errpos >= 0 ) {
	    dbg("line=%d, elem=%d : ERROR\n", nline, ++i );
	    break;
	}
	if ( flag_end ) break;
	// 次の行を読み込む
	n=ofmo_read_line(fp, MAXSTRLEN, MAXTOKEN, MAXSTRLEN, line, local);
	if ( n < 0 ) {
	    dbg("Unexpected EOF\n");
	    break;
	}
    }
    if ( !flag_end ) return -1;
	ofmo_data_put_vals("maxscf nprint itol icut method method_list method_dim_idx method_dim method_dim_n vqescr desc",
		                maxit, nprint,itol,icut,method,method_list,method_dim_idx,method_dim,method_dim_n,vqescr,desc);
   return 0;
}

/* ------------------------------------------------------
 * $SYSTEM部分の読み込み
 * 2011/03/02時点では、何もしない（$endまで流して読み込むだけ）
 * ------------------------------------------------------- */
static int input_system( FILE *fp, char **token, const int ntok ) {
    int i, nline = 0, n, errpos = -1, flag_end = false;
    if ( fp == NULL ) return 0;
    ofmo_strcpy_2d( token, local, ntok, MAXSTRLEN );
    n = ntok;
    while ( errpos < 0 && !flag_end ) {
	nline++;
	errpos = -1;
	for ( i=0; i<n; i++ ) {
	    if ( strcmp( local[i], "$system" ) == 0 ) continue;
	    if ( strcmp( local[i], "$end" ) == 0 ) {	// 終了
		flag_end = true;
		break;
	    }
	}
	if ( flag_end ) break;
	n=ofmo_read_line(fp, MAXSTRLEN, MAXTOKEN, MAXSTRLEN, line, local);
	// 次の行を読み込む
	if ( n < 0 ) {
	    dbg("Unexpected EOF\n");
	    break;
	}
    }
    if ( !flag_end ) return -1;

    return 0;
}

/* --------------------------------------------------------
 * $BASIS部分の読み込み
 * 以下の変数の値を読み込む
 * ngauss = (int)
 * diffsp = 重原子に対してSP型のdiffuse関数を加える(1)か否(0)か
 * diffs  = 水素原子に対してS型のdiffuse関数を加える(1)か否(0)か
 * gbasis = 関数のタイプ(STO, N21, N31, N311など）
 * ndfunc = 重原子に対するd型分極関数の個数
 * nffunc = 重原子に対するf型分極関数の個数
 * npfunc = 水素原子に対するp型分極関数の個数
 * -------------------------------------------------------- */
static int input_basis( FILE *fp, char **token, const int ntok ) {
    int ngauss=3, diffsp=false, diffs=false, npfunc=0, ndfunc=0, nffunc=0;
    char *gbasis;
    // local variables
    int i, nline = 0, n, ne, errpos = -1, flag_end = false;
    gbasis = (char*)malloc( sizeof(char) * MAXSTRLEN );
    gbasis[0] = '\0';
    if ( fp == NULL ) {
	ofmo_data_put_vals(
		"ngauss diffsp diffs ndfunc nffunc npfunc gbasis",
		ngauss, diffsp, diffs, ndfunc, nffunc, npfunc, gbasis );
	return 0;
    }
    ofmo_strcpy_2d( token, local, ntok, MAXSTRLEN );
    n = ntok;
    while ( errpos < 0 && !flag_end ) {
	nline++;
	errpos = -1;
	for ( i=0; i<n; i++ ) {
	    if ( strcmp( local[i], "$basis" ) == 0 ) continue;
	    if ( strcmp( local[i], "$end" ) == 0 ) {	// 終了
		flag_end = true;
		break;
	    }
	    ne = ofmo_split( local[i], "=", elems, MAXTOKEN, MAXSTRLEN );
	    if ( ne<2 ) {
		errpos = i;
		break;
	    }
	    if ( strcmp( elems[0], "ngauss" ) == 0 )
		ngauss = atoi( elems[1] );
	    else if ( strcmp( elems[0], "diffsp" ) == 0 )
		diffsp = (strcmp(elems[1], ".f.") == 0 ? false : true );
	    else if ( strcmp( elems[0], "diffs" ) == 0 )
		diffs = (strcmp(elems[1], ".f.") == 0 ? false : true );
	    else if ( strcmp( elems[0], "gbasis" ) == 0 )
		strcpy( gbasis, elems[1] );
	    else if ( strcmp( elems[0], "ndfunc" ) == 0 )
		ndfunc = atoi( elems[1] );
	    else if ( strcmp( elems[0], "nffunc" ) == 0 )
		nffunc = atoi( elems[1] );
	    else if ( strcmp( elems[0], "npfunc" ) == 0 )
		npfunc = atoi( elems[1] );
	}
	if ( errpos >= 0 ) {
	    dbg("line=%d, elem=%d : ERROR\n", nline, ++i );
	    break;
	}
	if ( flag_end ) break;
	// 次の行を読み込む
	n=ofmo_read_line(fp, MAXSTRLEN, MAXTOKEN, MAXSTRLEN, line, local);
	if ( n < 0 ) {
	    dbg("Unexpected EOF\n");
	    break;
	}
    }
    if ( !flag_end ) return -1;
    ofmo_data_put_vals("ngauss diffsp diffs ndfunc nffunc npfunc gbasis",
	    ngauss, diffsp, diffs, ndfunc, nffunc, npfunc, gbasis );

    return 0;
}

static int input_intgrl( FILE *fp, char **token, const int ntok ) {
    //long nintmx=15000, nintic=512*1024*1024;
    long nintmx=15000, nintic=512;
    // local variables
    int i, nline = 0, n, ne, errpos = -1, flag_end = false;
    if ( fp == NULL ) {
	ofmo_data_put_vals("nintmx nintic", nintmx, nintic );
	return 0;
    }
    ofmo_strcpy_2d( token, local, ntok, MAXSTRLEN );
    n = ntok;
    while ( errpos < 0 && !flag_end ) {
	nline++;
	errpos = -1;
	for ( i=0; i<n; i++ ) {
	    if ( strcmp( local[i], "$intgrl" ) == 0 ) continue;
	    if ( strcmp( local[i], "$end" ) == 0 ) {	// 終了
		flag_end = true;
		break;
	    }
	    ne = ofmo_split( local[i], "=", elems, MAXTOKEN, MAXSTRLEN );
	    if ( ne<2 ) {
		errpos = i;
		break;
	    }
	    if ( strcmp( elems[0], "nintmx" ) == 0 )
		nintmx = atol( elems[1] );
	    else if ( strcmp( elems[0], "nintic" ) == 0 )
		nintic = atol( elems[1] );
	}
	if ( errpos >= 0 ) {
	    dbg("line=%d, elem=%d : ERROR\n", nline, ++i );
	    break;
	}
	if ( flag_end ) break;
	// 次の行を読み込む
	n=ofmo_read_line(fp, MAXSTRLEN, MAXTOKEN, MAXSTRLEN, line, local);
	if ( n < 0 ) {
	    dbg("Unexpected EOF\n");
	    break;
	}
    }
    if ( !flag_end ) return -1;
    ofmo_data_put_vals("nintmx nintic", nintmx, nintic );

    return 0;
}

static int input_gddi( FILE *fp, char **token, const int ntok ) {
#ifdef USE_FALANX
    int ngroup=1, niogroup=2, nioprocs=1, master=0, iopos=-1;
#else
    int ngroup=1, niogroup=1, nioprocs=1, master=0, iopos=-1;
#endif
    // local variables
    int i, nline = 0, n, ne, errpos = -1, flag_end = false;
    if ( fp == NULL ) {
	ofmo_data_put_vals("ngroup niogroup nioprocs master iopos",
		ngroup, niogroup, nioprocs, master, iopos );
	return 0;
    }
    ofmo_strcpy_2d( token, local, ntok, MAXSTRLEN );
    n = ntok;
    while ( errpos < 0 && !flag_end ) {
	nline++;
	errpos = -1;
	for ( i=0; i<n; i++ ) {
	    if ( strcmp( local[i], "$gddi" ) == 0 ) continue;
	    if ( strcmp( local[i], "$end" ) == 0 ) {	// 終了
		flag_end = true;
		break;
	    }
	    if ( strncmp( local[i], "ngroup", 6) == 0 ) {
		ne=ofmo_split(local[i], "=", elems, MAXTOKEN, MAXSTRLEN);
		if ( ne<2 ) {
		    errpos = i;
		    break;
		}
		ngroup = atoi( elems[1] );
		if ( ngroup < 1 ) ngroup = 1;
	    } else if ( strncmp( local[i], "niogroup", 8) == 0 ) {
		ne=ofmo_split(local[i], "=", elems, MAXTOKEN, MAXSTRLEN);
		if ( ne<2 ) {
		    errpos = i;
		    break;
		}
		niogroup = atoi( elems[1] );
		if ( niogroup < 1 ) niogroup = 1;
	    } else if ( strncmp( local[i], "nioprocs", 8) == 0 ) {
		ne=ofmo_split(local[i], "=", elems, MAXTOKEN, MAXSTRLEN);
		if ( ne<2 ) {
		    errpos = i;
		    break;
		}
		nioprocs = atoi( elems[1] );
		if ( nioprocs < 1 ) nioprocs = 1;
	    } else if ( strncmp( local[i], "master", 6) == 0 ) {
		ne=ofmo_split(local[i], "=", elems, MAXTOKEN, MAXSTRLEN);
		if ( ne<2 ) {
		    errpos = i;
		    break;
		}
		master = atoi( elems[1] );
		if ( master < 0 ) master = 0;
	    } else if ( strncmp( local[i], "iopos", 5) == 0 ) {
		ne=ofmo_split(local[i], "=", elems, MAXTOKEN, MAXSTRLEN);
		if ( ne<2 ) {
		    errpos = i;
		    break;
		}
		iopos = atoi( elems[1] );
	    }
	}
	if ( errpos >= 0 ) {
	    dbg("line=%d, elem=%d : ERROR\n", nline, ++i );
	    break;
	}
	if ( flag_end ) break;
	// 次の行を読み込む
	n=ofmo_read_line(fp, MAXSTRLEN, MAXTOKEN, MAXSTRLEN, line, local);
	if ( n < 0 ) {
	    dbg("Unexpected EOF\n");
	    break;
	}
    }
    if ( !flag_end ) return -1;
    ofmo_data_put_vals("ngroup niogroup nioprocs master iopos",
	    ngroup, niogroup, nioprocs, master, iopos );

    return 0;
}

static int input_scf( FILE *fp, char **token, const int ntok ) {
    int diis=true, npunch=0, conv=7;
    // local variables
    int i, nline = 0, n, ne, errpos = -1, flag_end = false;
    if ( fp == NULL ) {
	ofmo_data_put_vals("diis npunch scfconv", diis, npunch, conv );
	return 0;
    }
    ofmo_strcpy_2d( token, local, ntok, MAXSTRLEN );
    n = ntok;
    while ( errpos < 0 && !flag_end ) {
	nline++;
	errpos = -1;
	for ( i=0; i<n; i++ ) {
	    if ( strcmp( local[i], "$scf" ) == 0 ) continue;
	    if ( strcmp( local[i], "$end" ) == 0 ) {	// 終了
		flag_end = true;
		break;
	    }
	    ne = ofmo_split( local[i], "=", elems, MAXTOKEN, MAXSTRLEN );
	    if ( ne<2 ) {
		errpos = i;
		break;
	    }
	    if ( strcmp( elems[0], "diis" ) == 0 )
		diis = ( strcmp(elems[1], ".f.")==0 ? false : true );
	    else if ( strcmp( elems[0], "npunch" ) == 0 )
		npunch = atoi( elems[1] );
	    else if ( strcmp( elems[0], "conv") == 0 ) {
		//conv = atoi( elems[1] );
		//if ( conv < 1 ) conv = 5;
                char *endptr=NULL;
                double convd;
		convd = strtod( elems[1], &endptr );
                if (convd <=0 ) {
                  dbg("CONV in $SCF must be a positive value.\n");
                  break;
                } else if ( convd < 1 ) conv = -log10(convd);
                else conv = convd;
	    }
	}
	if ( errpos >= 0 ) {
	    dbg("line=%d, elem=%d : ERROR\n", nline, ++i );
	    break;
	}
	if ( flag_end ) break;
	// 次の行を読み込む
	n=ofmo_read_line(fp, MAXSTRLEN, MAXTOKEN, MAXSTRLEN, line, local);
	if ( n < 0 ) {
	    dbg("Unexpected EOF\n");
	    break;
	}
    }
    if ( !flag_end ) return -1;
    ofmo_data_put_vals("diis npunch scfconv", diis, npunch, conv );

    return 0;
}

/* ------------------------------------------------------
 * $FMOPRP部分の読み込み
 * 2012/01/04時点では、SCCCONVとMAXSCCを読み込む
 * ------------------------------------------------------- */
static int input_fmoprp( FILE *fp, char **token, const int ntok ) {
    int i, nline = 0, n, ne, errpos = -1, flag_end = false;
    int maxscc=30, sccconv=7;
    if ( fp == NULL ) {
	ofmo_data_put_vals("sccconv maxscc", sccconv, maxscc );
	return 0;
    }
    ofmo_strcpy_2d( token, local, ntok, MAXSTRLEN );
    n = ntok;
    while ( errpos < 0 && !flag_end ) {
	nline++;
	errpos = -1;
	for ( i=0; i<n; i++ ) {
	    if ( strcmp( local[i], "$fmoprp" ) == 0 ) continue;
	    if ( strcmp( local[i], "$end" ) == 0 ) {	// 終了
		flag_end = true;
		break;
	    }
	    ne = ofmo_split( local[i], "=", elems, MAXTOKEN, MAXSTRLEN );
#if 0
	    if ( ne<2 ) {
		errpos = i;
		break;
	    }
#endif
	    if ( strcmp( elems[0], "maxit" ) == 0 ) {
		maxscc = atoi( elems[1] );
		if ( maxscc < 1 ) maxscc = 30;
	    } else if ( strcmp( elems[0], "conv" ) == 0 ) {
		//sccconv = atoi( elems[1] );
		//if ( sccconv < 1 ) sccconv = 5;
                char *endptr=NULL;
                double sccconvd;
		sccconvd = strtod( elems[1], &endptr );
                if (sccconvd <=0 ) {
                  dbg("CONV in $FMOPROP must be a positive value.\n");
                  break;
                } else if ( sccconvd < 1 ) sccconv = -log10(sccconvd);
                else sccconv = sccconvd;
	    }
	}
	if ( errpos >= 0 ) {
	    dbg("line=%d, elem=%d : ERROR\n", nline, ++i );
	    break;
	}
	if ( flag_end ) break;
	// 次の行を読み込む
	n=ofmo_read_line(fp, MAXSTRLEN, MAXTOKEN, MAXSTRLEN, line, local);
	if ( n < 0 ) {
	    dbg("Unexpected EOF\n");
	    break;
	}
    }
    ofmo_data_put_vals("sccconv maxscc", sccconv, maxscc );
    if ( !flag_end ) exit(1);
    if ( !flag_end ) return -1;

    return 0;
}

/*
   原子座標の入力
   （絶対に呼ばれる必要がある）
   */
static int input_fmoxyz( FILE *fp, char **token, const int ntok ) {
    double *atom_x, *atom_y, *atom_z;
    int *atomic_number;
    // local variables
    int nline = 0, n, flag_end = false;
    int nat = 0, pos=-1;

    // memory allocation
    if ( NATOM < 1 ) {
	dbg("# of atom hasn\'t be read, yet\n");
	return -1;
    }
    atom_x = (double*)malloc( sizeof(double) * NATOM );
    atom_y = (double*)malloc( sizeof(double) * NATOM );
    atom_z = (double*)malloc( sizeof(double) * NATOM );
    atomic_number = (int*)malloc( sizeof(int) * NATOM );
    /*if ( atom_x == NULL || atom_y == NULL ||
	    atom_z == NULL || atomic_number == NULL ) {
	//errpos = 101;
    }*/
    ofmo_strcpy_2d( token, local, ntok, MAXSTRLEN );
    n = ntok;
    while (1) {
	nline++;
	if ( strcmp( local[0], "$fmoxyz" ) == 0 ) {
	    // 次の行を読み込む
	    n=ofmo_read_line(fp, MAXSTRLEN, MAXTOKEN, MAXSTRLEN,
		    line, local);
	    if ( n < 0 ) {
		dbg("Unexpected EOF\n");
		break;
	    }
	    continue;
	}
	if ( strcmp( local[0], "$end" ) == 0 ) {
	    flag_end = true;
	    break;
	}
	if ( n < 5 ) {	// 要素が足りないとき
	    dbg("Illegal number of elements (line=%d)\n", nline );
	    break;
	}
	if ( nat == 0 ) {
	    if      ( isalpha(local[0][0]) ) pos = 0;
	    else if ( isalpha(local[1][0]) ) pos = 1;
	}

	// 読み込み
	//atomic_number[nat] = atoi( local[1] );
	atomic_number[nat] = AtomicNumber( local[pos] );
	atom_x[nat] = atof( local[2] );
	atom_y[nat] = atof( local[3] );
	atom_z[nat] = atof( local[4] );
	nat++;
	// 次の行を読み込む
	n=ofmo_read_line(fp, MAXSTRLEN, MAXTOKEN, MAXSTRLEN, line, local);
	if ( n < 0 ) {
	    dbg("Unexpected EOF\n");
	    break;
	}
    }

    if ( !flag_end ) return -1;
    if ( nat != NATOM ) {
	dbg("Inconsisten between nat(%d) vs NATOM(%d)\n", nat, NATOM );
	return -1;
    }

    /* change unit of length */
    double rb;
    rb = 1.e0 / BOHR_RADIUS;
    for ( int iat=0; iat<nat; iat++ ) {
	atom_x[iat] *= rb;
	atom_y[iat] *= rb;
	atom_z[iat] *= rb;
    }

    ofmo_data_put_vals("natom atn atx aty atz",
	    nat, atomic_number, atom_x, atom_y, atom_z );


    return 0;
}

/* 各モノマーの電荷を入力する
   （必ず呼ばれる必要がある）
   */
static int input_fmo_icharg( FILE *fp, char **token, const int ntok ) {
    int *icharg;
    int flag = false;	// icharg=true 読み込み中, =false 読み込んでない
    int valid = false;	// ichargの行に、要素があるか？
    // local variables
    int i, nline = 0, n, nr, ne, errpos = -1, flag_end = false;
    if ( NFRAG < 1 ) {
	dbg("# of fragment hasn\'t be read, yet\n");
	return -1;
    }
    icharg = (int*)malloc( sizeof(int) * NFRAG );
    ofmo_strcpy_2d( token, local, ntok, MAXSTRLEN );
    n = ntok;
    while ( errpos < 0 && !flag_end ) {
	nline++;
	if ( flag == false ) {
	    for ( i=0; i<n; i++ ) {
		if (strncmp( local[i], "icharg", 6) == 0 ) {// 発見！
		    strcpy( line, local[i] );// この行の残りをコピー
		    for ( int j=i+1; j<n; j++ ) strcat( line, local[j] );
		    ofmo_delete_all_space( line );
		    ne=ofmo_split(line, "=", elems, MAXTOKEN, MAXSTRLEN);
		    if ( ne > 1 ) {	// '='以降をコピー
			strcpy( line, elems[1] );
			valid = true;
		    }
		    nr = 0;
		    flag = true;
		    break;
		}
	    }
	} else {
	    //ne = ofmo_split( line, ",", elems, MAXTOKEN, MAXSTRLEN );
	    ofmo_replace_char( line, ',', ' ' );
	    ofmo_delete_prepost_space( line );
	    ne = ofmo_split( line, " ", elems, MAXTOKEN, MAXSTRLEN );
	    for ( i=0; i<ne; i++ ) icharg[nr++] = atoi( elems[i] );
	    if ( nr == NFRAG ) flag_end = true;
	}
	if ( flag_end ) break;
	// （最初の行以外で）次の行を読み込む
	if ( valid == false ) {
	    if ( fgets( line, MAXSTRLEN, fp ) == NULL ) {
		dbg("Unexpected EOF\n");
		break;
	    }
	    ofmo_tolower( line );
	    ofmo_delete_all_space( line );
	}
	valid = false;
    }
    if ( !flag_end ) {
	dbg(" %d elements has been read, but not reach nfrag(%d)\n",
		nr, NFRAG );
	return -1;
    }
    ofmo_data_put_vals("icharg", icharg );
    return 0;
}

/* 各原子がどのモノマーに含まれているか？のデータを入力する
   （必ず呼ばれる必要がある）
   */
static int input_fmo_indat( FILE *fp, char **token, const int ntok ) {
    int *atom2frg;
    int flag = false;	// indat=true 読み込み中, =false 読み込んでない
    int valid = false;	// indatの行に、要素があるか？
    // local variables
    int i, nline = 0, n, nr, ne, errpos = -1, flag_end = false;
    int mode = -1; /* 0 = 通常、1 = 短縮形 */
    int nf, is, prev;
    if ( NATOM < 1 ) {
	dbg("# of atom hasn\'t be read, yet\n");
	return -1;
    }
    atom2frg = (int*)malloc( sizeof(int) * NATOM );
    ofmo_strcpy_2d( token, local, ntok, MAXSTRLEN );
    n = ntok;
    while ( errpos < 0 && !flag_end ) {
	nline++;
	if ( flag == false ) {
	    for ( i=0; i<n; i++ ) {
		if (strncmp( local[i], "indat", 5) == 0 ) {// 発見！
		    strcpy( line, local[i] );// この行の残りをコピー
		    for ( int j=i+1; j<n; j++ ) strcat( line, local[j] );
		    //ofmo_delete_all_space( line );
		    ne=ofmo_split( line, "=", elems, MAXTOKEN, MAXSTRLEN);
		    if ( ne > 1 ) {	// '='以降をコピー
			strcpy( line, elems[1] );
			valid = true;
		    }
		    nr = 0;
		    flag = true;
		    break;
		}
	    }
	} else {
	    int k, tmp;
	    //ne = ofmo_split( line, ",", elems, MAXTOKEN, MAXSTRLEN );
	    ofmo_replace_char( line, ',', ' ' );
            ofmo_delete_multiple_space( line );
	    ofmo_delete_prepost_space( line );
	    ne = ofmo_split( line, " ", elems, MAXTOKEN, MAXSTRLEN );
	    is = 0;
	    if ( mode == -1 ) {
		mode = ( line[0] == '0' ? 1 : 0 );
		if ( mode == 1 ) {
		    nf = 0;
		    prev = 0;
		    is = 1;
		} else is = 0;
	    }
	    if ( mode == 0 ) {
		for ( i=is; i<ne; i++ ) atom2frg[nr++] = atoi(elems[i])-1;
	    } else if ( mode == 1 ) {
		for ( i=is; i<ne; i++ ) {
		    tmp = atoi( elems[i] );
		    //printf("elems[%d]=%s\n", i, elems[i] );
		    //fflush(stdout);

		    if ( tmp > 0 ) {
			if ( prev > 0 ) {	// ++ = 前の要素を追加
			    // debug
			    /*printf("prev=%d, tmp=%d\n", prev, tmp );
			    fflush(stdout);*/
			    atom2frg[prev-1] = nf;
			    nr++;
			}
		    } else if ( tmp < 0 ) {
			if ( prev > 0 ) {	// -+ = まとめて追加
			    for ( k=prev; k<=(-tmp); k++ ) {
				atom2frg[k-1] = nf;
				nr++;
			    }
			} else {	// 0-, -- = エラー
			    errpos = i;
			    break;
			}
		    } else {	// 今回が０
			if ( prev == 0 ) {	// 00 = エラー
			    errpos = i;
			    break;
			}
			if ( prev > 0 ) {	// +0
			    atom2frg[prev-1] = nf;
			    nr++;
			}
			nf++;
		    }

		    prev = tmp;
		}	// for ( i=is )
	    }	// if ( mode == )
	    if ( nr == NATOM ) flag_end = true;
	}
	if ( flag_end ) break;
	// （最初の行以外で）次の行を読み込む
	if ( valid == false ) {
	    if ( fgets( line, MAXSTRLEN, fp ) == NULL ) {
		dbg("Unexpected EOF\n");
		break;
	    }
	    ofmo_tolower( line );
	    //ofmo_delete_all_space( line );
	}
	valid = false;
    }
    if ( !flag_end ) {
	dbg(" %d elements has been read, but not reach nfrag(%d)\n",
		nr, NATOM );
	return -1;
    }
    ofmo_data_put_vals("at2frg", atom2frg );
    return 0;
}

/* ----------------------------------------------------
 * FMO関連の入力（必ず呼ばれる必要がある）
 * 以下の変数を読む
 * nfrag = フラグメント数
 * nbody = N-body (2011/0721現在は、呼んでいない）
 * icharg[] = 各フラグメントの電荷
 * indat[]  = 各原子が属するフラグメントの番号
 * ---------------------------------------------------- */
static int input_fmo( FILE *fp, char **token, const int ntok ) {
    int nfrag=1, nbody=2;
    double lptc=2.e0, laop=1.e0, ldim=2.e0;
    int flag = false;	/* nfragを呼んだ=true, 読んでいない =false */
    // local variables
    long pos;
    int i, nline = 0, n, ne, errpos = -1, flag_end = false;
    int flag_ichrg = false, flag_indat = false;
    pos = ftell( fp );
    ofmo_strcpy_2d( token, local, ntok, MAXSTRLEN );
    n = ntok;
    while ( errpos < 0 && !flag_end ) {
	nline++;
	errpos = -1;
	for ( i=0; i<n; i++ ) {
	    if ( strcmp( local[i], "$fmo" ) == 0 ) continue;
	    if ( strcmp( local[i], "$end" ) == 0 ) {	// 終了
		flag_end = true;
		break;
	    }
	    if ( strncmp( local[i], "nfrag", 5)==0 && flag == false ) {
		ne=ofmo_split(local[i], "=", elems, MAXTOKEN, MAXSTRLEN);
		if ( ne < 2 ) {
		    errpos = i;
		    break;
		}
		NFRAG = nfrag = atoi( elems[1] );
		flag = true;
		fseek( fp, pos, SEEK_SET );
	    } else if ( strncmp( local[i], "respap", 6)==0 ) {
		ne=ofmo_split(local[i], "=", elems, MAXTOKEN, MAXSTRLEN);
		if ( ne < 2 ) {
		    errpos = i;
		    break;
		}
		laop = atof( elems[1] );
		if ( laop < 0.e0 ) laop = 1.e0;
	    } else if ( strncmp( local[i], "resppc", 6)==0 ) {
		ne=ofmo_split(local[i], "=", elems, MAXTOKEN, MAXSTRLEN);
		if ( ne < 2 ) {
		    errpos = i;
		    break;
		}
		lptc = atof( elems[1] );
		if ( lptc < 0.e0 ) lptc = 2.e0;
	    } else if ( strncmp( local[i], "resdim", 6)==0 ) {
		ne=ofmo_split(local[i], "=", elems, MAXTOKEN, MAXSTRLEN);
		if ( ne < 2 ) {
		    errpos = i;
		    break;
		}
		ldim = atof( elems[1] );
		if ( ldim < 0.e0 ) ldim = 2.e0;
	    } else if ( strncmp( local[i], "nbody", 5)==0 ) {
		ne=ofmo_split(local[i], "=", elems, MAXTOKEN, MAXSTRLEN);
		if ( ne < 2 ) {
		    errpos = i;
		    break;
		}
		nbody = atoi( elems[1] );
		if ( nbody < 0 ) nbody = 2;
	    } else if ( strncmp( local[i], "icharg", 6)==0 ) {
		if ( flag == true && flag_ichrg == false ) {
		    input_fmo_icharg( fp, local, n );
		    flag_ichrg = true;
		}
	    } else if ( strncmp( local[i], "indat", 5)==0 ) {
		if ( flag == true && flag_indat == false ) {
		    input_fmo_indat( fp, local, n );
		    flag_indat = true;
		}
	    }
	}
	if ( flag_end == true ) break;
	if ( errpos >= 0 ) dbg("line=%d, elem=%d : ERROR\n", nline, ++i );
	// 次の行を読み込む
	n=ofmo_read_line(fp, MAXSTRLEN, MAXTOKEN, MAXSTRLEN, line, local);
	if ( n < 0 ) {
	    dbg("Unexpected EOF\n");
	    break;
	}
    }
    if ( !flag_end ) return -1;
    ofmo_data_put_vals("nfrag nbody lptc laop ldim",
	    nfrag, nbody, lptc, laop, ldim );
    return 0;
}

/*
 * LMO関連のパラメタを読み取る（必ず呼ばれる必要がある）
 * 以下の制限がある
 * 1. １種類のLMOパラメタしか読まない
 * 2. 各行の2番目のパラメタは無視する
 * */
static int input_fmolmo( FILE *fp, char **token, const int ntok ) {
    int nlmo=0, nfunc=0, *lea = NULL;
    double *clmo = NULL;
    // local variables
    int i, n, flag_end = false, reading, iao=-1;
    int lmo = -1;
    char *bs_name;
    int nbs=0;
    bs_name = (char*)malloc( sizeof(char) * MAXSTRLEN );
    bs_name[0] = '\0';
    ofmo_strcpy_2d( token, local, ntok, MAXSTRLEN );
    n = ntok;
    while (1) {
	n=ofmo_read_line(fp, MAXSTRLEN, MAXTOKEN, MAXSTRLEN, line, local);
	if ( n < 0 ) {
	    dbg("Unexpected EOF\n");
	    break;
	}
	if ( strcmp( local[0], "$end" ) == 0 ) {
	    if ( iao != nfunc ) {
		dbg("error: Illegal number of AO(%d)\n", iao );
	    }
	    flag_end = true;
	    break;
	}
        if (nbs>0) continue;
      /*
	if ( strcmp( local[0], "$fmolmo") == 0 ||
		strcmp( local[0], "$fmohyb" ) == 0 ) {
	} else if ( strcmp( local[0], "$end" ) == 0 ) {
	    if ( iao != nfunc ) {
		dbg("error: Illegal number of AO(%d)\n", iao );
	    }
	    flag_end = true;
	    break;
	} else if ( n == 3 ) {
            */
	if ( n == 3 ) {
	    strcpy( bs_name, local[0] );
	    ofmo_toupper( bs_name );
	    nfunc = atoi( local[1] );
	    nlmo  = atoi( local[2] );
	    clmo = (double*)malloc( sizeof(double)*nfunc*nlmo );
	    lea  = (int*)malloc( sizeof(int)*nlmo );
	    lmo = -1;
	    reading=false;
	} else {
	    if ( reading == false ) {
		if ( n > 5 ) {// first line of each LMO
		    iao=0;
		    lmo++;
		    reading = true;
		    lea[lmo] = atoi( local[0] );
		    for ( i=2; i<n; i++, iao++ )
			clmo[lmo*nfunc+iao] = atof( local[i] );
		    if ( iao == nfunc ) reading = false;
		}
	    } else {
		for ( i=0; i<n; i++, iao++ )
		    clmo[lmo*nfunc+iao] = atof( local[i] );
		if ( iao == nfunc ) reading = false;
	    }
	    if ( lmo == nlmo-1 && reading == false ) nbs++;
	}
        /*
	n=ofmo_read_line(fp, MAXSTRLEN, MAXTOKEN, MAXSTRLEN, line, local);
	if ( n < 0 ) {
	    dbg("Unexpected EOF\n");
	    break;
	}
        */
    }
    if ( !flag_end ) return -1;
    ofmo_data_put_vals("nlmo nfunc clmo ban lmobs",
	    nlmo, nfunc, clmo, lea, bs_name );
    return 0;
}

/* 結合原子の情報を入力する（必ず呼ばれる必要がある）
   */
static int input_fmobnd( FILE *fp, char **token, const int ntok ) {
    int nbond=0, *welec=NULL, *woelec=NULL;
    // local variables
    int n, flag_end = false, ibond=-1;
    int counted_nbond = false;
    long pos=0;

    ofmo_strcpy_2d( token, local, ntok, MAXSTRLEN );
    n = ntok;
    while (1) {
	if ( strcmp(local[0], "$fmobnd") == 0 ) {
	    ibond = 0;
	    pos = ftell( fp );
	} else if ( strcmp( local[0], "$end" ) == 0 ) {
	    if ( counted_nbond == false ) {
		counted_nbond = true;
		nbond = ibond;
		welec  = (int*)malloc( sizeof(int) * nbond );
		woelec = (int*)malloc( sizeof(int) * nbond );
		ibond = 0;
		fseek( fp, pos, SEEK_SET );
	    } else {
		flag_end = true;
		break;
	    }
	} else if ( counted_nbond ) {
	    woelec[ibond] = (-atoi( local[0])) - 1  ;/* 0-offset */
	    welec[ibond]  =  atoi( local[1] ) - 1 ;/* 0-offset */
	    ibond++;
	} else ibond++;
	n=ofmo_read_line(fp, MAXSTRLEN, MAXTOKEN, MAXSTRLEN, line, local);
	if ( n < 0 ) {
	    dbg("Unexpected EOF\n");
	    break;
	}
    }

    if ( !flag_end ) return -1;
    ofmo_data_put_vals("nbond welec woelec", nbond, welec, woelec );
    return 0;
}

/* ------------------------------------------------------
 * $DATA部分の読み込み
 * 2011/03/02時点では、何もしない（$endまで流して読み込むだけ）
 * ------------------------------------------------------- */
static int input_data( FILE *fp, char **token, const int ntok ) {
    int i, nline = 0, n, flag_end = false;
    char **now;
    if ( fp == NULL ) return 0;
    while (1) {
	if ( nline == 0 ) {
	    now = token;
	    n   = ntok;
	    nline++;
	} else {
	    if ( fgets( line, MAXSTRLEN, fp ) == NULL ) {
		dbg("Unexpected EOF\n");
		return -1;
	    }
	    nline++;
	    ofmo_delete_prepost_space( line );
	    n = ofmo_split( line, " \t\n\r", local, MAXTOKEN, MAXSTRLEN );
	    if ( n<1 ) continue;
	    now = local;
	}
	for ( i=0; i<n; i++ ) {
	    if ( strcasecmp( now[i], "$data" ) == 0 ) continue;
	    if ( strcasecmp( now[i], "$end" ) == 0 ) {
		flag_end = true;
		break;
	    }
	}
	if ( flag_end ) break;
    }
    if ( !flag_end ) return -1;
    return 0;
}

static int input_nat( FILE* fp ) {
    int is_reading_nat = false, nat = 0;
    // $fmoxyzのサーチ
    while ( fgets( line, MAXSTRLEN, fp ) != NULL ) {
	ofmo_delete_prepost_space( line );
	if ( strncasecmp( line, "$fmoxyz", 7 ) == 0 ) {
	    is_reading_nat = true;
	} else if ( is_reading_nat == true ) {
	    if ( strncasecmp( line, "$end", 4 ) == 0 ) {
		return nat;
	    } else nat++;
	}
    }
    return nat;
}

static int input_vqeprp( FILE *fp, char **token, const int ntok ){
	int monhomo = -1, monlumo = -1, dimhomo = -1, dimlumo = -1;
	int monent = -1, diment = -1;

	int i, nline = 0, n, ne, errpos = -1, flag_end = false;
	if ( fp == NULL ) {
	ofmo_data_put_vals("monhomo monlumo dimhomo dimlumo monent diment",
		monhomo, monlumo, dimhomo, dimlumo, monent, diment);
	return 0;
    }

	ofmo_strcpy_2d( token, local, ntok, MAXSTRLEN );
	n = ntok;
	while( errpos < 0 && !flag_end ){
	nline++;
	errpos = -1;
	for ( i=0; i<n; i++ ){
		if ( strcmp( local[i], "$vqeprp" ) == 0 ) continue;
		if ( strcmp( local[i], "$end" ) == 0 ) {
			flag_end = true;
			break;
		}
		ne = ofmo_split( local[i], "=", elems, MAXTOKEN, MAXSTRLEN );
		if ( ne<2 ){
			errpos = i;
			break;
		}
		if ( strcmp( elems[0], "monhomo" ) == 0)
			monhomo = atoi( elems[1] );
		if ( strcmp( elems[0], "monlumo" ) == 0)
			monlumo = atoi( elems[1] );
		if ( strcmp( elems[0], "dimhomo" ) == 0)
			dimhomo = atoi( elems[1] );
		if ( strcmp( elems[0], "dimlumo" ) == 0)
			dimlumo = atoi( elems[1] );
		if ( strcmp( elems[0], "monent" ) == 0)
			monent = atoi( elems[1] );
		if ( strcmp( elems[0], "diment" ) == 0)
			diment = atoi( elems[1] );
	}
	if ( errpos >= 0 ){
		dbg("line=%d, elem=%d : ERROR\n", nline, ++i );
	    break;
	}
	if ( flag_end ) break;
	n=ofmo_read_line(fp, MAXSTRLEN, MAXTOKEN, MAXSTRLEN, line, local);
	if ( n < 0 ) {
	    dbg("Unexpected EOF\n");
	    break;
	}
	}
	if ( !flag_end ) return -1;
	ofmo_data_put_vals("monhomo monlumo dimhomo dimlumo monent diment",
		monhomo, monlumo, dimhomo, dimlumo, monent, diment);
    return 0;

}

/** Function to read input data in GAMESS format
 *
 * Drives a function that reads the data for each section
 *
 * @param[in] fp File descriptor of input file
 *
 * @note The input file must be open when this function is called
 *
 * @ingroup ofmo-input
 * */
int ofmo_read_gamess_input( const char *filename ) {
    int ntok, i, ierr=0, nat;
    char **token;
    int ctrl, sys, basis, intgrl, gddi, scf, fmo, fmoprp, fmoxyz;
    int fmolmo, fmobnd, data;
    FILE *fp;
    if ( (fp=fopen( filename, "r")) == NULL ) {
	dbg("Failure in opening file (%s)\n", filename );
	return -1;
    }

    token = ofmo_alloc_char_2d( MAXTOKEN, MAXSTRLEN );
    local = ofmo_alloc_char_2d( MAXTOKEN, MAXSTRLEN );
    elems = ofmo_alloc_char_2d( MAXTOKEN, MAXSTRLEN );

    if ( token == NULL || local == NULL || elems == NULL ) {
	dbg("Failure in allocation temporary space\n");
	return -1;
    }

    if ( (nat=input_nat( fp )) < 1 ) {
	dbg("Illegal number of atom (%d)\n", nat );
	return -1;
    }
    NATOM = nat;

    rewind( fp );
    ctrl = sys = basis = intgrl = gddi = scf = fmo = fmoprp
	= fmoxyz = fmolmo = fmobnd = data = false;
    while ( fgets( line, MAXSTRLEN, fp ) != NULL ) {
	ofmo_delete_prepost_space( line );
	if ( line[0] == '!' ) continue;
	ntok = ofmo_split( line, " \t\n\r", token, MAXTOKEN, MAXSTRLEN );
	if ( ntok < 1 ) continue;
	for ( i=0; i<ntok; i++ ) ofmo_tolower( token[i] );
	if ( token[0][0] == '$' ) {
	    if      ( strcmp( &token[0][1], "contrl" ) == 0 ) {
		ierr=input_contrl( fp, token, ntok );
		//printf("== 1 ==\n");
		ctrl = true;
	    } else if ( strcmp( &token[0][1], "system" ) == 0 ) {
		ierr=input_system( fp, token, ntok );
		//printf("== 2 ==\n");
		sys = true;
	    } else if ( strcmp( &token[0][1], "basis" ) == 0 ) {
		ierr=input_basis( fp, token, ntok );
		//printf("== 3 ==\n");
		basis = true;
	    } else if ( strcmp( &token[0][1], "intgrl" ) == 0 ) {
		ierr=input_intgrl( fp, token, ntok );
		//printf("== 4 ==\n");
		intgrl = true;
	    } else if ( strcmp( &token[0][1], "gddi" ) == 0 ) {
		ierr=input_gddi( fp, token, ntok );
		//printf("== 5 ==\n");
		gddi = true;
	    } else if ( strcmp( &token[0][1], "scf" ) == 0 ) {
		ierr=input_scf( fp, token, ntok );
		//printf("== 6 ==\n");
		scf = true;
	    } else if ( strcmp( &token[0][1], "fmo" ) == 0 ) {
		ierr=input_fmo( fp, token, ntok );
		//printf("== 7 ==\n");
		fmo = true;
	    } else if ( strcmp( &token[0][1], "fmoprp" ) == 0 ) {
		ierr=input_fmoprp( fp, token, ntok );
		//printf("== 8 ==\n");
		fmoprp = true;
	    } else if ( strcmp( &token[0][1], "fmoxyz" ) == 0 ) {
		ierr=input_fmoxyz( fp, token, ntok );
		//printf("== 9 ==\n");
		fmoxyz = true;
	    } else if ( strcmp( &token[0][1], "fmolmo" ) == 0 ||
		    strcmp( &token[0][1], "fmohyb" ) == 0 ) {
		ierr=input_fmolmo( fp, token, ntok );
		//printf("== 10 ==\n");
		fmolmo = true;
	    } else if ( strcmp( &token[0][1], "fmobnd" ) == 0 ) {
		ierr=input_fmobnd( fp, token, ntok );
		//printf("== 11 ==\n");
		fmobnd = true;
	    } else if ( strcmp( &token[0][1], "data" ) == 0 ) {
		ierr=input_data( fp, token, ntok );
		//printf("== 12 ==\n");
		data = true;
	    } else if ( strcmp( &token[0][1], "vqeprp" ) == 0 ) {
		ierr=input_vqeprp( fp, token, ntok);
		data = true;
		} else {
		printf("Illegal section name (%s)\n", &token[0][1] );
		ierr = -2;
	    }
	}
    }
    int ierr0 = 0;
    if ( !ctrl ) input_contrl( NULL, token, ntok );
    if ( !sys  ) input_system( NULL, token, ntok );
    if ( !basis ) input_basis( NULL, token, ntok );
    if ( !intgrl ) input_intgrl( NULL, token, ntok );
    if ( !gddi ) input_gddi( NULL, token, ntok );
    if ( !scf ) input_scf( NULL, token, ntok );
    if ( !fmoprp ) input_fmoprp( NULL, token, ntok );
    if ( !data ) input_data( NULL, token, ntok );
    if ( !fmo ) { ierr0++; dbg("error: no $fmo section\n"); }
    if ( !fmoxyz ) { ierr0++; dbg("error: no $fmoxyz section\n"); }
    if ( !fmolmo ) { ierr0++; dbg("error: no $fmolmo section\n"); }
    if ( !fmobnd ) { ierr0++; dbg("error: no $fmobnd section\n"); }

    ofmo_free_char_2d( token );
    ofmo_free_char_2d( local );
    ofmo_free_char_2d( elems );
    if ( ierr0 > 0 ) ierr = -1;
    return ierr;
}
