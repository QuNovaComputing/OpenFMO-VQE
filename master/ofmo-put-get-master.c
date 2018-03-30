#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "ofmo-def.h"
#include "ofmo-data-struct.h"
#include "ofmo-data.h"

#include "ofmo-misc.h"
#include "ofmo-prof.h"
/** メモリサーバーからデータを取得する関数

  設定されたデフォルトのメモリサーバーからデータを取得する
  */
int ofmo_master_get( MPI_Comm comm_mserv, int data_id,
	int ifrag, double D[] ) {
    static int nfrag, *nfao, *nfatom, nfao_total, nfatom_total;
    static int called = false;
    if ( !called ) {
	int ierr;
	ierr = ofmo_data_get_vals("nfrag nfao nfatom",
		&nfrag, &nfao, &nfatom );
	if ( ierr != 0 ) return -1;
	nfao_total   = ofmo_isum2( nfrag, nfao );
	nfatom_total = ofmo_isum2( nfrag, nfatom );
	called = true;
    }

    int nelem, target,count, ierr = 0;
    int imsg[OFMO_IMSG_SZ], imsg_sz=OFMO_IMSG_SZ;
    int tag_ret=OFMO_TAG_RET, tag_data=OFMO_TAG_DAT;
    MPI_Status status;
    if ( data_id == OFMO_TOTAL_AOPOP || data_id == OFMO_TOTAL_ATPOP ) {
	nelem  = ofmo_get_data_nelems( data_id, 0 );
	target = ofmo_get_data_target( data_id, 0 );
    } else if ( ifrag < 0 ) {
	if ( data_id == OFMO_AOPOP1 || data_id == OFMO_AOPOP2 ) {
	    nelem  = nfao_total;
	    target = ofmo_get_data_target( data_id, 0 );
	} else if ( data_id == OFMO_ATPOP1 || data_id == OFMO_ATPOP2 ) {
	    nelem  = nfatom_total;
	    target = ofmo_get_data_target( data_id, 0 );
	} else if ( data_id == OFMO_ENERGY || data_id == OFMO_ENERGY0 ) {
	    nelem  = nfrag;
	    target = ofmo_get_data_target( data_id, 0 );
	} else {
	    return -1;
	}
    } else {
	nelem  = ofmo_get_data_nelems( data_id, ifrag );
	target = ofmo_get_data_target( data_id, ifrag );
    }
    imsg[OFMO_I_CMD]    = OFMO_GET;
    imsg[OFMO_I_METHOD] = data_id;
    imsg[OFMO_I_SCC]    = ifrag;
    MPI_Send( imsg, imsg_sz, MPI_INT, target, tag_ret, comm_mserv );
    MPI_Recv( D, nelem, MPI_DOUBLE, target, tag_data, comm_mserv,
	    &status );
    MPI_Get_count( &status, MPI_DOUBLE, &count );
    if ( count != nelem ) {
	fprintf( stderr, "ERROR: there are incomsistency\n"
		"  dataid = %d, ifrag = %d\n"
		"  %d elems must be recv, but %d recv\n",
		data_id, ifrag, nelem, count );
	ierr = -1;
    }
    return ierr;
}

/** １つのメモリサーバーにデータを送る（データを更新する）関数

  メモリサーバーのコミュニケータを指定する
 */
int ofmo_master_put( MPI_Comm comm_mserv, int data_id,
	int ifrag, double *src ) {
    static int nfrag, *nfao, *nfatom, nfao_total, nfatom_total;
    static int nao, natom;
    static int called = false;
    if ( !called ) {
	ofmo_data_get_vals("nfrag nfao nfatom nao natom",
		&nfrag, &nfao, &nfatom, &nao, &natom );
	nfao_total   = ofmo_isum2( nfrag, nfao );
	nfatom_total = ofmo_isum2( nfrag, nfatom );
	called = true;
    }
    int nelem, target;
    int imsg[OFMO_IMSG_SZ], imsg_sz=OFMO_IMSG_SZ;
    int tag_ret=OFMO_TAG_RET, tag_data=OFMO_TAG_DAT;
    if ( ifrag >= nfrag ) {
	if ( fp_prof ) {
	    fdbg( fp_prof, "ERROR: Illegal monomer number (%d)\n", ifrag);
	    fflush( fp_prof );
	}
	return -1;
    }
    nelem = 0;
    if ( data_id == OFMO_DENS1 || data_id == OFMO_DENS2 ) {
	if ( ifrag >= 0 ) {
	    int n;
	    n = nfao[ifrag];
	    nelem = (n*n+n)>>1;
	    target = ofmo_get_data_target( data_id, ifrag );
	}
    } else if ( data_id == OFMO_AOPOP1 || data_id == OFMO_AOPOP2 ) {
	nelem = ( ifrag<0 ? nfao_total : nfao[ifrag] );
	target = ofmo_get_data_target( data_id, 0 );
    } else if ( data_id == OFMO_ATPOP1 || data_id == OFMO_ATPOP2 ) {
	nelem = ( ifrag<0 ? nfatom_total : nfatom[ifrag] );
	target = ofmo_get_data_target( data_id, 0 );
    } else if ( data_id == OFMO_DISTA ) {
	nelem = nfrag;
	target = ofmo_get_data_target( data_id, ifrag );
    } else if ( data_id == OFMO_ENERGY || data_id == OFMO_ENERGY0 ) {
	nelem = ( ifrag<0 ? nfrag : 1 );
	target = ofmo_get_data_target( data_id, 0 );
    } else if ( data_id == OFMO_TOTAL_AOPOP ) {
	nelem  = nao;
	target = ofmo_get_data_target( data_id, 0 );
    } else if ( data_id == OFMO_TOTAL_ATPOP ) {
	nelem  = natom;
	target = ofmo_get_data_target( data_id, 0 );
    }
    if ( nelem <= 0 ) {
	if ( fp_prof ) {
	    fdbg( fp_prof,
		    "ERROR: Illegal parameter (ifrag=%d, data ID=%d)\n",
		    ifrag, data_id );
	    fflush( fp_prof );
	}
	return -1;
    }
    imsg[OFMO_I_CMD]    = OFMO_PUT;
    imsg[OFMO_I_METHOD] = data_id;
    imsg[OFMO_I_SCC]    = ifrag;

    MPI_Send( imsg, imsg_sz, MPI_INT, target, tag_ret, comm_mserv );
    MPI_Send( src, nelem, MPI_DOUBLE, target, tag_data, comm_mserv );
    MPI_Recv( MPI_BOTTOM, 0, MPI_BYTE, target, tag_ret, comm_mserv,
	    MPI_STATUS_IGNORE );
    return 0;
}

int ofmo_master_put_all( int nmservs, MPI_Comm comm_mservs[],
	int data_id, int ifrag, double *src ) {
    static int nfrag, *nfao, *nfatom, nfao_total, nfatom_total;
    static int nao, natom;
    static int called = false;
    if ( !called ) {
	ofmo_data_get_vals("nfrag nfao nfatom nao natom",
		&nfrag, &nfao, &nfatom, &nao, &natom );
	nfao_total   = ofmo_isum2( nfrag, nfao );
	nfatom_total = ofmo_isum2( nfrag, nfatom );
	called = true;
    }
    int nelem, target;
    int imsg[OFMO_IMSG_SZ], imsg_sz=OFMO_IMSG_SZ;
    int tag_ret=OFMO_TAG_RET, tag_data=OFMO_TAG_DAT;
    if ( ifrag >= nfrag ) {
	if ( fp_prof ) {
	    fdbg( fp_prof, "ERROR: Illegal monomer number (%d)\n", ifrag);
	    fflush( fp_prof );
	}
	return -1;
    }
    nelem = 0;
    if ( data_id == OFMO_DENS1 || data_id == OFMO_DENS2 ) {
	if ( ifrag >= 0 ) {
	    int n;
	    n = nfao[ifrag];
	    nelem = (n*n+n)>>1;
	    target = ofmo_get_data_target( data_id, ifrag );
	}
    } else if ( data_id == OFMO_AOPOP1 || data_id == OFMO_AOPOP2 ) {
	nelem = ( ifrag<0 ? nfao_total : nfao[ifrag] );
	target = ofmo_get_data_target( data_id, 0 );
    } else if ( data_id == OFMO_ATPOP1 || data_id == OFMO_ATPOP2 ) {
	nelem = ( ifrag<0 ? nfatom_total : nfatom[ifrag] );
	target = ofmo_get_data_target( data_id, 0 );
    } else if ( data_id == OFMO_DISTA ) {
	nelem = nfrag;
	target = ofmo_get_data_target( data_id, ifrag );
    } else if ( data_id == OFMO_ENERGY || data_id == OFMO_ENERGY0 ) {
	nelem = ( ifrag<0 ? nfrag : 1 );
	target = ofmo_get_data_target( data_id, 0 );
    } else if ( data_id == OFMO_TOTAL_AOPOP ) {
	nelem  = nao;
	target = ofmo_get_data_target( data_id, 0 );
    } else if ( data_id == OFMO_TOTAL_ATPOP ) {
	nelem  = natom;
	target = ofmo_get_data_target( data_id, 0 );
    }
    if ( nelem <= 0 ) {
	if ( fp_prof ) {
	    fdbg( fp_prof,
		    "ERROR: Illegal parameter (ifrag=%d, data ID=%d)\n",
		    ifrag, data_id );
	    fflush( fp_prof );
	}
	return -1;
    }
    imsg[OFMO_I_CMD]    = OFMO_PUT;
    imsg[OFMO_I_METHOD] = data_id;
    imsg[OFMO_I_SCC]    = ifrag;

    for ( int mserv_id=0; mserv_id<nmservs; mserv_id++ ) {
	MPI_Send( imsg, imsg_sz, MPI_INT, target, tag_ret,
		comm_mservs[mserv_id] );
	MPI_Send( src, nelem, MPI_DOUBLE, target, tag_data,
		comm_mservs[mserv_id] );
	MPI_Recv( MPI_BOTTOM, 0, MPI_BYTE, target, tag_ret,
		comm_mservs[mserv_id], MPI_STATUS_IGNORE );
    }
    return 0;
}
/* ******** ******* */

static int dold = OFMO_DENS1;
static int dnew = OFMO_DENS2;
static int aopold = OFMO_AOPOP1;
static int aopnew = OFMO_AOPOP2;
static int atpold = OFMO_ATPOP1;
static int atpnew = OFMO_ATPOP2;

int ofmo_get_monomer_density_master( MPI_Comm comm_mserv,
	int ifrag, double D[] ) {
    ofmo_master_get( comm_mserv, dold, ifrag, D );
    return 0;
}

int ofmo_get_monomer_aopop_master( MPI_Comm comm_mserv,
	int ifrag, double aopop[] ) {
    ofmo_master_get( comm_mserv, aopold, ifrag, aopop );
    return 0;
}

int ofmo_get_monomer_atpop_master( MPI_Comm comm_mserv,
	int ifrag, double atpop[] ) {
    ofmo_master_get( comm_mserv, atpold, ifrag, atpop );
    return 0;
}

int ofmo_update_monomer_data_master() {
    int itmp;
    itmp = dold;     dold = dnew;     dnew = itmp;
    itmp = aopold; aopold = aopnew; aopnew = itmp;
    itmp = atpold; atpold = atpnew; atpnew = itmp;
    return 0;
}

void ofmo_reset_monomer_data_master() {
    dold = OFMO_DENS1;
    dnew = OFMO_DENS2;
    aopold = OFMO_AOPOP1;
    aopnew = OFMO_AOPOP2;
    atpold = OFMO_ATPOP1;
    atpnew = OFMO_ATPOP2;
}

/** ファイルからデータを読み込む関数
  */
static int *IBUF = NULL;
static double *DBUF = NULL;
static int NFRAG = 0;

static void dealloc() {
    Free( IBUF );
    Free( DBUF );
    IBUF = NULL;
    DBUF = NULL;
    NFRAG = 0;
}

static int alloc() {
    static int called = false;
    if ( called ) return 0;
    int maxnfao, nfrag, n2, nao, n;
    ofmo_data_get_vals("nfrag maxnfao nao", &nfrag, &maxnfao, &nao );
    n2 = (maxnfao*maxnfao+maxnfao)>>1;
    n = n2;
    if ( nfrag > n ) n = nfrag;
    if ( nao > n )   n = nao;
    IBUF = (int*)malloc( sizeof(int) * nfrag );
    DBUF = (double*)malloc( sizeof(double) * n );
    NFRAG = nfrag;
    atexit( dealloc );
    called = true;
    return 0;
}

int ofmo_read_and_put_all( int nmservs, MPI_Comm comm_mservs[],
	char *fileheader ) {
    int nfrag, *nelems, nread, type, size, i, ifrag;
    FILE *fp;
    char filename[MAXSTRLEN];
    char suffix[6][5] = { "dens", "aop", "atp", "ene", "ene0", "dist" };
    int data_ids[6];
    alloc();
    nelems = IBUF;
    data_ids[0] = OFMO_DENS1;
    data_ids[1] = OFMO_AOPOP1;
    data_ids[2] = OFMO_ATPOP1;
    data_ids[3] = OFMO_ENERGY;
    data_ids[4] = OFMO_ENERGY0;
    data_ids[5] = OFMO_DISTA;
    // 必要なデータのファイルがあるかチェック
    for ( i=0; i<6; i++ ) {
	snprintf( filename, MAXSTRLEN, "%s-%s.dat", fileheader, suffix[i] );
	if ( (fp=fopen( filename, "r" )) == NULL ) return -1;
	fclose( fp );
    }
    for ( i=0; i<6; i++ ) {
	snprintf( filename, MAXSTRLEN, "%s-%s.dat", fileheader, suffix[i] );
	fp=fopen( filename, "r" );
	fread( &type, sizeof(int), 1, fp );
	fread( &nfrag, sizeof(int), 1, fp );
	if ( nfrag != NFRAG ) {
	    dbg("ERROR: Illgel value of nfrag (given=%d, file=%d)\n",
		    NFRAG, nfrag );
	    fclose( fp );
	    return -1;
	}
	if ( type == 1 ) {
	    fread( &size, sizeof(int), 1, fp );
	    for ( ifrag=0; ifrag<nfrag; ifrag++ ) nelems[ifrag] = size;
	} else if ( type == 2 ) {
	    nread = fread( nelems, sizeof(int), nfrag, fp );
	    if ( nread != nfrag ) {
		dbg("ERROR: Unexected EOF in reading nelems"
			" nfrag=%d, nread=%d\n", nfrag, nread );
		fclose( fp );
		return -1;
	    }
	}
	if ( type == 2 || type == 1 ) {
	    for ( ifrag=0; ifrag<nfrag; ifrag++ ) {
		nread = fread( DBUF, sizeof(double), nelems[ifrag],
			fp );
		if ( nread != nelems[ifrag] ) {
		    dbg("ERROR: Bad number of elemnts (%d vs %d )\n",
			    nread, nelems[ifrag] );
		    fclose( fp );
		    return -1;
		}
		ofmo_master_put_all( nmservs, comm_mservs,
			data_ids[i], ifrag, DBUF );
	    }
	} else if ( type == 0 ) {
	    nread = fread( DBUF, sizeof(double), nfrag, fp );
	    if ( nread != nfrag ) {
		dbg("ERROR: Bad number of elemnts (%d vs %d )\n",
			nread, nelems[ifrag] );
		fclose( fp );
		return -1;
	    }
	    ofmo_master_put_all( nmservs, comm_mservs,
		    data_ids[i], -1, DBUF );
	} else {
	    fclose( fp );
	    return -1;
	}
	fclose( fp );
    }	// for ( int i );
    ofmo_reset_monomer_data_master();
    return 0;
}

int ofmo_read_and_put( MPI_Comm comm_mserv, char *fileheader ) {
    int nfrag, *nelems, nread, type, size, ifrag, i;
    FILE *fp;
    char filename[MAXSTRLEN];
    char suffix[6][5] = { "dens", "aop", "atp", "ene", "ene0", "dist" };
    int data_ids[6];
    alloc();
    nelems = IBUF;
    data_ids[0] = OFMO_DENS1;
    data_ids[1] = OFMO_AOPOP1;
    data_ids[2] = OFMO_ATPOP1;
    data_ids[3] = OFMO_ENERGY;
    data_ids[4] = OFMO_ENERGY0;
    data_ids[5] = OFMO_DISTA;
    // 必要なデータのファイルがあるかチェック
    for ( i=0; i<6; i++ ) {
	snprintf( filename, MAXSTRLEN, "%s-%s.dat", fileheader, suffix[i] );
	if ( (fp=fopen( filename, "r" )) == NULL ) return -1;
	fclose( fp );
    }
    for ( i=0; i<6; i++ ) {
	snprintf( filename, MAXSTRLEN, "%s-%s.dat", fileheader, suffix[i] );
	fp=fopen( filename, "r" );
	fread( &type, sizeof(int), 1, fp );
	fread( &nfrag, sizeof(int), 1, fp );
	if ( nfrag != NFRAG ) {
	    dbg("ERROR: Illgel value of nfrag (given=%d, file=%d)\n",
		    NFRAG, nfrag );
	    fclose( fp );
	    return -1;
	}
	if ( type == 1 ) {
	    fread( &size, sizeof(int), 1, fp );
	    for ( ifrag=0; ifrag<nfrag; ifrag++ ) nelems[ifrag] = size;
	} else if ( type == 2 ) {
	    nread = fread( nelems, sizeof(int), nfrag, fp );
	    if ( nread != nfrag ) {
		dbg("ERROR: Unexected EOF in reading nelems"
			" nfrag=%d, nread=%d\n", nfrag, nread );
		fclose( fp );
		return -1;
	    }
	}
	if ( type == 2 || type == 1 ) {
	    for ( ifrag=0; ifrag<nfrag; ifrag++ ) {
		nread = fread( DBUF, sizeof(double), nelems[ifrag],
			fp );
		if ( nread != nelems[ifrag] ) {
		    dbg("ERROR: Bad number of elemnts (%d vs %d )\n",
			    nread, nelems[ifrag] );
		    fclose( fp );
		    return -1;
		}
		ofmo_master_put( comm_mserv, data_ids[i], ifrag, DBUF );
	    }
	} else if ( type == 0 ) {
	    nread = fread( DBUF, sizeof(double), nfrag, fp );
	    if ( nread != nfrag ) {
		dbg("ERROR: Bad number of elemnts (%d vs %d )\n",
			nread, nelems[ifrag] );
		fclose( fp );
		return -1;
	    }
	    ofmo_master_put( comm_mserv, data_ids[i], -1, DBUF );
	} else {
	    fclose( fp );
	    return -1;
	}
	fclose( fp );
    }	// for ( int i );
    ofmo_reset_monomer_data_master();
    return 0;
}

/** 計算で得られた（FMO計算続行に必要な）モノマーに間するデータを
  ファイルに出力する関数
  */
int ofmo_get_and_write( MPI_Comm comm_mserv, char *fileheader ) {
    int nfrag, *nelems, type, ifrag, i;
    FILE *fp;
    char filename[MAXSTRLEN];
    char suffix[6][5] = { "dens", "aop", "atp", "ene", "ene0", "dist" };
    int types[6] = { 2, 2, 2, 0, 0, 1 };
    int data_ids[6];
    alloc();
    data_ids[0] = dold;
    data_ids[1] = aopold;
    data_ids[2] = atpold;
    data_ids[3] = OFMO_ENERGY;
    data_ids[4] = OFMO_ENERGY0;
    data_ids[5] = OFMO_DISTA;
    nfrag = ofmo_get_nfrag();
    for ( i=0; i<6; i++ ) {
	snprintf( filename, MAXSTRLEN, "%s-%s.dat", fileheader, suffix[i] );
	if( (fp=fopen( filename, "w" )) == NULL ) {
	    dbg("ERROR: Failure in open file %s\n", filename );
	    return -1;
	}
	type = types[i];
	nelems = ofmo_getadd_data_nelems( data_ids[i] );
	fwrite( &type, sizeof(int), 1, fp );
	fwrite( &nfrag, sizeof(int), 1, fp );
	if      ( type == 1 ) fwrite( &nfrag, sizeof(int), 1, fp );
	else if ( type == 2 ) fwrite( nelems, sizeof(int), nfrag, fp );

	if ( type == 2 || type == 1 ) {
	    for ( ifrag=0; ifrag<nfrag; ifrag++ ) {
		ofmo_master_get( comm_mserv, data_ids[i], ifrag, DBUF );
		fwrite( DBUF, sizeof(double), nelems[ifrag], fp );
	    }
	} else if ( type == 0 ) {
	    ofmo_master_get( comm_mserv, data_ids[i], -1, DBUF );
	    fwrite( DBUF, sizeof(double), nfrag, fp );
	} else {
	    fclose( fp );
	    return -1;
	}
	fclose( fp );
    }	// for ( int i );
    return 0;
}

/** あるメモリサーバーのデータを別のメモリサーバーにコピーする関数

  マスタープロセスを経由して送る
*/
int ofmo_copy_mserv_data( MPI_Comm comm_mserv_src,
	MPI_Comm comm_mserv_dst ) {
    int nfrag, ndata, data_id, ifrag;
    alloc();
    ndata = ofmo_get_ndata();
    nfrag = ofmo_get_nfrag();
    for ( data_id=0; data_id<ndata; data_id++ ) {
	if ( data_id == OFMO_DENS1 || data_id == OFMO_DENS2 ||
		data_id == OFMO_AOPOP1 || data_id == OFMO_AOPOP2 ||
		data_id == OFMO_ATPOP1 || data_id == OFMO_ATPOP2 ||
		data_id == OFMO_DISTA ) {
	    for ( ifrag=0; ifrag<nfrag; ifrag++ ) {
		ofmo_master_get( comm_mserv_src, data_id, ifrag, DBUF );
		ofmo_master_put( comm_mserv_dst, data_id, ifrag, DBUF );
	    }
	} else {
	    ofmo_master_get( comm_mserv_src, data_id, -1, DBUF );
	    ofmo_master_put( comm_mserv_dst, data_id, -1, DBUF );
	}
    }
    return 0;
}
