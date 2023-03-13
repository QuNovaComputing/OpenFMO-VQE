#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ofmo-def.h"

/* INPUT_BUFFER_SIZE = buffer size in input (in MB) */
#include <zlib.h>

#ifndef INPUT_BUFFER_SIZE
#define INPUT_BUFFER_SIZE 8
#endif

#ifndef OUTPUT_BUFFER_SIZE
#define OUTPUT_BUFFER_SIZE 1
#endif

static int tmp_data_size      = 0;
static unsigned char *tmp_data    = NULL;

static z_stream z;

static void finalize_buf() {
    Free( tmp_data );
    tmp_data_size = 0;
}

static int init_buf() {
    static int called = false;
    char *p;
    if ( called ) return 0;
    if ( (p=getenv("OFMO_TMP_DATA_SIZE")) != NULL) {
	tmp_data_size = atoi( p );
    }
    if ( tmp_data_size < 1 ) tmp_data_size = OUTPUT_BUFFER_SIZE;
    tmp_data_size *= 1024 * 1024;
    tmp_data    = (unsigned char*)malloc( tmp_data_size );
    if ( tmp_data == NULL ) return -1;
    atexit( finalize_buf );
    called = true;
    return 0;
}

/** 圧縮したデータを保存する領域を確保する
  (別ファイルでもよいかも・・・）
   */
unsigned char * ofmo_alloc_compressed_data( int *max_comped_data_size ) {
    int MAX_COMP_DATA_SIZE=0;
    char *p;
    unsigned char *comped_data;
    if ( (p=getenv("OFMO_COMP_DATA_SIZE")) != NULL ) {
	MAX_COMP_DATA_SIZE = atoi( p );
    }
    if ( MAX_COMP_DATA_SIZE < 1 ) MAX_COMP_DATA_SIZE = INPUT_BUFFER_SIZE;
    MAX_COMP_DATA_SIZE *= 1024*1024;	// convert unit from MB to byte
    comped_data = (unsigned char*)malloc( MAX_COMP_DATA_SIZE );
    *max_comped_data_size = MAX_COMP_DATA_SIZE;
    if ( comped_data == NULL ) *max_comped_data_size = 0;
    return comped_data;
}

/** ファイルを読み込んで圧縮し、結果をバッファに書き込む
 * */
int ofmo_compress_data( const char filename[], unsigned char comped_data[],
	const int max_comped_data_size ) {
    int count, flush, status;
    FILE *fp;

    if ( (fp=fopen( filename, "r" )) == NULL ) {
	dbg("Failure in opening file (%s)\n", filename );
	return -1;
    }
    if ( init_buf() != 0 ) {
	dbg("Failure in allocation of temporary"
		" array used in data compress\n");
	return -1;
    }

    z.zalloc = Z_NULL;
    z.zfree = Z_NULL;
    z.opaque = Z_NULL;

    if (deflateInit(&z, Z_DEFAULT_COMPRESSION) != Z_OK) {
	dbg( "deflateInit: %s\n", (z.msg) ? z.msg : "???");
	return -1;
    }
    z.avail_in  = 0;
    z.next_out  = comped_data;
    z.avail_out = max_comped_data_size;

    flush = Z_NO_FLUSH;
    while (1) {
	if (z.avail_in == 0) {
	    z.next_in = tmp_data;
	    z.avail_in = fread( tmp_data, 1, tmp_data_size, fp );
	    if ( z.avail_in < tmp_data_size ) flush = Z_FINISH;
	}
	status = deflate(&z, flush);
	if (status == Z_STREAM_END) break;
	if (status != Z_OK) {
	    dbg( "deflate: %s\n", (z.msg) ? z.msg : "???");
	    return -1;
	}
	if (z.avail_out == 0) {
	    dbg("Too large compressed data size ( > %d byte )\n",
		    max_comped_data_size );
	}
    }
    count = max_comped_data_size - z.avail_out;

    if (deflateEnd(&z) != Z_OK) {
	dbg( "deflateEnd: %s\n", (z.msg) ? z.msg : "???");
	return -1;
    }
    fclose( fp );
    return count;
}

/** バッファに格納された圧縮ファイルを展開して、出力ファイルに
 * 書き込む
 * */
int ofmo_decompress_data( const char* filename,
	unsigned char* comped_data, int comped_data_size ) {
    FILE *fp;
    int count, status;

    if ( (fp=fopen( filename, "w")) == NULL ) {
	dbg("Failure in opening file (%s)\n", filename );
	return -1;
    }
    if ( init_buf() != 0 ) {
	dbg("Failure in allocation of temporary"
		" array used in data compress\n");
	return -1;
    }
    z.zalloc = Z_NULL;
    z.zfree = Z_NULL;
    z.opaque = Z_NULL;

    z.next_in = Z_NULL;
    z.avail_in = 0;
    if (inflateInit(&z) != Z_OK) {
	fprintf(stderr, "inflateInit: %s\n", (z.msg) ? z.msg : "???");
	return -1;
    }
    z.next_out = tmp_data;
    z.avail_out = tmp_data_size;
    status = Z_OK;

    z.avail_in = comped_data_size;
    z.next_in  = comped_data;
    while ( status != Z_STREAM_END ) {
	status = inflate(&z, Z_NO_FLUSH);
	if (status == Z_STREAM_END) break;
	if (status != Z_OK) {
	    fprintf(stderr, "inflate: %s\n", (z.msg) ? z.msg : "???");
	    return -1;
	}
	if (z.avail_out == 0) {
	    if (fwrite( tmp_data, 1, tmp_data_size, fp)
		    != tmp_data_size ) {
		fprintf(stderr, "Write error\n");
		return -1;
	    }
	    z.next_out = tmp_data;
	    z.avail_out = tmp_data_size;
	}
    }
    if ( (count = tmp_data_size - z.avail_out) != 0) {
	if (fwrite(tmp_data, 1, count, fp) != count) {
	    fprintf(stderr, "Write error\n");
	    return -1;
	}
    }

    fclose( fp );

    if (inflateEnd(&z) != Z_OK) {
	fprintf(stderr, "inflateEnd: %s\n", (z.msg) ? z.msg : "???");
	return -1;
    }
    return 0;
}
