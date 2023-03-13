
#ifndef _OFMO_MASTER_MAIN_H_
#define _OFMO_MASTER_MAIN_H_


#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdarg.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <getopt.h>
#include <assert.h>
#include <libgen.h>
#include <limits.h>
#include <errno.h>

#include <math.h>

#include <mpi.h>

#ifdef _OPENMP
#include <omp.h>
#else
#include "omp-dummy.h"
#endif


#include "ofmo-def.h"
#include "ofmo-prof.h"
#include "ofmo-data.h"

#include "ofmo-string.h"
#include "ofmo-misc.h"

#include "ofmo-scf.h"
#include "ofmo-integ.h"
#include "ofmo-data-struct.h"

#include <falanx.h>
#include <falanx_log.h>
#include <datastore.h>


#ifdef FJ_MAIN
#define MASTER_main MAIN__
#else
#define MASTER_main main
#endif

#define MAXTOKLEN MAXSTRLEN


struct ofmo_master_config_st {

    int ngroup;
    int ngroup_dm; /* #groups for dimer calculation */
    int nioprocs;
    int niogroup;
    int nmaxprocs;
    int group_size;
    int nfrag;
    int nbody;
    bool need_init_dens;
    char input[MAXSTRLEN];
    char header[MAXSTRLEN];
    char port_prefix[MAXSTRLEN];
    char local_input_dir[MAXSTRLEN];
    char dens[MAXSTRLEN];

    int ndev;
    long eribfsz;

    int bufsz_aop;
    int bufsz_atp;
    double* atpop_total;
    double* aopop_total;
    double* aopop_frg;
    double* atpop_frg;
    int* fsao2tuao;
    int* fatom2tatom;

    int** issued_job_list;
};
typedef struct ofmo_master_config_st ofmo_master_config_t;


/* 外部で定義された関数 */
extern int ofmo_init(const char* filename, MPI_Comm comm);
extern unsigned char* ofmo_alloc_compressed_data(int*);
extern int ofmo_compress_data(const char filename[], unsigned char comped_data[], const int max_comped_data_size);
extern int ofmo_frag_init();
extern int ofmo_create_worker_prof_name(const char*, const char*, const int, char*, size_t); 

#include "ofmo-task-util.h"
#include "ofmo-worker-task.h"
#include "ofmo-datastore.h"


#endif

