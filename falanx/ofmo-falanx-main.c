/*
 * falanx版OpenFMO マスタープロセス
 */

#include "ofmo-falanx-main.h"
#ifdef USE_CUDA
#include "cuda/cuda-drv.h"
#include "cuda/cuda-ifc4c-calc.h"
#endif
#include "ofmo-tlog.h"

#define IGNORE_DDV             (NULL)
#define IGNORE_MONOMER_LIST    (NULL)
#define REQUIRE_RECONFIG(CONF) ((CONF).ngroup_dm > 0)
#define CHECK_STRTOL(L) ((errno == ERANGE && ((L) == LONG_MAX || (L) == LONG_MIN)) || (errno != 0 && (L) == 0))
#define RANGE_IN_INT(L) (((L) <= INT_MAX) && ((L) >= INT_MIN))



/* -----------------------------------------------------
 * モノマーをAO数の大きい順に並び替えたリストの作成
 * ----------------------------------------------------- */
static int* frag_order = NULL;

/** ２つの整数を比較する関数 **/
static int comp2(const void* p1, const void* p2) {
    return (*(int*) p2) - (*(int*) p1);
}

static void dealloc_frag_order() {
    if (frag_order != NULL) { free(frag_order); }
    frag_order = NULL;
}

static int
ofmo_init_frag_order(void)
{
    static int called = false;
    int i, i2, nfrag;
    int* nfao;

    if (called) {
        return 0;
    }

    if (ofmo_data_get_vals("nfrag nfao", &nfrag, &nfao) != 0) {
        dbg("error\n");
        return -1;
    }
    frag_order = (int*) malloc(sizeof(int) * nfrag * 2);
    if (frag_order == NULL) {
        return -1;
    }
    for (i = 0, i2 = 0; i < nfrag; i++, i2 += 2) {
        frag_order[i2 + 0] = nfao[i];
        frag_order[i2 + 1] = i;
    }
    qsort(frag_order, nfrag, sizeof(int) * 2, comp2);
    for (i = 0, i2 = 0; i < nfrag; i++, i2 += 2) {
        frag_order[i] = frag_order[i2 + 1];
    }
    atexit(dealloc_frag_order);
    called = true;
    return 0;
}


/* ----- SCC収束条件を決める関数 ----- */
static int
ofmo_get_monomer_scfconv(double dE)
{
    static int scfconv_now = 2, called = false, scfconv;
    static double sccconvd;
    static int nit=0;
    int scfconv_old;
    double     de;
    if (!called) {
        int sccconv;
        double scfconvd;
        ofmo_data_get_vals("scfconv sccconv", &scfconv, &sccconv );
        called = true;
        scfconvd = pow(10,-scfconv);
        sccconvd = pow(10,-sccconv);
        if (scfconvd>sccconvd)
            printf("Warning: CONV in $SCF (%8.2e) > CONV in $FMOPRP (%8.2e)!\n",
                   scfconvd, sccconvd);
    }
    de = fabs(dE);
    scfconv_old = scfconv_now;
    nit++;
#if 0
    if      (scfconv_now <= 4 && de < 1.e-4) { scfconv_now = scfconv; }
    else if (scfconv_now <= 3 && de < 1.e-3) { scfconv_now = 4; }
    else if (scfconv_now <= 2 && de < 1.e-1) { scfconv_now = 3; }
#else
    if      (scfconv_now == 2 && (nit>5 || de < 1.e-1)) scfconv_now+=2;
    else if (scfconv_now >= 3 && (nit>3 || de < pow(10,-scfconv_now))) scfconv_now+=2;
#endif
    if (scfconv_now > scfconv || de < 10*sccconvd) { scfconv_now = scfconv; }
    if (scfconv_now != scfconv_old) nit=1;
    return scfconv_now;
}

/*
 * 配列に対するaccumulate処理を行う関数
 */
static void
ofmo_accumulate(
    const int n, const double data[], const int index[],
    double target[])
{
    int i;
    for (i = 0; i < n; i++) {
        target[index[i]] += data[i];
    }
}

/*
 * populationデータを取り出す.
 */
static int
ofmo_get_pop_data(
    falanx_context_t* context, char* key, size_t keysize, const int mode,
    int* fnao, double aopop_frg[], int fsao2tuao[],
    int* fnat, double atpop_frg[], int fatom2tatom[])
{
    int rc = -1;
    char outpopkey[MAX_KEY_LENGTH];
    data_store_t* ds = NULL;
    ofmo_scf_output_ao_t aopop;
    ofmo_scf_output_at_t atpop;

    *fnao = *fnat = 0;

    ds = falanx_context_get_data_store(context);
    assert(ds != NULL);

    if (!(mode == OFMO_BOTH_POPS || mode == OFMO_AOPOP_ONLY || mode == OFMO_ATPOP_ONLY)) {
        return -1; // invalid mode
    }

    if (mode == OFMO_BOTH_POPS || mode == OFMO_AOPOP_ONLY) {
        /* write directly */
        aopop.daopop = aopop_frg;
        aopop.fsao2tuao = fsao2tuao;

        ofmo_scf_get_ao_pop_key(key, outpopkey, MAX_KEY_LENGTH);
        rc = ofmo_master_get_scf_ao_pop(ds, outpopkey, strlen(outpopkey)+1, &aopop);
        if (rc == 0) {
            *fnao = aopop.fnao;
        }
    }
    if (mode == OFMO_BOTH_POPS || mode == OFMO_ATPOP_ONLY) {
        /* write directly */
        atpop.datpop = atpop_frg;
        atpop.fatom2tatom = fatom2tatom;

        ofmo_scf_get_at_pop_key(key, outpopkey, MAX_KEY_LENGTH);
        rc = ofmo_master_get_scf_at_pop(ds, outpopkey, strlen(outpopkey)+1, &atpop);
        if (rc == 0) {
            *fnat = atpop.fnat;
        }
    }
    return rc;
}

/*
  Copy string from src to dst if size fit (strlen(src) < size).
  The rest part is filled by null.  *dst is returned.
  Otherwise, NULL is returned with dst[] = "\0...".
 */
static char *safe_strncpy(char *dst, const char *src, int size)
{
    int len;

    memset(dst, 0, size);
    if (NULL == src) return NULL;
    len = strlen(src);
    if (len >= size) return NULL;
    memcpy(dst, src, len);

    return dst;
}



static void
ofmo_init_population_data(
    ofmo_master_config_t* config,
    int nao, int natom, int maxnfao, int maxnfatom)
{
    config->bufsz_aop   = maxnfao * config->nbody;
    config->bufsz_atp   = maxnfatom * config->nbody;
    config->atpop_total = (double*) calloc(natom, sizeof(double));
    config->aopop_total = (double*) calloc(nao, sizeof(double));
    config->aopop_frg   = (double*) malloc(sizeof(double) * config->bufsz_aop);
    config->atpop_frg   = (double*) malloc(sizeof(double) * config->bufsz_atp);
    config->fsao2tuao   = (int*) malloc(sizeof(int) * config->bufsz_aop);
    config->fatom2tatom = (int*) malloc(sizeof(int) * config->bufsz_atp);
}

static void
ofmo_finalize_config(ofmo_master_config_t* config)
{
    Free(config->atpop_total);
    Free(config->aopop_total);
    Free(config->aopop_frg);
    Free(config->atpop_frg);
    Free(config->fsao2tuao);
    Free(config->fatom2tatom);
}

static int
ofmo_resolve_group_size(int nmaxprocs, int niogroup, int ngroup)
{
    return (nmaxprocs - niogroup - 2) / ngroup;
}

/*
 * ワーカーグループのサイズを求める.
 */
static int
ofmo_master_resolve_group_size(ofmo_master_config_t* config)
{
    return ofmo_resolve_group_size(config->nmaxprocs, config->niogroup, config->ngroup);
}

/*
 * dimer計算用ワーカーグループのサイズを求める.
 */
static int
ofmo_master_resolve_group_size_dimer(ofmo_master_config_t* config)
{
    if (config->ngroup_dm > 0) {
        return ofmo_resolve_group_size(config->nmaxprocs, config->niogroup, config->ngroup_dm);
    } else {
        return ofmo_resolve_group_size(config->nmaxprocs, config->niogroup, config->ngroup);
    }
}

static void
ofmo_master_timestamp(const char *message)
{
    static double ET0, et0;
    double t;

    t = MPI_Wtime();

    if (NULL == message) {
        et0 = ET0 = t;
        return;
    }

    printf("##T %38s : lap time= %10.6f ( total etime= %10.6f )\n", 
           message, t - et0, t - ET0);
    fflush(stdout);
    et0 = t;

    return;
}

static int
ofmo_master_resolve_data_struct(ofmo_master_config_t* config)
{
    const int NDATA = 11;

    /* FMO計算で用いるデータ構造の決定 */
    ofmo_data_get_vals("nfrag", &config->nfrag);

    ofmo_make_data_struct(NDATA, config->nfrag, config->niogroup);

    printf(" # of data = %d\n", ofmo_get_ndata());
    printf("nfrag = %d\n", config->nfrag);

    ofmo_master_timestamp("Finish Initialization");

    return 0;
}

void show_help(const char *myname, const int verbose)
{
  fprintf(stderr,"Usage: %s [options] [input [InitDens]]\n", myname);
  fprintf(stderr,"  -ng #[:#]: # groups\n");
  fprintf(stderr,"  -np #: # total MPI procs\n");
  fprintf(stderr,"  -B #: buffer size / proc (MB, default: %d)\n", 512);
  fprintf(stderr,"  -v: verbose\n");
  fprintf(stderr,"  -h: show this help\n");
#ifdef USE_CUDA
  fprintf(stderr," Options for GPGPU:\n");
  fprintf(stderr,"  -d #: # devices (default:1)\n");
#endif
}


/*
 * ngroupオプションをパースする.
 * --ngroup/--ng monomer[:dimer]
 */
static int
parse_ngroup(const char* opt_str, int* outMonomer, int* outDimer)
{
    long l;
    char* endptr;

    l = strtol(opt_str, &endptr, 10);
    if (CHECK_STRTOL(l) || (!RANGE_IN_INT(l)) || (opt_str == endptr)) {
        return -1;
    }
    *outMonomer = l;

    if (*endptr == ':') {
        opt_str = endptr+1;
        l = strtol(opt_str, &endptr, 10);
        if (CHECK_STRTOL(l) || (!RANGE_IN_INT(l)) || (opt_str == endptr)) {
            return -1;
        }
        *outDimer = l;

        return 0;
    } else if (*endptr == '\0') {
        *outDimer = -1;
        return 0;
    }
    return -1;
}

static void
ofmo_master_show_info(ofmo_master_config_t* config)
{
    int  resultlen;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    int total_procs = 2 + config->niogroup * config->nioprocs
                        + config->ngroup * config->group_size;
    /* 2 = #client + #task_queue */
    /* nioprocs = 1 */

    printf("nmaxprocs    = %d\n",          config->nmaxprocs);
    printf("niogroup     = %d\n",          config->niogroup);   /* DS数 */
    printf("ngroup       = %d\n",          config->ngroup);     /* グループ数 */
    printf("group size   = %d\n",          config->group_size); /* グループサイズ */
    if (config->ngroup_dm > 0) {
        int dm_grp_sz = ofmo_master_resolve_group_size_dimer(config);
        printf("in Dimer Calculation.\n");
        printf("ngroup       = %d\n", config->ngroup_dm);
        printf("group size   = %d\n", dm_grp_sz);
    }
    printf("# of threads = %d\n",          omp_get_max_threads());
    printf("nintic       = %ld\n",         config->eribfsz );
#ifdef USE_CUDA
    printf("# GPGPU      = %d\n",          config->ndev );
#endif
    printf("input file name       = %s\n", config->input);
    printf("prefix of output file = %s\n", config->header);
    // printf("port name prefix      = %s\n", config->port_prefix);
    printf("local input dir.      = %s\n", config->local_input_dir);
    printf("# of total invoked procs = %d (=%diog*%diop+2+%dgr*%dpr)\n", total_procs,
           config->niogroup, config->nioprocs, config->ngroup, config->group_size);


    // debug
    MPI_Get_processor_name(hostname, &resultlen);
    printf("master host name = %s\n", hostname);
    printf("-----------------------------------------------------\n");
}

static int
ofmo_init_config(
    ofmo_master_config_t* config,
    MPI_Comm comm,
    int      argc,
    char*    argv[])
{
    int ierr;
    int rank;
    int c;
    char *endptr;
    extern char *optarg;
    extern int optind, opterr;
    int help = 0;
    int np = -1;
    char myname[MAXSTRLEN];
    int ngroupA, nioprocsA, niogroupA, iopos, master, verbose;
    static struct option long_options[] = {
        {"ngroup",    1, 0, 0},
        {"ng",        1, 0, 0},
        {"niogroup",  1, 0, 1},			//TI: used as # of DS ranks??
        {"nioprocs",  1, 0, 2},
        {"master",    1, 0, 3},
        {"iopos",     1, 0, 4},
        {"nmaxprocs", 1, 0, 'n'},
        {"np",        1, 0, 'n'},
        {"buffer",    1, 0, 'B'},
        // {"bindir",    1, 0, 5},
        {"scrdir",    1, 0, 6},
        {"ndev",      1, 0, 'd'},
        {"help",      0, 0, 'h'},
        {"verbose",   0, 0, 'v'},
        {0, 0, 0, 0}
    };


    MPI_Comm_rank(comm, &rank);

    ofmo_master_timestamp(NULL);
    config->atpop_total = NULL;
    config->aopop_total = NULL;
    config->aopop_frg = NULL;
    config->atpop_frg = NULL;
    config->fsao2tuao = NULL;
    config->fatom2tatom = NULL;
    config->issued_job_list = NULL;
    
    /* 入力データ、環境変数などの読み込み
     * マスタープロセスだけが、直接、これらのデータにアクセスできると仮定
     * 入力ファイル名の取得
     */
    safe_strncpy(config->input,            getenv("OFMO_INPUT"),           MAXSTRLEN);
    safe_strncpy(config->dens,             getenv("OFMO_FILE_NAME"),       MAXSTRLEN);
    safe_strncpy(config->local_input_dir,  getenv("OFMO_LOCAL_INPUT_DIR"), MAXSTRLEN);
    safe_strncpy(config->header,           getenv("OFMO_HEADER"),          MAXSTRLEN);

    /* 起動したプロセス数の取得. */
    MPI_Comm_size(MPI_COMM_WORLD, &config->nmaxprocs);

    verbose = 0;
    help = 0;
    np = -1;
    config->ngroup = config->niogroup = config->nioprocs = -1;
    ngroupA = niogroupA = nioprocsA = -1;
    master = -1;
    iopos = -2;

    config->eribfsz = -1;
    config->ndev    = 0;

    safe_strncpy(myname, basename(argv[0]), MAXSTRLEN);
    config->local_input_dir[0] = '\0';

    /* コマンドラインオプションを解析する */
    while ((c=getopt_long_only(argc, argv, "B:d:n:hv", long_options, NULL))!=-1) {
      switch(c) {
        case 0:
          if (parse_ngroup(optarg, &config->ngroup, &config->ngroup_dm) < 0)
            help = -1;
          ngroupA = config->ngroup;
          break;
        case 1:
          config->niogroup = strtol(optarg, &endptr, 10);
          if (endptr == optarg) help = -1;
          //if (config->niogroup <= 1) help = -1;
          niogroupA = config->niogroup;
          break;
        case 2:
          config->nioprocs = strtol(optarg, &endptr, 10);
          if (endptr == optarg) help = -1;
          //if (config->nioprocs <= 0) help = -1;
          if (config->nioprocs != 1) help = -1;
          nioprocsA = config->nioprocs;
          break;
        case 3:
          master = strtol(optarg, &endptr, 10);
          if (endptr == optarg) help = -1;
          if (master < 0) help = -1;
          break;
        case 4:
          iopos = strtol(optarg, &endptr, 10);
          if (endptr == optarg) help = -1;
          if (iopos < -1) help = -1;
          break;
	/*
        case 5:
          safe_strncpy(bindir, optarg, MAXSTRLEN);
          if (bindir[0]=='\0') help = -1;
          break;
	*/
        case 6:
          safe_strncpy(config->local_input_dir, optarg, MAXSTRLEN);
          if (config->local_input_dir[0] == '\0') help = -1;
          break;
        case 'n':
          np = strtol(optarg, &endptr, 10);
          if (endptr == optarg) help = -1;
          if (np <= 0) help = -1;
          break;
        case 'B':
          config->eribfsz = strtol(optarg, &endptr, 10);
          if (endptr == optarg) help = -1;
          if (config->eribfsz < 0) help = -1;
          break;
        case 'd':
          config->ndev = strtol(optarg, &endptr, 10);
#ifndef USE_CUDA
          config->ndev = 0;
#endif
          if (endptr == optarg) help = -1;
          if (config->ndev < 0) help = -1;
          break;
        case 'v':
          verbose = 1;
          break;
        case 'h':
        default:
          help = 1;
          break;
      }
    }
    argc -= optind;
    argv += optind;

    if (help != 0) {
        if (rank == 0) {
            show_help(myname, verbose);
        }
        MPI_Finalize();
        return 1;
    }

    if (argc > 0)
        safe_strncpy(config->input, argv[0], MAXSTRLEN);

    if (argc > 1)
        safe_strncpy(config->dens,  argv[1], MAXSTRLEN);

    if ( config->input[0]=='\0' ) {
        if (rank == 0) {
            dbg("Input file is not specified.\n");
            show_help(myname, verbose);
        }
        MPI_Finalize();
        return 2;
    }

    // if (np > 0) nmaxprocs = np;		//TI: this is not relevant in falanx; nmaxprocs == # of ranks
    if ( config->nmaxprocs < 1 ) {
        if (rank == 0) {
            dbg("Illegal number of procs. (%d)\n", config->nmaxprocs );
            show_help(myname, verbose);
        }
        MPI_Finalize();
        return 4;
    }

#if 0
    /* メモリサーバーサービス名のprefix名取得 */
    if (NULL == safe_strncpy( config->port_prefix, getenv("OFMO_PORT_NAME_PREFIX"), MAXSTRLEN )) {
        dbg("OFMO_PORT_NAME_PREFIX is not set\n");
        MPI_Abort(MPI_COMM_WORLD, 5);
    }
#endif

    /* プロファイルデータファイル名のヘッダ取得 */
    if ( config->header[0] == '\0') {
        int pos = strcspn(config->input, ".");
        strncpy(config->header, config->input, pos);
        config->header[pos] = '\0';
    }
    /* ローカルな入力ファイルの置き場所 */
    if ( config->local_input_dir[0] == '\0' ) {
        strncpy(config->local_input_dir, "./", MAXSTRLEN);
    }

    /* 入力データを読み込む */
    if ( ofmo_init(config->input, comm) != 0 ) {
        dbg("failure in input from file (%s)\n", config->input);
        MPI_Finalize();
        return 6;
    }
    ierr = ofmo_data_get_vals("ngroup niogroup nioprocs nbody",
                                         &config->ngroup,
                                       &config->niogroup,
                                      &config->nioprocs,
                                      &config->nbody);
    if ( ierr != 0 ) {
        dbg("failure in get from database\n");
        MPI_Finalize();
        return 66;
    }

    /* コマンドラインオプションが設定されていれば、それを優先する. */
    if (ngroupA > 0)
        config->ngroup = ngroupA;

    if (niogroupA > 0)
        config->niogroup = niogroupA;

    if (nioprocsA > 0)
        config->nioprocs = nioprocsA;


    /* 値を更新する */
    ierr = ofmo_data_put_vals("ngroup niogroup nioprocs",
                                         &config->ngroup,
                                       &config->niogroup,
                                      &config->nioprocs);
    if (master >= 0)
        ofmo_data_put_vals("master", &master);

    if (iopos >= -1)
        ofmo_data_put_vals("iopos", &iopos);

    if ( config->eribfsz < 0 ) {
        ierr = ofmo_data_get_vals("nintic", &config->eribfsz );
        if (ierr != 0)
            config->eribfsz = 0;
    }

    /* 各ワーカーのプロセス数（グループサイズ）の計算 */
    config->group_size = ofmo_master_resolve_group_size(config);

    if (rank == 0) {
        ofmo_master_show_info(config);
    }
    
    if (config->niogroup<=1) {
        if (rank==0) {dbg("niogroup must be >=2\n");};
        MPI_Finalize();
        return 67;
    }
    if (config->nioprocs!=1) {
        if (rank==0) {dbg("nioprocs must be 1\n");};
        MPI_Finalize();
        return 68;
    }
    if (config->group_size < 1) {
        if (rank == 0) { dbg("Illgel group size (%d)\n", config->group_size ); }
        MPI_Finalize();
        return 7;
    }
    if (config->ngroup_dm > 0) {
      int dm_grp_sz = ofmo_master_resolve_group_size_dimer(config);
      if (dm_grp_sz < 1) {
        if (rank == 0) { dbg("Illgel group size for dimer (%d)\n", dm_grp_sz ); }
        MPI_Finalize();
        return 8;
      }
    }

    return 0;
}

static void
ofmo_init_job_list(ofmo_master_config_t* config)
{
    int wid;
    int nworkers = config->ngroup;

    config->issued_job_list = ofmo_alloc_imatrix(nworkers, (MAXNJOB * 2));
    for (wid = 0; wid < nworkers; wid++) {
        for (int i = 0; i < (MAXNJOB * 2); i++) {
            config->issued_job_list[wid][i] = -1;
        }
    }
    // モノマーをAO数の大きい順にソート
    ofmo_init_frag_order();
}

/*
 * タスクの結果を待つ.
 *
 * @return タスクが正常終了の場合、FLX_SUCCESS、かつ outKey に出力キーが書き込まれる.
 *         タスクが無い場合 FLX_NO_WAIT_TASK が返る.
 *         タスクが失敗した場合、FLX_ERR_PROC_FAILED かつ、outKeyに入力キーが書き込まれる.
 *         それ以外はエラーが返る.
 */
static int
ofmo_falanx_wait_response(
    falanx_context_t* context,
    char* outKey,
    size_t outKeyLim,
    size_t* outKeySize,
    falanx_id_t* outId)
{
    int rc;
    falanx_id_t id;

    rc = falanx_task_wait_any(context, &id);
    if (rc == FLX_SUCCESS) {
        rc = falanx_get_result_value(context, id, outKey, outKeyLim, outKeySize);
        if (outId != NULL) {
            *outId = id;
        }
    }
    return rc;
}


/* --------------------------------------------------
 * 初期密度行列計算
 * --------------------------------------------------
 */
static void
ofmo_master_init_dens(
    ofmo_master_config_t* config,
    falanx_context_t*     context)
{
    int rc = 0;
    int data_ids[6];

    config->need_init_dens = true;

    // ファイルからデータを読み込む
    // if ((p = getenv("OFMO_FILE_NAME")) != NULL) {
    if (config->dens[0] != '\0') {
        rc = ofmo_read_and_put_all(config, context, config->dens);
        if (rc == 0) {
            config->need_init_dens = false;
        }
    }

    if (config->need_init_dens) {
        int ifrag = 0;
        char inputkey[MAX_KEY_LENGTH], outKey[MAX_KEY_LENGTH];
        size_t outKeySize = 0;

        /* モノマーデータ(タイプ)を取り出す */
        ofmo_get_monomer_data_master(data_ids);

        /* まずはじめに、すべてのワーカーに一通り、ジョブを投げる */
        for (ifrag = 0; ifrag < config->nfrag; ifrag += config->group_size) {

            ofmo_init_dens_inkey_make(inputkey, MAX_KEY_LENGTH, ifrag, data_ids);
            rc = falanx_task_submit(context,
                                    OFMO_CMD_INIT_DENS,
                                    inputkey,
                                    strlen(inputkey)+1,
                                    0, NULL);
            assert(rc == FLX_SUCCESS);
        }

        /* ワーカーでの密度行列計算が終了するのを待つ */
        while (1) {
            rc = ofmo_falanx_wait_response(context,
                                           outKey, MAX_KEY_LENGTH,
                                           &outKeySize,
                                           NULL);
            if (rc == FLX_NO_WAIT_TASK) {
                break;
            } else if (rc != FLX_SUCCESS) {
                dbg("ofmo_falanx_wait_response failed:%d\n", rc);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }

        /* モノマー間距離の計算を依頼する */
        rc = falanx_task_submit(context,
                                OFMO_CMD_DISTANCE,
                                "", 1,
                                0, NULL);
        assert(rc == FLX_SUCCESS);
        printf("calc. inter-frag distance list \n");
        fflush(stdout);

        rc = ofmo_falanx_wait_response(context,
                                       outKey, MAX_KEY_LENGTH, &outKeySize,
                                       NULL);
        assert(rc == FLX_SUCCESS);

        /* モノマーデータの交換（更新）を行う */
        ofmo_update_monomer_data_master();
    }

    ofmo_master_timestamp("Finish to Calculate Initial Density");
}


/*
 * SCF計算のoutput key から結果を取り出す.
 */
static int
ofmo_falanx_get_scf_output(
    falanx_context_t* context,
    char*   key,
    size_t  size,
    double* outEnergy,
    double* outEnergy0,
    double* outDdv,
    int*    outMonList)
{
    int rc;
    data_store_t* ds = NULL;
    ofmo_scf_output_en_t out;

    ds = falanx_context_get_data_store(context);
    assert(ds != NULL);

    rc = data_store_get_value(ds, key, size, (char*)&out, sizeof(out));
    if (rc != FLX_SUCCESS) {
        dbg("data_store_get_value failed: %d\n", rc);
        return -1;
    }
    if (outEnergy)  { *outEnergy  = out.energy;  }
    if (outEnergy0) { *outEnergy0 = out.energy0; }
    if (outDdv)     { *outDdv     = out.ddv;     }
    if (outMonList) {
        ofmo_scf_outkey_parse(key, size, NULL, NULL, NULL, outMonList);
    }
    return 0;
}


/* --------------------------------------------------
 * モノマーSCC計算
 * --------------------------------------------------
 */
static void
ofmo_master_monomer_scc(
    ofmo_master_config_t* config,
    falanx_context_t*     context,
    double*               outTotalEnergy0)
{
    int    rc;
    int    i;
    int    itera, maxscc, sccconv, scfconv_now;
    double o_energy, o_energy0;
    double total_energy, total_energy0;
    double total_energy0_prev = 0.e0, dE = HUGE_VAL;
    double scc_tol;
    char inputkey[MAX_KEY_LENGTH];
    char outputkey[MAX_KEY_LENGTH];
    size_t outputkeySize = 0;
    int data_ids[6];

    ofmo_data_get_vals("maxscc sccconv", &maxscc, &sccconv);
    scc_tol = pow(0.1e0, (double) sccconv);

    if (!config->need_init_dens) {
        /* ファイルから読み込んだ場合, total_energy0_prev を計算する */
        data_store_t* ds = NULL;
        double* e0s;
        ds = falanx_context_get_data_store(context);
        assert(ds != NULL);
        e0s = (double*) malloc(sizeof(double) * config->nfrag);
        assert(e0s != NULL);
        ofmo_master_get_energy_all(ds, OFMO_ENERGY0, e0s, config->nfrag);
        total_energy0_prev = 0.e0;
        for (i = 0; i < config->nfrag; i++) {
            total_energy0_prev += e0s[i];
        }
        free(e0s);
    }

    printf("---- start SCC (maxscc=%d, sccconv=%d) ----\n", maxscc, sccconv);
    fflush(stdout);

    for (itera = 1; itera <= maxscc; itera++) {
        double step_time;
        step_time = MPI_Wtime();

        scfconv_now = ofmo_get_monomer_scfconv(dE);
        total_energy = total_energy0 = 0.e0;

        /* モノマーデータ(タイプ)を取り出す */
        ofmo_get_monomer_data_master(data_ids);

        /* すべてのモノマー計算を投げる */
        for (i = 0; i < config->nfrag; i++) {
            ofmo_scf_inkey_make(inputkey, MAX_KEY_LENGTH,
                                OFMO_RHF,      /* I_METHO */
                                1,             /* NMON */
                                itera,         /* SCC */
                                scfconv_now,   /* CONV */
                                frag_order[i], -1, -1, /* MON1, (MON2, MON3) */
                                data_ids       /* (new) 現在のモノマーデータ */
                                );

            rc = falanx_task_submit(context,
                                    OFMO_CMD_SCF,
                                    inputkey, strlen(inputkey)+1,
                                    0, NULL);
            assert(rc == FLX_SUCCESS);
        }

        /* タスクの終了を待つ */
        while (1) {        
            rc = ofmo_falanx_wait_response(context,
                                           outputkey, MAX_KEY_LENGTH,
                                           &outputkeySize,
                                           NULL);
            if (rc == FLX_SUCCESS) {
                rc = ofmo_falanx_get_scf_output(context,
                                                outputkey, outputkeySize,
                                                &o_energy, &o_energy0,
                                                IGNORE_DDV, IGNORE_MONOMER_LIST);
                assert(rc == 0);

                total_energy  += o_energy;
                total_energy0 += o_energy0;
            } else if (rc == FLX_ERR_PROC_FAILED) {
                dbg("re-submit task: %s\n", outputkey);
                rc = falanx_task_submit(context,
                                        OFMO_CMD_SCF,
                                        outputkey, outputkeySize,
                                        0, NULL);
                assert(rc == FLX_SUCCESS);
            } else if (rc == FLX_NO_WAIT_TASK) {
                break;
            } else {
                dbg("failed:%d\n", rc);
                MPI_Abort(MPI_COMM_WORLD, 10003);
            }
        }

        step_time = MPI_Wtime() - step_time;
        printf("##SCC itera=%2d  conv=%d ", itera, scfconv_now);
        printf(" energy = %16.8f  energy0 = %16.8f  ( %10.6f )\n",
               total_energy, total_energy0, step_time);
        fflush(stdout);

        /* モノマーデータの交換（更新）を行う */
        ofmo_update_monomer_data_master();

        dE = total_energy0 - total_energy0_prev;
        if (fabs(dE) < scc_tol) {
            printf("==== SCC converged ====\n");
            break;
        }
        total_energy0_prev = total_energy0;
    }

    ofmo_master_timestamp("---- Finish SCC procedure ----");

    *outTotalEnergy0 = total_energy0;

}

static void
ofmo_master_write_monomer_data(
    ofmo_master_config_t* config,
    falanx_context_t*     context)
{
    char* p = getenv("OFMO_FILE_NAME");
    if (p != NULL) {
        ofmo_get_and_write(context, p, config->nfrag);
    }
}

/*
 * dimer タスクの投入
 * @return
 */
static int
ofmo_master_submit_dimer_scf(
    ofmo_master_config_t* config,
    falanx_context_t*     context)
{
    int rc, i, j;
    int maxscc, scfconv;
    int ifrag, jfrag;
    int nscf_dimer = 0;
    int data_ids[6];
    char inputkey[MAX_KEY_LENGTH];
    double ldim;
    double* dist;
    data_store_t* ds = NULL;

    ofmo_data_get_vals("maxscc ldim scfconv", &maxscc, &ldim, &scfconv);

    dist = (double*) malloc(sizeof(double) * config->nfrag);
    assert(dist != NULL);

    ds = falanx_context_get_data_store(context);
    assert(ds != NULL);

    ofmo_get_monomer_data_master(data_ids);

    for (i = 1; i < config->nfrag; i++) {
        ifrag = frag_order[i];
        ofmo_master_get(ds, OFMO_DISTA, ifrag, dist);
        for (j = 0; j < i; j++) {
            jfrag = frag_order[j];
            if (dist[jfrag] < ldim) {
                ofmo_scf_inkey_make(inputkey, MAX_KEY_LENGTH,
                                    OFMO_RHF,
                                    2,
                                    maxscc+100,
                                    scfconv,
                                    ifrag, jfrag, -1,
                                    data_ids
                                    );

                rc = falanx_task_submit(context,
                                        OFMO_CMD_SCF,
                                        inputkey, strlen(inputkey)+1,
                                        0, NULL);
                assert(rc == FLX_SUCCESS);
                nscf_dimer++;
            }
        }
    }

    free(dist);

    return nscf_dimer;
}

static inline void
ofmo_master_set_job(int* joblist, int i, int ifrag, int jfrag)
{
    if (ifrag > jfrag) {
        joblist[i]   = ifrag;
        joblist[i+1] = jfrag;
    } else {
        joblist[i]   = jfrag;
        joblist[i+1] = ifrag;
    }
}

/*
 * Approxタスクの投入
 *
 * @return nes_dimer
 */
static int
ofmo_master_submit_dimer_approx(
    ofmo_master_config_t* config,
    falanx_context_t*     context,
    int*                  outNumEsdimer)
{
    int rc;
    int ii, i, j, k, ie;
    int njob, seq, nes_dimer;
    int * ifrags, jfrag;
    int data_ids[6];
    int joblist[MAXNJOB*2];
    char inputkey[MAX_KEY_LENGTH];
    double ldim;
    double** dists = NULL;
    data_store_t* ds = NULL;

    ds = falanx_context_get_data_store(context);
    assert(ds != NULL);

    ofmo_data_get_vals("ldim", &ldim);

    dists  = ofmo_alloc_dmatrix(NLINE, config->nfrag);
    ifrags = (int*) malloc(sizeof(int) * NLINE);

    ofmo_get_monomer_data_master(data_ids);

    seq = njob = nes_dimer = 0;
    for (ii = 0; ii < config->nfrag; ii += NLINE) {
        ie = ii + NLINE;
        if (ie > config->nfrag) {
            ie = config->nfrag;
        }
        /* フラグメント番号と距離の保存 */
        for (i = ii, k = 0; i < ie; i++, k++) {
            ifrags[k] = frag_order[i];
            ofmo_master_get(ds, OFMO_DISTA, ifrags[k], dists[k]);
        }
        for (j = 0; j < ie; j++) {
            jfrag = frag_order[j];
            for (i = ii, k = 0; i < ie; i++, k++) {
                if (j >= i) {
                    continue;
                }
                if (dists[k][jfrag] >= ldim) {
                    ofmo_master_set_job(joblist, njob*2, ifrags[k], jfrag);
                    njob++;
                    nes_dimer++;
                    if (njob == MAXNJOB) {
                        /* ジョブが溜まったのでタスクを投入する */
                        ofmo_approx_inkey_make(inputkey, MAX_KEY_LENGTH,
                                               data_ids[0], njob, seq);
                        data_store_set(ds, inputkey, strlen(inputkey)+1,
                                      (char*)joblist, sizeof(int)*MAXNJOB*2);
                        rc = falanx_task_submit(context,
                                                OFMO_CMD_APPROX,
                                                inputkey,
                                                strlen(inputkey)+1,
                                                0, NULL);
                        assert(rc == FLX_SUCCESS);
                        njob = 0;
                        seq++;
                    }
                }
            }
        }
    }

    /* 未送信ジョブを投入する */
    if (njob != 0) {
        ofmo_approx_inkey_make(inputkey, MAX_KEY_LENGTH,
                               data_ids[0], njob, seq);
        data_store_set(ds, inputkey, strlen(inputkey)+1,
                      (char*)joblist, sizeof(int)*MAXNJOB*2);
        rc = falanx_task_submit(context,
                                OFMO_CMD_APPROX,
                                inputkey,
                                strlen(inputkey)+1,
                                0, NULL);

        assert(rc == FLX_SUCCESS);
        seq++;
    }

    ofmo_free_dmatrix(dists);
    free(ifrags);

    *outNumEsdimer = nes_dimer;

    return seq;
}

/*
 * outkeyに該当するタスク名を返す.
 * @return タスク名. 該当するタスクがない場合はNULLを返す.
 */
static const char*
ofmo_select_task_name(char* outkey, size_t size)
{
    if (ofmo_is_dimer_scf_input_key(outkey, size)) {
        return OFMO_CMD_SCF;
    } else if (ofmo_is_approx_input_key(outkey, size)) {
        return OFMO_CMD_APPROX;
    } else {
        return NULL;
    }
}

/*
 * dimer/esdimer(approxy)タスクの結果を待つ.
 */
static int
ofmo_master_wait_dimer_esdimer(
    ofmo_master_config_t* config,
    falanx_context_t*     context,
    int                   nscf_dimer,
    int                   nesdimer_tasks,
    double*               outDe0scf,
    double*               outTotalEsdimer0)
{
    int rc;
    double  de0scf;
    double* menergy0;
    char    outputkey[MAX_KEY_LENGTH];
    size_t  outputkeySize;
    data_store_t* ds = NULL;
    double o_energy, o_energy0, o_ddv;
    double total_esdimer0 = 0.0, eesdimer0 = 0.0;
    int monomer_list[3];
    int fnao, fnat;
    int cnt_nscf_dimer = nscf_dimer;


    menergy0 = (double*) malloc(sizeof(double) * config->nfrag);
    assert(menergy0 != NULL);

    ds = falanx_context_get_data_store(context);
    assert(ds != NULL);

    ofmo_master_get_energy_all(ds, OFMO_ENERGY0, menergy0, config->nfrag);

    de0scf = 0.e0;

    while (1) {
        rc = ofmo_falanx_wait_response(context,
                                       outputkey, MAX_KEY_LENGTH,
                                       &outputkeySize,
                                       NULL);
        if (rc == FLX_NO_WAIT_TASK) {
            break;
        } else if (rc == FLX_ERR_PROC_FAILED) {
            const char* cmd = ofmo_select_task_name(outputkey, outputkeySize);
            if (cmd == NULL) {
                dbg("Unknown task output returned. %s\n", outputkey);
                MPI_Abort(MPI_COMM_WORLD, 3);
            }
            dbg("re-submit task: %s\n", outputkey);
            rc = falanx_task_submit(context,
                                    cmd,
                                    outputkey,
                                    outputkeySize,
                                    0, NULL);
            assert(rc == FLX_SUCCESS);
            continue;
        } else if (rc != FLX_SUCCESS) {
            dbg("failed:%d\n", rc);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if (ofmo_is_dimer_scf_output_key(outputkey, outputkeySize)) {
            rc = ofmo_falanx_get_scf_output(context, outputkey, outputkeySize,
                                            &o_energy, &o_energy0, &o_ddv, monomer_list);
            assert(rc == 0);
            assert((monomer_list[0] >= 0) && (monomer_list[1] >= 0));

            de0scf += (o_energy0 - menergy0[monomer_list[0]] - menergy0[monomer_list[1]]);
            de0scf += o_ddv;
            ofmo_get_pop_data(context, outputkey, outputkeySize, OFMO_BOTH_POPS,
                              &fnao, config->aopop_frg, config->fsao2tuao,
                              &fnat, config->atpop_frg, config->fatom2tatom);
            ofmo_accumulate(fnao, config->aopop_frg, config->fsao2tuao, config->aopop_total);
            ofmo_accumulate(fnat, config->atpop_frg, config->fatom2tatom, config->atpop_total);

            cnt_nscf_dimer--;
            if (cnt_nscf_dimer == 0){
                ofmo_master_timestamp("---- Finish to Calc. Dimer SCF ----");
                printf("       ( # of scf dimer = %d )\n", nscf_dimer);
                printf("        dE0(SCF dimer) = %16.8f\n", de0scf);
                fflush(stdout);
            }
        } else if (ofmo_is_approx_output_key(outputkey, outputkeySize)) {
            rc = data_store_get_value(ds, outputkey, outputkeySize,
                                      (char*)&eesdimer0, sizeof(double));
            assert(rc == FLX_SUCCESS);

            total_esdimer0 += eesdimer0;

            nesdimer_tasks--;
            if (nesdimer_tasks == 0) {
                ofmo_master_timestamp("---- Finish to Calc. ES Dimer ----");
            }
        } else {
            dbg("Unknown task output returned. %s\n", outputkey);
            MPI_Abort(MPI_COMM_WORLD, 3);
        }
    }

    *outDe0scf = de0scf;
    *outTotalEsdimer0 = total_esdimer0;

    free(menergy0);

    return 0;
}

static void
ofmo_master_dimer_calc(
    ofmo_master_config_t* config,
    falanx_context_t*     context,
    double*               outDe0scf,
    double*               outTotalEsdimer0,
    int*                  outNumEsdimer)
{
    int nscf_dimer, nesdimer_tasks;

    nscf_dimer = ofmo_master_submit_dimer_scf(config, context);
    nesdimer_tasks = ofmo_master_submit_dimer_approx(config,
                                                    context,
                                             outNumEsdimer);
    ofmo_master_wait_dimer_esdimer(config,
                                  context,
                               nscf_dimer,
                           nesdimer_tasks,
                                outDe0scf,
                        outTotalEsdimer0);
}


static void
ofmo_master_print_energy(
    double total_energy0,
    double de0scf,
    double total_esdimer0,
    int    nes_dimer)
{
    double total_fmo_energy = total_energy0 + de0scf + total_esdimer0;
    printf("    (# of ES dimer = %d)\n", nes_dimer);
    printf("    dE0(ES dimer) = %16.8f\n", total_esdimer0);
    printf("    total energy(FMO) = %16.8f\n", total_fmo_energy);
    fflush(stdout);
}

/* ------------------------------------------------------------
 * 全体のAO population、ならびに、atomic populationの計算
 * ------------------------------------------------------------
 */
static void
ofmo_master_total_population(
    ofmo_master_config_t* config,
    falanx_context_t*     context)
{
    int natom;
    int iao;
    int iat;
    int ifrag;
    int* atomic_number, * nfao, * nfatom, ** ifatom;
    int* atm_lcs, * ushel_lqn, * ushel_ini;
    int** msao2tuao;
    double charge;
    data_store_t* ds = NULL;


    ds = falanx_context_get_data_store(context);
    assert(ds != NULL);

    ofmo_data_get_vals("natom msao2tuao atn atm_lcs ushel_lqn ushel_ini nfao nfatom ifatom",
        &natom, &msao2tuao, &atomic_number, &atm_lcs, &ushel_lqn,
        &ushel_ini, &nfao, &nfatom, &ifatom);

    for (ifrag = 0; ifrag < config->nfrag; ifrag++) {
        ofmo_get_monomer_aopop_master(ds, ifrag, config->aopop_frg);

        ofmo_get_monomer_atpop_master(ds, ifrag, config->atpop_frg);

        for (iao = 0; iao < nfao[ifrag]; iao++) {
            config->aopop_total[msao2tuao[ifrag][iao]] += config->aopop_frg[iao];
        }
        for (iat = 0; iat < nfatom[ifrag]; iat++) {
            config->atpop_total[ifatom[ifrag][iat]] += config->atpop_frg[iat];
        }
    }
    printf("---- ATOMIC POPULATION DATA -----\n"
           "   SN AN  ATOMIC POP  NET CHARGE\n"
           "---------------------------------\n");
    charge = 0.e0;
    for (iat = 0; iat < natom; iat++) {
        charge += ((double) atomic_number[iat] - config->atpop_total[iat]);
        printf(" %4d %2d  %10.3f  %10.3f\n",
               (iat + 1), atomic_number[iat], config->atpop_total[iat],
               ((double) atomic_number[iat] - config->atpop_total[iat]));
    }
    printf("---------------------------------\n"
           "total charge = %10.3f\n"
           "---------------------------------\n", charge);

    ofmo_master_timestamp("-- calc AO and atomic population --");
}


static int
ofmo_reconfig_worker(
    MPI_Comm basecomm,
    falanx_context_t* context,
    int cl_size,
    int ds_size,
    int wk_size,
    int* outNumWorkers)
{
    int rank, size, color, n;
    int q_base, ds_base, wk_base;

    assert(basecomm != MPI_COMM_NULL);
    MPI_Comm_rank(basecomm, &rank);
    MPI_Comm_size(basecomm, &size);

    q_base  = cl_size;
    ds_base = q_base + 1;
    wk_base = ds_base + ds_size;

    if (rank < wk_base) {
        color = MPI_UNDEFINED;
    } else {
        color = FLX_CL_WORKER + (rank - wk_base) / wk_size;
    }

    if (FLX_SUCCESS != falanx_split(context, color, rank)) {
        return FLX_ERROR;
    }

    n = (size - wk_base) / wk_size;
    if ((size - wk_base) % wk_size) {
        n++;
    }
    if (outNumWorkers != NULL) {
        *outNumWorkers = n;
    }

    return FLX_SUCCESS;
}

/*
 * --------------------------------------------------
 * マスターのメインプログラム
 * --------------------------------------------------
 */
int
MASTER_main(int argc, char* argv[])
{
    int rc;
    int rank, size;
    int num_workers;
    int color;
    falanx_context_t* context = NULL;
    ofmo_master_config_t config = { 0 };
    int init_provided;
    int nao, natom, maxnfao, maxnfatom;
    /* dimer SCFタスクの結果 */
    double de0scf;
    /* monomer SCF タスクの結果 */
    double total_energy0;
    /* approx タスクの結果 */
    double total_esdimer0;
    int    nes_dimer;


    MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &init_provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    assert(size >= 4); /* Cl, Tq, Ds, Wk */

    if (MPI_THREAD_SERIALIZED != init_provided) {
        if (rank == 0) {
            printf("MPI_THREAD_SERIALIZED is not supported\n");
            MPI_Abort(MPI_COMM_WORLD, 10001);
        }
    }
    TLOG_INITIALIZE();

    /* 初期化 */
    rc = ofmo_init_config(&config, MPI_COMM_WORLD, argc, argv);
    if (rc != 0)
        return rc;

    context = falanx_init_simple(MPI_COMM_WORLD,
                                 1,                 /* #Client      */
                                 config.niogroup,   /* #DS          */
                                 config.group_size, /* Worker size  */
                                 &num_workers);
    assert(context != NULL);

    falanx_register_function(context, OFMO_CMD_INIT_DENS, ofmo_init_dens_task);
    falanx_register_function(context, OFMO_CMD_DISTANCE, ofmo_distance_task);
    falanx_register_function(context, OFMO_CMD_SCF, ofmo_scf_task);
    falanx_register_function(context, OFMO_CMD_APPROX, ofmo_approx_task);


    if (falanx_context_is_worker(context)) {
      TLOG_LOG_IN(2);
        /* temporary */
        /* workerの初期化処理 */
        char prof_name[MAXSTRLEN];
        int maxnfao, maxlqn;
        MPI_Comm wcomm;

        ofmo_data_put_vals("nintic", config.eribfsz);

        color = falanx_context_get_color(context);
        wcomm = falanx_context_get_mpicomm_by_color(context, color);
        assert(wcomm != MPI_COMM_NULL);
        ofmo_data_get_vals("maxnfao maxlqn", &maxnfao, &maxlqn);
        ofmo_create_worker_prof_name(config.input, config.header,
                                     rank, prof_name, MAXSTRLEN);

        ofmo_prof_init(prof_name, wcomm);

#ifdef USE_CUDA
        {
          int ierr;
          int myrank, nprocs;
          int ndev = config.ndev;
          if ( fp_prof ) fprintf( fp_prof,"ndev = %d\n", ndev);
          MPI_Comm_rank(wcomm, &myrank);
          MPI_Comm_size(wcomm, &nprocs);
          ierr = cuda_Init(ndev, myrank, nprocs, wcomm);
          assert(ierr>=0);
        }
#endif
        ofmo_frag_init();
        ofmo_integ_init(maxlqn);
        ofmo_scf_init(config.nbody * maxnfao);
      TLOG_LOG_OUT(2);
    }


    rc = falanx_start(context);
    if (falanx_context_is_client(context)) {
      TLOG_LOG_IN(2);

        /* ***** FMO計算で用いるデータ構造の決定 ***** */
        ofmo_master_resolve_data_struct(&config);

        /* ***** populationデータ作成のための初期化 ***** */
        ofmo_data_get_vals("nao natom maxnfao maxnfatom",
                           &nao, &natom, &maxnfao, &maxnfatom);
        ofmo_init_population_data(&config, nao, natom, maxnfao, maxnfatom);

        /* ****** ワーカーに指示したジョブのリスト ***** */
        ofmo_init_job_list(&config);

        /* ****** 初期密度行列計算 ***** */
        ofmo_master_init_dens(&config, context);
        if (config.nbody == 0) {
            falanx_stop(context);
            goto finish;
        }

      TLOG_LOG_OUT(2);
      TLOG_LOG_IN(3);
        /* ****** モノマーSCC計算 ****** */
        ofmo_master_monomer_scc(&config, context, &total_energy0);

        /* ****** モノマーのデータ（密度行列、AO populationなど）をファイルに出力 */
        ofmo_master_write_monomer_data(&config, context);
      TLOG_LOG_OUT(3);

    }
#if 0
//    falanx_stop(context);
//    goto finish;
#ifdef USE_CUDA
    falanx_stop(context);
    TLOG_LOG_EVENT(1);
    cuda_print_wifc4c();
    falanx_start(context);
#endif
#endif

    if (REQUIRE_RECONFIG(config)) {
        int group_size;
        falanx_stop(context);

        if (rank == 0) { ofmo_master_timestamp("Start Worker Reconfig"); }
        group_size = ofmo_master_resolve_group_size_dimer(&config);
        rc = ofmo_reconfig_worker(MPI_COMM_WORLD,
                                  context,
                                  1,                 /* #Client      */
                                  config.niogroup,   /* #DS          */
                                  group_size,        /* Worker size  */
                                  NULL);
        if (rc != FLX_SUCCESS) {
            dbg("Could'nt reconfig worker size\n");
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
        config.group_size = group_size;

        if (falanx_context_is_worker(context)) {
            char prof_name[MAXSTRLEN];
            int maxnfao, maxlqn;
            MPI_Comm wcomm;

            color = falanx_context_get_color(context);
            wcomm = falanx_context_get_mpicomm_by_color(context, color);
            assert(wcomm != MPI_COMM_NULL);
            ofmo_data_get_vals("maxnfao maxlqn", &maxnfao, &maxlqn);
            ofmo_create_worker_prof_name(config.input, config.header,
                                         rank, prof_name, MAXSTRLEN);

            ofmo_prof_reinit(prof_name, wcomm);

#ifdef USE_CUDA
            {
              int ierr;
              int myrank, nprocs;
              MPI_Comm_rank(wcomm, &myrank);
              MPI_Comm_size(wcomm, &nprocs);
              ierr = cuda_Reconfig(myrank, nprocs, wcomm);
              assert(ierr>=0);
            }
#endif
        }

        if (rank == 0) { ofmo_master_timestamp("Finish Worker Reconfig"); }
        falanx_start(context);
    }

    if (falanx_context_is_client(context)) {
      TLOG_LOG_IN(4);
        /* ****** dimerの計算 ****** */
        if (config.nbody >= 2) {
            ofmo_master_dimer_calc(&config, context, &de0scf, &total_esdimer0, &nes_dimer);
            ofmo_master_print_energy(total_energy0, de0scf, total_esdimer0, nes_dimer);
        }

        /* ***** 全体のAO population、ならびに、atomic populationの計算 */
        if (config.nbody > 0) {
            ofmo_master_total_population(&config, context);
        }

      TLOG_LOG_OUT(4);
        falanx_stop(context);
    }


finish:
#ifdef USE_CUDA
    if (color >= FLX_CL_WORKER) {
      int ierr = cuda_Finalize();
      assert(ierr>=0);
      ofmo_show_cache_prof();
    }
#endif
    falanx_finalize(context);
    ofmo_finalize_config(&config);

    TLOG_LOG_EVENT(1);
    TLOG_FINALIZE();
    MPI_Finalize();

    return 0;
}

