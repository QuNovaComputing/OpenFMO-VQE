/*
 * Skelton program for RHF calculation
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <libgen.h>
#include <assert.h>

#include "ofmo-parallel.h"
#include "ofmo-prof.h"
#include "ofmo-scf.h"
#include "ofmo-twoint.h"
#include "skel-rhf-main.h"
#include "skel-rhf-data.h"
#include "skel-rhf-calc.h"
#include "skel-w2e.h"

#ifdef USE_CUDA
#include "cuda/cuda-drv.h"
#include "cuda/cudalib.h"
#include "cuda/cuda-integ.h"
#endif

#define DEFAULT_BUFFER_SIZE 0
#define DEFAULT_MAX_SCF_CYCLE 30;

// ----------------

static void skel_rhf_show_help(const char *myname)
{
  fprintf(stderr,"Usage: %s [-snvh][-B buffer] input [density]\n", myname);
  fprintf(stderr,"  -B buf: # buffer size (MB, default: %d)\n", DEFAULT_BUFFER_SIZE);
  fprintf(stderr,"  -s: sync\n");
  fprintf(stderr,"  -n: dryrun\n");
  fprintf(stderr,"  -v: verbose\n");
  fprintf(stderr,"  -h: show this help\n");
#ifdef USE_CUDA
  fprintf(stderr," Options for GPGPU:\n");
  fprintf(stderr,"  -d ndev: # devices (default:1)\n");
//  fprintf(stderr,"  -b nblk: # blocks\n");
//  fprintf(stderr,"  -t nthb: # threads\n");
//  fprintf(stderr,"  -a type: int. algorithm (0:OS_X, 1:OS, 2:RYS_X, 3:RYS)\n");
#endif
}

// ----------------

static void skel_rhf_show_config(skel_rhf_config_t* config)
{
  printf("-------- parallel information --------\n");
  printf(" # of process = %d\n", config->nprocs );
  printf(" # of threads = %d\n", config->nthreads );
  printf(" Buffer size  = %ld (MB)\n", config->eribfsz );
#ifdef USE_CUDA
  printf("-------- GPGPU information -----------\n");
  printf(" # of devices = %d\n", config->ndev );
  cuda_Print_DEFS();
#endif
  printf("-------- input information -----------\n");
  printf(" Input   = %s\n", config->input );
  if ( config->density[0] != '\0' ) printf(" Density = %s\n", config->density );
  printf("--------------------------------------\n");
  fflush(stdout);
}

// ----------------

static int skel_rhf_init_config(
    skel_rhf_config_t* config,
    MPI_Comm comm,
    int argc,
    char* argv[])
{
  int c;
  char myname[MAXSTRLEN];
  char *endptr;
  extern char *optarg;
  extern int optind, opterr;
  int optsync = false;
  int help = 0;
  int dryrun = false;
  int verbose = false;
  int nprocs, myrank, maxthreads;
  long buffer_size=DEFAULT_BUFFER_SIZE;
  int maxscfcyc = DEFAULT_MAX_SCF_CYCLE;
  int ndev = 0;
#ifdef USE_CUDA
  int type =INTTYPE_QUERY;
  int nblk = -1;
  int nthb = -1;
#endif

  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
  maxthreads = omp_get_max_threads();

  strncpy(myname, basename(argv[0]), MAXSTRLEN);

//#ifdef USE_CUDA
  while ((c=getopt(argc, argv, "B:C:d:b:t:c:a:snvh"))!=-1) {
//#else
//  while ((c=getopt(argc, argv, "B:C:snvh"))!=-1) {
//#endif
    switch(c) {
      case 'B':
        buffer_size = strtol(optarg, &endptr, 10);
        if (endptr == optarg) buffer_size=DEFAULT_BUFFER_SIZE;
        if (buffer_size<0) help=-1;
        break;
      case 'C':
        maxscfcyc = strtol(optarg, &endptr, 10);
        if (endptr == optarg) maxscfcyc = DEFAULT_MAX_SCF_CYCLE;
        if (maxscfcyc<=0) help=-1;
        break;
      case 'd':
        ndev = strtol(optarg, &endptr, 10);
        if (endptr == optarg) ndev=-1;
        if (ndev<0) help=-1;
        break;
#ifdef USE_CUDA
      case 'b':
        nblk = strtol(optarg, &endptr, 10);
        if (endptr == optarg) nblk=-1;
        if (nblk<=0) help=-1;
        break;
      case 't':
        nthb = strtol(optarg, &endptr, 10);
        if (endptr == optarg) nthb=-1;
        if (nthb<=0) help=-1;
        break;
      case 'a':
        type = strtol(optarg, &endptr, 10);
        if (type<0 || type>3) help=-1;
        break;
#endif
      case 's':
        optsync=true;
        break;
      case 'n':
        dryrun=true;
        break;
      case 'v':
        verbose=true;
        break;
      case '?':
        help=-1;
        break;
      case 'h':
      default:
        help=1;
        break;
    }
  }
  argc -= optind;
  argv += optind;

  if (help!=0 || argc < 1) {
    if (myrank==0) { skel_rhf_show_help(myname); };
    MPI_Finalize();
    return 2;
  }

  config->comm = comm;
  config->nprocs = nprocs;
  config->myrank = myrank;
  config->nthreads = maxthreads;
  config->eribfsz = buffer_size;
  config->maxscfcyc = maxscfcyc;
  config->verbose = verbose;
  config->dryrun = dryrun;
  config->optsync = optsync;
#ifdef USE_CUDA
  config->ndev = ndev;
  config->nblk = nblk;
  config->nthb = nthb;
  config->type = type;
#endif
  config->scfe = 1.e-8;
  config->scfd = 1.e-5;
  config->eps_ps4 = 1.e-20F;
  config->eps_eri = 1.e-15F;
  config->eps_sch = 1.e-12F;

  strncpy(config->input, argv[0], MAXSTRLEN);
  config->density[0] = '\0';
  if ( argc > 1 ) { strncpy(config->density, argv[1], MAXSTRLEN); };

  if (myrank==0) { skel_rhf_show_config(config); };

  return 0;
}

// ----------------

#ifdef FJ_MAIN
int MAIN__( int argc, char* argv[] ) {
#else
int main( int argc, char* argv[] ) {
#endif
  int rc;
  int maxthreads;
  int nprocs, myrank, provided;
  skel_rhf_config_t config = { 0 };
  skel_rhf_data_t data = { 0 };
  skel_rhf_array_t array = { 0 };
  MPI_Comm wcomm = MPI_COMM_WORLD;
  int master = false;
  double wt0, wt1;
  double t0, t1;
  double energy;

  maxthreads = omp_get_max_threads();

#ifdef _OPENMP
  //MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &provided );
  MPI_Init_thread( &argc, &argv, MPI_THREAD_SINGLE, &provided );
#else
  MPI_Init( &argc, &argv);
#endif
  MPI_Comm_size( wcomm, &nprocs );
  MPI_Comm_rank( wcomm, &myrank );
  master = (myrank == 0);
  wt0 = MPI_Wtime();

#ifdef _OPENMP
  //if ( provided != MPI_THREAD_FUNNELED ) {
  if ( provided != MPI_THREAD_SINGLE ) {
    if ( master ) {
      //printf("This program requires MPI_THREAD_FUNNELED\n");
      printf("This program requires MPI_THREAD_SINGLE\n");
      MPI_Abort( wcomm, 1 );
    }
  }
#endif

  // parse args
  rc = skel_rhf_init_config(&config, wcomm, argc, argv);
  if (rc!=0) return rc;
  setup_w2e(wcomm, config.optsync);

  // setup profile
  ofmo_prof_init(NULL, wcomm);

  // read input file
  rc = skel_rhf_init_data(&config, &data);
  assert(rc>=0);

  // init cuda device
#ifdef USE_CUDA
  t0 = MPI_Wtime();
  {
    int ndev = config.ndev;
    if ( fp_prof ) {fprintf( fp_prof,"ndev = %d\n", ndev);}
    rc = cuda_Init(ndev, myrank, nprocs, wcomm);
    assert(rc>=0);
    cuda_set_optsync(config.optsync);
  }
  t1 = MPI_Wtime();
  if ( master ) { printf( "etime (%s) = %9.6f\n", "CUDA Init", t1-t0 ); }
#endif

  ofmo_scf_init(data.nao);

  // make cutoff table
  rc = skel_rhf_cutoff_make_table(&config, &data);

  // allocate arrays
  rc = skel_rhf_init_array(&config, &data, &array);

  // 1-e integrals
  rc = skel_rhf_oneint(&config, &data, &array);

  // construct initial density
  rc = skel_rhf_init_density(&config, &data, &array);

  // buffer 2-e integrals
  t0 = MPI_Wtime();
  rc = skel_rhf_twoint_first(&config, &data, &array);
  t1 = MPI_Wtime();
  if ( master && config.eribfsz>0 ) { printf( "etime (%s) = %9.6f\n", "buffer int2e", t1-t0 ); }

  // scf calc.
  t0 = MPI_Wtime();
  energy = skel_rhf_scf(&config, &data, &array);
  t1 = MPI_Wtime();
  if ( master ) { printf( "etime (%s) = %9.6f\n", "SCF", t1-t0 ); }

  fflush(stdout);
  if ( master ) {
    printf( "Energy = %17.10f\n", energy );
  }

  skel_rhf_fin_array(&array);
  skel_rhf_fin_data(&data);

  MPI_Finalize();

#ifdef USE_CUDA
  if ( myrank == 0 ) {
    int ndev = config.ndev;
    int nblocks, nthreads;
    nblocks = cuda_get_numBlocks();
    nthreads = cuda_get_numThreads();
    printf("NDEV (NBLKS, NTHB): %d (%3d,%3d)\n", ndev, nblocks, nthreads);
  }
  if (myrank==0) print_w2e();
  rc = cuda_Finalize();
  if (rc<0) exit(rc);
#else
  if (myrank==0) print_w2e();
#endif

  wt1 = MPI_Wtime();
  if ( master ) { printf( "etime = %9.6f\n", wt1-wt0 ); }

  return 0;
}

// ----------------
