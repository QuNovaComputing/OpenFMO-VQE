/**
 * @file ofmo-twoint.c
 * @brief A file that describes a group of functions that perform miscellaneous
 * processing other than integral calculations related to 2-electron integration.
 *
 * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ofmo-def.h"

//#define Free(p) if ((p)!=NULL) free( (p) )

#ifdef _OPENMP
#include <omp.h>
#else
static int omp_get_thread_num() { return 0; }
static int omp_get_max_threads() { return 1; }
#endif

// Memory related to the buffer method
static int	*LAST_IJCS         = NULL;
static int	*LAST_KLCS         = NULL;
static int	*LAST_ERI_TYPE     = NULL;
static size_t	*EBUF_NON_ZERO_ERI = NULL;

static double	**EBUF_VAL         = NULL;
static short int **EBUF_IND4       = NULL;
static size_t	*EBUF_MAX_NZERI    = NULL;

// Integral temporary storage area used in the direct method
#define TMP_INTEG_SIZE	64	// MB
static size_t	*ETMP_MAX_NZERI    = NULL;
static double	**ETMP_VAL	   = NULL;
static short    **ETMP_IND4	   = NULL;
static size_t	*ETMP_NZERI        = NULL;

// Memory related to Fock matrix calculation when threads are parallel
static int	*MAX_NAO       = NULL;
static double	**G_THREADS    = NULL;

static double	*DENS_SQUARE   = NULL;
static int	BUF_NAO        = 0;

static float    *DENS_CS   = NULL;
static int      BUF_NCS        = 0;
static int      TWOINT_NCS = -1;

static int	NTHREADS       = 1;

double _twoint_inv2_;
double _twoint_inv3_;
double _twoint_spi2_;

static void ofmo_twoint_finalize() {
    int i;
    Free( LAST_IJCS );
    Free( LAST_KLCS );
    Free( LAST_ERI_TYPE );
    Free( EBUF_NON_ZERO_ERI );
    Free( DENS_SQUARE );
    Free( EBUF_MAX_NZERI );
    Free( MAX_NAO );
    if ( EBUF_VAL != NULL ) {
	for ( i=0; i<NTHREADS; i++ ) Free( EBUF_VAL[i] );
	Free( EBUF_VAL );
    }
    if ( EBUF_IND4 != NULL ) {
	for ( i=0; i<NTHREADS; i++ ) Free( EBUF_IND4[i] );
	Free( EBUF_IND4 );
    }
    if ( G_THREADS != NULL ) {
	for ( i=0; i<NTHREADS; i++ ) Free( G_THREADS[i] );
	Free( G_THREADS );
    }
    MAX_NAO        = NULL;
    //
    Free( ETMP_MAX_NZERI );
    Free( ETMP_NZERI );
    if ( ETMP_VAL != NULL ) {
	for ( i=0; i<NTHREADS; i++ ) Free( ETMP_VAL[i] );
	Free( ETMP_VAL );
    }
    if ( ETMP_IND4 != NULL ) {
	for ( i=0; i<NTHREADS; i++ ) Free( ETMP_IND4[i] );
	Free( ETMP_IND4 );
    }
}

/** Initialization function of 2-electron integral function
 * @ingroup integ-misc
 *
 * buffered direct SCF法をMPI/OpenMP hybrid並列で実行するために必要な
 * 各種変数領域の確保と、その初期化を行う
 *
 * @attention スレッド並列ではない部分から呼び出す
 *
 * */
int ofmo_twoint_init() {
    static int called = false;
    int nthreads, i;

    if ( called ) return 0;
    nthreads = omp_get_max_threads();
    LAST_IJCS        = (int*)malloc( sizeof(int) * nthreads );
    LAST_KLCS        = (int*)malloc( sizeof(int) * nthreads );
    LAST_ERI_TYPE    = (int*)malloc( sizeof(int) * nthreads );
    EBUF_NON_ZERO_ERI= (size_t*)malloc( sizeof(size_t) * nthreads );
    EBUF_VAL         = (double**)malloc( sizeof(double*) * nthreads );
    EBUF_IND4        = (short int**)malloc( sizeof(short int*) * nthreads );
    EBUF_MAX_NZERI   = (size_t*)malloc( sizeof(size_t) * nthreads );
    G_THREADS        = (double**)malloc( sizeof(double*) * nthreads );
    MAX_NAO          = (int*)malloc( sizeof(int) * nthreads );
    // for direct method
    ETMP_MAX_NZERI = (size_t*)malloc( sizeof(size_t) * nthreads );
    ETMP_NZERI     = (size_t*)malloc( sizeof(size_t) * nthreads );
    ETMP_VAL       = (double**)malloc( sizeof(double*) * nthreads );
    ETMP_IND4      = (short**)malloc( sizeof(short*) * nthreads );
    for ( i=0; i<nthreads; i++ ) {
	ETMP_VAL[i] = NULL;
	ETMP_IND4[i] = NULL;
    }
#pragma omp parallel
    {
	int mythread;
	size_t max_nzeri;
	mythread = omp_get_thread_num();
	max_nzeri = TMP_INTEG_SIZE*1024*1024/
	    ( sizeof(double)*1 + sizeof(short)*4 );
	ETMP_VAL[mythread] = (double*)malloc(sizeof(double)*max_nzeri );
	ETMP_IND4[mythread] = (short*)malloc( sizeof(short)*max_nzeri*4 );
	ETMP_NZERI[mythread] = 0;
	ETMP_MAX_NZERI[mythread] = max_nzeri;
    }

    for ( i=0; i<nthreads; i++ )
	LAST_IJCS[i] = LAST_KLCS[i] = -1;
    //for ( i=0; i<nthreads; i++ ) LAST_ERI_TYPE[i] = 1000000;
    for ( i=0; i<nthreads; i++ ) LAST_ERI_TYPE[i] = -1;
    for ( i=0; i<nthreads; i++ )
	EBUF_NON_ZERO_ERI[i] = EBUF_MAX_NZERI[i] = MAX_NAO[i] = 0;
    for ( i=0; i<nthreads; i++ ) {
	EBUF_VAL[i]  = NULL;
	EBUF_IND4[i] = NULL;
	G_THREADS[i] = NULL;
    }

    NTHREADS = nthreads;
    // added
    double pi2;
    pi2 = 2.e0 * atan(1.e0);
    _twoint_inv2_ = 1.e0 / 2.e0;
    _twoint_inv3_ = 1.e0 / 3.e0;
    _twoint_spi2_ = sqrt(pi2);

    atexit( ofmo_twoint_finalize );
    called = true;
    return 0;
}

/** Buffer size setting
 * 
 * buffered direct SCFで計算した２電子積分を格納するための領域（バッファ）
 * を確保する関数。
 * MPI/OpenMP hybrid並列での利用を考慮してある。
 *
 * @attention hybrid並列時にはスレッド並列部分から呼び出す
 *
 * @param[in] mythread スレッド番号
 * @param[in] ebuf_buffer_size_mb MB単位でのバッファサイズ
 *
 * @ingroup integ-misc
 *
 * */
size_t ofmo_twoint_set_buffer_size( const int mythread,
	const size_t ebuf_buffer_size_mb ) {
    size_t ebuf_max_nzeri;
    ebuf_max_nzeri = (size_t) (ebuf_buffer_size_mb * 1024 * 1024 /
	(double)( sizeof(double)+4*sizeof(short int) ) );

    if ( ebuf_max_nzeri > EBUF_MAX_NZERI[mythread] ) {
	if ( EBUF_VAL[mythread] != NULL ) Free( EBUF_VAL[mythread] );
	if ( EBUF_IND4[mythread] != NULL ) Free( EBUF_IND4[mythread] );
	EBUF_VAL[mythread] =
	    (double*)malloc( sizeof(double)*ebuf_max_nzeri );
	EBUF_IND4[mythread] =
	    (short int*)malloc( sizeof(short int)*ebuf_max_nzeri*4 );
	EBUF_MAX_NZERI[mythread] = ebuf_max_nzeri;
    }
    return EBUF_MAX_NZERI[mythread];
}

size_t ofmo_twoint_get_max_nzeri( const int mythread ) {
    return EBUF_MAX_NZERI[mythread];
}
/** Function that secures the temporary area used in G matrix calculation
 *
 * スレッド並列（hybrid並列を含む）でG行列（Fock行列の２電子積分部分）を
 * 生成する場合には、各スレッドは、独自の領域に部分G行列を計算した後、
 * リダクション処理を行う、というステップを踏む。
 * スレッド並列化した場合には、すべてのスレッドが共通のG行列に更新を
 * かける、という手段も使えるが、そうすると排他制御にコストが
 * かかってしまい、性能向上が得られない。
 * したがって、G行列作成時にスレッド毎に部分G行列を作成して、
 * 後でまとめる、という手法をとっている。
 * 
 * @attention スレッド並列部分で呼び出す必要がある
 *
 * @param[in] mythread スレッド番号
 * @param[in] maxnao 確保したいG行列のAO数。考えられるAO数の最大値を
 * 最初に与えておくと、無駄なfree, mallocの呼び出しが減らせる。
 *
 * @ingroup integ-misc
 *
 * */
double* ofmo_twoint_alloc_local_gmat( const int mythread,
	const int maxnao ) {
    int nnao;
    if ( maxnao > MAX_NAO[mythread] ) {
	nnao = maxnao * maxnao;
	if ( G_THREADS[mythread] != NULL ) Free( G_THREADS[mythread] );
	G_THREADS[mythread] = (double*)malloc( sizeof(double) * nnao );
	MAX_NAO[mythread] = maxnao;
    }
    return G_THREADS[mythread];
}

/** Allocate a buffer to store the density matrix as a square matrix
 *
 * G行列作成時には対称行列である密度行列を参照するが、そのとき、圧縮形式
 * のまま参照するよりも、正方行列として参照したほうがインデックスなどの
 * 計算量が少なくて済む。そのため、圧縮形式の密度行列を正方行列に
 * 展開するための領域が必要となる。
 * この関数では、そのための領域を確保して、確保した並列のポインタを
 * 返す。
 *
 * @attention スレッド並列領域外から呼ぶ
 *
 * @param[in] maxnao 確保したい密度行列のAO数。考えられるAO数の最大値を
 * 最初に与えておくと、無駄なfree, mallocの呼び出しが減らせる。
 *
 * @ingroup integ-misc
 *
 * */
double* ofmo_twoint_alloc_square_density( const int maxnao ) {
    if ( maxnao > BUF_NAO ) {
	if ( DENS_SQUARE != NULL ) Free( DENS_SQUARE );
	DENS_SQUARE = (double*)malloc( sizeof(double) * maxnao * maxnao );
	BUF_NAO = maxnao;
    }
    return DENS_SQUARE;
}

/** A function that registers the CS pair number when the buffer is full
 *
 * buffered direct SCF法では、バッファが一杯になった時点で計算していた
 * ２電子積分の２つのCSペア番号を記録する必要がある。
 * この関数では、そのうち、外側のCSペア番号を登録する。
 * スレッド並列を考慮してある。
 *
 * @param[in] mythread スレッド番号
 * @param[in] last_ijcs バッファが一杯になった時点での、外側CSペア番号
 *
 * @ingroup integ-misc
 *
 * */
void ofmo_twoint_set_last_ijcs( const int mythread, const int last_ijcs ) {
    LAST_IJCS[mythread] = last_ijcs;
}

/** A function that registers the CS pair number when the buffer is full
 *
 * buffered direct SCF法では、バッファが一杯になった時点で計算していた
 * ２電子積分の２つのCSペア番号を記録する必要がある。
 * この関数では、そのうち、内側のCSペア番号を登録する。
 * スレッド並列を考慮してある。
 *
 * @param[in] mythread スレッド番号
 * @param[in] last_ijcs バッファが一杯になった時点での、内側CSペア番号
 *
 * @ingroup integ-misc
 *
 * */
void ofmo_twoint_set_last_klcs( const int mythread, const int last_klcs ) {
    LAST_KLCS[mythread] = last_klcs;
}

/** A function that registers the integral type that was calculated when the buffer was full
 *
 * buffered direct SCF法では、バッファが一杯になった時点で計算していた
 * ２電子積分の積分タイプを記録する必要がある。
 * この関数では、その登録処理を行う。
 *
 * @param[in] mythread スレッド番号
 * @param[in] last_eri_type バッファが一杯になった時点で計算していた
 * ２電子積分のタイプ
 *
 * @ingroup integ-misc
 *
 * */
void ofmo_twoint_set_last_eri_type( const int mythread,
	const int last_eri_type ) {
    LAST_ERI_TYPE[mythread] = last_eri_type;
}

/** A function that registers the number of contracted two-electron integrals stored in the buffer
 *
 * buffered direct SCF法では、繰り返し計算のたびに行うG行列計算時に、
 * バッファに保存している２電子積分と、再計算する２電子積分を用いる。
 * その際に必要となるバッファに保存している２電子積分の個数を登録する
 * 関数である。
 *
 * @param[in] mythread スレッド番号
 * @param[in] ebuf_non_zero_eri
 *
 * @ingroup integ-misc
 * */
void ofmo_twoint_set_stored_nzeri( const int mythread,
	const size_t ebuf_non_zero_eri ) {
    EBUF_NON_ZERO_ERI[mythread] = ebuf_non_zero_eri;
}

/** A function that obtains the integral type that was calculated when the buffer was full
 *
 * In the buffered direct SCF method, information on the integral type of
 * the two-electron integral calculated when the buffer is full is
 * required at the start of the calculation of the two-electron integral not stored in the buffer.
 * This function is a function that returns its integral type.
 *
 * @param[in] mythread Thread number
 * @return Integral type of two-electron integral that was calculated when the buffer was full
 *
 * @ingroup integ-misc
 *
 * */
int ofmo_twoint_get_last_eri_type( const int mythread ) {
    return LAST_ERI_TYPE[mythread];
}

/** A function that gets the CS pair number when the buffer is full
 *
 * buffered direct SCF法では、バッファに保存されていない
 * ２電子積分の計算開始時に、バッファが一杯になった時点で計算していた
 * ２電子積分の２つのCSペア番号の情報が必要となる。
 * この関数では、そのうち、外側のCSペア番号を得る関数である。
 * スレッド並列を考慮してある。
 *
 * @param[in] mythread スレッド番号
 * @return バッファが一杯になった時点での外側CSペア番号
 *
 * @ingroup integ-misc
 *
 * */
int ofmo_twoint_get_last_ijcs( const int mythread ) {
    return LAST_IJCS[mythread];
}

/** バッファが一杯になって時点でのCSペア番号を得る関数
 *
 * buffered direct SCF法では、バッファに保存されていない
 * ２電子積分の計算開始時に、バッファが一杯になった時点で計算していた
 * ２電子積分の２つのCSペア番号の情報が必要となる。
 * この関数では、そのうち、内側のCSペア番号を得る関数である。
 * スレッド並列を考慮してある。
 *
 * @param[in] mythread スレッド番号
 * @return バッファが一杯になった時点での内側CSペア番号
 *
 * @ingroup integ-misc
 *
 * */
int ofmo_twoint_get_last_klcs( const int mythread ) {
    return LAST_KLCS[mythread];
}

/** A function that returns the start address of the buffer that stores the calculated two-electron integral
 *
 * buffered direct SCF法で用いる計算した縮約２電子積分を保存する
 * バッファ領域の先頭アドレスを返す。
 *
 * @param[in] mythread スレッド番号
 * @return バッファ領域の先頭アドレス
 *
 * @ingroup integ-misc
 *
 * */
double* ofmo_twoint_get_ebuf_eri( const int mythread ) {
    return EBUF_VAL[mythread];
}

/** Stores the four AO subscripts related to the two-electron integral stored in the buffer
 * 領域の先頭アドレスを返す関数
 *
 * @param[in] mythread スレッド番号
 * @return AO添字格納用の配列の先頭アドレス
 *
 * @ingroup integ-misc
 * */
short int* ofmo_twoint_get_ebuf_ind4( const int mythread ) {
    return EBUF_IND4[mythread];
}

/** A function that returns the number of reduced two-electron integrals stored in the buffer
 *
 * @param[in] mythread スレッド番号
 * @return Number of reduced two-electron integrals stored in the buffer
 *
 * @ingroup integ-misc
 *
 * */
size_t ofmo_twoint_get_stored_nzeri( const int mythread ) {
    return EBUF_NON_ZERO_ERI[mythread];
}

/* Functions related to the temporary domain of the integral used in the direct method and IFC4C */
void ofmo_twoint_set_stored_integ( const int mythread,
	const size_t ninteg ) {
    ETMP_NZERI[mythread] = ninteg;
}

size_t ofmo_twoint_get_max_stored_integ( const int mythread ) {
    return ETMP_MAX_NZERI[mythread];
}

size_t ofmo_twoint_get_stored_integ( const int mythread ) {
    return ETMP_NZERI[mythread];
}

double* ofmo_twoint_getadd_integ_val( const int mythread ) {
    return ETMP_VAL[mythread];
}

short* ofmo_twoint_getadd_integ_ind4( const int mythread ) {
    return ETMP_IND4[mythread];
}

// --------------------------------
// largest last_eri_type among all workers
static int g_LAST_ERI_TYPE = -1;

void ofmo_twoint_set_global_last_eri_type( const int last_eri_type ) {
    g_LAST_ERI_TYPE = last_eri_type;
}

int ofmo_twoint_get_global_last_eri_type( void ) {
    return g_LAST_ERI_TYPE;
}

// --------------------------------
// Dcs

float* ofmo_twoint_get_Dcs( void ) { return DENS_CS; }

void ofmo_twoint_free_Dcs( void ) {
  Free( DENS_CS );
  DENS_CS = NULL;
  BUF_NCS = 0;
  TWOINT_NCS = -1;
}

float* ofmo_twoint_alloc_Dcs( const int maxncs ) {
  if ( maxncs > BUF_NCS ) {
    ofmo_twoint_free_Dcs();
    DENS_CS = (float*)malloc( sizeof(float) * maxncs * maxncs );
    BUF_NCS = maxncs;
  }
  return DENS_CS;
}


#ifndef MIN2
#define MIN2(a, b) ((a)<(b))? (a): (b)
#endif
#ifndef MAX2
#define MAX2(a, b) ((a)>(b))? (a): (b)
#endif

float* ofmo_twoint_gen_Dcs( const int maxlqn, const int nao,
    const int leading_cs[], const double D[])
{
  int ncs = leading_cs[maxlqn+1];
  int nL[] = {1, 3, 6, 10};
  float *Dcs = NULL;
#pragma omp barrier
#pragma omp single
  Dcs = ofmo_twoint_alloc_Dcs(ncs);

  Dcs = DENS_CS;

//#pragma omp parallel
  {
    int iao00 = 0;
    for (int La=0; La<=maxlqn; La++) {
      int inao = nL[La];
      int ics0 = leading_cs[La];
      int ics1 = leading_cs[La+1];
#pragma omp for schedule(guided) nowait
      for (int ics=ics0; ics<ics1; ics++) {
        int iao0=iao00+inao*(ics-ics0);
        int iao1=iao0+inao;
        int jao0 = 0;
        for (int Lb=0; Lb<=La; Lb++) {
          int jnao = nL[Lb];
          int jcs0 = leading_cs[Lb];
          int jcs1 = leading_cs[Lb+1]-1;
          jcs1 = MIN2(jcs1, ics);
          for (int jcs=jcs0; jcs<=jcs1; jcs++) {
            double dmax = 0.0e0;
            for (int iao=iao0; iao<iao1; iao++) {
              int jao1=MIN2(jao0+jnao-1, iao);
              int iao2 = iao*(iao+1)/2;
              for (int jao=jao0; jao<=jao1; jao++) {
                dmax = MAX2(dmax, fabs(D[iao2+jao]));
              }
            }
            Dcs[ics*ncs+jcs] = (float)dmax;
            Dcs[jcs*ncs+ics] = (float)dmax;
            jao0 += jnao;
          }
        }
        //      iao0 += inao;
      }
      iao00 += inao*(ics1-ics0);
    }
  }
#pragma omp single
  TWOINT_NCS = ncs;
  return Dcs;
}

float ofmo_twoint_dmax6(const int i, const int j, const int k, const int l)
{
  int ncs = TWOINT_NCS;
  float *Dcs = DENS_CS;

  float dmax = 1.0;
  int i2 = i * ncs;
  int j2 = j * ncs;
  int k2 = k * ncs;

  dmax = MAX2(Dcs[i2+j], Dcs[k2+l]) * 4;
  dmax = MAX2(Dcs[i2+k], dmax);
  dmax = MAX2(Dcs[i2+l], dmax);
  dmax = MAX2(Dcs[j2+k], dmax);
  dmax = MAX2(Dcs[j2+l], dmax);

  return dmax;
}

float ofmo_twoint_dmax2(const int k, const int l)
{
  int ncs = TWOINT_NCS;
  float *Dcs = DENS_CS;

  return 2 * Dcs[k*ncs+l];
}

// --------------------------------
// EPS

static float TWOINT_EPS_PS4 = 1e-20F;
static float TWOINT_EPS_ERI = 1e-15F;
static float TWOINT_EPS_SCH = 1e-12F;

float ofmo_twoint_eps_eri(float eps)
{
  if (eps>0) TWOINT_EPS_ERI = eps;
  return TWOINT_EPS_ERI;
}

float ofmo_twoint_eps_ps4(float eps)
{
  if (eps>0) TWOINT_EPS_PS4 = eps;
  return TWOINT_EPS_PS4;
}

float ofmo_twoint_eps_sch(float eps)
{
  if (eps>0) TWOINT_EPS_SCH = eps;
  return TWOINT_EPS_SCH;
}

// --------------------------------

