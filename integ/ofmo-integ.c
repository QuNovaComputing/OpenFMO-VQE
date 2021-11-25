/**
 * @file ofmo-integ.c
 * A file that defines the top-level functions for performing
 * various molecular integrals required by the Hartree-Fock molecular
 * orbital method and the FMO method based on it.
 * */

/**
 * @defgroup integ 分子積分クラス
 * 通常のHartree-Fock(HF)分子軌道計算で用いる１電子積分（運動エネルギー
 * 積分、核-引力積分、重なり積分）や２電子積分だけでなく、FMO法に出現する
 * ４中心（３中心、２中心）の各種クーロン積分、および、カットオフテーブル
 * を作成する関数などを定義している。
 *
 * @ingroup ofmo
 *
 * */

/**
 * @defgroup integ-top Top class of molecular integral
 * @brief 積分計算を必要とする関数から呼ばれる最上位関数
 * @ingroup integ
 * */

/**
 * @defgroup integ-med 同じタイプの積分をまとめて行う関数群
 * @brief 同じタイプの積分をまとめて計算するための関数クラス
 * @ingroup integ
 * */

/**
 * @defgroup integ-core １つの縮約積分を計算する関数クラス
 * @brief １つの縮約積分を計算するための関数クラス
 * @ingroup integ
 * */

/**
 * @defgroup integ-misc 雑多な処理を行う関数クラス
 * @brief 積分計算以外の処理を行う関数クラス
 * @ingroup integ
 * */

/**
 * @defgroup integ-fmt 誤差関数計算を行う関数クラス
 * @brief 分子積分で必要となる誤差関数の計算を行う関数クラス
 * @ingroup integ
 * */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "ofmo-cutoff.h"
#include "ofmo-ifc2c.h"
#include "ofmo-oneint.h"
#include "ofmo-twoint.h"
#include "ofmo-twoint-buffer.h"
#include "ofmo-twoint-direct.h"
#include "ofmo-ifc4c.h"
#include "ofmo-ifc3c.h"
#include "fmt.h"
#include "fmt-m.h"

#include "ofmo-rys-xxxx.h"
#include "ofmo-os-xxxx.h"

#include "ofmo-def.h"
#include "ofmo-prof.h"
#include "ofmo-parallel.h"

#include "ofmo-tlog.h"

#ifdef OFMO_SKELETON
#include "rhf/skel-w2e.h"
#else
#define start_w2e()
#define set_w2e(Labcd)
#endif

double x_coef;    //DFT; HF exchange coef.

//#define Free(a) if ( a != NULL ) free( a ); a = NULL

extern int ofmo_twoint_xxxx(
	const int *pnworkers, const int *pworkerid,
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// for partially direct SCF
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs );
extern int ofmo_twoint_direct_xxxx(
	const int *pnworkers, const int *pworkerid,
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// for direct SCF
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_twoint_rys_xxxx(
	const int *pnworkers, const int *pworkerid,
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// for partially direct SCF
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs );
extern int ofmo_twoint_direct_rys_xxxx(
	const int *pnworkers, const int *pworkerid,
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// for direct SCF
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] );

extern int ofmo_OS_integ_init( const int maxlqn );
extern int ofmo_oneint_gen_init( const int maxlqn );

extern int ofmo_Rys_integ_init( const int maxlqn );

/*// for Fortran code
extern void fmt4_initialize_();
extern void fmt4_gen_initialize_();
extern void fort_const_();*/

#ifdef USE_CUDA
#include "cuda/cudalib.h"
#include "cuda/cuda-integ.h"
#include "cuda/cuda-ifc4c.h"
#include "cuda/cuda-ifc4c-calc.h"
extern int cuda_fmt_initialize();
#endif

/* ====================================================
   負荷分散などのための制御変数に関係する関数群
   ==================================================== */
// to control load-balancing
static int *target_type = NULL;
static size_t *loop_offset = NULL;
static int OFMO_MAX_THREADS  = 1;

static void finalize_ctrl() {
    Free( target_type );
    Free( loop_offset );
    OFMO_MAX_THREADS  = 1;
}

static int init_ctrl() {
    static int called = false;
    int maxthreads, i;
    if ( !called ) {
	maxthreads = omp_get_max_threads();
	target_type = (int*)malloc( sizeof(int) * maxthreads );
	loop_offset = (size_t*)malloc( sizeof(size_t) * maxthreads );
	OFMO_MAX_THREADS = maxthreads;
	for ( i=0; i<maxthreads; i++ ) {
	    target_type[i] = -1;
	    loop_offset[i] = 0;
	}

	atexit( finalize_ctrl );
	called = true;
    }
    return 0;
}

size_t ofmo_integ_get_loop_offset( const int mythread ) {
    return loop_offset[mythread];
}

void ofmo_integ_set_loop_offset( const int mythread,
	const size_t offset ) {
    loop_offset[mythread] = offset;
}

void ofmo_integ_set_target_type( const int mythread,
	const int ttype ) {
    target_type[mythread] = ttype;
}

static int ofmo_integ_get_target_type( const int mythread ) {
    return target_type[mythread];
}

int ofmo_integ_init( int maxlqn ) {
    static int called = false;
    if ( called == false )  {
	init_ctrl();	// added for load-balancing
	fmt_initialize( maxlqn );
        fmt_m_init();
#ifdef USE_CUDA
        int ret = cuda_fmt_initialize();
        if (ret<0) exit(1);
#endif
	ofmo_twoint_init();
	ofmo_ifc3c_os_init();
	ofmo_ifc3c_rys_init();
	ofmo_ifc2c_init();
	ofmo_oneint_init();
	ofmo_cutoff_init();
	//
	ofmo_OS_integ_init( maxlqn );
	ofmo_oneint_gen_init( maxlqn );
	//
	ofmo_Rys_integ_init( maxlqn );
	called = true;
    }
    return 0;
}

extern int ofmo_cutoff_xx(
	// input arguments
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[], const int shel_add[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	// output arguments
	int leading_cs_pair[],
	double csp_schwarz[], int csp_ics[], int csp_jcs[],
	int csp_leading_ps_pair[],
	double psp_zeta[], double psp_dkps[], double psp_xiza[] );

/* top level code for making cutoff table */
static int (*calc_schwarz[]) ( 
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[], const int shel_add[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	int leading_cs_pair[],
	double csp_schwarz[], int csp_ics[], int csp_jcs[],
	int csp_leading_ps_pair[],
	double psp_zeta[], double psp_dkps[], double psp_xiza[]) = {
    ofmo_cutoff_ss_, /*ofmo_cutoff_ps_, ofmo_cutoff_pp_,
    ofmo_cutoff_ds_, ofmo_cutoff_dp_, ofmo_cutoff_dd_,*/
    ofmo_cutoff_xx, ofmo_cutoff_xx,
    ofmo_cutoff_xx, ofmo_cutoff_xx, ofmo_cutoff_xx
};

/** Cutoff table creation function for using Schwarz inequalities
 * @ingroup integ-top
 *
 * Create a cutoff table using Schwarz's inequality.
 * Given the sort basis function data, a cutoff table for using
 * Schwarz's inequalities is calculated and returned.
 *
 * @attention 
 * @li Various arrays for output need to be reserved at the time of calling.
 * @li It is thread-safe if the various output arrays
 * (such as \ c csp_schwarz [] and \ c csp_ics []) are not shared by multiple threads.
 *
 * @param[in] maxlqn Maximum orbital quantum number
 * @param[in] leading_cs[lqn] First CS number of orbital quantum number \c lqn
 * @param[in] shel_tem[ics] CS reduction length of CS number \c ics
 * @param[in] shel_atm[ics] CS number \c ics CS number of the atom to which CS belongs
 * @param[in] shel_add[ics] CS number \c ics CS contains the first PS number of PS
 * @param[in] atom_x[iat] Atom number \c iat x coordinate (au unit)
 * @param[in] atom_y[iat] Atom number \c iat y coordinate (au unit)
 * @param[in] atom_z[iat] Atom number \c iat z coordinate (au unit)
 * @param[in] prim_exp[ips] PS orbit index of PS number \c ips
 * @param[in] prim_coe[ips] PS number \c ips PS standardization constant-included reduction factor
 *
 * @param[out] leading_cs_pair[itype] CS pair type number First CS pair number of \c itype
 * @param[out] csp_schwarz[icsp] Schwarz integral of CS pair number \c icsp
 * @param[out] csp_ics[icsp] CS pair number \c icsp's first CS number
 * @param[out] csp_jcs[icsp] CS pair number The second CS number for \c icsp.
 * However, it is \ f $ \ tt {csp \ _ics [icsp]} \ ge \ tt {csp \ _jcs [icsp]} \ f $.
 * @param[out] csp_leading_ps_pair[icsp]  CSペア番号 \c icsp に含まれる
 *     PSペアの先頭PSペア番号
 * @param[out] psp_zeta[ipsp] PSペア番号 \c ipsp の軌道指数和
 *     \f$ \zeta = \zeta_a + \zeta_b \f$
 * @param[out] psp_dkps[ipsp] PSペア番号 \c ipsp の線型結合定数
 *     \f[ K_{ab} = \sqrt2 \pi^{5/4} \frac1{\zeta_a+\zeta_b}
 *     \exp\left[ -\frac{\zeta_a \zeta_b}{\zeta_a + \zeta_b}
 *     ( \boldmath A \unboldmath - \boldmath B \unboldmath )^2
 *     \right]\f]
 * @param[out] psp_xiza[ipsp] PSペア番号 \c ipsp の
 *     \f$ \frac{\xi}{\zeta_a} = \frac{\zeta_b}{\zeta_a+\zeta_b} \f$
 *
 * @retval  0 正常終了
 * @retval -1 異常終了（サポートしていない積分タイプがあったなど）
 *
 * */
int ofmo_cutoff_make_table(
	// input arguments
	const int maxlqn, const int leading_cs[],
	const int shel_tem[], const int shel_atm[], const int shel_add[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	// output arguments
	int leading_cs_pair[],
	double csp_schwarz[], int csp_ics[], int csp_jcs[],
	int csp_leading_ps_pair[],
	double psp_zeta[], double psp_dkps[], double psp_xiza[] ) {
    int La, Lb, Lab;

    leading_cs_pair[0] = 0;
    csp_leading_ps_pair[0] = 0;
    for ( La=0; La<=maxlqn; La++ ) {
	for ( Lb=0; Lb<=La; Lb++ ) {
	    Lab = La*(La+1)/2 + Lb;
	    if ( Lab == 0 ) {
		ofmo_cutoff_ss_(
			&La, &Lb, leading_cs,
			shel_tem, shel_atm, shel_add,
			atom_x, atom_y, atom_z,
			prim_exp, prim_coe,

			leading_cs_pair,
			csp_schwarz, csp_ics, csp_jcs,
			csp_leading_ps_pair,
			psp_zeta, psp_dkps, psp_xiza );
	    } else {
		ofmo_cutoff_xx(
			&La, &Lb, leading_cs,
			shel_tem, shel_atm, shel_add,
			atom_x, atom_y, atom_z,
			prim_exp, prim_coe,

			leading_cs_pair,
			csp_schwarz, csp_ics, csp_jcs,
			csp_leading_ps_pair,
			psp_zeta, psp_dkps, psp_xiza );
	    }
#ifdef SORT_CSP
            ofmo_cutoff_sort_( La, Lb, leading_cs,
                shel_tem, shel_atm, shel_add,
                atom_x, atom_y, atom_z,
                prim_exp, prim_coe,
                leading_cs_pair,
                csp_schwarz, csp_ics, csp_jcs,
                csp_leading_ps_pair,
                psp_zeta, psp_dkps, psp_xiza );
#endif
	}	// Lb
    }		// La
    return 0;
}

static int isum( const int n, const int a[] ) {
    int i, sum;
    for ( i=0, sum=0; i<n; i++ ) sum += a[i];
    return sum;
}

extern int ofmo_oneint_xx(
	const int *pnworkers, const int *pworkerid,
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[],
	const int shel_add[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	const int *pnat, const int atomic_number[],
	double S[], double H[] );

static int (*oneint_func[]) (
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb,
	// basis set data for fragment
	const int leading_cs[],
	const int shel_tem[], const int shel_atm[],
	const int shel_add[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	const int *pnat, const int atomic_number[],
	double S[], double H[] ) = {
    /*ofmo_oneint_ss_, ofmo_oneint_ps_, ofmo_oneint_pp_,
    ofmo_oneint_ds_, ofmo_oneint_dp_, ofmo_oneint_dd_,*/
    ofmo_oneint_ss__, ofmo_oneint_ps__, ofmo_oneint_pp__,
    /*ofmo_oneint_xx, ofmo_oneint_xx,
    ofmo_oneint_xx, ofmo_oneint_xx, ofmo_oneint_xx,*/
    ofmo_oneint_ds__, ofmo_oneint_dp__, ofmo_oneint_dd__,
};

/** Sorted one-electron integral calculation function
 * @ingroup integ-top
 *
 * Given a sorted basis set, it calculates the one-electron integral
 * and returns the sorted overlap integral, the one-electron Hamilton matrix.
 *
 * @attention
 * @li 出力用の各種配列は、呼び出し時には確保されている必要がある
 * @li スレッド並列時の関数呼び出しは、スレッド並列領域内から行う必要
 *     がある
 * @li １プロセスで実行する場合には、関数終了時点で（スレッド
 *     並列化時には全スレッドが関数から返った時点で）完全な重なり行列や
 *     一電子ハミルトン行列が得られる。
 * @li 複数プロセスで実行する場合には、関数終了時点では部分の結果しか
 *     得られていない。完全な結果を得るためには、\c MPI_Allreduce などの
 *     関数を用いたリダクション処理が必要である。
 *
 * @param[in] nworkers 計算に用いるワーカプロセス（スレッド）数
 * @param[in] workerid 各ワーカプロセス（スレッド）のID。
 *     \f$ 0\le\tt{workerid}<\tt{nworkers} \f$である。
 * @param[in] maxlqn 最大軌道量子数
 * @param[in] leading_cs[lqn] 軌道量子数 \c lqn の先頭CS番号
 * @param[in] shel_tem[ics] CS番号 \c ics のCSの縮約長
 * @param[in] shel_atm[ics] CS番号 \c ics のCSが属する原子の番号
 * @param[in] shel_add[ics] CS番号 \c ics のCSに含まれるPSの先頭PS番号
 * @param[in] shel_ini[ics] CS番号 \c ics のCSに含まれるAOの先頭AO番号
 * @param[in] atom_x[iat] 原子の番号 \c iat のx座標（au単位）
 * @param[in] atom_y[iat] 原子の番号 \c iat のy座標（au単位）
 * @param[in] atom_z[iat] 原子の番号 \c iat のz座標（au単位）
 * @param[in] prim_exp[ips] PS番号 \c ips のPSの軌道指数
 * @param[in] prim_coe[ips] PS番号 \c ips のPSの規格化定数込みの縮約係数
 *
 * @param[in] nat 原子数
 * @param[in] atomic_number[iat] 原子の番号 \c iat の原子番号
 *
 * @param[out] S[] 重なり行列(圧縮"U"形式)。
 * @param[out] H[] 一電子ハミルトン行列(圧縮"U"形式)。
 *
 * @retval  0 正常終了
 * @retval -1 異常終了（いま(2011/06/13)のところ考えていない）
 * */
int ofmo_integ_oneint_sorted(
	const int nworkers, const int workerid,
	const int maxlqn, const int leading_cs[],
	const int shel_tem[], const int shel_atm[],
	const int shel_add[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	const int nat, const int atomic_number[],
	double S[], double H[]) {
    int La, Lb, Lab, sum;
    //double dsum;

    sum = isum( nat, atomic_number );
    ofmo_oneint_set_sum_atomic_numbers( sum );

    for ( La=0; La<=maxlqn; La++ ) {
	for ( Lb=0; Lb<=La; Lb++ ) {
	    Lab = La*(La+1)/2+Lb;
	    if ( Lab == 0 ) {
		ofmo_oneint_ss__(
			&nworkers, &workerid,
			&La, &Lb, leading_cs,
			shel_tem, shel_atm, shel_add, shel_ini,
			atom_x, atom_y, atom_z, prim_exp, prim_coe,
			&nat, atomic_number, S, H );

	    } else {
		ofmo_oneint_xx(
			&nworkers, &workerid,
			&La, &Lb, leading_cs,
			shel_tem, shel_atm, shel_add, shel_ini,
			atom_x, atom_y, atom_z, prim_exp, prim_coe,
			&nat, atomic_number, S, H );
	    }
	}
    }
    return 0;
}

static int (*calc_twoint_buffer[])(
	const int *pnworkers, const int *pworkerid,
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// for partially direct SCF
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs ) = {
    ofmo_twoint_buffer_ssss__,	// (SS,SS)はすべて同じコードを使用
    /*// original (OS, 個別、d関数は元の並び）
    ofmo_twoint_buffer_psss__, ofmo_twoint_buffer_psps__,
    ofmo_twoint_buffer_ppss__, ofmo_twoint_buffer_ppps__,
    ofmo_twoint_buffer_pppp__, ofmo_twoint_buffer_dsss__,
    ofmo_twoint_buffer_dsps__, ofmo_twoint_buffer_dspp__,
    ofmo_twoint_buffer_dsds__, ofmo_twoint_buffer_dpss__,
    ofmo_twoint_buffer_dpps__, ofmo_twoint_buffer_dppp__,
    ofmo_twoint_buffer_dpds__, ofmo_twoint_buffer_dpdp__,
    ofmo_twoint_buffer_ddss__, ofmo_twoint_buffer_ddps__,
    ofmo_twoint_buffer_ddpp__, ofmo_twoint_buffer_ddds__,
    ofmo_twoint_buffer_dddp__, ofmo_twoint_buffer_dddd__,*/
    /*// Obara-Saika(一般式、C言語)
    ofmo_twoint_xxxx, ofmo_twoint_xxxx,
    ofmo_twoint_xxxx, ofmo_twoint_xxxx, ofmo_twoint_xxxx,
    ofmo_twoint_xxxx, ofmo_twoint_xxxx, ofmo_twoint_xxxx,
    ofmo_twoint_xxxx, ofmo_twoint_xxxx, ofmo_twoint_xxxx,
    ofmo_twoint_xxxx, ofmo_twoint_xxxx, ofmo_twoint_xxxx,
    ofmo_twoint_xxxx, ofmo_twoint_xxxx, ofmo_twoint_xxxx,
    ofmo_twoint_xxxx, ofmo_twoint_xxxx, ofmo_twoint_xxxx,*/
    // Obara-Saika（個別、C言語）
    ofmo_twoint_os_psss, ofmo_twoint_os_psps,
    ofmo_twoint_os_ppss, ofmo_twoint_os_ppps, ofmo_twoint_os_pppp,
    ofmo_twoint_os_dsss, ofmo_twoint_os_dsps, ofmo_twoint_os_dspp,
    ofmo_twoint_os_dsds, ofmo_twoint_os_dpss, ofmo_twoint_os_dpps,
    ofmo_twoint_os_dppp, ofmo_twoint_os_dpds, ofmo_twoint_os_dpdp,
    ofmo_twoint_os_ddss, ofmo_twoint_os_ddps, ofmo_twoint_os_ddpp,
    ofmo_twoint_os_ddds, ofmo_twoint_os_dddp, ofmo_twoint_os_dddd,
    /*// Rys求積法(一般式、C言語)
    ofmo_twoint_rys_xxxx, ofmo_twoint_rys_xxxx,
    ofmo_twoint_rys_xxxx, ofmo_twoint_rys_xxxx, ofmo_twoint_rys_xxxx,
    ofmo_twoint_rys_xxxx, ofmo_twoint_rys_xxxx, ofmo_twoint_rys_xxxx,
    ofmo_twoint_rys_xxxx, ofmo_twoint_rys_xxxx, ofmo_twoint_rys_xxxx,
    ofmo_twoint_rys_xxxx, ofmo_twoint_rys_xxxx, ofmo_twoint_rys_xxxx,
    ofmo_twoint_rys_xxxx, ofmo_twoint_rys_xxxx, ofmo_twoint_rys_xxxx,
    ofmo_twoint_rys_xxxx, ofmo_twoint_rys_xxxx, ofmo_twoint_rys_xxxx,*/
    /* // Rsy求積法（個別、C言語）
    ofmo_twoint_rys_psss, ofmo_twoint_rys_psps,
    ofmo_twoint_rys_ppss, ofmo_twoint_rys_ppps,
    ofmo_twoint_rys_pppp, ofmo_twoint_rys_dsss,
    ofmo_twoint_rys_dsps, ofmo_twoint_rys_dspp,
    ofmo_twoint_rys_dsds, ofmo_twoint_rys_dpss,
    ofmo_twoint_rys_dpps, ofmo_twoint_rys_dppp,
    ofmo_twoint_rys_dpds, ofmo_twoint_rys_dpdp,
    ofmo_twoint_rys_ddss, ofmo_twoint_rys_ddps,
    ofmo_twoint_rys_ddpp, ofmo_twoint_rys_ddds,
    ofmo_twoint_rys_dddp, ofmo_twoint_rys_dddd,*/
};

static int (*calc_twoint_direct[])(
	// paralleization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis set & cutoff table data
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// concerned about buffered direct method
	const long *petmp_max_nzeri, long *petmp_non_zero_eri,
	double etmp_val[], short int etmp_ind4[],
	const int *plast_ijcs, const int *plast_klcs,
	// density matrix & G-matrix data
	const int *pnao, const double Ds[], double G[] ) = {
    ofmo_twoint_direct_ssss__,
    /*// Obara-Saika式（一般式、C言語）
    ofmo_twoint_direct_xxxx, ofmo_twoint_direct_xxxx,
    ofmo_twoint_direct_xxxx, ofmo_twoint_direct_xxxx,
    ofmo_twoint_direct_xxxx, ofmo_twoint_direct_xxxx,
    ofmo_twoint_direct_xxxx, ofmo_twoint_direct_xxxx,
    ofmo_twoint_direct_xxxx, ofmo_twoint_direct_xxxx,
    ofmo_twoint_direct_xxxx, ofmo_twoint_direct_xxxx,
    ofmo_twoint_direct_xxxx, ofmo_twoint_direct_xxxx,
    ofmo_twoint_direct_xxxx, ofmo_twoint_direct_xxxx,
    ofmo_twoint_direct_xxxx, ofmo_twoint_direct_xxxx,
    ofmo_twoint_direct_xxxx, ofmo_twoint_direct_xxxx,*/
    // Obara-Saika式（個別、C言語）
    ofmo_twoint_direct_os_psss, ofmo_twoint_direct_os_psps,
    ofmo_twoint_direct_os_ppss, ofmo_twoint_direct_os_ppps,
    ofmo_twoint_direct_os_pppp, ofmo_twoint_direct_os_dsss,
    ofmo_twoint_direct_os_dsps, ofmo_twoint_direct_os_dspp,
    ofmo_twoint_direct_os_dsds, ofmo_twoint_direct_os_dpss,
    ofmo_twoint_direct_os_dpps, ofmo_twoint_direct_os_dppp,
    ofmo_twoint_direct_os_dpds, ofmo_twoint_direct_os_dpdp,
    ofmo_twoint_direct_os_ddss, ofmo_twoint_direct_os_ddps,
    ofmo_twoint_direct_os_ddpp, ofmo_twoint_direct_os_ddds,
    ofmo_twoint_direct_os_dddp, ofmo_twoint_direct_os_dddd,
    /*// Rys求積法（一般式、C言語）
    ofmo_twoint_direct_rys_xxxx, ofmo_twoint_direct_rys_xxxx,
    ofmo_twoint_direct_rys_xxxx, ofmo_twoint_direct_rys_xxxx,
    ofmo_twoint_direct_rys_xxxx, ofmo_twoint_direct_rys_xxxx,
    ofmo_twoint_direct_rys_xxxx, ofmo_twoint_direct_rys_xxxx,
    ofmo_twoint_direct_rys_xxxx, ofmo_twoint_direct_rys_xxxx,
    ofmo_twoint_direct_rys_xxxx, ofmo_twoint_direct_rys_xxxx,
    ofmo_twoint_direct_rys_xxxx, ofmo_twoint_direct_rys_xxxx,
    ofmo_twoint_direct_rys_xxxx, ofmo_twoint_direct_rys_xxxx,
    ofmo_twoint_direct_rys_xxxx, ofmo_twoint_direct_rys_xxxx,
    ofmo_twoint_direct_rys_xxxx, ofmo_twoint_direct_rys_xxxx,*/
    /*// Rys求積法（個別、C言語）
    ofmo_twoint_direct_rys_psss, ofmo_twoint_direct_rys_psps,
    ofmo_twoint_direct_rys_ppss, ofmo_twoint_direct_rys_ppps,
    ofmo_twoint_direct_rys_pppp, ofmo_twoint_direct_rys_dsss,
    ofmo_twoint_direct_rys_dsps, ofmo_twoint_direct_rys_dspp,
    ofmo_twoint_direct_rys_dsds, ofmo_twoint_direct_rys_dpss,
    ofmo_twoint_direct_rys_dpps, ofmo_twoint_direct_rys_dppp,
    ofmo_twoint_direct_rys_dpds, ofmo_twoint_direct_rys_dpdp,
    ofmo_twoint_direct_rys_ddss, ofmo_twoint_direct_rys_ddps,
    ofmo_twoint_direct_rys_ddpp, ofmo_twoint_direct_rys_ddds,
    ofmo_twoint_direct_rys_dddp, ofmo_twoint_direct_rys_dddd,*/
};

/** @example oneint-serial.c
 * １電子積分関数のシリアル実行時のコード例。
 * */

/** @example oneint-mt.c
 * １電子積分関数のスレッド並列実行時のコード例。
 *
 * 関数は、スレッド並列領域内から呼ばれている。
 * 関数が終了して全スレッドで同期がとれた時点で、完全な行列が得られる。
 * */

/** @example oneint-mpi.c
 * １電子積分関数のフラットMPI並列時のコード例。
 *
 * 関数呼び出し後にリダクション処理を行うことで、完全な行列が得られる。
 * */

/** @example oneint-hybrid.c
 * １電子積分関数のOpenMPとMPIによるハイブリッド並列化の例。
 *
 * フラットMPIの場合と同様に、関数呼び出し後にリダクション処理を行う
 * ことで、完全な行列が得られる。
 *
 * */

/** Electron repulsion integral (two-electron integral) calculation function（１）
 * @ingroup integ-top
 *
 * Calculates and stores two-electron integrals that fit in a buffer of the specified size.
 * When the buffer is full, no further integral calculations are performed.
 * Also, the two-electron Hamilton matrix (G matrix) is not calculated.
 *
 * @attention
 * @li It must be called once before calling the G matrix generation function
 * 	   \c ofmo_integ_gen_gmat . Calling the \c ofmo_integ_gen_gmat function
 *     without calling this function may result in strange results.
 * @li When performing SCF calculation multiple times in one execution like
 *     FMO calculation, it is necessary to call it once each time SCF calculation is performed.
 * @li Thread parallelization is performed by OpenMP. In order to execute threads
 *     in parallel, this function must be called in the thread parallel area.
 * @li If you set nworkers and \c workerid appropriately, it also supports hybrid
 *     parallelization in combination with MPI.
 * @li A buffer for storing the two-electron integral is reserved inside the
 *     function. In addition, the secured buffer is released at the end of the
 *     program.
 *
 * @param[in] nworkers Number of worker processes (threads) used for calculation
 * @param[in] workerid ID of each worker process (thread).
 *     \f$ 0\le\tt{workerid}<\tt{nworkers} \f$である。
 * @param[in] ebuf_buffer_size_mb Buffer size (in MB) for storing the value of
 *     the two-electron integral. This buffer stores not only the integral value
 *     but also the four subscripts.
 * @param[in] maxlqn Maximum orbital quantum number
 * @param[in] shel_atm[ics] CS番号 \c ics のCSが属する原子の番号
 * @param[in] shel_ini[ics] CS番号 \c ics のCSに含まれるAOの先頭AO番号
 * @param[in] atom_x[iat] 原子の番号 \c iat のx座標（au単位）
 * @param[in] atom_y[iat] 原子の番号 \c iat のy座標（au単位）
 * @param[in] atom_z[iat] 原子の番号 \c iat のz座標（au単位）
 * @param[in] leading_cs_pair[itype] CSペアタイプ番号 \c itype の
 *     先頭CSペア番号
 * @param[in] csp_schwarz[icsp] CSペア番号 \c icsp のSchwarz積分
 * @param[in] csp_ics[icsp] CSペア番号 \c icsp の1つ目のCS番号
 * @param[in] csp_jcs[icsp] CSペア番号 \c icsp の2つめのCS番号。ただし、
 *     \f$ \tt{csp\_ics[icsp]} \ge \tt{csp\_jcs[icsp]} \f$ である。
 * @param[in] csp_leading_ps_pair[icsp]  CSペア番号 \c icsp に含まれる
 *     PSペアの先頭PSペア番号
 * @param[in] psp_zeta[ipsp] PSペア番号 \c ipsp の軌道指数和
 *     \f$ \zeta = \zeta_a + \zeta_b \f$
 * @param[in] psp_dkps[ipsp] PSペア番号 \c ipsp の線型結合定数
 *     \f[ K_{ab} = \sqrt2 \pi^{5/4} \frac1{\zeta_a+\zeta_b}
 *     \exp\left[ -\frac{\zeta_a \zeta_b}{\zeta_a + \zeta_b}
 *     ( \boldmath A \unboldmath - \boldmath B \unboldmath )^2
 *     \right]\f]
 * @param[in] psp_xiza[ipsp] PSペア番号 \c ipsp の
 *     \f$ \frac{\xi}{\zeta_a} = \frac{\zeta_b}{\zeta_a+\zeta_b} \f$
 *
 * 
 * @retval  0 All integrals are stored in the buffer and exit（in core SCF）
 * @retval  1 End when the buffer is full（partially direct SCF）
 *
 * */
int ofmo_integ_twoint_first(
	// parallelization
	const int nworkers, const int workerid,
	// buffer size
	const size_t ebuf_buffer_size_mb,
	// basis set & cutoff table data
	const int maxlqn, const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const int leading_cs_pair[], const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[] ) {
    int La, Lb, Lc, Ld, Lab, Lcd, Labcd;
    int mythread = 0, ret_code;
    int last_ijcs, last_klcs;
    long ebuf_non_zero_eri, ebuf_max_nzeri;
    short int *ebuf_ind4;
    double *ebuf_eri;
    // debug
    char *CS = "spdfg";
    // added for load-balancing
    int local_id=workerid;
    size_t offset;

    mythread = omp_get_thread_num();
    ebuf_max_nzeri = (long)
	ofmo_twoint_set_buffer_size( mythread, ebuf_buffer_size_mb );
    ebuf_ind4 = ofmo_twoint_get_ebuf_ind4( mythread );
    ebuf_eri  = ofmo_twoint_get_ebuf_eri( mythread );
    last_ijcs = last_klcs = -1;
    ebuf_non_zero_eri = 0;

    // added for load-balancing
    //offset   = ofmo_integ_get_loop_offset( mythread );
    offset = 0;

    // initialize
    ofmo_twoint_set_last_eri_type( mythread, -1 );
    if (ebuf_max_nzeri<=0) return -1;

    for ( La=0; La<=maxlqn; La++ ) {
	for ( Lb=0; Lb<=La; Lb++ ) {
	    Lab = La*(La+1)/2 + Lb;
	    for ( Lc=0; Lc<=La; Lc++ ) {
		for ( Ld=0; Ld<=(Lc==La? Lb : Lc ); Ld++ ) {
		    Lcd = Lc*(Lc+1)/2 + Ld;
		    Labcd = Lab*(Lab+1)/2 + Lcd;
#ifdef USE_CUDA
                    if (cuda_use_Device() && cuda_get_optCPU(Labcd)!=0) {
                      offset+=(leading_cs_pair[Lab+1]-leading_cs_pair[Lab]);
                      continue;
                    }
#endif
		    local_id =
			(int)((offset+(size_t)workerid)%(size_t)nworkers);

		    ret_code = calc_twoint_buffer[Labcd](
			&nworkers, &local_id,
			//&nworkers, &workerid,
			&La, &Lb, &Lc, &Ld,
			shel_atm, shel_ini, atom_x, atom_y, atom_z,
			leading_cs_pair, csp_schwarz, csp_ics, csp_jcs,
			csp_leading_ps_pair,
			psp_zeta, psp_dkps, psp_xiza,
			&ebuf_max_nzeri, &ebuf_non_zero_eri,
			ebuf_eri, ebuf_ind4,
			&last_ijcs, &last_klcs );
		    if ( ret_code == OFMO_EBUF_FULL ) {
			ofmo_twoint_set_last_eri_type( mythread, Labcd );
			ofmo_twoint_set_last_ijcs( mythread, last_ijcs );
			ofmo_twoint_set_last_klcs( mythread, last_klcs );
			ofmo_twoint_set_stored_nzeri( mythread,
				(size_t)ebuf_non_zero_eri );
			//ofmo_integ_set_loop_offset( mythread, offset );
			//return 0;
			return Labcd;
		    }
		    // added for load-balancing
		    offset+=(leading_cs_pair[Lab+1]-leading_cs_pair[Lab]);
		}// Ld
	    }	// Lc
	}	// Lb
    }		// La
    ofmo_twoint_set_last_eri_type( mythread, 100000 );
    ofmo_twoint_set_last_ijcs( mythread, last_ijcs );
    ofmo_twoint_set_last_klcs( mythread, last_klcs );
    ofmo_twoint_set_stored_nzeri( mythread, (size_t)ebuf_non_zero_eri );
    // added for load-balancing
    //ofmo_integ_set_loop_offset( mythread, offset );
    //return 0;
    return 100000;
}

int ofmo_integ_add_fock( const int nao, const size_t nstored_eri,
	const double eri_val[], const short int eri_ind4[],
	const double D[], double G[] ) {
    int i, j, k, l, ij, kl, ik, il, jk, jl, i0, j0;
    size_t ix, ix4;
    double x, x4;

    /*// debug
    {
	int mythread;
	mythread = omp_get_thread_num();
	printf("thrd= %d,  nzeri= %lld\n", mythread,
		(long long)nstored_eri );
	fflush( stdout );
    }*/
    for ( ix=0, ix4=0; ix<nstored_eri; ix++, ix4+=4 ) {
	x = eri_val[ix];
	i = (int)eri_ind4[ix4+0];
	j = (int)eri_ind4[ix4+1];
	k = (int)eri_ind4[ix4+2];
	l = (int)eri_ind4[ix4+3];

	x4 = x * 4.e0;
	x  = x * x_coef;   //DFT
	i0 = i*nao;
	j0 = j*nao;
	ij = i0 + j;
	ik = i0 + k;
	il = i0 + l;
	jk = j0 + k;
	jl = j0 + l;
	kl = k*nao + l;
	G[ij] += D[kl]*x4;
	G[kl] += D[ij]*x4;
	G[ik] -= D[jl]*x;
	G[il] -= D[jk]*x;
	G[jk] -= D[il]*x;
	G[jl] -= D[ik]*x;
    }
    return 0;
}
/*
 * Compute the Fock matrix using the two-electron integration stored in the buffer
 * However, the density matrix D and Fock matrix G of the arguments that are not accurate
 * unless the square matrix is ​​folded into a compressed format are the thread parallel parts
 * that are treated as square matrices.
 * Calling G expects to be initialized The following variables should be different
 * for each thread.
 *  non_zero_eri
 *  ebuf_eri[]
 *  ebuf_ind4[]
 *  G[]
 * */
static int ofmo_twoint_fock_incore_partial(
	const int mythread, const int nao, const double D[],
	double G[] ) {
    short int *ebuf_ind4;
    double *ebuf_eri;
    size_t non_zero_eri;
    ebuf_eri     = ofmo_twoint_get_ebuf_eri( mythread );
    ebuf_ind4    = ofmo_twoint_get_ebuf_ind4( mythread );
    non_zero_eri = ofmo_twoint_get_stored_nzeri( mythread );
    ofmo_integ_add_fock( nao, non_zero_eri, ebuf_eri, ebuf_ind4, D, G );

    return 0;
}

static void unpack_matrix( const int n, const double SRC[], double DST[] )
{
    int i, j, ij;
    ij = 0;
    for ( i=0; i<n; i++ ) {
	for ( j=0; j<=i; j++ ) {
	    DST[i*n+j] = DST[j*n+i] = SRC[ij];
	    ij++;
	}
    }
}

/** Electron repulsion integral (two-electron integral) calculation function (2)
 * @ingroup integ-top
 *
 * The two-electron Hamilton matrix (G matrix) is calculated using the two-electron integral
 * calculated and saved by calling the \c ofmo_integ_twoint function.
 * Integrals that are not stored in the buffer are calculated and added to the G matrix.
 *
 * @attention
 * @li This function performs thread parallelization using OpenMP.
 *     In order to execute threads in parallel, it is necessary to call
 *     this function from within the thread parallel area.
 * @li If \c nworkers and worker ID \c workerid are set appropriately, hybrid
 *     parallel execution of OpenMP and MPI is possible. In order to obtain
 *     a complete G matrix at the time of MPI parallel, it is necessary to
 *     perform reduction processing using the \c MPI_Allreduce function etc.
 *     after the end of this function.
 * @li \c ofmo_integ_twoint_first must be called in advance.
 * 	   Otherwise, the results may be strange.
 * @li The obtained G matrix is ​​sorted by the magnitude of the orbital
 *     quantum number. If you want the original G matrix, you need to sort
 *     the elements.
 *
 * @param[in] nworkers 計算に用いるワーカプロセス（スレッド）数
 * @param[in] workerid 各ワーカプロセス（スレッド）のID。
 *     \f$ 0\le\tt{workerid}<\tt{nworkers} \f$である。
 * @param[in] maxlqn 最大軌道量子数
 * @param[in] shel_atm[ics] CS番号 \c ics のCSが属する原子の番号
 * @param[in] shel_ini[ics] CS番号 \c ics のCSに含まれるAOの先頭AO番号
 * @param[in] atom_x[iat] 原子の番号 \c iat のx座標（au単位）
 * @param[in] atom_y[iat] 原子の番号 \c iat のy座標（au単位）
 * @param[in] atom_z[iat] 原子の番号 \c iat のz座標（au単位）
 * @param[in] leading_cs_pair[itype] CSペアタイプ番号 \c itype の
 *     先頭CSペア番号
 * @param[in] csp_schwarz[icsp] CSペア番号 \c icsp のSchwarz積分
 * @param[in] csp_ics[icsp] CSペア番号 \c icsp の1つ目のCS番号
 * @param[in] csp_jcs[icsp] CSペア番号 \c icsp の2つめのCS番号。ただし、
 *     \f$ \tt{csp\_ics[icsp]} \ge \tt{csp\_jcs[icsp]} \f$ である。
 * @param[in] csp_leading_ps_pair[icsp]  CSペア番号 \c icsp に含まれる
 *     PSペアの先頭PSペア番号
 * @param[in] psp_zeta[ipsp] PSペア番号 \c ipsp の軌道指数和
 *     \f$ \zeta = \zeta_a + \zeta_b \f$
 * @param[in] psp_dkps[ipsp] PSペア番号 \c ipsp の線型結合定数
 *     \f[ K_{ab} = \sqrt2 \pi^{5/4} \frac1{\zeta_a+\zeta_b}
 *     \exp\left[ -\frac{\zeta_a \zeta_b}{\zeta_a + \zeta_b}
 *     ( \boldmath A \unboldmath - \boldmath B \unboldmath )^2
 *     \right]\f]
 * @param[in] psp_xiza[ipsp] PSペア番号 \c ipsp の
 *     \f$ \frac{\xi}{\zeta_a} = \frac{\zeta_b}{\zeta_a+\zeta_b} \f$
 * @param[in] nao AO数
 * @param[in] D[] 密度行列（圧縮"U"形式）
 *
 * @param[out] G[] 二電子ハミルトン行列（G行列、圧縮"U"形式）
 *
 * @retval  0 正常終了（すべての積分が保存されても、バッファサイズの不足で
 *     保存されていない積分があっても、正常終了である）
 * @retval -1 異常終了（2011/0613現在では考えていない）
 *
 * */
int ofmo_integ_gen_gmat(
	// parallelization
	const int nworkers, const int workerid,
	// basis set & cutoff table data
	const int maxlqn, const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
        const int leading_cs[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// density matrix data & G-matrix (output)
	const int nao, const double D[], double G[] ) {
	// DFT flag
	//const double cHFx ) {
    int nao2;
    double *D_SQ=NULL;
    int mythread = 0;
    int La, Lb, Lc, Ld, Lab, Lcd, Labcd;
    int last_eri_type, last_ijcs, last_klcs;
    int nnao;
    double *Gtmp;
    long nzeri, max_nzeri;
    short *etmp_ind4;
    double *etmp_val;
    // added
    int local_id;
    size_t offset;

    //DFT_Bgn
    //x_coef = cHFx;
    x_coef = 1.0;
    //DFT_End
    // debug
    char *CS = "spdfg";
    int g_last_eri_type = ofmo_twoint_get_global_last_eri_type();
#pragma omp critical
    D_SQ = ofmo_twoint_alloc_square_density( nao );
#pragma omp single
    {
	unpack_matrix( nao, D, D_SQ );
	nao2 = nao*(nao+1)/2;
	memset( G, '\0', sizeof(double)*nao2 );
    }

    float *Dcs;
    Dcs = ofmo_twoint_gen_Dcs(maxlqn, nao, leading_cs, D);
#ifdef USE_CUDA
#pragma omp master
    {
      int ret = 0;
      int ncs = leading_cs[maxlqn+1];
      ret = cuda_genGmat_Init(ncs, nao, Dcs, D_SQ, x_coef);
      if (ret<0) exit(1);
    }
#endif

    mythread = omp_get_thread_num();
#pragma omp barrier

    Gtmp = ofmo_twoint_alloc_local_gmat( mythread, nao );
    last_eri_type = ofmo_twoint_get_last_eri_type( mythread );
    nnao = nao*nao;
    memset( Gtmp, '\0', sizeof(double)*nnao );
    //
    etmp_ind4 = ofmo_twoint_getadd_integ_ind4( mythread );
    etmp_val  = ofmo_twoint_getadd_integ_val( mythread );
    max_nzeri = (long)ofmo_twoint_get_max_stored_integ( mythread );
    ofmo_twoint_set_stored_integ( mythread, 0 );
    // added for load-balancing
    //offset   = ofmo_integ_get_loop_offset( mythread );
    offset = 0;

#ifndef USE_CUDA
    nzeri = 0;
    for ( La=0; La<=maxlqn; La++ ) {
	for ( Lb=0; Lb<=La; Lb++ ) {
	    Lab = La*(La+1)/2 + Lb;
	    for ( Lc=0; Lc<=La; Lc++ ) {
		for ( Ld=0; Ld<=(Lc==La? Lb : Lc ); Ld++ ) {
		    Lcd = Lc*(Lc+1)/2 + Ld;
		    Labcd = Lab*(Lab+1)/2 + Lcd;
						/*// debug*/
			/*printf("#D thd=%d, (%c%c|%c%c) ijcs=%d, klcs=%d\n",
				mythread, CS[La], CS[Lb], CS[Lc], CS[Ld],
				last_ijcs, last_klcs );
			fflush(stdout);*/
		    //if ( Labcd < last_eri_type ) continue;
		    if ( Labcd < last_eri_type ) {
                      offset+=(leading_cs_pair[Lab+1]-leading_cs_pair[Lab]);
                      continue;
                    }
		    if ( Labcd == last_eri_type ) {
			last_ijcs =
			    ofmo_twoint_get_last_ijcs( mythread );
			last_klcs =
			    ofmo_twoint_get_last_klcs( mythread );
		    } else {
			last_ijcs = last_klcs = -1;
		    }
		    local_id =
			(int)((offset+(size_t)workerid)%(size_t)nworkers);
                    start_w2e();
		    calc_twoint_direct[Labcd](
				&nworkers, &local_id,
				//&nworkers, &workerid,
				&La, &Lb, &Lc, &Ld,
				shel_atm, shel_ini,
				atom_x, atom_y, atom_z,
				leading_cs_pair,
				csp_schwarz, csp_ics, csp_jcs,
				csp_leading_ps_pair,
				psp_zeta, psp_dkps, psp_xiza,
				&max_nzeri, &nzeri,
				etmp_val, etmp_ind4,
				&last_ijcs, &last_klcs,
				&nao, D_SQ, Gtmp );
                    set_w2e(Labcd);
		    // added for load-balancing
		    offset+=(leading_cs_pair[Lab+1]-leading_cs_pair[Lab]);
		}// Ld
	    }    // Lc
	}        // Lb
    }	         // La
    if ( nzeri > 0 )
	ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4, D_SQ, Gtmp );
    start_w2e();
    ofmo_twoint_fock_incore_partial( mythread, nao, D_SQ, Gtmp );
    set_w2e(-1);
#else /* USE_CUDA */
    nzeri = 0;
    // idev=1 for GPU, 0 for CPU
    for (int idev=1; idev>=0; idev--) {
    offset = 0;
    for ( La=0; La<=maxlqn; La++ ) {
      for ( Lb=0; Lb<=La; Lb++ ) {
          Lab = La*(La+1)/2 + Lb;
          for ( Lc=0; Lc<=La; Lc++ ) {
              for ( Ld=0; Ld<=(Lc==La? Lb : Lc ); Ld++ ) {
                  Lcd = Lc*(Lc+1)/2 + Ld;
                  Labcd = Lab*(Lab+1)/2 + Lcd;
                  int gpu = (cuda_use_Device() && cuda_get_optCPU(Labcd)!=0);
                  //if ( Labcd < last_eri_type ) continue;
                  if (!gpu && Labcd < last_eri_type ) {
                    offset+=(leading_cs_pair[Lab+1]-leading_cs_pair[Lab]);
                    continue;
                  }
                  if ( Labcd == last_eri_type ) {
                      last_ijcs =
                          ofmo_twoint_get_last_ijcs( mythread );
                      last_klcs =
                          ofmo_twoint_get_last_klcs( mythread );
                  } else {
                      last_ijcs = last_klcs = -1;
                  }
		  local_id =
			(int)((offset+(size_t)workerid)%(size_t)nworkers);
                  start_w2e();
                    //if (idev==0 && Labcd<=g_last_eri_type) {
                  if (idev==0 && !gpu) {
                  calc_twoint_direct[Labcd](
                          &nworkers, &local_id,
                          //&nworkers, &workerid,
                          &La, &Lb, &Lc, &Ld,
                          shel_atm, shel_ini,
                          atom_x, atom_y, atom_z,
                          leading_cs_pair,
                          csp_schwarz, csp_ics, csp_jcs,
                          csp_leading_ps_pair,
                          psp_zeta, psp_dkps, psp_xiza,
                          &max_nzeri, &nzeri,
                          etmp_val, etmp_ind4,
                          &last_ijcs, &last_klcs,
                          &nao, D_SQ, Gtmp );
                    //} else if (idev>0 && Labcd>g_last_eri_type) {
                  } else if (idev>0 && gpu) {
                  cuda_calc_twoint_direct(Labcd,
                          nworkers, local_id,
                          //nworkers, workerid,
                          La, Lb, Lc, Ld,
                          shel_atm, shel_ini,
                          atom_x, atom_y, atom_z,
                          leading_cs_pair,
                          csp_schwarz, csp_ics, csp_jcs,
                          csp_leading_ps_pair,
                          psp_zeta, psp_dkps, psp_xiza,
                          max_nzeri, &nzeri,
                          etmp_val, etmp_ind4,
                          last_ijcs, last_klcs,
                          nao, D_SQ, Gtmp );
                  }
                  set_w2e(Labcd);
		  // added for load-balancing
		  offset+=(leading_cs_pair[Lab+1]-leading_cs_pair[Lab]);
              }// Ld
          }    // Lc
      }        // Lb
    }          // La
    } // idev
    if ( nzeri > 0 )
      ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4, D_SQ, Gtmp );
    start_w2e();
    ofmo_twoint_fock_incore_partial( mythread, nao, D_SQ, Gtmp );
    set_w2e(-1);

#pragma omp master
    {
      int ret;
      ret = cuda_genGmat_Add(nao, Gtmp);
      if (ret<0) exit(1);
    }
#endif /* USE_CUDA */

    // Add the result of the calculated thread to G sequentially
#pragma omp critical
    {
	int i, j, ij;
	ij = 0;
	for ( i=0; i<nao; i++ ) {
	    for ( j=0; j<i; j++ ) {
		G[ij] += Gtmp[i*nao+j] + Gtmp[j*nao+i];
		ij++;
	    }
	    G[ij] += Gtmp[i*nao+i];
	    ij++;
	}
    }
    return 0;
}

//static int (*calc_twoint_direct[])(
static int (*calc_ifc4c[])(
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc, const int *pLd,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) = {
    // original
    /*ofmo_ifc4c_ssss__, ofmo_ifc4c_ssps__, ofmo_ifc4c_sspp__,
    ofmo_ifc4c_ssds__, ofmo_ifc4c_ssdp__, ofmo_ifc4c_ssdd__,
    ofmo_ifc4c_psss__, ofmo_ifc4c_psps__, ofmo_ifc4c_pspp__,
    ofmo_ifc4c_psds__, ofmo_ifc4c_psdp__, ofmo_ifc4c_psdd__,
    ofmo_ifc4c_ppss__, ofmo_ifc4c_ppps__, ofmo_ifc4c_pppp__,
    ofmo_ifc4c_ppds__, ofmo_ifc4c_ppdp__, ofmo_ifc4c_ppdd__,
    ofmo_ifc4c_dsss__, ofmo_ifc4c_dsps__, ofmo_ifc4c_dspp__,
    ofmo_ifc4c_dsds__, ofmo_ifc4c_dsdp__, ofmo_ifc4c_dsdd__,
    ofmo_ifc4c_dpss__, ofmo_ifc4c_dpps__, ofmo_ifc4c_dppp__,
    ofmo_ifc4c_dpds__, ofmo_ifc4c_dpdp__, ofmo_ifc4c_dpdd__,
    ofmo_ifc4c_ddss__, ofmo_ifc4c_ddps__, ofmo_ifc4c_ddpp__,
    ofmo_ifc4c_ddds__, ofmo_ifc4c_dddp__, ofmo_ifc4c_dddd__,*/
    // OS
    ofmo_ifc4c_os_ssss, ofmo_ifc4c_os_ssps, ofmo_ifc4c_os_sspp,
    ofmo_ifc4c_os_ssds, ofmo_ifc4c_os_ssdp, ofmo_ifc4c_os_ssdd,
    ofmo_ifc4c_os_psss, ofmo_ifc4c_os_psps, ofmo_ifc4c_os_pspp,
    ofmo_ifc4c_os_psds, ofmo_ifc4c_os_psdp, ofmo_ifc4c_os_psdd,
    ofmo_ifc4c_os_ppss, ofmo_ifc4c_os_ppps, ofmo_ifc4c_os_pppp,
    ofmo_ifc4c_os_ppds, ofmo_ifc4c_os_ppdp, ofmo_ifc4c_os_ppdd,
    ofmo_ifc4c_os_dsss, ofmo_ifc4c_os_dsps, ofmo_ifc4c_os_dspp,
    ofmo_ifc4c_os_dsds, ofmo_ifc4c_os_dsdp, ofmo_ifc4c_os_dsdd,
    ofmo_ifc4c_os_dpss, ofmo_ifc4c_os_dpps, ofmo_ifc4c_os_dppp,
    ofmo_ifc4c_os_dpds, ofmo_ifc4c_os_dpdp, ofmo_ifc4c_os_dpdd,
    ofmo_ifc4c_os_ddss, ofmo_ifc4c_os_ddps, ofmo_ifc4c_os_ddpp,
    ofmo_ifc4c_os_ddds, ofmo_ifc4c_os_dddp, ofmo_ifc4c_os_dddd,
    // Rys
    /*ofmo_ifc4c_rys_ssss, ofmo_ifc4c_rys_ssps, ofmo_ifc4c_rys_sspp,
    ofmo_ifc4c_rys_ssds, ofmo_ifc4c_rys_ssdp, ofmo_ifc4c_rys_ssdd,
    ofmo_ifc4c_rys_psss, ofmo_ifc4c_rys_psps, ofmo_ifc4c_rys_pspp,
    ofmo_ifc4c_rys_psds, ofmo_ifc4c_rys_psdp, ofmo_ifc4c_rys_psdd,
    ofmo_ifc4c_rys_ppss, ofmo_ifc4c_rys_ppps, ofmo_ifc4c_rys_pppp,
    ofmo_ifc4c_rys_ppds, ofmo_ifc4c_rys_ppdp, ofmo_ifc4c_rys_ppdd,
    ofmo_ifc4c_rys_dsss, ofmo_ifc4c_rys_dsps, ofmo_ifc4c_rys_dspp,
    ofmo_ifc4c_rys_dsds, ofmo_ifc4c_rys_dsdp, ofmo_ifc4c_rys_dsdd,
    ofmo_ifc4c_rys_dpss, ofmo_ifc4c_rys_dpps, ofmo_ifc4c_rys_dppp,
    ofmo_ifc4c_rys_dpds, ofmo_ifc4c_rys_dpdp, ofmo_ifc4c_rys_dpdd,
    ofmo_ifc4c_rys_ddss, ofmo_ifc4c_rys_ddps, ofmo_ifc4c_rys_ddpp,
    ofmo_ifc4c_rys_ddds, ofmo_ifc4c_rys_dddp, ofmo_ifc4c_rys_dddd,*/
};

/** A function that calculates the 4-center Coulomb interaction term
 * @ingroup integ-top
 *
 * Calculate the 4-center Coulomb interaction term between the two monomers
 * that appears in the FMO calculation.
 * The 4-center Coulomb integral is calculated to obtain the 4-center Coulomb
 * interaction term based on the integral and the given density matrix.
 *
 * @attention
 * @li Thread parallelization is performed using OpenMP. At the time of thread parallelism,
 * 	   it is necessary to call this function in the thread parallel area.
 * @li \c nworkers If you set the worker ID \c workerid appropriately, hybrid parallel
 * 	   execution of OpenMP and MPI is possible. In order to obtain a complete Coulomb term
 *     during MPI parallelism, it is necessary to perform reduction processing using the
 *     \c MPI_Allreduce function after the end of this function.
 * @li The resulting Coulomb term \c V_frg[] is sorted by the magnitude of the orbital
 *     quantum number. If you want the Coulomb term in the original sequence, you need to
 *     reorder the elements.
 *
 * @param[in] nworkers Number of worker processes (threads) used for calculation
 * @param[in] workerid Worker ID.
 *     \f$ 0\le \tt{workerid} < \tt{nworkers} \f$
 *
 * @param[in] maxlqn Maximum orbital quantum number
 * @param[in] shel_atm_frg[ics] CS number of the target fragment
 *     The number of the atom to which \c ics of cics belongs
 * @param[in] shel_ini_frg[ics] CS number of the target fragment
 *     The first AO number of the AO contained in the CS of \c ics
 * @param[in] atom_x_frg[iat] X-coordinate (au unit) of atom number \c iat of the target fragment
 * @param[in] atom_y_frg[iat] Y-coordinate (au unit) of atom number \c iat of the target fragment
 * @param[in] atom_z_frg[iat] Z-coordinate (au unit) of atom number \c iat of the target fragment
 * @param[in] leading_cs_pair_frg[itype] The first CS pair number of the CS pair type number
 *     \c itype of the target fragment
 * @param[in] csp_schwarz_frg[icsp] Schwarz integral of CS pair number \c icsp of target fragment
 * @param[in] csp_ics_frg[icsp] The first CS number of the CS pair number \c icsp of the target fragment
 * @param[in] csp_jcs_frg[icsp] The second CS number of the target fragment, CS pair number
 *     \c icsp. However, it is \ f $ \ tt {csp \ _ics [icsp]} \ ge \ tt {csp \ _jcs [icsp]} \ f $.
 * @param[in] csp_leading_ps_pair_frg[icsp] The first PS pair number of the PS pair included
 *     in the CS pair number \c icsp of the target fragment
 * @param[in] psp_zeta_frg[ipsp] 対象フラグメントの、
 *     PSペア番号 \c ipsp の軌道指数和
 *     \f$ \zeta = \zeta_a + \zeta_b \f$
 * @param[in] psp_dkps_frg[ipsp] 対象フラグメントの、
 *     PSペア番号 \c ipsp の線型結合定数
 *     \f[ K_{ab} = \sqrt2 \pi^{5/4} \frac1{\zeta_a+\zeta_b}
 *     \exp\left[ -\frac{\zeta_a \zeta_b}{\zeta_a + \zeta_b}
 *     ( \boldmath A \unboldmath - \boldmath B \unboldmath )^2
 *     \right]\f]
 * @param[in] psp_xiza_frg[ipsp] 対象フラグメントの、
 *     PSペア番号 \c ipsp の、
 *     \f$ \frac{\xi}{\zeta_a} = \frac{\zeta_b}{\zeta_a+\zeta_b} \f$
 *
 * @param[in] shel_atm_mon[ics] Of the partner monomer
 *     CS番号 \c ics のCSが属する原子の番号
 * @param[in] shel_ini_mon[ics] 相手モノマーの、
 *     CS番号 \c ics のCSに含まれるAOの先頭AO番号
 * @param[in] atom_x_mon[iat] X-coordinate (au unit) of atom number \c iat of the mating monomer
 * @param[in] atom_y_mon[iat] 相手モノマーの、
 *     原子の番号 \c iat のy座標（au単位）
 * @param[in] atom_z_mon[iat] 相手モノマーの、
 *     原子の番号 \c iat のz座標（au単位）
 * @param[in] leading_cs_pair_mon[itype] 相手モノマーの、
 *     CSペアタイプ番号 \c itype の先頭CSペア番号
 * @param[in] csp_schwarz_mon[icsp] 相手モノマーの、
 *     CSペア番号 \c icsp のSchwarz積分
 * @param[in] csp_ics_mon[icsp] 相手モノマーの、
 * CSペア番号 \c icsp の1つ目のCS番号
 * @param[in] csp_jcs_mon[icsp] 相手モノマーの、
 *     CSペア番号 \c icsp の2つめのCS番号。ただし、
 *     \f$ \tt{csp\_ics[icsp]} \ge \tt{csp\_jcs[icsp]} \f$ である。
 * @param[in] csp_leading_ps_pair_mon[icsp] ターゲットフラグメントの、
 *     CSペア番号 \c icsp に含まれるPSペアの先頭PSペア番号
 * @param[in] psp_zeta_mon[ipsp] 相手モノマーの、
 *     PSペア番号 \c ipsp の軌道指数和
 *     \f$ \zeta = \zeta_a + \zeta_b \f$
 * @param[in] psp_dkps_mon[ipsp] 相手モノマーの、
 *     PSペア番号 \c ipsp の線型結合定数
 *     \f[ K_{ab} = \sqrt2 \pi^{5/4} \frac1{\zeta_a+\zeta_b}
 *     \exp\left[ -\frac{\zeta_a \zeta_b}{\zeta_a + \zeta_b}
 *     ( \boldmath A \unboldmath - \boldmath B \unboldmath )^2
 *     \right]\f]
 * @param[in] psp_xiza_mon[ipsp] 相手モノマーの、
 *     PSペア番号 \c ipsp の
 *     \f$ \frac{\xi}{\zeta_a} = \frac{\zeta_b}{\zeta_a+\zeta_b} \f$
 *
 * @param[in] nao_mon Number of AOs of partner monomer
 * @param[in] D_mon[] Density matrix of mating monomer (compressed "U" format)
 *
 * @param[out] V_frg[] 4-centric Coulomb interaction term(G matrix, compressed "U" form)
 * with the mating monomer in the fragment of interest.
 * This array needs to give a separate area for each thread.
 * 
 * @retval  0 正常終了（すべての積分が保存されても、バッファサイズの不足で
 *     保存されていない積分があっても、正常終了である）
 * @retval -1 異常終了（2011/06/14現在では考えていない）
 *
 * */
int ofmo_integ_ifc4c_sorted_partial(
	// parallelization
	const int nworkers, const int workerid,
	// basis and cutoff table data for fragment
	const int maxlqn,
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis and cutoff table data for monomer
	const int shel_atm_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[], const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
        //
        const int leading_cs_mon[],
	// density matrix of monomer
	const int nao_mon, const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
    int La, Lb, Lc, Ld, Lab, Lcd, Labcd;
    // added for load-balancing
    int offset=0, local_id=workerid, mythread;
    //
    int type, pos;
    int ma[] = {0, 1, 1, 2, 2, 2,};
    int mb[] = {0, 0, 1, 0, 1, 2,};
    // added for load-balancing
    mythread = omp_get_thread_num();
    //type     = ofmo_integ_get_target_type( mythread );
    type = -1;	// 動的負荷分散しない

#pragma omp master
    TLOG_LOG_IN(5);

    float *Dcs;
    Dcs = ofmo_twoint_gen_Dcs(maxlqn, nao_mon, leading_cs_mon, D_mon);
#ifdef USE_CUDA
#pragma omp master
    {
      int ret = 0;
      int ncs_mon = leading_cs_mon[maxlqn+1];
      ret = cuda_ifc4c_SetDcs(ncs_mon, Dcs);
      if (ret<0) exit(1);
      ret = cuda_ifc4c_calc_Init();
      if (ret<0) exit(1);
    }
#endif

#ifndef USE_CUDA
    offset   = ofmo_integ_get_loop_offset( mythread );
    offset = 0;
    local_id = (offset+workerid)%nworkers;
    for ( La=0; La<=maxlqn; La++ ) {
	for ( Lb=0; Lb<=La; Lb++ ) {
	    Lab = La*(La+1)/2 + Lb;
	    for ( Lc=0; Lc<=maxlqn; Lc++ ) {
		for ( Ld=0; Ld<=Lc; Ld++ ) {
		    Lcd = Lc*(Lc+1)/2 + Ld;
		    //Labcd = Lab*maxlqn2 + Lcd;
		    Labcd = Lab*6 + Lcd; // Here 6 = {1s, 2s, 2p, 3s, 3p, 3d}
		    calc_ifc4c[Labcd](
			    &nworkers, &local_id,
			    &La, &Lb, &Lc, &Ld,
			    shel_atm_frg, shel_ini_frg,
			    atom_x_frg, atom_y_frg, atom_z_frg,
			    leading_cs_pair_frg,
			    csp_schwarz_frg, csp_ics_frg, csp_jcs_frg,
			    csp_leading_ps_pair_frg,
			    psp_zeta_frg, psp_dkps_frg, psp_xiza_frg,
			    shel_atm_mon, shel_ini_mon,
			    atom_x_mon, atom_y_mon, atom_z_mon,
			    leading_cs_pair_mon,
			    csp_schwarz_mon, csp_ics_mon, csp_jcs_mon,
			    csp_leading_ps_pair_mon,
			    psp_zeta_mon, psp_dkps_mon, psp_xiza_mon,
			    D_mon, V_frg );
		    // added for load-balancing
		    offset +=
			(leading_cs_pair_frg[Lab+1]-leading_cs_pair_frg[Lab]);
		    local_id = (offset+workerid)%nworkers;

		}
	    }
	}
    }
#else /* USE_CUDA */
    offset   = ofmo_integ_get_loop_offset( mythread );
    offset = 0;
    local_id = (offset+workerid)%nworkers;
    // idev=1 for GPU, 0 for CPU
    for (int idev=1; idev>=0; idev--) {
    for ( La=0; La<=maxlqn; La++ ) {
	for ( Lb=0; Lb<=La; Lb++ ) {
	    Lab = La*(La+1)/2 + Lb;
	    for ( Lc=0; Lc<=maxlqn; Lc++ ) {
		for ( Ld=0; Ld<=Lc; Ld++ ) {
		    Lcd = Lc*(Lc+1)/2 + Ld;
		    //Labcd = Lab*maxlqn2 + Lcd;
		    Labcd = Lab*6 + Lcd;
		    cuda_ifc4c_calc(idev,
			    &nworkers, &local_id,
			    &La, &Lb, &Lc, &Ld,
			    shel_atm_frg, shel_ini_frg,
			    atom_x_frg, atom_y_frg, atom_z_frg,
			    leading_cs_pair_frg,
			    csp_schwarz_frg, csp_ics_frg, csp_jcs_frg,
			    csp_leading_ps_pair_frg,
			    psp_zeta_frg, psp_dkps_frg, psp_xiza_frg,
			    shel_atm_mon, shel_ini_mon,
			    atom_x_mon, atom_y_mon, atom_z_mon,
			    leading_cs_pair_mon,
			    csp_schwarz_mon, csp_ics_mon, csp_jcs_mon,
			    csp_leading_ps_pair_mon,
			    psp_zeta_mon, psp_dkps_mon, psp_xiza_mon,
			    D_mon, V_frg );
		    // added for load-balancing
		    offset +=
			(leading_cs_pair_frg[Lab+1]-leading_cs_pair_frg[Lab]);
		    local_id = (offset+workerid)%nworkers;

		}
	    }
	}
    }
    }
#endif /* USE_CUDA */
    // added for load-balancing
    ofmo_integ_set_loop_offset( mythread, offset );
#pragma omp master
    TLOG_LOG_OUT(5);
    return 0;
}


static int (*calc_ifc3c[]) (
	// parallelization
	const int *pnworkers, const int *pworkerid,
	// integral type data
	const int *pLa, const int *pLb, const int *pLc,
	// basis and cutoff table data for fragment
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] ) = {
    // OS
    ofmo_ifc3c_os_ssss, ofmo_ifc3c_os_sspp, ofmo_ifc3c_os_ssdd,
    ofmo_ifc3c_os_psss, ofmo_ifc3c_os_pspp, ofmo_ifc3c_os_psdd,
    ofmo_ifc3c_os_ppss, ofmo_ifc3c_os_pppp, ofmo_ifc3c_os_ppdd,
    ofmo_ifc3c_os_dsss, ofmo_ifc3c_os_dspp, ofmo_ifc3c_os_dsdd,
    ofmo_ifc3c_os_dpss, ofmo_ifc3c_os_dppp, ofmo_ifc3c_os_dpdd,
    ofmo_ifc3c_os_ddss, ofmo_ifc3c_os_ddpp, ofmo_ifc3c_os_dddd,
    // Rys
    /*ofmo_ifc3c_rys_ssss, ofmo_ifc3c_rys_sspp, ofmo_ifc3c_rys_ssdd,
    ofmo_ifc3c_rys_psss, ofmo_ifc3c_rys_pspp, ofmo_ifc3c_rys_psdd,
    ofmo_ifc3c_rys_ppss, ofmo_ifc3c_rys_pppp, ofmo_ifc3c_rys_ppdd,
    ofmo_ifc3c_rys_dsss, ofmo_ifc3c_rys_dspp, ofmo_ifc3c_rys_dsdd,
    ofmo_ifc3c_rys_dpss, ofmo_ifc3c_rys_dppp, ofmo_ifc3c_rys_dpdd,
    ofmo_ifc3c_rys_ddss, ofmo_ifc3c_rys_ddpp, ofmo_ifc3c_rys_dddd,*/
};

/** A function that calculates the 3-center Coulomb interaction term
 * @ingroup integ-top
 *
 * FMO計算に現れる、２つのモノマー間の３中心クーロン相互作用項を
 * 計算する。
 * ３中心クーロン積分を計算して、積分と与えられた密度行列を元に、
 * ３中心クーロン相互作用項を求める。
 *
 * @attention
 * @li この関数はOpenMPを用いたスレッド並列化が行われている。スレッド並列
 *     実行のためには、スレッド並列領域内からこの関数を呼び出す
 *     必要がある。
 * @li \c nworkers と \c workerid を適切に設定すると、OpenMPとMPIの
 *     ハイブリッド並列実行が可能である。MPI並列を利用する際には、
 *     関数終了後に、\c MPI_Allreduce 関数などを用いたリダクション処理
 *     を行うことで、完全なクーロン項が得られる。
 * @li 得られるクーロン項 \c V_frg[] は、軌道量子数の大きさで
 *     ソートされたものである。元の並びのクーロン項が欲しい場合には、
 *     要素の並べ替えが必要である。
 *
 * @param[in] nworkers 計算に用いるワーカプロセス（スレッド）数
 * @param[in] workerid ワーカID。
 *     \f$ 0\le \tt{workerid} < \tt{nworkers} \f$
 *
 * @param[in] maxlqn 最大軌道量子数
 * @param[in] shel_atm_frg[ics] 対象フラグメントの、
 *     CS番号 \c ics のCSが属する原子の番号
 * @param[in] shel_ini_frg[ics] 対象フラグメントの、
 *     CS番号 \c ics のCSに含まれるAOの先頭AO番号
 * @param[in] atom_x_frg[iat] 対象フラグメントの、
 *     原子の番号 \c iat のx座標（au単位）
 * @param[in] atom_y_frg[iat] 対象フラグメントの、
 *     原子の番号 \c iat のy座標（au単位）
 * @param[in] atom_z_frg[iat] 対象フラグメントの、
 *     原子の番号 \c iat のz座標（au単位）
 * @param[in] leading_cs_pair_frg[itype] 対象フラグメントの、
 *     CSペアタイプ番号 \c itype の先頭CSペア番号
 * @param[in] csp_schwarz_frg[icsp] 対象フラグメントの、
 *     CSペア番号 \c icsp のSchwarz積分
 * @param[in] csp_ics_frg[icsp] 対象フラグメントの、
 * CSペア番号 \c icsp の1つ目のCS番号
 * @param[in] csp_jcs_frg[icsp] 対象フラグメントの、
 *     CSペア番号 \c icsp の2つめのCS番号。ただし、
 *     \f$ \tt{csp\_ics[icsp]} \ge \tt{csp\_jcs[icsp]} \f$ である。
 * @param[in] csp_leading_ps_pair_frg[icsp] 対象フラグメントの、
 *     CSペア番号 \c icsp に含まれるPSペアの先頭PSペア番号
 * @param[in] psp_zeta_frg[ipsp] 対象フラグメントの、
 *     PSペア番号 \c ipsp の軌道指数和
 *     \f$ \zeta = \zeta_a + \zeta_b \f$
 * @param[in] psp_dkps_frg[ipsp] 対象フラグメントの、
 *     PSペア番号 \c ipsp の線型結合定数
 *     \f[ K_{ab} = \sqrt2 \pi^{5/4} \frac1{\zeta_a+\zeta_b}
 *     \exp\left[ -\frac{\zeta_a \zeta_b}{\zeta_a + \zeta_b}
 *     ( \boldmath A \unboldmath - \boldmath B \unboldmath )^2
 *     \right]\f]
 * @param[in] psp_xiza_frg[ipsp] 対象フラグメントの、
 *     PSペア番号 \c ipsp の
 *     \f$ \frac{\xi}{\zeta_a} = \frac{\zeta_b}{\zeta_a+\zeta_b} \f$
 *
 * @param[in] leading_cs_mon[lqn] 相手モノマーの、
 *     軌道量子数 \c lqn の先頭CS番号
 * @param[in] shel_tem_mon[ics] 相手モノマーの、CS番号 \c ics のCSの縮約長
 * @param[in] shel_atm_mon[ics] 相手モノマーの、
 *     CS番号 \c ics のCSが属する原子の番号
 * @param[in] shel_add_mon[ics] 相手モノマーの、CS番号 \c ics のCSに属する
 *     PSの先頭PS番号
 * @param[in] shel_ini_mon[ics] 相手モノマーの、CS番号 \c ics のCSに
 *     含まれるAOの先頭AO番号
 * @param[in] atom_x_mon[iat] 相手モノマーの、
 *     原子の番号 \c iat のx座標（au単位）
 * @param[in] atom_y_mon[iat] 相手モノマーの、
 *     原子の番号 \c iat のy座標（au単位）
 * @param[in] atom_z_mon[iat] 相手モノマーの、
 *     原子の番号 \c iat のz座標（au単位）
 * @param[in] prim_exp_mon[ips] 相手モノマーの、PS番号 \c ips のPSの
 *     軌道指数
 * @param[in] prim_coe_mon[ips] 相手モノマーの、PS番号 \c ips のPSの
 *     規格化定数込みの縮約係数
 *
 * @param[in] ao_pop_mon[] AO population of the partner monomer
 *
 * @param[out] V_frg[] 対象フラグメントにおける
 *     相手モノマーとの間の３中心クーロン相互作用項
 *     （G行列、圧縮"U"形式）。この配列は、スレッドごとに別領域
 *     を与える必要がある。
 * 
 * @retval  0 正常終了（すべての積分が保存されても、バッファサイズの不足で
 *     保存されていない積分があっても、正常終了である）
 * @retval -1 異常終了（2011/06/14現在では考えていない）
 *
 * */
int ofmo_integ_ifc3c_sorted_partial(
	// parallelization
	const int nworkers, const int workerid,
	// basis and cutoff table data for fragment
	const int maxlqn,
	const int shel_atm_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[], const int leading_cs_pair_frg[],
	//const double csp_schwarz_frg[],
	const int csp_ics_frg[], const int csp_jcs_frg[],
	const int csp_leading_ps_pair_frg[],
	const double psp_zeta_frg[], const double psp_dkps_frg[],
	const double psp_xiza_frg[],
	// basis set data for monomer
	const int leading_cs_mon[],
	const int shel_tem_mon[], const int shel_atm_mon[],
	const int shel_add_mon[], const int shel_ini_mon[],
	const double atom_x_mon[], const double atom_y_mon[],
	const double atom_z_mon[],
	const double prim_exp_mon[], const double prim_coe_mon[],
	// monomer AO population
	const double ao_pop_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
    int La, Lb, Lc, Lab, Labc;
    // added for load-balancing
    int offset, local_id=workerid, mythread;
    //
    int type, pos;
    int ma[] = {0, 1, 1, 2, 2, 2,};
    int mb[] = {0, 0, 1, 0, 1, 2,};
    // added for load-balancing
    mythread = omp_get_thread_num();
    //type = ofmo_integ_get_target_type( mythread );	// 動的負荷分散しない
    type = -1;

    offset   = ofmo_integ_get_loop_offset( mythread );
    local_id = (offset+workerid)%nworkers;
    for ( La=0; La<=maxlqn; La++ ) {
	for ( Lb=0; Lb<=La; Lb++ ) {
	    Lab = La*(La+1)/2 + Lb;
	    for ( Lc=0; Lc<=maxlqn; Lc++ ) {
		Labc = 3*Lab + Lc;
		calc_ifc3c[Labc](
			&nworkers, &local_id,
			&La, &Lb, &Lc,
			shel_atm_frg, shel_ini_frg,
			atom_x_frg, atom_y_frg, atom_z_frg,
			leading_cs_pair_frg,
			csp_ics_frg, csp_jcs_frg,
			csp_leading_ps_pair_frg,
			psp_zeta_frg, psp_dkps_frg, psp_xiza_frg,
			leading_cs_mon,
			shel_tem_mon, shel_atm_mon, shel_add_mon,
			shel_ini_mon, atom_x_mon, atom_y_mon,
			atom_z_mon,
			prim_exp_mon, prim_coe_mon,
			ao_pop_mon, V_frg );
		// added for load-balancing
		offset += (leading_cs_pair_frg[Lab+1]
			- leading_cs_pair_frg[Lab]);
		local_id = (offset+workerid)%nworkers;
	    }
	}
    }
    ofmo_integ_set_loop_offset( mythread, offset );
    return 0;
}

static int (*calc_ifc2c[]) (
	// parallelization
	const int *nworkers, const int *workerid,
	// integral type data
	const int *pLa, const int *pLb,
	// basis set data for fragment
	const int leading_cs_frg[],
	const int shel_tem_frg[], const int shel_atm_frg[],
	const int shel_add_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[],
	const double prim_exp_frg[], const double prim_coe_frg[],
	// atomic charge and atomic coordinate data for monomer
	const int *nat_mon, const double atom_x_mon[],
	const double atom_y_mon[], const double atom_z_mon[],
	const double atm_pop_mon[],
	// output
	double V_frg[] ) = {
    ofmo_ifc2c_ss__, ofmo_ifc2c_ps__, ofmo_ifc2c_pp__,
    ofmo_ifc2c_ds__, ofmo_ifc2c_dp__, ofmo_ifc2c_dd__,
};

/** ２中心クーロン相互作用項の計算を行う関数
 * @ingroup integ-top
 *
 * FMO計算で現れる２中心クーロン相互作用項を計算する関数。
 * 計算対象のフラグメントのソート基底関数と、相手モノマーのatomic
 * populationを与えると、２中心クーロン作用項が計算される。
 *
 * @attention
 * @li スレッド並列実行を行う場合には、スレッド並列領域内から
 *     この関数を呼び出す必要がある。
 * @li \c nworkers と \c workerid を適切に設定すると、OpenMPとMPIの
 *     ハイブリッド並列実行が可能である。ただし、MPI並列利用時には、
 *     関数終了後に、\c MPI_Allreduce などを用いたリダクション処理を
 *     行うことで、完全なクーロン項が得られる
 * @li 得られるクーロン項 \c V_frg[] は、軌道量子数の大きさで
 *     ソートされたものである。元の並びのクーロン項が欲しい場合には、
 *     要素の並べ替えが必要である。
 *
 * @param[in] nworkers 計算に用いるワーカプロセス（スレッド）数
 * @param[in] workerid 各ワーカプロセス（スレッド）のID。
 *     \f$ 0\le\tt{workerid}<\tt{nworkers} \f$である。
 * @param[in] maxlqn 最大軌道量子数
 * @param[in] leading_cs_frg[lqn] 対象フラグメントの、
 *     軌道量子数 \c lqn の先頭CS番号
 * @param[in] shel_tem_frg[ics] 対象フラグメントの、
 *     CS番号 \c ics のCSの縮約長
 * @param[in] shel_atm_frg[ics] 対象フラグメントの、
 *     CS番号 \c ics のCSが属する原子の番号
 * @param[in] shel_add_frg[ics] 対象フラグメントの、
 *     CS番号 \c ics のCSに含まれるPSの先頭PS番号
 * @param[in] shel_ini_frg[ics] 対象フラグメントの、
 *     CS番号 \c ics のCSに含まれるAOの先頭AO番号
 * @param[in] atom_x_frg[iat] 対象フラグメントの、
 *     原子の番号 \c iat のx座標（au単位）
 * @param[in] atom_y_frg[iat] 対象フラグメントの、
 *     原子の番号 \c iat のy座標（au単位）
 * @param[in] atom_z_frg[iat] 対象フラグメントの、
 *     原子の番号 \c iat のz座標（au単位）
 * @param[in] prim_exp_frg[ips] 対象フラグメントの、
 *     PS番号 \c ips のPSの軌道指数
 * @param[in] prim_coe_frg[ips] 対象フラグメントの、
 *     PS番号 \c ips のPSの規格化定数込みの縮約係数
 *
 * @param[in] nat_mon 相手モノマーの、原子数
 * @param[in] atom_x_mon[iat] 相手モノマーの、
 *     原子の番号 \c iat のx座標（au単位）
 * @param[in] atom_y_mon[iat] 相手モノマーの、
 *     原子の番号 \c iat のy座標（au単位）
 * @param[in] atom_z_mon[iat] 相手モノマーの、
 *     原子の番号 \c iat のz座標（au単位）
 * @param[in] atm_pop_mon[iat] 相手モノマーの、
 *     原子の番号 \c iat の原子番号
 *
 * @param[out] V_frg[] 対象フラグメントにおける
 *     相手モノマーとの間の２中心クーロン相互作用項
 *     （G行列、圧縮"U"形式）。この配列は、同一プロセス内のスレッド
 *     間で共有である。
 *
 * @retval  0 正常終了
 * @retval -1 異常終了（いま(2011/06/13)のところ考えていない）
 *
 * */
int ofmo_integ_ifc2c_sorted_partial(
	// parallelization
	const int nworkers, const int workerid,
	// input data of fragment
	const int maxlqn, const int leading_cs_frg[],
	const int shel_tem_frg[], const int shel_atm_frg[],
	const int shel_add_frg[], const int shel_ini_frg[],
	const double atom_x_frg[], const double atom_y_frg[],
	const double atom_z_frg[],
	const double prim_exp_frg[], const double prim_coe_frg[],
	// input data of counter monomer
	const int nat_mon, const double atom_x_mon[],
	const double atom_y_mon[], const double atom_z_mon[],
	const double atm_pop_mon[],
	// output data
	double V_frg[] ) {
    int La, Lb, Lab;
    int type, mythread;
    int ma[] = { 0, 1, 1, 2, 2, 2, };
    int mb[] = { 0, 0, 1, 0, 1, 2, };
    mythread = omp_get_thread_num();
    //type = ofmo_integ_get_target_type( mythread );
    type = -1;
    for ( La=0; La<=maxlqn; La++ ) {
	for ( Lb=0; Lb<=La; Lb++ ) {
	    Lab = La*(La+1)/2 + Lb;
	    calc_ifc2c[Lab](
		    &nworkers, &workerid,
		    &La, &Lb, leading_cs_frg,
		    shel_tem_frg, shel_atm_frg, shel_add_frg,
		    shel_ini_frg,
		    atom_x_frg, atom_y_frg, atom_z_frg,
		    prim_exp_frg, prim_coe_frg,
		    &nat_mon, atom_x_mon, atom_y_mon, atom_z_mon,
		    atm_pop_mon, V_frg );
	}
    }
    return 0;
}

static void print_eri(const double eri_val[], const short int eri_idx4[],
  const int nzeri){
    int ix, ix4;
    int i,j,k,l;
	double v;
    for(ix=0, ix4=0; ix<nzeri; ix++, ix4+=4){
        i = (int) eri_idx4[ix4+0];
        j = (int) eri_idx4[ix4+1];
        k = (int) eri_idx4[ix4+2];
        l = (int) eri_idx4[ix4+3];
		v = eri_val[ix];
		if(i == j){
			v *= 2;
		}
		if(k == l){
			v *= 2;
		}
		if((i == k && j == l) ||
		   (i == l && j == k)){
			v *= 2;
		}
        printf("%d %d %d %d %.7f\n", i, j, k, l, v);
    }
}

int ofmo_integ_export_eri(
	// parallelization
	const int nworkers, const int workerid,
	// basis set & cutoff table data
	const int maxlqn, const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[], const int leading_cs_pair[],
        const int leading_cs[],
	const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// density matrix data & G-matrix (output)
	const int nao, const double Ct[], double mo_tei[] ) {

	int nao2;
    double *D_SQ=NULL, *Gtmp;
	int mythread=0;
    int La, Lb, Lc, Ld, Lab, Lcd, Labcd;
	int ix, ix4;
	int iao,jao,kao,lao;
	int imo,jmo,kmo,lmo,mo_idx;
	int imo4, jmo3, kmo2;
	int ic, jc, kc, lc;
    int last_eri_type, last_ijcs, last_klcs;
    int nnao, nao_2, nao_3, nao_4;
    long nzeri, max_nzeri;
    short *etmp_ind4;
    double *etmp_val;
    int local_id;
	int ret;
    size_t offset;

	double eval;

	nao2 = nao*(nao+1)/2;
	nao_2 = nao * nao;
	nao_3 = nao_2 * nao;
	nao_4 = nao_3 * nao;
	
	mythread = omp_get_thread_num();
	Gtmp = ofmo_twoint_alloc_local_gmat( mythread, nao );
	nnao = nao*nao;
    memset( Gtmp, '\0', sizeof(double)*nnao );
	etmp_ind4 = ofmo_twoint_getadd_integ_ind4( mythread );
    etmp_val  = ofmo_twoint_getadd_integ_val( mythread );
    max_nzeri = (long)ofmo_twoint_get_max_stored_integ( mythread );
    offset = 0;

	/* TODO : Parallelize here.
	 * short-term goal : use enough nzeris.
	 * long-term goal  : duplicate calc_twoint_directs that doesn't erase buffer.
	 */
	ret = 0;
#pragma omp critical
{
    for ( La=0; La<=maxlqn; La++ ) {
	for ( Lb=0; Lb<=La; Lb++ ) {
	    Lab = La*(La+1)/2 + Lb;
	    for ( Lc=0; Lc<=La; Lc++ ) {
		for ( Ld=0; Ld<=(Lc==La? Lb : Lc ); Ld++ ) {
			if(ret == OFMO_EBUF_FULL){
				// printf("ERROR : Not enough buffer\n");
				break;
			}
		    ofmo_twoint_set_stored_integ( mythread, 0 );
			last_eri_type = -1;
			nzeri = 0;
		    Lcd = Lc*(Lc+1)/2 + Ld;
		    Labcd = Lab*(Lab+1)/2 + Lcd;
						/*// debug*/
			/*printf("#D thd=%d, (%c%c|%c%c) ijcs=%d, klcs=%d\n",
				mythread, CS[La], CS[Lb], CS[Lc], CS[Ld],
				last_ijcs, last_klcs );
			fflush(stdout);*/
		    //if ( Labcd < last_eri_type ) continue;
		    /*if ( Labcd < last_eri_type ) {
                      offset+=(leading_cs_pair[Lab+1]-leading_cs_pair[Lab]);
                      continue;
                    }
		    if ( Labcd == last_eri_type ) {
			last_ijcs =
			    ofmo_twoint_get_last_ijcs( mythread );
			last_klcs =
			    ofmo_twoint_get_last_klcs( mythread );
		    } else {
			last_ijcs = last_klcs = -1;
		    }
			*/
			last_ijcs = last_klcs = -1;
		    local_id =
			(int)((offset+(size_t)workerid)%(size_t)nworkers);
                    start_w2e();
		    ret = calc_twoint_direct[Labcd](
				&nworkers, &local_id,
				//&nworkers, &workerid,
				&La, &Lb, &Lc, &Ld,
				shel_atm, shel_ini,
				atom_x, atom_y, atom_z,
				leading_cs_pair,
				csp_schwarz, csp_ics, csp_jcs,
				csp_leading_ps_pair,
				psp_zeta, psp_dkps, psp_xiza,
				&max_nzeri, &nzeri,
				etmp_val, etmp_ind4,
				&last_ijcs, &last_klcs,
				&nao, D_SQ, Gtmp );
                    set_w2e(Labcd);
		    // added for load-balancing
		    offset+=(leading_cs_pair[Lab+1]-leading_cs_pair[Lab]);
			if(ret == OFMO_EBUF_FULL){
				// printf("ERROR : Not enough buffer\n");
				break;
			}
			else{
				for(ix=0, ix4=0; ix<nzeri; ix++, ix4+=4){
					eval = etmp_val[ix];
					iao = etmp_ind4[ix4+0];
					jao = etmp_ind4[ix4+1];
					kao = etmp_ind4[ix4+2];
					lao = etmp_ind4[ix4+3];
					if(iao == jao) eval *= 2;
					if(kao == lao) eval *= 2;
					if((iao == kao && jao == lao) ||
					   (iao == lao && jao == kao)) eval *= 2;
					printf("AO : %d %d %d %d | %f\n", iao, jao, kao, lao, eval);
					fflush(stdout);
					for(imo=0, imo4=0, ic=iao; imo<nao;   imo++, imo4+=nao_3, ic+=nao){
					for(jmo=0, jmo3=0, jc=jao; jmo<imo+1; jmo++, jmo3+=nao_2, jc+=nao){
					for(kmo=0, kmo2=0, kc=kao; kmo<imo+1; kmo++, kmo2+=nao,   kc+=nao){
					for(lmo=0, 		   lc=lao; lmo<kmo+1; lmo++, 			  lc+=nao){
						/* mo_idx = imo * nao * nao * nao
								+jmo * nao * nao
								+kmo * nao
								+lmo; */
						mo_idx = imo4 + jmo3 + kmo2 + lmo;
						mo_tei[mo_idx] += eval * Ct[ic] * Ct[jc] * Ct[lc] * Ct[kc];
					}}}}
				}
				//print_eri(etmp_val, etmp_ind4, nzeri);
			}
		}// Ld
	    }    // Lc
	}        // Lb
    }	         // La
} // OMP CRITICAL
	if(ret != 0) return -1;
	return 0;
}
