/**
 * @file ofmo-ifc3c.c 同じタイプの３中心クーロン積分を行う関数群
 * 2011/11/02現在、４中心積分用のコードを呼び出している（無駄が多い）
 * */
/**
 * @defgroup integ-ifc3c ３中心クーロン積分を行う関数群
 *
 * 同じタイプの３中心クーロン積分を行い、環境ポテンシャル項に加算する
 * 関数群。
 *
 * 初期化関数（\c ofmo_ifc3c_init）以外は同じ引数をもっているので、
 * 以下にその内容を示す。
 *
 * @param[in] nworkers 計算に用いるワーカプロセス（スレッド）数
 * @param[in] workerid 各ワーカプロセス（スレッド）のID
 *     （\f$ 0\le\tt{workerid}<\tt{nworkers} \f$）
 * @param[in] La １つ目の軌道量子数
 * @param[in] Lb ２つ目の軌道量子数
 *     （\f$ \tt{La} \ge \tt{Lb} \f$）
 * @param[in] Lc ３つ目の軌道量子数
 * @param[in] shel_atm_frg[ics] フラグメントのCS番号 \c ics のCSが属する
 * 原子の番号
 * @param[in] shel_ini_frg[ics] フラグメントのCS番号 \c ics のCSに
 * 含まれるAOの先頭AO番号
 * @param[in] atom_x_frg[iat] フラグメントの原子の番号 \c iat のx座標
 * （au単位）
 * @param[in] atom_y_frg[iat] フラグメントの原子の番号 \c iat のy座標
 * （au単位）
 * @param[in] atom_z_frg[iat] フラグメントの原子の番号 \c iat のz座標
 * （au単位）
 * @param[in] leading_cs_pair_frg[itype] フラグメントのCSペアタイプ番号
 * \c itype の先頭CSペア番号
 * @param[in] csp_ics_frg[icsp] フラグメントのCSペア番号 \c icsp の
 * 1つ目のCS番号
 * @param[in] csp_jcs_frg[icsp] フラグメントのCSペア番号 \c icsp の
 * 2つめのCS番号。ただし、
 * \f$ \tt{csp\_ics[icsp]} \ge \tt{csp\_jcs[icsp]} \f$ である。
 * @param[in] csp_leading_ps_pair_frg[icsp] フラグメントのCSペア番号
 * \c icsp に含まれるPSペアの先頭PSペア番号
 * @param[in] psp_zeta_frg[ipsp] フラグメントのPSペア番号 \c ipsp の
 * 軌道指数和\f$ \zeta = \zeta_a + \zeta_b \f$
 * @param[in] psp_dkps_frg[ipsp] フラグメントのPSペア番号 \c ipsp
 * の線型結合定数
 *     \f[ K_{ab} = \sqrt2 \pi^{5/4} \frac1{\zeta_a+\zeta_b}
 *     \exp\left[ -\frac{\zeta_a \zeta_b}{\zeta_a + \zeta_b}
 *     ( \boldmath A \unboldmath - \boldmath B \unboldmath )^2
 *     \right]\f]
 * @param[in] psp_xiza_frg[ipsp] フラグメントのPSペア番号 \c ipsp の
 *     \f$ \frac{\xi}{\zeta_a} = \frac{\zeta_b}{\zeta_a+\zeta_b} \f$
 *
 * @param[in] leading_cs_mon[lqn] 相手モノマーの
 *     軌道量子数 /c lqn の先頭CS番号
 * @param[in] shel_tem_mon[ics] 相手モノマーの
 *     CS番号 \c ics のCSの縮約長
 * @param[in] shel_atm_mon[ics] 相手モノマーの
 *     CS番号 \c ics のCSが属する原子の番号
 * @param[in] shel_add_mon[ics] 相手モノマーの
 *     CS番号 \c ics のCSに含まれるPSの先頭PS番号
 * @param[in] shel_ini_mon[ics] 相手モノマーの
 *     CS番号 \c ics のCSの先頭AO番号
 * @param[in] atom_x_mon[iat] 相手モノマーの
 *     原子の番号 \c iat のx座標（au単位）
 * @param[in] atom_y_mon[iat] 相手モノマーの
 *     原子の番号 \c iat のy座標（au単位）
 * @param[in] atom_z_mon[iat] 相手モノマーの
 *     原子の番号 \c iat のz座標（au単位）
 * @param[in] prim_exp_mon[ips] 相手モノマーの
 *     PS番号 \c ips のPSの軌道指数
 * @param[in] prim_coe_mon[ips] 相手モノマーの
 *     PS番号 \c ips のPSの規格化定数込みの縮約係数
 * @param[in] ao_pop_mon[iao] 相手モノマーのAO番号\c iao のAO population
 *
 * @param[out] V_frg[] この関数で計算した３中心クーロンポテンシャルが
 * 加算された環境ポテンシャル項（圧縮U形式）
 *
 * @attention
 * 引数のV_frg[]がスレッド毎に異なる領域をさしている場合には、
 * そのままスレッド並列に対応できる。
 *
 * @ingroup integ-med
 * */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAXTEM 10
#define MAXTEM2 MAXTEM*MAXTEM

#define ZERO    0.e0
#define ONE     1.e0
#define TWO     2.e0
#define FOUR    4.e0
#define HALF    0.5e0

#define EPS_ERI 1.e-15
#define EPS_PS4 1.e-30

#include "ofmo-twoint-core.h"

#ifndef false
#define false 0
#endif
#ifndef true
#define true 1
#endif

extern void ofmo_twoint_core_os_ssss(
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *nijps, const double vzeta[], const double vdkab[],
        const double vxiza[], const double BA[3],
        const int *nklps, const double veta[], const double vdkcd[],
        const double vxizc[], const double DC[3], const double AC[3],
        double *DINT );
extern void ofmo_twoint_core_os_psss(
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *nijps, const double vzeta[], const double vdkab[],
        const double vxiza[], const double BA[3],
        const int *nklps, const double veta[], const double vdkcd[],
        const double vxizc[], const double DC[3], const double AC[3],
        double *DINT );
extern void ofmo_twoint_core_os_psps(
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *nijps, const double vzeta[], const double vdkab[],
        const double vxiza[], const double BA[3],
        const int *nklps, const double veta[], const double vdkcd[],
        const double vxizc[], const double DC[3], const double AC[3],
        double *DINT );
extern void ofmo_twoint_core_os_ppss(
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *nijps, const double vzeta[], const double vdkab[],
        const double vxiza[], const double BA[3],
        const int *nklps, const double veta[], const double vdkcd[],
        const double vxizc[], const double DC[3], const double AC[3],
        double *DINT );
extern void ofmo_twoint_core_os_ppps(
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *nijps, const double vzeta[], const double vdkab[],
        const double vxiza[], const double BA[3],
        const int *nklps, const double veta[], const double vdkcd[],
        const double vxizc[], const double DC[3], const double AC[3],
        double *DINT );
extern void ofmo_twoint_core_os_pppp(
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *nijps, const double vzeta[], const double vdkab[],
        const double vxiza[], const double BA[3],
        const int *nklps, const double veta[], const double vdkcd[],
        const double vxizc[], const double DC[3], const double AC[3],
        double *DINT );
extern void ofmo_twoint_core_os_dsss(
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *nijps, const double vzeta[], const double vdkab[],
        const double vxiza[], const double BA[3],
        const int *nklps, const double veta[], const double vdkcd[],
        const double vxizc[], const double DC[3], const double AC[3],
        double *DINT );
extern void ofmo_twoint_core_os_dsps(
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *nijps, const double vzeta[], const double vdkab[],
        const double vxiza[], const double BA[3],
        const int *nklps, const double veta[], const double vdkcd[],
        const double vxizc[], const double DC[3], const double AC[3],
        double *DINT );
extern void ofmo_twoint_core_os_dspp(
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *nijps, const double vzeta[], const double vdkab[],
        const double vxiza[], const double BA[3],
        const int *nklps, const double veta[], const double vdkcd[],
        const double vxizc[], const double DC[3], const double AC[3],
        double *DINT );
extern void ofmo_twoint_core_os_dsds(
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *nijps, const double vzeta[], const double vdkab[],
        const double vxiza[], const double BA[3],
        const int *nklps, const double veta[], const double vdkcd[],
        const double vxizc[], const double DC[3], const double AC[3],
        double *DINT );
extern void ofmo_twoint_core_os_dpss(
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *nijps, const double vzeta[], const double vdkab[],
        const double vxiza[], const double BA[3],
        const int *nklps, const double veta[], const double vdkcd[],
        const double vxizc[], const double DC[3], const double AC[3],
        double *DINT );
extern void ofmo_twoint_core_os_dpps(
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *nijps, const double vzeta[], const double vdkab[],
        const double vxiza[], const double BA[3],
        const int *nklps, const double veta[], const double vdkcd[],
        const double vxizc[], const double DC[3], const double AC[3],
        double *DINT );
extern void ofmo_twoint_core_os_dppp(
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *nijps, const double vzeta[], const double vdkab[],
        const double vxiza[], const double BA[3],
        const int *nklps, const double veta[], const double vdkcd[],
        const double vxizc[], const double DC[3], const double AC[3],
        double *DINT );
extern void ofmo_twoint_core_os_dpds(
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *nijps, const double vzeta[], const double vdkab[],
        const double vxiza[], const double BA[3],
        const int *nklps, const double veta[], const double vdkcd[],
        const double vxizc[], const double DC[3], const double AC[3],
        double *DINT );
extern void ofmo_twoint_core_os_dpdp(
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *nijps, const double vzeta[], const double vdkab[],
        const double vxiza[], const double BA[3],
        const int *nklps, const double veta[], const double vdkcd[],
        const double vxizc[], const double DC[3], const double AC[3],
        double *DINT );
extern void ofmo_twoint_core_os_ddss(
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *nijps, const double vzeta[], const double vdkab[],
        const double vxiza[], const double BA[3],
        const int *nklps, const double veta[], const double vdkcd[],
        const double vxizc[], const double DC[3], const double AC[3],
        double *DINT );
extern void ofmo_twoint_core_os_ddps(
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *nijps, const double vzeta[], const double vdkab[],
        const double vxiza[], const double BA[3],
        const int *nklps, const double veta[], const double vdkcd[],
        const double vxizc[], const double DC[3], const double AC[3],
        double *DINT );
extern void ofmo_twoint_core_os_ddpp(
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *nijps, const double vzeta[], const double vdkab[],
        const double vxiza[], const double BA[3],
        const int *nklps, const double veta[], const double vdkcd[],
        const double vxizc[], const double DC[3], const double AC[3],
        double *DINT );
extern void ofmo_twoint_core_os_ddds(
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *nijps, const double vzeta[], const double vdkab[],
        const double vxiza[], const double BA[3],
        const int *nklps, const double veta[], const double vdkcd[],
        const double vxizc[], const double DC[3], const double AC[3],
        double *DINT );
extern void ofmo_twoint_core_os_dddp(
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *nijps, const double vzeta[], const double vdkab[],
        const double vxiza[], const double BA[3],
        const int *nklps, const double veta[], const double vdkcd[],
        const double vxizc[], const double DC[3], const double AC[3],
        double *DINT );
extern void ofmo_twoint_core_os_dddd(
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *nijps, const double vzeta[], const double vdkab[],
        const double vxiza[], const double BA[3],
        const int *nklps, const double veta[], const double vdkcd[],
        const double vxizc[], const double DC[3], const double AC[3],
        double *DINT );

static double _CK_;
static double vxizc[MAXTEM2];
static double DC[3];

int ofmo_ifc3c_os_init() {
    int i;
    double pi;
    static int called = false;
    if ( called ) return 0;
    pi     = FOUR * atan( ONE );
    _CK_   = sqrt( TWO * sqrt(pi) ) * pi;
    for ( i=0; i<MAXTEM2; i++ ) vxizc[i] = ZERO;
    for ( i=0; i<3; i++ ) DC[i] = ZERO;
    called = true;
    return 0;
}

static int ifc3c_calc_cd_pair_params(
	const double prim_exp[], const double prim_coe[],
	const int kps0, const int kps1,
	double veta[], double vdkcd[] ) {
    int kps, lps, nklps;
    double zeta_c, zeta_d, coef_c, coef_d;
    double sqr_eta, eta;
    for ( kps=kps0, nklps=0; kps<kps1; kps++ ) {
	zeta_c = prim_exp[kps];
	coef_c = prim_coe[kps];
	for ( lps=kps0; lps<kps1; lps++, nklps++ ) {
	    zeta_d = prim_exp[lps];
	    coef_d = prim_coe[lps];

	    sqr_eta      = sqrt( ONE / (zeta_c+zeta_d) );
	    eta          = sqr_eta*sqr_eta;
	    veta[nklps]  = eta;
	    vdkcd[nklps] = _CK_ * coef_c * coef_d * sqr_eta * eta;
	}
    }
    return nklps;
}

/** (ss,ss)タイプの３中心クーロン積分をまとめて行う関数
 * @ingroup integ-ifc3c
 * */
int ofmo_ifc3c_os_ssss(
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
	double V_frg[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc;
    int Lab, IJ, i;
    int ijcs, ijcs0, ijcs1;
    int kat, kao, kcs, kcs0, kcs1, kps0, kps1;
    int ijps0, nijps, nklps;
    int ics, iat, iao, jcs, jat, jao;
    double A[3], B[3], C[3], BA[3], AC[3];
    double SSSS[1];
    double veta[MAXTEM2], vdkcd[MAXTEM2];
    //
    int ixx, ncs0, dx, res, pos;
    Lab   = La * (La+1)/2 + Lb;
    ijcs0 = leading_cs_pair_frg[Lab];
    ijcs1 = leading_cs_pair_frg[Lab+1];
    if ( nworkers < 0 ) {
	pos  = workerid;
	ncs0 = leading_cs_mon[Lc+1] - leading_cs_mon[Lc];
	dx   = (ncs0>>4);	// >>3 means /8
	res  = (ncs0 & 0x000f );	// residues in division by 8
	kcs0 = leading_cs_mon[Lc] + ( pos<res ? pos*(dx+1) : pos*dx+res );
	kcs1 = ( pos<res ? kcs0+(dx+1) : kcs0 + dx );
	ixx  = 1;
    } else {
	kcs0  = leading_cs_mon[Lc] + workerid;
	kcs1  = leading_cs_mon[Lc+1];
	ixx   = nworkers;
    }

    for ( kcs=kcs0; kcs<kcs1; kcs+=ixx ) {
	kps0 = shel_add_mon[kcs];
	kps1 = kps0 + shel_tem_mon[kcs];
	kat  = shel_atm_mon[kcs];
	kao  = shel_ini_mon[kcs];
	C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	nklps = ifc3c_calc_cd_pair_params(
		prim_exp_mon, prim_coe_mon, kps0, kps1, veta, vdkcd );

	for ( ijcs=ijcs0; ijcs<ijcs1; ijcs++ ) {
	    //val_ab = csp_schwarz_frg[ijcs];
	    ics    = csp_ics_frg[ijcs];
	    jcs    = csp_jcs_frg[ijcs];
	    ijps0  = csp_leading_ps_pair_frg[ijcs];
	    nijps  = csp_leading_ps_pair_frg[ijcs+1]-ijps0;
	    iat    = shel_atm_frg[ics];
	    jat    = shel_atm_frg[jcs];
	    iao    = shel_ini_frg[ics];
	    jao    = shel_ini_frg[jcs];
	    IJ     = ((iao*iao+iao)>>1) + jao;
	    A[0]=atom_x_frg[iat]; A[1]=atom_y_frg[iat]; A[2]=atom_z_frg[iat];
	    B[0]=atom_x_frg[jat]; B[1]=atom_y_frg[jat]; B[2]=atom_z_frg[jat];
	    for ( i=0; i<3; i++ ) {
		BA[i] = B[i] - A[i];
		AC[i] = A[i] - C[i];
	    }
	    //ofmo_twoint_core_os_ssss( &La, &Lb, &Lc, &Lc,
	    twoint_core_ssss__(
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, veta, vdkcd, vxizc, DC, AC,     SSSS );
	    if ( fabs(SSSS[0]) > EPS_ERI ) {
		V_frg[IJ] += ao_pop_mon[kao] * SSSS[0];
	    }
	}	// ijcs
    }	// kcs
    return 0;
}

/** (ss,pp)タイプの３中心クーロン積分をまとめて行う関数
 * @ingroup integ-ifc3c
 * */
int ofmo_ifc3c_os_sspp(
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
	double V_frg[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc;
    int Lab, IJ, i, k, kk;
    int ijcs, ijcs0, ijcs1;
    int kat, kao, kao0, kcs, kcs0, kcs1, kps0, kps1;
    int ijps0, nijps, nklps;
    int ics, iat, iao, jcs, jat, jao;
    double A[3], B[3], C[3], BA[3], CA[3];
    double SSPP[3*3];	// 本当なら対角成分のみ
    double veta[MAXTEM2], vdkcd[MAXTEM2];
    //
    int ixx, ncs0, dx, res, pos;
    Lab   = La * (La+1)/2 + Lb;
    ijcs0 = leading_cs_pair_frg[Lab];
    ijcs1 = leading_cs_pair_frg[Lab+1];
    if ( nworkers < 0 ) {
	pos  = workerid;
	ncs0 = leading_cs_mon[Lc+1] - leading_cs_mon[Lc];
	dx   = (ncs0>>4);	// >>3 means /8
	res  = (ncs0 & 0x000f );	// residues in division by 8
	kcs0 = leading_cs_mon[Lc] + ( pos<res ? pos*(dx+1) : pos*dx+res );
	kcs1 = ( pos<res ? kcs0+(dx+1) : kcs0 + dx );
	ixx  = 1;
    } else {
	kcs0  = leading_cs_mon[Lc] + workerid;
	kcs1  = leading_cs_mon[Lc+1];
	ixx   = nworkers;
    }

    for ( kcs=kcs0; kcs<kcs1; kcs+=ixx ) {
	kps0 = shel_add_mon[kcs];
	kps1 = kps0 + shel_tem_mon[kcs];
	kat  = shel_atm_mon[kcs];
	kao0 = shel_ini_mon[kcs];
	C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	nklps = ifc3c_calc_cd_pair_params(
		prim_exp_mon, prim_coe_mon, kps0, kps1, veta, vdkcd );
	for ( ijcs=ijcs0; ijcs<ijcs1; ijcs++ ) {
	    //val_ab = csp_schwarz_frg[ijcs];
	    ics    = csp_ics_frg[ijcs];
	    jcs    = csp_jcs_frg[ijcs];
	    ijps0  = csp_leading_ps_pair_frg[ijcs];
	    nijps  = csp_leading_ps_pair_frg[ijcs+1]-ijps0;
	    iat    = shel_atm_frg[ics];
	    jat    = shel_atm_frg[jcs];
	    iao    = shel_ini_frg[ics];
	    jao    = shel_ini_frg[jcs];
	    IJ     = ((iao*iao+iao)>>1) + jao;
	    A[0]=atom_x_frg[iat]; A[1]=atom_y_frg[iat]; A[2]=atom_z_frg[iat];
	    B[0]=atom_x_frg[jat]; B[1]=atom_y_frg[jat]; B[2]=atom_z_frg[jat];
	    for ( i=0; i<3; i++ ) {
		BA[i] = B[i] - A[i];
		CA[i] = C[i] - A[i];
	    }
	    // abとcdをひっくり返す
	    ofmo_twoint_core_os_ppss( &La, &Lb, &Lc, &Lc,
		    &nklps, veta, vdkcd, vxizc, DC,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA, CA,    SSPP );

	    for ( k=0, kao=kao0; k<3; k++, kao++ ) {
		kk = k*3+k;
		if ( fabs( SSPP[kk] ) > EPS_ERI ) {
		    V_frg[IJ] += ao_pop_mon[kao] * SSPP[kk];
		}
	    }
	}	// ijcs
    }	// kcs
    return 0;
}

/** (ss,dd)タイプの３中心クーロン積分をまとめて行う関数
 * @ingroup integ-ifc3c
 * */
int ofmo_ifc3c_os_ssdd(
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
	double V_frg[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc;
    int Lab, IJ, i, k, kk;
    int ijcs, ijcs0, ijcs1;
    int kat, kao, kao0, kcs, kcs0, kcs1, kps0, kps1;
    int ijps0, nijps, nklps;
    int ics, iat, iao, jcs, jat, jao;
    double A[3], B[3], C[3], BA[3], CA[3];
    double SSDD[6*6];	// 本当なら対角成分のみ
    double veta[MAXTEM2], vdkcd[MAXTEM2];
    //
    int ixx, ncs0, dx, res, pos;
    Lab   = La * (La+1)/2 + Lb;
    ijcs0 = leading_cs_pair_frg[Lab];
    ijcs1 = leading_cs_pair_frg[Lab+1];
    if ( nworkers < 0 ) {
	pos  = workerid;
	ncs0 = leading_cs_mon[Lc+1] - leading_cs_mon[Lc];
	dx   = (ncs0>>4);	// >>3 means /8
	res  = (ncs0 & 0x000f );	// residues in division by 8
	kcs0 = leading_cs_mon[Lc] + ( pos<res ? pos*(dx+1) : pos*dx+res );
	kcs1 = ( pos<res ? kcs0+(dx+1) : kcs0 + dx );
	ixx  = 1;
    } else {
	kcs0  = leading_cs_mon[Lc] + workerid;
	kcs1  = leading_cs_mon[Lc+1];
	ixx   = nworkers;
    }

    for ( kcs=kcs0; kcs<kcs1; kcs+=ixx ) {
	kps0 = shel_add_mon[kcs];
	kps1 = kps0 + shel_tem_mon[kcs];
	kat  = shel_atm_mon[kcs];
	kao0 = shel_ini_mon[kcs];
	C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	nklps = ifc3c_calc_cd_pair_params(
		prim_exp_mon, prim_coe_mon, kps0, kps1, veta, vdkcd );
	for ( ijcs=ijcs0; ijcs<ijcs1; ijcs++ ) {
	    //val_ab = csp_schwarz_frg[ijcs];
	    ics    = csp_ics_frg[ijcs];
	    jcs    = csp_jcs_frg[ijcs];
	    ijps0  = csp_leading_ps_pair_frg[ijcs];
	    nijps  = csp_leading_ps_pair_frg[ijcs+1]-ijps0;
	    iat    = shel_atm_frg[ics];
	    jat    = shel_atm_frg[jcs];
	    iao    = shel_ini_frg[ics];
	    jao    = shel_ini_frg[jcs];
	    IJ     = ((iao*iao+iao)>>1) + jao;
	    A[0]=atom_x_frg[iat]; A[1]=atom_y_frg[iat]; A[2]=atom_z_frg[iat];
	    B[0]=atom_x_frg[jat]; B[1]=atom_y_frg[jat]; B[2]=atom_z_frg[jat];
	    for ( i=0; i<3; i++ ) {
		BA[i] = B[i] - A[i];
		CA[i] = C[i] - A[i];
	    }
	    // abとcdをひっくり返す
	    ofmo_twoint_core_os_ddss( &La, &Lb, &Lc, &Lc,
		    &nklps, veta, vdkcd, vxizc, DC,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA, CA,    SSDD );

	    for ( k=0, kao=kao0; k<6; k++, kao++ ) {
		kk = k*6+k;
		if ( fabs( SSDD[kk] ) > EPS_ERI ) {
		    V_frg[IJ] += ao_pop_mon[kao] * SSDD[kk];
		}
	    }
	}	// ijcs
    }	// kcs
    return 0;
}

/** (ps,ss)タイプの３中心クーロン積分をまとめて行う関数
 * @ingroup integ-ifc3c
 * */
int ofmo_ifc3c_os_psss(
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
	double V_frg[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc;
    int Lab, IJ, i;
    int ijcs, ijcs0, ijcs1;
    int kat, kao, kcs, kcs0, kcs1, kps0, kps1;
    int ijps0, nijps, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    double A[3], B[3], C[3], BA[3], AC[3];
    double PSSS[3];
    double veta[MAXTEM2], vdkcd[MAXTEM2];
    //
    int ixx, ncs0, dx, res, pos;
    Lab   = La * (La+1)/2 + Lb;
    ijcs0 = leading_cs_pair_frg[Lab];
    ijcs1 = leading_cs_pair_frg[Lab+1];
    if ( nworkers < 0 ) {
	pos  = workerid;
	ncs0 = leading_cs_mon[Lc+1] - leading_cs_mon[Lc];
	dx   = (ncs0>>4);	// >>3 means /8
	res  = (ncs0 & 0x000f );	// residues in division by 8
	kcs0 = leading_cs_mon[Lc] + ( pos<res ? pos*(dx+1) : pos*dx+res );
	kcs1 = ( pos<res ? kcs0+(dx+1) : kcs0 + dx );
	ixx  = 1;
    } else {
	kcs0  = leading_cs_mon[Lc] + workerid;
	kcs1  = leading_cs_mon[Lc+1];
	ixx   = nworkers;
    }

    for ( kcs=kcs0; kcs<kcs1; kcs+=ixx ) {
	kps0 = shel_add_mon[kcs];
	kps1 = kps0 + shel_tem_mon[kcs];
	kat  = shel_atm_mon[kcs];
	kao  = shel_ini_mon[kcs];
	C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	nklps = ifc3c_calc_cd_pair_params(
		prim_exp_mon, prim_coe_mon, kps0, kps1, veta, vdkcd );
	for ( ijcs=ijcs0; ijcs<ijcs1; ijcs++ ) {
	    //val_ab = csp_schwarz_frg[ijcs];
	    ics    = csp_ics_frg[ijcs];
	    jcs    = csp_jcs_frg[ijcs];
	    ijps0  = csp_leading_ps_pair_frg[ijcs];
	    nijps  = csp_leading_ps_pair_frg[ijcs+1]-ijps0;
	    iat    = shel_atm_frg[ics];
	    jat    = shel_atm_frg[jcs];
	    iao0   = shel_ini_frg[ics];
	    jao    = shel_ini_frg[jcs];
	    A[0]=atom_x_frg[iat]; A[1]=atom_y_frg[iat]; A[2]=atom_z_frg[iat];
	    B[0]=atom_x_frg[jat]; B[1]=atom_y_frg[jat]; B[2]=atom_z_frg[jat];
	    for ( i=0; i<3; i++ ) {
		BA[i] = B[i] - A[i];
		AC[i] = A[i] - C[i];
	    }
	    ofmo_twoint_core_os_psss( &La, &Lb, &Lc, &Lc,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, veta, vdkcd, vxizc, DC, AC,     PSSS );
	    for ( i=0, iao=iao0; i<3; i++, iao++ ) {
		IJ = ( (iao*iao+iao)>>1 ) + jao;
		if ( fabs(PSSS[i]) > EPS_ERI ) {
		    V_frg[IJ] += ao_pop_mon[kao] * PSSS[i];
		}
	    }
	}	// ijcs
    }	// kcs
    return 0;
}

/** (ps,pp)タイプの３中心クーロン積分をまとめて行う関数
 * @ingroup integ-ifc3c
 * */
int ofmo_ifc3c_os_pspp(
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
	double V_frg[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc;
    int Lab, IJ, i, k, kk, ix;
    int ijcs, ijcs0, ijcs1;
    int kat, kao, kao0, kcs, kcs0, kcs1, kps0, kps1;
    int ijps0, nijps, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    double A[3], B[3], C[3], BA[3], CA[3];
    double PSPP[3*3*3];	// 本来は対角成分だけ
    double veta[MAXTEM2], vdkcd[MAXTEM2];
    //
    int ixx, ncs0, dx, res, pos;
    Lab   = La * (La+1)/2 + Lb;
    ijcs0 = leading_cs_pair_frg[Lab];
    ijcs1 = leading_cs_pair_frg[Lab+1];
    if ( nworkers < 0 ) {
	pos  = workerid;
	ncs0 = leading_cs_mon[Lc+1] - leading_cs_mon[Lc];
	dx   = (ncs0>>4);	// >>3 means /8
	res  = (ncs0 & 0x000f );	// residues in division by 8
	kcs0 = leading_cs_mon[Lc] + ( pos<res ? pos*(dx+1) : pos*dx+res );
	kcs1 = ( pos<res ? kcs0+(dx+1) : kcs0 + dx );
	ixx  = 1;
    } else {
	kcs0  = leading_cs_mon[Lc] + workerid;
	kcs1  = leading_cs_mon[Lc+1];
	ixx   = nworkers;
    }

    for ( kcs=kcs0; kcs<kcs1; kcs+=ixx ) {
	kps0 = shel_add_mon[kcs];
	kps1 = kps0 + shel_tem_mon[kcs];
	kat  = shel_atm_mon[kcs];
	kao0 = shel_ini_mon[kcs];
	C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	nklps = ifc3c_calc_cd_pair_params(
		prim_exp_mon, prim_coe_mon, kps0, kps1, veta, vdkcd );
	for ( ijcs=ijcs0; ijcs<ijcs1; ijcs++ ) {
	    //val_ab = csp_schwarz_frg[ijcs];
	    ics    = csp_ics_frg[ijcs];
	    jcs    = csp_jcs_frg[ijcs];
	    ijps0  = csp_leading_ps_pair_frg[ijcs];
	    nijps  = csp_leading_ps_pair_frg[ijcs+1]-ijps0;
	    iat    = shel_atm_frg[ics];
	    jat    = shel_atm_frg[jcs];
	    iao0   = shel_ini_frg[ics];
	    jao    = shel_ini_frg[jcs];
	    A[0]=atom_x_frg[iat]; A[1]=atom_y_frg[iat]; A[2]=atom_z_frg[iat];
	    B[0]=atom_x_frg[jat]; B[1]=atom_y_frg[jat]; B[2]=atom_z_frg[jat];
	    for ( i=0; i<3; i++ ) {
		BA[i] = B[i] - A[i];
		CA[i] = C[i] - A[i];
	    }
	    // abとcdをひっくり返す
	    ofmo_twoint_core_os_ppps( &La, &Lb, &Lc, &Lc,
		    &nklps, veta, vdkcd, vxizc, DC,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA, CA,    PSPP );
	    for ( k=0, kao=kao0; k<3; k++, kao++ ) {
		kk = k*3+k;
		for ( i=0, iao=iao0, ix=kk*3; i<3; i++, iao++, ix++ ) {
		    IJ = ( (iao*iao+iao)>>1 ) + jao;
		    if ( fabs(PSPP[ix]) > EPS_ERI ) {
			V_frg[IJ] += ao_pop_mon[kao] * PSPP[ix];
		    }
		}
	    }
	}	// ijcs
    }	// kcs
    return 0;
}

/** (ps,dd)タイプの３中心クーロン積分をまとめて行う関数
 * @ingroup integ-ifc3c
 * */
int ofmo_ifc3c_os_psdd(
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
	double V_frg[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc;
    int Lab, IJ, i, k, kk, ix;
    int ijcs, ijcs0, ijcs1;
    int kat, kao, kao0, kcs, kcs0, kcs1, kps0, kps1;
    int ijps0, nijps, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    double A[3], B[3], C[3], BA[3], CA[3];
    double PSDD[3*6*6];	// 本当なら対角成分のみ
    double veta[MAXTEM2], vdkcd[MAXTEM2];
    //
    int ixx, ncs0, dx, res, pos;
    Lab   = La * (La+1)/2 + Lb;
    ijcs0 = leading_cs_pair_frg[Lab];
    ijcs1 = leading_cs_pair_frg[Lab+1];
    if ( nworkers < 0 ) {
	pos  = workerid;
	ncs0 = leading_cs_mon[Lc+1] - leading_cs_mon[Lc];
	dx   = (ncs0>>4);	// >>3 means /8
	res  = (ncs0 & 0x000f );	// residues in division by 8
	kcs0 = leading_cs_mon[Lc] + ( pos<res ? pos*(dx+1) : pos*dx+res );
	kcs1 = ( pos<res ? kcs0+(dx+1) : kcs0 + dx );
	ixx  = 1;
    } else {
	kcs0  = leading_cs_mon[Lc] + workerid;
	kcs1  = leading_cs_mon[Lc+1];
	ixx   = nworkers;
    }

    for ( kcs=kcs0; kcs<kcs1; kcs+=ixx ) {
	kps0 = shel_add_mon[kcs];
	kps1 = kps0 + shel_tem_mon[kcs];
	kat  = shel_atm_mon[kcs];
	kao0 = shel_ini_mon[kcs];
	C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	nklps = ifc3c_calc_cd_pair_params(
		prim_exp_mon, prim_coe_mon, kps0, kps1, veta, vdkcd );
	for ( ijcs=ijcs0; ijcs<ijcs1; ijcs++ ) {
	    //val_ab = csp_schwarz_frg[ijcs];
	    ics    = csp_ics_frg[ijcs];
	    jcs    = csp_jcs_frg[ijcs];
	    ijps0  = csp_leading_ps_pair_frg[ijcs];
	    nijps  = csp_leading_ps_pair_frg[ijcs+1]-ijps0;
	    iat    = shel_atm_frg[ics];
	    jat    = shel_atm_frg[jcs];
	    iao0   = shel_ini_frg[ics];
	    jao    = shel_ini_frg[jcs];
	    A[0]=atom_x_frg[iat]; A[1]=atom_y_frg[iat]; A[2]=atom_z_frg[iat];
	    B[0]=atom_x_frg[jat]; B[1]=atom_y_frg[jat]; B[2]=atom_z_frg[jat];
	    for ( i=0; i<3; i++ ) {
		BA[i] = B[i] - A[i];
		CA[i] = C[i] - A[i];
	    }
	    // abとcdをひっくり返す
	    ofmo_twoint_core_os_ddps( &La, &Lb, &Lc, &Lc,
		    &nklps, veta, vdkcd, vxizc, DC,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA, CA,    PSDD );

	    for ( k=0, kao=kao0; k<6; k++, kao++ ) {
		kk = k*6+k;
		for ( i=0, iao=iao0, ix=kk*3; i<3; i++, iao++, ix++ ) {
		    if ( fabs( PSDD[ix] ) > EPS_ERI ) {
			IJ = ((iao*iao+iao)>>1) + jao;
			V_frg[IJ] += ao_pop_mon[kao] * PSDD[ix];
		    }
		}
	    }
	}	// ijcs
    }	// kcs
    return 0;
}

/** (pp,ss)タイプの３中心クーロン積分をまとめて行う関数
 * @ingroup integ-ifc3c
 * */
int ofmo_ifc3c_os_ppss(
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
	double V_frg[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc;
    int Lab, IJ, I2, i, j, ix;
    int ijcs, ijcs0, ijcs1;
    int kat, kao, kcs, kcs0, kcs1, kps0, kps1;
    int ijps0, nijps, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    double A[3], B[3], C[3], BA[3], AC[3];
    double PPSS[3*3];
    double veta[MAXTEM2], vdkcd[MAXTEM2];
    //
    int ixx, ncs0, dx, res, pos;
    Lab   = La * (La+1)/2 + Lb;
    ijcs0 = leading_cs_pair_frg[Lab];
    ijcs1 = leading_cs_pair_frg[Lab+1];
    if ( nworkers < 0 ) {
	pos  = workerid;
	ncs0 = leading_cs_mon[Lc+1] - leading_cs_mon[Lc];
	dx   = (ncs0>>4);	// >>3 means /8
	res  = (ncs0 & 0x000f );	// residues in division by 8
	kcs0 = leading_cs_mon[Lc] + ( pos<res ? pos*(dx+1) : pos*dx+res );
	kcs1 = ( pos<res ? kcs0+(dx+1) : kcs0 + dx );
	ixx  = 1;
    } else {
	kcs0  = leading_cs_mon[Lc] + workerid;
	kcs1  = leading_cs_mon[Lc+1];
	ixx   = nworkers;
    }

    for ( kcs=kcs0; kcs<kcs1; kcs+=ixx ) {
	kps0 = shel_add_mon[kcs];
	kps1 = kps0 + shel_tem_mon[kcs];
	kat  = shel_atm_mon[kcs];
	kao  = shel_ini_mon[kcs];
	C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	nklps = ifc3c_calc_cd_pair_params(
		prim_exp_mon, prim_coe_mon, kps0, kps1, veta, vdkcd );
	for ( ijcs=ijcs0; ijcs<ijcs1; ijcs++ ) {
	    //val_ab = csp_schwarz_frg[ijcs];
	    ics    = csp_ics_frg[ijcs];
	    jcs    = csp_jcs_frg[ijcs];
	    ijps0  = csp_leading_ps_pair_frg[ijcs];
	    nijps  = csp_leading_ps_pair_frg[ijcs+1]-ijps0;
	    iat    = shel_atm_frg[ics];
	    jat    = shel_atm_frg[jcs];
	    iao0   = shel_ini_frg[ics];
	    jao0   = shel_ini_frg[jcs];
	    A[0]=atom_x_frg[iat]; A[1]=atom_y_frg[iat]; A[2]=atom_z_frg[iat];
	    B[0]=atom_x_frg[jat]; B[1]=atom_y_frg[jat]; B[2]=atom_z_frg[jat];
	    for ( i=0; i<3; i++ ) {
		BA[i] = B[i] - A[i];
		AC[i] = A[i] - C[i];
	    }
	    ofmo_twoint_core_os_ppss( &La, &Lb, &Lc, &Lc,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, veta, vdkcd, vxizc, DC, AC,     PPSS );
	    for ( i=0, iao=iao0, ix=0; i<3; i++, iao++ ) {
		I2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<3; j++, jao++, ix++ ) {
		    if ( jao > iao ) continue;
		    IJ = I2 + jao;
		    if ( fabs(PPSS[ix]) > EPS_ERI ) {
			V_frg[IJ] += ao_pop_mon[kao] * PPSS[ix];
		    }
		}
	    }
	}	// ijcs
    }	// kcs
    return 0;
}

/** (pp,pp)タイプの３中心クーロン積分をまとめて行う関数
 * @ingroup integ-ifc3c
 * */
int ofmo_ifc3c_os_pppp(
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
	double V_frg[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc;
    int Lab, IJ, I2, i, j, k, kk, ix;
    int ijcs, ijcs0, ijcs1;
    int kat, kao, kao0, kcs, kcs0, kcs1, kps0, kps1;
    int ijps0, nijps, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    double A[3], B[3], C[3], BA[3], AC[3];
    double PPPP[3*3*3*3];	// 本来は対角成分だけ
    double veta[MAXTEM2], vdkcd[MAXTEM2];
    //
    int ixx, ncs0, dx, res, pos;
    Lab   = La * (La+1)/2 + Lb;
    ijcs0 = leading_cs_pair_frg[Lab];
    ijcs1 = leading_cs_pair_frg[Lab+1];
    if ( nworkers < 0 ) {
	pos  = workerid;
	ncs0 = leading_cs_mon[Lc+1] - leading_cs_mon[Lc];
	dx   = (ncs0>>4);	// >>3 means /8
	res  = (ncs0 & 0x000f );	// residues in division by 8
	kcs0 = leading_cs_mon[Lc] + ( pos<res ? pos*(dx+1) : pos*dx+res );
	kcs1 = ( pos<res ? kcs0+(dx+1) : kcs0 + dx );
	ixx  = 1;
    } else {
	kcs0  = leading_cs_mon[Lc] + workerid;
	kcs1  = leading_cs_mon[Lc+1];
	ixx   = nworkers;
    }

    for ( kcs=kcs0; kcs<kcs1; kcs+=ixx ) {
	kps0 = shel_add_mon[kcs];
	kps1 = kps0 + shel_tem_mon[kcs];
	kat  = shel_atm_mon[kcs];
	kao0 = shel_ini_mon[kcs];
	C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	nklps = ifc3c_calc_cd_pair_params(
		prim_exp_mon, prim_coe_mon, kps0, kps1, veta, vdkcd );
	for ( ijcs=ijcs0; ijcs<ijcs1; ijcs++ ) {
	    //val_ab = csp_schwarz_frg[ijcs];
	    ics    = csp_ics_frg[ijcs];
	    jcs    = csp_jcs_frg[ijcs];
	    ijps0  = csp_leading_ps_pair_frg[ijcs];
	    nijps  = csp_leading_ps_pair_frg[ijcs+1]-ijps0;
	    iat    = shel_atm_frg[ics];
	    jat    = shel_atm_frg[jcs];
	    iao0   = shel_ini_frg[ics];
	    jao0   = shel_ini_frg[jcs];
	    A[0]=atom_x_frg[iat]; A[1]=atom_y_frg[iat]; A[2]=atom_z_frg[iat];
	    B[0]=atom_x_frg[jat]; B[1]=atom_y_frg[jat]; B[2]=atom_z_frg[jat];
	    for ( i=0; i<3; i++ ) {
		BA[i] = B[i] - A[i];
		AC[i] = A[i] - C[i];
	    }
	    ofmo_twoint_core_os_pppp( &La, &Lb, &Lc, &Lc,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, veta, vdkcd, vxizc, DC, AC,     PPPP );
	    for ( i=0, iao=iao0, ix=0; i<3; i++, iao++ ) {
		I2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<3; j++, jao++, ix+=9 ) {
		    if ( jao > iao ) continue;
		    IJ = I2 + jao;
		    for ( k=0, kao=kao0; k<3; k++, kao++ ) {
			kk = k*3+k;
			if ( fabs(PPPP[ix+kk]) > EPS_ERI ) {
			    V_frg[IJ] += ao_pop_mon[kao] * PPPP[ix+kk];
			}
		    }
		}
	    }
	}	// ijcs
    }	// kcs
    return 0;
}

/** (pp,dd)タイプの３中心クーロン積分をまとめて行う関数
 * @ingroup integ-ifc3c
 * */
int ofmo_ifc3c_os_ppdd(
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
	double V_frg[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc;
    int Lab, IJ, I2, i, j, k, kk, ix;
    int ijcs, ijcs0, ijcs1;
    int kat, kao, kao0, kcs, kcs0, kcs1, kps0, kps1;
    int ijps0, nijps, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    double A[3], B[3], C[3], BA[3], CA[3];
    double PPDD[3*3*6*6];	// 本当なら対角成分のみ
    double veta[MAXTEM2], vdkcd[MAXTEM2];
    //
    int ixx, ncs0, dx, res, pos;
    Lab   = La * (La+1)/2 + Lb;
    ijcs0 = leading_cs_pair_frg[Lab];
    ijcs1 = leading_cs_pair_frg[Lab+1];
    if ( nworkers < 0 ) {
	pos  = workerid;
	ncs0 = leading_cs_mon[Lc+1] - leading_cs_mon[Lc];
	dx   = (ncs0>>4);	// >>3 means /8
	res  = (ncs0 & 0x000f );	// residues in division by 8
	kcs0 = leading_cs_mon[Lc] + ( pos<res ? pos*(dx+1) : pos*dx+res );
	kcs1 = ( pos<res ? kcs0+(dx+1) : kcs0 + dx );
	ixx  = 1;
    } else {
	kcs0  = leading_cs_mon[Lc] + workerid;
	kcs1  = leading_cs_mon[Lc+1];
	ixx   = nworkers;
    }

    for ( kcs=kcs0; kcs<kcs1; kcs+=ixx ) {
	kps0 = shel_add_mon[kcs];
	kps1 = kps0 + shel_tem_mon[kcs];
	kat  = shel_atm_mon[kcs];
	kao0 = shel_ini_mon[kcs];
	C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	nklps = ifc3c_calc_cd_pair_params(
		prim_exp_mon, prim_coe_mon, kps0, kps1, veta, vdkcd );
	for ( ijcs=ijcs0; ijcs<ijcs1; ijcs++ ) {
	    //val_ab = csp_schwarz_frg[ijcs];
	    ics    = csp_ics_frg[ijcs];
	    jcs    = csp_jcs_frg[ijcs];
	    ijps0  = csp_leading_ps_pair_frg[ijcs];
	    nijps  = csp_leading_ps_pair_frg[ijcs+1]-ijps0;
	    iat    = shel_atm_frg[ics];
	    jat    = shel_atm_frg[jcs];
	    iao0   = shel_ini_frg[ics];
	    jao0   = shel_ini_frg[jcs];
	    //IJ     = ((iao*iao+iao)>>1) + jao;
	    A[0]=atom_x_frg[iat]; A[1]=atom_y_frg[iat]; A[2]=atom_z_frg[iat];
	    B[0]=atom_x_frg[jat]; B[1]=atom_y_frg[jat]; B[2]=atom_z_frg[jat];
	    for ( i=0; i<3; i++ ) {
		BA[i] = B[i] - A[i];
		CA[i] = C[i] - A[i];
	    }
	    // abとcdをひっくり返す
	    ofmo_twoint_core_os_ddpp( &La, &Lb, &Lc, &Lc,
		    &nklps, veta, vdkcd, vxizc, DC,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA, CA,    PPDD );
	    for ( k=0, kao=kao0; k<6; k++, kao++ ) {
		kk = k*6+k;
		for ( i=0, iao=iao0, ix=kk*9; i<3; i++, iao++ ) {
		    I2 = (iao*iao+iao)>>1;
		    for ( j=0, jao=jao0; j<3; j++, jao++, ix++ ) {
			if ( jao > iao ) continue;
			IJ = I2 + jao;
			if ( fabs( PPDD[ix] ) > EPS_ERI ) {
			    V_frg[IJ] += ao_pop_mon[kao] * PPDD[ix];
			}
		    }
		}
	    }
	}	// ijcs
    }	// kcs
    return 0;
}

/** (ds,ss)タイプの３中心クーロン積分をまとめて行う関数
 * @ingroup integ-ifc3c
 * */
int ofmo_ifc3c_os_dsss(
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
	double V_frg[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc;
    int Lab, IJ, i;
    int ijcs, ijcs0, ijcs1;
    int kat, kao, kcs, kcs0, kcs1, kps0, kps1;
    int ijps0, nijps, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    double A[3], B[3], C[3], BA[3], AC[3];
    double DSSS[6];
    double veta[MAXTEM2], vdkcd[MAXTEM2];
    //
    int ixx, ncs0, dx, res, pos;
    Lab   = La * (La+1)/2 + Lb;
    ijcs0 = leading_cs_pair_frg[Lab];
    ijcs1 = leading_cs_pair_frg[Lab+1];
    if ( nworkers < 0 ) {
	pos  = workerid;
	ncs0 = leading_cs_mon[Lc+1] - leading_cs_mon[Lc];
	dx   = (ncs0>>4);	// >>3 means /8
	res  = (ncs0 & 0x000f );	// residues in division by 8
	kcs0 = leading_cs_mon[Lc] + ( pos<res ? pos*(dx+1) : pos*dx+res );
	kcs1 = ( pos<res ? kcs0+(dx+1) : kcs0 + dx );
	ixx  = 1;
    } else {
	kcs0  = leading_cs_mon[Lc] + workerid;
	kcs1  = leading_cs_mon[Lc+1];
	ixx   = nworkers;
    }

    for ( kcs=kcs0; kcs<kcs1; kcs+=ixx ) {
	kps0 = shel_add_mon[kcs];
	kps1 = kps0 + shel_tem_mon[kcs];
	kat  = shel_atm_mon[kcs];
	kao  = shel_ini_mon[kcs];
	C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	nklps = ifc3c_calc_cd_pair_params(
		prim_exp_mon, prim_coe_mon, kps0, kps1, veta, vdkcd );
	for ( ijcs=ijcs0; ijcs<ijcs1; ijcs++ ) {
	    //val_ab = csp_schwarz_frg[ijcs];
	    ics    = csp_ics_frg[ijcs];
	    jcs    = csp_jcs_frg[ijcs];
	    ijps0  = csp_leading_ps_pair_frg[ijcs];
	    nijps  = csp_leading_ps_pair_frg[ijcs+1]-ijps0;
	    iat    = shel_atm_frg[ics];
	    jat    = shel_atm_frg[jcs];
	    iao0   = shel_ini_frg[ics];
	    jao    = shel_ini_frg[jcs];
	    A[0]=atom_x_frg[iat]; A[1]=atom_y_frg[iat]; A[2]=atom_z_frg[iat];
	    B[0]=atom_x_frg[jat]; B[1]=atom_y_frg[jat]; B[2]=atom_z_frg[jat];
	    for ( i=0; i<3; i++ ) {
		BA[i] = B[i] - A[i];
		AC[i] = A[i] - C[i];
	    }
	    ofmo_twoint_core_os_dsss( &La, &Lb, &Lc, &Lc,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, veta, vdkcd, vxizc, DC, AC,     DSSS );
	    for ( i=0, iao=iao0; i<6; i++, iao++ ) {
		IJ = ( (iao*iao+iao)>>1 ) + jao;
		if ( fabs(DSSS[i]) > EPS_ERI ) {
		    V_frg[IJ] += ao_pop_mon[kao] * DSSS[i];
		}
	    }
	}	// ijcs
    }	// kcs
    return 0;
}

/** (ds,pp)タイプの３中心クーロン積分をまとめて行う関数
 * @ingroup integ-ifc3c
 * */
int ofmo_ifc3c_os_dspp(
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
	double V_frg[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc;
    int Lab, IJ, i, k, kk;
    int ijcs, ijcs0, ijcs1;
    int kat, kao, kao0, kcs, kcs0, kcs1, kps0, kps1;
    int ijps0, nijps, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    double A[3], B[3], C[3], BA[3], AC[3];
    double DSPP[6*3*3];	// 本来は対角成分だけ
    double veta[MAXTEM2], vdkcd[MAXTEM2];
    //
    int ixx, ncs0, dx, res, pos;
    Lab   = La * (La+1)/2 + Lb;
    ijcs0 = leading_cs_pair_frg[Lab];
    ijcs1 = leading_cs_pair_frg[Lab+1];
    if ( nworkers < 0 ) {
	pos  = workerid;
	ncs0 = leading_cs_mon[Lc+1] - leading_cs_mon[Lc];
	dx   = (ncs0>>4);	// >>3 means /8
	res  = (ncs0 & 0x000f );	// residues in division by 8
	kcs0 = leading_cs_mon[Lc] + ( pos<res ? pos*(dx+1) : pos*dx+res );
	kcs1 = ( pos<res ? kcs0+(dx+1) : kcs0 + dx );
	ixx  = 1;
    } else {
	kcs0  = leading_cs_mon[Lc] + workerid;
	kcs1  = leading_cs_mon[Lc+1];
	ixx   = nworkers;
    }

    for ( kcs=kcs0; kcs<kcs1; kcs+=ixx ) {
	kps0 = shel_add_mon[kcs];
	kps1 = kps0 + shel_tem_mon[kcs];
	kat  = shel_atm_mon[kcs];
	kao0 = shel_ini_mon[kcs];
	C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	nklps = ifc3c_calc_cd_pair_params(
		prim_exp_mon, prim_coe_mon, kps0, kps1, veta, vdkcd );
	for ( ijcs=ijcs0; ijcs<ijcs1; ijcs++ ) {
	    //val_ab = csp_schwarz_frg[ijcs];
	    ics    = csp_ics_frg[ijcs];
	    jcs    = csp_jcs_frg[ijcs];
	    ijps0  = csp_leading_ps_pair_frg[ijcs];
	    nijps  = csp_leading_ps_pair_frg[ijcs+1]-ijps0;
	    iat    = shel_atm_frg[ics];
	    jat    = shel_atm_frg[jcs];
	    iao0   = shel_ini_frg[ics];
	    jao    = shel_ini_frg[jcs];
	    A[0]=atom_x_frg[iat]; A[1]=atom_y_frg[iat]; A[2]=atom_z_frg[iat];
	    B[0]=atom_x_frg[jat]; B[1]=atom_y_frg[jat]; B[2]=atom_z_frg[jat];
	    for ( i=0; i<3; i++ ) {
		BA[i] = B[i] - A[i];
		AC[i] = A[i] - C[i];
	    }
	    ofmo_twoint_core_os_dspp( &La, &Lb, &Lc, &Lc,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, veta, vdkcd, vxizc, DC, AC,     DSPP );
	    for ( i=0, iao=iao0; i<6; i++, iao++ ) {
		IJ = ( (iao*iao+iao)>>1 ) + jao;
		for ( k=0, kao=kao0; k<3; k++, kao++ ) {
		    kk = k*3+k;
		    if ( fabs(DSPP[i*9+kk]) > EPS_ERI ) {
			V_frg[IJ] += ao_pop_mon[kao] * DSPP[i*9+kk];
		    }
		}
	    }
	}	// ijcs
    }	// kcs
    return 0;
}

/** (ds,dd)タイプの３中心クーロン積分をまとめて行う関数
 * @ingroup integ-ifc3c
 * */
int ofmo_ifc3c_os_dsdd(
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
	double V_frg[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc;
    int Lab, IJ, i, k, kk, ix;
    int ijcs, ijcs0, ijcs1;
    int kat, kao, kao0, kcs, kcs0, kcs1, kps0, kps1;
    int ijps0, nijps, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    double A[3], B[3], C[3], BA[3], CA[3];
    double DSDD[6*6*6];	// 本当なら対角成分のみ
    double veta[MAXTEM2], vdkcd[MAXTEM2];
    //
    int ixx, ncs0, dx, res, pos;
    Lab   = La * (La+1)/2 + Lb;
    ijcs0 = leading_cs_pair_frg[Lab];
    ijcs1 = leading_cs_pair_frg[Lab+1];
    if ( nworkers < 0 ) {
	pos  = workerid;
	ncs0 = leading_cs_mon[Lc+1] - leading_cs_mon[Lc];
	dx   = (ncs0>>4);	// >>3 means /8
	res  = (ncs0 & 0x000f );	// residues in division by 8
	kcs0 = leading_cs_mon[Lc] + ( pos<res ? pos*(dx+1) : pos*dx+res );
	kcs1 = ( pos<res ? kcs0+(dx+1) : kcs0 + dx );
	ixx  = 1;
    } else {
	kcs0  = leading_cs_mon[Lc] + workerid;
	kcs1  = leading_cs_mon[Lc+1];
	ixx   = nworkers;
    }

    for ( kcs=kcs0; kcs<kcs1; kcs+=ixx ) {
	kps0 = shel_add_mon[kcs];
	kps1 = kps0 + shel_tem_mon[kcs];
	kat  = shel_atm_mon[kcs];
	kao0 = shel_ini_mon[kcs];
	C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	nklps = ifc3c_calc_cd_pair_params(
		prim_exp_mon, prim_coe_mon, kps0, kps1, veta, vdkcd );
	for ( ijcs=ijcs0; ijcs<ijcs1; ijcs++ ) {
	    //val_ab = csp_schwarz_frg[ijcs];
	    ics    = csp_ics_frg[ijcs];
	    jcs    = csp_jcs_frg[ijcs];
	    ijps0  = csp_leading_ps_pair_frg[ijcs];
	    nijps  = csp_leading_ps_pair_frg[ijcs+1]-ijps0;
	    iat    = shel_atm_frg[ics];
	    jat    = shel_atm_frg[jcs];
	    iao0   = shel_ini_frg[ics];
	    jao    = shel_ini_frg[jcs];
	    IJ     = ((iao*iao+iao)>>1) + jao;
	    A[0]=atom_x_frg[iat]; A[1]=atom_y_frg[iat]; A[2]=atom_z_frg[iat];
	    B[0]=atom_x_frg[jat]; B[1]=atom_y_frg[jat]; B[2]=atom_z_frg[jat];
	    for ( i=0; i<3; i++ ) {
		BA[i] = B[i] - A[i];
		CA[i] = C[i] - A[i];
	    }
	    // abとcdをひっくり返す
	    ofmo_twoint_core_os_ddds( &La, &Lb, &Lc, &Lc,
		    &nklps, veta, vdkcd, vxizc, DC,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA, CA,    DSDD );

	    for ( k=0, kao=kao0; k<6; k++, kao++ ) {
		kk = k*6+k;
		for ( i=0, iao=iao0, ix=kk*6; i<6; i++, iao++, ix++ ) {
		    IJ = ( (iao*iao+iao)>>1 ) + jao;
		    if ( fabs( DSDD[ix] ) > EPS_ERI ) {
			V_frg[IJ] += ao_pop_mon[kao] * DSDD[ix];
		    }
		}
	    }
	}	// ijcs
    }	// kcs
    return 0;
}

/** (dp,ss)タイプの３中心クーロン積分をまとめて行う関数
 * @ingroup integ-ifc3c
 * */
int ofmo_ifc3c_os_dpss(
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
	double V_frg[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc;
    int Lab, IJ, I2, i, j, ix;
    int ijcs, ijcs0, ijcs1;
    int kat, kao, kcs, kcs0, kcs1, kps0, kps1;
    int ijps0, nijps, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    double A[3], B[3], C[3], BA[3], AC[3];
    double DPSS[6*3];
    double veta[MAXTEM2], vdkcd[MAXTEM2];
    //
    int ixx, ncs0, dx, res, pos;
    Lab   = La * (La+1)/2 + Lb;
    ijcs0 = leading_cs_pair_frg[Lab];
    ijcs1 = leading_cs_pair_frg[Lab+1];
    if ( nworkers < 0 ) {
	pos  = workerid;
	ncs0 = leading_cs_mon[Lc+1] - leading_cs_mon[Lc];
	dx   = (ncs0>>4);	// >>3 means /8
	res  = (ncs0 & 0x000f );	// residues in division by 8
	kcs0 = leading_cs_mon[Lc] + ( pos<res ? pos*(dx+1) : pos*dx+res );
	kcs1 = ( pos<res ? kcs0+(dx+1) : kcs0 + dx );
	ixx  = 1;
    } else {
	kcs0  = leading_cs_mon[Lc] + workerid;
	kcs1  = leading_cs_mon[Lc+1];
	ixx   = nworkers;
    }

    for ( kcs=kcs0; kcs<kcs1; kcs+=ixx ) {
	kps0 = shel_add_mon[kcs];
	kps1 = kps0 + shel_tem_mon[kcs];
	kat  = shel_atm_mon[kcs];
	kao  = shel_ini_mon[kcs];
	C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	nklps = ifc3c_calc_cd_pair_params(
		prim_exp_mon, prim_coe_mon, kps0, kps1, veta, vdkcd );
	for ( ijcs=ijcs0; ijcs<ijcs1; ijcs++ ) {
	    //val_ab = csp_schwarz_frg[ijcs];
	    ics    = csp_ics_frg[ijcs];
	    jcs    = csp_jcs_frg[ijcs];
	    ijps0  = csp_leading_ps_pair_frg[ijcs];
	    nijps  = csp_leading_ps_pair_frg[ijcs+1]-ijps0;
	    iat    = shel_atm_frg[ics];
	    jat    = shel_atm_frg[jcs];
	    iao0   = shel_ini_frg[ics];
	    jao0   = shel_ini_frg[jcs];
	    A[0]=atom_x_frg[iat]; A[1]=atom_y_frg[iat]; A[2]=atom_z_frg[iat];
	    B[0]=atom_x_frg[jat]; B[1]=atom_y_frg[jat]; B[2]=atom_z_frg[jat];
	    for ( i=0; i<3; i++ ) {
		BA[i] = B[i] - A[i];
		AC[i] = A[i] - C[i];
	    }
	    ofmo_twoint_core_os_dpss( &La, &Lb, &Lc, &Lc,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, veta, vdkcd, vxizc, DC, AC,     DPSS );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		I2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<3; j++, jao++, ix++ ) {
		    IJ = I2 + jao;
		    if ( fabs(DPSS[ix]) > EPS_ERI ) {
			V_frg[IJ] += ao_pop_mon[kao] * DPSS[ix];
		    }
		}
	    }
	}	// ijcs
    }	// kcs
    return 0;
}

/** (dp,pp)タイプの３中心クーロン積分をまとめて行う関数
 * @ingroup integ-ifc3c
 * */
int ofmo_ifc3c_os_dppp(
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
	double V_frg[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc;
    int Lab, IJ, I2, i, j, k, kk, ix;
    int ijcs, ijcs0, ijcs1;
    int kat, kao, kao0, kcs, kcs0, kcs1, kps0, kps1;
    int ijps0, nijps, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    double A[3], B[3], C[3], BA[3], AC[3];
    double DPPP[6*3*3*3];	// 本来は対角成分だけ
    double veta[MAXTEM2], vdkcd[MAXTEM2];
    //
    int ixx, ncs0, dx, res, pos;
    Lab   = La * (La+1)/2 + Lb;
    ijcs0 = leading_cs_pair_frg[Lab];
    ijcs1 = leading_cs_pair_frg[Lab+1];
    if ( nworkers < 0 ) {
	pos  = workerid;
	ncs0 = leading_cs_mon[Lc+1] - leading_cs_mon[Lc];
	dx   = (ncs0>>4);	// >>3 means /8
	res  = (ncs0 & 0x000f );	// residues in division by 8
	kcs0 = leading_cs_mon[Lc] + ( pos<res ? pos*(dx+1) : pos*dx+res );
	kcs1 = ( pos<res ? kcs0+(dx+1) : kcs0 + dx );
	ixx  = 1;
    } else {
	kcs0  = leading_cs_mon[Lc] + workerid;
	kcs1  = leading_cs_mon[Lc+1];
	ixx   = nworkers;
    }

    for ( kcs=kcs0; kcs<kcs1; kcs+=ixx ) {
	kps0 = shel_add_mon[kcs];
	kps1 = kps0 + shel_tem_mon[kcs];
	kat  = shel_atm_mon[kcs];
	kao0 = shel_ini_mon[kcs];
	C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	nklps = ifc3c_calc_cd_pair_params(
		prim_exp_mon, prim_coe_mon, kps0, kps1, veta, vdkcd );
	for ( ijcs=ijcs0; ijcs<ijcs1; ijcs++ ) {
	    //val_ab = csp_schwarz_frg[ijcs];
	    ics    = csp_ics_frg[ijcs];
	    jcs    = csp_jcs_frg[ijcs];
	    ijps0  = csp_leading_ps_pair_frg[ijcs];
	    nijps  = csp_leading_ps_pair_frg[ijcs+1]-ijps0;
	    iat    = shel_atm_frg[ics];
	    jat    = shel_atm_frg[jcs];
	    iao0   = shel_ini_frg[ics];
	    jao0   = shel_ini_frg[jcs];
	    A[0]=atom_x_frg[iat]; A[1]=atom_y_frg[iat]; A[2]=atom_z_frg[iat];
	    B[0]=atom_x_frg[jat]; B[1]=atom_y_frg[jat]; B[2]=atom_z_frg[jat];
	    for ( i=0; i<3; i++ ) {
		BA[i] = B[i] - A[i];
		AC[i] = A[i] - C[i];
	    }
	    ofmo_twoint_core_os_dppp( &La, &Lb, &Lc, &Lc,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, veta, vdkcd, vxizc, DC, AC,     DPPP );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		I2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<3; j++, jao++, ix+=9 ) {
		    IJ = I2 + jao;
		    for ( k=0, kao=kao0; k<3; k++, kao++ ) {
			kk = k*3+k;
			if ( fabs(DPPP[ix+kk]) > EPS_ERI ) {
			    V_frg[IJ] += ao_pop_mon[kao] * DPPP[ix+kk];
			}
		    }
		}
	    }
	}	// ijcs
    }	// kcs
    return 0;
}

/** (dp,dd)タイプの３中心クーロン積分をまとめて行う関数
 * @ingroup integ-ifc3c
 * */
int ofmo_ifc3c_os_dpdd(
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
	double V_frg[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc;
    int Lab, IJ, I2, i, j, k, kk, ix;
    int ijcs, ijcs0, ijcs1;
    int kat, kao, kao0, kcs, kcs0, kcs1, kps0, kps1;
    int ijps0, nijps, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    double A[3], B[3], C[3], BA[3], CA[3];
    double DPDD[6*3*6*6];	// 本当なら対角成分のみ
    double veta[MAXTEM2], vdkcd[MAXTEM2];
    //
    int ixx, ncs0, dx, res, pos;
    Lab   = La * (La+1)/2 + Lb;
    ijcs0 = leading_cs_pair_frg[Lab];
    ijcs1 = leading_cs_pair_frg[Lab+1];
    if ( nworkers < 0 ) {
	pos  = workerid;
	ncs0 = leading_cs_mon[Lc+1] - leading_cs_mon[Lc];
	dx   = (ncs0>>4);	// >>3 means /8
	res  = (ncs0 & 0x000f );	// residues in division by 8
	kcs0 = leading_cs_mon[Lc] + ( pos<res ? pos*(dx+1) : pos*dx+res );
	kcs1 = ( pos<res ? kcs0+(dx+1) : kcs0 + dx );
	ixx  = 1;
    } else {
	kcs0  = leading_cs_mon[Lc] + workerid;
	kcs1  = leading_cs_mon[Lc+1];
	ixx   = nworkers;
    }

    for ( kcs=kcs0; kcs<kcs1; kcs+=ixx ) {
	kps0 = shel_add_mon[kcs];
	kps1 = kps0 + shel_tem_mon[kcs];
	kat  = shel_atm_mon[kcs];
	kao0 = shel_ini_mon[kcs];
	C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	nklps = ifc3c_calc_cd_pair_params(
		prim_exp_mon, prim_coe_mon, kps0, kps1, veta, vdkcd );
	for ( ijcs=ijcs0; ijcs<ijcs1; ijcs++ ) {
	    //val_ab = csp_schwarz_frg[ijcs];
	    ics    = csp_ics_frg[ijcs];
	    jcs    = csp_jcs_frg[ijcs];
	    ijps0  = csp_leading_ps_pair_frg[ijcs];
	    nijps  = csp_leading_ps_pair_frg[ijcs+1]-ijps0;
	    iat    = shel_atm_frg[ics];
	    jat    = shel_atm_frg[jcs];
	    iao0   = shel_ini_frg[ics];
	    jao0   = shel_ini_frg[jcs];
	    IJ     = ((iao*iao+iao)>>1) + jao;
	    A[0]=atom_x_frg[iat]; A[1]=atom_y_frg[iat]; A[2]=atom_z_frg[iat];
	    B[0]=atom_x_frg[jat]; B[1]=atom_y_frg[jat]; B[2]=atom_z_frg[jat];
	    for ( i=0; i<3; i++ ) {
		BA[i] = B[i] - A[i];
		CA[i] = C[i] - A[i];
	    }
	    // abとcdをひっくり返す
	    ofmo_twoint_core_os_dddp( &La, &Lb, &Lc, &Lc,
		    &nklps, veta, vdkcd, vxizc, DC,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA, CA,    DPDD );
	    for ( k=0, kao=kao0; k<6; k++, kao++ ) {
		kk = k*6+k;
		for ( i=0, iao=iao0, ix=kk*6*3; i<6; i++, iao++ ) {
		    I2 = (iao*iao+iao)>>1;
		    for ( j=0, jao=jao0; j<3; j++, jao++, ix++ ) {
			IJ = I2 + jao;
			if ( fabs( DPDD[ix] ) > EPS_ERI ) {
			    V_frg[IJ] += ao_pop_mon[kao] * DPDD[ix];
			}
		    }
		}
	    }
	}	// ijcs
    }	// kcs
    return 0;
}

/** (dd,ss)タイプの３中心クーロン積分をまとめて行う関数
 * @ingroup integ-ifc3c
 * */
int ofmo_ifc3c_os_ddss(
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
	double V_frg[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc;
    int Lab, IJ, I2, i, j, ix;
    int ijcs, ijcs0, ijcs1;
    int kat, kao, kcs, kcs0, kcs1, kps0, kps1;
    int ijps0, nijps, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    double A[3], B[3], C[3], BA[3], AC[3];
    double DDSS[6*6];
    double veta[MAXTEM2], vdkcd[MAXTEM2];
    //
    int ixx, ncs0, dx, res, pos;
    Lab   = La * (La+1)/2 + Lb;
    ijcs0 = leading_cs_pair_frg[Lab];
    ijcs1 = leading_cs_pair_frg[Lab+1];
    if ( nworkers < 0 ) {
	pos  = workerid;
	ncs0 = leading_cs_mon[Lc+1] - leading_cs_mon[Lc];
	dx   = (ncs0>>4);	// >>3 means /8
	res  = (ncs0 & 0x000f );	// residues in division by 8
	kcs0 = leading_cs_mon[Lc] + ( pos<res ? pos*(dx+1) : pos*dx+res );
	kcs1 = ( pos<res ? kcs0+(dx+1) : kcs0 + dx );
	ixx  = 1;
    } else {
	kcs0  = leading_cs_mon[Lc] + workerid;
	kcs1  = leading_cs_mon[Lc+1];
	ixx   = nworkers;
    }

    for ( kcs=kcs0; kcs<kcs1; kcs+=ixx ) {
	kps0 = shel_add_mon[kcs];
	kps1 = kps0 + shel_tem_mon[kcs];
	kat  = shel_atm_mon[kcs];
	kao  = shel_ini_mon[kcs];
	C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	nklps = ifc3c_calc_cd_pair_params(
		prim_exp_mon, prim_coe_mon, kps0, kps1, veta, vdkcd );
	for ( ijcs=ijcs0; ijcs<ijcs1; ijcs++ ) {
	    //val_ab = csp_schwarz_frg[ijcs];
	    ics    = csp_ics_frg[ijcs];
	    jcs    = csp_jcs_frg[ijcs];
	    ijps0  = csp_leading_ps_pair_frg[ijcs];
	    nijps  = csp_leading_ps_pair_frg[ijcs+1]-ijps0;
	    iat    = shel_atm_frg[ics];
	    jat    = shel_atm_frg[jcs];
	    iao0   = shel_ini_frg[ics];
	    jao0   = shel_ini_frg[jcs];
	    A[0]=atom_x_frg[iat]; A[1]=atom_y_frg[iat]; A[2]=atom_z_frg[iat];
	    B[0]=atom_x_frg[jat]; B[1]=atom_y_frg[jat]; B[2]=atom_z_frg[jat];
	    for ( i=0; i<3; i++ ) {
		BA[i] = B[i] - A[i];
		AC[i] = A[i] - C[i];
	    }
	    ofmo_twoint_core_os_ddss( &La, &Lb, &Lc, &Lc,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, veta, vdkcd, vxizc, DC, AC,     DDSS );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		I2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<6; j++, jao++, ix++ ) {
		    if ( jao > iao ) continue;
		    IJ = I2 + jao;
		    if ( fabs(DDSS[ix]) > EPS_ERI ) {
			V_frg[IJ] += ao_pop_mon[kao] * DDSS[ix];
		    }
		}
	    }
	}	// ijcs
    }	// kcs
    return 0;
}

/** (dd,pp)タイプの３中心クーロン積分をまとめて行う関数
 * @ingroup integ-ifc3c
 * */
int ofmo_ifc3c_os_ddpp(
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
	double V_frg[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc;
    int Lab, IJ, I2, i, j, k, kk, ix;
    int ijcs, ijcs0, ijcs1;
    int kat, kao, kao0, kcs, kcs0, kcs1, kps0, kps1;
    int ijps0, nijps, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    double A[3], B[3], C[3], BA[3], AC[3];
    double DDPP[6*6*3*3];	// 本来は対角成分だけ
    double veta[MAXTEM2], vdkcd[MAXTEM2];
    //
    int ixx, ncs0, dx, res, pos;
    Lab   = La * (La+1)/2 + Lb;
    ijcs0 = leading_cs_pair_frg[Lab];
    ijcs1 = leading_cs_pair_frg[Lab+1];
    if ( nworkers < 0 ) {
	pos  = workerid;
	ncs0 = leading_cs_mon[Lc+1] - leading_cs_mon[Lc];
	dx   = (ncs0>>4);	// >>3 means /8
	res  = (ncs0 & 0x000f );	// residues in division by 8
	kcs0 = leading_cs_mon[Lc] + ( pos<res ? pos*(dx+1) : pos*dx+res );
	kcs1 = ( pos<res ? kcs0+(dx+1) : kcs0 + dx );
	ixx  = 1;
    } else {
	kcs0  = leading_cs_mon[Lc] + workerid;
	kcs1  = leading_cs_mon[Lc+1];
	ixx   = nworkers;
    }

    for ( kcs=kcs0; kcs<kcs1; kcs+=ixx ) {
	kps0 = shel_add_mon[kcs];
	kps1 = kps0 + shel_tem_mon[kcs];
	kat  = shel_atm_mon[kcs];
	kao0 = shel_ini_mon[kcs];
	C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	nklps = ifc3c_calc_cd_pair_params(
		prim_exp_mon, prim_coe_mon, kps0, kps1, veta, vdkcd );
	for ( ijcs=ijcs0; ijcs<ijcs1; ijcs++ ) {
	    //val_ab = csp_schwarz_frg[ijcs];
	    ics    = csp_ics_frg[ijcs];
	    jcs    = csp_jcs_frg[ijcs];
	    ijps0  = csp_leading_ps_pair_frg[ijcs];
	    nijps  = csp_leading_ps_pair_frg[ijcs+1]-ijps0;
	    iat    = shel_atm_frg[ics];
	    jat    = shel_atm_frg[jcs];
	    iao0   = shel_ini_frg[ics];
	    jao0   = shel_ini_frg[jcs];
	    A[0]=atom_x_frg[iat]; A[1]=atom_y_frg[iat]; A[2]=atom_z_frg[iat];
	    B[0]=atom_x_frg[jat]; B[1]=atom_y_frg[jat]; B[2]=atom_z_frg[jat];
	    for ( i=0; i<3; i++ ) {
		BA[i] = B[i] - A[i];
		AC[i] = A[i] - C[i];
	    }
	    ofmo_twoint_core_os_ddpp( &La, &Lb, &Lc, &Lc,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, veta, vdkcd, vxizc, DC, AC,     DDPP );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		I2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<6; j++, jao++, ix+=9 ) {
		    if ( jao > iao ) continue;
		    IJ = I2 + jao;
		    for ( k=0, kao=kao0; k<3; k++, kao++ ) {
			kk = k*3+k;
			if ( fabs(DDPP[ix+kk]) > EPS_ERI ) {
			    V_frg[IJ] += ao_pop_mon[kao] * DDPP[ix+kk];
			}
		    }
		}
	    }
	}	// ijcs
    }	// kcs
    return 0;
}

/** (dd,dd)タイプの３中心クーロン積分をまとめて行う関数
 * @ingroup integ-ifc3c
 * */
int ofmo_ifc3c_os_dddd(
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
	double V_frg[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc;
    int Lab, IJ, I2, i, j, k, kk, ix;
    int ijcs, ijcs0, ijcs1;
    int kat, kao, kao0, kcs, kcs0, kcs1, kps0, kps1;
    int ijps0, nijps, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    double A[3], B[3], C[3], BA[3], AC[3];
    double DDDD[6*6*6*6];	// 本当なら対角成分のみ
    double veta[MAXTEM2], vdkcd[MAXTEM2];
    //
    int ixx, ncs0, dx, res, pos;
    Lab   = La * (La+1)/2 + Lb;
    ijcs0 = leading_cs_pair_frg[Lab];
    ijcs1 = leading_cs_pair_frg[Lab+1];
    if ( nworkers < 0 ) {
	pos  = workerid;
	ncs0 = leading_cs_mon[Lc+1] - leading_cs_mon[Lc];
	dx   = (ncs0>>4);	// >>3 means /8
	res  = (ncs0 & 0x000f );	// residues in division by 8
	kcs0 = leading_cs_mon[Lc] + ( pos<res ? pos*(dx+1) : pos*dx+res );
	kcs1 = ( pos<res ? kcs0+(dx+1) : kcs0 + dx );
	ixx  = 1;
    } else {
	kcs0  = leading_cs_mon[Lc] + workerid;
	kcs1  = leading_cs_mon[Lc+1];
	ixx   = nworkers;
    }

    for ( kcs=kcs0; kcs<kcs1; kcs+=ixx ) {
	kps0 = shel_add_mon[kcs];
	kps1 = kps0 + shel_tem_mon[kcs];
	kat  = shel_atm_mon[kcs];
	kao0 = shel_ini_mon[kcs];
	C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	nklps = ifc3c_calc_cd_pair_params(
		prim_exp_mon, prim_coe_mon, kps0, kps1, veta, vdkcd );
	for ( ijcs=ijcs0; ijcs<ijcs1; ijcs++ ) {
	    //val_ab = csp_schwarz_frg[ijcs];
	    ics    = csp_ics_frg[ijcs];
	    jcs    = csp_jcs_frg[ijcs];
	    ijps0  = csp_leading_ps_pair_frg[ijcs];
	    nijps  = csp_leading_ps_pair_frg[ijcs+1]-ijps0;
	    iat    = shel_atm_frg[ics];
	    jat    = shel_atm_frg[jcs];
	    iao0   = shel_ini_frg[ics];
	    jao0   = shel_ini_frg[jcs];
	    A[0]=atom_x_frg[iat]; A[1]=atom_y_frg[iat]; A[2]=atom_z_frg[iat];
	    B[0]=atom_x_frg[jat]; B[1]=atom_y_frg[jat]; B[2]=atom_z_frg[jat];
	    for ( i=0; i<3; i++ ) {
		BA[i] = B[i] - A[i];
		AC[i] = A[i] - C[i];
	    }
	    ofmo_twoint_core_os_dddd( &La, &Lb, &Lc, &Lc,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, veta, vdkcd, vxizc, DC, AC,     DDDD );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		I2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<6; j++, jao++, ix+=36 ) {
		    if ( jao > iao ) continue;
		    IJ = I2 + jao;
		    for ( k=0, kao=kao0; k<6; k++, kao++ ) {
			kk = k*6+k;
			if ( fabs(DDDD[ix+kk]) > EPS_ERI ) {
			    V_frg[IJ] += ao_pop_mon[kao] * DDDD[ix+kk];
			}
		    }
		}
	    }
	}	// ijcs
    }	// kcs
    return 0;
}

