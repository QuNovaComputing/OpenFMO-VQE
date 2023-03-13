/**
 * @file ofmo-ifc4c.c 同じタイプの３中心クーロン積分を行う関数群
 * */
/**
 * @defgroup integ-ifc4c ４中心クーロン積分を行う関数群
 *
 * 同じタイプの4中心クーロン積分を行い、環境ポテンシャル項に加算する
 * 関数群。
 *
 * すべての関数が同じ引数をもっているので、以下にその内容を示す。
 *
 * @param[in] nworkers 計算に用いるワーカプロセス（スレッド）数
 * @param[in] workerid 各ワーカプロセス（スレッド）のID
 *     （\f$ 0\le\tt{workerid}<\tt{nworkers} \f$）
 * @param[in] La １つ目の軌道量子数
 * @param[in] Lb ２つ目の軌道量子数（\f$ \tt{La} \ge \tt{Lb} \f$）
 * @param[in] Lc ３つ目の軌道量子数
 * @param[in] Ld ４つ目の軌道量子数（\f$ \tt{Lc} \ge \tt{Ld} \f$）
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
 * @param[int] csp_schwarz_frg[icsp] フラグメントのCSペア番号 \c icsp の
 * Schwarz積分の値
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
 * @param[in] shel_atm_mon[ics] 相手モノマーのCS番号 \c ics のCSが属する
 * 原子の番号
 * @param[in] shel_ini_mon[ics] 相手モノマーのCS番号 \c ics のCSに
 * 含まれるAOの先頭AO番号
 * @param[in] atom_x_mon[iat] 相手モノマーの原子の番号 \c iat のx座標
 * （au単位）
 * @param[in] atom_y_mon[iat] 相手モノマーの原子の番号 \c iat のy座標
 * （au単位）
 * @param[in] atom_z_mon[iat] 相手モノマーの原子の番号 \c iat のz座標
 * （au単位）
 * @param[in] leading_cs_pair_mon[itype] 相手モノマーのCSペアタイプ番号
 * \c itype の先頭CSペア番号
 * @param[int] csp_schwarz_mon[icsp] 相手モノマーのCSペア番号 \c icsp の
 * Schwarz積分の値
 * @param[in] csp_ics_mon[icsp] 相手モノマーのCSペア番号 \c icsp の
 * 1つ目のCS番号
 * @param[in] csp_jcs_mon[icsp] 相手モノマーのCSペア番号 \c icsp の
 * 2つめのCS番号。ただし、
 * \f$ \tt{csp\_ics[icsp]} \ge \tt{csp\_jcs[icsp]} \f$ である。
 * @param[in] csp_leading_ps_pair_mon[icsp] 相手モノマーのCSペア番号
 * \c icsp に含まれるPSペアの先頭PSペア番号
 * @param[in] psp_zeta_mon[ipsp] 相手モノマーのPSペア番号 \c ipsp の
 * 軌道指数和\f$ \zeta = \zeta_a + \zeta_b \f$
 * @param[in] psp_dkps_mon[ipsp] 相手モノマーのPSペア番号 \c ipsp
 * の線型結合定数
 *     \f[ K_{ab} = \sqrt2 \pi^{5/4} \frac1{\zeta_a+\zeta_b}
 *     \exp\left[ -\frac{\zeta_a \zeta_b}{\zeta_a + \zeta_b}
 *     ( \boldmath A \unboldmath - \boldmath B \unboldmath )^2
 *     \right]\f]
 * @param[in] psp_xiza_mon[ipsp] 相手モノマーのPSペア番号 \c ipsp の
 *     \f$ \frac{\xi}{\zeta_a} = \frac{\zeta_b}{\zeta_a+\zeta_b} \f$
 *
 * @param[in] D_mon[] 相手モノマーの密度行列（正方行列形式）
 *
 * @param[out] V_frg[] この関数で計算した４中心クーロンポテンシャルが
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

#define MAXPSPAIR       100

#define ZERO    0.e0
#define ONE     1.e0
#define TWO     2.e0
#define FOUR    4.e0
#define HALF    0.5e0

#define EPS_ERI 1.e-15
#define EPS_PS4 1.e-30

#include "ofmo-twoint-core.h"

/** (ss,ss)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_ssss__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, jcs, jat, jao;
    int kcs, kat, kao, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    double SSSS[1];
    //
    int ixx, ncsp, dx, res, pos;

    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao    = shel_ini_mon[kcs];
	    lao    = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    KL    = ((kao*kao+kao)>>1) + lao;
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_ssss__(
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,   AC,      SSSS );
	    coe = (kao==lao? ONE : TWO );
	    if ( fabs(SSSS[0]) > EPS_ERI ) {
		V_frg[IJ] += coe * D_mon[KL] * SSSS[0];
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_ssss_


/** (ss,ps)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_ssps__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, k;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, jcs, jat, jao;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], CA[3];
    double val_ab, val_cd;
    double SSPS[3];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao    = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		CA[i] = C[i] - A[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_psss__(
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,   CA,      SSPS );
	    for ( k=0, kao=kao0; k<3; k++, kao++ ) {
		KL = ((kao*kao+kao)>>1) + lao;
		if ( fabs(SSPS[k]) > EPS_ERI ) {
		    V_frg[IJ] += TWO * D_mon[KL] * SSPS[k];
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_ssps_


/** (ss,pp)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_sspp__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, k, l, K2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, jcs, jat, jao;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], CA[3];
    double val_ab, val_cd, coe;
    double SSPP[3*3];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao0   = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		CA[i] = C[i] - A[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_ppss__(
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,   CA,      SSPP );
	    for ( k=0, kao=kao0, ix=0; k<3; k++, kao++ ) {
		K2 = (kao*kao+kao)>>1;
		for ( l=0, lao=lao0; l<3; l++, lao++, ix++ ) {
		    if ( lao>kao ) continue;
		    KL = K2 + lao;
		    coe = (kao==lao? ONE : TWO );
		    if ( fabs(SSPP[ix]) > EPS_ERI ) {
			V_frg[IJ] += coe * D_mon[KL] * SSPP[ix];
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_sspp_


/** (ss,ds)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_ssds__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, k;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, jcs, jat, jao;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], CA[3];
    double val_ab, val_cd;
    double SSDS[6];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao    = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		CA[i] = C[i] - A[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dsss__(
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,   CA,      SSDS );
	    for ( k=0, kao=kao0; k<6; k++, kao++ ) {
		KL = ((kao*kao+kao)>>1) + lao;
		if ( fabs(SSDS[k]) > EPS_ERI ) {
		    V_frg[IJ] += TWO * D_mon[KL] * SSDS[k];
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_ssds_


/** (ss,dp)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_ssdp__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, k, l, K2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, jcs, jat, jao;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], CA[3];
    double val_ab, val_cd;
    double SSDP[6*3];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao0   = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		CA[i] = C[i] - A[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dpss__(
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,   CA,      SSDP );
	    for ( k=0, kao=kao0, ix=0; k<6; k++, kao++ ) {
		K2 = (kao*kao+kao)>>1;
		for ( l=0, lao=lao0; l<3; l++, lao++, ix++ ) {
		    KL = K2 + lao;
		    if ( fabs(SSDP[ix]) > EPS_ERI ) {
			V_frg[IJ] += TWO * D_mon[KL] * SSDP[ix];
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_ssdp_


/** (ss,dd)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_ssdd__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, k, l, K2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, jcs, jat, jao;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], CA[3];
    double val_ab, val_cd, coe;
    double SSDD[6*6];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao0   = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		CA[i] = C[i] - A[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_ddss__(
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,   CA,      SSDD );
	    for ( k=0, kao=kao0, ix=0; k<6; k++, kao++ ) {
		K2 = (kao*kao+kao)>>1;
		for ( l=0, lao=lao0; l<6; l++, lao++, ix++ ) {
		    if ( lao>kao ) continue;
		    KL = K2 + lao;
		    coe = (kao==lao? ONE : TWO );
		    if ( fabs(SSDD[ix]) > EPS_ERI ) {
			V_frg[IJ] += coe * D_mon[KL] * SSDD[ix];
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_ssdd_


/** (ps,ss)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_psss__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    int kcs, kat, kao, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    double PSSS[3];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao    = shel_ini_mon[kcs];
	    lao    = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    KL    = ((kao*kao+kao)>>1) + lao;
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_psss__(
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,   AC,      PSSS );
	    coe = (kao==lao? ONE : TWO );
	    for ( i=0, iao=iao0; i<3; i++, iao++ ) {
		IJ = ((iao*iao+iao)>>1) + jao;
		if ( fabs(PSSS[i]) > EPS_ERI ) {
		    V_frg[IJ] += coe * D_mon[KL] * PSSS[i];
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_psss_


/** (ps,ps)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_psps__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, k, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd;
    double PSPS[3*3];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao    = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_psps__(
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,   AC,      PSPS );
	    for ( i=0, iao=iao0, ix=0; i<3; i++, iao++ ) {
		IJ = ((iao*iao+iao)>>1) + jao;
		for ( k=0, kao=kao0; k<3; k++, kao++, ix++ ) {
		    KL = ((kao*kao+kao)>>1) + lao;
		    if ( fabs(PSPS[ix]) > EPS_ERI ) {
			V_frg[IJ] += TWO * D_mon[KL] * PSPS[ix];
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_psps_


/** (ps,pp)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_pspp__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, k, l, K2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], CA[3];
    double val_ab, val_cd, coe;
    double PSPP[3*3*3];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao0   = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		CA[i] = C[i] - A[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_ppps__(
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,   CA,      PSPP );
	    for ( k=0, kao=kao0, ix=0; k<3; k++, kao++ ) {
		K2 = (kao*kao+kao)>>1;
		for ( l=0, lao=lao0; l<3; l++, lao++ ) {
		    if ( lao>kao ) { ix+=3; continue; }
		    KL = K2 + lao;
		    coe = (kao==lao? ONE : TWO );
		    for ( i=0, iao=iao0; i<3; i++, iao++, ix++ ) {
			IJ = ((iao*iao+iao)>>1) + jao;
			if ( fabs(PSPP[ix]) > EPS_ERI ) {
			    V_frg[IJ] += coe * D_mon[KL] * PSPP[ix];
			}
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_pspp_


/** (ps,ds)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_psds__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, k, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], CA[3];
    double val_ab, val_cd;
    double PSDS[3*6];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao    = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		CA[i] = C[i] - A[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dsps__(
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,   CA,      PSDS );
	    for ( k=0, kao=kao0, ix=0; k<6; k++, kao++ ) {
		KL = ((kao*kao+kao)>>1) + lao;
		for ( i=0, iao=iao0; i<3; i++, iao++, ix++ ) {
		    IJ = ((iao*iao+iao)>>1) + jao;
		    if ( fabs(PSDS[ix]) > EPS_ERI ) {
			V_frg[IJ] += TWO * D_mon[KL] * PSDS[ix];
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_psds_


/** (ps,dp)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_psdp__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, k, l, K2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], CA[3];
    double val_ab, val_cd;
    double PSDP[3*6*3];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao0   = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		CA[i] = C[i] - A[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dpps__(
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,   CA,      PSDP );
	    for ( k=0, kao=kao0, ix=0; k<6; k++, kao++ ) {
		K2 = (kao*kao+kao)>>1;
		for ( l=0, lao=lao0; l<3; l++, lao++ ) {
		    KL = K2 + lao;
		    for ( i=0, iao=iao0; i<3; i++, iao++, ix++ ) {
			IJ = ((iao*iao+iao)>>1) + jao;
			if ( fabs(PSDP[ix]) > EPS_ERI ) {
			    V_frg[IJ] += TWO * D_mon[KL] * PSDP[ix];
			}
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_psdp_


/** (ps,dd)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_psdd__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, k, l, K2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], CA[3];
    double val_ab, val_cd, coe;
    double PSDD[3*6*6];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao0   = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		CA[i] = C[i] - A[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_ddps__(
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,   CA,      PSDD );
	    for ( k=0, kao=kao0, ix=0; k<6; k++, kao++ ) {
		K2 = (kao*kao+kao)>>1;
		for ( l=0, lao=lao0; l<6; l++, lao++ ) {
		    if ( lao>kao ) { ix+=3; continue; }
		    KL = K2 + lao;
		    coe = (kao==lao? ONE : TWO );
		    for ( i=0, iao=iao0; i<3; i++, iao++, ix++ ) {
			IJ = ((iao*iao+iao)>>1) + jao;
			if ( fabs(PSDD[ix]) > EPS_ERI ) {
			    V_frg[IJ] += coe * D_mon[KL] * PSDD[ix];
			}
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_psdd_


/** (pp,ss)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_ppss__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, j, I2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    double PPSS[3*3];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao    = shel_ini_mon[kcs];
	    lao    = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    KL    = ((kao*kao+kao)>>1) + lao;
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_ppss__(
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,   AC,      PPSS );
	    coe = (kao==lao? ONE : TWO );
	    for ( i=0, iao=iao0, ix=0; i<3; i++, iao++ ) {
		I2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<3; j++, jao++, ix++ ) {
		    if ( jao>iao ) continue;
		    IJ = I2 + jao;
		    if ( fabs(PPSS[ix]) > EPS_ERI ) {
			V_frg[IJ] += coe * D_mon[KL] * PPSS[ix];
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_ppss_


/** (pp,ps)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_ppps__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, j, k, I2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd;
    double PPPS[3*3*3];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao    = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_ppps__(
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,   AC,      PPPS );
	    for ( i=0, iao=iao0, ix=0; i<3; i++, iao++ ) {
		I2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<3; j++, jao++ ) {
		    if ( jao>iao ) { ix+=3; continue; }
		    IJ = I2 + jao;
		    for ( k=0, kao=kao0; k<3; k++, kao++, ix++ ) {
			KL = ((kao*kao+kao)>>1) + lao;
			if ( fabs(PPPS[ix]) > EPS_ERI ) {
			    V_frg[IJ] += TWO * D_mon[KL] * PPPS[ix];
			}
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_ppps_


/** (pp,pp)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_pppp__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, j, k, l, I2, K2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    double PPPP[3*3*3*3];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao0   = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_pppp__(
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,   AC,      PPPP );
	    for ( i=0, iao=iao0, ix=0; i<3; i++, iao++ ) {
		I2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<3; j++, jao++ ) {
		    if ( jao>iao ) { ix+=3*3; continue; }
		    IJ = I2 + jao;
		    for ( k=0, kao=kao0; k<3; k++, kao++ ) {
			K2 = (kao*kao+kao)>>1;
			for ( l=0, lao=lao0; l<3; l++, lao++, ix++ ) {
			    if ( lao>kao ) continue;
			    coe = (kao==lao? ONE : TWO );
			    KL = K2 + lao;
			    if ( fabs(PPPP[ix]) > EPS_ERI ) {
				V_frg[IJ] += coe * D_mon[KL] * PPPP[ix];
			    }
			}
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_pppp_


/** (pp,ds)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_ppds__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, j, k, I2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], CA[3];
    double val_ab, val_cd;
    double PPDS[3*3*6];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao    = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		CA[i] = C[i] - A[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dspp__(
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,   CA,      PPDS );
	    for ( k=0, kao=kao0, ix=0; k<6; k++, kao++ ) {
		KL = ((kao*kao+kao)>>1) + lao;
		for ( i=0, iao=iao0; i<3; i++, iao++ ) {
		    I2 = (iao*iao+iao)>>1;
		    for ( j=0, jao=jao0; j<3; j++, jao++, ix++ ) {
			if ( jao>iao ) continue;
			IJ = I2 + jao;
			if ( fabs(PPDS[ix]) > EPS_ERI ) {
			    V_frg[IJ] += TWO * D_mon[KL] * PPDS[ix];
			}
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_ppds_


/** (pp,dp)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_ppdp__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, j, k, l, I2, K2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], CA[3];
    double val_ab, val_cd;
    double PPDP[3*3*6*3];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao0   = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		CA[i] = C[i] - A[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dppp__(
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,   CA,      PPDP );
	    for ( k=0, kao=kao0, ix=0; k<6; k++, kao++ ) {
		K2 = (kao*kao+kao)>>1;
		for ( l=0, lao=lao0; l<3; l++, lao++ ) {
		    KL = K2 + lao;
		    for ( i=0, iao=iao0; i<3; i++, iao++ ) {
			I2 = (iao*iao+iao)>>1;
			for ( j=0, jao=jao0; j<3; j++, jao++, ix++ ) {
			    if ( jao>iao ) continue;
			    IJ = I2 + jao;
			    if ( fabs(PPDP[ix]) > EPS_ERI ) {
				V_frg[IJ] += TWO * D_mon[KL] * PPDP[ix];
			    }
			}
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_ppdp_


/** (pp,dd)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_ppdd__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, j, k, l, I2, K2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], CA[3];
    double val_ab, val_cd, coe;
    double PPDD[3*3*6*6];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao0   = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		CA[i] = C[i] - A[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_ddpp__(
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,   CA,      PPDD );
	    for ( k=0, kao=kao0, ix=0; k<6; k++, kao++ ) {
		K2 = (kao*kao+kao)>>1;
		for ( l=0, lao=lao0; l<6; l++, lao++ ) {
		    if ( lao>kao ) { ix+=3*3; continue; }
		    KL = K2 + lao;
		    coe = (kao==lao? ONE : TWO );
		    for ( i=0, iao=iao0; i<3; i++, iao++ ) {
			I2 = (iao*iao+iao)>>1;
			for ( j=0, jao=jao0; j<3; j++, jao++, ix++ ) {
			    if ( jao>iao ) continue;
			    IJ = I2 + jao;
			    if ( fabs(PPDD[ix]) > EPS_ERI ) {
				V_frg[IJ] += coe * D_mon[KL] * PPDD[ix];
			    }
			}
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_ppdd_


/** (ds,ss)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_dsss__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    int kcs, kat, kao, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    double DSSS[6];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao    = shel_ini_mon[kcs];
	    lao    = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    KL    = ((kao*kao+kao)>>1) + lao;
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dsss__(
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,   AC,      DSSS );
	    coe = (kao==lao? ONE : TWO );
	    for ( i=0, iao=iao0; i<6; i++, iao++ ) {
		IJ = ((iao*iao+iao)>>1) + jao;
		if ( fabs(DSSS[i]) > EPS_ERI ) {
		    V_frg[IJ] += coe * D_mon[KL] * DSSS[i];
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_dsss_


/** (ds,ps)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_dsps__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, k, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd;
    double DSPS[6*3];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao    = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dsps__(
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,   AC,      DSPS );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		IJ = ((iao*iao+iao)>>1) + jao;
		for ( k=0, kao=kao0; k<3; k++, kao++, ix++ ) {
		    KL = ((kao*kao+kao)>>1) + lao;
		    if ( fabs(DSPS[ix]) > EPS_ERI ) {
			V_frg[IJ] += TWO * D_mon[KL] * DSPS[ix];
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_dsps_


/** (ds,pp)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_dspp__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, k, l, K2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    double DSPP[6*3*3];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao0   = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dspp__(
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,   AC,      DSPP );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		IJ = ((iao*iao+iao)>>1) + jao;
		for ( k=0, kao=kao0; k<3; k++, kao++ ) {
		    K2 = (kao*kao+kao)>>1;
		    for ( l=0, lao=lao0; l<3; l++, lao++, ix++ ) {
			if ( lao>kao ) continue;
			coe = (kao==lao? ONE : TWO );
			KL = K2 + lao;
			if ( fabs(DSPP[ix]) > EPS_ERI ) {
			    V_frg[IJ] += coe * D_mon[KL] * DSPP[ix];
			}
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_dspp_


/** (ds,ds)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_dsds__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, k, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd;
    double DSDS[6*6];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao    = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dsds__(
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,   AC,      DSDS );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		IJ = ((iao*iao+iao)>>1) + jao;
		for ( k=0, kao=kao0; k<6; k++, kao++, ix++ ) {
		    KL = ((kao*kao+kao)>>1) + lao;
		    if ( fabs(DSDS[ix]) > EPS_ERI ) {
			V_frg[IJ] += TWO * D_mon[KL] * DSDS[ix];
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_dsds_


/** (ds,dp)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_dsdp__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, k, l, K2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], CA[3];
    double val_ab, val_cd;
    double DSDP[6*6*3];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao0   = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		CA[i] = C[i] - A[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dpds__(
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,   CA,      DSDP );
	    for ( k=0, kao=kao0, ix=0; k<6; k++, kao++ ) {
		K2 = (kao*kao+kao)>>1;
		for ( l=0, lao=lao0; l<3; l++, lao++ ) {
		    KL = K2 + lao;
		    for ( i=0, iao=iao0; i<6; i++, iao++, ix++ ) {
			IJ = ((iao*iao+iao)>>1) + jao;
			if ( fabs(DSDP[ix]) > EPS_ERI ) {
			    V_frg[IJ] += TWO * D_mon[KL] * DSDP[ix];
			}
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_dsdp_


/** (ds,dd)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_dsdd__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, k, l, K2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], CA[3];
    double val_ab, val_cd, coe;
    double DSDD[6*6*6];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao0   = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		CA[i] = C[i] - A[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_ddds__(
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,   CA,      DSDD );
	    for ( k=0, kao=kao0, ix=0; k<6; k++, kao++ ) {
		K2 = (kao*kao+kao)>>1;
		for ( l=0, lao=lao0; l<6; l++, lao++ ) {
		    if ( lao>kao ) { ix+=6; continue; }
		    KL = K2 + lao;
		    coe = (kao==lao? ONE : TWO );
		    for ( i=0, iao=iao0; i<6; i++, iao++, ix++ ) {
			IJ = ((iao*iao+iao)>>1) + jao;
			if ( fabs(DSDD[ix]) > EPS_ERI ) {
			    V_frg[IJ] += coe * D_mon[KL] * DSDD[ix];
			}
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_dsdd_


/** (dp,ss)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_dpss__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, j, I2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    double DPSS[6*3];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao    = shel_ini_mon[kcs];
	    lao    = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    KL    = ((kao*kao+kao)>>1) + lao;
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dpss__(
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,   AC,      DPSS );
	    coe = (kao==lao? ONE : TWO );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		I2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<3; j++, jao++, ix++ ) {
		    IJ = I2 + jao;
		    if ( fabs(DPSS[ix]) > EPS_ERI ) {
			V_frg[IJ] += coe * D_mon[KL] * DPSS[ix];
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_dpss_


/** (dp,ps)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_dpps__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, j, k, I2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd;
    double DPPS[6*3*3];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao    = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dpps__(
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,   AC,      DPPS );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		I2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<3; j++, jao++ ) {
		    IJ = I2 + jao;
		    for ( k=0, kao=kao0; k<3; k++, kao++, ix++ ) {
			KL = ((kao*kao+kao)>>1) + lao;
			if ( fabs(DPPS[ix]) > EPS_ERI ) {
			    V_frg[IJ] += TWO * D_mon[KL] * DPPS[ix];
			}
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_dpps_


/** (dp,pp)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_dppp__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, j, k, l, I2, K2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    double DPPP[6*3*3*3];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao0   = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dppp__(
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,   AC,      DPPP );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		I2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<3; j++, jao++ ) {
		    IJ = I2 + jao;
		    for ( k=0, kao=kao0; k<3; k++, kao++ ) {
			K2 = (kao*kao+kao)>>1;
			for ( l=0, lao=lao0; l<3; l++, lao++, ix++ ) {
			    if ( lao>kao ) continue;
			    coe = (kao==lao? ONE : TWO );
			    KL = K2 + lao;
			    if ( fabs(DPPP[ix]) > EPS_ERI ) {
				V_frg[IJ] += coe * D_mon[KL] * DPPP[ix];
			    }
			}
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_dppp_


/** (dp,ds)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_dpds__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, j, k, I2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd;
    double DPDS[6*3*6];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao    = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dpds__(
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,   AC,      DPDS );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		I2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<3; j++, jao++ ) {
		    IJ = I2 + jao;
		    for ( k=0, kao=kao0; k<6; k++, kao++, ix++ ) {
			KL = ((kao*kao+kao)>>1) + lao;
			if ( fabs(DPDS[ix]) > EPS_ERI ) {
			    V_frg[IJ] += TWO * D_mon[KL] * DPDS[ix];
			}
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_dpds_


/** (dp,dp)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_dpdp__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, j, k, l, I2, K2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd;
    double DPDP[6*3*6*3];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao0   = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dpdp__(
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,   AC,      DPDP );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		I2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<3; j++, jao++ ) {
		    IJ = I2 + jao;
		    for ( k=0, kao=kao0; k<6; k++, kao++ ) {
			K2 = (kao*kao+kao)>>1;
			for ( l=0, lao=lao0; l<3; l++, lao++, ix++ ) {
			    KL = K2 + lao;
			    if ( fabs(DPDP[ix]) > EPS_ERI ) {
				V_frg[IJ] += TWO * D_mon[KL] * DPDP[ix];
			    }
			}
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_dpdp_


/** (dp,dd)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_dpdd__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, j, k, l, I2, K2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], CA[3];
    double val_ab, val_cd, coe;
    double DPDD[6*3*6*6];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao0   = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		CA[i] = C[i] - A[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dddp__(
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,   CA,      DPDD );
	    for ( k=0, kao=kao0, ix=0; k<6; k++, kao++ ) {
		K2 = (kao*kao+kao)>>1;
		for ( l=0, lao=lao0; l<6; l++, lao++ ) {
		    if ( lao>kao ) { ix+=6*3; continue; }
		    KL = K2 + lao;
		    coe = (kao==lao? ONE : TWO );
		    for ( i=0, iao=iao0; i<6; i++, iao++ ) {
			I2 = (iao*iao+iao)>>1;
			for ( j=0, jao=jao0; j<3; j++, jao++, ix++ ) {
			    IJ = I2 + jao;
			    if ( fabs(DPDD[ix]) > EPS_ERI ) {
				V_frg[IJ] += coe * D_mon[KL] * DPDD[ix];
			    }
			}
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_dpdd_


/** (dd,ss)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_ddss__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, j, I2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    double DDSS[6*6];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao    = shel_ini_mon[kcs];
	    lao    = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    KL    = ((kao*kao+kao)>>1) + lao;
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_ddss__(
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,   AC,      DDSS );
	    coe = (kao==lao? ONE : TWO );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		I2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<6; j++, jao++, ix++ ) {
		    if ( jao>iao ) continue;
		    IJ = I2 + jao;
		    if ( fabs(DDSS[ix]) > EPS_ERI ) {
			V_frg[IJ] += coe * D_mon[KL] * DDSS[ix];
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_ddss_


/** (dd,ps)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_ddps__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, j, k, I2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd;
    double DDPS[6*6*3];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao    = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_ddps__(
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,   AC,      DDPS );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		I2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<6; j++, jao++ ) {
		    if ( jao>iao ) { ix+=3; continue; }
		    IJ = I2 + jao;
		    for ( k=0, kao=kao0; k<3; k++, kao++, ix++ ) {
			KL = ((kao*kao+kao)>>1) + lao;
			if ( fabs(DDPS[ix]) > EPS_ERI ) {
			    V_frg[IJ] += TWO * D_mon[KL] * DDPS[ix];
			}
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_ddps_


/** (dd,pp)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_ddpp__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, j, k, l, I2, K2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    double DDPP[6*6*3*3];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao0   = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_ddpp__(
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,   AC,      DDPP );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		I2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<6; j++, jao++ ) {
		    if ( jao>iao ) { ix+=3*3; continue; }
		    IJ = I2 + jao;
		    for ( k=0, kao=kao0; k<3; k++, kao++ ) {
			K2 = (kao*kao+kao)>>1;
			for ( l=0, lao=lao0; l<3; l++, lao++, ix++ ) {
			    if ( lao>kao ) continue;
			    coe = (kao==lao? ONE : TWO );
			    KL = K2 + lao;
			    if ( fabs(DDPP[ix]) > EPS_ERI ) {
				V_frg[IJ] += coe * D_mon[KL] * DDPP[ix];
			    }
			}
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_ddpp_


/** (dd,ds)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_ddds__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, j, k, I2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd;
    double DDDS[6*6*6];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao    = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_ddds__(
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,   AC,      DDDS );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		I2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<6; j++, jao++ ) {
		    if ( jao>iao ) { ix+=6; continue; }
		    IJ = I2 + jao;
		    for ( k=0, kao=kao0; k<6; k++, kao++, ix++ ) {
			KL = ((kao*kao+kao)>>1) + lao;
			if ( fabs(DDDS[ix]) > EPS_ERI ) {
			    V_frg[IJ] += TWO * D_mon[KL] * DDDS[ix];
			}
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_ddds_


/** (dd,dp)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_dddp__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, j, k, l, I2, K2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd;
    double DDDP[6*6*6*3];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao0   = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dddp__(
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,   AC,      DDDP );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		I2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<6; j++, jao++ ) {
		    if ( jao>iao ) { ix+=6*3; continue; }
		    IJ = I2 + jao;
		    for ( k=0, kao=kao0; k<6; k++, kao++ ) {
			K2 = (kao*kao+kao)>>1;
			for ( l=0, lao=lao0; l<3; l++, lao++, ix++ ) {
			    KL = K2 + lao;
			    if ( fabs(DDDP[ix]) > EPS_ERI ) {
				V_frg[IJ] += TWO * D_mon[KL] * DDDP[ix];
			    }
			}
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_dddp_


/** (dd,dd)タイプの４中心クーロンポテンシャル項をまとめて計算する
 * @ingroup integ-ifc4c
 * */
int ofmo_ifc4c_dddd__(
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
	const double atom_z_mon[],
	const int leading_cs_pair_mon[],
	const double csp_schwarz_mon[],
	const int csp_ics_mon[], const int csp_jcs_mon[],
	const int csp_leading_ps_pair_mon[],
	const double psp_zeta_mon[], const double psp_dkps_mon[],
	const double psp_xiza_mon[],
	// density matrix of monomer
	const double D_mon[],
	// (output) Coulomb potential
	double V_frg[] ) {
	    int nworkers=*pnworkers, workerid=*pworkerid;
	    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int Lab, Lcd, IJ, KL, i, j, k, l, I2, K2, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    double DDDD[6*6*6*6];
    //
    int ixx, ncsp, dx, res, pos;
    Lab = La * (La+1)/2 + Lb;
    Lcd = Lc * (Lc+1)/2 + Ld;
    if ( nworkers < 0 ) {
	pos   = workerid;
	ncsp  = leading_cs_pair_frg[Lab+1] - leading_cs_pair_frg[Lab];
	dx    = (ncsp>>7);	// >>5 means /32
	res   = (ncsp & 0x007f);// residues in division by 32
	ijcs0 = leading_cs_pair_frg[Lab]
	    + (pos<res ? pos*(dx+1) : pos*dx+res );
	ijcs1 = ijcs0 + ( pos<res ? dx+1 : dx );
	ixx   = 1;
    } else {
	ijcs0 = leading_cs_pair_frg[Lab] + workerid;
	ijcs1 = leading_cs_pair_frg[Lab+1];
	ixx   = nworkers;
    }
    klcs0 = leading_cs_pair_mon[Lcd];
    klcs1 = leading_cs_pair_mon[Lcd+1];
    for ( ijcs=ijcs0; ijcs<ijcs1; ijcs+=ixx ) {
	val_ab = csp_schwarz_frg[ijcs];
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
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];

//#pragma omp for schedule(guided)
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz_mon[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics_mon[klcs];
	    lcs    = csp_jcs_mon[klcs];
	    klps0  = csp_leading_ps_pair_mon[klcs];
	    nklps  = csp_leading_ps_pair_mon[klcs+1]-klps0;
	    kat    = shel_atm_mon[kcs];
	    lat    = shel_atm_mon[lcs];
	    kao0   = shel_ini_mon[kcs];
	    lao0   = shel_ini_mon[lcs];
	    C[0]=atom_x_mon[kat]; C[1]=atom_y_mon[kat]; C[2]=atom_z_mon[kat];
	    D[0]=atom_x_mon[lat]; D[1]=atom_y_mon[lat]; D[2]=atom_z_mon[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dddd__(
		    &nijps, &psp_zeta_frg[ijps0], &psp_dkps_frg[ijps0],
		    &psp_xiza_frg[ijps0], BA,
		    &nklps, &psp_zeta_mon[klps0], &psp_dkps_mon[klps0],
		    &psp_xiza_mon[klps0], DC,   AC,      DDDD );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		I2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<6; j++, jao++ ) {
		    if ( jao>iao ) { ix+=6*6; continue; }
		    IJ = I2 + jao;
		    for ( k=0, kao=kao0; k<6; k++, kao++ ) {
			K2 = (kao*kao+kao)>>1;
			for ( l=0, lao=lao0; l<6; l++, lao++, ix++ ) {
			    if ( lao>kao ) continue;
			    coe = (kao==lao? ONE : TWO );
			    KL = K2 + lao;
			    if ( fabs(DDDD[ix]) > EPS_ERI ) {
				V_frg[IJ] += coe * D_mon[KL] * DDDD[ix];
			    }
			}
		    }
		}
	    }
	}	// klcs
    }		// ijcs
    return 0;
}	// end of ofmo_ifc4c_dddd_

