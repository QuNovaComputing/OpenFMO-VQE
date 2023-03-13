/**
 * @file ofmo-ifc2c.c
 * ２中心クーロン相互作用項を計算する中位の関数群を
 * 記述している
 */

/**
 * @defgroup integ-ifc2c ２中心クーロン積分を行う関数群
 *
 * ２中心クーロン相互作用項を計算する中位の関数群を
 * 記述している。
 * 中位の関数では、各タイプの積分をまとめて計算する。
 *
 * すべての関数は、同じ引数をとるので、以下にそれらの内容を示す。
 *
 * @param[in] nworkers 計算に用いるワーカプロセス（スレッド）数
 * @param[in] workerid 各ワーカプロセス（スレッド）のID。
 *     \f$ 0\le\tt{workerid}<\tt{nworkers} \f$である。
 * @param[in] La CSペアのうち、１つめのCSの軌道量子数
 * @param[in] Lb CSペアのうち、２つ目のCSの軌道量子数
 *     （\f$ \tt{La} \ge \tt{Lb} \ge 0 \f$
 * @param[in] leading_cs_frg[lqn] 対象フラグメントの、
 *     軌道量子数 /c lqn の先頭CS番号
 * @param[in] shel_tem_frg[ics] 対象フラグメントの、
 *     CS番号 \c ics のCSの縮約長
 * @param[in] shel_atm_frg[ics] 対象フラグメントの、
 *     CS番号 \c ics のCSが属する原子の番号
 * @param[in] shel_add_frg[ics] 対象フラグメントの、
 *     CS番号 \c ics のCSに含まれるPSの先頭PS番号
 * @param[in] shel_ini_frg[ics] 対象フラグメントの、
 *     CS番号 \c ics のCSの先頭AO番号
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
 * @param[in] nat_mon 相手モノマーの、原子数
 * @param[in] atom_x_mon[iat] 相手モノマーの、
 *     原子の番号 \c iat のx座標（au単位）
 * @param[in] atom_y_mon[iat] 相手モノマーの、
 *     原子の番号 \c iat のy座標（au単位）
 * @param[in] atom_z_mon[iat] 相手モノマーの、
 *     原子の番号 \c iat のz座標（au単位）
 * @param[in] atm_pop_mon[iat] 相手モノマーの、
 *     原子の番号 \c iat の原子の原子番号
 *
 * @param[out] V_frg[] 計算された２中心クーロン相互作用項の要素を
 *     代入する配列（圧縮"U"形式）
 *
 * @ingroup integ-med
 * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ofmo-ifc2c-core.h"
/** SSタイプの２中心相互作用項を計算する関数
 * @ingroup integ-ifc2c
 * */
int ofmo_ifc2c_ss__(
	// parallelization
	const int *pnworkers, const int *pworkerid,
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
	const int *pnat_mon, const double atom_x_mon[],
	const double atom_y_mon[], const double atom_z_mon[],
	const double atm_pop_mon[],
	// output
	double V_frg[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, nat_mon=*pnat_mon;
    int iao2, ijao;
    int ips0, ics, ics0, ics1, iat, iao;
    int jps0, jcs, jcs0      , jat, jao;
    int nps_i, nps_j;
    double A[3], B[3];
    double SS[1];
    
    ics0 = leading_cs_frg[La];
    ics1 = leading_cs_frg[La+1];
    jcs0 = leading_cs_frg[Lb];
    for ( ics=ics0+workerid; ics<ics1; ics+=nworkers ) {
	ips0  = shel_add_frg[ics];
	iat   = shel_atm_frg[ics];
	iao   = shel_ini_frg[ics];
	nps_i = shel_tem_frg[ics];
	A[0] = atom_x_frg[ iat ];
	A[1] = atom_y_frg[ iat ];
	A[2] = atom_z_frg[ iat ];
	iao2 = iao*(iao+1)/2;
	for ( jcs=jcs0; jcs<=ics; jcs++ ) {
	    jps0  = shel_add_frg[jcs];
	    jat   = shel_atm_frg[jcs];
	    jao   = shel_ini_frg[jcs];
	    nps_j = shel_tem_frg[jcs];
	    B[0] = atom_x_frg[ jat ];
	    B[1] = atom_y_frg[ jat ];
	    B[2] = atom_z_frg[ jat ];
	    ijao  = iao2+jao;
	    ifc2c_core_ss__(
		    &ips0, &nps_i, A,  &jps0, &nps_j, B,
		    prim_exp_frg, prim_coe_frg,
		    &nat_mon, atom_x_mon, atom_y_mon, atom_z_mon,
		    atm_pop_mon,    SS );
	    V_frg[ijao] += SS[0];
	}
    }
    return 0;
}

/** PSタイプの２中心相互作用項を計算する関数
 * @ingroup integ-ifc2c
 * */
int ofmo_ifc2c_ps__(
	// parallelization
	const int *pnworkers, const int *pworkerid,
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
	const int *pnat_mon, const double atom_x_mon[],
	const double atom_y_mon[], const double atom_z_mon[],
	const double atm_pop_mon[],
	// output
	double V_frg[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, nat_mon=*pnat_mon;
    int i, ijao;
    int ips0, ics, ics0, ics1, iat, iao, iao0;
    int jps0, jcs, jcs0, jcs1, jat, jao;
    int nps_i, nps_j;
    double A[3], B[3];
    double PS[3];
    
    ics0 = leading_cs_frg[La];
    ics1 = leading_cs_frg[La+1];
    jcs0 = leading_cs_frg[Lb];
    jcs1 = leading_cs_frg[Lb+1];
    for ( ics=ics0+workerid; ics<ics1; ics+=nworkers ) {
	ips0  = shel_add_frg[ics];
	iat   = shel_atm_frg[ics];
	iao0  = shel_ini_frg[ics];
	nps_i = shel_tem_frg[ics];
	A[0] = atom_x_frg[ iat ];
	A[1] = atom_y_frg[ iat ];
	A[2] = atom_z_frg[ iat ];
	for ( jcs=jcs0; jcs<jcs1; jcs++ ) {
	    jps0  = shel_add_frg[jcs];
	    jat   = shel_atm_frg[jcs];
	    jao   = shel_ini_frg[jcs];
	    nps_j = shel_tem_frg[jcs];
	    B[0] = atom_x_frg[ jat ];
	    B[1] = atom_y_frg[ jat ];
	    B[2] = atom_z_frg[ jat ];
	    ifc2c_core_ps__(
		    &ips0, &nps_i, A,  &jps0, &nps_j, B,
		    prim_exp_frg, prim_coe_frg,
		    &nat_mon, atom_x_mon, atom_y_mon, atom_z_mon,
		    atm_pop_mon,    PS );
	    for ( i=0, iao=iao0; i<3; i++, iao++ ) {
		ijao = ((iao*iao+iao)>>1) + jao;
		V_frg[ijao] += PS[i];
	    }
	}
    }
    return 0;
}

/** PPタイプの２中心相互作用項を計算する関数
 * @ingroup integ-ifc2c
 * */
int ofmo_ifc2c_pp__(
	// parallelization
	const int *pnworkers, const int *pworkerid,
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
	const int *pnat_mon, const double atom_x_mon[],
	const double atom_y_mon[], const double atom_z_mon[],
	const double atm_pop_mon[],
	// output
	double V_frg[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, nat_mon=*pnat_mon;
    int i, j, iao2, ijao, ix;
    int ips0, ics, ics0, ics1, iat, iao, iao0;
    int jps0, jcs, jcs0      , jat, jao, jao0;
    int nps_i, nps_j;
    double A[3], B[3];
    double PP[3*3];
    
    ics0 = leading_cs_frg[La];
    ics1 = leading_cs_frg[La+1];
    jcs0 = leading_cs_frg[Lb];
    for ( ics=ics0+workerid; ics<ics1; ics+=nworkers ) {
	ips0  = shel_add_frg[ics];
	iat   = shel_atm_frg[ics];
	iao0  = shel_ini_frg[ics];
	nps_i = shel_tem_frg[ics];
	A[0] = atom_x_frg[ iat ];
	A[1] = atom_y_frg[ iat ];
	A[2] = atom_z_frg[ iat ];
	for ( jcs=jcs0; jcs<=ics; jcs++ ) {
	    jps0  = shel_add_frg[jcs];
	    jat   = shel_atm_frg[jcs];
	    jao0  = shel_ini_frg[jcs];
	    nps_j = shel_tem_frg[jcs];
	    B[0] = atom_x_frg[ jat ];
	    B[1] = atom_y_frg[ jat ];
	    B[2] = atom_z_frg[ jat ];
	    ifc2c_core_pp__(
		    &ips0, &nps_i, A,  &jps0, &nps_j, B,
		    prim_exp_frg, prim_coe_frg,
		    &nat_mon, atom_x_mon, atom_y_mon, atom_z_mon,
		    atm_pop_mon,    PP );
	    for ( i=0, iao=iao0, ix=0; i<3; i++, iao++ ) {
		iao2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<3; j++, jao++, ix++ ) {
		    if ( jao>iao ) continue;
		    ijao = iao2 + jao;
		    V_frg[ijao] += PP[ix];
		}
	    }
	}
    }
    return 0;
}

/** DSタイプの２中心相互作用項を計算する関数
 * @ingroup integ-ifc2c
 * */
int ofmo_ifc2c_ds__(
	// parallelization
	const int *pnworkers, const int *pworkerid,
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
	const int *pnat_mon, const double atom_x_mon[],
	const double atom_y_mon[], const double atom_z_mon[],
	const double atm_pop_mon[],
	// output
	double V_frg[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, nat_mon=*pnat_mon;
    int i, ijao;
    int ips0, ics, ics0, ics1, iat, iao, iao0;
    int jps0, jcs, jcs0, jcs1, jat, jao;
    int nps_i, nps_j;
    double A[3], B[3];
    double DS[6];
    
    ics0 = leading_cs_frg[La];
    ics1 = leading_cs_frg[La+1];
    jcs0 = leading_cs_frg[Lb];
    jcs1 = leading_cs_frg[Lb+1];
    for ( ics=ics0+workerid; ics<ics1; ics+=nworkers ) {
	ips0  = shel_add_frg[ics];
	iat   = shel_atm_frg[ics];
	iao0  = shel_ini_frg[ics];
	nps_i = shel_tem_frg[ics];
	A[0] = atom_x_frg[ iat ];
	A[1] = atom_y_frg[ iat ];
	A[2] = atom_z_frg[ iat ];
	for ( jcs=jcs0; jcs<jcs1; jcs++ ) {
	    jps0  = shel_add_frg[jcs];
	    jat   = shel_atm_frg[jcs];
	    jao   = shel_ini_frg[jcs];
	    nps_j = shel_tem_frg[jcs];
	    B[0] = atom_x_frg[ jat ];
	    B[1] = atom_y_frg[ jat ];
	    B[2] = atom_z_frg[ jat ];
	    ifc2c_core_ds__(
		    &ips0, &nps_i, A,  &jps0, &nps_j, B,
		    prim_exp_frg, prim_coe_frg,
		    &nat_mon, atom_x_mon, atom_y_mon, atom_z_mon,
		    atm_pop_mon,    DS );
	    for ( i=0, iao=iao0; i<6; i++, iao++ ) {
		ijao = ((iao*iao+iao)>>1) + jao;
		V_frg[ijao] += DS[i];
	    }
	}
    }
    return 0;
}

/** DPタイプの２中心相互作用項を計算する関数
 * @ingroup integ-ifc2c
 * */
int ofmo_ifc2c_dp__(
	// parallelization
	const int *pnworkers, const int *pworkerid,
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
	const int *pnat_mon, const double atom_x_mon[],
	const double atom_y_mon[], const double atom_z_mon[],
	const double atm_pop_mon[],
	// output
	double V_frg[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, nat_mon=*pnat_mon;
    int i, j, iao2, ijao, ix;
    int ips0, ics, ics0, ics1, iat, iao, iao0;
    int jps0, jcs, jcs0, jcs1, jat, jao, jao0;
    int nps_i, nps_j;
    double A[3], B[3];
    double DP[6*3];
    
    ics0 = leading_cs_frg[La];
    ics1 = leading_cs_frg[La+1];
    jcs0 = leading_cs_frg[Lb];
    jcs1 = leading_cs_frg[Lb+1];
    for ( ics=ics0+workerid; ics<ics1; ics+=nworkers ) {
	ips0  = shel_add_frg[ics];
	iat   = shel_atm_frg[ics];
	iao0  = shel_ini_frg[ics];
	nps_i = shel_tem_frg[ics];
	A[0] = atom_x_frg[ iat ];
	A[1] = atom_y_frg[ iat ];
	A[2] = atom_z_frg[ iat ];
	for ( jcs=jcs0; jcs<jcs1; jcs++ ) {
	    jps0  = shel_add_frg[jcs];
	    jat   = shel_atm_frg[jcs];
	    jao0  = shel_ini_frg[jcs];
	    nps_j = shel_tem_frg[jcs];
	    B[0] = atom_x_frg[ jat ];
	    B[1] = atom_y_frg[ jat ];
	    B[2] = atom_z_frg[ jat ];
	    ifc2c_core_dp__(
		    &ips0, &nps_i, A,  &jps0, &nps_j, B,
		    prim_exp_frg, prim_coe_frg,
		    &nat_mon, atom_x_mon, atom_y_mon, atom_z_mon,
		    atm_pop_mon,    DP );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		iao2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<3; j++, jao++, ix++ ) {
		    ijao = iao2 + jao;
		    V_frg[ijao] += DP[ix];
		}
	    }
	}
    }
    return 0;
}

/** DDタイプの２中心相互作用項を計算する関数
 * @ingroup integ-ifc2c
 * */
int ofmo_ifc2c_dd__(
	// parallelization
	const int *pnworkers, const int *pworkerid,
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
	const int *pnat_mon, const double atom_x_mon[],
	const double atom_y_mon[], const double atom_z_mon[],
	const double atm_pop_mon[],
	// output
	double V_frg[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, nat_mon=*pnat_mon;
    int i, j, iao2, ijao, ix;
    int ips0, ics, ics0, ics1, iat, iao, iao0;
    int jps0, jcs, jcs0      , jat, jao, jao0;
    int nps_i, nps_j;
    double A[3], B[3];
    double DD[6*6];
    
    ics0 = leading_cs_frg[La];
    ics1 = leading_cs_frg[La+1];
    jcs0 = leading_cs_frg[Lb];
    for ( ics=ics0+workerid; ics<ics1; ics+=nworkers ) {
	ips0  = shel_add_frg[ics];
	iat   = shel_atm_frg[ics];
	iao0  = shel_ini_frg[ics];
	nps_i = shel_tem_frg[ics];
	A[0] = atom_x_frg[ iat ];
	A[1] = atom_y_frg[ iat ];
	A[2] = atom_z_frg[ iat ];
	for ( jcs=jcs0; jcs<=ics; jcs++ ) {
	    jps0  = shel_add_frg[jcs];
	    jat   = shel_atm_frg[jcs];
	    jao0  = shel_ini_frg[jcs];
	    nps_j = shel_tem_frg[jcs];
	    B[0] = atom_x_frg[ jat ];
	    B[1] = atom_y_frg[ jat ];
	    B[2] = atom_z_frg[ jat ];
	    ifc2c_core_dd__(
		    &ips0, &nps_i, A,  &jps0, &nps_j, B,
		    prim_exp_frg, prim_coe_frg,
		    &nat_mon, atom_x_mon, atom_y_mon, atom_z_mon,
		    atm_pop_mon,    DD );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		iao2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<6; j++, jao++, ix++ ) {
		    if ( jao>iao ) continue;
		    ijao = iao2 + jao;
		    V_frg[ijao] += DD[ix];
		}
	    }
	}
    }
    return 0;
}
