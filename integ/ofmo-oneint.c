/**
 * @file ofmo-oneint.c
 * 通常のHartree-Fock SCF計算で用いる１電子積分計算を行う関数群
 * */
/**
 * @defgroup integ-oneint １電子積分行う関数群
 *
 * 通常のHartree-Fock計算で用いる、重なり積分、および、１電子ハミルトン行列
 * の計算を行う関数群。各積分タイプの１電子積分をまとめて行う。
 *
 * すべての関数は、同じ引数をとるので、以下にそれらの内容を示す。
 *
 * @param[in] La CSペアのうち、１つめのCSの軌道量子数
 * @param[in] Lb CSペアのうち、２つ目のCSの軌道量子数
 *     （\f$ \tt{La} \ge \tt{Lb} \ge 0 \f$
 * @param[in] leading_cs[lqn] 軌道量子数 /c lqn の先頭CS番号
 * @param[in] shel_tem[ics] CS番号 \c ics のCSの縮約長
 * @param[in] shel_atm[ics] CS番号 \c ics のCSが属する原子の番号
 * @param[in] shel_add[ics] CS番号 \c ics のCSに含まれるPSの先頭PS番号
 * @param[in] atom_x[iat] 原子の番号 \c iat のx座標（au単位）
 * @param[in] atom_y[iat] 原子の番号 \c iat のy座標（au単位）
 * @param[in] atom_z[iat] 原子の番号 \c iat のz座標（au単位）
 * @param[in] prim_exp[ips] PS番号 \c ips のPSの軌道指数
 * @param[in] prim_coe[ips] PS番号 \c ips のPSの規格化定数込みの縮約係数
 * @param[in] nat 原子数
 * @param[in] atomic_number[iat] 原子の番号 \c iat の原子の原子番号
 * @param[out] S[] 計算された重なり行列の要素を格納する配列
 *     （圧縮"U"形式）
 * @param[out] H[] 計算された１電子ハミルトン行列の要素を格納する配列
 *     （圧縮"U"形式）
 *
 * @ingroup integ-med
 *
 * */
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ofmo-oneint-core.h"

/** SSタイプの１電子積分を行う関数
 * @ingroup integ-oneint
 * */
int ofmo_oneint_ss__(
	const int *pnworkers, const int *pworkerid,
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[],
	const int shel_add[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	const int *pnat, const int atomic_number[],
	double S[], double H[] ) {
    int iao2, ijao;
    int ips0, ics, ics0, ics1, iat, iao;
    int jps0, jcs, jcs0      , jat, jao;
    int nps_i, nps_j;
    double A[3], B[3];
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, nat=*pnat;

    double oviss, hcoress;

/*#pragma omp master
    {
	printf("SS\n");
	fflush(stdout);
    }*/

    ics0 = jcs0 = leading_cs[La];
    ics1 = leading_cs[La+1];
    for ( ics=ics0+workerid; ics<ics1; ics+=nworkers ) {
	ips0  = shel_add[ics];
	iat   = shel_atm[ics];
	iao   = shel_ini[ics];
	nps_i = shel_tem[ics];
	A[0]=atom_x[ iat ]; A[1]=atom_y[ iat ]; A[2]=atom_z[ iat ];
	iao2 = iao*(iao+1)/2;
	for ( jcs=jcs0; jcs<=ics; jcs++ ) {
	    jps0  = shel_add[jcs];
	    jat   = shel_atm[jcs];
	    jao   = shel_ini[jcs];
	    nps_j = shel_tem[jcs];
	    B[0]=atom_x[ jat ]; B[1]=atom_y[ jat ]; B[2]=atom_z[ jat ];
	    ijao  = iao2+jao;
	    /*// debug
#pragma omp master
	    {
		printf("iao, jao= %d, %d\n", iao, jao );
		fflush(stdout);
	    }*/
	    oneint_core_ss__(
		    &ips0, &nps_i, A,  &jps0, &nps_j, B,  prim_exp, prim_coe,
		    &nat, atom_x, atom_y, atom_z, atomic_number,
		    &oviss, &hcoress );
	    S[ijao] = oviss;
	    H[ijao] = hcoress;
	}
    }
    return 0;
}
/** PSタイプの１電子積分を行う関数
 * @ingroup integ-oneint
 * */
int ofmo_oneint_ps__(
	const int *pnworkers, const int *pworkerid,
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[],
	const int shel_add[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	const int *pnat, const int atomic_number[],
	double S[], double H[] ) {
    int i, iao2, ijao;
    int ips0, ics, ics0, ics1, iat, iao, iao0;
    int jps0, jcs, jcs0, jcs1, jat, jao;
    int nps_i, nps_j;
    double A[3], B[3];
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, nat=*pnat;

    double ovips[3], hcoreps[3];

    ics0 = leading_cs[La];
    ics1 = leading_cs[La+1];
    jcs0 = leading_cs[Lb];
    jcs1 = leading_cs[Lb+1];
    for ( ics=ics0+workerid; ics<ics1; ics+=nworkers ) {
	ips0  = shel_add[ics];
	iat   = shel_atm[ics];
	iao0  = shel_ini[ics];
	nps_i = shel_tem[ics];
	A[0]=atom_x[ iat ]; A[1]=atom_y[ iat ]; A[2]=atom_z[ iat ];
	for ( jcs=jcs0; jcs<jcs1; jcs++ ) {
	    jps0  = shel_add[jcs];
	    jat   = shel_atm[jcs];
	    jao   = shel_ini[jcs];
	    nps_j = shel_tem[jcs];
	    B[0]=atom_x[ jat ]; B[1]=atom_y[ jat ]; B[2]=atom_z[ jat ];
	    oneint_core_ps__(
		    &ips0, &nps_i, A,  &jps0, &nps_j, B,  prim_exp, prim_coe,
		    &nat, atom_x, atom_y, atom_z, atomic_number,
		    ovips, hcoreps );
	    for ( i=0, iao=iao0; i<3; i++, iao++ ) {
		iao2 = (iao*iao+iao)>>1;
		ijao = iao2 + jao;
		S[ijao] = ovips[i];
		H[ijao] = hcoreps[i];
	    }
	}
    }
    return 0;
}
/** PPタイプの１電子積分を行う関数
 * @ingroup integ-oneint
 * */
int ofmo_oneint_pp__(
	const int *pnworkers, const int *pworkerid,
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[],
	const int shel_add[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	const int *pnat, const int atomic_number[],
	double S[], double H[] ) {
    int i, j, ij, ii, iao2, ijao;
    int ips0, ics, ics0, ics1, iat, iao, iao0;
    int jps0, jcs, jcs0      , jat, jao, jao0;
    int nps_i, nps_j;
    double A[3], B[3];
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, nat=*pnat;

    double ovipp[3*3], hcorepp[3*3];

    ics0 = jcs0 = leading_cs[La];
    ics1 = leading_cs[La+1];
    for ( ics=ics0+workerid; ics<ics1; ics+=nworkers ) {
	ips0  = shel_add[ics];
	iat   = shel_atm[ics];
	iao0  = shel_ini[ics];
	nps_i = shel_tem[ics];
	A[0]=atom_x[ iat ]; A[1]=atom_y[ iat ]; A[2]=atom_z[ iat ];
	for ( jcs=jcs0; jcs<=ics; jcs++ ) {
	    jps0  = shel_add[jcs];
	    jat   = shel_atm[jcs];
	    jao0  = shel_ini[jcs];
	    nps_j = shel_tem[jcs];
	    B[0]=atom_x[ jat ]; B[1]=atom_y[ jat ]; B[2]=atom_z[ jat ];
	    oneint_core_pp__(
		    &ips0, &nps_i, A,  &jps0, &nps_j, B,  prim_exp, prim_coe,
		    &nat, atom_x, atom_y, atom_z, atomic_number,
		    ovipp, hcorepp );
	    if ( ics != jcs ) {
		ij = 0;
		for ( i=0, iao=iao0; i<3; i++, iao++ ) {
		    iao2 = (iao*iao+iao)>>1;
		    for ( j=0, jao=jao0; j<3; j++,jao++ ) {
			ijao = iao2 + jao;
			S[ijao] = ovipp[ij];
			H[ijao] = hcorepp[ij];
			ij++;
		    }
		}
	    } else {
		for ( i=0, ii=0, iao=iao0; i<3; i++, iao++, ii+=3 ) {
		    iao2 = (iao*iao+iao)>>1;
		    for ( j=0, jao=jao0, ij=ii; j<=i; j++, jao++, ij++ ) {
			ijao = iao2 + jao;
			S[ijao] = ovipp[ij];
			H[ijao] = hcorepp[ij];
		    }
		}
	    }
	}
    }
    return 0;
}
/** DSタイプの１電子積分を行う関数
 * @ingroup integ-oneint
 * */
int ofmo_oneint_ds__(
	const int *pnworkers, const int *pworkerid,
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[],
	const int shel_add[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	const int *pnat, const int atomic_number[],
	double S[], double H[] ) {
    int i, iao2, ijao;
    int ips0, ics, ics0, ics1, iat, iao, iao0;
    int jps0, jcs, jcs0, jcs1, jat, jao;
    int nps_i, nps_j;
    double A[3], B[3];
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, nat=*pnat;

    double ovids[6], hcoreds[6];

    ics0 = leading_cs[La];
    ics1 = leading_cs[La+1];
    jcs0 = leading_cs[Lb];
    jcs1 = leading_cs[Lb+1];
    for ( ics=ics0+workerid; ics<ics1; ics+=nworkers ) {
	ips0  = shel_add[ics];
	iat   = shel_atm[ics];
	iao0  = shel_ini[ics];
	nps_i = shel_tem[ics];
	A[0]=atom_x[ iat ]; A[1]=atom_y[ iat ]; A[2]=atom_z[ iat ];
	for ( jcs=jcs0; jcs<jcs1; jcs++ ) {
	    jps0  = shel_add[jcs];
	    jat   = shel_atm[jcs];
	    jao   = shel_ini[jcs];
	    nps_j = shel_tem[jcs];
	    B[0]=atom_x[ jat ]; B[1]=atom_y[ jat ]; B[2]=atom_z[ jat ];
	    oneint_core_ds__(
		    &ips0, &nps_i, A,  &jps0, &nps_j, B,  prim_exp, prim_coe,
		    &nat, atom_x, atom_y, atom_z, atomic_number,
		    ovids, hcoreds );
	    for ( i=0, iao=iao0; i<6; i++, iao++ ) {
		iao2 = (iao*iao+iao)>>1;
		ijao = iao2 + jao;
		S[ijao] = ovids[i];
		H[ijao] = hcoreds[i];
	    }
	}
    }
    return 0;
}
/** DPタイプの１電子積分を行う関数
 * @ingroup integ-oneint
 * */
int ofmo_oneint_dp__(
	const int *pnworkers, const int *pworkerid,
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[],
	const int shel_add[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	const int *pnat, const int atomic_number[],
	double S[], double H[] ) {
    int i, j, ij, iao2, ijao;
    int ips0, ics, ics0, ics1, iat, iao, iao0;
    int jps0, jcs, jcs0, jcs1, jat, jao, jao0;
    int nps_i, nps_j;
    double A[3], B[3];
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, nat=*pnat;

    double ovidp[6*3], hcoredp[6*3];

    ics0 = leading_cs[La];
    ics1 = leading_cs[La+1];
    jcs0 = leading_cs[Lb];
    jcs1 = leading_cs[Lb+1];
    for ( ics=ics0+workerid; ics<ics1; ics+=nworkers ) {
	ips0  = shel_add[ics];
	iat   = shel_atm[ics];
	iao0  = shel_ini[ics];
	nps_i = shel_tem[ics];
	A[0]=atom_x[ iat ]; A[1]=atom_y[ iat ]; A[2]=atom_z[ iat ];
	for ( jcs=jcs0; jcs<jcs1; jcs++ ) {
	    jps0  = shel_add[jcs];
	    jat   = shel_atm[jcs];
	    jao0  = shel_ini[jcs];
	    nps_j = shel_tem[jcs];
	    B[0]=atom_x[ jat ]; B[1]=atom_y[ jat ]; B[2]=atom_z[ jat ];
	    oneint_core_dp__(
		    &ips0, &nps_i, A,  &jps0, &nps_j, B,  prim_exp, prim_coe,
		    &nat, atom_x, atom_y, atom_z, atomic_number,
		    ovidp, hcoredp );
	    ij = 0;
	    for ( i=0, iao=iao0; i<6; i++, iao++ ) {
		iao2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<3; j++, jao++ ) {
		    ijao = iao2 + jao;
		    S[ijao] = ovidp[ij];
		    H[ijao] = hcoredp[ij];
		    ij++;
		}
	    }
	}
    }
    return 0;
}

/** DDタイプの１電子積分を行う関数
 * @ingroup integ-oneint
 * */
int ofmo_oneint_dd__(
	const int *pnworkers, const int *pworkerid,
	const int *pLa, const int *pLb, const int leading_cs[],
	const int shel_tem[], const int shel_atm[],
	const int shel_add[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	const double prim_exp[], const double prim_coe[],
	const int *pnat, const int atomic_number[],
	double S[], double H[] ) {
    int i, j, ij, ii, iao2, ijao;
    int ips0, ics, ics0, ics1, iat, iao, iao0;
    int jps0, jcs, jcs0      , jat, jao, jao0;
    int nps_i, nps_j;
    double A[3], B[3];
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, nat=*pnat;

    double ovidd[6*6], hcoredd[6*6];

    ics0 = jcs0 = leading_cs[La];
    ics1 = leading_cs[La+1];
    for ( ics=ics0+workerid; ics<ics1; ics+=nworkers ) {
	ips0  = shel_add[ics];
	iat   = shel_atm[ics];
	iao0  = shel_ini[ics];
	nps_i = shel_tem[ics];
	A[0]=atom_x[ iat ]; A[1]=atom_y[ iat ]; A[2]=atom_z[ iat ];
	for ( jcs=jcs0; jcs<=ics; jcs++ ) {
	    jps0  = shel_add[jcs];
	    jat   = shel_atm[jcs];
	    jao0  = shel_ini[jcs];
	    nps_j = shel_tem[jcs];
	    B[0]=atom_x[ jat ]; B[1]=atom_y[ jat ]; B[2]=atom_z[ jat ];
	    oneint_core_dd__(
		    &ips0, &nps_i, A,  &jps0, &nps_j, B,  prim_exp, prim_coe,
		    &nat, atom_x, atom_y, atom_z, atomic_number,
		    ovidd, hcoredd );
	    if ( ics != jcs ) {
		ij = 0;
		for ( i=0, iao=iao0; i<6; i++, iao++ ) {
		    iao2 = (iao*iao+iao)>>1;
		    for ( j=0, jao=jao0; j<6; j++,jao++ ) {
			ijao = iao2 + jao;
			S[ijao] = ovidd[ij];
			H[ijao] = hcoredd[ij];
			ij++;
		    }
		}
	    } else {
		for ( i=0, ii=0, iao=iao0; i<6; i++, iao++, ii+=6 ) {
		    iao2 = (iao*iao+iao)>>1;
		    for ( j=0, jao=jao0, ij=ii; j<=i; j++, jao++, ij++ ) {
			ijao = iao2 + jao;
			S[ijao] = ovidd[ij];
			H[ijao] = hcoredd[ij];
		    }
		}
	    }
	}
    }
    return 0;
}

