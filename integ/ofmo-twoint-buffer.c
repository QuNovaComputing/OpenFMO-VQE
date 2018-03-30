/**
 * @file ofmo-twoint-buffer.c
 * 各タイプの二電子積分を計算して、バッファに保存する関数群。
 * */
/**
 * @defgroup integ-twoint-buffer ２電子積分を計算してバッファに保存する
 * buffered direct SCF法において、決められた容量を持つバッファに
 * 貯められるだけ２電子積分を計算して、保存する関数群。
 *
 * 各タイプ（2011/06/14現在、(ss,ss)〜(dd,dd)の21種類）の二電子積分
 * を計算して、計算した結果を対応する添字とともにバッファに保存する。
 * 基本的には、バッファがいっぱいになるまで計算するだけで、G行列生成は
 * 行わない。
 *
 * 各関数は、同じ引数を持っているので、以下にその内容を記す。
 *
 * @param[in] nworkers 計算に用いるワーカプロセス（スレッド）数
 * @param[in] workerid 各ワーカプロセス（スレッド）のID
 *     （\f$ 0\le\tt{workerid}<\tt{nworkers} \f$）
 * @param[in] ebuf_max_nzeri バッファに保存できる最大積分数
 * @param[in] La １つ目の軌道量子数
 * @param[in] Lb １つ目の軌道量子数
 *     （\f$ \tt{La} \ge \tt{Lb} \f$）
 * @param[in] Lc １つ目の軌道量子数
 *     （\f$ \tt{La} \ge \tt{Lc} \f$）
 * @param[in] Ld １つ目の軌道量子数
 *     （\f$ \tt{Lc} \ge \tt{Ld} \f$、かつ、
 *     \f$ \frac{\tt{La}(\tt{La}+1)}2+\tt{Lb} \ge
 *     \frac{\tt{Lc}(\tt{Lc}+1)}2+\tt{Ld} \f$）
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
 * @param[out] ebuf_val[] 計算して十分な大きさを持つSchwarz積分の値
 * @param[out] ebuf_ind4[] 保存したSchwarz積分に対応する4つの添字。
 *     <tt>short int</tt>の値であり、4つの添字がメモリ上の連続する領域に
 *     保存される。
 *
 * @param[in,out] *ebuf_non_zero_eri
 *     （入力時）この関数呼び出しまでに保存されたSchwarz積分
 *     総数。（出力時）この関数終了時点までに保存したSchwarz積分総数。
 * @param[out] *last_ijcs
 *     バッファが一杯になった時点での外側CSペア番号 \c ijcs の値。
 *     この関数での計算中にバッファが一杯になった場合のみ意味を持つ
 * @param[out] *last_klcs
 *     バッファが一杯になった時点での内側CSペア番号 \c klcs の値。
 *     この関数での計算中にバッファが一杯になった場合のみ意味を持つ。
 *
 * @retval OFMO_EBUF_FULL ２電子積分保存用のバッファが一杯になった場合
 * @retval OFMO_EBUF_NOFULL バッファが一杯にならなかった場合
 *
 * @ingroup integ-med
 * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define OFMO_EBUF_FULL		1
#define OFMO_EBUF_NOFULL	0

#define EPS_PS4 1.e-30
#define EPS_ERI	1.e-15
//#define EPS_PS4 1.e-15
//#define EPS_ERI	1.e-11
//#define EPS_PS4 1.e-14
//#define EPS_ERI	1.e-10

#define true  1
#define false 0

#define ZERO	0.e0
#define ONE	1.e0
#define HALF	.5e0

#include "ofmo-twoint-core.h"

/** (ss,ss)タイプの２電子積分関数
 * @ingroup integ-twoint-buffer
 * */
int ofmo_twoint_buffer_ssss__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs
	) {
    int Lab, Lcd, i;
    int ijcs, ijcs0, ijcs1, ijao;
    int klcs, klcs0,        klao;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, jcs, jat, jao, kcs, kat, kao, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    double SSSS[1];
    long nzeri, max_nzeri, nzeri4;
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    long ebuf_max_nzeri = *pebuf_max_nzeri;

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs0 = leading_cs_pair[Lab];
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];

    nzeri     = *ebuf_non_zero_eri;
    max_nzeri = ebuf_max_nzeri - (1-1);
    nzeri4    = nzeri*4;
    if ( nzeri >= max_nzeri ) { 
	*last_ijcs = ijcs0+workerid;
	*last_klcs = klcs0 - 1;
	*ebuf_non_zero_eri = nzeri;
	return OFMO_EBUF_FULL;
    }

    for ( ijcs=ijcs0+workerid; ijcs<ijcs1; ijcs+=nworkers ) {
	val_ab = csp_schwarz[ijcs];
	ics    = csp_ics[ijcs];
	jcs    = csp_jcs[ijcs];
	ijps0  = csp_leading_ps_pair[ijcs];
	nijps  = csp_leading_ps_pair[ijcs+1]-ijps0;
	iat    = shel_atm[ics];
	jat    = shel_atm[jcs];
	iao    = shel_ini[ics];
	jao    = shel_ini[jcs];
	ijao   = ((iao*iao+iao)>>1) + jao;
	A[0]=atom_x[iat]; A[1]=atom_y[iat]; A[2]=atom_z[iat];
	B[0]=atom_x[jat]; B[1]=atom_y[jat]; B[2]=atom_z[jat];
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];
	
	for ( klcs=klcs0; klcs<=ijcs; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
	    klps0  = csp_leading_ps_pair[klcs];
	    nklps  = csp_leading_ps_pair[klcs+1]-klps0;
	    kat    = shel_atm[kcs];
	    lat    = shel_atm[lcs];
	    kao    = shel_ini[kcs];
	    lao    = shel_ini[lcs];
	    klao   = ( (kao*kao+kao)>>1 ) + lao;
	    C[0]=atom_x[kat]; C[1]=atom_y[kat]; C[2]=atom_z[kat];
	    D[0]=atom_x[lat]; D[1]=atom_y[lat]; D[2]=atom_z[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_ssss__(
		    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
		    &psp_xiza[ijps0], BA,
		    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
		    &psp_xiza[klps0], DC,   AC,      SSSS );
	    if ( fabs(SSSS[0]) > EPS_ERI ) {
		coe = ONE;
		if ( iao == jao ) coe = HALF;
		if ( kao == lao ) coe *= HALF;
		if ( ijao == klao ) coe *= HALF;
		ebuf_val[nzeri]     = coe*SSSS[0];
		/*
		// for check sum
		lsum += ebuf_val[nzeri];
		*/

		ebuf_ind4[nzeri4+0] = (short int)iao;
		ebuf_ind4[nzeri4+1] = (short int)jao;
		ebuf_ind4[nzeri4+2] = (short int)kao;
		ebuf_ind4[nzeri4+3] = (short int)lao;
		nzeri++;
		nzeri4+=4;
		if ( nzeri >= max_nzeri ) {
		    *last_ijcs = ijcs;
		    *last_klcs = klcs;
		    *ebuf_non_zero_eri = nzeri;
		    return OFMO_EBUF_FULL;
		}
	    }	// if ( fabs(SSSS[0]) > EPS_ERI )
	}	// for ( klcs );
    }		// for ( ijcs );
    *ebuf_non_zero_eri = nzeri;
    /*
    // for check sum
#pragma omp atomic
    check_sum[0] += lsum;
    */

    return OFMO_EBUF_NOFULL;
}

/** (ps,ss)タイプの２電子積分関数
 * @ingroup integ-twoint-buffer
 * */
int ofmo_twoint_buffer_psss__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs
	) {
    int Lab, Lcd, i;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, kcs, kat, kao, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    double PSSS[3];
    long nzeri, max_nzeri, nzeri4;
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    long ebuf_max_nzeri = *pebuf_max_nzeri;
    /*
    // for check sum
    double lsum = ZERO;
    */

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs0 = leading_cs_pair[Lab];
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];

    nzeri     = *ebuf_non_zero_eri;
    max_nzeri = ebuf_max_nzeri - 3;
    nzeri4    = nzeri*4;
    if ( nzeri >= max_nzeri ) { 
	*last_ijcs = ijcs0+workerid;
	*last_klcs = klcs0 - 1;
	*ebuf_non_zero_eri = nzeri;
	return OFMO_EBUF_FULL;
    }
    /*
    // debug
    printf("(ps,ss) :  nwks=%d, wkid=%d, max nzeri=%lld\n",
	    nworkers, workerid, max_nzeri );
    fflush(stdout);
    */

    for ( ijcs=ijcs0+workerid; ijcs<ijcs1; ijcs+=nworkers ) {
	val_ab = csp_schwarz[ijcs];
	ics    = csp_ics[ijcs];
	jcs    = csp_jcs[ijcs];
	ijps0  = csp_leading_ps_pair[ijcs];
	nijps  = csp_leading_ps_pair[ijcs+1]-ijps0;
	iat    = shel_atm[ics];
	jat    = shel_atm[jcs];
	iao0   = shel_ini[ics];
	jao    = shel_ini[jcs];
	A[0]=atom_x[iat]; A[1]=atom_y[iat]; A[2]=atom_z[iat];
	B[0]=atom_x[jat]; B[1]=atom_y[jat]; B[2]=atom_z[jat];
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];
	
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
	    klps0  = csp_leading_ps_pair[klcs];
	    nklps  = csp_leading_ps_pair[klcs+1]-klps0;
	    kat    = shel_atm[kcs];
	    lat    = shel_atm[lcs];
	    kao    = shel_ini[kcs];
	    lao    = shel_ini[lcs];
	    C[0]=atom_x[kat]; C[1]=atom_y[kat]; C[2]=atom_z[kat];
	    D[0]=atom_x[lat]; D[1]=atom_y[lat]; D[2]=atom_z[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_psss__(
		    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
		    &psp_xiza[ijps0], BA,
		    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
		    &psp_xiza[klps0], DC,   AC,      PSSS );
	    coe = ( kao == lao ? HALF : ONE );
	    for ( i=0, iao=iao0; i<3; i++, iao++ ) {
		if ( fabs(PSSS[i]) > EPS_ERI ) {
		    ebuf_val[nzeri]     = coe*PSSS[i];
		    /*
		// for check sum
		lsum += ebuf_val[nzeri];
		*/

		    ebuf_ind4[nzeri4+0] = (short int)iao;
		    ebuf_ind4[nzeri4+1] = (short int)jao;
		    ebuf_ind4[nzeri4+2] = (short int)kao;
		    ebuf_ind4[nzeri4+3] = (short int)lao;
		    nzeri++;
		    nzeri4+=4;
		}	// if ( fabs(PSSS[i]) > EPS_ERI )
	    }
	    if ( nzeri >= max_nzeri ) {
		*last_ijcs = ijcs;
		*last_klcs = klcs;
		*ebuf_non_zero_eri = nzeri;
		return OFMO_EBUF_FULL;
	    }
	}	// for ( klcs );
    }		// for ( ijcs );
    *ebuf_non_zero_eri = nzeri;
    /*
    // for check sum
#pragma omp atomic
    check_sum[1] += lsum;
    */

    return OFMO_EBUF_NOFULL;
}

/** (ps,ps)タイプの２電子積分関数
 * @ingroup integ-twoint-buffer
 * */
int ofmo_twoint_buffer_psps__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs
	) {
    int Lab, Lcd, i, k, ix;
    int ipat, IJ, KL;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    double PSPS[3*3];
    long nzeri, max_nzeri, nzeri4;
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    long ebuf_max_nzeri = *pebuf_max_nzeri;
    /*
    // for check sum
    double lsum = ZERO;
    */

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs0 = leading_cs_pair[Lab];
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];

    nzeri     = *ebuf_non_zero_eri;
    max_nzeri = ebuf_max_nzeri - 3*3;
    nzeri4    = nzeri*4;
    if ( nzeri >= max_nzeri ) { 
	*last_ijcs = ijcs0+workerid;
	*last_klcs = klcs0 - 1;
	*ebuf_non_zero_eri = nzeri;
	return OFMO_EBUF_FULL;
    }
    /*
    // debug
    printf("(ps,ps) :  nwks=%d, wkid=%d, max nzeri=%lld\n",
	    nworkers, workerid, max_nzeri );
    fflush(stdout);
    */

    for ( ijcs=ijcs0+workerid; ijcs<ijcs1; ijcs+=nworkers ) {
	val_ab = csp_schwarz[ijcs];
	ics    = csp_ics[ijcs];
	jcs    = csp_jcs[ijcs];
	ijps0  = csp_leading_ps_pair[ijcs];
	nijps  = csp_leading_ps_pair[ijcs+1]-ijps0;
	iat    = shel_atm[ics];
	jat    = shel_atm[jcs];
	iao0   = shel_ini[ics];
	jao    = shel_ini[jcs];
	A[0]=atom_x[iat]; A[1]=atom_y[iat]; A[2]=atom_z[iat];
	B[0]=atom_x[jat]; B[1]=atom_y[jat]; B[2]=atom_z[jat];
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];
	
	for ( klcs=klcs0; klcs<=ijcs; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
	    klps0  = csp_leading_ps_pair[klcs];
	    nklps  = csp_leading_ps_pair[klcs+1]-klps0;
	    kat    = shel_atm[kcs];
	    lat    = shel_atm[lcs];
	    kao0   = shel_ini[kcs];
	    lao    = shel_ini[lcs];
	    C[0]=atom_x[kat]; C[1]=atom_y[kat]; C[2]=atom_z[kat];
	    D[0]=atom_x[lat]; D[1]=atom_y[lat]; D[2]=atom_z[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_psps__(
		    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
		    &psp_xiza[ijps0], BA,
		    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
		    &psp_xiza[klps0], DC,   AC,      PSPS );
	    /* だめな加算方法
	    for ( i=0, iao=iao0, ix=0; i<3; i++, iao++ ) {
		for ( k=0, kao=kao0; k<3; k++, kao++, ix++ ) {
		    if ( kao > iao ) continue;
		    if ( fabs(PSPS[ix]) > EPS_ERI ) {
			coe = ONE;
			if ( iao==kao && jao==lao ) coe = HALF;
			ebuf_val[nzeri]     = coe*PSPS[ix];
		// for check sum
		lsum += ebuf_val[nzeri];

			ebuf_ind4[nzeri4+0] = (short int)iao;
			ebuf_ind4[nzeri4+1] = (short int)jao;
			ebuf_ind4[nzeri4+2] = (short int)kao;
			ebuf_ind4[nzeri4+3] = (short int)lao;
			nzeri++;
			nzeri4+=4;
		    }
		}
	    }
	    */
	    ipat = ( (ics==kcs && jcs>lcs) ? true : false);
#ifdef SORT_CSP
            int ijgekl = (ics>kcs);
            if (ics==kcs) ijgekl = (jcs>=lcs);
            if (!ijgekl) ipat = ( (ics==kcs && jcs<lcs) ? true : false);
#endif
	    for ( i=0, iao=iao0, ix=0; i<3; i++, iao++ ) {
		IJ = ((iao*iao+iao)>>1) + jao;
		for ( k=0, kao=kao0; k<3; k++, kao++, ix++ ) {
		    KL = ((kao*kao+kao)>>1) + lao ;
		    if ( fabs(PSPS[ix]) <= EPS_ERI ) continue;
#ifndef SORT_CSP
		    if ( IJ >= KL ) {
#else
                    if ( (ijgekl&&IJ>=KL) || (!ijgekl&&KL>=IJ) ) {
#endif
			coe = ONE;
			if ( iao==kao && jao==lao ) coe = HALF;
			ebuf_val[nzeri]     = coe*PSPS[ix];
			/*
		// for check sum
		lsum += ebuf_val[nzeri];
		*/

			ebuf_ind4[nzeri4+0] = (short int)iao;
			ebuf_ind4[nzeri4+1] = (short int)jao;
			ebuf_ind4[nzeri4+2] = (short int)kao;
			ebuf_ind4[nzeri4+3] = (short int)lao;
			nzeri++;
			nzeri4+=4;
		    } else if ( ipat ) {
			ebuf_val[nzeri]     = PSPS[ix];
			/*
		// for check sum
		lsum += ebuf_val[nzeri];
		*/

			ebuf_ind4[nzeri4+0] = (short int)kao;
			ebuf_ind4[nzeri4+1] = (short int)lao;
			ebuf_ind4[nzeri4+2] = (short int)iao;
			ebuf_ind4[nzeri4+3] = (short int)jao;
			nzeri++;
			nzeri4+=4;
		    }
		}
	    }
	    if ( nzeri >= max_nzeri ) {
		*last_ijcs = ijcs;
		*last_klcs = klcs;
		*ebuf_non_zero_eri = nzeri;
		return OFMO_EBUF_FULL;
	    }
	}	// for ( klcs );
    }		// for ( ijcs );
    *ebuf_non_zero_eri = nzeri;
    /*
    // for check sum
#pragma omp atomic
    check_sum[2] += lsum;
    */

    return OFMO_EBUF_NOFULL;
}

/** (pp,ss)タイプの２電子積分関数
 * @ingroup integ-twoint-buffer
 * */
int ofmo_twoint_buffer_ppss__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs
	) {
    int Lab, Lcd, i, j, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe, coe0;
    double PPSS[3*3];
    long nzeri, max_nzeri, nzeri4;
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    long ebuf_max_nzeri = *pebuf_max_nzeri;
    /*
    // for check sum
    double lsum = ZERO;
    */

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs0 = leading_cs_pair[Lab];
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];

    nzeri     = *ebuf_non_zero_eri;
    max_nzeri = ebuf_max_nzeri - 3*3;
    nzeri4    = nzeri*4;
    if ( nzeri >= max_nzeri ) { 
	*last_ijcs = ijcs0+workerid;
	*last_klcs = klcs0 - 1;
	*ebuf_non_zero_eri = nzeri;
	return OFMO_EBUF_FULL;
    }
    /*
    // debug
    printf("(pp,ss) :  nwks=%d, wkid=%d, max nzeri=%lld\n",
	    nworkers, workerid, max_nzeri );
    fflush(stdout);
    */

    for ( ijcs=ijcs0+workerid; ijcs<ijcs1; ijcs+=nworkers ) {
	val_ab = csp_schwarz[ijcs];
	ics    = csp_ics[ijcs];
	jcs    = csp_jcs[ijcs];
	ijps0  = csp_leading_ps_pair[ijcs];
	nijps  = csp_leading_ps_pair[ijcs+1]-ijps0;
	iat    = shel_atm[ics];
	jat    = shel_atm[jcs];
	iao0   = shel_ini[ics];
	jao0   = shel_ini[jcs];
	A[0]=atom_x[iat]; A[1]=atom_y[iat]; A[2]=atom_z[iat];
	B[0]=atom_x[jat]; B[1]=atom_y[jat]; B[2]=atom_z[jat];
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];
	
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
	    klps0  = csp_leading_ps_pair[klcs];
	    nklps  = csp_leading_ps_pair[klcs+1]-klps0;
	    kat    = shel_atm[kcs];
	    lat    = shel_atm[lcs];
	    kao    = shel_ini[kcs];
	    lao    = shel_ini[lcs];
	    C[0]=atom_x[kat]; C[1]=atom_y[kat]; C[2]=atom_z[kat];
	    D[0]=atom_x[lat]; D[1]=atom_y[lat]; D[2]=atom_z[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_ppss__(
		    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
		    &psp_xiza[ijps0], BA,
		    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
		    &psp_xiza[klps0], DC,   AC,      PPSS );
	    coe0 = ( kao==lao ? HALF : ONE );
	    for ( i=0, iao=iao0, ix=0; i<3; i++, iao++ ) {
		for ( j=0, jao=jao0; j<3; j++, jao++, ix++ ) {
		    if ( jao > iao ) continue;
		    if ( fabs(PPSS[ix]) > EPS_ERI ) {
			coe = coe0;
			if ( iao == jao ) coe *= HALF;
			ebuf_val[nzeri]     = coe*PPSS[ix];
			/*
		// for check sum
		lsum += ebuf_val[nzeri];
		*/

			ebuf_ind4[nzeri4+0] = (short int)iao;
			ebuf_ind4[nzeri4+1] = (short int)jao;
			ebuf_ind4[nzeri4+2] = (short int)kao;
			ebuf_ind4[nzeri4+3] = (short int)lao;
			nzeri++;
			nzeri4+=4;
		    }	// if ( fabs(PPSS[ix]) > EPS_ERI )
		}
	    }
	    if ( nzeri >= max_nzeri ) {
		*last_ijcs = ijcs;
		*last_klcs = klcs;
		*ebuf_non_zero_eri = nzeri;
		return OFMO_EBUF_FULL;
	    }
	}	// for ( klcs );
    }		// for ( ijcs );
    *ebuf_non_zero_eri = nzeri;
    /*
    // for check sum
#pragma omp atomic
    check_sum[3] += lsum;
    */

    return OFMO_EBUF_NOFULL;
}

/** (pp,ps)タイプの２電子積分関数
 * @ingroup integ-twoint-buffer
 * */
int ofmo_twoint_buffer_ppps__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs
	) {
    int Lab, Lcd, i, j, k, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    double PPPS[3*3*3];
    long nzeri, max_nzeri, nzeri4;
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    long ebuf_max_nzeri = *pebuf_max_nzeri;
    /*
    // for check sum
    double lsum = ZERO;
    */

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs0 = leading_cs_pair[Lab];
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];

    nzeri     = *ebuf_non_zero_eri;
    max_nzeri = ebuf_max_nzeri - 3*3*3;
    nzeri4    = nzeri*4;
    if ( nzeri >= max_nzeri ) { 
	*last_ijcs = ijcs0+workerid;
	*last_klcs = klcs0 - 1;
	*ebuf_non_zero_eri = nzeri;
	return OFMO_EBUF_FULL;
    }
    /*
    // debug
    printf("(pp,ps) :  nwks=%d, wkid=%d, max nzeri=%lld\n",
	    nworkers, workerid, max_nzeri );
    fflush(stdout);
    */

    for ( ijcs=ijcs0+workerid; ijcs<ijcs1; ijcs+=nworkers ) {
	val_ab = csp_schwarz[ijcs];
	ics    = csp_ics[ijcs];
	jcs    = csp_jcs[ijcs];
	ijps0  = csp_leading_ps_pair[ijcs];
	nijps  = csp_leading_ps_pair[ijcs+1]-ijps0;
	iat    = shel_atm[ics];
	jat    = shel_atm[jcs];
	iao0   = shel_ini[ics];
	jao0   = shel_ini[jcs];
	A[0]=atom_x[iat]; A[1]=atom_y[iat]; A[2]=atom_z[iat];
	B[0]=atom_x[jat]; B[1]=atom_y[jat]; B[2]=atom_z[jat];
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];
	
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
	    klps0  = csp_leading_ps_pair[klcs];
	    nklps  = csp_leading_ps_pair[klcs+1]-klps0;
	    kat    = shel_atm[kcs];
	    lat    = shel_atm[lcs];
	    kao0   = shel_ini[kcs];
	    lao    = shel_ini[lcs];
	    C[0]=atom_x[kat]; C[1]=atom_y[kat]; C[2]=atom_z[kat];
	    D[0]=atom_x[lat]; D[1]=atom_y[lat]; D[2]=atom_z[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_ppps__(
		    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
		    &psp_xiza[ijps0], BA,
		    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
		    &psp_xiza[klps0], DC,   AC,      PPPS );
	    for ( i=0, iao=iao0, ix=0; i<3; i++, iao++ ) {
		for ( j=0, jao=jao0; j<3; j++, jao++ ) {
		    if ( jao>iao ) { ix += 3; continue; }
		    coe = ( iao==jao ? HALF : ONE );
		    for ( k=0, kao=kao0; k<3; k++, kao++, ix++ ) {
			if ( fabs(PPPS[ix]) > EPS_ERI ) {
			    ebuf_val[nzeri]     = coe*PPPS[ix];
			    /*
		// for check sum
		lsum += ebuf_val[nzeri];
		*/

			    if ( iao >= kao ) {
				ebuf_ind4[nzeri4+0] = (short int)iao;
				ebuf_ind4[nzeri4+1] = (short int)jao;
				ebuf_ind4[nzeri4+2] = (short int)kao;
				ebuf_ind4[nzeri4+3] = (short int)lao;
			    } else {
				ebuf_ind4[nzeri4+0] = (short int)kao;
				ebuf_ind4[nzeri4+1] = (short int)lao;
				ebuf_ind4[nzeri4+2] = (short int)iao;
				ebuf_ind4[nzeri4+3] = (short int)jao;
			    }
			    nzeri++;
			    nzeri4+=4;
			}	// if ( fabs(PPSS[ix]) > EPS_ERI )
		    }
		}
	    }
	    if ( nzeri >= max_nzeri ) {
		*last_ijcs = ijcs;
		*last_klcs = klcs;
		*ebuf_non_zero_eri = nzeri;
		return OFMO_EBUF_FULL;
	    }
	}	// for ( klcs );
    }		// for ( ijcs );
    *ebuf_non_zero_eri = nzeri;
    /*
    // for check sum
#pragma omp atomic
    check_sum[4] += lsum;
    */

    return OFMO_EBUF_NOFULL;
}

/** (pp,pp)タイプの２電子積分関数
 * @ingroup integ-twoint-buffer
 * */
int ofmo_twoint_buffer_pppp__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs
	) {
    int Lab, Lcd, i, j, k, l, ipat;
    int I2, IJ, K2, KL;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0 ;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe, coe0;
    double PPPP[3*3*3*3], pppp;
    long nzeri, max_nzeri, nzeri4;
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    long ebuf_max_nzeri = *pebuf_max_nzeri;
    /*
    // for check sum
    double lsum = ZERO;
    */

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs0 = leading_cs_pair[Lab];
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];

    nzeri     = *ebuf_non_zero_eri;
    max_nzeri = ebuf_max_nzeri - 3*3*3*3;
    nzeri4    = nzeri*4;
    if ( nzeri >= max_nzeri ) { 
	*last_ijcs = ijcs0+workerid;
	*last_klcs = klcs0 - 1;
	*ebuf_non_zero_eri = nzeri;
	return OFMO_EBUF_FULL;
    }
    /*
    // debug
    printf("(pp,pp) :  nwks=%d, wkid=%d, max nzeri=%lld\n",
	    nworkers, workerid, max_nzeri );
    fflush(stdout);
    */

    for ( ijcs=ijcs0+workerid; ijcs<ijcs1; ijcs+=nworkers ) {
	val_ab = csp_schwarz[ijcs];
	ics    = csp_ics[ijcs];
	jcs    = csp_jcs[ijcs];
	ijps0  = csp_leading_ps_pair[ijcs];
	nijps  = csp_leading_ps_pair[ijcs+1]-ijps0;
	iat    = shel_atm[ics];
	jat    = shel_atm[jcs];
	iao0   = shel_ini[ics];
	jao0   = shel_ini[jcs];
	A[0]=atom_x[iat]; A[1]=atom_y[iat]; A[2]=atom_z[iat];
	B[0]=atom_x[jat]; B[1]=atom_y[jat]; B[2]=atom_z[jat];
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];
	
	for ( klcs=klcs0; klcs<=ijcs; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
	    klps0  = csp_leading_ps_pair[klcs];
	    nklps  = csp_leading_ps_pair[klcs+1]-klps0;
	    kat    = shel_atm[kcs];
	    lat    = shel_atm[lcs];
	    kao0   = shel_ini[kcs];
	    lao0   = shel_ini[lcs];
	    C[0]=atom_x[kat]; C[1]=atom_y[kat]; C[2]=atom_z[kat];
	    D[0]=atom_x[lat]; D[1]=atom_y[lat]; D[2]=atom_z[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_pppp__(
		    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
		    &psp_xiza[ijps0], BA,
		    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
		    &psp_xiza[klps0], DC,   AC,      PPPP );
	    ipat = ( (ics==kcs && jcs>lcs) ? true : false );
#ifdef SORT_CSP
            int ijgekl = (ics>kcs);
            if (ics==kcs) ijgekl = (jcs>=lcs);
            if (!ijgekl) ipat = ( (ics==kcs && jcs<lcs) ? true : false);
#endif
	    for ( i=0, iao=iao0; i<3; i++, iao++ ) {
		I2 = iao*(iao+1)/2;
		for ( j=0, jao=jao0; j<3; j++, jao++ ) {
		    if ( jao>iao ) continue;
		    IJ = I2 + jao;
		    coe0 = ( iao==jao ? HALF : ONE );
		    for ( k=0, kao=kao0; k<3; k++, kao++ ) {
			K2 = kao*(kao+1)/2;
			for ( l=0, lao=lao0; l<3; l++, lao++ ) {
			    if ( lao>kao ) continue;
			    pppp = PPPP[i*27+j*9+k*3+l];
			    if ( fabs(pppp) > EPS_ERI ) {
				KL = K2 + lao;
#ifndef SORT_CSP
				if ( IJ >= KL ) {
#else
                    if ( (ijgekl&&IJ>=KL) || (!ijgekl&&KL>=IJ) ) {
#endif
				    coe = coe0;
				    if ( kao==lao ) coe *= HALF;
				    if ( KL == IJ ) coe *= HALF;
				    ebuf_val[nzeri]     = coe*pppp;
				    /*
		// for check sum
		lsum += ebuf_val[nzeri];
		*/

				    ebuf_ind4[nzeri4+0] = (short int)iao;
				    ebuf_ind4[nzeri4+1] = (short int)jao;
				    ebuf_ind4[nzeri4+2] = (short int)kao;
				    ebuf_ind4[nzeri4+3] = (short int)lao;
				    nzeri++;
				    nzeri4+=4;
				} else if ( ipat ) {
				    coe = coe0;
				    if ( kao==lao ) coe*=HALF;
				    ebuf_val[nzeri]     = coe*pppp;
				    /*
		// for check sum
		lsum += ebuf_val[nzeri];
		*/

				    ebuf_ind4[nzeri4+0] = (short int)kao;
				    ebuf_ind4[nzeri4+1] = (short int)lao;
				    ebuf_ind4[nzeri4+2] = (short int)iao;
				    ebuf_ind4[nzeri4+3] = (short int)jao;
				    nzeri++;
				    nzeri4+=4;
				}
			    }
			}	// l
		    }		// k
		}		// j
	    }			// i
	    if ( nzeri >= max_nzeri ) {
		*last_ijcs = ijcs;
		*last_klcs = klcs;
		*ebuf_non_zero_eri = nzeri;
		return OFMO_EBUF_FULL;
	    }
	}	// for ( klcs );
    }		// for ( ijcs );
    *ebuf_non_zero_eri = nzeri;
    /*
    // for check sum
#pragma omp atomic
    check_sum[5] += lsum;
    */

    return OFMO_EBUF_NOFULL;
}

/** (ds,ss)タイプの２電子積分関数
 * @ingroup integ-twoint-buffer
 * */
int ofmo_twoint_buffer_dsss__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs
	) {
    int Lab, Lcd, i;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, kcs, kat, kao, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    double DSSS[6];
    long nzeri, max_nzeri, nzeri4;
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    long ebuf_max_nzeri = *pebuf_max_nzeri;
    /*
    // for check sum
    double lsum = ZERO;
    */

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs0 = leading_cs_pair[Lab];
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];

    nzeri     = *ebuf_non_zero_eri;
    max_nzeri = ebuf_max_nzeri - 6;
    nzeri4    = nzeri*4;
    if ( nzeri >= max_nzeri ) { 
	*last_ijcs = ijcs0+workerid;
	*last_klcs = klcs0 - 1;
	*ebuf_non_zero_eri = nzeri;
	return OFMO_EBUF_FULL;
    }
    /*
    // debug
    printf("(ds,ss) :  nwks=%d, wkid=%d, max nzeri=%lld\n",
	    nworkers, workerid, max_nzeri );
    fflush(stdout);
    */

    for ( ijcs=ijcs0+workerid; ijcs<ijcs1; ijcs+=nworkers ) {
	val_ab = csp_schwarz[ijcs];
	ics    = csp_ics[ijcs];
	jcs    = csp_jcs[ijcs];
	ijps0  = csp_leading_ps_pair[ijcs];
	nijps  = csp_leading_ps_pair[ijcs+1]-ijps0;
	iat    = shel_atm[ics];
	jat    = shel_atm[jcs];
	iao0   = shel_ini[ics];
	jao    = shel_ini[jcs];
	A[0]=atom_x[iat]; A[1]=atom_y[iat]; A[2]=atom_z[iat];
	B[0]=atom_x[jat]; B[1]=atom_y[jat]; B[2]=atom_z[jat];
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];
	
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
	    klps0  = csp_leading_ps_pair[klcs];
	    nklps  = csp_leading_ps_pair[klcs+1]-klps0;
	    kat    = shel_atm[kcs];
	    lat    = shel_atm[lcs];
	    kao    = shel_ini[kcs];
	    lao    = shel_ini[lcs];
	    C[0]=atom_x[kat]; C[1]=atom_y[kat]; C[2]=atom_z[kat];
	    D[0]=atom_x[lat]; D[1]=atom_y[lat]; D[2]=atom_z[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dsss__(
		    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
		    &psp_xiza[ijps0], BA,
		    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
		    &psp_xiza[klps0], DC,   AC,      DSSS );
	    // debug
	    //if ( fabs(DSSS[0]) > 1.e-4 )
	    //printf("ijcs, klcs, dsss(0) = %4d %4d %16.12f\n", ijcs, klcs, DSSS[0] );

	    coe = ( kao == lao ? HALF : ONE );
	    for ( i=0, iao=iao0; i<6; i++, iao++ ) {
		if ( fabs(DSSS[i]) > EPS_ERI ) {
		    ebuf_val[nzeri]     = coe*DSSS[i];
		    /*
		// for check sum
		lsum += ebuf_val[nzeri];
		*/

		    ebuf_ind4[nzeri4+0] = (short int)iao;
		    ebuf_ind4[nzeri4+1] = (short int)jao;
		    ebuf_ind4[nzeri4+2] = (short int)kao;
		    ebuf_ind4[nzeri4+3] = (short int)lao;
		    nzeri++;
		    nzeri4+=4;
		}	// if ( fabs(DSSS[i]) > EPS_ERI )
	    }
	    if ( nzeri >= max_nzeri ) {
		*last_ijcs = ijcs;
		*last_klcs = klcs;
		*ebuf_non_zero_eri = nzeri;
		return OFMO_EBUF_FULL;
	    }
	}	// for ( klcs );
    }		// for ( ijcs );
    *ebuf_non_zero_eri = nzeri;
    /*
    // for check sum
#pragma omp atomic
    check_sum[6] += lsum;
    */

    return OFMO_EBUF_NOFULL;
}

/** (ds,ps)タイプの２電子積分関数
 * @ingroup integ-twoint-buffer
 * */
int ofmo_twoint_buffer_dsps__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs
	) {
    int Lab, Lcd, i, k, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd;
    double DSPS[6*3];
    long nzeri, max_nzeri, nzeri4;
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    long ebuf_max_nzeri = *pebuf_max_nzeri;
    /*
    // for check sum
    double lsum = ZERO;
    */

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs0 = leading_cs_pair[Lab];
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];

    nzeri     = *ebuf_non_zero_eri;
    max_nzeri = ebuf_max_nzeri - 6*3;
    nzeri4    = nzeri*4;
    if ( nzeri >= max_nzeri ) { 
	*last_ijcs = ijcs0+workerid;
	*last_klcs = klcs0 - 1;
	*ebuf_non_zero_eri = nzeri;
	return OFMO_EBUF_FULL;
    }
    /*
    // debug
    printf("(ds,ps) :  nwks=%d, wkid=%d, max nzeri=%lld\n",
	    nworkers, workerid, max_nzeri );
    fflush(stdout);
    */

    for ( ijcs=ijcs0+workerid; ijcs<ijcs1; ijcs+=nworkers ) {
	val_ab = csp_schwarz[ijcs];
	ics    = csp_ics[ijcs];
	jcs    = csp_jcs[ijcs];
	ijps0  = csp_leading_ps_pair[ijcs];
	nijps  = csp_leading_ps_pair[ijcs+1]-ijps0;
	iat    = shel_atm[ics];
	jat    = shel_atm[jcs];
	iao0   = shel_ini[ics];
	jao    = shel_ini[jcs];
	A[0]=atom_x[iat]; A[1]=atom_y[iat]; A[2]=atom_z[iat];
	B[0]=atom_x[jat]; B[1]=atom_y[jat]; B[2]=atom_z[jat];
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];
	
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
	    klps0  = csp_leading_ps_pair[klcs];
	    nklps  = csp_leading_ps_pair[klcs+1]-klps0;
	    kat    = shel_atm[kcs];
	    lat    = shel_atm[lcs];
	    kao0   = shel_ini[kcs];
	    lao    = shel_ini[lcs];
	    C[0]=atom_x[kat]; C[1]=atom_y[kat]; C[2]=atom_z[kat];
	    D[0]=atom_x[lat]; D[1]=atom_y[lat]; D[2]=atom_z[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dsps__(
		    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
		    &psp_xiza[ijps0], BA,
		    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
		    &psp_xiza[klps0], DC,   AC,      DSPS );

	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		for ( k=0, kao=kao0; k<3; k++, kao++, ix++ ) {
		    if ( fabs(DSPS[ix]) > EPS_ERI ) {
			ebuf_val[nzeri]     = DSPS[ix];
			/*
		// for check sum
		lsum += ebuf_val[nzeri];
		*/

			ebuf_ind4[nzeri4+0] = (short int)iao;
			ebuf_ind4[nzeri4+1] = (short int)jao;
			ebuf_ind4[nzeri4+2] = (short int)kao;
			ebuf_ind4[nzeri4+3] = (short int)lao;
			nzeri++;
			nzeri4+=4;
		    }	// if ( fabs(DSPS[i]) > EPS_ERI )
		}
	    }
	    if ( nzeri >= max_nzeri ) {
		*last_ijcs = ijcs;
		*last_klcs = klcs;
		*ebuf_non_zero_eri = nzeri;
		return OFMO_EBUF_FULL;
	    }
	}	// for ( klcs );
    }		// for ( ijcs );
    *ebuf_non_zero_eri = nzeri;
    /*
    // for check sum
#pragma omp atomic
    check_sum[7] += lsum;
    */

    return OFMO_EBUF_NOFULL;
}

/** (ds,pp)タイプの２電子積分関数
 * @ingroup integ-twoint-buffer
 * */
int ofmo_twoint_buffer_dspp__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs
	) {
    int Lab, Lcd, i, k, l, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    double DSPP[6*3*3];
    long nzeri, max_nzeri, nzeri4;
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    long ebuf_max_nzeri = *pebuf_max_nzeri;
    Lab = La*(La+1)/2+Lb;
    /*
    // for check sum
    double lsum = ZERO;
    */

    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs0 = leading_cs_pair[Lab];
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];

    nzeri     = *ebuf_non_zero_eri;
    max_nzeri = ebuf_max_nzeri - 6*3*3;
    nzeri4    = nzeri*4;
    if ( nzeri >= max_nzeri ) { 
	*last_ijcs = ijcs0+workerid;
	*last_klcs = klcs0 - 1;
	*ebuf_non_zero_eri = nzeri;
	return OFMO_EBUF_FULL;
    }
    /*
    // debug
    printf("(ds,pp) :  nwks=%d, wkid=%d, max nzeri=%lld\n",
	    nworkers, workerid, max_nzeri );
    fflush(stdout);
    */

    for ( ijcs=ijcs0+workerid; ijcs<ijcs1; ijcs+=nworkers ) {
	val_ab = csp_schwarz[ijcs];
	ics    = csp_ics[ijcs];
	jcs    = csp_jcs[ijcs];
	ijps0  = csp_leading_ps_pair[ijcs];
	nijps  = csp_leading_ps_pair[ijcs+1]-ijps0;
	iat    = shel_atm[ics];
	jat    = shel_atm[jcs];
	iao0   = shel_ini[ics];
	jao    = shel_ini[jcs];
	A[0]=atom_x[iat]; A[1]=atom_y[iat]; A[2]=atom_z[iat];
	B[0]=atom_x[jat]; B[1]=atom_y[jat]; B[2]=atom_z[jat];
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];
	
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
	    klps0  = csp_leading_ps_pair[klcs];
	    nklps  = csp_leading_ps_pair[klcs+1]-klps0;
	    kat    = shel_atm[kcs];
	    lat    = shel_atm[lcs];
	    kao0   = shel_ini[kcs];
	    lao0   = shel_ini[lcs];
	    C[0]=atom_x[kat]; C[1]=atom_y[kat]; C[2]=atom_z[kat];
	    D[0]=atom_x[lat]; D[1]=atom_y[lat]; D[2]=atom_z[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dspp__(
		    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
		    &psp_xiza[ijps0], BA,
		    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
		    &psp_xiza[klps0], DC,   AC,      DSPP );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		for (k=0, kao=kao0; k<3; k++, kao++ ) {
		    for ( l=0, lao=lao0; l<3; l++, lao++, ix++ ) {
			if ( lao > kao ) continue;
			coe = ( kao == lao ? HALF : ONE );
			if ( fabs(DSPP[ix]) > EPS_ERI ) {
			    ebuf_val[nzeri]     = coe*DSPP[ix];
			    /*
		// for check sum
		lsum += ebuf_val[nzeri];
		*/

			    ebuf_ind4[nzeri4+0] = (short int)iao;
			    ebuf_ind4[nzeri4+1] = (short int)jao;
			    ebuf_ind4[nzeri4+2] = (short int)kao;
			    ebuf_ind4[nzeri4+3] = (short int)lao;
			    nzeri++;
			    nzeri4+=4;
			}	// if ( fabs(DSPP[i]) > EPS_ERI )
		    }
		}
	    }
	    if ( nzeri >= max_nzeri ) {
		*last_ijcs = ijcs;
		*last_klcs = klcs;
		*ebuf_non_zero_eri = nzeri;
		return OFMO_EBUF_FULL;
	    }
	}	// for ( klcs );
    }		// for ( ijcs );
    *ebuf_non_zero_eri = nzeri;
    /*
    // for check sum
#pragma omp atomic
    check_sum[8] += lsum;
    */

    return OFMO_EBUF_NOFULL;
}

/** (ds,ds)タイプの２電子積分関数
 * @ingroup integ-twoint-buffer
 * */
int ofmo_twoint_buffer_dsds__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs
	) {
    int Lab, Lcd, i, k, ix;
    int IJ, KL, ipat;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    double DSDS[6*6];
    long nzeri, max_nzeri, nzeri4;
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    long ebuf_max_nzeri = *pebuf_max_nzeri;
    /*
    // for check sum
    double lsum = ZERO;
    */

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs0 = leading_cs_pair[Lab];
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];

    nzeri     = *ebuf_non_zero_eri;
    max_nzeri = ebuf_max_nzeri - 6*6;
    nzeri4    = nzeri*4;
    if ( nzeri >= max_nzeri ) { 
	*last_ijcs = ijcs0+workerid;
	*last_klcs = klcs0 - 1;
	*ebuf_non_zero_eri = nzeri;
	return OFMO_EBUF_FULL;
    }
    /*
    // debug
    printf("(ds,ds) :  nwks=%d, wkid=%d, max nzeri=%lld\n",
	    nworkers, workerid, max_nzeri );
    fflush(stdout);
    */

    for ( ijcs=ijcs0+workerid; ijcs<ijcs1; ijcs+=nworkers ) {
	val_ab = csp_schwarz[ijcs];
	ics    = csp_ics[ijcs];
	jcs    = csp_jcs[ijcs];
	ijps0  = csp_leading_ps_pair[ijcs];
	nijps  = csp_leading_ps_pair[ijcs+1]-ijps0;
	iat    = shel_atm[ics];
	jat    = shel_atm[jcs];
	iao0   = shel_ini[ics];
	jao    = shel_ini[jcs];
	A[0]=atom_x[iat]; A[1]=atom_y[iat]; A[2]=atom_z[iat];
	B[0]=atom_x[jat]; B[1]=atom_y[jat]; B[2]=atom_z[jat];
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];
	
	for ( klcs=klcs0; klcs<=ijcs; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
	    klps0  = csp_leading_ps_pair[klcs];
	    nklps  = csp_leading_ps_pair[klcs+1]-klps0;
	    kat    = shel_atm[kcs];
	    lat    = shel_atm[lcs];
	    kao0   = shel_ini[kcs];
	    lao    = shel_ini[lcs];
	    C[0]=atom_x[kat]; C[1]=atom_y[kat]; C[2]=atom_z[kat];
	    D[0]=atom_x[lat]; D[1]=atom_y[lat]; D[2]=atom_z[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dsds__(
		    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
		    &psp_xiza[ijps0], BA,
		    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
		    &psp_xiza[klps0], DC,   AC,      DSDS );
	    /* ダメな加算方法
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		for ( k=0, kao=kao0; k<6; k++, kao++, ix++ ) {
		    if ( kao > iao ) continue;
		    if ( fabs(DSDS[ix]) > EPS_ERI ) {
			coe = ONE;
			if ( iao==kao && jao==lao ) coe = HALF;
			ebuf_val[nzeri]     = coe*DSDS[ix];
		// for check sum
		lsum += ebuf_val[nzeri];

			ebuf_ind4[nzeri4+0] = (short int)iao;
			ebuf_ind4[nzeri4+1] = (short int)jao;
			ebuf_ind4[nzeri4+2] = (short int)kao;
			ebuf_ind4[nzeri4+3] = (short int)lao;
			nzeri++;
			nzeri4+=4;
		    }
		}
	    }
	    */
	    ipat = ( (ics==kcs && jcs>lcs) ? true : false);
#ifdef SORT_CSP
            int ijgekl = (ics>kcs);
            if (ics==kcs) ijgekl = (jcs>=lcs);
            if (!ijgekl) ipat = ( (ics==kcs && jcs<lcs) ? true : false);
#endif
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		IJ = ( (iao*iao+iao)>>1) + jao;
		for ( k=0, kao=kao0; k<6; k++, kao++, ix++ ) {
		    KL = ( (kao*kao+kao)>>1) + lao;
		    if ( fabs(DSDS[ix]) <= EPS_ERI )  continue;
#ifndef SORT_CSP
		    if ( IJ >= KL ) {
#else
                    if ( (ijgekl&&IJ>=KL) || (!ijgekl&&KL>=IJ) ) {
#endif
			coe = ( IJ==KL? HALF : ONE);
			ebuf_val[nzeri]     = coe*DSDS[ix];
			/*
		// for check sum
		lsum += ebuf_val[nzeri];
		*/

			ebuf_ind4[nzeri4+0] = (short int)iao;
			ebuf_ind4[nzeri4+1] = (short int)jao;
			ebuf_ind4[nzeri4+2] = (short int)kao;
			ebuf_ind4[nzeri4+3] = (short int)lao;
			nzeri++;
			nzeri4+=4;
		    } else if ( ipat ) {
			ebuf_val[nzeri]     = DSDS[ix];
			/*
		// for check sum
		lsum += ebuf_val[nzeri];
		*/

			ebuf_ind4[nzeri4+0] = (short int)kao;
			ebuf_ind4[nzeri4+1] = (short int)lao;
			ebuf_ind4[nzeri4+2] = (short int)iao;
			ebuf_ind4[nzeri4+3] = (short int)jao;
			nzeri++;
			nzeri4+=4;
		    }
		}
	    }
	    if ( nzeri >= max_nzeri ) {
		*last_ijcs = ijcs;
		*last_klcs = klcs;
		*ebuf_non_zero_eri = nzeri;
		return OFMO_EBUF_FULL;
	    }
	}	// for ( klcs );
    }		// for ( ijcs );
    *ebuf_non_zero_eri = nzeri;
    /*
    // for check sum
#pragma omp atomic
    check_sum[9] += lsum;
    */

    return OFMO_EBUF_NOFULL;
}

/** (dp,ss)タイプの２電子積分関数
 * @ingroup integ-twoint-buffer
 * */
int ofmo_twoint_buffer_dpss__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs
	) {
    int Lab, Lcd, i, j, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    double DPSS[6*3];
    long nzeri, max_nzeri, nzeri4;
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    long ebuf_max_nzeri = *pebuf_max_nzeri;
    /*
    // for check sum
    double lsum = ZERO;
    */

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs0 = leading_cs_pair[Lab];
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];

    nzeri     = *ebuf_non_zero_eri;
    max_nzeri = ebuf_max_nzeri - 6*3;
    nzeri4    = nzeri*4;
    if ( nzeri >= max_nzeri ) { 
	*last_ijcs = ijcs0+workerid;
	*last_klcs = klcs0 - 1;
	*ebuf_non_zero_eri = nzeri;
	return OFMO_EBUF_FULL;
    }
    /*
    // debug
    printf("(dp,ss) :  nwks=%d, wkid=%d, max nzeri=%lld\n",
	    nworkers, workerid, max_nzeri );
    fflush(stdout);
    */

    for ( ijcs=ijcs0+workerid; ijcs<ijcs1; ijcs+=nworkers ) {
	val_ab = csp_schwarz[ijcs];
	ics    = csp_ics[ijcs];
	jcs    = csp_jcs[ijcs];
	ijps0  = csp_leading_ps_pair[ijcs];
	nijps  = csp_leading_ps_pair[ijcs+1]-ijps0;
	iat    = shel_atm[ics];
	jat    = shel_atm[jcs];
	iao0   = shel_ini[ics];
	jao0   = shel_ini[jcs];
	A[0]=atom_x[iat]; A[1]=atom_y[iat]; A[2]=atom_z[iat];
	B[0]=atom_x[jat]; B[1]=atom_y[jat]; B[2]=atom_z[jat];
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];
	
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
	    klps0  = csp_leading_ps_pair[klcs];
	    nklps  = csp_leading_ps_pair[klcs+1]-klps0;
	    kat    = shel_atm[kcs];
	    lat    = shel_atm[lcs];
	    kao    = shel_ini[kcs];
	    lao    = shel_ini[lcs];
	    C[0]=atom_x[kat]; C[1]=atom_y[kat]; C[2]=atom_z[kat];
	    D[0]=atom_x[lat]; D[1]=atom_y[lat]; D[2]=atom_z[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dpss__(
		    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
		    &psp_xiza[ijps0], BA,
		    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
		    &psp_xiza[klps0], DC,   AC,      DPSS );
	    coe = ( kao == lao ? HALF : ONE );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		for ( j=0, jao=jao0; j<3; j++, jao++, ix++ ) {
		    if ( fabs(DPSS[ix]) > EPS_ERI ) {
			ebuf_val[nzeri]     = coe*DPSS[ix];
			/*
		// for check sum
		lsum += ebuf_val[nzeri];
		*/

			ebuf_ind4[nzeri4+0] = (short int)iao;
			ebuf_ind4[nzeri4+1] = (short int)jao;
			ebuf_ind4[nzeri4+2] = (short int)kao;
			ebuf_ind4[nzeri4+3] = (short int)lao;
			nzeri++;
			nzeri4+=4;
		    }	// if ( fabs(PSSS[i]) > EPS_ERI )
		}
	    }
	    if ( nzeri >= max_nzeri ) {
		*last_ijcs = ijcs;
		*last_klcs = klcs;
		*ebuf_non_zero_eri = nzeri;
		return OFMO_EBUF_FULL;
	    }
	}	// for ( klcs );
    }		// for ( ijcs );
    *ebuf_non_zero_eri = nzeri;
    /*
    // for check sum
#pragma omp atomic
    check_sum[10] += lsum;
    */

    return OFMO_EBUF_NOFULL;
}

/** (dp,ps)タイプの２電子積分関数
 * @ingroup integ-twoint-buffer
 * */
int ofmo_twoint_buffer_dpps__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs
	) {
    int Lab, Lcd, i, j, k, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd;
    double DPPS[6*3*3];
    long nzeri, max_nzeri, nzeri4;
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    long ebuf_max_nzeri = *pebuf_max_nzeri;
    /*
    // for check sum
    double lsum = ZERO;
    */

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs0 = leading_cs_pair[Lab];
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];

    nzeri     = *ebuf_non_zero_eri;
    max_nzeri = ebuf_max_nzeri - 6*3*3;
    nzeri4    = nzeri*4;
    if ( nzeri >= max_nzeri ) { 
	*last_ijcs = ijcs0+workerid;
	*last_klcs = klcs0 - 1;
	*ebuf_non_zero_eri = nzeri;
	return OFMO_EBUF_FULL;
    }
    /*
    // debug
    printf("(dp,ps) :  nwks=%d, wkid=%d, max nzeri=%lld\n",
	    nworkers, workerid, max_nzeri );
    fflush(stdout);
    */

    for ( ijcs=ijcs0+workerid; ijcs<ijcs1; ijcs+=nworkers ) {
	val_ab = csp_schwarz[ijcs];
	ics    = csp_ics[ijcs];
	jcs    = csp_jcs[ijcs];
	ijps0  = csp_leading_ps_pair[ijcs];
	nijps  = csp_leading_ps_pair[ijcs+1]-ijps0;
	iat    = shel_atm[ics];
	jat    = shel_atm[jcs];
	iao0   = shel_ini[ics];
	jao0   = shel_ini[jcs];
	A[0]=atom_x[iat]; A[1]=atom_y[iat]; A[2]=atom_z[iat];
	B[0]=atom_x[jat]; B[1]=atom_y[jat]; B[2]=atom_z[jat];
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];
	
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
	    klps0  = csp_leading_ps_pair[klcs];
	    nklps  = csp_leading_ps_pair[klcs+1]-klps0;
	    kat    = shel_atm[kcs];
	    lat    = shel_atm[lcs];
	    kao0   = shel_ini[kcs];
	    lao    = shel_ini[lcs];
	    C[0]=atom_x[kat]; C[1]=atom_y[kat]; C[2]=atom_z[kat];
	    D[0]=atom_x[lat]; D[1]=atom_y[lat]; D[2]=atom_z[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dpps__(
		    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
		    &psp_xiza[ijps0], BA,
		    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
		    &psp_xiza[klps0], DC,   AC,      DPPS );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		for ( j=0, jao=jao0; j<3; j++, jao++ ) {
		    for ( k=0, kao=kao0; k<3; k++, kao++, ix++ ) {
			if ( fabs(DPPS[ix]) > EPS_ERI ) {
			    ebuf_val[nzeri]     = DPPS[ix];
			    /*
		// for check sum
		lsum += ebuf_val[nzeri];
		*/

			    ebuf_ind4[nzeri4+0] = (short int)iao;
			    ebuf_ind4[nzeri4+1] = (short int)jao;
			    ebuf_ind4[nzeri4+2] = (short int)kao;
			    ebuf_ind4[nzeri4+3] = (short int)lao;
			    nzeri++;
			    nzeri4+=4;
			}	// if ( fabs(PSSS[i]) > EPS_ERI )
		    }
		}
	    }
	    if ( nzeri >= max_nzeri ) {
		*last_ijcs = ijcs;
		*last_klcs = klcs;
		*ebuf_non_zero_eri = nzeri;
		return OFMO_EBUF_FULL;
	    }
	}	// for ( klcs );
    }		// for ( ijcs );
    *ebuf_non_zero_eri = nzeri;
    /*
    // for check sum
#pragma omp atomic
    check_sum[11] += lsum;
    */

    return OFMO_EBUF_NOFULL;
}

/** (dp,pp)タイプの２電子積分関数
 * @ingroup integ-twoint-buffer
 * */
int ofmo_twoint_buffer_dppp__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs
	) {
    int Lab, Lcd, i, j, k, l, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    double DPPP[6*3*3*3];
    long nzeri, max_nzeri, nzeri4;
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    long ebuf_max_nzeri = *pebuf_max_nzeri;
    /*
    // for check sum
    double lsum = ZERO;
    */

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs0 = leading_cs_pair[Lab];
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];

    nzeri     = *ebuf_non_zero_eri;
    max_nzeri = ebuf_max_nzeri - 6*3*3*3;
    nzeri4    = nzeri*4;
    if ( nzeri >= max_nzeri ) { 
	*last_ijcs = ijcs0+workerid;
	*last_klcs = klcs0 - 1;
	*ebuf_non_zero_eri = nzeri;
	return OFMO_EBUF_FULL;
    }
    /*
    // debug
    printf("(dp,pp) :  nwks=%d, wkid=%d, max nzeri=%lld\n",
	    nworkers, workerid, max_nzeri );
    fflush(stdout);
    */

    for ( ijcs=ijcs0+workerid; ijcs<ijcs1; ijcs+=nworkers ) {
	val_ab = csp_schwarz[ijcs];
	ics    = csp_ics[ijcs];
	jcs    = csp_jcs[ijcs];
	ijps0  = csp_leading_ps_pair[ijcs];
	nijps  = csp_leading_ps_pair[ijcs+1]-ijps0;
	iat    = shel_atm[ics];
	jat    = shel_atm[jcs];
	iao0   = shel_ini[ics];
	jao0   = shel_ini[jcs];
	A[0]=atom_x[iat]; A[1]=atom_y[iat]; A[2]=atom_z[iat];
	B[0]=atom_x[jat]; B[1]=atom_y[jat]; B[2]=atom_z[jat];
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];
	
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
	    klps0  = csp_leading_ps_pair[klcs];
	    nklps  = csp_leading_ps_pair[klcs+1]-klps0;
	    kat    = shel_atm[kcs];
	    lat    = shel_atm[lcs];
	    kao0   = shel_ini[kcs];
	    lao0   = shel_ini[lcs];
	    C[0]=atom_x[kat]; C[1]=atom_y[kat]; C[2]=atom_z[kat];
	    D[0]=atom_x[lat]; D[1]=atom_y[lat]; D[2]=atom_z[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dppp__(
		    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
		    &psp_xiza[ijps0], BA,
		    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
		    &psp_xiza[klps0], DC,   AC,      DPPP );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		for ( j=0, jao=jao0; j<3; j++, jao++ ) {
		    for (k=0, kao=kao0; k<3; k++, kao++ ) {
			for ( l=0, lao=lao0; l<3; l++, lao++, ix++ ) {
			    if ( lao > kao ) continue;
			    coe = ( kao == lao ? HALF : ONE );
			    if ( fabs(DPPP[ix]) > EPS_ERI ) {
				ebuf_val[nzeri]     = coe*DPPP[ix];
				/*
		// for check sum
		lsum += ebuf_val[nzeri];
		*/

				ebuf_ind4[nzeri4+0] = (short int)iao;
				ebuf_ind4[nzeri4+1] = (short int)jao;
				ebuf_ind4[nzeri4+2] = (short int)kao;
				ebuf_ind4[nzeri4+3] = (short int)lao;
				nzeri++;
				nzeri4+=4;
			    }	// if ( fabs(DSPP[i]) > EPS_ERI )
			}
		    }
		}
	    }
	    if ( nzeri >= max_nzeri ) {
		*last_ijcs = ijcs;
		*last_klcs = klcs;
		*ebuf_non_zero_eri = nzeri;
		return OFMO_EBUF_FULL;
	    }
	}	// for ( klcs );
    }		// for ( ijcs );
    *ebuf_non_zero_eri = nzeri;
    /*
    // for check sum
#pragma omp atomic
    check_sum[12] += lsum;
    */

    return OFMO_EBUF_NOFULL;
}

/** (dp,ds)タイプの２電子積分関数
 * @ingroup integ-twoint-buffer
 * */
int ofmo_twoint_buffer_dpds__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs
	) {
    int Lab, Lcd, i, j, k, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd;
    double DPDS[6*3*6];
    long nzeri, max_nzeri, nzeri4;
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    long ebuf_max_nzeri = *pebuf_max_nzeri;
    /*
    // for check sum
    double lsum = ZERO;
    */

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs0 = leading_cs_pair[Lab];
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];

    nzeri     = *ebuf_non_zero_eri;
    max_nzeri = ebuf_max_nzeri - 6*3*6;
    nzeri4    = nzeri*4;
    if ( nzeri >= max_nzeri ) { 
	*last_ijcs = ijcs0+workerid;
	*last_klcs = klcs0 - 1;
	*ebuf_non_zero_eri = nzeri;
	return OFMO_EBUF_FULL;
    }
    /*
    // debug
    printf("(dp,ds) :  nwks=%d, wkid=%d, max nzeri=%lld\n",
	    nworkers, workerid, max_nzeri );
    fflush(stdout);
    */

    for ( ijcs=ijcs0+workerid; ijcs<ijcs1; ijcs+=nworkers ) {
	val_ab = csp_schwarz[ijcs];
	ics    = csp_ics[ijcs];
	jcs    = csp_jcs[ijcs];
	ijps0  = csp_leading_ps_pair[ijcs];
	nijps  = csp_leading_ps_pair[ijcs+1]-ijps0;
	iat    = shel_atm[ics];
	jat    = shel_atm[jcs];
	iao0   = shel_ini[ics];
	jao0   = shel_ini[jcs];
	A[0]=atom_x[iat]; A[1]=atom_y[iat]; A[2]=atom_z[iat];
	B[0]=atom_x[jat]; B[1]=atom_y[jat]; B[2]=atom_z[jat];
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];
	
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
	    klps0  = csp_leading_ps_pair[klcs];
	    nklps  = csp_leading_ps_pair[klcs+1]-klps0;
	    kat    = shel_atm[kcs];
	    lat    = shel_atm[lcs];
	    kao0   = shel_ini[kcs];
	    lao    = shel_ini[lcs];
	    C[0]=atom_x[kat]; C[1]=atom_y[kat]; C[2]=atom_z[kat];
	    D[0]=atom_x[lat]; D[1]=atom_y[lat]; D[2]=atom_z[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dpds__(
		    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
		    &psp_xiza[ijps0], BA,
		    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
		    &psp_xiza[klps0], DC,   AC,      DPDS );

	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		for ( j=0, jao=jao0; j<3; j++, jao++ ) {
		    for ( k=0, kao=kao0; k<6; k++, kao++, ix++ ) {
			if ( fabs(DPDS[ix]) > EPS_ERI ) {
			    ebuf_val[nzeri]     = DPDS[ix];
			    /*
		// for check sum
		lsum += ebuf_val[nzeri];
		*/

			    ebuf_ind4[nzeri4+0] = (short int)iao;
			    ebuf_ind4[nzeri4+1] = (short int)jao;
			    ebuf_ind4[nzeri4+2] = (short int)kao;
			    ebuf_ind4[nzeri4+3] = (short int)lao;
			    nzeri++;
			    nzeri4+=4;
			}	// if ( fabs(DPDS[i]) > EPS_ERI )
		    }
		}
	    }
	    if ( nzeri >= max_nzeri ) {
		*last_ijcs = ijcs;
		*last_klcs = klcs;
		*ebuf_non_zero_eri = nzeri;
		return OFMO_EBUF_FULL;
	    }
	}	// for ( klcs );
    }		// for ( ijcs );
    *ebuf_non_zero_eri = nzeri;
    /*
    // for check sum
#pragma omp atomic
    check_sum[13] += lsum;
    */

    return OFMO_EBUF_NOFULL;
}

/** (dp,dp)タイプの２電子積分関数
 * @ingroup integ-twoint-buffer
 * */
int ofmo_twoint_buffer_dpdp__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs
	) {
    int Lab, Lcd, i, j, k, l, ix, ipat;
    int IJ, KL;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    double DPDP[6*3*6*3];
    long nzeri, max_nzeri, nzeri4;
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    long ebuf_max_nzeri = *pebuf_max_nzeri;
    /*
    // for check sum
    double lsum = ZERO;
    */
    //int nnn = 0;

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs0 = leading_cs_pair[Lab];
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];

    nzeri     = *ebuf_non_zero_eri;
    max_nzeri = ebuf_max_nzeri - 6*3*6*3;
    nzeri4    = nzeri*4;
    if ( nzeri >= max_nzeri ) { 
	*last_ijcs = ijcs0+workerid;
	*last_klcs = klcs0 - 1;
	*ebuf_non_zero_eri = nzeri;
	return OFMO_EBUF_FULL;
    }
    /*
    // debug
    printf("(dp,dp) :  nwks=%d, wkid=%d, max nzeri=%lld\n",
	    nworkers, workerid, max_nzeri );
    fflush(stdout);
    */

    for ( ijcs=ijcs0+workerid; ijcs<ijcs1; ijcs+=nworkers ) {
	val_ab = csp_schwarz[ijcs];
	ics    = csp_ics[ijcs];
	jcs    = csp_jcs[ijcs];
	ijps0  = csp_leading_ps_pair[ijcs];
	nijps  = csp_leading_ps_pair[ijcs+1]-ijps0;
	iat    = shel_atm[ics];
	jat    = shel_atm[jcs];
	iao0   = shel_ini[ics];
	jao0   = shel_ini[jcs];
	A[0]=atom_x[iat]; A[1]=atom_y[iat]; A[2]=atom_z[iat];
	B[0]=atom_x[jat]; B[1]=atom_y[jat]; B[2]=atom_z[jat];
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];
	
	for ( klcs=klcs0; klcs<=ijcs; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
	    klps0  = csp_leading_ps_pair[klcs];
	    nklps  = csp_leading_ps_pair[klcs+1]-klps0;
	    kat    = shel_atm[kcs];
	    lat    = shel_atm[lcs];
	    kao0   = shel_ini[kcs];
	    lao0   = shel_ini[lcs];
	    C[0]=atom_x[kat]; C[1]=atom_y[kat]; C[2]=atom_z[kat];
	    D[0]=atom_x[lat]; D[1]=atom_y[lat]; D[2]=atom_z[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dpdp__(
		    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
		    &psp_xiza[ijps0], BA,
		    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
		    &psp_xiza[klps0], DC,   AC,      DPDP );
	    ipat = ( (ics==kcs && jcs>lcs) ? true : false );
#ifdef SORT_CSP
            int ijgekl = (ics>kcs);
            if (ics==kcs) ijgekl = (jcs>=lcs);
            if (!ijgekl) ipat = ( (ics==kcs && jcs<lcs) ? true : false);
#endif
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		for ( j=0, jao=jao0; j<3; j++, jao++ ) {
		    IJ = ((iao*iao+iao)>>1) + jao;
		    for ( k=0, kao=kao0; k<6; k++, kao++ ) {
			for ( l=0, lao=lao0; l<3; l++, lao++, ix++ ) {
			    KL = ((kao*kao+kao)>>1) + lao;
			    if ( fabs(DPDP[ix]) > EPS_ERI ) {
#ifndef SORT_CSP
				if ( IJ >= KL ) {
#else
                    if ( (ijgekl&&IJ>=KL) || (!ijgekl&&KL>=IJ) ) {
#endif
				    coe = ( IJ==KL ? HALF : ONE );
				    ebuf_val[nzeri]     = coe*DPDP[ix];
				    /*
		// for check sum
		lsum += ebuf_val[nzeri];
		*/

				    ebuf_ind4[nzeri4+0] = (short int)iao;
				    ebuf_ind4[nzeri4+1] = (short int)jao;
				    ebuf_ind4[nzeri4+2] = (short int)kao;
				    ebuf_ind4[nzeri4+3] = (short int)lao;
				    nzeri++;
				    nzeri4+=4;
				    /*
				    //debug
				    nnn++;
				    if ( nnn < 1000 ) {
					printf("%3d %3d %3d %3d",
						iao, jao, kao, lao);
					printf("  %23.15f\n", 
						ebuf_val[nzeri-1] );
				    }
				    */
				} else if ( ipat ) {
				    ebuf_val[nzeri]     = DPDP[ix];
				    /*
		// for check sum
		lsum += ebuf_val[nzeri];
		*/

				    ebuf_ind4[nzeri4+0] = (short int)kao;
				    ebuf_ind4[nzeri4+1] = (short int)lao;
				    ebuf_ind4[nzeri4+2] = (short int)iao;
				    ebuf_ind4[nzeri4+3] = (short int)jao;
				    nzeri++;
				    nzeri4+=4;
				    /*
				    // debug
				    nnn++;
				    if ( nnn < 1000 ) {
					printf("%3d %3d %3d %3d",
						iao, jao, kao, lao);
					printf("  %23.15f\n", 
						ebuf_val[nzeri-1] );
				    }
				    */
				}
			    }
			}
		    }
		}
	    }
	    if ( nzeri >= max_nzeri ) {
		*last_ijcs = ijcs;
		*last_klcs = klcs;
		*ebuf_non_zero_eri = nzeri;
		return OFMO_EBUF_FULL;
	    }
	}	// for ( klcs );
    }		// for ( ijcs );
    *ebuf_non_zero_eri = nzeri;
    /*
    // for check sum
#pragma omp atomic
    check_sum[14] += lsum;
    */

    //printf("dpdp  = %d\n", nnn );

    return OFMO_EBUF_NOFULL;
}

/** (dd,ss)タイプの２電子積分関数
 * @ingroup integ-twoint-buffer
 * */
int ofmo_twoint_buffer_ddss__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs
	) {
    int Lab, Lcd, i, j, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe, coe0;
    double DDSS[6*6];
    long nzeri, max_nzeri, nzeri4;
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    long ebuf_max_nzeri = *pebuf_max_nzeri;
    /*
    // for check sum
    double lsum = ZERO;
    */

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs0 = leading_cs_pair[Lab];
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];

    nzeri     = *ebuf_non_zero_eri;
    max_nzeri = ebuf_max_nzeri - 6*6;
    nzeri4    = nzeri*4;
    if ( nzeri >= max_nzeri ) { 
	*last_ijcs = ijcs0+workerid;
	*last_klcs = klcs0 - 1;
	*ebuf_non_zero_eri = nzeri;
	return OFMO_EBUF_FULL;
    }
    /*
    // debug
    printf("(dd,ss) :  nwks=%d, wkid=%d, max nzeri=%lld, nzeri=%lld\n",
	    nworkers, workerid, max_nzeri, nzeri );
    fflush(stdout);
    */

    for ( ijcs=ijcs0+workerid; ijcs<ijcs1; ijcs+=nworkers ) {
	val_ab = csp_schwarz[ijcs];
	ics    = csp_ics[ijcs];
	jcs    = csp_jcs[ijcs];
	ijps0  = csp_leading_ps_pair[ijcs];
	nijps  = csp_leading_ps_pair[ijcs+1]-ijps0;
	iat    = shel_atm[ics];
	jat    = shel_atm[jcs];
	iao0   = shel_ini[ics];
	jao0   = shel_ini[jcs];
	A[0]=atom_x[iat]; A[1]=atom_y[iat]; A[2]=atom_z[iat];
	B[0]=atom_x[jat]; B[1]=atom_y[jat]; B[2]=atom_z[jat];
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];
	
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
	    klps0  = csp_leading_ps_pair[klcs];
	    nklps  = csp_leading_ps_pair[klcs+1]-klps0;
	    kat    = shel_atm[kcs];
	    lat    = shel_atm[lcs];
	    kao    = shel_ini[kcs];
	    lao    = shel_ini[lcs];
	    C[0]=atom_x[kat]; C[1]=atom_y[kat]; C[2]=atom_z[kat];
	    D[0]=atom_x[lat]; D[1]=atom_y[lat]; D[2]=atom_z[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_ddss__(
		    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
		    &psp_xiza[ijps0], BA,
		    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
		    &psp_xiza[klps0], DC,   AC,      DDSS );
	    coe0 = ( kao==lao ? HALF : ONE );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		for ( j=0, jao=jao0; j<6; j++, jao++, ix++ ) {
		    if ( jao > iao ) continue;
		    if ( fabs(DDSS[ix]) > EPS_ERI ) {
			coe = coe0;
			if ( iao == jao ) coe *= HALF;
			ebuf_val[nzeri]     = coe*DDSS[ix];
			/*
		// for check sum
		lsum += ebuf_val[nzeri];
		*/

			ebuf_ind4[nzeri4+0] = (short int)iao;
			ebuf_ind4[nzeri4+1] = (short int)jao;
			ebuf_ind4[nzeri4+2] = (short int)kao;
			ebuf_ind4[nzeri4+3] = (short int)lao;
			nzeri++;
			nzeri4+=4;
		    }	// if ( fabs(PPSS[ix]) > EPS_ERI )
		}
	    }
	    if ( nzeri >= max_nzeri ) {
		*last_ijcs = ijcs;
		*last_klcs = klcs;
		*ebuf_non_zero_eri = nzeri;
		return OFMO_EBUF_FULL;
	    }
	}	// for ( klcs );
    }		// for ( ijcs );
    *ebuf_non_zero_eri = nzeri;
    /*
    // for check sum
#pragma omp atomic
    check_sum[15] += lsum;
    */

    return OFMO_EBUF_NOFULL;
}

/** (dd,ps)タイプの２電子積分関数
 * */
int ofmo_twoint_buffer_ddps__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs
	) {
    int Lab, Lcd, i, j, k, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    double DDPS[6*6*3];
    long nzeri, max_nzeri, nzeri4;
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    long ebuf_max_nzeri = *pebuf_max_nzeri;
    /*
    // for check sum
    double lsum = ZERO;
    */

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs0 = leading_cs_pair[Lab];
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];

    nzeri     = *ebuf_non_zero_eri;
    max_nzeri = ebuf_max_nzeri - 6*6*3;
    nzeri4    = nzeri*4;
    if ( nzeri >= max_nzeri ) { 
	*last_ijcs = ijcs0+workerid;
	*last_klcs = klcs0 - 1;
	*ebuf_non_zero_eri = nzeri;
	return OFMO_EBUF_FULL;
    }
    /*
    // debug
    printf("(dd,ps) :  nwks=%d, wkid=%d, max nzeri=%lld, nzeri=%lld\n",
	    nworkers, workerid, max_nzeri, nzeri );
    fflush(stdout);
    */

    for ( ijcs=ijcs0+workerid; ijcs<ijcs1; ijcs+=nworkers ) {
	val_ab = csp_schwarz[ijcs];
	ics    = csp_ics[ijcs];
	jcs    = csp_jcs[ijcs];
	ijps0  = csp_leading_ps_pair[ijcs];
	nijps  = csp_leading_ps_pair[ijcs+1]-ijps0;
	iat    = shel_atm[ics];
	jat    = shel_atm[jcs];
	iao0   = shel_ini[ics];
	jao0   = shel_ini[jcs];
	A[0]=atom_x[iat]; A[1]=atom_y[iat]; A[2]=atom_z[iat];
	B[0]=atom_x[jat]; B[1]=atom_y[jat]; B[2]=atom_z[jat];
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];
	
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
	    klps0  = csp_leading_ps_pair[klcs];
	    nklps  = csp_leading_ps_pair[klcs+1]-klps0;
	    kat    = shel_atm[kcs];
	    lat    = shel_atm[lcs];
	    kao0   = shel_ini[kcs];
	    lao    = shel_ini[lcs];
	    C[0]=atom_x[kat]; C[1]=atom_y[kat]; C[2]=atom_z[kat];
	    D[0]=atom_x[lat]; D[1]=atom_y[lat]; D[2]=atom_z[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_ddps__(
		    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
		    &psp_xiza[ijps0], BA,
		    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
		    &psp_xiza[klps0], DC,   AC,      DDPS );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		for ( j=0, jao=jao0; j<6; j++, jao++ ) {
		    if ( jao>iao ) { ix += 3; continue; }
		    coe = ( iao==jao ? HALF : ONE );
		    for ( k=0, kao=kao0; k<3; k++, kao++, ix++ ) {
			if ( fabs(DDPS[ix]) > EPS_ERI ) {
			    ebuf_val[nzeri]     = coe*DDPS[ix];
			    /*
		// for check sum
		lsum += ebuf_val[nzeri];
		*/

			    ebuf_ind4[nzeri4+0] = (short int)iao;
			    ebuf_ind4[nzeri4+1] = (short int)jao;
			    ebuf_ind4[nzeri4+2] = (short int)kao;
			    ebuf_ind4[nzeri4+3] = (short int)lao;
			    nzeri++;
			    nzeri4+=4;
			}	// if ( fabs(PPSS[ix]) > EPS_ERI )
		    }
		}
	    }
	    if ( nzeri >= max_nzeri ) {
		*last_ijcs = ijcs;
		*last_klcs = klcs;
		*ebuf_non_zero_eri = nzeri;
		return OFMO_EBUF_FULL;
	    }
	}	// for ( klcs );
    }		// for ( ijcs );
    *ebuf_non_zero_eri = nzeri;
    /*
    // for check sum
#pragma omp atomic
    check_sum[16] += lsum;
    */

    return OFMO_EBUF_NOFULL;
}

/** (dd,pp)タイプの２電子積分関数
 * */
int ofmo_twoint_buffer_ddpp__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs
	) {
    int Lab, Lcd, i, j, k, l, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe, coe0;
    double DDPP[6*6*3*3];
    long nzeri, max_nzeri, nzeri4;
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    long ebuf_max_nzeri = *pebuf_max_nzeri;
    /*
    // for check sum
    double lsum = ZERO;
    */

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs0 = leading_cs_pair[Lab];
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];

    nzeri     = *ebuf_non_zero_eri;
    max_nzeri = ebuf_max_nzeri - 6*6*3*3;
    nzeri4    = nzeri*4;
    if ( nzeri >= max_nzeri ) { 
	*last_ijcs = ijcs0+workerid;
	*last_klcs = klcs0 - 1;
	*ebuf_non_zero_eri = nzeri;
	return OFMO_EBUF_FULL;
    }
    /*
    // debug
    printf("(dd,pp) :  nwks=%d, wkid=%d, max nzeri=%lld, nzeri=%lld\n",
	    nworkers, workerid, max_nzeri, nzeri );
    fflush(stdout);
    */

    for ( ijcs=ijcs0+workerid; ijcs<ijcs1; ijcs+=nworkers ) {
	val_ab = csp_schwarz[ijcs];
	ics    = csp_ics[ijcs];
	jcs    = csp_jcs[ijcs];
	ijps0  = csp_leading_ps_pair[ijcs];
	nijps  = csp_leading_ps_pair[ijcs+1]-ijps0;
	iat    = shel_atm[ics];
	jat    = shel_atm[jcs];
	iao0   = shel_ini[ics];
	jao0   = shel_ini[jcs];
	A[0]=atom_x[iat]; A[1]=atom_y[iat]; A[2]=atom_z[iat];
	B[0]=atom_x[jat]; B[1]=atom_y[jat]; B[2]=atom_z[jat];
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];
	
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
	    klps0  = csp_leading_ps_pair[klcs];
	    nklps  = csp_leading_ps_pair[klcs+1]-klps0;
	    kat    = shel_atm[kcs];
	    lat    = shel_atm[lcs];
	    kao0   = shel_ini[kcs];
	    lao0   = shel_ini[lcs];
	    C[0]=atom_x[kat]; C[1]=atom_y[kat]; C[2]=atom_z[kat];
	    D[0]=atom_x[lat]; D[1]=atom_y[lat]; D[2]=atom_z[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_ddpp__(
		    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
		    &psp_xiza[ijps0], BA,
		    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
		    &psp_xiza[klps0], DC,   AC,      DDPP );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		for ( j=0, jao=jao0; j<6; j++, jao++ ) {
		    if ( jao>iao ) { ix += 3*3; continue; }
		    coe0 = ( iao==jao ? HALF : ONE );
		    for ( k=0, kao=kao0; k<3; k++, kao++ ) {
			for ( l=0, lao=lao0; l<3; l++, lao++, ix++ ) {
			    if ( lao > kao ) continue;
			    coe = coe0 * ( kao==lao ? HALF : ONE );
			    if ( fabs(DDPP[ix]) > EPS_ERI ) {
				ebuf_val[nzeri]     = coe*DDPP[ix];
				/*
		// for check sum
		lsum += ebuf_val[nzeri];
		*/

				ebuf_ind4[nzeri4+0] = (short int)iao;
				ebuf_ind4[nzeri4+1] = (short int)jao;
				ebuf_ind4[nzeri4+2] = (short int)kao;
				ebuf_ind4[nzeri4+3] = (short int)lao;
				nzeri++;
				nzeri4+=4;
			    }	// if ( fabs(DDPP[ix]) > EPS_ERI )
			}
		    }
		}
	    }
	    if ( nzeri >= max_nzeri ) {
		*last_ijcs = ijcs;
		*last_klcs = klcs;
		*ebuf_non_zero_eri = nzeri;
		return OFMO_EBUF_FULL;
	    }
	}	// for ( klcs );
    }		// for ( ijcs );
    *ebuf_non_zero_eri = nzeri;
    /*
    // for check sum
#pragma omp atomic
    check_sum[17] += lsum;
    */

    return OFMO_EBUF_NOFULL;
}

/** (dd,ds)タイプの２電子積分関数
 * */
int ofmo_twoint_buffer_ddds__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs
	) {
    int Lab, Lcd, i, j, k, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    double DDDS[6*6*6];
    long nzeri, max_nzeri, nzeri4;
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    long ebuf_max_nzeri = *pebuf_max_nzeri;
    /*
    // for check sum
    double lsum = ZERO;
    */

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs0 = leading_cs_pair[Lab];
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];

    nzeri     = *ebuf_non_zero_eri;
    max_nzeri = ebuf_max_nzeri - 6*6*6;
    nzeri4    = nzeri*4;
    if ( nzeri >= max_nzeri ) { 
	*last_ijcs = ijcs0+workerid;
	*last_klcs = klcs0 - 1;
	*ebuf_non_zero_eri = nzeri;
	return OFMO_EBUF_FULL;
    }
    /*
    // debug
    printf("(dd,ds) :  nwks=%d, wkid=%d, max nzeri=%lld\n",
	    nworkers, workerid, max_nzeri );
    fflush(stdout);
    */

    for ( ijcs=ijcs0+workerid; ijcs<ijcs1; ijcs+=nworkers ) {
	val_ab = csp_schwarz[ijcs];
	ics    = csp_ics[ijcs];
	jcs    = csp_jcs[ijcs];
	ijps0  = csp_leading_ps_pair[ijcs];
	nijps  = csp_leading_ps_pair[ijcs+1]-ijps0;
	iat    = shel_atm[ics];
	jat    = shel_atm[jcs];
	iao0   = shel_ini[ics];
	jao0   = shel_ini[jcs];
	A[0]=atom_x[iat]; A[1]=atom_y[iat]; A[2]=atom_z[iat];
	B[0]=atom_x[jat]; B[1]=atom_y[jat]; B[2]=atom_z[jat];
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];
	
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
	    klps0  = csp_leading_ps_pair[klcs];
	    nklps  = csp_leading_ps_pair[klcs+1]-klps0;
	    kat    = shel_atm[kcs];
	    lat    = shel_atm[lcs];
	    kao0   = shel_ini[kcs];
	    lao    = shel_ini[lcs];
	    C[0]=atom_x[kat]; C[1]=atom_y[kat]; C[2]=atom_z[kat];
	    D[0]=atom_x[lat]; D[1]=atom_y[lat]; D[2]=atom_z[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_ddds__(
		    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
		    &psp_xiza[ijps0], BA,
		    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
		    &psp_xiza[klps0], DC,   AC,      DDDS );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		for ( j=0, jao=jao0; j<6; j++, jao++ ) {
		    if ( jao>iao ) { ix += 6; continue; }
		    coe = ( iao==jao ? HALF : ONE );
		    for ( k=0, kao=kao0; k<6; k++, kao++, ix++ ) {
			if ( fabs(DDDS[ix]) > EPS_ERI ) {
			    ebuf_val[nzeri]     = coe*DDDS[ix];
			    /*
		// for check sum
		lsum += ebuf_val[nzeri];
		*/

			    if ( iao >= kao ) {
				ebuf_ind4[nzeri4+0] = (short int)iao;
				ebuf_ind4[nzeri4+1] = (short int)jao;
				ebuf_ind4[nzeri4+2] = (short int)kao;
				ebuf_ind4[nzeri4+3] = (short int)lao;
			    } else {
				ebuf_ind4[nzeri4+0] = (short int)kao;
				ebuf_ind4[nzeri4+1] = (short int)lao;
				ebuf_ind4[nzeri4+2] = (short int)iao;
				ebuf_ind4[nzeri4+3] = (short int)jao;
			    }
			    nzeri++;
			    nzeri4+=4;
			}	// if ( fabs(PPSS[ix]) > EPS_ERI )
		    }
		}
	    }
	    if ( nzeri >= max_nzeri ) {
		*last_ijcs = ijcs;
		*last_klcs = klcs;
		*ebuf_non_zero_eri = nzeri;
		return OFMO_EBUF_FULL;
	    }
	}	// for ( klcs );
    }		// for ( ijcs );
    *ebuf_non_zero_eri = nzeri;
    /*
    // for check sum
#pragma omp atomic
    check_sum[18] += lsum;
    */

    return OFMO_EBUF_NOFULL;
}

/** (dd,dp)タイプの２電子積分関数
 * */
int ofmo_twoint_buffer_dddp__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs
	) {
    int Lab, Lcd, i, j, k, l, ix;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    double DDDP[6*6*6*3];
    long nzeri, max_nzeri, nzeri4;
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    long ebuf_max_nzeri = *pebuf_max_nzeri;
    /*
    // for check sum
    double lsum = ZERO;
    */

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs0 = leading_cs_pair[Lab];
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];

    nzeri     = *ebuf_non_zero_eri;
    max_nzeri = ebuf_max_nzeri - 6*6*6*3;
    nzeri4    = nzeri*4;
    if ( nzeri >= max_nzeri ) { 
	*last_ijcs = ijcs0+workerid;
	*last_klcs = klcs0 - 1;
	*ebuf_non_zero_eri = nzeri;
	return OFMO_EBUF_FULL;
    }
    /*
    // debug
    printf("(dd,dp) :  nwks=%d, wkid=%d, max nzeri=%lld\n",
	    nworkers, workerid, max_nzeri );
    fflush(stdout);
    */

    for ( ijcs=ijcs0+workerid; ijcs<ijcs1; ijcs+=nworkers ) {
	val_ab = csp_schwarz[ijcs];
	ics    = csp_ics[ijcs];
	jcs    = csp_jcs[ijcs];
	ijps0  = csp_leading_ps_pair[ijcs];
	nijps  = csp_leading_ps_pair[ijcs+1]-ijps0;
	iat    = shel_atm[ics];
	jat    = shel_atm[jcs];
	iao0   = shel_ini[ics];
	jao0   = shel_ini[jcs];
	A[0]=atom_x[iat]; A[1]=atom_y[iat]; A[2]=atom_z[iat];
	B[0]=atom_x[jat]; B[1]=atom_y[jat]; B[2]=atom_z[jat];
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];
	
	for ( klcs=klcs0; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
	    klps0  = csp_leading_ps_pair[klcs];
	    nklps  = csp_leading_ps_pair[klcs+1]-klps0;
	    kat    = shel_atm[kcs];
	    lat    = shel_atm[lcs];
	    kao0   = shel_ini[kcs];
	    lao0   = shel_ini[lcs];
	    C[0]=atom_x[kat]; C[1]=atom_y[kat]; C[2]=atom_z[kat];
	    D[0]=atom_x[lat]; D[1]=atom_y[lat]; D[2]=atom_z[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    /*// debug
	    if ( ics==36 && jcs==35 && kcs == 35 && lcs==25 ) {
		DDDP[0] = 1.e0;
	    } else {
		DDDP[0] = 0.e0;
	    }*/
	    twoint_core_dddp__(
		    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
		    &psp_xiza[ijps0], BA,
		    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
		    &psp_xiza[klps0], DC,   AC,      DDDP );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		for ( j=0, jao=jao0; j<6; j++, jao++ ) {
		    if ( jao>iao ) { ix += 6*3; continue; }
		    coe = ( iao==jao ? HALF : ONE );
		    for ( k=0, kao=kao0; k<6; k++, kao++ ) {
			for ( l=0, lao=lao0; l<3; l++, lao++, ix++ ) {
			    if ( fabs(DDDP[ix]) > EPS_ERI ) {
				ebuf_val[nzeri]     = coe*DDDP[ix];
				/*
		// for check sum
		lsum += ebuf_val[nzeri];
		*/

				if ( iao >= kao ) {
				    ebuf_ind4[nzeri4+0] = (short int)iao;
				    ebuf_ind4[nzeri4+1] = (short int)jao;
				    ebuf_ind4[nzeri4+2] = (short int)kao;
				    ebuf_ind4[nzeri4+3] = (short int)lao;
				    /*if ( ics==36 && jcs==35 && kcs == 35 && lcs==25 ) {

				    // debug
				    printf(" %2d %2d %2d %2d : %1d %1d %1d %1d : %2d %2d %2d %2d : %16.12f\n",
					    iao, jao, kao, lao, (i+1), (j+1), (k+1), (l+1),
					    iat, jat, kat, lat, DDDP[ix]);
				    }*/

				} else {
				    ebuf_ind4[nzeri4+0] = (short int)kao;
				    ebuf_ind4[nzeri4+1] = (short int)lao;
				    ebuf_ind4[nzeri4+2] = (short int)iao;
				    ebuf_ind4[nzeri4+3] = (short int)jao;
				    /*// debug
				    printf(" %2d %2d %2d %2d : %1d %1d %1d %1d : %2d %2d %2d %2d : %16.12f\n",
					    kao, lao, iao, jao, (k+1), (l+1), (i+1), (j+1),
					    kat, lat, iat, jat, DDDP[ix]);
					    */
				}
				nzeri++;
				nzeri4+=4;
			    }	// if ( fabs(PPSS[ix]) > EPS_ERI )
			}
		    }
		}
	    }
	    if ( nzeri >= max_nzeri ) {
		*last_ijcs = ijcs;
		*last_klcs = klcs;
		*ebuf_non_zero_eri = nzeri;
		return OFMO_EBUF_FULL;
	    }
	}	// for ( klcs );
    }		// for ( ijcs );
    *ebuf_non_zero_eri = nzeri;
    /*
    // for check sum
#pragma omp atomic
    check_sum[19] += lsum;
    */

    return OFMO_EBUF_NOFULL;
}

/** (dd,dd)タイプの２電子積分関数
 * */
int ofmo_twoint_buffer_dddd__(
	// parallelization
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
	// output
	const long *pebuf_max_nzeri, long *ebuf_non_zero_eri,
	double ebuf_val[], short int ebuf_ind4[],
	int *last_ijcs, int *last_klcs
	) {
    int Lab, Lcd, i, j, k, l, ipat, ix;
    int I2, IJ, K2, KL;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0 ;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe, coe0;
    double DDDD[6*6*6*6];
    long nzeri, max_nzeri, nzeri4;
    //
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    long ebuf_max_nzeri = *pebuf_max_nzeri;
    /*
    // for check sum
    double lsum = ZERO;
    */

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs0 = leading_cs_pair[Lab];
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];

    nzeri     = *ebuf_non_zero_eri;
    max_nzeri = ebuf_max_nzeri - 6*6*6*6;
    nzeri4    = nzeri*4;
    if ( nzeri >= max_nzeri ) { 
	*last_ijcs = ijcs0+workerid;
	*last_klcs = klcs0 - 1;
	*ebuf_non_zero_eri = nzeri;
	return OFMO_EBUF_FULL;
    }

    for ( ijcs=ijcs0+workerid; ijcs<ijcs1; ijcs+=nworkers ) {
	val_ab = csp_schwarz[ijcs];
	ics    = csp_ics[ijcs];
	jcs    = csp_jcs[ijcs];
	ijps0  = csp_leading_ps_pair[ijcs];
	nijps  = csp_leading_ps_pair[ijcs+1]-ijps0;
	iat    = shel_atm[ics];
	jat    = shel_atm[jcs];
	iao0   = shel_ini[ics];
	jao0   = shel_ini[jcs];
	A[0]=atom_x[iat]; A[1]=atom_y[iat]; A[2]=atom_z[iat];
	B[0]=atom_x[jat]; B[1]=atom_y[jat]; B[2]=atom_z[jat];
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];
	
	for ( klcs=klcs0; klcs<=ijcs; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < EPS_PS4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
	    klps0  = csp_leading_ps_pair[klcs];
	    nklps  = csp_leading_ps_pair[klcs+1]-klps0;
	    kat    = shel_atm[kcs];
	    lat    = shel_atm[lcs];
	    kao0   = shel_ini[kcs];
	    lao0   = shel_ini[lcs];
	    C[0]=atom_x[kat]; C[1]=atom_y[kat]; C[2]=atom_z[kat];
	    D[0]=atom_x[lat]; D[1]=atom_y[lat]; D[2]=atom_z[lat];
	    for ( i=0; i<3; i++ ) {
		AC[i] = A[i] - C[i];
		DC[i] = D[i] - C[i];
	    }
	    twoint_core_dddd__(
		    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
		    &psp_xiza[ijps0], BA,
		    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
		    &psp_xiza[klps0], DC,   AC,      DDDD );
	    ipat = ( (ics==kcs && jcs>lcs) ? true : false );
#ifdef SORT_CSP
            int ijgekl = (ics>kcs);
            if (ics==kcs) ijgekl = (jcs>=lcs);
            if (!ijgekl) ipat = ( (ics==kcs && jcs<lcs) ? true : false);
#endif
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		I2 = iao*(iao+1)/2;
		for ( j=0, jao=jao0; j<6; j++, jao++ ) {
		    if ( jao>iao ) { ix+=6*6; continue; }
		    IJ = I2 + jao;
		    coe0 = ( iao==jao ? HALF : ONE );
		    for ( k=0, kao=kao0; k<6; k++, kao++ ) {
			K2 = kao*(kao+1)/2;
			for ( l=0, lao=lao0; l<6; l++, lao++, ix++ ) {
			    if ( lao>kao ) continue;
			    if ( fabs(DDDD[ix]) > EPS_ERI ) {
				KL = K2 + lao;
#ifndef SORT_CSP
				if ( IJ >= KL ) {
#else
                    if ( (ijgekl&&IJ>=KL) || (!ijgekl&&KL>=IJ) ) {
#endif
				    coe = coe0;
				    if ( kao==lao ) coe *= HALF;
				    if ( KL == IJ ) coe *= HALF;
				    ebuf_val[nzeri]     = coe*DDDD[ix];
				    /*
		// for check sum
		lsum += ebuf_val[nzeri];
		*/

				    ebuf_ind4[nzeri4+0] = (short int)iao;
				    ebuf_ind4[nzeri4+1] = (short int)jao;
				    ebuf_ind4[nzeri4+2] = (short int)kao;
				    ebuf_ind4[nzeri4+3] = (short int)lao;
				    nzeri++;
				    nzeri4+=4;
				} else if ( ipat ) {
				    coe = coe0;
				    if ( kao==lao ) coe*=HALF;
				    ebuf_val[nzeri]     = coe*DDDD[ix];
				    /*
		// for check sum
		lsum += ebuf_val[nzeri];
		*/

				    ebuf_ind4[nzeri4+0] = (short int)kao;
				    ebuf_ind4[nzeri4+1] = (short int)lao;
				    ebuf_ind4[nzeri4+2] = (short int)iao;
				    ebuf_ind4[nzeri4+3] = (short int)jao;
				    nzeri++;
				    nzeri4+=4;
				}
			    }
			}	// l
		    }		// k
		}		// j
	    }			// i
	    if ( nzeri >= max_nzeri ) {
		*last_ijcs = ijcs;
		*last_klcs = klcs;
		*ebuf_non_zero_eri = nzeri;
		return OFMO_EBUF_FULL;
	    }
	}	// for ( klcs );
    }		// for ( ijcs );
    *ebuf_non_zero_eri = nzeri;
    /*
    // for check sum
#pragma omp atomic
    check_sum[20] += lsum;
    */

    return OFMO_EBUF_NOFULL;
}
