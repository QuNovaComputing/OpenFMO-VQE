/**
 * @file ofmo-twoint-direct.c
 * 各タイプの二電子積分を計算して２電子ハミルトン行列（G行列）を
 * 計算する
 * この関数では、積分計算とG行列計算の両方を行っている
 * */
/**
 * @defgroup integ-twoint-direct direct法によるG行列作成を行う関数群
 * buffered direct SCF法において、バッファに保存出来なかった２電子積分を
 * 計算しながらG行列生成を行う関数群
 *
 * 各関数は、同じ引数を持っているので、以下にその内容を記す。
 *
 * @param[in] nworkers 計算に用いるワーカプロセス（スレッド）数
 * @param[in] workerid 各ワーカプロセス（スレッド）のID
 *     （\f$ 0\le\tt{workerid}<\tt{nworkers} \f$）
 * @param[in] ebuf_max_nzeri バッファに保存できる最大積分数
 * @param[in] La １つ目の軌道量子数
 * @param[in] Lb ２つ目の軌道量子数
 *     （\f$ \tt{La} \ge \tt{Lb} \f$）
 * @param[in] Lc ３つ目の軌道量子数
 *     （\f$ \tt{La} \ge \tt{Lc} \f$）
 * @param[in] Ld ４つ目の軌道量子数
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
 * @param[in] last_ijcs 積分計算を始める際の、外側CSペアループ変数の初期値
 * @param[in] last_klcs 積分計算を始める際の、外側CSペアループ変数の初期値
 *
 * @param[in] nao AO数
 * @param[in] Ds[] 密度行列（正方行列形式）
 * @param[out] G[] 計算されたG行列要素を格納する配列（正方行列形式）
 *
 * @ingroup integ-med
 * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ofmo-def.h"

#define OFMO_EBUF_FULL		-1
#define OFMO_EBUF_NOFULL	0

#define EPS_PS4 1.e-30
#define EPS_ERI	1.e-15

#define ZERO	0.e0
#define ONE	1.e0
#define HALF	.5e0

#include "ofmo-twoint-core.h"
#include "ofmo-twoint.h"

extern int ofmo_integ_add_fock( const int nao, const size_t nstored_eri,
	const double eri_val[], const short int eri_ind4[],
	const double D[], double G[] );
extern float ofmo_twoint_dmax6(const int i, const int j, const int k, const int l);

/** (ss,ss)タイプの積分計算、および、G行列計算を行う関数
 *
 * @ingroup integ-twoint-direct
 * */
int ofmo_twoint_direct_ssss__(
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
	const int *pnao, const double Ds[], double G[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int last_ijcs=*plast_ijcs, last_klcs=*plast_klcs, nao=*pnao;
    long max_nzeri=*petmp_max_nzeri;
    long nzeri4, nzeri=*petmp_non_zero_eri;
    int Lab, Lcd, i;
    int ijcs,        ijcs1, ijao;
    int klcs, klcs0,        klao;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, jcs, jat, jao, kcs, kat, kao, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);
    float eps_sch = ofmo_twoint_eps_sch(0);
    double SSSS[1];

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    if ( last_ijcs != -1 ) {
	ijcs = last_ijcs;
	klcs = last_klcs+1;
    } else {
	ijcs = leading_cs_pair[Lab] + workerid;
	klcs = klcs0;
    }
    max_nzeri-=0;
    if ( nzeri>= max_nzeri ) {
	ofmo_integ_add_fock( nao, (size_t)nzeri, etmp_val, etmp_ind4, Ds, G );
	nzeri = 0;
    }
    nzeri4 = nzeri<<2;

    for (  ; ijcs<ijcs1; ijcs+=nworkers ) {
	val_ab = csp_schwarz[ijcs];
	ics    = csp_ics[ijcs];
	jcs    = csp_jcs[ijcs];
	ijps0  = csp_leading_ps_pair[ijcs];
	nijps  = csp_leading_ps_pair[ijcs+1]-ijps0;
	iat    = shel_atm[ics];
	jat    = shel_atm[jcs];
	iao    = shel_ini[ics];
	jao    = shel_ini[jcs];
	ijao   = iao*(iao+1)/2+jao;
	A[0]=atom_x[iat]; A[1]=atom_y[iat]; A[2]=atom_z[iat];
	B[0]=atom_x[jat]; B[1]=atom_y[jat]; B[2]=atom_z[jat];
	for ( i=0; i<3; i++ ) BA[i] = B[i] - A[i];
	
	for ( ; klcs<=ijcs; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < eps_ps4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
            float dmax = ofmo_twoint_dmax6(ics,jcs,kcs,lcs);
            if ( dmax*val_ab*val_cd < eps_sch ) continue;
	    klps0  = csp_leading_ps_pair[klcs];
	    nklps  = csp_leading_ps_pair[klcs+1]-klps0;
	    kat    = shel_atm[kcs];
	    lat    = shel_atm[lcs];
	    kao    = shel_ini[kcs];
	    lao    = shel_ini[lcs];
	    klao   = kao*(kao+1)/2+lao;
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
	    if ( fabs(SSSS[0]) > eps_eri ) {
		double x;
		coe = ONE;
		if ( iao == jao ) coe = HALF;
		if ( kao == lao ) coe *= HALF;
		if ( ijao == klao ) coe *= HALF;
		x  = coe * SSSS[0];
		etmp_val[nzeri] = x;
		etmp_ind4[nzeri4+0] = iao;
		etmp_ind4[nzeri4+1] = jao;
		etmp_ind4[nzeri4+2] = kao;
		etmp_ind4[nzeri4+3] = lao;
		nzeri++;
		nzeri4+=4;
	    }	// if ( fabs(SSSS[0]) > eps_eri )
	    if ( nzeri >= max_nzeri ) {
		ofmo_integ_add_fock( nao, (size_t)nzeri, etmp_val, etmp_ind4,
			Ds, G );
		nzeri = nzeri4= 0;
	    }
	}	// for ( klcs );
	klcs = klcs0;
    }		// for ( ijcs );
    *petmp_non_zero_eri = nzeri;
    return 0;
}

/** (ps,ss)タイプの積分計算、および、G行列計算を行う関数
 *
 * @ingroup integ-twoint-direct
 * */
int ofmo_twoint_direct_psss__(
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
	const int *pnao, const double Ds[], double G[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int last_ijcs=*plast_ijcs, last_klcs=*plast_klcs, nao=*pnao;
    long max_nzeri=*petmp_max_nzeri;
    long nzeri4, nzeri=*petmp_non_zero_eri;
    int Lab, Lcd, i;
    int ijcs,        ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    int kcs, kat, kao, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);
    float eps_sch = ofmo_twoint_eps_sch(0);
    double PSSS[3];

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];
    if ( last_ijcs != -1 ) {
	ijcs = last_ijcs;
	klcs = last_klcs+1;
    } else {
	ijcs = leading_cs_pair[Lab] + workerid;
	klcs = klcs0;
    }
    max_nzeri-=3;
    if ( nzeri>= max_nzeri ) {
	ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4, Ds, G );
	nzeri = 0;
    }
    nzeri4 = nzeri<<2;

    for (  ; ijcs<ijcs1; ijcs+=nworkers ) {
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
	
	for ( ; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < eps_ps4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
            float dmax = ofmo_twoint_dmax6(ics,jcs,kcs,lcs);
            if ( dmax*val_ab*val_cd < eps_sch ) continue;
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
		if ( fabs(PSSS[i]) > eps_eri ) {
		    double x;
		    x  = coe * PSSS[i];
		    etmp_val[nzeri] = x;
		    etmp_ind4[nzeri4+0] = iao;
		    etmp_ind4[nzeri4+1] = jao;
		    etmp_ind4[nzeri4+2] = kao;
		    etmp_ind4[nzeri4+3] = lao;
		    nzeri++;
		    nzeri4+=4;
		}
	    }	// for ( i, iao);
	    if ( nzeri >= max_nzeri ) {
		ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4,
			Ds, G );
		nzeri = nzeri4= 0;
	    }
	}	// for ( klcs );
	klcs = klcs0;
    }		// for ( ijcs );
    *petmp_non_zero_eri = nzeri;
    return 0;
}

/** (ps,ps)タイプの積分計算、および、G行列計算を行う関数
 *
 * @ingroup integ-twoint-direct
 * */
int ofmo_twoint_direct_psps__(
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
	const int *pnao, const double Ds[], double G[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int last_ijcs=*plast_ijcs, last_klcs=*plast_klcs, nao=*pnao;
    long max_nzeri=*petmp_max_nzeri;
    long nzeri4, nzeri=*petmp_non_zero_eri;
    int Lab, Lcd, i, k, ix;
    int ipat, IJ, KL;
    int ijcs,        ijcs1;
    int klcs, klcs0;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);
    float eps_sch = ofmo_twoint_eps_sch(0);
    double PSPS[3*3];

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    if ( last_ijcs != -1 ) {
	ijcs = last_ijcs;
	klcs = last_klcs+1;
    } else {
	ijcs = leading_cs_pair[Lab] + workerid;
	klcs = klcs0;
    }
    max_nzeri-=3*3;
    if ( nzeri>= max_nzeri ) {
	ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4, Ds, G );
	nzeri = 0;
    }
    nzeri4 = nzeri<<2;

    for (  ; ijcs<ijcs1; ijcs+=nworkers ) {
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
	
	for ( ; klcs<=ijcs; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < eps_ps4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
            float dmax = ofmo_twoint_dmax6(ics,jcs,kcs,lcs);
            if ( dmax*val_ab*val_cd < eps_sch ) continue;
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
	    ipat = ( (ics==kcs && jcs>lcs) ? true : false);
	    for ( i=0, iao=iao0, ix=0; i<3; i++, iao++ ) {
		IJ = ((iao*iao+iao)>>1) + jao;
		for ( k=0, kao=kao0; k<3; k++, kao++, ix++ ) {
		    KL = ((kao*kao+kao)>>1) + lao ;
		    if ( fabs(PSPS[ix]) <= eps_eri ) continue;
		    if ( IJ>=KL || ipat ) {
			double x;
			coe = ONE;
			if ( IJ == KL ) coe = HALF;
			x  = coe * PSPS[ix];
			etmp_val[nzeri] = x;
			etmp_ind4[nzeri4+0] = iao;
			etmp_ind4[nzeri4+1] = jao;
			etmp_ind4[nzeri4+2] = kao;
			etmp_ind4[nzeri4+3] = lao;
			nzeri++;
			nzeri4+=4;
		    }
		}	// for ( kao )
	    }		// for ( iao )
	    if ( nzeri >= max_nzeri ) {
		ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4,
			Ds, G );
		nzeri = nzeri4= 0;
	    }
	}	// for ( klcs );
	klcs = klcs0;
    }		// for ( ijcs );
    *petmp_non_zero_eri = nzeri;
    return 0;
}

/** (pp,ss)タイプの積分計算、および、G行列計算を行う関数
 *
 * @ingroup integ-twoint-direct
 * */
int ofmo_twoint_direct_ppss__(
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
	const int *pnao, const double Ds[], double G[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int last_ijcs=*plast_ijcs, last_klcs=*plast_klcs, nao=*pnao;
    long max_nzeri=*petmp_max_nzeri;
    long nzeri4, nzeri=*petmp_non_zero_eri;
    int Lab, Lcd, i, j, ix;
    int ijcs,        ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe, coe0;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);
    float eps_sch = ofmo_twoint_eps_sch(0);
    double PPSS[3*3];

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];
    if ( last_ijcs != -1 ) {
	ijcs = last_ijcs;
	klcs = last_klcs+1;
    } else {
	ijcs = leading_cs_pair[Lab] + workerid;
	klcs = klcs0;
    }
    max_nzeri-=3*3;
    if ( nzeri>= max_nzeri ) {
	ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4, Ds, G );
	nzeri = 0;
    }
    nzeri4 = nzeri<<2;

    for (  ; ijcs<ijcs1; ijcs+=nworkers ) {
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
	
	for ( ; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < eps_ps4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
            float dmax = ofmo_twoint_dmax6(ics,jcs,kcs,lcs);
            if ( dmax*val_ab*val_cd < eps_sch ) continue;
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
	    coe0 = ( kao == lao ? HALF : ONE );
	    for ( i=0, iao=iao0, ix=0; i<3; i++, iao++ ) {
		for ( j=0, jao=jao0; j<3; j++, jao++, ix++ ) {
		    if ( jao > iao ) continue;
		    if ( fabs(PPSS[ix]) > eps_eri ) {
			double x;
			coe = coe0;
			if ( iao == jao ) coe *= HALF;
			x  = coe * PPSS[ix];
			etmp_val[nzeri] = x;
			etmp_ind4[nzeri4+0] = iao;
			etmp_ind4[nzeri4+1] = jao;
			etmp_ind4[nzeri4+2] = kao;
			etmp_ind4[nzeri4+3] = lao;
			nzeri++;
			nzeri4+=4;
		    }
		} // for ( jao );
	    }	// for ( iao );
	    if ( nzeri >= max_nzeri ) {
		ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4,
			Ds, G );
		nzeri = nzeri4= 0;
	    }
	}	// for ( klcs );
	klcs = klcs0;
    }		// for ( ijcs );
    *petmp_non_zero_eri = nzeri;
    return 0;
}

/** (pp,ps)タイプの積分計算、および、G行列計算を行う関数
 *
 * @ingroup integ-twoint-direct
 * */
int ofmo_twoint_direct_ppps__(
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
	const int *pnao, const double Ds[], double G[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int last_ijcs=*plast_ijcs, last_klcs=*plast_klcs, nao=*pnao;
    long max_nzeri=*petmp_max_nzeri;
    long nzeri4, nzeri=*petmp_non_zero_eri;
    int Lab, Lcd, i, j, k, ix;
    int ijcs,        ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);
    float eps_sch = ofmo_twoint_eps_sch(0);
    double PPPS[3*3*3];

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];
    if ( last_ijcs != -1 ) {
	ijcs = last_ijcs;
	klcs = last_klcs+1;
    } else {
	ijcs = leading_cs_pair[Lab] + workerid;
	klcs = klcs0;
    }
    max_nzeri-=3*3*3;
    if ( nzeri>= max_nzeri ) {
	ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4, Ds, G );
	nzeri = 0;
    }
    nzeri4 = nzeri<<2;

    for (  ; ijcs<ijcs1; ijcs+=nworkers ) {
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
	
	for ( ; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < eps_ps4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
            float dmax = ofmo_twoint_dmax6(ics,jcs,kcs,lcs);
            if ( dmax*val_ab*val_cd < eps_sch ) continue;
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
			if ( fabs(PPPS[ix]) > eps_eri ) {
			    double x;
			    x  = coe * PPPS[ix];
			    etmp_val[nzeri] = x;
			    etmp_ind4[nzeri4+0] = iao;
			    etmp_ind4[nzeri4+1] = jao;
			    etmp_ind4[nzeri4+2] = kao;
			    etmp_ind4[nzeri4+3] = lao;
			    nzeri++;
			    nzeri4+=4;
			}	// if ( fabs );
		    }	// for ( kao )
		}	// for ( jao )
	    }	// for ( iao )
	    if ( nzeri >= max_nzeri ) {
		ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4,
			Ds, G );
		nzeri = nzeri4= 0;
	    }
	}	// for ( klcs );
	klcs = klcs0;
    }		// for ( ijcs );
    *petmp_non_zero_eri = nzeri;
    return 0;
}

/** (pp,pp)タイプの積分計算、および、G行列計算を行う関数
 *
 * @ingroup integ-twoint-direct
 * */
int ofmo_twoint_direct_pppp__(
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
	const int *pnao, const double Ds[], double G[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int last_ijcs=*plast_ijcs, last_klcs=*plast_klcs, nao=*pnao;
    long max_nzeri=*petmp_max_nzeri;
    long nzeri4, nzeri=*petmp_non_zero_eri;
    int Lab, Lcd, i, j, k, l, ipat;
    int I2, IJ, K2, KL;
    int ijcs,        ijcs1;
    int klcs, klcs0;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe, coe0;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);
    float eps_sch = ofmo_twoint_eps_sch(0);
    double PPPP[3*3*3*3], pppp;

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    if ( last_ijcs != -1 ) {
	ijcs = last_ijcs;
	klcs = last_klcs+1;
    } else {
	ijcs = leading_cs_pair[Lab] + workerid;
	klcs = klcs0;
    }
    max_nzeri-=3*3*3*3;
    if ( nzeri>= max_nzeri ) {
	ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4, Ds, G );
	nzeri = 0;
    }
    nzeri4 = nzeri<<2;

    for (  ; ijcs<ijcs1; ijcs+=nworkers ) {
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
	
	for ( ; klcs<=ijcs; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < eps_ps4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
            float dmax = ofmo_twoint_dmax6(ics,jcs,kcs,lcs);
            if ( dmax*val_ab*val_cd < eps_sch ) continue;
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
			    if ( fabs(pppp) > eps_eri ) {
				double x;
				KL = K2 + lao;
				if ( IJ >= KL  || ipat ) {
				    coe = coe0;
				    if ( kao==lao ) coe *= HALF;
				    if ( KL == IJ ) coe *= HALF;
				    x  = coe * pppp;
				    etmp_val[nzeri] = x;
				    etmp_ind4[nzeri4+0] = iao;
				    etmp_ind4[nzeri4+1] = jao;
				    etmp_ind4[nzeri4+2] = kao;
				    etmp_ind4[nzeri4+3] = lao;
				    nzeri++;
				    nzeri4+=4;
				}
			    } // if ( fabs )
			} // for ( lao )
		    } // for ( kao )
		} // for ( jao )
	    }	// for ( iao )
	    if ( nzeri >= max_nzeri ) {
		ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4,
			Ds, G );
		nzeri = nzeri4= 0;
	    }
	}	// for ( klcs );
	klcs = klcs0;
    }		// for ( ijcs );
    *petmp_non_zero_eri = nzeri;
    return 0;
}

/** (ds,ss)タイプの積分計算、および、G行列計算を行う関数
 *
 * @ingroup integ-twoint-direct
 * */
int ofmo_twoint_direct_dsss__(
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
	const int *pnao, const double Ds[], double G[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int last_ijcs=*plast_ijcs, last_klcs=*plast_klcs, nao=*pnao;
    long max_nzeri=*petmp_max_nzeri;
    long nzeri4, nzeri=*petmp_non_zero_eri;
    int Lab, Lcd, i;
    int ijcs,        ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    int kcs, kat, kao, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);
    float eps_sch = ofmo_twoint_eps_sch(0);
    double DSSS[6];

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];
    if ( last_ijcs != -1 ) {
	ijcs = last_ijcs;
	klcs = last_klcs+1;
    } else {
	ijcs = leading_cs_pair[Lab] + workerid;
	klcs = klcs0;
    }
    max_nzeri-=6;
    if ( nzeri>= max_nzeri ) {
	ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4, Ds, G );
	nzeri = 0;
    }
    nzeri4 = nzeri<<2;

    for (  ; ijcs<ijcs1; ijcs+=nworkers ) {
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
	
	for ( ; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < eps_ps4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
            float dmax = ofmo_twoint_dmax6(ics,jcs,kcs,lcs);
            if ( dmax*val_ab*val_cd < eps_sch ) continue;
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
	    coe = ( kao == lao ? HALF : ONE );
	    for ( i=0, iao=iao0; i<6; i++, iao++ ) {
		if ( fabs(DSSS[i]) > eps_eri ) {
		    double x;
		    x  = coe * DSSS[i];
		    etmp_val[nzeri] = x;
		    etmp_ind4[nzeri4+0] = iao;
		    etmp_ind4[nzeri4+1] = jao;
		    etmp_ind4[nzeri4+2] = kao;
		    etmp_ind4[nzeri4+3] = lao;
		    nzeri++;
		    nzeri4+=4;
		}
	    }	// for ( i, iao);
	    if ( nzeri >= max_nzeri ) {
		ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4,
			Ds, G );
		nzeri = nzeri4= 0;
	    }
	}	// for ( klcs );
	klcs = klcs0;
    }		 // for ( ijcs );
    /*
    // for check sum
#pragma omp atomic
    check_sum[6] += lsum;
    */

    *petmp_non_zero_eri = nzeri;
    return 0;
}

/** (ds,ps)タイプの積分計算、および、G行列計算を行う関数
 *
 * @ingroup integ-twoint-direct
 * */
int ofmo_twoint_direct_dsps__(
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
	const int *pnao, const double Ds[], double G[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int last_ijcs=*plast_ijcs, last_klcs=*plast_klcs, nao=*pnao;
    long max_nzeri=*petmp_max_nzeri;
    long nzeri4, nzeri=*petmp_non_zero_eri;
    int Lab, Lcd, i, k, ix;
    int ijcs,        ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);
    float eps_sch = ofmo_twoint_eps_sch(0);
    double DSPS[6*3];

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];
    if ( last_ijcs != -1 ) {
	ijcs = last_ijcs;
	klcs = last_klcs+1;
    } else {
	ijcs = leading_cs_pair[Lab] + workerid;
	klcs = klcs0;
    }
    max_nzeri-=6*3;
    if ( nzeri>= max_nzeri ) {
	ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4, Ds, G );
	nzeri = 0;
    }
    nzeri4 = nzeri<<2;

    for (  ; ijcs<ijcs1; ijcs+=nworkers ) {
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
	
	for ( ; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < eps_ps4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
            float dmax = ofmo_twoint_dmax6(ics,jcs,kcs,lcs);
            if ( dmax*val_ab*val_cd < eps_sch ) continue;
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
		for (k=0, kao=kao0; k<3; k++, kao++, ix++ ) {
		    if ( fabs(DSPS[ix]) > eps_eri ) {
			double x;
			x  = DSPS[ix];
			etmp_val[nzeri] = x;
			etmp_ind4[nzeri4+0] = iao;
			etmp_ind4[nzeri4+1] = jao;
			etmp_ind4[nzeri4+2] = kao;
			etmp_ind4[nzeri4+3] = lao;
			nzeri++;
			nzeri4+=4;
		    }
		}
	    }	// for ( i, iao);
	    if ( nzeri >= max_nzeri ) {
		ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4,
			Ds, G );
		nzeri = nzeri4= 0;
	    }
	}	// for ( klcs );
	klcs = klcs0;
    }		// for ( ijcs );
    *petmp_non_zero_eri = nzeri;
    return 0;
}

/** (ds,pp)タイプの積分計算、および、G行列計算を行う関数
 *
 * @ingroup integ-twoint-direct
 * */
int ofmo_twoint_direct_dspp__(
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
	const int *pnao, const double Ds[], double G[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int last_ijcs=*plast_ijcs, last_klcs=*plast_klcs, nao=*pnao;
    long max_nzeri=*petmp_max_nzeri;
    long nzeri4, nzeri=*petmp_non_zero_eri;
    int Lab, Lcd, i, k, l, ix;
    int ijcs,        ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);
    float eps_sch = ofmo_twoint_eps_sch(0);
    double DSPP[6*3*3];

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];
    if ( last_ijcs != -1 ) {
	ijcs = last_ijcs;
	klcs = last_klcs+1;
    } else {
	ijcs = leading_cs_pair[Lab] + workerid;
	klcs = klcs0;
    }
    max_nzeri-=6*3*3;
    if ( nzeri>= max_nzeri ) {
	ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4, Ds, G );
	nzeri = 0;
    }
    nzeri4 = nzeri<<2;

    for (  ; ijcs<ijcs1; ijcs+=nworkers ) {
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
	
	for ( ; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < eps_ps4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
            float dmax = ofmo_twoint_dmax6(ics,jcs,kcs,lcs);
            if ( dmax*val_ab*val_cd < eps_sch ) continue;
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
			if ( fabs(DSPP[ix]) > eps_eri ) {
			    double x;
			    x  = coe * DSPP[ix];
			    etmp_val[nzeri] = x;
			    etmp_ind4[nzeri4+0] = iao;
			    etmp_ind4[nzeri4+1] = jao;
			    etmp_ind4[nzeri4+2] = kao;
			    etmp_ind4[nzeri4+3] = lao;
			    nzeri++;
			    nzeri4+=4;
			}
		    }
		}
	    }	// for ( i, iao);
	    if ( nzeri >= max_nzeri ) {
		ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4,
			Ds, G );
		nzeri = nzeri4= 0;
	    }
	}	// for ( klcs );
	klcs = klcs0;
    }		// for ( ijcs );
    *petmp_non_zero_eri = nzeri;
    return 0;
}

/** (ds,ds)タイプの積分計算、および、G行列計算を行う関数
 *
 * @ingroup integ-twoint-direct
 * */
int ofmo_twoint_direct_dsds__(
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
	const int *pnao, const double Ds[], double G[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int last_ijcs=*plast_ijcs, last_klcs=*plast_klcs, nao=*pnao;
    long max_nzeri=*petmp_max_nzeri;
    long nzeri4, nzeri=*petmp_non_zero_eri;
    int Lab, Lcd, i, k, ix;
    int IJ, KL, ipat;
    int ijcs,        ijcs1;
    int klcs, klcs0;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);
    float eps_sch = ofmo_twoint_eps_sch(0);
    double DSDS[6*6];

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    if ( last_ijcs != -1 ) {
	ijcs = last_ijcs;
	klcs = last_klcs+1;
    } else {
	ijcs = leading_cs_pair[Lab] + workerid;
	klcs = klcs0;
    }
    max_nzeri-=6*6;
    if ( nzeri>= max_nzeri ) {
	ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4, Ds, G );
	nzeri = 0;
    }
    nzeri4 = nzeri<<2;

    for (  ; ijcs<ijcs1; ijcs+=nworkers ) {
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
	
	for ( ; klcs<=ijcs; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < eps_ps4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
            float dmax = ofmo_twoint_dmax6(ics,jcs,kcs,lcs);
            if ( dmax*val_ab*val_cd < eps_sch ) continue;
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
	    ipat = ( (ics==kcs && jcs>lcs) ? true : false);
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		IJ = ( (iao*iao+iao)>>1) + jao;
		for ( k=0, kao=kao0; k<6; k++, kao++, ix++ ) {
		    KL = ( (kao*kao+kao)>>1) + lao;
		    if ( fabs(DSDS[ix]) <= eps_eri ) continue;
		    if ( IJ >= KL || ipat ) {
			double x;
			coe = ONE;
			if ( iao==kao && jao==lao ) coe = HALF;
			x  = coe * DSDS[ix];
			etmp_val[nzeri] = x;
			etmp_ind4[nzeri4+0] = iao;
			etmp_ind4[nzeri4+1] = jao;
			etmp_ind4[nzeri4+2] = kao;
			etmp_ind4[nzeri4+3] = lao;
			nzeri++;
			nzeri4+=4;
		    }
		}	// for ( kao )
	    }		// for ( iao )
	    if ( nzeri >= max_nzeri ) {
		ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4,
			Ds, G );
		nzeri = nzeri4= 0;
	    }
	}	// for ( klcs );
	klcs = klcs0;
    }		// for ( ijcs );
    *petmp_non_zero_eri = nzeri;
    return 0;
}

/** (dp,ss)タイプの積分計算、および、G行列計算を行う関数
 *
 * @ingroup integ-twoint-direct
 * */
int ofmo_twoint_direct_dpss__(
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
	const int *pnao, const double Ds[], double G[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int last_ijcs=*plast_ijcs, last_klcs=*plast_klcs, nao=*pnao;
    long max_nzeri=*petmp_max_nzeri;
    long nzeri4, nzeri=*petmp_non_zero_eri;
    int Lab, Lcd, i, j, ix;
    int ijcs,        ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);
    float eps_sch = ofmo_twoint_eps_sch(0);
    double DPSS[6*3];

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];
    if ( last_ijcs != -1 ) {
	ijcs = last_ijcs;
	klcs = last_klcs+1;
    } else {
	ijcs = leading_cs_pair[Lab] + workerid;
	klcs = klcs0;
    }
    max_nzeri-=6*3;
    if ( nzeri>= max_nzeri ) {
	ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4, Ds, G );
	nzeri = 0;
    }
    nzeri4 = nzeri<<2;

    for (  ; ijcs<ijcs1; ijcs+=nworkers ) {
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
	
	for ( ; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < eps_ps4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
            float dmax = ofmo_twoint_dmax6(ics,jcs,kcs,lcs);
            if ( dmax*val_ab*val_cd < eps_sch ) continue;
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
		    if ( fabs(DPSS[ix]) > eps_eri ) {
			double x;
			x  = coe * DPSS[ix];
			etmp_val[nzeri] = x;
			etmp_ind4[nzeri4+0] = iao;
			etmp_ind4[nzeri4+1] = jao;
			etmp_ind4[nzeri4+2] = kao;
			etmp_ind4[nzeri4+3] = lao;
			nzeri++;
			nzeri4+=4;
		    }
		}
	    }	// for ( i, iao);
	    if ( nzeri >= max_nzeri ) {
		ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4,
			Ds, G );
		nzeri = nzeri4= 0;
	    }
	}	// for ( klcs );
	klcs = klcs0;
    }		// for ( ijcs );
    *petmp_non_zero_eri = nzeri;
    return 0;
}

/** (dp,ps)タイプの積分計算、および、G行列計算を行う関数
 *
 * @ingroup integ-twoint-direct
 * */
int ofmo_twoint_direct_dpps__(
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
	const int *pnao, const double Ds[], double G[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int last_ijcs=*plast_ijcs, last_klcs=*plast_klcs, nao=*pnao;
    long max_nzeri=*petmp_max_nzeri;
    long nzeri4, nzeri=*petmp_non_zero_eri;
    int Lab, Lcd, i, j, k, ix;
    int ijcs,        ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);
    float eps_sch = ofmo_twoint_eps_sch(0);
    double DPPS[6*3*3];

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];
    if ( last_ijcs != -1 ) {
	ijcs = last_ijcs;
	klcs = last_klcs+1;
    } else {
	ijcs = leading_cs_pair[Lab] + workerid;
	klcs = klcs0;
    }
    max_nzeri-=6*3*3;
    if ( nzeri>= max_nzeri ) {
	ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4, Ds, G );
	nzeri = 0;
    }
    nzeri4 = nzeri<<2;

    for (  ; ijcs<ijcs1; ijcs+=nworkers ) {
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
	
	for ( ; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < eps_ps4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
            float dmax = ofmo_twoint_dmax6(ics,jcs,kcs,lcs);
            if ( dmax*val_ab*val_cd < eps_sch ) continue;
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
		    for (k=0, kao=kao0; k<3; k++, kao++, ix++ ) {
			if ( fabs(DPPS[ix]) > eps_eri ) {
			    double x;
			    x  = DPPS[ix];
			    etmp_val[nzeri] = x;
			    etmp_ind4[nzeri4+0] = iao;
			    etmp_ind4[nzeri4+1] = jao;
			    etmp_ind4[nzeri4+2] = kao;
			    etmp_ind4[nzeri4+3] = lao;
			    nzeri++;
			    nzeri4+=4;
			}
		    }
		}
	    }	// for ( i, iao);
	    if ( nzeri >= max_nzeri ) {
		ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4,
			Ds, G );
		nzeri = nzeri4= 0;
	    }
	}	// for ( klcs );
	klcs = klcs0;
    }		// for ( ijcs );
    *petmp_non_zero_eri = nzeri;
    return 0;
}

/** (dp,pp)タイプの積分計算、および、G行列計算を行う関数
 *
 * @ingroup integ-twoint-direct
 * */
int ofmo_twoint_direct_dppp__(
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
	const int *pnao, const double Ds[], double G[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int last_ijcs=*plast_ijcs, last_klcs=*plast_klcs, nao=*pnao;
    long max_nzeri=*petmp_max_nzeri;
    long nzeri4, nzeri=*petmp_non_zero_eri;
    int Lab, Lcd, i, j, k, l, ix;
    int ijcs,        ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);
    float eps_sch = ofmo_twoint_eps_sch(0);
    double DPPP[6*3*3*3];

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];
    if ( last_ijcs != -1 ) {
	ijcs = last_ijcs;
	klcs = last_klcs+1;
    } else {
	ijcs = leading_cs_pair[Lab] + workerid;
	klcs = klcs0;
    }
    max_nzeri-=6*3*3*3;
    if ( nzeri>= max_nzeri ) {
	ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4, Ds, G );
	nzeri = 0;
    }
    nzeri4 = nzeri<<2;

    for (  ; ijcs<ijcs1; ijcs+=nworkers ) {
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
	
	for ( ; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < eps_ps4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
            float dmax = ofmo_twoint_dmax6(ics,jcs,kcs,lcs);
            if ( dmax*val_ab*val_cd < eps_sch ) continue;
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
			    if ( fabs(DPPP[ix]) > eps_eri ) {
				double x;
				x  = coe * DPPP[ix];
				etmp_val[nzeri] = x;
				etmp_ind4[nzeri4+0] = iao;
				etmp_ind4[nzeri4+1] = jao;
				etmp_ind4[nzeri4+2] = kao;
				etmp_ind4[nzeri4+3] = lao;
				nzeri++;
				nzeri4+=4;
			    }
			}
		    }
		}
	    }	// for ( i, iao);
	    if ( nzeri >= max_nzeri ) {
		ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4,
			Ds, G );
		nzeri = nzeri4= 0;
	    }
	}	// for ( klcs );
	klcs = klcs0;
    }		// for ( ijcs );
    *petmp_non_zero_eri = nzeri;
    return 0;
}

/** (dp,ds)タイプの積分計算、および、G行列計算を行う関数
 *
 * @ingroup integ-twoint-direct
 * */
int ofmo_twoint_direct_dpds__(
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
	const int *pnao, const double Ds[], double G[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int last_ijcs=*plast_ijcs, last_klcs=*plast_klcs, nao=*pnao;
    long max_nzeri=*petmp_max_nzeri;
    long nzeri4, nzeri=*petmp_non_zero_eri;
    int Lab, Lcd, i, j, k, ix;
    int ijcs,        ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);
    float eps_sch = ofmo_twoint_eps_sch(0);
    double DPDS[6*3*6];

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];
    if ( last_ijcs != -1 ) {
	ijcs = last_ijcs;
	klcs = last_klcs+1;
    } else {
	ijcs = leading_cs_pair[Lab] + workerid;
	klcs = klcs0;
    }
    max_nzeri-=6*3*6;
    if ( nzeri>= max_nzeri ) {
	ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4, Ds, G );
	nzeri = 0;
    }
    nzeri4 = nzeri<<2;

    for (  ; ijcs<ijcs1; ijcs+=nworkers ) {
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
	
	for ( ; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < eps_ps4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
            float dmax = ofmo_twoint_dmax6(ics,jcs,kcs,lcs);
            if ( dmax*val_ab*val_cd < eps_sch ) continue;
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
		    for (k=0, kao=kao0; k<6; k++, kao++, ix++ ) {
			if ( fabs(DPDS[ix]) > eps_eri ) {
			    double x;
			    x  = DPDS[ix];
			    etmp_val[nzeri] = x;
			    etmp_ind4[nzeri4+0] = iao;
			    etmp_ind4[nzeri4+1] = jao;
			    etmp_ind4[nzeri4+2] = kao;
			    etmp_ind4[nzeri4+3] = lao;
			    nzeri++;
			    nzeri4+=4;
			}
		    }
		}
	    }	// for ( i, iao);
	    if ( nzeri >= max_nzeri ) {
		ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4,
			Ds, G );
		nzeri = nzeri4= 0;
	    }
	}	// for ( klcs );
	klcs = klcs0;
    }		// for ( ijcs );
    *petmp_non_zero_eri = nzeri;
    return 0;
}

/** (dp,dp)タイプの積分計算、および、G行列計算を行う関数
 *
 * @ingroup integ-twoint-direct
 * */
int ofmo_twoint_direct_dpdp__(
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
	const int *pnao, const double Ds[], double G[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int last_ijcs=*plast_ijcs, last_klcs=*plast_klcs, nao=*pnao;
    long max_nzeri=*petmp_max_nzeri;
    long nzeri4, nzeri=*petmp_non_zero_eri;
    int Lab, Lcd, i, j, k, l, ix, ipat;
    int IJ, KL;
    int ijcs,        ijcs1;
    int klcs, klcs0;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);
    float eps_sch = ofmo_twoint_eps_sch(0);
    double DPDP[6*3*6*3];

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    if ( last_ijcs != -1 ) {
	ijcs = last_ijcs;
	klcs = last_klcs+1;
    } else {
	ijcs = leading_cs_pair[Lab] + workerid;
	klcs = klcs0;
    }
    max_nzeri-=6*3*6*3;
    if ( nzeri>= max_nzeri ) {
	ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4, Ds, G );
	nzeri = 0;
    }
    nzeri4 = nzeri<<2;

    for (  ; ijcs<ijcs1; ijcs+=nworkers ) {
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
	
	for ( ; klcs<=ijcs; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < eps_ps4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
            float dmax = ofmo_twoint_dmax6(ics,jcs,kcs,lcs);
            if ( dmax*val_ab*val_cd < eps_sch ) continue;
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
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		for ( j=0, jao=jao0; j<3; j++, jao++ ) {
		    IJ = ((iao*iao+iao)>>1) + jao;
		    for ( k=0, kao=kao0; k<6; k++, kao++ ) {
			for ( l=0, lao=lao0; l<3; l++, lao++, ix++ ) {
			    KL = ((kao*kao+kao)>>1) + lao;
			    if ( fabs(DPDP[ix]) > eps_eri ) {
				if ( IJ>=KL || ipat ) {
				    double x;
				    coe = ( IJ==KL ? HALF : ONE );
				    x  = coe * DPDP[ix];
				    etmp_val[nzeri] = x;
				    etmp_ind4[nzeri4+0] = iao;
				    etmp_ind4[nzeri4+1] = jao;
				    etmp_ind4[nzeri4+2] = kao;
				    etmp_ind4[nzeri4+3] = lao;
				    nzeri++;
				    nzeri4+=4;
				}
			    }
			}
		    }	// for ( kao )
		}
	    }		// for ( iao )
	    if ( nzeri >= max_nzeri ) {
		ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4,
			Ds, G );
		nzeri = nzeri4= 0;
	    }
	}	// for ( klcs );
	klcs = klcs0;
    }		// for ( ijcs );
    *petmp_non_zero_eri = nzeri;
    return 0;
}

/** (dd,ss)タイプの積分計算、および、G行列計算を行う関数
 *
 * @ingroup integ-twoint-direct
 * */
int ofmo_twoint_direct_ddss__(
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
	const int *pnao, const double Ds[], double G[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int last_ijcs=*plast_ijcs, last_klcs=*plast_klcs, nao=*pnao;
    long max_nzeri=*petmp_max_nzeri;
    long nzeri4, nzeri=*petmp_non_zero_eri;
    int Lab, Lcd, i, j, ix;
    int ijcs,        ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe, coe0;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);
    float eps_sch = ofmo_twoint_eps_sch(0);
    double DDSS[6*6];

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];
    if ( last_ijcs != -1 ) {
	ijcs = last_ijcs;
	klcs = last_klcs+1;
    } else {
	ijcs = leading_cs_pair[Lab] + workerid;
	klcs = klcs0;
    }
    max_nzeri-=6*6;
    if ( nzeri>= max_nzeri ) {
	ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4, Ds, G );
	nzeri = 0;
    }
    nzeri4 = nzeri<<2;

    for (  ; ijcs<ijcs1; ijcs+=nworkers ) {
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
	
	for ( ; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < eps_ps4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
            float dmax = ofmo_twoint_dmax6(ics,jcs,kcs,lcs);
            if ( dmax*val_ab*val_cd < eps_sch ) continue;
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
	    coe0 = ( kao == lao ? HALF : ONE );
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		for ( j=0, jao=jao0; j<6; j++, jao++, ix++ ) {
		    if ( jao > iao ) continue;
		    if ( fabs(DDSS[ix]) > eps_eri ) {
			double x;
			coe = coe0;
			if ( iao == jao ) coe *= HALF;
			x  = coe * DDSS[ix];
			etmp_val[nzeri] = x;
			etmp_ind4[nzeri4+0] = iao;
			etmp_ind4[nzeri4+1] = jao;
			etmp_ind4[nzeri4+2] = kao;
			etmp_ind4[nzeri4+3] = lao;
			nzeri++;
			nzeri4+=4;
		    }
		} // for ( jao );
	    }	// for ( iao );
	    if ( nzeri >= max_nzeri ) {
		ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4,
			Ds, G );
		nzeri = nzeri4= 0;
	    }
	}	// for ( klcs );
	klcs = klcs0;
    }		// for ( ijcs );
    *petmp_non_zero_eri = nzeri;
    return 0;
}

/** (dd,ps)タイプの積分計算、および、G行列計算を行う関数
 *
 * @ingroup integ-twoint-direct
 * */
int ofmo_twoint_direct_ddps__(
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
	const int *pnao, const double Ds[], double G[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int last_ijcs=*plast_ijcs, last_klcs=*plast_klcs, nao=*pnao;
    long max_nzeri=*petmp_max_nzeri;
    long nzeri4, nzeri=*petmp_non_zero_eri;
    int Lab, Lcd, i, j, k, ix;
    int ijcs,        ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);
    float eps_sch = ofmo_twoint_eps_sch(0);
    double DDPS[6*6*3];

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];
    if ( last_ijcs != -1 ) {
	ijcs = last_ijcs;
	klcs = last_klcs+1;
    } else {
	ijcs = leading_cs_pair[Lab] + workerid;
	klcs = klcs0;
    }
    max_nzeri-=6*6*3;
    if ( nzeri>= max_nzeri ) {
	ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4, Ds, G );
	nzeri = 0;
    }
    nzeri4 = nzeri<<2;

    for (  ; ijcs<ijcs1; ijcs+=nworkers ) {
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
	
	for ( ; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < eps_ps4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
            float dmax = ofmo_twoint_dmax6(ics,jcs,kcs,lcs);
            if ( dmax*val_ab*val_cd < eps_sch ) continue;
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
			if ( fabs(DDPS[ix]) > eps_eri ) {
			    double x;
			    x  = coe * DDPS[ix];
			    etmp_val[nzeri] = x;
			    etmp_ind4[nzeri4+0] = iao;
			    etmp_ind4[nzeri4+1] = jao;
			    etmp_ind4[nzeri4+2] = kao;
			    etmp_ind4[nzeri4+3] = lao;
			    nzeri++;
			    nzeri4+=4;
			}	// if ( fabs );
		    }	// for ( kao )
		}	// for ( jao )
	    }	// for ( iao )
	    if ( nzeri >= max_nzeri ) {
		ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4,
			Ds, G );
		nzeri = nzeri4= 0;
	    }
	}	// for ( klcs );
	klcs = klcs0;
    }		// for ( ijcs );
    *petmp_non_zero_eri = nzeri;
    return 0;
}

/** (dd,pp)タイプの積分計算、および、G行列計算を行う関数
 *
 * @ingroup integ-twoint-direct
 * */
int ofmo_twoint_direct_ddpp__(
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
	const int *pnao, const double Ds[], double G[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int last_ijcs=*plast_ijcs, last_klcs=*plast_klcs, nao=*pnao;
    long max_nzeri=*petmp_max_nzeri;
    long nzeri4, nzeri=*petmp_non_zero_eri;
    int Lab, Lcd, i, j, k, l, ix;
    int ijcs,        ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe, coe0;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);
    float eps_sch = ofmo_twoint_eps_sch(0);
    double DDPP[6*6*3*3];

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];
    if ( last_ijcs != -1 ) {
	ijcs = last_ijcs;
	klcs = last_klcs+1;
    } else {
	ijcs = leading_cs_pair[Lab] + workerid;
	klcs = klcs0;
    }
    max_nzeri-=6*6*3*3;
    if ( nzeri>= max_nzeri ) {
	ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4, Ds, G );
	nzeri = 0;
    }
    nzeri4 = nzeri<<2;

    for (  ; ijcs<ijcs1; ijcs+=nworkers ) {
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
	
	for ( ; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < eps_ps4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
            float dmax = ofmo_twoint_dmax6(ics,jcs,kcs,lcs);
            if ( dmax*val_ab*val_cd < eps_sch ) continue;
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
		    if ( jao>iao ) { ix+=3*3; continue; }
		    coe0 = (iao==jao ? HALF : ONE );
		    for (k=0, kao=kao0; k<3; k++, kao++ ) {
			for ( l=0, lao=lao0; l<3; l++, lao++, ix++ ) {
			    if ( lao > kao ) continue;
			    coe = coe0 * ( kao==lao ? HALF : ONE );
			    if ( fabs(DDPP[ix]) > eps_eri ) {
				double x;
				x  = coe * DDPP[ix];
				etmp_val[nzeri] = x;
				etmp_ind4[nzeri4+0] = iao;
				etmp_ind4[nzeri4+1] = jao;
				etmp_ind4[nzeri4+2] = kao;
				etmp_ind4[nzeri4+3] = lao;
				nzeri++;
				nzeri4+=4;
			    }
			}
		    }
		}
	    }	// for ( i, iao);
	    if ( nzeri >= max_nzeri ) {
		ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4,
			Ds, G );
		nzeri = nzeri4= 0;
	    }
	}	// for ( klcs );
	klcs = klcs0;
    }		// for ( ijcs );
    *petmp_non_zero_eri = nzeri;
    return 0;
}

/** (dd,ds)タイプの積分計算、および、G行列計算を行う関数
 *
 * @ingroup integ-twoint-direct
 * */
int ofmo_twoint_direct_ddds__(
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
	const int *pnao, const double Ds[], double G[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int last_ijcs=*plast_ijcs, last_klcs=*plast_klcs, nao=*pnao;
    long max_nzeri=*petmp_max_nzeri;
    long nzeri4, nzeri=*petmp_non_zero_eri;
    int Lab, Lcd, i, j, k, ix;
    int ijcs,        ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);
    float eps_sch = ofmo_twoint_eps_sch(0);
    double DDDS[6*6*6];

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];
    if ( last_ijcs != -1 ) {
	ijcs = last_ijcs;
	klcs = last_klcs+1;
    } else {
	ijcs = leading_cs_pair[Lab] + workerid;
	klcs = klcs0;
    }
    max_nzeri-=6*6*6;
    if ( nzeri>= max_nzeri ) {
	ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4, Ds, G );
	nzeri = 0;
    }
    nzeri4 = nzeri<<2;

    for (  ; ijcs<ijcs1; ijcs+=nworkers ) {
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
	
	for ( ; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < eps_ps4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
            float dmax = ofmo_twoint_dmax6(ics,jcs,kcs,lcs);
            if ( dmax*val_ab*val_cd < eps_sch ) continue;
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
			if ( fabs(DDDS[ix]) > eps_eri ) {
			    double x;
			    x  = coe * DDDS[ix];
			    etmp_val[nzeri] = x;
			    etmp_ind4[nzeri4+0] = iao;
			    etmp_ind4[nzeri4+1] = jao;
			    etmp_ind4[nzeri4+2] = kao;
			    etmp_ind4[nzeri4+3] = lao;
			    nzeri++;
			    nzeri4+=4;
			}	// if ( fabs );
		    }	// for ( kao )
		}	// for ( jao )
	    }	// for ( iao )
	    if ( nzeri >= max_nzeri ) {
		ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4,
			Ds, G );
		nzeri = nzeri4= 0;
	    }
	}	// for ( klcs );
	klcs = klcs0;
    }		// for ( ijcs );
    *petmp_non_zero_eri = nzeri;
    return 0;
}

/** (dd,dp)タイプの積分計算、および、G行列計算を行う関数
 *
 * @ingroup integ-twoint-direct
 * */
int ofmo_twoint_direct_dddp__(
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
	const int *pnao, const double Ds[], double G[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int last_ijcs=*plast_ijcs, last_klcs=*plast_klcs, nao=*pnao;
    long max_nzeri=*petmp_max_nzeri;
    long nzeri4, nzeri=*petmp_non_zero_eri;
    int Lab, Lcd, i, j, k, l, ix;
    int ijcs,        ijcs1;
    int klcs, klcs0, klcs1;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);
    float eps_sch = ofmo_twoint_eps_sch(0);
    double DDDP[6*6*6*3];

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];
    if ( last_ijcs != -1 ) {
	ijcs = last_ijcs;
	klcs = last_klcs+1;
    } else {
	ijcs = leading_cs_pair[Lab] + workerid;
	klcs = klcs0;
    }
    max_nzeri-=6*6*6*3;
    if ( nzeri>= max_nzeri ) {
	ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4, Ds, G );
	nzeri = 0;
    }
    nzeri4 = nzeri<<2;

    for (  ; ijcs<ijcs1; ijcs+=nworkers ) {
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
	
	for ( ; klcs<klcs1; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < eps_ps4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
            float dmax = ofmo_twoint_dmax6(ics,jcs,kcs,lcs);
            if ( dmax*val_ab*val_cd < eps_sch ) continue;
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
			    if ( fabs(DDDP[ix]) > eps_eri ) {
				double x;
				x  = coe * DDDP[ix];
				etmp_val[nzeri] = x;
				etmp_ind4[nzeri4+0] = iao;
				etmp_ind4[nzeri4+1] = jao;
				etmp_ind4[nzeri4+2] = kao;
				etmp_ind4[nzeri4+3] = lao;
				nzeri++;
				nzeri4+=4;
			    }	// if ( fabs );
			}
		    }	// for ( kao )
		}	// for ( jao )
	    }	// for ( iao )
	    if ( nzeri>= max_nzeri ) {
		ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4, Ds, G );
		nzeri = nzeri4 = 0;
	    }
	}	// for ( klcs );
	klcs = klcs0;
    }		// for ( ijcs );
        *petmp_non_zero_eri = nzeri;
    return 0;
}

/** (dd,dd)タイプの積分計算、および、G行列計算を行う関数
 *
 * @ingroup integ-twoint-direct
 * */
int ofmo_twoint_direct_dddd__(
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
	const int *pnao, const double Ds[], double G[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int last_ijcs=*plast_ijcs, last_klcs=*plast_klcs, nao=*pnao;
    long max_nzeri=*petmp_max_nzeri;
    long nzeri4, nzeri=*petmp_non_zero_eri;
    int Lab, Lcd, i, j, k, l, ix, ipat;
    int I2, IJ, K2, KL;
    int ijcs,        ijcs1;
    int klcs, klcs0;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe, coe0;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);
    float eps_sch = ofmo_twoint_eps_sch(0);
    double DDDD[6*6*6*6];

    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    if ( last_ijcs != -1 ) {
	ijcs = last_ijcs;
	klcs = last_klcs+1;
    } else {
	ijcs = leading_cs_pair[Lab] + workerid;
	klcs = klcs0;
    }
    max_nzeri-=6*6*6*6;
    if ( nzeri>= max_nzeri ) {
	ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4, Ds, G );
	nzeri = 0;
    }
    nzeri4 = nzeri<<2;

    for (  ; ijcs<ijcs1; ijcs+=nworkers ) {
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
	
	for ( ; klcs<=ijcs; klcs++ ) {
	    val_cd = csp_schwarz[klcs];
	    if ( val_ab*val_cd < eps_ps4 ) continue;
	    kcs    = csp_ics[klcs];
	    lcs    = csp_jcs[klcs];
            float dmax = ofmo_twoint_dmax6(ics,jcs,kcs,lcs);
            if ( dmax*val_ab*val_cd < eps_sch ) continue;
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
	    for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
		I2 = (iao*iao+iao)>>1;
		for ( j=0, jao=jao0; j<6; j++, jao++ ) {
		    if ( jao>iao ) { ix+=6*6; continue; }
		    IJ = I2 + jao;
		    coe0 = ( iao==jao ? HALF : ONE );
		    for ( k=0, kao=kao0; k<6; k++, kao++ ) {
			K2 = (kao*kao+kao)>>1;
			for ( l=0, lao=lao0; l<6; l++, lao++, ix++ ) {
			    if ( lao>kao ) continue;
			    if ( fabs(DDDD[ix]) > eps_eri ) {
				double x;
				KL = K2 + lao;
				if ( IJ >= KL || ipat ) {
				    coe = coe0;
				    if ( kao==lao ) coe *= HALF;
				    if ( KL == IJ ) coe *= HALF;
				    x  = coe * DDDD[ix];
				    etmp_val[nzeri] = x;
				    etmp_ind4[nzeri4+0] = iao;
				    etmp_ind4[nzeri4+1] = jao;
				    etmp_ind4[nzeri4+2] = kao;
				    etmp_ind4[nzeri4+3] = lao;
				    nzeri++;
				    nzeri4+=4;
				}
			    } // if ( fabs )
			} // for ( lao )
		    } // for ( kao )
		} // for ( jao )
	    }	// for ( iao )
	    if ( nzeri >= max_nzeri ) {
		ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4,
			Ds, G );
		nzeri = nzeri4= 0;
	    }
	}	// for ( klcs );
	klcs = klcs0;
    }		// for ( ijcs );
    *petmp_non_zero_eri = nzeri;
    return 0;
}
