/**
 * @file ofmo-projection.c
 *
 * 射影演算子を作成する関数を定義したファイル
 * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ofmo-def.h"
#include "ofmo-data.h"
#include "ofmo-prof.h"

#include "ofmo-mat.h"
#include "ofmo-basis.h"

static int NLMO  = 0;
static int NFUNC = 0;
//static int NAO   = 0;
static int *tuao2fsao = NULL;
static double *C0 = NULL;	/* 回転前の混成軌道係数 */
static double *CR = NULL;	/* 回転後の混成軌道係数 */
static double *CW = NULL;	/* 回転後の混成軌道係数（フル）*/
static double *CS = NULL;  /* ソート済みの回転後の混成軌道係数（フル）*/
static double *CTS = NULL; /* C'Sの代入先 */
static double *S = NULL; /* 正方行列形式のソート済み重なり行列 */
static double *P = NULL; /* 正方行列形式のソート済み射影演算子 */

static void dealloc() {
    Free( C0 );
    Free( CR );
    Free( CW );
    Free( CS );
    Free( CTS );
    Free( S );
    Free( P );
    Free( tuao2fsao );
    C0 = NULL;
    CR = NULL;
    CW = NULL;
    CS = NULL;
    CTS = NULL;
    S = NULL;
    P = NULL;
    tuao2fsao = NULL;
    NLMO = 0;
    NFUNC = 0;
    //NAO = 0;
}

static int alloc( const int nlmo, const int nfunc ) {
    static int called = false;
    if ( nlmo <= NLMO && nfunc <= NFUNC ) return 0;
    dealloc();
    int nn, maxnfao, nbody, nao, nao_total, ierr;
    ierr = ofmo_data_get_vals("maxnfao nbody nao",
	    &maxnfao, &nbody, &nao_total );
    if ( ierr != 0 ) {
	if ( fp_prof ) fdbg(fp_prof, "ERROR\n");
	return -1;
    }
    nao = nbody * maxnfao;
    nn = nfunc * nlmo;
    C0 = (double*)malloc( sizeof(double) * nn );
    CR = (double*)malloc( sizeof(double) * nn );
    CW = (double*)malloc( sizeof(double) * nlmo * nao );
    CS = (double*)malloc( sizeof(double) * nlmo * nao );
    CTS = (double*)malloc( sizeof(double) * nlmo * nao );
    S = (double*)malloc( sizeof(double) * nao * nao );
    P = (double*)malloc( sizeof(double) * nao * nao );
    tuao2fsao = (int*)malloc( sizeof(int) * nao_total );
    NLMO  = nlmo;
    NFUNC = nfunc;
    //NAO   = nao;
    if ( !called ) {
	atexit( dealloc );
	called = true;
    }
    return 0;
}

/** Projection operator A function that creates a term
 *
 * A function that calculates the projection operator term applied to the part
 * where the bond of the monomer is broken.
 *
 * @param[in] nmonomer フラグメントを構成するモノマー数
 * @param[in] monomer_list[] フラグメントを構成するモノマー番号のリスト
 * @param[in] nao フラグメントのAO数
 * @param[in] fsao2tuao[] ソートされたフラグメントAOの全体の非ソートAO
 * 番号
 * @param[in] Ss[] ソートされた重なり行列（圧縮U形式）
 * @param[out] Ps[] Sorted projection operator terms (compressed U format)
 *
 * @ingroup ofmo-calc
 * */
int ofmo_projection_operator(
	const int nmonomer, const int monomer_list[],
	const int nao, const int fsao2tuao[], const double Ss[],
	double Ps[]) {
    /* local variables */
    int imon, ifrag, ibd, i, ibda, flag, kmon, tfrag1, tfrag2, jbda;
    double rx, ry, rz, rxy2, rxy, r2, r;
    double cph, sph, cth, sth;
    double *rot, rmu[3*3];
    double *cr, *c0, *cw, *c, *s, *cts, bk=1.e6;
    int ics0, ics1, ics, lmo, ix, lqn, iao, iao0, jao;
    int plmo, k;
    double t;
    /* 外部で定義されているはずの変数 */
    static int called = false;
    static double rsfact2[MAXLQN+1], sfact2[MAXLQN+1];
    static double *CLMO;
    static int nlmo, nfunc, *flg_z, nao_total;
    static int *nfbond, **fbondsn1/*, **fbda*/;
    static int *at2frg, *welec, *woelec;
    static double *atom_x, *atom_y, *atom_z;
    static int *atm_lcs, *ushel_ini, *ushel_lqn;
    /* 初期設定（最初の1回だけ）*/
    if ( !called ) {
	int ifact2[MAXLQN+1], ierr;
	ifact2[0] = 1;
	for ( lqn=1; lqn<=MAXLQN; lqn++ )
	    ifact2[lqn] = ifact2[lqn-1]*(2*lqn-1);
	for ( lqn=0; lqn<=MAXLQN; lqn++ ) {
	    sfact2[lqn]  = sqrt( (double)ifact2[lqn] );
	    rsfact2[lqn] = 1.e0 / sfact2[lqn];
	}
	ierr = ofmo_data_get_vals(
		"nlmo nfunc clmo ban nfbond fbsn1 woelec welec "
		"at2frg atx aty atz atm_lcs ushel_ini ushel_lqn nao",
		&nlmo, &nfunc, &CLMO, &flg_z, &nfbond, &fbondsn1,
		&woelec, &welec, &at2frg, &atom_x, &atom_y, &atom_z,
		&atm_lcs, &ushel_ini, &ushel_lqn, &nao_total );
	if ( ierr != 0 ) return -1;
	if ( alloc( nlmo, nfunc ) != 0 ) return -1;
	called = true;
    }
    for ( int iao=0; iao<nao_total; iao++ ) tuao2fsao[iao] = -1;
    for ( int iao=0; iao<nao; iao++ ) tuao2fsao[ fsao2tuao[iao] ] = iao;
    ofmo_unpack_matrix( nao, Ss, S );
    memset( P, '\0', sizeof(double)*nao*nao );

    for ( imon=0; imon<nmonomer; imon++ ) {
	ifrag = monomer_list[imon];
	for ( ibd=0; ibd<nfbond[ifrag]; ibd++ ) {
	    i = fbondsn1[ifrag][ibd];
	    i = ( i>0 ? i : -i );
	    ibda = woelec[i-1];	/* 電子を与える結合原子の全体でのSN */
	    jbda = welec[i-1];	/* 電子をもらう結合原子の全体でのSN */
	    tfrag1 = at2frg[ibda];
	    tfrag2 = at2frg[jbda];
	    /* 電子を与える結合原子がフラグメントを構成する
	     * 他のモノマーの原子の場合には、何もしない */
	    flag = false;
	    for ( kmon=0; kmon<nmonomer; kmon++ ) {
		if ( kmon == imon ) continue;
		if ( tfrag1 == monomer_list[kmon] ) flag = true;
		if ( tfrag2 == monomer_list[kmon] ) flag = true;
	    }
	    if ( flag ) continue;
	    /* z軸を結合方向に回転させるための角度に対応する余弦と
	     * 正弦の値を算出 */
	    rx = atom_x[jbda] - atom_x[ibda];
	    ry = atom_y[jbda] - atom_y[ibda];
	    rz = atom_z[jbda] - atom_z[ibda];
	    rxy2 = rx*rx + ry*ry;
	    if ( rxy2 < 1.e-8 ) {
		cph = 1.e0;
		sph = 0.e0;
		rxy = 0.e0;
		cth = (rz > 0.e0 ? 1.e0 : -1.e0 );
		sth = 0.e0;
	    } else {
		rxy = sqrt( rxy2 );
		r2  = rxy2 + rz*rz;
		r   = sqrt( r2 );
		cph = rx / rxy;
		sph = ry / rxy;
		cth = rz / r;
		sth = rxy / r;
	    }
	    /* z軸を結合方向に回転させる行列の要素 */
	    /*rmt[0+0*3] = cph*cth; rmt[0+1*3] = -sph; rmt[0+2*3] = cph*sth;
	    rmt[1+0*3] = sph*cth; rmt[1+1*3] =  cph; rmt[1+2*3] = sph*sth;
	    rmt[2+0*3] = -sth;    rmt[2+1*3] = 0.e0; rmt[2+2*3] = cth;*/
	    //
            rmu[0]=cph*cph*cth+sph*sph; rmu[3]=    cph*sph*(cth-1.e0);
	    rmu[6]=cph*sth;
            rmu[1]=    cph*sph*(cth-1.e0); rmu[4]=sph*sph*cth+cph*cph;
	    rmu[7]=sph*sth;
            rmu[2]=           -cph*sth; rmu[5]=           -sph*sth;
	    rmu[8]=    cth;
	    rot = rmu;

	    /* 混成軌道の係数行列をコピー */
	    for ( i=0; i<nlmo*nfunc; i++ ) C0[i] = CLMO[i];
	    /* */
	    ics0 = atm_lcs[ibda];
	    ics1 = atm_lcs[ibda+1];
	    for ( lmo=0, cr=CR, c0=C0; lmo<nlmo;
		    lmo++, cr+=nfunc, c0+=nfunc ) {
		ix = 0;
		for ( ics=ics0; ics<ics1; ics++ ) {
		    lqn = ushel_lqn[ics];
		    if ( lqn == 0 ) {	/* s-type */
			cr[ix+0] = c0[ix+0];
			ix++;
		    } else if ( lqn == 1 ) {	/* p-type */
			cr[ix+0] = rot[0+0*3]*c0[0+ix]
			         + rot[0+1*3]*c0[1+ix]
				 + rot[0+2*3]*c0[2+ix];
			cr[ix+1] = rot[1+0*3]*c0[0+ix]
			         + rot[1+1*3]*c0[1+ix]
				 + rot[1+2*3]*c0[2+ix];
			cr[ix+2] = rot[2+0*3]*c0[0+ix]
			         + rot[2+1*3]*c0[1+ix]
				 + rot[2+2*3]*c0[2+ix];
			ix+=3;
		    } else if ( lqn == 2 ) {	/* d-type */
			c0[ix+0]*=rsfact2[2];	/* 6dを仮定している */
			c0[ix+1]*=rsfact2[2];
			c0[ix+2]*=rsfact2[2];
			cr[ix+0] = rot[0+0*3]*rot[0+0*3]*c0[ix+0]
			      + rot[0+1*3]*rot[0+1*3]*c0[ix+1]
			      + rot[0+2*3]*rot[0+2*3]*c0[ix+2]
			      + rot[0+0*3]*rot[0+1*3]*c0[ix+3]
			      + rot[0+1*3]*rot[0+2*3]*c0[ix+4]
			      + rot[0+2*3]*rot[0+0*3]*c0[ix+5];

			cr[ix+1] = rot[1+0*3]*rot[1+0*3]*c0[ix+0]
			      + rot[1+1*3]*rot[1+1*3]*c0[ix+1]
			      + rot[1+2*3]*rot[1+2*3]*c0[ix+2]
			      + rot[1+0*3]*rot[1+1*3]*c0[ix+3]
			      + rot[1+1*3]*rot[1+2*3]*c0[ix+4]
			      + rot[1+2*3]*rot[1+0*3]*c0[ix+5];

			cr[ix+2] = rot[2+0*3]*rot[2+0*3]*c0[ix+0]
			      + rot[2+1*3]*rot[2+1*3]*c0[ix+1]
			      + rot[2+2*3]*rot[2+2*3]*c0[ix+2]
			      + rot[2+0*3]*rot[2+1*3]*c0[ix+3]
			      + rot[2+1*3]*rot[2+2*3]*c0[ix+4]
			      + rot[2+2*3]*rot[2+0*3]*c0[ix+5];

			cr[ix+3] = rot[0+0*3]*rot[1+0*3]*c0[ix+0]*2.e0
			      + rot[0+1*3]*rot[1+1*3]*c0[ix+1]*2.e0
			      + rot[0+2*3]*rot[1+2*3]*c0[ix+2]*2.e0
			      + (rot[0+0*3]*rot[1+1*3]+rot[0+3*1]*rot[1+3*0])*c0[ix+3]
			      + (rot[0+1*3]*rot[1+2*3]+rot[0+3*2]*rot[1+3*1])*c0[ix+4]
			      + (rot[0+2*3]*rot[1+0*3]+rot[0+3*0]*rot[1+3*2])*c0[ix+5];


			cr[ix+4] = rot[2+0*3]*rot[0+0*3]*c0[ix+0]*2.e0
			      + rot[2+1*3]*rot[0+1*3]*c0[ix+1]*2.e0
			      + rot[2+2*3]*rot[0+2*3]*c0[ix+2]*2.e0
			      + (rot[2+0*3]*rot[0+1*3]+rot[2+1*3]*rot[0+0*3])*c0[ix+3]
			      + (rot[2+1*3]*rot[0+2*3]+rot[2+2*3]*rot[0+1*3])*c0[ix+4]
			      + (rot[2+2*3]*rot[0+0*3]+rot[2+0*3]*rot[0+2*3])*c0[ix+5];

			cr[ix+5] = rot[1+0*3]*rot[2+0*3]*c0[ix+0]*2.e0
			      + rot[1+1*3]*rot[2+1*3]*c0[ix+1]*2.e0
			      + rot[1+2*3]*rot[2+2*3]*c0[ix+2]*2.e0
			      + (rot[1+0*3]*rot[2+1*3]+rot[1+1*3]*rot[2+0*3])*c0[ix+3]
			      + (rot[1+1*3]*rot[2+2*3]+rot[1+2*3]*rot[2+1*3])*c0[ix+4]
			      + (rot[1+2*3]*rot[2+0*3]+rot[1+0*3]*rot[2+2*3])*c0[ix+5];
			/*
			cr[ix+4] = rot[1+0*3]*rot[2+0*3]*c0[ix+0]*2.e0
			      + rot[1+1*3]*rot[2+1*3]*c0[ix+1]*2.e0
			      + rot[1+2*3]*rot[2+2*3]*c0[ix+2]*2.e0
			      + (rot[1+0*3]*rot[2+1*3]+rot[1+1*3]*rot[2+0*3])*c0[ix+3]
			      + (rot[1+1*3]*rot[2+2*3]+rot[1+2*3]*rot[2+1*3])*c0[ix+4]
			      + (rot[1+2*3]*rot[2+0*3]+rot[1+0*3]*rot[2+2*3])*c0[ix+5];

			cr[ix+5] = rot[2+0*3]*rot[0+0*3]*c0[ix+0]*2.e0
			      + rot[2+1*3]*rot[0+1*3]*c0[ix+1]*2.e0
			      + rot[2+2*3]*rot[0+2*3]*c0[ix+2]*2.e0
			      + (rot[2+0*3]*rot[0+1*3]+rot[2+1*3]*rot[0+0*3])*c0[ix+3]
			      + (rot[2+1*3]*rot[0+2*3]+rot[2+2*3]*rot[0+1*3])*c0[ix+4]
			      + (rot[2+2*3]*rot[0+0*3]+rot[2+0*3]*rot[0+2*3])*c0[ix+5];
			      */
			cr[ix+0] *= sfact2[2];
			cr[ix+1] *= sfact2[2];
			cr[ix+2] *= sfact2[2];
			c0[ix+0]*=sfact2[2];
			c0[ix+1]*=sfact2[2];
			c0[ix+2]*=sfact2[2];
			ix+=6;
		    }
		}	// for ( ics )
	    }	// for ( lmo )
	    /* 射影演算子用のMO係数行列（非ソート）に代入 */
	    int I0;
	    for ( i=0; i<nlmo*nao; i++ ) CW[i] = 0.e0;

	    for ( lmo=0, c=CR, cw=CW; lmo<nlmo;
		    lmo++,c+=nfunc, cw+=nao ) {
		iao=0;
		for ( ics=ics0; ics<ics1; ics++ ) {
		    iao0 = ushel_ini[ics];
		    if ( (I0=tuao2fsao[iao0]) < 0 ) {
			if( fp_prof ) fdbg(fp_prof, "error\n");
			return -1;
		    }
		    lqn  = ushel_lqn[ics];
		    if ( lqn == 0 ) {
			cw[I0+0] = c[iao];
			iao++;
		    } else if ( lqn == 1 ) {
			cw[I0+0] = c[iao+0];
			cw[I0+1] = c[iao+1];
			cw[I0+2] = c[iao+2];
			iao+=3;
		    } else if ( lqn == 2 ) {
			cw[I0+0] = c[iao+0];
			cw[I0+1] = c[iao+1];
			cw[I0+2] = c[iao+2];
			cw[I0+3] = c[iao+3];
			cw[I0+4] = c[iao+4];
			cw[I0+5] = c[iao+5];
			iao+=6;
		    }
		}	// for ( ics )
		// MO係数行列をソートする
	    }
	    cts = CTS;
	    /* P=Bk*(C'S)'*(C'S)を計算する */
	    /* 電子を与える結合原子が、このモノマーの原子の場合には、
	     * 結合性の混成軌道を除外する */
	    if ( at2frg[ibda] == ifrag ) {
		//plmo=1; c0=cs;
		plmo=1; c0=CW;
	    } else {
		//plmo=4; c0=cs+nfunc;
		plmo=4; c0=CW+nao;
	    }
	    /* C'S */
	    for ( lmo=0; lmo<plmo; lmo++, c0+=nao ) {
		for ( iao=0, s=S; iao<nao; iao++, s+=nao ) {
		    t = 0.e0;
		    for ( k=0; k<nao; k++ ) t += c0[k]*s[k];
		    cts[lmo*nao+iao] = t;
		}
	    }
	    /* P=bk*(C'S)'*(C'S) */
	    for ( iao=0; iao<nao; iao++ ) {
		for ( jao=0; jao<nao; jao++ ) {
		    t = 0.e0;
		    for ( k=0; k<plmo; k++ )
			t += cts[k*nao+iao]*cts[k*nao+jao];
		    P[iao*nao+jao] += t;
		}
	    }
	}	// for ( iat )
    }		// for ( imon )
    for ( i=0; i<nao*nao; i++ ) P[i] *= bk;
    ofmo_pack_matrix( nao, P, Ps );
    return 0;
}
