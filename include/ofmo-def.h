#ifndef _OFMO_DEF_H_
#define _OFMO_DEF_H_

#ifndef false
#define false 0
#endif

#ifndef true
#define true 1
#endif

#ifndef MAXLQN
#define MAXLQN 12
#endif

#ifndef dbg
#define dbg(...) \
    ( printf("%s %u @%s:",__FILE__, __LINE__, __func__), \
      printf(" "__VA_ARGS__), fflush(stdout) )
#endif

// added in 2011/10/12
#ifndef fdbg
#define fdbg(fp, ...) \
    ( fprintf(fp, "%s %u @%s:",__FILE__, __LINE__, __func__), \
      fprintf(fp, " "__VA_ARGS__), fflush(fp) )
#endif

//#define Free(p) if ((p)!=NULL) free( (p) )
#define Free(a) do { if ( a != NULL ) { free( a ); a = NULL;} } while (0)

#ifndef MAXSTRLEN
#define MAXSTRLEN 256
#endif

#ifndef MAXTOKEN
#define MAXTOKEN 64
#endif

#ifndef MAXATOMICNUMBER
#define MAXATOMICNUMBER 103
#endif

#define BOHR_RADIUS 0.52917724924e0

// データ塊数の最大値
#ifndef OFMO_MAX_DATA
#define OFMO_MAX_DATA 16
#endif

// マスタからの依頼のID
#define OFMO_FINALIZE		0
#define OFMO_INIT_DENS		1
#define OFMO_SCF		2
#define OFMO_APPROX		3
#define OFMO_DISTANCE		4
#define OFMO_BARRIER		5
#define OFMO_ACCEPT		6
#define OFMO_REJECT		7
#define OFMO_UPDATE_DATA	8
#define OFMO_CHANGE_MSERV	9
#define OFMO_DISCONNECT		10
#define OFMO_RESET_MON_DATA	11
#define OFMO_GET_POP_DATA	12

// 計算手法
#define OFMO_UNDEF	0
#define OFMO_RHF	1
#define OFMO_RIMP2	2
#define OFMO_DFT	3

// ワーカからマスタに転送するpopulationデータの種類
#define	OFMO_BOTH_POPS	0	// AO pop. Atomic pop. の両方
#define OFMO_AOPOP_ONLY	1	// AO pop. のみ
#define OFMO_ATPOP_ONLY	2	// atomic pop. のみ

// ES dimerの計算や4中心項の数
// 注意点：MAXNJOB <= MAXNIFC4C/2
#define MAXNIFC4C	64
#define MAXNJOB		16
#define NLINE		4

// マスター、ワーカー、メモリサーバー間の信号伝達の要素数

#define OFMO_I_CMD	0	// 計算の種類、処理内容
#define OFMO_I_METHOD	1	// 計算手法
#define OFMO_I_SCC	2	// 現在のSCC繰り返し回数
#define OFMO_I_CONV	3	// 収束条件
#define OFMO_I_NMON	5	// 関係するモノマー数
				// ES dimer計算時は、計算する近似ダイマー数
#define OFMO_I_MON1	6	// モノマー１
#define OFMO_I_MON2	7	// モノマー２
#define OFMO_I_MON3	8	// モノマー３

#define OFMO_IMSG_SZ	(MAXNJOB*2+OFMO_I_MON1)
#define OFMO_DMSG_SZ	6

#define OFMO_D_ENERGY	0
#define OFMO_D_ENERGY0	1
#define OFMO_D_DDV	2

#define OFMO_TAG_SIG	11
#define OFMO_TAG_CMD	12
#define OFMO_TAG_RET	13
#define OFMO_TAG_END	14
#define OFMO_TAG_DAT	15

// シグナル用変数のおおきさ
#define NBUF	10

// メモリサーバーへのシグナル
#define OFMO_PUT	10
#define OFMO_GET	11
#define OFMO_ACC	12
#define OFMO_ZCR	13

// メモリサーバーが保持しているデータのID
#define OFMO_DENS1	0
#define OFMO_DENS2	1
#define OFMO_AOPOP1	2
#define OFMO_AOPOP2	3
#define OFMO_ATPOP1	4
#define OFMO_ATPOP2	5
#define OFMO_DISTA	6
#define OFMO_ENERGY	7
#define OFMO_ENERGY0	8
#define OFMO_TOTAL_AOPOP	9
#define OFMO_TOTAL_ATPOP	10

// 環境ポテンシャルの計算レベル
#define OFMO_IFC0C	0
#define OFMO_IFC4C	1
#define OFMO_IFC3C	2
#define OFMO_IFC2C	3


#if 0
// defined but not used
static double check_sum( int n, double d[] ) {
    double sum=0.e0;
    for ( int i=0; i<n; i++ ) sum += d[i];
    return sum;
}
#endif

#endif
