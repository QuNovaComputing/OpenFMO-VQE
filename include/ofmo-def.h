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

// Maximum number of data chunks
#ifndef OFMO_MAX_DATA
#define OFMO_MAX_DATA 16
#endif

// ID of the request from the master
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

// Calculation method
#define OFMO_UNDEF	0
#define OFMO_RHF	1
#define OFMO_RIMP2	2
#define OFMO_DFT	3
#define OFMO_VQE  4
#define OFMO_RHF_VQE  5
#define OFMO_VQE_RHF  6

// Types of population data to transfer from the worker to the master
#define	OFMO_BOTH_POPS	0	// AO pop. Atomic pop. Both
#define OFMO_AOPOP_ONLY	1	// AO pop. Only
#define OFMO_ATPOP_ONLY	2	// atomic pop. Only

// ES dimer calculation and number of 4 central terms
// Critical Point：MAXNJOB <= MAXNIFC4C/2
#define MAXNIFC4C	64
#define MAXNJOB		16
#define NLINE		4

// Number of elements of signal transmission between master, worker and memory server

#define OFMO_I_CMD	0	// Calculation type, processing content
#define OFMO_I_METHOD	1	// Calculation method
#define OFMO_I_SCC	2	// Current SCC repeat count
#define OFMO_I_CONV	3	// Convergence condition
#define OFMO_I_NMON	5	// Number of monomers involved
				// When calculating ES dimer, the approximate number of dimers to be calculated
#define OFMO_I_MON1	6	// モノマー１
#define OFMO_I_MON2	7	// Monomer 2
#define OFMO_I_MON3	8	// Monomer 3

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
