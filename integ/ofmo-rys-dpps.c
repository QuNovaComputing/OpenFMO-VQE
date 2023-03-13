#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "ofmo-twoint.h"

#ifndef false
#define false 0
#endif
#ifndef true
#define true 1
#endif

#define HALF    0.5e0
#define ONE     1.e0
#define ZERO    0.e0

#define EPS_PS4 1.e-30
#define EPS_ERI 1.e-15

#define OFMO_EBUF_FULL          1
#define OFMO_EBUF_NOFULL        0

#define EPS_PS_PAIR     1.e-32
#define EPS_CS_PAIR2    1.e-30
#define MAXNPSPAIR 100

extern double* ofmo_integ_getadd_xint( const int mythread );
extern double* ofmo_integ_getadd_yint( const int mythread );
extern double* ofmo_integ_getadd_zint( const int mythread );
extern double* ofmo_integ_getadd_eh( const int mythread );
extern double* ofmo_integ_getadd_eri( const int mythread );

extern void calc_root( const int nroot, const double T,
        double *U, double *W );
extern int ofmo_integ_add_fock( const int nao, const size_t nstored_eri,
        const double eri_val[], const short int eri_ind4[],
        const double D[], double G[] );

extern double* ofmo_getadd_dfact();
static double *DFACT;

static void ofmo_hrr_clear_dpps( double *eh ) {
    int i;
    // (DS|PS)
    for ( i=0; i<(0+18); i++ ) eh[i] = 0.e0;
    // (FS|PS)
    for ( i=18; i<(18+30); i++ ) eh[i] = 0.e0;
}

static void ofmo_hrr_coef_dpps(
        double *eh, double *DINT ) {
    int i, j, k, l, iao, jao, kao, lao, ix;
    double coef_a, coef_ab, coef_abc;
    double *th;
    th = &eh[48];
    ix = 0;
    for ( i=0, iao=4; i<6; i++, iao++ ) {
        coef_a = DFACT[iao];
        for ( j=0, jao=1; j<3; j++, jao++ ) {
            coef_ab = coef_a * DFACT[jao];
            for ( k=0, kao=1; k<3; k++, kao++ ) {
                coef_abc = coef_ab * DFACT[kao];
                for ( l=0, lao=0; l<1; l++, lao++ ) {
                    DINT[ix] = coef_abc * DFACT[lao] * th[ix];
                    ix++;
                }
            }
        }
    }
}

static void ofmo_hrr_calc_dpps( double *eh,
        const double BA[3], const double DC[3] ) {
    // HRR for (XX|XS)-type integral (center AB)
    // (DP,PS)
    eh[  48] = eh[  18] - BA[0]*eh[   0];
    eh[  49] = eh[  19] - BA[0]*eh[   1];
    eh[  50] = eh[  20] - BA[0]*eh[   2];
    eh[  51] = eh[  27] - BA[1]*eh[   0];
    eh[  52] = eh[  28] - BA[1]*eh[   1];
    eh[  53] = eh[  29] - BA[1]*eh[   2];
    eh[  54] = eh[  30] - BA[2]*eh[   0];
    eh[  55] = eh[  31] - BA[2]*eh[   1];
    eh[  56] = eh[  32] - BA[2]*eh[   2];
    eh[  57] = eh[  33] - BA[0]*eh[   3];
    eh[  58] = eh[  34] - BA[0]*eh[   4];
    eh[  59] = eh[  35] - BA[0]*eh[   5];
    eh[  60] = eh[  21] - BA[1]*eh[   3];
    eh[  61] = eh[  22] - BA[1]*eh[   4];
    eh[  62] = eh[  23] - BA[1]*eh[   5];
    eh[  63] = eh[  42] - BA[2]*eh[   3];
    eh[  64] = eh[  43] - BA[2]*eh[   4];
    eh[  65] = eh[  44] - BA[2]*eh[   5];
    eh[  66] = eh[  36] - BA[0]*eh[   6];
    eh[  67] = eh[  37] - BA[0]*eh[   7];
    eh[  68] = eh[  38] - BA[0]*eh[   8];
    eh[  69] = eh[  45] - BA[1]*eh[   6];
    eh[  70] = eh[  46] - BA[1]*eh[   7];
    eh[  71] = eh[  47] - BA[1]*eh[   8];
    eh[  72] = eh[  24] - BA[2]*eh[   6];
    eh[  73] = eh[  25] - BA[2]*eh[   7];
    eh[  74] = eh[  26] - BA[2]*eh[   8];
    eh[  75] = eh[  27] - BA[0]*eh[   9];
    eh[  76] = eh[  28] - BA[0]*eh[  10];
    eh[  77] = eh[  29] - BA[0]*eh[  11];
    eh[  78] = eh[  33] - BA[1]*eh[   9];
    eh[  79] = eh[  34] - BA[1]*eh[  10];
    eh[  80] = eh[  35] - BA[1]*eh[  11];
    eh[  81] = eh[  39] - BA[2]*eh[   9];
    eh[  82] = eh[  40] - BA[2]*eh[  10];
    eh[  83] = eh[  41] - BA[2]*eh[  11];
    eh[  84] = eh[  30] - BA[0]*eh[  12];
    eh[  85] = eh[  31] - BA[0]*eh[  13];
    eh[  86] = eh[  32] - BA[0]*eh[  14];
    eh[  87] = eh[  39] - BA[1]*eh[  12];
    eh[  88] = eh[  40] - BA[1]*eh[  13];
    eh[  89] = eh[  41] - BA[1]*eh[  14];
    eh[  90] = eh[  36] - BA[2]*eh[  12];
    eh[  91] = eh[  37] - BA[2]*eh[  13];
    eh[  92] = eh[  38] - BA[2]*eh[  14];
    eh[  93] = eh[  39] - BA[0]*eh[  15];
    eh[  94] = eh[  40] - BA[0]*eh[  16];
    eh[  95] = eh[  41] - BA[0]*eh[  17];
    eh[  96] = eh[  42] - BA[1]*eh[  15];
    eh[  97] = eh[  43] - BA[1]*eh[  16];
    eh[  98] = eh[  44] - BA[1]*eh[  17];
    eh[  99] = eh[  45] - BA[2]*eh[  15];
    eh[ 100] = eh[  46] - BA[2]*eh[  16];
    eh[ 101] = eh[  47] - BA[2]*eh[  17];
}

static void ofmo_xyzint_dpps(
        const double *F00, const double *B00, const double *B10,
        const double *B01, const double *C00, const double *CP00,
        double *xint, double *yint, double *zint ) {
    int Lab, Lcd;
    int m, m3, N, M, ix3, ix2, ix1, ix0, nroot;
    double C10[3], CP10[3], CP01[3], C01[3];
    // (0,0)
    xint[  0]=1.e0;
    yint[  0]=1.e0;
    zint[  0]=F00[0];
    xint[  1]=1.e0;
    yint[  1]=1.e0;
    zint[  1]=F00[1];
    xint[  2]=1.e0;
    yint[  2]=1.e0;
    zint[  2]=F00[2];
    // (1,0)
    xint[  6]=C00[ 0];
    yint[  6]=C00[ 1];
    zint[  6]=C00[ 2]*F00[0];
    xint[  7]=C00[ 3];
    yint[  7]=C00[ 4];
    zint[  7]=C00[ 5]*F00[1];
    xint[  8]=C00[ 6];
    yint[  8]=C00[ 7];
    zint[  8]=C00[ 8]*F00[2];
    // (0,1)
    xint[  3]=CP00[ 0];
    yint[  3]=CP00[ 1];
    zint[  3]=CP00[ 2]*F00[0];
    xint[  4]=CP00[ 3];
    yint[  4]=CP00[ 4];
    zint[  4]=CP00[ 5]*F00[1];
    xint[  5]=CP00[ 6];
    yint[  5]=CP00[ 7];
    zint[  5]=CP00[ 8]*F00[2];
    // (1,1)
    xint[  9]=CP00[ 0]*xint[  6]+B00[0];
    yint[  9]=CP00[ 1]*yint[  6]+B00[0];
    zint[  9]=CP00[ 2]*zint[  6]+B00[0]*F00[0];
    xint[ 10]=CP00[ 3]*xint[  7]+B00[1];
    yint[ 10]=CP00[ 4]*yint[  7]+B00[1];
    zint[ 10]=CP00[ 5]*zint[  7]+B00[1]*F00[1];
    xint[ 11]=CP00[ 6]*xint[  8]+B00[2];
    yint[ 11]=CP00[ 7]*yint[  8]+B00[2];
    zint[ 11]=CP00[ 8]*zint[  8]+B00[2]*F00[2];
    // (N,0) and (N,1)
    for ( m=0; m<3; m++ ) {
        C10[m]  = 0.e0;
        CP10[m] = B00[m];
    }
    // (2,0)
    C10[0] += B10[0];
    xint[ 12]=C00[ 0]*xint[  6]+C10[0]*xint[  0];
    yint[ 12]=C00[ 1]*yint[  6]+C10[0]*yint[  0];
    zint[ 12]=C00[ 2]*zint[  6]+C10[0]*zint[  0];
    C10[1] += B10[1];
    xint[ 13]=C00[ 3]*xint[  7]+C10[1]*xint[  1];
    yint[ 13]=C00[ 4]*yint[  7]+C10[1]*yint[  1];
    zint[ 13]=C00[ 5]*zint[  7]+C10[1]*zint[  1];
    C10[2] += B10[2];
    xint[ 14]=C00[ 6]*xint[  8]+C10[2]*xint[  2];
    yint[ 14]=C00[ 7]*yint[  8]+C10[2]*yint[  2];
    zint[ 14]=C00[ 8]*zint[  8]+C10[2]*zint[  2];
    // (2,1)
    CP10[0] += B00[0];
    xint[ 15]=CP00[ 0]*xint[ 12]+CP10[0]*xint[  6];
    yint[ 15]=CP00[ 1]*yint[ 12]+CP10[0]*yint[  6];
    zint[ 15]=CP00[ 2]*zint[ 12]+CP10[0]*zint[  6];
    CP10[1] += B00[1];
    xint[ 16]=CP00[ 3]*xint[ 13]+CP10[1]*xint[  7];
    yint[ 16]=CP00[ 4]*yint[ 13]+CP10[1]*yint[  7];
    zint[ 16]=CP00[ 5]*zint[ 13]+CP10[1]*zint[  7];
    CP10[2] += B00[2];
    xint[ 17]=CP00[ 6]*xint[ 14]+CP10[2]*xint[  8];
    yint[ 17]=CP00[ 7]*yint[ 14]+CP10[2]*yint[  8];
    zint[ 17]=CP00[ 8]*zint[ 14]+CP10[2]*zint[  8];
    // (3,0)
    C10[0] += B10[0];
    xint[ 18]=C00[ 0]*xint[ 12]+C10[0]*xint[  6];
    yint[ 18]=C00[ 1]*yint[ 12]+C10[0]*yint[  6];
    zint[ 18]=C00[ 2]*zint[ 12]+C10[0]*zint[  6];
    C10[1] += B10[1];
    xint[ 19]=C00[ 3]*xint[ 13]+C10[1]*xint[  7];
    yint[ 19]=C00[ 4]*yint[ 13]+C10[1]*yint[  7];
    zint[ 19]=C00[ 5]*zint[ 13]+C10[1]*zint[  7];
    C10[2] += B10[2];
    xint[ 20]=C00[ 6]*xint[ 14]+C10[2]*xint[  8];
    yint[ 20]=C00[ 7]*yint[ 14]+C10[2]*yint[  8];
    zint[ 20]=C00[ 8]*zint[ 14]+C10[2]*zint[  8];
    // (3,1)
    CP10[0] += B00[0];
    xint[ 21]=CP00[ 0]*xint[ 18]+CP10[0]*xint[ 12];
    yint[ 21]=CP00[ 1]*yint[ 18]+CP10[0]*yint[ 12];
    zint[ 21]=CP00[ 2]*zint[ 18]+CP10[0]*zint[ 12];
    CP10[1] += B00[1];
    xint[ 22]=CP00[ 3]*xint[ 19]+CP10[1]*xint[ 13];
    yint[ 22]=CP00[ 4]*yint[ 19]+CP10[1]*yint[ 13];
    zint[ 22]=CP00[ 5]*zint[ 19]+CP10[1]*zint[ 13];
    CP10[2] += B00[2];
    xint[ 23]=CP00[ 6]*xint[ 20]+CP10[2]*xint[ 14];
    yint[ 23]=CP00[ 7]*yint[ 20]+CP10[2]*yint[ 14];
    zint[ 23]=CP00[ 8]*zint[ 20]+CP10[2]*zint[ 14];
}

static void ofmo_form_dpps(
        const double *xint, const double *yint, const double *zint,
        double *eh ) {
    // (DS|PS)
    eh[   0] += xint[ 15]*yint[  0]*zint[  0];
    eh[   0] += xint[ 16]*yint[  1]*zint[  1];
    eh[   0] += xint[ 17]*yint[  2]*zint[  2];
    eh[   1] += xint[ 12]*yint[  3]*zint[  0];
    eh[   1] += xint[ 13]*yint[  4]*zint[  1];
    eh[   1] += xint[ 14]*yint[  5]*zint[  2];
    eh[   2] += xint[ 12]*yint[  0]*zint[  3];
    eh[   2] += xint[ 13]*yint[  1]*zint[  4];
    eh[   2] += xint[ 14]*yint[  2]*zint[  5];
    eh[   3] += xint[  3]*yint[ 12]*zint[  0];
    eh[   3] += xint[  4]*yint[ 13]*zint[  1];
    eh[   3] += xint[  5]*yint[ 14]*zint[  2];
    eh[   4] += xint[  0]*yint[ 15]*zint[  0];
    eh[   4] += xint[  1]*yint[ 16]*zint[  1];
    eh[   4] += xint[  2]*yint[ 17]*zint[  2];
    eh[   5] += xint[  0]*yint[ 12]*zint[  3];
    eh[   5] += xint[  1]*yint[ 13]*zint[  4];
    eh[   5] += xint[  2]*yint[ 14]*zint[  5];
    eh[   6] += xint[  3]*yint[  0]*zint[ 12];
    eh[   6] += xint[  4]*yint[  1]*zint[ 13];
    eh[   6] += xint[  5]*yint[  2]*zint[ 14];
    eh[   7] += xint[  0]*yint[  3]*zint[ 12];
    eh[   7] += xint[  1]*yint[  4]*zint[ 13];
    eh[   7] += xint[  2]*yint[  5]*zint[ 14];
    eh[   8] += xint[  0]*yint[  0]*zint[ 15];
    eh[   8] += xint[  1]*yint[  1]*zint[ 16];
    eh[   8] += xint[  2]*yint[  2]*zint[ 17];
    eh[   9] += xint[  9]*yint[  6]*zint[  0];
    eh[   9] += xint[ 10]*yint[  7]*zint[  1];
    eh[   9] += xint[ 11]*yint[  8]*zint[  2];
    eh[  10] += xint[  6]*yint[  9]*zint[  0];
    eh[  10] += xint[  7]*yint[ 10]*zint[  1];
    eh[  10] += xint[  8]*yint[ 11]*zint[  2];
    eh[  11] += xint[  6]*yint[  6]*zint[  3];
    eh[  11] += xint[  7]*yint[  7]*zint[  4];
    eh[  11] += xint[  8]*yint[  8]*zint[  5];
    eh[  12] += xint[  9]*yint[  0]*zint[  6];
    eh[  12] += xint[ 10]*yint[  1]*zint[  7];
    eh[  12] += xint[ 11]*yint[  2]*zint[  8];
    eh[  13] += xint[  6]*yint[  3]*zint[  6];
    eh[  13] += xint[  7]*yint[  4]*zint[  7];
    eh[  13] += xint[  8]*yint[  5]*zint[  8];
    eh[  14] += xint[  6]*yint[  0]*zint[  9];
    eh[  14] += xint[  7]*yint[  1]*zint[ 10];
    eh[  14] += xint[  8]*yint[  2]*zint[ 11];
    eh[  15] += xint[  3]*yint[  6]*zint[  6];
    eh[  15] += xint[  4]*yint[  7]*zint[  7];
    eh[  15] += xint[  5]*yint[  8]*zint[  8];
    eh[  16] += xint[  0]*yint[  9]*zint[  6];
    eh[  16] += xint[  1]*yint[ 10]*zint[  7];
    eh[  16] += xint[  2]*yint[ 11]*zint[  8];
    eh[  17] += xint[  0]*yint[  6]*zint[  9];
    eh[  17] += xint[  1]*yint[  7]*zint[ 10];
    eh[  17] += xint[  2]*yint[  8]*zint[ 11];
    // (FS|PS)
    eh[  18] += xint[ 21]*yint[  0]*zint[  0];
    eh[  18] += xint[ 22]*yint[  1]*zint[  1];
    eh[  18] += xint[ 23]*yint[  2]*zint[  2];
    eh[  19] += xint[ 18]*yint[  3]*zint[  0];
    eh[  19] += xint[ 19]*yint[  4]*zint[  1];
    eh[  19] += xint[ 20]*yint[  5]*zint[  2];
    eh[  20] += xint[ 18]*yint[  0]*zint[  3];
    eh[  20] += xint[ 19]*yint[  1]*zint[  4];
    eh[  20] += xint[ 20]*yint[  2]*zint[  5];
    eh[  21] += xint[  3]*yint[ 18]*zint[  0];
    eh[  21] += xint[  4]*yint[ 19]*zint[  1];
    eh[  21] += xint[  5]*yint[ 20]*zint[  2];
    eh[  22] += xint[  0]*yint[ 21]*zint[  0];
    eh[  22] += xint[  1]*yint[ 22]*zint[  1];
    eh[  22] += xint[  2]*yint[ 23]*zint[  2];
    eh[  23] += xint[  0]*yint[ 18]*zint[  3];
    eh[  23] += xint[  1]*yint[ 19]*zint[  4];
    eh[  23] += xint[  2]*yint[ 20]*zint[  5];
    eh[  24] += xint[  3]*yint[  0]*zint[ 18];
    eh[  24] += xint[  4]*yint[  1]*zint[ 19];
    eh[  24] += xint[  5]*yint[  2]*zint[ 20];
    eh[  25] += xint[  0]*yint[  3]*zint[ 18];
    eh[  25] += xint[  1]*yint[  4]*zint[ 19];
    eh[  25] += xint[  2]*yint[  5]*zint[ 20];
    eh[  26] += xint[  0]*yint[  0]*zint[ 21];
    eh[  26] += xint[  1]*yint[  1]*zint[ 22];
    eh[  26] += xint[  2]*yint[  2]*zint[ 23];
    eh[  27] += xint[ 15]*yint[  6]*zint[  0];
    eh[  27] += xint[ 16]*yint[  7]*zint[  1];
    eh[  27] += xint[ 17]*yint[  8]*zint[  2];
    eh[  28] += xint[ 12]*yint[  9]*zint[  0];
    eh[  28] += xint[ 13]*yint[ 10]*zint[  1];
    eh[  28] += xint[ 14]*yint[ 11]*zint[  2];
    eh[  29] += xint[ 12]*yint[  6]*zint[  3];
    eh[  29] += xint[ 13]*yint[  7]*zint[  4];
    eh[  29] += xint[ 14]*yint[  8]*zint[  5];
    eh[  30] += xint[ 15]*yint[  0]*zint[  6];
    eh[  30] += xint[ 16]*yint[  1]*zint[  7];
    eh[  30] += xint[ 17]*yint[  2]*zint[  8];
    eh[  31] += xint[ 12]*yint[  3]*zint[  6];
    eh[  31] += xint[ 13]*yint[  4]*zint[  7];
    eh[  31] += xint[ 14]*yint[  5]*zint[  8];
    eh[  32] += xint[ 12]*yint[  0]*zint[  9];
    eh[  32] += xint[ 13]*yint[  1]*zint[ 10];
    eh[  32] += xint[ 14]*yint[  2]*zint[ 11];
    eh[  33] += xint[  9]*yint[ 12]*zint[  0];
    eh[  33] += xint[ 10]*yint[ 13]*zint[  1];
    eh[  33] += xint[ 11]*yint[ 14]*zint[  2];
    eh[  34] += xint[  6]*yint[ 15]*zint[  0];
    eh[  34] += xint[  7]*yint[ 16]*zint[  1];
    eh[  34] += xint[  8]*yint[ 17]*zint[  2];
    eh[  35] += xint[  6]*yint[ 12]*zint[  3];
    eh[  35] += xint[  7]*yint[ 13]*zint[  4];
    eh[  35] += xint[  8]*yint[ 14]*zint[  5];
    eh[  36] += xint[  9]*yint[  0]*zint[ 12];
    eh[  36] += xint[ 10]*yint[  1]*zint[ 13];
    eh[  36] += xint[ 11]*yint[  2]*zint[ 14];
    eh[  37] += xint[  6]*yint[  3]*zint[ 12];
    eh[  37] += xint[  7]*yint[  4]*zint[ 13];
    eh[  37] += xint[  8]*yint[  5]*zint[ 14];
    eh[  38] += xint[  6]*yint[  0]*zint[ 15];
    eh[  38] += xint[  7]*yint[  1]*zint[ 16];
    eh[  38] += xint[  8]*yint[  2]*zint[ 17];
    eh[  39] += xint[  9]*yint[  6]*zint[  6];
    eh[  39] += xint[ 10]*yint[  7]*zint[  7];
    eh[  39] += xint[ 11]*yint[  8]*zint[  8];
    eh[  40] += xint[  6]*yint[  9]*zint[  6];
    eh[  40] += xint[  7]*yint[ 10]*zint[  7];
    eh[  40] += xint[  8]*yint[ 11]*zint[  8];
    eh[  41] += xint[  6]*yint[  6]*zint[  9];
    eh[  41] += xint[  7]*yint[  7]*zint[ 10];
    eh[  41] += xint[  8]*yint[  8]*zint[ 11];
    eh[  42] += xint[  3]*yint[ 12]*zint[  6];
    eh[  42] += xint[  4]*yint[ 13]*zint[  7];
    eh[  42] += xint[  5]*yint[ 14]*zint[  8];
    eh[  43] += xint[  0]*yint[ 15]*zint[  6];
    eh[  43] += xint[  1]*yint[ 16]*zint[  7];
    eh[  43] += xint[  2]*yint[ 17]*zint[  8];
    eh[  44] += xint[  0]*yint[ 12]*zint[  9];
    eh[  44] += xint[  1]*yint[ 13]*zint[ 10];
    eh[  44] += xint[  2]*yint[ 14]*zint[ 11];
    eh[  45] += xint[  3]*yint[  6]*zint[ 12];
    eh[  45] += xint[  4]*yint[  7]*zint[ 13];
    eh[  45] += xint[  5]*yint[  8]*zint[ 14];
    eh[  46] += xint[  0]*yint[  9]*zint[ 12];
    eh[  46] += xint[  1]*yint[ 10]*zint[ 13];
    eh[  46] += xint[  2]*yint[ 11]*zint[ 14];
    eh[  47] += xint[  0]*yint[  6]*zint[ 15];
    eh[  47] += xint[  1]*yint[  7]*zint[ 16];
    eh[  47] += xint[  2]*yint[  8]*zint[ 17];
}

void ofmo_twoint_core_rys_dpps( const int mythread,
        const int *nijps, const double *vzeta, const double *vdkab,
        const double vxiza[], const double BA[3],
        const int *nklps, const double *veta, const double *vdkcd,
        const double *vxizc, const double DC[3], const double AC[3],
        double *DINT ) {
    int ijps, klps, i;
    double cssss, zeta, dkab, xiza, eta, xizc, dk, T;
    double zeta2, eta2, rz, PA[3], QC[3];
    double PQ2, sqrho, rho, PC[3], QP[3];
    double C00[9], CP00[9], B00[3], B10[3], B01[3], F00[3];
    double rrho, rze, W[13], U[13];
    double u2, duminv, dm2inv, dum;
    int m, m3;
    double *xint, *yint, *zint, *eh;
    xint = ofmo_integ_getadd_xint( mythread );
    yint = ofmo_integ_getadd_yint( mythread );
    zint = ofmo_integ_getadd_zint( mythread );
    eh   = ofmo_integ_getadd_eh( mythread );
    DFACT = ofmo_getadd_dfact();
    ofmo_hrr_clear_dpps( eh );
    for ( ijps=0; ijps<(*nijps); ijps++ ) {
        zeta  = vzeta[ijps];
        dkab  = vdkab[ijps];
        xiza  = vxiza[ijps];
        zeta2 = HALF * zeta;
        for ( i=0; i<3; i++ ) {
            PC[i] = AC[i] + xiza*BA[i];
            PA[i] = xiza * BA[i];
        }
        for ( klps=0; klps<(*nklps); klps++ ) {
            eta  = veta[klps];
            dk   = dkab * vdkcd[klps];
            xizc = vxizc[klps];
            eta2 = HALF * eta;
            PQ2  = ZERO;
            for ( i=0; i<3; i++ ) {
                QC[i] = xizc*DC[i];
                QP[i] = xizc*DC[i] - PC[i];
                PQ2  += QP[i]*QP[i];
            }
            rrho  = zeta + eta;
            rze   = zeta * eta;
            sqrho = sqrt(1.e0/rrho);
            rho   = sqrho * sqrho;
            rz    = rho * zeta;
            T     = rho * PQ2;
            cssss = sqrho * dk;
            calc_root( 3, T, U, W );
            for ( m=m3=0; m<3; m++, m3+=3 ) {
                u2     = rho * U[m];
                F00[m] = cssss * W[m];
                duminv = 1.e0 / ( 1.e0 + rrho * u2 );
                dm2inv = 0.5e0 * duminv;
                B00[m] = dm2inv * rze * u2;
                B10[m] = dm2inv * ( zeta + rze*u2 );
                B01[m] = dm2inv * ( eta  + rze*u2 );
                dum    = zeta * u2 * duminv;
                for ( i=0; i<3; i++ ) C00[m3+i]  = PA[i] + dum * QP[i];
                dum    = eta * u2 * duminv;
                for ( i=0; i<3; i++ ) CP00[m3+i] = QC[i] - dum * QP[i];
            }
            ofmo_xyzint_dpps( F00, B00, B10, B01, C00, CP00,
                     xint, yint, zint );
            ofmo_form_dpps( xint, yint, zint, eh );
        }
    }
    ofmo_hrr_calc_dpps( eh, BA, DC );
    ofmo_hrr_coef_dpps( eh, DINT );
}

int ofmo_twoint_rys_dpps(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
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
        int *last_ijcs, int *last_klcs ) {
    int Lab, Lcd, i, j, k, l, ipat, ix;
    int I2, IJ, K2, KL;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1, max_klcs;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe, coe0;
    double *DINTEG;
    long nzeri, max_nzeri, nzeri4;
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    long ebuf_max_nzeri = *pebuf_max_nzeri;
    int mythread;
    DFACT = ofmo_getadd_dfact();
    mythread = omp_get_thread_num();
    DINTEG = ofmo_integ_getadd_eri( mythread );
    Lab = La*(La+1)/2+Lb;
    Lcd = Lc*(Lc+1)/2+Ld;
    ijcs0 = leading_cs_pair[Lab];
    ijcs1 = leading_cs_pair[Lab+1];
    klcs0 = leading_cs_pair[Lcd];
    klcs1 = leading_cs_pair[Lcd+1];
    nzeri     = *ebuf_non_zero_eri;
    max_nzeri = ebuf_max_nzeri - 6*3*3*1;
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
        max_klcs = ( Lab == Lcd ? ijcs+1 : klcs1 );
        for ( klcs=klcs0; klcs<max_klcs; klcs++ ) {
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
            ofmo_twoint_core_rys_dpps( mythread,
                    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
                    &psp_xiza[ijps0], BA,
                    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
                    &psp_xiza[klps0], DC,   AC,      DINTEG );
            ipat=((Lab != Lcd) || (ics==kcs && jcs>lcs) ? true : false);
            for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
                I2 = (iao*iao+iao)>>1;
                for ( j=0, jao=jao0; j<3; j++, jao++ ) {
                    if ( jao>iao ) { ix+=3*1; continue; }
                    IJ = I2 + jao;
                    coe0 = ( iao==jao ? HALF : ONE );
                    for ( k=0, kao=kao0; k<3; k++, kao++ ) {
                        K2 = (kao*kao+kao)>>1;
                        for ( l=0, lao=lao0; l<1; l++, lao++, ix++ ) {
                            if ( lao>kao ) continue;
                            if ( fabs(DINTEG[ix]) > EPS_ERI ) {
                                KL = K2 + lao;
                                if ( IJ >= KL ) {
                                    coe = coe0;
                                    if ( kao==lao ) coe *= HALF;
                                    if ( KL == IJ ) coe *= HALF;
                                    ebuf_val[nzeri]     = coe*DINTEG[ix];
                                    ebuf_ind4[nzeri4+0] = (short int)iao;
                                    ebuf_ind4[nzeri4+1] = (short int)jao;
                                    ebuf_ind4[nzeri4+2] = (short int)kao;
                                    ebuf_ind4[nzeri4+3] = (short int)lao;
                                    nzeri++;
                                    nzeri4+=4;
                                } else if ( ipat ) {
                                    coe = coe0;
                                    if ( kao==lao ) coe*=HALF;
                                    ebuf_val[nzeri]     = coe*DINTEG[ix];
                                    ebuf_ind4[nzeri4+0] = (short int)kao;
                                    ebuf_ind4[nzeri4+1] = (short int)lao;
                                    ebuf_ind4[nzeri4+2] = (short int)iao;
                                    ebuf_ind4[nzeri4+3] = (short int)jao;
                                    nzeri++;
                                    nzeri4+=4;
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
        }
    }
    *ebuf_non_zero_eri = nzeri;
    return OFMO_EBUF_NOFULL;
}

int ofmo_twoint_direct_rys_dpps(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *shel_atm, const int *shel_ini,
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
        const int *pnao, const double Ds[], double G[] ) {
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    int last_ijcs=*plast_ijcs, last_klcs=*plast_klcs, nao=*pnao;
    long max_nzeri=*petmp_max_nzeri;
    long nzeri4, nzeri=*petmp_non_zero_eri;
    //
    int Lab, Lcd, i, j, k, l, ipat, ix;
    int I2, IJ, K2, KL;
    int ijcs, ijcs0, ijcs1;
    int klcs, klcs0, klcs1, max_klcs;
    int ijps0, nijps, klps0, nklps;
    int ics, iat, iao, iao0, jcs, jat, jao, jao0;
    int kcs, kat, kao, kao0, lcs, lat, lao, lao0;
    double A[3], B[3], C[3], D[3], BA[3], DC[3], AC[3];
    double val_ab, val_cd, coe, coe0;
    double *DINTEG;
    int mythread;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);
    float eps_sch = ofmo_twoint_eps_sch(0);

    DFACT = ofmo_getadd_dfact();
    mythread = omp_get_thread_num();
    DINTEG = ofmo_integ_getadd_eri( mythread );
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
    max_nzeri -= 6*3*3*1;
    nzeri4    = nzeri*4;
    if ( nzeri >= max_nzeri ) {
        ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4, Ds, G );
        nzeri = nzeri4 = 0;
    }
    for ( ; ijcs<ijcs1; ijcs+=nworkers ) {
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
        max_klcs = ( Lab == Lcd ? ijcs+1 : klcs1 );
        for ( ; klcs<max_klcs; klcs++ ) {
            val_cd = csp_schwarz[klcs];
            if ( val_ab*val_cd < eps_ps4 ) continue;
            kcs    = csp_ics[klcs];
            lcs    = csp_jcs[klcs];
            if ( val_ab*val_cd*ofmo_twoint_dmax6(ics,jcs,kcs,lcs) < eps_sch ) continue;
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
            ofmo_twoint_core_rys_dpps( mythread,
                    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
                    &psp_xiza[ijps0], BA,
                    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
                    &psp_xiza[klps0], DC,   AC,      DINTEG );
            ipat=((Lab != Lcd) || (ics==kcs && jcs>lcs) ? true : false);
            for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
                I2 = (iao*iao+iao)>>1;
                for ( j=0, jao=jao0; j<3; j++, jao++ ) {
                    if ( jao>iao ) { ix+=3*1; continue; }
                    IJ = I2 + jao;
                    coe0 = ( iao==jao ? HALF : ONE );
                    for ( k=0, kao=kao0; k<3; k++, kao++ ) {
                        K2 = (kao*kao+kao)>>1;
                        for ( l=0, lao=lao0; l<1; l++, lao++, ix++ ) {
                            if ( lao>kao ) continue;
                            if ( fabs(DINTEG[ix]) > eps_eri ) {
                                KL = K2 + lao;
                                if ( IJ >= KL ) {
                                    coe = coe0;
                                    if ( kao==lao ) coe *= HALF;
                                    if ( KL == IJ ) coe *= HALF;
                                    etmp_val[nzeri]     = coe*DINTEG[ix];
                                    etmp_ind4[nzeri4+0] = (short int)iao;
                                    etmp_ind4[nzeri4+1] = (short int)jao;
                                    etmp_ind4[nzeri4+2] = (short int)kao;
                                    etmp_ind4[nzeri4+3] = (short int)lao;
                                    nzeri++;
                                    nzeri4+=4;
                                } else if ( ipat ) {
                                    coe = coe0;
                                    if ( kao==lao ) coe*=HALF;
                                    etmp_val[nzeri]     = coe*DINTEG[ix];
                                    etmp_ind4[nzeri4+0] = (short int)kao;
                                    etmp_ind4[nzeri4+1] = (short int)lao;
                                    etmp_ind4[nzeri4+2] = (short int)iao;
                                    etmp_ind4[nzeri4+3] = (short int)jao;
                                    nzeri++;
                                    nzeri4+=4;
                                }
                            }
                        }
                    }
                }
            }
            if ( nzeri >= max_nzeri ) {
                ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4,
                        Ds, G );
                nzeri = nzeri4= 0;
            }
        }
        klcs = klcs0;
    }
    *petmp_non_zero_eri = nzeri;
    return 0;
}
