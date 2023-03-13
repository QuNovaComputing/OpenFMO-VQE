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

static void ofmo_hrr_clear_dpds( double *eh ) {
    int i;
    // (DS|DS)
    for ( i=0; i<(0+36); i++ ) eh[i] = 0.e0;
    // (FS|DS)
    for ( i=36; i<(36+60); i++ ) eh[i] = 0.e0;
}

static void ofmo_hrr_coef_dpds(
        double *eh, double *DINT ) {
    int i, j, k, l, iao, jao, kao, lao, ix;
    double coef_a, coef_ab, coef_abc;
    double *th;
    th = &eh[96];
    ix = 0;
    for ( i=0, iao=4; i<6; i++, iao++ ) {
        coef_a = DFACT[iao];
        for ( j=0, jao=1; j<3; j++, jao++ ) {
            coef_ab = coef_a * DFACT[jao];
            for ( k=0, kao=4; k<6; k++, kao++ ) {
                coef_abc = coef_ab * DFACT[kao];
                for ( l=0, lao=0; l<1; l++, lao++ ) {
                    DINT[ix] = coef_abc * DFACT[lao] * th[ix];
                    ix++;
                }
            }
        }
    }
}

static void ofmo_hrr_calc_dpds( double *eh,
        const double BA[3], const double DC[3] ) {
    // HRR for (XX|XS)-type integral (center AB)
    // (DP,DS)
    eh[  96] = eh[  36] - BA[0]*eh[   0];
    eh[  97] = eh[  37] - BA[0]*eh[   1];
    eh[  98] = eh[  38] - BA[0]*eh[   2];
    eh[  99] = eh[  39] - BA[0]*eh[   3];
    eh[ 100] = eh[  40] - BA[0]*eh[   4];
    eh[ 101] = eh[  41] - BA[0]*eh[   5];
    eh[ 102] = eh[  54] - BA[1]*eh[   0];
    eh[ 103] = eh[  55] - BA[1]*eh[   1];
    eh[ 104] = eh[  56] - BA[1]*eh[   2];
    eh[ 105] = eh[  57] - BA[1]*eh[   3];
    eh[ 106] = eh[  58] - BA[1]*eh[   4];
    eh[ 107] = eh[  59] - BA[1]*eh[   5];
    eh[ 108] = eh[  60] - BA[2]*eh[   0];
    eh[ 109] = eh[  61] - BA[2]*eh[   1];
    eh[ 110] = eh[  62] - BA[2]*eh[   2];
    eh[ 111] = eh[  63] - BA[2]*eh[   3];
    eh[ 112] = eh[  64] - BA[2]*eh[   4];
    eh[ 113] = eh[  65] - BA[2]*eh[   5];
    eh[ 114] = eh[  66] - BA[0]*eh[   6];
    eh[ 115] = eh[  67] - BA[0]*eh[   7];
    eh[ 116] = eh[  68] - BA[0]*eh[   8];
    eh[ 117] = eh[  69] - BA[0]*eh[   9];
    eh[ 118] = eh[  70] - BA[0]*eh[  10];
    eh[ 119] = eh[  71] - BA[0]*eh[  11];
    eh[ 120] = eh[  42] - BA[1]*eh[   6];
    eh[ 121] = eh[  43] - BA[1]*eh[   7];
    eh[ 122] = eh[  44] - BA[1]*eh[   8];
    eh[ 123] = eh[  45] - BA[1]*eh[   9];
    eh[ 124] = eh[  46] - BA[1]*eh[  10];
    eh[ 125] = eh[  47] - BA[1]*eh[  11];
    eh[ 126] = eh[  84] - BA[2]*eh[   6];
    eh[ 127] = eh[  85] - BA[2]*eh[   7];
    eh[ 128] = eh[  86] - BA[2]*eh[   8];
    eh[ 129] = eh[  87] - BA[2]*eh[   9];
    eh[ 130] = eh[  88] - BA[2]*eh[  10];
    eh[ 131] = eh[  89] - BA[2]*eh[  11];
    eh[ 132] = eh[  72] - BA[0]*eh[  12];
    eh[ 133] = eh[  73] - BA[0]*eh[  13];
    eh[ 134] = eh[  74] - BA[0]*eh[  14];
    eh[ 135] = eh[  75] - BA[0]*eh[  15];
    eh[ 136] = eh[  76] - BA[0]*eh[  16];
    eh[ 137] = eh[  77] - BA[0]*eh[  17];
    eh[ 138] = eh[  90] - BA[1]*eh[  12];
    eh[ 139] = eh[  91] - BA[1]*eh[  13];
    eh[ 140] = eh[  92] - BA[1]*eh[  14];
    eh[ 141] = eh[  93] - BA[1]*eh[  15];
    eh[ 142] = eh[  94] - BA[1]*eh[  16];
    eh[ 143] = eh[  95] - BA[1]*eh[  17];
    eh[ 144] = eh[  48] - BA[2]*eh[  12];
    eh[ 145] = eh[  49] - BA[2]*eh[  13];
    eh[ 146] = eh[  50] - BA[2]*eh[  14];
    eh[ 147] = eh[  51] - BA[2]*eh[  15];
    eh[ 148] = eh[  52] - BA[2]*eh[  16];
    eh[ 149] = eh[  53] - BA[2]*eh[  17];
    eh[ 150] = eh[  54] - BA[0]*eh[  18];
    eh[ 151] = eh[  55] - BA[0]*eh[  19];
    eh[ 152] = eh[  56] - BA[0]*eh[  20];
    eh[ 153] = eh[  57] - BA[0]*eh[  21];
    eh[ 154] = eh[  58] - BA[0]*eh[  22];
    eh[ 155] = eh[  59] - BA[0]*eh[  23];
    eh[ 156] = eh[  66] - BA[1]*eh[  18];
    eh[ 157] = eh[  67] - BA[1]*eh[  19];
    eh[ 158] = eh[  68] - BA[1]*eh[  20];
    eh[ 159] = eh[  69] - BA[1]*eh[  21];
    eh[ 160] = eh[  70] - BA[1]*eh[  22];
    eh[ 161] = eh[  71] - BA[1]*eh[  23];
    eh[ 162] = eh[  78] - BA[2]*eh[  18];
    eh[ 163] = eh[  79] - BA[2]*eh[  19];
    eh[ 164] = eh[  80] - BA[2]*eh[  20];
    eh[ 165] = eh[  81] - BA[2]*eh[  21];
    eh[ 166] = eh[  82] - BA[2]*eh[  22];
    eh[ 167] = eh[  83] - BA[2]*eh[  23];
    eh[ 168] = eh[  60] - BA[0]*eh[  24];
    eh[ 169] = eh[  61] - BA[0]*eh[  25];
    eh[ 170] = eh[  62] - BA[0]*eh[  26];
    eh[ 171] = eh[  63] - BA[0]*eh[  27];
    eh[ 172] = eh[  64] - BA[0]*eh[  28];
    eh[ 173] = eh[  65] - BA[0]*eh[  29];
    eh[ 174] = eh[  78] - BA[1]*eh[  24];
    eh[ 175] = eh[  79] - BA[1]*eh[  25];
    eh[ 176] = eh[  80] - BA[1]*eh[  26];
    eh[ 177] = eh[  81] - BA[1]*eh[  27];
    eh[ 178] = eh[  82] - BA[1]*eh[  28];
    eh[ 179] = eh[  83] - BA[1]*eh[  29];
    eh[ 180] = eh[  72] - BA[2]*eh[  24];
    eh[ 181] = eh[  73] - BA[2]*eh[  25];
    eh[ 182] = eh[  74] - BA[2]*eh[  26];
    eh[ 183] = eh[  75] - BA[2]*eh[  27];
    eh[ 184] = eh[  76] - BA[2]*eh[  28];
    eh[ 185] = eh[  77] - BA[2]*eh[  29];
    eh[ 186] = eh[  78] - BA[0]*eh[  30];
    eh[ 187] = eh[  79] - BA[0]*eh[  31];
    eh[ 188] = eh[  80] - BA[0]*eh[  32];
    eh[ 189] = eh[  81] - BA[0]*eh[  33];
    eh[ 190] = eh[  82] - BA[0]*eh[  34];
    eh[ 191] = eh[  83] - BA[0]*eh[  35];
    eh[ 192] = eh[  84] - BA[1]*eh[  30];
    eh[ 193] = eh[  85] - BA[1]*eh[  31];
    eh[ 194] = eh[  86] - BA[1]*eh[  32];
    eh[ 195] = eh[  87] - BA[1]*eh[  33];
    eh[ 196] = eh[  88] - BA[1]*eh[  34];
    eh[ 197] = eh[  89] - BA[1]*eh[  35];
    eh[ 198] = eh[  90] - BA[2]*eh[  30];
    eh[ 199] = eh[  91] - BA[2]*eh[  31];
    eh[ 200] = eh[  92] - BA[2]*eh[  32];
    eh[ 201] = eh[  93] - BA[2]*eh[  33];
    eh[ 202] = eh[  94] - BA[2]*eh[  34];
    eh[ 203] = eh[  95] - BA[2]*eh[  35];
}

static void ofmo_xyzint_dpds(
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
    xint[  9]=C00[ 0];
    yint[  9]=C00[ 1];
    zint[  9]=C00[ 2]*F00[0];
    xint[ 10]=C00[ 3];
    yint[ 10]=C00[ 4];
    zint[ 10]=C00[ 5]*F00[1];
    xint[ 11]=C00[ 6];
    yint[ 11]=C00[ 7];
    zint[ 11]=C00[ 8]*F00[2];
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
    xint[ 12]=CP00[ 0]*xint[  9]+B00[0];
    yint[ 12]=CP00[ 1]*yint[  9]+B00[0];
    zint[ 12]=CP00[ 2]*zint[  9]+B00[0]*F00[0];
    xint[ 13]=CP00[ 3]*xint[ 10]+B00[1];
    yint[ 13]=CP00[ 4]*yint[ 10]+B00[1];
    zint[ 13]=CP00[ 5]*zint[ 10]+B00[1]*F00[1];
    xint[ 14]=CP00[ 6]*xint[ 11]+B00[2];
    yint[ 14]=CP00[ 7]*yint[ 11]+B00[2];
    zint[ 14]=CP00[ 8]*zint[ 11]+B00[2]*F00[2];
    // (N,0) and (N,1)
    for ( m=0; m<3; m++ ) {
        C10[m]  = 0.e0;
        CP10[m] = B00[m];
    }
    // (2,0)
    C10[0] += B10[0];
    xint[ 18]=C00[ 0]*xint[  9]+C10[0]*xint[  0];
    yint[ 18]=C00[ 1]*yint[  9]+C10[0]*yint[  0];
    zint[ 18]=C00[ 2]*zint[  9]+C10[0]*zint[  0];
    C10[1] += B10[1];
    xint[ 19]=C00[ 3]*xint[ 10]+C10[1]*xint[  1];
    yint[ 19]=C00[ 4]*yint[ 10]+C10[1]*yint[  1];
    zint[ 19]=C00[ 5]*zint[ 10]+C10[1]*zint[  1];
    C10[2] += B10[2];
    xint[ 20]=C00[ 6]*xint[ 11]+C10[2]*xint[  2];
    yint[ 20]=C00[ 7]*yint[ 11]+C10[2]*yint[  2];
    zint[ 20]=C00[ 8]*zint[ 11]+C10[2]*zint[  2];
    // (2,1)
    CP10[0] += B00[0];
    xint[ 21]=CP00[ 0]*xint[ 18]+CP10[0]*xint[  9];
    yint[ 21]=CP00[ 1]*yint[ 18]+CP10[0]*yint[  9];
    zint[ 21]=CP00[ 2]*zint[ 18]+CP10[0]*zint[  9];
    CP10[1] += B00[1];
    xint[ 22]=CP00[ 3]*xint[ 19]+CP10[1]*xint[ 10];
    yint[ 22]=CP00[ 4]*yint[ 19]+CP10[1]*yint[ 10];
    zint[ 22]=CP00[ 5]*zint[ 19]+CP10[1]*zint[ 10];
    CP10[2] += B00[2];
    xint[ 23]=CP00[ 6]*xint[ 20]+CP10[2]*xint[ 11];
    yint[ 23]=CP00[ 7]*yint[ 20]+CP10[2]*yint[ 11];
    zint[ 23]=CP00[ 8]*zint[ 20]+CP10[2]*zint[ 11];
    // (3,0)
    C10[0] += B10[0];
    xint[ 27]=C00[ 0]*xint[ 18]+C10[0]*xint[  9];
    yint[ 27]=C00[ 1]*yint[ 18]+C10[0]*yint[  9];
    zint[ 27]=C00[ 2]*zint[ 18]+C10[0]*zint[  9];
    C10[1] += B10[1];
    xint[ 28]=C00[ 3]*xint[ 19]+C10[1]*xint[ 10];
    yint[ 28]=C00[ 4]*yint[ 19]+C10[1]*yint[ 10];
    zint[ 28]=C00[ 5]*zint[ 19]+C10[1]*zint[ 10];
    C10[2] += B10[2];
    xint[ 29]=C00[ 6]*xint[ 20]+C10[2]*xint[ 11];
    yint[ 29]=C00[ 7]*yint[ 20]+C10[2]*yint[ 11];
    zint[ 29]=C00[ 8]*zint[ 20]+C10[2]*zint[ 11];
    // (3,1)
    CP10[0] += B00[0];
    xint[ 30]=CP00[ 0]*xint[ 27]+CP10[0]*xint[ 18];
    yint[ 30]=CP00[ 1]*yint[ 27]+CP10[0]*yint[ 18];
    zint[ 30]=CP00[ 2]*zint[ 27]+CP10[0]*zint[ 18];
    CP10[1] += B00[1];
    xint[ 31]=CP00[ 3]*xint[ 28]+CP10[1]*xint[ 19];
    yint[ 31]=CP00[ 4]*yint[ 28]+CP10[1]*yint[ 19];
    zint[ 31]=CP00[ 5]*zint[ 28]+CP10[1]*zint[ 19];
    CP10[2] += B00[2];
    xint[ 32]=CP00[ 6]*xint[ 29]+CP10[2]*xint[ 20];
    yint[ 32]=CP00[ 7]*yint[ 29]+CP10[2]*yint[ 20];
    zint[ 32]=CP00[ 8]*zint[ 29]+CP10[2]*zint[ 20];
    // (0,M) and (1,M)
    for ( m=0; m<3; m++ ) {
        CP01[m] = 0.e0;
        C01[m]  = B00[m];
    }

    // (0,2)
    CP01[0] += B01[0];
    xint[  6]=CP00[ 0]*xint[  3]+CP01[0]*xint[  0];
    yint[  6]=CP00[ 1]*yint[  3]+CP01[0]*yint[  0];
    zint[  6]=CP00[ 2]*zint[  3]+CP01[0]*zint[  0];
    CP01[1] += B01[1];
    xint[  7]=CP00[ 3]*xint[  4]+CP01[1]*xint[  1];
    yint[  7]=CP00[ 4]*yint[  4]+CP01[1]*yint[  1];
    zint[  7]=CP00[ 5]*zint[  4]+CP01[1]*zint[  1];
    CP01[2] += B01[2];
    xint[  8]=CP00[ 6]*xint[  5]+CP01[2]*xint[  2];
    yint[  8]=CP00[ 7]*yint[  5]+CP01[2]*yint[  2];
    zint[  8]=CP00[ 8]*zint[  5]+CP01[2]*zint[  2];
    // (1,2)
    C01[0] += B00[0];
    xint[ 15]=C00[ 0]*xint[  6]+C01[0]*xint[  3];
    yint[ 15]=C00[ 1]*yint[  6]+C01[0]*yint[  3];
    zint[ 15]=C00[ 2]*zint[  6]+C01[0]*zint[  3];
    C01[1] += B00[1];
    xint[ 16]=C00[ 3]*xint[  7]+C01[1]*xint[  4];
    yint[ 16]=C00[ 4]*yint[  7]+C01[1]*yint[  4];
    zint[ 16]=C00[ 5]*zint[  7]+C01[1]*zint[  4];
    C01[2] += B00[2];
    xint[ 17]=C00[ 6]*xint[  8]+C01[2]*xint[  5];
    yint[ 17]=C00[ 7]*yint[  8]+C01[2]*yint[  5];
    zint[ 17]=C00[ 8]*zint[  8]+C01[2]*zint[  5];
    // (N,M)
    for ( m=0; m<3; m++ ) C01[m] = B00[m];
    for ( m=0; m<3; m++ ) {
        C01[m] += B00[m];
        C10[m]  = B10[m];
    }
    // (2,2)
    xint[ 24]=C00[ 0]*xint[ 15]+C10[0]*xint[  6]+C01[0]*xint[ 12];
    yint[ 24]=C00[ 1]*yint[ 15]+C10[0]*yint[  6]+C01[0]*yint[ 12];
    zint[ 24]=C00[ 2]*zint[ 15]+C10[0]*zint[  6]+C01[0]*zint[ 12];
    C10[0] += B10[0];
    xint[ 25]=C00[ 3]*xint[ 16]+C10[1]*xint[  7]+C01[1]*xint[ 13];
    yint[ 25]=C00[ 4]*yint[ 16]+C10[1]*yint[  7]+C01[1]*yint[ 13];
    zint[ 25]=C00[ 5]*zint[ 16]+C10[1]*zint[  7]+C01[1]*zint[ 13];
    C10[1] += B10[1];
    xint[ 26]=C00[ 6]*xint[ 17]+C10[2]*xint[  8]+C01[2]*xint[ 14];
    yint[ 26]=C00[ 7]*yint[ 17]+C10[2]*yint[  8]+C01[2]*yint[ 14];
    zint[ 26]=C00[ 8]*zint[ 17]+C10[2]*zint[  8]+C01[2]*zint[ 14];
    C10[2] += B10[2];
    // (3,2)
    xint[ 33]=C00[ 0]*xint[ 24]+C10[0]*xint[ 15]+C01[0]*xint[ 21];
    yint[ 33]=C00[ 1]*yint[ 24]+C10[0]*yint[ 15]+C01[0]*yint[ 21];
    zint[ 33]=C00[ 2]*zint[ 24]+C10[0]*zint[ 15]+C01[0]*zint[ 21];
    C10[0] += B10[0];
    xint[ 34]=C00[ 3]*xint[ 25]+C10[1]*xint[ 16]+C01[1]*xint[ 22];
    yint[ 34]=C00[ 4]*yint[ 25]+C10[1]*yint[ 16]+C01[1]*yint[ 22];
    zint[ 34]=C00[ 5]*zint[ 25]+C10[1]*zint[ 16]+C01[1]*zint[ 22];
    C10[1] += B10[1];
    xint[ 35]=C00[ 6]*xint[ 26]+C10[2]*xint[ 17]+C01[2]*xint[ 23];
    yint[ 35]=C00[ 7]*yint[ 26]+C10[2]*yint[ 17]+C01[2]*yint[ 23];
    zint[ 35]=C00[ 8]*zint[ 26]+C10[2]*zint[ 17]+C01[2]*zint[ 23];
    C10[2] += B10[2];
}

static void ofmo_form_dpds(
        const double *xint, const double *yint, const double *zint,
        double *eh ) {
    // (DS|DS)
    eh[   0] += xint[ 24]*yint[  0]*zint[  0];
    eh[   0] += xint[ 25]*yint[  1]*zint[  1];
    eh[   0] += xint[ 26]*yint[  2]*zint[  2];
    eh[   1] += xint[ 18]*yint[  6]*zint[  0];
    eh[   1] += xint[ 19]*yint[  7]*zint[  1];
    eh[   1] += xint[ 20]*yint[  8]*zint[  2];
    eh[   2] += xint[ 18]*yint[  0]*zint[  6];
    eh[   2] += xint[ 19]*yint[  1]*zint[  7];
    eh[   2] += xint[ 20]*yint[  2]*zint[  8];
    eh[   3] += xint[ 21]*yint[  3]*zint[  0];
    eh[   3] += xint[ 22]*yint[  4]*zint[  1];
    eh[   3] += xint[ 23]*yint[  5]*zint[  2];
    eh[   4] += xint[ 21]*yint[  0]*zint[  3];
    eh[   4] += xint[ 22]*yint[  1]*zint[  4];
    eh[   4] += xint[ 23]*yint[  2]*zint[  5];
    eh[   5] += xint[ 18]*yint[  3]*zint[  3];
    eh[   5] += xint[ 19]*yint[  4]*zint[  4];
    eh[   5] += xint[ 20]*yint[  5]*zint[  5];
    eh[   6] += xint[  6]*yint[ 18]*zint[  0];
    eh[   6] += xint[  7]*yint[ 19]*zint[  1];
    eh[   6] += xint[  8]*yint[ 20]*zint[  2];
    eh[   7] += xint[  0]*yint[ 24]*zint[  0];
    eh[   7] += xint[  1]*yint[ 25]*zint[  1];
    eh[   7] += xint[  2]*yint[ 26]*zint[  2];
    eh[   8] += xint[  0]*yint[ 18]*zint[  6];
    eh[   8] += xint[  1]*yint[ 19]*zint[  7];
    eh[   8] += xint[  2]*yint[ 20]*zint[  8];
    eh[   9] += xint[  3]*yint[ 21]*zint[  0];
    eh[   9] += xint[  4]*yint[ 22]*zint[  1];
    eh[   9] += xint[  5]*yint[ 23]*zint[  2];
    eh[  10] += xint[  3]*yint[ 18]*zint[  3];
    eh[  10] += xint[  4]*yint[ 19]*zint[  4];
    eh[  10] += xint[  5]*yint[ 20]*zint[  5];
    eh[  11] += xint[  0]*yint[ 21]*zint[  3];
    eh[  11] += xint[  1]*yint[ 22]*zint[  4];
    eh[  11] += xint[  2]*yint[ 23]*zint[  5];
    eh[  12] += xint[  6]*yint[  0]*zint[ 18];
    eh[  12] += xint[  7]*yint[  1]*zint[ 19];
    eh[  12] += xint[  8]*yint[  2]*zint[ 20];
    eh[  13] += xint[  0]*yint[  6]*zint[ 18];
    eh[  13] += xint[  1]*yint[  7]*zint[ 19];
    eh[  13] += xint[  2]*yint[  8]*zint[ 20];
    eh[  14] += xint[  0]*yint[  0]*zint[ 24];
    eh[  14] += xint[  1]*yint[  1]*zint[ 25];
    eh[  14] += xint[  2]*yint[  2]*zint[ 26];
    eh[  15] += xint[  3]*yint[  3]*zint[ 18];
    eh[  15] += xint[  4]*yint[  4]*zint[ 19];
    eh[  15] += xint[  5]*yint[  5]*zint[ 20];
    eh[  16] += xint[  3]*yint[  0]*zint[ 21];
    eh[  16] += xint[  4]*yint[  1]*zint[ 22];
    eh[  16] += xint[  5]*yint[  2]*zint[ 23];
    eh[  17] += xint[  0]*yint[  3]*zint[ 21];
    eh[  17] += xint[  1]*yint[  4]*zint[ 22];
    eh[  17] += xint[  2]*yint[  5]*zint[ 23];
    eh[  18] += xint[ 15]*yint[  9]*zint[  0];
    eh[  18] += xint[ 16]*yint[ 10]*zint[  1];
    eh[  18] += xint[ 17]*yint[ 11]*zint[  2];
    eh[  19] += xint[  9]*yint[ 15]*zint[  0];
    eh[  19] += xint[ 10]*yint[ 16]*zint[  1];
    eh[  19] += xint[ 11]*yint[ 17]*zint[  2];
    eh[  20] += xint[  9]*yint[  9]*zint[  6];
    eh[  20] += xint[ 10]*yint[ 10]*zint[  7];
    eh[  20] += xint[ 11]*yint[ 11]*zint[  8];
    eh[  21] += xint[ 12]*yint[ 12]*zint[  0];
    eh[  21] += xint[ 13]*yint[ 13]*zint[  1];
    eh[  21] += xint[ 14]*yint[ 14]*zint[  2];
    eh[  22] += xint[ 12]*yint[  9]*zint[  3];
    eh[  22] += xint[ 13]*yint[ 10]*zint[  4];
    eh[  22] += xint[ 14]*yint[ 11]*zint[  5];
    eh[  23] += xint[  9]*yint[ 12]*zint[  3];
    eh[  23] += xint[ 10]*yint[ 13]*zint[  4];
    eh[  23] += xint[ 11]*yint[ 14]*zint[  5];
    eh[  24] += xint[ 15]*yint[  0]*zint[  9];
    eh[  24] += xint[ 16]*yint[  1]*zint[ 10];
    eh[  24] += xint[ 17]*yint[  2]*zint[ 11];
    eh[  25] += xint[  9]*yint[  6]*zint[  9];
    eh[  25] += xint[ 10]*yint[  7]*zint[ 10];
    eh[  25] += xint[ 11]*yint[  8]*zint[ 11];
    eh[  26] += xint[  9]*yint[  0]*zint[ 15];
    eh[  26] += xint[ 10]*yint[  1]*zint[ 16];
    eh[  26] += xint[ 11]*yint[  2]*zint[ 17];
    eh[  27] += xint[ 12]*yint[  3]*zint[  9];
    eh[  27] += xint[ 13]*yint[  4]*zint[ 10];
    eh[  27] += xint[ 14]*yint[  5]*zint[ 11];
    eh[  28] += xint[ 12]*yint[  0]*zint[ 12];
    eh[  28] += xint[ 13]*yint[  1]*zint[ 13];
    eh[  28] += xint[ 14]*yint[  2]*zint[ 14];
    eh[  29] += xint[  9]*yint[  3]*zint[ 12];
    eh[  29] += xint[ 10]*yint[  4]*zint[ 13];
    eh[  29] += xint[ 11]*yint[  5]*zint[ 14];
    eh[  30] += xint[  6]*yint[  9]*zint[  9];
    eh[  30] += xint[  7]*yint[ 10]*zint[ 10];
    eh[  30] += xint[  8]*yint[ 11]*zint[ 11];
    eh[  31] += xint[  0]*yint[ 15]*zint[  9];
    eh[  31] += xint[  1]*yint[ 16]*zint[ 10];
    eh[  31] += xint[  2]*yint[ 17]*zint[ 11];
    eh[  32] += xint[  0]*yint[  9]*zint[ 15];
    eh[  32] += xint[  1]*yint[ 10]*zint[ 16];
    eh[  32] += xint[  2]*yint[ 11]*zint[ 17];
    eh[  33] += xint[  3]*yint[ 12]*zint[  9];
    eh[  33] += xint[  4]*yint[ 13]*zint[ 10];
    eh[  33] += xint[  5]*yint[ 14]*zint[ 11];
    eh[  34] += xint[  3]*yint[  9]*zint[ 12];
    eh[  34] += xint[  4]*yint[ 10]*zint[ 13];
    eh[  34] += xint[  5]*yint[ 11]*zint[ 14];
    eh[  35] += xint[  0]*yint[ 12]*zint[ 12];
    eh[  35] += xint[  1]*yint[ 13]*zint[ 13];
    eh[  35] += xint[  2]*yint[ 14]*zint[ 14];
    // (FS|DS)
    eh[  36] += xint[ 33]*yint[  0]*zint[  0];
    eh[  36] += xint[ 34]*yint[  1]*zint[  1];
    eh[  36] += xint[ 35]*yint[  2]*zint[  2];
    eh[  37] += xint[ 27]*yint[  6]*zint[  0];
    eh[  37] += xint[ 28]*yint[  7]*zint[  1];
    eh[  37] += xint[ 29]*yint[  8]*zint[  2];
    eh[  38] += xint[ 27]*yint[  0]*zint[  6];
    eh[  38] += xint[ 28]*yint[  1]*zint[  7];
    eh[  38] += xint[ 29]*yint[  2]*zint[  8];
    eh[  39] += xint[ 30]*yint[  3]*zint[  0];
    eh[  39] += xint[ 31]*yint[  4]*zint[  1];
    eh[  39] += xint[ 32]*yint[  5]*zint[  2];
    eh[  40] += xint[ 30]*yint[  0]*zint[  3];
    eh[  40] += xint[ 31]*yint[  1]*zint[  4];
    eh[  40] += xint[ 32]*yint[  2]*zint[  5];
    eh[  41] += xint[ 27]*yint[  3]*zint[  3];
    eh[  41] += xint[ 28]*yint[  4]*zint[  4];
    eh[  41] += xint[ 29]*yint[  5]*zint[  5];
    eh[  42] += xint[  6]*yint[ 27]*zint[  0];
    eh[  42] += xint[  7]*yint[ 28]*zint[  1];
    eh[  42] += xint[  8]*yint[ 29]*zint[  2];
    eh[  43] += xint[  0]*yint[ 33]*zint[  0];
    eh[  43] += xint[  1]*yint[ 34]*zint[  1];
    eh[  43] += xint[  2]*yint[ 35]*zint[  2];
    eh[  44] += xint[  0]*yint[ 27]*zint[  6];
    eh[  44] += xint[  1]*yint[ 28]*zint[  7];
    eh[  44] += xint[  2]*yint[ 29]*zint[  8];
    eh[  45] += xint[  3]*yint[ 30]*zint[  0];
    eh[  45] += xint[  4]*yint[ 31]*zint[  1];
    eh[  45] += xint[  5]*yint[ 32]*zint[  2];
    eh[  46] += xint[  3]*yint[ 27]*zint[  3];
    eh[  46] += xint[  4]*yint[ 28]*zint[  4];
    eh[  46] += xint[  5]*yint[ 29]*zint[  5];
    eh[  47] += xint[  0]*yint[ 30]*zint[  3];
    eh[  47] += xint[  1]*yint[ 31]*zint[  4];
    eh[  47] += xint[  2]*yint[ 32]*zint[  5];
    eh[  48] += xint[  6]*yint[  0]*zint[ 27];
    eh[  48] += xint[  7]*yint[  1]*zint[ 28];
    eh[  48] += xint[  8]*yint[  2]*zint[ 29];
    eh[  49] += xint[  0]*yint[  6]*zint[ 27];
    eh[  49] += xint[  1]*yint[  7]*zint[ 28];
    eh[  49] += xint[  2]*yint[  8]*zint[ 29];
    eh[  50] += xint[  0]*yint[  0]*zint[ 33];
    eh[  50] += xint[  1]*yint[  1]*zint[ 34];
    eh[  50] += xint[  2]*yint[  2]*zint[ 35];
    eh[  51] += xint[  3]*yint[  3]*zint[ 27];
    eh[  51] += xint[  4]*yint[  4]*zint[ 28];
    eh[  51] += xint[  5]*yint[  5]*zint[ 29];
    eh[  52] += xint[  3]*yint[  0]*zint[ 30];
    eh[  52] += xint[  4]*yint[  1]*zint[ 31];
    eh[  52] += xint[  5]*yint[  2]*zint[ 32];
    eh[  53] += xint[  0]*yint[  3]*zint[ 30];
    eh[  53] += xint[  1]*yint[  4]*zint[ 31];
    eh[  53] += xint[  2]*yint[  5]*zint[ 32];
    eh[  54] += xint[ 24]*yint[  9]*zint[  0];
    eh[  54] += xint[ 25]*yint[ 10]*zint[  1];
    eh[  54] += xint[ 26]*yint[ 11]*zint[  2];
    eh[  55] += xint[ 18]*yint[ 15]*zint[  0];
    eh[  55] += xint[ 19]*yint[ 16]*zint[  1];
    eh[  55] += xint[ 20]*yint[ 17]*zint[  2];
    eh[  56] += xint[ 18]*yint[  9]*zint[  6];
    eh[  56] += xint[ 19]*yint[ 10]*zint[  7];
    eh[  56] += xint[ 20]*yint[ 11]*zint[  8];
    eh[  57] += xint[ 21]*yint[ 12]*zint[  0];
    eh[  57] += xint[ 22]*yint[ 13]*zint[  1];
    eh[  57] += xint[ 23]*yint[ 14]*zint[  2];
    eh[  58] += xint[ 21]*yint[  9]*zint[  3];
    eh[  58] += xint[ 22]*yint[ 10]*zint[  4];
    eh[  58] += xint[ 23]*yint[ 11]*zint[  5];
    eh[  59] += xint[ 18]*yint[ 12]*zint[  3];
    eh[  59] += xint[ 19]*yint[ 13]*zint[  4];
    eh[  59] += xint[ 20]*yint[ 14]*zint[  5];
    eh[  60] += xint[ 24]*yint[  0]*zint[  9];
    eh[  60] += xint[ 25]*yint[  1]*zint[ 10];
    eh[  60] += xint[ 26]*yint[  2]*zint[ 11];
    eh[  61] += xint[ 18]*yint[  6]*zint[  9];
    eh[  61] += xint[ 19]*yint[  7]*zint[ 10];
    eh[  61] += xint[ 20]*yint[  8]*zint[ 11];
    eh[  62] += xint[ 18]*yint[  0]*zint[ 15];
    eh[  62] += xint[ 19]*yint[  1]*zint[ 16];
    eh[  62] += xint[ 20]*yint[  2]*zint[ 17];
    eh[  63] += xint[ 21]*yint[  3]*zint[  9];
    eh[  63] += xint[ 22]*yint[  4]*zint[ 10];
    eh[  63] += xint[ 23]*yint[  5]*zint[ 11];
    eh[  64] += xint[ 21]*yint[  0]*zint[ 12];
    eh[  64] += xint[ 22]*yint[  1]*zint[ 13];
    eh[  64] += xint[ 23]*yint[  2]*zint[ 14];
    eh[  65] += xint[ 18]*yint[  3]*zint[ 12];
    eh[  65] += xint[ 19]*yint[  4]*zint[ 13];
    eh[  65] += xint[ 20]*yint[  5]*zint[ 14];
    eh[  66] += xint[ 15]*yint[ 18]*zint[  0];
    eh[  66] += xint[ 16]*yint[ 19]*zint[  1];
    eh[  66] += xint[ 17]*yint[ 20]*zint[  2];
    eh[  67] += xint[  9]*yint[ 24]*zint[  0];
    eh[  67] += xint[ 10]*yint[ 25]*zint[  1];
    eh[  67] += xint[ 11]*yint[ 26]*zint[  2];
    eh[  68] += xint[  9]*yint[ 18]*zint[  6];
    eh[  68] += xint[ 10]*yint[ 19]*zint[  7];
    eh[  68] += xint[ 11]*yint[ 20]*zint[  8];
    eh[  69] += xint[ 12]*yint[ 21]*zint[  0];
    eh[  69] += xint[ 13]*yint[ 22]*zint[  1];
    eh[  69] += xint[ 14]*yint[ 23]*zint[  2];
    eh[  70] += xint[ 12]*yint[ 18]*zint[  3];
    eh[  70] += xint[ 13]*yint[ 19]*zint[  4];
    eh[  70] += xint[ 14]*yint[ 20]*zint[  5];
    eh[  71] += xint[  9]*yint[ 21]*zint[  3];
    eh[  71] += xint[ 10]*yint[ 22]*zint[  4];
    eh[  71] += xint[ 11]*yint[ 23]*zint[  5];
    eh[  72] += xint[ 15]*yint[  0]*zint[ 18];
    eh[  72] += xint[ 16]*yint[  1]*zint[ 19];
    eh[  72] += xint[ 17]*yint[  2]*zint[ 20];
    eh[  73] += xint[  9]*yint[  6]*zint[ 18];
    eh[  73] += xint[ 10]*yint[  7]*zint[ 19];
    eh[  73] += xint[ 11]*yint[  8]*zint[ 20];
    eh[  74] += xint[  9]*yint[  0]*zint[ 24];
    eh[  74] += xint[ 10]*yint[  1]*zint[ 25];
    eh[  74] += xint[ 11]*yint[  2]*zint[ 26];
    eh[  75] += xint[ 12]*yint[  3]*zint[ 18];
    eh[  75] += xint[ 13]*yint[  4]*zint[ 19];
    eh[  75] += xint[ 14]*yint[  5]*zint[ 20];
    eh[  76] += xint[ 12]*yint[  0]*zint[ 21];
    eh[  76] += xint[ 13]*yint[  1]*zint[ 22];
    eh[  76] += xint[ 14]*yint[  2]*zint[ 23];
    eh[  77] += xint[  9]*yint[  3]*zint[ 21];
    eh[  77] += xint[ 10]*yint[  4]*zint[ 22];
    eh[  77] += xint[ 11]*yint[  5]*zint[ 23];
    eh[  78] += xint[ 15]*yint[  9]*zint[  9];
    eh[  78] += xint[ 16]*yint[ 10]*zint[ 10];
    eh[  78] += xint[ 17]*yint[ 11]*zint[ 11];
    eh[  79] += xint[  9]*yint[ 15]*zint[  9];
    eh[  79] += xint[ 10]*yint[ 16]*zint[ 10];
    eh[  79] += xint[ 11]*yint[ 17]*zint[ 11];
    eh[  80] += xint[  9]*yint[  9]*zint[ 15];
    eh[  80] += xint[ 10]*yint[ 10]*zint[ 16];
    eh[  80] += xint[ 11]*yint[ 11]*zint[ 17];
    eh[  81] += xint[ 12]*yint[ 12]*zint[  9];
    eh[  81] += xint[ 13]*yint[ 13]*zint[ 10];
    eh[  81] += xint[ 14]*yint[ 14]*zint[ 11];
    eh[  82] += xint[ 12]*yint[  9]*zint[ 12];
    eh[  82] += xint[ 13]*yint[ 10]*zint[ 13];
    eh[  82] += xint[ 14]*yint[ 11]*zint[ 14];
    eh[  83] += xint[  9]*yint[ 12]*zint[ 12];
    eh[  83] += xint[ 10]*yint[ 13]*zint[ 13];
    eh[  83] += xint[ 11]*yint[ 14]*zint[ 14];
    eh[  84] += xint[  6]*yint[ 18]*zint[  9];
    eh[  84] += xint[  7]*yint[ 19]*zint[ 10];
    eh[  84] += xint[  8]*yint[ 20]*zint[ 11];
    eh[  85] += xint[  0]*yint[ 24]*zint[  9];
    eh[  85] += xint[  1]*yint[ 25]*zint[ 10];
    eh[  85] += xint[  2]*yint[ 26]*zint[ 11];
    eh[  86] += xint[  0]*yint[ 18]*zint[ 15];
    eh[  86] += xint[  1]*yint[ 19]*zint[ 16];
    eh[  86] += xint[  2]*yint[ 20]*zint[ 17];
    eh[  87] += xint[  3]*yint[ 21]*zint[  9];
    eh[  87] += xint[  4]*yint[ 22]*zint[ 10];
    eh[  87] += xint[  5]*yint[ 23]*zint[ 11];
    eh[  88] += xint[  3]*yint[ 18]*zint[ 12];
    eh[  88] += xint[  4]*yint[ 19]*zint[ 13];
    eh[  88] += xint[  5]*yint[ 20]*zint[ 14];
    eh[  89] += xint[  0]*yint[ 21]*zint[ 12];
    eh[  89] += xint[  1]*yint[ 22]*zint[ 13];
    eh[  89] += xint[  2]*yint[ 23]*zint[ 14];
    eh[  90] += xint[  6]*yint[  9]*zint[ 18];
    eh[  90] += xint[  7]*yint[ 10]*zint[ 19];
    eh[  90] += xint[  8]*yint[ 11]*zint[ 20];
    eh[  91] += xint[  0]*yint[ 15]*zint[ 18];
    eh[  91] += xint[  1]*yint[ 16]*zint[ 19];
    eh[  91] += xint[  2]*yint[ 17]*zint[ 20];
    eh[  92] += xint[  0]*yint[  9]*zint[ 24];
    eh[  92] += xint[  1]*yint[ 10]*zint[ 25];
    eh[  92] += xint[  2]*yint[ 11]*zint[ 26];
    eh[  93] += xint[  3]*yint[ 12]*zint[ 18];
    eh[  93] += xint[  4]*yint[ 13]*zint[ 19];
    eh[  93] += xint[  5]*yint[ 14]*zint[ 20];
    eh[  94] += xint[  3]*yint[  9]*zint[ 21];
    eh[  94] += xint[  4]*yint[ 10]*zint[ 22];
    eh[  94] += xint[  5]*yint[ 11]*zint[ 23];
    eh[  95] += xint[  0]*yint[ 12]*zint[ 21];
    eh[  95] += xint[  1]*yint[ 13]*zint[ 22];
    eh[  95] += xint[  2]*yint[ 14]*zint[ 23];
}

void ofmo_twoint_core_rys_dpds( const int mythread,
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
    ofmo_hrr_clear_dpds( eh );
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
            ofmo_xyzint_dpds( F00, B00, B10, B01, C00, CP00,
                     xint, yint, zint );
            ofmo_form_dpds( xint, yint, zint, eh );
        }
    }
    ofmo_hrr_calc_dpds( eh, BA, DC );
    ofmo_hrr_coef_dpds( eh, DINT );
}

int ofmo_twoint_rys_dpds(
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
    max_nzeri = ebuf_max_nzeri - 6*3*6*1;
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
            ofmo_twoint_core_rys_dpds( mythread,
                    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
                    &psp_xiza[ijps0], BA,
                    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
                    &psp_xiza[klps0], DC,   AC,      DINTEG );
            ipat=((Lab != Lcd) || (ics==kcs && jcs>lcs) ? true : false);
            for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
                I2 = (iao*iao+iao)>>1;
                for ( j=0, jao=jao0; j<3; j++, jao++ ) {
                    if ( jao>iao ) { ix+=6*1; continue; }
                    IJ = I2 + jao;
                    coe0 = ( iao==jao ? HALF : ONE );
                    for ( k=0, kao=kao0; k<6; k++, kao++ ) {
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

int ofmo_twoint_direct_rys_dpds(
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
    max_nzeri -= 6*3*6*1;
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
            ofmo_twoint_core_rys_dpds( mythread,
                    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
                    &psp_xiza[ijps0], BA,
                    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
                    &psp_xiza[klps0], DC,   AC,      DINTEG );
            ipat=((Lab != Lcd) || (ics==kcs && jcs>lcs) ? true : false);
            for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
                I2 = (iao*iao+iao)>>1;
                for ( j=0, jao=jao0; j<3; j++, jao++ ) {
                    if ( jao>iao ) { ix+=6*1; continue; }
                    IJ = I2 + jao;
                    coe0 = ( iao==jao ? HALF : ONE );
                    for ( k=0, kao=kao0; k<6; k++, kao++ ) {
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
