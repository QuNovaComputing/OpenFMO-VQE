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

static void ofmo_hrr_clear_ddps( double *eh ) {
    int i;
    // (DS|PS)
    for ( i=0; i<(0+18); i++ ) eh[i] = 0.e0;
    // (FS|PS)
    for ( i=18; i<(18+30); i++ ) eh[i] = 0.e0;
    // (GS|PS)
    for ( i=48; i<(48+45); i++ ) eh[i] = 0.e0;
}

static void ofmo_hrr_coef_ddps(
        double *eh, double *DINT ) {
    int i, j, k, l, iao, jao, kao, lao, ix;
    double coef_a, coef_ab, coef_abc;
    double *th;
    th = &eh[237];
    ix = 0;
    for ( i=0, iao=4; i<6; i++, iao++ ) {
        coef_a = DFACT[iao];
        for ( j=0, jao=4; j<6; j++, jao++ ) {
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

static void ofmo_hrr_calc_ddps( double *eh,
        const double BA[3], const double DC[3] ) {
    // HRR for (XX|XS)-type integral (center AB)
    // (DP,PS)
    eh[  93] = eh[  18] - BA[0]*eh[   0];
    eh[  94] = eh[  19] - BA[0]*eh[   1];
    eh[  95] = eh[  20] - BA[0]*eh[   2];
    eh[  96] = eh[  27] - BA[1]*eh[   0];
    eh[  97] = eh[  28] - BA[1]*eh[   1];
    eh[  98] = eh[  29] - BA[1]*eh[   2];
    eh[  99] = eh[  30] - BA[2]*eh[   0];
    eh[ 100] = eh[  31] - BA[2]*eh[   1];
    eh[ 101] = eh[  32] - BA[2]*eh[   2];
    eh[ 102] = eh[  33] - BA[0]*eh[   3];
    eh[ 103] = eh[  34] - BA[0]*eh[   4];
    eh[ 104] = eh[  35] - BA[0]*eh[   5];
    eh[ 105] = eh[  21] - BA[1]*eh[   3];
    eh[ 106] = eh[  22] - BA[1]*eh[   4];
    eh[ 107] = eh[  23] - BA[1]*eh[   5];
    eh[ 108] = eh[  42] - BA[2]*eh[   3];
    eh[ 109] = eh[  43] - BA[2]*eh[   4];
    eh[ 110] = eh[  44] - BA[2]*eh[   5];
    eh[ 111] = eh[  36] - BA[0]*eh[   6];
    eh[ 112] = eh[  37] - BA[0]*eh[   7];
    eh[ 113] = eh[  38] - BA[0]*eh[   8];
    eh[ 114] = eh[  45] - BA[1]*eh[   6];
    eh[ 115] = eh[  46] - BA[1]*eh[   7];
    eh[ 116] = eh[  47] - BA[1]*eh[   8];
    eh[ 117] = eh[  24] - BA[2]*eh[   6];
    eh[ 118] = eh[  25] - BA[2]*eh[   7];
    eh[ 119] = eh[  26] - BA[2]*eh[   8];
    eh[ 120] = eh[  27] - BA[0]*eh[   9];
    eh[ 121] = eh[  28] - BA[0]*eh[  10];
    eh[ 122] = eh[  29] - BA[0]*eh[  11];
    eh[ 123] = eh[  33] - BA[1]*eh[   9];
    eh[ 124] = eh[  34] - BA[1]*eh[  10];
    eh[ 125] = eh[  35] - BA[1]*eh[  11];
    eh[ 126] = eh[  39] - BA[2]*eh[   9];
    eh[ 127] = eh[  40] - BA[2]*eh[  10];
    eh[ 128] = eh[  41] - BA[2]*eh[  11];
    eh[ 129] = eh[  30] - BA[0]*eh[  12];
    eh[ 130] = eh[  31] - BA[0]*eh[  13];
    eh[ 131] = eh[  32] - BA[0]*eh[  14];
    eh[ 132] = eh[  39] - BA[1]*eh[  12];
    eh[ 133] = eh[  40] - BA[1]*eh[  13];
    eh[ 134] = eh[  41] - BA[1]*eh[  14];
    eh[ 135] = eh[  36] - BA[2]*eh[  12];
    eh[ 136] = eh[  37] - BA[2]*eh[  13];
    eh[ 137] = eh[  38] - BA[2]*eh[  14];
    eh[ 138] = eh[  39] - BA[0]*eh[  15];
    eh[ 139] = eh[  40] - BA[0]*eh[  16];
    eh[ 140] = eh[  41] - BA[0]*eh[  17];
    eh[ 141] = eh[  42] - BA[1]*eh[  15];
    eh[ 142] = eh[  43] - BA[1]*eh[  16];
    eh[ 143] = eh[  44] - BA[1]*eh[  17];
    eh[ 144] = eh[  45] - BA[2]*eh[  15];
    eh[ 145] = eh[  46] - BA[2]*eh[  16];
    eh[ 146] = eh[  47] - BA[2]*eh[  17];
    // (FP,PS)
    eh[ 147] = eh[  48] - BA[0]*eh[  18];
    eh[ 148] = eh[  49] - BA[0]*eh[  19];
    eh[ 149] = eh[  50] - BA[0]*eh[  20];
    eh[ 150] = eh[  57] - BA[1]*eh[  18];
    eh[ 151] = eh[  58] - BA[1]*eh[  19];
    eh[ 152] = eh[  59] - BA[1]*eh[  20];
    eh[ 153] = eh[  60] - BA[2]*eh[  18];
    eh[ 154] = eh[  61] - BA[2]*eh[  19];
    eh[ 155] = eh[  62] - BA[2]*eh[  20];
    eh[ 156] = eh[  72] - BA[0]*eh[  21];
    eh[ 157] = eh[  73] - BA[0]*eh[  22];
    eh[ 158] = eh[  74] - BA[0]*eh[  23];
    eh[ 159] = eh[  51] - BA[1]*eh[  21];
    eh[ 160] = eh[  52] - BA[1]*eh[  22];
    eh[ 161] = eh[  53] - BA[1]*eh[  23];
    eh[ 162] = eh[  84] - BA[2]*eh[  21];
    eh[ 163] = eh[  85] - BA[2]*eh[  22];
    eh[ 164] = eh[  86] - BA[2]*eh[  23];
    eh[ 165] = eh[  75] - BA[0]*eh[  24];
    eh[ 166] = eh[  76] - BA[0]*eh[  25];
    eh[ 167] = eh[  77] - BA[0]*eh[  26];
    eh[ 168] = eh[  87] - BA[1]*eh[  24];
    eh[ 169] = eh[  88] - BA[1]*eh[  25];
    eh[ 170] = eh[  89] - BA[1]*eh[  26];
    eh[ 171] = eh[  54] - BA[2]*eh[  24];
    eh[ 172] = eh[  55] - BA[2]*eh[  25];
    eh[ 173] = eh[  56] - BA[2]*eh[  26];
    eh[ 174] = eh[  57] - BA[0]*eh[  27];
    eh[ 175] = eh[  58] - BA[0]*eh[  28];
    eh[ 176] = eh[  59] - BA[0]*eh[  29];
    eh[ 177] = eh[  63] - BA[1]*eh[  27];
    eh[ 178] = eh[  64] - BA[1]*eh[  28];
    eh[ 179] = eh[  65] - BA[1]*eh[  29];
    eh[ 180] = eh[  69] - BA[2]*eh[  27];
    eh[ 181] = eh[  70] - BA[2]*eh[  28];
    eh[ 182] = eh[  71] - BA[2]*eh[  29];
    eh[ 183] = eh[  60] - BA[0]*eh[  30];
    eh[ 184] = eh[  61] - BA[0]*eh[  31];
    eh[ 185] = eh[  62] - BA[0]*eh[  32];
    eh[ 186] = eh[  69] - BA[1]*eh[  30];
    eh[ 187] = eh[  70] - BA[1]*eh[  31];
    eh[ 188] = eh[  71] - BA[1]*eh[  32];
    eh[ 189] = eh[  66] - BA[2]*eh[  30];
    eh[ 190] = eh[  67] - BA[2]*eh[  31];
    eh[ 191] = eh[  68] - BA[2]*eh[  32];
    eh[ 192] = eh[  63] - BA[0]*eh[  33];
    eh[ 193] = eh[  64] - BA[0]*eh[  34];
    eh[ 194] = eh[  65] - BA[0]*eh[  35];
    eh[ 195] = eh[  72] - BA[1]*eh[  33];
    eh[ 196] = eh[  73] - BA[1]*eh[  34];
    eh[ 197] = eh[  74] - BA[1]*eh[  35];
    eh[ 198] = eh[  78] - BA[2]*eh[  33];
    eh[ 199] = eh[  79] - BA[2]*eh[  34];
    eh[ 200] = eh[  80] - BA[2]*eh[  35];
    eh[ 201] = eh[  66] - BA[0]*eh[  36];
    eh[ 202] = eh[  67] - BA[0]*eh[  37];
    eh[ 203] = eh[  68] - BA[0]*eh[  38];
    eh[ 204] = eh[  81] - BA[1]*eh[  36];
    eh[ 205] = eh[  82] - BA[1]*eh[  37];
    eh[ 206] = eh[  83] - BA[1]*eh[  38];
    eh[ 207] = eh[  75] - BA[2]*eh[  36];
    eh[ 208] = eh[  76] - BA[2]*eh[  37];
    eh[ 209] = eh[  77] - BA[2]*eh[  38];
    eh[ 210] = eh[  69] - BA[0]*eh[  39];
    eh[ 211] = eh[  70] - BA[0]*eh[  40];
    eh[ 212] = eh[  71] - BA[0]*eh[  41];
    eh[ 213] = eh[  78] - BA[1]*eh[  39];
    eh[ 214] = eh[  79] - BA[1]*eh[  40];
    eh[ 215] = eh[  80] - BA[1]*eh[  41];
    eh[ 216] = eh[  81] - BA[2]*eh[  39];
    eh[ 217] = eh[  82] - BA[2]*eh[  40];
    eh[ 218] = eh[  83] - BA[2]*eh[  41];
    eh[ 219] = eh[  78] - BA[0]*eh[  42];
    eh[ 220] = eh[  79] - BA[0]*eh[  43];
    eh[ 221] = eh[  80] - BA[0]*eh[  44];
    eh[ 222] = eh[  84] - BA[1]*eh[  42];
    eh[ 223] = eh[  85] - BA[1]*eh[  43];
    eh[ 224] = eh[  86] - BA[1]*eh[  44];
    eh[ 225] = eh[  90] - BA[2]*eh[  42];
    eh[ 226] = eh[  91] - BA[2]*eh[  43];
    eh[ 227] = eh[  92] - BA[2]*eh[  44];
    eh[ 228] = eh[  81] - BA[0]*eh[  45];
    eh[ 229] = eh[  82] - BA[0]*eh[  46];
    eh[ 230] = eh[  83] - BA[0]*eh[  47];
    eh[ 231] = eh[  90] - BA[1]*eh[  45];
    eh[ 232] = eh[  91] - BA[1]*eh[  46];
    eh[ 233] = eh[  92] - BA[1]*eh[  47];
    eh[ 234] = eh[  87] - BA[2]*eh[  45];
    eh[ 235] = eh[  88] - BA[2]*eh[  46];
    eh[ 236] = eh[  89] - BA[2]*eh[  47];
    // (DD,PS)
    eh[ 237] = eh[ 147] - BA[0]*eh[  93];
    eh[ 238] = eh[ 148] - BA[0]*eh[  94];
    eh[ 239] = eh[ 149] - BA[0]*eh[  95];
    eh[ 240] = eh[ 177] - BA[1]*eh[  96];
    eh[ 241] = eh[ 178] - BA[1]*eh[  97];
    eh[ 242] = eh[ 179] - BA[1]*eh[  98];
    eh[ 243] = eh[ 189] - BA[2]*eh[  99];
    eh[ 244] = eh[ 190] - BA[2]*eh[ 100];
    eh[ 245] = eh[ 191] - BA[2]*eh[ 101];
    eh[ 246] = eh[ 150] - BA[0]*eh[  96];
    eh[ 247] = eh[ 151] - BA[0]*eh[  97];
    eh[ 248] = eh[ 152] - BA[0]*eh[  98];
    eh[ 249] = eh[ 153] - BA[0]*eh[  99];
    eh[ 250] = eh[ 154] - BA[0]*eh[ 100];
    eh[ 251] = eh[ 155] - BA[0]*eh[ 101];
    eh[ 252] = eh[ 180] - BA[1]*eh[  99];
    eh[ 253] = eh[ 181] - BA[1]*eh[ 100];
    eh[ 254] = eh[ 182] - BA[1]*eh[ 101];
    eh[ 255] = eh[ 192] - BA[0]*eh[ 102];
    eh[ 256] = eh[ 193] - BA[0]*eh[ 103];
    eh[ 257] = eh[ 194] - BA[0]*eh[ 104];
    eh[ 258] = eh[ 159] - BA[1]*eh[ 105];
    eh[ 259] = eh[ 160] - BA[1]*eh[ 106];
    eh[ 260] = eh[ 161] - BA[1]*eh[ 107];
    eh[ 261] = eh[ 225] - BA[2]*eh[ 108];
    eh[ 262] = eh[ 226] - BA[2]*eh[ 109];
    eh[ 263] = eh[ 227] - BA[2]*eh[ 110];
    eh[ 264] = eh[ 195] - BA[0]*eh[ 105];
    eh[ 265] = eh[ 196] - BA[0]*eh[ 106];
    eh[ 266] = eh[ 197] - BA[0]*eh[ 107];
    eh[ 267] = eh[ 198] - BA[0]*eh[ 108];
    eh[ 268] = eh[ 199] - BA[0]*eh[ 109];
    eh[ 269] = eh[ 200] - BA[0]*eh[ 110];
    eh[ 270] = eh[ 162] - BA[1]*eh[ 108];
    eh[ 271] = eh[ 163] - BA[1]*eh[ 109];
    eh[ 272] = eh[ 164] - BA[1]*eh[ 110];
    eh[ 273] = eh[ 201] - BA[0]*eh[ 111];
    eh[ 274] = eh[ 202] - BA[0]*eh[ 112];
    eh[ 275] = eh[ 203] - BA[0]*eh[ 113];
    eh[ 276] = eh[ 231] - BA[1]*eh[ 114];
    eh[ 277] = eh[ 232] - BA[1]*eh[ 115];
    eh[ 278] = eh[ 233] - BA[1]*eh[ 116];
    eh[ 279] = eh[ 171] - BA[2]*eh[ 117];
    eh[ 280] = eh[ 172] - BA[2]*eh[ 118];
    eh[ 281] = eh[ 173] - BA[2]*eh[ 119];
    eh[ 282] = eh[ 204] - BA[0]*eh[ 114];
    eh[ 283] = eh[ 205] - BA[0]*eh[ 115];
    eh[ 284] = eh[ 206] - BA[0]*eh[ 116];
    eh[ 285] = eh[ 207] - BA[0]*eh[ 117];
    eh[ 286] = eh[ 208] - BA[0]*eh[ 118];
    eh[ 287] = eh[ 209] - BA[0]*eh[ 119];
    eh[ 288] = eh[ 234] - BA[1]*eh[ 117];
    eh[ 289] = eh[ 235] - BA[1]*eh[ 118];
    eh[ 290] = eh[ 236] - BA[1]*eh[ 119];
    eh[ 291] = eh[ 174] - BA[0]*eh[ 120];
    eh[ 292] = eh[ 175] - BA[0]*eh[ 121];
    eh[ 293] = eh[ 176] - BA[0]*eh[ 122];
    eh[ 294] = eh[ 195] - BA[1]*eh[ 123];
    eh[ 295] = eh[ 196] - BA[1]*eh[ 124];
    eh[ 296] = eh[ 197] - BA[1]*eh[ 125];
    eh[ 297] = eh[ 216] - BA[2]*eh[ 126];
    eh[ 298] = eh[ 217] - BA[2]*eh[ 127];
    eh[ 299] = eh[ 218] - BA[2]*eh[ 128];
    eh[ 300] = eh[ 177] - BA[0]*eh[ 123];
    eh[ 301] = eh[ 178] - BA[0]*eh[ 124];
    eh[ 302] = eh[ 179] - BA[0]*eh[ 125];
    eh[ 303] = eh[ 180] - BA[0]*eh[ 126];
    eh[ 304] = eh[ 181] - BA[0]*eh[ 127];
    eh[ 305] = eh[ 182] - BA[0]*eh[ 128];
    eh[ 306] = eh[ 198] - BA[1]*eh[ 126];
    eh[ 307] = eh[ 199] - BA[1]*eh[ 127];
    eh[ 308] = eh[ 200] - BA[1]*eh[ 128];
    eh[ 309] = eh[ 183] - BA[0]*eh[ 129];
    eh[ 310] = eh[ 184] - BA[0]*eh[ 130];
    eh[ 311] = eh[ 185] - BA[0]*eh[ 131];
    eh[ 312] = eh[ 213] - BA[1]*eh[ 132];
    eh[ 313] = eh[ 214] - BA[1]*eh[ 133];
    eh[ 314] = eh[ 215] - BA[1]*eh[ 134];
    eh[ 315] = eh[ 207] - BA[2]*eh[ 135];
    eh[ 316] = eh[ 208] - BA[2]*eh[ 136];
    eh[ 317] = eh[ 209] - BA[2]*eh[ 137];
    eh[ 318] = eh[ 186] - BA[0]*eh[ 132];
    eh[ 319] = eh[ 187] - BA[0]*eh[ 133];
    eh[ 320] = eh[ 188] - BA[0]*eh[ 134];
    eh[ 321] = eh[ 189] - BA[0]*eh[ 135];
    eh[ 322] = eh[ 190] - BA[0]*eh[ 136];
    eh[ 323] = eh[ 191] - BA[0]*eh[ 137];
    eh[ 324] = eh[ 216] - BA[1]*eh[ 135];
    eh[ 325] = eh[ 217] - BA[1]*eh[ 136];
    eh[ 326] = eh[ 218] - BA[1]*eh[ 137];
    eh[ 327] = eh[ 210] - BA[0]*eh[ 138];
    eh[ 328] = eh[ 211] - BA[0]*eh[ 139];
    eh[ 329] = eh[ 212] - BA[0]*eh[ 140];
    eh[ 330] = eh[ 222] - BA[1]*eh[ 141];
    eh[ 331] = eh[ 223] - BA[1]*eh[ 142];
    eh[ 332] = eh[ 224] - BA[1]*eh[ 143];
    eh[ 333] = eh[ 234] - BA[2]*eh[ 144];
    eh[ 334] = eh[ 235] - BA[2]*eh[ 145];
    eh[ 335] = eh[ 236] - BA[2]*eh[ 146];
    eh[ 336] = eh[ 213] - BA[0]*eh[ 141];
    eh[ 337] = eh[ 214] - BA[0]*eh[ 142];
    eh[ 338] = eh[ 215] - BA[0]*eh[ 143];
    eh[ 339] = eh[ 216] - BA[0]*eh[ 144];
    eh[ 340] = eh[ 217] - BA[0]*eh[ 145];
    eh[ 341] = eh[ 218] - BA[0]*eh[ 146];
    eh[ 342] = eh[ 225] - BA[1]*eh[ 144];
    eh[ 343] = eh[ 226] - BA[1]*eh[ 145];
    eh[ 344] = eh[ 227] - BA[1]*eh[ 146];
}

static void ofmo_xyzint_ddps(
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
    // (4,0)
    C10[0] += B10[0];
    xint[ 24]=C00[ 0]*xint[ 18]+C10[0]*xint[ 12];
    yint[ 24]=C00[ 1]*yint[ 18]+C10[0]*yint[ 12];
    zint[ 24]=C00[ 2]*zint[ 18]+C10[0]*zint[ 12];
    C10[1] += B10[1];
    xint[ 25]=C00[ 3]*xint[ 19]+C10[1]*xint[ 13];
    yint[ 25]=C00[ 4]*yint[ 19]+C10[1]*yint[ 13];
    zint[ 25]=C00[ 5]*zint[ 19]+C10[1]*zint[ 13];
    C10[2] += B10[2];
    xint[ 26]=C00[ 6]*xint[ 20]+C10[2]*xint[ 14];
    yint[ 26]=C00[ 7]*yint[ 20]+C10[2]*yint[ 14];
    zint[ 26]=C00[ 8]*zint[ 20]+C10[2]*zint[ 14];
    // (4,1)
    CP10[0] += B00[0];
    xint[ 27]=CP00[ 0]*xint[ 24]+CP10[0]*xint[ 18];
    yint[ 27]=CP00[ 1]*yint[ 24]+CP10[0]*yint[ 18];
    zint[ 27]=CP00[ 2]*zint[ 24]+CP10[0]*zint[ 18];
    CP10[1] += B00[1];
    xint[ 28]=CP00[ 3]*xint[ 25]+CP10[1]*xint[ 19];
    yint[ 28]=CP00[ 4]*yint[ 25]+CP10[1]*yint[ 19];
    zint[ 28]=CP00[ 5]*zint[ 25]+CP10[1]*zint[ 19];
    CP10[2] += B00[2];
    xint[ 29]=CP00[ 6]*xint[ 26]+CP10[2]*xint[ 20];
    yint[ 29]=CP00[ 7]*yint[ 26]+CP10[2]*yint[ 20];
    zint[ 29]=CP00[ 8]*zint[ 26]+CP10[2]*zint[ 20];
}

static void ofmo_form_ddps(
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
    // (GS|PS)
    eh[  48] += xint[ 27]*yint[  0]*zint[  0];
    eh[  48] += xint[ 28]*yint[  1]*zint[  1];
    eh[  48] += xint[ 29]*yint[  2]*zint[  2];
    eh[  49] += xint[ 24]*yint[  3]*zint[  0];
    eh[  49] += xint[ 25]*yint[  4]*zint[  1];
    eh[  49] += xint[ 26]*yint[  5]*zint[  2];
    eh[  50] += xint[ 24]*yint[  0]*zint[  3];
    eh[  50] += xint[ 25]*yint[  1]*zint[  4];
    eh[  50] += xint[ 26]*yint[  2]*zint[  5];
    eh[  51] += xint[  3]*yint[ 24]*zint[  0];
    eh[  51] += xint[  4]*yint[ 25]*zint[  1];
    eh[  51] += xint[  5]*yint[ 26]*zint[  2];
    eh[  52] += xint[  0]*yint[ 27]*zint[  0];
    eh[  52] += xint[  1]*yint[ 28]*zint[  1];
    eh[  52] += xint[  2]*yint[ 29]*zint[  2];
    eh[  53] += xint[  0]*yint[ 24]*zint[  3];
    eh[  53] += xint[  1]*yint[ 25]*zint[  4];
    eh[  53] += xint[  2]*yint[ 26]*zint[  5];
    eh[  54] += xint[  3]*yint[  0]*zint[ 24];
    eh[  54] += xint[  4]*yint[  1]*zint[ 25];
    eh[  54] += xint[  5]*yint[  2]*zint[ 26];
    eh[  55] += xint[  0]*yint[  3]*zint[ 24];
    eh[  55] += xint[  1]*yint[  4]*zint[ 25];
    eh[  55] += xint[  2]*yint[  5]*zint[ 26];
    eh[  56] += xint[  0]*yint[  0]*zint[ 27];
    eh[  56] += xint[  1]*yint[  1]*zint[ 28];
    eh[  56] += xint[  2]*yint[  2]*zint[ 29];
    eh[  57] += xint[ 21]*yint[  6]*zint[  0];
    eh[  57] += xint[ 22]*yint[  7]*zint[  1];
    eh[  57] += xint[ 23]*yint[  8]*zint[  2];
    eh[  58] += xint[ 18]*yint[  9]*zint[  0];
    eh[  58] += xint[ 19]*yint[ 10]*zint[  1];
    eh[  58] += xint[ 20]*yint[ 11]*zint[  2];
    eh[  59] += xint[ 18]*yint[  6]*zint[  3];
    eh[  59] += xint[ 19]*yint[  7]*zint[  4];
    eh[  59] += xint[ 20]*yint[  8]*zint[  5];
    eh[  60] += xint[ 21]*yint[  0]*zint[  6];
    eh[  60] += xint[ 22]*yint[  1]*zint[  7];
    eh[  60] += xint[ 23]*yint[  2]*zint[  8];
    eh[  61] += xint[ 18]*yint[  3]*zint[  6];
    eh[  61] += xint[ 19]*yint[  4]*zint[  7];
    eh[  61] += xint[ 20]*yint[  5]*zint[  8];
    eh[  62] += xint[ 18]*yint[  0]*zint[  9];
    eh[  62] += xint[ 19]*yint[  1]*zint[ 10];
    eh[  62] += xint[ 20]*yint[  2]*zint[ 11];
    eh[  63] += xint[ 15]*yint[ 12]*zint[  0];
    eh[  63] += xint[ 16]*yint[ 13]*zint[  1];
    eh[  63] += xint[ 17]*yint[ 14]*zint[  2];
    eh[  64] += xint[ 12]*yint[ 15]*zint[  0];
    eh[  64] += xint[ 13]*yint[ 16]*zint[  1];
    eh[  64] += xint[ 14]*yint[ 17]*zint[  2];
    eh[  65] += xint[ 12]*yint[ 12]*zint[  3];
    eh[  65] += xint[ 13]*yint[ 13]*zint[  4];
    eh[  65] += xint[ 14]*yint[ 14]*zint[  5];
    eh[  66] += xint[ 15]*yint[  0]*zint[ 12];
    eh[  66] += xint[ 16]*yint[  1]*zint[ 13];
    eh[  66] += xint[ 17]*yint[  2]*zint[ 14];
    eh[  67] += xint[ 12]*yint[  3]*zint[ 12];
    eh[  67] += xint[ 13]*yint[  4]*zint[ 13];
    eh[  67] += xint[ 14]*yint[  5]*zint[ 14];
    eh[  68] += xint[ 12]*yint[  0]*zint[ 15];
    eh[  68] += xint[ 13]*yint[  1]*zint[ 16];
    eh[  68] += xint[ 14]*yint[  2]*zint[ 17];
    eh[  69] += xint[ 15]*yint[  6]*zint[  6];
    eh[  69] += xint[ 16]*yint[  7]*zint[  7];
    eh[  69] += xint[ 17]*yint[  8]*zint[  8];
    eh[  70] += xint[ 12]*yint[  9]*zint[  6];
    eh[  70] += xint[ 13]*yint[ 10]*zint[  7];
    eh[  70] += xint[ 14]*yint[ 11]*zint[  8];
    eh[  71] += xint[ 12]*yint[  6]*zint[  9];
    eh[  71] += xint[ 13]*yint[  7]*zint[ 10];
    eh[  71] += xint[ 14]*yint[  8]*zint[ 11];
    eh[  72] += xint[  9]*yint[ 18]*zint[  0];
    eh[  72] += xint[ 10]*yint[ 19]*zint[  1];
    eh[  72] += xint[ 11]*yint[ 20]*zint[  2];
    eh[  73] += xint[  6]*yint[ 21]*zint[  0];
    eh[  73] += xint[  7]*yint[ 22]*zint[  1];
    eh[  73] += xint[  8]*yint[ 23]*zint[  2];
    eh[  74] += xint[  6]*yint[ 18]*zint[  3];
    eh[  74] += xint[  7]*yint[ 19]*zint[  4];
    eh[  74] += xint[  8]*yint[ 20]*zint[  5];
    eh[  75] += xint[  9]*yint[  0]*zint[ 18];
    eh[  75] += xint[ 10]*yint[  1]*zint[ 19];
    eh[  75] += xint[ 11]*yint[  2]*zint[ 20];
    eh[  76] += xint[  6]*yint[  3]*zint[ 18];
    eh[  76] += xint[  7]*yint[  4]*zint[ 19];
    eh[  76] += xint[  8]*yint[  5]*zint[ 20];
    eh[  77] += xint[  6]*yint[  0]*zint[ 21];
    eh[  77] += xint[  7]*yint[  1]*zint[ 22];
    eh[  77] += xint[  8]*yint[  2]*zint[ 23];
    eh[  78] += xint[  9]*yint[ 12]*zint[  6];
    eh[  78] += xint[ 10]*yint[ 13]*zint[  7];
    eh[  78] += xint[ 11]*yint[ 14]*zint[  8];
    eh[  79] += xint[  6]*yint[ 15]*zint[  6];
    eh[  79] += xint[  7]*yint[ 16]*zint[  7];
    eh[  79] += xint[  8]*yint[ 17]*zint[  8];
    eh[  80] += xint[  6]*yint[ 12]*zint[  9];
    eh[  80] += xint[  7]*yint[ 13]*zint[ 10];
    eh[  80] += xint[  8]*yint[ 14]*zint[ 11];
    eh[  81] += xint[  9]*yint[  6]*zint[ 12];
    eh[  81] += xint[ 10]*yint[  7]*zint[ 13];
    eh[  81] += xint[ 11]*yint[  8]*zint[ 14];
    eh[  82] += xint[  6]*yint[  9]*zint[ 12];
    eh[  82] += xint[  7]*yint[ 10]*zint[ 13];
    eh[  82] += xint[  8]*yint[ 11]*zint[ 14];
    eh[  83] += xint[  6]*yint[  6]*zint[ 15];
    eh[  83] += xint[  7]*yint[  7]*zint[ 16];
    eh[  83] += xint[  8]*yint[  8]*zint[ 17];
    eh[  84] += xint[  3]*yint[ 18]*zint[  6];
    eh[  84] += xint[  4]*yint[ 19]*zint[  7];
    eh[  84] += xint[  5]*yint[ 20]*zint[  8];
    eh[  85] += xint[  0]*yint[ 21]*zint[  6];
    eh[  85] += xint[  1]*yint[ 22]*zint[  7];
    eh[  85] += xint[  2]*yint[ 23]*zint[  8];
    eh[  86] += xint[  0]*yint[ 18]*zint[  9];
    eh[  86] += xint[  1]*yint[ 19]*zint[ 10];
    eh[  86] += xint[  2]*yint[ 20]*zint[ 11];
    eh[  87] += xint[  3]*yint[  6]*zint[ 18];
    eh[  87] += xint[  4]*yint[  7]*zint[ 19];
    eh[  87] += xint[  5]*yint[  8]*zint[ 20];
    eh[  88] += xint[  0]*yint[  9]*zint[ 18];
    eh[  88] += xint[  1]*yint[ 10]*zint[ 19];
    eh[  88] += xint[  2]*yint[ 11]*zint[ 20];
    eh[  89] += xint[  0]*yint[  6]*zint[ 21];
    eh[  89] += xint[  1]*yint[  7]*zint[ 22];
    eh[  89] += xint[  2]*yint[  8]*zint[ 23];
    eh[  90] += xint[  3]*yint[ 12]*zint[ 12];
    eh[  90] += xint[  4]*yint[ 13]*zint[ 13];
    eh[  90] += xint[  5]*yint[ 14]*zint[ 14];
    eh[  91] += xint[  0]*yint[ 15]*zint[ 12];
    eh[  91] += xint[  1]*yint[ 16]*zint[ 13];
    eh[  91] += xint[  2]*yint[ 17]*zint[ 14];
    eh[  92] += xint[  0]*yint[ 12]*zint[ 15];
    eh[  92] += xint[  1]*yint[ 13]*zint[ 16];
    eh[  92] += xint[  2]*yint[ 14]*zint[ 17];
}

void ofmo_twoint_core_rys_ddps( const int mythread,
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
    ofmo_hrr_clear_ddps( eh );
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
            ofmo_xyzint_ddps( F00, B00, B10, B01, C00, CP00,
                     xint, yint, zint );
            ofmo_form_ddps( xint, yint, zint, eh );
        }
    }
    ofmo_hrr_calc_ddps( eh, BA, DC );
    ofmo_hrr_coef_ddps( eh, DINT );
}

int ofmo_twoint_rys_ddps(
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
    max_nzeri = ebuf_max_nzeri - 6*6*3*1;
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
            ofmo_twoint_core_rys_ddps( mythread,
                    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
                    &psp_xiza[ijps0], BA,
                    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
                    &psp_xiza[klps0], DC,   AC,      DINTEG );
            ipat=((Lab != Lcd) || (ics==kcs && jcs>lcs) ? true : false);
            for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
                I2 = (iao*iao+iao)>>1;
                for ( j=0, jao=jao0; j<6; j++, jao++ ) {
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

int ofmo_twoint_direct_rys_ddps(
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
    max_nzeri -= 6*6*3*1;
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
            ofmo_twoint_core_rys_ddps( mythread,
                    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
                    &psp_xiza[ijps0], BA,
                    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
                    &psp_xiza[klps0], DC,   AC,      DINTEG );
            ipat=((Lab != Lcd) || (ics==kcs && jcs>lcs) ? true : false);
            for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
                I2 = (iao*iao+iao)>>1;
                for ( j=0, jao=jao0; j<6; j++, jao++ ) {
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
