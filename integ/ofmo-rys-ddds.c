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

static void ofmo_hrr_clear_ddds( double *eh ) {
    int i;
    // (DS|DS)
    for ( i=0; i<(0+36); i++ ) eh[i] = 0.e0;
    // (FS|DS)
    for ( i=36; i<(36+60); i++ ) eh[i] = 0.e0;
    // (GS|DS)
    for ( i=96; i<(96+90); i++ ) eh[i] = 0.e0;
}

static void ofmo_hrr_coef_ddds(
        double *eh, double *DINT ) {
    int i, j, k, l, iao, jao, kao, lao, ix;
    double coef_a, coef_ab, coef_abc;
    double *th;
    th = &eh[474];
    ix = 0;
    for ( i=0, iao=4; i<6; i++, iao++ ) {
        coef_a = DFACT[iao];
        for ( j=0, jao=4; j<6; j++, jao++ ) {
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

static void ofmo_hrr_calc_ddds( double *eh,
        const double BA[3], const double DC[3] ) {
    // HRR for (XX|XS)-type integral (center AB)
    // (DP,DS)
    eh[ 186] = eh[  36] - BA[0]*eh[   0];
    eh[ 187] = eh[  37] - BA[0]*eh[   1];
    eh[ 188] = eh[  38] - BA[0]*eh[   2];
    eh[ 189] = eh[  39] - BA[0]*eh[   3];
    eh[ 190] = eh[  40] - BA[0]*eh[   4];
    eh[ 191] = eh[  41] - BA[0]*eh[   5];
    eh[ 192] = eh[  54] - BA[1]*eh[   0];
    eh[ 193] = eh[  55] - BA[1]*eh[   1];
    eh[ 194] = eh[  56] - BA[1]*eh[   2];
    eh[ 195] = eh[  57] - BA[1]*eh[   3];
    eh[ 196] = eh[  58] - BA[1]*eh[   4];
    eh[ 197] = eh[  59] - BA[1]*eh[   5];
    eh[ 198] = eh[  60] - BA[2]*eh[   0];
    eh[ 199] = eh[  61] - BA[2]*eh[   1];
    eh[ 200] = eh[  62] - BA[2]*eh[   2];
    eh[ 201] = eh[  63] - BA[2]*eh[   3];
    eh[ 202] = eh[  64] - BA[2]*eh[   4];
    eh[ 203] = eh[  65] - BA[2]*eh[   5];
    eh[ 204] = eh[  66] - BA[0]*eh[   6];
    eh[ 205] = eh[  67] - BA[0]*eh[   7];
    eh[ 206] = eh[  68] - BA[0]*eh[   8];
    eh[ 207] = eh[  69] - BA[0]*eh[   9];
    eh[ 208] = eh[  70] - BA[0]*eh[  10];
    eh[ 209] = eh[  71] - BA[0]*eh[  11];
    eh[ 210] = eh[  42] - BA[1]*eh[   6];
    eh[ 211] = eh[  43] - BA[1]*eh[   7];
    eh[ 212] = eh[  44] - BA[1]*eh[   8];
    eh[ 213] = eh[  45] - BA[1]*eh[   9];
    eh[ 214] = eh[  46] - BA[1]*eh[  10];
    eh[ 215] = eh[  47] - BA[1]*eh[  11];
    eh[ 216] = eh[  84] - BA[2]*eh[   6];
    eh[ 217] = eh[  85] - BA[2]*eh[   7];
    eh[ 218] = eh[  86] - BA[2]*eh[   8];
    eh[ 219] = eh[  87] - BA[2]*eh[   9];
    eh[ 220] = eh[  88] - BA[2]*eh[  10];
    eh[ 221] = eh[  89] - BA[2]*eh[  11];
    eh[ 222] = eh[  72] - BA[0]*eh[  12];
    eh[ 223] = eh[  73] - BA[0]*eh[  13];
    eh[ 224] = eh[  74] - BA[0]*eh[  14];
    eh[ 225] = eh[  75] - BA[0]*eh[  15];
    eh[ 226] = eh[  76] - BA[0]*eh[  16];
    eh[ 227] = eh[  77] - BA[0]*eh[  17];
    eh[ 228] = eh[  90] - BA[1]*eh[  12];
    eh[ 229] = eh[  91] - BA[1]*eh[  13];
    eh[ 230] = eh[  92] - BA[1]*eh[  14];
    eh[ 231] = eh[  93] - BA[1]*eh[  15];
    eh[ 232] = eh[  94] - BA[1]*eh[  16];
    eh[ 233] = eh[  95] - BA[1]*eh[  17];
    eh[ 234] = eh[  48] - BA[2]*eh[  12];
    eh[ 235] = eh[  49] - BA[2]*eh[  13];
    eh[ 236] = eh[  50] - BA[2]*eh[  14];
    eh[ 237] = eh[  51] - BA[2]*eh[  15];
    eh[ 238] = eh[  52] - BA[2]*eh[  16];
    eh[ 239] = eh[  53] - BA[2]*eh[  17];
    eh[ 240] = eh[  54] - BA[0]*eh[  18];
    eh[ 241] = eh[  55] - BA[0]*eh[  19];
    eh[ 242] = eh[  56] - BA[0]*eh[  20];
    eh[ 243] = eh[  57] - BA[0]*eh[  21];
    eh[ 244] = eh[  58] - BA[0]*eh[  22];
    eh[ 245] = eh[  59] - BA[0]*eh[  23];
    eh[ 246] = eh[  66] - BA[1]*eh[  18];
    eh[ 247] = eh[  67] - BA[1]*eh[  19];
    eh[ 248] = eh[  68] - BA[1]*eh[  20];
    eh[ 249] = eh[  69] - BA[1]*eh[  21];
    eh[ 250] = eh[  70] - BA[1]*eh[  22];
    eh[ 251] = eh[  71] - BA[1]*eh[  23];
    eh[ 252] = eh[  78] - BA[2]*eh[  18];
    eh[ 253] = eh[  79] - BA[2]*eh[  19];
    eh[ 254] = eh[  80] - BA[2]*eh[  20];
    eh[ 255] = eh[  81] - BA[2]*eh[  21];
    eh[ 256] = eh[  82] - BA[2]*eh[  22];
    eh[ 257] = eh[  83] - BA[2]*eh[  23];
    eh[ 258] = eh[  60] - BA[0]*eh[  24];
    eh[ 259] = eh[  61] - BA[0]*eh[  25];
    eh[ 260] = eh[  62] - BA[0]*eh[  26];
    eh[ 261] = eh[  63] - BA[0]*eh[  27];
    eh[ 262] = eh[  64] - BA[0]*eh[  28];
    eh[ 263] = eh[  65] - BA[0]*eh[  29];
    eh[ 264] = eh[  78] - BA[1]*eh[  24];
    eh[ 265] = eh[  79] - BA[1]*eh[  25];
    eh[ 266] = eh[  80] - BA[1]*eh[  26];
    eh[ 267] = eh[  81] - BA[1]*eh[  27];
    eh[ 268] = eh[  82] - BA[1]*eh[  28];
    eh[ 269] = eh[  83] - BA[1]*eh[  29];
    eh[ 270] = eh[  72] - BA[2]*eh[  24];
    eh[ 271] = eh[  73] - BA[2]*eh[  25];
    eh[ 272] = eh[  74] - BA[2]*eh[  26];
    eh[ 273] = eh[  75] - BA[2]*eh[  27];
    eh[ 274] = eh[  76] - BA[2]*eh[  28];
    eh[ 275] = eh[  77] - BA[2]*eh[  29];
    eh[ 276] = eh[  78] - BA[0]*eh[  30];
    eh[ 277] = eh[  79] - BA[0]*eh[  31];
    eh[ 278] = eh[  80] - BA[0]*eh[  32];
    eh[ 279] = eh[  81] - BA[0]*eh[  33];
    eh[ 280] = eh[  82] - BA[0]*eh[  34];
    eh[ 281] = eh[  83] - BA[0]*eh[  35];
    eh[ 282] = eh[  84] - BA[1]*eh[  30];
    eh[ 283] = eh[  85] - BA[1]*eh[  31];
    eh[ 284] = eh[  86] - BA[1]*eh[  32];
    eh[ 285] = eh[  87] - BA[1]*eh[  33];
    eh[ 286] = eh[  88] - BA[1]*eh[  34];
    eh[ 287] = eh[  89] - BA[1]*eh[  35];
    eh[ 288] = eh[  90] - BA[2]*eh[  30];
    eh[ 289] = eh[  91] - BA[2]*eh[  31];
    eh[ 290] = eh[  92] - BA[2]*eh[  32];
    eh[ 291] = eh[  93] - BA[2]*eh[  33];
    eh[ 292] = eh[  94] - BA[2]*eh[  34];
    eh[ 293] = eh[  95] - BA[2]*eh[  35];
    // (FP,DS)
    eh[ 294] = eh[  96] - BA[0]*eh[  36];
    eh[ 295] = eh[  97] - BA[0]*eh[  37];
    eh[ 296] = eh[  98] - BA[0]*eh[  38];
    eh[ 297] = eh[  99] - BA[0]*eh[  39];
    eh[ 298] = eh[ 100] - BA[0]*eh[  40];
    eh[ 299] = eh[ 101] - BA[0]*eh[  41];
    eh[ 300] = eh[ 114] - BA[1]*eh[  36];
    eh[ 301] = eh[ 115] - BA[1]*eh[  37];
    eh[ 302] = eh[ 116] - BA[1]*eh[  38];
    eh[ 303] = eh[ 117] - BA[1]*eh[  39];
    eh[ 304] = eh[ 118] - BA[1]*eh[  40];
    eh[ 305] = eh[ 119] - BA[1]*eh[  41];
    eh[ 306] = eh[ 120] - BA[2]*eh[  36];
    eh[ 307] = eh[ 121] - BA[2]*eh[  37];
    eh[ 308] = eh[ 122] - BA[2]*eh[  38];
    eh[ 309] = eh[ 123] - BA[2]*eh[  39];
    eh[ 310] = eh[ 124] - BA[2]*eh[  40];
    eh[ 311] = eh[ 125] - BA[2]*eh[  41];
    eh[ 312] = eh[ 144] - BA[0]*eh[  42];
    eh[ 313] = eh[ 145] - BA[0]*eh[  43];
    eh[ 314] = eh[ 146] - BA[0]*eh[  44];
    eh[ 315] = eh[ 147] - BA[0]*eh[  45];
    eh[ 316] = eh[ 148] - BA[0]*eh[  46];
    eh[ 317] = eh[ 149] - BA[0]*eh[  47];
    eh[ 318] = eh[ 102] - BA[1]*eh[  42];
    eh[ 319] = eh[ 103] - BA[1]*eh[  43];
    eh[ 320] = eh[ 104] - BA[1]*eh[  44];
    eh[ 321] = eh[ 105] - BA[1]*eh[  45];
    eh[ 322] = eh[ 106] - BA[1]*eh[  46];
    eh[ 323] = eh[ 107] - BA[1]*eh[  47];
    eh[ 324] = eh[ 168] - BA[2]*eh[  42];
    eh[ 325] = eh[ 169] - BA[2]*eh[  43];
    eh[ 326] = eh[ 170] - BA[2]*eh[  44];
    eh[ 327] = eh[ 171] - BA[2]*eh[  45];
    eh[ 328] = eh[ 172] - BA[2]*eh[  46];
    eh[ 329] = eh[ 173] - BA[2]*eh[  47];
    eh[ 330] = eh[ 150] - BA[0]*eh[  48];
    eh[ 331] = eh[ 151] - BA[0]*eh[  49];
    eh[ 332] = eh[ 152] - BA[0]*eh[  50];
    eh[ 333] = eh[ 153] - BA[0]*eh[  51];
    eh[ 334] = eh[ 154] - BA[0]*eh[  52];
    eh[ 335] = eh[ 155] - BA[0]*eh[  53];
    eh[ 336] = eh[ 174] - BA[1]*eh[  48];
    eh[ 337] = eh[ 175] - BA[1]*eh[  49];
    eh[ 338] = eh[ 176] - BA[1]*eh[  50];
    eh[ 339] = eh[ 177] - BA[1]*eh[  51];
    eh[ 340] = eh[ 178] - BA[1]*eh[  52];
    eh[ 341] = eh[ 179] - BA[1]*eh[  53];
    eh[ 342] = eh[ 108] - BA[2]*eh[  48];
    eh[ 343] = eh[ 109] - BA[2]*eh[  49];
    eh[ 344] = eh[ 110] - BA[2]*eh[  50];
    eh[ 345] = eh[ 111] - BA[2]*eh[  51];
    eh[ 346] = eh[ 112] - BA[2]*eh[  52];
    eh[ 347] = eh[ 113] - BA[2]*eh[  53];
    eh[ 348] = eh[ 114] - BA[0]*eh[  54];
    eh[ 349] = eh[ 115] - BA[0]*eh[  55];
    eh[ 350] = eh[ 116] - BA[0]*eh[  56];
    eh[ 351] = eh[ 117] - BA[0]*eh[  57];
    eh[ 352] = eh[ 118] - BA[0]*eh[  58];
    eh[ 353] = eh[ 119] - BA[0]*eh[  59];
    eh[ 354] = eh[ 126] - BA[1]*eh[  54];
    eh[ 355] = eh[ 127] - BA[1]*eh[  55];
    eh[ 356] = eh[ 128] - BA[1]*eh[  56];
    eh[ 357] = eh[ 129] - BA[1]*eh[  57];
    eh[ 358] = eh[ 130] - BA[1]*eh[  58];
    eh[ 359] = eh[ 131] - BA[1]*eh[  59];
    eh[ 360] = eh[ 138] - BA[2]*eh[  54];
    eh[ 361] = eh[ 139] - BA[2]*eh[  55];
    eh[ 362] = eh[ 140] - BA[2]*eh[  56];
    eh[ 363] = eh[ 141] - BA[2]*eh[  57];
    eh[ 364] = eh[ 142] - BA[2]*eh[  58];
    eh[ 365] = eh[ 143] - BA[2]*eh[  59];
    eh[ 366] = eh[ 120] - BA[0]*eh[  60];
    eh[ 367] = eh[ 121] - BA[0]*eh[  61];
    eh[ 368] = eh[ 122] - BA[0]*eh[  62];
    eh[ 369] = eh[ 123] - BA[0]*eh[  63];
    eh[ 370] = eh[ 124] - BA[0]*eh[  64];
    eh[ 371] = eh[ 125] - BA[0]*eh[  65];
    eh[ 372] = eh[ 138] - BA[1]*eh[  60];
    eh[ 373] = eh[ 139] - BA[1]*eh[  61];
    eh[ 374] = eh[ 140] - BA[1]*eh[  62];
    eh[ 375] = eh[ 141] - BA[1]*eh[  63];
    eh[ 376] = eh[ 142] - BA[1]*eh[  64];
    eh[ 377] = eh[ 143] - BA[1]*eh[  65];
    eh[ 378] = eh[ 132] - BA[2]*eh[  60];
    eh[ 379] = eh[ 133] - BA[2]*eh[  61];
    eh[ 380] = eh[ 134] - BA[2]*eh[  62];
    eh[ 381] = eh[ 135] - BA[2]*eh[  63];
    eh[ 382] = eh[ 136] - BA[2]*eh[  64];
    eh[ 383] = eh[ 137] - BA[2]*eh[  65];
    eh[ 384] = eh[ 126] - BA[0]*eh[  66];
    eh[ 385] = eh[ 127] - BA[0]*eh[  67];
    eh[ 386] = eh[ 128] - BA[0]*eh[  68];
    eh[ 387] = eh[ 129] - BA[0]*eh[  69];
    eh[ 388] = eh[ 130] - BA[0]*eh[  70];
    eh[ 389] = eh[ 131] - BA[0]*eh[  71];
    eh[ 390] = eh[ 144] - BA[1]*eh[  66];
    eh[ 391] = eh[ 145] - BA[1]*eh[  67];
    eh[ 392] = eh[ 146] - BA[1]*eh[  68];
    eh[ 393] = eh[ 147] - BA[1]*eh[  69];
    eh[ 394] = eh[ 148] - BA[1]*eh[  70];
    eh[ 395] = eh[ 149] - BA[1]*eh[  71];
    eh[ 396] = eh[ 156] - BA[2]*eh[  66];
    eh[ 397] = eh[ 157] - BA[2]*eh[  67];
    eh[ 398] = eh[ 158] - BA[2]*eh[  68];
    eh[ 399] = eh[ 159] - BA[2]*eh[  69];
    eh[ 400] = eh[ 160] - BA[2]*eh[  70];
    eh[ 401] = eh[ 161] - BA[2]*eh[  71];
    eh[ 402] = eh[ 132] - BA[0]*eh[  72];
    eh[ 403] = eh[ 133] - BA[0]*eh[  73];
    eh[ 404] = eh[ 134] - BA[0]*eh[  74];
    eh[ 405] = eh[ 135] - BA[0]*eh[  75];
    eh[ 406] = eh[ 136] - BA[0]*eh[  76];
    eh[ 407] = eh[ 137] - BA[0]*eh[  77];
    eh[ 408] = eh[ 162] - BA[1]*eh[  72];
    eh[ 409] = eh[ 163] - BA[1]*eh[  73];
    eh[ 410] = eh[ 164] - BA[1]*eh[  74];
    eh[ 411] = eh[ 165] - BA[1]*eh[  75];
    eh[ 412] = eh[ 166] - BA[1]*eh[  76];
    eh[ 413] = eh[ 167] - BA[1]*eh[  77];
    eh[ 414] = eh[ 150] - BA[2]*eh[  72];
    eh[ 415] = eh[ 151] - BA[2]*eh[  73];
    eh[ 416] = eh[ 152] - BA[2]*eh[  74];
    eh[ 417] = eh[ 153] - BA[2]*eh[  75];
    eh[ 418] = eh[ 154] - BA[2]*eh[  76];
    eh[ 419] = eh[ 155] - BA[2]*eh[  77];
    eh[ 420] = eh[ 138] - BA[0]*eh[  78];
    eh[ 421] = eh[ 139] - BA[0]*eh[  79];
    eh[ 422] = eh[ 140] - BA[0]*eh[  80];
    eh[ 423] = eh[ 141] - BA[0]*eh[  81];
    eh[ 424] = eh[ 142] - BA[0]*eh[  82];
    eh[ 425] = eh[ 143] - BA[0]*eh[  83];
    eh[ 426] = eh[ 156] - BA[1]*eh[  78];
    eh[ 427] = eh[ 157] - BA[1]*eh[  79];
    eh[ 428] = eh[ 158] - BA[1]*eh[  80];
    eh[ 429] = eh[ 159] - BA[1]*eh[  81];
    eh[ 430] = eh[ 160] - BA[1]*eh[  82];
    eh[ 431] = eh[ 161] - BA[1]*eh[  83];
    eh[ 432] = eh[ 162] - BA[2]*eh[  78];
    eh[ 433] = eh[ 163] - BA[2]*eh[  79];
    eh[ 434] = eh[ 164] - BA[2]*eh[  80];
    eh[ 435] = eh[ 165] - BA[2]*eh[  81];
    eh[ 436] = eh[ 166] - BA[2]*eh[  82];
    eh[ 437] = eh[ 167] - BA[2]*eh[  83];
    eh[ 438] = eh[ 156] - BA[0]*eh[  84];
    eh[ 439] = eh[ 157] - BA[0]*eh[  85];
    eh[ 440] = eh[ 158] - BA[0]*eh[  86];
    eh[ 441] = eh[ 159] - BA[0]*eh[  87];
    eh[ 442] = eh[ 160] - BA[0]*eh[  88];
    eh[ 443] = eh[ 161] - BA[0]*eh[  89];
    eh[ 444] = eh[ 168] - BA[1]*eh[  84];
    eh[ 445] = eh[ 169] - BA[1]*eh[  85];
    eh[ 446] = eh[ 170] - BA[1]*eh[  86];
    eh[ 447] = eh[ 171] - BA[1]*eh[  87];
    eh[ 448] = eh[ 172] - BA[1]*eh[  88];
    eh[ 449] = eh[ 173] - BA[1]*eh[  89];
    eh[ 450] = eh[ 180] - BA[2]*eh[  84];
    eh[ 451] = eh[ 181] - BA[2]*eh[  85];
    eh[ 452] = eh[ 182] - BA[2]*eh[  86];
    eh[ 453] = eh[ 183] - BA[2]*eh[  87];
    eh[ 454] = eh[ 184] - BA[2]*eh[  88];
    eh[ 455] = eh[ 185] - BA[2]*eh[  89];
    eh[ 456] = eh[ 162] - BA[0]*eh[  90];
    eh[ 457] = eh[ 163] - BA[0]*eh[  91];
    eh[ 458] = eh[ 164] - BA[0]*eh[  92];
    eh[ 459] = eh[ 165] - BA[0]*eh[  93];
    eh[ 460] = eh[ 166] - BA[0]*eh[  94];
    eh[ 461] = eh[ 167] - BA[0]*eh[  95];
    eh[ 462] = eh[ 180] - BA[1]*eh[  90];
    eh[ 463] = eh[ 181] - BA[1]*eh[  91];
    eh[ 464] = eh[ 182] - BA[1]*eh[  92];
    eh[ 465] = eh[ 183] - BA[1]*eh[  93];
    eh[ 466] = eh[ 184] - BA[1]*eh[  94];
    eh[ 467] = eh[ 185] - BA[1]*eh[  95];
    eh[ 468] = eh[ 174] - BA[2]*eh[  90];
    eh[ 469] = eh[ 175] - BA[2]*eh[  91];
    eh[ 470] = eh[ 176] - BA[2]*eh[  92];
    eh[ 471] = eh[ 177] - BA[2]*eh[  93];
    eh[ 472] = eh[ 178] - BA[2]*eh[  94];
    eh[ 473] = eh[ 179] - BA[2]*eh[  95];
    // (DD,DS)
    eh[ 474] = eh[ 294] - BA[0]*eh[ 186];
    eh[ 475] = eh[ 295] - BA[0]*eh[ 187];
    eh[ 476] = eh[ 296] - BA[0]*eh[ 188];
    eh[ 477] = eh[ 297] - BA[0]*eh[ 189];
    eh[ 478] = eh[ 298] - BA[0]*eh[ 190];
    eh[ 479] = eh[ 299] - BA[0]*eh[ 191];
    eh[ 480] = eh[ 354] - BA[1]*eh[ 192];
    eh[ 481] = eh[ 355] - BA[1]*eh[ 193];
    eh[ 482] = eh[ 356] - BA[1]*eh[ 194];
    eh[ 483] = eh[ 357] - BA[1]*eh[ 195];
    eh[ 484] = eh[ 358] - BA[1]*eh[ 196];
    eh[ 485] = eh[ 359] - BA[1]*eh[ 197];
    eh[ 486] = eh[ 378] - BA[2]*eh[ 198];
    eh[ 487] = eh[ 379] - BA[2]*eh[ 199];
    eh[ 488] = eh[ 380] - BA[2]*eh[ 200];
    eh[ 489] = eh[ 381] - BA[2]*eh[ 201];
    eh[ 490] = eh[ 382] - BA[2]*eh[ 202];
    eh[ 491] = eh[ 383] - BA[2]*eh[ 203];
    eh[ 492] = eh[ 300] - BA[0]*eh[ 192];
    eh[ 493] = eh[ 301] - BA[0]*eh[ 193];
    eh[ 494] = eh[ 302] - BA[0]*eh[ 194];
    eh[ 495] = eh[ 303] - BA[0]*eh[ 195];
    eh[ 496] = eh[ 304] - BA[0]*eh[ 196];
    eh[ 497] = eh[ 305] - BA[0]*eh[ 197];
    eh[ 498] = eh[ 306] - BA[0]*eh[ 198];
    eh[ 499] = eh[ 307] - BA[0]*eh[ 199];
    eh[ 500] = eh[ 308] - BA[0]*eh[ 200];
    eh[ 501] = eh[ 309] - BA[0]*eh[ 201];
    eh[ 502] = eh[ 310] - BA[0]*eh[ 202];
    eh[ 503] = eh[ 311] - BA[0]*eh[ 203];
    eh[ 504] = eh[ 360] - BA[1]*eh[ 198];
    eh[ 505] = eh[ 361] - BA[1]*eh[ 199];
    eh[ 506] = eh[ 362] - BA[1]*eh[ 200];
    eh[ 507] = eh[ 363] - BA[1]*eh[ 201];
    eh[ 508] = eh[ 364] - BA[1]*eh[ 202];
    eh[ 509] = eh[ 365] - BA[1]*eh[ 203];
    eh[ 510] = eh[ 384] - BA[0]*eh[ 204];
    eh[ 511] = eh[ 385] - BA[0]*eh[ 205];
    eh[ 512] = eh[ 386] - BA[0]*eh[ 206];
    eh[ 513] = eh[ 387] - BA[0]*eh[ 207];
    eh[ 514] = eh[ 388] - BA[0]*eh[ 208];
    eh[ 515] = eh[ 389] - BA[0]*eh[ 209];
    eh[ 516] = eh[ 318] - BA[1]*eh[ 210];
    eh[ 517] = eh[ 319] - BA[1]*eh[ 211];
    eh[ 518] = eh[ 320] - BA[1]*eh[ 212];
    eh[ 519] = eh[ 321] - BA[1]*eh[ 213];
    eh[ 520] = eh[ 322] - BA[1]*eh[ 214];
    eh[ 521] = eh[ 323] - BA[1]*eh[ 215];
    eh[ 522] = eh[ 450] - BA[2]*eh[ 216];
    eh[ 523] = eh[ 451] - BA[2]*eh[ 217];
    eh[ 524] = eh[ 452] - BA[2]*eh[ 218];
    eh[ 525] = eh[ 453] - BA[2]*eh[ 219];
    eh[ 526] = eh[ 454] - BA[2]*eh[ 220];
    eh[ 527] = eh[ 455] - BA[2]*eh[ 221];
    eh[ 528] = eh[ 390] - BA[0]*eh[ 210];
    eh[ 529] = eh[ 391] - BA[0]*eh[ 211];
    eh[ 530] = eh[ 392] - BA[0]*eh[ 212];
    eh[ 531] = eh[ 393] - BA[0]*eh[ 213];
    eh[ 532] = eh[ 394] - BA[0]*eh[ 214];
    eh[ 533] = eh[ 395] - BA[0]*eh[ 215];
    eh[ 534] = eh[ 396] - BA[0]*eh[ 216];
    eh[ 535] = eh[ 397] - BA[0]*eh[ 217];
    eh[ 536] = eh[ 398] - BA[0]*eh[ 218];
    eh[ 537] = eh[ 399] - BA[0]*eh[ 219];
    eh[ 538] = eh[ 400] - BA[0]*eh[ 220];
    eh[ 539] = eh[ 401] - BA[0]*eh[ 221];
    eh[ 540] = eh[ 324] - BA[1]*eh[ 216];
    eh[ 541] = eh[ 325] - BA[1]*eh[ 217];
    eh[ 542] = eh[ 326] - BA[1]*eh[ 218];
    eh[ 543] = eh[ 327] - BA[1]*eh[ 219];
    eh[ 544] = eh[ 328] - BA[1]*eh[ 220];
    eh[ 545] = eh[ 329] - BA[1]*eh[ 221];
    eh[ 546] = eh[ 402] - BA[0]*eh[ 222];
    eh[ 547] = eh[ 403] - BA[0]*eh[ 223];
    eh[ 548] = eh[ 404] - BA[0]*eh[ 224];
    eh[ 549] = eh[ 405] - BA[0]*eh[ 225];
    eh[ 550] = eh[ 406] - BA[0]*eh[ 226];
    eh[ 551] = eh[ 407] - BA[0]*eh[ 227];
    eh[ 552] = eh[ 462] - BA[1]*eh[ 228];
    eh[ 553] = eh[ 463] - BA[1]*eh[ 229];
    eh[ 554] = eh[ 464] - BA[1]*eh[ 230];
    eh[ 555] = eh[ 465] - BA[1]*eh[ 231];
    eh[ 556] = eh[ 466] - BA[1]*eh[ 232];
    eh[ 557] = eh[ 467] - BA[1]*eh[ 233];
    eh[ 558] = eh[ 342] - BA[2]*eh[ 234];
    eh[ 559] = eh[ 343] - BA[2]*eh[ 235];
    eh[ 560] = eh[ 344] - BA[2]*eh[ 236];
    eh[ 561] = eh[ 345] - BA[2]*eh[ 237];
    eh[ 562] = eh[ 346] - BA[2]*eh[ 238];
    eh[ 563] = eh[ 347] - BA[2]*eh[ 239];
    eh[ 564] = eh[ 408] - BA[0]*eh[ 228];
    eh[ 565] = eh[ 409] - BA[0]*eh[ 229];
    eh[ 566] = eh[ 410] - BA[0]*eh[ 230];
    eh[ 567] = eh[ 411] - BA[0]*eh[ 231];
    eh[ 568] = eh[ 412] - BA[0]*eh[ 232];
    eh[ 569] = eh[ 413] - BA[0]*eh[ 233];
    eh[ 570] = eh[ 414] - BA[0]*eh[ 234];
    eh[ 571] = eh[ 415] - BA[0]*eh[ 235];
    eh[ 572] = eh[ 416] - BA[0]*eh[ 236];
    eh[ 573] = eh[ 417] - BA[0]*eh[ 237];
    eh[ 574] = eh[ 418] - BA[0]*eh[ 238];
    eh[ 575] = eh[ 419] - BA[0]*eh[ 239];
    eh[ 576] = eh[ 468] - BA[1]*eh[ 234];
    eh[ 577] = eh[ 469] - BA[1]*eh[ 235];
    eh[ 578] = eh[ 470] - BA[1]*eh[ 236];
    eh[ 579] = eh[ 471] - BA[1]*eh[ 237];
    eh[ 580] = eh[ 472] - BA[1]*eh[ 238];
    eh[ 581] = eh[ 473] - BA[1]*eh[ 239];
    eh[ 582] = eh[ 348] - BA[0]*eh[ 240];
    eh[ 583] = eh[ 349] - BA[0]*eh[ 241];
    eh[ 584] = eh[ 350] - BA[0]*eh[ 242];
    eh[ 585] = eh[ 351] - BA[0]*eh[ 243];
    eh[ 586] = eh[ 352] - BA[0]*eh[ 244];
    eh[ 587] = eh[ 353] - BA[0]*eh[ 245];
    eh[ 588] = eh[ 390] - BA[1]*eh[ 246];
    eh[ 589] = eh[ 391] - BA[1]*eh[ 247];
    eh[ 590] = eh[ 392] - BA[1]*eh[ 248];
    eh[ 591] = eh[ 393] - BA[1]*eh[ 249];
    eh[ 592] = eh[ 394] - BA[1]*eh[ 250];
    eh[ 593] = eh[ 395] - BA[1]*eh[ 251];
    eh[ 594] = eh[ 432] - BA[2]*eh[ 252];
    eh[ 595] = eh[ 433] - BA[2]*eh[ 253];
    eh[ 596] = eh[ 434] - BA[2]*eh[ 254];
    eh[ 597] = eh[ 435] - BA[2]*eh[ 255];
    eh[ 598] = eh[ 436] - BA[2]*eh[ 256];
    eh[ 599] = eh[ 437] - BA[2]*eh[ 257];
    eh[ 600] = eh[ 354] - BA[0]*eh[ 246];
    eh[ 601] = eh[ 355] - BA[0]*eh[ 247];
    eh[ 602] = eh[ 356] - BA[0]*eh[ 248];
    eh[ 603] = eh[ 357] - BA[0]*eh[ 249];
    eh[ 604] = eh[ 358] - BA[0]*eh[ 250];
    eh[ 605] = eh[ 359] - BA[0]*eh[ 251];
    eh[ 606] = eh[ 360] - BA[0]*eh[ 252];
    eh[ 607] = eh[ 361] - BA[0]*eh[ 253];
    eh[ 608] = eh[ 362] - BA[0]*eh[ 254];
    eh[ 609] = eh[ 363] - BA[0]*eh[ 255];
    eh[ 610] = eh[ 364] - BA[0]*eh[ 256];
    eh[ 611] = eh[ 365] - BA[0]*eh[ 257];
    eh[ 612] = eh[ 396] - BA[1]*eh[ 252];
    eh[ 613] = eh[ 397] - BA[1]*eh[ 253];
    eh[ 614] = eh[ 398] - BA[1]*eh[ 254];
    eh[ 615] = eh[ 399] - BA[1]*eh[ 255];
    eh[ 616] = eh[ 400] - BA[1]*eh[ 256];
    eh[ 617] = eh[ 401] - BA[1]*eh[ 257];
    eh[ 618] = eh[ 366] - BA[0]*eh[ 258];
    eh[ 619] = eh[ 367] - BA[0]*eh[ 259];
    eh[ 620] = eh[ 368] - BA[0]*eh[ 260];
    eh[ 621] = eh[ 369] - BA[0]*eh[ 261];
    eh[ 622] = eh[ 370] - BA[0]*eh[ 262];
    eh[ 623] = eh[ 371] - BA[0]*eh[ 263];
    eh[ 624] = eh[ 426] - BA[1]*eh[ 264];
    eh[ 625] = eh[ 427] - BA[1]*eh[ 265];
    eh[ 626] = eh[ 428] - BA[1]*eh[ 266];
    eh[ 627] = eh[ 429] - BA[1]*eh[ 267];
    eh[ 628] = eh[ 430] - BA[1]*eh[ 268];
    eh[ 629] = eh[ 431] - BA[1]*eh[ 269];
    eh[ 630] = eh[ 414] - BA[2]*eh[ 270];
    eh[ 631] = eh[ 415] - BA[2]*eh[ 271];
    eh[ 632] = eh[ 416] - BA[2]*eh[ 272];
    eh[ 633] = eh[ 417] - BA[2]*eh[ 273];
    eh[ 634] = eh[ 418] - BA[2]*eh[ 274];
    eh[ 635] = eh[ 419] - BA[2]*eh[ 275];
    eh[ 636] = eh[ 372] - BA[0]*eh[ 264];
    eh[ 637] = eh[ 373] - BA[0]*eh[ 265];
    eh[ 638] = eh[ 374] - BA[0]*eh[ 266];
    eh[ 639] = eh[ 375] - BA[0]*eh[ 267];
    eh[ 640] = eh[ 376] - BA[0]*eh[ 268];
    eh[ 641] = eh[ 377] - BA[0]*eh[ 269];
    eh[ 642] = eh[ 378] - BA[0]*eh[ 270];
    eh[ 643] = eh[ 379] - BA[0]*eh[ 271];
    eh[ 644] = eh[ 380] - BA[0]*eh[ 272];
    eh[ 645] = eh[ 381] - BA[0]*eh[ 273];
    eh[ 646] = eh[ 382] - BA[0]*eh[ 274];
    eh[ 647] = eh[ 383] - BA[0]*eh[ 275];
    eh[ 648] = eh[ 432] - BA[1]*eh[ 270];
    eh[ 649] = eh[ 433] - BA[1]*eh[ 271];
    eh[ 650] = eh[ 434] - BA[1]*eh[ 272];
    eh[ 651] = eh[ 435] - BA[1]*eh[ 273];
    eh[ 652] = eh[ 436] - BA[1]*eh[ 274];
    eh[ 653] = eh[ 437] - BA[1]*eh[ 275];
    eh[ 654] = eh[ 420] - BA[0]*eh[ 276];
    eh[ 655] = eh[ 421] - BA[0]*eh[ 277];
    eh[ 656] = eh[ 422] - BA[0]*eh[ 278];
    eh[ 657] = eh[ 423] - BA[0]*eh[ 279];
    eh[ 658] = eh[ 424] - BA[0]*eh[ 280];
    eh[ 659] = eh[ 425] - BA[0]*eh[ 281];
    eh[ 660] = eh[ 444] - BA[1]*eh[ 282];
    eh[ 661] = eh[ 445] - BA[1]*eh[ 283];
    eh[ 662] = eh[ 446] - BA[1]*eh[ 284];
    eh[ 663] = eh[ 447] - BA[1]*eh[ 285];
    eh[ 664] = eh[ 448] - BA[1]*eh[ 286];
    eh[ 665] = eh[ 449] - BA[1]*eh[ 287];
    eh[ 666] = eh[ 468] - BA[2]*eh[ 288];
    eh[ 667] = eh[ 469] - BA[2]*eh[ 289];
    eh[ 668] = eh[ 470] - BA[2]*eh[ 290];
    eh[ 669] = eh[ 471] - BA[2]*eh[ 291];
    eh[ 670] = eh[ 472] - BA[2]*eh[ 292];
    eh[ 671] = eh[ 473] - BA[2]*eh[ 293];
    eh[ 672] = eh[ 426] - BA[0]*eh[ 282];
    eh[ 673] = eh[ 427] - BA[0]*eh[ 283];
    eh[ 674] = eh[ 428] - BA[0]*eh[ 284];
    eh[ 675] = eh[ 429] - BA[0]*eh[ 285];
    eh[ 676] = eh[ 430] - BA[0]*eh[ 286];
    eh[ 677] = eh[ 431] - BA[0]*eh[ 287];
    eh[ 678] = eh[ 432] - BA[0]*eh[ 288];
    eh[ 679] = eh[ 433] - BA[0]*eh[ 289];
    eh[ 680] = eh[ 434] - BA[0]*eh[ 290];
    eh[ 681] = eh[ 435] - BA[0]*eh[ 291];
    eh[ 682] = eh[ 436] - BA[0]*eh[ 292];
    eh[ 683] = eh[ 437] - BA[0]*eh[ 293];
    eh[ 684] = eh[ 450] - BA[1]*eh[ 288];
    eh[ 685] = eh[ 451] - BA[1]*eh[ 289];
    eh[ 686] = eh[ 452] - BA[1]*eh[ 290];
    eh[ 687] = eh[ 453] - BA[1]*eh[ 291];
    eh[ 688] = eh[ 454] - BA[1]*eh[ 292];
    eh[ 689] = eh[ 455] - BA[1]*eh[ 293];
}

static void ofmo_xyzint_ddds(
        const double *F00, const double *B00, const double *B10,
        const double *B01, const double *C00, const double *CP00,
        double *xint, double *yint, double *zint ) {
    int Lab, Lcd;
    int m, m3, N, M, ix3, ix2, ix1, ix0, nroot;
    double C10[4], CP10[4], CP01[4], C01[4];
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
    xint[  3]=1.e0;
    yint[  3]=1.e0;
    zint[  3]=F00[3];
    // (1,0)
    xint[ 12]=C00[ 0];
    yint[ 12]=C00[ 1];
    zint[ 12]=C00[ 2]*F00[0];
    xint[ 13]=C00[ 3];
    yint[ 13]=C00[ 4];
    zint[ 13]=C00[ 5]*F00[1];
    xint[ 14]=C00[ 6];
    yint[ 14]=C00[ 7];
    zint[ 14]=C00[ 8]*F00[2];
    xint[ 15]=C00[ 9];
    yint[ 15]=C00[10];
    zint[ 15]=C00[11]*F00[3];
    // (0,1)
    xint[  4]=CP00[ 0];
    yint[  4]=CP00[ 1];
    zint[  4]=CP00[ 2]*F00[0];
    xint[  5]=CP00[ 3];
    yint[  5]=CP00[ 4];
    zint[  5]=CP00[ 5]*F00[1];
    xint[  6]=CP00[ 6];
    yint[  6]=CP00[ 7];
    zint[  6]=CP00[ 8]*F00[2];
    xint[  7]=CP00[ 9];
    yint[  7]=CP00[10];
    zint[  7]=CP00[11]*F00[3];
    // (1,1)
    xint[ 16]=CP00[ 0]*xint[ 12]+B00[0];
    yint[ 16]=CP00[ 1]*yint[ 12]+B00[0];
    zint[ 16]=CP00[ 2]*zint[ 12]+B00[0]*F00[0];
    xint[ 17]=CP00[ 3]*xint[ 13]+B00[1];
    yint[ 17]=CP00[ 4]*yint[ 13]+B00[1];
    zint[ 17]=CP00[ 5]*zint[ 13]+B00[1]*F00[1];
    xint[ 18]=CP00[ 6]*xint[ 14]+B00[2];
    yint[ 18]=CP00[ 7]*yint[ 14]+B00[2];
    zint[ 18]=CP00[ 8]*zint[ 14]+B00[2]*F00[2];
    xint[ 19]=CP00[ 9]*xint[ 15]+B00[3];
    yint[ 19]=CP00[10]*yint[ 15]+B00[3];
    zint[ 19]=CP00[11]*zint[ 15]+B00[3]*F00[3];
    // (N,0) and (N,1)
    for ( m=0; m<4; m++ ) {
        C10[m]  = 0.e0;
        CP10[m] = B00[m];
    }
    // (2,0)
    C10[0] += B10[0];
    xint[ 24]=C00[ 0]*xint[ 12]+C10[0]*xint[  0];
    yint[ 24]=C00[ 1]*yint[ 12]+C10[0]*yint[  0];
    zint[ 24]=C00[ 2]*zint[ 12]+C10[0]*zint[  0];
    C10[1] += B10[1];
    xint[ 25]=C00[ 3]*xint[ 13]+C10[1]*xint[  1];
    yint[ 25]=C00[ 4]*yint[ 13]+C10[1]*yint[  1];
    zint[ 25]=C00[ 5]*zint[ 13]+C10[1]*zint[  1];
    C10[2] += B10[2];
    xint[ 26]=C00[ 6]*xint[ 14]+C10[2]*xint[  2];
    yint[ 26]=C00[ 7]*yint[ 14]+C10[2]*yint[  2];
    zint[ 26]=C00[ 8]*zint[ 14]+C10[2]*zint[  2];
    C10[3] += B10[3];
    xint[ 27]=C00[ 9]*xint[ 15]+C10[3]*xint[  3];
    yint[ 27]=C00[10]*yint[ 15]+C10[3]*yint[  3];
    zint[ 27]=C00[11]*zint[ 15]+C10[3]*zint[  3];
    // (2,1)
    CP10[0] += B00[0];
    xint[ 28]=CP00[ 0]*xint[ 24]+CP10[0]*xint[ 12];
    yint[ 28]=CP00[ 1]*yint[ 24]+CP10[0]*yint[ 12];
    zint[ 28]=CP00[ 2]*zint[ 24]+CP10[0]*zint[ 12];
    CP10[1] += B00[1];
    xint[ 29]=CP00[ 3]*xint[ 25]+CP10[1]*xint[ 13];
    yint[ 29]=CP00[ 4]*yint[ 25]+CP10[1]*yint[ 13];
    zint[ 29]=CP00[ 5]*zint[ 25]+CP10[1]*zint[ 13];
    CP10[2] += B00[2];
    xint[ 30]=CP00[ 6]*xint[ 26]+CP10[2]*xint[ 14];
    yint[ 30]=CP00[ 7]*yint[ 26]+CP10[2]*yint[ 14];
    zint[ 30]=CP00[ 8]*zint[ 26]+CP10[2]*zint[ 14];
    CP10[3] += B00[3];
    xint[ 31]=CP00[ 9]*xint[ 27]+CP10[3]*xint[ 15];
    yint[ 31]=CP00[10]*yint[ 27]+CP10[3]*yint[ 15];
    zint[ 31]=CP00[11]*zint[ 27]+CP10[3]*zint[ 15];
    // (3,0)
    C10[0] += B10[0];
    xint[ 36]=C00[ 0]*xint[ 24]+C10[0]*xint[ 12];
    yint[ 36]=C00[ 1]*yint[ 24]+C10[0]*yint[ 12];
    zint[ 36]=C00[ 2]*zint[ 24]+C10[0]*zint[ 12];
    C10[1] += B10[1];
    xint[ 37]=C00[ 3]*xint[ 25]+C10[1]*xint[ 13];
    yint[ 37]=C00[ 4]*yint[ 25]+C10[1]*yint[ 13];
    zint[ 37]=C00[ 5]*zint[ 25]+C10[1]*zint[ 13];
    C10[2] += B10[2];
    xint[ 38]=C00[ 6]*xint[ 26]+C10[2]*xint[ 14];
    yint[ 38]=C00[ 7]*yint[ 26]+C10[2]*yint[ 14];
    zint[ 38]=C00[ 8]*zint[ 26]+C10[2]*zint[ 14];
    C10[3] += B10[3];
    xint[ 39]=C00[ 9]*xint[ 27]+C10[3]*xint[ 15];
    yint[ 39]=C00[10]*yint[ 27]+C10[3]*yint[ 15];
    zint[ 39]=C00[11]*zint[ 27]+C10[3]*zint[ 15];
    // (3,1)
    CP10[0] += B00[0];
    xint[ 40]=CP00[ 0]*xint[ 36]+CP10[0]*xint[ 24];
    yint[ 40]=CP00[ 1]*yint[ 36]+CP10[0]*yint[ 24];
    zint[ 40]=CP00[ 2]*zint[ 36]+CP10[0]*zint[ 24];
    CP10[1] += B00[1];
    xint[ 41]=CP00[ 3]*xint[ 37]+CP10[1]*xint[ 25];
    yint[ 41]=CP00[ 4]*yint[ 37]+CP10[1]*yint[ 25];
    zint[ 41]=CP00[ 5]*zint[ 37]+CP10[1]*zint[ 25];
    CP10[2] += B00[2];
    xint[ 42]=CP00[ 6]*xint[ 38]+CP10[2]*xint[ 26];
    yint[ 42]=CP00[ 7]*yint[ 38]+CP10[2]*yint[ 26];
    zint[ 42]=CP00[ 8]*zint[ 38]+CP10[2]*zint[ 26];
    CP10[3] += B00[3];
    xint[ 43]=CP00[ 9]*xint[ 39]+CP10[3]*xint[ 27];
    yint[ 43]=CP00[10]*yint[ 39]+CP10[3]*yint[ 27];
    zint[ 43]=CP00[11]*zint[ 39]+CP10[3]*zint[ 27];
    // (4,0)
    C10[0] += B10[0];
    xint[ 48]=C00[ 0]*xint[ 36]+C10[0]*xint[ 24];
    yint[ 48]=C00[ 1]*yint[ 36]+C10[0]*yint[ 24];
    zint[ 48]=C00[ 2]*zint[ 36]+C10[0]*zint[ 24];
    C10[1] += B10[1];
    xint[ 49]=C00[ 3]*xint[ 37]+C10[1]*xint[ 25];
    yint[ 49]=C00[ 4]*yint[ 37]+C10[1]*yint[ 25];
    zint[ 49]=C00[ 5]*zint[ 37]+C10[1]*zint[ 25];
    C10[2] += B10[2];
    xint[ 50]=C00[ 6]*xint[ 38]+C10[2]*xint[ 26];
    yint[ 50]=C00[ 7]*yint[ 38]+C10[2]*yint[ 26];
    zint[ 50]=C00[ 8]*zint[ 38]+C10[2]*zint[ 26];
    C10[3] += B10[3];
    xint[ 51]=C00[ 9]*xint[ 39]+C10[3]*xint[ 27];
    yint[ 51]=C00[10]*yint[ 39]+C10[3]*yint[ 27];
    zint[ 51]=C00[11]*zint[ 39]+C10[3]*zint[ 27];
    // (4,1)
    CP10[0] += B00[0];
    xint[ 52]=CP00[ 0]*xint[ 48]+CP10[0]*xint[ 36];
    yint[ 52]=CP00[ 1]*yint[ 48]+CP10[0]*yint[ 36];
    zint[ 52]=CP00[ 2]*zint[ 48]+CP10[0]*zint[ 36];
    CP10[1] += B00[1];
    xint[ 53]=CP00[ 3]*xint[ 49]+CP10[1]*xint[ 37];
    yint[ 53]=CP00[ 4]*yint[ 49]+CP10[1]*yint[ 37];
    zint[ 53]=CP00[ 5]*zint[ 49]+CP10[1]*zint[ 37];
    CP10[2] += B00[2];
    xint[ 54]=CP00[ 6]*xint[ 50]+CP10[2]*xint[ 38];
    yint[ 54]=CP00[ 7]*yint[ 50]+CP10[2]*yint[ 38];
    zint[ 54]=CP00[ 8]*zint[ 50]+CP10[2]*zint[ 38];
    CP10[3] += B00[3];
    xint[ 55]=CP00[ 9]*xint[ 51]+CP10[3]*xint[ 39];
    yint[ 55]=CP00[10]*yint[ 51]+CP10[3]*yint[ 39];
    zint[ 55]=CP00[11]*zint[ 51]+CP10[3]*zint[ 39];
    // (0,M) and (1,M)
    for ( m=0; m<4; m++ ) {
        CP01[m] = 0.e0;
        C01[m]  = B00[m];
    }

    // (0,2)
    CP01[0] += B01[0];
    xint[  8]=CP00[ 0]*xint[  4]+CP01[0]*xint[  0];
    yint[  8]=CP00[ 1]*yint[  4]+CP01[0]*yint[  0];
    zint[  8]=CP00[ 2]*zint[  4]+CP01[0]*zint[  0];
    CP01[1] += B01[1];
    xint[  9]=CP00[ 3]*xint[  5]+CP01[1]*xint[  1];
    yint[  9]=CP00[ 4]*yint[  5]+CP01[1]*yint[  1];
    zint[  9]=CP00[ 5]*zint[  5]+CP01[1]*zint[  1];
    CP01[2] += B01[2];
    xint[ 10]=CP00[ 6]*xint[  6]+CP01[2]*xint[  2];
    yint[ 10]=CP00[ 7]*yint[  6]+CP01[2]*yint[  2];
    zint[ 10]=CP00[ 8]*zint[  6]+CP01[2]*zint[  2];
    CP01[3] += B01[3];
    xint[ 11]=CP00[ 9]*xint[  7]+CP01[3]*xint[  3];
    yint[ 11]=CP00[10]*yint[  7]+CP01[3]*yint[  3];
    zint[ 11]=CP00[11]*zint[  7]+CP01[3]*zint[  3];
    // (1,2)
    C01[0] += B00[0];
    xint[ 20]=C00[ 0]*xint[  8]+C01[0]*xint[  4];
    yint[ 20]=C00[ 1]*yint[  8]+C01[0]*yint[  4];
    zint[ 20]=C00[ 2]*zint[  8]+C01[0]*zint[  4];
    C01[1] += B00[1];
    xint[ 21]=C00[ 3]*xint[  9]+C01[1]*xint[  5];
    yint[ 21]=C00[ 4]*yint[  9]+C01[1]*yint[  5];
    zint[ 21]=C00[ 5]*zint[  9]+C01[1]*zint[  5];
    C01[2] += B00[2];
    xint[ 22]=C00[ 6]*xint[ 10]+C01[2]*xint[  6];
    yint[ 22]=C00[ 7]*yint[ 10]+C01[2]*yint[  6];
    zint[ 22]=C00[ 8]*zint[ 10]+C01[2]*zint[  6];
    C01[3] += B00[3];
    xint[ 23]=C00[ 9]*xint[ 11]+C01[3]*xint[  7];
    yint[ 23]=C00[10]*yint[ 11]+C01[3]*yint[  7];
    zint[ 23]=C00[11]*zint[ 11]+C01[3]*zint[  7];
    // (N,M)
    for ( m=0; m<4; m++ ) C01[m] = B00[m];
    for ( m=0; m<4; m++ ) {
        C01[m] += B00[m];
        C10[m]  = B10[m];
    }
    // (2,2)
    xint[ 32]=C00[ 0]*xint[ 20]+C10[0]*xint[  8]+C01[0]*xint[ 16];
    yint[ 32]=C00[ 1]*yint[ 20]+C10[0]*yint[  8]+C01[0]*yint[ 16];
    zint[ 32]=C00[ 2]*zint[ 20]+C10[0]*zint[  8]+C01[0]*zint[ 16];
    C10[0] += B10[0];
    xint[ 33]=C00[ 3]*xint[ 21]+C10[1]*xint[  9]+C01[1]*xint[ 17];
    yint[ 33]=C00[ 4]*yint[ 21]+C10[1]*yint[  9]+C01[1]*yint[ 17];
    zint[ 33]=C00[ 5]*zint[ 21]+C10[1]*zint[  9]+C01[1]*zint[ 17];
    C10[1] += B10[1];
    xint[ 34]=C00[ 6]*xint[ 22]+C10[2]*xint[ 10]+C01[2]*xint[ 18];
    yint[ 34]=C00[ 7]*yint[ 22]+C10[2]*yint[ 10]+C01[2]*yint[ 18];
    zint[ 34]=C00[ 8]*zint[ 22]+C10[2]*zint[ 10]+C01[2]*zint[ 18];
    C10[2] += B10[2];
    xint[ 35]=C00[ 9]*xint[ 23]+C10[3]*xint[ 11]+C01[3]*xint[ 19];
    yint[ 35]=C00[10]*yint[ 23]+C10[3]*yint[ 11]+C01[3]*yint[ 19];
    zint[ 35]=C00[11]*zint[ 23]+C10[3]*zint[ 11]+C01[3]*zint[ 19];
    C10[3] += B10[3];
    // (3,2)
    xint[ 44]=C00[ 0]*xint[ 32]+C10[0]*xint[ 20]+C01[0]*xint[ 28];
    yint[ 44]=C00[ 1]*yint[ 32]+C10[0]*yint[ 20]+C01[0]*yint[ 28];
    zint[ 44]=C00[ 2]*zint[ 32]+C10[0]*zint[ 20]+C01[0]*zint[ 28];
    C10[0] += B10[0];
    xint[ 45]=C00[ 3]*xint[ 33]+C10[1]*xint[ 21]+C01[1]*xint[ 29];
    yint[ 45]=C00[ 4]*yint[ 33]+C10[1]*yint[ 21]+C01[1]*yint[ 29];
    zint[ 45]=C00[ 5]*zint[ 33]+C10[1]*zint[ 21]+C01[1]*zint[ 29];
    C10[1] += B10[1];
    xint[ 46]=C00[ 6]*xint[ 34]+C10[2]*xint[ 22]+C01[2]*xint[ 30];
    yint[ 46]=C00[ 7]*yint[ 34]+C10[2]*yint[ 22]+C01[2]*yint[ 30];
    zint[ 46]=C00[ 8]*zint[ 34]+C10[2]*zint[ 22]+C01[2]*zint[ 30];
    C10[2] += B10[2];
    xint[ 47]=C00[ 9]*xint[ 35]+C10[3]*xint[ 23]+C01[3]*xint[ 31];
    yint[ 47]=C00[10]*yint[ 35]+C10[3]*yint[ 23]+C01[3]*yint[ 31];
    zint[ 47]=C00[11]*zint[ 35]+C10[3]*zint[ 23]+C01[3]*zint[ 31];
    C10[3] += B10[3];
    // (4,2)
    xint[ 56]=C00[ 0]*xint[ 44]+C10[0]*xint[ 32]+C01[0]*xint[ 40];
    yint[ 56]=C00[ 1]*yint[ 44]+C10[0]*yint[ 32]+C01[0]*yint[ 40];
    zint[ 56]=C00[ 2]*zint[ 44]+C10[0]*zint[ 32]+C01[0]*zint[ 40];
    C10[0] += B10[0];
    xint[ 57]=C00[ 3]*xint[ 45]+C10[1]*xint[ 33]+C01[1]*xint[ 41];
    yint[ 57]=C00[ 4]*yint[ 45]+C10[1]*yint[ 33]+C01[1]*yint[ 41];
    zint[ 57]=C00[ 5]*zint[ 45]+C10[1]*zint[ 33]+C01[1]*zint[ 41];
    C10[1] += B10[1];
    xint[ 58]=C00[ 6]*xint[ 46]+C10[2]*xint[ 34]+C01[2]*xint[ 42];
    yint[ 58]=C00[ 7]*yint[ 46]+C10[2]*yint[ 34]+C01[2]*yint[ 42];
    zint[ 58]=C00[ 8]*zint[ 46]+C10[2]*zint[ 34]+C01[2]*zint[ 42];
    C10[2] += B10[2];
    xint[ 59]=C00[ 9]*xint[ 47]+C10[3]*xint[ 35]+C01[3]*xint[ 43];
    yint[ 59]=C00[10]*yint[ 47]+C10[3]*yint[ 35]+C01[3]*yint[ 43];
    zint[ 59]=C00[11]*zint[ 47]+C10[3]*zint[ 35]+C01[3]*zint[ 43];
    C10[3] += B10[3];
}

static void ofmo_form_ddds(
        const double *xint, const double *yint, const double *zint,
        double *eh ) {
    // (DS|DS)
    eh[   0] += xint[ 32]*yint[  0]*zint[  0];
    eh[   0] += xint[ 33]*yint[  1]*zint[  1];
    eh[   0] += xint[ 34]*yint[  2]*zint[  2];
    eh[   0] += xint[ 35]*yint[  3]*zint[  3];
    eh[   1] += xint[ 24]*yint[  8]*zint[  0];
    eh[   1] += xint[ 25]*yint[  9]*zint[  1];
    eh[   1] += xint[ 26]*yint[ 10]*zint[  2];
    eh[   1] += xint[ 27]*yint[ 11]*zint[  3];
    eh[   2] += xint[ 24]*yint[  0]*zint[  8];
    eh[   2] += xint[ 25]*yint[  1]*zint[  9];
    eh[   2] += xint[ 26]*yint[  2]*zint[ 10];
    eh[   2] += xint[ 27]*yint[  3]*zint[ 11];
    eh[   3] += xint[ 28]*yint[  4]*zint[  0];
    eh[   3] += xint[ 29]*yint[  5]*zint[  1];
    eh[   3] += xint[ 30]*yint[  6]*zint[  2];
    eh[   3] += xint[ 31]*yint[  7]*zint[  3];
    eh[   4] += xint[ 28]*yint[  0]*zint[  4];
    eh[   4] += xint[ 29]*yint[  1]*zint[  5];
    eh[   4] += xint[ 30]*yint[  2]*zint[  6];
    eh[   4] += xint[ 31]*yint[  3]*zint[  7];
    eh[   5] += xint[ 24]*yint[  4]*zint[  4];
    eh[   5] += xint[ 25]*yint[  5]*zint[  5];
    eh[   5] += xint[ 26]*yint[  6]*zint[  6];
    eh[   5] += xint[ 27]*yint[  7]*zint[  7];
    eh[   6] += xint[  8]*yint[ 24]*zint[  0];
    eh[   6] += xint[  9]*yint[ 25]*zint[  1];
    eh[   6] += xint[ 10]*yint[ 26]*zint[  2];
    eh[   6] += xint[ 11]*yint[ 27]*zint[  3];
    eh[   7] += xint[  0]*yint[ 32]*zint[  0];
    eh[   7] += xint[  1]*yint[ 33]*zint[  1];
    eh[   7] += xint[  2]*yint[ 34]*zint[  2];
    eh[   7] += xint[  3]*yint[ 35]*zint[  3];
    eh[   8] += xint[  0]*yint[ 24]*zint[  8];
    eh[   8] += xint[  1]*yint[ 25]*zint[  9];
    eh[   8] += xint[  2]*yint[ 26]*zint[ 10];
    eh[   8] += xint[  3]*yint[ 27]*zint[ 11];
    eh[   9] += xint[  4]*yint[ 28]*zint[  0];
    eh[   9] += xint[  5]*yint[ 29]*zint[  1];
    eh[   9] += xint[  6]*yint[ 30]*zint[  2];
    eh[   9] += xint[  7]*yint[ 31]*zint[  3];
    eh[  10] += xint[  4]*yint[ 24]*zint[  4];
    eh[  10] += xint[  5]*yint[ 25]*zint[  5];
    eh[  10] += xint[  6]*yint[ 26]*zint[  6];
    eh[  10] += xint[  7]*yint[ 27]*zint[  7];
    eh[  11] += xint[  0]*yint[ 28]*zint[  4];
    eh[  11] += xint[  1]*yint[ 29]*zint[  5];
    eh[  11] += xint[  2]*yint[ 30]*zint[  6];
    eh[  11] += xint[  3]*yint[ 31]*zint[  7];
    eh[  12] += xint[  8]*yint[  0]*zint[ 24];
    eh[  12] += xint[  9]*yint[  1]*zint[ 25];
    eh[  12] += xint[ 10]*yint[  2]*zint[ 26];
    eh[  12] += xint[ 11]*yint[  3]*zint[ 27];
    eh[  13] += xint[  0]*yint[  8]*zint[ 24];
    eh[  13] += xint[  1]*yint[  9]*zint[ 25];
    eh[  13] += xint[  2]*yint[ 10]*zint[ 26];
    eh[  13] += xint[  3]*yint[ 11]*zint[ 27];
    eh[  14] += xint[  0]*yint[  0]*zint[ 32];
    eh[  14] += xint[  1]*yint[  1]*zint[ 33];
    eh[  14] += xint[  2]*yint[  2]*zint[ 34];
    eh[  14] += xint[  3]*yint[  3]*zint[ 35];
    eh[  15] += xint[  4]*yint[  4]*zint[ 24];
    eh[  15] += xint[  5]*yint[  5]*zint[ 25];
    eh[  15] += xint[  6]*yint[  6]*zint[ 26];
    eh[  15] += xint[  7]*yint[  7]*zint[ 27];
    eh[  16] += xint[  4]*yint[  0]*zint[ 28];
    eh[  16] += xint[  5]*yint[  1]*zint[ 29];
    eh[  16] += xint[  6]*yint[  2]*zint[ 30];
    eh[  16] += xint[  7]*yint[  3]*zint[ 31];
    eh[  17] += xint[  0]*yint[  4]*zint[ 28];
    eh[  17] += xint[  1]*yint[  5]*zint[ 29];
    eh[  17] += xint[  2]*yint[  6]*zint[ 30];
    eh[  17] += xint[  3]*yint[  7]*zint[ 31];
    eh[  18] += xint[ 20]*yint[ 12]*zint[  0];
    eh[  18] += xint[ 21]*yint[ 13]*zint[  1];
    eh[  18] += xint[ 22]*yint[ 14]*zint[  2];
    eh[  18] += xint[ 23]*yint[ 15]*zint[  3];
    eh[  19] += xint[ 12]*yint[ 20]*zint[  0];
    eh[  19] += xint[ 13]*yint[ 21]*zint[  1];
    eh[  19] += xint[ 14]*yint[ 22]*zint[  2];
    eh[  19] += xint[ 15]*yint[ 23]*zint[  3];
    eh[  20] += xint[ 12]*yint[ 12]*zint[  8];
    eh[  20] += xint[ 13]*yint[ 13]*zint[  9];
    eh[  20] += xint[ 14]*yint[ 14]*zint[ 10];
    eh[  20] += xint[ 15]*yint[ 15]*zint[ 11];
    eh[  21] += xint[ 16]*yint[ 16]*zint[  0];
    eh[  21] += xint[ 17]*yint[ 17]*zint[  1];
    eh[  21] += xint[ 18]*yint[ 18]*zint[  2];
    eh[  21] += xint[ 19]*yint[ 19]*zint[  3];
    eh[  22] += xint[ 16]*yint[ 12]*zint[  4];
    eh[  22] += xint[ 17]*yint[ 13]*zint[  5];
    eh[  22] += xint[ 18]*yint[ 14]*zint[  6];
    eh[  22] += xint[ 19]*yint[ 15]*zint[  7];
    eh[  23] += xint[ 12]*yint[ 16]*zint[  4];
    eh[  23] += xint[ 13]*yint[ 17]*zint[  5];
    eh[  23] += xint[ 14]*yint[ 18]*zint[  6];
    eh[  23] += xint[ 15]*yint[ 19]*zint[  7];
    eh[  24] += xint[ 20]*yint[  0]*zint[ 12];
    eh[  24] += xint[ 21]*yint[  1]*zint[ 13];
    eh[  24] += xint[ 22]*yint[  2]*zint[ 14];
    eh[  24] += xint[ 23]*yint[  3]*zint[ 15];
    eh[  25] += xint[ 12]*yint[  8]*zint[ 12];
    eh[  25] += xint[ 13]*yint[  9]*zint[ 13];
    eh[  25] += xint[ 14]*yint[ 10]*zint[ 14];
    eh[  25] += xint[ 15]*yint[ 11]*zint[ 15];
    eh[  26] += xint[ 12]*yint[  0]*zint[ 20];
    eh[  26] += xint[ 13]*yint[  1]*zint[ 21];
    eh[  26] += xint[ 14]*yint[  2]*zint[ 22];
    eh[  26] += xint[ 15]*yint[  3]*zint[ 23];
    eh[  27] += xint[ 16]*yint[  4]*zint[ 12];
    eh[  27] += xint[ 17]*yint[  5]*zint[ 13];
    eh[  27] += xint[ 18]*yint[  6]*zint[ 14];
    eh[  27] += xint[ 19]*yint[  7]*zint[ 15];
    eh[  28] += xint[ 16]*yint[  0]*zint[ 16];
    eh[  28] += xint[ 17]*yint[  1]*zint[ 17];
    eh[  28] += xint[ 18]*yint[  2]*zint[ 18];
    eh[  28] += xint[ 19]*yint[  3]*zint[ 19];
    eh[  29] += xint[ 12]*yint[  4]*zint[ 16];
    eh[  29] += xint[ 13]*yint[  5]*zint[ 17];
    eh[  29] += xint[ 14]*yint[  6]*zint[ 18];
    eh[  29] += xint[ 15]*yint[  7]*zint[ 19];
    eh[  30] += xint[  8]*yint[ 12]*zint[ 12];
    eh[  30] += xint[  9]*yint[ 13]*zint[ 13];
    eh[  30] += xint[ 10]*yint[ 14]*zint[ 14];
    eh[  30] += xint[ 11]*yint[ 15]*zint[ 15];
    eh[  31] += xint[  0]*yint[ 20]*zint[ 12];
    eh[  31] += xint[  1]*yint[ 21]*zint[ 13];
    eh[  31] += xint[  2]*yint[ 22]*zint[ 14];
    eh[  31] += xint[  3]*yint[ 23]*zint[ 15];
    eh[  32] += xint[  0]*yint[ 12]*zint[ 20];
    eh[  32] += xint[  1]*yint[ 13]*zint[ 21];
    eh[  32] += xint[  2]*yint[ 14]*zint[ 22];
    eh[  32] += xint[  3]*yint[ 15]*zint[ 23];
    eh[  33] += xint[  4]*yint[ 16]*zint[ 12];
    eh[  33] += xint[  5]*yint[ 17]*zint[ 13];
    eh[  33] += xint[  6]*yint[ 18]*zint[ 14];
    eh[  33] += xint[  7]*yint[ 19]*zint[ 15];
    eh[  34] += xint[  4]*yint[ 12]*zint[ 16];
    eh[  34] += xint[  5]*yint[ 13]*zint[ 17];
    eh[  34] += xint[  6]*yint[ 14]*zint[ 18];
    eh[  34] += xint[  7]*yint[ 15]*zint[ 19];
    eh[  35] += xint[  0]*yint[ 16]*zint[ 16];
    eh[  35] += xint[  1]*yint[ 17]*zint[ 17];
    eh[  35] += xint[  2]*yint[ 18]*zint[ 18];
    eh[  35] += xint[  3]*yint[ 19]*zint[ 19];
    // (FS|DS)
    eh[  36] += xint[ 44]*yint[  0]*zint[  0];
    eh[  36] += xint[ 45]*yint[  1]*zint[  1];
    eh[  36] += xint[ 46]*yint[  2]*zint[  2];
    eh[  36] += xint[ 47]*yint[  3]*zint[  3];
    eh[  37] += xint[ 36]*yint[  8]*zint[  0];
    eh[  37] += xint[ 37]*yint[  9]*zint[  1];
    eh[  37] += xint[ 38]*yint[ 10]*zint[  2];
    eh[  37] += xint[ 39]*yint[ 11]*zint[  3];
    eh[  38] += xint[ 36]*yint[  0]*zint[  8];
    eh[  38] += xint[ 37]*yint[  1]*zint[  9];
    eh[  38] += xint[ 38]*yint[  2]*zint[ 10];
    eh[  38] += xint[ 39]*yint[  3]*zint[ 11];
    eh[  39] += xint[ 40]*yint[  4]*zint[  0];
    eh[  39] += xint[ 41]*yint[  5]*zint[  1];
    eh[  39] += xint[ 42]*yint[  6]*zint[  2];
    eh[  39] += xint[ 43]*yint[  7]*zint[  3];
    eh[  40] += xint[ 40]*yint[  0]*zint[  4];
    eh[  40] += xint[ 41]*yint[  1]*zint[  5];
    eh[  40] += xint[ 42]*yint[  2]*zint[  6];
    eh[  40] += xint[ 43]*yint[  3]*zint[  7];
    eh[  41] += xint[ 36]*yint[  4]*zint[  4];
    eh[  41] += xint[ 37]*yint[  5]*zint[  5];
    eh[  41] += xint[ 38]*yint[  6]*zint[  6];
    eh[  41] += xint[ 39]*yint[  7]*zint[  7];
    eh[  42] += xint[  8]*yint[ 36]*zint[  0];
    eh[  42] += xint[  9]*yint[ 37]*zint[  1];
    eh[  42] += xint[ 10]*yint[ 38]*zint[  2];
    eh[  42] += xint[ 11]*yint[ 39]*zint[  3];
    eh[  43] += xint[  0]*yint[ 44]*zint[  0];
    eh[  43] += xint[  1]*yint[ 45]*zint[  1];
    eh[  43] += xint[  2]*yint[ 46]*zint[  2];
    eh[  43] += xint[  3]*yint[ 47]*zint[  3];
    eh[  44] += xint[  0]*yint[ 36]*zint[  8];
    eh[  44] += xint[  1]*yint[ 37]*zint[  9];
    eh[  44] += xint[  2]*yint[ 38]*zint[ 10];
    eh[  44] += xint[  3]*yint[ 39]*zint[ 11];
    eh[  45] += xint[  4]*yint[ 40]*zint[  0];
    eh[  45] += xint[  5]*yint[ 41]*zint[  1];
    eh[  45] += xint[  6]*yint[ 42]*zint[  2];
    eh[  45] += xint[  7]*yint[ 43]*zint[  3];
    eh[  46] += xint[  4]*yint[ 36]*zint[  4];
    eh[  46] += xint[  5]*yint[ 37]*zint[  5];
    eh[  46] += xint[  6]*yint[ 38]*zint[  6];
    eh[  46] += xint[  7]*yint[ 39]*zint[  7];
    eh[  47] += xint[  0]*yint[ 40]*zint[  4];
    eh[  47] += xint[  1]*yint[ 41]*zint[  5];
    eh[  47] += xint[  2]*yint[ 42]*zint[  6];
    eh[  47] += xint[  3]*yint[ 43]*zint[  7];
    eh[  48] += xint[  8]*yint[  0]*zint[ 36];
    eh[  48] += xint[  9]*yint[  1]*zint[ 37];
    eh[  48] += xint[ 10]*yint[  2]*zint[ 38];
    eh[  48] += xint[ 11]*yint[  3]*zint[ 39];
    eh[  49] += xint[  0]*yint[  8]*zint[ 36];
    eh[  49] += xint[  1]*yint[  9]*zint[ 37];
    eh[  49] += xint[  2]*yint[ 10]*zint[ 38];
    eh[  49] += xint[  3]*yint[ 11]*zint[ 39];
    eh[  50] += xint[  0]*yint[  0]*zint[ 44];
    eh[  50] += xint[  1]*yint[  1]*zint[ 45];
    eh[  50] += xint[  2]*yint[  2]*zint[ 46];
    eh[  50] += xint[  3]*yint[  3]*zint[ 47];
    eh[  51] += xint[  4]*yint[  4]*zint[ 36];
    eh[  51] += xint[  5]*yint[  5]*zint[ 37];
    eh[  51] += xint[  6]*yint[  6]*zint[ 38];
    eh[  51] += xint[  7]*yint[  7]*zint[ 39];
    eh[  52] += xint[  4]*yint[  0]*zint[ 40];
    eh[  52] += xint[  5]*yint[  1]*zint[ 41];
    eh[  52] += xint[  6]*yint[  2]*zint[ 42];
    eh[  52] += xint[  7]*yint[  3]*zint[ 43];
    eh[  53] += xint[  0]*yint[  4]*zint[ 40];
    eh[  53] += xint[  1]*yint[  5]*zint[ 41];
    eh[  53] += xint[  2]*yint[  6]*zint[ 42];
    eh[  53] += xint[  3]*yint[  7]*zint[ 43];
    eh[  54] += xint[ 32]*yint[ 12]*zint[  0];
    eh[  54] += xint[ 33]*yint[ 13]*zint[  1];
    eh[  54] += xint[ 34]*yint[ 14]*zint[  2];
    eh[  54] += xint[ 35]*yint[ 15]*zint[  3];
    eh[  55] += xint[ 24]*yint[ 20]*zint[  0];
    eh[  55] += xint[ 25]*yint[ 21]*zint[  1];
    eh[  55] += xint[ 26]*yint[ 22]*zint[  2];
    eh[  55] += xint[ 27]*yint[ 23]*zint[  3];
    eh[  56] += xint[ 24]*yint[ 12]*zint[  8];
    eh[  56] += xint[ 25]*yint[ 13]*zint[  9];
    eh[  56] += xint[ 26]*yint[ 14]*zint[ 10];
    eh[  56] += xint[ 27]*yint[ 15]*zint[ 11];
    eh[  57] += xint[ 28]*yint[ 16]*zint[  0];
    eh[  57] += xint[ 29]*yint[ 17]*zint[  1];
    eh[  57] += xint[ 30]*yint[ 18]*zint[  2];
    eh[  57] += xint[ 31]*yint[ 19]*zint[  3];
    eh[  58] += xint[ 28]*yint[ 12]*zint[  4];
    eh[  58] += xint[ 29]*yint[ 13]*zint[  5];
    eh[  58] += xint[ 30]*yint[ 14]*zint[  6];
    eh[  58] += xint[ 31]*yint[ 15]*zint[  7];
    eh[  59] += xint[ 24]*yint[ 16]*zint[  4];
    eh[  59] += xint[ 25]*yint[ 17]*zint[  5];
    eh[  59] += xint[ 26]*yint[ 18]*zint[  6];
    eh[  59] += xint[ 27]*yint[ 19]*zint[  7];
    eh[  60] += xint[ 32]*yint[  0]*zint[ 12];
    eh[  60] += xint[ 33]*yint[  1]*zint[ 13];
    eh[  60] += xint[ 34]*yint[  2]*zint[ 14];
    eh[  60] += xint[ 35]*yint[  3]*zint[ 15];
    eh[  61] += xint[ 24]*yint[  8]*zint[ 12];
    eh[  61] += xint[ 25]*yint[  9]*zint[ 13];
    eh[  61] += xint[ 26]*yint[ 10]*zint[ 14];
    eh[  61] += xint[ 27]*yint[ 11]*zint[ 15];
    eh[  62] += xint[ 24]*yint[  0]*zint[ 20];
    eh[  62] += xint[ 25]*yint[  1]*zint[ 21];
    eh[  62] += xint[ 26]*yint[  2]*zint[ 22];
    eh[  62] += xint[ 27]*yint[  3]*zint[ 23];
    eh[  63] += xint[ 28]*yint[  4]*zint[ 12];
    eh[  63] += xint[ 29]*yint[  5]*zint[ 13];
    eh[  63] += xint[ 30]*yint[  6]*zint[ 14];
    eh[  63] += xint[ 31]*yint[  7]*zint[ 15];
    eh[  64] += xint[ 28]*yint[  0]*zint[ 16];
    eh[  64] += xint[ 29]*yint[  1]*zint[ 17];
    eh[  64] += xint[ 30]*yint[  2]*zint[ 18];
    eh[  64] += xint[ 31]*yint[  3]*zint[ 19];
    eh[  65] += xint[ 24]*yint[  4]*zint[ 16];
    eh[  65] += xint[ 25]*yint[  5]*zint[ 17];
    eh[  65] += xint[ 26]*yint[  6]*zint[ 18];
    eh[  65] += xint[ 27]*yint[  7]*zint[ 19];
    eh[  66] += xint[ 20]*yint[ 24]*zint[  0];
    eh[  66] += xint[ 21]*yint[ 25]*zint[  1];
    eh[  66] += xint[ 22]*yint[ 26]*zint[  2];
    eh[  66] += xint[ 23]*yint[ 27]*zint[  3];
    eh[  67] += xint[ 12]*yint[ 32]*zint[  0];
    eh[  67] += xint[ 13]*yint[ 33]*zint[  1];
    eh[  67] += xint[ 14]*yint[ 34]*zint[  2];
    eh[  67] += xint[ 15]*yint[ 35]*zint[  3];
    eh[  68] += xint[ 12]*yint[ 24]*zint[  8];
    eh[  68] += xint[ 13]*yint[ 25]*zint[  9];
    eh[  68] += xint[ 14]*yint[ 26]*zint[ 10];
    eh[  68] += xint[ 15]*yint[ 27]*zint[ 11];
    eh[  69] += xint[ 16]*yint[ 28]*zint[  0];
    eh[  69] += xint[ 17]*yint[ 29]*zint[  1];
    eh[  69] += xint[ 18]*yint[ 30]*zint[  2];
    eh[  69] += xint[ 19]*yint[ 31]*zint[  3];
    eh[  70] += xint[ 16]*yint[ 24]*zint[  4];
    eh[  70] += xint[ 17]*yint[ 25]*zint[  5];
    eh[  70] += xint[ 18]*yint[ 26]*zint[  6];
    eh[  70] += xint[ 19]*yint[ 27]*zint[  7];
    eh[  71] += xint[ 12]*yint[ 28]*zint[  4];
    eh[  71] += xint[ 13]*yint[ 29]*zint[  5];
    eh[  71] += xint[ 14]*yint[ 30]*zint[  6];
    eh[  71] += xint[ 15]*yint[ 31]*zint[  7];
    eh[  72] += xint[ 20]*yint[  0]*zint[ 24];
    eh[  72] += xint[ 21]*yint[  1]*zint[ 25];
    eh[  72] += xint[ 22]*yint[  2]*zint[ 26];
    eh[  72] += xint[ 23]*yint[  3]*zint[ 27];
    eh[  73] += xint[ 12]*yint[  8]*zint[ 24];
    eh[  73] += xint[ 13]*yint[  9]*zint[ 25];
    eh[  73] += xint[ 14]*yint[ 10]*zint[ 26];
    eh[  73] += xint[ 15]*yint[ 11]*zint[ 27];
    eh[  74] += xint[ 12]*yint[  0]*zint[ 32];
    eh[  74] += xint[ 13]*yint[  1]*zint[ 33];
    eh[  74] += xint[ 14]*yint[  2]*zint[ 34];
    eh[  74] += xint[ 15]*yint[  3]*zint[ 35];
    eh[  75] += xint[ 16]*yint[  4]*zint[ 24];
    eh[  75] += xint[ 17]*yint[  5]*zint[ 25];
    eh[  75] += xint[ 18]*yint[  6]*zint[ 26];
    eh[  75] += xint[ 19]*yint[  7]*zint[ 27];
    eh[  76] += xint[ 16]*yint[  0]*zint[ 28];
    eh[  76] += xint[ 17]*yint[  1]*zint[ 29];
    eh[  76] += xint[ 18]*yint[  2]*zint[ 30];
    eh[  76] += xint[ 19]*yint[  3]*zint[ 31];
    eh[  77] += xint[ 12]*yint[  4]*zint[ 28];
    eh[  77] += xint[ 13]*yint[  5]*zint[ 29];
    eh[  77] += xint[ 14]*yint[  6]*zint[ 30];
    eh[  77] += xint[ 15]*yint[  7]*zint[ 31];
    eh[  78] += xint[ 20]*yint[ 12]*zint[ 12];
    eh[  78] += xint[ 21]*yint[ 13]*zint[ 13];
    eh[  78] += xint[ 22]*yint[ 14]*zint[ 14];
    eh[  78] += xint[ 23]*yint[ 15]*zint[ 15];
    eh[  79] += xint[ 12]*yint[ 20]*zint[ 12];
    eh[  79] += xint[ 13]*yint[ 21]*zint[ 13];
    eh[  79] += xint[ 14]*yint[ 22]*zint[ 14];
    eh[  79] += xint[ 15]*yint[ 23]*zint[ 15];
    eh[  80] += xint[ 12]*yint[ 12]*zint[ 20];
    eh[  80] += xint[ 13]*yint[ 13]*zint[ 21];
    eh[  80] += xint[ 14]*yint[ 14]*zint[ 22];
    eh[  80] += xint[ 15]*yint[ 15]*zint[ 23];
    eh[  81] += xint[ 16]*yint[ 16]*zint[ 12];
    eh[  81] += xint[ 17]*yint[ 17]*zint[ 13];
    eh[  81] += xint[ 18]*yint[ 18]*zint[ 14];
    eh[  81] += xint[ 19]*yint[ 19]*zint[ 15];
    eh[  82] += xint[ 16]*yint[ 12]*zint[ 16];
    eh[  82] += xint[ 17]*yint[ 13]*zint[ 17];
    eh[  82] += xint[ 18]*yint[ 14]*zint[ 18];
    eh[  82] += xint[ 19]*yint[ 15]*zint[ 19];
    eh[  83] += xint[ 12]*yint[ 16]*zint[ 16];
    eh[  83] += xint[ 13]*yint[ 17]*zint[ 17];
    eh[  83] += xint[ 14]*yint[ 18]*zint[ 18];
    eh[  83] += xint[ 15]*yint[ 19]*zint[ 19];
    eh[  84] += xint[  8]*yint[ 24]*zint[ 12];
    eh[  84] += xint[  9]*yint[ 25]*zint[ 13];
    eh[  84] += xint[ 10]*yint[ 26]*zint[ 14];
    eh[  84] += xint[ 11]*yint[ 27]*zint[ 15];
    eh[  85] += xint[  0]*yint[ 32]*zint[ 12];
    eh[  85] += xint[  1]*yint[ 33]*zint[ 13];
    eh[  85] += xint[  2]*yint[ 34]*zint[ 14];
    eh[  85] += xint[  3]*yint[ 35]*zint[ 15];
    eh[  86] += xint[  0]*yint[ 24]*zint[ 20];
    eh[  86] += xint[  1]*yint[ 25]*zint[ 21];
    eh[  86] += xint[  2]*yint[ 26]*zint[ 22];
    eh[  86] += xint[  3]*yint[ 27]*zint[ 23];
    eh[  87] += xint[  4]*yint[ 28]*zint[ 12];
    eh[  87] += xint[  5]*yint[ 29]*zint[ 13];
    eh[  87] += xint[  6]*yint[ 30]*zint[ 14];
    eh[  87] += xint[  7]*yint[ 31]*zint[ 15];
    eh[  88] += xint[  4]*yint[ 24]*zint[ 16];
    eh[  88] += xint[  5]*yint[ 25]*zint[ 17];
    eh[  88] += xint[  6]*yint[ 26]*zint[ 18];
    eh[  88] += xint[  7]*yint[ 27]*zint[ 19];
    eh[  89] += xint[  0]*yint[ 28]*zint[ 16];
    eh[  89] += xint[  1]*yint[ 29]*zint[ 17];
    eh[  89] += xint[  2]*yint[ 30]*zint[ 18];
    eh[  89] += xint[  3]*yint[ 31]*zint[ 19];
    eh[  90] += xint[  8]*yint[ 12]*zint[ 24];
    eh[  90] += xint[  9]*yint[ 13]*zint[ 25];
    eh[  90] += xint[ 10]*yint[ 14]*zint[ 26];
    eh[  90] += xint[ 11]*yint[ 15]*zint[ 27];
    eh[  91] += xint[  0]*yint[ 20]*zint[ 24];
    eh[  91] += xint[  1]*yint[ 21]*zint[ 25];
    eh[  91] += xint[  2]*yint[ 22]*zint[ 26];
    eh[  91] += xint[  3]*yint[ 23]*zint[ 27];
    eh[  92] += xint[  0]*yint[ 12]*zint[ 32];
    eh[  92] += xint[  1]*yint[ 13]*zint[ 33];
    eh[  92] += xint[  2]*yint[ 14]*zint[ 34];
    eh[  92] += xint[  3]*yint[ 15]*zint[ 35];
    eh[  93] += xint[  4]*yint[ 16]*zint[ 24];
    eh[  93] += xint[  5]*yint[ 17]*zint[ 25];
    eh[  93] += xint[  6]*yint[ 18]*zint[ 26];
    eh[  93] += xint[  7]*yint[ 19]*zint[ 27];
    eh[  94] += xint[  4]*yint[ 12]*zint[ 28];
    eh[  94] += xint[  5]*yint[ 13]*zint[ 29];
    eh[  94] += xint[  6]*yint[ 14]*zint[ 30];
    eh[  94] += xint[  7]*yint[ 15]*zint[ 31];
    eh[  95] += xint[  0]*yint[ 16]*zint[ 28];
    eh[  95] += xint[  1]*yint[ 17]*zint[ 29];
    eh[  95] += xint[  2]*yint[ 18]*zint[ 30];
    eh[  95] += xint[  3]*yint[ 19]*zint[ 31];
    // (GS|DS)
    eh[  96] += xint[ 56]*yint[  0]*zint[  0];
    eh[  96] += xint[ 57]*yint[  1]*zint[  1];
    eh[  96] += xint[ 58]*yint[  2]*zint[  2];
    eh[  96] += xint[ 59]*yint[  3]*zint[  3];
    eh[  97] += xint[ 48]*yint[  8]*zint[  0];
    eh[  97] += xint[ 49]*yint[  9]*zint[  1];
    eh[  97] += xint[ 50]*yint[ 10]*zint[  2];
    eh[  97] += xint[ 51]*yint[ 11]*zint[  3];
    eh[  98] += xint[ 48]*yint[  0]*zint[  8];
    eh[  98] += xint[ 49]*yint[  1]*zint[  9];
    eh[  98] += xint[ 50]*yint[  2]*zint[ 10];
    eh[  98] += xint[ 51]*yint[  3]*zint[ 11];
    eh[  99] += xint[ 52]*yint[  4]*zint[  0];
    eh[  99] += xint[ 53]*yint[  5]*zint[  1];
    eh[  99] += xint[ 54]*yint[  6]*zint[  2];
    eh[  99] += xint[ 55]*yint[  7]*zint[  3];
    eh[ 100] += xint[ 52]*yint[  0]*zint[  4];
    eh[ 100] += xint[ 53]*yint[  1]*zint[  5];
    eh[ 100] += xint[ 54]*yint[  2]*zint[  6];
    eh[ 100] += xint[ 55]*yint[  3]*zint[  7];
    eh[ 101] += xint[ 48]*yint[  4]*zint[  4];
    eh[ 101] += xint[ 49]*yint[  5]*zint[  5];
    eh[ 101] += xint[ 50]*yint[  6]*zint[  6];
    eh[ 101] += xint[ 51]*yint[  7]*zint[  7];
    eh[ 102] += xint[  8]*yint[ 48]*zint[  0];
    eh[ 102] += xint[  9]*yint[ 49]*zint[  1];
    eh[ 102] += xint[ 10]*yint[ 50]*zint[  2];
    eh[ 102] += xint[ 11]*yint[ 51]*zint[  3];
    eh[ 103] += xint[  0]*yint[ 56]*zint[  0];
    eh[ 103] += xint[  1]*yint[ 57]*zint[  1];
    eh[ 103] += xint[  2]*yint[ 58]*zint[  2];
    eh[ 103] += xint[  3]*yint[ 59]*zint[  3];
    eh[ 104] += xint[  0]*yint[ 48]*zint[  8];
    eh[ 104] += xint[  1]*yint[ 49]*zint[  9];
    eh[ 104] += xint[  2]*yint[ 50]*zint[ 10];
    eh[ 104] += xint[  3]*yint[ 51]*zint[ 11];
    eh[ 105] += xint[  4]*yint[ 52]*zint[  0];
    eh[ 105] += xint[  5]*yint[ 53]*zint[  1];
    eh[ 105] += xint[  6]*yint[ 54]*zint[  2];
    eh[ 105] += xint[  7]*yint[ 55]*zint[  3];
    eh[ 106] += xint[  4]*yint[ 48]*zint[  4];
    eh[ 106] += xint[  5]*yint[ 49]*zint[  5];
    eh[ 106] += xint[  6]*yint[ 50]*zint[  6];
    eh[ 106] += xint[  7]*yint[ 51]*zint[  7];
    eh[ 107] += xint[  0]*yint[ 52]*zint[  4];
    eh[ 107] += xint[  1]*yint[ 53]*zint[  5];
    eh[ 107] += xint[  2]*yint[ 54]*zint[  6];
    eh[ 107] += xint[  3]*yint[ 55]*zint[  7];
    eh[ 108] += xint[  8]*yint[  0]*zint[ 48];
    eh[ 108] += xint[  9]*yint[  1]*zint[ 49];
    eh[ 108] += xint[ 10]*yint[  2]*zint[ 50];
    eh[ 108] += xint[ 11]*yint[  3]*zint[ 51];
    eh[ 109] += xint[  0]*yint[  8]*zint[ 48];
    eh[ 109] += xint[  1]*yint[  9]*zint[ 49];
    eh[ 109] += xint[  2]*yint[ 10]*zint[ 50];
    eh[ 109] += xint[  3]*yint[ 11]*zint[ 51];
    eh[ 110] += xint[  0]*yint[  0]*zint[ 56];
    eh[ 110] += xint[  1]*yint[  1]*zint[ 57];
    eh[ 110] += xint[  2]*yint[  2]*zint[ 58];
    eh[ 110] += xint[  3]*yint[  3]*zint[ 59];
    eh[ 111] += xint[  4]*yint[  4]*zint[ 48];
    eh[ 111] += xint[  5]*yint[  5]*zint[ 49];
    eh[ 111] += xint[  6]*yint[  6]*zint[ 50];
    eh[ 111] += xint[  7]*yint[  7]*zint[ 51];
    eh[ 112] += xint[  4]*yint[  0]*zint[ 52];
    eh[ 112] += xint[  5]*yint[  1]*zint[ 53];
    eh[ 112] += xint[  6]*yint[  2]*zint[ 54];
    eh[ 112] += xint[  7]*yint[  3]*zint[ 55];
    eh[ 113] += xint[  0]*yint[  4]*zint[ 52];
    eh[ 113] += xint[  1]*yint[  5]*zint[ 53];
    eh[ 113] += xint[  2]*yint[  6]*zint[ 54];
    eh[ 113] += xint[  3]*yint[  7]*zint[ 55];
    eh[ 114] += xint[ 44]*yint[ 12]*zint[  0];
    eh[ 114] += xint[ 45]*yint[ 13]*zint[  1];
    eh[ 114] += xint[ 46]*yint[ 14]*zint[  2];
    eh[ 114] += xint[ 47]*yint[ 15]*zint[  3];
    eh[ 115] += xint[ 36]*yint[ 20]*zint[  0];
    eh[ 115] += xint[ 37]*yint[ 21]*zint[  1];
    eh[ 115] += xint[ 38]*yint[ 22]*zint[  2];
    eh[ 115] += xint[ 39]*yint[ 23]*zint[  3];
    eh[ 116] += xint[ 36]*yint[ 12]*zint[  8];
    eh[ 116] += xint[ 37]*yint[ 13]*zint[  9];
    eh[ 116] += xint[ 38]*yint[ 14]*zint[ 10];
    eh[ 116] += xint[ 39]*yint[ 15]*zint[ 11];
    eh[ 117] += xint[ 40]*yint[ 16]*zint[  0];
    eh[ 117] += xint[ 41]*yint[ 17]*zint[  1];
    eh[ 117] += xint[ 42]*yint[ 18]*zint[  2];
    eh[ 117] += xint[ 43]*yint[ 19]*zint[  3];
    eh[ 118] += xint[ 40]*yint[ 12]*zint[  4];
    eh[ 118] += xint[ 41]*yint[ 13]*zint[  5];
    eh[ 118] += xint[ 42]*yint[ 14]*zint[  6];
    eh[ 118] += xint[ 43]*yint[ 15]*zint[  7];
    eh[ 119] += xint[ 36]*yint[ 16]*zint[  4];
    eh[ 119] += xint[ 37]*yint[ 17]*zint[  5];
    eh[ 119] += xint[ 38]*yint[ 18]*zint[  6];
    eh[ 119] += xint[ 39]*yint[ 19]*zint[  7];
    eh[ 120] += xint[ 44]*yint[  0]*zint[ 12];
    eh[ 120] += xint[ 45]*yint[  1]*zint[ 13];
    eh[ 120] += xint[ 46]*yint[  2]*zint[ 14];
    eh[ 120] += xint[ 47]*yint[  3]*zint[ 15];
    eh[ 121] += xint[ 36]*yint[  8]*zint[ 12];
    eh[ 121] += xint[ 37]*yint[  9]*zint[ 13];
    eh[ 121] += xint[ 38]*yint[ 10]*zint[ 14];
    eh[ 121] += xint[ 39]*yint[ 11]*zint[ 15];
    eh[ 122] += xint[ 36]*yint[  0]*zint[ 20];
    eh[ 122] += xint[ 37]*yint[  1]*zint[ 21];
    eh[ 122] += xint[ 38]*yint[  2]*zint[ 22];
    eh[ 122] += xint[ 39]*yint[  3]*zint[ 23];
    eh[ 123] += xint[ 40]*yint[  4]*zint[ 12];
    eh[ 123] += xint[ 41]*yint[  5]*zint[ 13];
    eh[ 123] += xint[ 42]*yint[  6]*zint[ 14];
    eh[ 123] += xint[ 43]*yint[  7]*zint[ 15];
    eh[ 124] += xint[ 40]*yint[  0]*zint[ 16];
    eh[ 124] += xint[ 41]*yint[  1]*zint[ 17];
    eh[ 124] += xint[ 42]*yint[  2]*zint[ 18];
    eh[ 124] += xint[ 43]*yint[  3]*zint[ 19];
    eh[ 125] += xint[ 36]*yint[  4]*zint[ 16];
    eh[ 125] += xint[ 37]*yint[  5]*zint[ 17];
    eh[ 125] += xint[ 38]*yint[  6]*zint[ 18];
    eh[ 125] += xint[ 39]*yint[  7]*zint[ 19];
    eh[ 126] += xint[ 32]*yint[ 24]*zint[  0];
    eh[ 126] += xint[ 33]*yint[ 25]*zint[  1];
    eh[ 126] += xint[ 34]*yint[ 26]*zint[  2];
    eh[ 126] += xint[ 35]*yint[ 27]*zint[  3];
    eh[ 127] += xint[ 24]*yint[ 32]*zint[  0];
    eh[ 127] += xint[ 25]*yint[ 33]*zint[  1];
    eh[ 127] += xint[ 26]*yint[ 34]*zint[  2];
    eh[ 127] += xint[ 27]*yint[ 35]*zint[  3];
    eh[ 128] += xint[ 24]*yint[ 24]*zint[  8];
    eh[ 128] += xint[ 25]*yint[ 25]*zint[  9];
    eh[ 128] += xint[ 26]*yint[ 26]*zint[ 10];
    eh[ 128] += xint[ 27]*yint[ 27]*zint[ 11];
    eh[ 129] += xint[ 28]*yint[ 28]*zint[  0];
    eh[ 129] += xint[ 29]*yint[ 29]*zint[  1];
    eh[ 129] += xint[ 30]*yint[ 30]*zint[  2];
    eh[ 129] += xint[ 31]*yint[ 31]*zint[  3];
    eh[ 130] += xint[ 28]*yint[ 24]*zint[  4];
    eh[ 130] += xint[ 29]*yint[ 25]*zint[  5];
    eh[ 130] += xint[ 30]*yint[ 26]*zint[  6];
    eh[ 130] += xint[ 31]*yint[ 27]*zint[  7];
    eh[ 131] += xint[ 24]*yint[ 28]*zint[  4];
    eh[ 131] += xint[ 25]*yint[ 29]*zint[  5];
    eh[ 131] += xint[ 26]*yint[ 30]*zint[  6];
    eh[ 131] += xint[ 27]*yint[ 31]*zint[  7];
    eh[ 132] += xint[ 32]*yint[  0]*zint[ 24];
    eh[ 132] += xint[ 33]*yint[  1]*zint[ 25];
    eh[ 132] += xint[ 34]*yint[  2]*zint[ 26];
    eh[ 132] += xint[ 35]*yint[  3]*zint[ 27];
    eh[ 133] += xint[ 24]*yint[  8]*zint[ 24];
    eh[ 133] += xint[ 25]*yint[  9]*zint[ 25];
    eh[ 133] += xint[ 26]*yint[ 10]*zint[ 26];
    eh[ 133] += xint[ 27]*yint[ 11]*zint[ 27];
    eh[ 134] += xint[ 24]*yint[  0]*zint[ 32];
    eh[ 134] += xint[ 25]*yint[  1]*zint[ 33];
    eh[ 134] += xint[ 26]*yint[  2]*zint[ 34];
    eh[ 134] += xint[ 27]*yint[  3]*zint[ 35];
    eh[ 135] += xint[ 28]*yint[  4]*zint[ 24];
    eh[ 135] += xint[ 29]*yint[  5]*zint[ 25];
    eh[ 135] += xint[ 30]*yint[  6]*zint[ 26];
    eh[ 135] += xint[ 31]*yint[  7]*zint[ 27];
    eh[ 136] += xint[ 28]*yint[  0]*zint[ 28];
    eh[ 136] += xint[ 29]*yint[  1]*zint[ 29];
    eh[ 136] += xint[ 30]*yint[  2]*zint[ 30];
    eh[ 136] += xint[ 31]*yint[  3]*zint[ 31];
    eh[ 137] += xint[ 24]*yint[  4]*zint[ 28];
    eh[ 137] += xint[ 25]*yint[  5]*zint[ 29];
    eh[ 137] += xint[ 26]*yint[  6]*zint[ 30];
    eh[ 137] += xint[ 27]*yint[  7]*zint[ 31];
    eh[ 138] += xint[ 32]*yint[ 12]*zint[ 12];
    eh[ 138] += xint[ 33]*yint[ 13]*zint[ 13];
    eh[ 138] += xint[ 34]*yint[ 14]*zint[ 14];
    eh[ 138] += xint[ 35]*yint[ 15]*zint[ 15];
    eh[ 139] += xint[ 24]*yint[ 20]*zint[ 12];
    eh[ 139] += xint[ 25]*yint[ 21]*zint[ 13];
    eh[ 139] += xint[ 26]*yint[ 22]*zint[ 14];
    eh[ 139] += xint[ 27]*yint[ 23]*zint[ 15];
    eh[ 140] += xint[ 24]*yint[ 12]*zint[ 20];
    eh[ 140] += xint[ 25]*yint[ 13]*zint[ 21];
    eh[ 140] += xint[ 26]*yint[ 14]*zint[ 22];
    eh[ 140] += xint[ 27]*yint[ 15]*zint[ 23];
    eh[ 141] += xint[ 28]*yint[ 16]*zint[ 12];
    eh[ 141] += xint[ 29]*yint[ 17]*zint[ 13];
    eh[ 141] += xint[ 30]*yint[ 18]*zint[ 14];
    eh[ 141] += xint[ 31]*yint[ 19]*zint[ 15];
    eh[ 142] += xint[ 28]*yint[ 12]*zint[ 16];
    eh[ 142] += xint[ 29]*yint[ 13]*zint[ 17];
    eh[ 142] += xint[ 30]*yint[ 14]*zint[ 18];
    eh[ 142] += xint[ 31]*yint[ 15]*zint[ 19];
    eh[ 143] += xint[ 24]*yint[ 16]*zint[ 16];
    eh[ 143] += xint[ 25]*yint[ 17]*zint[ 17];
    eh[ 143] += xint[ 26]*yint[ 18]*zint[ 18];
    eh[ 143] += xint[ 27]*yint[ 19]*zint[ 19];
    eh[ 144] += xint[ 20]*yint[ 36]*zint[  0];
    eh[ 144] += xint[ 21]*yint[ 37]*zint[  1];
    eh[ 144] += xint[ 22]*yint[ 38]*zint[  2];
    eh[ 144] += xint[ 23]*yint[ 39]*zint[  3];
    eh[ 145] += xint[ 12]*yint[ 44]*zint[  0];
    eh[ 145] += xint[ 13]*yint[ 45]*zint[  1];
    eh[ 145] += xint[ 14]*yint[ 46]*zint[  2];
    eh[ 145] += xint[ 15]*yint[ 47]*zint[  3];
    eh[ 146] += xint[ 12]*yint[ 36]*zint[  8];
    eh[ 146] += xint[ 13]*yint[ 37]*zint[  9];
    eh[ 146] += xint[ 14]*yint[ 38]*zint[ 10];
    eh[ 146] += xint[ 15]*yint[ 39]*zint[ 11];
    eh[ 147] += xint[ 16]*yint[ 40]*zint[  0];
    eh[ 147] += xint[ 17]*yint[ 41]*zint[  1];
    eh[ 147] += xint[ 18]*yint[ 42]*zint[  2];
    eh[ 147] += xint[ 19]*yint[ 43]*zint[  3];
    eh[ 148] += xint[ 16]*yint[ 36]*zint[  4];
    eh[ 148] += xint[ 17]*yint[ 37]*zint[  5];
    eh[ 148] += xint[ 18]*yint[ 38]*zint[  6];
    eh[ 148] += xint[ 19]*yint[ 39]*zint[  7];
    eh[ 149] += xint[ 12]*yint[ 40]*zint[  4];
    eh[ 149] += xint[ 13]*yint[ 41]*zint[  5];
    eh[ 149] += xint[ 14]*yint[ 42]*zint[  6];
    eh[ 149] += xint[ 15]*yint[ 43]*zint[  7];
    eh[ 150] += xint[ 20]*yint[  0]*zint[ 36];
    eh[ 150] += xint[ 21]*yint[  1]*zint[ 37];
    eh[ 150] += xint[ 22]*yint[  2]*zint[ 38];
    eh[ 150] += xint[ 23]*yint[  3]*zint[ 39];
    eh[ 151] += xint[ 12]*yint[  8]*zint[ 36];
    eh[ 151] += xint[ 13]*yint[  9]*zint[ 37];
    eh[ 151] += xint[ 14]*yint[ 10]*zint[ 38];
    eh[ 151] += xint[ 15]*yint[ 11]*zint[ 39];
    eh[ 152] += xint[ 12]*yint[  0]*zint[ 44];
    eh[ 152] += xint[ 13]*yint[  1]*zint[ 45];
    eh[ 152] += xint[ 14]*yint[  2]*zint[ 46];
    eh[ 152] += xint[ 15]*yint[  3]*zint[ 47];
    eh[ 153] += xint[ 16]*yint[  4]*zint[ 36];
    eh[ 153] += xint[ 17]*yint[  5]*zint[ 37];
    eh[ 153] += xint[ 18]*yint[  6]*zint[ 38];
    eh[ 153] += xint[ 19]*yint[  7]*zint[ 39];
    eh[ 154] += xint[ 16]*yint[  0]*zint[ 40];
    eh[ 154] += xint[ 17]*yint[  1]*zint[ 41];
    eh[ 154] += xint[ 18]*yint[  2]*zint[ 42];
    eh[ 154] += xint[ 19]*yint[  3]*zint[ 43];
    eh[ 155] += xint[ 12]*yint[  4]*zint[ 40];
    eh[ 155] += xint[ 13]*yint[  5]*zint[ 41];
    eh[ 155] += xint[ 14]*yint[  6]*zint[ 42];
    eh[ 155] += xint[ 15]*yint[  7]*zint[ 43];
    eh[ 156] += xint[ 20]*yint[ 24]*zint[ 12];
    eh[ 156] += xint[ 21]*yint[ 25]*zint[ 13];
    eh[ 156] += xint[ 22]*yint[ 26]*zint[ 14];
    eh[ 156] += xint[ 23]*yint[ 27]*zint[ 15];
    eh[ 157] += xint[ 12]*yint[ 32]*zint[ 12];
    eh[ 157] += xint[ 13]*yint[ 33]*zint[ 13];
    eh[ 157] += xint[ 14]*yint[ 34]*zint[ 14];
    eh[ 157] += xint[ 15]*yint[ 35]*zint[ 15];
    eh[ 158] += xint[ 12]*yint[ 24]*zint[ 20];
    eh[ 158] += xint[ 13]*yint[ 25]*zint[ 21];
    eh[ 158] += xint[ 14]*yint[ 26]*zint[ 22];
    eh[ 158] += xint[ 15]*yint[ 27]*zint[ 23];
    eh[ 159] += xint[ 16]*yint[ 28]*zint[ 12];
    eh[ 159] += xint[ 17]*yint[ 29]*zint[ 13];
    eh[ 159] += xint[ 18]*yint[ 30]*zint[ 14];
    eh[ 159] += xint[ 19]*yint[ 31]*zint[ 15];
    eh[ 160] += xint[ 16]*yint[ 24]*zint[ 16];
    eh[ 160] += xint[ 17]*yint[ 25]*zint[ 17];
    eh[ 160] += xint[ 18]*yint[ 26]*zint[ 18];
    eh[ 160] += xint[ 19]*yint[ 27]*zint[ 19];
    eh[ 161] += xint[ 12]*yint[ 28]*zint[ 16];
    eh[ 161] += xint[ 13]*yint[ 29]*zint[ 17];
    eh[ 161] += xint[ 14]*yint[ 30]*zint[ 18];
    eh[ 161] += xint[ 15]*yint[ 31]*zint[ 19];
    eh[ 162] += xint[ 20]*yint[ 12]*zint[ 24];
    eh[ 162] += xint[ 21]*yint[ 13]*zint[ 25];
    eh[ 162] += xint[ 22]*yint[ 14]*zint[ 26];
    eh[ 162] += xint[ 23]*yint[ 15]*zint[ 27];
    eh[ 163] += xint[ 12]*yint[ 20]*zint[ 24];
    eh[ 163] += xint[ 13]*yint[ 21]*zint[ 25];
    eh[ 163] += xint[ 14]*yint[ 22]*zint[ 26];
    eh[ 163] += xint[ 15]*yint[ 23]*zint[ 27];
    eh[ 164] += xint[ 12]*yint[ 12]*zint[ 32];
    eh[ 164] += xint[ 13]*yint[ 13]*zint[ 33];
    eh[ 164] += xint[ 14]*yint[ 14]*zint[ 34];
    eh[ 164] += xint[ 15]*yint[ 15]*zint[ 35];
    eh[ 165] += xint[ 16]*yint[ 16]*zint[ 24];
    eh[ 165] += xint[ 17]*yint[ 17]*zint[ 25];
    eh[ 165] += xint[ 18]*yint[ 18]*zint[ 26];
    eh[ 165] += xint[ 19]*yint[ 19]*zint[ 27];
    eh[ 166] += xint[ 16]*yint[ 12]*zint[ 28];
    eh[ 166] += xint[ 17]*yint[ 13]*zint[ 29];
    eh[ 166] += xint[ 18]*yint[ 14]*zint[ 30];
    eh[ 166] += xint[ 19]*yint[ 15]*zint[ 31];
    eh[ 167] += xint[ 12]*yint[ 16]*zint[ 28];
    eh[ 167] += xint[ 13]*yint[ 17]*zint[ 29];
    eh[ 167] += xint[ 14]*yint[ 18]*zint[ 30];
    eh[ 167] += xint[ 15]*yint[ 19]*zint[ 31];
    eh[ 168] += xint[  8]*yint[ 36]*zint[ 12];
    eh[ 168] += xint[  9]*yint[ 37]*zint[ 13];
    eh[ 168] += xint[ 10]*yint[ 38]*zint[ 14];
    eh[ 168] += xint[ 11]*yint[ 39]*zint[ 15];
    eh[ 169] += xint[  0]*yint[ 44]*zint[ 12];
    eh[ 169] += xint[  1]*yint[ 45]*zint[ 13];
    eh[ 169] += xint[  2]*yint[ 46]*zint[ 14];
    eh[ 169] += xint[  3]*yint[ 47]*zint[ 15];
    eh[ 170] += xint[  0]*yint[ 36]*zint[ 20];
    eh[ 170] += xint[  1]*yint[ 37]*zint[ 21];
    eh[ 170] += xint[  2]*yint[ 38]*zint[ 22];
    eh[ 170] += xint[  3]*yint[ 39]*zint[ 23];
    eh[ 171] += xint[  4]*yint[ 40]*zint[ 12];
    eh[ 171] += xint[  5]*yint[ 41]*zint[ 13];
    eh[ 171] += xint[  6]*yint[ 42]*zint[ 14];
    eh[ 171] += xint[  7]*yint[ 43]*zint[ 15];
    eh[ 172] += xint[  4]*yint[ 36]*zint[ 16];
    eh[ 172] += xint[  5]*yint[ 37]*zint[ 17];
    eh[ 172] += xint[  6]*yint[ 38]*zint[ 18];
    eh[ 172] += xint[  7]*yint[ 39]*zint[ 19];
    eh[ 173] += xint[  0]*yint[ 40]*zint[ 16];
    eh[ 173] += xint[  1]*yint[ 41]*zint[ 17];
    eh[ 173] += xint[  2]*yint[ 42]*zint[ 18];
    eh[ 173] += xint[  3]*yint[ 43]*zint[ 19];
    eh[ 174] += xint[  8]*yint[ 12]*zint[ 36];
    eh[ 174] += xint[  9]*yint[ 13]*zint[ 37];
    eh[ 174] += xint[ 10]*yint[ 14]*zint[ 38];
    eh[ 174] += xint[ 11]*yint[ 15]*zint[ 39];
    eh[ 175] += xint[  0]*yint[ 20]*zint[ 36];
    eh[ 175] += xint[  1]*yint[ 21]*zint[ 37];
    eh[ 175] += xint[  2]*yint[ 22]*zint[ 38];
    eh[ 175] += xint[  3]*yint[ 23]*zint[ 39];
    eh[ 176] += xint[  0]*yint[ 12]*zint[ 44];
    eh[ 176] += xint[  1]*yint[ 13]*zint[ 45];
    eh[ 176] += xint[  2]*yint[ 14]*zint[ 46];
    eh[ 176] += xint[  3]*yint[ 15]*zint[ 47];
    eh[ 177] += xint[  4]*yint[ 16]*zint[ 36];
    eh[ 177] += xint[  5]*yint[ 17]*zint[ 37];
    eh[ 177] += xint[  6]*yint[ 18]*zint[ 38];
    eh[ 177] += xint[  7]*yint[ 19]*zint[ 39];
    eh[ 178] += xint[  4]*yint[ 12]*zint[ 40];
    eh[ 178] += xint[  5]*yint[ 13]*zint[ 41];
    eh[ 178] += xint[  6]*yint[ 14]*zint[ 42];
    eh[ 178] += xint[  7]*yint[ 15]*zint[ 43];
    eh[ 179] += xint[  0]*yint[ 16]*zint[ 40];
    eh[ 179] += xint[  1]*yint[ 17]*zint[ 41];
    eh[ 179] += xint[  2]*yint[ 18]*zint[ 42];
    eh[ 179] += xint[  3]*yint[ 19]*zint[ 43];
    eh[ 180] += xint[  8]*yint[ 24]*zint[ 24];
    eh[ 180] += xint[  9]*yint[ 25]*zint[ 25];
    eh[ 180] += xint[ 10]*yint[ 26]*zint[ 26];
    eh[ 180] += xint[ 11]*yint[ 27]*zint[ 27];
    eh[ 181] += xint[  0]*yint[ 32]*zint[ 24];
    eh[ 181] += xint[  1]*yint[ 33]*zint[ 25];
    eh[ 181] += xint[  2]*yint[ 34]*zint[ 26];
    eh[ 181] += xint[  3]*yint[ 35]*zint[ 27];
    eh[ 182] += xint[  0]*yint[ 24]*zint[ 32];
    eh[ 182] += xint[  1]*yint[ 25]*zint[ 33];
    eh[ 182] += xint[  2]*yint[ 26]*zint[ 34];
    eh[ 182] += xint[  3]*yint[ 27]*zint[ 35];
    eh[ 183] += xint[  4]*yint[ 28]*zint[ 24];
    eh[ 183] += xint[  5]*yint[ 29]*zint[ 25];
    eh[ 183] += xint[  6]*yint[ 30]*zint[ 26];
    eh[ 183] += xint[  7]*yint[ 31]*zint[ 27];
    eh[ 184] += xint[  4]*yint[ 24]*zint[ 28];
    eh[ 184] += xint[  5]*yint[ 25]*zint[ 29];
    eh[ 184] += xint[  6]*yint[ 26]*zint[ 30];
    eh[ 184] += xint[  7]*yint[ 27]*zint[ 31];
    eh[ 185] += xint[  0]*yint[ 28]*zint[ 28];
    eh[ 185] += xint[  1]*yint[ 29]*zint[ 29];
    eh[ 185] += xint[  2]*yint[ 30]*zint[ 30];
    eh[ 185] += xint[  3]*yint[ 31]*zint[ 31];
}

void ofmo_twoint_core_rys_ddds( const int mythread,
        const int *nijps, const double *vzeta, const double *vdkab,
        const double vxiza[], const double BA[3],
        const int *nklps, const double *veta, const double *vdkcd,
        const double *vxizc, const double DC[3], const double AC[3],
        double *DINT ) {
    int ijps, klps, i;
    double cssss, zeta, dkab, xiza, eta, xizc, dk, T;
    double zeta2, eta2, rz, PA[3], QC[3];
    double PQ2, sqrho, rho, PC[3], QP[3];
    double C00[12], CP00[12], B00[4], B10[4], B01[4], F00[4];
    double rrho, rze, W[13], U[13];
    double u2, duminv, dm2inv, dum;
    int m, m3;
    double *xint, *yint, *zint, *eh;
    xint = ofmo_integ_getadd_xint( mythread );
    yint = ofmo_integ_getadd_yint( mythread );
    zint = ofmo_integ_getadd_zint( mythread );
    eh   = ofmo_integ_getadd_eh( mythread );
    DFACT = ofmo_getadd_dfact();
    ofmo_hrr_clear_ddds( eh );
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
            calc_root( 4, T, U, W );
            for ( m=m3=0; m<4; m++, m3+=3 ) {
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
            ofmo_xyzint_ddds( F00, B00, B10, B01, C00, CP00,
                     xint, yint, zint );
            ofmo_form_ddds( xint, yint, zint, eh );
        }
    }
    ofmo_hrr_calc_ddds( eh, BA, DC );
    ofmo_hrr_coef_ddds( eh, DINT );
}

int ofmo_twoint_rys_ddds(
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
    max_nzeri = ebuf_max_nzeri - 6*6*6*1;
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
            ofmo_twoint_core_rys_ddds( mythread,
                    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
                    &psp_xiza[ijps0], BA,
                    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
                    &psp_xiza[klps0], DC,   AC,      DINTEG );
            ipat=((Lab != Lcd) || (ics==kcs && jcs>lcs) ? true : false);
            for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
                I2 = (iao*iao+iao)>>1;
                for ( j=0, jao=jao0; j<6; j++, jao++ ) {
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

int ofmo_twoint_direct_rys_ddds(
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
    max_nzeri -= 6*6*6*1;
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
            ofmo_twoint_core_rys_ddds( mythread,
                    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
                    &psp_xiza[ijps0], BA,
                    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
                    &psp_xiza[klps0], DC,   AC,      DINTEG );
            ipat=((Lab != Lcd) || (ics==kcs && jcs>lcs) ? true : false);
            for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
                I2 = (iao*iao+iao)>>1;
                for ( j=0, jao=jao0; j<6; j++, jao++ ) {
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
