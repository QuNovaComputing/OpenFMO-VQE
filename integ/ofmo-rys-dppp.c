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

static void ofmo_hrr_clear_dppp( double *eh ) {
    int i;
    // (DS|PS)
    for ( i=0; i<(0+18); i++ ) eh[i] = 0.e0;
    // (DS|DS)
    for ( i=18; i<(18+36); i++ ) eh[i] = 0.e0;
    // (FS|PS)
    for ( i=54; i<(54+30); i++ ) eh[i] = 0.e0;
    // (FS|DS)
    for ( i=84; i<(84+60); i++ ) eh[i] = 0.e0;
}

static void ofmo_hrr_coef_dppp(
        double *eh, double *DINT ) {
    int i, j, k, l, iao, jao, kao, lao, ix;
    double coef_a, coef_ab, coef_abc;
    double *th;
    th = &eh[306];
    ix = 0;
    for ( i=0, iao=4; i<6; i++, iao++ ) {
        coef_a = DFACT[iao];
        for ( j=0, jao=1; j<3; j++, jao++ ) {
            coef_ab = coef_a * DFACT[jao];
            for ( k=0, kao=1; k<3; k++, kao++ ) {
                coef_abc = coef_ab * DFACT[kao];
                for ( l=0, lao=1; l<3; l++, lao++ ) {
                    DINT[ix] = coef_abc * DFACT[lao] * th[ix];
                    ix++;
                }
            }
        }
    }
}

static void ofmo_hrr_calc_dppp( double *eh,
        const double BA[3], const double DC[3] ) {
    // HRR for (XX|XS)-type integral (center AB)
    // (DP,PS)
    eh[ 144] = eh[  54] - BA[0]*eh[   0];
    eh[ 145] = eh[  55] - BA[0]*eh[   1];
    eh[ 146] = eh[  56] - BA[0]*eh[   2];
    eh[ 147] = eh[  63] - BA[1]*eh[   0];
    eh[ 148] = eh[  64] - BA[1]*eh[   1];
    eh[ 149] = eh[  65] - BA[1]*eh[   2];
    eh[ 150] = eh[  66] - BA[2]*eh[   0];
    eh[ 151] = eh[  67] - BA[2]*eh[   1];
    eh[ 152] = eh[  68] - BA[2]*eh[   2];
    eh[ 153] = eh[  69] - BA[0]*eh[   3];
    eh[ 154] = eh[  70] - BA[0]*eh[   4];
    eh[ 155] = eh[  71] - BA[0]*eh[   5];
    eh[ 156] = eh[  57] - BA[1]*eh[   3];
    eh[ 157] = eh[  58] - BA[1]*eh[   4];
    eh[ 158] = eh[  59] - BA[1]*eh[   5];
    eh[ 159] = eh[  78] - BA[2]*eh[   3];
    eh[ 160] = eh[  79] - BA[2]*eh[   4];
    eh[ 161] = eh[  80] - BA[2]*eh[   5];
    eh[ 162] = eh[  72] - BA[0]*eh[   6];
    eh[ 163] = eh[  73] - BA[0]*eh[   7];
    eh[ 164] = eh[  74] - BA[0]*eh[   8];
    eh[ 165] = eh[  81] - BA[1]*eh[   6];
    eh[ 166] = eh[  82] - BA[1]*eh[   7];
    eh[ 167] = eh[  83] - BA[1]*eh[   8];
    eh[ 168] = eh[  60] - BA[2]*eh[   6];
    eh[ 169] = eh[  61] - BA[2]*eh[   7];
    eh[ 170] = eh[  62] - BA[2]*eh[   8];
    eh[ 171] = eh[  63] - BA[0]*eh[   9];
    eh[ 172] = eh[  64] - BA[0]*eh[  10];
    eh[ 173] = eh[  65] - BA[0]*eh[  11];
    eh[ 174] = eh[  69] - BA[1]*eh[   9];
    eh[ 175] = eh[  70] - BA[1]*eh[  10];
    eh[ 176] = eh[  71] - BA[1]*eh[  11];
    eh[ 177] = eh[  75] - BA[2]*eh[   9];
    eh[ 178] = eh[  76] - BA[2]*eh[  10];
    eh[ 179] = eh[  77] - BA[2]*eh[  11];
    eh[ 180] = eh[  66] - BA[0]*eh[  12];
    eh[ 181] = eh[  67] - BA[0]*eh[  13];
    eh[ 182] = eh[  68] - BA[0]*eh[  14];
    eh[ 183] = eh[  75] - BA[1]*eh[  12];
    eh[ 184] = eh[  76] - BA[1]*eh[  13];
    eh[ 185] = eh[  77] - BA[1]*eh[  14];
    eh[ 186] = eh[  72] - BA[2]*eh[  12];
    eh[ 187] = eh[  73] - BA[2]*eh[  13];
    eh[ 188] = eh[  74] - BA[2]*eh[  14];
    eh[ 189] = eh[  75] - BA[0]*eh[  15];
    eh[ 190] = eh[  76] - BA[0]*eh[  16];
    eh[ 191] = eh[  77] - BA[0]*eh[  17];
    eh[ 192] = eh[  78] - BA[1]*eh[  15];
    eh[ 193] = eh[  79] - BA[1]*eh[  16];
    eh[ 194] = eh[  80] - BA[1]*eh[  17];
    eh[ 195] = eh[  81] - BA[2]*eh[  15];
    eh[ 196] = eh[  82] - BA[2]*eh[  16];
    eh[ 197] = eh[  83] - BA[2]*eh[  17];
    // (DP,DS)
    eh[ 198] = eh[  84] - BA[0]*eh[  18];
    eh[ 199] = eh[  85] - BA[0]*eh[  19];
    eh[ 200] = eh[  86] - BA[0]*eh[  20];
    eh[ 201] = eh[  87] - BA[0]*eh[  21];
    eh[ 202] = eh[  88] - BA[0]*eh[  22];
    eh[ 203] = eh[  89] - BA[0]*eh[  23];
    eh[ 204] = eh[ 102] - BA[1]*eh[  18];
    eh[ 205] = eh[ 103] - BA[1]*eh[  19];
    eh[ 206] = eh[ 104] - BA[1]*eh[  20];
    eh[ 207] = eh[ 105] - BA[1]*eh[  21];
    eh[ 208] = eh[ 106] - BA[1]*eh[  22];
    eh[ 209] = eh[ 107] - BA[1]*eh[  23];
    eh[ 210] = eh[ 108] - BA[2]*eh[  18];
    eh[ 211] = eh[ 109] - BA[2]*eh[  19];
    eh[ 212] = eh[ 110] - BA[2]*eh[  20];
    eh[ 213] = eh[ 111] - BA[2]*eh[  21];
    eh[ 214] = eh[ 112] - BA[2]*eh[  22];
    eh[ 215] = eh[ 113] - BA[2]*eh[  23];
    eh[ 216] = eh[ 114] - BA[0]*eh[  24];
    eh[ 217] = eh[ 115] - BA[0]*eh[  25];
    eh[ 218] = eh[ 116] - BA[0]*eh[  26];
    eh[ 219] = eh[ 117] - BA[0]*eh[  27];
    eh[ 220] = eh[ 118] - BA[0]*eh[  28];
    eh[ 221] = eh[ 119] - BA[0]*eh[  29];
    eh[ 222] = eh[  90] - BA[1]*eh[  24];
    eh[ 223] = eh[  91] - BA[1]*eh[  25];
    eh[ 224] = eh[  92] - BA[1]*eh[  26];
    eh[ 225] = eh[  93] - BA[1]*eh[  27];
    eh[ 226] = eh[  94] - BA[1]*eh[  28];
    eh[ 227] = eh[  95] - BA[1]*eh[  29];
    eh[ 228] = eh[ 132] - BA[2]*eh[  24];
    eh[ 229] = eh[ 133] - BA[2]*eh[  25];
    eh[ 230] = eh[ 134] - BA[2]*eh[  26];
    eh[ 231] = eh[ 135] - BA[2]*eh[  27];
    eh[ 232] = eh[ 136] - BA[2]*eh[  28];
    eh[ 233] = eh[ 137] - BA[2]*eh[  29];
    eh[ 234] = eh[ 120] - BA[0]*eh[  30];
    eh[ 235] = eh[ 121] - BA[0]*eh[  31];
    eh[ 236] = eh[ 122] - BA[0]*eh[  32];
    eh[ 237] = eh[ 123] - BA[0]*eh[  33];
    eh[ 238] = eh[ 124] - BA[0]*eh[  34];
    eh[ 239] = eh[ 125] - BA[0]*eh[  35];
    eh[ 240] = eh[ 138] - BA[1]*eh[  30];
    eh[ 241] = eh[ 139] - BA[1]*eh[  31];
    eh[ 242] = eh[ 140] - BA[1]*eh[  32];
    eh[ 243] = eh[ 141] - BA[1]*eh[  33];
    eh[ 244] = eh[ 142] - BA[1]*eh[  34];
    eh[ 245] = eh[ 143] - BA[1]*eh[  35];
    eh[ 246] = eh[  96] - BA[2]*eh[  30];
    eh[ 247] = eh[  97] - BA[2]*eh[  31];
    eh[ 248] = eh[  98] - BA[2]*eh[  32];
    eh[ 249] = eh[  99] - BA[2]*eh[  33];
    eh[ 250] = eh[ 100] - BA[2]*eh[  34];
    eh[ 251] = eh[ 101] - BA[2]*eh[  35];
    eh[ 252] = eh[ 102] - BA[0]*eh[  36];
    eh[ 253] = eh[ 103] - BA[0]*eh[  37];
    eh[ 254] = eh[ 104] - BA[0]*eh[  38];
    eh[ 255] = eh[ 105] - BA[0]*eh[  39];
    eh[ 256] = eh[ 106] - BA[0]*eh[  40];
    eh[ 257] = eh[ 107] - BA[0]*eh[  41];
    eh[ 258] = eh[ 114] - BA[1]*eh[  36];
    eh[ 259] = eh[ 115] - BA[1]*eh[  37];
    eh[ 260] = eh[ 116] - BA[1]*eh[  38];
    eh[ 261] = eh[ 117] - BA[1]*eh[  39];
    eh[ 262] = eh[ 118] - BA[1]*eh[  40];
    eh[ 263] = eh[ 119] - BA[1]*eh[  41];
    eh[ 264] = eh[ 126] - BA[2]*eh[  36];
    eh[ 265] = eh[ 127] - BA[2]*eh[  37];
    eh[ 266] = eh[ 128] - BA[2]*eh[  38];
    eh[ 267] = eh[ 129] - BA[2]*eh[  39];
    eh[ 268] = eh[ 130] - BA[2]*eh[  40];
    eh[ 269] = eh[ 131] - BA[2]*eh[  41];
    eh[ 270] = eh[ 108] - BA[0]*eh[  42];
    eh[ 271] = eh[ 109] - BA[0]*eh[  43];
    eh[ 272] = eh[ 110] - BA[0]*eh[  44];
    eh[ 273] = eh[ 111] - BA[0]*eh[  45];
    eh[ 274] = eh[ 112] - BA[0]*eh[  46];
    eh[ 275] = eh[ 113] - BA[0]*eh[  47];
    eh[ 276] = eh[ 126] - BA[1]*eh[  42];
    eh[ 277] = eh[ 127] - BA[1]*eh[  43];
    eh[ 278] = eh[ 128] - BA[1]*eh[  44];
    eh[ 279] = eh[ 129] - BA[1]*eh[  45];
    eh[ 280] = eh[ 130] - BA[1]*eh[  46];
    eh[ 281] = eh[ 131] - BA[1]*eh[  47];
    eh[ 282] = eh[ 120] - BA[2]*eh[  42];
    eh[ 283] = eh[ 121] - BA[2]*eh[  43];
    eh[ 284] = eh[ 122] - BA[2]*eh[  44];
    eh[ 285] = eh[ 123] - BA[2]*eh[  45];
    eh[ 286] = eh[ 124] - BA[2]*eh[  46];
    eh[ 287] = eh[ 125] - BA[2]*eh[  47];
    eh[ 288] = eh[ 126] - BA[0]*eh[  48];
    eh[ 289] = eh[ 127] - BA[0]*eh[  49];
    eh[ 290] = eh[ 128] - BA[0]*eh[  50];
    eh[ 291] = eh[ 129] - BA[0]*eh[  51];
    eh[ 292] = eh[ 130] - BA[0]*eh[  52];
    eh[ 293] = eh[ 131] - BA[0]*eh[  53];
    eh[ 294] = eh[ 132] - BA[1]*eh[  48];
    eh[ 295] = eh[ 133] - BA[1]*eh[  49];
    eh[ 296] = eh[ 134] - BA[1]*eh[  50];
    eh[ 297] = eh[ 135] - BA[1]*eh[  51];
    eh[ 298] = eh[ 136] - BA[1]*eh[  52];
    eh[ 299] = eh[ 137] - BA[1]*eh[  53];
    eh[ 300] = eh[ 138] - BA[2]*eh[  48];
    eh[ 301] = eh[ 139] - BA[2]*eh[  49];
    eh[ 302] = eh[ 140] - BA[2]*eh[  50];
    eh[ 303] = eh[ 141] - BA[2]*eh[  51];
    eh[ 304] = eh[ 142] - BA[2]*eh[  52];
    eh[ 305] = eh[ 143] - BA[2]*eh[  53];
    // (DP,PP)
    eh[ 306] = eh[ 198] - DC[0]*eh[ 144];
    eh[ 307] = eh[ 201] - DC[1]*eh[ 144];
    eh[ 308] = eh[ 202] - DC[2]*eh[ 144];
    eh[ 309] = eh[ 201] - DC[0]*eh[ 145];
    eh[ 310] = eh[ 199] - DC[1]*eh[ 145];
    eh[ 311] = eh[ 203] - DC[2]*eh[ 145];
    eh[ 312] = eh[ 202] - DC[0]*eh[ 146];
    eh[ 313] = eh[ 203] - DC[1]*eh[ 146];
    eh[ 314] = eh[ 200] - DC[2]*eh[ 146];
    eh[ 315] = eh[ 204] - DC[0]*eh[ 147];
    eh[ 316] = eh[ 207] - DC[1]*eh[ 147];
    eh[ 317] = eh[ 208] - DC[2]*eh[ 147];
    eh[ 318] = eh[ 207] - DC[0]*eh[ 148];
    eh[ 319] = eh[ 205] - DC[1]*eh[ 148];
    eh[ 320] = eh[ 209] - DC[2]*eh[ 148];
    eh[ 321] = eh[ 208] - DC[0]*eh[ 149];
    eh[ 322] = eh[ 209] - DC[1]*eh[ 149];
    eh[ 323] = eh[ 206] - DC[2]*eh[ 149];
    eh[ 324] = eh[ 210] - DC[0]*eh[ 150];
    eh[ 325] = eh[ 213] - DC[1]*eh[ 150];
    eh[ 326] = eh[ 214] - DC[2]*eh[ 150];
    eh[ 327] = eh[ 213] - DC[0]*eh[ 151];
    eh[ 328] = eh[ 211] - DC[1]*eh[ 151];
    eh[ 329] = eh[ 215] - DC[2]*eh[ 151];
    eh[ 330] = eh[ 214] - DC[0]*eh[ 152];
    eh[ 331] = eh[ 215] - DC[1]*eh[ 152];
    eh[ 332] = eh[ 212] - DC[2]*eh[ 152];
    eh[ 333] = eh[ 216] - DC[0]*eh[ 153];
    eh[ 334] = eh[ 219] - DC[1]*eh[ 153];
    eh[ 335] = eh[ 220] - DC[2]*eh[ 153];
    eh[ 336] = eh[ 219] - DC[0]*eh[ 154];
    eh[ 337] = eh[ 217] - DC[1]*eh[ 154];
    eh[ 338] = eh[ 221] - DC[2]*eh[ 154];
    eh[ 339] = eh[ 220] - DC[0]*eh[ 155];
    eh[ 340] = eh[ 221] - DC[1]*eh[ 155];
    eh[ 341] = eh[ 218] - DC[2]*eh[ 155];
    eh[ 342] = eh[ 222] - DC[0]*eh[ 156];
    eh[ 343] = eh[ 225] - DC[1]*eh[ 156];
    eh[ 344] = eh[ 226] - DC[2]*eh[ 156];
    eh[ 345] = eh[ 225] - DC[0]*eh[ 157];
    eh[ 346] = eh[ 223] - DC[1]*eh[ 157];
    eh[ 347] = eh[ 227] - DC[2]*eh[ 157];
    eh[ 348] = eh[ 226] - DC[0]*eh[ 158];
    eh[ 349] = eh[ 227] - DC[1]*eh[ 158];
    eh[ 350] = eh[ 224] - DC[2]*eh[ 158];
    eh[ 351] = eh[ 228] - DC[0]*eh[ 159];
    eh[ 352] = eh[ 231] - DC[1]*eh[ 159];
    eh[ 353] = eh[ 232] - DC[2]*eh[ 159];
    eh[ 354] = eh[ 231] - DC[0]*eh[ 160];
    eh[ 355] = eh[ 229] - DC[1]*eh[ 160];
    eh[ 356] = eh[ 233] - DC[2]*eh[ 160];
    eh[ 357] = eh[ 232] - DC[0]*eh[ 161];
    eh[ 358] = eh[ 233] - DC[1]*eh[ 161];
    eh[ 359] = eh[ 230] - DC[2]*eh[ 161];
    eh[ 360] = eh[ 234] - DC[0]*eh[ 162];
    eh[ 361] = eh[ 237] - DC[1]*eh[ 162];
    eh[ 362] = eh[ 238] - DC[2]*eh[ 162];
    eh[ 363] = eh[ 237] - DC[0]*eh[ 163];
    eh[ 364] = eh[ 235] - DC[1]*eh[ 163];
    eh[ 365] = eh[ 239] - DC[2]*eh[ 163];
    eh[ 366] = eh[ 238] - DC[0]*eh[ 164];
    eh[ 367] = eh[ 239] - DC[1]*eh[ 164];
    eh[ 368] = eh[ 236] - DC[2]*eh[ 164];
    eh[ 369] = eh[ 240] - DC[0]*eh[ 165];
    eh[ 370] = eh[ 243] - DC[1]*eh[ 165];
    eh[ 371] = eh[ 244] - DC[2]*eh[ 165];
    eh[ 372] = eh[ 243] - DC[0]*eh[ 166];
    eh[ 373] = eh[ 241] - DC[1]*eh[ 166];
    eh[ 374] = eh[ 245] - DC[2]*eh[ 166];
    eh[ 375] = eh[ 244] - DC[0]*eh[ 167];
    eh[ 376] = eh[ 245] - DC[1]*eh[ 167];
    eh[ 377] = eh[ 242] - DC[2]*eh[ 167];
    eh[ 378] = eh[ 246] - DC[0]*eh[ 168];
    eh[ 379] = eh[ 249] - DC[1]*eh[ 168];
    eh[ 380] = eh[ 250] - DC[2]*eh[ 168];
    eh[ 381] = eh[ 249] - DC[0]*eh[ 169];
    eh[ 382] = eh[ 247] - DC[1]*eh[ 169];
    eh[ 383] = eh[ 251] - DC[2]*eh[ 169];
    eh[ 384] = eh[ 250] - DC[0]*eh[ 170];
    eh[ 385] = eh[ 251] - DC[1]*eh[ 170];
    eh[ 386] = eh[ 248] - DC[2]*eh[ 170];
    eh[ 387] = eh[ 252] - DC[0]*eh[ 171];
    eh[ 388] = eh[ 255] - DC[1]*eh[ 171];
    eh[ 389] = eh[ 256] - DC[2]*eh[ 171];
    eh[ 390] = eh[ 255] - DC[0]*eh[ 172];
    eh[ 391] = eh[ 253] - DC[1]*eh[ 172];
    eh[ 392] = eh[ 257] - DC[2]*eh[ 172];
    eh[ 393] = eh[ 256] - DC[0]*eh[ 173];
    eh[ 394] = eh[ 257] - DC[1]*eh[ 173];
    eh[ 395] = eh[ 254] - DC[2]*eh[ 173];
    eh[ 396] = eh[ 258] - DC[0]*eh[ 174];
    eh[ 397] = eh[ 261] - DC[1]*eh[ 174];
    eh[ 398] = eh[ 262] - DC[2]*eh[ 174];
    eh[ 399] = eh[ 261] - DC[0]*eh[ 175];
    eh[ 400] = eh[ 259] - DC[1]*eh[ 175];
    eh[ 401] = eh[ 263] - DC[2]*eh[ 175];
    eh[ 402] = eh[ 262] - DC[0]*eh[ 176];
    eh[ 403] = eh[ 263] - DC[1]*eh[ 176];
    eh[ 404] = eh[ 260] - DC[2]*eh[ 176];
    eh[ 405] = eh[ 264] - DC[0]*eh[ 177];
    eh[ 406] = eh[ 267] - DC[1]*eh[ 177];
    eh[ 407] = eh[ 268] - DC[2]*eh[ 177];
    eh[ 408] = eh[ 267] - DC[0]*eh[ 178];
    eh[ 409] = eh[ 265] - DC[1]*eh[ 178];
    eh[ 410] = eh[ 269] - DC[2]*eh[ 178];
    eh[ 411] = eh[ 268] - DC[0]*eh[ 179];
    eh[ 412] = eh[ 269] - DC[1]*eh[ 179];
    eh[ 413] = eh[ 266] - DC[2]*eh[ 179];
    eh[ 414] = eh[ 270] - DC[0]*eh[ 180];
    eh[ 415] = eh[ 273] - DC[1]*eh[ 180];
    eh[ 416] = eh[ 274] - DC[2]*eh[ 180];
    eh[ 417] = eh[ 273] - DC[0]*eh[ 181];
    eh[ 418] = eh[ 271] - DC[1]*eh[ 181];
    eh[ 419] = eh[ 275] - DC[2]*eh[ 181];
    eh[ 420] = eh[ 274] - DC[0]*eh[ 182];
    eh[ 421] = eh[ 275] - DC[1]*eh[ 182];
    eh[ 422] = eh[ 272] - DC[2]*eh[ 182];
    eh[ 423] = eh[ 276] - DC[0]*eh[ 183];
    eh[ 424] = eh[ 279] - DC[1]*eh[ 183];
    eh[ 425] = eh[ 280] - DC[2]*eh[ 183];
    eh[ 426] = eh[ 279] - DC[0]*eh[ 184];
    eh[ 427] = eh[ 277] - DC[1]*eh[ 184];
    eh[ 428] = eh[ 281] - DC[2]*eh[ 184];
    eh[ 429] = eh[ 280] - DC[0]*eh[ 185];
    eh[ 430] = eh[ 281] - DC[1]*eh[ 185];
    eh[ 431] = eh[ 278] - DC[2]*eh[ 185];
    eh[ 432] = eh[ 282] - DC[0]*eh[ 186];
    eh[ 433] = eh[ 285] - DC[1]*eh[ 186];
    eh[ 434] = eh[ 286] - DC[2]*eh[ 186];
    eh[ 435] = eh[ 285] - DC[0]*eh[ 187];
    eh[ 436] = eh[ 283] - DC[1]*eh[ 187];
    eh[ 437] = eh[ 287] - DC[2]*eh[ 187];
    eh[ 438] = eh[ 286] - DC[0]*eh[ 188];
    eh[ 439] = eh[ 287] - DC[1]*eh[ 188];
    eh[ 440] = eh[ 284] - DC[2]*eh[ 188];
    eh[ 441] = eh[ 288] - DC[0]*eh[ 189];
    eh[ 442] = eh[ 291] - DC[1]*eh[ 189];
    eh[ 443] = eh[ 292] - DC[2]*eh[ 189];
    eh[ 444] = eh[ 291] - DC[0]*eh[ 190];
    eh[ 445] = eh[ 289] - DC[1]*eh[ 190];
    eh[ 446] = eh[ 293] - DC[2]*eh[ 190];
    eh[ 447] = eh[ 292] - DC[0]*eh[ 191];
    eh[ 448] = eh[ 293] - DC[1]*eh[ 191];
    eh[ 449] = eh[ 290] - DC[2]*eh[ 191];
    eh[ 450] = eh[ 294] - DC[0]*eh[ 192];
    eh[ 451] = eh[ 297] - DC[1]*eh[ 192];
    eh[ 452] = eh[ 298] - DC[2]*eh[ 192];
    eh[ 453] = eh[ 297] - DC[0]*eh[ 193];
    eh[ 454] = eh[ 295] - DC[1]*eh[ 193];
    eh[ 455] = eh[ 299] - DC[2]*eh[ 193];
    eh[ 456] = eh[ 298] - DC[0]*eh[ 194];
    eh[ 457] = eh[ 299] - DC[1]*eh[ 194];
    eh[ 458] = eh[ 296] - DC[2]*eh[ 194];
    eh[ 459] = eh[ 300] - DC[0]*eh[ 195];
    eh[ 460] = eh[ 303] - DC[1]*eh[ 195];
    eh[ 461] = eh[ 304] - DC[2]*eh[ 195];
    eh[ 462] = eh[ 303] - DC[0]*eh[ 196];
    eh[ 463] = eh[ 301] - DC[1]*eh[ 196];
    eh[ 464] = eh[ 305] - DC[2]*eh[ 196];
    eh[ 465] = eh[ 304] - DC[0]*eh[ 197];
    eh[ 466] = eh[ 305] - DC[1]*eh[ 197];
    eh[ 467] = eh[ 302] - DC[2]*eh[ 197];
}

static void ofmo_xyzint_dppp(
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

static void ofmo_form_dppp(
        const double *xint, const double *yint, const double *zint,
        double *eh ) {
    // (DS|PS)
    eh[   0] += xint[ 21]*yint[  0]*zint[  0];
    eh[   0] += xint[ 22]*yint[  1]*zint[  1];
    eh[   0] += xint[ 23]*yint[  2]*zint[  2];
    eh[   1] += xint[ 18]*yint[  3]*zint[  0];
    eh[   1] += xint[ 19]*yint[  4]*zint[  1];
    eh[   1] += xint[ 20]*yint[  5]*zint[  2];
    eh[   2] += xint[ 18]*yint[  0]*zint[  3];
    eh[   2] += xint[ 19]*yint[  1]*zint[  4];
    eh[   2] += xint[ 20]*yint[  2]*zint[  5];
    eh[   3] += xint[  3]*yint[ 18]*zint[  0];
    eh[   3] += xint[  4]*yint[ 19]*zint[  1];
    eh[   3] += xint[  5]*yint[ 20]*zint[  2];
    eh[   4] += xint[  0]*yint[ 21]*zint[  0];
    eh[   4] += xint[  1]*yint[ 22]*zint[  1];
    eh[   4] += xint[  2]*yint[ 23]*zint[  2];
    eh[   5] += xint[  0]*yint[ 18]*zint[  3];
    eh[   5] += xint[  1]*yint[ 19]*zint[  4];
    eh[   5] += xint[  2]*yint[ 20]*zint[  5];
    eh[   6] += xint[  3]*yint[  0]*zint[ 18];
    eh[   6] += xint[  4]*yint[  1]*zint[ 19];
    eh[   6] += xint[  5]*yint[  2]*zint[ 20];
    eh[   7] += xint[  0]*yint[  3]*zint[ 18];
    eh[   7] += xint[  1]*yint[  4]*zint[ 19];
    eh[   7] += xint[  2]*yint[  5]*zint[ 20];
    eh[   8] += xint[  0]*yint[  0]*zint[ 21];
    eh[   8] += xint[  1]*yint[  1]*zint[ 22];
    eh[   8] += xint[  2]*yint[  2]*zint[ 23];
    eh[   9] += xint[ 12]*yint[  9]*zint[  0];
    eh[   9] += xint[ 13]*yint[ 10]*zint[  1];
    eh[   9] += xint[ 14]*yint[ 11]*zint[  2];
    eh[  10] += xint[  9]*yint[ 12]*zint[  0];
    eh[  10] += xint[ 10]*yint[ 13]*zint[  1];
    eh[  10] += xint[ 11]*yint[ 14]*zint[  2];
    eh[  11] += xint[  9]*yint[  9]*zint[  3];
    eh[  11] += xint[ 10]*yint[ 10]*zint[  4];
    eh[  11] += xint[ 11]*yint[ 11]*zint[  5];
    eh[  12] += xint[ 12]*yint[  0]*zint[  9];
    eh[  12] += xint[ 13]*yint[  1]*zint[ 10];
    eh[  12] += xint[ 14]*yint[  2]*zint[ 11];
    eh[  13] += xint[  9]*yint[  3]*zint[  9];
    eh[  13] += xint[ 10]*yint[  4]*zint[ 10];
    eh[  13] += xint[ 11]*yint[  5]*zint[ 11];
    eh[  14] += xint[  9]*yint[  0]*zint[ 12];
    eh[  14] += xint[ 10]*yint[  1]*zint[ 13];
    eh[  14] += xint[ 11]*yint[  2]*zint[ 14];
    eh[  15] += xint[  3]*yint[  9]*zint[  9];
    eh[  15] += xint[  4]*yint[ 10]*zint[ 10];
    eh[  15] += xint[  5]*yint[ 11]*zint[ 11];
    eh[  16] += xint[  0]*yint[ 12]*zint[  9];
    eh[  16] += xint[  1]*yint[ 13]*zint[ 10];
    eh[  16] += xint[  2]*yint[ 14]*zint[ 11];
    eh[  17] += xint[  0]*yint[  9]*zint[ 12];
    eh[  17] += xint[  1]*yint[ 10]*zint[ 13];
    eh[  17] += xint[  2]*yint[ 11]*zint[ 14];
    // (DS|DS)
    eh[  18] += xint[ 24]*yint[  0]*zint[  0];
    eh[  18] += xint[ 25]*yint[  1]*zint[  1];
    eh[  18] += xint[ 26]*yint[  2]*zint[  2];
    eh[  19] += xint[ 18]*yint[  6]*zint[  0];
    eh[  19] += xint[ 19]*yint[  7]*zint[  1];
    eh[  19] += xint[ 20]*yint[  8]*zint[  2];
    eh[  20] += xint[ 18]*yint[  0]*zint[  6];
    eh[  20] += xint[ 19]*yint[  1]*zint[  7];
    eh[  20] += xint[ 20]*yint[  2]*zint[  8];
    eh[  21] += xint[ 21]*yint[  3]*zint[  0];
    eh[  21] += xint[ 22]*yint[  4]*zint[  1];
    eh[  21] += xint[ 23]*yint[  5]*zint[  2];
    eh[  22] += xint[ 21]*yint[  0]*zint[  3];
    eh[  22] += xint[ 22]*yint[  1]*zint[  4];
    eh[  22] += xint[ 23]*yint[  2]*zint[  5];
    eh[  23] += xint[ 18]*yint[  3]*zint[  3];
    eh[  23] += xint[ 19]*yint[  4]*zint[  4];
    eh[  23] += xint[ 20]*yint[  5]*zint[  5];
    eh[  24] += xint[  6]*yint[ 18]*zint[  0];
    eh[  24] += xint[  7]*yint[ 19]*zint[  1];
    eh[  24] += xint[  8]*yint[ 20]*zint[  2];
    eh[  25] += xint[  0]*yint[ 24]*zint[  0];
    eh[  25] += xint[  1]*yint[ 25]*zint[  1];
    eh[  25] += xint[  2]*yint[ 26]*zint[  2];
    eh[  26] += xint[  0]*yint[ 18]*zint[  6];
    eh[  26] += xint[  1]*yint[ 19]*zint[  7];
    eh[  26] += xint[  2]*yint[ 20]*zint[  8];
    eh[  27] += xint[  3]*yint[ 21]*zint[  0];
    eh[  27] += xint[  4]*yint[ 22]*zint[  1];
    eh[  27] += xint[  5]*yint[ 23]*zint[  2];
    eh[  28] += xint[  3]*yint[ 18]*zint[  3];
    eh[  28] += xint[  4]*yint[ 19]*zint[  4];
    eh[  28] += xint[  5]*yint[ 20]*zint[  5];
    eh[  29] += xint[  0]*yint[ 21]*zint[  3];
    eh[  29] += xint[  1]*yint[ 22]*zint[  4];
    eh[  29] += xint[  2]*yint[ 23]*zint[  5];
    eh[  30] += xint[  6]*yint[  0]*zint[ 18];
    eh[  30] += xint[  7]*yint[  1]*zint[ 19];
    eh[  30] += xint[  8]*yint[  2]*zint[ 20];
    eh[  31] += xint[  0]*yint[  6]*zint[ 18];
    eh[  31] += xint[  1]*yint[  7]*zint[ 19];
    eh[  31] += xint[  2]*yint[  8]*zint[ 20];
    eh[  32] += xint[  0]*yint[  0]*zint[ 24];
    eh[  32] += xint[  1]*yint[  1]*zint[ 25];
    eh[  32] += xint[  2]*yint[  2]*zint[ 26];
    eh[  33] += xint[  3]*yint[  3]*zint[ 18];
    eh[  33] += xint[  4]*yint[  4]*zint[ 19];
    eh[  33] += xint[  5]*yint[  5]*zint[ 20];
    eh[  34] += xint[  3]*yint[  0]*zint[ 21];
    eh[  34] += xint[  4]*yint[  1]*zint[ 22];
    eh[  34] += xint[  5]*yint[  2]*zint[ 23];
    eh[  35] += xint[  0]*yint[  3]*zint[ 21];
    eh[  35] += xint[  1]*yint[  4]*zint[ 22];
    eh[  35] += xint[  2]*yint[  5]*zint[ 23];
    eh[  36] += xint[ 15]*yint[  9]*zint[  0];
    eh[  36] += xint[ 16]*yint[ 10]*zint[  1];
    eh[  36] += xint[ 17]*yint[ 11]*zint[  2];
    eh[  37] += xint[  9]*yint[ 15]*zint[  0];
    eh[  37] += xint[ 10]*yint[ 16]*zint[  1];
    eh[  37] += xint[ 11]*yint[ 17]*zint[  2];
    eh[  38] += xint[  9]*yint[  9]*zint[  6];
    eh[  38] += xint[ 10]*yint[ 10]*zint[  7];
    eh[  38] += xint[ 11]*yint[ 11]*zint[  8];
    eh[  39] += xint[ 12]*yint[ 12]*zint[  0];
    eh[  39] += xint[ 13]*yint[ 13]*zint[  1];
    eh[  39] += xint[ 14]*yint[ 14]*zint[  2];
    eh[  40] += xint[ 12]*yint[  9]*zint[  3];
    eh[  40] += xint[ 13]*yint[ 10]*zint[  4];
    eh[  40] += xint[ 14]*yint[ 11]*zint[  5];
    eh[  41] += xint[  9]*yint[ 12]*zint[  3];
    eh[  41] += xint[ 10]*yint[ 13]*zint[  4];
    eh[  41] += xint[ 11]*yint[ 14]*zint[  5];
    eh[  42] += xint[ 15]*yint[  0]*zint[  9];
    eh[  42] += xint[ 16]*yint[  1]*zint[ 10];
    eh[  42] += xint[ 17]*yint[  2]*zint[ 11];
    eh[  43] += xint[  9]*yint[  6]*zint[  9];
    eh[  43] += xint[ 10]*yint[  7]*zint[ 10];
    eh[  43] += xint[ 11]*yint[  8]*zint[ 11];
    eh[  44] += xint[  9]*yint[  0]*zint[ 15];
    eh[  44] += xint[ 10]*yint[  1]*zint[ 16];
    eh[  44] += xint[ 11]*yint[  2]*zint[ 17];
    eh[  45] += xint[ 12]*yint[  3]*zint[  9];
    eh[  45] += xint[ 13]*yint[  4]*zint[ 10];
    eh[  45] += xint[ 14]*yint[  5]*zint[ 11];
    eh[  46] += xint[ 12]*yint[  0]*zint[ 12];
    eh[  46] += xint[ 13]*yint[  1]*zint[ 13];
    eh[  46] += xint[ 14]*yint[  2]*zint[ 14];
    eh[  47] += xint[  9]*yint[  3]*zint[ 12];
    eh[  47] += xint[ 10]*yint[  4]*zint[ 13];
    eh[  47] += xint[ 11]*yint[  5]*zint[ 14];
    eh[  48] += xint[  6]*yint[  9]*zint[  9];
    eh[  48] += xint[  7]*yint[ 10]*zint[ 10];
    eh[  48] += xint[  8]*yint[ 11]*zint[ 11];
    eh[  49] += xint[  0]*yint[ 15]*zint[  9];
    eh[  49] += xint[  1]*yint[ 16]*zint[ 10];
    eh[  49] += xint[  2]*yint[ 17]*zint[ 11];
    eh[  50] += xint[  0]*yint[  9]*zint[ 15];
    eh[  50] += xint[  1]*yint[ 10]*zint[ 16];
    eh[  50] += xint[  2]*yint[ 11]*zint[ 17];
    eh[  51] += xint[  3]*yint[ 12]*zint[  9];
    eh[  51] += xint[  4]*yint[ 13]*zint[ 10];
    eh[  51] += xint[  5]*yint[ 14]*zint[ 11];
    eh[  52] += xint[  3]*yint[  9]*zint[ 12];
    eh[  52] += xint[  4]*yint[ 10]*zint[ 13];
    eh[  52] += xint[  5]*yint[ 11]*zint[ 14];
    eh[  53] += xint[  0]*yint[ 12]*zint[ 12];
    eh[  53] += xint[  1]*yint[ 13]*zint[ 13];
    eh[  53] += xint[  2]*yint[ 14]*zint[ 14];
    // (FS|PS)
    eh[  54] += xint[ 30]*yint[  0]*zint[  0];
    eh[  54] += xint[ 31]*yint[  1]*zint[  1];
    eh[  54] += xint[ 32]*yint[  2]*zint[  2];
    eh[  55] += xint[ 27]*yint[  3]*zint[  0];
    eh[  55] += xint[ 28]*yint[  4]*zint[  1];
    eh[  55] += xint[ 29]*yint[  5]*zint[  2];
    eh[  56] += xint[ 27]*yint[  0]*zint[  3];
    eh[  56] += xint[ 28]*yint[  1]*zint[  4];
    eh[  56] += xint[ 29]*yint[  2]*zint[  5];
    eh[  57] += xint[  3]*yint[ 27]*zint[  0];
    eh[  57] += xint[  4]*yint[ 28]*zint[  1];
    eh[  57] += xint[  5]*yint[ 29]*zint[  2];
    eh[  58] += xint[  0]*yint[ 30]*zint[  0];
    eh[  58] += xint[  1]*yint[ 31]*zint[  1];
    eh[  58] += xint[  2]*yint[ 32]*zint[  2];
    eh[  59] += xint[  0]*yint[ 27]*zint[  3];
    eh[  59] += xint[  1]*yint[ 28]*zint[  4];
    eh[  59] += xint[  2]*yint[ 29]*zint[  5];
    eh[  60] += xint[  3]*yint[  0]*zint[ 27];
    eh[  60] += xint[  4]*yint[  1]*zint[ 28];
    eh[  60] += xint[  5]*yint[  2]*zint[ 29];
    eh[  61] += xint[  0]*yint[  3]*zint[ 27];
    eh[  61] += xint[  1]*yint[  4]*zint[ 28];
    eh[  61] += xint[  2]*yint[  5]*zint[ 29];
    eh[  62] += xint[  0]*yint[  0]*zint[ 30];
    eh[  62] += xint[  1]*yint[  1]*zint[ 31];
    eh[  62] += xint[  2]*yint[  2]*zint[ 32];
    eh[  63] += xint[ 21]*yint[  9]*zint[  0];
    eh[  63] += xint[ 22]*yint[ 10]*zint[  1];
    eh[  63] += xint[ 23]*yint[ 11]*zint[  2];
    eh[  64] += xint[ 18]*yint[ 12]*zint[  0];
    eh[  64] += xint[ 19]*yint[ 13]*zint[  1];
    eh[  64] += xint[ 20]*yint[ 14]*zint[  2];
    eh[  65] += xint[ 18]*yint[  9]*zint[  3];
    eh[  65] += xint[ 19]*yint[ 10]*zint[  4];
    eh[  65] += xint[ 20]*yint[ 11]*zint[  5];
    eh[  66] += xint[ 21]*yint[  0]*zint[  9];
    eh[  66] += xint[ 22]*yint[  1]*zint[ 10];
    eh[  66] += xint[ 23]*yint[  2]*zint[ 11];
    eh[  67] += xint[ 18]*yint[  3]*zint[  9];
    eh[  67] += xint[ 19]*yint[  4]*zint[ 10];
    eh[  67] += xint[ 20]*yint[  5]*zint[ 11];
    eh[  68] += xint[ 18]*yint[  0]*zint[ 12];
    eh[  68] += xint[ 19]*yint[  1]*zint[ 13];
    eh[  68] += xint[ 20]*yint[  2]*zint[ 14];
    eh[  69] += xint[ 12]*yint[ 18]*zint[  0];
    eh[  69] += xint[ 13]*yint[ 19]*zint[  1];
    eh[  69] += xint[ 14]*yint[ 20]*zint[  2];
    eh[  70] += xint[  9]*yint[ 21]*zint[  0];
    eh[  70] += xint[ 10]*yint[ 22]*zint[  1];
    eh[  70] += xint[ 11]*yint[ 23]*zint[  2];
    eh[  71] += xint[  9]*yint[ 18]*zint[  3];
    eh[  71] += xint[ 10]*yint[ 19]*zint[  4];
    eh[  71] += xint[ 11]*yint[ 20]*zint[  5];
    eh[  72] += xint[ 12]*yint[  0]*zint[ 18];
    eh[  72] += xint[ 13]*yint[  1]*zint[ 19];
    eh[  72] += xint[ 14]*yint[  2]*zint[ 20];
    eh[  73] += xint[  9]*yint[  3]*zint[ 18];
    eh[  73] += xint[ 10]*yint[  4]*zint[ 19];
    eh[  73] += xint[ 11]*yint[  5]*zint[ 20];
    eh[  74] += xint[  9]*yint[  0]*zint[ 21];
    eh[  74] += xint[ 10]*yint[  1]*zint[ 22];
    eh[  74] += xint[ 11]*yint[  2]*zint[ 23];
    eh[  75] += xint[ 12]*yint[  9]*zint[  9];
    eh[  75] += xint[ 13]*yint[ 10]*zint[ 10];
    eh[  75] += xint[ 14]*yint[ 11]*zint[ 11];
    eh[  76] += xint[  9]*yint[ 12]*zint[  9];
    eh[  76] += xint[ 10]*yint[ 13]*zint[ 10];
    eh[  76] += xint[ 11]*yint[ 14]*zint[ 11];
    eh[  77] += xint[  9]*yint[  9]*zint[ 12];
    eh[  77] += xint[ 10]*yint[ 10]*zint[ 13];
    eh[  77] += xint[ 11]*yint[ 11]*zint[ 14];
    eh[  78] += xint[  3]*yint[ 18]*zint[  9];
    eh[  78] += xint[  4]*yint[ 19]*zint[ 10];
    eh[  78] += xint[  5]*yint[ 20]*zint[ 11];
    eh[  79] += xint[  0]*yint[ 21]*zint[  9];
    eh[  79] += xint[  1]*yint[ 22]*zint[ 10];
    eh[  79] += xint[  2]*yint[ 23]*zint[ 11];
    eh[  80] += xint[  0]*yint[ 18]*zint[ 12];
    eh[  80] += xint[  1]*yint[ 19]*zint[ 13];
    eh[  80] += xint[  2]*yint[ 20]*zint[ 14];
    eh[  81] += xint[  3]*yint[  9]*zint[ 18];
    eh[  81] += xint[  4]*yint[ 10]*zint[ 19];
    eh[  81] += xint[  5]*yint[ 11]*zint[ 20];
    eh[  82] += xint[  0]*yint[ 12]*zint[ 18];
    eh[  82] += xint[  1]*yint[ 13]*zint[ 19];
    eh[  82] += xint[  2]*yint[ 14]*zint[ 20];
    eh[  83] += xint[  0]*yint[  9]*zint[ 21];
    eh[  83] += xint[  1]*yint[ 10]*zint[ 22];
    eh[  83] += xint[  2]*yint[ 11]*zint[ 23];
    // (FS|DS)
    eh[  84] += xint[ 33]*yint[  0]*zint[  0];
    eh[  84] += xint[ 34]*yint[  1]*zint[  1];
    eh[  84] += xint[ 35]*yint[  2]*zint[  2];
    eh[  85] += xint[ 27]*yint[  6]*zint[  0];
    eh[  85] += xint[ 28]*yint[  7]*zint[  1];
    eh[  85] += xint[ 29]*yint[  8]*zint[  2];
    eh[  86] += xint[ 27]*yint[  0]*zint[  6];
    eh[  86] += xint[ 28]*yint[  1]*zint[  7];
    eh[  86] += xint[ 29]*yint[  2]*zint[  8];
    eh[  87] += xint[ 30]*yint[  3]*zint[  0];
    eh[  87] += xint[ 31]*yint[  4]*zint[  1];
    eh[  87] += xint[ 32]*yint[  5]*zint[  2];
    eh[  88] += xint[ 30]*yint[  0]*zint[  3];
    eh[  88] += xint[ 31]*yint[  1]*zint[  4];
    eh[  88] += xint[ 32]*yint[  2]*zint[  5];
    eh[  89] += xint[ 27]*yint[  3]*zint[  3];
    eh[  89] += xint[ 28]*yint[  4]*zint[  4];
    eh[  89] += xint[ 29]*yint[  5]*zint[  5];
    eh[  90] += xint[  6]*yint[ 27]*zint[  0];
    eh[  90] += xint[  7]*yint[ 28]*zint[  1];
    eh[  90] += xint[  8]*yint[ 29]*zint[  2];
    eh[  91] += xint[  0]*yint[ 33]*zint[  0];
    eh[  91] += xint[  1]*yint[ 34]*zint[  1];
    eh[  91] += xint[  2]*yint[ 35]*zint[  2];
    eh[  92] += xint[  0]*yint[ 27]*zint[  6];
    eh[  92] += xint[  1]*yint[ 28]*zint[  7];
    eh[  92] += xint[  2]*yint[ 29]*zint[  8];
    eh[  93] += xint[  3]*yint[ 30]*zint[  0];
    eh[  93] += xint[  4]*yint[ 31]*zint[  1];
    eh[  93] += xint[  5]*yint[ 32]*zint[  2];
    eh[  94] += xint[  3]*yint[ 27]*zint[  3];
    eh[  94] += xint[  4]*yint[ 28]*zint[  4];
    eh[  94] += xint[  5]*yint[ 29]*zint[  5];
    eh[  95] += xint[  0]*yint[ 30]*zint[  3];
    eh[  95] += xint[  1]*yint[ 31]*zint[  4];
    eh[  95] += xint[  2]*yint[ 32]*zint[  5];
    eh[  96] += xint[  6]*yint[  0]*zint[ 27];
    eh[  96] += xint[  7]*yint[  1]*zint[ 28];
    eh[  96] += xint[  8]*yint[  2]*zint[ 29];
    eh[  97] += xint[  0]*yint[  6]*zint[ 27];
    eh[  97] += xint[  1]*yint[  7]*zint[ 28];
    eh[  97] += xint[  2]*yint[  8]*zint[ 29];
    eh[  98] += xint[  0]*yint[  0]*zint[ 33];
    eh[  98] += xint[  1]*yint[  1]*zint[ 34];
    eh[  98] += xint[  2]*yint[  2]*zint[ 35];
    eh[  99] += xint[  3]*yint[  3]*zint[ 27];
    eh[  99] += xint[  4]*yint[  4]*zint[ 28];
    eh[  99] += xint[  5]*yint[  5]*zint[ 29];
    eh[ 100] += xint[  3]*yint[  0]*zint[ 30];
    eh[ 100] += xint[  4]*yint[  1]*zint[ 31];
    eh[ 100] += xint[  5]*yint[  2]*zint[ 32];
    eh[ 101] += xint[  0]*yint[  3]*zint[ 30];
    eh[ 101] += xint[  1]*yint[  4]*zint[ 31];
    eh[ 101] += xint[  2]*yint[  5]*zint[ 32];
    eh[ 102] += xint[ 24]*yint[  9]*zint[  0];
    eh[ 102] += xint[ 25]*yint[ 10]*zint[  1];
    eh[ 102] += xint[ 26]*yint[ 11]*zint[  2];
    eh[ 103] += xint[ 18]*yint[ 15]*zint[  0];
    eh[ 103] += xint[ 19]*yint[ 16]*zint[  1];
    eh[ 103] += xint[ 20]*yint[ 17]*zint[  2];
    eh[ 104] += xint[ 18]*yint[  9]*zint[  6];
    eh[ 104] += xint[ 19]*yint[ 10]*zint[  7];
    eh[ 104] += xint[ 20]*yint[ 11]*zint[  8];
    eh[ 105] += xint[ 21]*yint[ 12]*zint[  0];
    eh[ 105] += xint[ 22]*yint[ 13]*zint[  1];
    eh[ 105] += xint[ 23]*yint[ 14]*zint[  2];
    eh[ 106] += xint[ 21]*yint[  9]*zint[  3];
    eh[ 106] += xint[ 22]*yint[ 10]*zint[  4];
    eh[ 106] += xint[ 23]*yint[ 11]*zint[  5];
    eh[ 107] += xint[ 18]*yint[ 12]*zint[  3];
    eh[ 107] += xint[ 19]*yint[ 13]*zint[  4];
    eh[ 107] += xint[ 20]*yint[ 14]*zint[  5];
    eh[ 108] += xint[ 24]*yint[  0]*zint[  9];
    eh[ 108] += xint[ 25]*yint[  1]*zint[ 10];
    eh[ 108] += xint[ 26]*yint[  2]*zint[ 11];
    eh[ 109] += xint[ 18]*yint[  6]*zint[  9];
    eh[ 109] += xint[ 19]*yint[  7]*zint[ 10];
    eh[ 109] += xint[ 20]*yint[  8]*zint[ 11];
    eh[ 110] += xint[ 18]*yint[  0]*zint[ 15];
    eh[ 110] += xint[ 19]*yint[  1]*zint[ 16];
    eh[ 110] += xint[ 20]*yint[  2]*zint[ 17];
    eh[ 111] += xint[ 21]*yint[  3]*zint[  9];
    eh[ 111] += xint[ 22]*yint[  4]*zint[ 10];
    eh[ 111] += xint[ 23]*yint[  5]*zint[ 11];
    eh[ 112] += xint[ 21]*yint[  0]*zint[ 12];
    eh[ 112] += xint[ 22]*yint[  1]*zint[ 13];
    eh[ 112] += xint[ 23]*yint[  2]*zint[ 14];
    eh[ 113] += xint[ 18]*yint[  3]*zint[ 12];
    eh[ 113] += xint[ 19]*yint[  4]*zint[ 13];
    eh[ 113] += xint[ 20]*yint[  5]*zint[ 14];
    eh[ 114] += xint[ 15]*yint[ 18]*zint[  0];
    eh[ 114] += xint[ 16]*yint[ 19]*zint[  1];
    eh[ 114] += xint[ 17]*yint[ 20]*zint[  2];
    eh[ 115] += xint[  9]*yint[ 24]*zint[  0];
    eh[ 115] += xint[ 10]*yint[ 25]*zint[  1];
    eh[ 115] += xint[ 11]*yint[ 26]*zint[  2];
    eh[ 116] += xint[  9]*yint[ 18]*zint[  6];
    eh[ 116] += xint[ 10]*yint[ 19]*zint[  7];
    eh[ 116] += xint[ 11]*yint[ 20]*zint[  8];
    eh[ 117] += xint[ 12]*yint[ 21]*zint[  0];
    eh[ 117] += xint[ 13]*yint[ 22]*zint[  1];
    eh[ 117] += xint[ 14]*yint[ 23]*zint[  2];
    eh[ 118] += xint[ 12]*yint[ 18]*zint[  3];
    eh[ 118] += xint[ 13]*yint[ 19]*zint[  4];
    eh[ 118] += xint[ 14]*yint[ 20]*zint[  5];
    eh[ 119] += xint[  9]*yint[ 21]*zint[  3];
    eh[ 119] += xint[ 10]*yint[ 22]*zint[  4];
    eh[ 119] += xint[ 11]*yint[ 23]*zint[  5];
    eh[ 120] += xint[ 15]*yint[  0]*zint[ 18];
    eh[ 120] += xint[ 16]*yint[  1]*zint[ 19];
    eh[ 120] += xint[ 17]*yint[  2]*zint[ 20];
    eh[ 121] += xint[  9]*yint[  6]*zint[ 18];
    eh[ 121] += xint[ 10]*yint[  7]*zint[ 19];
    eh[ 121] += xint[ 11]*yint[  8]*zint[ 20];
    eh[ 122] += xint[  9]*yint[  0]*zint[ 24];
    eh[ 122] += xint[ 10]*yint[  1]*zint[ 25];
    eh[ 122] += xint[ 11]*yint[  2]*zint[ 26];
    eh[ 123] += xint[ 12]*yint[  3]*zint[ 18];
    eh[ 123] += xint[ 13]*yint[  4]*zint[ 19];
    eh[ 123] += xint[ 14]*yint[  5]*zint[ 20];
    eh[ 124] += xint[ 12]*yint[  0]*zint[ 21];
    eh[ 124] += xint[ 13]*yint[  1]*zint[ 22];
    eh[ 124] += xint[ 14]*yint[  2]*zint[ 23];
    eh[ 125] += xint[  9]*yint[  3]*zint[ 21];
    eh[ 125] += xint[ 10]*yint[  4]*zint[ 22];
    eh[ 125] += xint[ 11]*yint[  5]*zint[ 23];
    eh[ 126] += xint[ 15]*yint[  9]*zint[  9];
    eh[ 126] += xint[ 16]*yint[ 10]*zint[ 10];
    eh[ 126] += xint[ 17]*yint[ 11]*zint[ 11];
    eh[ 127] += xint[  9]*yint[ 15]*zint[  9];
    eh[ 127] += xint[ 10]*yint[ 16]*zint[ 10];
    eh[ 127] += xint[ 11]*yint[ 17]*zint[ 11];
    eh[ 128] += xint[  9]*yint[  9]*zint[ 15];
    eh[ 128] += xint[ 10]*yint[ 10]*zint[ 16];
    eh[ 128] += xint[ 11]*yint[ 11]*zint[ 17];
    eh[ 129] += xint[ 12]*yint[ 12]*zint[  9];
    eh[ 129] += xint[ 13]*yint[ 13]*zint[ 10];
    eh[ 129] += xint[ 14]*yint[ 14]*zint[ 11];
    eh[ 130] += xint[ 12]*yint[  9]*zint[ 12];
    eh[ 130] += xint[ 13]*yint[ 10]*zint[ 13];
    eh[ 130] += xint[ 14]*yint[ 11]*zint[ 14];
    eh[ 131] += xint[  9]*yint[ 12]*zint[ 12];
    eh[ 131] += xint[ 10]*yint[ 13]*zint[ 13];
    eh[ 131] += xint[ 11]*yint[ 14]*zint[ 14];
    eh[ 132] += xint[  6]*yint[ 18]*zint[  9];
    eh[ 132] += xint[  7]*yint[ 19]*zint[ 10];
    eh[ 132] += xint[  8]*yint[ 20]*zint[ 11];
    eh[ 133] += xint[  0]*yint[ 24]*zint[  9];
    eh[ 133] += xint[  1]*yint[ 25]*zint[ 10];
    eh[ 133] += xint[  2]*yint[ 26]*zint[ 11];
    eh[ 134] += xint[  0]*yint[ 18]*zint[ 15];
    eh[ 134] += xint[  1]*yint[ 19]*zint[ 16];
    eh[ 134] += xint[  2]*yint[ 20]*zint[ 17];
    eh[ 135] += xint[  3]*yint[ 21]*zint[  9];
    eh[ 135] += xint[  4]*yint[ 22]*zint[ 10];
    eh[ 135] += xint[  5]*yint[ 23]*zint[ 11];
    eh[ 136] += xint[  3]*yint[ 18]*zint[ 12];
    eh[ 136] += xint[  4]*yint[ 19]*zint[ 13];
    eh[ 136] += xint[  5]*yint[ 20]*zint[ 14];
    eh[ 137] += xint[  0]*yint[ 21]*zint[ 12];
    eh[ 137] += xint[  1]*yint[ 22]*zint[ 13];
    eh[ 137] += xint[  2]*yint[ 23]*zint[ 14];
    eh[ 138] += xint[  6]*yint[  9]*zint[ 18];
    eh[ 138] += xint[  7]*yint[ 10]*zint[ 19];
    eh[ 138] += xint[  8]*yint[ 11]*zint[ 20];
    eh[ 139] += xint[  0]*yint[ 15]*zint[ 18];
    eh[ 139] += xint[  1]*yint[ 16]*zint[ 19];
    eh[ 139] += xint[  2]*yint[ 17]*zint[ 20];
    eh[ 140] += xint[  0]*yint[  9]*zint[ 24];
    eh[ 140] += xint[  1]*yint[ 10]*zint[ 25];
    eh[ 140] += xint[  2]*yint[ 11]*zint[ 26];
    eh[ 141] += xint[  3]*yint[ 12]*zint[ 18];
    eh[ 141] += xint[  4]*yint[ 13]*zint[ 19];
    eh[ 141] += xint[  5]*yint[ 14]*zint[ 20];
    eh[ 142] += xint[  3]*yint[  9]*zint[ 21];
    eh[ 142] += xint[  4]*yint[ 10]*zint[ 22];
    eh[ 142] += xint[  5]*yint[ 11]*zint[ 23];
    eh[ 143] += xint[  0]*yint[ 12]*zint[ 21];
    eh[ 143] += xint[  1]*yint[ 13]*zint[ 22];
    eh[ 143] += xint[  2]*yint[ 14]*zint[ 23];
}

void ofmo_twoint_core_rys_dppp( const int mythread,
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
    ofmo_hrr_clear_dppp( eh );
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
            ofmo_xyzint_dppp( F00, B00, B10, B01, C00, CP00,
                     xint, yint, zint );
            ofmo_form_dppp( xint, yint, zint, eh );
        }
    }
    ofmo_hrr_calc_dppp( eh, BA, DC );
    ofmo_hrr_coef_dppp( eh, DINT );
}

int ofmo_twoint_rys_dppp(
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
    max_nzeri = ebuf_max_nzeri - 6*3*3*3;
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
            ofmo_twoint_core_rys_dppp( mythread,
                    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
                    &psp_xiza[ijps0], BA,
                    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
                    &psp_xiza[klps0], DC,   AC,      DINTEG );
            ipat=((Lab != Lcd) || (ics==kcs && jcs>lcs) ? true : false);
            for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
                I2 = (iao*iao+iao)>>1;
                for ( j=0, jao=jao0; j<3; j++, jao++ ) {
                    if ( jao>iao ) { ix+=3*3; continue; }
                    IJ = I2 + jao;
                    coe0 = ( iao==jao ? HALF : ONE );
                    for ( k=0, kao=kao0; k<3; k++, kao++ ) {
                        K2 = (kao*kao+kao)>>1;
                        for ( l=0, lao=lao0; l<3; l++, lao++, ix++ ) {
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

int ofmo_twoint_direct_rys_dppp(
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
    max_nzeri -= 6*3*3*3;
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
            ofmo_twoint_core_rys_dppp( mythread,
                    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
                    &psp_xiza[ijps0], BA,
                    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
                    &psp_xiza[klps0], DC,   AC,      DINTEG );
            ipat=((Lab != Lcd) || (ics==kcs && jcs>lcs) ? true : false);
            for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
                I2 = (iao*iao+iao)>>1;
                for ( j=0, jao=jao0; j<3; j++, jao++ ) {
                    if ( jao>iao ) { ix+=3*3; continue; }
                    IJ = I2 + jao;
                    coe0 = ( iao==jao ? HALF : ONE );
                    for ( k=0, kao=kao0; k<3; k++, kao++ ) {
                        K2 = (kao*kao+kao)>>1;
                        for ( l=0, lao=lao0; l<3; l++, lao++, ix++ ) {
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
