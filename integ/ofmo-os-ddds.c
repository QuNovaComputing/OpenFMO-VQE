#include <stdio.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "ofmo-twoint.h"
#include "fmt-m.h"

#ifndef false
#define false 0
#endif
#ifndef true
#define true 1
#endif

#define HALF    0.5e0
#define ONE     1.e0
#define ZERO    0.e0


#define OFMO_EBUF_FULL          1
#define OFMO_EBUF_NOFULL        0

#define EPS_PS_PAIR     1.e-32
#define EPS_CS_PAIR2    1.e-30
#define MAXNPSPAIR 100

extern double* ofmo_os_getadd_vrr( const int mythread );
extern double* ofmo_os_getadd_hrr( const int mythread );
extern double* ofmo_os_getadd_eri( const int mythread );

extern void fmt( double F[],
        const int m, const double T, const double cssss );
extern int ofmo_integ_add_fock( const int nao, const size_t nstored_eri,
        const double eri_val[], const short int eri_ind4[],
        const double D[], double G[] );

extern double* ofmo_getadd_dfact();
static double *DFACT = NULL;

static void ofmo_hrr_clear_ddds( double *eh ) {
    int i;
    // (DS|DS)
    for ( i=0; i<(0+36); i++ ) eh[i] = 0.e0;
    // (FS|DS)
    for ( i=36; i<(36+60); i++ ) eh[i] = 0.e0;
    // (GS|DS)
    for ( i=96; i<(96+90); i++ ) eh[i] = 0.e0;
}

static void ofmo_hrr_calc_ddds(
        const double BA[3], const double DC[3], double *eh ) {
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
    // HRR for (XX|XX)-type integral (center CD)
}

static void ofmo_hrr_coef_ddds(
        const double *eh, double *DINT ) {
    int i, j, k, l, iao, jao, kao, lao, ix, iy;
    double coef_a, coef_ab, coef_abc;
    ix = 474;
    iy = 0;

    for ( i=0, iao=4; i<6; i++, iao++ ) {
        coef_a = DFACT[iao];
        for ( j=0, jao=4; j<6; j++, jao++ ) {
            coef_ab = coef_a * DFACT[jao];
            for ( k=0, kao=4; k<6; k++, kao++ ) {
                coef_abc = coef_ab * DFACT[kao];
                DINT[iy] = coef_abc * eh[ix];
                iy++;
                ix++;
            }
        }
    }
}

static void ofmo_vrr_calc_ddds(
        const double T, const double cssss,
        const double zeta2, const double eta2, const double ze2,
        const double rz, const double re,
        const double PA[3], const double WP[3],
        const double QC[3], const double WQ[3],
        double *ev ) {
    // (ss|ss) m=0,6
    //fmt( &ev[0], 6, T, cssss );
    OFMO_FMT( &ev[0], 6, T, cssss );
    // (ps|ss) m=0,5
    ev[ 7]=PA[0]*ev[0]+WP[0]*ev[1];
    ev[ 8]=PA[1]*ev[0]+WP[1]*ev[1];
    ev[ 9]=PA[2]*ev[0]+WP[2]*ev[1];
    ev[10]=PA[0]*ev[1]+WP[0]*ev[2];
    ev[11]=PA[1]*ev[1]+WP[1]*ev[2];
    ev[12]=PA[2]*ev[1]+WP[2]*ev[2];
    ev[13]=PA[0]*ev[2]+WP[0]*ev[3];
    ev[14]=PA[1]*ev[2]+WP[1]*ev[3];
    ev[15]=PA[2]*ev[2]+WP[2]*ev[3];
    ev[16]=PA[0]*ev[3]+WP[0]*ev[4];
    ev[17]=PA[1]*ev[3]+WP[1]*ev[4];
    ev[18]=PA[2]*ev[3]+WP[2]*ev[4];
    ev[19]=PA[0]*ev[4]+WP[0]*ev[5];
    ev[20]=PA[1]*ev[4]+WP[1]*ev[5];
    ev[21]=PA[2]*ev[4]+WP[2]*ev[5];
    ev[22]=PA[0]*ev[5]+WP[0]*ev[6];
    ev[23]=PA[1]*ev[5]+WP[1]*ev[6];
    ev[24]=PA[2]*ev[5]+WP[2]*ev[6];
    // (ds|ss) m=0,4
    ev[25]=PA[0]*ev[ 7]+WP[0]*ev[10]+zeta2*(ev[0]-rz*ev[1]);
    ev[26]=PA[1]*ev[ 8]+WP[1]*ev[11]+zeta2*(ev[0]-rz*ev[1]);
    ev[27]=PA[2]*ev[ 9]+WP[2]*ev[12]+zeta2*(ev[0]-rz*ev[1]);
    ev[28]=PA[0]*ev[ 8]+WP[0]*ev[11];
    ev[29]=PA[0]*ev[ 9]+WP[0]*ev[12];
    ev[30]=PA[1]*ev[ 9]+WP[1]*ev[12];
    ev[31]=PA[0]*ev[10]+WP[0]*ev[13]+zeta2*(ev[1]-rz*ev[2]);
    ev[32]=PA[1]*ev[11]+WP[1]*ev[14]+zeta2*(ev[1]-rz*ev[2]);
    ev[33]=PA[2]*ev[12]+WP[2]*ev[15]+zeta2*(ev[1]-rz*ev[2]);
    ev[34]=PA[0]*ev[11]+WP[0]*ev[14];
    ev[35]=PA[0]*ev[12]+WP[0]*ev[15];
    ev[36]=PA[1]*ev[12]+WP[1]*ev[15];
    ev[37]=PA[0]*ev[13]+WP[0]*ev[16]+zeta2*(ev[2]-rz*ev[3]);
    ev[38]=PA[1]*ev[14]+WP[1]*ev[17]+zeta2*(ev[2]-rz*ev[3]);
    ev[39]=PA[2]*ev[15]+WP[2]*ev[18]+zeta2*(ev[2]-rz*ev[3]);
    ev[40]=PA[0]*ev[14]+WP[0]*ev[17];
    ev[41]=PA[0]*ev[15]+WP[0]*ev[18];
    ev[42]=PA[1]*ev[15]+WP[1]*ev[18];
    ev[43]=PA[0]*ev[16]+WP[0]*ev[19]+zeta2*(ev[3]-rz*ev[4]);
    ev[44]=PA[1]*ev[17]+WP[1]*ev[20]+zeta2*(ev[3]-rz*ev[4]);
    ev[45]=PA[2]*ev[18]+WP[2]*ev[21]+zeta2*(ev[3]-rz*ev[4]);
    ev[46]=PA[0]*ev[17]+WP[0]*ev[20];
    ev[47]=PA[0]*ev[18]+WP[0]*ev[21];
    ev[48]=PA[1]*ev[18]+WP[1]*ev[21];
    ev[49]=PA[0]*ev[19]+WP[0]*ev[22]+zeta2*(ev[4]-rz*ev[5]);
    ev[50]=PA[1]*ev[20]+WP[1]*ev[23]+zeta2*(ev[4]-rz*ev[5]);
    ev[51]=PA[2]*ev[21]+WP[2]*ev[24]+zeta2*(ev[4]-rz*ev[5]);
    ev[52]=PA[0]*ev[20]+WP[0]*ev[23];
    ev[53]=PA[0]*ev[21]+WP[0]*ev[24];
    ev[54]=PA[1]*ev[21]+WP[1]*ev[24];
    // (fs|ss) m=0,3
    ev[55]=PA[0]*ev[25]+WP[0]*ev[31]+2.e0*zeta2*(ev[ 7]-rz*ev[10]);
    ev[56]=PA[1]*ev[26]+WP[1]*ev[32]+2.e0*zeta2*(ev[ 8]-rz*ev[11]);
    ev[57]=PA[2]*ev[27]+WP[2]*ev[33]+2.e0*zeta2*(ev[ 9]-rz*ev[12]);
    ev[58]=PA[1]*ev[25]+WP[1]*ev[31];
    ev[59]=PA[2]*ev[25]+WP[2]*ev[31];
    ev[60]=PA[0]*ev[26]+WP[0]*ev[32];
    ev[61]=PA[0]*ev[27]+WP[0]*ev[33];
    ev[62]=PA[0]*ev[30]+WP[0]*ev[36];
    ev[63]=PA[2]*ev[26]+WP[2]*ev[32];
    ev[64]=PA[1]*ev[27]+WP[1]*ev[33];
    ev[65]=PA[0]*ev[31]+WP[0]*ev[37]+2.e0*zeta2*(ev[10]-rz*ev[13]);
    ev[66]=PA[1]*ev[32]+WP[1]*ev[38]+2.e0*zeta2*(ev[11]-rz*ev[14]);
    ev[67]=PA[2]*ev[33]+WP[2]*ev[39]+2.e0*zeta2*(ev[12]-rz*ev[15]);
    ev[68]=PA[1]*ev[31]+WP[1]*ev[37];
    ev[69]=PA[2]*ev[31]+WP[2]*ev[37];
    ev[70]=PA[0]*ev[32]+WP[0]*ev[38];
    ev[71]=PA[0]*ev[33]+WP[0]*ev[39];
    ev[72]=PA[0]*ev[36]+WP[0]*ev[42];
    ev[73]=PA[2]*ev[32]+WP[2]*ev[38];
    ev[74]=PA[1]*ev[33]+WP[1]*ev[39];
    ev[75]=PA[0]*ev[37]+WP[0]*ev[43]+2.e0*zeta2*(ev[13]-rz*ev[16]);
    ev[76]=PA[1]*ev[38]+WP[1]*ev[44]+2.e0*zeta2*(ev[14]-rz*ev[17]);
    ev[77]=PA[2]*ev[39]+WP[2]*ev[45]+2.e0*zeta2*(ev[15]-rz*ev[18]);
    ev[78]=PA[1]*ev[37]+WP[1]*ev[43];
    ev[79]=PA[2]*ev[37]+WP[2]*ev[43];
    ev[80]=PA[0]*ev[38]+WP[0]*ev[44];
    ev[81]=PA[0]*ev[39]+WP[0]*ev[45];
    ev[82]=PA[0]*ev[42]+WP[0]*ev[48];
    ev[83]=PA[2]*ev[38]+WP[2]*ev[44];
    ev[84]=PA[1]*ev[39]+WP[1]*ev[45];
    ev[85]=PA[0]*ev[43]+WP[0]*ev[49]+2.e0*zeta2*(ev[16]-rz*ev[19]);
    ev[86]=PA[1]*ev[44]+WP[1]*ev[50]+2.e0*zeta2*(ev[17]-rz*ev[20]);
    ev[87]=PA[2]*ev[45]+WP[2]*ev[51]+2.e0*zeta2*(ev[18]-rz*ev[21]);
    ev[88]=PA[1]*ev[43]+WP[1]*ev[49];
    ev[89]=PA[2]*ev[43]+WP[2]*ev[49];
    ev[90]=PA[0]*ev[44]+WP[0]*ev[50];
    ev[91]=PA[0]*ev[45]+WP[0]*ev[51];
    ev[92]=PA[0]*ev[48]+WP[0]*ev[54];
    ev[93]=PA[2]*ev[44]+WP[2]*ev[50];
    ev[94]=PA[1]*ev[45]+WP[1]*ev[51];
    // (gs|ss) m=0,2
    ev[ 95]=PA[0]*ev[55]+WP[0]*ev[65]+3.e0*zeta2*(ev[25]-rz*ev[31]);
    ev[ 96]=PA[1]*ev[56]+WP[1]*ev[66]+3.e0*zeta2*(ev[26]-rz*ev[32]);
    ev[ 97]=PA[2]*ev[57]+WP[2]*ev[67]+3.e0*zeta2*(ev[27]-rz*ev[33]);
    ev[ 98]=PA[1]*ev[55]+WP[1]*ev[65];
    ev[ 99]=PA[2]*ev[55]+WP[2]*ev[65];
    ev[100]=PA[0]*ev[60]+WP[0]*ev[70]+     zeta2*(ev[26]-rz*ev[32]);
    ev[101]=PA[0]*ev[61]+WP[0]*ev[71]+     zeta2*(ev[27]-rz*ev[33]);
    ev[102]=PA[1]*ev[59]+WP[1]*ev[69];
    ev[103]=PA[0]*ev[56]+WP[0]*ev[66];
    ev[104]=PA[0]*ev[57]+WP[0]*ev[67];
    ev[105]=PA[0]*ev[63]+WP[0]*ev[73];
    ev[106]=PA[0]*ev[64]+WP[0]*ev[74];
    ev[107]=PA[2]*ev[56]+WP[2]*ev[66];
    ev[108]=PA[1]*ev[57]+WP[1]*ev[67];
    ev[109]=PA[1]*ev[64]+WP[1]*ev[74]+     zeta2*(ev[27]-rz*ev[33]);
    ev[110]=PA[0]*ev[65]+WP[0]*ev[75]+3.e0*zeta2*(ev[31]-rz*ev[37]);
    ev[111]=PA[1]*ev[66]+WP[1]*ev[76]+3.e0*zeta2*(ev[32]-rz*ev[38]);
    ev[112]=PA[2]*ev[67]+WP[2]*ev[77]+3.e0*zeta2*(ev[33]-rz*ev[39]);
    ev[113]=PA[1]*ev[65]+WP[1]*ev[75];
    ev[114]=PA[2]*ev[65]+WP[2]*ev[75];
    ev[115]=PA[0]*ev[70]+WP[0]*ev[80]+     zeta2*(ev[32]-rz*ev[38]);
    ev[116]=PA[0]*ev[71]+WP[0]*ev[81]+     zeta2*(ev[33]-rz*ev[39]);
    ev[117]=PA[1]*ev[69]+WP[1]*ev[79];
    ev[118]=PA[0]*ev[66]+WP[0]*ev[76];
    ev[119]=PA[0]*ev[67]+WP[0]*ev[77];
    ev[120]=PA[0]*ev[73]+WP[0]*ev[83];
    ev[121]=PA[0]*ev[74]+WP[0]*ev[84];
    ev[122]=PA[2]*ev[66]+WP[2]*ev[76];
    ev[123]=PA[1]*ev[67]+WP[1]*ev[77];
    ev[124]=PA[1]*ev[74]+WP[1]*ev[84]+     zeta2*(ev[33]-rz*ev[39]);
    ev[125]=PA[0]*ev[75]+WP[0]*ev[85]+3.e0*zeta2*(ev[37]-rz*ev[43]);
    ev[126]=PA[1]*ev[76]+WP[1]*ev[86]+3.e0*zeta2*(ev[38]-rz*ev[44]);
    ev[127]=PA[2]*ev[77]+WP[2]*ev[87]+3.e0*zeta2*(ev[39]-rz*ev[45]);
    ev[128]=PA[1]*ev[75]+WP[1]*ev[85];
    ev[129]=PA[2]*ev[75]+WP[2]*ev[85];
    ev[130]=PA[0]*ev[80]+WP[0]*ev[90]+     zeta2*(ev[38]-rz*ev[44]);
    ev[131]=PA[0]*ev[81]+WP[0]*ev[91]+     zeta2*(ev[39]-rz*ev[45]);
    ev[132]=PA[1]*ev[79]+WP[1]*ev[89];
    ev[133]=PA[0]*ev[76]+WP[0]*ev[86];
    ev[134]=PA[0]*ev[77]+WP[0]*ev[87];
    ev[135]=PA[0]*ev[83]+WP[0]*ev[93];
    ev[136]=PA[0]*ev[84]+WP[0]*ev[94];
    ev[137]=PA[2]*ev[76]+WP[2]*ev[86];
    ev[138]=PA[1]*ev[77]+WP[1]*ev[87];
    ev[139]=PA[1]*ev[84]+WP[1]*ev[94]+     zeta2*(ev[39]-rz*ev[45]);
    // (ps|ps) m=[1,1]
    ev[140]=QC[0]*ev[10]+WQ[0]*ev[13]+ze2*ev[2];
    ev[141]=QC[1]*ev[10]+WQ[1]*ev[13];
    ev[142]=QC[2]*ev[10]+WQ[2]*ev[13];
    ev[143]=QC[0]*ev[11]+WQ[0]*ev[14];
    ev[144]=QC[1]*ev[11]+WQ[1]*ev[14]+ze2*ev[2];
    ev[145]=QC[2]*ev[11]+WQ[2]*ev[14];
    ev[146]=QC[0]*ev[12]+WQ[0]*ev[15];
    ev[147]=QC[1]*ev[12]+WQ[1]*ev[15];
    ev[148]=QC[2]*ev[12]+WQ[2]*ev[15]+ze2*ev[2];
    // (ds|ps) m=[0,1]
    ev[149]=QC[0]*ev[25]+WQ[0]*ev[31]+2.e0*ze2*ev[10];
    ev[150]=QC[1]*ev[25]+WQ[1]*ev[31];
    ev[151]=QC[2]*ev[25]+WQ[2]*ev[31];
    ev[152]=QC[0]*ev[26]+WQ[0]*ev[32];
    ev[153]=QC[1]*ev[26]+WQ[1]*ev[32]+2.e0*ze2*ev[11];
    ev[154]=QC[2]*ev[26]+WQ[2]*ev[32];
    ev[155]=QC[0]*ev[27]+WQ[0]*ev[33];
    ev[156]=QC[1]*ev[27]+WQ[1]*ev[33];
    ev[157]=QC[2]*ev[27]+WQ[2]*ev[33]+2.e0*ze2*ev[12];
    ev[158]=QC[0]*ev[28]+WQ[0]*ev[34]+     ze2*ev[11];
    ev[159]=QC[1]*ev[28]+WQ[1]*ev[34]+     ze2*ev[10];
    ev[160]=QC[2]*ev[28]+WQ[2]*ev[34];
    ev[161]=QC[0]*ev[29]+WQ[0]*ev[35]+     ze2*ev[12];
    ev[162]=QC[1]*ev[29]+WQ[1]*ev[35];
    ev[163]=QC[2]*ev[29]+WQ[2]*ev[35]+     ze2*ev[10];
    ev[164]=QC[0]*ev[30]+WQ[0]*ev[36];
    ev[165]=QC[1]*ev[30]+WQ[1]*ev[36]+     ze2*ev[12];
    ev[166]=QC[2]*ev[30]+WQ[2]*ev[36]+     ze2*ev[11];
    ev[167]=QC[0]*ev[31]+WQ[0]*ev[37]+2.e0*ze2*ev[13];
    ev[168]=QC[1]*ev[31]+WQ[1]*ev[37];
    ev[169]=QC[2]*ev[31]+WQ[2]*ev[37];
    ev[170]=QC[0]*ev[32]+WQ[0]*ev[38];
    ev[171]=QC[1]*ev[32]+WQ[1]*ev[38]+2.e0*ze2*ev[14];
    ev[172]=QC[2]*ev[32]+WQ[2]*ev[38];
    ev[173]=QC[0]*ev[33]+WQ[0]*ev[39];
    ev[174]=QC[1]*ev[33]+WQ[1]*ev[39];
    ev[175]=QC[2]*ev[33]+WQ[2]*ev[39]+2.e0*ze2*ev[15];
    ev[176]=QC[0]*ev[34]+WQ[0]*ev[40]+     ze2*ev[14];
    ev[177]=QC[1]*ev[34]+WQ[1]*ev[40]+     ze2*ev[13];
    ev[178]=QC[2]*ev[34]+WQ[2]*ev[40];
    ev[179]=QC[0]*ev[35]+WQ[0]*ev[41]+     ze2*ev[15];
    ev[180]=QC[1]*ev[35]+WQ[1]*ev[41];
    ev[181]=QC[2]*ev[35]+WQ[2]*ev[41]+     ze2*ev[13];
    ev[182]=QC[0]*ev[36]+WQ[0]*ev[42];
    ev[183]=QC[1]*ev[36]+WQ[1]*ev[42]+     ze2*ev[15];
    ev[184]=QC[2]*ev[36]+WQ[2]*ev[42]+     ze2*ev[14];
    // (fs|ps) m=[0,1]
    ev[185]=QC[0]*ev[55]+WQ[0]*ev[65]+3.e0*ze2*ev[31];
    ev[186]=QC[1]*ev[55]+WQ[1]*ev[65];
    ev[187]=QC[2]*ev[55]+WQ[2]*ev[65];
    ev[188]=QC[0]*ev[56]+WQ[0]*ev[66];
    ev[189]=QC[1]*ev[56]+WQ[1]*ev[66]+3.e0*ze2*ev[32];
    ev[190]=QC[2]*ev[56]+WQ[2]*ev[66];
    ev[191]=QC[0]*ev[57]+WQ[0]*ev[67];
    ev[192]=QC[1]*ev[57]+WQ[1]*ev[67];
    ev[193]=QC[2]*ev[57]+WQ[2]*ev[67]+3.e0*ze2*ev[33];
    ev[194]=QC[0]*ev[58]+WQ[0]*ev[68]+2.e0*ze2*ev[34];
    ev[195]=QC[1]*ev[58]+WQ[1]*ev[68]+     ze2*ev[31];
    ev[196]=QC[2]*ev[58]+WQ[2]*ev[68];
    ev[197]=QC[0]*ev[59]+WQ[0]*ev[69]+2.e0*ze2*ev[35];
    ev[198]=QC[1]*ev[59]+WQ[1]*ev[69];
    ev[199]=QC[2]*ev[59]+WQ[2]*ev[69]+     ze2*ev[31];
    ev[200]=QC[0]*ev[60]+WQ[0]*ev[70]+     ze2*ev[32];
    ev[201]=QC[1]*ev[60]+WQ[1]*ev[70]+2.e0*ze2*ev[34];
    ev[202]=QC[2]*ev[60]+WQ[2]*ev[70];
    ev[203]=QC[0]*ev[61]+WQ[0]*ev[71]+     ze2*ev[33];
    ev[204]=QC[1]*ev[61]+WQ[1]*ev[71];
    ev[205]=QC[2]*ev[61]+WQ[2]*ev[71]+2.e0*ze2*ev[35];
    ev[206]=QC[0]*ev[62]+WQ[0]*ev[72]+     ze2*ev[36];
    ev[207]=QC[1]*ev[62]+WQ[1]*ev[72]+     ze2*ev[35];
    ev[208]=QC[2]*ev[62]+WQ[2]*ev[72]+     ze2*ev[34];
    ev[209]=QC[0]*ev[63]+WQ[0]*ev[73];
    ev[210]=QC[1]*ev[63]+WQ[1]*ev[73]+2.e0*ze2*ev[36];
    ev[211]=QC[2]*ev[63]+WQ[2]*ev[73]+     ze2*ev[32];
    ev[212]=QC[0]*ev[64]+WQ[0]*ev[74];
    ev[213]=QC[1]*ev[64]+WQ[1]*ev[74]+     ze2*ev[33];
    ev[214]=QC[2]*ev[64]+WQ[2]*ev[74]+2.e0*ze2*ev[36];
    ev[215]=QC[0]*ev[65]+WQ[0]*ev[75]+3.e0*ze2*ev[37];
    ev[216]=QC[1]*ev[65]+WQ[1]*ev[75];
    ev[217]=QC[2]*ev[65]+WQ[2]*ev[75];
    ev[218]=QC[0]*ev[66]+WQ[0]*ev[76];
    ev[219]=QC[1]*ev[66]+WQ[1]*ev[76]+3.e0*ze2*ev[38];
    ev[220]=QC[2]*ev[66]+WQ[2]*ev[76];
    ev[221]=QC[0]*ev[67]+WQ[0]*ev[77];
    ev[222]=QC[1]*ev[67]+WQ[1]*ev[77];
    ev[223]=QC[2]*ev[67]+WQ[2]*ev[77]+3.e0*ze2*ev[39];
    ev[224]=QC[0]*ev[68]+WQ[0]*ev[78]+2.e0*ze2*ev[40];
    ev[225]=QC[1]*ev[68]+WQ[1]*ev[78]+     ze2*ev[37];
    ev[226]=QC[2]*ev[68]+WQ[2]*ev[78];
    ev[227]=QC[0]*ev[69]+WQ[0]*ev[79]+2.e0*ze2*ev[41];
    ev[228]=QC[1]*ev[69]+WQ[1]*ev[79];
    ev[229]=QC[2]*ev[69]+WQ[2]*ev[79]+     ze2*ev[37];
    ev[230]=QC[0]*ev[70]+WQ[0]*ev[80]+     ze2*ev[38];
    ev[231]=QC[1]*ev[70]+WQ[1]*ev[80]+2.e0*ze2*ev[40];
    ev[232]=QC[2]*ev[70]+WQ[2]*ev[80];
    ev[233]=QC[0]*ev[71]+WQ[0]*ev[81]+     ze2*ev[39];
    ev[234]=QC[1]*ev[71]+WQ[1]*ev[81];
    ev[235]=QC[2]*ev[71]+WQ[2]*ev[81]+2.e0*ze2*ev[41];
    ev[236]=QC[0]*ev[72]+WQ[0]*ev[82]+     ze2*ev[42];
    ev[237]=QC[1]*ev[72]+WQ[1]*ev[82]+     ze2*ev[41];
    ev[238]=QC[2]*ev[72]+WQ[2]*ev[82]+     ze2*ev[40];
    ev[239]=QC[0]*ev[73]+WQ[0]*ev[83];
    ev[240]=QC[1]*ev[73]+WQ[1]*ev[83]+2.e0*ze2*ev[42];
    ev[241]=QC[2]*ev[73]+WQ[2]*ev[83]+     ze2*ev[38];
    ev[242]=QC[0]*ev[74]+WQ[0]*ev[84];
    ev[243]=QC[1]*ev[74]+WQ[1]*ev[84]+     ze2*ev[39];
    ev[244]=QC[2]*ev[74]+WQ[2]*ev[84]+2.e0*ze2*ev[42];
    // (gs|ps) m=[0,1]
    ev[245]=QC[0]*ev[ 95]+WQ[0]*ev[110]+4.e0*ze2*ev[65];
    ev[246]=QC[1]*ev[ 95]+WQ[1]*ev[110];
    ev[247]=QC[2]*ev[ 95]+WQ[2]*ev[110];
    ev[248]=QC[0]*ev[ 96]+WQ[0]*ev[111];
    ev[249]=QC[1]*ev[ 96]+WQ[1]*ev[111]+4.e0*ze2*ev[66];
    ev[250]=QC[2]*ev[ 96]+WQ[2]*ev[111];
    ev[251]=QC[0]*ev[ 97]+WQ[0]*ev[112];
    ev[252]=QC[1]*ev[ 97]+WQ[1]*ev[112];
    ev[253]=QC[2]*ev[ 97]+WQ[2]*ev[112]+4.e0*ze2*ev[67];
    ev[254]=QC[0]*ev[ 98]+WQ[0]*ev[113]+3.e0*ze2*ev[68];
    ev[255]=QC[1]*ev[ 98]+WQ[1]*ev[113]+     ze2*ev[65];
    ev[256]=QC[2]*ev[ 98]+WQ[2]*ev[113];
    ev[257]=QC[0]*ev[ 99]+WQ[0]*ev[114]+3.e0*ze2*ev[69];
    ev[258]=QC[1]*ev[ 99]+WQ[1]*ev[114];
    ev[259]=QC[2]*ev[ 99]+WQ[2]*ev[114]+     ze2*ev[65];
    ev[260]=QC[0]*ev[100]+WQ[0]*ev[115]+2.e0*ze2*ev[70];
    ev[261]=QC[1]*ev[100]+WQ[1]*ev[115]+2.e0*ze2*ev[68];
    ev[262]=QC[2]*ev[100]+WQ[2]*ev[115];
    ev[263]=QC[0]*ev[101]+WQ[0]*ev[116]+2.e0*ze2*ev[71];
    ev[264]=QC[1]*ev[101]+WQ[1]*ev[116];
    ev[265]=QC[2]*ev[101]+WQ[2]*ev[116]+2.e0*ze2*ev[69];
    ev[266]=QC[0]*ev[102]+WQ[0]*ev[117]+2.e0*ze2*ev[72];
    ev[267]=QC[1]*ev[102]+WQ[1]*ev[117]+     ze2*ev[69];
    ev[268]=QC[2]*ev[102]+WQ[2]*ev[117]+     ze2*ev[68];
    ev[269]=QC[0]*ev[103]+WQ[0]*ev[118]+     ze2*ev[66];
    ev[270]=QC[1]*ev[103]+WQ[1]*ev[118]+3.e0*ze2*ev[70];
    ev[271]=QC[2]*ev[103]+WQ[2]*ev[118];
    ev[272]=QC[0]*ev[104]+WQ[0]*ev[119]+     ze2*ev[67];
    ev[273]=QC[1]*ev[104]+WQ[1]*ev[119];
    ev[274]=QC[2]*ev[104]+WQ[2]*ev[119]+3.e0*ze2*ev[71];
    ev[275]=QC[0]*ev[105]+WQ[0]*ev[120]+     ze2*ev[73];
    ev[276]=QC[1]*ev[105]+WQ[1]*ev[120]+2.e0*ze2*ev[72];
    ev[277]=QC[2]*ev[105]+WQ[2]*ev[120]+     ze2*ev[70];
    ev[278]=QC[0]*ev[106]+WQ[0]*ev[121]+     ze2*ev[74];
    ev[279]=QC[1]*ev[106]+WQ[1]*ev[121]+     ze2*ev[71];
    ev[280]=QC[2]*ev[106]+WQ[2]*ev[121]+2.e0*ze2*ev[72];
    ev[281]=QC[0]*ev[107]+WQ[0]*ev[122];
    ev[282]=QC[1]*ev[107]+WQ[1]*ev[122]+3.e0*ze2*ev[73];
    ev[283]=QC[2]*ev[107]+WQ[2]*ev[122]+     ze2*ev[66];
    ev[284]=QC[0]*ev[108]+WQ[0]*ev[123];
    ev[285]=QC[1]*ev[108]+WQ[1]*ev[123]+     ze2*ev[67];
    ev[286]=QC[2]*ev[108]+WQ[2]*ev[123]+3.e0*ze2*ev[74];
    ev[287]=QC[0]*ev[109]+WQ[0]*ev[124];
    ev[288]=QC[1]*ev[109]+WQ[1]*ev[124]+2.e0*ze2*ev[74];
    ev[289]=QC[2]*ev[109]+WQ[2]*ev[124]+2.e0*ze2*ev[73];
    ev[290]=QC[0]*ev[110]+WQ[0]*ev[125]+4.e0*ze2*ev[75];
    ev[291]=QC[1]*ev[110]+WQ[1]*ev[125];
    ev[292]=QC[2]*ev[110]+WQ[2]*ev[125];
    ev[293]=QC[0]*ev[111]+WQ[0]*ev[126];
    ev[294]=QC[1]*ev[111]+WQ[1]*ev[126]+4.e0*ze2*ev[76];
    ev[295]=QC[2]*ev[111]+WQ[2]*ev[126];
    ev[296]=QC[0]*ev[112]+WQ[0]*ev[127];
    ev[297]=QC[1]*ev[112]+WQ[1]*ev[127];
    ev[298]=QC[2]*ev[112]+WQ[2]*ev[127]+4.e0*ze2*ev[77];
    ev[299]=QC[0]*ev[113]+WQ[0]*ev[128]+3.e0*ze2*ev[78];
    ev[300]=QC[1]*ev[113]+WQ[1]*ev[128]+     ze2*ev[75];
    ev[301]=QC[2]*ev[113]+WQ[2]*ev[128];
    ev[302]=QC[0]*ev[114]+WQ[0]*ev[129]+3.e0*ze2*ev[79];
    ev[303]=QC[1]*ev[114]+WQ[1]*ev[129];
    ev[304]=QC[2]*ev[114]+WQ[2]*ev[129]+     ze2*ev[75];
    ev[305]=QC[0]*ev[115]+WQ[0]*ev[130]+2.e0*ze2*ev[80];
    ev[306]=QC[1]*ev[115]+WQ[1]*ev[130]+2.e0*ze2*ev[78];
    ev[307]=QC[2]*ev[115]+WQ[2]*ev[130];
    ev[308]=QC[0]*ev[116]+WQ[0]*ev[131]+2.e0*ze2*ev[81];
    ev[309]=QC[1]*ev[116]+WQ[1]*ev[131];
    ev[310]=QC[2]*ev[116]+WQ[2]*ev[131]+2.e0*ze2*ev[79];
    ev[311]=QC[0]*ev[117]+WQ[0]*ev[132]+2.e0*ze2*ev[82];
    ev[312]=QC[1]*ev[117]+WQ[1]*ev[132]+     ze2*ev[79];
    ev[313]=QC[2]*ev[117]+WQ[2]*ev[132]+     ze2*ev[78];
    ev[314]=QC[0]*ev[118]+WQ[0]*ev[133]+     ze2*ev[76];
    ev[315]=QC[1]*ev[118]+WQ[1]*ev[133]+3.e0*ze2*ev[80];
    ev[316]=QC[2]*ev[118]+WQ[2]*ev[133];
    ev[317]=QC[0]*ev[119]+WQ[0]*ev[134]+     ze2*ev[77];
    ev[318]=QC[1]*ev[119]+WQ[1]*ev[134];
    ev[319]=QC[2]*ev[119]+WQ[2]*ev[134]+3.e0*ze2*ev[81];
    ev[320]=QC[0]*ev[120]+WQ[0]*ev[135]+     ze2*ev[83];
    ev[321]=QC[1]*ev[120]+WQ[1]*ev[135]+2.e0*ze2*ev[82];
    ev[322]=QC[2]*ev[120]+WQ[2]*ev[135]+     ze2*ev[80];
    ev[323]=QC[0]*ev[121]+WQ[0]*ev[136]+     ze2*ev[84];
    ev[324]=QC[1]*ev[121]+WQ[1]*ev[136]+     ze2*ev[81];
    ev[325]=QC[2]*ev[121]+WQ[2]*ev[136]+2.e0*ze2*ev[82];
    ev[326]=QC[0]*ev[122]+WQ[0]*ev[137];
    ev[327]=QC[1]*ev[122]+WQ[1]*ev[137]+3.e0*ze2*ev[83];
    ev[328]=QC[2]*ev[122]+WQ[2]*ev[137]+     ze2*ev[76];
    ev[329]=QC[0]*ev[123]+WQ[0]*ev[138];
    ev[330]=QC[1]*ev[123]+WQ[1]*ev[138]+     ze2*ev[77];
    ev[331]=QC[2]*ev[123]+WQ[2]*ev[138]+3.e0*ze2*ev[84];
    ev[332]=QC[0]*ev[124]+WQ[0]*ev[139];
    ev[333]=QC[1]*ev[124]+WQ[1]*ev[139]+2.e0*ze2*ev[84];
    ev[334]=QC[2]*ev[124]+WQ[2]*ev[139]+2.e0*ze2*ev[83];
    // (ds|ds) m=[0,0]
    ev[335]=QC[0]*ev[149]+WQ[0]*ev[167]+eta2*(ev[25]-re*ev[31])
            +2.e0*ze2*ev[140];
    ev[336]=QC[1]*ev[150]+WQ[1]*ev[168]+eta2*(ev[25]-re*ev[31]);
    ev[337]=QC[2]*ev[151]+WQ[2]*ev[169]+eta2*(ev[25]-re*ev[31]);
    ev[338]=QC[0]*ev[150]+WQ[0]*ev[168]+2.e0*ze2*ev[141];
    ev[339]=QC[0]*ev[151]+WQ[0]*ev[169]+2.e0*ze2*ev[142];
    ev[340]=QC[1]*ev[151]+WQ[1]*ev[169];
    ev[341]=QC[0]*ev[152]+WQ[0]*ev[170]+eta2*(ev[26]-re*ev[32]);
    ev[342]=QC[1]*ev[153]+WQ[1]*ev[171]+eta2*(ev[26]-re*ev[32])
            +2.e0*ze2*ev[144];
    ev[343]=QC[2]*ev[154]+WQ[2]*ev[172]+eta2*(ev[26]-re*ev[32]);
    ev[344]=QC[0]*ev[153]+WQ[0]*ev[171];
    ev[345]=QC[0]*ev[154]+WQ[0]*ev[172];
    ev[346]=QC[1]*ev[154]+WQ[1]*ev[172]+2.e0*ze2*ev[145];
    ev[347]=QC[0]*ev[155]+WQ[0]*ev[173]+eta2*(ev[27]-re*ev[33]);
    ev[348]=QC[1]*ev[156]+WQ[1]*ev[174]+eta2*(ev[27]-re*ev[33]);
    ev[349]=QC[2]*ev[157]+WQ[2]*ev[175]+eta2*(ev[27]-re*ev[33])
            +2.e0*ze2*ev[148];
    ev[350]=QC[0]*ev[156]+WQ[0]*ev[174];
    ev[351]=QC[0]*ev[157]+WQ[0]*ev[175];
    ev[352]=QC[1]*ev[157]+WQ[1]*ev[175];
    ev[353]=QC[0]*ev[158]+WQ[0]*ev[176]+eta2*(ev[28]-re*ev[34])
            +     ze2*ev[143];
    ev[354]=QC[1]*ev[159]+WQ[1]*ev[177]+eta2*(ev[28]-re*ev[34])
            +     ze2*ev[141];
    ev[355]=QC[2]*ev[160]+WQ[2]*ev[178]+eta2*(ev[28]-re*ev[34]);
    ev[356]=QC[0]*ev[159]+WQ[0]*ev[177]+     ze2*ev[144];
    ev[357]=QC[0]*ev[160]+WQ[0]*ev[178]+     ze2*ev[145];
    ev[358]=QC[1]*ev[160]+WQ[1]*ev[178]+     ze2*ev[142];
    ev[359]=QC[0]*ev[161]+WQ[0]*ev[179]+eta2*(ev[29]-re*ev[35])
            +     ze2*ev[146];
    ev[360]=QC[1]*ev[162]+WQ[1]*ev[180]+eta2*(ev[29]-re*ev[35]);
    ev[361]=QC[2]*ev[163]+WQ[2]*ev[181]+eta2*(ev[29]-re*ev[35])
            +     ze2*ev[142];
    ev[362]=QC[0]*ev[162]+WQ[0]*ev[180]+     ze2*ev[147];
    ev[363]=QC[0]*ev[163]+WQ[0]*ev[181]+     ze2*ev[148];
    ev[364]=QC[1]*ev[163]+WQ[1]*ev[181];
    ev[365]=QC[0]*ev[164]+WQ[0]*ev[182]+eta2*(ev[30]-re*ev[36]);
    ev[366]=QC[1]*ev[165]+WQ[1]*ev[183]+eta2*(ev[30]-re*ev[36])
            +     ze2*ev[147];
    ev[367]=QC[2]*ev[166]+WQ[2]*ev[184]+eta2*(ev[30]-re*ev[36])
            +     ze2*ev[145];
    ev[368]=QC[0]*ev[165]+WQ[0]*ev[183];
    ev[369]=QC[0]*ev[166]+WQ[0]*ev[184];
    ev[370]=QC[1]*ev[166]+WQ[1]*ev[184]+     ze2*ev[148];
    // (fs|ds) m=[0,0]
    ev[371]=QC[0]*ev[185]+WQ[0]*ev[215]+eta2*(ev[55]-re*ev[65])
            +3.e0*ze2*ev[167];
    ev[372]=QC[1]*ev[186]+WQ[1]*ev[216]+eta2*(ev[55]-re*ev[65]);
    ev[373]=QC[2]*ev[187]+WQ[2]*ev[217]+eta2*(ev[55]-re*ev[65]);
    ev[374]=QC[0]*ev[186]+WQ[0]*ev[216]+3.e0*ze2*ev[168];
    ev[375]=QC[0]*ev[187]+WQ[0]*ev[217]+3.e0*ze2*ev[169];
    ev[376]=QC[1]*ev[187]+WQ[1]*ev[217];
    ev[377]=QC[0]*ev[188]+WQ[0]*ev[218]+eta2*(ev[56]-re*ev[66]);
    ev[378]=QC[1]*ev[189]+WQ[1]*ev[219]+eta2*(ev[56]-re*ev[66])
            +3.e0*ze2*ev[171];
    ev[379]=QC[2]*ev[190]+WQ[2]*ev[220]+eta2*(ev[56]-re*ev[66]);
    ev[380]=QC[0]*ev[189]+WQ[0]*ev[219];
    ev[381]=QC[0]*ev[190]+WQ[0]*ev[220];
    ev[382]=QC[1]*ev[190]+WQ[1]*ev[220]+3.e0*ze2*ev[172];
    ev[383]=QC[0]*ev[191]+WQ[0]*ev[221]+eta2*(ev[57]-re*ev[67]);
    ev[384]=QC[1]*ev[192]+WQ[1]*ev[222]+eta2*(ev[57]-re*ev[67]);
    ev[385]=QC[2]*ev[193]+WQ[2]*ev[223]+eta2*(ev[57]-re*ev[67])
            +3.e0*ze2*ev[175];
    ev[386]=QC[0]*ev[192]+WQ[0]*ev[222];
    ev[387]=QC[0]*ev[193]+WQ[0]*ev[223];
    ev[388]=QC[1]*ev[193]+WQ[1]*ev[223];
    ev[389]=QC[0]*ev[194]+WQ[0]*ev[224]+eta2*(ev[58]-re*ev[68])
            +2.e0*ze2*ev[176];
    ev[390]=QC[1]*ev[195]+WQ[1]*ev[225]+eta2*(ev[58]-re*ev[68])
            +     ze2*ev[168];
    ev[391]=QC[2]*ev[196]+WQ[2]*ev[226]+eta2*(ev[58]-re*ev[68]);
    ev[392]=QC[0]*ev[195]+WQ[0]*ev[225]+2.e0*ze2*ev[177];
    ev[393]=QC[0]*ev[196]+WQ[0]*ev[226]+2.e0*ze2*ev[178];
    ev[394]=QC[1]*ev[196]+WQ[1]*ev[226]+     ze2*ev[169];
    ev[395]=QC[0]*ev[197]+WQ[0]*ev[227]+eta2*(ev[59]-re*ev[69])
            +2.e0*ze2*ev[179];
    ev[396]=QC[1]*ev[198]+WQ[1]*ev[228]+eta2*(ev[59]-re*ev[69]);
    ev[397]=QC[2]*ev[199]+WQ[2]*ev[229]+eta2*(ev[59]-re*ev[69])
            +     ze2*ev[169];
    ev[398]=QC[0]*ev[198]+WQ[0]*ev[228]+2.e0*ze2*ev[180];
    ev[399]=QC[0]*ev[199]+WQ[0]*ev[229]+2.e0*ze2*ev[181];
    ev[400]=QC[1]*ev[199]+WQ[1]*ev[229];
    ev[401]=QC[0]*ev[200]+WQ[0]*ev[230]+eta2*(ev[60]-re*ev[70])
            +     ze2*ev[170];
    ev[402]=QC[1]*ev[201]+WQ[1]*ev[231]+eta2*(ev[60]-re*ev[70])
            +2.e0*ze2*ev[177];
    ev[403]=QC[2]*ev[202]+WQ[2]*ev[232]+eta2*(ev[60]-re*ev[70]);
    ev[404]=QC[0]*ev[201]+WQ[0]*ev[231]+     ze2*ev[171];
    ev[405]=QC[0]*ev[202]+WQ[0]*ev[232]+     ze2*ev[172];
    ev[406]=QC[1]*ev[202]+WQ[1]*ev[232]+2.e0*ze2*ev[178];
    ev[407]=QC[0]*ev[203]+WQ[0]*ev[233]+eta2*(ev[61]-re*ev[71])
            +     ze2*ev[173];
    ev[408]=QC[1]*ev[204]+WQ[1]*ev[234]+eta2*(ev[61]-re*ev[71]);
    ev[409]=QC[2]*ev[205]+WQ[2]*ev[235]+eta2*(ev[61]-re*ev[71])
            +2.e0*ze2*ev[181];
    ev[410]=QC[0]*ev[204]+WQ[0]*ev[234]+     ze2*ev[174];
    ev[411]=QC[0]*ev[205]+WQ[0]*ev[235]+     ze2*ev[175];
    ev[412]=QC[1]*ev[205]+WQ[1]*ev[235];
    ev[413]=QC[0]*ev[206]+WQ[0]*ev[236]+eta2*(ev[62]-re*ev[72])
            +     ze2*ev[182];
    ev[414]=QC[1]*ev[207]+WQ[1]*ev[237]+eta2*(ev[62]-re*ev[72])
            +     ze2*ev[180];
    ev[415]=QC[2]*ev[208]+WQ[2]*ev[238]+eta2*(ev[62]-re*ev[72])
            +     ze2*ev[178];
    ev[416]=QC[0]*ev[207]+WQ[0]*ev[237]+     ze2*ev[183];
    ev[417]=QC[0]*ev[208]+WQ[0]*ev[238]+     ze2*ev[184];
    ev[418]=QC[1]*ev[208]+WQ[1]*ev[238]+     ze2*ev[181];
    ev[419]=QC[0]*ev[209]+WQ[0]*ev[239]+eta2*(ev[63]-re*ev[73]);
    ev[420]=QC[1]*ev[210]+WQ[1]*ev[240]+eta2*(ev[63]-re*ev[73])
            +2.e0*ze2*ev[183];
    ev[421]=QC[2]*ev[211]+WQ[2]*ev[241]+eta2*(ev[63]-re*ev[73])
            +     ze2*ev[172];
    ev[422]=QC[0]*ev[210]+WQ[0]*ev[240];
    ev[423]=QC[0]*ev[211]+WQ[0]*ev[241];
    ev[424]=QC[1]*ev[211]+WQ[1]*ev[241]+2.e0*ze2*ev[184];
    ev[425]=QC[0]*ev[212]+WQ[0]*ev[242]+eta2*(ev[64]-re*ev[74]);
    ev[426]=QC[1]*ev[213]+WQ[1]*ev[243]+eta2*(ev[64]-re*ev[74])
            +     ze2*ev[174];
    ev[427]=QC[2]*ev[214]+WQ[2]*ev[244]+eta2*(ev[64]-re*ev[74])
            +2.e0*ze2*ev[184];
    ev[428]=QC[0]*ev[213]+WQ[0]*ev[243];
    ev[429]=QC[0]*ev[214]+WQ[0]*ev[244];
    ev[430]=QC[1]*ev[214]+WQ[1]*ev[244]+     ze2*ev[175];
    // (gs|ds) m=[0,0]
    ev[431]=QC[0]*ev[245]+WQ[0]*ev[290]+eta2*(ev[ 95]-re*ev[110])
            +4.e0*ze2*ev[215];
    ev[432]=QC[1]*ev[246]+WQ[1]*ev[291]+eta2*(ev[ 95]-re*ev[110]);
    ev[433]=QC[2]*ev[247]+WQ[2]*ev[292]+eta2*(ev[ 95]-re*ev[110]);
    ev[434]=QC[0]*ev[246]+WQ[0]*ev[291]+4.e0*ze2*ev[216];
    ev[435]=QC[0]*ev[247]+WQ[0]*ev[292]+4.e0*ze2*ev[217];
    ev[436]=QC[1]*ev[247]+WQ[1]*ev[292];
    ev[437]=QC[0]*ev[248]+WQ[0]*ev[293]+eta2*(ev[ 96]-re*ev[111]);
    ev[438]=QC[1]*ev[249]+WQ[1]*ev[294]+eta2*(ev[ 96]-re*ev[111])
            +4.e0*ze2*ev[219];
    ev[439]=QC[2]*ev[250]+WQ[2]*ev[295]+eta2*(ev[ 96]-re*ev[111]);
    ev[440]=QC[0]*ev[249]+WQ[0]*ev[294];
    ev[441]=QC[0]*ev[250]+WQ[0]*ev[295];
    ev[442]=QC[1]*ev[250]+WQ[1]*ev[295]+4.e0*ze2*ev[220];
    ev[443]=QC[0]*ev[251]+WQ[0]*ev[296]+eta2*(ev[ 97]-re*ev[112]);
    ev[444]=QC[1]*ev[252]+WQ[1]*ev[297]+eta2*(ev[ 97]-re*ev[112]);
    ev[445]=QC[2]*ev[253]+WQ[2]*ev[298]+eta2*(ev[ 97]-re*ev[112])
            +4.e0*ze2*ev[223];
    ev[446]=QC[0]*ev[252]+WQ[0]*ev[297];
    ev[447]=QC[0]*ev[253]+WQ[0]*ev[298];
    ev[448]=QC[1]*ev[253]+WQ[1]*ev[298];
    ev[449]=QC[0]*ev[254]+WQ[0]*ev[299]+eta2*(ev[ 98]-re*ev[113])
            +3.e0*ze2*ev[224];
    ev[450]=QC[1]*ev[255]+WQ[1]*ev[300]+eta2*(ev[ 98]-re*ev[113])
            +     ze2*ev[216];
    ev[451]=QC[2]*ev[256]+WQ[2]*ev[301]+eta2*(ev[ 98]-re*ev[113]);
    ev[452]=QC[0]*ev[255]+WQ[0]*ev[300]+3.e0*ze2*ev[225];
    ev[453]=QC[0]*ev[256]+WQ[0]*ev[301]+3.e0*ze2*ev[226];
    ev[454]=QC[1]*ev[256]+WQ[1]*ev[301]+     ze2*ev[217];
    ev[455]=QC[0]*ev[257]+WQ[0]*ev[302]+eta2*(ev[ 99]-re*ev[114])
            +3.e0*ze2*ev[227];
    ev[456]=QC[1]*ev[258]+WQ[1]*ev[303]+eta2*(ev[ 99]-re*ev[114]);
    ev[457]=QC[2]*ev[259]+WQ[2]*ev[304]+eta2*(ev[ 99]-re*ev[114])
            +     ze2*ev[217];
    ev[458]=QC[0]*ev[258]+WQ[0]*ev[303]+3.e0*ze2*ev[228];
    ev[459]=QC[0]*ev[259]+WQ[0]*ev[304]+3.e0*ze2*ev[229];
    ev[460]=QC[1]*ev[259]+WQ[1]*ev[304];
    ev[461]=QC[0]*ev[260]+WQ[0]*ev[305]+eta2*(ev[100]-re*ev[115])
            +2.e0*ze2*ev[230];
    ev[462]=QC[1]*ev[261]+WQ[1]*ev[306]+eta2*(ev[100]-re*ev[115])
            +2.e0*ze2*ev[225];
    ev[463]=QC[2]*ev[262]+WQ[2]*ev[307]+eta2*(ev[100]-re*ev[115]);
    ev[464]=QC[0]*ev[261]+WQ[0]*ev[306]+2.e0*ze2*ev[231];
    ev[465]=QC[0]*ev[262]+WQ[0]*ev[307]+2.e0*ze2*ev[232];
    ev[466]=QC[1]*ev[262]+WQ[1]*ev[307]+2.e0*ze2*ev[226];
    ev[467]=QC[0]*ev[263]+WQ[0]*ev[308]+eta2*(ev[101]-re*ev[116])
            +2.e0*ze2*ev[233];
    ev[468]=QC[1]*ev[264]+WQ[1]*ev[309]+eta2*(ev[101]-re*ev[116]);
    ev[469]=QC[2]*ev[265]+WQ[2]*ev[310]+eta2*(ev[101]-re*ev[116])
            +2.e0*ze2*ev[229];
    ev[470]=QC[0]*ev[264]+WQ[0]*ev[309]+2.e0*ze2*ev[234];
    ev[471]=QC[0]*ev[265]+WQ[0]*ev[310]+2.e0*ze2*ev[235];
    ev[472]=QC[1]*ev[265]+WQ[1]*ev[310];
    ev[473]=QC[0]*ev[266]+WQ[0]*ev[311]+eta2*(ev[102]-re*ev[117])
            +2.e0*ze2*ev[236];
    ev[474]=QC[1]*ev[267]+WQ[1]*ev[312]+eta2*(ev[102]-re*ev[117])
            +     ze2*ev[228];
    ev[475]=QC[2]*ev[268]+WQ[2]*ev[313]+eta2*(ev[102]-re*ev[117])
            +     ze2*ev[226];
    ev[476]=QC[0]*ev[267]+WQ[0]*ev[312]+2.e0*ze2*ev[237];
    ev[477]=QC[0]*ev[268]+WQ[0]*ev[313]+2.e0*ze2*ev[238];
    ev[478]=QC[1]*ev[268]+WQ[1]*ev[313]+     ze2*ev[229];
    ev[479]=QC[0]*ev[269]+WQ[0]*ev[314]+eta2*(ev[103]-re*ev[118])
            +     ze2*ev[218];
    ev[480]=QC[1]*ev[270]+WQ[1]*ev[315]+eta2*(ev[103]-re*ev[118])
            +3.e0*ze2*ev[231];
    ev[481]=QC[2]*ev[271]+WQ[2]*ev[316]+eta2*(ev[103]-re*ev[118]);
    ev[482]=QC[0]*ev[270]+WQ[0]*ev[315]+     ze2*ev[219];
    ev[483]=QC[0]*ev[271]+WQ[0]*ev[316]+     ze2*ev[220];
    ev[484]=QC[1]*ev[271]+WQ[1]*ev[316]+3.e0*ze2*ev[232];
    ev[485]=QC[0]*ev[272]+WQ[0]*ev[317]+eta2*(ev[104]-re*ev[119])
            +     ze2*ev[221];
    ev[486]=QC[1]*ev[273]+WQ[1]*ev[318]+eta2*(ev[104]-re*ev[119]);
    ev[487]=QC[2]*ev[274]+WQ[2]*ev[319]+eta2*(ev[104]-re*ev[119])
            +3.e0*ze2*ev[235];
    ev[488]=QC[0]*ev[273]+WQ[0]*ev[318]+     ze2*ev[222];
    ev[489]=QC[0]*ev[274]+WQ[0]*ev[319]+     ze2*ev[223];
    ev[490]=QC[1]*ev[274]+WQ[1]*ev[319];
    ev[491]=QC[0]*ev[275]+WQ[0]*ev[320]+eta2*(ev[105]-re*ev[120])
            +     ze2*ev[239];
    ev[492]=QC[1]*ev[276]+WQ[1]*ev[321]+eta2*(ev[105]-re*ev[120])
            +2.e0*ze2*ev[237];
    ev[493]=QC[2]*ev[277]+WQ[2]*ev[322]+eta2*(ev[105]-re*ev[120])
            +     ze2*ev[232];
    ev[494]=QC[0]*ev[276]+WQ[0]*ev[321]+     ze2*ev[240];
    ev[495]=QC[0]*ev[277]+WQ[0]*ev[322]+     ze2*ev[241];
    ev[496]=QC[1]*ev[277]+WQ[1]*ev[322]+2.e0*ze2*ev[238];
    ev[497]=QC[0]*ev[278]+WQ[0]*ev[323]+eta2*(ev[106]-re*ev[121])
            +     ze2*ev[242];
    ev[498]=QC[1]*ev[279]+WQ[1]*ev[324]+eta2*(ev[106]-re*ev[121])
            +     ze2*ev[234];
    ev[499]=QC[2]*ev[280]+WQ[2]*ev[325]+eta2*(ev[106]-re*ev[121])
            +2.e0*ze2*ev[238];
    ev[500]=QC[0]*ev[279]+WQ[0]*ev[324]+     ze2*ev[243];
    ev[501]=QC[0]*ev[280]+WQ[0]*ev[325]+     ze2*ev[244];
    ev[502]=QC[1]*ev[280]+WQ[1]*ev[325]+     ze2*ev[235];
    ev[503]=QC[0]*ev[281]+WQ[0]*ev[326]+eta2*(ev[107]-re*ev[122]);
    ev[504]=QC[1]*ev[282]+WQ[1]*ev[327]+eta2*(ev[107]-re*ev[122])
            +3.e0*ze2*ev[240];
    ev[505]=QC[2]*ev[283]+WQ[2]*ev[328]+eta2*(ev[107]-re*ev[122])
            +     ze2*ev[220];
    ev[506]=QC[0]*ev[282]+WQ[0]*ev[327];
    ev[507]=QC[0]*ev[283]+WQ[0]*ev[328];
    ev[508]=QC[1]*ev[283]+WQ[1]*ev[328]+3.e0*ze2*ev[241];
    ev[509]=QC[0]*ev[284]+WQ[0]*ev[329]+eta2*(ev[108]-re*ev[123]);
    ev[510]=QC[1]*ev[285]+WQ[1]*ev[330]+eta2*(ev[108]-re*ev[123])
            +     ze2*ev[222];
    ev[511]=QC[2]*ev[286]+WQ[2]*ev[331]+eta2*(ev[108]-re*ev[123])
            +3.e0*ze2*ev[244];
    ev[512]=QC[0]*ev[285]+WQ[0]*ev[330];
    ev[513]=QC[0]*ev[286]+WQ[0]*ev[331];
    ev[514]=QC[1]*ev[286]+WQ[1]*ev[331]+     ze2*ev[223];
    ev[515]=QC[0]*ev[287]+WQ[0]*ev[332]+eta2*(ev[109]-re*ev[124]);
    ev[516]=QC[1]*ev[288]+WQ[1]*ev[333]+eta2*(ev[109]-re*ev[124])
            +2.e0*ze2*ev[243];
    ev[517]=QC[2]*ev[289]+WQ[2]*ev[334]+eta2*(ev[109]-re*ev[124])
            +2.e0*ze2*ev[241];
    ev[518]=QC[0]*ev[288]+WQ[0]*ev[333];
    ev[519]=QC[0]*ev[289]+WQ[0]*ev[334];
    ev[520]=QC[1]*ev[289]+WQ[1]*ev[334]+2.e0*ze2*ev[244];
}

static void ofmo_vrr_cint_ddds( const double *ev, double *eh ) {
    int La=2, Lb=2, Lc=2, Ld=0;
    int i, ih, iv;
    // (DS|DS)
    for ( i=0, iv=335, ih=0; i<36; i++, iv++, ih++ ) eh[ih]+=ev[iv];
    // (FS|DS)
    for ( i=0, iv=371, ih=36; i<60; i++, iv++, ih++ ) eh[ih]+=ev[iv];
    // (GS|DS)
    for ( i=0, iv=431, ih=96; i<90; i++, iv++, ih++ ) eh[ih]+=ev[iv];
}

void ofmo_twoint_core_os_ddds(
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *nijps, const double vzeta[], const double vdkab[],
        const double vxiza[], const double BA[3],
        const int *nklps, const double veta[], const double vdkcd[],
        const double vxizc[], const double DC[3], const double AC[3],
        double *DINT ) {
    int ijps, klps, i;
    double cssss, zeta, dkab, xiza, eta, xizc, dk, T;
    double zeta2, eta2, ze2, rz, re, PA[3], WP[3], QC[3], WQ[3];
    double PQ2, sqrho, rho, PC[3], QP[3];
    double ev[521], eh[690];
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;

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
            sqrho = sqrt(1.e0/(zeta+eta));
            rho   = sqrho*sqrho;
            rz    = rho * zeta;
            re    = rho * eta;
            ze2 = rz * eta2;
            for ( i=0; i<3; i++ ) {
                WP[i] = rz*QP[i];
                WQ[i] = rz*QP[i] - QP[i];
            }
            T     = rho * PQ2;
            cssss = sqrho * dk;
            ofmo_vrr_calc_ddds(
                    T, cssss, zeta2, eta2, ze2, rz, re, PA, WP, QC, WQ,
                    ev );
            ofmo_vrr_cint_ddds( ev, eh );
        }	// for (klps)
    }	// for (ijps)
    ofmo_hrr_calc_ddds( BA, DC, eh );
    ofmo_hrr_coef_ddds( eh, DINT );
}

int ofmo_twoint_os_ddds(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int shel_atm[], const int shel_ini[],
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
    double DINTEG[6*6*6*1];
    long nzeri, max_nzeri, nzeri4;
    int nworkers=*pnworkers, workerid=*pworkerid;
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;
    long ebuf_max_nzeri = *pebuf_max_nzeri;
    int mythread;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);

    
    mythread = omp_get_thread_num();
    if ( DFACT == NULL ) DFACT = ofmo_getadd_dfact();
    
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
            if ( val_ab*val_cd < eps_ps4 ) continue;
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
            ofmo_twoint_core_os_ddds(
                    &La, &Lb, &Lc, &Ld,
                    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
                    &psp_xiza[ijps0], BA,
                    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
                    &psp_xiza[klps0], DC,   AC,      DINTEG );
            ipat=((Lab != Lcd)||(ics==kcs && jcs>lcs) ? true : false );
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
                        }	// l
                    }	// k
                }	// j
            }	// i
            if ( nzeri >= max_nzeri ) {
                *last_ijcs = ijcs;
                *last_klcs = klcs;
                *ebuf_non_zero_eri = nzeri;
                return OFMO_EBUF_FULL;
            }
        }	// for (klcs)
    }	// for (ijcs)
    *ebuf_non_zero_eri = nzeri;
    return OFMO_EBUF_NOFULL;
}

int ofmo_twoint_direct_os_ddds(
        const int *pnworkers, const int *pworkerid,
        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int shel_atm[], const int shel_ini[],
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
    double DINTEG[6*6*6*1];
    int mythread;
    float eps_eri = ofmo_twoint_eps_eri(0);
    float eps_ps4 = ofmo_twoint_eps_ps4(0);
    float eps_sch = ofmo_twoint_eps_sch(0);

    mythread = omp_get_thread_num();
    if ( DFACT == NULL ) DFACT = ofmo_getadd_dfact();
    
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
            ofmo_twoint_core_os_ddds(
                    &La, &Lb, &Lc, &Ld,
                    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
                    &psp_xiza[ijps0], BA,
                    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
                    &psp_xiza[klps0], DC,   AC,      DINTEG );
            ipat=((Lab != Lcd)||(ics==kcs && jcs>lcs) ? true : false );
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
                        }	// l
                    }	// k
                }	// j
            }	// i
            if ( nzeri >= max_nzeri ) {
                ofmo_integ_add_fock( nao, nzeri, etmp_val, etmp_ind4,
                        Ds, G );
                nzeri = nzeri4= 0;
            }
        }	// for (klcs)
        klcs = klcs0;
    }	// for (ijcs)
    *petmp_non_zero_eri = nzeri;
    return 0;
}
