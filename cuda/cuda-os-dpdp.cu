// dpdp
#include "cuda-twoint-core-os.h"

__device__ void gpu_hrr_clear_dpdp( double *eh ) {
    int i;
    // (DS|DS)
#pragma unroll
    for ( i=0; i<(0+36); i++ ) eh[i] = 0.e0;
    // (DS|FS)
#pragma unroll
    for ( i=36; i<(36+60); i++ ) eh[i] = 0.e0;
    // (FS|DS)
#pragma unroll
    for ( i=96; i<(96+60); i++ ) eh[i] = 0.e0;
    // (FS|FS)
#pragma unroll
    for ( i=156; i<(156+100); i++ ) eh[i] = 0.e0;
}

__device__ void gpu_hrr_calc_dpdp(
        const double BA[3], const double DC[3], double *eh ) {
    // (DP,DS)
    eh[ 256] = eh[  96] - BA[0]*eh[   0];
    eh[ 257] = eh[  97] - BA[0]*eh[   1];
    eh[ 258] = eh[  98] - BA[0]*eh[   2];
    eh[ 259] = eh[  99] - BA[0]*eh[   3];
    eh[ 260] = eh[ 100] - BA[0]*eh[   4];
    eh[ 261] = eh[ 101] - BA[0]*eh[   5];
    eh[ 262] = eh[ 114] - BA[1]*eh[   0];
    eh[ 263] = eh[ 115] - BA[1]*eh[   1];
    eh[ 264] = eh[ 116] - BA[1]*eh[   2];
    eh[ 265] = eh[ 117] - BA[1]*eh[   3];
    eh[ 266] = eh[ 118] - BA[1]*eh[   4];
    eh[ 267] = eh[ 119] - BA[1]*eh[   5];
    eh[ 268] = eh[ 120] - BA[2]*eh[   0];
    eh[ 269] = eh[ 121] - BA[2]*eh[   1];
    eh[ 270] = eh[ 122] - BA[2]*eh[   2];
    eh[ 271] = eh[ 123] - BA[2]*eh[   3];
    eh[ 272] = eh[ 124] - BA[2]*eh[   4];
    eh[ 273] = eh[ 125] - BA[2]*eh[   5];
    eh[ 274] = eh[ 126] - BA[0]*eh[   6];
    eh[ 275] = eh[ 127] - BA[0]*eh[   7];
    eh[ 276] = eh[ 128] - BA[0]*eh[   8];
    eh[ 277] = eh[ 129] - BA[0]*eh[   9];
    eh[ 278] = eh[ 130] - BA[0]*eh[  10];
    eh[ 279] = eh[ 131] - BA[0]*eh[  11];
    eh[ 280] = eh[ 102] - BA[1]*eh[   6];
    eh[ 281] = eh[ 103] - BA[1]*eh[   7];
    eh[ 282] = eh[ 104] - BA[1]*eh[   8];
    eh[ 283] = eh[ 105] - BA[1]*eh[   9];
    eh[ 284] = eh[ 106] - BA[1]*eh[  10];
    eh[ 285] = eh[ 107] - BA[1]*eh[  11];
    eh[ 286] = eh[ 144] - BA[2]*eh[   6];
    eh[ 287] = eh[ 145] - BA[2]*eh[   7];
    eh[ 288] = eh[ 146] - BA[2]*eh[   8];
    eh[ 289] = eh[ 147] - BA[2]*eh[   9];
    eh[ 290] = eh[ 148] - BA[2]*eh[  10];
    eh[ 291] = eh[ 149] - BA[2]*eh[  11];
    eh[ 292] = eh[ 132] - BA[0]*eh[  12];
    eh[ 293] = eh[ 133] - BA[0]*eh[  13];
    eh[ 294] = eh[ 134] - BA[0]*eh[  14];
    eh[ 295] = eh[ 135] - BA[0]*eh[  15];
    eh[ 296] = eh[ 136] - BA[0]*eh[  16];
    eh[ 297] = eh[ 137] - BA[0]*eh[  17];
    eh[ 298] = eh[ 150] - BA[1]*eh[  12];
    eh[ 299] = eh[ 151] - BA[1]*eh[  13];
    eh[ 300] = eh[ 152] - BA[1]*eh[  14];
    eh[ 301] = eh[ 153] - BA[1]*eh[  15];
    eh[ 302] = eh[ 154] - BA[1]*eh[  16];
    eh[ 303] = eh[ 155] - BA[1]*eh[  17];
    eh[ 304] = eh[ 108] - BA[2]*eh[  12];
    eh[ 305] = eh[ 109] - BA[2]*eh[  13];
    eh[ 306] = eh[ 110] - BA[2]*eh[  14];
    eh[ 307] = eh[ 111] - BA[2]*eh[  15];
    eh[ 308] = eh[ 112] - BA[2]*eh[  16];
    eh[ 309] = eh[ 113] - BA[2]*eh[  17];
    eh[ 310] = eh[ 114] - BA[0]*eh[  18];
    eh[ 311] = eh[ 115] - BA[0]*eh[  19];
    eh[ 312] = eh[ 116] - BA[0]*eh[  20];
    eh[ 313] = eh[ 117] - BA[0]*eh[  21];
    eh[ 314] = eh[ 118] - BA[0]*eh[  22];
    eh[ 315] = eh[ 119] - BA[0]*eh[  23];
    eh[ 316] = eh[ 126] - BA[1]*eh[  18];
    eh[ 317] = eh[ 127] - BA[1]*eh[  19];
    eh[ 318] = eh[ 128] - BA[1]*eh[  20];
    eh[ 319] = eh[ 129] - BA[1]*eh[  21];
    eh[ 320] = eh[ 130] - BA[1]*eh[  22];
    eh[ 321] = eh[ 131] - BA[1]*eh[  23];
    eh[ 322] = eh[ 138] - BA[2]*eh[  18];
    eh[ 323] = eh[ 139] - BA[2]*eh[  19];
    eh[ 324] = eh[ 140] - BA[2]*eh[  20];
    eh[ 325] = eh[ 141] - BA[2]*eh[  21];
    eh[ 326] = eh[ 142] - BA[2]*eh[  22];
    eh[ 327] = eh[ 143] - BA[2]*eh[  23];
    eh[ 328] = eh[ 120] - BA[0]*eh[  24];
    eh[ 329] = eh[ 121] - BA[0]*eh[  25];
    eh[ 330] = eh[ 122] - BA[0]*eh[  26];
    eh[ 331] = eh[ 123] - BA[0]*eh[  27];
    eh[ 332] = eh[ 124] - BA[0]*eh[  28];
    eh[ 333] = eh[ 125] - BA[0]*eh[  29];
    eh[ 334] = eh[ 138] - BA[1]*eh[  24];
    eh[ 335] = eh[ 139] - BA[1]*eh[  25];
    eh[ 336] = eh[ 140] - BA[1]*eh[  26];
    eh[ 337] = eh[ 141] - BA[1]*eh[  27];
    eh[ 338] = eh[ 142] - BA[1]*eh[  28];
    eh[ 339] = eh[ 143] - BA[1]*eh[  29];
    eh[ 340] = eh[ 132] - BA[2]*eh[  24];
    eh[ 341] = eh[ 133] - BA[2]*eh[  25];
    eh[ 342] = eh[ 134] - BA[2]*eh[  26];
    eh[ 343] = eh[ 135] - BA[2]*eh[  27];
    eh[ 344] = eh[ 136] - BA[2]*eh[  28];
    eh[ 345] = eh[ 137] - BA[2]*eh[  29];
    eh[ 346] = eh[ 138] - BA[0]*eh[  30];
    eh[ 347] = eh[ 139] - BA[0]*eh[  31];
    eh[ 348] = eh[ 140] - BA[0]*eh[  32];
    eh[ 349] = eh[ 141] - BA[0]*eh[  33];
    eh[ 350] = eh[ 142] - BA[0]*eh[  34];
    eh[ 351] = eh[ 143] - BA[0]*eh[  35];
    eh[ 352] = eh[ 144] - BA[1]*eh[  30];
    eh[ 353] = eh[ 145] - BA[1]*eh[  31];
    eh[ 354] = eh[ 146] - BA[1]*eh[  32];
    eh[ 355] = eh[ 147] - BA[1]*eh[  33];
    eh[ 356] = eh[ 148] - BA[1]*eh[  34];
    eh[ 357] = eh[ 149] - BA[1]*eh[  35];
    eh[ 358] = eh[ 150] - BA[2]*eh[  30];
    eh[ 359] = eh[ 151] - BA[2]*eh[  31];
    eh[ 360] = eh[ 152] - BA[2]*eh[  32];
    eh[ 361] = eh[ 153] - BA[2]*eh[  33];
    eh[ 362] = eh[ 154] - BA[2]*eh[  34];
    eh[ 363] = eh[ 155] - BA[2]*eh[  35];
    // (DP,FS)
    eh[ 364] = eh[ 156] - BA[0]*eh[  36];
    eh[ 365] = eh[ 157] - BA[0]*eh[  37];
    eh[ 366] = eh[ 158] - BA[0]*eh[  38];
    eh[ 367] = eh[ 159] - BA[0]*eh[  39];
    eh[ 368] = eh[ 160] - BA[0]*eh[  40];
    eh[ 369] = eh[ 161] - BA[0]*eh[  41];
    eh[ 370] = eh[ 162] - BA[0]*eh[  42];
    eh[ 371] = eh[ 163] - BA[0]*eh[  43];
    eh[ 372] = eh[ 164] - BA[0]*eh[  44];
    eh[ 373] = eh[ 165] - BA[0]*eh[  45];
    eh[ 374] = eh[ 186] - BA[1]*eh[  36];
    eh[ 375] = eh[ 187] - BA[1]*eh[  37];
    eh[ 376] = eh[ 188] - BA[1]*eh[  38];
    eh[ 377] = eh[ 189] - BA[1]*eh[  39];
    eh[ 378] = eh[ 190] - BA[1]*eh[  40];
    eh[ 379] = eh[ 191] - BA[1]*eh[  41];
    eh[ 380] = eh[ 192] - BA[1]*eh[  42];
    eh[ 381] = eh[ 193] - BA[1]*eh[  43];
    eh[ 382] = eh[ 194] - BA[1]*eh[  44];
    eh[ 383] = eh[ 195] - BA[1]*eh[  45];
    eh[ 384] = eh[ 196] - BA[2]*eh[  36];
    eh[ 385] = eh[ 197] - BA[2]*eh[  37];
    eh[ 386] = eh[ 198] - BA[2]*eh[  38];
    eh[ 387] = eh[ 199] - BA[2]*eh[  39];
    eh[ 388] = eh[ 200] - BA[2]*eh[  40];
    eh[ 389] = eh[ 201] - BA[2]*eh[  41];
    eh[ 390] = eh[ 202] - BA[2]*eh[  42];
    eh[ 391] = eh[ 203] - BA[2]*eh[  43];
    eh[ 392] = eh[ 204] - BA[2]*eh[  44];
    eh[ 393] = eh[ 205] - BA[2]*eh[  45];
    eh[ 394] = eh[ 206] - BA[0]*eh[  46];
    eh[ 395] = eh[ 207] - BA[0]*eh[  47];
    eh[ 396] = eh[ 208] - BA[0]*eh[  48];
    eh[ 397] = eh[ 209] - BA[0]*eh[  49];
    eh[ 398] = eh[ 210] - BA[0]*eh[  50];
    eh[ 399] = eh[ 211] - BA[0]*eh[  51];
    eh[ 400] = eh[ 212] - BA[0]*eh[  52];
    eh[ 401] = eh[ 213] - BA[0]*eh[  53];
    eh[ 402] = eh[ 214] - BA[0]*eh[  54];
    eh[ 403] = eh[ 215] - BA[0]*eh[  55];
    eh[ 404] = eh[ 166] - BA[1]*eh[  46];
    eh[ 405] = eh[ 167] - BA[1]*eh[  47];
    eh[ 406] = eh[ 168] - BA[1]*eh[  48];
    eh[ 407] = eh[ 169] - BA[1]*eh[  49];
    eh[ 408] = eh[ 170] - BA[1]*eh[  50];
    eh[ 409] = eh[ 171] - BA[1]*eh[  51];
    eh[ 410] = eh[ 172] - BA[1]*eh[  52];
    eh[ 411] = eh[ 173] - BA[1]*eh[  53];
    eh[ 412] = eh[ 174] - BA[1]*eh[  54];
    eh[ 413] = eh[ 175] - BA[1]*eh[  55];
    eh[ 414] = eh[ 236] - BA[2]*eh[  46];
    eh[ 415] = eh[ 237] - BA[2]*eh[  47];
    eh[ 416] = eh[ 238] - BA[2]*eh[  48];
    eh[ 417] = eh[ 239] - BA[2]*eh[  49];
    eh[ 418] = eh[ 240] - BA[2]*eh[  50];
    eh[ 419] = eh[ 241] - BA[2]*eh[  51];
    eh[ 420] = eh[ 242] - BA[2]*eh[  52];
    eh[ 421] = eh[ 243] - BA[2]*eh[  53];
    eh[ 422] = eh[ 244] - BA[2]*eh[  54];
    eh[ 423] = eh[ 245] - BA[2]*eh[  55];
    eh[ 424] = eh[ 216] - BA[0]*eh[  56];
    eh[ 425] = eh[ 217] - BA[0]*eh[  57];
    eh[ 426] = eh[ 218] - BA[0]*eh[  58];
    eh[ 427] = eh[ 219] - BA[0]*eh[  59];
    eh[ 428] = eh[ 220] - BA[0]*eh[  60];
    eh[ 429] = eh[ 221] - BA[0]*eh[  61];
    eh[ 430] = eh[ 222] - BA[0]*eh[  62];
    eh[ 431] = eh[ 223] - BA[0]*eh[  63];
    eh[ 432] = eh[ 224] - BA[0]*eh[  64];
    eh[ 433] = eh[ 225] - BA[0]*eh[  65];
    eh[ 434] = eh[ 246] - BA[1]*eh[  56];
    eh[ 435] = eh[ 247] - BA[1]*eh[  57];
    eh[ 436] = eh[ 248] - BA[1]*eh[  58];
    eh[ 437] = eh[ 249] - BA[1]*eh[  59];
    eh[ 438] = eh[ 250] - BA[1]*eh[  60];
    eh[ 439] = eh[ 251] - BA[1]*eh[  61];
    eh[ 440] = eh[ 252] - BA[1]*eh[  62];
    eh[ 441] = eh[ 253] - BA[1]*eh[  63];
    eh[ 442] = eh[ 254] - BA[1]*eh[  64];
    eh[ 443] = eh[ 255] - BA[1]*eh[  65];
    eh[ 444] = eh[ 176] - BA[2]*eh[  56];
    eh[ 445] = eh[ 177] - BA[2]*eh[  57];
    eh[ 446] = eh[ 178] - BA[2]*eh[  58];
    eh[ 447] = eh[ 179] - BA[2]*eh[  59];
    eh[ 448] = eh[ 180] - BA[2]*eh[  60];
    eh[ 449] = eh[ 181] - BA[2]*eh[  61];
    eh[ 450] = eh[ 182] - BA[2]*eh[  62];
    eh[ 451] = eh[ 183] - BA[2]*eh[  63];
    eh[ 452] = eh[ 184] - BA[2]*eh[  64];
    eh[ 453] = eh[ 185] - BA[2]*eh[  65];
    eh[ 454] = eh[ 186] - BA[0]*eh[  66];
    eh[ 455] = eh[ 187] - BA[0]*eh[  67];
    eh[ 456] = eh[ 188] - BA[0]*eh[  68];
    eh[ 457] = eh[ 189] - BA[0]*eh[  69];
    eh[ 458] = eh[ 190] - BA[0]*eh[  70];
    eh[ 459] = eh[ 191] - BA[0]*eh[  71];
    eh[ 460] = eh[ 192] - BA[0]*eh[  72];
    eh[ 461] = eh[ 193] - BA[0]*eh[  73];
    eh[ 462] = eh[ 194] - BA[0]*eh[  74];
    eh[ 463] = eh[ 195] - BA[0]*eh[  75];
    eh[ 464] = eh[ 206] - BA[1]*eh[  66];
    eh[ 465] = eh[ 207] - BA[1]*eh[  67];
    eh[ 466] = eh[ 208] - BA[1]*eh[  68];
    eh[ 467] = eh[ 209] - BA[1]*eh[  69];
    eh[ 468] = eh[ 210] - BA[1]*eh[  70];
    eh[ 469] = eh[ 211] - BA[1]*eh[  71];
    eh[ 470] = eh[ 212] - BA[1]*eh[  72];
    eh[ 471] = eh[ 213] - BA[1]*eh[  73];
    eh[ 472] = eh[ 214] - BA[1]*eh[  74];
    eh[ 473] = eh[ 215] - BA[1]*eh[  75];
    eh[ 474] = eh[ 226] - BA[2]*eh[  66];
    eh[ 475] = eh[ 227] - BA[2]*eh[  67];
    eh[ 476] = eh[ 228] - BA[2]*eh[  68];
    eh[ 477] = eh[ 229] - BA[2]*eh[  69];
    eh[ 478] = eh[ 230] - BA[2]*eh[  70];
    eh[ 479] = eh[ 231] - BA[2]*eh[  71];
    eh[ 480] = eh[ 232] - BA[2]*eh[  72];
    eh[ 481] = eh[ 233] - BA[2]*eh[  73];
    eh[ 482] = eh[ 234] - BA[2]*eh[  74];
    eh[ 483] = eh[ 235] - BA[2]*eh[  75];
    eh[ 484] = eh[ 196] - BA[0]*eh[  76];
    eh[ 485] = eh[ 197] - BA[0]*eh[  77];
    eh[ 486] = eh[ 198] - BA[0]*eh[  78];
    eh[ 487] = eh[ 199] - BA[0]*eh[  79];
    eh[ 488] = eh[ 200] - BA[0]*eh[  80];
    eh[ 489] = eh[ 201] - BA[0]*eh[  81];
    eh[ 490] = eh[ 202] - BA[0]*eh[  82];
    eh[ 491] = eh[ 203] - BA[0]*eh[  83];
    eh[ 492] = eh[ 204] - BA[0]*eh[  84];
    eh[ 493] = eh[ 205] - BA[0]*eh[  85];
    eh[ 494] = eh[ 226] - BA[1]*eh[  76];
    eh[ 495] = eh[ 227] - BA[1]*eh[  77];
    eh[ 496] = eh[ 228] - BA[1]*eh[  78];
    eh[ 497] = eh[ 229] - BA[1]*eh[  79];
    eh[ 498] = eh[ 230] - BA[1]*eh[  80];
    eh[ 499] = eh[ 231] - BA[1]*eh[  81];
    eh[ 500] = eh[ 232] - BA[1]*eh[  82];
    eh[ 501] = eh[ 233] - BA[1]*eh[  83];
    eh[ 502] = eh[ 234] - BA[1]*eh[  84];
    eh[ 503] = eh[ 235] - BA[1]*eh[  85];
    eh[ 504] = eh[ 216] - BA[2]*eh[  76];
    eh[ 505] = eh[ 217] - BA[2]*eh[  77];
    eh[ 506] = eh[ 218] - BA[2]*eh[  78];
    eh[ 507] = eh[ 219] - BA[2]*eh[  79];
    eh[ 508] = eh[ 220] - BA[2]*eh[  80];
    eh[ 509] = eh[ 221] - BA[2]*eh[  81];
    eh[ 510] = eh[ 222] - BA[2]*eh[  82];
    eh[ 511] = eh[ 223] - BA[2]*eh[  83];
    eh[ 512] = eh[ 224] - BA[2]*eh[  84];
    eh[ 513] = eh[ 225] - BA[2]*eh[  85];
    eh[ 514] = eh[ 226] - BA[0]*eh[  86];
    eh[ 515] = eh[ 227] - BA[0]*eh[  87];
    eh[ 516] = eh[ 228] - BA[0]*eh[  88];
    eh[ 517] = eh[ 229] - BA[0]*eh[  89];
    eh[ 518] = eh[ 230] - BA[0]*eh[  90];
    eh[ 519] = eh[ 231] - BA[0]*eh[  91];
    eh[ 520] = eh[ 232] - BA[0]*eh[  92];
    eh[ 521] = eh[ 233] - BA[0]*eh[  93];
    eh[ 522] = eh[ 234] - BA[0]*eh[  94];
    eh[ 523] = eh[ 235] - BA[0]*eh[  95];
    eh[ 524] = eh[ 236] - BA[1]*eh[  86];
    eh[ 525] = eh[ 237] - BA[1]*eh[  87];
    eh[ 526] = eh[ 238] - BA[1]*eh[  88];
    eh[ 527] = eh[ 239] - BA[1]*eh[  89];
    eh[ 528] = eh[ 240] - BA[1]*eh[  90];
    eh[ 529] = eh[ 241] - BA[1]*eh[  91];
    eh[ 530] = eh[ 242] - BA[1]*eh[  92];
    eh[ 531] = eh[ 243] - BA[1]*eh[  93];
    eh[ 532] = eh[ 244] - BA[1]*eh[  94];
    eh[ 533] = eh[ 245] - BA[1]*eh[  95];
    eh[ 534] = eh[ 246] - BA[2]*eh[  86];
    eh[ 535] = eh[ 247] - BA[2]*eh[  87];
    eh[ 536] = eh[ 248] - BA[2]*eh[  88];
    eh[ 537] = eh[ 249] - BA[2]*eh[  89];
    eh[ 538] = eh[ 250] - BA[2]*eh[  90];
    eh[ 539] = eh[ 251] - BA[2]*eh[  91];
    eh[ 540] = eh[ 252] - BA[2]*eh[  92];
    eh[ 541] = eh[ 253] - BA[2]*eh[  93];
    eh[ 542] = eh[ 254] - BA[2]*eh[  94];
    eh[ 543] = eh[ 255] - BA[2]*eh[  95];
    // HRR for (XX|XX)-type integral (center CD)
    // (DP,DP)
    eh[ 544] = eh[ 364] - DC[0]*eh[ 256];
    eh[ 545] = eh[ 367] - DC[1]*eh[ 256];
    eh[ 546] = eh[ 368] - DC[2]*eh[ 256];
    eh[ 547] = eh[ 369] - DC[0]*eh[ 257];
    eh[ 548] = eh[ 365] - DC[1]*eh[ 257];
    eh[ 549] = eh[ 372] - DC[2]*eh[ 257];
    eh[ 550] = eh[ 370] - DC[0]*eh[ 258];
    eh[ 551] = eh[ 373] - DC[1]*eh[ 258];
    eh[ 552] = eh[ 366] - DC[2]*eh[ 258];
    eh[ 553] = eh[ 367] - DC[0]*eh[ 259];
    eh[ 554] = eh[ 369] - DC[1]*eh[ 259];
    eh[ 555] = eh[ 371] - DC[2]*eh[ 259];
    eh[ 556] = eh[ 368] - DC[0]*eh[ 260];
    eh[ 557] = eh[ 371] - DC[1]*eh[ 260];
    eh[ 558] = eh[ 370] - DC[2]*eh[ 260];
    eh[ 559] = eh[ 371] - DC[0]*eh[ 261];
    eh[ 560] = eh[ 372] - DC[1]*eh[ 261];
    eh[ 561] = eh[ 373] - DC[2]*eh[ 261];
    eh[ 562] = eh[ 374] - DC[0]*eh[ 262];
    eh[ 563] = eh[ 377] - DC[1]*eh[ 262];
    eh[ 564] = eh[ 378] - DC[2]*eh[ 262];
    eh[ 565] = eh[ 379] - DC[0]*eh[ 263];
    eh[ 566] = eh[ 375] - DC[1]*eh[ 263];
    eh[ 567] = eh[ 382] - DC[2]*eh[ 263];
    eh[ 568] = eh[ 380] - DC[0]*eh[ 264];
    eh[ 569] = eh[ 383] - DC[1]*eh[ 264];
    eh[ 570] = eh[ 376] - DC[2]*eh[ 264];
    eh[ 571] = eh[ 377] - DC[0]*eh[ 265];
    eh[ 572] = eh[ 379] - DC[1]*eh[ 265];
    eh[ 573] = eh[ 381] - DC[2]*eh[ 265];
    eh[ 574] = eh[ 378] - DC[0]*eh[ 266];
    eh[ 575] = eh[ 381] - DC[1]*eh[ 266];
    eh[ 576] = eh[ 380] - DC[2]*eh[ 266];
    eh[ 577] = eh[ 381] - DC[0]*eh[ 267];
    eh[ 578] = eh[ 382] - DC[1]*eh[ 267];
    eh[ 579] = eh[ 383] - DC[2]*eh[ 267];
    eh[ 580] = eh[ 384] - DC[0]*eh[ 268];
    eh[ 581] = eh[ 387] - DC[1]*eh[ 268];
    eh[ 582] = eh[ 388] - DC[2]*eh[ 268];
    eh[ 583] = eh[ 389] - DC[0]*eh[ 269];
    eh[ 584] = eh[ 385] - DC[1]*eh[ 269];
    eh[ 585] = eh[ 392] - DC[2]*eh[ 269];
    eh[ 586] = eh[ 390] - DC[0]*eh[ 270];
    eh[ 587] = eh[ 393] - DC[1]*eh[ 270];
    eh[ 588] = eh[ 386] - DC[2]*eh[ 270];
    eh[ 589] = eh[ 387] - DC[0]*eh[ 271];
    eh[ 590] = eh[ 389] - DC[1]*eh[ 271];
    eh[ 591] = eh[ 391] - DC[2]*eh[ 271];
    eh[ 592] = eh[ 388] - DC[0]*eh[ 272];
    eh[ 593] = eh[ 391] - DC[1]*eh[ 272];
    eh[ 594] = eh[ 390] - DC[2]*eh[ 272];
    eh[ 595] = eh[ 391] - DC[0]*eh[ 273];
    eh[ 596] = eh[ 392] - DC[1]*eh[ 273];
    eh[ 597] = eh[ 393] - DC[2]*eh[ 273];
    eh[ 598] = eh[ 394] - DC[0]*eh[ 274];
    eh[ 599] = eh[ 397] - DC[1]*eh[ 274];
    eh[ 600] = eh[ 398] - DC[2]*eh[ 274];
    eh[ 601] = eh[ 399] - DC[0]*eh[ 275];
    eh[ 602] = eh[ 395] - DC[1]*eh[ 275];
    eh[ 603] = eh[ 402] - DC[2]*eh[ 275];
    eh[ 604] = eh[ 400] - DC[0]*eh[ 276];
    eh[ 605] = eh[ 403] - DC[1]*eh[ 276];
    eh[ 606] = eh[ 396] - DC[2]*eh[ 276];
    eh[ 607] = eh[ 397] - DC[0]*eh[ 277];
    eh[ 608] = eh[ 399] - DC[1]*eh[ 277];
    eh[ 609] = eh[ 401] - DC[2]*eh[ 277];
    eh[ 610] = eh[ 398] - DC[0]*eh[ 278];
    eh[ 611] = eh[ 401] - DC[1]*eh[ 278];
    eh[ 612] = eh[ 400] - DC[2]*eh[ 278];
    eh[ 613] = eh[ 401] - DC[0]*eh[ 279];
    eh[ 614] = eh[ 402] - DC[1]*eh[ 279];
    eh[ 615] = eh[ 403] - DC[2]*eh[ 279];
    eh[ 616] = eh[ 404] - DC[0]*eh[ 280];
    eh[ 617] = eh[ 407] - DC[1]*eh[ 280];
    eh[ 618] = eh[ 408] - DC[2]*eh[ 280];
    eh[ 619] = eh[ 409] - DC[0]*eh[ 281];
    eh[ 620] = eh[ 405] - DC[1]*eh[ 281];
    eh[ 621] = eh[ 412] - DC[2]*eh[ 281];
    eh[ 622] = eh[ 410] - DC[0]*eh[ 282];
    eh[ 623] = eh[ 413] - DC[1]*eh[ 282];
    eh[ 624] = eh[ 406] - DC[2]*eh[ 282];
    eh[ 625] = eh[ 407] - DC[0]*eh[ 283];
    eh[ 626] = eh[ 409] - DC[1]*eh[ 283];
    eh[ 627] = eh[ 411] - DC[2]*eh[ 283];
    eh[ 628] = eh[ 408] - DC[0]*eh[ 284];
    eh[ 629] = eh[ 411] - DC[1]*eh[ 284];
    eh[ 630] = eh[ 410] - DC[2]*eh[ 284];
    eh[ 631] = eh[ 411] - DC[0]*eh[ 285];
    eh[ 632] = eh[ 412] - DC[1]*eh[ 285];
    eh[ 633] = eh[ 413] - DC[2]*eh[ 285];
    eh[ 634] = eh[ 414] - DC[0]*eh[ 286];
    eh[ 635] = eh[ 417] - DC[1]*eh[ 286];
    eh[ 636] = eh[ 418] - DC[2]*eh[ 286];
    eh[ 637] = eh[ 419] - DC[0]*eh[ 287];
    eh[ 638] = eh[ 415] - DC[1]*eh[ 287];
    eh[ 639] = eh[ 422] - DC[2]*eh[ 287];
    eh[ 640] = eh[ 420] - DC[0]*eh[ 288];
    eh[ 641] = eh[ 423] - DC[1]*eh[ 288];
    eh[ 642] = eh[ 416] - DC[2]*eh[ 288];
    eh[ 643] = eh[ 417] - DC[0]*eh[ 289];
    eh[ 644] = eh[ 419] - DC[1]*eh[ 289];
    eh[ 645] = eh[ 421] - DC[2]*eh[ 289];
    eh[ 646] = eh[ 418] - DC[0]*eh[ 290];
    eh[ 647] = eh[ 421] - DC[1]*eh[ 290];
    eh[ 648] = eh[ 420] - DC[2]*eh[ 290];
    eh[ 649] = eh[ 421] - DC[0]*eh[ 291];
    eh[ 650] = eh[ 422] - DC[1]*eh[ 291];
    eh[ 651] = eh[ 423] - DC[2]*eh[ 291];
    eh[ 652] = eh[ 424] - DC[0]*eh[ 292];
    eh[ 653] = eh[ 427] - DC[1]*eh[ 292];
    eh[ 654] = eh[ 428] - DC[2]*eh[ 292];
    eh[ 655] = eh[ 429] - DC[0]*eh[ 293];
    eh[ 656] = eh[ 425] - DC[1]*eh[ 293];
    eh[ 657] = eh[ 432] - DC[2]*eh[ 293];
    eh[ 658] = eh[ 430] - DC[0]*eh[ 294];
    eh[ 659] = eh[ 433] - DC[1]*eh[ 294];
    eh[ 660] = eh[ 426] - DC[2]*eh[ 294];
    eh[ 661] = eh[ 427] - DC[0]*eh[ 295];
    eh[ 662] = eh[ 429] - DC[1]*eh[ 295];
    eh[ 663] = eh[ 431] - DC[2]*eh[ 295];
    eh[ 664] = eh[ 428] - DC[0]*eh[ 296];
    eh[ 665] = eh[ 431] - DC[1]*eh[ 296];
    eh[ 666] = eh[ 430] - DC[2]*eh[ 296];
    eh[ 667] = eh[ 431] - DC[0]*eh[ 297];
    eh[ 668] = eh[ 432] - DC[1]*eh[ 297];
    eh[ 669] = eh[ 433] - DC[2]*eh[ 297];
    eh[ 670] = eh[ 434] - DC[0]*eh[ 298];
    eh[ 671] = eh[ 437] - DC[1]*eh[ 298];
    eh[ 672] = eh[ 438] - DC[2]*eh[ 298];
    eh[ 673] = eh[ 439] - DC[0]*eh[ 299];
    eh[ 674] = eh[ 435] - DC[1]*eh[ 299];
    eh[ 675] = eh[ 442] - DC[2]*eh[ 299];
    eh[ 676] = eh[ 440] - DC[0]*eh[ 300];
    eh[ 677] = eh[ 443] - DC[1]*eh[ 300];
    eh[ 678] = eh[ 436] - DC[2]*eh[ 300];
    eh[ 679] = eh[ 437] - DC[0]*eh[ 301];
    eh[ 680] = eh[ 439] - DC[1]*eh[ 301];
    eh[ 681] = eh[ 441] - DC[2]*eh[ 301];
    eh[ 682] = eh[ 438] - DC[0]*eh[ 302];
    eh[ 683] = eh[ 441] - DC[1]*eh[ 302];
    eh[ 684] = eh[ 440] - DC[2]*eh[ 302];
    eh[ 685] = eh[ 441] - DC[0]*eh[ 303];
    eh[ 686] = eh[ 442] - DC[1]*eh[ 303];
    eh[ 687] = eh[ 443] - DC[2]*eh[ 303];
    eh[ 688] = eh[ 444] - DC[0]*eh[ 304];
    eh[ 689] = eh[ 447] - DC[1]*eh[ 304];
    eh[ 690] = eh[ 448] - DC[2]*eh[ 304];
    eh[ 691] = eh[ 449] - DC[0]*eh[ 305];
    eh[ 692] = eh[ 445] - DC[1]*eh[ 305];
    eh[ 693] = eh[ 452] - DC[2]*eh[ 305];
    eh[ 694] = eh[ 450] - DC[0]*eh[ 306];
    eh[ 695] = eh[ 453] - DC[1]*eh[ 306];
    eh[ 696] = eh[ 446] - DC[2]*eh[ 306];
    eh[ 697] = eh[ 447] - DC[0]*eh[ 307];
    eh[ 698] = eh[ 449] - DC[1]*eh[ 307];
    eh[ 699] = eh[ 451] - DC[2]*eh[ 307];
    eh[ 700] = eh[ 448] - DC[0]*eh[ 308];
    eh[ 701] = eh[ 451] - DC[1]*eh[ 308];
    eh[ 702] = eh[ 450] - DC[2]*eh[ 308];
    eh[ 703] = eh[ 451] - DC[0]*eh[ 309];
    eh[ 704] = eh[ 452] - DC[1]*eh[ 309];
    eh[ 705] = eh[ 453] - DC[2]*eh[ 309];
    eh[ 706] = eh[ 454] - DC[0]*eh[ 310];
    eh[ 707] = eh[ 457] - DC[1]*eh[ 310];
    eh[ 708] = eh[ 458] - DC[2]*eh[ 310];
    eh[ 709] = eh[ 459] - DC[0]*eh[ 311];
    eh[ 710] = eh[ 455] - DC[1]*eh[ 311];
    eh[ 711] = eh[ 462] - DC[2]*eh[ 311];
    eh[ 712] = eh[ 460] - DC[0]*eh[ 312];
    eh[ 713] = eh[ 463] - DC[1]*eh[ 312];
    eh[ 714] = eh[ 456] - DC[2]*eh[ 312];
    eh[ 715] = eh[ 457] - DC[0]*eh[ 313];
    eh[ 716] = eh[ 459] - DC[1]*eh[ 313];
    eh[ 717] = eh[ 461] - DC[2]*eh[ 313];
    eh[ 718] = eh[ 458] - DC[0]*eh[ 314];
    eh[ 719] = eh[ 461] - DC[1]*eh[ 314];
    eh[ 720] = eh[ 460] - DC[2]*eh[ 314];
    eh[ 721] = eh[ 461] - DC[0]*eh[ 315];
    eh[ 722] = eh[ 462] - DC[1]*eh[ 315];
    eh[ 723] = eh[ 463] - DC[2]*eh[ 315];
    eh[ 724] = eh[ 464] - DC[0]*eh[ 316];
    eh[ 725] = eh[ 467] - DC[1]*eh[ 316];
    eh[ 726] = eh[ 468] - DC[2]*eh[ 316];
    eh[ 727] = eh[ 469] - DC[0]*eh[ 317];
    eh[ 728] = eh[ 465] - DC[1]*eh[ 317];
    eh[ 729] = eh[ 472] - DC[2]*eh[ 317];
    eh[ 730] = eh[ 470] - DC[0]*eh[ 318];
    eh[ 731] = eh[ 473] - DC[1]*eh[ 318];
    eh[ 732] = eh[ 466] - DC[2]*eh[ 318];
    eh[ 733] = eh[ 467] - DC[0]*eh[ 319];
    eh[ 734] = eh[ 469] - DC[1]*eh[ 319];
    eh[ 735] = eh[ 471] - DC[2]*eh[ 319];
    eh[ 736] = eh[ 468] - DC[0]*eh[ 320];
    eh[ 737] = eh[ 471] - DC[1]*eh[ 320];
    eh[ 738] = eh[ 470] - DC[2]*eh[ 320];
    eh[ 739] = eh[ 471] - DC[0]*eh[ 321];
    eh[ 740] = eh[ 472] - DC[1]*eh[ 321];
    eh[ 741] = eh[ 473] - DC[2]*eh[ 321];
    eh[ 742] = eh[ 474] - DC[0]*eh[ 322];
    eh[ 743] = eh[ 477] - DC[1]*eh[ 322];
    eh[ 744] = eh[ 478] - DC[2]*eh[ 322];
    eh[ 745] = eh[ 479] - DC[0]*eh[ 323];
    eh[ 746] = eh[ 475] - DC[1]*eh[ 323];
    eh[ 747] = eh[ 482] - DC[2]*eh[ 323];
    eh[ 748] = eh[ 480] - DC[0]*eh[ 324];
    eh[ 749] = eh[ 483] - DC[1]*eh[ 324];
    eh[ 750] = eh[ 476] - DC[2]*eh[ 324];
    eh[ 751] = eh[ 477] - DC[0]*eh[ 325];
    eh[ 752] = eh[ 479] - DC[1]*eh[ 325];
    eh[ 753] = eh[ 481] - DC[2]*eh[ 325];
    eh[ 754] = eh[ 478] - DC[0]*eh[ 326];
    eh[ 755] = eh[ 481] - DC[1]*eh[ 326];
    eh[ 756] = eh[ 480] - DC[2]*eh[ 326];
    eh[ 757] = eh[ 481] - DC[0]*eh[ 327];
    eh[ 758] = eh[ 482] - DC[1]*eh[ 327];
    eh[ 759] = eh[ 483] - DC[2]*eh[ 327];
    eh[ 760] = eh[ 484] - DC[0]*eh[ 328];
    eh[ 761] = eh[ 487] - DC[1]*eh[ 328];
    eh[ 762] = eh[ 488] - DC[2]*eh[ 328];
    eh[ 763] = eh[ 489] - DC[0]*eh[ 329];
    eh[ 764] = eh[ 485] - DC[1]*eh[ 329];
    eh[ 765] = eh[ 492] - DC[2]*eh[ 329];
    eh[ 766] = eh[ 490] - DC[0]*eh[ 330];
    eh[ 767] = eh[ 493] - DC[1]*eh[ 330];
    eh[ 768] = eh[ 486] - DC[2]*eh[ 330];
    eh[ 769] = eh[ 487] - DC[0]*eh[ 331];
    eh[ 770] = eh[ 489] - DC[1]*eh[ 331];
    eh[ 771] = eh[ 491] - DC[2]*eh[ 331];
    eh[ 772] = eh[ 488] - DC[0]*eh[ 332];
    eh[ 773] = eh[ 491] - DC[1]*eh[ 332];
    eh[ 774] = eh[ 490] - DC[2]*eh[ 332];
    eh[ 775] = eh[ 491] - DC[0]*eh[ 333];
    eh[ 776] = eh[ 492] - DC[1]*eh[ 333];
    eh[ 777] = eh[ 493] - DC[2]*eh[ 333];
    eh[ 778] = eh[ 494] - DC[0]*eh[ 334];
    eh[ 779] = eh[ 497] - DC[1]*eh[ 334];
    eh[ 780] = eh[ 498] - DC[2]*eh[ 334];
    eh[ 781] = eh[ 499] - DC[0]*eh[ 335];
    eh[ 782] = eh[ 495] - DC[1]*eh[ 335];
    eh[ 783] = eh[ 502] - DC[2]*eh[ 335];
    eh[ 784] = eh[ 500] - DC[0]*eh[ 336];
    eh[ 785] = eh[ 503] - DC[1]*eh[ 336];
    eh[ 786] = eh[ 496] - DC[2]*eh[ 336];
    eh[ 787] = eh[ 497] - DC[0]*eh[ 337];
    eh[ 788] = eh[ 499] - DC[1]*eh[ 337];
    eh[ 789] = eh[ 501] - DC[2]*eh[ 337];
    eh[ 790] = eh[ 498] - DC[0]*eh[ 338];
    eh[ 791] = eh[ 501] - DC[1]*eh[ 338];
    eh[ 792] = eh[ 500] - DC[2]*eh[ 338];
    eh[ 793] = eh[ 501] - DC[0]*eh[ 339];
    eh[ 794] = eh[ 502] - DC[1]*eh[ 339];
    eh[ 795] = eh[ 503] - DC[2]*eh[ 339];
    eh[ 796] = eh[ 504] - DC[0]*eh[ 340];
    eh[ 797] = eh[ 507] - DC[1]*eh[ 340];
    eh[ 798] = eh[ 508] - DC[2]*eh[ 340];
    eh[ 799] = eh[ 509] - DC[0]*eh[ 341];
    eh[ 800] = eh[ 505] - DC[1]*eh[ 341];
    eh[ 801] = eh[ 512] - DC[2]*eh[ 341];
    eh[ 802] = eh[ 510] - DC[0]*eh[ 342];
    eh[ 803] = eh[ 513] - DC[1]*eh[ 342];
    eh[ 804] = eh[ 506] - DC[2]*eh[ 342];
    eh[ 805] = eh[ 507] - DC[0]*eh[ 343];
    eh[ 806] = eh[ 509] - DC[1]*eh[ 343];
    eh[ 807] = eh[ 511] - DC[2]*eh[ 343];
    eh[ 808] = eh[ 508] - DC[0]*eh[ 344];
    eh[ 809] = eh[ 511] - DC[1]*eh[ 344];
    eh[ 810] = eh[ 510] - DC[2]*eh[ 344];
    eh[ 811] = eh[ 511] - DC[0]*eh[ 345];
    eh[ 812] = eh[ 512] - DC[1]*eh[ 345];
    eh[ 813] = eh[ 513] - DC[2]*eh[ 345];
    eh[ 814] = eh[ 514] - DC[0]*eh[ 346];
    eh[ 815] = eh[ 517] - DC[1]*eh[ 346];
    eh[ 816] = eh[ 518] - DC[2]*eh[ 346];
    eh[ 817] = eh[ 519] - DC[0]*eh[ 347];
    eh[ 818] = eh[ 515] - DC[1]*eh[ 347];
    eh[ 819] = eh[ 522] - DC[2]*eh[ 347];
    eh[ 820] = eh[ 520] - DC[0]*eh[ 348];
    eh[ 821] = eh[ 523] - DC[1]*eh[ 348];
    eh[ 822] = eh[ 516] - DC[2]*eh[ 348];
    eh[ 823] = eh[ 517] - DC[0]*eh[ 349];
    eh[ 824] = eh[ 519] - DC[1]*eh[ 349];
    eh[ 825] = eh[ 521] - DC[2]*eh[ 349];
    eh[ 826] = eh[ 518] - DC[0]*eh[ 350];
    eh[ 827] = eh[ 521] - DC[1]*eh[ 350];
    eh[ 828] = eh[ 520] - DC[2]*eh[ 350];
    eh[ 829] = eh[ 521] - DC[0]*eh[ 351];
    eh[ 830] = eh[ 522] - DC[1]*eh[ 351];
    eh[ 831] = eh[ 523] - DC[2]*eh[ 351];
    eh[ 832] = eh[ 524] - DC[0]*eh[ 352];
    eh[ 833] = eh[ 527] - DC[1]*eh[ 352];
    eh[ 834] = eh[ 528] - DC[2]*eh[ 352];
    eh[ 835] = eh[ 529] - DC[0]*eh[ 353];
    eh[ 836] = eh[ 525] - DC[1]*eh[ 353];
    eh[ 837] = eh[ 532] - DC[2]*eh[ 353];
    eh[ 838] = eh[ 530] - DC[0]*eh[ 354];
    eh[ 839] = eh[ 533] - DC[1]*eh[ 354];
    eh[ 840] = eh[ 526] - DC[2]*eh[ 354];
    eh[ 841] = eh[ 527] - DC[0]*eh[ 355];
    eh[ 842] = eh[ 529] - DC[1]*eh[ 355];
    eh[ 843] = eh[ 531] - DC[2]*eh[ 355];
    eh[ 844] = eh[ 528] - DC[0]*eh[ 356];
    eh[ 845] = eh[ 531] - DC[1]*eh[ 356];
    eh[ 846] = eh[ 530] - DC[2]*eh[ 356];
    eh[ 847] = eh[ 531] - DC[0]*eh[ 357];
    eh[ 848] = eh[ 532] - DC[1]*eh[ 357];
    eh[ 849] = eh[ 533] - DC[2]*eh[ 357];
    eh[ 850] = eh[ 534] - DC[0]*eh[ 358];
    eh[ 851] = eh[ 537] - DC[1]*eh[ 358];
    eh[ 852] = eh[ 538] - DC[2]*eh[ 358];
    eh[ 853] = eh[ 539] - DC[0]*eh[ 359];
    eh[ 854] = eh[ 535] - DC[1]*eh[ 359];
    eh[ 855] = eh[ 542] - DC[2]*eh[ 359];
    eh[ 856] = eh[ 540] - DC[0]*eh[ 360];
    eh[ 857] = eh[ 543] - DC[1]*eh[ 360];
    eh[ 858] = eh[ 536] - DC[2]*eh[ 360];
    eh[ 859] = eh[ 537] - DC[0]*eh[ 361];
    eh[ 860] = eh[ 539] - DC[1]*eh[ 361];
    eh[ 861] = eh[ 541] - DC[2]*eh[ 361];
    eh[ 862] = eh[ 538] - DC[0]*eh[ 362];
    eh[ 863] = eh[ 541] - DC[1]*eh[ 362];
    eh[ 864] = eh[ 540] - DC[2]*eh[ 362];
    eh[ 865] = eh[ 541] - DC[0]*eh[ 363];
    eh[ 866] = eh[ 542] - DC[1]*eh[ 363];
    eh[ 867] = eh[ 543] - DC[2]*eh[ 363];
}

__device__ void gpu_hrr_coef_dpdp(
        const double *eh, double *DINT ) {
    int i, j, k, l, iao, jao, kao, lao, ix, iy;
    double coef_a, coef_ab, coef_abc;
    ix = 544;
    iy = 0;

#pragma unroll
    for ( i=0, iao=4; i<6; i++, iao++ ) {
        coef_a = LDG(DFACT[iao]);
#pragma unroll
        for ( j=0, jao=1; j<3; j++, jao++ ) {
            coef_ab = coef_a * LDG(DFACT[jao]);
#pragma unroll
            for ( k=0, kao=4; k<6; k++, kao++ ) {
                coef_abc = coef_ab * LDG(DFACT[kao]);
#pragma unroll
                for ( l=0, lao=1; l<3; l++, lao++ ) {
                    DINT[iy] = coef_abc * LDG(DFACT[lao]) * eh[ix];
                    iy++;
                    ix++;
                }
            }
        }
    }
}

__device__ void gpu_vrr_calc_dpdp(
        const double T, const double cssss,
        const double zeta2, const double eta2, const double ze2,
        const double rz, const double re,
        const double PA[3], const double WP[3],
        const double QC[3], const double WQ[3],
        double *ev ) {
    // (ss|ss) m=0,6
    //fmt( &ev[0], 6, T, cssss );
    //OFMO_FMT( &ev[0], 6, T, cssss );
#if   CUDA_FMT_M == 3
    gpu_fmt6_method3( T, cssss, ev );
#elif CUDA_FMT_M == 2
    gpu_fmt6_method2( T, cssss, ev );
#elif CUDA_FMT_M == 1
    gpu_fmt6_method1( T, cssss, ev );
#else
    gpu_fmt6( ev, T, cssss );
#endif
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
    // (ss|ps) m=[2,2]
    ev[95]=QC[0]*ev[2]+WQ[0]*ev[3];
    ev[96]=QC[1]*ev[2]+WQ[1]*ev[3];
    ev[97]=QC[2]*ev[2]+WQ[2]*ev[3];
    // (ps|ps) m=[1,2]
    ev[ 98]=QC[0]*ev[10]+WQ[0]*ev[13]+ze2*ev[2];
    ev[ 99]=QC[1]*ev[10]+WQ[1]*ev[13];
    ev[100]=QC[2]*ev[10]+WQ[2]*ev[13];
    ev[101]=QC[0]*ev[11]+WQ[0]*ev[14];
    ev[102]=QC[1]*ev[11]+WQ[1]*ev[14]+ze2*ev[2];
    ev[103]=QC[2]*ev[11]+WQ[2]*ev[14];
    ev[104]=QC[0]*ev[12]+WQ[0]*ev[15];
    ev[105]=QC[1]*ev[12]+WQ[1]*ev[15];
    ev[106]=QC[2]*ev[12]+WQ[2]*ev[15]+ze2*ev[2];
    ev[107]=QC[0]*ev[13]+WQ[0]*ev[16]+ze2*ev[3];
    ev[108]=QC[1]*ev[13]+WQ[1]*ev[16];
    ev[109]=QC[2]*ev[13]+WQ[2]*ev[16];
    ev[110]=QC[0]*ev[14]+WQ[0]*ev[17];
    ev[111]=QC[1]*ev[14]+WQ[1]*ev[17]+ze2*ev[3];
    ev[112]=QC[2]*ev[14]+WQ[2]*ev[17];
    ev[113]=QC[0]*ev[15]+WQ[0]*ev[18];
    ev[114]=QC[1]*ev[15]+WQ[1]*ev[18];
    ev[115]=QC[2]*ev[15]+WQ[2]*ev[18]+ze2*ev[3];
    // (ds|ps) m=[0,2]
    ev[116]=QC[0]*ev[25]+WQ[0]*ev[31]+2.e0*ze2*ev[10];
    ev[117]=QC[1]*ev[25]+WQ[1]*ev[31];
    ev[118]=QC[2]*ev[25]+WQ[2]*ev[31];
    ev[119]=QC[0]*ev[26]+WQ[0]*ev[32];
    ev[120]=QC[1]*ev[26]+WQ[1]*ev[32]+2.e0*ze2*ev[11];
    ev[121]=QC[2]*ev[26]+WQ[2]*ev[32];
    ev[122]=QC[0]*ev[27]+WQ[0]*ev[33];
    ev[123]=QC[1]*ev[27]+WQ[1]*ev[33];
    ev[124]=QC[2]*ev[27]+WQ[2]*ev[33]+2.e0*ze2*ev[12];
    ev[125]=QC[0]*ev[28]+WQ[0]*ev[34]+     ze2*ev[11];
    ev[126]=QC[1]*ev[28]+WQ[1]*ev[34]+     ze2*ev[10];
    ev[127]=QC[2]*ev[28]+WQ[2]*ev[34];
    ev[128]=QC[0]*ev[29]+WQ[0]*ev[35]+     ze2*ev[12];
    ev[129]=QC[1]*ev[29]+WQ[1]*ev[35];
    ev[130]=QC[2]*ev[29]+WQ[2]*ev[35]+     ze2*ev[10];
    ev[131]=QC[0]*ev[30]+WQ[0]*ev[36];
    ev[132]=QC[1]*ev[30]+WQ[1]*ev[36]+     ze2*ev[12];
    ev[133]=QC[2]*ev[30]+WQ[2]*ev[36]+     ze2*ev[11];
    ev[134]=QC[0]*ev[31]+WQ[0]*ev[37]+2.e0*ze2*ev[13];
    ev[135]=QC[1]*ev[31]+WQ[1]*ev[37];
    ev[136]=QC[2]*ev[31]+WQ[2]*ev[37];
    ev[137]=QC[0]*ev[32]+WQ[0]*ev[38];
    ev[138]=QC[1]*ev[32]+WQ[1]*ev[38]+2.e0*ze2*ev[14];
    ev[139]=QC[2]*ev[32]+WQ[2]*ev[38];
    ev[140]=QC[0]*ev[33]+WQ[0]*ev[39];
    ev[141]=QC[1]*ev[33]+WQ[1]*ev[39];
    ev[142]=QC[2]*ev[33]+WQ[2]*ev[39]+2.e0*ze2*ev[15];
    ev[143]=QC[0]*ev[34]+WQ[0]*ev[40]+     ze2*ev[14];
    ev[144]=QC[1]*ev[34]+WQ[1]*ev[40]+     ze2*ev[13];
    ev[145]=QC[2]*ev[34]+WQ[2]*ev[40];
    ev[146]=QC[0]*ev[35]+WQ[0]*ev[41]+     ze2*ev[15];
    ev[147]=QC[1]*ev[35]+WQ[1]*ev[41];
    ev[148]=QC[2]*ev[35]+WQ[2]*ev[41]+     ze2*ev[13];
    ev[149]=QC[0]*ev[36]+WQ[0]*ev[42];
    ev[150]=QC[1]*ev[36]+WQ[1]*ev[42]+     ze2*ev[15];
    ev[151]=QC[2]*ev[36]+WQ[2]*ev[42]+     ze2*ev[14];
    ev[152]=QC[0]*ev[37]+WQ[0]*ev[43]+2.e0*ze2*ev[16];
    ev[153]=QC[1]*ev[37]+WQ[1]*ev[43];
    ev[154]=QC[2]*ev[37]+WQ[2]*ev[43];
    ev[155]=QC[0]*ev[38]+WQ[0]*ev[44];
    ev[156]=QC[1]*ev[38]+WQ[1]*ev[44]+2.e0*ze2*ev[17];
    ev[157]=QC[2]*ev[38]+WQ[2]*ev[44];
    ev[158]=QC[0]*ev[39]+WQ[0]*ev[45];
    ev[159]=QC[1]*ev[39]+WQ[1]*ev[45];
    ev[160]=QC[2]*ev[39]+WQ[2]*ev[45]+2.e0*ze2*ev[18];
    ev[161]=QC[0]*ev[40]+WQ[0]*ev[46]+     ze2*ev[17];
    ev[162]=QC[1]*ev[40]+WQ[1]*ev[46]+     ze2*ev[16];
    ev[163]=QC[2]*ev[40]+WQ[2]*ev[46];
    ev[164]=QC[0]*ev[41]+WQ[0]*ev[47]+     ze2*ev[18];
    ev[165]=QC[1]*ev[41]+WQ[1]*ev[47];
    ev[166]=QC[2]*ev[41]+WQ[2]*ev[47]+     ze2*ev[16];
    ev[167]=QC[0]*ev[42]+WQ[0]*ev[48];
    ev[168]=QC[1]*ev[42]+WQ[1]*ev[48]+     ze2*ev[18];
    ev[169]=QC[2]*ev[42]+WQ[2]*ev[48]+     ze2*ev[17];
    // (fs|ps) m=[0,2]
    ev[170]=QC[0]*ev[55]+WQ[0]*ev[65]+3.e0*ze2*ev[31];
    ev[171]=QC[1]*ev[55]+WQ[1]*ev[65];
    ev[172]=QC[2]*ev[55]+WQ[2]*ev[65];
    ev[173]=QC[0]*ev[56]+WQ[0]*ev[66];
    ev[174]=QC[1]*ev[56]+WQ[1]*ev[66]+3.e0*ze2*ev[32];
    ev[175]=QC[2]*ev[56]+WQ[2]*ev[66];
    ev[176]=QC[0]*ev[57]+WQ[0]*ev[67];
    ev[177]=QC[1]*ev[57]+WQ[1]*ev[67];
    ev[178]=QC[2]*ev[57]+WQ[2]*ev[67]+3.e0*ze2*ev[33];
    ev[179]=QC[0]*ev[58]+WQ[0]*ev[68]+2.e0*ze2*ev[34];
    ev[180]=QC[1]*ev[58]+WQ[1]*ev[68]+     ze2*ev[31];
    ev[181]=QC[2]*ev[58]+WQ[2]*ev[68];
    ev[182]=QC[0]*ev[59]+WQ[0]*ev[69]+2.e0*ze2*ev[35];
    ev[183]=QC[1]*ev[59]+WQ[1]*ev[69];
    ev[184]=QC[2]*ev[59]+WQ[2]*ev[69]+     ze2*ev[31];
    ev[185]=QC[0]*ev[60]+WQ[0]*ev[70]+     ze2*ev[32];
    ev[186]=QC[1]*ev[60]+WQ[1]*ev[70]+2.e0*ze2*ev[34];
    ev[187]=QC[2]*ev[60]+WQ[2]*ev[70];
    ev[188]=QC[0]*ev[61]+WQ[0]*ev[71]+     ze2*ev[33];
    ev[189]=QC[1]*ev[61]+WQ[1]*ev[71];
    ev[190]=QC[2]*ev[61]+WQ[2]*ev[71]+2.e0*ze2*ev[35];
    ev[191]=QC[0]*ev[62]+WQ[0]*ev[72]+     ze2*ev[36];
    ev[192]=QC[1]*ev[62]+WQ[1]*ev[72]+     ze2*ev[35];
    ev[193]=QC[2]*ev[62]+WQ[2]*ev[72]+     ze2*ev[34];
    ev[194]=QC[0]*ev[63]+WQ[0]*ev[73];
    ev[195]=QC[1]*ev[63]+WQ[1]*ev[73]+2.e0*ze2*ev[36];
    ev[196]=QC[2]*ev[63]+WQ[2]*ev[73]+     ze2*ev[32];
    ev[197]=QC[0]*ev[64]+WQ[0]*ev[74];
    ev[198]=QC[1]*ev[64]+WQ[1]*ev[74]+     ze2*ev[33];
    ev[199]=QC[2]*ev[64]+WQ[2]*ev[74]+2.e0*ze2*ev[36];
    ev[200]=QC[0]*ev[65]+WQ[0]*ev[75]+3.e0*ze2*ev[37];
    ev[201]=QC[1]*ev[65]+WQ[1]*ev[75];
    ev[202]=QC[2]*ev[65]+WQ[2]*ev[75];
    ev[203]=QC[0]*ev[66]+WQ[0]*ev[76];
    ev[204]=QC[1]*ev[66]+WQ[1]*ev[76]+3.e0*ze2*ev[38];
    ev[205]=QC[2]*ev[66]+WQ[2]*ev[76];
    ev[206]=QC[0]*ev[67]+WQ[0]*ev[77];
    ev[207]=QC[1]*ev[67]+WQ[1]*ev[77];
    ev[208]=QC[2]*ev[67]+WQ[2]*ev[77]+3.e0*ze2*ev[39];
    ev[209]=QC[0]*ev[68]+WQ[0]*ev[78]+2.e0*ze2*ev[40];
    ev[210]=QC[1]*ev[68]+WQ[1]*ev[78]+     ze2*ev[37];
    ev[211]=QC[2]*ev[68]+WQ[2]*ev[78];
    ev[212]=QC[0]*ev[69]+WQ[0]*ev[79]+2.e0*ze2*ev[41];
    ev[213]=QC[1]*ev[69]+WQ[1]*ev[79];
    ev[214]=QC[2]*ev[69]+WQ[2]*ev[79]+     ze2*ev[37];
    ev[215]=QC[0]*ev[70]+WQ[0]*ev[80]+     ze2*ev[38];
    ev[216]=QC[1]*ev[70]+WQ[1]*ev[80]+2.e0*ze2*ev[40];
    ev[217]=QC[2]*ev[70]+WQ[2]*ev[80];
    ev[218]=QC[0]*ev[71]+WQ[0]*ev[81]+     ze2*ev[39];
    ev[219]=QC[1]*ev[71]+WQ[1]*ev[81];
    ev[220]=QC[2]*ev[71]+WQ[2]*ev[81]+2.e0*ze2*ev[41];
    ev[221]=QC[0]*ev[72]+WQ[0]*ev[82]+     ze2*ev[42];
    ev[222]=QC[1]*ev[72]+WQ[1]*ev[82]+     ze2*ev[41];
    ev[223]=QC[2]*ev[72]+WQ[2]*ev[82]+     ze2*ev[40];
    ev[224]=QC[0]*ev[73]+WQ[0]*ev[83];
    ev[225]=QC[1]*ev[73]+WQ[1]*ev[83]+2.e0*ze2*ev[42];
    ev[226]=QC[2]*ev[73]+WQ[2]*ev[83]+     ze2*ev[38];
    ev[227]=QC[0]*ev[74]+WQ[0]*ev[84];
    ev[228]=QC[1]*ev[74]+WQ[1]*ev[84]+     ze2*ev[39];
    ev[229]=QC[2]*ev[74]+WQ[2]*ev[84]+2.e0*ze2*ev[42];
    ev[230]=QC[0]*ev[75]+WQ[0]*ev[85]+3.e0*ze2*ev[43];
    ev[231]=QC[1]*ev[75]+WQ[1]*ev[85];
    ev[232]=QC[2]*ev[75]+WQ[2]*ev[85];
    ev[233]=QC[0]*ev[76]+WQ[0]*ev[86];
    ev[234]=QC[1]*ev[76]+WQ[1]*ev[86]+3.e0*ze2*ev[44];
    ev[235]=QC[2]*ev[76]+WQ[2]*ev[86];
    ev[236]=QC[0]*ev[77]+WQ[0]*ev[87];
    ev[237]=QC[1]*ev[77]+WQ[1]*ev[87];
    ev[238]=QC[2]*ev[77]+WQ[2]*ev[87]+3.e0*ze2*ev[45];
    ev[239]=QC[0]*ev[78]+WQ[0]*ev[88]+2.e0*ze2*ev[46];
    ev[240]=QC[1]*ev[78]+WQ[1]*ev[88]+     ze2*ev[43];
    ev[241]=QC[2]*ev[78]+WQ[2]*ev[88];
    ev[242]=QC[0]*ev[79]+WQ[0]*ev[89]+2.e0*ze2*ev[47];
    ev[243]=QC[1]*ev[79]+WQ[1]*ev[89];
    ev[244]=QC[2]*ev[79]+WQ[2]*ev[89]+     ze2*ev[43];
    ev[245]=QC[0]*ev[80]+WQ[0]*ev[90]+     ze2*ev[44];
    ev[246]=QC[1]*ev[80]+WQ[1]*ev[90]+2.e0*ze2*ev[46];
    ev[247]=QC[2]*ev[80]+WQ[2]*ev[90];
    ev[248]=QC[0]*ev[81]+WQ[0]*ev[91]+     ze2*ev[45];
    ev[249]=QC[1]*ev[81]+WQ[1]*ev[91];
    ev[250]=QC[2]*ev[81]+WQ[2]*ev[91]+2.e0*ze2*ev[47];
    ev[251]=QC[0]*ev[82]+WQ[0]*ev[92]+     ze2*ev[48];
    ev[252]=QC[1]*ev[82]+WQ[1]*ev[92]+     ze2*ev[47];
    ev[253]=QC[2]*ev[82]+WQ[2]*ev[92]+     ze2*ev[46];
    ev[254]=QC[0]*ev[83]+WQ[0]*ev[93];
    ev[255]=QC[1]*ev[83]+WQ[1]*ev[93]+2.e0*ze2*ev[48];
    ev[256]=QC[2]*ev[83]+WQ[2]*ev[93]+     ze2*ev[44];
    ev[257]=QC[0]*ev[84]+WQ[0]*ev[94];
    ev[258]=QC[1]*ev[84]+WQ[1]*ev[94]+     ze2*ev[45];
    ev[259]=QC[2]*ev[84]+WQ[2]*ev[94]+2.e0*ze2*ev[48];
    // (ps|ds) m=[1,1]
    ev[260]=QC[0]*ev[ 98]+WQ[0]*ev[107]+eta2*(ev[10]-re*ev[13])+ze2*ev[95]
            ;
    ev[261]=QC[1]*ev[ 99]+WQ[1]*ev[108]+eta2*(ev[10]-re*ev[13]);
    ev[262]=QC[2]*ev[100]+WQ[2]*ev[109]+eta2*(ev[10]-re*ev[13]);
    ev[263]=QC[0]*ev[ 99]+WQ[0]*ev[108]+ze2*ev[96];
    ev[264]=QC[0]*ev[100]+WQ[0]*ev[109]+ze2*ev[97];
    ev[265]=QC[1]*ev[100]+WQ[1]*ev[109];
    ev[266]=QC[0]*ev[101]+WQ[0]*ev[110]+eta2*(ev[11]-re*ev[14]);
    ev[267]=QC[1]*ev[102]+WQ[1]*ev[111]+eta2*(ev[11]-re*ev[14])+ze2*ev[96]
            ;
    ev[268]=QC[2]*ev[103]+WQ[2]*ev[112]+eta2*(ev[11]-re*ev[14]);
    ev[269]=QC[0]*ev[102]+WQ[0]*ev[111];
    ev[270]=QC[0]*ev[103]+WQ[0]*ev[112];
    ev[271]=QC[1]*ev[103]+WQ[1]*ev[112]+ze2*ev[97];
    ev[272]=QC[0]*ev[104]+WQ[0]*ev[113]+eta2*(ev[12]-re*ev[15]);
    ev[273]=QC[1]*ev[105]+WQ[1]*ev[114]+eta2*(ev[12]-re*ev[15]);
    ev[274]=QC[2]*ev[106]+WQ[2]*ev[115]+eta2*(ev[12]-re*ev[15])+ze2*ev[97]
            ;
    ev[275]=QC[0]*ev[105]+WQ[0]*ev[114];
    ev[276]=QC[0]*ev[106]+WQ[0]*ev[115];
    ev[277]=QC[1]*ev[106]+WQ[1]*ev[115];
    // (ds|ds) m=[0,1]
    ev[278]=QC[0]*ev[116]+WQ[0]*ev[134]+eta2*(ev[25]-re*ev[31])
            +2.e0*ze2*ev[ 98];
    ev[279]=QC[1]*ev[117]+WQ[1]*ev[135]+eta2*(ev[25]-re*ev[31]);
    ev[280]=QC[2]*ev[118]+WQ[2]*ev[136]+eta2*(ev[25]-re*ev[31]);
    ev[281]=QC[0]*ev[117]+WQ[0]*ev[135]+2.e0*ze2*ev[ 99];
    ev[282]=QC[0]*ev[118]+WQ[0]*ev[136]+2.e0*ze2*ev[100];
    ev[283]=QC[1]*ev[118]+WQ[1]*ev[136];
    ev[284]=QC[0]*ev[119]+WQ[0]*ev[137]+eta2*(ev[26]-re*ev[32]);
    ev[285]=QC[1]*ev[120]+WQ[1]*ev[138]+eta2*(ev[26]-re*ev[32])
            +2.e0*ze2*ev[102];
    ev[286]=QC[2]*ev[121]+WQ[2]*ev[139]+eta2*(ev[26]-re*ev[32]);
    ev[287]=QC[0]*ev[120]+WQ[0]*ev[138];
    ev[288]=QC[0]*ev[121]+WQ[0]*ev[139];
    ev[289]=QC[1]*ev[121]+WQ[1]*ev[139]+2.e0*ze2*ev[103];
    ev[290]=QC[0]*ev[122]+WQ[0]*ev[140]+eta2*(ev[27]-re*ev[33]);
    ev[291]=QC[1]*ev[123]+WQ[1]*ev[141]+eta2*(ev[27]-re*ev[33]);
    ev[292]=QC[2]*ev[124]+WQ[2]*ev[142]+eta2*(ev[27]-re*ev[33])
            +2.e0*ze2*ev[106];
    ev[293]=QC[0]*ev[123]+WQ[0]*ev[141];
    ev[294]=QC[0]*ev[124]+WQ[0]*ev[142];
    ev[295]=QC[1]*ev[124]+WQ[1]*ev[142];
    ev[296]=QC[0]*ev[125]+WQ[0]*ev[143]+eta2*(ev[28]-re*ev[34])
            +     ze2*ev[101];
    ev[297]=QC[1]*ev[126]+WQ[1]*ev[144]+eta2*(ev[28]-re*ev[34])
            +     ze2*ev[ 99];
    ev[298]=QC[2]*ev[127]+WQ[2]*ev[145]+eta2*(ev[28]-re*ev[34]);
    ev[299]=QC[0]*ev[126]+WQ[0]*ev[144]+     ze2*ev[102];
    ev[300]=QC[0]*ev[127]+WQ[0]*ev[145]+     ze2*ev[103];
    ev[301]=QC[1]*ev[127]+WQ[1]*ev[145]+     ze2*ev[100];
    ev[302]=QC[0]*ev[128]+WQ[0]*ev[146]+eta2*(ev[29]-re*ev[35])
            +     ze2*ev[104];
    ev[303]=QC[1]*ev[129]+WQ[1]*ev[147]+eta2*(ev[29]-re*ev[35]);
    ev[304]=QC[2]*ev[130]+WQ[2]*ev[148]+eta2*(ev[29]-re*ev[35])
            +     ze2*ev[100];
    ev[305]=QC[0]*ev[129]+WQ[0]*ev[147]+     ze2*ev[105];
    ev[306]=QC[0]*ev[130]+WQ[0]*ev[148]+     ze2*ev[106];
    ev[307]=QC[1]*ev[130]+WQ[1]*ev[148];
    ev[308]=QC[0]*ev[131]+WQ[0]*ev[149]+eta2*(ev[30]-re*ev[36]);
    ev[309]=QC[1]*ev[132]+WQ[1]*ev[150]+eta2*(ev[30]-re*ev[36])
            +     ze2*ev[105];
    ev[310]=QC[2]*ev[133]+WQ[2]*ev[151]+eta2*(ev[30]-re*ev[36])
            +     ze2*ev[103];
    ev[311]=QC[0]*ev[132]+WQ[0]*ev[150];
    ev[312]=QC[0]*ev[133]+WQ[0]*ev[151];
    ev[313]=QC[1]*ev[133]+WQ[1]*ev[151]+     ze2*ev[106];
    ev[314]=QC[0]*ev[134]+WQ[0]*ev[152]+eta2*(ev[31]-re*ev[37])
            +2.e0*ze2*ev[107];
    ev[315]=QC[1]*ev[135]+WQ[1]*ev[153]+eta2*(ev[31]-re*ev[37]);
    ev[316]=QC[2]*ev[136]+WQ[2]*ev[154]+eta2*(ev[31]-re*ev[37]);
    ev[317]=QC[0]*ev[135]+WQ[0]*ev[153]+2.e0*ze2*ev[108];
    ev[318]=QC[0]*ev[136]+WQ[0]*ev[154]+2.e0*ze2*ev[109];
    ev[319]=QC[1]*ev[136]+WQ[1]*ev[154];
    ev[320]=QC[0]*ev[137]+WQ[0]*ev[155]+eta2*(ev[32]-re*ev[38]);
    ev[321]=QC[1]*ev[138]+WQ[1]*ev[156]+eta2*(ev[32]-re*ev[38])
            +2.e0*ze2*ev[111];
    ev[322]=QC[2]*ev[139]+WQ[2]*ev[157]+eta2*(ev[32]-re*ev[38]);
    ev[323]=QC[0]*ev[138]+WQ[0]*ev[156];
    ev[324]=QC[0]*ev[139]+WQ[0]*ev[157];
    ev[325]=QC[1]*ev[139]+WQ[1]*ev[157]+2.e0*ze2*ev[112];
    ev[326]=QC[0]*ev[140]+WQ[0]*ev[158]+eta2*(ev[33]-re*ev[39]);
    ev[327]=QC[1]*ev[141]+WQ[1]*ev[159]+eta2*(ev[33]-re*ev[39]);
    ev[328]=QC[2]*ev[142]+WQ[2]*ev[160]+eta2*(ev[33]-re*ev[39])
            +2.e0*ze2*ev[115];
    ev[329]=QC[0]*ev[141]+WQ[0]*ev[159];
    ev[330]=QC[0]*ev[142]+WQ[0]*ev[160];
    ev[331]=QC[1]*ev[142]+WQ[1]*ev[160];
    ev[332]=QC[0]*ev[143]+WQ[0]*ev[161]+eta2*(ev[34]-re*ev[40])
            +     ze2*ev[110];
    ev[333]=QC[1]*ev[144]+WQ[1]*ev[162]+eta2*(ev[34]-re*ev[40])
            +     ze2*ev[108];
    ev[334]=QC[2]*ev[145]+WQ[2]*ev[163]+eta2*(ev[34]-re*ev[40]);
    ev[335]=QC[0]*ev[144]+WQ[0]*ev[162]+     ze2*ev[111];
    ev[336]=QC[0]*ev[145]+WQ[0]*ev[163]+     ze2*ev[112];
    ev[337]=QC[1]*ev[145]+WQ[1]*ev[163]+     ze2*ev[109];
    ev[338]=QC[0]*ev[146]+WQ[0]*ev[164]+eta2*(ev[35]-re*ev[41])
            +     ze2*ev[113];
    ev[339]=QC[1]*ev[147]+WQ[1]*ev[165]+eta2*(ev[35]-re*ev[41]);
    ev[340]=QC[2]*ev[148]+WQ[2]*ev[166]+eta2*(ev[35]-re*ev[41])
            +     ze2*ev[109];
    ev[341]=QC[0]*ev[147]+WQ[0]*ev[165]+     ze2*ev[114];
    ev[342]=QC[0]*ev[148]+WQ[0]*ev[166]+     ze2*ev[115];
    ev[343]=QC[1]*ev[148]+WQ[1]*ev[166];
    ev[344]=QC[0]*ev[149]+WQ[0]*ev[167]+eta2*(ev[36]-re*ev[42]);
    ev[345]=QC[1]*ev[150]+WQ[1]*ev[168]+eta2*(ev[36]-re*ev[42])
            +     ze2*ev[114];
    ev[346]=QC[2]*ev[151]+WQ[2]*ev[169]+eta2*(ev[36]-re*ev[42])
            +     ze2*ev[112];
    ev[347]=QC[0]*ev[150]+WQ[0]*ev[168];
    ev[348]=QC[0]*ev[151]+WQ[0]*ev[169];
    ev[349]=QC[1]*ev[151]+WQ[1]*ev[169]+     ze2*ev[115];
    // (fs|ds) m=[0,1]
    ev[350]=QC[0]*ev[170]+WQ[0]*ev[200]+eta2*(ev[55]-re*ev[65])
            +3.e0*ze2*ev[134];
    ev[351]=QC[1]*ev[171]+WQ[1]*ev[201]+eta2*(ev[55]-re*ev[65]);
    ev[352]=QC[2]*ev[172]+WQ[2]*ev[202]+eta2*(ev[55]-re*ev[65]);
    ev[353]=QC[0]*ev[171]+WQ[0]*ev[201]+3.e0*ze2*ev[135];
    ev[354]=QC[0]*ev[172]+WQ[0]*ev[202]+3.e0*ze2*ev[136];
    ev[355]=QC[1]*ev[172]+WQ[1]*ev[202];
    ev[356]=QC[0]*ev[173]+WQ[0]*ev[203]+eta2*(ev[56]-re*ev[66]);
    ev[357]=QC[1]*ev[174]+WQ[1]*ev[204]+eta2*(ev[56]-re*ev[66])
            +3.e0*ze2*ev[138];
    ev[358]=QC[2]*ev[175]+WQ[2]*ev[205]+eta2*(ev[56]-re*ev[66]);
    ev[359]=QC[0]*ev[174]+WQ[0]*ev[204];
    ev[360]=QC[0]*ev[175]+WQ[0]*ev[205];
    ev[361]=QC[1]*ev[175]+WQ[1]*ev[205]+3.e0*ze2*ev[139];
    ev[362]=QC[0]*ev[176]+WQ[0]*ev[206]+eta2*(ev[57]-re*ev[67]);
    ev[363]=QC[1]*ev[177]+WQ[1]*ev[207]+eta2*(ev[57]-re*ev[67]);
    ev[364]=QC[2]*ev[178]+WQ[2]*ev[208]+eta2*(ev[57]-re*ev[67])
            +3.e0*ze2*ev[142];
    ev[365]=QC[0]*ev[177]+WQ[0]*ev[207];
    ev[366]=QC[0]*ev[178]+WQ[0]*ev[208];
    ev[367]=QC[1]*ev[178]+WQ[1]*ev[208];
    ev[368]=QC[0]*ev[179]+WQ[0]*ev[209]+eta2*(ev[58]-re*ev[68])
            +2.e0*ze2*ev[143];
    ev[369]=QC[1]*ev[180]+WQ[1]*ev[210]+eta2*(ev[58]-re*ev[68])
            +     ze2*ev[135];
    ev[370]=QC[2]*ev[181]+WQ[2]*ev[211]+eta2*(ev[58]-re*ev[68]);
    ev[371]=QC[0]*ev[180]+WQ[0]*ev[210]+2.e0*ze2*ev[144];
    ev[372]=QC[0]*ev[181]+WQ[0]*ev[211]+2.e0*ze2*ev[145];
    ev[373]=QC[1]*ev[181]+WQ[1]*ev[211]+     ze2*ev[136];
    ev[374]=QC[0]*ev[182]+WQ[0]*ev[212]+eta2*(ev[59]-re*ev[69])
            +2.e0*ze2*ev[146];
    ev[375]=QC[1]*ev[183]+WQ[1]*ev[213]+eta2*(ev[59]-re*ev[69]);
    ev[376]=QC[2]*ev[184]+WQ[2]*ev[214]+eta2*(ev[59]-re*ev[69])
            +     ze2*ev[136];
    ev[377]=QC[0]*ev[183]+WQ[0]*ev[213]+2.e0*ze2*ev[147];
    ev[378]=QC[0]*ev[184]+WQ[0]*ev[214]+2.e0*ze2*ev[148];
    ev[379]=QC[1]*ev[184]+WQ[1]*ev[214];
    ev[380]=QC[0]*ev[185]+WQ[0]*ev[215]+eta2*(ev[60]-re*ev[70])
            +     ze2*ev[137];
    ev[381]=QC[1]*ev[186]+WQ[1]*ev[216]+eta2*(ev[60]-re*ev[70])
            +2.e0*ze2*ev[144];
    ev[382]=QC[2]*ev[187]+WQ[2]*ev[217]+eta2*(ev[60]-re*ev[70]);
    ev[383]=QC[0]*ev[186]+WQ[0]*ev[216]+     ze2*ev[138];
    ev[384]=QC[0]*ev[187]+WQ[0]*ev[217]+     ze2*ev[139];
    ev[385]=QC[1]*ev[187]+WQ[1]*ev[217]+2.e0*ze2*ev[145];
    ev[386]=QC[0]*ev[188]+WQ[0]*ev[218]+eta2*(ev[61]-re*ev[71])
            +     ze2*ev[140];
    ev[387]=QC[1]*ev[189]+WQ[1]*ev[219]+eta2*(ev[61]-re*ev[71]);
    ev[388]=QC[2]*ev[190]+WQ[2]*ev[220]+eta2*(ev[61]-re*ev[71])
            +2.e0*ze2*ev[148];
    ev[389]=QC[0]*ev[189]+WQ[0]*ev[219]+     ze2*ev[141];
    ev[390]=QC[0]*ev[190]+WQ[0]*ev[220]+     ze2*ev[142];
    ev[391]=QC[1]*ev[190]+WQ[1]*ev[220];
    ev[392]=QC[0]*ev[191]+WQ[0]*ev[221]+eta2*(ev[62]-re*ev[72])
            +     ze2*ev[149];
    ev[393]=QC[1]*ev[192]+WQ[1]*ev[222]+eta2*(ev[62]-re*ev[72])
            +     ze2*ev[147];
    ev[394]=QC[2]*ev[193]+WQ[2]*ev[223]+eta2*(ev[62]-re*ev[72])
            +     ze2*ev[145];
    ev[395]=QC[0]*ev[192]+WQ[0]*ev[222]+     ze2*ev[150];
    ev[396]=QC[0]*ev[193]+WQ[0]*ev[223]+     ze2*ev[151];
    ev[397]=QC[1]*ev[193]+WQ[1]*ev[223]+     ze2*ev[148];
    ev[398]=QC[0]*ev[194]+WQ[0]*ev[224]+eta2*(ev[63]-re*ev[73]);
    ev[399]=QC[1]*ev[195]+WQ[1]*ev[225]+eta2*(ev[63]-re*ev[73])
            +2.e0*ze2*ev[150];
    ev[400]=QC[2]*ev[196]+WQ[2]*ev[226]+eta2*(ev[63]-re*ev[73])
            +     ze2*ev[139];
    ev[401]=QC[0]*ev[195]+WQ[0]*ev[225];
    ev[402]=QC[0]*ev[196]+WQ[0]*ev[226];
    ev[403]=QC[1]*ev[196]+WQ[1]*ev[226]+2.e0*ze2*ev[151];
    ev[404]=QC[0]*ev[197]+WQ[0]*ev[227]+eta2*(ev[64]-re*ev[74]);
    ev[405]=QC[1]*ev[198]+WQ[1]*ev[228]+eta2*(ev[64]-re*ev[74])
            +     ze2*ev[141];
    ev[406]=QC[2]*ev[199]+WQ[2]*ev[229]+eta2*(ev[64]-re*ev[74])
            +2.e0*ze2*ev[151];
    ev[407]=QC[0]*ev[198]+WQ[0]*ev[228];
    ev[408]=QC[0]*ev[199]+WQ[0]*ev[229];
    ev[409]=QC[1]*ev[199]+WQ[1]*ev[229]+     ze2*ev[142];
    ev[410]=QC[0]*ev[200]+WQ[0]*ev[230]+eta2*(ev[65]-re*ev[75])
            +3.e0*ze2*ev[152];
    ev[411]=QC[1]*ev[201]+WQ[1]*ev[231]+eta2*(ev[65]-re*ev[75]);
    ev[412]=QC[2]*ev[202]+WQ[2]*ev[232]+eta2*(ev[65]-re*ev[75]);
    ev[413]=QC[0]*ev[201]+WQ[0]*ev[231]+3.e0*ze2*ev[153];
    ev[414]=QC[0]*ev[202]+WQ[0]*ev[232]+3.e0*ze2*ev[154];
    ev[415]=QC[1]*ev[202]+WQ[1]*ev[232];
    ev[416]=QC[0]*ev[203]+WQ[0]*ev[233]+eta2*(ev[66]-re*ev[76]);
    ev[417]=QC[1]*ev[204]+WQ[1]*ev[234]+eta2*(ev[66]-re*ev[76])
            +3.e0*ze2*ev[156];
    ev[418]=QC[2]*ev[205]+WQ[2]*ev[235]+eta2*(ev[66]-re*ev[76]);
    ev[419]=QC[0]*ev[204]+WQ[0]*ev[234];
    ev[420]=QC[0]*ev[205]+WQ[0]*ev[235];
    ev[421]=QC[1]*ev[205]+WQ[1]*ev[235]+3.e0*ze2*ev[157];
    ev[422]=QC[0]*ev[206]+WQ[0]*ev[236]+eta2*(ev[67]-re*ev[77]);
    ev[423]=QC[1]*ev[207]+WQ[1]*ev[237]+eta2*(ev[67]-re*ev[77]);
    ev[424]=QC[2]*ev[208]+WQ[2]*ev[238]+eta2*(ev[67]-re*ev[77])
            +3.e0*ze2*ev[160];
    ev[425]=QC[0]*ev[207]+WQ[0]*ev[237];
    ev[426]=QC[0]*ev[208]+WQ[0]*ev[238];
    ev[427]=QC[1]*ev[208]+WQ[1]*ev[238];
    ev[428]=QC[0]*ev[209]+WQ[0]*ev[239]+eta2*(ev[68]-re*ev[78])
            +2.e0*ze2*ev[161];
    ev[429]=QC[1]*ev[210]+WQ[1]*ev[240]+eta2*(ev[68]-re*ev[78])
            +     ze2*ev[153];
    ev[430]=QC[2]*ev[211]+WQ[2]*ev[241]+eta2*(ev[68]-re*ev[78]);
    ev[431]=QC[0]*ev[210]+WQ[0]*ev[240]+2.e0*ze2*ev[162];
    ev[432]=QC[0]*ev[211]+WQ[0]*ev[241]+2.e0*ze2*ev[163];
    ev[433]=QC[1]*ev[211]+WQ[1]*ev[241]+     ze2*ev[154];
    ev[434]=QC[0]*ev[212]+WQ[0]*ev[242]+eta2*(ev[69]-re*ev[79])
            +2.e0*ze2*ev[164];
    ev[435]=QC[1]*ev[213]+WQ[1]*ev[243]+eta2*(ev[69]-re*ev[79]);
    ev[436]=QC[2]*ev[214]+WQ[2]*ev[244]+eta2*(ev[69]-re*ev[79])
            +     ze2*ev[154];
    ev[437]=QC[0]*ev[213]+WQ[0]*ev[243]+2.e0*ze2*ev[165];
    ev[438]=QC[0]*ev[214]+WQ[0]*ev[244]+2.e0*ze2*ev[166];
    ev[439]=QC[1]*ev[214]+WQ[1]*ev[244];
    ev[440]=QC[0]*ev[215]+WQ[0]*ev[245]+eta2*(ev[70]-re*ev[80])
            +     ze2*ev[155];
    ev[441]=QC[1]*ev[216]+WQ[1]*ev[246]+eta2*(ev[70]-re*ev[80])
            +2.e0*ze2*ev[162];
    ev[442]=QC[2]*ev[217]+WQ[2]*ev[247]+eta2*(ev[70]-re*ev[80]);
    ev[443]=QC[0]*ev[216]+WQ[0]*ev[246]+     ze2*ev[156];
    ev[444]=QC[0]*ev[217]+WQ[0]*ev[247]+     ze2*ev[157];
    ev[445]=QC[1]*ev[217]+WQ[1]*ev[247]+2.e0*ze2*ev[163];
    ev[446]=QC[0]*ev[218]+WQ[0]*ev[248]+eta2*(ev[71]-re*ev[81])
            +     ze2*ev[158];
    ev[447]=QC[1]*ev[219]+WQ[1]*ev[249]+eta2*(ev[71]-re*ev[81]);
    ev[448]=QC[2]*ev[220]+WQ[2]*ev[250]+eta2*(ev[71]-re*ev[81])
            +2.e0*ze2*ev[166];
    ev[449]=QC[0]*ev[219]+WQ[0]*ev[249]+     ze2*ev[159];
    ev[450]=QC[0]*ev[220]+WQ[0]*ev[250]+     ze2*ev[160];
    ev[451]=QC[1]*ev[220]+WQ[1]*ev[250];
    ev[452]=QC[0]*ev[221]+WQ[0]*ev[251]+eta2*(ev[72]-re*ev[82])
            +     ze2*ev[167];
    ev[453]=QC[1]*ev[222]+WQ[1]*ev[252]+eta2*(ev[72]-re*ev[82])
            +     ze2*ev[165];
    ev[454]=QC[2]*ev[223]+WQ[2]*ev[253]+eta2*(ev[72]-re*ev[82])
            +     ze2*ev[163];
    ev[455]=QC[0]*ev[222]+WQ[0]*ev[252]+     ze2*ev[168];
    ev[456]=QC[0]*ev[223]+WQ[0]*ev[253]+     ze2*ev[169];
    ev[457]=QC[1]*ev[223]+WQ[1]*ev[253]+     ze2*ev[166];
    ev[458]=QC[0]*ev[224]+WQ[0]*ev[254]+eta2*(ev[73]-re*ev[83]);
    ev[459]=QC[1]*ev[225]+WQ[1]*ev[255]+eta2*(ev[73]-re*ev[83])
            +2.e0*ze2*ev[168];
    ev[460]=QC[2]*ev[226]+WQ[2]*ev[256]+eta2*(ev[73]-re*ev[83])
            +     ze2*ev[157];
    ev[461]=QC[0]*ev[225]+WQ[0]*ev[255];
    ev[462]=QC[0]*ev[226]+WQ[0]*ev[256];
    ev[463]=QC[1]*ev[226]+WQ[1]*ev[256]+2.e0*ze2*ev[169];
    ev[464]=QC[0]*ev[227]+WQ[0]*ev[257]+eta2*(ev[74]-re*ev[84]);
    ev[465]=QC[1]*ev[228]+WQ[1]*ev[258]+eta2*(ev[74]-re*ev[84])
            +     ze2*ev[159];
    ev[466]=QC[2]*ev[229]+WQ[2]*ev[259]+eta2*(ev[74]-re*ev[84])
            +2.e0*ze2*ev[169];
    ev[467]=QC[0]*ev[228]+WQ[0]*ev[258];
    ev[468]=QC[0]*ev[229]+WQ[0]*ev[259];
    ev[469]=QC[1]*ev[229]+WQ[1]*ev[259]+     ze2*ev[160];
    // (ds|fs) m=[0,0]
    ev[470]=QC[0]*ev[278]+WQ[0]*ev[314]+2.e0*eta2*(ev[116]-re*ev[134])
            +2.e0*ze2*ev[260];
    ev[471]=QC[1]*ev[279]+WQ[1]*ev[315]+2.e0*eta2*(ev[117]-re*ev[135]);
    ev[472]=QC[2]*ev[280]+WQ[2]*ev[316]+2.e0*eta2*(ev[118]-re*ev[136]);
    ev[473]=QC[1]*ev[278]+WQ[1]*ev[314];
    ev[474]=QC[2]*ev[278]+WQ[2]*ev[314];
    ev[475]=QC[0]*ev[279]+WQ[0]*ev[315]+2.e0*ze2*ev[261];
    ev[476]=QC[0]*ev[280]+WQ[0]*ev[316]+2.e0*ze2*ev[262];
    ev[477]=QC[0]*ev[283]+WQ[0]*ev[319]+2.e0*ze2*ev[265];
    ev[478]=QC[2]*ev[279]+WQ[2]*ev[315];
    ev[479]=QC[1]*ev[280]+WQ[1]*ev[316];
    ev[480]=QC[0]*ev[284]+WQ[0]*ev[320]+2.e0*eta2*(ev[119]-re*ev[137]);
    ev[481]=QC[1]*ev[285]+WQ[1]*ev[321]+2.e0*eta2*(ev[120]-re*ev[138])
            +2.e0*ze2*ev[267];
    ev[482]=QC[2]*ev[286]+WQ[2]*ev[322]+2.e0*eta2*(ev[121]-re*ev[139]);
    ev[483]=QC[1]*ev[284]+WQ[1]*ev[320]+2.e0*ze2*ev[266];
    ev[484]=QC[2]*ev[284]+WQ[2]*ev[320];
    ev[485]=QC[0]*ev[285]+WQ[0]*ev[321];
    ev[486]=QC[0]*ev[286]+WQ[0]*ev[322];
    ev[487]=QC[0]*ev[289]+WQ[0]*ev[325];
    ev[488]=QC[2]*ev[285]+WQ[2]*ev[321];
    ev[489]=QC[1]*ev[286]+WQ[1]*ev[322]+2.e0*ze2*ev[268];
    ev[490]=QC[0]*ev[290]+WQ[0]*ev[326]+2.e0*eta2*(ev[122]-re*ev[140]);
    ev[491]=QC[1]*ev[291]+WQ[1]*ev[327]+2.e0*eta2*(ev[123]-re*ev[141]);
    ev[492]=QC[2]*ev[292]+WQ[2]*ev[328]+2.e0*eta2*(ev[124]-re*ev[142])
            +2.e0*ze2*ev[274];
    ev[493]=QC[1]*ev[290]+WQ[1]*ev[326];
    ev[494]=QC[2]*ev[290]+WQ[2]*ev[326]+2.e0*ze2*ev[272];
    ev[495]=QC[0]*ev[291]+WQ[0]*ev[327];
    ev[496]=QC[0]*ev[292]+WQ[0]*ev[328];
    ev[497]=QC[0]*ev[295]+WQ[0]*ev[331];
    ev[498]=QC[2]*ev[291]+WQ[2]*ev[327]+2.e0*ze2*ev[273];
    ev[499]=QC[1]*ev[292]+WQ[1]*ev[328];
    ev[500]=QC[0]*ev[296]+WQ[0]*ev[332]+2.e0*eta2*(ev[125]-re*ev[143])
            +     ze2*ev[266];
    ev[501]=QC[1]*ev[297]+WQ[1]*ev[333]+2.e0*eta2*(ev[126]-re*ev[144])
            +     ze2*ev[261];
    ev[502]=QC[2]*ev[298]+WQ[2]*ev[334]+2.e0*eta2*(ev[127]-re*ev[145]);
    ev[503]=QC[1]*ev[296]+WQ[1]*ev[332]+     ze2*ev[260];
    ev[504]=QC[2]*ev[296]+WQ[2]*ev[332];
    ev[505]=QC[0]*ev[297]+WQ[0]*ev[333]+     ze2*ev[267];
    ev[506]=QC[0]*ev[298]+WQ[0]*ev[334]+     ze2*ev[268];
    ev[507]=QC[0]*ev[301]+WQ[0]*ev[337]+     ze2*ev[271];
    ev[508]=QC[2]*ev[297]+WQ[2]*ev[333];
    ev[509]=QC[1]*ev[298]+WQ[1]*ev[334]+     ze2*ev[262];
    ev[510]=QC[0]*ev[302]+WQ[0]*ev[338]+2.e0*eta2*(ev[128]-re*ev[146])
            +     ze2*ev[272];
    ev[511]=QC[1]*ev[303]+WQ[1]*ev[339]+2.e0*eta2*(ev[129]-re*ev[147]);
    ev[512]=QC[2]*ev[304]+WQ[2]*ev[340]+2.e0*eta2*(ev[130]-re*ev[148])
            +     ze2*ev[262];
    ev[513]=QC[1]*ev[302]+WQ[1]*ev[338];
    ev[514]=QC[2]*ev[302]+WQ[2]*ev[338]+     ze2*ev[260];
    ev[515]=QC[0]*ev[303]+WQ[0]*ev[339]+     ze2*ev[273];
    ev[516]=QC[0]*ev[304]+WQ[0]*ev[340]+     ze2*ev[274];
    ev[517]=QC[0]*ev[307]+WQ[0]*ev[343]+     ze2*ev[277];
    ev[518]=QC[2]*ev[303]+WQ[2]*ev[339]+     ze2*ev[261];
    ev[519]=QC[1]*ev[304]+WQ[1]*ev[340];
    ev[520]=QC[0]*ev[308]+WQ[0]*ev[344]+2.e0*eta2*(ev[131]-re*ev[149]);
    ev[521]=QC[1]*ev[309]+WQ[1]*ev[345]+2.e0*eta2*(ev[132]-re*ev[150])
            +     ze2*ev[273];
    ev[522]=QC[2]*ev[310]+WQ[2]*ev[346]+2.e0*eta2*(ev[133]-re*ev[151])
            +     ze2*ev[268];
    ev[523]=QC[1]*ev[308]+WQ[1]*ev[344]+     ze2*ev[272];
    ev[524]=QC[2]*ev[308]+WQ[2]*ev[344]+     ze2*ev[266];
    ev[525]=QC[0]*ev[309]+WQ[0]*ev[345];
    ev[526]=QC[0]*ev[310]+WQ[0]*ev[346];
    ev[527]=QC[0]*ev[313]+WQ[0]*ev[349];
    ev[528]=QC[2]*ev[309]+WQ[2]*ev[345]+     ze2*ev[267];
    ev[529]=QC[1]*ev[310]+WQ[1]*ev[346]+     ze2*ev[274];
    // (fs|fs) m=[0,0]
    ev[530]=QC[0]*ev[350]+WQ[0]*ev[410]+2.e0*eta2*(ev[170]-re*ev[200])
            +3.e0*ze2*ev[314];
    ev[531]=QC[1]*ev[351]+WQ[1]*ev[411]+2.e0*eta2*(ev[171]-re*ev[201]);
    ev[532]=QC[2]*ev[352]+WQ[2]*ev[412]+2.e0*eta2*(ev[172]-re*ev[202]);
    ev[533]=QC[1]*ev[350]+WQ[1]*ev[410];
    ev[534]=QC[2]*ev[350]+WQ[2]*ev[410];
    ev[535]=QC[0]*ev[351]+WQ[0]*ev[411]+3.e0*ze2*ev[315];
    ev[536]=QC[0]*ev[352]+WQ[0]*ev[412]+3.e0*ze2*ev[316];
    ev[537]=QC[0]*ev[355]+WQ[0]*ev[415]+3.e0*ze2*ev[319];
    ev[538]=QC[2]*ev[351]+WQ[2]*ev[411];
    ev[539]=QC[1]*ev[352]+WQ[1]*ev[412];
    ev[540]=QC[0]*ev[356]+WQ[0]*ev[416]+2.e0*eta2*(ev[173]-re*ev[203]);
    ev[541]=QC[1]*ev[357]+WQ[1]*ev[417]+2.e0*eta2*(ev[174]-re*ev[204])
            +3.e0*ze2*ev[321];
    ev[542]=QC[2]*ev[358]+WQ[2]*ev[418]+2.e0*eta2*(ev[175]-re*ev[205]);
    ev[543]=QC[1]*ev[356]+WQ[1]*ev[416]+3.e0*ze2*ev[320];
    ev[544]=QC[2]*ev[356]+WQ[2]*ev[416];
    ev[545]=QC[0]*ev[357]+WQ[0]*ev[417];
    ev[546]=QC[0]*ev[358]+WQ[0]*ev[418];
    ev[547]=QC[0]*ev[361]+WQ[0]*ev[421];
    ev[548]=QC[2]*ev[357]+WQ[2]*ev[417];
    ev[549]=QC[1]*ev[358]+WQ[1]*ev[418]+3.e0*ze2*ev[322];
    ev[550]=QC[0]*ev[362]+WQ[0]*ev[422]+2.e0*eta2*(ev[176]-re*ev[206]);
    ev[551]=QC[1]*ev[363]+WQ[1]*ev[423]+2.e0*eta2*(ev[177]-re*ev[207]);
    ev[552]=QC[2]*ev[364]+WQ[2]*ev[424]+2.e0*eta2*(ev[178]-re*ev[208])
            +3.e0*ze2*ev[328];
    ev[553]=QC[1]*ev[362]+WQ[1]*ev[422];
    ev[554]=QC[2]*ev[362]+WQ[2]*ev[422]+3.e0*ze2*ev[326];
    ev[555]=QC[0]*ev[363]+WQ[0]*ev[423];
    ev[556]=QC[0]*ev[364]+WQ[0]*ev[424];
    ev[557]=QC[0]*ev[367]+WQ[0]*ev[427];
    ev[558]=QC[2]*ev[363]+WQ[2]*ev[423]+3.e0*ze2*ev[327];
    ev[559]=QC[1]*ev[364]+WQ[1]*ev[424];
    ev[560]=QC[0]*ev[368]+WQ[0]*ev[428]+2.e0*eta2*(ev[179]-re*ev[209])
            +2.e0*ze2*ev[332];
    ev[561]=QC[1]*ev[369]+WQ[1]*ev[429]+2.e0*eta2*(ev[180]-re*ev[210])
            +     ze2*ev[315];
    ev[562]=QC[2]*ev[370]+WQ[2]*ev[430]+2.e0*eta2*(ev[181]-re*ev[211]);
    ev[563]=QC[1]*ev[368]+WQ[1]*ev[428]+     ze2*ev[314];
    ev[564]=QC[2]*ev[368]+WQ[2]*ev[428];
    ev[565]=QC[0]*ev[369]+WQ[0]*ev[429]+2.e0*ze2*ev[333];
    ev[566]=QC[0]*ev[370]+WQ[0]*ev[430]+2.e0*ze2*ev[334];
    ev[567]=QC[0]*ev[373]+WQ[0]*ev[433]+2.e0*ze2*ev[337];
    ev[568]=QC[2]*ev[369]+WQ[2]*ev[429];
    ev[569]=QC[1]*ev[370]+WQ[1]*ev[430]+     ze2*ev[316];
    ev[570]=QC[0]*ev[374]+WQ[0]*ev[434]+2.e0*eta2*(ev[182]-re*ev[212])
            +2.e0*ze2*ev[338];
    ev[571]=QC[1]*ev[375]+WQ[1]*ev[435]+2.e0*eta2*(ev[183]-re*ev[213]);
    ev[572]=QC[2]*ev[376]+WQ[2]*ev[436]+2.e0*eta2*(ev[184]-re*ev[214])
            +     ze2*ev[316];
    ev[573]=QC[1]*ev[374]+WQ[1]*ev[434];
    ev[574]=QC[2]*ev[374]+WQ[2]*ev[434]+     ze2*ev[314];
    ev[575]=QC[0]*ev[375]+WQ[0]*ev[435]+2.e0*ze2*ev[339];
    ev[576]=QC[0]*ev[376]+WQ[0]*ev[436]+2.e0*ze2*ev[340];
    ev[577]=QC[0]*ev[379]+WQ[0]*ev[439]+2.e0*ze2*ev[343];
    ev[578]=QC[2]*ev[375]+WQ[2]*ev[435]+     ze2*ev[315];
    ev[579]=QC[1]*ev[376]+WQ[1]*ev[436];
    ev[580]=QC[0]*ev[380]+WQ[0]*ev[440]+2.e0*eta2*(ev[185]-re*ev[215])
            +     ze2*ev[320];
    ev[581]=QC[1]*ev[381]+WQ[1]*ev[441]+2.e0*eta2*(ev[186]-re*ev[216])
            +2.e0*ze2*ev[333];
    ev[582]=QC[2]*ev[382]+WQ[2]*ev[442]+2.e0*eta2*(ev[187]-re*ev[217]);
    ev[583]=QC[1]*ev[380]+WQ[1]*ev[440]+2.e0*ze2*ev[332];
    ev[584]=QC[2]*ev[380]+WQ[2]*ev[440];
    ev[585]=QC[0]*ev[381]+WQ[0]*ev[441]+     ze2*ev[321];
    ev[586]=QC[0]*ev[382]+WQ[0]*ev[442]+     ze2*ev[322];
    ev[587]=QC[0]*ev[385]+WQ[0]*ev[445]+     ze2*ev[325];
    ev[588]=QC[2]*ev[381]+WQ[2]*ev[441];
    ev[589]=QC[1]*ev[382]+WQ[1]*ev[442]+2.e0*ze2*ev[334];
    ev[590]=QC[0]*ev[386]+WQ[0]*ev[446]+2.e0*eta2*(ev[188]-re*ev[218])
            +     ze2*ev[326];
    ev[591]=QC[1]*ev[387]+WQ[1]*ev[447]+2.e0*eta2*(ev[189]-re*ev[219]);
    ev[592]=QC[2]*ev[388]+WQ[2]*ev[448]+2.e0*eta2*(ev[190]-re*ev[220])
            +2.e0*ze2*ev[340];
    ev[593]=QC[1]*ev[386]+WQ[1]*ev[446];
    ev[594]=QC[2]*ev[386]+WQ[2]*ev[446]+2.e0*ze2*ev[338];
    ev[595]=QC[0]*ev[387]+WQ[0]*ev[447]+     ze2*ev[327];
    ev[596]=QC[0]*ev[388]+WQ[0]*ev[448]+     ze2*ev[328];
    ev[597]=QC[0]*ev[391]+WQ[0]*ev[451]+     ze2*ev[331];
    ev[598]=QC[2]*ev[387]+WQ[2]*ev[447]+2.e0*ze2*ev[339];
    ev[599]=QC[1]*ev[388]+WQ[1]*ev[448];
    ev[600]=QC[0]*ev[392]+WQ[0]*ev[452]+2.e0*eta2*(ev[191]-re*ev[221])
            +     ze2*ev[344];
    ev[601]=QC[1]*ev[393]+WQ[1]*ev[453]+2.e0*eta2*(ev[192]-re*ev[222])
            +     ze2*ev[339];
    ev[602]=QC[2]*ev[394]+WQ[2]*ev[454]+2.e0*eta2*(ev[193]-re*ev[223])
            +     ze2*ev[334];
    ev[603]=QC[1]*ev[392]+WQ[1]*ev[452]+     ze2*ev[338];
    ev[604]=QC[2]*ev[392]+WQ[2]*ev[452]+     ze2*ev[332];
    ev[605]=QC[0]*ev[393]+WQ[0]*ev[453]+     ze2*ev[345];
    ev[606]=QC[0]*ev[394]+WQ[0]*ev[454]+     ze2*ev[346];
    ev[607]=QC[0]*ev[397]+WQ[0]*ev[457]+     ze2*ev[349];
    ev[608]=QC[2]*ev[393]+WQ[2]*ev[453]+     ze2*ev[333];
    ev[609]=QC[1]*ev[394]+WQ[1]*ev[454]+     ze2*ev[340];
    ev[610]=QC[0]*ev[398]+WQ[0]*ev[458]+2.e0*eta2*(ev[194]-re*ev[224]);
    ev[611]=QC[1]*ev[399]+WQ[1]*ev[459]+2.e0*eta2*(ev[195]-re*ev[225])
            +2.e0*ze2*ev[345];
    ev[612]=QC[2]*ev[400]+WQ[2]*ev[460]+2.e0*eta2*(ev[196]-re*ev[226])
            +     ze2*ev[322];
    ev[613]=QC[1]*ev[398]+WQ[1]*ev[458]+2.e0*ze2*ev[344];
    ev[614]=QC[2]*ev[398]+WQ[2]*ev[458]+     ze2*ev[320];
    ev[615]=QC[0]*ev[399]+WQ[0]*ev[459];
    ev[616]=QC[0]*ev[400]+WQ[0]*ev[460];
    ev[617]=QC[0]*ev[403]+WQ[0]*ev[463];
    ev[618]=QC[2]*ev[399]+WQ[2]*ev[459]+     ze2*ev[321];
    ev[619]=QC[1]*ev[400]+WQ[1]*ev[460]+2.e0*ze2*ev[346];
    ev[620]=QC[0]*ev[404]+WQ[0]*ev[464]+2.e0*eta2*(ev[197]-re*ev[227]);
    ev[621]=QC[1]*ev[405]+WQ[1]*ev[465]+2.e0*eta2*(ev[198]-re*ev[228])
            +     ze2*ev[327];
    ev[622]=QC[2]*ev[406]+WQ[2]*ev[466]+2.e0*eta2*(ev[199]-re*ev[229])
            +2.e0*ze2*ev[346];
    ev[623]=QC[1]*ev[404]+WQ[1]*ev[464]+     ze2*ev[326];
    ev[624]=QC[2]*ev[404]+WQ[2]*ev[464]+2.e0*ze2*ev[344];
    ev[625]=QC[0]*ev[405]+WQ[0]*ev[465];
    ev[626]=QC[0]*ev[406]+WQ[0]*ev[466];
    ev[627]=QC[0]*ev[409]+WQ[0]*ev[469];
    ev[628]=QC[2]*ev[405]+WQ[2]*ev[465]+2.e0*ze2*ev[345];
    ev[629]=QC[1]*ev[406]+WQ[1]*ev[466]+     ze2*ev[328];
}

__device__ void gpu_vrr_cint_dpdp( const double *ev, double *eh ) {
    int La=2, Lb=1, Lc=2, Ld=1;
    int i, ih, iv;
    // (DS|DS)
#pragma unroll
    for ( i=0, iv=278, ih=0; i<36; i++, iv++, ih++ ) eh[ih]+=ev[iv];
    // (DS|FS)
#pragma unroll
    for ( i=0, iv=470, ih=36; i<60; i++, iv++, ih++ ) eh[ih]+=ev[iv];
    // (FS|DS)
#pragma unroll
    for ( i=0, iv=350, ih=96; i<60; i++, iv++, ih++ ) eh[ih]+=ev[iv];
    // (FS|FS)
#pragma unroll
    for ( i=0, iv=530, ih=156; i<100; i++, iv++, ih++ ) eh[ih]+=ev[iv];
}

__device__ void gpu_twoint_core_os_dpdp(
//        const int *pLa, const int *pLb, const int *pLc, const int *pLd,
        const int *nijps, const double vzeta[], const double vdkab[],
        const double vxiza[], const double BA[3],
        const int *nklps, const double veta[], const double vdkcd[],
        const double vxizc[], const double DC[3], const double AC[3],
        double *DINT ) {
    int ijps, klps, i;
    double cssss, zeta, dkab, xiza, eta, xizc, dk, T;
    double zeta2, eta2, ze2, rz, re, PA[3], WP[3], QC[3], WQ[3];
    double PQ2, sqrho, rho, PC[3], QP[3];
    double ev[630], eh[868];
//    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;

//    DFACT = ofmo_getadd_dfact();
    gpu_hrr_clear_dpdp( eh );
    for ( ijps=0; ijps<(*nijps); ijps++ ) {
        zeta  = LDG(vzeta[ijps]);
        dkab  = LDG(vdkab[ijps]);
        xiza  = LDG(vxiza[ijps]);
        zeta2 = HALF * zeta;
#pragma unroll
        for ( i=0; i<3; i++ ) {
            PC[i] = AC[i] + xiza*BA[i];
            PA[i] = xiza * BA[i];
        }
        for ( klps=0; klps<(*nklps); klps++ ) {
            eta  = LDG(veta[klps]);
            dk   = dkab * LDG(vdkcd[klps]);
            xizc = LDG(vxizc[klps]);
            eta2 = HALF * eta;
            PQ2  = ZERO;
#pragma unroll
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
#pragma unroll
            for ( i=0; i<3; i++ ) {
                WP[i] = rz*QP[i];
                WQ[i] = rz*QP[i] - QP[i];
            }
            T     = rho * PQ2;
            cssss = sqrho * dk;
            gpu_vrr_calc_dpdp(
                    T, cssss, zeta2, eta2, ze2, rz, re, PA, WP, QC, WQ,
                    ev );
            gpu_vrr_cint_dpdp( ev, eh );
        }	// for (klps)
    }	// for (ijps)
    gpu_hrr_calc_dpdp( BA, DC, eh );
    gpu_hrr_coef_dpdp( eh, DINT );
}

#if 0
int ofmo_twoint_os_dpdp(
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
    double DINTEG[6*3*6*3];
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
    max_nzeri = ebuf_max_nzeri - 6*3*6*3;
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
            ofmo_twoint_core_os_dpdp(
                    &La, &Lb, &Lc, &Ld,
                    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
                    &psp_xiza[ijps0], BA,
                    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
                    &psp_xiza[klps0], DC,   AC,      DINTEG );
            ipat=((Lab != Lcd)||(ics==kcs && jcs>lcs) ? true : false );
#ifdef SORT_CSP
            int ijgekl = (ics>kcs);
            if (ics==kcs) ijgekl = (jcs>=lcs);
            if (!ijgekl) ipat = ( (ics==kcs && jcs<lcs) ? true : false);
#endif
            for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
                I2 = (iao*iao+iao)>>1;
                for ( j=0, jao=jao0; j<3; j++, jao++ ) {
                    if ( jao>iao ) { ix+=6*3; continue; }
                    IJ = I2 + jao;
                    coe0 = ( iao==jao ? HALF : ONE );
                    for ( k=0, kao=kao0; k<6; k++, kao++ ) {
                        K2 = (kao*kao+kao)>>1;
                        for ( l=0, lao=lao0; l<3; l++, lao++, ix++ ) {
                            if ( lao>kao ) continue;
                            if ( fabs(DINTEG[ix]) > eps_eri ) {
                                KL = K2 + lao;
#ifndef SORT_CSP
                                if ( IJ >= KL ) {
#else
                                if ((ijgekl&&IJ>=KL) || (!ijgekl&&KL>=IJ)) {
#endif
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

int ofmo_twoint_direct_os_dpdp(
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
    double DINTEG[6*3*6*3];
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

    max_nzeri -= 6*3*6*3;
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
            ofmo_twoint_core_os_dpdp(
                    &La, &Lb, &Lc, &Ld,
                    &nijps, &psp_zeta[ijps0], &psp_dkps[ijps0],
                    &psp_xiza[ijps0], BA,
                    &nklps, &psp_zeta[klps0], &psp_dkps[klps0],
                    &psp_xiza[klps0], DC,   AC,      DINTEG );
            ipat=((Lab != Lcd)||(ics==kcs && jcs>lcs) ? true : false );
#ifdef SORT_CSP
            int ijgekl = (ics>kcs);
            if (ics==kcs) ijgekl = (jcs>=lcs);
            if (!ijgekl) ipat = ( (ics==kcs && jcs<lcs) ? true : false);
#endif
            for ( i=0, iao=iao0, ix=0; i<6; i++, iao++ ) {
                I2 = (iao*iao+iao)>>1;
                for ( j=0, jao=jao0; j<3; j++, jao++ ) {
                    if ( jao>iao ) { ix+=6*3; continue; }
                    IJ = I2 + jao;
                    coe0 = ( iao==jao ? HALF : ONE );
                    for ( k=0, kao=kao0; k<6; k++, kao++ ) {
                        K2 = (kao*kao+kao)>>1;
                        for ( l=0, lao=lao0; l<3; l++, lao++, ix++ ) {
                            if ( lao>kao ) continue;
                            if ( fabs(DINTEG[ix]) > eps_eri ) {
                                KL = K2 + lao;
#ifndef SORT_CSP
                                if ( IJ >= KL ) {
#else
                                if ((ijgekl&&IJ>=KL) || (!ijgekl&&KL>=IJ)) {
#endif
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
#endif // 0
