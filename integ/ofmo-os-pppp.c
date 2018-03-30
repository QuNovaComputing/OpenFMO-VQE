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

static void ofmo_hrr_clear_pppp( double *eh ) {
    int i;
    // (PS|PS)
    for ( i=0; i<(0+9); i++ ) eh[i] = 0.e0;
    // (PS|DS)
    for ( i=9; i<(9+18); i++ ) eh[i] = 0.e0;
    // (DS|PS)
    for ( i=27; i<(27+18); i++ ) eh[i] = 0.e0;
    // (DS|DS)
    for ( i=45; i<(45+36); i++ ) eh[i] = 0.e0;
}

static void ofmo_hrr_calc_pppp(
        const double BA[3], const double DC[3], double *eh ) {
    // (PP,PS)
    eh[  81] = eh[  27] - BA[0]*eh[   0];
    eh[  82] = eh[  28] - BA[0]*eh[   1];
    eh[  83] = eh[  29] - BA[0]*eh[   2];
    eh[  84] = eh[  36] - BA[1]*eh[   0];
    eh[  85] = eh[  37] - BA[1]*eh[   1];
    eh[  86] = eh[  38] - BA[1]*eh[   2];
    eh[  87] = eh[  39] - BA[2]*eh[   0];
    eh[  88] = eh[  40] - BA[2]*eh[   1];
    eh[  89] = eh[  41] - BA[2]*eh[   2];
    eh[  90] = eh[  36] - BA[0]*eh[   3];
    eh[  91] = eh[  37] - BA[0]*eh[   4];
    eh[  92] = eh[  38] - BA[0]*eh[   5];
    eh[  93] = eh[  30] - BA[1]*eh[   3];
    eh[  94] = eh[  31] - BA[1]*eh[   4];
    eh[  95] = eh[  32] - BA[1]*eh[   5];
    eh[  96] = eh[  42] - BA[2]*eh[   3];
    eh[  97] = eh[  43] - BA[2]*eh[   4];
    eh[  98] = eh[  44] - BA[2]*eh[   5];
    eh[  99] = eh[  39] - BA[0]*eh[   6];
    eh[ 100] = eh[  40] - BA[0]*eh[   7];
    eh[ 101] = eh[  41] - BA[0]*eh[   8];
    eh[ 102] = eh[  42] - BA[1]*eh[   6];
    eh[ 103] = eh[  43] - BA[1]*eh[   7];
    eh[ 104] = eh[  44] - BA[1]*eh[   8];
    eh[ 105] = eh[  33] - BA[2]*eh[   6];
    eh[ 106] = eh[  34] - BA[2]*eh[   7];
    eh[ 107] = eh[  35] - BA[2]*eh[   8];
    // (PP,DS)
    eh[ 108] = eh[  45] - BA[0]*eh[   9];
    eh[ 109] = eh[  46] - BA[0]*eh[  10];
    eh[ 110] = eh[  47] - BA[0]*eh[  11];
    eh[ 111] = eh[  48] - BA[0]*eh[  12];
    eh[ 112] = eh[  49] - BA[0]*eh[  13];
    eh[ 113] = eh[  50] - BA[0]*eh[  14];
    eh[ 114] = eh[  63] - BA[1]*eh[   9];
    eh[ 115] = eh[  64] - BA[1]*eh[  10];
    eh[ 116] = eh[  65] - BA[1]*eh[  11];
    eh[ 117] = eh[  66] - BA[1]*eh[  12];
    eh[ 118] = eh[  67] - BA[1]*eh[  13];
    eh[ 119] = eh[  68] - BA[1]*eh[  14];
    eh[ 120] = eh[  69] - BA[2]*eh[   9];
    eh[ 121] = eh[  70] - BA[2]*eh[  10];
    eh[ 122] = eh[  71] - BA[2]*eh[  11];
    eh[ 123] = eh[  72] - BA[2]*eh[  12];
    eh[ 124] = eh[  73] - BA[2]*eh[  13];
    eh[ 125] = eh[  74] - BA[2]*eh[  14];
    eh[ 126] = eh[  63] - BA[0]*eh[  15];
    eh[ 127] = eh[  64] - BA[0]*eh[  16];
    eh[ 128] = eh[  65] - BA[0]*eh[  17];
    eh[ 129] = eh[  66] - BA[0]*eh[  18];
    eh[ 130] = eh[  67] - BA[0]*eh[  19];
    eh[ 131] = eh[  68] - BA[0]*eh[  20];
    eh[ 132] = eh[  51] - BA[1]*eh[  15];
    eh[ 133] = eh[  52] - BA[1]*eh[  16];
    eh[ 134] = eh[  53] - BA[1]*eh[  17];
    eh[ 135] = eh[  54] - BA[1]*eh[  18];
    eh[ 136] = eh[  55] - BA[1]*eh[  19];
    eh[ 137] = eh[  56] - BA[1]*eh[  20];
    eh[ 138] = eh[  75] - BA[2]*eh[  15];
    eh[ 139] = eh[  76] - BA[2]*eh[  16];
    eh[ 140] = eh[  77] - BA[2]*eh[  17];
    eh[ 141] = eh[  78] - BA[2]*eh[  18];
    eh[ 142] = eh[  79] - BA[2]*eh[  19];
    eh[ 143] = eh[  80] - BA[2]*eh[  20];
    eh[ 144] = eh[  69] - BA[0]*eh[  21];
    eh[ 145] = eh[  70] - BA[0]*eh[  22];
    eh[ 146] = eh[  71] - BA[0]*eh[  23];
    eh[ 147] = eh[  72] - BA[0]*eh[  24];
    eh[ 148] = eh[  73] - BA[0]*eh[  25];
    eh[ 149] = eh[  74] - BA[0]*eh[  26];
    eh[ 150] = eh[  75] - BA[1]*eh[  21];
    eh[ 151] = eh[  76] - BA[1]*eh[  22];
    eh[ 152] = eh[  77] - BA[1]*eh[  23];
    eh[ 153] = eh[  78] - BA[1]*eh[  24];
    eh[ 154] = eh[  79] - BA[1]*eh[  25];
    eh[ 155] = eh[  80] - BA[1]*eh[  26];
    eh[ 156] = eh[  57] - BA[2]*eh[  21];
    eh[ 157] = eh[  58] - BA[2]*eh[  22];
    eh[ 158] = eh[  59] - BA[2]*eh[  23];
    eh[ 159] = eh[  60] - BA[2]*eh[  24];
    eh[ 160] = eh[  61] - BA[2]*eh[  25];
    eh[ 161] = eh[  62] - BA[2]*eh[  26];
    // HRR for (XX|XX)-type integral (center CD)
    // (PP,PP)
    eh[ 162] = eh[ 108] - DC[0]*eh[  81];
    eh[ 163] = eh[ 111] - DC[1]*eh[  81];
    eh[ 164] = eh[ 112] - DC[2]*eh[  81];
    eh[ 165] = eh[ 111] - DC[0]*eh[  82];
    eh[ 166] = eh[ 109] - DC[1]*eh[  82];
    eh[ 167] = eh[ 113] - DC[2]*eh[  82];
    eh[ 168] = eh[ 112] - DC[0]*eh[  83];
    eh[ 169] = eh[ 113] - DC[1]*eh[  83];
    eh[ 170] = eh[ 110] - DC[2]*eh[  83];
    eh[ 171] = eh[ 114] - DC[0]*eh[  84];
    eh[ 172] = eh[ 117] - DC[1]*eh[  84];
    eh[ 173] = eh[ 118] - DC[2]*eh[  84];
    eh[ 174] = eh[ 117] - DC[0]*eh[  85];
    eh[ 175] = eh[ 115] - DC[1]*eh[  85];
    eh[ 176] = eh[ 119] - DC[2]*eh[  85];
    eh[ 177] = eh[ 118] - DC[0]*eh[  86];
    eh[ 178] = eh[ 119] - DC[1]*eh[  86];
    eh[ 179] = eh[ 116] - DC[2]*eh[  86];
    eh[ 180] = eh[ 120] - DC[0]*eh[  87];
    eh[ 181] = eh[ 123] - DC[1]*eh[  87];
    eh[ 182] = eh[ 124] - DC[2]*eh[  87];
    eh[ 183] = eh[ 123] - DC[0]*eh[  88];
    eh[ 184] = eh[ 121] - DC[1]*eh[  88];
    eh[ 185] = eh[ 125] - DC[2]*eh[  88];
    eh[ 186] = eh[ 124] - DC[0]*eh[  89];
    eh[ 187] = eh[ 125] - DC[1]*eh[  89];
    eh[ 188] = eh[ 122] - DC[2]*eh[  89];
    eh[ 189] = eh[ 126] - DC[0]*eh[  90];
    eh[ 190] = eh[ 129] - DC[1]*eh[  90];
    eh[ 191] = eh[ 130] - DC[2]*eh[  90];
    eh[ 192] = eh[ 129] - DC[0]*eh[  91];
    eh[ 193] = eh[ 127] - DC[1]*eh[  91];
    eh[ 194] = eh[ 131] - DC[2]*eh[  91];
    eh[ 195] = eh[ 130] - DC[0]*eh[  92];
    eh[ 196] = eh[ 131] - DC[1]*eh[  92];
    eh[ 197] = eh[ 128] - DC[2]*eh[  92];
    eh[ 198] = eh[ 132] - DC[0]*eh[  93];
    eh[ 199] = eh[ 135] - DC[1]*eh[  93];
    eh[ 200] = eh[ 136] - DC[2]*eh[  93];
    eh[ 201] = eh[ 135] - DC[0]*eh[  94];
    eh[ 202] = eh[ 133] - DC[1]*eh[  94];
    eh[ 203] = eh[ 137] - DC[2]*eh[  94];
    eh[ 204] = eh[ 136] - DC[0]*eh[  95];
    eh[ 205] = eh[ 137] - DC[1]*eh[  95];
    eh[ 206] = eh[ 134] - DC[2]*eh[  95];
    eh[ 207] = eh[ 138] - DC[0]*eh[  96];
    eh[ 208] = eh[ 141] - DC[1]*eh[  96];
    eh[ 209] = eh[ 142] - DC[2]*eh[  96];
    eh[ 210] = eh[ 141] - DC[0]*eh[  97];
    eh[ 211] = eh[ 139] - DC[1]*eh[  97];
    eh[ 212] = eh[ 143] - DC[2]*eh[  97];
    eh[ 213] = eh[ 142] - DC[0]*eh[  98];
    eh[ 214] = eh[ 143] - DC[1]*eh[  98];
    eh[ 215] = eh[ 140] - DC[2]*eh[  98];
    eh[ 216] = eh[ 144] - DC[0]*eh[  99];
    eh[ 217] = eh[ 147] - DC[1]*eh[  99];
    eh[ 218] = eh[ 148] - DC[2]*eh[  99];
    eh[ 219] = eh[ 147] - DC[0]*eh[ 100];
    eh[ 220] = eh[ 145] - DC[1]*eh[ 100];
    eh[ 221] = eh[ 149] - DC[2]*eh[ 100];
    eh[ 222] = eh[ 148] - DC[0]*eh[ 101];
    eh[ 223] = eh[ 149] - DC[1]*eh[ 101];
    eh[ 224] = eh[ 146] - DC[2]*eh[ 101];
    eh[ 225] = eh[ 150] - DC[0]*eh[ 102];
    eh[ 226] = eh[ 153] - DC[1]*eh[ 102];
    eh[ 227] = eh[ 154] - DC[2]*eh[ 102];
    eh[ 228] = eh[ 153] - DC[0]*eh[ 103];
    eh[ 229] = eh[ 151] - DC[1]*eh[ 103];
    eh[ 230] = eh[ 155] - DC[2]*eh[ 103];
    eh[ 231] = eh[ 154] - DC[0]*eh[ 104];
    eh[ 232] = eh[ 155] - DC[1]*eh[ 104];
    eh[ 233] = eh[ 152] - DC[2]*eh[ 104];
    eh[ 234] = eh[ 156] - DC[0]*eh[ 105];
    eh[ 235] = eh[ 159] - DC[1]*eh[ 105];
    eh[ 236] = eh[ 160] - DC[2]*eh[ 105];
    eh[ 237] = eh[ 159] - DC[0]*eh[ 106];
    eh[ 238] = eh[ 157] - DC[1]*eh[ 106];
    eh[ 239] = eh[ 161] - DC[2]*eh[ 106];
    eh[ 240] = eh[ 160] - DC[0]*eh[ 107];
    eh[ 241] = eh[ 161] - DC[1]*eh[ 107];
    eh[ 242] = eh[ 158] - DC[2]*eh[ 107];
}

static void ofmo_hrr_coef_pppp(
        const double *eh, double *DINT ) {
    memcpy( DINT, &eh[162], sizeof(double)*3*3*3*3 );
}

static void ofmo_vrr_calc_pppp(
        const double T, const double cssss,
        const double zeta2, const double eta2, const double ze2,
        const double rz, const double re,
        const double PA[3], const double WP[3],
        const double QC[3], const double WQ[3],
        double *ev ) {
    // (ss|ss) m=0,4
    //fmt( &ev[0], 4, T, cssss );
    OFMO_FMT( &ev[0], 4, T, cssss );
    // (ps|ss) m=0,3
    ev[ 5]=PA[0]*ev[0]+WP[0]*ev[1];
    ev[ 6]=PA[1]*ev[0]+WP[1]*ev[1];
    ev[ 7]=PA[2]*ev[0]+WP[2]*ev[1];
    ev[ 8]=PA[0]*ev[1]+WP[0]*ev[2];
    ev[ 9]=PA[1]*ev[1]+WP[1]*ev[2];
    ev[10]=PA[2]*ev[1]+WP[2]*ev[2];
    ev[11]=PA[0]*ev[2]+WP[0]*ev[3];
    ev[12]=PA[1]*ev[2]+WP[1]*ev[3];
    ev[13]=PA[2]*ev[2]+WP[2]*ev[3];
    ev[14]=PA[0]*ev[3]+WP[0]*ev[4];
    ev[15]=PA[1]*ev[3]+WP[1]*ev[4];
    ev[16]=PA[2]*ev[3]+WP[2]*ev[4];
    // (ds|ss) m=0,2
    ev[17]=PA[0]*ev[ 5]+WP[0]*ev[ 8]+zeta2*(ev[0]-rz*ev[1]);
    ev[18]=PA[1]*ev[ 6]+WP[1]*ev[ 9]+zeta2*(ev[0]-rz*ev[1]);
    ev[19]=PA[2]*ev[ 7]+WP[2]*ev[10]+zeta2*(ev[0]-rz*ev[1]);
    ev[20]=PA[0]*ev[ 6]+WP[0]*ev[ 9];
    ev[21]=PA[0]*ev[ 7]+WP[0]*ev[10];
    ev[22]=PA[1]*ev[ 7]+WP[1]*ev[10];
    ev[23]=PA[0]*ev[ 8]+WP[0]*ev[11]+zeta2*(ev[1]-rz*ev[2]);
    ev[24]=PA[1]*ev[ 9]+WP[1]*ev[12]+zeta2*(ev[1]-rz*ev[2]);
    ev[25]=PA[2]*ev[10]+WP[2]*ev[13]+zeta2*(ev[1]-rz*ev[2]);
    ev[26]=PA[0]*ev[ 9]+WP[0]*ev[12];
    ev[27]=PA[0]*ev[10]+WP[0]*ev[13];
    ev[28]=PA[1]*ev[10]+WP[1]*ev[13];
    ev[29]=PA[0]*ev[11]+WP[0]*ev[14]+zeta2*(ev[2]-rz*ev[3]);
    ev[30]=PA[1]*ev[12]+WP[1]*ev[15]+zeta2*(ev[2]-rz*ev[3]);
    ev[31]=PA[2]*ev[13]+WP[2]*ev[16]+zeta2*(ev[2]-rz*ev[3]);
    ev[32]=PA[0]*ev[12]+WP[0]*ev[15];
    ev[33]=PA[0]*ev[13]+WP[0]*ev[16];
    ev[34]=PA[1]*ev[13]+WP[1]*ev[16];
    // (ss|ps) m=[1,1]
    ev[35]=QC[0]*ev[1]+WQ[0]*ev[2];
    ev[36]=QC[1]*ev[1]+WQ[1]*ev[2];
    ev[37]=QC[2]*ev[1]+WQ[2]*ev[2];
    // (ps|ps) m=[0,1]
    ev[38]=QC[0]*ev[ 5]+WQ[0]*ev[ 8]+ze2*ev[1];
    ev[39]=QC[1]*ev[ 5]+WQ[1]*ev[ 8];
    ev[40]=QC[2]*ev[ 5]+WQ[2]*ev[ 8];
    ev[41]=QC[0]*ev[ 6]+WQ[0]*ev[ 9];
    ev[42]=QC[1]*ev[ 6]+WQ[1]*ev[ 9]+ze2*ev[1];
    ev[43]=QC[2]*ev[ 6]+WQ[2]*ev[ 9];
    ev[44]=QC[0]*ev[ 7]+WQ[0]*ev[10];
    ev[45]=QC[1]*ev[ 7]+WQ[1]*ev[10];
    ev[46]=QC[2]*ev[ 7]+WQ[2]*ev[10]+ze2*ev[1];
    ev[47]=QC[0]*ev[ 8]+WQ[0]*ev[11]+ze2*ev[2];
    ev[48]=QC[1]*ev[ 8]+WQ[1]*ev[11];
    ev[49]=QC[2]*ev[ 8]+WQ[2]*ev[11];
    ev[50]=QC[0]*ev[ 9]+WQ[0]*ev[12];
    ev[51]=QC[1]*ev[ 9]+WQ[1]*ev[12]+ze2*ev[2];
    ev[52]=QC[2]*ev[ 9]+WQ[2]*ev[12];
    ev[53]=QC[0]*ev[10]+WQ[0]*ev[13];
    ev[54]=QC[1]*ev[10]+WQ[1]*ev[13];
    ev[55]=QC[2]*ev[10]+WQ[2]*ev[13]+ze2*ev[2];
    // (ds|ps) m=[0,1]
    ev[56]=QC[0]*ev[17]+WQ[0]*ev[23]+2.e0*ze2*ev[ 8];
    ev[57]=QC[1]*ev[17]+WQ[1]*ev[23];
    ev[58]=QC[2]*ev[17]+WQ[2]*ev[23];
    ev[59]=QC[0]*ev[18]+WQ[0]*ev[24];
    ev[60]=QC[1]*ev[18]+WQ[1]*ev[24]+2.e0*ze2*ev[ 9];
    ev[61]=QC[2]*ev[18]+WQ[2]*ev[24];
    ev[62]=QC[0]*ev[19]+WQ[0]*ev[25];
    ev[63]=QC[1]*ev[19]+WQ[1]*ev[25];
    ev[64]=QC[2]*ev[19]+WQ[2]*ev[25]+2.e0*ze2*ev[10];
    ev[65]=QC[0]*ev[20]+WQ[0]*ev[26]+     ze2*ev[ 9];
    ev[66]=QC[1]*ev[20]+WQ[1]*ev[26]+     ze2*ev[ 8];
    ev[67]=QC[2]*ev[20]+WQ[2]*ev[26];
    ev[68]=QC[0]*ev[21]+WQ[0]*ev[27]+     ze2*ev[10];
    ev[69]=QC[1]*ev[21]+WQ[1]*ev[27];
    ev[70]=QC[2]*ev[21]+WQ[2]*ev[27]+     ze2*ev[ 8];
    ev[71]=QC[0]*ev[22]+WQ[0]*ev[28];
    ev[72]=QC[1]*ev[22]+WQ[1]*ev[28]+     ze2*ev[10];
    ev[73]=QC[2]*ev[22]+WQ[2]*ev[28]+     ze2*ev[ 9];
    ev[74]=QC[0]*ev[23]+WQ[0]*ev[29]+2.e0*ze2*ev[11];
    ev[75]=QC[1]*ev[23]+WQ[1]*ev[29];
    ev[76]=QC[2]*ev[23]+WQ[2]*ev[29];
    ev[77]=QC[0]*ev[24]+WQ[0]*ev[30];
    ev[78]=QC[1]*ev[24]+WQ[1]*ev[30]+2.e0*ze2*ev[12];
    ev[79]=QC[2]*ev[24]+WQ[2]*ev[30];
    ev[80]=QC[0]*ev[25]+WQ[0]*ev[31];
    ev[81]=QC[1]*ev[25]+WQ[1]*ev[31];
    ev[82]=QC[2]*ev[25]+WQ[2]*ev[31]+2.e0*ze2*ev[13];
    ev[83]=QC[0]*ev[26]+WQ[0]*ev[32]+     ze2*ev[12];
    ev[84]=QC[1]*ev[26]+WQ[1]*ev[32]+     ze2*ev[11];
    ev[85]=QC[2]*ev[26]+WQ[2]*ev[32];
    ev[86]=QC[0]*ev[27]+WQ[0]*ev[33]+     ze2*ev[13];
    ev[87]=QC[1]*ev[27]+WQ[1]*ev[33];
    ev[88]=QC[2]*ev[27]+WQ[2]*ev[33]+     ze2*ev[11];
    ev[89]=QC[0]*ev[28]+WQ[0]*ev[34];
    ev[90]=QC[1]*ev[28]+WQ[1]*ev[34]+     ze2*ev[13];
    ev[91]=QC[2]*ev[28]+WQ[2]*ev[34]+     ze2*ev[12];
    // (ps|ds) m=[0,0]
    ev[ 92]=QC[0]*ev[38]+WQ[0]*ev[47]+eta2*(ev[5]-re*ev[ 8])+ze2*ev[35];
    ev[ 93]=QC[1]*ev[39]+WQ[1]*ev[48]+eta2*(ev[5]-re*ev[ 8]);
    ev[ 94]=QC[2]*ev[40]+WQ[2]*ev[49]+eta2*(ev[5]-re*ev[ 8]);
    ev[ 95]=QC[0]*ev[39]+WQ[0]*ev[48]+ze2*ev[36];
    ev[ 96]=QC[0]*ev[40]+WQ[0]*ev[49]+ze2*ev[37];
    ev[ 97]=QC[1]*ev[40]+WQ[1]*ev[49];
    ev[ 98]=QC[0]*ev[41]+WQ[0]*ev[50]+eta2*(ev[6]-re*ev[ 9]);
    ev[ 99]=QC[1]*ev[42]+WQ[1]*ev[51]+eta2*(ev[6]-re*ev[ 9])+ze2*ev[36];
    ev[100]=QC[2]*ev[43]+WQ[2]*ev[52]+eta2*(ev[6]-re*ev[ 9]);
    ev[101]=QC[0]*ev[42]+WQ[0]*ev[51];
    ev[102]=QC[0]*ev[43]+WQ[0]*ev[52];
    ev[103]=QC[1]*ev[43]+WQ[1]*ev[52]+ze2*ev[37];
    ev[104]=QC[0]*ev[44]+WQ[0]*ev[53]+eta2*(ev[7]-re*ev[10]);
    ev[105]=QC[1]*ev[45]+WQ[1]*ev[54]+eta2*(ev[7]-re*ev[10]);
    ev[106]=QC[2]*ev[46]+WQ[2]*ev[55]+eta2*(ev[7]-re*ev[10])+ze2*ev[37];
    ev[107]=QC[0]*ev[45]+WQ[0]*ev[54];
    ev[108]=QC[0]*ev[46]+WQ[0]*ev[55];
    ev[109]=QC[1]*ev[46]+WQ[1]*ev[55];
    // (ds|ds) m=[0,0]
    ev[110]=QC[0]*ev[56]+WQ[0]*ev[74]+eta2*(ev[17]-re*ev[23])
            +2.e0*ze2*ev[47];
    ev[111]=QC[1]*ev[57]+WQ[1]*ev[75]+eta2*(ev[17]-re*ev[23]);
    ev[112]=QC[2]*ev[58]+WQ[2]*ev[76]+eta2*(ev[17]-re*ev[23]);
    ev[113]=QC[0]*ev[57]+WQ[0]*ev[75]+2.e0*ze2*ev[48];
    ev[114]=QC[0]*ev[58]+WQ[0]*ev[76]+2.e0*ze2*ev[49];
    ev[115]=QC[1]*ev[58]+WQ[1]*ev[76];
    ev[116]=QC[0]*ev[59]+WQ[0]*ev[77]+eta2*(ev[18]-re*ev[24]);
    ev[117]=QC[1]*ev[60]+WQ[1]*ev[78]+eta2*(ev[18]-re*ev[24])
            +2.e0*ze2*ev[51];
    ev[118]=QC[2]*ev[61]+WQ[2]*ev[79]+eta2*(ev[18]-re*ev[24]);
    ev[119]=QC[0]*ev[60]+WQ[0]*ev[78];
    ev[120]=QC[0]*ev[61]+WQ[0]*ev[79];
    ev[121]=QC[1]*ev[61]+WQ[1]*ev[79]+2.e0*ze2*ev[52];
    ev[122]=QC[0]*ev[62]+WQ[0]*ev[80]+eta2*(ev[19]-re*ev[25]);
    ev[123]=QC[1]*ev[63]+WQ[1]*ev[81]+eta2*(ev[19]-re*ev[25]);
    ev[124]=QC[2]*ev[64]+WQ[2]*ev[82]+eta2*(ev[19]-re*ev[25])
            +2.e0*ze2*ev[55];
    ev[125]=QC[0]*ev[63]+WQ[0]*ev[81];
    ev[126]=QC[0]*ev[64]+WQ[0]*ev[82];
    ev[127]=QC[1]*ev[64]+WQ[1]*ev[82];
    ev[128]=QC[0]*ev[65]+WQ[0]*ev[83]+eta2*(ev[20]-re*ev[26])
            +     ze2*ev[50];
    ev[129]=QC[1]*ev[66]+WQ[1]*ev[84]+eta2*(ev[20]-re*ev[26])
            +     ze2*ev[48];
    ev[130]=QC[2]*ev[67]+WQ[2]*ev[85]+eta2*(ev[20]-re*ev[26]);
    ev[131]=QC[0]*ev[66]+WQ[0]*ev[84]+     ze2*ev[51];
    ev[132]=QC[0]*ev[67]+WQ[0]*ev[85]+     ze2*ev[52];
    ev[133]=QC[1]*ev[67]+WQ[1]*ev[85]+     ze2*ev[49];
    ev[134]=QC[0]*ev[68]+WQ[0]*ev[86]+eta2*(ev[21]-re*ev[27])
            +     ze2*ev[53];
    ev[135]=QC[1]*ev[69]+WQ[1]*ev[87]+eta2*(ev[21]-re*ev[27]);
    ev[136]=QC[2]*ev[70]+WQ[2]*ev[88]+eta2*(ev[21]-re*ev[27])
            +     ze2*ev[49];
    ev[137]=QC[0]*ev[69]+WQ[0]*ev[87]+     ze2*ev[54];
    ev[138]=QC[0]*ev[70]+WQ[0]*ev[88]+     ze2*ev[55];
    ev[139]=QC[1]*ev[70]+WQ[1]*ev[88];
    ev[140]=QC[0]*ev[71]+WQ[0]*ev[89]+eta2*(ev[22]-re*ev[28]);
    ev[141]=QC[1]*ev[72]+WQ[1]*ev[90]+eta2*(ev[22]-re*ev[28])
            +     ze2*ev[54];
    ev[142]=QC[2]*ev[73]+WQ[2]*ev[91]+eta2*(ev[22]-re*ev[28])
            +     ze2*ev[52];
    ev[143]=QC[0]*ev[72]+WQ[0]*ev[90];
    ev[144]=QC[0]*ev[73]+WQ[0]*ev[91];
    ev[145]=QC[1]*ev[73]+WQ[1]*ev[91]+     ze2*ev[55];
}

static void ofmo_vrr_cint_pppp( const double *ev, double *eh ) {
    int La=1, Lb=1, Lc=1, Ld=1;
    int i, ih, iv;
    // (PS|PS)
    for ( i=0, iv=38, ih=0; i<9; i++, iv++, ih++ ) eh[ih]+=ev[iv];
    // (PS|DS)
    for ( i=0, iv=92, ih=9; i<18; i++, iv++, ih++ ) eh[ih]+=ev[iv];
    // (DS|PS)
    for ( i=0, iv=56, ih=27; i<18; i++, iv++, ih++ ) eh[ih]+=ev[iv];
    // (DS|DS)
    for ( i=0, iv=110, ih=45; i<36; i++, iv++, ih++ ) eh[ih]+=ev[iv];
}

void ofmo_twoint_core_os_pppp(
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
    double ev[146], eh[243];
    int La=*pLa, Lb=*pLb, Lc=*pLc, Ld=*pLd;

    ofmo_hrr_clear_pppp( eh );
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
            ofmo_vrr_calc_pppp(
                    T, cssss, zeta2, eta2, ze2, rz, re, PA, WP, QC, WQ,
                    ev );
            ofmo_vrr_cint_pppp( ev, eh );
        }	// for (klps)
    }	// for (ijps)
    ofmo_hrr_calc_pppp( BA, DC, eh );
    ofmo_hrr_coef_pppp( eh, DINT );
}

int ofmo_twoint_os_pppp(
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
    double DINTEG[3*3*3*3];
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
    max_nzeri = ebuf_max_nzeri - 3*3*3*3;
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
            ofmo_twoint_core_os_pppp(
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
            for ( i=0, iao=iao0, ix=0; i<3; i++, iao++ ) {
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

int ofmo_twoint_direct_os_pppp(
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
    double DINTEG[3*3*3*3];
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
    
    max_nzeri -= 3*3*3*3;
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
            ofmo_twoint_core_os_pppp(
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
            for ( i=0, iao=iao0, ix=0; i<3; i++, iao++ ) {
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
