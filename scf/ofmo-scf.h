/**
 * @file ofmo-scf.h
 *
 * RHF-SCF計算のサンプル関数のプロトタイプ宣言を行っている
 * ヘッダファイル
 *
 * */
#ifndef _OFMO_SCF_H_
#define _OFMO_SCF_H_

enum ofmo_scf_convType { scf0, scc, scf };

void ofmo_scf_set_convType(enum ofmo_scf_convType type);
void ofmo_scf_set_minimum_scf(int mincyc);

extern int ofmo_scf_init( const int nao );

extern double ofmo_calc_nuclear_repulsion(
	const int nat, const int atomic_number[],
	const double atom_x[], const double atom_y[],
	const double atom_z[] );

extern int ofmo_scf_make_rhf_density(const int nao, const int nocc,
	const double C[], double D[]);

extern double ofmo_scf_rhf_energy(const int nao, const double D[],
	const double H[], const double F[]);

extern int ofmo_scf_init_density_ehuckel(
	const int nat, const int ncs, const int nao, const int maxlqn,
	const int nocc, const int atomic_number[],
	const int leading_cs[], const int shel_atm[], const int shel_ini[],
	const double SP[], double DP[], double aopop[], double atpop[] );

extern int ofmo_scf_mulliken_population(
	const int nat, const int nao, const int maxlqn,
	const int leading_cs[], const int shel_atm[], const int shel_ini[],
	const double SP[], const double DP[],
	double aopop[], double atpop[] );

extern int ofmo_scf_rhf(
	MPI_Comm comm, const int maxlqn, const double Enuc,
	const int ncs, const int nao,
	const int leading_cs[],
	// 基底関数データ
	const int shel_atm[], const int shel_ini[],
	const double atom_x[], const double atom_y[],
	const double atom_z[],
	// カットオフテーブルデータ
	const int leading_cs_pair[], const double csp_schwarz[],
	const int csp_ics[], const int csp_jcs[],
	const int csp_leading_ps_pair[],
	const double psp_zeta[], const double psp_dkps[],
	const double psp_xiza[],
	// 分子データ
	const int nat,
	const int nocc,
	// 積分データ
	double S[], double H[],
	// 制御用データ
	const int maxscfcyc, const double scfe, const double scfd,
	// 結果代入用データ
	double D[], double C[], double moe[], double *Eelec);
#endif
