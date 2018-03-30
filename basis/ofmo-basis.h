/**
 * @file ofmo-basis.h
 *
 * 非ソート、ソート基底関数データを取得するための
 * 関数群のプロトタイプ宣言を行っているヘッダファイル
 * */
#ifndef _OFMO_BASIS_H_
#define _OFMO_BASIS_H_

extern int ofmo_get_basis_size( const int nat, const int nsbs,
	char **basis_name, const int atomic_number[],
	const int atom_basis[],
	int *maxlqn, int *ncs, int *nao, int *nps );

extern int ofmo_assign_basis(
	const int nat, const int nsbs, char *basis_name[],
	const int atomic_number[], const int atom_basis[],
	int ushel_lqn[], int ushel_tem[], int ushel_atm[],
	int ushel_add[], int ushel_ini[],
	double uprim_exp[], double uprim_coe[] );

extern int ofmo_sort_basis(
	const int maxlqn, const int ncs,
	const int ushel_lqn[], const int ushel_tem[],
	const int ushel_atm[], const int ushel_add[],
	const int ushel_ini[],
	const double uprim_exp[], const double uprim_coe[],
	int leading_cs[],
	int shel_tem[], int shel_atm[], int shel_add[], int shel_ini[],
	double prim_exp[], double prim_coe[], int s2u[]);

extern int ofmo_alloc_unsorted_basis( const int ncs, const int nps,
	int **ushel_lqn, int **ushel_tem, int **ushel_atm,
	int **ushel_add, int **ushel_ini,
	double **uprim_exp, double **uprim_coe );

extern int ofmo_alloc_sorted_basis(
	const int maxlqn, const int ncs, const int nao, const int nps,
	int **leading_cs, int **shel_tem, int **shel_atm, int **shel_add,
	int **shel_ini, double **prim_exp, double **prim_coe,
	int **s2u );

extern void ofmo_show_unsorted_basis_params( FILE *fp, const int ncs,
	const int ushel_lqn[], const int ushel_tem[],
	const int ushel_atm[], const int ushel_add[],
	const int ushel_ini[],
	const double uprim_exp[], const double uprim_coe[] );

extern void ofmo_show_sorted_basis_params( FILE *fp,
	const int maxlqn, const int nao, const int leading_cs[],
	const int shel_tem[], const int shel_atm[],
	const int shel_add[], const int shel_ini[],
	const double prim_exp[], const double prim_coe[],
	const int s2u[] );

extern void ofmo_unsort_packed_matrix(const int nao, const int s2u[],
	const double src[], double dist[]);

extern void ofmo_sort_packed_matrix(const int nao, const int s2u[],
	const double src[], double dest[]);

extern void ofmo_sort_vector( const int nao, const int sao2uao[],
	const double v0[], double vs[] );

#endif
