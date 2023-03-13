/**
 * @file ofmo-basis-database.h
 *
 * 基底関数の元データを管理する関数群のプロトタイプ宣言を
 * 行っているファイル
 *
 * */
#ifndef _FMO_BASIS_DATABASE_H_
#define _FMO_BASIS_DATABASE_H_

#include <stdarg.h>
#include "ofmo-def.h"

extern void ofmo_show_available_basis();
extern int ofmo_get_basis_data( const char basis_name[],
	char* format, ...);

#endif
