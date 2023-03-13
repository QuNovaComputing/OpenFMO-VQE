#ifndef _OFMO_DATA_H_
#define _OFMO_DATA_H_

#include <stdio.h>
#include <stdarg.h>

extern int ofmo_data_show_all();
extern int ofmo_data_get_vals( const char* format, ...);
extern int ofmo_data_put_vals( const char* format, ... );

#endif
