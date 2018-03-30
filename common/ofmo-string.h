#ifndef _OFMO_STRING_H_
#define _OFMO_STRING_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern char** ofmo_alloc_char_2d( const int ntok, const int maxtoklen );
extern void ofmo_free_char_2d( char **s );

extern void ofmo_strcpy_2d( char** src, char** dest,
	const int ntok, const int maxtoklen );

extern void ofmo_tolower( char s[] );
extern void ofmo_delete_prepost_space( char s[] );
extern int ofmo_split( const char s[], const char sep[], char **token,
	const int maxtoken, const int maxtoklen );
extern void ofmo_delete_all_space( char s[] );
extern void ofmo_delete_multiple_space( char s[] );

extern void ofmo_toupper( char s[] );
extern void ofmo_replace_char( char s[], const char org, const char rep );
extern int ofmo_extract_alpha(char s[] );
extern int AtomicNumber( const char name[] );
extern char* AtomicSymbol( const int atn );

extern int ofmo_read_line( FILE *fp, const int maxstrlen,
	const int maxtoken, const int maxtoklen, char *line,
	char **token );

#endif
