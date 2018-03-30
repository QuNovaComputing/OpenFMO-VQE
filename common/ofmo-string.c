#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "ofmo-def.h"

static char *data[] = { "x",
    "h", "he", "li", "be", "b", "c", "n", "o", "f", "ne", 
    "na", "mg", "al", "si", "p", "s", "cl", "ar", "k", "ca", 
    "sc", "ti", "v", "cr", "mn", "fe", "co", "ni", "cu", "zn", 
    "ga", "ge", "as", "se", "br", "kr", "rb", "sr", "y", "zr", 
    "nb", "mo", "tc", "ru", "rh", "pd", "ag", "cd", "in", "sn", 
    "sb", "te", "i", "xe", "cs", "ba", "la", "ce", "pr", "nd", 
    "pm", "sm", "eu", "gd", "tb", "dy", "ho", "er", "tm", "yb", 
    "lu", "hf", "ta", "w", "re", "os", "ir", "pt", "au", "hg", 
    "tl", "pb", "bi", "po", "at", "rn", "fr", "ra", "ac", "th", 
    "pa", "u", "np", "pu", "am", "cm", "bk", "cf", "es", "fm", 
    "md", "no", "lr"
};

// local variables used in functions defined in this file
static char strtmp[MAXSTRLEN*MAXTOKEN];

/* 2次元文字列配列を確保する（トークンへの分割で用いるため） */
char** ofmo_alloc_char_2d( const int ntok, const int maxtoklen ) {
    int i;
    char **s;
    // check
    if ( ntok > MAXTOKEN ) {
	dbg("Illegal number of token (%d):"
		" must be less equal MAXTOKEN(%d)\n" , ntok, MAXTOKEN );
	return NULL;
    }
    if ( maxtoklen > MAXSTRLEN ) {
	dbg("Illegal length of token (%d):"
		"must be less equal MAXSTRLEN (%d)\n",
		maxtoklen, MAXSTRLEN );
	return NULL;
    }
    s = (char**)malloc( sizeof(char*) * ntok );
    s[0] = (char*)malloc( sizeof(char) * ntok * maxtoklen );
    for ( i=1; i<ntok; i++ ) s[i] = s[i-1] + maxtoklen;
    for ( i=0; i<ntok; i++ ) s[i][0] = '\0';
    return s;
}

/* 2次元文字列配列を解放する */
void ofmo_free_char_2d( char **s ) {
    if ( s != NULL ) {
	if ( s[0] != NULL ) free( s[0] );
	free( s );
    }
}

/* 文字列配列をコピー */
void ofmo_strcpy_2d( char** src, char** dest,
	const int ntok, const int maxtoklen ) {
    int i, len;
    /*// debug
    printf("src=%p, dest=%p, ntok=%d, maxtoklen=%d\n",
	    src[0], dest[0], ntok, maxtoklen );
    fflush(stdout);*/

    for ( i=0; i<ntok; i++ ) {
	len = strlen( src[i] );
	if ( len >= (maxtoklen-1) ) len = maxtoklen-1;
	strncpy( dest[i], src[i], len );
	dest[i][len] = '\0';
    }
}

/* 与えられた文字列を小文字に変換する */
void ofmo_tolower( char s[] ) {
    int i, n;
    n = strlen( s );
    for ( i=0; i<n; i++ ) s[i] = tolower( s[i] );
}

/* 文字列の先頭、及び、末尾の表示不可能文字（スペースを含む）を
 * 削除する */
void ofmo_delete_prepost_space( char s[] ) {
    int n, c, start, end;
    strcpy( strtmp, s );
    n = strlen( s );
    // search start & end points
    for ( start=0; start<n; start++ ) if ( isgraph(strtmp[start]) ) break;
    for ( end=n-1; end>=0; end-- ) if ( isgraph(strtmp[end]) ) break;
    c = end - start + 1;
    if ( c > 0 ) {
	strncpy( s, &strtmp[start], c );
	s[c] = '\0';
    } else {
	s[0] = '\0';
    }
}

/*
  関数名：split
  機能  ：文字列をトークンに分割する
  引数：
   const char* s = (in)分割する文字列
   const char* sep = (in)区切り文字
   char **token = (out)トークンが代入される配列
	token[maxtoken][maxtoklen]
   int maxtoken = (in)トークンの最大値
   int maxtoklen = (in)トークンあたりの最大文字数
  戻り値：
   トークン数
*/
int ofmo_split( const char s[], const char sep[], char **token,
	const int maxtoken, const int maxtoklen ) {
    char *p;
    int count=0, len;
    strcpy( strtmp, s );
    // トークンの切り分け
    p = strtok( strtmp, sep );
    if ( p == NULL ) return 0;
    len = strlen( p );
    if ( len >= maxtoklen ) len = maxtoklen - 1;
    strncpy( token[count], p, len );
    token[count][len] = '\0';
    count++;
    while ( (p=strtok( NULL, sep )) != NULL && count < maxtoken ) {
	len = strlen( p );
	if ( len >= maxtoklen ) len = maxtoklen - 1;
	strncpy( token[count], p, len );
	token[count][len] = '\0';
	count++;
    }
    return count;
}

/* 文字列のすべての表示不可能文字（スペースを含む）を削除する */
void ofmo_delete_all_space( char s[] ) {
    int n, c, i;
    strcpy( strtmp, s );
    n = strlen( s );
    c=0;
    for ( i=0; i<n; i++ ) {
	if ( isgraph(strtmp[i]) ) {
	    s[c] = strtmp[i];
	    c++;
	}
    }
    s[c] = '\0';
}

void ofmo_delete_multiple_space( char s[] ) {
    int n, c, i;
    int flag = 0;
    strcpy( strtmp, s );
    n = strlen( s );
    c=0;
    for ( i=0; i<n; i++ ) {
	if ( isgraph(strtmp[i]) ) {
	    s[c] = strtmp[i];
	    c++;
            flag=0;
	} else {
          if (flag==0) {
	    s[c] = ' ';
	    c++;
          }
          flag=1;
        }
    }
    s[c] = '\0';
}

/* 与えられた文字列を大文字に変換する */
void ofmo_toupper( char s[] ) {
    int i, n;
    n = strlen( s );
    for ( i=0; i<n; i++ ) s[i] = toupper( s[i] );
}

// 与えられた文字列の特定の文字を別の文字に置き換える
void ofmo_replace_char( char s[], const char org, const char rep ) {
    int i, n;
    n = strlen( s );
    for ( i=0; i<n; i++ ) if ( s[i] == org ) s[i] = rep;
}

/** 文字列からアルファベットのみを抜き出す
  戻り値：アルファベット以外の文字数
 * */
int ofmo_extract_alpha( char s[] ) {
    int i,n, c, k;
    char tmp[MAXSTRLEN];
    n = strlen( s );
    strcpy( tmp, s );
    for ( i=0, c=0, k=0; i<n; i++ ) {
	if ( isalpha( tmp[i] ) ) s[c++] = tmp[i];
	else                   k++;
    }
    s[c] = '\0';
    return k;
}

void ofmo_delete_multiple_space_and_comma_space( char s[] )
{
    int n, c, i;
    int flag = 0;
    strcpy( strtmp, s );
    n = strlen( s );
    c=0;
    for ( i=0; i<n; i++ ) {
	if ( isgraph(strtmp[i]) ) {
	    s[c] = strtmp[i];
	    c++;
            flag=0;
            if ( strtmp[i]==',' ) flag=1;
	} else {
          if (flag==0) {
	    s[c] = ' ';
	    c++;
          }
          flag=1;
        }
    }
    s[c] = '\0';
}

/** １行読み込んで、小文字に変換したのち、空白文字で分割する
  */
int ofmo_read_line( FILE *fp, const int maxstrlen,
	const int maxtoken, const int maxtoklen, char line[],
	char **token ) {
    if ( fgets( line, maxstrlen, fp ) == NULL ) return -1;
    ofmo_delete_prepost_space( line );
    ofmo_tolower( line );
    //ofmo_delete_multiple_space_and_comma_space( line )
    return ofmo_split( line, " \t\n\r", token, maxtoken, maxtoklen );
}

/* 与えられた文字列（元素記号）に対応する原子番号を返す */
int AtomicNumber( const char *name ) {
    int maxatn = MAXATOMICNUMBER;
    int i;
    char tmp[MAXSTRLEN];

    if (strlen(name) >= MAXSTRLEN || strlen(name) < 1) return -2;
    strcpy( tmp, name );
    ofmo_extract_alpha( tmp );
    ofmo_tolower( tmp );
    
    for ( i=0; i<=maxatn; i++ ) {
	if ( strcmp(tmp, data[i]) == 0 ) return i;
    }
    dbg("Illegal atomic symbol (%s)\n", tmp);
    return -1;
}

char* AtomicSymbol( const int atn )
{
  int atom=0;
  if (atn>0 && atn<=MAXATOMICNUMBER) atom=atn;
  strncpy(strtmp, data[atn], 3);
  strtmp[0]=toupper(strtmp[0]);
  return strtmp;
}

/* AOのタイプを表す文字列を返す */
/* 作成途中 */
/*static char ***AO_TYPE = NULL;
static int _MAXLQN_;
static int ofmo_ao_type_init() {
    int nnao[] = { 1, 3, 6, 10, 15 };
    int tot_clen, tot_tlen, ic, it;
    _MAXLQN_ = sizeof(nnao)/sizeof(int) - 1;
    tot_clen = tot_tlen = 0;
    for ( int i=0; i<=_MAXLQN_; i++ ) {
	tot_clen += (nnao[i]*(i+1));
	tot_tlen += nnao[i];
    }
    AO_TYPE = (char***)malloc( sizeof(char**) * (_MAXLQN_+1) );
    AT_TYPE[0] = (char**)malloc( sizeof(char*) * tot_tlen );
    AT_TYPE[0][0] = (char*)malloc( sizeof(char) * tot_clen );
    for ( int lqn=1; lqn<=_MAXLQN_; lqn++ )
	AT_TYPE[lqn] = AT_TYPE[lqn-1] + nnao[lqn-1];
    ic = it = 0;
    for ( int lqn=0; lqn<=_MAXLQN_; lqn++ ) {
	for ( int iao=0; iao<nnao[lqn]; iao++ )
}*/
