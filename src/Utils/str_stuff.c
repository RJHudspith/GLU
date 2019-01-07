/**
   @file str_stuff.c
   @brief common functions that involve chars
 */
#include "Mainfile.h"

// append a char
void
append_char( char **str ,
	     const char *tmp )
{
  *str = realloc( *str , (1+strlen(tmp)+strlen(*str))*sizeof(char) ) ;
  *str = strcat( *str , tmp ) ;
}

// strcmp defaults to 0 if they are equal which is contrary to standard if statements
int
are_equal( const char *str_1 , const char *str_2 ) {
  return !strcmp( str_1 , str_2 ) ;
}
