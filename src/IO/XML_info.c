/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (XML_info.c) is part of GLU.

    GLU is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GLU is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GLU.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   @file XML_info.c
   @brief dirty parser for the xml information in Scidac and ILDG config files
 */
#include "Mainfile.h"

#include <errno.h>

#include "str_stuff.h" // are_equal

/**
  These functions came from ahmidas, they claim to have had trouble parsing ILDG
  xml data due to some poor choices by ETMC, just to be safe I include them as
  my code for understanding the xml is similar to theirs
 */
inline static char *realFront(char *string) { return string + strspn(string, " \t"); }

// returns 0 if there is nothing else returns the value we
// are looking for
static int
get_int_tag( char *pch , const char *tag )
{
  const int length = strlen( tag ) ;
  char tmp[ length + 2 ] , *endptr ;
  sprintf( tmp , "/%s" , tag ) ;
  // loop tokens
  int value = 0 ;
  errno = 0 ;
  while( ( pch = strtok( 0 , "<>" ) ) ) {
    if( !strncmp( pch , tmp , length+1 ) ) break ;
    value = (int)strtol( realFront( pch ) , &endptr , 10 ) ;
    if( errno == ERANGE || endptr == realFront( pch ) ) {
      return GLU_FAILURE ;
    }
  }
  return value ;
}

// the result is in pch
static int
get_prec_tag( char *pch , const char *tag )
{
  const int length = strlen( tag ) ;
  char tmp[ length + 2 ] , *endptr ;
  sprintf( tmp , "/%s" , tag ) ;
  char prec[ 8 ] = "00" ;
  // loop tokens
  while( ( pch = strtok( 0 , "<>" ) ) ) {
    if( !strncmp( pch , tmp , length+1 ) ) break ;
    sprintf( prec , "%s" , pch ) ;
  }
  if( !strcmp( prec , "F" ) ) return FLOAT_PREC ;
  if( !strcmp( prec , "D" ) ) return DOUBLE_PREC ;
  const int sd = (int)strtol( realFront( prec ) , &endptr , 10 ) ;
  if( sd == 32 ) return FLOAT_PREC ;
  if( sd == 64 ) return DOUBLE_PREC ;
  return GLU_FAILURE ;
}

// parse the information from the c-str of xml info
int
parse_and_set_xml_SCIDAC( char *xml_info ,
			  struct head_data *HEAD_DATA ) 
{
  // scidac is always big endian, which is nice
  HEAD_DATA -> endianess = B_ENDIAN ;

  // We use the C tokenize capabilities to parse this string
  char *pch = strtok( xml_info , "<>" ) ;
  // this is not an error
  if (strncmp( pch , "?xml" , 4 ) ) {
    return 0 ;
  }
  int dimensions = 0 ;

  // set up the search for extra dimensions, extra space for terminating null
  const char open[4][3] = { "lx" , "ly" , "lz" , "lt" } ;
  char search[ND][3] ;
  size_t mu ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    if( mu < 4 ) {
      sprintf( search[mu] , "%s" , open[mu] ) ;
    } else {
      sprintf( search[mu] , "l%zu" , mu ) ;
    }
  }
  sprintf( search[ND-1] , "%s" , open[3] ) ;

  // We've removed the XML header, now we can set up a state 
  // machine to parse the file
  while( ( pch = strtok( 0 , "<>" ) ) ) {
    while( ( pch = strtok( 0 , "<>" ) ) ) {
      
      // break if we are at the end of a scidac file or record
      if( are_equal( pch , "/scidacRecord" ) ) break ;
      if( are_equal( pch , "/scidacFile" ) ) break ;
      if( are_equal( pch , "/ildgFormat" ) ) break ;

      // get number of spacetime dimensions
      if( are_equal( pch , "spacetime" ) ) {
	if( ( dimensions = get_int_tag( pch , "spacetime" ) ) != 0 ) {
	  if( dimensions != ND ) {
	    fprintf( stderr , "[IO] ND mismatch compiled %d vs. read %d \n" ,
		     ND , dimensions ) ;
	    return GLU_FAILURE ;
	  }
	}
        #ifdef DEBUG
	printf( "[IO] spacetime :: %d \n" , dimensions ) ;
        #endif
      }

      // have a look at the colors
      if( are_equal( pch , "colors" ) ) {
	if( ( dimensions = get_int_tag( pch , "colors" ) ) != 0 ) {
	  if( dimensions != NC ) {
	    printf( "[IO] NC mismatch compiled %d vs. read %d \n" ,
		    NC , dimensions ) ;
	    return GLU_FAILURE ;
	  }
	}
      }

      // have a look at the storage type size
      if( are_equal( pch , "typesize" ) ) {
	if( ( dimensions = get_int_tag( pch , "typesize" ) ) != 0 ) {
	  HEAD_DATA -> config_type = OUTPUT_NCxNC ; // binary file output
	  if( dimensions == 2 * ND * NCNC ) {
	    HEAD_DATA -> precision = FLOAT_PREC ;
	  } else if( dimensions == ( 4 * ND * NCNC ) ) {
	    HEAD_DATA -> precision = DOUBLE_PREC ;
	  } else {
	    fprintf( stderr , "[IO] Storage type not recognised \n" ) ;
	    return GLU_FAILURE ;
	  }
	}
      }

      // get the floating point output precision ...
      if( are_equal( pch , "precision" ) ) {
	if( ( dimensions = get_prec_tag( pch , "precision" ) ) != -1 ) {
	  if( dimensions == FLOAT_PREC ) {
	    #ifdef verbose
	    fprintf( stdout , "[IO] Attempting to read "
		     "single precision data\n" ) ;
	    #endif
	    HEAD_DATA -> precision = FLOAT_PREC ;
	  } else if( dimensions == DOUBLE_PREC ) {
	    #ifdef verbose
	    fprintf( stdout , "[IO] Attempting to read double "
		     "precision data\n" ) ;
	    #endif
	    HEAD_DATA -> precision = DOUBLE_PREC ;
	  } else {
	    fprintf( stderr , "[IO] Precision %d not understood \n" , 
		     dimensions ) ;
	    return GLU_FAILURE ;
	  }
	}
      }

      // grok the lattice dimensions
      if ( are_equal( pch, "dims" ) ) {
	while ( ( pch = strtok( 0 , "<>" ) ) ) {
	  // break up the dimensions ...
	  char *token = strtok( pch , " " ) ;
	  if( token == NULL ) return GLU_FAILURE ;
	  int idx = 0 ;
	  char *pEnd ;
	  Latt.dims[ idx++ ] = strtol( token , &pEnd , 10 ) ;
	  while( ( token = strtok( NULL , " " ) ) != NULL ) {
	    if(  ( are_equal( token , "/dims" ) ) ) break ;
	    Latt.dims[ idx ] = strtol( token , &pEnd , 10 ) ;
	    idx++ ;
	  } 
	  // and break the xml
	  break;
	}
	continue ;
      }

      // I just take the first checksum ...
      if( are_equal( pch , "suma" ) ) {
	while( ( pch = strtok( 0 , "<>" ) ) ) {
	  if( are_equal( pch , "/suma" ) ) break ;
	  sscanf( pch , "%x" , &(HEAD_DATA -> checksum) ) ;
	}
	continue ;
      }

      ////////////////////////  ILDG SPECIFIC TAGS ////////////////////
      // If I were to make up a format it would be exactly the       //
      // same as another one apart from tiny, pointless, differences //

      // geometry, not sure about making this ND-generic much prefer the
      // scidac "dims" array
      int length = 0 ;
      for( mu = 0 ; mu < ND ; mu++ ) {
	if( are_equal( pch , search[mu] ) ) {
	  if( ( length = get_int_tag( pch , search[mu] ) ) != 0 ) {
	    Latt.dims[ mu ] = length ;
	  }
	  continue ;
	}
      }

      // for getting the field ... su3gauge is actually NCxNC by the look of it
      // loop tokens
      if( are_equal( pch , "field" ) ) {
	while( ( pch = strtok( 0 , "<>" ) ) ) {
	  if( are_equal( pch , "/field" ) ) break ;
	  // sometimes there is whitespace around these for some unknown reason
	  char *token = realFront( pch ) ;
	  char compare[10] ;
	  sprintf( compare , "su%dgauge" , NC ) ;
	  if( !strncmp( compare , token , 8 ) ) {
	    HEAD_DATA -> config_type = OUTPUT_NCxNC ;
	  } else {
	    fprintf( stderr , "[ILDG] Expected %s, got \"%s\" \n" , 
		     compare , token ) ;
	    return GLU_FAILURE ;
	  }
	}
	continue ;
      }
      // +any others that may come up ...
    }
    break;
  }
  return GLU_SUCCESS ;
}
// YUCK

#ifdef DEBUG
  #undef DEBUG
#endif
