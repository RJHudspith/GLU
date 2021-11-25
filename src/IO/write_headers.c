/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (write_headers.c) is part of GLU.

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
   @file write_headers.c
   @brief allows for the writing out of header files   
 */
#include "Mainfile.h"

#include "GLU_bswap.h"   // byteswaps
#include "GLU_timer.h"   // for the date
#include "plaqs_links.h" // for the plaquette

// little helper for the MILC header
static void
swap_and_write_int( FILE *out ,
		    void *arr ,
		    const int N )
{
#ifdef OUT_BIG
  if( !WORDS_BIGENDIAN ) { bswap_32( N , arr ) ; }
#else
  if( WORDS_BIGENDIAN ) { bswap_32( N , arr ) ; }
#endif
  fwrite( arr , sizeof( uint32_t ) , N , out ) ;
}

// for writing out MILC configurations ...
void
write_header_MILC( FILE *__restrict out ,
		   const uint32_t milc_cksum29 ,
		   const uint32_t milc_cksum31 )
{
  // dirty, only write out correct format when compiled to single precision
#ifdef SINGLE_PREC
  int magic[ 1 ] = { 20103 } ;
#else
  int magic[ 1 ] = { 20104 } ;
#endif
  swap_and_write_int( out , magic , 1 ) ;

  int dims[ ND ] ;
  for (int i = 0; i < ND; i ++) dims[i] = Latt.dims[i];
  swap_and_write_int( out , dims , ND ) ;

  char str[ 64 ] ;
  sprintf( str , "%s GMT" , get_date( ) ) ;
  fprintf( stdout , "[IO] MILC time %s\n" , str ) ;
  fwrite( str , sizeof( char ) , 64 , out ) ;

  int output_type[ 1 ] = { 0 } ; //naturally-ordered type we support
  swap_and_write_int( out , output_type , 1 ) ;

  int checksums[ 2 ] = { milc_cksum29 , milc_cksum31 } ;
  swap_and_write_int( out , checksums , 2 ) ;
  
  return ;
}

// writes the usual nersc header ...
void
write_header_NERSC( FILE *__restrict out ,
		    const double tr ,
		    const double plaq ,
		    const uint32_t chksum ,
		    const char *details ,
		    const GLU_output type )
{
 // fill in the normal header crap //
 fprintf( out , "BEGIN_HEADER\n" ) ; 
 fprintf( out , "HDR_VERSION = 1.0\n" ) ;

 switch( type ) {
 case OUTPUT_SMALL : 
   fprintf( out , "DATATYPE = %dD_SU%d_GAUGE_SMALL\n" , 
	    ND , NC ) ; 
   break ;
 case OUTPUT_GAUGE : 
   fprintf( out , "DATATYPE = %dD_SU%d_GAUGE\n" , 
	    ND , NC ) ; 
   break ;
 default :
   fprintf( out , "DATATYPE = %dD_SU%d_GAUGE_%dx%d\n" , 
	    ND , NC , NC , NC ) ;  
   break ;
 }
 fprintf( out , "STORAGE_FORMAT = 1.0\n" ) ; 
 size_t mu ;
 for( mu = 0 ; mu < ND ;  mu++ ) {
   fprintf( out , "DIMENSION_%zu = %zu\n" , mu + 1 , Latt.dims[mu] ) ; 
 }
 fprintf( out , "CHECKSUM = %x\n" , chksum ) ; 
 fprintf( out , "LINK_TRACE = %1.15f\n" , tr ) ; 
 fprintf( out , "PLAQUETTE = %1.15f\n" , plaq ) ; 
 // only use PERIODIC
 for( mu = 0 ; mu < ND ; mu++ ) {
   fprintf( out , "BOUNDARY_%zu = PERIODIC\n" , mu + 1 ) ; 
 }
 fprintf( out , "ENSEMBLE_ID = ukqcd\n" ) ; 
 fprintf( out , "SEQUENCE_NUMBER = %zu\n" , Latt.flow ) ; 
 fprintf( out , "ENSEMBLE_LABEL = %s\n" , details ) ; 
 fprintf( out , "CREATOR_MACHINE = %s\n" , getenv( "USER" ) ) ; 
#ifdef HAVE_TIME_H
 fprintf( out , "CREATION_DATE = %s\n" , get_date( ) ) ;
#endif
 //this bit's a pain -> going to have to fight the machine endianness
#ifdef SINGLE_PREC
   #ifdef OUT_BIG
      fprintf( out , "FLOATING_POINT = IEEE32BIG\n" ) ; 
   #else
      fprintf( out , "FLOATING_POINT = IEEE32LITTLE\n" ) ; 
   #endif
#else
   #ifdef OUT_BIG
      fprintf( out , "FLOATING_POINT = IEEE64BIG\n" ) ; 
   #else
      fprintf( out , "FLOATING_POINT = IEEE64LITTLE\n" ) ; 
   #endif
#endif
 fprintf( out , "END_HEADER\n" ) ; 

 return ;
}
