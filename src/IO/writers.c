/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (writers.c) is part of GLU.

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
   @file writers.c
   @brief writes out a NERSC configuration
 */
#include "Mainfile.h"

#include "crc.h"           // MILC, SCIDAC and ILDG formats use a crc
#include "GLU_bswap.h"     // byteswapping
#include "GLU_memcheck.h"  // can we use the faster/mem expensive code?
#include "HIREP.h"         // write out a HIREP config
#include "plaqs_links.h"   // compute plaquette and links again
#include "Scidac.h"        // write out a scidac file
#include "write_headers.h" // write the config header

// avoid type-punning issues with unions
#ifdef SINGLE_PREC
typedef union {
  float val ;
  uint32_t chk ;
} U32 ;
#else
typedef union {
  double val ;
  uint32_t chk[2] ;
} U32 ;
#endif

// compute the checksum
static uint32_t
checksum( const GLU_real *__restrict p ,
	  const size_t idx )
{
  #ifdef SINGLE_PREC
  U32 in ;
  register uint32_t hi ;
  in.val = p[ idx ] ;
  hi = in.chk ;
  return (uint32_t)( hi ) ;
  #else
  U32 in ;
  register uint32_t hi , lo ;
  // point to the value
  in.val = p[ idx ] ;
  hi = in.chk[0] ;
  lo = in.chk[1] ; 
  return (uint32_t)( hi + lo ) ;
  #endif
}

// construct the loop variables ...
static GLU_output
construct_loop_variables( size_t *LATT_LOOP , 
			  size_t *LOOP_VAR ,
			  const GLU_output type )
{
  // set up our loop variables, SMALL is "special"
  switch( type ) {
  case OUTPUT_SMALL :
    fprintf( stdout , "[IO] Writing %dD_SU%d_GAUGE_SMALL binary data\n", 
	     ND , NC ) ;
    *LATT_LOOP = ND * LOOP_SMALL * LVOLUME ;
    *LOOP_VAR = LOOP_SMALL ;
    return OUTPUT_SMALL ;
  case OUTPUT_GAUGE : 
    fprintf( stdout , "[IO] Writing %dD_SU%d_GAUGE binary data\n", 
	     ND , NC ) ;
    *LATT_LOOP = ND * LOOP_GAUGE * LVOLUME ;
    *LOOP_VAR = LOOP_GAUGE ;
    return OUTPUT_GAUGE ;
  case OUTPUT_MILC :
  case OUTPUT_SCIDAC :
  case OUTPUT_ILDG :
  case OUTPUT_NCxNC :
    fprintf( stdout , "[IO] Writing %dD_SU%d_GAUGE_%dx%d binary data\n", 
	     ND , NC , NC , NC ) ;
    *LATT_LOOP = ND * LOOP_NCxNC * LVOLUME ;
    *LOOP_VAR = LOOP_NCxNC ;
    return OUTPUT_NCxNC ;
  default :
    fprintf( stderr , "[IO] Unrecognised output type .. Defaulting "
	     "to NERSC's %dD_SU%d_GAUGE binary data\n", ND , NC ) ;
    *LATT_LOOP = ND * LOOP_GAUGE * LVOLUME ;
    *LOOP_VAR = LOOP_GAUGE ;
    return OUTPUT_GAUGE ;
  }
}

// produce an array with the site's data stored
static void
grab_sitedata( GLU_real *__restrict chunk , 
	       const struct site lat ,
	       const size_t LOOP_VAR , 
	       const GLU_output checktype )
{
  GLU_real tmp[ LOOP_VAR ] ;
  size_t j , val = 0 , mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    switch( checktype ) {
    case OUTPUT_SMALL :
      shorten( tmp , lat.O[mu] ) ;
      for( j = 0 ; j < LOOP_VAR ; j++ ) {
	chunk[ val++ ] = tmp[ j ] ;
      }
      break ;
    default :
      for( j = 0 ; j < LOOP_VAR/2 ; j++ ) {
	chunk[ val++ ] = creal( lat.O[ mu ][ j ] ) ;
	chunk[ val++ ] = cimag( lat.O[ mu ][ j ] ) ;
      }
    }
  }
  return ;
}

// byteswap for output
static void
swap_for_output( GLU_real *uout , 
		 const size_t SIZE )
{
  // byte swaps
#ifdef OUT_BIG
  if( !WORDS_BIGENDIAN ) { 
    #ifdef SINGLE_PREC
    bswap_32( SIZE , uout ) ; 
    #else
    bswap_64( SIZE , uout ) ; 
    #endif
  }
#else // write a in little endian format
  if( WORDS_BIGENDIAN ) {
    #ifdef SINGLE_PREC
    bswap_32( SIZE , uout ) ; 
    #else
    bswap_64( SIZE , uout ) ; 
    #endif
  }
#endif
}

// compute checksums
static void
compute_checksum( uint32_t *nersc_cksum , 
		  uint32_t *milc_cksum29 , 
		  uint32_t *milc_cksum31 , 
		  uint32_t *scidac_cksum29 ,
		  uint32_t *scidac_cksum31 ,
		  const struct site *__restrict lat ,
		  const size_t LOOP_VAR ,
		  const GLU_output checktype )
{
  // initialise
  int i ;
  uint32_t k = 0 , sum29 = 0 , sum31 = 0 , CRCsum29 = 0 , CRCsum31 = 0 ;

  // loop some arbitrary length
  #pragma omp parallel for private(i) reduction(+:k) reduction(^:sum29) reduction(^:sum31) reduction(^:CRCsum29) reduction(^:CRCsum31) 
  for( i = 0 ; i < LVOLUME ; i++ ) { 

    // usual allocations
    GLU_real site[ ND * LOOP_VAR ] ;
    size_t rank29 = ( i * ND * LOOP_VAR ) % 29 ;
    size_t rank31 = ( i * ND * LOOP_VAR ) % 31 ;
    register uint32_t k_loc = 0 ;
    register uint32_t sum29_loc = 0 , sum31_loc = 0 ;
    uint32_t work ;

    // all of the site data is in GLU_real array site
    grab_sitedata( site , lat[i] , LOOP_VAR , checktype ) ;

    size_t j ;
    for( j = 0 ; j < ND * LOOP_VAR ; j++ ) {
      // compute the actual checksum
      work = (uint32_t)checksum( site , j ) ;

      // nersc checksum is just a sum of the data
      k_loc += work;
	
      // the first milc checksums
      sum29_loc ^= ( work << rank29 | work >> ( 32 - rank29 ) ) ;
      sum31_loc ^= ( work << rank31 | work >> ( 32 - rank31 ) ) ;
      rank29 = ( rank29 < 28 ) ? rank29 + 1 : 0 ;
      rank31 = ( rank31 < 30 ) ? rank31 + 1 : 0 ;
    }
    // and compute the CRC for the outputted data
    swap_for_output( site , ND * LOOP_VAR ) ;
    const size_t r29 = i % 29 ;
    const size_t r31 = i % 31 ;
    work = (uint32_t)crc32( 0 , (unsigned char*)site , 
			    ND * LOOP_VAR * sizeof( GLU_real ) ) ;
    CRCsum29 = CRCsum29 ^ (uint32_t) ( work << r29 | work >> ( 32 - r29 ) ) ;
    CRCsum31 = CRCsum31 ^ (uint32_t) ( work << r31 | work >> ( 32 - r31 ) ) ;

    // and the local values reductions
    // nersc
    k = k + (uint32_t)k_loc ;

    // milc
    sum29 = sum29 ^ sum29_loc ;
    sum31 = sum31 ^ sum31_loc ;
  }
  *nersc_cksum = k ;
  *milc_cksum29 = sum29 ;
  *milc_cksum31 = sum31 ;
  *scidac_cksum29 = CRCsum29 ;
  *scidac_cksum31 = CRCsum31 ;
  return ;
}

// can write in big or little endian, overwrites uout!
static void
byteswap_and_write( FILE *__restrict outfile ,
		    GLU_real *__restrict uout ,
		    const size_t SIZE )
{
  swap_for_output( uout , SIZE ) ; // do the byteswaps if necessary
  fwrite( uout , sizeof( GLU_real ) , SIZE , outfile ) ; 
  return ;
}

// data copying
static void
copy_data( GLU_real *uout , 
	   const struct site *lat ,
	   const size_t START ,
	   const size_t END ,
	   const size_t LOOP_VAR ,
	   const GLU_output checktype )
{
  size_t i ;
#pragma omp parallel for private(i)
  for( i = START ; i < END ; i++ ) { 
    // temporary
    GLU_real site[ ND * LOOP_VAR ] ;

    // set this up here
    size_t val = LOOP_VAR * ( i * ND ) ;

    // puts the matrix in glu_real temporary
    grab_sitedata( site , lat[i] , LOOP_VAR , checktype ) ;

    // loop AntiHermitian_projised matrix-direction index
    size_t j ;
    for( j = 0 ; j < ND * LOOP_VAR ; j++ ) {
      #ifdef OUT_DOUBLE
      uout[ val++ ] = (double)site[ j ] ; 
      #else
      uout[ val++ ] = (float)site[ j ] ; 
      #endif
    }
  }
  return ;
}

// writes the fields as binary data
static void
write_binary_data( const struct site *__restrict lat ,
		   FILE *__restrict outfile ,
		   const char *__restrict details ,
		   const GLU_output checktype ,
		   const size_t LATT_LOOP ,
		   const size_t LOOP_VAR )
{ 
  // write out in working precision
  GLU_real *uout = ( GLU_real* ) malloc( LATT_LOOP  * sizeof( GLU_real ) ) ; 

  // copy lattice data into uout 
  copy_data( uout , lat , 0 , LVOLUME , LOOP_VAR , checktype ) ;

  // and write it out
  byteswap_and_write( outfile , uout , LATT_LOOP ) ;

  free( uout ) ;
  return ;
}

// writes the fields as binary data -> SLICE by SLICE version
static void
write_binary_data_cheap( const struct site *__restrict lat ,
			 FILE *__restrict outfile ,
			 const char *__restrict details ,
			 const GLU_output checktype ,
			 const size_t LOOP_VAR )
{ 
  // write out in working precision
  GLU_real *uout = ( GLU_real* ) malloc( LCU * LOOP_VAR * 
					 ND * sizeof( GLU_real ) ) ; 

  size_t t ;
  for( t = 0 ; t < Latt.dims[ ND - 1 ] ; t++ ) {
    
    // copy lattice data into uout 
    copy_data( uout , lat , LCU*t , LCU*(t+1) , LOOP_VAR , checktype ) ;
    
    // and write it out
    byteswap_and_write( outfile , uout , LCU * LOOP_VAR * ND ) ;
  }

  free( uout ) ;
  return ;
}

// This is just the wrapping function for the writers
int
write_lat( struct site *__restrict lat , 
	   FILE *__restrict out , 
	   const GLU_output type , 
	   char *__restrict details )
{
  // HIREP stuff is all handled in HIREP.c
  if( type == OUTPUT_HIREP ) {
    printf( "[IO] Writing out in HIREP format\n" ) ;
    write_gauge_field( lat , out ) ;
    return GLU_SUCCESS ;
  }

  // loop variables
  size_t LOOP_VAR , LATT_LOOP ;
  const GLU_output checktype = construct_loop_variables( 
            &LATT_LOOP , &LOOP_VAR , type ) ;

  // compute the checksums ...
  uint32_t nersc_cksum , milc_cksum29 , milc_cksum31 ;
  uint32_t scidac_cksum29 , scidac_cksum31 ;

  compute_checksum( &nersc_cksum , 
		    &milc_cksum29 , &milc_cksum31 , 
		    &scidac_cksum29 , &scidac_cksum31 , 
		    lat , LOOP_VAR , checktype ) ;

  // and begin writing the headers
  switch( type ) {
  case OUTPUT_HIREP : return GLU_FAILURE ;
  case OUTPUT_SMALL : 
  case OUTPUT_GAUGE :
  case OUTPUT_NCxNC :
    write_header_NERSC( out , links(lat) , av_plaquette(lat) , 
			nersc_cksum , details , checktype ) ;
    break ;
  case OUTPUT_MILC :
    write_header_MILC( out , milc_cksum29 , milc_cksum31 ) ;
    break ;
  case OUTPUT_SCIDAC :
    write_header_SCIDAC( out ) ;
    break ;
  case OUTPUT_ILDG :
    write_header_ILDG( out ) ;
    break ;
  default :
    return GLU_FAILURE ;
  }

  // binary output is all the same for these types ...
  const int SAFETY = have_memory_readers_writers( type ) ;
  if( SAFETY != FAST ) {
    fprintf( stdout , "[IO] Using the slower/memory cheaper code \n" ) ;
    write_binary_data_cheap( lat , out , details , type , LOOP_VAR ) ;
  } else {
    fprintf( stdout , "[IO] Using the faster/memory expensive code \n" ) ;
    write_binary_data( lat , out , details , type , LATT_LOOP , LOOP_VAR ) ;
  }

  // need to write out the scidac checksums here
  switch( type ) {
  case OUTPUT_SCIDAC :
  case OUTPUT_ILDG :
    write_trailing_header_SCIDAC( out , scidac_cksum29 , scidac_cksum31 ) ;
    break ;
  default :
    break ;
  }
  return GLU_SUCCESS ;
}

