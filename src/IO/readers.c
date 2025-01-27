/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (readers.c) is part of GLU.

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
   @file readers.c
   @brief binary data file reader supports #config_size outputs
 */
#include "Mainfile.h"

#include "crc.h"           // for the scidac circular checksum
#include "LU.h"            // LU_det_overwrite()
#include "GLU_bswap.h"     // bswap_*()
#include "gramschmidt.h"   // reunit2()
#include "random_config.h" // latt_reunitU()

// read a full matrix, should we just read in all 
// of the complex values?
static inline void 
complete_NCxNC( GLU_complex *__restrict O ,
		const GLU_real *__restrict uout )
{
  memcpy( O , uout , NCNC * sizeof( GLU_complex ) ) ;
  return ;
}

// recombines O from the top NC-1 rows of matrix
static void 
complete_top( GLU_complex *__restrict O , 
	      const GLU_real *__restrict uout ) 
{
  memcpy( O , uout , ( NCNC-NC ) * sizeof( GLU_complex ) ) ;
#if NC == 3
  O[6] = conj( O[1] * O[5] - O[2] * O[4] ) ;
  O[7] = conj( O[2] * O[3] - O[0] * O[5] ) ;
  O[8] = conj( O[0] * O[4] - O[1] * O[3] ) ;
#elif NC == 2
  O[2] = -conj( O[1] ) ;
  O[3] =  conj( O[0] ) ;
#else
  // and complete, taken from the gramschmidt code, should consider a minors function ?
  size_t i ;
  GLU_complex array[ ( NC - 1 ) * ( NC - 1 ) ] ;
  for( i = (NCNC-NC) ; i < NCNC ; i++ ) { // our bona-fide minor index
    size_t idx = 0 , j ;
    for( j = 0 ; j < ( NCNC - NC ) ; j++ ) {
      if( ( j%NC != i%NC ) ) { // remove columns and rows
	// pack array
	array[idx] = O[j] ;
	idx ++ ;
      } 
    }
    // compute the minors of the bottom row
    #if ( NC%2 == 0 )
    register const GLU_real mulfact = ( i % 2 == 0 ) ? -1.0 : 1.0 ; 
    #else
    register const GLU_real mulfact = ( i % 2 == 0 ) ? 1.0 : -1.0 ; 
    #endif
    O[i]= conj( mulfact * (GLU_complex)LU_det_overwrite( NC-1 , array ) ) ;
  }
#endif
  return ;
}

// construct some common LOOP variables
static int
construct_loop_variables( size_t *LATT_LOOP , 
			  size_t *LOOP_VAR ,
			  const GLU_output type )
{
  // dump it all in memory...
  switch( type ) {
  case OUTPUT_SMALL :
    *LATT_LOOP = ND * LOOP_SMALL * Latt.Volume ;
    *LOOP_VAR = LOOP_SMALL ;
    return GLU_SUCCESS ;
  case OUTPUT_GAUGE : 
    *LATT_LOOP = ND * LOOP_GAUGE * Latt.Volume ;
    *LOOP_VAR = LOOP_GAUGE ;
    return GLU_SUCCESS ;
  case OUTPUT_NCxNC :
    *LATT_LOOP = ND * LOOP_NCxNC * Latt.Volume ;
    *LOOP_VAR = LOOP_NCxNC ;
    return GLU_SUCCESS ;
  default :
    // actually should try and read an NCxNC config perhaps?
    fprintf( stderr , "[IO] Unrecognised input type .. "
	     "Leaving in disgust\n" ) ;
    return GLU_FAILURE ;
  }
}

// CRC checksum calculator
// name and idea come from ETMC
static void 
DML_checksum_accum( uint32_t *checksuma , 
		    uint32_t *checksumb , 
		    const uint32_t rank, 
		    char *buf, 
		    const size_t size )
{
  const uint32_t rank29 = rank % 29 ;
  const uint32_t rank31 = rank % 31 ;
  const uint32_t work = (uint32_t)crc32(0, (unsigned char*)buf, size);
  *checksuma ^= work<<rank29 | work>>(32-rank29);
  *checksumb ^= work<<rank31 | work>>(32-rank31);
  return ;
}

#if NC<4
// recreate O from the "short" definition
static void 
exhume_O( GLU_complex *__restrict S , 
	  const GLU_real *__restrict uout ) 
{
  rebuild( S , uout ) ;
  gram_reunit( S ) ;
  return ;
}
#endif

// rebuilt using the switch for the different allowd formats
static void
rebuild_lat( GLU_complex *__restrict link ,
	     const GLU_real *__restrict utemp ,
	     const GLU_output type )
{
  // smash all the read values into lat
  switch( type ) { 
  case OUTPUT_SMALL :
    #if NC<4
    exhume_O( link , utemp ) ;
    #else
    exit(1) ;
    #endif
    return ;
  case OUTPUT_GAUGE :
    complete_top( link , utemp ) ;
    return ;
  case OUTPUT_NCxNC :
    complete_NCxNC( link , utemp ) ;
    return ;
  default : return ;
  }  
  return ;
}

//////////// MEMSPENSIVE VERSION ///////////////
uint32_t
lattice_reader_suNC( struct site *lat , 
		     FILE *__restrict in , 
		     const struct head_data HEAD_DATA )
{
  // this is checked previously, nice to be certain though
  if( in == NULL ) {
    fprintf( stderr , "[IO] Error opening config file!!\n" ) ; 
    return GLU_FAILURE ; 
  }

  // loop variables
  size_t LOOP_VAR , LATT_LOOP ;
  if( construct_loop_variables( &LATT_LOOP , &LOOP_VAR , 
				HEAD_DATA.config_type ) == GLU_FAILURE ) {
    return GLU_FAILURE ;
  }

  static double *uind , *p ; 
  static float *uin , *q ; 
  uint32_t CRCsum29 = 0 , CRCsum31 = 0 , CRC_BQCD = 0 ;

  if( HEAD_DATA.precision == DOUBLE_PREC ) {
    uind = ( double* )malloc( LATT_LOOP * sizeof( double ) ) ; 
    if( fread( uind , sizeof( double ) , LATT_LOOP , in ) != LATT_LOOP ) {
      fprintf( stderr , "[IO] Configuration File reading failure "
	       ".. Leaving \n" ) ;
      free( uind ) ;
      return GLU_FAILURE ;
    }
    // scidac checksum is on the RAW binary data, not the byteswapped
    if( Latt.head == SCIDAC_HEADER ||
	Latt.head == ILDG_SCIDAC_HEADER ) {
      size_t i ;
      #pragma omp parallel for private(i) reduction(^:CRCsum29) reduction(^:CRCsum31)
      for( i = 0 ; i < LVOLUME ; i++ ) {
	const uint32_t rank29 = i % 29 ;
	const uint32_t rank31 = i % 31 ;
	const uint32_t work =
	  (uint32_t)crc32(0, (unsigned char*)( uind + i*ND*LOOP_VAR ) ,
			  sizeof(double) * ND * LOOP_VAR );
	CRCsum29 = CRCsum29 ^ ( work<<rank29 | work>>(32-rank29));
	CRCsum31 = CRCsum31 ^ ( work<<rank31 | work>>(32-rank31) );
	// TDOD -> an unthreaded bswap here ?
      }
    }
    // BQCD checksum is just the crc of the whole thing and is not thread safe
    // for the moment possibly for ever
    if( Latt.head == ILDG_BQCD_HEADER ) {
      size_t i ;
      for( i = 0 ; i < LVOLUME ; i++ ) {
	// BQCD's
	CKSUM_ADD( ( uin + ( i * ND * LOOP_VAR ) ) , 
		   sizeof(float) * ND * LOOP_VAR ) ;
      }
      uint32_t nbytes ;
      CKSUM_GET( &CRC_BQCD , &nbytes ) ;
    }
    // and then we do the byte swap
    if( HEAD_DATA.endianess != WORDS_BIGENDIAN ) {
      bswap_64( LATT_LOOP , uind ) ;  
    }
    p = uind ; 
  } else {
    uin = ( float* )malloc( LATT_LOOP * sizeof( float ) ) ; 
    if( fread( uin , sizeof( float ) , LATT_LOOP  , in ) != LATT_LOOP ) {
      fprintf( stderr , "[IO] Configuration File reading failure "
	       ".. Leaving \n" ) ;
      free( uin ) ;
      return GLU_FAILURE ;
    }
    // ILDG checksum is on the raw data
    if( Latt.head == SCIDAC_HEADER ||
	Latt.head == ILDG_SCIDAC_HEADER ) {
      size_t i ;
      #pragma omp parallel for private(i) reduction(^:CRCsum29) reduction(^:CRCsum31)
      for( i = 0 ; i < LVOLUME ; i++ ) {
	const uint32_t rank29 = i % 29 ;
	const uint32_t rank31 = i % 31 ;
	const uint32_t work =
	  (uint32_t)crc32(0, (unsigned char*)( uin + i*ND*LOOP_VAR ) ,
			  sizeof(float) * ND * LOOP_VAR );
	CRCsum29 = CRCsum29 ^ ( work<<rank29 | work>>(32-rank29));
	CRCsum31 = CRCsum31 ^ ( work<<rank31 | work>>(32-rank31) );
	// TDOD -> an unthreaded bswap here ?
      }
    }
    // BQCD checksum is just the crc of the whole thing and is kinda expensive
    if( Latt.head == ILDG_BQCD_HEADER ) {
      size_t i ;
      for( i = 0 ; i < LVOLUME ; i++ ) {
	// BQCD's
	CKSUM_ADD( ( uin + ( i * ND * LOOP_VAR ) ) , 
		   sizeof(float) * ND * LOOP_VAR ) ;
      }
      uint32_t nbytes ;
      CKSUM_GET( &CRC_BQCD , &nbytes ) ;
    }
    // and byteswap if necessary
    if( HEAD_DATA.endianess != WORDS_BIGENDIAN ) {
      bswap_32( LATT_LOOP , uin ) ;
    }
    q = uin ;
  }

  // compute the crcs and poke into gauge links
  uint32_t k = 0 , sum29 = 0 , sum31 = 0 ;
  size_t i ;
  #pragma omp parallel for private(i) reduction(+:k) reduction(^:sum29) reduction(^:sum31)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    // t is the AntiHermitian_projised idx
    size_t t = ND * LOOP_VAR * i ;
    size_t rank29 = t % 29 ;
    size_t rank31 = t % 31 ;
    register uint32_t k_loc = 0 ;
    register uint32_t sum29_loc = 0 , sum31_loc = 0 ;

    // general variables ...
    GLU_real utemp[ LOOP_VAR ] ;
    uint32_t res = 0 ;
    size_t mu , j , count = 0 ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      for( j = 0 ; j < LOOP_VAR ; j++ ) {
	// should I shuffle this around ?
	if( HEAD_DATA.precision == DOUBLE_PREC ) {
	  // compute the checksum ...
	  uint32_t *buf = ( uint32_t* )( p + t ) ; 
	  res = *buf + *( buf + 1 ) ; 
	  // put value into temporary
	  *( utemp + j ) = ( GLU_real )*( p + t ) ;
	} else {
	  // nersc checksum ...
	  res = *( uint32_t *)( q + t ) ;
	  // and put the value in the temporary
	  *( utemp + j ) = ( GLU_real )*( q + t ) ; 
	}
	count++ ;

	// milc checksums ...
	sum29_loc ^= (uint32_t)( res << rank29 | res >> ( 32 - rank29 ) ) ;
	sum31_loc ^= (uint32_t)( res << rank31 | res >> ( 32 - rank31 ) ) ;

	/// and perform the mods
	rank29 = ( rank29 < 28 ) ? rank29 + 1 : 0 ;
	rank31 = ( rank31 < 30 ) ? rank31 + 1 : 0 ;

	// local sum
	k_loc += res ; 

	t++ ;
      }
      // smash all the read values into lat
      rebuild_lat( lat[i].O[mu] , utemp , HEAD_DATA.config_type ) ;
    }
    // reductions go here ...
    // nersc
    k = k + (uint32_t)k_loc ;

    // milc
    sum29 = sum29 ^ (uint32_t)sum29_loc ;
    sum31 = sum31 ^ (uint32_t)sum31_loc ;
  }

  if( HEAD_DATA.precision == DOUBLE_PREC ) {
    free( uind ) ; 
  } else {
    // reunitarise up to working precision
    #ifndef SINGLE_PREC
    latt_reunitU( lat ) ;
    #endif
    free( uin ) ;
  }

  // if we are reading a MILC file we output the sum29 checksum
  switch( Latt.head ) {
  case NERSC_HEADER : return k ;
  case MILC_HEADER : return sum29 ;
  case ILDG_SCIDAC_HEADER :
  case SCIDAC_HEADER : return CRCsum29 ;
  case ILDG_BQCD_HEADER : return CRC_BQCD ;
  case LIME_HEADER : return GLU_SUCCESS ;
  case JLQCD_HEADER : return GLU_SUCCESS ;
  default : return GLU_FAILURE ; // should print an error here
  }
}

// MEMCHEAP READER
uint32_t
lattice_reader_suNC_cheaper( struct site *lat , 
			     FILE *__restrict in , 
			     const struct head_data HEAD_DATA )
{
  // this is checked previously, nice to be certain though
  if( in == NULL ) {
    fprintf( stderr , "[IO] Error opening config file!!..."
	     "Returning with error \n" ) ; 
    return GLU_FAILURE ; 
  }

  // loop variables
  size_t LOOP_VAR , LATT_LOOP ;
  if( construct_loop_variables( &LATT_LOOP , &LOOP_VAR , 
				HEAD_DATA.config_type ) == GLU_FAILURE ) {
    return GLU_FAILURE ;
  }

  static double *uind , *p ; 
  static float *uin , *q ; 

  if( HEAD_DATA.precision == DOUBLE_PREC ) {
    uind = ( double* )malloc( ND*LOOP_VAR * sizeof( double ) ) ; 
  } else {
    uin = ( float* )malloc( ND*LOOP_VAR * sizeof( float ) ) ; 
  }

  uint32_t k = 0 , sum29 = 0 , sum31 = 0 ; // NERSC & MILC checksums ...
  uint32_t CRCsum29 = 0 , CRCsum31 = 0 ;
  size_t i ;
  // do this in serial reads in site by site
  for( i = 0 ; i < LVOLUME ; i++ ) {

    // read in a site, could we not use fseek to allow this to be done in parallel?
    if( HEAD_DATA.precision == DOUBLE_PREC ) {
      if( fread( uind , sizeof( double ) , ND*LOOP_VAR , in ) != ND*LOOP_VAR ) {
	printf("[IO] Configuration File reading failure .. Leaving \n" ) ;
	free( uind ) ;
	return GLU_FAILURE ;
      }
      DML_checksum_accum( &CRCsum29 , &CRCsum31 , 
			  i , (char*)( uind ) ,
			  sizeof(double) * ND * LOOP_VAR ) ;
      // BQCD's
      CKSUM_ADD( uind , sizeof(double) * ND * LOOP_VAR ) ;
      if( HEAD_DATA.endianess != WORDS_BIGENDIAN ) {
	bswap_64( ND*LOOP_VAR , uind ) ;  
      }
      p = uind ;
    } else { 
      if( fread( uin , sizeof( float ) , ND*LOOP_VAR , in ) != ND*LOOP_VAR ) {
	printf("[IO] Configuration File reading failure .. Leaving \n" ) ;
	free( uin ) ;
	return GLU_FAILURE ;
      }
      DML_checksum_accum( &CRCsum29 , &CRCsum31 , 
			  i , (char*)( uin ) ,
			  sizeof(float) * ND * LOOP_VAR ) ;
      // BQCD's
      CKSUM_ADD( uin , sizeof(float) * ND * LOOP_VAR ) ;
      if( HEAD_DATA.endianess != WORDS_BIGENDIAN ) {
	bswap_32( ND*LOOP_VAR , uin ) ;
      }
      q = uin ;
    }
    size_t t = 0 ;

    register uint32_t k_loc = 0 , sum29_loc = 0 , sum31_loc = 0 ;
    size_t rank29 = ( ND * LOOP_VAR * i ) % 29 ;
    size_t rank31 = ( ND * LOOP_VAR * i ) % 31 ;

    // general variables ...
    GLU_real utemp[ LOOP_VAR ] ;
    uint32_t res = 0 ;
    size_t mu , j ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      for( j = 0 ; j < LOOP_VAR ; j++ ) {
	// should I shuffle this around ?
	if( HEAD_DATA.precision == DOUBLE_PREC ) {
	  // compute the checksum ...
	  uint32_t *buf = ( uint32_t* )( p + t ) ; 
	  res = *buf + *( buf + 1 ) ; 
	  // put value into temporary
	  *( utemp + j ) = ( GLU_real )*( p + t ) ;
	} else {
	  // nersc checksum ...
	  res = *( uint32_t *)( q + t ) ;
	  // and put the value in the temporary
	  *( utemp + j ) = ( GLU_real )*( q + t ) ; 
	} 
	// local computations
	sum29_loc ^= (uint32_t)( res << rank29 | res >> ( 32 - rank29 ) ) ;
	sum31_loc ^= (uint32_t)( res << rank31 | res >> ( 32 - rank31 ) ) ;

	/// and perform the mods
	rank29 = ( rank29 < 28 ) ? rank29 + 1 : 0 ;
	rank31 = ( rank31 < 30 ) ? rank31 + 1 : 0 ;

	// local sum
	k_loc += res ; 

	t++ ;
      }
      // smash all the read values into lat
      rebuild_lat( lat[i].O[mu] , utemp , HEAD_DATA.config_type ) ;
    }
    // nersc
    k = k + (uint32_t)k_loc ;

    // milc
    sum29 = sum29 ^ (uint32_t)sum29_loc ;
    sum31 = sum31 ^ (uint32_t)sum31_loc ;
  }

  // BQCD checksum is just the crc of the whol thing
  uint32_t CRC_BQCD , nbytes ;
  CKSUM_GET( &CRC_BQCD , &nbytes ) ;

  if( HEAD_DATA.precision == DOUBLE_PREC ) {
    free( uind ) ; 
  } else {
    // reunitarise up to working precision
    #ifndef SINGLE_PREC
    latt_reunitU( lat ) ;
    #endif
    free( uin ) ;
  }

  // if we are reading a MILC file we output the sum29 checksum
  switch( Latt.head ) {
  case NERSC_HEADER : return k ;
  case MILC_HEADER : return sum29 ;
  case ILDG_SCIDAC_HEADER :
  case SCIDAC_HEADER : return CRCsum29 ;
  case ILDG_BQCD_HEADER : return CRC_BQCD ;
  case LIME_HEADER : return GLU_SUCCESS ;
  case JLQCD_HEADER : return GLU_SUCCESS ;
  default : return GLU_FAILURE ; // should print an error here
  }
}

// clean this up for scope
#ifdef DEBUG
  #undef DEBUG
#endif
