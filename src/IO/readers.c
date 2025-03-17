/*
  Copyright 2013-2025 Renwick James Hudspith
  
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
   @brief binary data file reader fully POSIX'd
 */
#include "Mainfile.h"

// POSIX IO crap
#include <sys/mman.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

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
static inline void 
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
    fprintf( stderr , "[IO] Unrecognised input type .. Leaving in disgust\n" ) ;
    return GLU_FAILURE ;
  }
}

// CRC checksum calculator
// name and idea come from ETMC
static inline void 
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
static inline void 
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
    memcpy( link , utemp , NCNC*sizeof( GLU_complex ) ) ;
    return ;
  default : return ;
  }  
  return ;
}

uint32_t
lattice_reader_suNC_posix( struct site *lat ,
			   const char *config_in ,
			   FILE *__restrict in ,
			   const struct head_data HEAD_DATA )
{
  // get the current position in the file as we should have read the file header information
  const size_t bstart = ftell( in ) ;

  fprintf( stdout , "[IO] file offset %zu\n" , bstart ) ;

  // loop variables
  size_t LOOP_VAR , LATT_LOOP ;
  if( construct_loop_variables( &LATT_LOOP , &LOOP_VAR , 
				HEAD_DATA.config_type ) == GLU_FAILURE ) {
    return GLU_FAILURE ;
  }
  fprintf( stdout , "[IO] LOOP_VAR = %zu\n" , LOOP_VAR ) ;

  // memory map the file
  int fd = open( config_in , O_RDONLY , S_IRUSR ) ;
  struct stat sb ;
  if( fstat( fd , &sb ) == -1 ) {
    perror( "fstat" ) ;
    exit(1) ;
  }

  // create mmap instance and try and reduce pagefaults
  char *mm = mmap( NULL , sb.st_size , PROT_READ , MAP_SHARED | MAP_POPULATE , fd , 0 ) ;
  madvise( mm , sb.st_size , MADV_HUGEPAGE ) ;
  mm += bstart ;

  // pointers to the start of the binary data
  const double *uind = (const double*)mm ;
  const float  *uinf = (const float*)mm ;

  // don't know why I want to support this
  uint32_t CRC_BQCD = 0 ;
  if( Latt.head == ILDG_BQCD_HEADER ) {
    if( HEAD_DATA.precision == DOUBLE_PREC ) {
      for( size_t i = 0 ; i < LVOLUME ; i++ ) {
	CKSUM_ADD( ( uind + ( i*ND*LOOP_VAR ) ) , sizeof(double)*ND*LOOP_VAR ) ;
      }
    } else {
      for( size_t i = 0 ; i < LVOLUME ; i++ ) {
	CKSUM_ADD( ( uinf + ( i*ND*LOOP_VAR ) ) , sizeof(float)*ND*LOOP_VAR ) ;
      }
    }
    uint32_t nbytes ;
    CKSUM_GET( &CRC_BQCD , &nbytes ) ;
  }
  
  // reduction checksums only do the 29 guy
  uint32_t k = 0 , sum29 = 0 , CRCsum29 = 0 ;

  // idea is to fill a buffer and then copy to lat and do checksum reductions
  size_t i ;
#pragma omp parallel for private(i) reduction(+:k) reduction(^:sum29) reduction(^:CRCsum29)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    double dbl_temp[ ND*LOOP_VAR ] GLUalign ;
    float flt_temp[ ND*LOOP_VAR ] GLUalign ;

    // the one in our own precision
    GLU_real utemp[ LOOP_VAR] GLUalign ;

    // checksum shit
    size_t rank29 = (ND*LOOP_VAR*i)%29 , CRCrank29 = i%29 ;
    register uint32_t k_loc = 0 , sum29_loc = 0 , res = 0 ;
    rank29 = ( rank29 < 28 ) ? rank29 + 1 : 0 ;

    // fill the buffers and maybe do the crc32 if needed otherwise fuck it lol
    // these are on the un-byte-swapped data
    if( HEAD_DATA.precision == DOUBLE_PREC ) {
      memcpy( dbl_temp , uind+(LOOP_VAR*(ND*i)) , sizeof( double )*ND*LOOP_VAR ) ;
      if( Latt.head == SCIDAC_HEADER || Latt.head == ILDG_SCIDAC_HEADER ) {
	res = (uint32_t)crc32(0, (unsigned char*)dbl_temp , sizeof(double)*ND*LOOP_VAR );
      }
      if( HEAD_DATA.endianess != WORDS_BIGENDIAN ) {
	bswap_64( ND*LOOP_VAR , dbl_temp ) ;
      }
    } else {
      memcpy( flt_temp , uinf+(LOOP_VAR*(ND*i)) , sizeof( float )*LOOP_VAR*ND ) ;
      if( Latt.head == SCIDAC_HEADER || Latt.head == ILDG_SCIDAC_HEADER ) {
	res = (uint32_t)crc32(0, (unsigned char*)flt_temp , sizeof(float)*ND*LOOP_VAR );
      }
      if( HEAD_DATA.endianess != WORDS_BIGENDIAN ) {
	bswap_32( LOOP_VAR*ND , flt_temp ) ;
      }
    }
    CRCsum29 = CRCsum29^( res<<CRCrank29 | res>>(32-CRCrank29) );
    
    for( int mu = 0 ; mu < ND ; mu++ ) {
      if( HEAD_DATA.precision == DOUBLE_PREC ) {
	for( size_t j = 0 ; j < LOOP_VAR ; j++ ) {
	  utemp[j] = dbl_temp[j+LOOP_VAR*mu] ;
	  const uint32_t *buf = (uint32_t*)&dbl_temp[j+LOOP_VAR*mu] ;
	  res = buf[0]+buf[1] ;
	  sum29_loc ^= (uint32_t)( res << rank29 | res >> ( 32 - rank29 ) ) ;	  
	  k_loc += res ;
	}
      } else {
	for( size_t j = 0 ; j < LOOP_VAR ; j++ ) {
	  utemp[j] = flt_temp[j+LOOP_VAR*mu] ;
	  const uint32_t *buf = (uint32_t*)&flt_temp[j+LOOP_VAR*mu] ;
	  res = buf[0] ;
	  sum29_loc ^= (uint32_t)( res << rank29 | res >> ( 32 - rank29 ) ) ;
	  k_loc += res ;
	}
      }
      // smash all the read values into lat
      rebuild_lat( lat[i].O[mu] , utemp , HEAD_DATA.config_type ) ;
    }
    k = k + (uint32_t)k_loc ;

    // milc
    sum29 = sum29 ^ (uint32_t)sum29_loc ;
  }

  // close the file descriptor
  close( fd ) ;
  munmap( mm , sb.st_size ) ;

  // seek to the FILE to the end of the data for looking up stupid checksums
  if( HEAD_DATA.precision == DOUBLE_PREC ) {
    fseek( in , LATT_LOOP*sizeof(double) , SEEK_CUR ) ;
  } else {
    fseek( in , LATT_LOOP*sizeof(float) , SEEK_CUR ) ;
  }

  // if we are reading a MILC file we output the sum29 checksum
  switch( Latt.head ) {
  case NERSC_HEADER       : return k ;
  case MILC_HEADER        : return sum29 ;
  case ILDG_SCIDAC_HEADER : case SCIDAC_HEADER : return CRCsum29 ;
  case ILDG_BQCD_HEADER   : return CRC_BQCD ;
  case LIME_HEADER        : return GLU_SUCCESS ;
  case JLQCD_HEADER       : return GLU_SUCCESS ;
  default                 :
    fprintf( stderr , "[IO] unrecognised header\n" ) ;
    return GLU_FAILURE ;
  }
}
