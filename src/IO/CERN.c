/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (CERN.c) is part of GLU.

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
   @file CERN.c
   @brief code to read and write CERN configs generated by openQCD

   bizarrely their geometry is
   t,x,y,z which makes little sense to me. I imagine there is a very good reason
   for this.... Probably. Anyway, fully memory-mapped and chunked IO here which I will probably do
   for the NERSC files too

   @warning these aren't #ND generic although in principle this wouldn't be a big problem to extend
 */
#include "Mainfile.h"

#include <sys/mman.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include "geometry.h"      // gen_site()
#include "GLU_bswap.h"     // byte swapping arrays
#include "plaqs_links.h"   // compute the plaquette
#include "GLU_timer.h"

// converts CLS to GLU. Assumes CLS is t,x,y,z and GLU is t,z,y,x
// assume this continues to ND-1,0,1,2,3.. ND-1,ND-2,....,1,0
static inline size_t
CLS_to_GLU( int *is_odd,
	    const size_t i )
{
  int x[ ND ] , sum = 0 ;
  size_t subvol = 1 ;
  for( int mu = 0 ; mu < ND-1 ; mu++ ) {
    x[ ND-2-mu ] = ( ( i - i % subvol ) / subvol ) % Latt.dims[ mu ] ;
    subvol *= Latt.dims[ mu ] ;
    sum += x[ ND-2-mu ] ;
  }
  x[ ND-1 ] = ( ( i - i % subvol ) / subvol ) % Latt.dims[ ND-1 ] ;
  sum += x[ ND-1 ] ;
  *is_odd = sum&1 ;  
  return gen_site( x ) ;
}

// filling CLS single site fields with our gauge fields 
static inline void
field_to_out( double complex *uout ,
	      const size_t i ,
	      const struct site *lat ,
	      const size_t idx )
{
  double complex *u = (double complex*)(uout + 2*ND*NCNC*(i/2)) ;
  // create a lookup table with the hope that this gets cached
#ifndef SINGLE_PREC
  const struct site *lp[5] = { &lat[idx] ,
			       &lat[lat[idx].back[3]] ,
			       &lat[lat[idx].back[0]] ,
			       &lat[lat[idx].back[1]] ,
			       &lat[lat[idx].back[2]] } ;
  memcpy( u        , lp[0]->O[3] , NCNC*sizeof(GLU_complex) ) ;
  memcpy( u+NCNC   , lp[1]->O[3] , NCNC*sizeof(GLU_complex) ) ;
  memcpy( u+2*NCNC , lp[0]->O[0] , NCNC*sizeof(GLU_complex) ) ;
  memcpy( u+3*NCNC , lp[2]->O[0] , NCNC*sizeof(GLU_complex) ) ;
  memcpy( u+4*NCNC , lp[0]->O[1] , NCNC*sizeof(GLU_complex) ) ;
  memcpy( u+5*NCNC , lp[3]->O[1] , NCNC*sizeof(GLU_complex) ) ;
  memcpy( u+6*NCNC , lp[0]->O[2] , NCNC*sizeof(GLU_complex) ) ;
  memcpy( u+7*NCNC , lp[4]->O[2] , NCNC*sizeof(GLU_complex) ) ;
#else
  // t first
  size_t shift = lat[idx].back[ND-1] , j ;
  for( j = 0 ; j < NCNC ; j++ ) {
    u[j] = (double complex)lat[idx].O[ND-1][j] ;
  }
  u += NCNC ;
  for( j = 0 ; j < NCNC ; j++ ) {
    u[j] = (double complex)lat[shift].O[ND-1][j] ;
  }
  u += NCNC ;
  // and then the others
  for( int mu = 0 ; mu < ND-1 ; mu++ ) {
    size_t shift = lat[idx].back[mu] ;
    for( j = 0 ; j < NCNC ; j++ ) {
      u[j] = (double complex)lat[idx].O[mu][j] ;
    }
    u+= NCNC ;
    for( j = 0 ; j < NCNC ; j++ ) {
      u[j] = (double complex)lat[shift].O[mu][j] ;
    }
    u+= NCNC ;
  }
#endif
}

// read a CERN gauge field
int
read_CLS_field( struct site *__restrict lat , 
		const char *config_in ,
		uint32_t *chksum )
{
  if( ND != 4 ) return GLU_FAILURE ;
#ifdef WORDS_BIFENDIAN
  fprintf( stderr , "[IO] CERN we don't support having a different file endianess yet\n" ) ;
  return GLU_FAILURE ;
#endif

  // open the file
  int fd = open( config_in , O_RDONLY , S_IRUSR ) ;
  
  struct stat sb;
  if( fstat(fd, &sb) == -1 ){
    perror( "fstat" ) ;
    return GLU_FAILURE ;
  }
  
  // memory map the file so the OS can just deal with it
  char *mm = mmap( NULL , sb.st_size , PROT_READ , MAP_PRIVATE , fd , 0 ) ;

  // usual CERN junk at the top of the file : dimensions and plaquette
  uint32_t *dims = (uint32_t*)mm ;
  fprintf( stdout , "[IO] CLS dims xyzt (%dx%dx%dx%d)\n" ,
	   dims[1] , dims[2] , dims[3] , dims[0] ) ;
  mm += 4*sizeof(uint32_t) ;
  
  double *plq = (double*)mm ;
  fprintf( stdout , "[IO] CLS header plaquette %1.15f\n" , *plq/(double)NC ) ;
  mm += sizeof(double);

  // pun out into a nicer type for us
  const double complex *uin = (const double complex*)mm ;
  
  // loop the CLS checkerboard volume in parallel
  size_t i ;
  #pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    // CLS order
    int is_odd = 0 ;
    const size_t idx = CLS_to_GLU( &is_odd , i ) ;
    if( is_odd ) {
      // local pointer
      const double complex *u = (const double complex*)(uin + 2*ND*NCNC*(i/2)) ;
      #ifndef SINGLE_PREC
      // create a lookup table with the hope that this gets cached
      struct site *lp[5] = { &lat[ idx ] ,
			     &lat[ lat[idx].back[3] ] ,
			     &lat[ lat[idx].back[0] ] ,
			     &lat[ lat[idx].back[1] ] ,
			     &lat[ lat[idx].back[2] ] } ;      
      memcpy( lp[0]->O[3] , u        , NCNC*sizeof(GLU_complex) ) ;
      memcpy( lp[1]->O[3] , u+NCNC   , NCNC*sizeof(GLU_complex) ) ;
      memcpy( lp[0]->O[0] , u+2*NCNC , NCNC*sizeof(GLU_complex) ) ;
      memcpy( lp[2]->O[0] , u+3*NCNC , NCNC*sizeof(GLU_complex) ) ;
      memcpy( lp[0]->O[1] , u+4*NCNC , NCNC*sizeof(GLU_complex) ) ;
      memcpy( lp[3]->O[1] , u+5*NCNC , NCNC*sizeof(GLU_complex) ) ;
      memcpy( lp[0]->O[2] , u+6*NCNC , NCNC*sizeof(GLU_complex) ) ;
      memcpy( lp[4]->O[2] , u+7*NCNC , NCNC*sizeof(GLU_complex) ) ;
      #else
      size_t mu , j ;
      size_t shift = lat[idx].back[ND-1] ;
      // t first
      for( j = 0 ; j < NCNC ; j++ ) {
	lat[idx].O[ND-1][j] = (GLU_real)creal(u[j])+I*(GLU_real)cimag(u[j]) ;
      }
      u += NCNC ;
      for( j = 0 ; j < NCNC ; j++ ) {
	lat[shift].O[ND-1][j] = (GLU_real)creal(u[j])+I*(GLU_real)cimag(u[j]) ;
      }
      u += NCNC ;
      // then the others (xyz)
      for( mu = 0 ;  mu < ND - 1 ; mu++ ) {
	for( j = 0 ; j < NCNC ; j++ ) {
	  lat[idx].O[mu][j] = (GLU_real)creal(u[j])+I*(GLU_real)cimag(u[j]) ;
	}
	u += NCNC ;
	size_t shift = lat[idx].back[mu] ;
	for( j = 0 ; j < NCNC ; j++ ) {
	  lat[shift].O[mu][j] = (GLU_real)creal(u[j])+I*(GLU_real)cimag(u[j]) ;
	}
	u += NCNC ;
      }
      #endif
    }
  }
  // close the file descriptor
  close( fd ) ;
  munmap( mm , sb.st_size ) ;

  *chksum = 0 ;

  return GLU_SUCCESS ; 
}

// writes my gauge field out in the CERN format using mmapd posix IO. Sadly performance isn't great
void
write_CLS_field_mmap( const struct site *__restrict lat ,
		      const char *outfile )
{
  if( ND != 4 ) {
    fprintf( stderr , "No known CERN format for Nd != 4\n" ) ;
    return ;
  }
#ifdef WORDS_BIGENDIAN
  fprintf( stderr , "[IO] CERN we don't support having a different file endianess yet\n" ) ;
  return ;
#endif

  // full size of the gauge fields + dimensions + plaquette
  const size_t size = LVOLUME*sizeof(double complex)*ND*NCNC + 4*sizeof(int) + sizeof(double) ;
  
  // open a file
  int fd = open( outfile , O_RDWR | O_CREAT , (mode_t)0600 ) ;
  if( fd == -1 ) {
    perror("open") ;
    exit(1) ;
  }
  if( ftruncate( fd , size ) == -1 ) exit(1) ;
  
  // now we map the file
  char *data = mmap( NULL , size, PROT_WRITE , MAP_SHARED , fd , 0 ) ;
  if( data == MAP_FAILED ) {
    perror( "mmap" ) ;
    exit(1) ;
  }
  
  // write the dimensions
  uint32_t *pt = (uint32_t*)data ;
  pt[0] = (uint32_t)Latt.dims[3] ; pt[1] = (uint32_t)Latt.dims[0] ;
  pt[2] = (uint32_t)Latt.dims[1] ; pt[3] = (uint32_t)Latt.dims[2] ;
  data += 4*sizeof(uint32_t) ;

  // write the plaquette
  double *plaq = (double*)data ;
  plaq[0] = NC*av_plaquette( lat ) ;
  data += sizeof(double) ;

  // finally point to the fields we will write
  double complex *uout = (double complex*)data ;

  // loop entire volume again in CLS order
  size_t i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    // CLS order
    int is_odd = 0 ;
    const size_t idx = CLS_to_GLU( &is_odd , i ) ;
    if( is_odd ) {
      field_to_out( uout , i , lat , idx ) ;
    }
  }
  
  // sync it all at the end
  msync( data , size , MS_SYNC ) ;
  munmap( data , size ) ;
  close(fd) ;
  
  return ;
}

// writes my gauge field out in the CERN format 
void
write_CLS_field( const struct site *__restrict lat ,
		 FILE *__restrict outfile )
{
  if( ND != 4 ) {
    fprintf( stderr , "No known CERN format for Nd != 4\n" ) ;
    return ;
  }
    
  uint32_t NAV[ ND ] = { Latt.dims[3] , Latt.dims[0] , Latt.dims[1] , Latt.dims[2] } ;
  if( WORDS_BIGENDIAN ) {
    bswap_32( ND , NAV ) ;
  }
  
  // first of all set calculate the plaquette checksum in header
  double plaq[1] = { 0 } ;

  plaq[0] = NC * av_plaquette( lat ) ;
  if( WORDS_BIGENDIAN ) {
    bswap_64( 1 , plaq ) ;
  }
  
  // write the checksums
  fwrite( NAV , ( ND ) * sizeof( uint32_t ) , 1 , outfile ) ;
  fwrite( plaq , sizeof( double ) , 1 , outfile ) ;

  // create a copy
  double complex *uout = NULL ;
  if( GLU_malloc( (void**)&uout , ALIGNMENT, LVOLUME*ND*NCNC*sizeof(double complex) ) != 0 ) {
    return ;
  }
  
  // the idea is to do all of the translation stuff in parallel and then dump it to a file
  size_t i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    // CLS order
    int is_odd = 0 ;
    const size_t idx = CLS_to_GLU( &is_odd , i ) ;
    if( is_odd ) {
      field_to_out( uout , i , lat , idx ) ;
    }
  }

  // write it out
  if( WORDS_BIGENDIAN ) {
    bswap_64( 8*LVOLUME*NCNC , uout ) ;
  }
  fwrite( uout , sizeof( double ) , 8*LVOLUME*NCNC , outfile ) ;

  if( uout != NULL ) {
    free( uout ) ;
  }
  
  return ;
}
