/*
Copyright 2013-2025 Renwick James Hudspith

    This file (par_MWC_4096.c) is part of GLU.

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
   @file par_MWC_4096.c
   @brief stupidly long period generator
 */
#include "Mainfile.h"

#include "GLU_bswap.h" // bswap32 and bswap64

static uint32_t **table ;
static uint32_t *carry ;  // carry
static uint32_t *mcwc_i ; // index

// Returns a uint32_t
static uint32_t
par_MWC_4096( const uint32_t thread )
{
  mcwc_i[thread] = 4095 ;
  mcwc_i[thread] = ( mcwc_i[thread]+1 )&4095 ;

  const uint64_t t = 18782LL*table[thread][mcwc_i[thread]] + carry[thread] ;
  carry[thread] = ( t >> 32 ) ; 
  uint32_t x = t + carry[thread] ; 
  if( x < carry[thread] ){
    x++ ; 
    carry[thread]++ ;
  }
  return( table[thread][mcwc_i[thread]] = 0xfffffffe - x ) ;    
}

// free the table
void
free_par_MWC_4096( void )
{
  size_t i ;
  for( i = 0 ; i < Latt.Nthreads ; i++ ) {
    if( table[i] != NULL ) {
      free( table[i] ) ;
    }
  }
  if( mcwc_i != NULL ) {
    free( mcwc_i ) ;
  }
  if( carry != NULL ) {
    free( carry ) ;
  }
  if( table != NULL ) {
    free( table ) ;
  }
  return ;
}

// set the allocated table
void
GLU_set_par_MWC_4096_table( const uint32_t seed[ Latt.Nthreads ] )
{
  table = malloc( RNG_TABLE * sizeof( uint32_t* ) ) ;
  mcwc_i = malloc( Latt.Nthreads * sizeof( uint32_t ) ) ;
  carry  = malloc( Latt.Nthreads * sizeof( uint32_t ) ) ;
  size_t i , j ;
  for( i = 0 ; i < Latt.Nthreads ; i++ ) {
    table[i] = malloc( RNG_TABLE * sizeof( uint32_t ) ) ;
    table[i][0] = seed[i] ;
    for( j = 1 ; j < 4096 ; j++ ) {
      table[i][j] = ( 1812433253UL * ( table[i][j-1] ^ ( table[i][j-1] >> 30)) + j ) ;
    }
    // set the carry and the mwc index
    carry[i]  = table[i][RNG_TABLE-1] % 809430660UL ;
    mcwc_i[i] = RNG_TABLE - 1 ;
  }
  return ;
}

// generate a double between [ 0 , 1 )
double 
par_MWC_4096_dbl( const uint32_t thread )
{
#ifdef FULL_DOUBLE
  // to remove the type-punning error gcc always gives about strict aliasing
  union {
    uint32_t theInts[2] ;
    uint64_t ulint ;
    double x ;
  } lu2dbl ;
  // pack the bottom and top half of an unsigned long with uints
  lu2dbl.theInts[1] = par_MWC_4096( thread ) ;
  lu2dbl.theInts[0] = par_MWC_4096( thread ) ;
  // shift down so that only the mantissa has non-zero entries
  lu2dbl.ulint = ( lu2dbl.ulint >> 12 ) | 0x3FF0000000000000UL ;
  // have to subtract one because of the implicit 2^0+mantissa
  return lu2dbl.x - 1.0 ;
#else
  return (double)( par_MWC_4096( thread ) * 2.3283064365386963e-10 ) ;
#endif
}

// read the MWC table
int
read_par_MWC_4096_table( FILE *rng_file )
{
  table  = malloc( Latt.Nthreads * sizeof( uint32_t* ) ) ;
  mcwc_i = malloc( Latt.Nthreads * sizeof( uint32_t ) ) ;
  carry  = malloc( Latt.Nthreads * sizeof( uint32_t ) ) ;
  size_t i ;
  for( i = 0 ; i < Latt.Nthreads ; i++ ) {
    table[i] = malloc( RNG_TABLE * sizeof( uint32_t ) ) ;
    if( fread( table[i] , sizeof( uint32_t ) , RNG_TABLE , rng_file ) != RNG_TABLE ) {
      fprintf( stderr , "[PAR_RNG] file read failure\n" ) ;
      return GLU_FAILURE ;
    }
  }
  if( fread( mcwc_i , sizeof( uint32_t ) , Latt.Nthreads , rng_file ) != Latt.Nthreads ) {
    fprintf( stderr , "[PAR_RNG] file read failure\n" ) ;
    return GLU_FAILURE ;
  }
  if( fread( carry , sizeof( uint32_t ) , Latt.Nthreads , rng_file ) != Latt.Nthreads ) {
    fprintf( stderr , "[PAR_RNG] file read failure\n" ) ;
    return GLU_FAILURE ;
  }
  // swap them if we aren't big endian
  if( !WORDS_BIGENDIAN ) { 
    for( i = 0 ; i < Latt.Nthreads ; i++ ) {
      bswap_32( RNG_TABLE , table[i] ) ;
    }
    bswap_32( Latt.Nthreads , mcwc_i ) ;
    bswap_32( Latt.Nthreads , carry ) ;
  }
  return GLU_SUCCESS ;
}

// write the table and periferies
void
write_par_MWC_4096_table( FILE *rng_file )
{
  // write out big endian binary data
  size_t i ;
  for( i = 0 ; i < Latt.Nthreads ; i++ ) {
    // put into big endian
    if( !WORDS_BIGENDIAN ) { 
      bswap_32( RNG_TABLE , table[i] ) ; 
    }
    // write the table
    fwrite( table[i] , sizeof( uint32_t ) , RNG_TABLE , rng_file ) ;
    // reswap
    if( !WORDS_BIGENDIAN ) { 
      bswap_32( RNG_TABLE , table[i] ) ; 
    }
  }
  // write out the carries
  if( !WORDS_BIGENDIAN ) { 
    bswap_32( Latt.Nthreads , mcwc_i ) ;
    bswap_32( Latt.Nthreads , carry ) ; 
  }
  fwrite( mcwc_i , sizeof( uint32_t ) , Latt.Nthreads , rng_file ) ;
  fwrite( carry , sizeof( uint32_t ) , Latt.Nthreads , rng_file ) ;
  if( !WORDS_BIGENDIAN ) { 
    bswap_32( Latt.Nthreads , mcwc_i ) ;
    bswap_32( Latt.Nthreads , carry ) ; 
  }

  return ;
}
