/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (par_KISS.c) is part of GLU.

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
   @file par_KISS.c
   @brief parallel KISS rng
 */
#include "Mainfile.h"

#include "GLU_bswap.h"

// rng table
static uint32_t **table ; 

// unsigned integer version
static uint32_t
JKISS32( const uint32_t thread )
{
  int t ;
  table[thread][1] ^= ( table[thread][1] << 5 ) ;
  table[thread][1] ^= ( table[thread][1] >> 7 ) ;
  table[thread][1] ^= ( table[thread][1] << 22 ) ;
  t = table[thread][2] + table[thread][3] + table[thread][4] ;
  table[thread][2] = table[thread][3] ;
  table[thread][4] = t < 0 ;
  table[thread][3] = t&2147483647 ;
  table[thread][0] += 1411392427 ;
  return table[thread][0] + table[thread][1] + table[thread][3] ;
}

void
free_par_KISS( void )
{
  size_t i ;
  for( i = 0 ; i < Latt.Nthreads ; i++ ) {
    free( table[i] ) ;
  }
  free( table ) ;
}

void
GLU_set_par_KISS_table( const uint32_t seed[ Latt.Nthreads ] )
{
  table = malloc( Latt.Nthreads * sizeof( uint32_t* ) ) ;
  size_t i , j ;
  for( i = 0 ; i < Latt.Nthreads ; i++ ) {
    table[i] = malloc( RNG_TABLE * sizeof( uint32_t ) ) ;
    table[i][0] = seed[i] ;
    for( j = 1 ; j < RNG_TABLE ; j++ ) {
      table[i][j] = ( 1812433253UL * ( table[i][j-1] ^ ( table[i][j-1] >> 30)) + j ) ;
    }
  }
  return ;
}

// generates a random double number
double
par_KISS_dbl( const uint32_t thread ) 
{
#ifdef FULL_DOUBLE
  // to remove the type-punning error gcc always gives about strict aliasing
  union {
    uint32_t theInts[2] ;
    uint64_t ulint ;
    double x ;
  } lu2dbl ;
  // pack the bottom and top half of an unsigned long with uints
  lu2dbl.theInts[0] = JKISS32( thread ) ;
  lu2dbl.theInts[1] = JKISS32( thread ) ;
  // shift down so that only the mantissa has non-zero entries
  lu2dbl.ulint = ( lu2dbl.ulint >> 12 ) | 0x3FF0000000000000ULL ;
  // have to subtract one because of the implicit 2^0+mantissa
  return lu2dbl.x - 1.0 ;
#else
  return JKISS32( thread ) * 2.3283064365386963e-10 ;
#endif
}

// read the MWC table
int
read_par_KISS_table( FILE *rng_file )
{
  table  = malloc( Latt.Nthreads * sizeof( uint32_t* ) ) ;
  size_t i ;
  for( i = 0 ; i < Latt.Nthreads ; i++ ) {
    table[i] = malloc( RNG_TABLE * sizeof( uint32_t ) ) ;
    if( fread( table[i] , sizeof( uint32_t ) , RNG_TABLE , rng_file ) != RNG_TABLE ) {
      fprintf( stderr , "[PAR_RNG] file read failure\n" ) ;
      return GLU_FAILURE ;
    }
  }
  return GLU_SUCCESS ;
}

// write the table and periferies
int
write_par_KISS_table( FILE *rng_file )
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
  return GLU_SUCCESS ;
}
