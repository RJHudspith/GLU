/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (par_WELL_512.c) is part of GLU.

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
   @file par_WELL_512.c
   @brief well rng changed from public domain code http://lomont.org/Math/Papers/2008/Lomont_PRNG_2008.pdf
 */
#include "Mainfile.h"

#include "GLU_bswap.h" // bswap32

static uint32_t **table ;
static uint32_t *well_i ;

// returns a 32bit unsigned int using the Well
static uint32_t
par_WELL_512( const uint32_t thread )
{
  uint32_t *index = ( well_i + thread ) ;
  uint32_t *t = table[ thread ] ;

  register uint32_t a  = t[ *index ] ;
  register uint32_t c  = t[( *index + 13 )&15 ] ;
  const register uint32_t b  = a^c^(a<<16)^(c<<15) ;
  c  = t[ (*index+9 )&15 ] ;
  c ^= (c>>11) ;
  a  = t[ *index ] = b^c ; 
  const register uint32_t d  = a^((a<<5)&0xDA442D24UL) ;

  *index = ( *index + 15)&15;
  a = t[ *index ];
  return ( t[ *index ] = a^b^d^(a<<2)^(b<<18)^(c<<28) ) ;
}

// free the table
void
free_par_WELL_512( void )
{
  size_t i ;
  for( i = 0 ; i < Latt.Nthreads ; i++ ) {
    free( table[i] ) ;
  }
  free( table ) ;
  free( well_i ) ;
}

// set the table
void
GLU_set_par_WELL_512_table( const uint32_t seed[ Latt.Nthreads ] )
{
  table  = malloc( Latt.Nthreads * sizeof( uint32_t* ) ) ;
  well_i = malloc( Latt.Nthreads * sizeof( uint32_t ) ) ;
  size_t i , j ;
  for( i = 0 ; i < Latt.Nthreads ; i++ ) {
    table[i] = malloc( RNG_TABLE * sizeof( uint32_t ) ) ;
    table[i][0] = seed[i] ;
    for( j = 1 ; j < RNG_TABLE ; j++ ) {
      table[i][j] = ( 1812433253UL * ( table[i][j-1] ^ ( table[i][j-1] >> 30)) + j ) ;
    }
    // initialise the state to the first index
    well_i[i] = 0 ;
  }
  return ;
}

// generate a double between [ 0 , 1 )
double 
par_WELL_512_dbl( const uint32_t thread )
{
#ifdef FULL_DOUBLE
  // to remove the type-punning error gcc always gives about strict aliasing
  union {
    uint32_t theInts[2] ;
    uint64_t ulint ;
    double x ;
  } lu2dbl ;
  // pack the bottom and top half of an unsigned long with uints
  lu2dbl.theInts[1] = par_WELL_512( thread ) ;
  lu2dbl.theInts[0] = par_WELL_512( thread ) ;
  // shift down so that only the mantissa has non-zero entries
  lu2dbl.ulint = ( lu2dbl.ulint >> 12 ) | 0x3FF0000000000000UL ;
  // have to subtract one because of the implicit 2^0+mantissa
  return lu2dbl.x - 1.0 ;
#else
  return (double)( par_WELL_512( thread ) * 2.3283064365386963e-10 ) ;
#endif
}


// read the MWC table
int
read_par_WELL_512_table( FILE *rng_file )
{
  table  = malloc( Latt.Nthreads * sizeof( uint32_t* ) ) ;
  well_i = malloc( Latt.Nthreads * sizeof( uint32_t ) ) ;
  size_t i ;
  for( i = 0 ; i < Latt.Nthreads ; i++ ) {
    table[i] = malloc( RNG_TABLE * sizeof( uint32_t ) ) ;
    if( fread( table[i] , sizeof( uint32_t ) , RNG_TABLE , rng_file ) != RNG_TABLE ) {
      fprintf( stderr , "[PAR_RNG] file read failure\n" ) ;
      return GLU_FAILURE ;
    }
  }
  if( fread( well_i , sizeof( uint32_t ) , Latt.Nthreads , rng_file ) != Latt.Nthreads ) {
    fprintf( stderr , "[PAR_RNG] file read failure\n" ) ;
    return GLU_FAILURE ;
  }
  // swap them if we aren't big endian
  if( !WORDS_BIGENDIAN ) { 
    for( i = 0 ; i < Latt.Nthreads ; i++ ) {
      bswap_32( RNG_TABLE , table[i] ) ;
    }
    bswap_32( Latt.Nthreads , well_i ) ;
  }
  return GLU_SUCCESS ;
}

// write the table and periferies
void
write_par_WELL_512_table( FILE *rng_file )
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
    bswap_32( Latt.Nthreads , well_i ) ;
  }
  fwrite( well_i , sizeof( uint32_t ) , Latt.Nthreads , rng_file ) ;
  if( !WORDS_BIGENDIAN ) { 
    bswap_32( Latt.Nthreads , well_i ) ;
  }
  return ;
}
