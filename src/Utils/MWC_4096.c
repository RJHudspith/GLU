/*
    Copyright 2013 Renwick James Hudspith

    This file (MWC_4096.c) is part of GLU.

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
   @file MWC_4096.c
   @brief Marsaglia's ultra-long period complementary-multiply-with carry generator
 */

#include <stdint.h> // for uint32_t and uint64_t
#include <stdlib.h> // for malloc in a function

// uncomment for a better double
//#define FULL_DOUBLE 

static uint32_t *table ;
static uint64_t carry ; // carry
static int mcwc_i ;     // index

// set the allocated table
void
GLU_set_MWC4096_table( const uint32_t seed )
{
  table = ( uint32_t* )malloc( 4096 * sizeof( uint32_t ) ) ;
  table[0] = seed ;
  int i ;
  for( i = 1 ; i < 4096 ; i++ ) {
    table[i] = ( 1812433253UL * ( table[i-1] ^ ( table[i-1] >> 30)) + i ) ;
  }
  // set the carry and the mwc index
  carry = table[4095] % 809430660UL ;
  mcwc_i = 4096 ;
  return ;
}

// free the allocated table
void
GLU_free_MWC4096_table( void )
{
  free( table ) ;
}

// Returns a uint32_t
static uint32_t
MWC4096( void )
{
  //mcwc_i = 4096 ;
  mcwc_i = ( mcwc_i + 1 )&4095 ;
  const uint64_t t = 18782ULL*table[ mcwc_i ] + carry ; 
  uint32_t x = t + carry ;
  if( x < carry ) {
    x++ ;
    carry++ ;
  }
  return ( table[mcwc_i] = 0xFFFFFFFEU-x ) ;
}

// generate a double between [ 0 , 1 )
double 
mwc_4096_dbl( void )
{
#ifdef FULL_DOUBLE
  // to remove the type-punning error gcc always gives about strict aliasing
  union {
    uint32_t theInts[2] ;
    uint64_t ulint ;
    double x ;
  } lu2dbl ;
  // pack the bottom and top half of an unsigned long with uints
  lu2dbl.theInts[1] = MWC4096( ) ;
  lu2dbl.theInts[0] = MWC4096( ) ;
  // shift down so that only the mantissa has non-zero entries
  lu2dbl.ulint = ( lu2dbl.ulint >> 12 ) | 0x3FF0000000000000UL ;
  // have to subtract one because of the implicit 2^0+mantissa
  return lu2dbl.x - 1.0 ;
#else
  return (double)( MWC4096( ) * 2.3283064365386963e-10 ) ;
#endif
}

