/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (XOR1024.c) is part of GLU.

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

#include <stdint.h>
#include <stdlib.h>

static uint64_t *table ;

// hmmmm, 64-bit generator seeded with 32-bit?
void
GLU_set_XOR1024_table( const uint32_t seed )
{
  table = ( uint64_t* )malloc( 16 * sizeof( uint64_t ) ) ;
  table[0] = ( uint64_t )seed ;
  int i ;
  for( i = 1 ; i < 15 ; i++ ) {
    table[i] = ( 1812433253UL * ( table[i-1] ^ ( table[i-1] >> 30)) + i ) ;
  }
  return ;
}

// free the allocated table
void
GLU_free_XOR1024_table( void )
{
  free( table ) ;
}

int p;
static uint64_t next( void ) 
{
  uint64_t s0 = table[ p ];
  uint64_t s1 = table[ p = ( p + 1 ) & 15 ];
  s1 ^= s1 << 31; // a
  s1 ^= s1 >> 11; // b
  s0 ^= s0 >> 30; // c
  return ( table[ p ] = s0 ^ s1 ) * 1181783497276652981LL;
}

double
XOR1024_dbl( void ) 
{
  return next( ) * 5.421010862427522e-20 ;
}
