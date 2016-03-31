/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (KISS_rng.h) is part of GLU.

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
   @file KISS_rng.h
   @brief function prototypes for the KISS rng by Marsaglia

   Makes sure seed_x , seed_y , seed_w < 2^32 and seed_z < 698769068
 **/

#ifndef GLU_KISS_RNG_H
#define GLU_KISS_RNG_H

/**
   @fn void GLU_set_KISS_table( const uint32_t seed )
   @brief allocate the RNG table (is 5 uint32_t's)
   @param seed :: single RNG seed
 */
void
GLU_set_KISS_table( const uint32_t seed ) ;

/**
   @fn void GLU_free_KISS_table( void )
   @brief free the allocated RNG table
 */
void
GLU_free_KISS_table( void ) ;

/**
   @fn void seed_rand_JKISS32( const uint32_t table[5] )
   @brief seed the JKISS32 RNG
   @param table of 5 seed values
 */
void
seed_rand_JKISS32( const uint32_t table[5] ) ;

/**
   @fn double KISS_dbl( void ) 
   @brief randomly generated double
   @return a double precision random number in the interval [0,1)
 */
double
KISS_dbl( void ) ;

#endif 
