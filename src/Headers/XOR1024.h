/*
    Copyright 2013 Renwick James Hudspith

    This file (XOR1024_rng.h) is part of GLU.

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
   @file XOR1024.h
   @brief function prototypes for the 64-bit XOR1024 RNG
 **/

#ifndef GLU_XOR1024_RNG_H
#define GLU_XOR1024_RNG_H

/**
   @fn void GLU_set_XOR1024_table( const uint32_t seed )
   @brief allocate the RNG table (is 16 uint32_t's)
   @param seed :: single RNG seed
 */
void
GLU_set_XOR1024_table( const uint32_t seed ) ;

/**
   @fn void GLU_free_XOR1024_table( void )
   @brief free the allocated RNG table
 */
void
GLU_free_XOR1024_table( void ) ;

/**
   @fn double XOR1024_dbl( void ) 
   @brief randomly generated double
   @return a double precision random number in the interval [0,1)
 */
double
XOR1024_dbl( void ) ;

#endif 
