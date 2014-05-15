/*
    Copyright 2013 Renwick James Hudspith

    This file (MWC_4096.h) is part of GLU.

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
   @file MWC_4096.h
   @brief prototype functions for Marsaglia's massive period generator
 */

#ifndef GLU_MWC_4096_H
#define GLU_MWC_4096_H

/**
   @fn void GLU_set_MWC4096_table( const uint32_t seed )
   @brief allocate the MWC4096 table
   @param seed :: single seed value
   
   Uses the lagged fibonacci generator of Knuth to fill the table from the seed
 */
void
GLU_set_MWC4096_table( const uint32_t seed ) ;

/**
   @fn void GLU_free_MWC4096_table( void )
   @brief deallocate the MWC4096 table
 */
void
GLU_free_MWC4096_table( void ) ;

/**
   @fn double mwc_4096_dbl( )
   @brief generates a double-precision random number from the long unsigned int
 */
double
mwc_4096_dbl( void ) ;

#endif
