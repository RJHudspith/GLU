/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (well_rng_19937a.h) is part of GLU.

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
   @file well_rng_19937a.h
   @brief prototype functions for the WELL rng
   Is our default, because it is Jamie's favourite
   @warning does not play nice with openMP, be careful
 */

#ifndef GLU_WELL_RNG_19937A_H
#define GLU_WELL_RNG_19937A_H

/* ***************************************************************************** */
/* Copyright:      Francois Panneton and Pierre L'Ecuyer, University of Montreal */
/*                 Makoto Matsumoto, Hiroshima University                        */
/* Notice:         This code can be used freely for personal, academic,          */
/*                 or non-commercial purposes. For commercial purposes,          */
/*                 please contact P. L'Ecuyer at: lecuyer@iro.UMontreal.ca       */
/* ***************************************************************************** */

/**
   @fn void GLU_set_WELL19937_table( const uint32_t seed )
   @brief allocate the well RNG table
   @param seed :: initial seed value

   Uses the lagged fibonacci of Knuth to fill the table based on the seed
 */
void
GLU_set_WELL19937_table( const uint32_t seed ) ;

/**
   @fn void GLU_free_WELL19937_table( void )
   @brief free the table
 */
void
GLU_free_WELL19937_table( void ) ;

/**
   @fn extern double (*WELLRNG19937a) ( void )
   @returns a double precision random number in [ 0 , 1 )
 */
extern double
( *WELLRNG19937a ) ( void ) ;

#endif

