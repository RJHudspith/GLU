/*
Copyright 2013-2025 Renwick James Hudspith

    This file (par_rng.h) is part of GLU.

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
   @file par_rng.h
   @brief prototype functions for the parallel rng
 */
#ifndef PAR_RNG_H
#define PAR_RNG_H

/**
   @fn void free_par_rng( void )
   @brief free the rng table
 */
void
free_par_rng( void ) ;

/**
   @fn int initialise_par_rng( const char *rng_file )
   @brief initialise the rng
   @param rng_file :: file to read rng state or NULL if we want to generate one
   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int
initialise_par_rng( const char *rng_file ) ;

/**
   @fn void par_generate_NCxNC( GLU_complex U[ NCNC ] , const uint32_t thread )
   @brief generate a #NC x #NC matrix of gaussian random numbers
   @param U :: matrix of random numbers
   @param thread :: parallel thread index
 */
void 
par_generate_NCxNC( GLU_complex U[ NCNC ] , 
		    const uint32_t thread ) ;

/**
   @fn GLU_real par_polar( const uint32_t thread )
   @brief gaussian random numbers from the box-mueller algorithm
   @param thread :: parallel thread index
   @return a real gaussian random number with distribution sigma of 1
 */
GLU_real
par_polar( const uint32_t thread ) ;

/**
   @fn GLU_complex par_polar_box( const uint32_t thread )
   @brief gaussian random numbers from the polar box-mueller algorithm
   @param thread :: parallel thread index
   @return a gaussian distributed complex number
 */
GLU_complex
par_polar_box( const uint32_t thread ) ;

/**
   @fn double par_rng_dbl( const uint32_t thread )
   @brief a uniform double precision number distributed between 0 and 1
   @param thread :: paralllel thread index
   @return a random double
 */
double
par_rng_dbl( const uint32_t thread ) ;

/**
   @fn uint32_t par_rng_int( const uint32_t thread )
   @brief a random uint32_t between 0 and UINT_MAX
   @param thread :: parallel thread index
   @return a uint32_t
 */
uint32_t
par_rng_int( const uint32_t thread ) ;

/**
   @fn int read_par_rng_state( const char *infile )
   @brief read the rng table from infile
   @param infile :: file to be read
   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int
read_par_rng_state( const char *infile ) ;

/**
   @fn int read_par_rng_state( const char *outfile )
   @brief read the rng table to outfile
   @param outfile :: file to be written
*/
void
write_par_rng_state( const char *outfile ) ;

#endif
