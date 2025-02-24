/*
Copyright 2013-2025 Renwick James Hudspith

    This file (GLU_rng.h) is part of GLU.

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
   @file GLU_rng.h
   @brief prototype functions for the random number generation
 */
#ifndef GLU_RNG_H
#define GLU_RNG_H

/**
   @fn void generate_NCxNC( GLU_complex U[ NCNC ] )
   @brief generate a random NCxNC matrix calls polar()
 */
void 
generate_NCxNC( GLU_complex U[ NCNC ] ) ;

/**
   @fn int initialise_seed( void )
   @brief initialises the original seed for the RNG
   @warning attempts to read from /dev/urandom !!
   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int
initialise_seed( void ) ;

/**
   @fn GLU_real polar( void )
   @brief Box-Mueller algorithm for uniformly distributed real gaussian numbers 

   uses the slower cosine definition for the box-mueller calls rng_dbl()
 */
GLU_real 
polar( void ) ;

/**
   @fn GLU_complex polar_box( void )
   @brief Box-Mueller algorithm for uniformly distributed complex gaussian numbers 
   Uses the faster polar Box-Mueller algorithm for generating the numbers plus recursion calls rng_dbl()
 **/
GLU_complex
polar_box( void ) ;

/**
   @fn double rng_dbl( void )
   @brief return a double in [0,1)
 */
double
rng_dbl( void ) ;

/**
   @fn void rng_free( void )
   @brief free the RNG

   If the rng hasn't been set up it does nothing
 */
void 
rng_free( void ) ;

/**
   @fn void rng_init( void )
   @brief initialise the RNG

   This code uses the mersenne twister algorithm 
   <a href="linkURL">http://en.wikipedia.org/wiki/Mersenne_twister </a> if #GSL_RNG is defined
   Using gsl's implementation 
   <a href="linkURL">http://www.gnu.org/software/gsl/ </a>

   <br><br>
   Otherwise we use the WELL algorithm copied pretty much verbatim from <a href="linkURL">http://www.iro.umontreal.ca/~panneton/WELLRNG.html</a>

   <br><br>
   Or we use one of Marsaglia's algorithms. The KISS or the complementary multiply with carry MWC_1038 or MWC_4096.
 **/
void 
rng_init( void ) ;

/**
   @fn unsigned int rng_int( void )
   @brief return an unsigned int

   sometimes we just want an integer calls rng_dbl()
 */
unsigned int
rng_int( void ) ;

#endif
