/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (GLU_rng.c) is part of GLU.

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
   @file GLU_rng.c
   @brief functions used in the generation of random numbers
 **/

#include "Mainfile.h"

#include "GLU_rng.h" // include self for alphabetising

// have we inited the rng?
static GLU_bool INIT_RNG = GLU_FALSE ;

// all the new rngs
#include "KISS_rng.h"
#include "MWC_1038.h"
#include "MWC_4096.h"
#include "well_rng_19937a.h"
#include "XOR1024.h"

#ifdef GSL_RNG
  // include for the GSL_RNG is in "GLU_definitions.h"
  static gsl_rng *r ;
#endif

// generate a random NCxNC matrix ... Use gaussian numbers?
void 
generate_NCxNC( GLU_complex U[ NCNC ] )
{
  int i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    // gaussians
    U[i] = polar_box() ;
  }
  return ;
}

// initialise the seed of the RNG by reading from /dev/urandom
int
initialise_seed( void )
{
  // should be initialised in input_reader, but just in case have a call in gsl_init too!
  if( Latt.Seed[0] == 0 ) {
    // if we haven't specified a Seed we go to the Unix entropy pool urandom
    // urandom provides the entropy pool value from random or a pseudorandom number if that is not available
    if( Latt.Seed[0] == 0 ) {
      FILE *urandom = fopen( "/dev/urandom" , "r" ) ;

      if( urandom == NULL ) {
	fprintf( stderr , "[RNG] /dev/urandom not opened!! ... Exiting \n" ) ;
	return GLU_FAILURE ;
      }

      const int check = fread( Latt.Seed , sizeof( Latt.Seed ) , 1 , urandom ) ;
      if( unlikely( check != 1 ) ) { 
	fprintf( stderr , "[RNG] Entropy pool Seed not read properly ! "
		 "... Exiting \n" ) ; 
	return GLU_FAILURE ;
      }

      fclose( urandom ) ;
    }
  }
  return GLU_SUCCESS ;
}

// Gaussian distributed doubles ocassionally called
GLU_real
polar( void )
{
  const GLU_real u = rng_dbl( ) , v = rng_dbl( ) ;
  return sqrt( -2. * log( u ) ) * cos( TWOPI * v ) ;
}

// looks more complicated , is faster.
GLU_complex
polar_box( void )
{ 
  register const GLU_real u = (GLU_real)( 2. * rng_dbl( ) - 1. ) ;
  register const GLU_real v = (GLU_real)( 2. * rng_dbl( ) - 1. ) ;
  const GLU_real s = u * u + v * v ;
  return s < 1. ? sqrt( -log( s ) / s ) * ( u + I * v ) : polar_box() ;
}

// Free the memory allocated ....
void 
rng_free( void )
{
  if( INIT_RNG == GLU_TRUE ) {
    #ifdef GSL_RNG
    gsl_rng_free( r ) ;
    #elif defined KISS_RNG
    GLU_free_KISS_table( ) ;
    #elif defined MWC_1038_RNG
    GLU_free_MWC1038_table( ) ;
    #elif defined MWC_4096_RNG
    GLU_free_MWC4096_table( ) ;
    #elif defined XOR1024_RNG
    GLU_free_XOR1024_table( ) ;
    #else
    GLU_free_WELL19937_table( ) ;
    #endif
  }
  return ;
}

// accessors for doubles
double
rng_dbl( void )
{
#ifdef MWC_1038_RNG
  return mwc_1038_dbl( ) ;
#elif defined KISS_RNG
  return KISS_dbl( ) ;
#elif defined GSL_RNG
  return gsl_rng_uniform( r ) ;
#elif defined MWC_4096_RNG
  return mwc_4096_dbl( ) ;
#elif defined XOR1024_RNG
  return XOR1024_dbl( ) ;
#else
  return WELLRNG19937a( ) ;
#endif
}

/// initialise the rng
void 
rng_init( void )
{
  // if we are already inited then we don't do it again ...
  if( INIT_RNG == GLU_FALSE ) {

    #ifdef GSL_RNG
    // mt seeded from Unix entropy pool //
    r = gsl_rng_alloc( gsl_rng_mt19937 ) ;
    gsl_rng_set( r , Latt.Seed[0] ) ;
    #else
    // probably a good case for some callbacks here //
      #if defined KISS_RNG
      GLU_set_KISS_table( Latt.Seed[0] ) ;
      #elif defined MWC_1038_RNG
      GLU_set_MWC1038_table( Latt.Seed[0] ) ;
      #elif defined MWC_4096_RNG
      GLU_set_MWC4096_table( Latt.Seed[0] ) ;
      #elif defined XOR1024_RNG
      GLU_set_XOR1024_table( Latt.Seed[0] ) ;
      #else
      GLU_set_WELL19937_table( Latt.Seed[0] ) ;
      #endif
    #endif

    // should discard the first 2000 or so (warm up phase)
    int i ;
    for( i = 0 ; i < 2E3 ; i++ ) { 
      rng_dbl() ;
    }

    INIT_RNG = GLU_TRUE ;
  }
  return ;
}

// accessor for ints
uint32_t
rng_int( void )
{
  return (uint32_t)( UINT32_MAX * rng_dbl() ) ;
}
