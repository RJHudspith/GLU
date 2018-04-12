/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (plan_ffts.c) is part of GLU.

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
   @file plan_ffts.c
   @brief planner for the FFTS we use everywhere
 */
#include "Mainfile.h"

#include "GLU_timer.h" // tells us how long we spent planning FFTs

// guard the whole thing
#ifdef HAVE_FFTW3_H

// I have reason to believe that this gives better results
// FFTW_PATIENT OVERWRITES the data in the arrray given, be warned!
#define GLU_PLAN FFTW_PATIENT

// for ease of reading
enum{ NOPLAN = 0 } ;

// see if we have wisdom already
static char *
obtain_wisdom( int *planflag ,
	       const size_t dims[ ND ] ,
	       const int DIR , 
	       const char *type )
{
  *planflag = NOPLAN ;
  char *str = malloc( 256 * sizeof( char ) ) ;
  sprintf( str , "%s" , type ) ;
  if( DIR < 0 ) return str ;
#ifndef CONDOR_MODE
  FILE *wizzard ;
  size_t mu ;
  char prec_str[ 16 ] ; 
  *planflag = NOPLAN ; 
  #ifdef SINGLE_PREC
  sprintf( prec_str , "FLOAT" ) ;
  #else
  sprintf( prec_str , "DOUBLE" ) ;
  #endif
  sprintf( str , "%s/Local/Wisdom/%s_%sSU%d_" , 
	   HAVE_PREFIX , prec_str , type , NC ) ;
  for( mu = 0 ; mu < DIR - 1 ; mu++ ) {
    sprintf( str , "%s%zux" , str , dims[ mu ] ) ;
  }
  sprintf( str , "%s%zu.wisdom" , str , dims[ DIR - 1 ] ) ;
  if( ( wizzard = fopen( str , "r" ) ) == NULL ) {
    fprintf( stdout , "\n[FFTW] No wisdom to be obtained here ... planning\n" ) ; 
  } else {
    #ifdef verbose
    fprintf( stdout , "\n[FFTW] Successful wisdom attained\n" ) ;
    #endif
    *planflag = fftw_import_wisdom_from_file( wizzard ) ; 
    fclose( wizzard ) ; 
  }
  // condor mode ifdef
#else
  fprintf( stdout , "[FFTW] Creating plan on CONDOR host\n" ) ; 
#endif
  return str ;
}

// record both the forward and backward
void
create_plans_DFT( struct fftw_stuff *FFTW ,
		  const size_t dims[ ND ] ,
		  const int ARR_SIZE ,
		  const int DIR )
{
  // set up our fft
  size_t VOL = 1 , mu , i ;
  int dimes[ DIR ] , planflag ;
  
  // swap these defs around as FFTW and I disagree on orderings
  for( mu = 0 ; mu < DIR ; mu++ ) {
    dimes[ mu ] = dims[ DIR - 1 - mu ] ;
    VOL *= dims[ DIR - 1 - mu ] ;
  }

  #ifdef verbose
  start_timer( ) ;
  #endif

  FFTW -> forward  = malloc( ARR_SIZE * sizeof( fftw_plan ) ) ; 
  FFTW -> backward = malloc( ARR_SIZE * sizeof( fftw_plan ) ) ; 
  FFTW -> out      = fftw_malloc( ARR_SIZE * sizeof( GLU_complex* ) ) ; 
  FFTW -> in       = fftw_malloc( ARR_SIZE * sizeof( GLU_complex* ) ) ; 

#pragma omp parallel
  {
#pragma omp for private(i)
    for(  i = 0 ; i < ARR_SIZE ; i++  ) {
      FFTW -> out[i] = ( GLU_complex* )fftw_malloc( VOL * sizeof( GLU_complex ) ) ; 
      FFTW -> in[i] = ( GLU_complex* )fftw_malloc( VOL * sizeof( GLU_complex ) ) ; 
    }
  }
  
  char *str = obtain_wisdom( &planflag , dims , DIR , "" ) ;

  for( mu = 0 ; mu < ARR_SIZE ; mu++ ) {
    FFTW -> forward[mu] = fftw_plan_dft( DIR , dimes ,
					 FFTW -> in[mu] , FFTW -> out[mu] , 
					 FFTW_FORWARD , GLU_PLAN ) ; 
    FFTW -> backward[mu] = fftw_plan_dft( DIR , dimes ,
					  FFTW -> out[mu] , FFTW -> in[mu] , 
					  FFTW_BACKWARD , GLU_PLAN ) ;
  }

  // I want to know how long FFTW is taking to plan its FFTs
  #ifdef verbose
  print_time( ) ;
  fprintf( stdout , "[FFTW] plans finished\n\n" ) ;
  #endif

#ifndef CONDOR_MODE
  if( planflag == NOPLAN )  {
    FILE *wizzard = fopen( str , "w" ) ; 
    fftw_export_wisdom_to_file( wizzard ) ; 
    fclose( wizzard ) ; 
  }
#endif
  free( str ) ;

  return ;
}

/// Small plan does not care about the NC unlike the above
void
small_create_plans_DFT( struct fftw_small_stuff *FFTW ,
			const size_t dims[ ND ] ,
			const int DIR )
{
  // set up our fft
  size_t VOL = 1 , mu ;
  int dimes[ DIR ] , planflag ;
  // swap these defs around
  for( mu = 0 ; mu < DIR ; mu++ ) {
    dimes[ mu ] = dims[ DIR - 1 - mu ] ;
    VOL *= dims[ DIR - 1 - mu ] ;
  }

  // allocations
  FFTW -> in = fftw_malloc( VOL * sizeof( GLU_complex ) ) ;
  FFTW -> out = fftw_malloc( VOL * sizeof( GLU_complex ) ) ;

  // initialise the clock
  #ifdef verbose
  start_timer( ) ;
  #endif

  char *str = obtain_wisdom( &planflag , dims , DIR , "single_" ) ;

  FFTW -> forward = fftw_plan_dft( DIR , dimes ,
				   FFTW -> in , FFTW -> out , 
				   FFTW_FORWARD , GLU_PLAN ) ; 
  FFTW -> backward = fftw_plan_dft( DIR , dimes ,
				    FFTW -> out , FFTW -> in , 
				    FFTW_BACKWARD , GLU_PLAN ) ; 

  // I want to know how long FFTW is taking to plan its FFTs
  #ifdef verbose
  print_time( ) ;
  fprintf( stdout , "[FFTW] plans finished\n\n" ) ;
  #endif

#ifndef CONDOR_MODE
  if( planflag == NOPLAN ) {
    FILE *wizzard = fopen( str , "w" ) ; 
    fftw_export_wisdom_to_file( wizzard ) ; 
    fclose( wizzard ) ; 
  }
#endif
  free( str ) ;

  return ;
}

void
clean_up_fftw( struct fftw_stuff FFTW ,
	       const int ARR_SIZE )
{
  size_t i ;
  for( i = 0 ; i < ARR_SIZE; i++ ) {
    fftw_destroy_plan( FFTW.forward[i] ) ;   
    fftw_destroy_plan( FFTW.backward[i] ) ;  
    fftw_free( FFTW.out[i] ) ; 
    fftw_free( FFTW.in[i] ) ; 
  }
  fftw_free( FFTW.out ) ; 
  fftw_free( FFTW.in ) ; 
  free( FFTW.forward ) ; 
  free( FFTW.backward ) ; 
  fftw_cleanup( ) ;
  if( FFTW.psq != NULL ) {
    free( FFTW.psq ) ; 
  }
  return ;
}

void
small_clean_up_fftw( struct fftw_small_stuff FFTW )
{
  fftw_free( FFTW.out ) ; 
  fftw_free( FFTW.in ) ; 
  fftw_cleanup( ) ;
  if( FFTW.psq != NULL ) {
    free( FFTW.psq ) ; 
  }
  return ;
}

// and clean this up
#ifdef GLU_PLAN
  #undef GLU_PLAN
#endif

#endif
