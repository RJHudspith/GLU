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

// I have reason to believe that this gives better results //FFTW_PATIENT
// FFTW_PATIENT OVERWRITES the data in the arrray given, be warned!
#define GLU_PLAN FFTW_PATIENT

// for ease of reading
enum{ NOPLAN = 0 } ;

// have a look for parallel ffts
int
parallel_ffts( void )
{
  // initialise parallel fft ?
#ifdef OMP_FFTW
  if( fftw_init_threads() == 0 ) {
    return GLU_FAILURE ;
  } else {
    // in here I set the number of fftw threads to be the same
    // as the usual number of parallel threads ..
    fftw_plan_with_nthreads( Latt.Nthreads ) ;
    fprintf( stdout , "[PAR] FFTW using %d thread(s) \n" , Latt.Nthreads ) ;
  }
#endif
  return GLU_SUCCESS ;
}

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
  // openmp'd wisdom
  #ifdef OMP_FFTW
  sprintf( str , "%s/Local/Wisdom/%s_%sOMPnt%d_SU%d_" , 
	   HAVE_PREFIX , prec_str , type , nthreads , NC ) ;
  #else
  sprintf( str , "%s/Local/Wisdom/%s_%sSU%d_" , 
	   HAVE_PREFIX , prec_str , type , NC ) ;
  #endif
  for( mu = 0 ; mu < DIR - 1 ; mu++ ) {
    sprintf( str , "%s%zux" , str , dims[ mu ] ) ;
  }
  sprintf( str , "%s%zu.wisdom" , str , dims[ DIR - 1 ] ) ;
  if( ( wizzard = fopen( str , "r" ) ) == NULL ) {
    fprintf( stdout , "\n[FFTW] No wisdom to be obtained here ... planning" ) ; 
  } else {
    #ifdef verbose
    fprintf( stdout , "\n[FFTW] Successful wisdom attained" ) ;
    #endif
    *planflag = fftw_import_wisdom_from_file( wizzard ) ; 
    fclose( wizzard ) ; 
  }
  // condor mode ifdef
#else
  fprintf( stdout , "[FFTW] Creating plan on CONDOR host" ) ; 
#endif
  return str ;
}

// record both the forward and backward
void
create_plans_DFT( fftw_plan *__restrict forward , 
		  fftw_plan *__restrict backward ,
		  GLU_complex *__restrict *__restrict in , 
		  GLU_complex *__restrict *__restrict out ,
		  const size_t dims[ ND ] ,
		  const int ARR_SIZE ,
		  const int DIR )
{
  // set up our fft
  int dimes[ DIR ] , mu , planflag ;
  
  // swap these defs around as FFTW and I disagree on orderings
  for( mu = 0 ; mu < DIR ; mu++ ) {
    dimes[ mu ] = dims[ DIR - 1 - mu ] ;
  }

  #ifdef verbose
  start_timer( ) ;
  #endif
  
  char *str = obtain_wisdom( &planflag , dims , DIR , "" ) ;

  for( mu = 0 ; mu < ARR_SIZE ; mu++ ) {
    forward[mu] = fftw_plan_dft( DIR , dimes , in[mu] , out[mu] , 
				 FFTW_FORWARD , GLU_PLAN ) ; 
    backward[mu] = fftw_plan_dft( DIR , dimes , out[mu] , in[mu] , 
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

// plans for the DHT, only need one plan as it is its own inverse
void
create_plans_DHT( fftw_plan *__restrict plan , 
		  GLU_real *__restrict *__restrict in , 
		  GLU_real *__restrict *__restrict out ,
		  const size_t dims[ ND ] ,
		  const int ARR_SIZE ,
		  const int DIR )
{
  // set up our fft
  int dimes[ DIR ] , mu , planflag ;

  // swap these defs around 
  for( mu = 0 ; mu < DIR ; mu++ ) {
    dimes[ mu ] = dims[ DIR - 1 - mu ] ;
  }

  fftw_r2r_kind l[ND] ; 

  // we use a real to real dht, is is very much a halfcomplex to real trans .
  for( mu = 0 ; mu < ND ; mu++ ) {
    l[mu] = FFTW_DHT ;
  }

  // initialise the timer
  #ifdef verbose
  start_timer( ) ;
  #endif
  
  char *str = obtain_wisdom( &planflag , dims , DIR , "DHT_" ) ;

  // and organise the plans
  for( mu = 0 ; mu < ARR_SIZE ; mu++ ) {
    plan[mu] = fftw_plan_r2r( DIR , dimes , in[mu] , out[mu] , l , 
			      GLU_PLAN ) ; 
  }

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

/// Small plan does not care about the NC unlike the above
void
small_create_plans_DFT( fftw_plan *__restrict forward , 
			fftw_plan *__restrict backward ,
			GLU_complex *__restrict in , 
			GLU_complex *__restrict out ,
			const size_t dims[ ND ] ,
			const int DIR )
{
  // set up our fft
  int dimes[ DIR ] , mu , planflag ;
  // swap these defs around
  for( mu = 0 ; mu < DIR ; mu++ ) {
    dimes[ mu ] = dims[ DIR - 1 - mu ] ;
  }

  // initialise the clock
  #ifdef verbose
  start_timer( ) ;
  #endif

  char *str = obtain_wisdom( &planflag , dims , DIR , "single_" ) ;

  *forward = fftw_plan_dft( DIR , dimes , in , out , 
			    FFTW_FORWARD , GLU_PLAN ) ; 
  *backward = fftw_plan_dft( DIR , dimes , out , in , 
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

// and clean this up
#ifdef GLU_PLAN
  #undef GLU_PLAN
#endif

#endif
