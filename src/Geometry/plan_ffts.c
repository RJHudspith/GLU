/*
    Copyright 2013 Renwick James Hudspith

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

// if we have these 
#if ( defined OMP_FFTW ) && ( defined HAVE_OMP_H )
 #include <omp.h>
 static int nthreads = 1 ;
#endif

// for ease of reading
enum{ NOPLAN = 0 } ;

// have a look for parallel ffts
short int
parallel_ffts()
{
  // initialise parallel fft ?
#ifdef OMP_FFTW
  if( fftw_init_threads() == 0 ) {
    return GLU_FAILURE ;
  } else {
    // in here I set the number of fftw threads to be the same
    // as the usual number of parallel threads ..
    #pragma omp parallel
    { nthreads = omp_get_num_threads( ) ; } // set nthreads
    fftw_plan_with_nthreads( nthreads ) ;
    printf("[PAR] FFTW using %d thread(s) \n" , nthreads ) ;
  }
#endif
  return GLU_SUCCESS ;
}

// record both the forward and backward
void
create_plans_DFT( fftw_plan *__restrict forward , 
		  fftw_plan *__restrict backward ,
		  GLU_complex *__restrict *__restrict in , 
		  GLU_complex *__restrict *__restrict out , 
		  const int ARR_SIZE ,
		  const int DIR )
{
  // set up our fft
  int dimes[ DIR ] , mu ;
  // swap these defs around
  for( mu = 0 ; mu < DIR ; mu++ ) {
    dimes[ mu ] = Latt.dims[ DIR - 1 - mu ] ;
  }

  start_timer( ) ;

#ifndef CONDOR_MODE
  char str[512] , prec_str[ 16 ] ; 
  int check = NOPLAN ; 

  #ifdef SINGLE_PREC
  sprintf( prec_str , "FLOAT" ) ;
  #else
  sprintf( prec_str , "DOUBLE" ) ;
  #endif

  #ifdef OMP_FFTW
  sprintf( str , "%s/Local/Wisdom/%s_OMPnt%d_SU%d_" , 
	   HAVE_PREFIX , prec_str , nthreads , NC ) ;
  #else
  sprintf( str , "%s/Local/Wisdom/%s_SU%d_" , 
	   HAVE_PREFIX , prec_str , NC ) ;
  #endif
  for( mu = 0 ; mu < DIR - 1 ; mu++ ) {
    sprintf( str , "%s%dx" , str , Latt.dims[ mu ] ) ;
  }
  sprintf( str , "%s%d.wisdom" , str , Latt.dims[ DIR - 1 ] ) ;
  FILE *wizzard = fopen( str , "r" ) ; 
  if( wizzard == NULL ) {
    printf( "\n[FFTW] No wisdom to be obtained here ... planning" ) ; 
  } else {
    printf( "\n[FFTW] Successful wisdom attained" ) ; 
    check = fftw_import_wisdom_from_file( wizzard ) ; 
    fclose( wizzard ) ; 
  }
#else
  printf( "\n[FFTW] Creating plan on CONDOR host" ) ; 
#endif

  for( mu = 0 ; mu < ARR_SIZE ; mu++ ) {
    forward[mu] = fftw_plan_dft( DIR , dimes , in[mu] , out[mu] , 
				 FFTW_FORWARD , GLU_PLAN ) ; 
    backward[mu] = fftw_plan_dft( DIR , dimes , out[mu] , in[mu] , 
				  FFTW_BACKWARD , GLU_PLAN ) ;
  }

  // I want to know how long FFTW is taking to plan its FFTs
  print_time( ) ;
  printf( "[FFTW] plans finished\n\n" ) ;

#ifndef CONDOR_MODE
  if( check == NOPLAN )  {
    wizzard = fopen( str , "w" ) ; 
    fftw_export_wisdom_to_file( wizzard ) ; 
    fclose( wizzard ) ; 
  }
#endif
  return ;
}

// plans for the DHT, only need one plan as it is its own inverse
void
create_plans_DHT( fftw_plan *__restrict plan , 
		  GLU_real *__restrict *__restrict in , 
		  GLU_real *__restrict *__restrict out , 
		  const int ARR_SIZE ,
		  const int DIR )
{
  // set up our fft
  int dimes[ DIR ] ;

  // swap these defs around 
  int mu ;
  for( mu = 0 ; mu < DIR ; mu++ ) {
    dimes[ mu ] = Latt.dims[ DIR - 1 - mu ] ;
  }

  fftw_r2r_kind l[ND] ; 

  // we use a real to real dht, is is very much a halfcomplex to real trans .
  for( mu = 0 ; mu < ND ; mu++ ) {
    l[mu] = FFTW_DHT ;
  }

  // initialise the timer
  start_timer( ) ;

#ifndef CONDOR_MODE
  char str[512] , prec_str[16] ; 
  int check = 0 ;

#ifdef SINGLE_PREC
  sprintf( prec_str , "FLOAT" ) ;
#else
  sprintf( prec_str , "DOUBLE" ) ;
#endif

#ifdef OMP_FFTW
  sprintf( str , "%s/Local/Wisdom/%s_DHT_OMPnt%d_" , 
	   HAVE_PREFIX , prec_str , nthreads ) ;
#else
  sprintf( str , "%s/Local/Wisdom/%s_DHT_" , 
	   HAVE_PREFIX , prec_str ) ;
#endif
  for( mu = 0 ; mu < DIR - 1 ; mu++ ){
    sprintf( str , "%s%dx" , str , Latt.dims[ mu ] ) ;
  }
  sprintf( str , "%s%d.wisdom", str , Latt.dims[ DIR - 1 ] ) ;
  FILE *wizzard = fopen( str , "r" ) ; 
  if( wizzard == NULL ) {
    printf( "\n[FFTW] No wisdom to be obtained here" ) ; 
  } else {
    printf( "\n[FFTW] Successful wisdom attained" ) ; 
    check = fftw_import_wisdom_from_file( wizzard ) ; 
    fclose( wizzard ) ; 
  }
#else
  printf( "[FFTW] Creating plan on CONDOR host" ) ; 
#endif

  // and organise the plans
  for( mu = 0 ; mu < ARR_SIZE ; mu++ ) {
    plan[mu] = fftw_plan_r2r( DIR , dimes , in[mu] , out[mu] , l , 
			      GLU_PLAN ) ; 
  }

  // I want to know how long FFTW is taking to plan its FFTs
  print_time( ) ;
  printf( "[FFTW] plans finished\n\n" ) ;

#ifndef CONDOR_MODE
  if( check == NOPLAN ) {
    wizzard = fopen( str , "w" ) ; 
    fftw_export_wisdom_to_file( wizzard ) ; 
    fclose( wizzard ) ; 
  }
#endif
  return ;
}

/// Small plan does not care about the NC unlike the above
void
small_create_plans_DFT( fftw_plan *__restrict forward , 
			fftw_plan *__restrict backward ,
			GLU_complex *__restrict in , 
			GLU_complex *__restrict out ,
			const int DIR )
{
  // set up our fft
  int dimes[ DIR ] ;
  // swap these defs around
  int mu ;
  for( mu = 0 ; mu < DIR ; mu++ ) {
    dimes[ mu ] = Latt.dims[ DIR - 1 - mu ] ;
  }

  // initialise the clock
  start_timer( ) ; 

#ifndef CONDOR_MODE
  char str[512] , prec_str[16] ; 
  int check = NOPLAN ; 

#ifdef SINGLE_PREC
  sprintf( prec_str , "FLOAT" ) ;
#else
  sprintf( prec_str , "DOUBLE" ) ;
#endif

#ifdef OMP_FFTW
  sprintf( str , "%s/Local/Wisdom/%s_OMPnt%d_single_" , 
	   HAVE_PREFIX , prec_str , nthreads ) ;
#else
  sprintf( str , "%s/Local/Wisdom/%s_single_" , 
	   HAVE_PREFIX , prec_str ) ;
#endif
  for( mu = 0 ; mu < DIR - 1 ; mu++ ) {
    sprintf( str , "%s%dx" , str , Latt.dims[ mu ] ) ;
  }
  sprintf( str , "%s%d.wisdom", str , Latt.dims[ DIR - 1 ] ) ;
  FILE *wizzard = fopen( str , "r" ) ; 
  if( wizzard == NULL ) {
    printf( "\n[FFTW] No wisdom to be obtained here" ) ; 
  } else {
    printf( "\n[FFTW] Successful wisdom attained" ) ; 
    check = fftw_import_wisdom_from_file( wizzard ) ; 
    fclose( wizzard ) ; 
  }
#else
  printf( "\n[FFTW] Creating plan on CONDOR host" ) ; 
#endif

  *forward = fftw_plan_dft( DIR , dimes , in , out , 
			    FFTW_FORWARD , GLU_PLAN ) ; 
  *backward = fftw_plan_dft( DIR , dimes , out , in , 
			     FFTW_BACKWARD , GLU_PLAN ) ; 

  // I want to know how long FFTW is taking to plan its FFTs
  print_time( ) ;
  printf( "[FFTW] plans finished\n\n" ) ;

#ifndef CONDOR_MODE
  if( check == NOPLAN ) {
    wizzard = fopen( str , "w" ) ; 
    fftw_export_wisdom_to_file( wizzard ) ; 
    fclose( wizzard ) ; 
  }
#endif
  return ;
}

// and clean this up
#ifdef GLU_PLAN
  #undef GLU_PLAN
#endif

#endif
