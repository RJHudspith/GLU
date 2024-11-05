/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (KPHB.c) is part of GLU.

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
   @file KPHB.c
   @brief Heat bath over-relaxation algorithm
 */
#include "Mainfile.h"

#include "draughtboard.h"  // draughtboarding
#include "GLUlib_wrap.h"   // write out a configuration
#include "GLU_timer.h"     // print_time()
#include "gramschmidt.h"   // gram_reunit()
#include "hb.h"            // heat-bath
#include "par_rng.h"       // parallel rngs
#include "plaqs_links.h"   // av_plaquette()
#include "POLY.h"          // poly()
#include "random_config.h" // reunit_latt()
#include "relax.h"         // overrelaxation
#include "str_stuff.h"     // append_char

//#define POLY_MEAS

// measure plaquette and polyakov loop
static void
perform_measurements( double *PLAQred ,
		      double *POLYred ,
		      const struct site *lat ,
		      const size_t conf_idx )
{
  double pl = 0.0 ;
  size_t k ;
  
#pragma omp for private(k)
  for( k = 0 ; k < Latt.Nthreads*CLINE ; k++ ) {
    PLAQred[ k ] = 0.0 ;
  }
  
  av_plaquette_th( PLAQred , lat ) ;
  
  for( k = 0 ; k < Latt.Nthreads ; k++ ) {
    pl += PLAQred[ 3 + CLINE*k ] ;
  }
  
#pragma omp master
  {
    // write out the plaquette
    fprintf( stdout , "[UPDATE] %zu :: {P} %1.12f \n" , 
	     conf_idx , pl/( NC*(ND-1)*(ND-2)*LVOLUME ) ) ;
  }

  // icc really hates this for some reason
#if POLY_MEAS
  double re = 0 , im = 0 ;
  size_t mu ;
  
#pragma omp for private(k)
  for( k = 0 ; k < Latt.Nthreads*CLINE ; k++ ) {
    POLYred[ k ] = 0.0 ;
  }
  
  // write the polyakov loops, (re,im) |L|
  for( mu = 0 ; mu < ND ; mu++ ) {
    poly_th( POLYred , lat , mu ) ;
  }
  
#pragma omp master
  {
    for( mu = 0 ; mu < ND ; mu++ ) {
      re = im = 0.0 ;
      for( k = 0 ; k < Latt.Nthreads ; k++ ) {
	re += POLYred[ 2*mu + k*CLINE ] ;
	im += POLYred[ 2*mu + 1 + k*CLINE ] ;
      }
      fprintf( stdout , "[UPDATE] {L_%zu} ( %1.12e , %1.12e ) %1.12e \n" ,
	       mu , re , im , sqrt( re*re + im*im ) ) ;
    }
    fflush( stdout ) ;
  }
#endif
  return ;
}

// updates the lattice
static void
update_lattice( struct site *lat ,
		const double inverse_beta ,
		const struct draughtboard db ,
		const size_t Nor )
{
  // do a Nhb heat baths
  hb_lattice( lat , inverse_beta , db ) ;
  
  // and some number of over-relaxations
  // this keeps the action the same while increasing tunneling into 
  // different topological sectors
  size_t or ;
  for( or = 0 ; or < Nor ; or++ ) {
    // perform over-relaxation 
    OR_lattice( lat , db ) ;
  }
  return ;
}

// perform a heatbath
int
hb_update( struct site *lat ,
	   const struct hb_info HBINFO ,
	   const char *traj_name ,
	   const GLU_output storage , 
	   const char *output_details )
{
  // strip the number in the infile if it has one
  char *pch = strtok( (char*)traj_name , "." ) ;
  char str[ strlen(pch)+6+sizeof(size_t) ] ;
  sprintf( str , "%s.%zu.rand" , pch , Latt.flow ) ;
  
  if( HBINFO.continuation == GLU_FALSE ) {
    fprintf( stdout , "[UPDATE] initialising par_rng from pool\n" ) ;
    if( initialise_par_rng( NULL ) == GLU_FAILURE ) {
      return GLU_FAILURE ;
    }
  } else {
    fprintf( stdout , "[UPDATE] restarting from a previous trajectory\n" ) ;
    if( initialise_par_rng( str ) == GLU_FAILURE ) {
      return GLU_FAILURE ;
    }
  }
  
  // initialise the draughtboard
  struct draughtboard db ;
#ifdef IMPROVED_SMEARING
  if( init_improved_cb( &db ) == GLU_FAILURE ) {
    return GLU_FAILURE ;
  }
#else
  if( init_cb( &db , LVOLUME , ND ) == GLU_FAILURE ) {
    return GLU_FAILURE ;
  }
#endif

  // this guy appears throughout the HB algorithm
  const double inverse_beta = NC/(HBINFO.beta) ;

  // give us some information
  fprintf( stdout , "[UPDATE] Performing %zu HB-OR iterations\n" , 
	   HBINFO.iterations ) ;
  fprintf( stdout , "[UPDATE] NSTOCH updates %d\n" , NSTOCH ) ;
  fprintf( stdout , "[UPDATE] Themalising for %zu iterations\n" , 
	   HBINFO.therm ) ;
  fprintf( stdout , "[UPDATE] %zu over-relaxations per heatbath\n" ,
	   HBINFO.Nor ) ;
  fprintf( stdout , "[UPDATE] Saving every %zu iteration(s)\n" ,
	   HBINFO.Nsave ) ;
  fprintf( stdout , "[UPDATE] Using beta = %1.12f \n" , HBINFO.beta ) ;

  double *PLAQred = calloc( CLINE*Latt.Nthreads , sizeof( double ) ) ;
  double *POLYred = calloc( CLINE*Latt.Nthreads , sizeof( double ) ) ;
  
  // thermalise
  start_timer( ) ;
  
#pragma omp parallel
  {   
    size_t i ;
    for( i = 0 ; i < HBINFO.therm ; i++ ) {
      update_lattice( lat , inverse_beta , db , HBINFO.Nor ) ;
      if( !(i&15) ) {
        #pragma omp master
	{
	  fprintf( stdout , "\n[UPDATE] config %zu done \n" , i ) ;
	  print_time() ;
	}
      }
    }

    #pragma omp master
    {
      if( HBINFO.therm != 0 ) {
	fprintf( stdout , "\n[UPDATE] Trajectory Thermalised \n" ) ;
	print_time( ) ;
      }
    }
    
    // iterate the number of runs
    const size_t start = Latt.flow ;
    for( i = Latt.flow ; i < HBINFO.iterations ; i++ ) {

      // if we are saving the data print out the plaquette and write a file
      if( i%HBINFO.Nmeasure == 0 ) {
	perform_measurements( PLAQred , POLYred , lat , i ) ;
      }

      // if we hit a save point we write out the configuration
      if( i%HBINFO.Nsave == 0 && i != start ) {
	#pragma omp single
	{
	  // write a configuration
	  sprintf( str , "%s.%zu" , pch , i ) ;
	  if( write_configuration( lat , str , storage , output_details )
	      != GLU_FAILURE ) {
	    // write out the rng state
	    sprintf( str , "%s.%zu.rand" , pch , i ) ;
	    write_par_rng_state( str ) ;
	  }
	}
      }

      // perform a hb-OR step
      update_lattice( lat , inverse_beta , db , HBINFO.Nor ) ;

      #pragma omp master
      {
	Latt.flow = i + 1 ; 
      }
    }
  }

  // free the reduction array
  free( PLAQred ) ;
  free( POLYred ) ;

  // free the draughtboard
  free_cb( &db ) ;

  // tell us how long this generation took
  print_time() ;

  // free the rng
  free_par_rng( ) ;

  return GLU_SUCCESS ;
}

