/**
   @file KPHB.c
   @brief Kennedy-Pendleton heat bath algorithm

   TODO :: test this works/ compiles OK
           draughtboard somewhere else, Geometry? Own file?
 */
#include "Mainfile.h"

#include <omp.h>

#include "draughtboard.h"  // draughtboarding
#include "GLU_timer.h"     // print_time()
#include "hb.h"            // heat-bath
#include "par_rng.h"       // parallel KISS rng
#include "plaqs_links.h"   // av_plaquette()
#include "random_config.h" // reunit_latt()
#include "relax.h"         // overrelaxation
#include "GLUlib_wrap.h"   // write out a configuration

static void
update_lattice( struct site *lat ,
		struct site *staple ,
		const double inverse_beta ,
		const struct draughtboard db ,
		const size_t Nor )
{
  // do a heat bath
  hb_lattice( lat , staple , inverse_beta , db ) ;

  // and some number of over-relaxations
  size_t or ;
  for( or = 0 ; or < Nor ; or++ ) {
    OR_lattice( lat , staple , db ) ;
  }

  // reunitarise the gauge field? If NC gets large this can be a problem
  latt_reunitU( lat ) ;

  return ;
}

// perform a heatbath
int
hb_update( struct site *lat ,
	   const struct hb_info HBINFO ,
	   const char *traj_name ,
	   const GLU_output storage , 
	   const char *output_details ,
	   const GLU_bool continuation )
{
  // counters
  size_t i ;

  // seed the parallel RNG
  char str[ 256 ] ;

  // strip the number in the infile if it has one
  char *pch = strtok( (char*)traj_name , "." ) ;

  sprintf( str , "%s.%zu.rand" , pch , Latt.flow ) ;
  if( continuation == GLU_FALSE ) {
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
  init_cb( &db , LVOLUME , ND ) ;

  // temporary space for staples
  struct site *staple = NULL ; 
  GLU_malloc( (void**)&staple , ALIGNMENT , 
	      (db.Nred>db.Nblack?db.Nred:db.Nblack)*sizeof( struct site ) ) ;

  // this guy appears throughout the HB algorithm
  const double inverse_beta = NC/(HBINFO.beta) ;

  // give us some information
  fprintf( stdout , "[UPDATE] Performing %zu HB-OR iterations\n" , 
	   HBINFO.iterations ) ;
  fprintf( stdout , "[UPDATE] Themalising for %zu iterations\n" , 
	   HBINFO.therm ) ;
  fprintf( stdout , "[UPDATE] %zu over-relaxations per heatbath\n" ,
	   HBINFO.Nor ) ;
  fprintf( stdout , "[UPDATE] Saving every %zu iteration(s)\n" ,
	   HBINFO.Nsave ) ;

  // thermalise
  start_timer( ) ;
  for( i = 0 ; i < HBINFO.therm ; i++ ) {
    update_lattice( lat , staple , inverse_beta , db , HBINFO.Nor ) ;
    if( !(i&15) ) {
      fprintf( stdout , "\n[UPDATE] config %zu done \n" , i ) ;
      print_time() ;
    }
  }
  print_time( ) ;

  // iterate the number of runs
  const size_t start = Latt.flow ;
  for( i = Latt.flow ; i < HBINFO.iterations ; i++ ) {
    update_lattice( lat , staple , inverse_beta , db , HBINFO.Nor ) ;
    // set the lattice flow
    // if we are saving the data print out the plaquette and write a file
    if( i%HBINFO.Nmeasure == 0 ) {
      fprintf( stdout , "[UPDATE] %zu :: {P} %1.12f \n" , 
	       i , av_plaquette( lat ) ) ;
    }
    // if we hit a save point we write out the configuration
    if( i%HBINFO.Nsave == 0 && i != start ) {
      // write a configuration
      sprintf( str , "%s.%zu" , traj_name , i ) ;
      write_configuration( lat , str , storage , output_details ) ;
      // write out the rng state
      sprintf( str , "%s.rand" , str ) ;
      write_par_rng_state( str ) ;
    }
    Latt.flow = i + 1 ;
  }

  // free the draughtboard
  free_cb( &db ) ;

  // free the staples
  free( staple ) ;

  // free the rng
  free_par_rng( ) ;

  return GLU_SUCCESS ;
}

