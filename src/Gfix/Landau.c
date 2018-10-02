/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (Landau.c) is part of GLU.

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
   @file Landau.c
   @brief Landau gauge fixing routines

  Landau gauge fixing code, malloc's the ffts and
  calculates psq. Then calls the FA code in gtrans.c
  If we do not reach adequate convergence #GF_GLU_FAILURES restart(s)
  with a random gauge transform is(are) performed.
 */
#include "Mainfile.h"

#include "FACG.h"          // Fourier-Accelerated Conjugate Gradient
#include "geometry.h"      // lattice geometry, used for psq
#include "gftests.h"       // derivative evaluations
#include "GLU_malloc.h"    // malllocs
#include "MAG.h"           // for randomly restarting the MAG
#include "plan_ffts.h"     // FFTW wrappers
#include "plaqs_links.h"   // average plaquette and link trace
#include "read_headers.h"  // understands header formats
#include "read_config.h"   // configuration reader
#include "random_config.h" // random config and lattice reunit

// output the data, pass lat for the plaquette
static void
output_fixing_info( struct site *lat ,
		    const double theta ,
		    const size_t iters )
{
  // reunitarise just to limit the damage from round-off
  latt_reunitU( lat ) ;

  ////////// Print out the Gauge Fixing information /////////////

  printf( "[GF] Plaquette :: %1.15f \n[GF] Accuracy :: %1.4e\n" , 
	  av_plaquette( lat ) , theta ) ;
  GLU_real tr ;
  const double link = indivlinks( lat , &tr ) ;
  printf( "[GF] Iters :: %zu\n[GF] Link trace :: %1.15f || Maximum :: %1.15f\n" ,
	  iters , link , tr / NC ) ; 
  double lin , log ;
  const_time( lat , &lin , &log ) ; 
  printf( "[GF] Temporal constance || Lin %e || Log %e \n" , lin , log ) ;
  gauge_functional( lat ) ;
  printf( "[GF] Functional :: %1.15f\n" , gauge_functional( lat ) ) ;

  ///////////////////////////////////////////////////////////////
  return ;
}

// printing helper ...
static void
print_failure_info( const size_t failure , 
		    const size_t iters ,
		    const double theta ) 
{
  if( failure < GF_GLU_FAILURES ) {
    fprintf( stderr , "\n[GF] Non-convergence ... Randomly-restarting \n"
	     "\n[GF] Failure :: %zu || Iters :: %zu || ACC :: %e \n" 
	     , failure , iters , theta ) ;
  } else {
    fprintf( stderr , "\n[GF] Insufficient Convergence ......\n\n"
	    "[GF] Failures :: %zu || Total iterations :: %zu \n\n"
	    "[GF] -> Try reducing the tuning parameter OR/AND\n"
	    "increasing the maximum number of iterations -< \n\n" , 
	    failure , iters ) ; 
  }
  return ;
}

#if ( defined LUXURY_GAUGE )

// fast routine used here
static size_t
luxury_copy_fast( struct site *lat ,
		  struct fftw_stuff *FFTW ,
		  double *tr ,
		  const double acc ,
		  const size_t max_iters ) 
{
  size_t i , iters = 0 , copies ;

  struct site *lat_copy = NULL , *lat_best = NULL ;
  GLU_complex **gauge = NULL ;

  // allocate new gauge field
  if( GLU_malloc( (void**)&gauge , 16 , LVOLUME * sizeof( GLU_complex* ) ) != 0 ) {
    fprintf( stderr , "[GF] luxury copy fast failed to allocate "
	     "temporary gauge\n" ) ;
    goto memfree ;
  }
#pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    gauge[i] = ( GLU_complex* )malloc( NCNC * sizeof( GLU_complex ) ) ; 
    identity( gauge[i] ) ;
  }
  
  if( ( lat_copy = allocate_lat( ) ) == NULL ||
      ( lat_best = allocate_lat( ) ) == NULL ) {
    fprintf( stderr , "[GF] luxury gauge temporary allocation failure\n" ) ;
    goto memfree ;
  }

#ifdef BEST_COPY
  double maxlink = 1.0 , newlink ;
#else
  double maxlink = 0.0 , newlink ;
#endif
  // loop over the number of gauge copies !
  for( copies = 0 ; copies < LUXURY_GAUGE ; copies++ ) {
    
    // copy our lattice fields
#pragma omp parallel for private(i)
    for( i = 0 ; i < LVOLUME ; i++ ) {
      size_t mu ;
      for( mu = 0 ; mu < ND ; mu++ ) {
	equiv( lat_copy[i].O[mu] , lat[i].O[mu] ) ;
      }
    }

    // perform a random gauge transform
    random_transform( lat_copy , gauge ) ;
    // set the accuracy to be low
    const double tempacc = 1E-6 ;
    const double max = 1000 ;
    #ifdef GLU_GFIX_SD
    iters = FASD( lat_copy , FFTW , tr , tempacc , max ) ; 
    #else
    iters = FACG( lat_copy , FFTW , tr , tempacc , max ) ; 
    #endif
    
    // compute the link , wrap this to the functional?
    newlink = gauge_functional( lat_copy ) ;
    fprintf( stdout , "  [COPY] %zu [FUNCTIONAL] %1.15f [ITER] %zu " , 
	     copies , newlink , iters ) ; 
    #ifdef BEST_COPY 
    // the best copy is defined as the effective minimisation of the Gauge-functional 
    if( newlink < maxlink && iters != max && iters != 123456789 ) 
    #else
    if( newlink > maxlink && iters != max && iters != 123456789 ) 
    #endif
      {
	maxlink = newlink ;
	fprintf( stdout , " -> Copy accepted\n" ) ;
	// copy our lattice fields
        #pragma omp parallel for private(i)
	for( i = 0 ; i < LVOLUME ; i++ ) {
	  size_t mu ;
	  for( mu = 0 ; mu < ND ; mu++ ) {
	    equiv( lat_best[i].O[mu] , lat_copy[i].O[mu] ) ;
	  }
	}
      } else {
      double diff = newlink - maxlink ;
      fprintf( stdout , " -> Copy rejected %e\n" , diff ) ;
    }
  }

 
  // set our lattice to our chosen copy
#pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      equiv( lat[i].O[mu] , lat_best[i].O[mu] ) ;
    }
  }
  
  // final convergence run god I hope this one doesn't fail! Pretty unlikely
  #ifdef GLU_GFIX_SD
  iters += FASD( lat , FFTW , tr , acc , max_iters ) ; 
  #else
  iters += FACG( lat , FFTW , tr , acc , max_iters ) ; 
  #endif
  
  printf( "FINISHED\n" ) ;

 memfree :

  // free the gauge transformation matrices
  if( gauge != NULL ) {
#pragma omp parallel for private(i)
    for( i = 0 ; i < LVOLUME ; i++ ) {
      free( gauge[i] ) ;   
    }
    free( gauge ) ;
  }

  // free the lattice temporaries
  free_lat( lat_copy ) ;
  free_lat( lat_best ) ;
  
  return iters ;
}

#endif // luxury gauge

// cute little callback
static size_t
( *FA_callback ) ( struct site *lat ,
		   struct fftw_stuff *FFTW ,
		   double *tr ,
		   const double acc ,
		   const size_t max_iters ) ;

// callback selector
static void
select_callback( void ) 
{
#ifdef LUXURY_GAUGE
  FA_callback = luxury_copy_fast ;
#else
  #ifdef GLU_GFIX_SD
  FA_callback = FASD ;
  #else
  FA_callback = FACG ;
  #endif
#endif
  return ;
}

// should return failure if this really messes up
int
grab_file( struct site *lat , 
	   GLU_complex **gauge , 
	   const char *infile )
{
  struct head_data HEAD_DATA ;
  FILE *config = fopen( infile , "rb" ) ;
  // assume the chksum is correct
  if( read_header( config , &HEAD_DATA , GLU_FALSE ) == GLU_FAILURE ||
      get_config_SUNC( config , lat , HEAD_DATA ) == GLU_FAILURE ) {
    printf( "[IO] binary file read error \n" ) ;
    return GLU_FAILURE ;
  } 
  random_transform( lat , gauge ) ; 
  fclose( config ) ;
  // could reperform checksum calculations here ....
  return GLU_SUCCESS ;
}

// Landau gauge fixing routine uses callbacks
size_t 
Landau( struct site *lat ,
	const double accuracy ,
	const size_t max_iter ,
	const char *infile )
{
  double theta = 0. ; 
  size_t i ;
  
  struct fftw_stuff FFTW ;
  GLU_complex **gauge = NULL ;
  
#ifdef HAVE_FFTW3_H

  FFTW.psq = malloc( LVOLUME * sizeof( GLU_real ) ) ; 
  #pragma omp parallel
  {
    #pragma omp for private(i)
    for(  i = 0 ; i < LVOLUME ; i++  ) {
      FFTW.psq[i] = MAX_LANDAU / ( gen_p_sq( i , ND )  ) ; 
    }
  }

  /////////////// Look for Wisdom //////////////
  create_plans_DFT( &FFTW , Latt.dims , TRUE_HERM , ND ) ;

  ///////// End of the search for Wisdom //////
#else 
  FFTW.in = malloc( ( TRUE_HERM ) * sizeof( GLU_complex* ) ) ; 
  #pragma omp parallel for private(i)
  for(  i = 0 ; i < TRUE_HERM ; i++  ) {
    GLU_malloc( (void**)&FFTW.in[i] , ALIGNMENT , LVOLUME * sizeof( GLU_complex ) ) ;
  }
  // these are really dummy variables that don't get used in the SD
  FFTW.out = NULL ; FFTW.psq = NULL ;
  FFTW.forward = NULL , FFTW.backward = NULL ;
#endif

  // set up the FA method callback
  select_callback( )  ;

  size_t iters = FA_callback( lat , &FFTW , &theta , accuracy , max_iter ) ;

  // random restart portion of the code
  size_t failure = 0 ; 
  if( ( iters == 123456789 ) ) {

    // allocate new gauge field
    if( GLU_malloc( (void**)&gauge , 16 , LVOLUME * sizeof( GLU_complex* ) ) != 0 ) {
      fprintf( stderr , "[GF] GF_wrap_landau failed to allocate "
	       "temporary gauge\n" ) ;
      return GLU_FAILURE ;
    }
#pragma omp parallel for private(i)
    for( i = 0 ; i < LVOLUME ; i++ ) {
      gauge[i] = ( GLU_complex* )malloc( NCNC * sizeof( GLU_complex ) ) ; 
      identity( gauge[i] ) ;
    }

    print_failure_info( failure , max_iter , theta ) ;
    
    size_t iters_loc = 0 ;
    for( failure = 1 ; failure < GF_GLU_FAILURES ; failure++ ) {
      //repeat procedure have to read in the file :( 
      if( grab_file( lat , gauge , infile ) == GLU_FAILURE ) { 
	printf( "[IO] Something funky happened when trying to read in config again!\n" ) ;
	goto MemFree ; 
      }
      // and the callback
      iters_loc = FA_callback( lat , &FFTW , &theta , accuracy , max_iter ) ;

      // if we succeed we break the loop
      if( iters_loc != 123456789 ) {
	iters += iters_loc ; 
	break ;
      } else {// print information about the last failure
	iters += max_iter ;	
	print_failure_info( failure , iters_loc , theta ) ;
      }
    }
  }
  // end of random transform loop -> Usually not needed //

 MemFree :

  // free the gauge transformation matrices
  if( gauge != NULL ) {
#pragma omp parallel for private(i)
    for( i = 0 ; i < LVOLUME ; i++ ) {
      free( gauge[i] ) ;   
    }
    free( gauge ) ;
  }
  
#ifdef HAVE_FFTW3_H
  // free mallocs
  clean_up_fftw( FFTW , TRUE_HERM ) ;
#else
  #pragma omp parallel for private(i)
  for( i = 0 ; i < TRUE_HERM ; i++ ) {
    free( FFTW.in[i] ) ; 
  }
  free( FFTW.in ) ; 
#endif

  if( failure == GF_GLU_FAILURES ) {
    printf( "\n[GF] Failure to converge to a sufficient solution \n" 
	    "[GF] Please lower the GF_TUNE parameter from %g \n" , Latt.gf_alpha ) ;
    return GLU_FAILURE ;
  } else {
    output_fixing_info( lat , theta , iters ) ;
  }

  return iters ; 
}
