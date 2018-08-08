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
#include "MAG.h"           // for randomly restarting the MAG
#include "plan_ffts.h"     // FFTW wrappers
#include "plaqs_links.h"   // average plaquette and link trace
#include "read_headers.h"  // understands header formats
#include "read_config.h"   // configuration reader
#include "random_config.h" // random config and lattice reunit

// output the data, pass lat for the plaquette
static void
output_fixing_info( struct site *__restrict lat ,
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
		    const size_t iters ) 
{
  if( failure < GF_GLU_FAILURES ) {
    fprintf( stderr , "\n[GF] Non-convergence ... Randomly-restarting \n"
	     "\n[GF] Failure :: %zu || Iters :: %zu \n" 
	     , failure , iters ) ; 
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
luxury_copy_fast( struct site *__restrict lat ,
		  GLU_complex *__restrict *__restrict gauge ,
		  GLU_complex *__restrict *__restrict out ,
		  GLU_complex *__restrict *__restrict in ,
		  const void *forward ,
		  const void *backward ,
		  const GLU_real *__restrict psq ,
		  double *tr ,
		  const double acc ,
		  const size_t max_iters ) 
{
  size_t i , iters = 0 , copies ;

  struct site *lat_copy = NULL , *lat_best = NULL ;
  if( GLU_malloc( (void**)&lat_copy , ALIGNMENT , LVOLUME * sizeof( struct site ) ) != 0 ||
      GLU_malloc( (void**)&lat_best , ALIGNMENT , LVOLUME * sizeof( struct site ) ) != 0 ) {
    fprintf( stderr , "[GF] luxury_copy_fast temporary lattice allocation failure\n" ) ;
    return iters ;
  }

  init_navig( lat_copy ) ;
  init_navig( lat_best ) ;

#ifdef BEST_COPY
  double maxlink = 1.0 , newlink ;
#else
  double maxlink = 0.0 , newlink ;
#endif
  // loop over the number of gauge copies !
  for( copies = 0 ; copies < LUXURY_GAUGE ; copies++ ) {
    // copy our lattice fields
    #pragma omp parallel for private(i)
    PFOR(  i = 0 ; i < LVOLUME ; i++  ) {
      memcpy( &lat_copy[i] , &lat[i] , sizeof( struct site ) ) ;
    }
    // perform a random gauge transform
    random_transform( lat_copy , gauge ) ;
    // set the accuracy to be low
    const double tempacc = 1E-6 ;
    const double max = 1000 ;
    #ifdef GLU_GFIX_SD
    iters = FASD( lat_copy ,  
		  out , in ,
		  forward , backward , 
		  psq , tr , tempacc , max ) ; 
    #else
    iters = FACG( lat_copy ,
		  out , in ,
		  forward , backward , 
		  psq , tr , tempacc , max ) ; 
    #endif
    
    // compute the link , wrap this to the functional?
    newlink = gauge_functional( lat_copy ) ;
    fprintf( stdout , "  [COPY] %zu [FUNCTIONAL] %1.15f [ITER] %zu " , 
	     copies , newlink , iters ) ; 
    #ifdef BEST_COPY 
    // the best copy is defined as the effective minimisation of the Gauge-functional 
    if( newlink < maxlink && iters != max ) 
    #else
    if( newlink > maxlink && iters != max ) 
    #endif
      {
	fprintf( stdout , " -> Copy accepted \n" ) ;
	maxlink = newlink ;
        #pragma omp parallel for private(i)
	PFOR(  i = 0 ; i < LVOLUME ; i++  ) {
	  memcpy( &lat_best[i] , &lat_copy[i] , sizeof( struct site ) ) ;
	}
      } else { fprintf( stdout , " -> Copy rejected %e\n" , *tr ) ; }
  }
  // set our lattice to our chosen copy
  #pragma omp parallel for private(i)
  PFOR(  i = 0 ; i < LVOLUME ; i++  ) {
    memcpy( &lat[i] , &lat_best[i] , sizeof( struct site ) ) ;
  }
  // final convergence run god I hope this one doesn't fail! Pretty unlikely
  #ifdef GLU_GFIX_SD
  iters += FASD( lat , 
		 out , in , 
		 forward , backward , 
		 psq , 
		 tr , acc , max_iters ) ; 
  #else
  iters += FACG( lat , 
		 out , in , 
		 forward , backward , 
		 psq , 
		 tr , acc , max_iters ) ; 
  #endif

  // free the lattice temporaries
  free( lat_copy ) ;
  free( lat_best ) ;

  return iters ;
}

#endif // luxury gauge

// cute little callback
static size_t
( *FA_callback ) ( struct site *__restrict lat ,
		   struct fftw_stuff FFTW ,
		   double *tr ,
		   const double acc ,
		   const size_t max_iters ) ;

// callback selector
static void
select_callback( const int improvement ) 
{
  #ifdef LUXURY_GAUGE
  if( improvement != SMPREC_IMPROVE ) {
    FA_callback = luxury_copy_fast ; 
  } else { 
    #ifdef GLU_GFIX_SD
    FA_callback = FASD ; 
    #else
    FA_callback = FACG ;
    #endif
  }
  #else
  if( improvement == SMPREC_IMPROVE ) {
    // fix this!!!
    //FA_callback = FASD_SMEAR ;
    FA_callback = FASD ;
  } else {
    #ifdef GLU_GFIX_SD
    FA_callback = FASD ;
    #else
    FA_callback = FACG ;
    #endif
  } 
  #endif
  return ;
}

// should return failure if this really messes up
int
grab_file( struct site *__restrict lat , 
	   GLU_complex *__restrict *__restrict gauge , 
	   const char *__restrict infile )
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
Landau( struct site *__restrict lat ,
	GLU_complex *__restrict *__restrict gauge ,
	const double accuracy ,
	const size_t iter ,
	const char *__restrict infile ,
	const GF_improvements improvement )
{
  double theta = 0. ; 
  size_t i ;
  
  struct fftw_stuff FFTW ;

#ifdef HAVE_FFTW3_H

  FFTW.psq      = malloc( LVOLUME * sizeof( GLU_real ) ) ; 
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
  PFOR(  i = 0 ; i < TRUE_HERM ; i++  ) {
    GLU_malloc( (void**)&FFTW.in[i] , ALIGNMENT , LVOLUME * sizeof( GLU_complex ) ) ;
  }
  // these are really dummy variables that don't get used in the SD
  FFTW.out = NULL ; FFTW.psq = NULL ;
  FFTW.forward = NULL , FFTW.backward = NULL ;
#endif

  // set up the FA method callback
  select_callback( improvement )  ;

  size_t iters = FA_callback( lat , FFTW , &theta , accuracy , iter ) ;

  // random restart portion of the code
  size_t failure = 0 ; 
  if( ( iters == 123456789 ) && ( improvement != SMPREC_IMPROVE ) ) {

    print_failure_info( failure , iter ) ;
    
    size_t iters_loc = 0 ;
    for( failure = 1 ; failure < GF_GLU_FAILURES ; failure++ ) {
      //repeat procedure have to read in the file :( 
      if( grab_file( lat , gauge , infile ) == GLU_FAILURE ) { 
	printf( "[IO] Something funky happened when trying to read in config again!\n" ) ;
	goto MemFree ; 
      }
      if( improvement == MAG_IMPROVE ) { mag( lat , gauge ) ; }
      
      // and the callback
      iters_loc = FA_callback( lat , FFTW ,
			       &theta , accuracy , iter ) ;

      // if we succeed we break the loop
      if( iters_loc != 123456789 ) {
	iters += iters_loc ; 
	break ;
      } else {// print information about the last failure
	iters += iter ;	
	print_failure_info( failure , iters_loc ) ;
	printf( "\n[GF] Failure :: %zu || Accuracy %1.4e\n" , 
		failure , theta ) ;
      }
    }
  }
  // end of random transform loop -> Usually not needed //

 MemFree :
  
#ifdef HAVE_FFTW3_H
  // free mallocs
  clean_up_fftw( FFTW , TRUE_HERM ) ;
#else
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < TRUE_HERM ; i++ ) {
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
