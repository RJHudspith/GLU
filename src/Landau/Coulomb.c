/*
    Copyright 2013 Renwick James Hudspith

    This file (Coulomb.c) is part of GLU.

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
  @file Coulomb.c
  @brief coulomb gauge fixing code, fixes each time-slice in turn and restarts with a random gauge transform if there is a problem ...
 */

#include "Mainfile.h"

#include "CFACG.h"        // Coulomb Conjugate Gradient
#include "geometry.h"     // for the p^2 calc
#include "plan_ffts.h"    // fftw plan wrapper
#include "plaqs_links.h"  // plaqutte and link measurements

// cute little callback
static int 
( *FA_callback ) ( struct site *__restrict lat ,
		   GLU_complex *__restrict *__restrict out , 
		   GLU_complex *__restrict *__restrict in , 
		   const void *__restrict forward , 
		   const void *__restrict backward , 
		   const GLU_real *psq , 
		   const double acc , 
		   const int max_iters ) ;

// callback selector
static void
select_callback( void ) 
{
#ifdef GLU_GFIX_SD
  FA_callback = Coulomb_FASD ;
#else
  FA_callback = Coulomb_FACG ; query_probes_Coulomb( ) ;
#endif
  return ;
}

// Coulomb gauge fixing code
int 
Coulomb( struct site *__restrict lat , 
	 const double accuracy , 
	 const int iter )
{
  double splink , tlink ;
  all_links( lat , &splink , &tlink ) ;

  printf( "[GF] Initial Tlink :: %1.15f || Slink :: %1.15f \n"
	  "[GF] Plaquette :: %1.15f \n", 
	  tlink , splink , av_plaquette( lat ) ) ; 

  // put an ifdef guard here as SD requires none of this ...
#ifdef HAVE_FFTW3_H

  if( parallel_ffts( ) == GLU_FAILURE ) {
    printf( "Problem with initialising the OpenMP FFTW routines \n" ) ;
    // should clean up the memory here
    return GLU_FAILURE ;
  }

  fftw_plan *forward = malloc( ( TRUE_HERM ) * sizeof( fftw_plan ) ) ; 
  fftw_plan *backward = malloc( ( TRUE_HERM ) * sizeof( fftw_plan ) ) ; 
  GLU_complex **out = fftw_malloc( ( TRUE_HERM ) * sizeof( GLU_complex* ) ) ; 
  GLU_complex **in = fftw_malloc( ( TRUE_HERM ) * sizeof( GLU_complex* ) ) ; 

  int i ; 
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < TRUE_HERM ; i ++  ) {
    out[i] = ( GLU_complex* )fftw_malloc( LCU * sizeof( GLU_complex ) ) ; 
    in[i] = ( GLU_complex* )fftw_malloc( LCU * sizeof( GLU_complex ) ) ; 
  }

  GLU_real *psq = malloc( LCU * sizeof( GLU_real ) ) ; 

  // we calculate the lattice p-squared here and pass it to the FFT-accelerator
#pragma omp parallel for private(i) 
  PFOR( i = 0 ; i < LCU ; i ++ ) {
    psq[i] = MAX_COULOMB / ( gen_p_sq( i , ND - 1 ) ) ;
  }

  // create the fftw plans, or read them if they are stored
  create_plans_DFT( forward , backward , in , out , TRUE_HERM , ND - 1 ) ;

  /////  End of the search for Wisdom  ////
#else

 GLU_complex **in = malloc( ( TRUE_HERM ) * sizeof( GLU_complex* ) ) ; 
  int i ;
  #pragma omp parallel for private(i)
  PFOR(  i = 0 ; i < TRUE_HERM ; i++  ) {
    in[i] = ( GLU_complex* )malloc( LCU * sizeof( GLU_complex ) ) ; 
  }
  // these are dummy arrays
  GLU_complex **out = NULL ;
  GLU_real *psq = NULL ;
  int *forward = NULL , *backward = NULL ;

#endif
  select_callback( ) ;

  const int iters =  FA_callback( lat , out , in ,
				  forward , backward , 
				  psq , 
				  accuracy , iter ) ;
 
#ifdef HAVE_FFTW3_H

  // cleanup the mallocs and whatnot
#pragma omp parallel for private(i) 
  PFOR( i = 0 ; i < TRUE_HERM ; i ++  ) {
      fftw_destroy_plan( forward[i] ) ;  
      fftw_destroy_plan( backward[i] ) ;  
      fftw_free( out[i] ) ; 
      fftw_free( in[i] ) ; 
    }
  free( forward ) ;
  free( backward );
  fftw_cleanup(  ) ; 
  free( psq ) ; 
  free( out ) ; 
  free( in ) ; 

#else
  
  #pragma omp parallel for private(i) 
  PFOR( i = 0 ; i < TRUE_HERM ; i ++  ) {
    free( in[i] ) ;
  }
  free( in ) ;

#endif

  all_links( lat , &splink , &tlink ) ;
  printf( "[GF] Tuning :: %f || Iterations :: %d ||\n[GF] Final Tlink :: %1.15f "
	  "|| Slink :: %1.15f \n[GF] Plaquette :: %1.15f \n" , 
	  Latt.gf_alpha , iters , tlink , splink , av_plaquette( lat ) ) ; 

  return iters ; 
}

