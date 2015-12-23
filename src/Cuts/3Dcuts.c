/*
    Copyright 2013 Renwick James Hudspith

    This file (3Dcuts.c) is part of GLU.

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
   @file 3Dcuts.c
   @brief This code computes the insantaneous spatial and instantaneous
   temporal gluon propagators

   it writes them to a file once the usual momentum cuts have been performed.

   Write format ( BIG ENDIAN )  <br>
   Number of Momenta :: nmom    <br>
   mom_list :: ND-1 (int)a (int)b .... <br> 
   Number of Momenta :: nmom <br> 
   Instantaneous spatial glue-prop :: \f$ g_2^{spatial} \f$ <br> 
   Number of Momenta :: nmom <br> 
   Instantaneous temporal glue-prop :: \f$ g_2^{temporal} \f$ <br> 
*/

#include "Mainfile.h"

#include "cut_output.h"    // outputting
#include "cut_routines.h"  // momentum cuts
#include "plan_ffts.h"     // FFTW plan setting
#include "pspace_landau.h" // momentum space Coulomb condition check

#ifdef HAVE_FFTW3_H

static int
mom_gauge_spatial( struct site *__restrict A ,
		   const lie_field_def def )
{
  size_t i ; 

  // callback for the log definition
  void (*log)( GLU_complex Q[ NCNC ] ,
	       const GLU_complex U[ NCNC ] ) ;
  switch( def ) {
  case LINEAR_DEF :
    log = Hermitian_proj ;
    break ;
  case LOG_DEF : 
    log = exact_log_slow ; 
    break ;
  }

#pragma omp parallel for private(i) shared(A)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    GLU_complex temp[ NCNC ] ;
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      log( temp , A[i].O[mu] ) ;
      equiv( A[i].O[mu] , temp ) ;
    }
  }

  // FFTW routines
  GLU_complex *out = fftw_malloc( LCU * sizeof( GLU_complex ) ) ; 
  GLU_complex *in = fftw_malloc( LCU * sizeof( GLU_complex ) ) ; 

  fftw_plan forward , backward ;
  small_create_plans_DFT( &forward , &backward , in , out , ND - 1 ) ;

  //// End of the search for Wisdom /////
  
  //3D complex fourier transform each time-slice
  size_t t ;
  for( t = 0 ; t < Latt.dims[ ND - 1 ] ; t++ )  {
    const size_t slice = LCU *t ;
    size_t mu , i , j ; 
    // Loop directions and elements of the matrix
    for( mu = 0 ; mu < ND ; mu++ )
      for( j = 0 ; j < NCNC ; j++ ) {
	// FORWARDS
#ifdef CUT_FORWARD
        #pragma omp parallel for private(i)
	PFOR( i = 0 ; i < LCU ; i++ ) {
	  in[i] = A[ slice + i ].O[mu][j] ; 
	}
	fftw_execute( forward ) ; 
        #pragma omp parallel for private(i)
	PFOR( i = 0 ; i < LCU ; i++ ) {
	  A[ slice + i ].O[mu][j] = out[i] ; 
	}
#else
	// backwards
        #pragma omp parallel for private(i)
	PFOR( i = 0 ; i < LCU ; i++ ) {
	  out[i] = A[ slice + i ].O[mu][j] ; 
	}
	fftw_execute( backward ) ; 
        #pragma omp parallel for private(i)
	PFOR( i = 0 ; i < LCU ; i++ ) {
	  A[ slice + i ].O[mu][j] = in[i] ; 
	}
#endif
      }
  }

  //average into one sp_site size field
  // et voila! we have our fourier-transformed links in 0-2Pi BZ
  fftw_free( out ) ; 
  fftw_free( in ) ; 
  fftw_destroy_plan( forward ) ; 
  fftw_destroy_plan( backward ) ; 
#ifdef OMP_FFTW
  fftw_cleanup_threads( ) ;
#endif
  fftw_cleanup( ) ; 

  return GLU_SUCCESS ;
}

#else

static short int
mom_gauge_spatial( struct site *__restrict A ,
		   const lie_field_def def )
{
  fprintf( stderr , "[CUTS] WARNING! No FFT of the fields \n" ) ;
  return GLU_FAILURE ;
}

#endif

//will write our A field with cuts takes in our ( hopefully ) gauge-fixed links
int 
cuts_spatial ( struct site *__restrict A ,
	       const struct cut_info CUTINFO )
{
  fprintf( stdout , "\n[CUTS] Instantaneous Spatial/Temporal prop calculating ... \n" ) ; 
  // performs LOG and FFT of A, rewrites "A" to save space
  mom_gauge_spatial( A , CUTINFO.definition ) ; 

  if( check_psq( CUTINFO ) == GLU_FAILURE ) {
    return GLU_FAILURE ;
  }
  
  // calculate the momentum list
  int *num_mom = ( int* )malloc( sizeof( int ) ) ; 
  const struct veclist *list = compute_veclist( num_mom , CUTINFO ,
						ND-1 , GLU_FALSE ) ;

  // Landau correction required for Coulomb gauge too.
  correct_pspace_landau( A , list , num_mom , ND-1 ) ;

  // create some useful output and open the file
  char *str = output_str_struct( CUTINFO ) ;
  FILE *Aps = fopen( str , "wb" ) ; 

  write_mom_veclist( Aps , num_mom , list , ND - 1 ) ;

  //calculate the gluon prop and put into "g2_spatial" time-slice average it too
  /*          
    SPATIAL ::
			 Lt    ND - 1
                         ___     ___
                2        \       \
g( p ) =   ----------     >       >    Tr < A_( \mu )( p ) . A_( \mu )( -p ) >
          NCNC.(ND-2).V  /__     /__
                        t = 0   mu = 0

    TEMPORAL :: 
			 Lt  
                        ___   
                2       \    
g( p ) =   ----------    >    Tr < A_( ND - 1 )( p ) . A_( ND - 1 )( -p ) >
         NCNC.(ND-2).V  /__  
                       t = 0  

  */

  //this is the gluon propagator we will be writing
  double *g2_spatial = ( double* )malloc( num_mom[0] * sizeof( double ) ) ; 
  double *g2_temporal = ( double* )malloc( num_mom[0] * sizeof( double ) ) ; 

  // we average between the time slices to improve the signal some more ...
  const double g2_norm = ND > 2 ? 2.0 / ( ( NCNC - 1 ) * ( ND - 2 ) * LVOLUME ) : 2.0 / ( ( NCNC - 1 ) * LVOLUME ) ; 
  const double g0_norm = ND > 1 ? 2.0 / ( ( NCNC - 1 ) * ( ND - 1 ) * LVOLUME ) : 2.0 / ( ( NCNC - 1 ) * LVOLUME ) ; 

  // logically this is the zero momentum position as our list is symmetric
  const size_t midpoint = ( num_mom[0] - 1 ) >> 1 ;

  size_t i ;
#pragma omp parallel for private(i) 
  PFOR( i = 0  ;  i < num_mom[0]  ;  i++  ) {
    double sum_spatial = 0. ; 
    double sum_temporal = 0. ;
    GLU_complex tr ;
    size_t t , mu , mom , cnj ;
    for( t = 0 ;  t < Latt.dims[ ND - 1 ]  ;  t++ ) {
      mom = LCU * t + list[ i ].idx ;  
      cnj = LCU * t + list[ num_mom[0] - i - 1 ].idx ;  
      // prod all spatial directions
      for( mu = 0 ; mu < ND - 1 ; mu++ ) {      
	trace_ab( &tr , A[ mom ].O[mu] , A[ cnj ].O[mu] ) ;
	sum_spatial += (double)creal( tr ) ;
      }
      trace_ab( &tr , A[ mom ].O[ ND - 1 ] , A[ cnj ].O[ ND - 1 ] ) ;
      sum_temporal += (double)creal( tr ) ;
    }
      
    if( likely( i != midpoint ) ) {
      sum_spatial *= g2_norm ; 
      sum_temporal *= g2_norm ; 
    } else {	  
      sum_spatial *= g0_norm ; 
      sum_temporal *= g0_norm ; 
    }     
    
    g2_spatial[i] = sum_spatial ; 
    g2_temporal[i] = sum_temporal ; 
  }
  
  // write our fields to the file
  write_g2g3_to_list( Aps , g2_spatial , g2_temporal , num_mom ) ;
  
  // Tell us what we are outputting ...
  fprintf( stdout , "[CUTS] Writing finished...\n[CUTS] Outputting to %s \n" , str ) ;
 
  // FREE STUFF //
  fclose( Aps ) ;
  free( (void*)list ) ; 
  free( num_mom ) ; 
  free( str ) ;
  free( g2_spatial ) ;
  free( g2_temporal ) ;

  // Exit with a smile
  return GLU_SUCCESS ;
}


