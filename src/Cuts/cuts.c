/*
Copyright 2013-2025 Renwick James Hudspith

    This file (cuts.c) is part of GLU.

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
   @file cuts.c
   @brief file for calculating the lie fields from our gauge links, fourier transforming them and discarding the unwanted momenta from a cutting procedure
 */
#include "Mainfile.h"

#include "cut_output.h"    // output file
#include "cut_routines.h"  // momentum cut routines
#include "glueprop.h"      // gluon propagator calculation
#include "MOMgg.h"         // exceptional (BOUCAUD,CHETYRKIN et al)
#include "MOMggg.h"        // non-exceptional scheme
#include "plan_ffts.h"     // FFTW planning wrappers
#include "pspace_landau.h" // momentum space Landau correction

#ifdef HAVE_FFTW3_H

static int
mom_gauge( struct site *__restrict A ,
	   const lie_field_def def )
{
  // callback for the log definition
  void (*log)( GLU_complex Q[ NCNC ] ,
	       const GLU_complex U[ NCNC ] ) = Hermitian_proj ;
  switch( def ) {
  case LINEAR_DEF :
    break ;
  case LOG_DEF : 
    log = exact_log_slow ; 
    break ;
  }

  size_t i ;
  #pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) { 
    GLU_complex temp[ NCNC ] GLUalign ;
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      log( temp , A[i].O[mu] ) ;
      equiv( A[i].O[mu] , temp ) ;
    }
  }

  // FFTW routines
  struct fftw_small_stuff FFTW ;
  small_create_plans_DFT( &FFTW , Latt.dims , ND ) ;
  
  //forward transform
  size_t mu , j ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    for( j = 0 ; j < NCNC ; j++ ) {
      // FORWARD ONE
      #ifdef CUT_FORWARD
      #pragma omp parallel for private(i)
      for( i = 0 ; i < LVOLUME ; i++ ) {
	FFTW.in[i] = A[i].O[mu][j] ;
      }
      fftw_execute( FFTW.forward ) ;
      #pragma omp parallel for private(i)
      for( i = 0 ; i < LVOLUME ; i++ ) {
	A[i].O[mu][j] = FFTW.out[i] ;
      }
      #else
      // BACKWARD TRANSFORM
      #pragma omp parallel for private(i)
      for( i = 0 ; i < LVOLUME ; i++ ) {
	FFTW.out[i] = A[i].O[mu][j] ;
      }
      fftw_execute( FFTW.backward ) ;
      #pragma omp parallel for private(i)
      for( i = 0 ; i < LVOLUME ; i++ ) {
	A[i].O[mu][j] = FFTW.in[i] ;
      }
      #endif
    }
  }
  //et voila! we have our fourier-transformed links in 0-2Pi BZ
  small_clean_up_fftw( FFTW ) ;

  return GLU_SUCCESS ;
}

#else

static int
mom_gauge( struct site *__restrict A ,
	   const lie_field_def def )
{
  fprintf( stderr , "[CUTS] WARNING! No FFT taking place\n" ) ;
  return GLU_FAILURE ;
}

#endif

//will write our A field with cuts takes in our (hopefully) gauge-fixed links
int 
cuts_struct( struct site *__restrict A ,
	     const struct cut_info CUTINFO )
{
  // performs LOG and FFT , rewrites "A"
  if( mom_gauge( A , CUTINFO.definition ) == GLU_FAILURE ) {
    fprintf( stderr , "[CUTS] FFT of the fields failed\n" ) ;
    return GLU_FAILURE ;
  }

  if( check_psq( CUTINFO ) == GLU_FAILURE ) {
    return GLU_FAILURE ;
  }

  // this contains the momentum list
  size_t *in = malloc( sizeof(size_t) ) ;
 
  struct cut_info TEMP = CUTINFO ;
  // NONEXCEPTIONAL MUST USE PSQ TO CONSERVE MOMENTUM
  if( CUTINFO.dir == NONEXCEPTIONAL ) {
    TEMP.type = PSQ_CUT ;
  }
  
  const struct veclist *list = compute_veclist( in , TEMP ,
						ND , GLU_FALSE ) ;

  // correct for the Landau condition
  correct_pspace_landau( A , list , in , ND ) ;
  
  char *str = output_str_struct( CUTINFO ) ;  
  FILE *Ap = fopen( str , "wb" ) ; 

  int check = 0 ;
  switch( CUTINFO.dir ) {
  case FIELDS :
    //all of this work leads to us discarding (most of) our A fields
    write_mom_veclist( Ap , in , list , ND ) ;
    write_lattice_fields( Ap , A , list , in ) ;
    break ;
  case EXCEPTIONAL :      
    write_mom_veclist( Ap , in , list , ND ) ;
    check = write_exceptional_g2g3_MOMgg( Ap , A , list , in ) ;
    break ;
  case NONEXCEPTIONAL :
    // we write the mom-list internally in the g2g3 code as it can
    // be different from the cut-momentum one..
    write_nonexceptional_g2g3( Ap , A , list , in , CUTINFO.max_mom ) ;
    break ;
  case GLUON_PROPS :
    write_mom_veclist( Ap , in , list , ND ) ;
    compute_gluon_prop( Ap , A , list , in ) ;
    break ;
  default : // default behaviour to the Exceptional triple glue
    write_mom_veclist( Ap , in , list , ND ) ;
    check = write_exceptional_g2g3_MOMgg( Ap , A , list , in ) ;
    break ;
  }
  fprintf( stdout , "[CUTS] Writing finished...\n"
	            "[CUTS] Outputting to %s\n" , str ) ;

  // Free allocated memory
  fclose( Ap ) ;
  free( (void*)list ) ;
  free( in ) ;
  free( str ) ;

  if( check == GLU_FAILURE ) {
    return GLU_FAILURE ;
  } else {
    return GLU_SUCCESS ;
  }
}
