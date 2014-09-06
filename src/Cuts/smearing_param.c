/*
    Copyright 2013 Renwick James Hudspith

    This file (smearing_param.c) is part of GLU.

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
   @file smearing_param.c
   @brief computes the necessary gluon props for the smearing parameter f(q) 

   @ingroup Cuts
   writes the unsmeared gluon propagator and the smeared one
 */

#include "Mainfile.h"

#include "cut_routines.h"  // momentum cuts
#include "cut_output.h"    // output file
#include "geometry.h"      // general geometry for the p-calcs
#include "plan_ffts.h"     // FFTW plan wrappers
#include "plaqs_links.h"   // plaquettes and links calculations
#include "pspace_landau.h" // momentum space Landau correction
#include "SM_wrap.h"       // smearing operations wrapper

//#define WEAK_FIELD
//#define TADPOLE_IMPROVE

#ifdef HAVE_FFTW3_H

//def = 0 "lin" def = 1 "log"
static volatile int
mom_gauge( A , def )
     struct site *__restrict A ;
     const lie_field_def def ;
{
  if( parallel_ffts( ) == GLU_FAILURE ) {
    printf( "[PAR] Problem with initialising the OpenMP FFTW routines \n" ) ;
    return GLU_FAILURE ;
  }

  int i ;
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) { 
    GLU_complex temp[ NCNC ] ;
    int mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      switch( def )
	{
	case LINEAR_DEF :
	  Hermitian_proj( temp , A[i].O[mu] ) ;
	  break ;
	case LOG_DEF :
	  exact_log_slow( temp , A[i].O[mu] ) ;
	  break ;
	}
      equiv( A[i].O[mu] , temp ) ;
    }
  }
  
  /// FFTW routines
  GLU_complex *out = fftw_malloc( LVOLUME * sizeof( GLU_complex ) ) ;
  GLU_complex *in = fftw_malloc( LVOLUME * sizeof( GLU_complex ) ) ;

  fftw_plan forward , backward ;
  small_create_plans_DFT( &forward , &backward , in , out , ND ) ;

  ////// End of the search for Wisdom /////

  //forward transform
  int mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    int j ;
    for( j = 0 ; j < NCNC ; j++ ) {
      // forwards transform
      #ifdef CUT_FORWARD
      #pragma omp parallel for private(i)
      PFOR( i = 0 ; i < LVOLUME ; i++ ) {
	in[i] = A[i].O[mu][j] ;
      }
      fftw_execute( forward ) ;
      #pragma omp parallel for private(i)
      PFOR( i = 0 ; i < LVOLUME ; i++ ) {
	A[i].O[mu][j] = out[i] ;
      }
      #else
      // backwards transform
      #pragma omp parallel for private(i)
      PFOR( i = 0 ; i < LVOLUME ; i++ ) {
	out[i] = A[i].O[mu][j] ;
      }
      fftw_execute( backward ) ;
      #pragma omp parallel for private(i)
      PFOR( i = 0 ; i < LVOLUME ; i++ ) {
	A[i].O[mu][j] = in[i] ;
      }
      #endif
    }
  }

  //et voila! we have our fourier-transformed links in 0-2Pi BZ
  fftw_destroy_plan( forward ) ;
  fftw_destroy_plan( backward ) ;
  fftw_cleanup( ) ;
#ifdef OMP_FFTW
  fftw_cleanup_threads( ) ;
#endif
  fftw_free( out ) ;  
  fftw_free( in ) ;

  return GLU_SUCCESS ;
}

#else

static int
mom_gauge( A , def )
     struct site *__restrict A ;
     const lie_field_def def ;
{
  printf( "[CUTS] WARNING! Not performing an FFT \n" ) ;
  return GLU_FAILURE ;
}

#endif

// compute psq using the momentum lattice coordinates
static inline double
psq_calc( double mom[ ND ] ,
	  const int posit )
{
  int k[ ND ] , mu ;
  TwoPI_mpipi_momconv( k , posit , ND ) ;
  register double tspsq = 0. ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    #ifdef SIN_MOM
    mom[ mu ] = 2.0 * sin( k[ mu ] * Latt.twiddles[ mu ] * 0.5 ) ;
    #else
    mom[ mu ] = k[ mu ] * Latt.twiddles[ mu ] ;
    #endif
    tspsq += mom[ mu ] * mom[ mu ] ;
  }
  return ( tspsq < PREC_TOL ) ? 1.0 : tspsq ;
}

/**
   @fn static void project_trans_long( double *transverse , double *longitudinal , const struct site *__restrict A , const int posit ) 
  Transverse :
                
                  p_\mu p_\nu
  ( g_{\mu\nu} -  ----------- ) Tr[ G_{\mu\nu}(p^2) ] = G(p^2)
                      p^2

  Longitudinal :

  p_\mu p_\nu
  ----------- Tr[ G_{\mu\nu}(p^2) ] = F(p^2)
      p^2
 */
static void
project_trans_long( double *transverse ,
		    double *longitudinal ,
		    const struct site *__restrict A ,
		    const int posit ) 
{
  GLU_complex tr ;
  double mom[ ND ] , common ;
  // make this a constant
  const double spsq = 1.0 / psq_calc( mom , posit ) ;
  int mu , nu ;
  *transverse = *longitudinal = 0.0 ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    for( nu = 0 ; nu < ND ; nu++ ) {
      trace_ab_dag( &tr , A[posit].O[mu] , A[posit].O[nu] ) ;
      // common factor
      common = mom[mu] * mom[nu] * spsq * (double)creal(tr) ;
      // compute the longitudinal and transverse scalars
      *transverse += ( mu != nu ) ? -common : (double)creal(tr) - common ;
      *longitudinal += common ;
    }
  }
  return ;
}

#ifdef WEAK_FIELD
// used for comparison with pert-theory
static void
create_weak_field( struct site *__restrict A ) 
{
  int j ; 
#pragma omp parallel for private(j)
  for( j = 0 ; j < LVOLUME ; j++ ) {
    GLU_complex temp[ NCNC ] ;
    int mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      exact_log_slow( temp , A[j].O[mu] ) ;
      //Hermitian_proj( temp , A[j].O[mu] ) ;
      int k ;
      for( k = 0 ; k < NCNC ; k++ ) {
	temp[k] *= 0.01 ;
      }
      exponentiate( A[j].O[mu] , temp ) ;
    }
  }
  printf( "[CUTS] Weak field {p} %1.15f {l} %1.15f\n" , 
	  av_plaquette( A ) , links( A ) ) ; 
  return ;
}
#endif

// Computes the gluon propagators for smeared and unsmeared fields
int 
cuts_struct_smeared( struct site *__restrict A ,
		     const struct cut_info CUTINFO ,
		     const struct sm_info SMINFO )
{
  if( check_psq( CUTINFO ) == GLU_FAILURE ) {
    return GLU_FAILURE ;
  }

  // allocate room for the smeared gauge field
  struct site *SM_A = malloc( LVOLUME * sizeof ( struct site ) ) ;
  init_navig( SM_A ) ;

  // dirty, create a very weak field
#ifdef WEAK_FIELD
  create_weak_field( A ) ;
#endif

  int j ; 
#pragma omp parallel for private(j)
  for( j = 0 ; j < LVOLUME ; j++ ) {
    memcpy( &SM_A[j] , &A[j] , sizeof( struct site ) ) ; 
  }
  SM_wrap_struct( SM_A , SMINFO ) ;

#ifdef TADPOLE_IMPROVE
  const double u0 = 1.0 / sqrt( av_plaquette( A ) ) ;
  const double u0_SM = 1.0 / sqrt( av_plaquette( SM_A ) ) ;
  printf( "[CUTS] <WARNING> Tadpole improvement being used \n" ) ;
  printf( "[CUTS] Thin :: %f || Fat :: %f \n" , u0 , u0_SM ) ;
#endif

  // performs LOG and FFT , rewrites "A"
  mom_gauge( A , CUTINFO.definition ) ;
  mom_gauge( SM_A , CUTINFO.definition ) ;
  
  // this contains the length of the momentum list
  int *in = (int*)malloc( 1 * sizeof(int) ) ;
  const struct veclist *list = compute_veclist( in , CUTINFO , 
						ND , GLU_FALSE) ;

  // correct for the Landau condition, actually not necessary
  correct_pspace_landau( A , list , in , ND ) ;
  correct_pspace_landau( SM_A , list , in , ND ) ;

  // set up the outputs
  char *str = output_str_struct( CUTINFO ) ;  
  FILE *Ap = fopen( str , "wb" ) ; 

  // Here we write the two LANDAU gauge gluon propagators, 
  // unsmeared first and then smeared
  // traditional gluon 2pt function norms
  const double g2_norm = 2.0 / ( ( NCNC - 1 ) * ( ND - 1 ) * LVOLUME ) ;
  const double g0_norm = 2.0 / ( ( NCNC - 1 ) * ( ND ) * LVOLUME ) ;

  double *g2 = malloc( in[0] * sizeof( double ) ) ;
  double *g2_SM = malloc( in[0] * sizeof( double ) ) ;

  int i ;
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < in[0] ; i++ ) {
    // these require that this list is p -> -p symmetric, it is
    const int posit = list[i].idx ;
    const int conj  = list[ in[0] - i - 1 ].idx ;
    double trans , longitudinal ;

    // project out the transverse and longitudinal for the unsmeared case ...
    project_trans_long( &trans , &longitudinal , A , posit ) ;
    g2[i] = ( conj != posit ) ? trans * g2_norm : trans * g0_norm ;

    // projection for the smeared case ...
    project_trans_long( &trans , &longitudinal , SM_A , posit ) ;
    g2_SM[i] = ( conj != posit ) ? trans * g2_norm : trans * g0_norm ;

#ifdef TADPOLE_IMPROVE
    // our u0's are ( 1 / u0 ) ^ 2, one correction for each A-field 
    g2[i] *= ( u0 ) ;
    g2_SM[i] *= ( u0_SM ) ; // practically 1 
#endif
  }
  // create an appropriate output file
  write_mom_veclist( Ap , in , list , ND ) ;
  write_g2g3_to_list( Ap , g2 , g2_SM , in ) ;

  // tell us where the output file is ....
  printf( "[CUTS] Writing finished...\n[CUTS]Outputting to %s \n" , str ) ;

  // free the props
  free( g2 ) ;
  free( g2_SM ) ;

  // free everything else
  fclose( Ap ) ;
  free( (void*)list ) ;
  free( in ) ;
  free( str ) ;
  free( SM_A ) ;

  return GLU_SUCCESS ;
}

// undefine the weak field
#ifdef WEAK_FIELD
  #undef WEAK_FIELD
#endif

// undefine the tadpole improvement
#ifdef TADPOLE_IMPROVE
  #undef TADPOLE_IMPROVE
#endif
