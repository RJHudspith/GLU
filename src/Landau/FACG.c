/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (FACG.c) is part of GLU.

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
   @file FACG.c
   @brief conjugate gradient Fourier Accelerated Landau gauge fixing routines

   Is slightly different from the routines in sweet_FA.c although I could merge them. At the moment only set up for SU(3)
 */

#include "Mainfile.h"    // for all the definitions

#include "CG.h"          // routines used by both CG codes
#include "GLU_sums.h"    // round off resistant summations
#include "gftests.h"     // theta_test stopping condition
#include "gtrans.h"      // gauge transformations
#include "lin_derivs.h"  // linear approximation of lie matrices
#include "log_derivs.h"  // log-def of lie matrices 
#include "plaqs_links.h" // plaquette routine might be called

// these are the steps I choose, I need 5 points for the high-order cubic spline
// they DO NOT have to be evenly spaced, but DO have to be sorted
#if (NC > 3) && (defined deriv_full)
  #define LINE_NSTEPS 4
#else 
  #define LINE_NSTEPS 3
#endif

// if we are using the FACG the line search usually suggests a smaller set of probes
#ifdef HAVE_FFTW3_H
  #if (NC > 3) && (defined deriv_full)
  static const double al[LINE_NSTEPS] = { 0.0 , 0.1 , 0.2 , 0.3 } ;
  #else
  static const double al[LINE_NSTEPS] = { 0.0 , 0.15 , 0.3 } ;
  #endif
#else
static const double al[LINE_NSTEPS] = { 0.0 , 0.2 , 0.4 } ;
#endif

// zero alpha can be roughly computed in step
static double zero_alpha ;

// if we want to look at the fixing information
#ifdef CAREFUL
// helpers
static void
check_info1( const struct site *__restrict lat ,
	     double *link ,
	     const double newlink ,
	     const double theta , 
	     const size_t iters ) 
{
  if( iters % CAREFUL == 0 ) {
    *link = ( newlink == 0. ? 1. : newlink ) ; 
    fprintf( stdout , "[GF] ACC after %zu iterations :: %1.5e \n" , 
	     iters , theta ) ;
    fprintf( stdout , "[GF] Link :: %1.15f || Plaq :: %1.15f \n" , 
	     *link , av_plaquette( lat ) ) ; 
  }
  return ;
}

// and this one
static void
check_info2( const struct site *__restrict lat ,
	     const GLU_complex *__restrict *__restrict gauge ,
	     const double link ,
	     double *newlink ,
	     const double theta ,
	     const size_t iters )
{
  if( iters % CAREFUL == 0 ) {
    *newlink = links( lat ) ;
    fprintf( stdout , "[GF] CHROMA ACC :: %1.15e \n" , 1. - link/(*newlink) ) ;
    fprintf( stdout , "[GF] GAUGE ACC :: %1.15e \n" , gauge_test( gauge ) ) ;
    fprintf( stdout , "[GF] THETA ACC :: %1.15e \n" , theta ) ;
  }
  return ;
}
#endif

static int
FA_deriv(  GLU_complex *__restrict *__restrict in , 
	   struct site *__restrict lat , 
	   double *tr )
{
  size_t i ;
  // put the functional and the trace of the square of the deriv here
  double *alpha = NULL , *trAA ;
  if( GLU_malloc( (void**)&alpha , ALIGNMENT , LVOLUME * sizeof( double ) ) != 0 || 
      GLU_malloc( (void**)&trAA  , ALIGNMENT , LVOLUME * sizeof( double ) ) != 0 ) {
    return GLU_FAILURE ;
  }
  const double fact = 1.0 / (double)( NC * ND ) ;
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {

    GLU_complex sum[ HERMSIZE ] ;
    memset( &sum , 0 , HERMSIZE * sizeof( GLU_complex ) ) ;

    double functional = 0.0 ;
    #ifdef deriv_lin
    const double deriv = fast_deriv_AntiHermitian_proj( sum , &functional , lat , i ) ; 
    #elif defined deriv_linn
    const double deriv = fast_derivnn_AntiHermitian_proj( sum , lat , i ) ; 
    #elif defined deriv_full
    // OK, so if we do not use the approx log in the derivative we need to get the convergence going with the 
    // vandermonde definition, this idea was taken from arXiv:1010.5120v1
      #ifdef SINGLE_PREC
      const double deriv = *tr > 0.1 ?					\
	fast_deriv_AntiHermitian_proj( sum , &functional , lat , i )	\
        : log_deriv( sum , lat , i , ND ) ;
      #else
        #if NC < 4
        const double deriv = *tr > 0.1 ?				\
	  fast_deriv_AntiHermitian_proj( sum , &functional , lat , i )	\
	  : log_deriv( sum , &functional , lat , i , ND ) ; 
	#else
        const double deriv = log_deriv( sum , &functional , lat , i , ND ) ; 
        #endif
      #endif
    // nearest neighbour version
    #elif defined deriv_fulln    
    const double deriv = log_deriv_nn( sum , lat , i , ND ) ; 
    // next nearest neighbour version
    #elif defined deriv_fullnn
    const double deriv = *tr > 0.1 ? log_deriv_nnn( sum , lat , i , ND ) ;
    #endif

    // reductions
    trAA[i]  = deriv ;
    alpha[i] = functional * fact ;

    // make in anti-hermitian here!
    #if NC == 3
    // for SU(3) I pack the first FFT with the two explicitly real diagonal elements
    // to save on a Fourier transform ..
    in[0][i] = I * sum[0] - sum[3] ;
    in[1][i] = I * sum[1] ;
    in[2][i] = I * sum[2] ;
    in[3][i] = I * sum[4] ;
    #else
    size_t mu ;
    for( mu = 0 ; mu < HERMSIZE ; mu++ ) {
      in[mu][i] = I * sum[mu] ;
    }
    #endif
  }
#if ( defined deriv_lin ) || (defined deriv_linn )
  zero_alpha = 1.0 - kahan_summation( alpha , LVOLUME ) / (double)LVOLUME ;
#else
  #if NC < 4
  zero_alpha = *tr > 0.1 ? gauge_functional_fast( lat ) :\
    0.5 * kahan_summation( alpha , LVOLUME ) / (double)LVOLUME ;
  #else
  zero_alpha = 0.5 * kahan_summation( alpha , LVOLUME )  / (double)LVOLUME ;
  #endif
#endif
  *tr = kahan_summation( trAA , LVOLUME ) * ( GFNORM_LANDAU ) ; 
  free( alpha ) ; free( trAA ) ;
  return GLU_SUCCESS ;
}

/**
   @fn static void exponentiate_gauge( GLU_complex *__restrict *__restrict gauge , const GLU_complex *__restrict *__restrict in , const double alpha )
   @brief perform the exponentiation \f[g(x) = e^{alpha*\Delta_mu A_\mu}\f]
   @param gauge :: gauge transformation matrices generated by exponentiation
   @param in :: FFTW array carrying i.dA in shortened format
   @param alpha :: coupling
 */
static void
exponentiate_gauge( GLU_complex *__restrict *__restrict gauge , 
		    const GLU_complex *__restrict *__restrict in ,
		    const double alpha )
{
  size_t i ;
  // the derivative in this form is antihermitian i.e -> i.dA
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    set_gauge_matrix( gauge[i] , in , alpha , i ) ;
  }
  return ;
}

static double
lsearch( GLU_complex *__restrict *__restrict gauge , 
	 const struct site *__restrict lat ,
	 const GLU_complex *__restrict *__restrict in ,
	 const double *alphas ,
	 const size_t Nalphas )
{
  double val[ Nalphas ] ;
  size_t counter = 0 ;
  val[0] = zero_alpha ; // 0 is a freebie
#ifdef verbose
  fprintf( stdout , "[GF] Landau CG probe :: %f %1.15f \n" , alphas[0] , val[0] ) ;
#endif
  for( counter = 1 ; counter < Nalphas ; counter++ ) {
    exponentiate_gauge( gauge , (const GLU_complex**)in , alphas[counter] ) ;
    val[counter] = evaluate_alpha( (const GLU_complex**) gauge , 
				   lat , ND , LVOLUME , 0 ) ;
    // the last argument of this has to be 0 !!! 
#ifdef verbose
    fprintf( stdout ,  "[GF] Landau CG probe :: %f %1.15f \n" , 
	     alphas[counter] , val[counter] ) ;
#endif
  }
  // defined in CG.c
  return approx_minimum( Nalphas , alphas , val ) ;
}

/**
   @fn static double line_search( GLU_complex *__restrict *__restrict gauge , const struct site *__restrict lat , const GLU_complex *__restrict *__restrict in )
   @brief line search for the minimum of the gauge functional
 */
static double
line_search( GLU_complex *__restrict *__restrict gauge , 
	     const struct site *__restrict lat ,
	     const GLU_complex *__restrict *__restrict in )
{
  // for the early ones we have a bigger possible step to push us
  // deep down the initial gradient
  /*
  if( early_count < 3 ) {
    double al_early[ 5 ] = { 0.0 , 0.2 , 0.4 , 0.6 , 0.8 } ;
    early_count++ ;
    return lsearch( gauge , lat , in , al_early , 5 ) ;
  } else {
    return lsearch( gauge , lat , in , al , LINE_NSTEPS ) ;
  }
  */
  return lsearch( gauge , lat , in , al , LINE_NSTEPS ) ;
}

/**
   @fn static void steep_Landau_FA( GLU_complex *__restrict *__restrict gauge , struct site *__restrict lat , const void forward[ HERMSIZE ] , const void backward[ HERMSIZE ] , GLU_complex *__restrict *__restrict out , GLU_complex *__restrict *__restrict in , const GLU_real *psq , double *tr )
   @brief Fourier accelerated steepest descent method
 */
static void
steep_Landau_FA( GLU_complex *__restrict *__restrict gauge , 
		 struct site *__restrict lat ,
		 const void *__restrict forward , 
		 const void *__restrict backward , 
		 GLU_complex *__restrict *__restrict out , 
		 GLU_complex *__restrict *__restrict in , 
		 const GLU_real *psq , 
		 double *tr )
{
  // do a steepest-descents step with the result in "out"
  FA_deriv( in , lat , tr ) ;

  // and do the fourier acceleration
  FOURIER_ACCELERATE( in , out , forward , backward ,
		      psq , LVOLUME ) ;

  // and step length gf_alpha provided from the input file
#if (defined GLU_FR_CG) || (defined GLU_GFIX_SD)
  const double spline_min = Latt.gf_alpha ;
#else
  const double spline_min = line_search( gauge , lat , (const GLU_complex**)in ) ;
#endif

  // exponentiate
  exponentiate_gauge( gauge , (const GLU_complex**)in , spline_min ) ;

  // do the gauge transform here
  gtransform( lat , ( const GLU_complex **)gauge ) ;

  return ;
}

/**
   @fn static int steep_Landau_FACG( GLU_complex *__restrict *__restrict gauge , struct site *__restrict lat , const void forward[ HERMSIZE ] , const void backward[ HERMSIZE ] , GLU_complex *__restrict *__restrict out , GLU_complex *__restrict *__restrict in , const GLU_real *psq , double *tr , const double acc )
   @brief Fourier accelerated conjugate gradient method
 */
static int
steep_Landau_FACG( GLU_complex *__restrict *__restrict gauge , 
		   struct site *__restrict lat ,
		   const void *__restrict forward , 
		   const void *__restrict backward , 
		   GLU_complex *__restrict *__restrict out , 
		   GLU_complex *__restrict *__restrict in , 
		   GLU_complex *__restrict *__restrict in_old , 
		   GLU_complex *__restrict *__restrict sn ,
		   const GLU_real *psq , 
		   double *tr ,
		   const double acc ,
		   const size_t max_iters )
{
  // perform an SD start
  steep_Landau_FA( gauge , lat , forward , backward , out , in , psq , tr ) ;

  size_t i ;
  #pragma omp parallel for private(i)
  for( i = 0 ; i < TRUE_HERM ; i++ ) {
    // copy sn to be "in" the FA derivative
    size_t j ;
    for( j = 0 ; j < LVOLUME ; j++ ) {
      sn[i][j] = in_old[i][j] = in[i][j] ;
    }
  }

  // compute the quantity Tr( dA dA )
  double in_old_sum = sum_deriv( (const GLU_complex**)in , LVOLUME ) ;

  // loop a set number of CG-iterations
  size_t iters = 0 ;
  while( ( *tr > acc ) && ( iters < max_iters ) ) {

    // this ONLY works with out and in, make sure that is what we use
    FA_deriv( in , lat , tr ) ;

    // and FA
    FOURIER_ACCELERATE( in , out , forward , backward , psq , LVOLUME ) ;

    // summation of in * in, gets put in in_old_sum
    const double insum = sum_deriv( (const GLU_complex**)in , LVOLUME ) ;

    // is the sum of the trace of in * in_old
    const double sum_conj = sum_PR_numerator( (const GLU_complex**)in , 
					      (const GLU_complex**)in_old , 
					      LVOLUME ) ;

    // compute the beta value, who knows what value is best?
    double beta = PRfmax( 0.0 , ( sum_conj ) / in_old_sum ) ;

    // switch to the fletcher reeves
    if( *tr < CG_TOL ) {
      beta = insum / in_old_sum ;
    }
    in_old_sum = insum ;

    // compute sn, the conjugate matrix in FFTW's order
    #pragma omp parallel for private(i)
    for( i = 0 ; i < TRUE_HERM ; i++ ) {
      size_t j ;
      for( j = 0 ; j < LVOLUME ; j++ ) {
	sn[i][j] = in[i][j] + beta * ( sn[i][j] ) ;
	in_old[i][j] = in[i][j] ;
      }
    }

    // do a line search
    double spline_min = Latt.gf_alpha ;
    if( *tr > CG_TOL ) {
      spline_min = line_search( gauge , lat , (const GLU_complex**)sn ) ;
    }

    #ifdef verbose
    fprintf( stdout , "[GF] %zu BETA %e " , iters , beta ) ;
    fprintf( stdout , "[GF] SPLINE :: %e " , spline_min ) ;
    fprintf( stdout , "[GF] cgacc %e \n" , *tr ) ;
    #endif

    // spline_min == 0 is a special case where there looks like
    // there is no minimum
    if( spline_min == 0.0 ) continue ;

    // and step the optimal amount
    exponentiate_gauge( gauge , (const GLU_complex**)sn , spline_min ) ;

    // do the gauge transform here
    gtransform( lat , ( const GLU_complex **)gauge ) ;

    // increment
    iters++ ;
  }

  return iters ;
}

// overwrites the lattice links in lat to Landau gauge fixed links using CG
size_t
FACG( struct site *__restrict lat , 
      GLU_complex *__restrict *__restrict gauge , 
      GLU_complex *__restrict *__restrict out ,
      GLU_complex *__restrict *__restrict in ,
      const void *__restrict forward , 
      const void *__restrict backward , 
      const GLU_real *__restrict p_sq , 
      double *th ,
      const double acc ,
      const size_t max_iters )
{
  // set up the maximum and what have you
  GLU_real max ; 
  *th = theta_test_lin( lat , &max , ND ) ; 
  size_t iters = 0 ;
  
  // have the option to leave early before allocations
  if( *th < acc ) return iters ;
  
  // allocate the conjugate matrices
  GLU_complex **sn     = NULL ;
  GLU_complex **in_old = NULL ;

  if( GLU_malloc( (void**)&sn     , 16 , TRUE_HERM * sizeof( GLU_complex* ) ) != 0 ||
      GLU_malloc( (void**)&in_old , 16 , TRUE_HERM * sizeof( GLU_complex* ) ) != 0 ) {
    fprintf( stderr , "[GF] FACG failed to allocate CG temporaries\n" ) ;
    return GLU_FAILURE ;
  }

  // allocate the "traces" array
  allocate_traces( LVOLUME ) ;

  // allocate the conjugate directions
  size_t i ;
  #pragma omp parallel for private(i)
  for( i = 0 ; i < TRUE_HERM ; i++ ) {
    GLU_malloc( (void**)&sn[i]     , 16 , LVOLUME * sizeof( GLU_complex ) ) ;
    GLU_malloc( (void**)&in_old[i] , 16 , LVOLUME * sizeof( GLU_complex ) ) ;
  }

  #ifdef CAREFUL
  double newlink = 0. , link = 0. ; 
  #endif

  // I believe the algorithm wastes its time early on
  // these steps are better being SD'd
  // loop the CG
  while( ( *th > acc ) && ( iters < max_iters ) ) {
    #ifdef CAREFUL
    check_info1( lat , &link , newlink , *th , iters ) ;
    #endif
    size_t loc_iters = 1 ;
    
    // if the spline just keeps giving zeros we switch to the SD
    // stopping this method wasting its time doing meaningless line searches
    #ifdef GLU_FR_CG
    loc_iters = steep_Landau_FACG_FR( gauge , lat , forward , backward , 
				      out , in , in_old , sn , p_sq , th , acc ) ; 
    #else
    loc_iters = steep_Landau_FACG( gauge , lat , forward , backward , 
				   out , in , in_old , sn , p_sq , th , acc , 
				   max_iters ) ; 
    #endif

    // and increment the iteration counter
    iters += loc_iters ;

    #ifdef CAREFUL
    check_info2( lat , ( const GLU_complex **)gauge , link , &newlink , *th , iters ) ;
    #endif
  }

  /////////////////////////////////
  // Have this loop just in case //
  if(  unlikely( iters >= max_iters ) ) {
    if( ( *th < 1E3*acc ) ) {
      const int sumiters = iters ;
      iters = 0 ; 
      fprintf( stdout , "[GF] Continuation run \n" ) ;
      while( ( *th > acc ) && ( iters < max_iters ) ) {
	iters = steep_Landau_FACG( gauge , lat , forward , backward , 
				   out , in , in_old , sn , p_sq , th , acc ,
				   max_iters ) ; 
      }
      // if the continuation doesn't work we restart
      if( *th > acc ) return GLU_FAILURE ;
      // otherwise we sum up the total iters
      iters += sumiters ;
    } else {
      iters = 123456789 ;
    }
  } 
  
  // free the traces array
  free_traces( ) ;

  // free the conjugate-direction derivative
  #pragma omp parallel for private(i) 
  for( i = 0 ; i < TRUE_HERM ; i++ ) {
    free( in_old[i] ) ;
    free( sn[i] ) ;
  }
  free( in_old ) ;
  free( sn ) ;

  return iters ; 
}

//returns the global gauge transform on lat
size_t
FASD( struct site *__restrict lat ,
      GLU_complex *__restrict *__restrict gauge ,
      GLU_complex *__restrict *__restrict out ,
      GLU_complex *__restrict *__restrict in ,
      const void *__restrict forward ,
      const void *__restrict backward ,
      const GLU_real *__restrict p_sq ,
      double *th ,
      const double acc ,
      const size_t max_iters )
{
  GLU_real max ; 
  *th = theta_test_lin( lat , &max , ND ) ; 
  size_t iters = 0 ; 

  // have a reference for the plaquette to check convergence
  const double ref_plaquette = av_plaquette( lat ) ;

  // allocate the "traces" array
  allocate_traces( LVOLUME ) ;

  #ifdef CAREFUL
  double newlink = 0. , link = 0. ; 
  #endif

  while(  ( *th > acc ) && ( iters < max_iters ) ) {
    #ifdef CAREFUL
    check_info1( lat , &link , newlink , *th , iters ) ;
    #endif
    iters++ ; 
    steep_Landau_FA( gauge , lat , forward , backward , 
		     out , in , p_sq , th ) ; 
    #ifdef CAREFUL
    check_info2( lat , ( const GLU_complex **)gauge , link , &newlink , *th , iters ) ;
    #endif
    /////// have a check for ill-convergence, which is related to plaquette-drift /////
    if( ( iters&127 ) == 0 ) {
      // roughly, I will allow for a deviation around 0.1*PREC_TOL per iteration
      // for the average plaquette
      if( fabs( av_plaquette( lat ) - ref_plaquette ) > (0.1 * PREC_TOL * iters) ) {
	fprintf( stderr , "[GF] severe plaquette drift found\n" ) ;
	return 123456789 ;
      }
    }
  }
  /////////////////////////////////
  // Have this loop just in case //
  if(  unlikely( iters == max_iters ) ) {
    if( *th < ( 1E3 * acc ) ) {
      fprintf( stdout , "\n[GF] We will contunue this "
	       "instead of restarting ... \n" ) ;
      fprintf( stdout , "[GF] Link %1.15f || [GF_ ACC] %1.4e\n" , 
	       links( lat ) , *th ) ;
      // be wary, this could cause infinite tail recursion ...
      // in practice I don't think it will as we are probably over the hump
      return iters + FASD( lat , gauge , out , in , forward , backward , 
			   p_sq , th , acc , max_iters ) ;
    } else {
      return 123456789 ;
    }
  }

  //
  free_traces( ) ;

  //////////////////////////////////
  return iters ; 
}

// Fourier accelerated Steepest descents for the smeared fields
size_t
FASD_SMEAR( struct site *__restrict lat ,
	    GLU_complex *__restrict *__restrict gauge ,
	    GLU_complex *__restrict *__restrict out ,
	    GLU_complex *__restrict *__restrict in ,
	    const void *__restrict forward ,
	    const void *__restrict backward ,
	    const GLU_real *__restrict p_sq ,
	    double *th ,
	    const double acc ,
	    const size_t max_iters )
{
  GLU_real max = 0. ; 
  *th = theta_test_lin( lat , &max , ND ) ; 
  size_t iters = 0 , i ; 
  
  //malloc temporary gauge
  GLU_complex **gauge2 = NULL ;
  if( GLU_malloc( (void**)&gauge2 , 16 , LVOLUME * sizeof( GLU_complex* ) ) != 0 ) {
    fprintf( stderr , "[GF] FASD_SMEAR failed to allocate temporary gauge\n" ) ;
    return GLU_FAILURE ;
  }

  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    GLU_malloc( (void**)&gauge2[i] , 16 , NCNC * sizeof( GLU_complex ) ) ;
    identity( gauge2[i] ) ;
  }

  #ifdef CAREFUL
  double newlink = 0. , link = 0. ; 
  #endif
  while( *th > acc && iters < max_iters) {
    #ifdef CAREFUL
    check_info1( lat , &link , newlink , *th , iters ) ;
    #endif
    iters++ ; 
    steep_Landau_FA( gauge , lat , forward , backward , 
		     out , in , p_sq , th ) ; 
    // multiply through 
    #pragma omp parallel for private(i) 
    PFOR(  i = 0 ; i < LVOLUME ; i++  ) {
      GLU_complex temp[ NCNC ] GLUalign ;
      memcpy( temp , gauge2[i] , NCNC * sizeof( GLU_complex ) ) ;
      multab_suNC( gauge2[i] , gauge[i] , temp ) ; 
    }
    #ifdef CAREFUL
    check_info2( lat , ( const GLU_complex **)gauge , link , &newlink , *th , iters ) ;
    #endif
  }
  #pragma omp parallel for private(i)
  for(  i = 0 ; i < LVOLUME ; i++  ) {
    equiv( gauge[i] , gauge2[i] ) ; 
    free( gauge2[i] ) ; 
  }
  free( gauge2 ) ; 
  return iters ; 
}

// tells us the probes we are using
void
query_probes_Landau( void ) {
  fprintf( stdout , "[GF] Using the following probes for the CG \n" ) ;
  size_t mu ;
  for( mu = 0 ; mu < LINE_NSTEPS ; mu++ ) {
    fprintf( stdout , "[GF] probe-%zu %f \n" , mu , al[mu] ) ; 
  }
  return ;
}

// again keep local
#undef LINE_NSTEPS
