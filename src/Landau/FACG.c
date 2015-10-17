/*
    Copyright 2013 Renwick James Hudspith

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

#ifdef HAVE_FFTW3_H
  #if (NC > 3) && (defined deriv_full)
  static const double al[LINE_NSTEPS] = { 0.0 , 0.1 , 0.2 , 0.3 } ;
  #else
  static const double al[LINE_NSTEPS] = { 0.0 , 0.15 , 0.3 } ;
  #endif
#else
static const double al[LINE_NSTEPS] = { 0.0 , 0.2 , 0.4 } ;
#endif

#ifdef GLU_FR_CG
  // minimum conjugate gradient step allowed
  const static double spmin = 0.05 ;
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
	     const int iters ) 
{
  if( iters % CAREFUL == 0 ) {
    *link = ( newlink == 0. ? 1. : newlink ) ; 
    printf( "[GF] ACC after %d iterations :: %1.5e \n" , iters , theta ) ;
    printf( "[GF] Link :: %1.15f || Plaq :: %1.15f \n" , 
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
	     const int iters )
{
  if( iters % CAREFUL == 0 ) {
    *newlink = links( lat ) ;
    printf("[GF] CHROMA ACC :: %1.15e \n" , 1. - link/(*newlink) ) ;
    printf("[GF] GAUGE ACC :: %1.15e \n" , gauge_test( gauge ) ) ;
    printf("[GF] THETA ACC :: %1.15e \n" , theta ) ;
  }
  return ;
}
#endif

/**
   @fn static void FA_deriv( struct site *__restrict lat , const void forward[ HERMSIZE ] , const void backward[ HERMSIZE ] , GLU_complex *__restrict *__restrict out , GLU_complex *__restrict *__restrict in , const GLU_real *psq , double *tr )
   @brief Fourier Accelerated steepest descent Landau gauge fixing routine
 */
static void 
FA_deriv(  GLU_complex *__restrict *__restrict in , 
	   struct site *__restrict lat ,
	   const GLU_real *psq , 
	   double *tr )
{
  size_t i ;
  // put the functional and the trace of the square of the deriv here
  double *alpha = malloc( LVOLUME * sizeof( double ) ) ;
  double *trAA = malloc( LVOLUME * sizeof( double ) ) ;
  const double fact = 1.0 / (double)( NC * ND ) ;
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {

    #if ( defined deriv_lin || defined deriv_linn )
    GLU_complex sum[ HERMSIZE ] ;
    #else
    GLU_complex sum[ HERMSIZE ] = { } ;
    #endif

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
    const double deriv = *tr > 0.1 ? \
      approx_log_deriv_nn( sum , lat , i , ND ) : log_deriv_nn( sum , lat , i , ND ) ; 
    // next nearest neighbour version
    #elif defined deriv_fullnn
    const double deriv = *tr > 0.1 ? \
      approx_log_deriv_nnn( sum , lat , i , ND ) : log_deriv_nnn( sum , lat , i , ND ) ;
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
    int mu ;
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
  return ;
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
  int i ;
  // the derivative in this form is antihermitian i.e -> i.dA
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    set_gauge_matrix( gauge[i] , in , alpha , i ) ;
  }
  return ;
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
  int counter = 0 ;
  double val[LINE_NSTEPS] ;
  val[0] = zero_alpha ; // 0 is a freebie
#ifdef verbose
  printf( "[GF] Landau CG probe :: %f %e \n" , al[0] , val[0] ) ;
#endif
  for( counter = 1 ; counter < LINE_NSTEPS ; counter++ ) {
    exponentiate_gauge( gauge , (const GLU_complex**)in , al[counter] ) ;
    val[counter] = evaluate_alpha( (const GLU_complex**) gauge , 
				   lat , ND , LVOLUME , 0 ) ;
    // the last argument of this has to be 0 !!! 
#ifdef verbose
    printf( "[GF] Landau CG probe :: %f %e \n" , al[counter] , val[counter] ) ;
#endif
  }

  // defined in CG.c
  return approx_minimum( LINE_NSTEPS , al , val ) ;
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
  FA_deriv( in , lat , psq , tr ) ;

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

#ifdef GLU_FR_CG
/**
   @fn static int steep_Landau_FACG( GLU_complex *__restrict *__restrict gauge , struct site *__restrict lat , const void forward[ HERMSIZE ] , const void backward[ HERMSIZE ] , GLU_complex *__restrict *__restrict out , GLU_complex *__restrict *__restrict in , const GLU_real *psq , double *tr , const double acc )
   @brief Fourier accelerated conjugate gradient method (Fletcher Reeves)
 */
static int
steep_Landau_FACG_FR( GLU_complex *__restrict *__restrict gauge , 
		      struct site *__restrict lat ,
		      const void *__restrict forward , 
		      const void *__restrict backward , 
		      GLU_complex *__restrict *__restrict out , 
		      GLU_complex *__restrict *__restrict in , 
		      GLU_complex *__restrict *__restrict in_old , 
		      GLU_complex *__restrict *__restrict sn ,
		      const GLU_real *psq , 
		      double *tr ,
		      const double acc )
{
  // perform an SD start
  steep_Landau_FA( gauge , lat , forward , backward , out , in , psq , tr ) ;

  int i ;
  #pragma omp parallel for private(i)
  for( i = 0 ; i < TRUE_HERM ; i++ ) {
    // copy sn to be "in" the FA derivative
    memcpy( sn[i] , in[i] , LVOLUME * sizeof( GLU_complex ) ) ;
  }

  // compute the quantity Tr( dA dA )
  double in_old_sum = sum_deriv( (const GLU_complex**)in , LVOLUME ) ;

  // initialise this to the one provided in the input file
  double spline_min = Latt.gf_alpha ;

  // loop a set number of CG-iterations
  int iters = 0 ;
  while( ( *tr > acc ) && ( iters > CG_MAXITERS ) && ( spline_min > spmin ) ) {

    // this ONLY works with out and in, make sure that is what we use
    FA_deriv( in , lat , psq , tr ) ;

    // and FA
    FOURIER_ACCELERATE( in , out , forward , backward ,
			psq , LVOLUME ) ;

    // summations for computing the scaling parameter
    const double insum = sum_deriv( (const GLU_complex**)in , LVOLUME ) ;

    // like the Fletcher-Reeves I suppose
    const double beta = in_old_sum > 0.0 ? insum / in_old_sum : 0.0 ; 
    in_old_sum = insum ;

    // compute sn, the conjugate matrix in FFTW's order
    #pragma omp parallel for private(i)
    for( i = 0 ; i < TRUE_HERM ; i++ ) {
      int j ;
      for( j = 0 ; j < LVOLUME ; j++ ) {
	sn[i][j] = in[i][j] + beta * ( sn[i][j] ) ;
      }
    }

    // do a line search
    if( *tr > PREC_TOL ) {
      spline_min = line_search( gauge , lat , (const GLU_complex**)sn ) ;
    } else {
      spline_min = Latt.gf_alpha ;
    }

    #ifdef verbose
    printf( "[GF] %d BETA %f \n" , iters , beta ) ;
    printf( "[GF] %d SPLINE2 :: %f \n" , iters , spline_min ) ;
    printf( "[GF] cgacc %e \n" , *tr ) ;
    #endif

    // and step the optimal amount
    exponentiate_gauge( gauge , (const GLU_complex**)sn , spline_min ) ;

    // do the gauge transform here
    gtransform( lat , ( const GLU_complex **)gauge ) ;

    // increment
    iters++ ;
  }
  return iters ;
}
#endif

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
		   const int max_iters )
{
  // perform an SD start
  steep_Landau_FA( gauge , lat , forward , backward , out , in , psq , tr ) ;

  int i ;
  #pragma omp parallel for private(i)
  for( i = 0 ; i < TRUE_HERM ; i++ ) {
    // copy sn to be "in" the FA derivative
    int j ;
    for( j = 0 ; j < LVOLUME ; j++ ) {
      sn[i][j] = in_old[i][j] = in[i][j] ;
    }
  }

  // compute the quantity Tr( dA dA )
  double in_old_sum = sum_deriv( (const GLU_complex**)in , LVOLUME ) ;

  // initialise this to the one provided in the input file
  double spline_min = Latt.gf_alpha ;

  // loop a set number of CG-iterations
  int iters = 0 ;
  while( ( *tr > acc ) && ( iters < max_iters ) ) {

    // this ONLY works with out and in, make sure that is what we use
    FA_deriv( in , lat , psq , tr ) ;

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
      int j ;
      for( j = 0 ; j < LVOLUME ; j++ ) {
	sn[i][j] = in[i][j] + beta * ( sn[i][j] ) ;
	in_old[i][j] = in[i][j] ;
      }
    }

    // do a line search
    if( *tr > CG_TOL ) {
      spline_min = line_search( gauge , lat , (const GLU_complex**)sn ) ;
    } else {
      spline_min = Latt.gf_alpha ;
    }

    #ifdef verbose
    printf( "[GF] %d BETA %e " , iters , beta ) ;
    printf( "[GF] SPLINE :: %e " , spline_min ) ;
    printf( "[GF] cgacc %e \n" , *tr ) ;
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
int 
FACG( struct site *__restrict lat , 
      GLU_complex *__restrict *__restrict gauge , 
      GLU_complex *__restrict *__restrict out ,
      GLU_complex *__restrict *__restrict in ,
      const void *__restrict forward , 
      const void *__restrict backward , 
      const GLU_real *__restrict p_sq , 
      double *th ,
      const double acc ,
      const int max_iters )
{
  // set up the maximum and what have you
  GLU_real max ; 
  *th = theta_test_lin( lat , &max , ND ) ; 
  int iters = 0 ;

  // have the option to leave early before allocations
  if( *th < acc ) return iters ;

  // allocate the "traces" array
  allocate_traces( LVOLUME ) ;

  // allocate the conjugate matrix
  GLU_complex **sn = malloc( TRUE_HERM * sizeof(GLU_complex*) ) ;
  GLU_complex **in_old = malloc( TRUE_HERM * sizeof( GLU_complex* ) ) ;

  // allocate the conjugate directions
  int i ;
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
    int loc_iters = 1 ;
    
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
      printf( "[GF] Continuation run \n" ) ;
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
      iters = GLU_FAILURE ;
    }
  } 
  //////////////////////////////////


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
int 
FASD( struct site *__restrict lat , 
      GLU_complex *__restrict *__restrict gauge , 
      GLU_complex *__restrict *__restrict out ,
      GLU_complex *__restrict *__restrict in ,
      const void *__restrict forward , 
      const void *__restrict backward , 
      const GLU_real *__restrict p_sq , 
      double *th ,
      const double acc ,
      const int max_iters )
{
  GLU_real max ; 
  *th = theta_test_lin( lat , &max , ND ) ; 
  int iters = 0 ; 

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
	return GLU_FAILURE ;
      }
    }
  }
  /////////////////////////////////
  // Have this loop just in case //
  if(  unlikely( iters == max_iters ) ) {
    if( *th < ( 1E3 * acc ) ) {
      printf( "\n[GF] We will contunue this instead of restarting ... \n" ) ;
      printf( "[GF] Link %1.15f || [GF_ ACC] %1.4e\n" , links( lat ) , *th ) ;
      // be wary, this could cause infinite tail recursion ...
      // in practice I don't think it will as we are probably over the hump
      return iters + FASD( lat , gauge , out , in , forward , backward , 
			   p_sq , th , acc , max_iters ) ;
    } else {
      return GLU_FAILURE ;
    }
  }

  //
  free_traces( ) ;

  //////////////////////////////////
  return iters ; 
}

// Fourier accelerated Steepest descents for the smeared fields
int 
FASD_SMEAR( struct site *__restrict lat , 
	    GLU_complex *__restrict *__restrict gauge , 
	    GLU_complex *__restrict *__restrict out ,
	    GLU_complex *__restrict *__restrict in ,
	    const void *__restrict forward , 
	    const void *__restrict backward , 
	    const GLU_real *__restrict p_sq , 
	    double *th ,
	    const double acc ,
	    const int max_iters )
{
  GLU_real max = 0. ; 
  *th = theta_test_lin( lat , &max , ND ) ; 
  int iters = 0 , i ; 
  
  //malloc temporary gauge
  GLU_complex **gauge2 = malloc( LVOLUME * sizeof( GLU_complex* ) ) ;
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    //gauge2[i] = ( GLU_complex* )malloc( NCNC * sizeof( GLU_complex ) ) ; 
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
      GLU_complex temp[ NCNC ] ;
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
  printf( "[GF] Using the following probes for the CG \n" ) ;
  int mu ;
  for( mu = 0 ; mu < LINE_NSTEPS ; mu++ ) {
    printf( "[GF] probe-%d %f \n" , mu , al[mu] ) ; 
  }
  return ;
}

// again keep local
#undef LINE_NSTEPS
