#include "Mainfile.h"

#include "CG.h"          // routines used by both CG codes
#include "GLU_sums.h"    // round off resistant summations
#include "gftests.h"     // theta_test stopping condition
#include "gtrans.h"      // gauge transformations
#include "lin_derivs.h"  // linear approximation of lie matrices
#include "log_derivs.h"  // log-def of lie matrices 
#include "plaqs_links.h" // plaquette routine might be called
#include "invert.h"

static double zero_alpha = 0.0 ;

static int
FA_deriv2(  GLU_complex *__restrict *__restrict in , 
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

static void
exponentiate_gauge2( GLU_complex *__restrict *__restrict gauge , 
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

static void
steep_Landau_FALM( GLU_complex **gauge , 
		   struct site *lat ,
		   const void *forward , 
		   const void *backward , 
		   GLU_complex **out , 
		   GLU_complex **in , 
		   const GLU_real *psq , 
		   double *tr )
{ 
  // do a steepest-descents step with the result in "out"
  FA_deriv2( in , lat , tr ) ;

  // and do the fourier acceleration
  FOURIER_ACCELERATE( in , out , forward , backward ,
		      psq , LVOLUME ) ;

  // do stuff here when I work it out
  
  // exponentiate
  exponentiate_gauge2( gauge , (const GLU_complex**)in , 0.09 ) ;

  // do the gauge transform here
  gtransform( lat , ( const GLU_complex **)gauge ) ;

  printf( "THIS :: %f \n" , links( lat ) ) ;
  //exit(1) ;

  return ;
}

//returns the global gauge transform on lat
size_t
FALM( struct site *lat ,
      GLU_complex **gauge ,
      GLU_complex **out ,
      GLU_complex **in ,
      const void *forward ,
      const void *backward ,
      const GLU_real *p_sq ,
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
    steep_Landau_FALM( gauge , lat , forward , backward , 
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
      fprintf( stdout , "\n[GF] We will contunue this "
	       "instead of restarting ... \n" ) ;
      fprintf( stdout , "[GF] Link %1.15f || [GF_ ACC] %1.4e\n" , 
	       links( lat ) , *th ) ;
      // be wary, this could cause infinite tail recursion ...
      // in practice I don't think it will as we are probably over the hump
      return iters + FASD( lat , gauge , out , in , forward , backward , 
			   p_sq , th , acc , max_iters ) ;
    } else {
      return 0 ;
    }
  }

  //
  free_traces( ) ;

  //////////////////////////////////
  return iters ; 
}
