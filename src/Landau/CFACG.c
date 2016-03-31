/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (CFACG.c) is part of GLU.

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
   @file CFACG.c
   @brief Coulomb Gauge Fourier Accelerated Conjugate Gradient  
 */

#include "Mainfile.h"      // for all the definitions

#include "CG.h"            // for some of the shared routines
#include "GLU_sums.h"      // round-off resistant sums
#include "gftests.h"       // theta_test stopping condition
#include "gramschmidt.h"   // for reunit2, gram-schmidt reunitarisation
#include "gtrans.h"        // gauge transformations
#include "plaqs_links.h"   // plaquette routine might be called
#include "random_config.h" // for the random transformed gauge slices

// these are the steps I choose, I need 5 points for the high-order cubic spline
// they DO NOT have to be evenly spaced, but DO have to be sorted
#define LINE_NSTEPS 3

#ifdef HAVE_FFTW3_H
static const double alcg[LINE_NSTEPS] = { 0.0 , 0.17 , 0.34 } ;
#else
static const double alcg[LINE_NSTEPS] = { 0.0 , 0.2 , 0.4 } ;
#endif

// computed in-step when the derivative is taken
static double zero_alpha = 1.0 ;
// start with driving the SD and build up from there ...
static void
exponentiate_gauge_CG( GLU_complex *__restrict *__restrict gauge , 
		       const GLU_complex *__restrict *__restrict in ,
		       const double alpha )
{
  size_t i ;
  // the derivative in this form is antihermitian i.e -> i.dA
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LCU ; i++ ) {
    GLU_complex temp[ NCNC ] GLUalign , temp2[ NCNC ] GLUalign ;
    set_gauge_matrix( temp , in , alpha , i ) ;
    equiv( temp2 , gauge[i] ) ;
    multab_suNC( gauge[i] , temp , temp2 ) ;
  }
  return ;
}

// prints out some relevant information ...
static void 
get_info( const size_t t , const double tr , const size_t iters , 
	  const int control , const double accuracy , 
	  const GLU_bool switched )
{
  if( t == 0 ) {
    fprintf( stdout , "\n" ) ;
  } if( tr < accuracy ) {
    fprintf( stdout , "[GF] Slice :: %zu {Stopped by convergence} "
	     "\n[GF] Accuracy :: %1.5e"
	     " || Iterations :: %zu\n[GF] Failures :: %d\n" , 
	     t , tr , iters , control ) ; 
  } else {
    fprintf( stdout , "[GF] Slice :: %zu {Stopped by iterations too high} \n"
	     "[GF] Accuracy :: %1.5e || Iterations :: %zu\n"
	     "[GF] Failures :: %d \n" , 
	     t , tr , iters , control ) ; 
  }
  // if we perform an SD switch -> should we ever?
  if( switched == GLU_TRUE ) {
    fprintf( stdout , "[GF] Switched from CG to SD \n" ) ;
  }
  fprintf( stdout , "\n" ) ; 
  return ;
}

// perform a line search using GLUbic splines for approximately the best alpha
// which it returns
static double
line_search_Coulomb( GLU_complex *__restrict *__restrict gtransformed , 
		     const GLU_complex *__restrict *__restrict gauge , 
		     const struct site *__restrict lat ,
		     const GLU_complex *__restrict *__restrict in ,
		     const size_t t )
{
  // holds the probe evaluations
  double val[LINE_NSTEPS] ;

  // this has been computed with the derivative
  val[0] = zero_alpha ;

  // loop the evaluations
  size_t counter ;
  for( counter = 1 ; counter < LINE_NSTEPS ; counter++ ) {
    size_t i ;
    #pragma omp parallel for private(i)
    for( i = 0 ; i < LCU ; i++ ) {
      GLU_complex temp[ NCNC ] ;
      set_gauge_matrix( temp , in , alcg[counter] , i ) ;
      multab_suNC( gtransformed[i] , temp , gauge[i] ) ;
    }
    // and compute the functional
    val[counter] = evaluate_alpha( (const GLU_complex**)gtransformed , 
				   lat , ND-1 , LCU , t ) ;
  }
#ifdef verbose
  printf( "[SPLINE] MIN :: %f \n" , 
	  approx_minimum( LINE_NSTEPS , alcg , val ) ) ;
#endif
  // defined in CG.c
  return approx_minimum( LINE_NSTEPS , alcg , val ) ;
}

// bit that calculates the steepest ascent with fourier acceleration
static void
steep_deriv_CG( GLU_complex *__restrict *__restrict in ,
		struct sp_site_herm *__restrict rotato ,
		const struct site *__restrict lat , 
		const GLU_complex *__restrict *__restrict slice_gauge , 
		const size_t t ,
		double *tr )
{
  // sets the gtrans'd field into the temporary "rotato"
  zero_alpha = coul_gtrans_fields( rotato , lat , slice_gauge , t ) ;

  // gauge transform the whole slice
  double *trAA = malloc( LCU * sizeof( double ) ) ;
  size_t i ;
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LCU ; i ++ ) {
    // compute gauge rotated derivatives
    GLU_complex sum[ HERMSIZE ] ;
    size_t mu ;
    for( mu = 0 ; mu < HERMSIZE ; mu++ ) {
      sum[ mu ] = 0.0 ;
    }
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      const int back = lat[i].back[mu] ;
      #if NC == 3
      sum[0] += rotato[back].O[mu][0] - rotato[i].O[mu][0] ;
      sum[1] += rotato[back].O[mu][1] - rotato[i].O[mu][1] ;
      sum[2] += rotato[back].O[mu][2] - rotato[i].O[mu][2] ;
      sum[3] += rotato[back].O[mu][3] - rotato[i].O[mu][3] ;
      sum[4] += rotato[back].O[mu][4] - rotato[i].O[mu][4] ;
      #elif NC == 2
      sum[0] += rotato[back].O[mu][0] - rotato[i].O[mu][0] ;
      sum[1] += rotato[back].O[mu][1] - rotato[i].O[mu][1] ;
      #else 
      size_t nu ;
      for( nu = 0 ; nu < HERMSIZE ; nu++ ) {
	sum[ nu ] += rotato[back].O[mu][nu] - rotato[i].O[mu][nu] ;
      }
      #endif
      // nearest-neighbour derivatives
      #if ( defined deriv_linn ) || ( defined deriv_fulln ) 
      const size_t back2 = lat[back].back[mu] ;
      const size_t forw2 = lat[i].neighbor[mu] ;
        #if NC == 3
        sum[0] += rotato[back2].O[mu][0] - rotato[forw2].O[mu][0] ;
        sum[1] += rotato[back2].O[mu][1] - rotato[forw2].O[mu][1] ;
        sum[2] += rotato[back2].O[mu][2] - rotato[forw2].O[mu][2] ;
        sum[3] += rotato[back2].O[mu][3] - rotato[forw2].O[mu][3] ;
        sum[4] += rotato[back2].O[mu][4] - rotato[forw2].O[mu][4] ;
        #elif NC == 2 
	sum[0] += rotato[back2].O[mu][0] - rotato[forw2].O[mu][0] ;
        sum[1] += rotato[back2].O[mu][1] - rotato[forw2].O[mu][1] ;
        #else
        for( nu = 0 ; nu < HERMSIZE ; nu++ ) {
	  sum[ nu ] += rotato[back2].O[mu][nu] - rotato[forw2].O[mu][nu] ;
        }
        #endif
      #endif
    }
    register double der = 0.0 ;
    #if NC == 3
    der += creal( sum[0] * sum[0] ) ;
    der += creal( sum[1] * conj( sum[1] ) ) ;
    der += creal( sum[2] * conj( sum[2] ) ) ;
    der += creal( sum[3] * sum[3] ) ;
    der += creal( sum[4] * conj( sum[4] ) ) ;
    der += creal( sum[0] * sum[3] ) ;
    #elif NC == 2
    der += creal( sum[0] * sum[0] ) ;
    der += creal( sum[1] * conj( sum[1] ) ) ;
    #else
    GLU_real tr ;
    trace_ab_herm_short( &tr , sum , sum ) ;
    der = (double)tr ;
    #endif
    // The trace has a factor of 2, as does the def 2 tr(AA) = \delta_ab A^a A^b
    der *= 4.0 ;

    // make it antihermitian
    #if NC == 3
    // for SU(3) I pack the first FFT with the two explicitly
    // real diagonal elements
    // to save on a Fourier transform ..
    in[0][i] = I * sum[0] - sum[3] ;
    in[1][i] = I * sum[1] ;
    in[2][i] = I * sum[2] ;
    in[3][i] = I * sum[4] ;
    #elif NC == 2
    in[0][i] = I * sum[0] ;
    in[1][i] = I * sum[1] ;
    #else
    for( mu = 0 ; mu < HERMSIZE ; mu++ ) {
      in[mu][i] = I * sum[mu] ;
    }
    #endif
    trAA[i] = der ;
  }
 
  // normalise the measure
  *tr = kahan_summation( trAA , LCU ) * GFNORM_COULOMB ; 
  free( trAA ) ;
  return ;
}

// perform a Fourier-accelerated SD step
static void
steep_step_SD( const struct site *__restrict lat , 
	       GLU_complex *__restrict *__restrict slice_gauge , 
	       GLU_complex *__restrict *__restrict out , 
	       GLU_complex *__restrict *__restrict in , 
	       struct sp_site_herm *__restrict rotato ,
	       GLU_complex *__restrict *__restrict gtransformed , 
	       const void *__restrict forward , 
	       const void *__restrict backward , 
	       const GLU_real *__restrict psq , 
	       const size_t t ,
	       double *tr )
{
  // compute deriv and Fourier Accelerate
  steep_deriv_CG( in , rotato , lat , (const GLU_complex**)slice_gauge , t , tr ) ;

  // if we want to Fourier accelerate, we call this otherwise it is the SD
  FOURIER_ACCELERATE( in , out , forward , backward , psq , LCU ) ;

  double spline_min = Latt.gf_alpha ;
  if( gtransformed != NULL ) {
#if !(defined GLU_GFIX_SD) && !(defined GLU_FR_CG)
    spline_min = line_search_Coulomb( gtransformed , 
				      (const GLU_complex**)slice_gauge , lat , 
				      (const GLU_complex**)in , t ) ;
#endif
  }
  // exponentiate this
  exponentiate_gauge_CG( slice_gauge , (const GLU_complex**)in , spline_min ) ;
  return ;
}

// more up to date Polyak-Ribiere CG code
static size_t
steep_step_FACG( const struct site *__restrict lat , 
		 GLU_complex *__restrict *__restrict gauge , 
		 GLU_complex *__restrict *__restrict out , 
		 GLU_complex *__restrict *__restrict in , 
		 GLU_complex *__restrict *__restrict in_old ,
		 GLU_complex *__restrict *__restrict sn ,
		 struct sp_site_herm *__restrict rotato ,
		 GLU_complex *__restrict *__restrict gtransformed ,
		 const void *__restrict forward , 
		 const void *__restrict backward , 
		 const GLU_real *__restrict psq , 
		 const size_t t ,
		 double *tr ,
		 const double accuracy , 
		 const size_t max_iters )
{
  // SD start
  steep_step_SD( lat , gauge , out , in , rotato , gtransformed ,
		 forward , backward , psq , t , tr ) ;

  size_t i ;
  #pragma omp parallel for private(i)
  for( i = 0 ; i < TRUE_HERM ; i++ ) {
    // copy sn to be "in" the FA derivative
    size_t j ;
    for( j = 0 ; j < LCU ; j++ ) {
      sn[i][j] = in_old[i][j] = in[i][j] ;
    }
  }

  // compute the quantity Tr( dA dA )
  double in_old_sum = sum_deriv( (const GLU_complex**)in , LCU ) ;

  // loop a set number of CG-iterations
  size_t iters = 0 ;
  while( ( *tr > accuracy ) && ( iters < max_iters ) ) {

    // this ONLY works with in, make sure that is what we use
    steep_deriv_CG( in , rotato , lat , (const GLU_complex**)gauge , t , tr ) ;

    // if we want to Fourier accelerate, we call this otherwise it is the SD
    FOURIER_ACCELERATE( in , out , forward , backward , psq , LCU ) ;

    // summations for computing the scaling parameter
    const double insum = sum_deriv( (const GLU_complex**)in , LCU ) ;

    // is the sum of the trace of in * in_old
    const double sum_conj = sum_PR_numerator( (const GLU_complex**)in , 
					      (const GLU_complex**)in_old , 
					      LCU ) ;

    // compute the beta value, who knows what value is best?
    double beta = PRfmax( 0.0 , sum_conj / in_old_sum ) ;

    // switch to the fletcher reeves
    if( *tr < CG_TOL) {
      beta = insum / in_old_sum ;
    }
    in_old_sum = insum ;

    // compute sn, the conjugate matrix in FFTW's order
    #pragma omp parallel for private(i)
    for( i = 0 ; i < TRUE_HERM ; i++ ) {
      size_t j ;
      for( j = 0 ; j < LCU ; j++ ) {
	sn[i][j] = in[i][j] + beta * ( sn[i][j] ) ;
	in_old[i][j] = in[i][j] ;
      }
    }

    // line search up to a tolerance
    double spline_min = Latt.gf_alpha ;
    if( *tr > CG_TOL ) { 
      spline_min = line_search_Coulomb( gtransformed , 
					(const GLU_complex**)gauge , lat , 
					(const GLU_complex**)sn , t ) ;
    }

    #ifdef verbose
    fprintf( stdout , "[GF] %zu BETA %e \n" , iters , beta ) ;
    fprintf( stdout , "[GF] %zu SPLINE2 :: %e \n" , iters , spline_min ) ;
    fprintf( stdout , "[GF] cgacc %e \n" , *tr ) ;
    #endif

    // no point in performing this rotation, 0.0 is like a flag for a broken
    // interpolation
    if( spline_min == 0.0 ) continue ;

    // and step the optimal amount multiplying atomically into the gauge matrices
    exponentiate_gauge_CG( gauge , (const GLU_complex**)sn , spline_min ) ;

    // increment
    iters++ ;
  }
  return iters ;
}

// coulomb gauge fix on slice t
static size_t
steep_fix_FACG( const struct site *__restrict lat , 
		GLU_complex *__restrict *__restrict slice_gauge , 
		GLU_complex *__restrict *__restrict out , 
		GLU_complex *__restrict *__restrict in , 
		GLU_complex *__restrict *__restrict in_old ,
		GLU_complex *__restrict *__restrict sn , 
		struct sp_site_herm *__restrict rotato ,
		GLU_complex *__restrict *__restrict gtransformed , 
		const void *__restrict forward , 
		const void *__restrict backward , 
		const GLU_real *__restrict psq , 
		const size_t t ,
		const double accuracy ,
		const size_t max_iter )
{
  double tr = 1.0 ;
  size_t iters = 0 , temp_iters = 0 , control = 0 ;
  GLU_bool switched = GLU_FALSE , continuation = GLU_FALSE ;
  while ( ( tr > accuracy ) && ( iters < max_iter ) ) {

    // perform a Fourier accelerated step, this is where the CG can go ...
    size_t loc_iters = 
      steep_step_FACG( lat , slice_gauge , out , in , in_old , sn , 
		       rotato , gtransformed , forward , backward , 
		       psq , t , &tr , accuracy , max_iter ) ;

    // if we have gone over the max number of iters allowed
    // we randomly transform
    if( ( iters + loc_iters ) >= ( max_iter - 1 ) ) {
      // if we are close we continue, no use in wasting good work
      if ( tr < ( accuracy * 1E3 ) && continuation == GLU_FALSE ) {
	fprintf( stdout , "[GF] Continuation run %e \n" , tr ) ;
	temp_iters += ( iters + loc_iters ) ;
	loc_iters = iters = 0 ;
	continuation = GLU_TRUE ;
	// otherwise we randomly restart
      } else if( control < GF_GLU_FAILURES && tr > accuracy ) {
	fprintf( stdout , "[GF] Random transform \n" ) ;
	random_gtrans_slice( slice_gauge ) ;
	temp_iters += ( iters + loc_iters ) ;
	loc_iters = iters = 0 ;
	control++ ;
	tr = 1.0 ;
	// otherwise we have to leave
      } else {
	break ;
      }
      // and that is the end
    }

    // count some more iterations
    iters += loc_iters ;

    // end of the while loop
  }
  // and tell us about the slice's iteration
  get_info( t , tr , iters + temp_iters , control , accuracy , switched ) ; 

  return iters + temp_iters ;
}

// Fourier accelerated CG code is here
size_t
Coulomb_FACG( struct site  *__restrict lat , 
	      GLU_complex  *__restrict *__restrict out , 
	      GLU_complex  *__restrict *__restrict in , 
	      const void *__restrict forward , 
	      const void *__restrict backward , 
	      const GLU_real * __restrict p_sq ,
	      const double accuracy ,
	      const size_t max_iter )
{
  // allocations 
  GLU_complex **gt = NULL , **g = NULL , 
    **g_end = NULL , **g_up  = NULL ;
  GLU_complex **sn = NULL , **in_old = NULL ; // CG temporaries

  // allocate rotato
  struct sp_site_herm *rotato = NULL ; 

  // allocate traces
  allocate_traces( LCU ) ;

  // flag for if something goes wrong
  size_t tot_its = 0 ;

  // loop counter and timeslice number
  size_t i , t = 0 ;

  // temporary field allocations
  if( GLU_malloc( (void**)&gt     , ALIGNMENT , LCU * sizeof( GLU_complex* ) ) != 0 ||
      GLU_malloc( (void**)&g      , ALIGNMENT , LCU * sizeof( GLU_complex* ) ) != 0 ||
      GLU_malloc( (void**)&g_end  , ALIGNMENT , LCU * sizeof( GLU_complex* ) ) != 0 ||
      GLU_malloc( (void**)&g_up   , ALIGNMENT , LCU * sizeof( GLU_complex* ) ) != 0 ||
      GLU_malloc( (void**)&sn     , ALIGNMENT , TRUE_HERM * sizeof( GLU_complex* ) ) != 0 ||
      GLU_malloc( (void**)&in_old , ALIGNMENT , TRUE_HERM * sizeof( GLU_complex* ) ) != 0 ||
      GLU_malloc( (void**)&rotato , ALIGNMENT , LCU * sizeof( struct sp_site_herm ) ) != 0 ) {
    goto memfree ;
  }

#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < TRUE_HERM ; i ++  ) {
    GLU_malloc( (void**)&sn[i] , ALIGNMENT , LCU * sizeof( GLU_complex ) ) ;
    GLU_malloc( (void**)&in_old[i] , ALIGNMENT , LCU * sizeof( GLU_complex ) ) ; 
  }

  // allocate the transformation matrices
#pragma omp parallel for private(i) 
  PFOR( i = 0 ; i < LCU ; i++ ) {
    // allocations
    GLU_malloc( (void**)&g[i]     , ALIGNMENT , NCNC * sizeof( GLU_complex ) ) ; 
    GLU_malloc( (void**)&g_up[i]  , ALIGNMENT , NCNC * sizeof( GLU_complex ) ) ;
    GLU_malloc( (void**)&g_end[i] , ALIGNMENT , NCNC * sizeof( GLU_complex ) ) ;
    GLU_malloc( (void**)&gt[i]    , ALIGNMENT , NCNC * sizeof( GLU_complex ) ) ;
    identity( g[i] ) ;
    identity( g_end[i] ) ;
  }

  // OK so we have set up the gauge transformation matrices
  tot_its = steep_fix_FACG( lat , g_end , out , in , in_old , 
			    sn , rotato , gt , forward , backward , 
			    p_sq , t , accuracy , max_iter ) ; 
  
  // and t+1
  tot_its += steep_fix_FACG( lat , g , out , in , in_old , 
			     sn , rotato , gt , forward , backward , 
			     p_sq , t + 1 , accuracy , max_iter ) ; 

  // reunitarise the gauges to counteract the accumulated round-off error
  #pragma omp parallel for private(i) 
  PFOR( i = 0 ; i < LCU ; i++ ) { 
    gram_reunit( g_end[i] ) ; 
    gram_reunit( g[i] ) ; 
  }

  //gauge transform the links at x -> then set slice_gauge_up to be slice_gauge
  gtransform_slice( (const GLU_complex **)g_end , lat , 
		    (const GLU_complex **)g , t ) ; 

  //now we do the same for all time slices
  for( t = 2 ; t < Latt.dims[ ND - 1 ] ; t++ ) {

    // set the one above to be the identity so we can accumulate transforms
    #pragma omp parallel for private(i)
    PFOR( i = 0 ; i < LCU ; i++ ) { identity( g_up[i] ) ; }

    // gauge fix on this slice
    tot_its += steep_fix_FACG( lat , g_up , out , in , in_old , 
			       sn , rotato , gt , forward , backward , 
			       p_sq , t , accuracy , max_iter ) ; 

    // reunitarise to counteract the accumulated round-off error
    #pragma omp parallel for private(i) 
    PFOR( i = 0 ; i < LCU ; i++ ) { gram_reunit( g_up[i] ) ; }

    //gauge transform the links for this slice "g"
    gtransform_slice( (const GLU_complex **)g , lat , 
		      (const GLU_complex **)g_up , t - 1 ) ; 

    // and copy "g_up" (the working transformation matrices) into "g" 
    GLU_complex **ptr = g ;
    g = g_up ;
    g_up = ptr ;
  }

  // no need for a reunitarisation step as it has already been done
  // gauge transform the very final slice
  gtransform_slice( (const GLU_complex **)g , lat , 
		    (const GLU_complex **)g_end , t - 1 ) ; 

  // and free all of that memory, especially rotato
 memfree :

  if( sn != NULL ) {
    for( i = 0 ; i < TRUE_HERM ; i++ ) {
      free( sn[i]     ) ;
    }
  }
  if( in_old != NULL ) {
    for( i = 0 ; i < TRUE_HERM ; i++ ) {
      free( in_old[i] ) ;
    }
  }
  free( rotato ) ;
  free( in_old );
  free( sn     ) ;
  free_traces( ) ;

  // free temporary gauge transformation matrices
  if( gt != NULL ) {
#pragma omp parallel for private(i)
    PFOR( i = 0 ; i < LCU ; i ++  ) {
      free( gt[i]     ) ; 
    }
  }
  if( g != NULL ) {
#pragma omp parallel for private(i)
    PFOR( i = 0 ; i < LCU ; i ++  ) {
      free( g[i]     ) ; 
    }
  }
  if( g_up != NULL ) {
#pragma omp parallel for private(i)
    PFOR( i = 0 ; i < LCU ; i ++  ) {
      free( g_up[i]     ) ; 
    }
  }
  if( g_end != NULL ) {
#pragma omp parallel for private(i)
    PFOR( i = 0 ; i < LCU ; i ++  ) {
      free( g_end[i]     ) ; 
    }
  }
  free( gt    ) ;
  free( g     ) ; 
  free( g_up  ) ; 
  free( g_end ) ;
 
  // and return the total iterations
  return tot_its ;
}

// FASD coulomb gauge fix on slice t
static size_t
steep_fix_FA( const struct site *__restrict lat , 
	      GLU_complex *__restrict *__restrict slice_gauge , 
	      GLU_complex *__restrict *__restrict out , 
	      GLU_complex *__restrict *__restrict in , 
	      struct sp_site_herm *__restrict rotato ,
	      const void *__restrict forward , 
	      const void *__restrict backward , 
	      const GLU_real *__restrict psq , 
	      const size_t t ,
	      const double accuracy ,
	      const size_t max_iter )
{
  double tr = 1.0 ;
  size_t iters = 0 , temp_iters = 0 , control = 0 ;
  while ( ( tr > accuracy ) && ( iters < max_iter ) ) {
    //criteria:: We only randomly transform six times before complaining
    if( iters == ( max_iter - 1 ) ) {
      if( control < GF_GLU_FAILURES )  {
	random_gtrans_slice( slice_gauge ) ;
	iters = 0 ; 
	control++ ; 
	// continuation run
      } else if( tr < ( 1E3 * accuracy ) ) {
	iters = 0 ; 
	temp_iters += max_iter ;
      } else {
	temp_iters += max_iter ;
	break ;
      }
    }
    // perform a Fourier accelerated step
    steep_step_SD( lat , slice_gauge , out , in , rotato , NULL ,
		   forward , backward , psq , t , &tr ) ;

    // the log needs some help to stay on track
    #if ( defined deriv_full ) || ( defined deriv_fulln )
    if( (iters&127) == 0 ) {
      size_t i ;
      // reunitarise the computed gauge to counteract the accumulated round-off error
      #pragma omp parallel for private(i) 
      PFOR( i = 0 ; i < LCU ; i++ ) { gram_reunit( slice_gauge[i] ) ; }
    }
    #endif

    iters ++ ;
  }
  // and tell us about the slice's iteration
  get_info( t , tr , iters + temp_iters , control , accuracy , GLU_FALSE ) ; 
  return iters + temp_iters ;
}

// Coulomb gauge FASD driver
size_t
Coulomb_FASD( struct site  *__restrict lat , 
	      GLU_complex  *__restrict *__restrict out , 
	      GLU_complex  *__restrict *__restrict in , 
	      const void *__restrict forward , 
	      const void *__restrict backward , 
	      const GLU_real * __restrict p_sq ,
	      const double accuracy ,
	      const size_t max_iter )
{
  // allocations 
  GLU_complex **g = NULL , **g_end = NULL , **g_up  = NULL ;

  // allocate rotato
  struct sp_site_herm *rotato = NULL ;

  // flag for whether it succeeded or not
  size_t tot_its = 0 ;

  // allocate traces
  allocate_traces( LCU ) ;

  // initialise loop counter and timeslice index
  size_t i , t = 0 ;

  // allocate temporary gauge transformation matrices
  if( GLU_malloc( (void**)&g      , ALIGNMENT , LCU * sizeof( GLU_complex* ) ) != 0 ||
      GLU_malloc( (void**)&g_end  , ALIGNMENT , LCU * sizeof( GLU_complex* ) ) != 0 ||
      GLU_malloc( (void**)&g_up   , ALIGNMENT , LCU * sizeof( GLU_complex* ) ) != 0 ||
      GLU_malloc( (void**)&rotato , ALIGNMENT , LCU * sizeof( struct sp_site_herm ) ) != 0 ) {
    fprintf( stderr , "[GF] CFASD temporary gauge fields allocation failure\n" ) ;
    goto memfree ;
  }

  // and allocate the transformation matrices
#pragma omp parallel for private(i) 
  PFOR( i = 0 ; i < LCU ; i ++  ) {
    GLU_malloc( (void**)&g[i]     , ALIGNMENT , NCNC * sizeof( GLU_complex ) ) ;
    GLU_malloc( (void**)&g_up[i]  , ALIGNMENT , NCNC * sizeof( GLU_complex ) ) ;
    GLU_malloc( (void**)&g_end[i] , ALIGNMENT , NCNC * sizeof( GLU_complex ) ) ;
    identity( g[i] ) ;
    identity( g_end[i] ) ;
  }

  // OK so we have set up the gauge transformation matrices
  tot_its = steep_fix_FA( lat , g_end , out , in ,
			  rotato , forward , backward , p_sq , 
			  t , accuracy , max_iter ) ; 
  
  // and t+1
  tot_its += steep_fix_FA( lat , g , out , in , rotato ,
			   forward , backward , p_sq , 
			   t + 1 , accuracy , max_iter ) ; 

  // reunitarise the gauges to counteract the accumulated round-off error
  #pragma omp parallel for private(i) 
  PFOR( i = 0 ; i < LCU ; i++ ) { 
    gram_reunit( g_end[i] ) ; 
    gram_reunit( g[i] ) ; 
  }

  //gauge transform the links at x -> then set slice_gauge_up to be slice_gauge
  gtransform_slice( (const GLU_complex**)g_end , lat , 
		    (const GLU_complex**)g , t ) ; 

  //now we do the same for all time slices
  for( t = 2 ; t < Latt.dims[ ND - 1 ] ; t++ ) {

    // set the one above to be the identity so we can accumulate transforms
    #pragma omp parallel for private(i)
    PFOR( i = 0 ; i < LCU ; i++ ) { identity( g_up[i] ) ; }

    // gauge fix on this slice
    tot_its += steep_fix_FA( lat , g_up , out , in , 
			     rotato , forward , backward , p_sq , 
			     t , accuracy , max_iter ) ; 

    // reunitarise the working gauge trans matrices to counteract the accumulated round-off error
    #pragma omp parallel for private(i)
    PFOR( i = 0 ; i < LCU ; i++ ) { gram_reunit( g_up[i] ) ; }

    //gauge transform the links at x -> then set slice_gauge_up to be slice_gauge
    gtransform_slice( (const GLU_complex**)g , lat , 
		      (const GLU_complex**)g_up , t - 1 ) ; 

    // and copy "slice_up" (the working transformation matrices) into "slice" (the previous)
    GLU_complex **ptr = g ;
    g = g_up ;
    g_up = ptr ;
  }

  // gauge transform the very final slice, no need for reunit2
  gtransform_slice( (const GLU_complex**)g , lat , 
		    (const GLU_complex**)g_end , t - 1 ) ; 

 memfree :

  // free the temporary transformation matrix, that I called rotato
  free( rotato ) ;

  // and free the traces
  free_traces( ) ;

  // free temporary gauge transformation matrices
  if( g != NULL ) {
#pragma omp parallel for private(i)
    PFOR( i = 0 ; i < LCU ; i ++  ) {
      free( g[i]     ) ; 
    }
  }
  if( g_up != NULL ) {
#pragma omp parallel for private(i)
    PFOR( i = 0 ; i < LCU ; i ++  ) {
      free( g_up[i]     ) ; 
    }
  }
  if( g_end != NULL ) {
#pragma omp parallel for private(i)
    PFOR( i = 0 ; i < LCU ; i ++  ) {
      free( g_end[i]     ) ; 
    }
  }
  free( g     ) ; 
  free( g_up  ) ; 
  free( g_end ) ; 

  // and return the total iterations
  return tot_its ;
}

// tells us the probes we are using
void
query_probes_Coulomb( void ) {
  fprintf( stdout , "[GF] Using the following probes for the CG \n" ) ;
  size_t mu ;
  for( mu = 0 ; mu < LINE_NSTEPS ; mu++ ) {
    fprintf( stdout , "[GF] probe-%zu %f \n" , mu , alcg[mu] ) ; 
  }
  return ;
}

// keep this local
#undef LINE_NSTEPS
