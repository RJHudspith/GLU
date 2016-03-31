/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (CG.c) is part of GLU.

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
   @file CG.c
   @brief common code for the the CG routines
 */
#include "Mainfile.h"    // general definitions

#include "GLU_splines.h" // GLUbic spline interpolation code
#include "GLU_sums.h"    // round-off resistant summations
#include "gramschmidt.h" // for reunitarisation
#include "gtrans.h"      // gauge transformations

// this is the point where the gram-schmidt loses out
// to the n-ape det-rescaled projection
#if NC > 18
  #include "taylor_logs.h"
#endif

// some small memory for the stabler average
static double *traces = NULL ;

// gauge transformation and a log
#if (defined deriv_full) || (defined deriv_fulln)
static void
gtrans_log( GLU_complex A[ NCNC ] ,
	    const GLU_complex a[ NCNC ] ,
	    const GLU_complex b[ NCNC ] ,
	    const GLU_complex c[ NCNC ] )
{
  GLU_complex temp[ NCNC ] GLUalign ;
  equiv( temp , b ) ;
  gtransform_local( a , temp , c ) ;
  exact_log_slow( A , temp ) ;
  return ;
}
#endif

// could have several different searches here
double
approx_minimum( const size_t nmeas , 
		const double alphas[ nmeas ] ,
		const double functional[ nmeas ] )
{
  // compute the spline derivatives
  double derivative[ nmeas ] ;
  spline_derivative( derivative , alphas , functional , nmeas ) ;

  // compute the sum of the derivatives, if they are all small my thinking is
  // that we are at the limit of the precision of the functional
  register double sumder = 0.0 ;

  // best minimum index
  size_t bestmin = 0 ;
  // best alpha
  double func_min = 2.0 ;
  // number of derivatives with a minus sign
  size_t sumneg = 0 ;

  size_t i ;
  for( i = 0 ; i < nmeas ; i++ ) {

    // sum of the derivatives, if we are too flat
    // we exit returning the user-specified tuning alpha
    sumder += fabs( derivative[i] ) ;

    // we find a minimum using this dirty method
    if( derivative[i] < 0. ) {
      sumneg ++ ;
      // find the lowest minimum in case we have more than one
      // this is pretty unlikely unless we use a bunch of probes
      if( functional[i] < func_min ) {
	bestmin = i + 1 ;
	func_min = functional[i] ;
      }
    }

    #ifdef verbose
    fprintf( stdout , "[GF] der[%zu] %e \n" , i , derivative[i] ) ;
    #endif
  }

  #ifdef verbose
  fprintf( stdout , "[GF] sumneg  :: %zu \n" , sumneg ) ;
  fprintf( stdout , "[GF] bestmin :: %zu \n" , bestmin ) ;
  fprintf( stdout , "[GF] sumder  :: %e \n" , sumder ) ;
  #endif

  // if we are at the limit of precision we leave
  if( sumder < PREC_TOL ) {
    return Latt.gf_alpha ;
  }

  // at the moment this routine assumes a quadratic shape
  // if there are no negative terms in the derivative we return 0
  if( sumneg == 0 ) {
    return 0. ;
    // if it is all negative the best alpha is greater than our largest probe
    // we return the largest probe
  } else if( sumneg == nmeas ) {
    return alphas[nmeas-1] ;
    // otherwise we have bound the minimum and we solve for it 
  } else {
    const double result = cubic_min( alphas , functional , 
				     derivative , bestmin ) ;

    if( isnan( result ) ) { 
      return 0.0 ;
    } else {
      return result ;
    }
  }
  return 0.0 ;
}

// little wrapper for the traces array
void 
allocate_traces( const int LENGTH )
{
  traces = malloc( LENGTH * sizeof( double ) ) ;
  return ;
}

// Coulomb derivative term
double
coul_gtrans_fields( struct sp_site_herm *__restrict rotato ,
		    const struct site *__restrict lat ,
		    const GLU_complex *__restrict *__restrict slice_gauge ,
		    const size_t t ,
		    const double acc )
{
  const double fact = 1.0 / (double)( ( ND - 1 ) * NC ) ;
  size_t i ;
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LCU ; i ++ ) {
    GLU_complex temp[ NCNC ] GLUalign ;
    #ifdef deriv_lin
    double loc_sum = 0.0 ;
    #endif
    const size_t j = i + LCU * t ;
    size_t mu ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      const int it = lat[i].neighbor[mu] ;
      equiv( temp , lat[j].O[mu] ) ;
      gtransform_local( slice_gauge[i] , temp , slice_gauge[it] ) ;
      #if (defined deriv_lin) || (defined deriv_linn) 
      Hermitian_proj_short( rotato[i].O[mu] , temp ) ;
      #else
        #ifdef GLU_GFIX_SD
        if( acc < 0.1 ) {
          exact_log_slow_short( rotato[i].O[mu] , temp ) ;
        } else {
          Hermitian_proj_short( rotato[i].O[mu] , temp ) ;
        }
        #else
        // FACG: doing the above is detrimental?
        exact_log_slow_short( rotato[i].O[mu] , temp ) ;
        #endif
      #endif
      // compute val in here?
      #ifdef deriv_lin
      loc_sum += creal( trace( temp ) ) ;
      #endif
    }
    #ifdef deriv_lin
    traces[i] = (double)loc_sum * fact ;
    #endif
  }

  // this computes the \alpha == 0.0 contribution for the log def
#if defined deriv_full
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LCU ; i ++ ) {
    GLU_real tr ;
    register double loc_sum = 0.0 ;
    size_t mu ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      trace_ab_herm_short( &tr , rotato[i].O[mu] , rotato[i].O[mu] ) ;
      loc_sum += (double)tr ;
    }
    traces[i] = 0.5 * loc_sum * fact ;
  }
#endif

#if defined deriv_full
  return knuth_average( traces , LCU ) ;
#elif defined deriv_lin
  return 1.0 - knuth_average( traces , LCU ) ;
#endif
}

// for the polyak-ribiere
double
PRfmax( const double a , const double b )
{
  return ( b < a ? a : b ) ;
}

// gauge transformed functional evaluations
double
evaluate_alpha( const GLU_complex *__restrict *__restrict gauge , 
		const struct site *__restrict lat ,
		const int DIR ,
		const int LENGTH ,
		const int t ) 
{
  const double fact = 1.0 / (double)( NC * DIR ) ;

  // gauge transform to check the functional, this is the expensive bit, it is worth
  // really considering if there is something cheaper out there for us
  size_t i ;
  #pragma omp parallel for private(i)
  for( i = 0 ; i < LENGTH ; i++ ) {
    // some stacked allocations
    #if (defined deriv_full) || (defined deriv_fulln) || (defined deriv_fullnn)
    GLU_complex A[ NCNC ] GLUalign ;
    GLU_real trAA ;
    #endif
    // and compute the relevant slice
    const size_t j = LENGTH * t + i ;
    // gauge transform on site fields
    register double loc_sum = 0.0 ;
    size_t mu ;
    for( mu = 0 ; mu < DIR ; mu++ ) {
      // give it the right functional
      #if ( defined deriv_lin )
      // trace identity used here
      loc_sum += Re_trace_abc_dag_suNC( gauge[i] , lat[j].O[mu] , 
					gauge[lat[i].neighbor[mu]] ) ;
      #elif (defined deriv_linn )
      ///////////////////////////////////////////
      // TODO :: Think about the functionals being minimised
      #elif (defined deriv_fulln )
      ///////////////////////////////////////////
      // TODO :: Think about the functionals being minimised
      #else
      // gauge transform of U_\mu(x+\mu/2)
      gtrans_log( A , gauge[i] , lat[j].O[mu] , 
		  gauge[lat[i].neighbor[mu]] ) ;
      trace_ab_herm( &trAA , A , A ) ;
      loc_sum += 0.5 * (double)trAA ;
      #endif
    }
    // sum the trace or whatever
    traces[i] = (double)loc_sum * fact ;
  }

  // functional definitions are slightly different
#if ( defined deriv_lin ) || (defined deriv_linn )
  return 1.0 - knuth_average( traces , LENGTH ) ;
#else
  return knuth_average( traces , LENGTH ) ;
#endif
}

// is the same for Landau and Coulomb just with different LENGTHS
void
FOURIER_ACCELERATE( GLU_complex *__restrict *__restrict in ,
		    GLU_complex *__restrict *__restrict out ,
		    const void *__restrict forward ,
		    const void *__restrict backward ,
		    const GLU_real *__restrict psq ,
		    const size_t LENGTH ) 
{
#ifdef HAVE_FFTW3_H
  const fftw_plan *forw = ( const fftw_plan* )forward ;
  const fftw_plan *back = ( const fftw_plan* )backward ;
  ///// Fourier Acceleration //////
  #ifdef OMP_FFTW
  size_t mu , i ;
  for( mu = 0 ; mu < TRUE_HERM ; mu++ ) {
    PSPAWN fftw_execute( forw[mu] ) ; 
    #pragma omp parallel for private(i)
    PFOR( i = 0 ; i < LENGTH ; i++ ) {
      out[ mu ][ i ] *= psq[i] ;
    }
    PSPAWN fftw_execute( back[mu] ) ; 
  }
  PSYNC ;
  #else
  // single core FFT's
  size_t mu ;
  #pragma omp parallel for private(mu) schedule(dynamic)
  for( mu = 0 ; mu < TRUE_HERM ; mu++ ) {
    PSPAWN fftw_execute( forw[mu] ) ; 
    size_t i ;
    PFOR( i = 0 ; i < LENGTH ; i++ ) {
      out[ mu ][ i ] *= psq[i] ;
    }
    PSPAWN fftw_execute( back[mu] ) ; 
  }
  PSYNC ;
  #endif
#endif
  return ;
}

// and for freeing the traces array
void free_traces( void ) { free( traces ) ; }

// compute the functional quickly
double
gauge_functional_fast( const struct site *__restrict lat )
{
  const double fact = 1.0 / (double)( NC * ND ) ;
  size_t i ;
  #pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    #if (defined deriv_full) || (defined deriv_fulln) || (defined deriv_fullnn)
    GLU_complex A[ HERMSIZE ] ;
    GLU_real trAA ;
    #endif
    register double loc_sum = 0.0 ;
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      #if ( defined deriv_lin )
        #if NC == 3
        loc_sum += creal( lat[i].O[mu][0] ) ;
	loc_sum += creal( lat[i].O[mu][4] ) ;
	loc_sum += creal( lat[i].O[mu][8] ) ;
        #elif NC == 2
        loc_sum += creal( lat[i].O[mu][0] ) ;
	loc_sum += creal( lat[i].O[mu][3] ) ;
        #else
	loc_sum += creal( trace( lat[i].O[mu] ) ) ;
        #endif
      #elif (defined deriv_linn )
	GLU_complex temp[ NCNC ] GLUalign ;
	multab_suNC( temp , lat[i].O[mu] , lat[lat[i].neighbor[mu]].O[mu] ) ;
	loc_sum += creal( trace( temp ) ) ;
      #else
	// still need to think about these
      exact_log_slow_short( A , lat[i].O[mu] ) ;
      trace_ab_herm_short( &trAA , A , A ) ;
      loc_sum += 0.5 * (double)trAA ;
      #endif
    }
    traces[i] = loc_sum * fact ;
  }

#if ( defined deriv_lin ) || (defined deriv_linn )
  return 1.0 - knuth_average( traces , LVOLUME ) ;
#else
  return knuth_average( traces , LVOLUME ) ; 
#endif
}

// this is the same between Landau and Coulomb so I put it here 
void
set_gauge_matrix( GLU_complex *__restrict gx ,
		  const GLU_complex *__restrict *__restrict in ,
		  const double alpha ,
		  const size_t i )
{
#ifdef exp_exact
  GLU_complex short_gauge[ HERMSIZE ] ;
  #if NC==3
  short_gauge[ 0 ] = alpha * cimag( in[ 0 ][ i ] ) ;
  short_gauge[ 1 ] = -I * alpha * in[ 1 ][ i ] ;
  short_gauge[ 2 ] = -I * alpha * in[ 2 ][ i ] ;
  short_gauge[ 3 ] = -alpha * creal( in[ 0 ][ i ] ) ;
  short_gauge[ 4 ] = -I * alpha * in[ 3 ][ i ] ;
  #elif NC==2
  short_gauge[ 0 ] = alpha * cimag( in[ 0 ][ i ] ) ;
  short_gauge[ 1 ] = -I * alpha * in[ 1 ][ i ] ;
  #else 
  // this makes it hermitian, which is what the exponential expects
  int mu ;
  for( mu = 0 ; mu < HERMSIZE ; mu++ ) {
    short_gauge[mu] = -I * alpha * in[ mu ][ i ] ;
  }
  #endif
  exponentiate_short( gx , short_gauge ) ; 
#else
  #if NC==3
  gx[0] = 1.0 + I * alpha * cimag( in[0][i] ) ; 
  gx[1] = alpha * in[1][i] ;  
  gx[3] = -conj( gx[1] ) ; 
  gx[4] = 1.0 + -I * alpha * creal( in[0][i] ) ; 
  gx[6] = -alpha * conj( in[2][i] ) ;  
  gx[7] = -alpha * conj( in[3][i] ) ;
  #elif NC==2
  gx[0] = 1.0 + I * alpha * cimag( in[0][i] ) ; 
  gx[1] = alpha * in[1][i] ; 
  gx[2] = -conj( gx[1] ) ; 
  #else
  GLU_complex short_gauge[ HERMSIZE ] ;
  // this makes it antihermitian
  size_t mu ;
  for( mu = 0 ; mu < HERMSIZE ; mu++ ) {
    short_gauge[mu] = alpha * in[ mu ][ i ] ;
  }
  rebuild_antihermitian( gx , short_gauge ) ;
  add_constant( gx , 1.0 ) ; 
  #endif
  #if NC > 18
  nape_reunit( gx ) ; 
  #else
  gram_reunit( gx ) ; 
  #endif
#endif
}

// derivative
double
sum_deriv( const GLU_complex *__restrict *__restrict in , 
	   const size_t LENGTH )
{
  size_t i ;
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LENGTH ; i++ ) {
    register double loc_sum = 0.0 ;
    #if NC == 3
    loc_sum += creal(in[0][i]) * creal(in[0][i]) + cimag(in[0][i])*cimag(in[0][i]) 
      + creal(in[0][i]) * cimag( in[0][i] ) ;
    loc_sum += creal(in[1][i]) * creal(in[1][i]) + cimag(in[1][i])*cimag(in[1][i]) ; 
    loc_sum += creal(in[2][i]) * creal(in[2][i]) + cimag(in[2][i])*cimag(in[2][i]) ; 
    loc_sum += creal(in[3][i]) * creal(in[3][i]) + cimag(in[3][i])*cimag(in[3][i]) ; 
    #elif NC == 2
    loc_sum += creal(in[0][i]) * creal(in[0][i]) + cimag(in[0][i])*cimag(in[0][i]) ;
    loc_sum += creal(in[1][i]) * creal(in[1][i]) + cimag(in[1][i])*cimag(in[1][i]) ;
    #else
    GLU_complex temp[ HERMSIZE ] ;
    size_t mu ;
    for( mu = 0 ; mu < HERMSIZE ; mu++ ) {
      temp[mu] = in[mu][i] ;
    }
    GLU_real tr ;
    trace_ab_herm_short( &tr , temp , temp ) ;
    loc_sum = 0.5 * (double)tr ;
    #endif
    traces[i] = loc_sum ;
  }
  return 2.0 * kahan_summation( traces , LENGTH ) ; 
}

// Polak Ribiere Numerator
double
sum_PR_numerator( const GLU_complex *__restrict *__restrict in , 
		  const GLU_complex *__restrict *__restrict in_old ,
		  const size_t LENGTH ) 
{
  size_t i ;
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LENGTH ; i++ ) {
    register double loc_sum = 0.0 ;
#if NC == 3
    register GLU_complex temp = in[0][i] - in_old[0][i] ;
    loc_sum += 2.0 * ( creal( in[0][i] ) * creal( temp ) + cimag( in[0][i] ) * cimag( temp ) ) ;
    loc_sum += creal( in[0][i] ) * cimag( temp ) + cimag( in[0][i] ) * creal( temp ) ;
    temp = in[1][i] - in_old[1][i] ;
    loc_sum += 2.0 * ( creal( in[1][i] ) * creal( temp ) + cimag( in[1][i] ) * cimag( temp ) ) ;
    temp = in[2][i] - in_old[2][i] ;
    loc_sum += 2.0 * ( creal( in[2][i] ) * creal( temp ) + cimag( in[2][i] ) * cimag( temp ) ) ;
    temp = in[3][i] - in_old[3][i] ;
    loc_sum += 2.0 * ( creal( in[3][i] ) * creal( temp ) + cimag( in[3][i] ) * cimag( temp ) ) ;
#elif NC == 2
    register GLU_complex temp = in[0][i] - in_old[0][i] ;
    loc_sum += 2.0 * ( creal( in[0][i] ) * creal( temp ) + cimag( in[0][i] ) * cimag( temp ) ) ;
    temp = in[1][i] - in_old[1][i] ;
    loc_sum += 2.0 * ( creal( in[1][i] ) * creal( temp ) + cimag( in[1][i] ) * cimag( temp ) ) ;
#else
    GLU_complex temp[ HERMSIZE ] GLUalign , temp2[ NCNC ] GLUalign ;
    size_t mu ;
    for( mu = 0 ; mu < HERMSIZE ; mu++ ) {
      temp[mu] = in[mu][i] - in_old[mu][i] ;
      temp2[ mu ] = in[ mu ][ i ] ;
    }
    GLU_real tr ;
    trace_ab_herm_short( &tr , temp2 , temp ) ;
    loc_sum = (double)tr ;
#endif
    traces[i] = loc_sum ;
  }
  return kahan_summation( traces , LENGTH ) ;
}
