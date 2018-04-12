/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (line_search.c) is part of GLU.

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
   @file line_search.c
   @brief Line search routines
 */
#include "Mainfile.h"
#include "CG.h"          // set gauge matrix
#include "GLU_splines.h" // GLUbic spline interpolation code

// line search probes for the Coulomb CG-gauge fixing
#define PC1 (0.17)
#define PC2 (0.34)

#define PL1 (0.15)
#define PL2 (0.30)
static const double alphas[ 3 ] = { 0.0 , PL1 , PL2 } ;

void
exponentiate_gauge_CG( GLU_complex **gauge , 
		       const GLU_complex **in ,
		       const double alpha )
{
  size_t i ;
  // the derivative in this form is antihermitian i.e -> i.dA
#pragma omp for nowait private(i)
  PFOR( i = 0 ; i < LCU ; i++ ) {
    GLU_complex temp[ NCNC ] GLUalign , temp2[ NCNC ] GLUalign ;
    set_gauge_matrix( temp , in , alpha , i ) ;
    equiv( temp2 , gauge[i] ) ;
    multab_suNC( gauge[i] , temp , temp2 ) ;
  }
  return ;
}

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

    //if( derivative[i] > 0.0 ) break ;
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

// perform a line search using GLUbic splines for approximately the best alpha
void
line_search_Coulomb( double *red ,
		     GLU_complex **gauge ,
		     const struct s_site *rotato ,
		     const struct draughtboard db ,
		     const struct site *lat ,
		     const GLU_complex **in ,
		     const size_t t )
{
  size_t i , j ;
  #pragma omp for private(i)
  for( i = 0 ; i < db.Nsquare[0] ; i++ ) {
    GLU_complex A1[ NCNC ] GLUalign ;
    GLU_complex A2[ NCNC ] GLUalign ;
    GLU_complex B[ NCNC ] GLUalign ;
    double loc_v[ LINE_NSTEPS ] ;
    size_t mu , n ;
    for( mu = 0 ; mu < LINE_NSTEPS ; mu++ ) {
      loc_v[ mu ] = 0.0 ;
    }
    const size_t idx = db.square[0][i] ;
    size_t bck ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      bck = lat[idx].back[mu] ;
      loc_v[0] += creal( trace( rotato[idx].O[mu] ) ) ;
      loc_v[0] += creal( trace( rotato[bck].O[mu] ) ) ;
    }
    // evaluate for probes -> We are still calling set_gauge_matrix
    // wayyy too often!
    //    for( n = 1 ; n < LINE_NSTEPS ; n++ ) {
    set_gauge_matrix( A1 , in , PC1 , idx ) ;
    set_gauge_matrix( A2 , in , PC2 , idx ) ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      bck = lat[idx].back[mu] ;
      // positive ones
      set_gauge_matrix( B , in , PC1 , lat[idx].neighbor[mu] ) ;
      loc_v[1] += Re_trace_abc_dag_suNC( A1 , rotato[idx].O[mu] , B ) ;
      // negative ones
      set_gauge_matrix( B , in , PC1 , bck ) ;
      loc_v[1] += Re_trace_abc_dag_suNC( B , rotato[bck].O[mu] , A1 ) ;
      // positive ones
      set_gauge_matrix( B , in , PC2 , lat[idx].neighbor[mu] ) ;
      loc_v[2] += Re_trace_abc_dag_suNC( A2 , rotato[idx].O[mu] , B ) ;
      // negative ones
      set_gauge_matrix( B , in , PC2 , bck ) ;
      loc_v[2] += Re_trace_abc_dag_suNC( B , rotato[bck].O[mu] , A2 ) ;
    }
    const size_t th = get_GLU_thread() ;
    // reductions
    for( n = 0 ; n < LINE_NSTEPS ; n++ ) {
      red[ n + th * CLINE ] += loc_v[ n ] ;
    }
  }
      
  // compute the vals
  double val[ LINE_NSTEPS ] ;
  for( j = 0 ; j < LINE_NSTEPS ; j++ ) {
    val[j] = 0.0 ;
    for( i = 0 ; i < Latt.Nthreads ; i++ ) {
      val[ j ] += red[ j + CLINE*i ] ;
    }
    val[ j ] = 1 - val[j] / ( (ND-1)*NC*LCU ) ;
  }

  const double Calcg[ LINE_NSTEPS ] = { 0.0 , PC1 , PC2 } ;
  const double min = approx_minimum( LINE_NSTEPS , Calcg , val ) ;
 
  exponentiate_gauge_CG( gauge , in , min ) ;

  return ;
}

void
egauge_Landau( GLU_complex **gauge , 
	       const GLU_complex **in ,
	       const double alpha )
{
  size_t i ;
  // the derivative in this form is antihermitian i.e -> i.dA
#pragma omp for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    set_gauge_matrix( gauge[i] , in , alpha , i ) ;
  }
  return ;
}

void
line_search_Landau( double *red ,
		    GLU_complex **gauge , 
		    const struct site *lat ,
		    const GLU_complex **in )
{
  size_t k = 0 , i ; 
#pragma omp for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    register double loc_sum = 0.0 ;
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      #if NC == 3
      loc_sum += creal( lat[i].O[mu][0] )
	+ creal( lat[i].O[mu][4] )
	+ creal( lat[i].O[mu][8] ) ;
      #else
      loc_sum += creal( trace( lat[i].O[mu] ) ) ;
      #endif
    }
    const size_t th = get_GLU_thread() ;
    // reductions
    red[ th * CLINE ] += loc_sum ;
  }
  // loop probes
  for( k = 1 ; k < LINE_NSTEPS ; k++ ) {

#pragma omp for private(i)
    for( i = 0 ; i < LVOLUME ; i++ ) {
      set_gauge_matrix( gauge[i] , in , alphas[ k ] , i ) ;
    }
    
#pragma omp for private(i)
    for( i = 0 ; i < LVOLUME ; i++ ) {
      register double loc_sum = 0.0 ;
      size_t mu ;
      for( mu = 0 ; mu < ND ; mu++ ) {
	loc_sum += Re_trace_abc_dag_suNC( gauge[i] , lat[i].O[mu] , 
					  gauge[lat[i].neighbor[mu]] ) ;
      }
      
      const size_t th = get_GLU_thread() ;
      // reductions
      red[ k + th * CLINE ] += loc_sum ;
    }
  }
  
  // do the reduction
  double val[ LINE_NSTEPS ] ;
  for( k = 0 ; k < LINE_NSTEPS ; k++ ) {
    val[ k ] = 0.0 ;
    for( i = 0 ; i < Latt.Nthreads ; i++ ) {
      val[ k ] += red[ k + CLINE*i ] ;
    }
    val[ k ] = 1.0 - val[ k ] / ( ND*NC*LVOLUME ) ;
  }
  
  const double min = approx_minimum( LINE_NSTEPS , alphas , val ) ;
  
#pragma omp for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    set_gauge_matrix( gauge[i] , in , min , i ) ;
  }
  
  return ; 
}
