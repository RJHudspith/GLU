/*
v    Copyright 2013 Renwick James Hudspith

    This file (gftests.c) is part of GLU.

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
   @file gftests.c
   @brief has all of the gauge tests i.e. Plaquettes , Links , Gauge fixing tests
 */

#include "Mainfile.h"

#include "clover.h"     // field strength tensor
#include "lin_derivs.h" // Hermitian approximation derivatives
#include "log_derivs.h" // Logarithmic field definition derivatives

// the one with the openmp locks, just an excuse to try this out to be honest
#if ( defined HAVE_OMP_H ) && ( defined _OPENMP )
   #define GLU_OMP_MEAS
   #include <omp.h>
#endif

// look at constant-ness of our gauge fields in all directions
void
const_time( const struct site *__restrict lat ,
	    double *const_lin ,
	    double *const_log )
{
  GLU_complex sum1[ NCNC ] = { } ;
  GLU_complex sum1_log[ NCNC ] = { } ;

  // initialise the first layer
  int j ;
  for( j = 0 ; j < LCU ; j++ ) {
    GLU_complex A[ NCNC ] ;
    Hermitian_proj( A , lat[j].O[ ND - 1 ] ) ; 
    a_plus_b( sum1 , A ) ; 
      
    // log defs ..
    exact_log_slow( A , lat[j].O[ ND - 1 ] ) ; 
    a_plus_b( sum1_log , A ) ; 
  }

  //loop time slices
  double log_average = 0. , average = 0. ;
  int t ;
#pragma omp parallel for private(t) reduction(+:log_average) reduction(+:average)
  for( t = 1 ; t < Latt.dims[ ND - 1 ] ; t++ ) {  

      double av = 0. , avlog = 0. ; 

      GLU_complex sum2[ NCNC ] = { } ;
      GLU_complex sum2_log[ NCNC ] = { } ;
      //loop inner cube
      int i ;
      for( i = LCU * t ; i < LCU * t + LCU ; i++ ) {
	GLU_complex A[ NCNC ] ;
	Hermitian_proj( A , lat[i].O[ ND - 1 ] ) ; 
	a_plus_b( sum2 , A ) ; 
	
	exact_log_slow( A , lat[i].O[ ND - 1 ] ) ; 
	a_plus_b( sum2_log , A ) ; 
      }
      //find the difference from the first slice

      #if NC == 3

      register const double re0 = creal( sum1[0] ) - creal( sum2[0] ) ;
      av += (double)fabs( re0 ) ; 
      av += (double)cabs( sum1[1] - sum2[1] ) ;
      av += (double)cabs( sum1[2] - sum2[2] ) ;
      register const double re4 = creal( sum1[4] ) - creal( sum2[4] ) ;
      av += (double)fabs( re4 ); 
      av += (double)cabs( sum1[5] - sum2[5] ) ;
      // log defs
      register const double rel0 = creal( sum1_log[0] ) - creal( sum2_log[0] ) ;
      avlog += (double)fabs( rel0 ) ; 
      avlog += (double)cabs( sum1_log[1] - sum2_log[1] ) ;
      avlog += (double)cabs( sum1_log[2] - sum2_log[2] ) ;
      register const double rel4 = creal( sum1_log[4] ) - creal( sum2_log[4] ) ;
      avlog += (double)fabs( rel4 ); 
      avlog += (double)cabs( sum1_log[5] - sum2_log[5] ) ;

      #elif NC == 2

      register const double re0 = creal( sum1[0] ) - creal( sum1[0] ) ;
      av += (double)fabs( re0 ) ; 
      av += (double)cabs( sum1[1] - sum2[1] ) ;
      // log def
      register const double rel0 = creal( sum1[0] ) - creal( sum1[0] ) ;
      avlog += (double)fabs( rel0 ) ; 
      avlog += (double)cabs( sum1_log[1] - sum2_log[1] ) ;

      #else
  
      GLU_complex dif[ NCNC ] ; 
      GLU_complex dif_log[ NCNC ] ; 
      b_min_c( dif , sum1 , sum2 ) ; 
      b_min_c( dif_log , sum1_log , sum2_log ) ; 
      for( i = 0 ; i < NCNC ; i++ ) {
	av += (double)cabs( dif[i] ) ; 
	avlog += (double)cabs( dif_log[i] ) ;
      }

      #endif

      log_average = log_average + (double)avlog ;
      average = average + (double)av ;
  }

  // Write the average difference in time-constantness
  const double check = 1.0 / ( NCNC * ( Latt.dims[ ND - 1 ] - 1 ) ) ;

  *const_lin = average * check ;
  *const_log = log_average * check ;

  return ;
}

// gauge fixing functional calculator, dependent on the derivative type...
double
gauge_functional( const struct site *__restrict lat )
{
  // for the log defs the functional is actually Re(Tr(A*A)) so we do this
  int i ;
  double trAA = 0.0 ;
  #pragma omp parallel for private(i) reduction(+:trAA)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    // removes unused variable warning ...
    #if ( defined deriv_full ) || ( defined deriv_fulln ) || ( defined deriv_fullnn ) 
    GLU_complex A[ NCNC ] ;
    #endif
    GLU_real tr ;
    register double loc_tr = 0.0 ;
    int mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      #if ( defined deriv_linn ) || ( defined deriv_lin ) 
      tr = (GLU_real)creal( trace( lat[i].O[mu] ) ) ;
      #else
      exact_log_slow( A , lat[i].O[mu] ) ;
      trace_ab_herm( &tr , A , A ) ;
      // factor of 1/2 for consistency between the two defs
      tr *= 0.5 ;
      #endif
      loc_tr += (double)tr ;
    }
    trAA = trAA + (double)loc_tr ;
  }
  return 1.0 - trAA / (double)( ND * NC * LVOLUME ) ;
}

//test the gauge transform matrices should tend to one so 1-Re( tr( G( x ) ) )->0
double
gauge_test( const GLU_complex *__restrict *__restrict gauge ) 
{
  double sum = 0.0 ; 
  int i ; 
#pragma omp parallel for private(i) reduction(+:sum) 
  for( i = 0 ; i < LVOLUME ; i++ ) {
    double tr ; 
    speed_trace_Re( &tr , gauge[i] ) ;
    sum = sum + (double)tr ; 
  }
  return 1.0 - sum /(double)( NC * LVOLUME ) ; 
}
  

// gets the functional for a gauge transformation slice
double
gtrans_functional( const struct site *__restrict lat ,
		   const GLU_complex *__restrict *__restrict slice_gauge ,
		   const int t )
{
  const int slice_idx = LCU * t ;
  int i ;
  double tr = 0. ;
#if ( defined deriv_lin ) || ( defined deriv_linn )
  #pragma omp parallel for private(i) reduction(+:tr)
  for( i = 0 ; i < LCU ; i++ ) {
    const int j = slice_idx + i ;
    int mu ;
    double loc_tr = 0. ;
    GLU_complex trabc ;
    for( mu = 0 ; mu < ND - 1 ; mu ++  ){
      trace_abc_dag( &trabc , slice_gauge[i] , lat[j].O[mu] ,
		     slice_gauge[lat[i].neighbor[mu]] ) ;
      loc_tr += (double)creal( trabc ) ;
    }
    tr = tr + (double)loc_tr ;
  }
  return 1.0 - ( tr / (double)( ( ND - 1 ) * LCU * NC ) ) ;
#else
  // actually have to do the multiplications, yuck
  #pragma omp parallel for private(i) reduction(+:tr)
  for( i = 0 ; i < LCU ; i++ ) {
    GLU_complex A[ NCNC ] , temp[ NCNC ] , temp2[ NCNC ] ;
    GLU_real trAA ;
    double loc_tr = 0. ;
    const int j = slice_idx + i ;
    int mu ;
    for( mu = 0 ; mu < ND - 1 ; mu ++  ){
      const int it = lat[i].neighbor[mu] ;
      // gauge transform
      multab_dag_suNC( temp2 , lat[j].O[mu] , slice_gauge[it] ) ; 
      multab_suNC( temp , slice_gauge[i] , temp2 ) ; 
      exact_log_slow( A , temp ) ;
      trace_ab_herm( &trAA , A , A ) ;
      loc_tr += (double)trAA ;
    }
    tr = tr + (double)loc_tr ;
  }
  return tr / (double)( 2.0 * ( ND - 1 ) * LCU * NC ) ;
#endif
}

//test3 the standard theta test tr[Delta.Delta^{dagger}]
double
theta_test_lin( const struct site *__restrict lat , 
		GLU_real *max ,
		const int MAX_DIR ) 
{
  double tr = 0. ; 
  int i ; 
  *max = 0. ;
  #ifdef GLU_OMP_MEAS
  omp_lock_t writelock ;
  omp_init_lock( &writelock ) ;
  #endif
#pragma omp parallel for private(i) reduction(+:tr)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    GLU_complex temp[ HERMSIZE ] = { } ; 
    register double res = latt_deriv_AntiHermitian_proj( temp , lat , 
							 i , MAX_DIR ) ; 
    tr = tr + (double)res ; 
    #ifdef GLU_OMP_MEAS
    if( unlikely( res > *max ) ) { 
      omp_set_lock( &writelock ) ;
      *max = res ;
      omp_unset_lock( &writelock ) ;
    } 
    #else
    if( unlikely( res > *max ) ) { *max = res ; } 
    #endif
  }
  #ifdef GLU_OMP_MEAS
  omp_destroy_lock( &writelock ) ;
  #endif
  return 2. * tr / (double)( NC * LVOLUME ) ;
}

//test2 my test using log-def tr[Delta.Delta^{dagger}]
double
theta_test_log( const struct site *__restrict lat , 
		GLU_real *max ,
		const int MAX_DIR )
{
  double tr = 0. ; 
  int i ; 
  *max = 0. ; 
  #ifdef GLU_OMP_MEAS
  omp_lock_t writelock ;
  omp_init_lock( &writelock ) ;
  #endif
  #pragma omp parallel for private(i) reduction(+:tr) 
  for( i = 0 ; i < LVOLUME ; i++ ) {
    GLU_complex temp[ HERMSIZE ] = { } ; 
    double functional ;
    register const double res = log_deriv( temp , &functional , lat ,
					   i , MAX_DIR ) ; 
    tr = tr + (double)res ; 
    #ifdef GLU_OMP_MEAS
    if( unlikely( res > *max ) ) { 
      omp_set_lock( &writelock ) ;
      *max = res ;
      omp_unset_lock( &writelock ) ;
    } 
    #else
    if( unlikely( res > *max ) ) { *max = res ; } 
    #endif
  }
  #ifdef GLU_OMP_MEAS
  omp_destroy_lock( &writelock ) ;
  #endif
  return 2. * tr / ( NC * LVOLUME ) ; 
}

// and clean it up
#ifdef GLU_OMP_MEAS
  #undef GLU_OMP_MEAS
#endif
