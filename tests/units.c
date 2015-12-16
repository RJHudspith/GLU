/*
    Copyright 2013 Renwick James Hudspith

    This file (units.c) is part of GLU.

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
   @file units.c
   @brief a selection of unit tests for the code using minunit
 */

#include "Mainfile.h"

#include "exactQ.h"
#include "expMat.h"
#include "geometry.h"
#include "GLU_rng.h"
#include "GLUlib_wrap.h"
#include "gramschmidt.h"
#include "invert.h"
#include "LU.h"
#include "minunit.h"
#include "plaqs_links.h"
#include "POLY.h"
#include "random_config.h"

// these are set to their lexicographical order
static GLU_complex a[ NCNC ] , b[ NCNC ] , c[ NCNC ] ;

struct latt_info Latt ;
GLU_bool INIT_RNG ;
int tests_run = 0 ;
int tests_fail = 0 ;

// test the add constant routine
static char *add_constant_test( ) {
  GLU_bool is_ok = GLU_TRUE ;
  GLU_complex A[ NCNC ] ;
  equiv( A , a ) ;
  add_constant( A , 1.0 ) ;
  if( cabs( trace( A ) - trace( a ) - NC ) > PREC_TOL ) is_ok = GLU_FALSE ;
  mu_assert("[GLUnit] error : add constant broken", is_ok );
  return 0;
}

// test the a+b routine
static char *a_plus_b_test( ) {
  GLU_bool is_ok = GLU_TRUE ;
  GLU_complex d[ NCNC ] ;
  equiv( d , b ) ;
  a_plus_b( d , b ) ;
  int i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    if( cabs( d[i] - 2*i*(1.+I) ) > PREC_TOL ) is_ok = GLU_FALSE ;
  }  
  mu_assert("[GLUnit] error : add_plus_b broken", is_ok );
  return 0;
}

// test the add constant routine
static char *a_plus_CSxb_test( ) {
  GLU_bool is_ok = GLU_TRUE ;
  GLU_complex d[ NCNC ] ;
  equiv( d , a ) ;
  a_plus_CSxb( d , b , I ) ;
  int i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    if( fabs( creal( d[i] ) ) > PREC_TOL ) is_ok = GLU_FALSE ;
    if( fabs( cimag( d[i] ) - 2*cimag( b[i] ) ) > PREC_TOL ) is_ok = GLU_FALSE ;
  }  
  mu_assert("[GLUnit] error : a_plus_CSxb broken", is_ok );
  return 0;
}

// test the add constant routine
static char *a_plus_Sxb_test( ) {
  GLU_bool is_ok = GLU_TRUE ;
  GLU_complex d[ NCNC ] ;
  equiv( d , a ) ;
  a_plus_Sxb( d , b , -1.0 ) ;
  int i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    if( cabs( d[i] ) > PREC_TOL ) is_ok = GLU_FALSE ;
  }  
  mu_assert("[GLUnit] error : a_plus_Sxb broken", is_ok );
  return 0;
}

// test the add constant routine
static char *a_plus_Sxbminc_short_test( ) {
  GLU_complex hHa[ HERMSIZE ] , Ua[ NCNC ] ;

  generate_NCxNC( Ua ) ;
  reunit2( Ua ) ;

  Hermitian_proj_short( hHa , Ua ) ;

  // a += S * ( b - c ) 
  GLU_complex res[ HERMSIZE ] ;
  int i ;
  GLU_bool is_ok = GLU_TRUE ;
  for( i = 0 ; i < HERMSIZE ; i++ ) {
    res[ i ] = hHa[ i ] ;
  }
  a_plus_Sxbminc_short( hHa , -1.0 , hHa , hHa ) ;
  for( i = 0 ; i < HERMSIZE ; i++ ) {
    if( cabs( hHa[i] - res[i] ) > PREC_TOL ) is_ok = GLU_FALSE ;
  }
  mu_assert("[GLUnit] error : a_plus_Sxbminc broken", is_ok );
  return 0 ;
}

// test the a=b-c routine
static char *b_min_c_test( ) {
  GLU_bool is_ok = GLU_TRUE ;
  GLU_complex d[ NCNC ] ;
  b_min_c( d , b , c ) ;
  int i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    if( cabs( d[i] ) > PREC_TOL ) is_ok = GLU_FALSE ;
  } 
  mu_assert("[GLUnit] error : b_min_c broken", is_ok );
  return 0;
}

// test the cofactor transpose? -> probably dated
static char *cofactor_transpose_test( ) {
  // cofactor_transpose / determinant is the inverse!
  GLU_complex Ua[ NCNC ] ;

  // create random SU(N) matrix
  generate_NCxNC( Ua ) ;
  reunit2( Ua ) ;

  GLU_complex ctUa[ NCNC ] ;
  GLU_complex dt = cofactor_transpose( ctUa , Ua ) ;
  GLU_bool is_ok = GLU_TRUE ;
  if( cabs( dt - det( Ua ) ) > PREC_TOL ) is_ok = GLU_FALSE ;

  GLU_complex invUa[ NCNC ] ;
  inverse( invUa , Ua ) ;

  int i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    if( cabs( invUa[ i ] - ctUa[ i ] ) > PREC_TOL ) is_ok = GLU_FALSE ;
  }
  mu_assert("[GLUnit] error : cofactor_transpose broken", is_ok );
  return 0;
}

// test the equivalence of matrices
static char *equiv_test( ) {
  GLU_complex d[ NCNC ] ;
  identity( d ) ;
  equiv( d , a ) ;
  double sum = 0.0 ;
  int i ;
  for( i = 0 ; i < NCNC ; i++ ) sum += cabs( d[i] - a[i] ) ; 
  mu_assert("[GLUnit] error : equiv broken", !( sum > PREC_TOL ) );
  return 0;
}

// test the daggering
static char *dagger_test( ) {
  GLU_bool is_ok = GLU_TRUE ;
  GLU_complex d[ NCNC ] ;
  dagger( d , a ) ;
  int i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      if( cabs( conj( d[ j + NC * i ] ) - a[ i + NC * j ] ) > PREC_TOL ) 
	is_ok = GLU_FALSE ; 
    }
  }
  mu_assert("[GLUnit] error : dagger broken", is_ok );
  return 0;
}

// test the diagonal routine
static char *diag_test( ) {
  GLU_bool is_ok = GLU_TRUE ;
  GLU_complex d[ NCNC ] , sum = 0.0 ;
  diag( d , 0.5 ) ;
  int i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    sum += creal( d[ i ] ) ; 
  }
  if( cabs( sum - NC * 0.5 ) > PREC_TOL ) is_ok = GLU_FALSE ; 
  mu_assert("[GLUnit] error : diag broken", is_ok );
  return 0;
}

// test the add constant routine
static char *diag_vect_test( ) {
  GLU_bool is_ok = GLU_TRUE ;
  GLU_complex d[ NCNC ] , v[ NC ] , sum = 0.0 ;
  int i ;
  for( i = 0 ; i < NC ; i++ ) {
    v[ i ] = i ;
  }
  diag_vect( d , v ) ;
  for( i = 0 ; i < NCNC ; i++ ) {
    sum += d[ i ] ; 
  }
  if( cabs( sum - NC*(NC-1)/2 ) > PREC_TOL ) is_ok = GLU_FALSE ; 
  mu_assert("[GLUnit] error : diag_vect broken", is_ok );
  return 0;
}

// test the identity routine
static char *identity_test(  ) {
  GLU_bool is_ok = GLU_TRUE ;
  GLU_complex d[ NCNC ] ;
  identity( d ) ;
  // diagonal is 1, everything else is 0
  int i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      if( i != j ) {
	if( fabs( d[j+i*NC] ) > PREC_TOL ) {
	  is_ok = GLU_FALSE ;
	  break ;
	} 
      } else {
	if( fabs( d[j+i*NC] - 1.0 ) > PREC_TOL ) {
	  is_ok = GLU_FALSE ;
	  break ;
	} 
      }
    }
  }
  mu_assert("[GLUnit] error : identity is broken", is_ok );
  return 0;
}

// test the is unitary test
static char *is_unitary_test( ) {
  GLU_complex d[ NCNC ] ;
  identity( d ) ;
  mu_assert("[GLUnit] error : is_unitary broken", !is_unitary( d ) );
  return 0;
}

#if NC > 3
// test the LU_det comparing to the identity for now
static char *LU_det_test( ) {
  GLU_complex d[ NCNC ] ;
  identity( d ) ;
  mu_assert("[GLUnit] error : LU_det broken", 
	    cabs( LU_det( NC , d ) - 1.0 ) > PREC_TOL ? GLU_FALSE : GLU_TRUE  );
  return 0;
}
#endif

// test the determinant of the identity
static char *det_test( ) {
  GLU_complex d[ NCNC ] ;
  identity( d ) ;
  mu_assert("[GLUnit] error : det broken", 
	    cabs( det( d ) - 1.0 ) > PREC_TOL ? GLU_FALSE : GLU_TRUE  );
  return 0;
}

// test matrix multiply a vector
static char *mat_mult_vec_test(  ) {
  GLU_complex v[ NC ] , rv[ NC ] ;
  int i ;
  for( i = 0 ; i < NC ; i++ ) {
    v[ i ] = i ;
  }
  mat_mult_vec( rv , b , v ) ;
  mu_assert("[GLUnit] error : mat_mult_vec is broken",  1 );
  return 0 ;
}

// test the identity routine
static char *M_times_c_test(  ) {
  GLU_complex d[ NCNC ] ;
  int i ;
  GLU_bool is_ok = GLU_TRUE ;
  equiv( d , b ) ;
  M_times_c( d , 2.0 ) ;
  for( i = 0 ; i < NCNC ; i++ ) {
    if( cabs( d[i] - 2. * b[i] ) > PREC_TOL ) is_ok = GLU_FALSE ;
  }
  mu_assert("[GLUnit] error : M_times_c is broken", is_ok );
  return 0 ;
}

// matrix power
static char *matrix_power_test(  ) {
  GLU_complex Ua[ NCNC ] , Ub[ NCNC ] ;

  // generate a random matrix
  generate_NCxNC( Ua ) ;
  reunit2( Ua ) ;

  matrix_power( Ub , Ua , 11 ) ;

  GLU_complex slow[ NCNC ] ;
  equiv( slow , Ua ) ;
  int i ;
  for( i = 0 ; i < 10 ; i++ ) {
    multab_atomic_right( slow , Ua ) ;
  }

  GLU_bool is_ok = GLU_TRUE ;
  for( i = 0 ; i < NCNC ; i++ ) {
    if( cabs( slow[ i ] - Ub[ i ] ) > PREC_TOL ) is_ok = GLU_FALSE ; 
  }
  mu_assert("[GLUnit] error : matrix_power is broken", is_ok );
  return 0 ;
}

// test the outerproduct should all be ones!
static char *outerproduct_test(  ) {
  GLU_complex o[ NC ] , v[ NC ] ;
  int i ;
  GLU_bool is_ok = GLU_TRUE ;
  for( i = 0 ; i < NC ; i++ ) {
    o[ i ] = v[ i ] = 1.0 ;
  }
  GLU_complex d[ NCNC ] ;
  outerproduct( d , o , v ) ;
  for( i = 0 ; i < NCNC ; i++ ) {
    if( cabs( d[i] - 1.0 ) > PREC_TOL ) is_ok = GLU_FALSE ;
  }
  mu_assert("[GLUnit] error : outerproduct is broken", is_ok );
  return 0 ;
}

// test the hermitian packing
static char *pack_hermitian_test(  ) {
  GLU_complex Ha[ NCNC ] , sHa[ HERMSIZE ] ;
  Hermitian_proj( Ha , a ) ;
  pack_hermitian( sHa , Ha ) ;
  int i , j , idx = 0 ;
  GLU_bool is_ok = GLU_TRUE ;
  for( i = 0 ; i < NC-1 ; i++ ) {
    for( j = i ; j < NC ; j++ ) {
      if( cabs( Ha[j+i*NC] - sHa[idx] ) > PREC_TOL ) is_ok = GLU_FALSE ;
      idx++ ;
    }
  }
  mu_assert("[GLUnit] error : pack_hermitian is broken", is_ok );
  return 0 ;
}

// test the rebuilding antihermitian
static char *rebuild_antihermitian_test(  ) {
  GLU_complex Ha[ NCNC ] , sHa[ HERMSIZE ] , rHa[ NCNC ];
  Hermitian_proj( Ha , a ) ;
  pack_hermitian( sHa , Ha ) ;
  int i ;
  for( i = 0 ; i < HERMSIZE ; i++ ) {
    sHa[ i ] *= I ;
  }
  rebuild_antihermitian( rHa , sHa ) ;
  GLU_bool is_ok = GLU_TRUE ;
  for( i = 0 ; i < NCNC ; i++ ) {
    if( cabs( rHa[ i ] - I * Ha[ i ] ) > PREC_TOL ) is_ok = GLU_FALSE ;
  }
  mu_assert("[GLUnit] error : rebuild_hermitian is broken", is_ok );
  return 0 ;
}

// test the rebuilding hermitian
static char *rebuild_hermitian_test(  ) {
  GLU_complex Ha[ NCNC ] , sHa[ HERMSIZE ] , rHa[ NCNC ];
  Hermitian_proj( Ha , a ) ;
  pack_hermitian( sHa , Ha ) ;
  rebuild_hermitian( rHa , sHa ) ;
  int i ;
  GLU_bool is_ok = GLU_TRUE ;
  for( i = 0 ; i < NCNC ; i++ ) {
    if( cabs( rHa[ i ] - Ha[ i ] ) > PREC_TOL ) is_ok = GLU_FALSE ;
  }
  mu_assert("[GLUnit] error : rebuild_hermitian is broken", is_ok );
  return 0 ;
}

// test the shortening and rebuilding
static char *shorten_test(  ) {
  GLU_bool is_ok = GLU_TRUE ;
#if NC < 4
  GLU_real sA[ NCNC-1 ] ;
  GLU_complex Ua[ NCNC ] , rUa[ NCNC ] ;

  // random SU(N) matrix
  generate_NCxNC( Ua ) ;
  reunit2( Ua ) ;
  
  // test if shorten -> rebuild gives same answer
  shorten( sA , Ua ) ;
  rebuild( rUa , sA ) ;

  int i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    if( cabs( rUa[ i ] - Ua[ i ] ) > PREC_TOL ) is_ok = GLU_FALSE ;
  }
#endif
  mu_assert("[GLUnit] error : shorten is broken",  is_ok );
  return 0 ;
}

// test the identity routine
static char *speed_det_test(  ) {
  GLU_complex d[ NCNC ] , s ;
  identity( d ) ;
  speed_det( &s , d ) ;
  mu_assert("[GLUnit] error : det broken", 
	    cabs( s - 1.0 ) > PREC_TOL ? GLU_FALSE : GLU_TRUE  );
  return 0 ;
}

// test the identity routine
static char *speed_trace_test(  ) {
  GLU_complex tr ;
  int sum = 0 , i ;  
  speed_trace( &tr , b ) ;
  for( i = 0 ; i < NC ; i++ ) sum += i ;
  mu_assert("[GLUnit] error : trace is broken", !( creal( (int)tr - sum * ( NC+1 )) > PREC_TOL ) );

  double rtr ;
  speed_trace_Re( &rtr , b ) ;
  mu_assert("[GLUnit] error : speed_trace_Re is broken", !( creal( (int)tr - sum * ( NC+1 )) > PREC_TOL ) );
  return 0 ;
}

// test the identity routine
static char *trace_test(  ) {
  int sum = 0 ; 
  int i ;  
  for( i = 0 ; i < NC ; i++ ) sum += i ;
  mu_assert("[GLUnit] error : trace is broken", !( creal( trace( a ) - sum * ( NC+1 )) > PREC_TOL ) );
  return 0 ;
}

// test the identity routine
static char *trace_ab_test(  ) {
  // multiply
  GLU_complex d[ NCNC ] , tr ;
  GLU_bool is_ok = GLU_TRUE ;
  multab( d , a , b ) ;
  trace_ab( &tr , a , b ) ;
  if( cabs( trace( d ) - tr ) > PREC_TOL ) is_ok = GLU_FALSE ;
  mu_assert("[GLUnit] error : trace_ab is broken", is_ok );
  return 0 ;
}

// test the identity routine
static char *trace_abc_test(  ) {
  GLU_complex d[ NCNC ] , e[ NCNC ] , tr ;
  GLU_bool is_ok = GLU_TRUE ;
  multab( d , b , c ) ;
  multab( e , a , d ) ;
  trace_abc( &tr , a , b , c ) ;
  if( cabs( trace( e ) - tr ) > PREC_TOL ) is_ok = GLU_FALSE ;
  mu_assert("[GLUnit] error : trace_abc is broken", is_ok ) ;
  return 0 ;
}

// test the identity routine
static char *trace_abc_dag_test(  ) {
  GLU_complex d[ NCNC ] , e[ NCNC ] , tr ;
  GLU_bool is_ok = GLU_TRUE ;
  multab_dag( d , b , c ) ;
  multab( e , a , d ) ;
  trace_abc_dag( &tr , a , b , c ) ;
  if( cabs( trace( e ) - tr ) > PREC_TOL ) is_ok = GLU_FALSE ;
  mu_assert("[GLUnit] error : trace_abc_dag is broken", is_ok ) ;

  GLU_real rtr ;
  trace_abc_dag_Re( &rtr , a , b , c ) ;
  if( fabs( rtr - creal( tr ) ) > PREC_TOL ) is_ok = GLU_FALSE ;
  mu_assert("[GLUnit] error : trace_abc_dag_Re is broken", is_ok ) ;

  return 0 ;
}

// test the identity routine
static char *trace_ab_herm_test(  ) {
  // need two hermitian matrices and to multiply them
  GLU_complex Ua[ NCNC ] , Ub[ NCNC ] ;
  GLU_complex Ha[ NCNC ] , Hb[ NCNC ] ;
  GLU_bool is_ok = GLU_TRUE ;

  // random gaussians have had RNG checked
  generate_NCxNC( Ua ) ;
  generate_NCxNC( Ub ) ;

  // reunitarise
  reunit2( Ua ) ;
  reunit2( Ub ) ;

  // projections have been checked already
  Hermitian_proj( Ha , Ua ) ;
  Hermitian_proj( Hb , Ua ) ;

  // multiply
  GLU_complex HaHb[ NCNC ] ;
  GLU_real tr ;
  multab( HaHb , Ha , Hb ) ;

  trace_ab_herm( &tr , Ha , Hb ) ;

  if( cabs( trace( HaHb ) - tr ) > PREC_TOL ) is_ok = GLU_FALSE ;
  mu_assert("[GLUnit] error : trace_ab_herm is broken", is_ok ) ;
  return 0 ;
}

// test the trace of ab_dag
static char *trace_ab_dag_test(  ) {
  GLU_complex Ua[ NCNC ] , Ub[ NCNC ] , Uab[ NCNC ] ;
  // random gaussians have had RNG checked
  generate_NCxNC( Ua ) ;
  generate_NCxNC( Ub ) ;

  // reunitarise
  reunit2( Ua ) ;
  reunit2( Ub ) ;

  multab_dag( Uab , Ua , Ub ) ;

  // compute the trace
  GLU_complex tr ;
  GLU_bool is_ok = GLU_TRUE ;
  trace_ab_dag( &tr , Ua , Ub ) ;

  if( cabs( trace( Uab ) - tr ) > PREC_TOL ) is_ok = GLU_FALSE ;
  mu_assert("[GLUnit] error : trace_ab_dag is broken", is_ok ) ;

  // test the real part routine
  GLU_real rtr ;
  trace_ab_dag_Re( &rtr , Ua , Ub ) ;
  if( fabs( rtr - creal( tr ) ) > PREC_TOL ) is_ok = GLU_FALSE ;
  mu_assert("[GLUnit] error : trace_ab_dag_Re is broken", is_ok ) ;

  return 0 ;
}

// test the trace of two hermitian matrices
static char *trace_ab_herm_short_test(  ) {
  // need two hermitian matrices and to multiply them
  GLU_complex Ua[ NCNC ] , Ub[ NCNC ] ;
  GLU_complex Ha[ NCNC ] , Hb[ NCNC ] ;
  GLU_complex hHa[ HERMSIZE ] , hHb[ HERMSIZE ] ;
  GLU_bool is_ok = GLU_TRUE ;

  // random gaussians have had RNG checked
  generate_NCxNC( Ua ) ;
  generate_NCxNC( Ub ) ;

  // reunitarise
  reunit2( Ua ) ;
  reunit2( Ub ) ;

  // projections have been checked already
  Hermitian_proj( Ha , Ua ) ;
  Hermitian_proj( Hb , Ua ) ;
  Hermitian_proj_short( hHa , Ua ) ;
  Hermitian_proj_short( hHb , Ua ) ;

  // multiply
  GLU_complex HaHb[ NCNC ] ;
  GLU_real tr ;
  multab( HaHb , Ha , Hb ) ;

  trace_ab_herm_short( &tr , hHa , hHb ) ;

  if( cabs( trace( HaHb ) - tr ) > PREC_TOL ) is_ok = GLU_FALSE ;
  mu_assert("[GLUnit] error : trace_ab_herm_short is broken", is_ok ) ;
  return 0 ;
}

// test the identity routine
static char *transpose_test(  ) {
  GLU_complex d[ NCNC ] ;
  GLU_bool is_ok = GLU_TRUE ;
  transpose( d , a ) ;
  int i , j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      if( cabs( d[ j + i * NC ] - a[ i + j * NC ] ) > PREC_TOL )
	is_ok = GLU_FALSE ;
    }
  }
  mu_assert("[GLUnit] error : transpose is broken", is_ok ) ;
  return 0 ;
}

// test the identity routine
static char *zero_mat_test(  ) {
  GLU_complex d[ NCNC ] ;
  double sum = 0.0 ;
  int i ;
  zero_mat( d ) ;
  for( i = 0 ; i < NCNC ; i++ ) {
    sum += cabs( d[i] ) ;
  }
  mu_assert("[GLUnit] error : zero_mat is broken", sum == 0.0 );
  return 0 ;
}


// test the rng
static char *rng_test( ) {
  const double res = 
#ifdef GSL_RNG
    0.0
#elif defined KISS_RNG
    0.985860379878432
#elif defined MWC_1038_RNG
    0.367735110223293
#elif defined MWC_4096_RNG
    0.165379548445344
#else
    0.575690010329708
#endif
    ;
  mu_assert( "[GLUnit] error : rng has changed", !( fabs( rng_dbl() - res ) > 1E-15 ) );
  return 0;
}

// test our random numbers
static char *rng_tests( ) {
  mu_run_test( rng_test ) ;
  return 0 ;
}

// test the exponential e^{0} = Identity?
static char *exponentiate_test(  ) {
  GLU_complex A[ NCNC ] , Id[ NCNC ] ;
  zero_mat( A ) ;
  exponentiate( Id , A ) ;
  int i ;
  GLU_bool is_ok = GLU_TRUE ;
  for( i = 0 ; i < NCNC ; i++ ) {
    if( i % ( NC + 1 ) == 0 ) {
      if( cabs( Id[ i ] - 1.0 ) > PREC_TOL ) is_ok = GLU_FALSE ;
    } else {
      if( cabs( Id[ i ] ) > PREC_TOL ) is_ok = GLU_FALSE ;
    }
  }
  mu_assert( "[GLUnit] error : exponentiate is broken", is_ok );
  return 0 ;
}

static char *exponentiate_short_test( void ) 
{
  // Hermitian proj has been tested
  GLU_complex hHa[ HERMSIZE ] , Ha[ NCNC ] ;
  GLU_complex hUa[ NCNC ] , Ua[ NCNC ] ;

  // create a unitary matrix
  generate_NCxNC( Ua ) ;
  reunit2( Ua ) ;

  // hermitian proj it
  Hermitian_proj( Ha , Ua ) ;
  Hermitian_proj_short( hHa , Ua ) ;

  exponentiate( Ua , Ha ) ;
  exponentiate( hUa , Ha ) ;

  GLU_bool is_ok = GLU_TRUE ;
  int i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    if( cabs( Ua[i] - hUa[i] ) > PREC_TOL ) is_ok = GLU_FALSE ;
  }
  mu_assert( "[GLUnit] error : exponentiate_short is broken", is_ok );
  return 0 ;
}

// log and exponentiate should be equivalent
static char *exact_log_test(  ) {
  GLU_complex A[ NCNC ] , Ua[ NCNC ] , eA[ NCNC ] ;

  // create a unitary matrix
  generate_NCxNC( Ua ) ;
  reunit2( Ua ) ;

  exact_log_slow( A , Ua ) ;
  exponentiate( eA , A ) ;
  
  int i ;
  GLU_bool is_ok = GLU_TRUE ;
  for( i = 0 ; i < NCNC ; i++ ) {
    if( cabs( eA[i] - Ua[i] ) > PREC_TOL ) is_ok = GLU_FALSE ;
  }
  mu_assert( "[GLUnit] error : exact_log_slow is broken", is_ok ) ;
  return 0 ;
}

// log and exponentiate should be equivalent
static char *exact_log_short_test(  ) {
  GLU_complex A[ NCNC ] , Ua[ NCNC ] , eA[ NCNC ] ;
  GLU_complex hA[ NCNC ] , heA[ NCNC ] ;

  // create a unitary matrix
  generate_NCxNC( Ua ) ;
  reunit2( Ua ) ;

  exact_log_slow( A , Ua ) ;
  exact_log_slow_short( hA , Ua ) ;

  exponentiate( eA , A ) ;
  exponentiate_short( heA , hA ) ;

  int i ;
  GLU_bool is_ok = GLU_TRUE ;
  for( i = 0 ; i < NCNC ; i++ ) {
    if( cabs( eA[i] - heA[i] ) > PREC_TOL ) is_ok = GLU_FALSE ;
  }
  mu_assert( "[GLUnit] error : exact_log_slow_short is broken", is_ok ) ;
  return 0 ;
}

// log-exponential tests
static char *logexp_tests( ) {
  mu_run_test( exponentiate_test ) ;
  mu_run_test( exponentiate_short_test ) ;
  mu_run_test( exact_log_test ) ;
  mu_run_test( exact_log_short_test ) ;
  return 0 ;
}

// matrix inverse tests
static char *inverse_test( ) {
  GLU_complex Ua[ NCNC ] ;

  // have been checked already
  generate_NCxNC( Ua ) ;
  reunit2( Ua ) ;

  // perform numerical inverse
  GLU_complex iUa[ NCNC ] ;
  inverse( iUa , Ua ) ;

  // must be checked before too
  multab_atomic_right( iUa , Ua ) ;
  int i ;
  GLU_bool is_ok = GLU_TRUE ;
  for( i = 0 ; i < NCNC ; i++ ) {
    if( i%(NC+1) == 0 ) {
      if( cabs( iUa[ i ] - 1.0 ) > PREC_TOL ) is_ok = GLU_FALSE ;
    } else {
      if( cabs( iUa[ i ] ) > PREC_TOL ) is_ok = GLU_FALSE ;
    }
  }

  mu_assert( "[GLUnit] error : inverse is broken", is_ok ) ;
  return 0 ;
}

// matrix inverse tests
static char *inversion_tests( ) {
  mu_run_test( inverse_test ) ;
  return 0 ;
}

// matrix multiply tests
static char *multiply_tests( ) {
  GLU_complex Ua[ NCNC ] , Ub[ NCNC ] , res[ NCNC ] ;

  // random gaussians have had RNG checked
  generate_NCxNC( Ua ) ;
  generate_NCxNC( Ub ) ;

  // reunitarise
  reunit2( Ua ) ;
  mu_assert("[GLUnit] error : reunit2 is broken", !is_unitary( Ua ) );

  reunit2( Ub ) ;
  mu_assert("[GLUnit] error : reunit2 is broken", !is_unitary( Ub ) );

  // perform matrix multiplies and compare SU(N) version to standard
  // looped version
  printf( "multab \n" ) ;
  multab_suNC( res , Ua , Ub ) ;
  printf( "SU(N)\n" ) ;
  write_matrix( res ) ;
  multab( res , Ua , Ub ) ;
  printf( "General\n" ) ;
  write_matrix( res ) ;
  printf( "multab_atomic_right a.b\n" ) ;
  equiv( res , Ua ) ;
  multab_atomic_right( res , Ub ) ;
  write_matrix( res ) ;
  printf( "multab_atomic_left a.b\n" ) ;
  equiv( res , Ub ) ;
  multab_atomic_left( res , Ua ) ;
  write_matrix( res ) ;
  printf( "multabdag \n" ) ;
  multabdag_suNC( res , Ua , Ub ) ;
  printf( "SU(N)\n" ) ;
  write_matrix( res ) ;
  printf( "General\n" ) ;
  multabdag( res , Ua , Ub ) ;
  printf( "SU(N)\n" ) ;
  write_matrix( res ) ;
  printf( "multabdagdag \n" ) ;
  multab_dagdag_suNC( res , Ua , Ub ) ;
  printf( "SU(N)\n" ) ;
  write_matrix( res ) ;
  printf( "General\n" ) ;
  multab_dagdag( res , Ua , Ub ) ;
  printf( "SU(N)\n" ) ;
  write_matrix( res ) ;
  printf( "multab_dag \n" ) ;
  multab_dag_suNC( res , Ua , Ub ) ;
  printf( "SU(N)\n" ) ;
  write_matrix( res ) ;
  printf( "General\n" ) ;
  multab_dag( res , Ua , Ub ) ;
  write_matrix( res ) ;

  return 0 ;
}

// run all the tests
static char *all_tests( ) {
  // important tests
  mu_run_test( equiv_test ) ;
  mu_run_test( zero_mat_test ) ;
  mu_run_test( trace_test ) ;
  mu_run_test( dagger_test ) ;

  mu_run_test( add_constant_test ) ;
  mu_run_test( a_plus_b_test ) ;
  mu_run_test( a_plus_Sxbminc_short_test ) ;
  mu_run_test( a_plus_CSxb_test ) ;
  mu_run_test( a_plus_Sxb_test ) ;
  mu_run_test( b_min_c_test ) ;
  mu_run_test( diag_test ) ;
  mu_run_test( diag_vect_test ) ;
  mu_run_test( identity_test ) ;
  mu_run_test( is_unitary_test ) ;
#if NC>3
  mu_run_test( LU_det_test ) ;
#endif
  mu_run_test( det_test ) ;
  mu_run_test( cofactor_transpose_test ) ;
  mu_run_test( mat_mult_vec_test ) ;
  mu_run_test( matrix_power_test ) ;
  mu_run_test( M_times_c_test ) ;
  mu_run_test( pack_hermitian_test ) ;
  mu_run_test( rebuild_hermitian_test ) ;
  mu_run_test( rebuild_antihermitian_test ) ;
  mu_run_test( outerproduct_test ) ;
  mu_run_test( shorten_test ) ;
  mu_run_test( speed_det_test ) ;
  mu_run_test( speed_trace_test ) ;
  mu_run_test( trace_ab_test ) ;
  mu_run_test( trace_abc_test ) ;
  mu_run_test( trace_abc_dag_test ) ;
  mu_run_test( transpose_test ) ;
  mu_run_test( trace_ab_herm_test ) ;
  mu_run_test( trace_ab_herm_short_test ) ;
  mu_run_test( trace_ab_dag_test ) ;

  return 0 ;
}

struct site *lat , *glat ;

static char *av_plaquette_test( ) {
  GLU_bool is_ok = GLU_TRUE ;

  double sp_plaq , t_plaq , gsp_plaq , gt_plaq ;
  all_plaquettes( lat , &sp_plaq , &t_plaq ) ;
  all_plaquettes( glat , &gsp_plaq , &gt_plaq ) ;

  if( fabs( sp_plaq - gsp_plaq ) > PREC_TOL ) is_ok = GLU_FALSE ;
  mu_assert( "[GLUnit] error :spatial  plaquttes are broken", is_ok ) ;

  if( fabs( t_plaq - gt_plaq ) > PREC_TOL ) is_ok = GLU_FALSE ;
  mu_assert( "[GLUnit] error : temporal plaquttes are broken", is_ok ) ;

  return 0 ;
}

// test the gauge invariance of the polyakov loops in each direction
// after a gauge transformation
static char *polyakov_test( ) {
  GLU_bool is_ok = GLU_TRUE ;
  int mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    // soften this a little
    if( cabs( poly( lat , mu ) - poly( glat , mu ) ) > 10*PREC_TOL ) {
      is_ok = GLU_FALSE ;
    }
  }
  mu_assert( "[GLUnit] error : polyakov loops are broken", is_ok ) ;
  return 0 ;
}

// configuration-dependent startup
static char *config_tests( ) {
  // 4^ND periodic lattice
  int mu ; 
  for( mu = 0 ; mu < ND ; mu++ ) {
    Latt.dims[ mu ] = 4 ;
  }

  // initialise geometry so that we can use LVOLUME and stuff
  init_geom( ) ;

  // malloc our gauge field and initialise our lattice geometry
  lat = malloc( LVOLUME * sizeof ( struct site ) ) ;
  init_navig( lat ) ;

  // allocate gauge transformed fields
  glat = malloc( LVOLUME * sizeof ( struct site ) ) ;
  init_navig( glat ) ;

  // randomly generate an SU(NC) field
  random_suNC( lat ) ;

  int i ;
  for( i = 0 ; i < LVOLUME ; i++ ) {
    for( mu = 0 ; mu < ND ; mu++ ) {
      equiv( glat[i].O[mu] , lat[i].O[mu] ) ;
    }
  }

  // assumes gtrans works properly
  random_gtrans( glat ) ;

  // test gauge invariant stuff
  mu_run_test( av_plaquette_test ) ;
  mu_run_test( polyakov_test ) ;

  free( lat ) ;
  free( glat ) ;

  return 0 ;
}

// main
int main( const int argc , const char *argv[] )
{
  // inits
  attach_GLU( ) ;

  // set the seed to something I know answers for
  Latt.Seed[0] = 1 ;

  // initialise rng
  rng_init( ) ;

  // test rng
  char *test_results ;
  if( ( test_results = rng_tests( ) ) != 0 ) {
     goto TEST_FAILURE ;
  } 

  if( ( test_results = multiply_tests( ) ) != 0 ) {
    goto TEST_FAILURE ;
  } 

  if( ( test_results = logexp_tests( ) ) != 0 ) {
    goto TEST_FAILURE ;
  } 

  if( ( test_results = inversion_tests( ) ) != 0 ) {
    goto TEST_FAILURE ;
  } 

  // set the matrices' entries be their lexicographical index
  int i, j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      a[ j + i*NC ] = b[ j + i*NC ] = c[ j + i*NC ] = 
	(GLU_real)( j + i*NC ) + I * (GLU_real)( j + i*NC ) ;
    }
  }

  if( ( test_results = all_tests( ) ) != 0 ) {
    goto TEST_FAILURE ;
  } 

  if( ( test_results = config_tests( ) ) != 0 ) {
    goto TEST_FAILURE ;
  }

  // if we are successful we go there
  goto TEST_SUCCESS ;

 TEST_FAILURE :
  unstick_GLU( ) ;
  printf( "[GLUnit] Test failure \n%s \n" , test_results ) ;
  return GLU_FAILURE ;
 TEST_SUCCESS :
  unstick_GLU( ) ;
  printf( "[GLUnit] All %d tests passed \n" , tests_run ) ;
  return GLU_SUCCESS ;
}
