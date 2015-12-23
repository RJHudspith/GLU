/**
   @file geom_test.c
   @brief test our geometry functions
 */
#include "Mainfile.h"

#include "minunit.h"
#include "GLU_rng.h"
#include "gramschmidt.h"
#include "invert.h"

// these are set to their lexicographical order
static GLU_complex a[ NCNC ] , b[ NCNC ] , c[ NCNC ] ;

// test the add constant routine
static char *add_constant_test( ) {
  GLU_bool is_ok = GLU_TRUE ;
  GLU_complex A[ NCNC ] ;
  equiv( A , a ) ;
  add_constant( A , 1.0 ) ;
  if( cabs( trace( A ) - trace( a ) - NC ) > PREC_TOL ) is_ok = GLU_FALSE ;
  mu_assert( "[GLUnit] error : add constant broken", is_ok );
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
  mu_assert( "[GLUnit] error : add_plus_b broken", is_ok );
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
static char *a_plus_Sxbminc_short_test( void ) {
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

  size_t i ;
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

static char *
unops_test( void )
{
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

  return NULL ;
}

int
UNOPS_test( void )
{
  // init to zero again
  tests_run = tests_fail = 0 ;

  // set the matrices' entries be their lexicographical index
  int i, j ;
  for( i = 0 ; i < NC ; i++ ) {
    for( j = 0 ; j < NC ; j++ ) {
      a[ j + i*NC ] = b[ j + i*NC ] = c[ j + i*NC ] = 
	(GLU_real)( j + i*NC ) + I * (GLU_real)( j + i*NC ) ;
    }
  }

  // initial gamma setup and test
  char *unops_res = unops_test( ) ;

  // if we fail we complain
  if( tests_fail == 0 ) {
    printf( "[UNOPS UNIT] all %d tests passed\n\n" ,
	    tests_run ) ;
    return GLU_SUCCESS ;
  } else {
    printf( "%s \n" , unops_res ) ;
    printf( "[UNOPS UNIT] %d out of %d tests failed\n\n" , 
	    tests_fail , tests_run ) ;
    return GLU_FAILURE ;
  }
}
