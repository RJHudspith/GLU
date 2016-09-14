/**
   @file geom_test.c
   @brief test our geometry functions
 */
#include "Mainfile.h"

#include "minunit.h"
#include "par_rng.h"
#include "gramschmidt.h"
#include "invert.h"

// these are set to their lexicographical order
static GLU_complex a[ NCNC ] GLUalign , b[ NCNC ] GLUalign , 
  c[ NCNC ] GLUalign ;

// test the add constant routine
static char *add_constant_test( void ) {
  GLU_bool is_ok = GLU_TRUE ;
  GLU_complex A[ NCNC ] GLUalign ;
  equiv( A , a ) ;
  add_constant( A , 1.0 ) ;
  if( cabs( trace( A ) - trace( a ) - NC ) > PREC_TOL ) is_ok = GLU_FALSE ;
  mu_assert( "[GLUnit] error : add constant broken", is_ok );
  return 0;
}

// test the a+b routine
static char *a_plus_b_test( void ) {
  GLU_bool is_ok = GLU_TRUE ;
  GLU_complex d[ NCNC ] GLUalign ;
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
static char *a_plus_CSxb_test( void ) {
  GLU_bool is_ok = GLU_TRUE ;
  GLU_complex d[ NCNC ] GLUalign ;
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
static char *a_plus_Sxb_test( void ) {
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
  GLU_complex hHa[ HERMSIZE ] , Ua[ NCNC ] GLUalign ;

  par_generate_NCxNC( Ua , 0 ) ;
  gram_reunit( Ua ) ;

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
static char *b_min_c_test( void ) {
  GLU_bool is_ok = GLU_TRUE ;
  GLU_complex d[ NCNC ] GLUalign ;
  b_min_c( d , b , c ) ;
  int i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    if( cabs( d[i] ) > PREC_TOL ) is_ok = GLU_FALSE ;
  } 
  mu_assert("[GLUnit] error : b_min_c broken", is_ok );
  return 0;
}

// test the cofactor transpose? -> probably dated
static char *cofactor_transpose_test( void ) {
  // cofactor_transpose / determinant is the inverse!
  GLU_complex Ua[ NCNC ] GLUalign ;

  // create random SU(N) matrix
  Sunitary_gen( Ua , 0 ) ;

  GLU_complex ctUa[ NCNC ] GLUalign ;
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
static char *equiv_test( void ) {
  GLU_complex d[ NCNC ] GLUalign ;
  identity( d ) ;
  equiv( d , a ) ;
  double sum = 0.0 ;
  int i ;
  for( i = 0 ; i < NCNC ; i++ ) sum += cabs( d[i] - a[i] ) ; 
  mu_assert("[GLUnit] error : equiv broken", !( sum > PREC_TOL ) );
  return 0;
}

// test the daggering
static char *dagger_test( void ) {
  GLU_bool is_ok = GLU_TRUE ;
  GLU_complex d[ NCNC ] GLUalign ;
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
static char *diag_test( void ) {
  GLU_bool is_ok = GLU_TRUE ;
  GLU_complex d[ NCNC ] GLUalign , sum = 0.0 ;
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
static char *diag_vect_test( void ) {
  GLU_bool is_ok = GLU_TRUE ;
  GLU_complex d[ NCNC ] GLUalign , v[ NC ] , sum = 0.0 ;
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
static char *identity_test( void ) {
  GLU_bool is_ok = GLU_TRUE ;
  GLU_complex d[ NCNC ] GLUalign ;
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
static char *is_unitary_test( void ) {
  GLU_complex d[ NCNC ] ;
  identity( d ) ;
  mu_assert("[GLUnit] error : is_unitary broken", is_unitary( d ) );
  return 0;
}

#if NC > 3
// test the LU_det comparing to the identity for now
static char *LU_det_test( void ) {
  GLU_complex d[ NCNC ] ;
  identity( d ) ;
  mu_assert("[GLUnit] error : LU_det broken", 
	    cabs( det( d ) - 1.0 ) > PREC_TOL ? GLU_FALSE : GLU_TRUE  );
  return 0;
}
#endif

// test the determinant of the identity
static char *det_test( void ) {
  GLU_complex d[ NCNC ] GLUalign ;
  identity( d ) ;
  mu_assert("[GLUnit] error : det broken", 
	    cabs( det( d ) - 1.0 ) > PREC_TOL ? GLU_FALSE : GLU_TRUE  );
  return 0;
}

// test matrix multiply a vector
static char *mat_mult_vec_test( void ) {
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
static char *M_times_c_test( void ) {
  GLU_complex d[ NCNC ] GLUalign ;
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
static char *matrix_power_test( void ) {
  GLU_complex Ua[ NCNC ] GLUalign , Ub[ NCNC ] GLUalign ;

  // generate a random matrix
  Sunitary_gen( Ua , 0 ) ;

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
static char *outerproduct_test( void ) {
  GLU_complex o[ NC ] , v[ NC ] ;
  int i ;
  GLU_bool is_ok = GLU_TRUE ;
  for( i = 0 ; i < NC ; i++ ) {
    o[ i ] = v[ i ] = 1.0 ;
  }
  GLU_complex d[ NCNC ] GLUalign ;
  outerproduct( d , o , v ) ;
  for( i = 0 ; i < NCNC ; i++ ) {
    if( cabs( d[i] - 1.0 ) > PREC_TOL ) is_ok = GLU_FALSE ;
  }
  mu_assert("[GLUnit] error : outerproduct is broken", is_ok );
  return 0 ;
}

// test the hermitian packing
static char *pack_hermitian_test( void ) {
  GLU_complex Ha[ NCNC ] GLUalign , sHa[ HERMSIZE ] ;
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
static char *rebuild_antihermitian_test( void ) {
  GLU_complex Ha[ NCNC ] GLUalign , sHa[ HERMSIZE ] , 
    rHa[ NCNC ] GLUalign ;
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
static char *rebuild_hermitian_test( void ) {
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

// test the identity routine
static char *Re_trace_abc_dag_test( void ) {
  GLU_complex A[ NCNC ] GLUalign , B[ NCNC ] GLUalign ,
    C[ NCNC ] GLUalign ;
  Sunitary_gen( A , 0 ) ; Sunitary_gen( B , 0 ) ; Sunitary_gen( C , 0 ) ;
  GLU_real tr ;
  trace_abc_dag_Re( &tr , A , B , C ) ;
  mu_assert("[GLUnit] error : Re_trace_abc_dag is broken", 
	    !( fabs( tr - Re_trace_abc_dag_suNC( A , B , C ) ) 
	       > PREC_TOL ) ) ;
  return 0 ;
}

// test the shortening and rebuilding
static char *shorten_test( void ) {
  GLU_bool is_ok = GLU_TRUE ;
#if NC < 4
  GLU_real sA[ NCNC-1 ] ;
  GLU_complex Ua[ NCNC ] GLUalign , rUa[ NCNC ] GLUalign ;

  // random SU(N) matrix
  Sunitary_gen( Ua , 0 ) ;
  
  // test if shorten -> rebuild gives same answer
  shorten( sA , Ua ) ;
  rebuild( rUa , sA ) ;

  size_t i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    if( cabs( rUa[ i ] - Ua[ i ] ) > PREC_TOL ) is_ok = GLU_FALSE ;
  }
#endif
  mu_assert("[GLUnit] error : shorten is broken",  is_ok );
  return 0 ;
}

// test the identity routine
static char *speed_det_test( void ) {
  GLU_complex d[ NCNC ] GLUalign , s ;
  identity( d ) ;
  speed_det( &s , d ) ;
  mu_assert("[GLUnit] error : det broken", 
	    cabs( s - 1.0 ) > PREC_TOL ? GLU_FALSE : GLU_TRUE  );
  return 0 ;
}

// test the identity routine
static char *speed_trace_test( void ) {
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
static char *trace_test( void ) {
  int sum = 0 ; 
  int i ;  
  for( i = 0 ; i < NC ; i++ ) sum += i ;
  mu_assert("[GLUnit] error : trace is broken", !( creal( trace( a ) - sum * ( NC+1 )) > PREC_TOL ) );
  return 0 ;
}

// test the identity routine
static char *trace_ab_test( void ) {
  // multiply
  GLU_complex d[ NCNC ] GLUalign , tr ;
  GLU_bool is_ok = GLU_TRUE ;
  multab( d , a , b ) ;
  trace_ab( &tr , a , b ) ;
  if( cabs( trace( d ) - tr ) > PREC_TOL ) is_ok = GLU_FALSE ;
  mu_assert("[GLUnit] error : trace_ab is broken", is_ok );
  return 0 ;
}

// test the identity routine
static char *trace_abc_test( void ) {
  GLU_complex d[ NCNC ] GLUalign , e[ NCNC ] GLUalign , tr ;
  GLU_bool is_ok = GLU_TRUE ;
  multab( d , b , c ) ;
  multab( e , a , d ) ;
  trace_abc( &tr , a , b , c ) ;
  if( cabs( trace( e ) - tr ) > PREC_TOL ) is_ok = GLU_FALSE ;
  mu_assert("[GLUnit] error : trace_abc is broken", is_ok ) ;
  return 0 ;
}

// test the tr[ a.b.c^dag ] routine
static char *trace_abc_dag_test( void ) {
  GLU_complex d[ NCNC ] GLUalign , e[ NCNC ] GLUalign , tr ;
  GLU_bool is_ok = GLU_TRUE ;
  multab_dag( d , b , c ) ;
  multab( e , a , d ) ;
  trace_abc_dag( &tr , a , b , c ) ;

  if( cabs( trace( e ) - tr ) > PREC_TOL ) { 
    is_ok = GLU_FALSE ;
  }

  GLU_real rtr ;
  trace_abc_dag_Re( &rtr , a , b , c ) ;

  if( fabs( rtr - creal( tr ) ) > PREC_TOL ) { 
    is_ok = GLU_FALSE ;
  }
  mu_assert("[GLUnit] error : trace_abc_dag_Re is broken", is_ok ) ;

  return 0 ;
}

// test the identity routine
static char *trace_ab_herm_test( void ) {
  // need two hermitian matrices and to multiply them
  GLU_complex Ua[ NCNC ] GLUalign , Ub[ NCNC ] GLUalign ;
  GLU_complex Ha[ NCNC ] GLUalign , Hb[ NCNC ] GLUalign ;
  GLU_bool is_ok = GLU_TRUE ;

  // random gaussians have had RNG checked
  Sunitary_gen( Ua , 0 ) ;
  Sunitary_gen( Ub , 0 ) ;

  // projections have been checked already
  Hermitian_proj( Ha , Ua ) ;
  Hermitian_proj( Hb , Ua ) ;

  // multiply
  GLU_complex HaHb[ NCNC ] GLUalign ;
  GLU_real tr ;
  multab( HaHb , Ha , Hb ) ;

  trace_ab_herm( &tr , Ha , Hb ) ;

  if( cabs( trace( HaHb ) - tr ) > PREC_TOL ) is_ok = GLU_FALSE ;
  mu_assert("[GLUnit] error : trace_ab_herm is broken", is_ok ) ;
  return 0 ;
}

// test the trace of ab_dag
static char *trace_ab_dag_test( void ) {
  GLU_complex Ua[ NCNC ] GLUalign , Ub[ NCNC ] GLUalign , 
    Uab[ NCNC ] GLUalign ;

  // random gaussians have had RNG checked
  Sunitary_gen( Ua , 0 ) ;
  Sunitary_gen( Ub , 0 ) ;

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
static char *trace_ab_herm_short_test( void ) {
  // need two hermitian matrices and to multiply them
  GLU_complex Ua[ NCNC ] GLUalign , Ub[ NCNC ] GLUalign ;
  GLU_complex Ha[ NCNC ] GLUalign , Hb[ NCNC ] GLUalign ;
  GLU_complex hHa[ HERMSIZE ] , hHb[ HERMSIZE ] ;
  GLU_bool is_ok = GLU_TRUE ;

  // random gaussians have had RNG checked
  Sunitary_gen( Ua , 0 ) ;
  Sunitary_gen( Ub , 0 ) ;

  // projections have been checked already
  Hermitian_proj( Ha , Ua ) ;
  Hermitian_proj( Hb , Ua ) ;
  Hermitian_proj_short( hHa , Ua ) ;
  Hermitian_proj_short( hHb , Ua ) ;

  // multiply
  GLU_complex HaHb[ NCNC ] GLUalign ;
  GLU_real tr ;
  multab( HaHb , Ha , Hb ) ;

  trace_ab_herm_short( &tr , hHa , hHb ) ;

  if( cabs( trace( HaHb ) - tr ) > PREC_TOL ) is_ok = GLU_FALSE ;
  mu_assert("[GLUnit] error : trace_ab_herm_short is broken", is_ok ) ;
  return 0 ;
}

// test the identity routine
static char *transpose_test( void ) {
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
static char *zero_mat_test( void ) {
  GLU_complex d[ NCNC ] GLUalign ;
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
  mu_run_test( Re_trace_abc_dag_test ) ;
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
    fprintf( stdout , "[UNOPS UNIT] all %d tests passed\n\n" ,
	    tests_run ) ;
    return GLU_SUCCESS ;
  } else {
    fprintf( stderr , "%s \n" , unops_res ) ;
    fprintf( stderr , "[UNOPS UNIT] %d out of %d tests failed\n\n" , 
	    tests_fail , tests_run ) ;
    return GLU_FAILURE ;
  }
}
