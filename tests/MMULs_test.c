/**
   @file MMULs_test.c
   @brief test our matrix multiplies
 */
#include "Mainfile.h"

#include "par_rng.h"
#include "gramschmidt.h"
#include "minunit.h"

// temporary mats
static GLU_complex Ua[ NCNC ] GLUalign , Ub[ NCNC ] GLUalign , Uc[NCNC] ;
static GLU_complex res1[ NCNC ] GLUalign , res2[ NCNC ] GLUalign , res3[NCNC] ;

// equivalence test for matrices
static GLU_bool
are_equal( const GLU_complex a[ NCNC ] ,
	   const GLU_complex b[ NCNC ] ) 
{
  size_t i ;
  for( i = 0 ; i < NCNC ; i++ ) {
    if( cabs( a[i] - b[i] ) > PREC_TOL ) {
      return GLU_FALSE ;
    }
  }
  return GLU_TRUE ;
}

// check we can gram-schmidt a matrix 
static char *reunit_test( void )
{
  // random gaussians have had RNG checked
  Sunitary_gen( Ua , 0 ) ;
  Sunitary_gen( Ub , 0 ) ;
  Sunitary_gen( Uc , 0 ) ;

  mu_assert( "[GLUnit] error : Sunitary_gen is broken" , is_unitary( Ua ) ) ;
  mu_assert( "[GLUnit] error : Sunitary_gen is broken" , is_unitary( Ub ) ) ;
  mu_assert( "[GLUnit] error : Sunitary_gen is broken" , is_unitary( Uc ) ) ;
  return NULL ;
}

// perform matrix multiplies and compare SU(N) version to standard
// looped version
static char *multab_test( void )
{
  multab_suNC( res1 , Ua , Ub ) ;
  multab( res2 , Ua , Ub ) ;
  if( are_equal( res1 , res2 ) == GLU_FALSE ) {
    write_matrix( res1 ) ;
    write_matrix( res2 ) ; 
    mu_assert( "[GLUnit] error : MMULS differ" , 0 ) ;
  }
  return NULL ;
}

// test atomic multiply to the right
static char *multab_atomic_left_test( void )
{
  // multab atomic right
  multab( res2 , Ua , Ub ) ;
  equiv( res1 , Ub ) ;
  multab_atomic_left( res1 , Ua ) ;
  if( are_equal( res1 , res2 ) == GLU_FALSE ) {
    write_matrix( res1 ) ;
    write_matrix( res2 ) ;
    mu_assert( "[GLUnit] error : MMUL atomic left broken" , 0 ) ;
  }
  return NULL ;
}

// test atomic multiply to the right
static char *multab_atomic_right_test( void )
{
  // multab atomic right
  multab( res2 , Ua , Ub ) ;
  equiv( res1 , Ua ) ;
  multab_atomic_right( res1 , Ub ) ;
  if( are_equal( res1 , res2 ) == GLU_FALSE ) {
    write_matrix( res1 ) ;
    write_matrix( res2 ) ;
    mu_assert( "[GLUnit] error : MMUL atomic right broken" , 0 ) ;
  }
  return NULL ;
}

// test suNC multabdag 
static char *multabdag_test( void ) 
{
  multabdag_suNC( res1 , Ua , Ub ) ;
  multabdag( res2 , Ua , Ub ) ;
  if( are_equal( res1 , res2 ) == GLU_FALSE ) {
    write_matrix( res1 ) ;
    write_matrix( res2 ) ;
    mu_assert( "[GLUnit] error : multabdag broken" , 0 ) ;
  }
  return NULL ;
}

// test suNC multabdagdag 
static char *multab_dagdag_test( void ) 
{
  multab_dagdag_suNC( res1 , Ua , Ub ) ;
  multab_dagdag( res2 , Ua , Ub ) ;
  if( are_equal( res1 , res2 ) == GLU_FALSE ) {
    write_matrix( res1 ) ;
    write_matrix( res2 ) ;
    mu_assert( "[GLUnit] error : multab_dagdag broken" , 0 ) ;
  }
  // test that Ua^{\dagger}.Ub^{\dagger} = [Ub.Ua]^{\dagger} ?
  return NULL ;
}

// test suNC multab_dag 
static char *multab_dag_test( void ) 
{
  multab_dag_suNC( res1 , Ua , Ub ) ;
  multab_dag( res2 , Ua , Ub ) ;

  if( are_equal( res1 , res2 ) == GLU_FALSE ) {
    write_matrix( res1 ) ;
    write_matrix( res2 ) ;
    mu_assert( "[GLUnit] error : multab_dag" , 0 ) ;
  }
  // test that Ua^{\dagger}.Ub^{\dagger} = [Ub.Ua]^{\dagger} ?
  return NULL ;
}

// test bcd^dag
static char *multabcdag_test( void )
{
  for( int i = 0 ; i < 100000 ; i++ ) {
    multabcdag_suNC( res1 , Ua , Ub , Uc ) ;
    // slower guy    
    multab_suNC( res3 , Ua , Ub ) ;
    multab_dag_suNC( res2 , res3 , Uc ) ;
  }
  
  if( are_equal( res1 , res2 ) == GLU_FALSE ) {
    write_matrix( res1 ) ;
    write_matrix( res2 ) ;
    mu_assert( "[GLUnit] error : multabcdag" , 0 ) ;
  }
  return NULL ;
}

// test bcd^dag
static char *multadagbc_test( void )
{
  for( int i = 0 ; i < 100000 ; i++ ) {
    multadagbc_suNC( res1 , Ua , Ub , Uc ) ;
    // slower guy
    multab_suNC( res3 , Ub , Uc ) ;
    multabdag_suNC( res2 , Ua , res3 ) ;
  }
  
  if( are_equal( res1 , res2 ) == GLU_FALSE ) {
    write_matrix( res1 ) ;
    write_matrix( res2 ) ;
    mu_assert( "[GLUnit] error : multadagbc" , 0 ) ;
  }
  return NULL ;
}

// little runner for the tests
static char *
mmul_test( void )
{
  mu_run_test( reunit_test ) ;
  mu_run_test( multab_test ) ;
  mu_run_test( multab_atomic_left_test ) ;
  mu_run_test( multab_atomic_right_test ) ;
  mu_run_test( multabdag_test ) ;
  mu_run_test( multab_dagdag_test ) ;
  mu_run_test( multab_dag_test ) ;
  mu_run_test( multabcdag_test ) ;
  mu_run_test( multadagbc_test ) ;
  return NULL ;
}

int
MMUL_test( void )
{
  // init to zero again
  tests_run = tests_fail = 0 ;

  // initial gamma setup and test
  char *mmulres = mmul_test( ) ;

  if( tests_fail == 0 ) {
    printf( "[MMUL UNIT] all %d tests passed\n\n" ,
	    tests_run ) ;
    return GLU_SUCCESS ;
  } else {
    printf( "%s \n" , mmulres ) ;
    printf( "[MMUL UNIT] %d out of %d tests failed\n\n" , 
	    tests_fail , tests_run ) ;
    return GLU_FAILURE ;
  }
}
