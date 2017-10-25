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
#include "init.h"
#include "par_rng.h"
#include "GLUlib_wrap.h"
#include "gramschmidt.h"
#include "invert.h"
#include "LU.h"
#include "minunit.h"
#include "plaqs_links.h"
#include "POLY.h"
#include "random_config.h"

#include "geom_test.h"
#include "MMULs_test.h"
#include "U_Nops_test.h"

// these are set to their lexicographical order
struct site *lat , *glat ;

struct latt_info Latt ; // dimensions and stuff
GLU_bool INIT_RNG ;
int tests_run = 0 ;
int tests_fail = 0 ;

// test the rng
static char *rng_test( void ) {
  /*
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
  */
  return 0;
}

// test our random numbers
static char *rng_tests( void ) {
  mu_run_test( rng_test ) ;
  return 0 ;
}

// test the exponential e^{0} = Identity?
static char *exponentiate_test( void ) {
  GLU_complex A[ NCNC ] GLUalign , Id[ NCNC ] GLUalign ;
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
  GLU_complex hHa[ HERMSIZE ] , Ha[ NCNC ] GLUalign ;
  GLU_complex hUa[ NCNC ] GLUalign , Ua[ NCNC ] GLUalign ;

  // create a unitary matrix
  Sunitary_gen( Ua , 0 ) ;

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
static char *exact_log_test( void ) {
  GLU_complex A[ NCNC ] , Ua[ NCNC ] , eA[ NCNC ] ;

  // create a unitary matrix
  Sunitary_gen( Ua , 0 ) ;

  exact_log_slow( A , Ua ) ;
  exponentiate( eA , A ) ;
  
  int i ;
  GLU_bool is_ok = GLU_TRUE ;
  for( i = 0 ; i < NCNC ; i++ ) {
    if( cabs( eA[i] - Ua[i] ) > NC*PREC_TOL ) is_ok = GLU_FALSE ;
  }
#if NC < 15
  mu_assert( "[GLUnit] error : exact_log_slow is broken", is_ok ) ;
#endif
  return 0 ;
}

// log and exponentiate should be equivalent
static char *exact_log_short_test( void ) {
  GLU_complex A[ NCNC ] GLUalign , Ua[ NCNC ] GLUalign , 
    eA[ NCNC ] GLUalign ;
  GLU_complex hA[ NCNC ] GLUalign , heA[ NCNC ] GLUalign ;

  // create a unitary matrix
  Sunitary_gen( Ua , 0 ) ;

  exact_log_slow( A , Ua ) ;
  exact_log_slow_short( hA , Ua ) ;

  exponentiate( eA , A ) ;
  exponentiate_short( heA , hA ) ;

  size_t i ;
  GLU_bool is_ok = GLU_TRUE ;
  for( i = 0 ; i < NCNC ; i++ ) {
    if( cabs( eA[i] - heA[i] ) > NC*PREC_TOL ) is_ok = GLU_FALSE ;
  }
#if NC < 15
  mu_assert( "[GLUnit] error : exact_log_slow_short is broken", is_ok ) ;
#endif
  return 0 ;
}

// log-exponential tests
static char *logexp_tests( void ) {
  mu_run_test( exponentiate_test ) ;
  mu_run_test( exponentiate_short_test ) ;
  mu_run_test( exact_log_test ) ;
  mu_run_test( exact_log_short_test ) ;
  return 0 ;
}

// matrix inverse tests
static char *inverse_test( void ) {
  GLU_complex Ua[ NCNC ] GLUalign ;

  // have been checked already
  Sunitary_gen( Ua , 0 ) ;

  // perform numerical inverse
  GLU_complex iUa[ NCNC ] GLUalign ;
  inverse( iUa , Ua ) ;

  // must be checked before too
  multab_atomic_right( iUa , Ua ) ;

  size_t i ;
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
static char *inversion_tests( void ) {
  mu_run_test( inverse_test ) ;
  return 0 ;
}

static char *av_plaquette_test( void ) {
  GLU_bool is_ok = GLU_TRUE ;

  double sp_plaq , t_plaq , gsp_plaq , gt_plaq ;
  all_plaquettes( lat , &sp_plaq , &t_plaq ) ;
  all_plaquettes( glat , &gsp_plaq , &gt_plaq ) ;

  if( fabs( sp_plaq - gsp_plaq ) > PREC_TOL ) is_ok = GLU_FALSE ;
  mu_assert( "[GLUnit] error :spatial  plaquttes are broken" , is_ok ) ;

  if( fabs( t_plaq - gt_plaq ) > PREC_TOL ) is_ok = GLU_FALSE ;
  mu_assert( "[GLUnit] error : temporal plaquttes are broken" , is_ok ) ;

  return 0 ;
}

// test the gauge invariance of the polyakov loops in each direction
// after a gauge transformation
static char *polyakov_test( void ) {
  GLU_bool is_ok = GLU_TRUE ;
  int mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    // soften this a little
    if( cabs( poly( lat , mu ) - poly( glat , mu ) ) > 10*PREC_TOL ) {
      is_ok = GLU_FALSE ;
    }
  }
  mu_assert( "[GLUnit] error : polyakov loops are broken" , is_ok ) ;
  return 0 ;
}

// configuration-dependent startup
static char *config_tests( void ) 
{
  // malloc our gauge field and initialise our lattice geometry
  lat = NULL ; 
  if( ( lat = allocate_lat( ) ) == NULL ) {
    fprintf( stderr , "[CONFIG-UNIT] Gauge field allocation failure\n" ) ;
    return NULL ;
  }
  init_navig( lat ) ;

  // randomly generate an SU(NC) field
  random_suNC( lat ) ;

  // allocate gauge transformed fields
  glat = NULL ; 
  if( ( glat = allocate_lat( ) ) == NULL ) {
    fprintf( stderr , "[CONFIG-UNIT] Gauge field 2 allocation failure\n" ) ;
    return NULL ;
  }
  init_navig( glat ) ;

  size_t i , mu ;
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
int main( void )
{
  // inits
  attach_GLU( ) ;

  // 4^ND periodic lattice
  size_t mu ; 
  for( mu = 0 ; mu < ND ; mu++ ) {
    Latt.dims[ mu ] = 4 ;
  }
  
  // initialise geometry so that we can use LVOLUME and stuff
  init_latt( ) ;

  // set the seed to something I know answers for
  Latt.Seed[0] = 1 ;

  // initialise rng
  initialise_par_rng( NULL ) ;

  // a counter for all the tests we have done
  int full_tests_run = 0 ;

  // test rng
  char *test_results ;
  if( ( test_results = rng_tests( ) ) != 0 ) {
     goto TEST_FAILURE ;
  }
  full_tests_run += tests_run ;

  // matrix multiply tests
  if( MMUL_test() == GLU_FAILURE ) {
    goto TEST_FAILURE ;
  }
  full_tests_run += tests_run ;

  // first one is simple linear algebra tests
  if( UNOPS_test() == GLU_FAILURE ) {
    goto TEST_FAILURE ;
  }
  full_tests_run += tests_run ;

  // log-exponentiate tests
  if( ( test_results = logexp_tests( ) ) != 0 ) {
    goto TEST_FAILURE ;
  } 
  full_tests_run += tests_run ;

  // matrix inversion/determinant tests
  if( ( test_results = inversion_tests( ) ) != 0 ) {
    goto TEST_FAILURE ;
  }
  full_tests_run += tests_run ;

  // geometry tests
  if( geometry_test() == GLU_FAILURE ) {
    goto TEST_FAILURE ;
  }
  full_tests_run += tests_run ;

  // general lattice configuration tests
  if( ( test_results = config_tests( ) ) != 0 ) {
    goto TEST_FAILURE ;
  }
  full_tests_run += tests_run ;

  // if we are successful we go there
  goto TEST_SUCCESS ;

 TEST_FAILURE :
  unstick_GLU( ) ;
  fprintf( stderr , "[GLUnit] Test failure \n%s \n" , test_results ) ;
  return GLU_FAILURE ;
 TEST_SUCCESS :
  unstick_GLU( ) ;
  fprintf( stderr , "[GLUnit] All %d tests passed \n" , full_tests_run ) ;

  return GLU_SUCCESS ;
}
