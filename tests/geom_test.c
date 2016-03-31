/**
   @file geom_test.c
   @brief test our geometry functions
 */
#include "Mainfile.h"

#include "geometry.h"
#include "minunit.h"

static char*
gen_site_test( void )
{
  // test that {0,0,...,0} == 0
  int x[ ND ] , mu ;
  for( mu = 0 ; mu < ND ; mu++ ) x[mu] = 0 ;
  mu_assert( "[GEOM] gen_site broken" , gen_site(x)==0 ) ;
  // test that {L-1,L-1,L-1...,L-1} == LVOLUME-1
  for( mu = 0 ; mu < ND ; mu++ ) x[mu] = Latt.dims[mu]-1 ;
  mu_assert( "[GEOM] gen_site broken" , 
	     ( gen_site(x) == (LVOLUME - 1) ) ) ;
  return NULL ;
}

// test that (-l/2,-l/2,-l/2,-l/2) maps to site zero
static char*
get_site_2piBZ_test( void )
{
  int x[ ND ] , mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    x[mu] = -(int)Latt.dims[mu]/2 ;
  }
  mu_assert( "[GEOM] get_site_2piBZ broken" , 
	     ( get_site_2piBZ(x,ND) == 0 ) ) ;
  return NULL ;
}

// test (0,0,...,0) maps to zero
static char*
get_site_pipiBZ_test( void )
{
  int x[ ND ] , mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    x[mu] = 0 ;
  }
  mu_assert( "[GEOM] get_site_pipiBZ broken" , 
	     ( get_site_pipiBZ(x,ND) == 0 ) ) ;
  return NULL ;
}

// test that site 0 is (-l/2,-l/2..)
static char*
get_mom_pipi_test( void )
{
  int x[ ND ] , mu ;
  get_mom_pipi( x , 0 , ND ) ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    mu_assert( "[GEOM] get_mom_pipi broken" , 
	       ( x[mu] == -(int)Latt.dims[mu]/2 ) ) ;
  }
  return NULL ;
}

// test that site 0 is (0,0,0,...)
static char*
get_mom_2piBZ_test( void )
{
  int x[ ND ] , mu ;
  get_mom_2piBZ( x , 0 , ND ) ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    mu_assert( "[GEOM] get_mom_2piBZ broken" , 
	       ( x[mu] == 0 ) ) ;
  }
  return NULL ;
}

// test that site 0 is (-l/2,..)
static char*
get_mom_TwoPI_mpipi_momconv_test( void )
{
  int x[ ND ] , mu ;
  TwoPI_mpipi_momconv( x , 0 , ND ) ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    mu_assert( "[GEOM] get_mom_TwoPI_mpipi_momconv broken" , 
	       ( x[mu] == 0 ) ) ;
  }
  return NULL ;
}

// test for the HiRep geom has (t,x,y,z) z runs fastest
static char*
get_mom_2piBZ_hirep2_test( void )
{
  int x[ ND ] , mu ;
  // site 1 should be (0,0,1,0) in our order
  get_mom_2piBZ_hirep2( x , 1 ) ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    if( mu == ND-2 ) {
      mu_assert( "[GEOM] get_mom_2piBZ_hirep2 broken" , 
		 ( x[mu] == 1 ) ) ;
    } else {
      mu_assert( "[GEOM] get_mom_2piBZ_hirep2 broken" , 
		 ( x[mu] == 0 ) ) ;
    }
  }
  return NULL ;
}

// test should be zero
static char*
compute_p_test( void )
{
  GLU_real p[ ND ] ;
  int n[ ND ] , mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    n[ mu ] = 0 ;
  }
  compute_p( p , n , ND ) ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    mu_assert( "[GEOM] compute_p broken" , 
	       ( fabs( p[mu] ) < PREC_TOL ) ) ;
  }
  return NULL ;
}

// test the 0 site has 1 because of special case
static char*
gen_p_sq_test( void )
{
  // momentum at site zero on 2piBZ (fftw) lattice should be zero
  mu_assert( "[GEOM] gen_psq broken" , 
	     ( fabs( gen_p_sq( 0 , ND ) - 1 ) < PREC_TOL ) ) ;
  return NULL ;
}

// test U1 momentum condition at site 0 should get psq = 1 && flag =  1
static char*
gen_p_sq_feyn_test( void )
{
  int flag ;
  GLU_real psq = gen_p_sq_feyn( 0 , &flag ) ;
  mu_assert( "[GEOM] gen_p_sq_feyn broken psq" , 
	     ( fabs( psq - 1 ) < PREC_TOL ) ) ;
  mu_assert( "[GEOM] gen_p_sq_feyn broken flag" , 
	     ( flag == 1 ) ) ;
  return NULL ;
}

// lattice momentum is 0 at site midpoint in the -Pi->Pi lattice
static char*
gen_get_p_test( void )
{
  GLU_real p[ ND ] ;
  int mu , midpoint = Latt.dims[ ND-1 ]/2 ;
  for( mu = ND-1 ; mu > 0 ; mu-- ) {
    midpoint = Latt.dims[mu-1] * midpoint + Latt.dims[ mu - 1 ]/2 ; 
  }
  gen_get_p( p , midpoint , ND ) ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    mu_assert( "[GEOM] gen_get_p broken" , 
	       ( fabs( p[mu] ) < PREC_TOL ) ) ;
  }
  return NULL ;
}

static char*
get_vec_from_origin_test( void )
{
  // test that the site at (1,1,1,1) is (1,1 .. ) from origin
  int x[ ND ] , mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    x[ mu ] = 1 ;
  }
  get_vec_from_origin( x , gen_site( x ) , ND ) ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    mu_assert( "[GEOM] get_vec_from_origin broken" , x[mu] == 1 ) ;
  }
  return NULL ;
}

// compute rsq site 1 away from everything has rsq ND
static char*
compute_rsq_test( void )
{
  int x[ ND ] , mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    x[ mu ] = 1 ;
  }
  mu_assert( "[GEOM] compute_rsq broken" , 
	     compute_rsq( gen_site( x ) , ND ) == ND ) ;
  return NULL ;
}

// test by moving L away from origin and checking is zero
static char*
compute_spacing_test( void )
{
  int x[ ND ] , mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    x[ mu ] = Latt.dims[ mu ]  ;
  }
  mu_assert( "[GEOM] compute_spacing broken" , 
	     compute_spacing( x , 0 , ND ) == 0 ) ;
  return NULL ;
}

// shift in L directions and end up in same place
static char*
gen_shift_test( void )
{
  size_t mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    // shift forward Latt.dims times puts you back where you were
    size_t nu , i = 0 ;
    for( nu = 0 ; nu < Latt.dims[mu] ; nu++ ) {
      i = gen_shift( i , mu ) ;
    }
    mu_assert( "[GEOM] gen_shift broken forwards" , 
	       i == 0 ) ;
    i = 0 ;
    for( nu = 0 ; nu < Latt.dims[mu] ; nu++ ) {
      i = gen_shift( i , -mu-1 ) ;
    }
    mu_assert( "[GEOM] gen_shift broken backwards" , 
	       i == 0 ) ;
  }
  return NULL ;
}

// geometry tests
static char *
geom_test( void )
{
  mu_run_test( gen_site_test ) ;
  mu_run_test( get_site_2piBZ_test ) ;
  mu_run_test( get_site_pipiBZ_test ) ;
  mu_run_test( get_mom_pipi_test ) ;
  mu_run_test( get_mom_2piBZ_test ) ;
  mu_run_test( get_mom_TwoPI_mpipi_momconv_test ) ;
  mu_run_test( get_mom_2piBZ_hirep2_test ) ;
  mu_run_test( compute_p_test ) ;
  mu_run_test( gen_p_sq_test ) ;
  mu_run_test( gen_p_sq_feyn_test ) ;
  mu_run_test( gen_get_p_test ) ;
  mu_run_test( get_vec_from_origin_test ) ;
  mu_run_test( compute_rsq_test ) ;
  mu_run_test( compute_spacing_test ) ;
  mu_run_test( gen_shift_test ) ;

  return NULL ;
}

// test the geometries
int
geometry_test( void )
{
  // init to zero again
  tests_run = tests_fail = 0 ;

  // initial gamma setup and test
  char *geomres = geom_test( ) ;

  if( tests_fail == 0 ) {
    printf( "[GEOMETRY UNIT] all %d tests passed\n\n" ,
	    tests_run ) ;
    return GLU_SUCCESS ;
  } else {
    printf( "%s \n" , geomres ) ;
    printf( "[GEOMETRY UNIT] %d out of %d tests failed\n\n" , 
	    tests_fail , tests_run ) ;
    return GLU_FAILURE ;
  }
}
