/**
   @file HB.c
   @brief heatbath code

   Notes :: s0 and s1 are the top row of the SU(2) matrix, this is all we need
            due to the funny nature of links we have to loop over directions and 
	    then volume, this is a bit out of order for our AoS setup but fuck it

   Credit :: Lots of credit to T deGrand whose code (in MILC) I used as the basis 
             for this work and whose code I tested this against.
 */
#include "Mainfile.h"

#include "par_rng.h"    // parallel rng
#include "staples.h"    // all_staples()
#include "SU2_rotate.h" // rotation

// maximum number of heatbath updates
#define NHBMAX (25)

// Creutz heatbath algorithm
static inline GLU_bool
Creutz( double *__restrict d ,
	const double xl ,
	const double NORM ,
	const uint32_t thread ) 
{
  const register double x4 = par_rng_dbl( thread ) ;
  register const double a0 = 1.0 + \
    log( xl + ( 1.0 - xl ) * par_rng_dbl( thread ) ) * NORM ;
  *d = 1.0 - a0 ;
  return ( 1.0 - a0*a0 ) < x4*x4 ;
}

// Kennedy-Pendleton heatbath algorithm
static inline GLU_bool
KP( double *__restrict d ,
    const double xl ,
    const double NORM ,
    const uint32_t thread ) 
{  
  const register double x4 = par_rng_dbl( thread ) ;
  const register double x3 = cos( TWOPI * par_rng_dbl( thread ) ) ;
  *d = -( log( par_rng_dbl( thread ) ) + log( par_rng_dbl( thread ) ) * x3*x3 ) * NORM ;
  return ( 1. - 0.5*(*d) ) < x4*x4 ;
}

// generate SU(2) matrix proportional to boltzmann weight
static void
generate_SU2( GLU_complex *s0 ,
	      GLU_complex *s1 ,
	      const double NORM ,
	      const int thread )
{
  double xl = 1.0 , d = 0.0 ;
  size_t iters = 1 ;

  // function pointer for the algorithm we want
  GLU_bool (*K)( double *__restrict d ,
		 const double xl ,
		 const double NORM ,
		 const uint32_t thread ) ;

  // for small beta the Creutz algorithm is preferred
  switch( NORM > 0.5 ) {
  case 0 : K = KP ; break ;
  case 1 : K = Creutz ; xl = exp( -2./NORM ) ; break ;
  }

  // compute d which is (1-a0) = the leftmost su2 matrix element
  if( K( &d , xl , NORM , thread ) ) {
    // redo the algorithm
    while( K( &d , xl , NORM , thread ) && iters < NHBMAX ) {
      iters++ ;
    }
    // complain if we do too many hits?
    if( iters == NHBMAX ) {
      switch( NORM > 0.5 ) {
      case 0 :
	fprintf( stderr , "[KPHB] %d iterations of KPHB to no avail \n" , NHBMAX ) ;
	break ;
      case 1 :
	fprintf( stderr , "[KPHB] %d iterations of Creutz to no avail \n" , NHBMAX ) ;
	break ;
      }
    }
  }

  // generate SU(2) matrix
  const double a0  = 1.0 - d ;
  const double ar2 = fabs( 1.0 - a0*a0 ) ; 
  const double a3  = sqrt( ar2 ) * ( 2.0 * par_rng_dbl( thread ) - 1.0 ) ;
  const double rho = sqrt( fabs( ar2 - a3*a3 ) ) ;
  const double xr2 = TWOPI * par_rng_dbl( thread ) ;

  const double a2 = rho * sin( xr2 ) ;
  const double a1 = rho * cos( xr2 ) ;
  const double rS0 = creal( *s0 ) , iS0 = cimag( *s0 ) ;
  const double rS1 = creal( *s1 ) , iS1 = cimag( *s1 ) ;

  // unrolled matrix multiply
  *s0 = a0 * rS0 + a3 * iS0 + a2 * rS1 + a1 * iS1 + 
    I * ( -a0 * iS0 + rS0 * a3 - a2 * iS1 + rS1 * a1 ) ;
  *s1 = -a0 * rS1 + a3 * iS1 + a2 * rS0 - a1 * iS0 +
    I * ( -a0 * iS1 - a3 * rS1 + a2 * iS0 + a1 * rS0 ) ;

  return ;
}

// update matrices compute SU(2) subgroup of the staple
// compute heat-bath update and multiply through link
static void
hb( GLU_complex U[ NCNC ] , 
    const GLU_complex staple[ NCNC ] ,
    const double invbeta ,
    const int thread )
{
  GLU_complex s0 , s1 ;
  double scale ;
  size_t i ;
  for( i = 0 ; i < NSU2SUBGROUPS ; i++ ) {
    only_subgroup( &s0 , &s1 , &scale , U , staple , i ) ;
    generate_SU2( &s0 , &s1 , invbeta*scale , thread ) ;
    su2_rotate( U , s0 , s1 , i ) ;
  }
  return ;
}

/*
// playground for testing a macro definition
#if (defined _FOPENMP) && (defined HAVE_OMP_H)
#define STR(x) #x
#define STRINGIFY(x) STR(x) 
#define CONCATENATE(X,Y) X ( Y )
#define GLU_parallel_for(...) _Pragma( STRINGIFY(CONCATENATE(omp parallel for,__VA_ARGS__)))
#else
#define GLU_parallel_for(...) // do nothing
#endif
*/

// perform a heat-bath over the whole lattice
int
hb_lattice( struct site *lat ,
	    struct site *staple ,
	    const double invbeta ,
	    const struct draughtboard db )
{
  size_t i , mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    // compute staples surrounding red links
#pragma omp parallel for private(i)
    for( i = 0 ; i < db.Nred ; i++ ) {
      zero_mat( staple[i].O[mu] ) ;
      all_staples( staple[i].O[mu] , lat , db.red[i] , mu , ND , SM_APE ) ;
    }
    // heat bath all the red links
#pragma omp parallel for private(i)
    for( i = 0 ; i < db.Nred ; i++ ) {   
      hb( lat[ db.red[i] ].O[mu] , staple[i].O[mu] , invbeta , get_GLU_thread() ) ;
    }
    // compute the staples surrounding the black links
#pragma omp parallel for private(i)
    for( i = 0 ; i < db.Nblack ; i++ ) {
      zero_mat( staple[i].O[mu] ) ;
      all_staples( staple[i].O[mu] , lat , db.black[i] , mu , ND , SM_APE ) ;
    }
    // heat bath all black links
#pragma omp parallel for private(i)
    for( i = 0 ; i < db.Nblack ; i++ ) {
      hb( lat[ db.black[i] ].O[mu] , staple[i].O[mu] , invbeta , get_GLU_thread() ) ;
    }
  }
  return GLU_SUCCESS ;
}
