/**
   @brief balls.c
   @brief glueball calculator

   path is (confusingly) ordered from 1

   x,y,z,t = 1,2,3,4
   -x,-y,-z,-t = -1,-2,-3,-4
 */
#include "Mainfile.h"

#include "plan_ffts.h"
#include "SM_wrap.h"

int
link_path( GLU_complex result[ NCNC ] ,
	   const struct site *lat ,
	   const int *path ,
	   const size_t i ,
	   const size_t Npath )
{
  size_t j , k ;
  GLU_complex temp[ NCNC ] GLUalign ;

  size_t mu = (size_t)abs( path[0] ) - 1 ;
  if( path[0] < 0 ) {
    k = lat[i].back[ mu ] ;
    dagger( temp , lat[k].O[mu] ) ;
  } else {
    equiv( temp , lat[i].O[ path[0] - 1 ] ) ;
    k = lat[i].neighbor[ mu ] ;
  } 
  for( j = 1 ; j < Npath ; j++ ) {
    mu = (size_t)abs( path[j] ) - 1 ;
    // logic
    if( path[j] < 0 ) {
      k = lat[k].back[ mu ] ;
      multab_dag_suNC( result , temp , lat[k].O[mu] ) ;
    } else {
      multab_suNC( result , temp , lat[k].O[mu] ) ;
      k = lat[k].neighbor[ mu ] ;
    }
    equiv( temp , result ) ;
    //
  }
  return GLU_SUCCESS ;
}

int
GLUball_lattice( GLU_complex *result ,
		 const struct site *lat ,
		 const int *path ,
		 const size_t Npath )
{
  size_t i ;
  int sum = 0 ;
  for( i = 0 ; i < Npath ; i++ ) {
    sum += path[i] ;
  }
  // check that we are closing the loop
  if( sum != 0 ) {
    return GLU_FAILURE ;
  }
#pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    GLU_complex temp[ NCNC ] GLUalign ;
    link_path( temp , lat , path , i , Npath ) ;
    result[i] += trace( temp ) ;
  }
  return GLU_SUCCESS ;
}

// 
static int
correlator( double *gsp ,
	    const GLU_complex *result )
{
  // init parallel threads, maybe
  if( parallel_ffts( ) == GLU_FAILURE ) {
    fprintf( stderr , "[PAR] Problem with initialising "
	              "the OpenMP FFTW routines \n" ) ;
    // should clean up the memory here
    return GLU_FAILURE ;
  }

  // FFTW routines
  GLU_complex *out = fftw_malloc( LVOLUME * sizeof( GLU_complex ) ) ;
  GLU_complex *in = fftw_malloc( LVOLUME * sizeof( GLU_complex ) ) ;

  // create some plans
  fftw_plan forward , backward ;
  small_create_plans_DFT( &forward , &backward , in , out , Latt.dims , ND ) ;

  size_t i ;
  for( i = 0 ; i < LVOLUME ; i++ ) {
    in[i] = result[i] ;
  }

  fftw_execute( forward ) ;

  // convolve
  for( i = 0 ; i < LVOLUME ; i++ ) {
    out[i] *= conj( out[i] ) ;
  }

  fftw_execute( backward ) ;

  size_t t ;
  for( t = 0 ; t < Latt.dims[ ND-1 ] ; t++ ) {
    gsp[t] = 0.0 ;
    for( i = 0 ; i < LCU ; i++ ) {
      gsp[t] += creal( in[ i + LCU*t ] ) ;
    }
    //gsp[t] /= LVOLUME ;
  }

  // free the FFTs
  fftw_destroy_plan( backward ) ;
  fftw_destroy_plan( forward ) ;
  fftw_cleanup( ) ;
#ifdef OMP_FFTW
  fftw_cleanup_threads( ) ;
#endif
  fftw_free( out ) ;  
  fftw_free( in ) ;

  return GLU_SUCCESS ;
}

int
GLU_balls( struct site *lat ,
	   struct sm_info SMINFO )
{
  SM_wrap_struct( lat , SMINFO ) ;

  GLU_complex *result = NULL ;
  GLU_malloc( (void**)&result , 16 , LVOLUME*sizeof(GLU_complex) ) ;

  // create a plaquette correlator
  int path[ 4 ] ;

  size_t x ;
  for( x = 0 ; x < LVOLUME ; x++ ) {
    result[x] = 0.0 ;
  }

  // loop variations
  size_t i , j ;
  for( i = 0 ; i < ND-1 ; i++ ) {
    for( j = i+1 ; j < ND-1 ; j++ ) {
      // top right plaq
      path[0] =  i + 1 ; path[1] =  j + 1 ; 
      path[2] = -i - 1 ; path[3] = -j - 1 ;
      GLUball_lattice( result , lat , path , 4 ) ;
    }
  }

  // check average plaquette
  double sum = 0.0 ;
  for( i = 0 ; i < LVOLUME ; i++ ) {
    sum += creal( result[i] ) ;
  }
  printf( "Plaq :: %1.15f \n" , sum/(LVOLUME*18) ) ;

  // compute the correlator
  double C[ Latt.dims[ ND-1 ] ] ;
  correlator( C , result ) ;
      
  //
  size_t t ;
  for( t = 0 ; t < Latt.dims[ ND-1 ]/2 ; t++ ) {
    printf( "%zu %e \n" , t , C[t] ) ;
  }

  free( result ) ;

  return GLU_SUCCESS ;
}
