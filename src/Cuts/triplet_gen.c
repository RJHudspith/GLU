/*
    Copyright 2013 Renwick James Hudspith

    This file (triplet_gen.c) is part of GLU.

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
   @file triplet_gen.c
   @brief compute the momentum-conserving triplets

   @warning I have added the generic code for this using tail recursion. It is now ND-generic for what it is worth.
 */
#include "Mainfile.h"

#include "cut_routines.h" // for simorb ratios
#include "geometry.h"     // momentum calcs and what have you

// get a position of a possible momentum
static int
get_posit( const int *p ,
	   const int *__restrict *__restrict momentum ,
	   const int nmom )
{
  size_t check , nu , flag = 0 ;
  // loop momentum list
  for( check = 0 ; check < nmom ; check ++ ) {
    for( nu = 0 ; nu < ND ; nu++ ) {
      if( momentum[ check ][nu] != p[nu] ) {
	flag = 1 ;
	break ;
      }
    }
    if( flag == 0 ) { return check ; }
  }
  printf( "404 ! MOMENTUM NOT FOUND :( NMOM :: %d ( " , nmom ) ;
  for( nu = 0 ; nu < ND ; nu++ ) {
    printf( "%d " , p[nu] ) ;
  }
  printf( ")\n" ) ;
  return GLU_FAILURE ;
}

///// EXPERIMENTAL CUTS /////
#ifdef CUTS_EXPERIMENT
// put in a cutting routine ...
static size_t count_zeros( const int p[ ND ] ) 
{
  size_t mu , count = 0 ;
  for( mu = 0 ; mu < ND ; mu++ ) { if( p[mu] == 0 ) { count++ ; } } return count ;
}
// have a look to see if one of the momenta are over a value too .
static size_t over_val( const int p[ ND ] , const int val ) 
{
  size_t mu ;
  for( mu = 0 ; mu < ND ; mu++ ) { if( p[mu] > ( val*(double)Latt.dims[mu]/(double)Latt.dims[ND-1] ) ) return GLU_FAILURE ; }
  return GLU_SUCCESS ;   
}
#define NZEROS 1

static int
momcut( int p1[ ND ] , int p2[ ND ] )
{
  size_t mu ;
  int p3[ ND ] ;
  for( mu = 0 ; mu < ND ; mu++ ) { p3[mu] = -( p1[mu] + p2[mu] ) ; }
  if( count_zeros( p1 ) > NZEROS || over_val( p1 , Latt.dims[ND-1]/2 ) == GLU_FAILURE ) { return GLU_FAILURE ; }
  if( count_zeros( p2 ) > NZEROS || over_val( p2 , Latt.dims[ND-1]/2 ) == GLU_FAILURE ) { return GLU_FAILURE ; }
  if( count_zeros( p3 ) > NZEROS || over_val( p3 , Latt.dims[ND-1]/2 ) == GLU_FAILURE ) { return GLU_FAILURE ; }
  return GLU_SUCCESS ;
}
#undef NZEROS
#endif

// return the lattice momentum from the points in the BZ
static inline double
get_psq( const int *p )
{ 
  size_t nu ;
  register double psq = 0.0 ;
  for( nu = 0 ; nu < ND ; nu++ ) {
    const double cache = p[ nu ] * Latt.twiddles[ nu ] ;
    psq += cache * cache ; 
  }
  return psq ;
}

// generic delta function
static double **d ;
static void
init_delta( void )
{
  d = (double**)malloc( ND * sizeof( double* ) ) ;
  size_t i , j ;
  for( i = 0 ; i < ND ; i++ ) {
    d[i] = (double*)malloc( ND * sizeof( double ) ) ;
    for( j = 0 ; j < ND ; j++ ) {
      d[i][j] = ( i!=j ) ? 0.0 : 1.0 ;
    }
  }
  return ;
}

static void
free_delta( void )
{
  size_t i ;
  for( i = 0 ; i < ND ; i++ ) {
    free( d[i] ) ;
  }
  free( d ) ;
  return ;
}

//  Recurse the second looking for again equivalent psq's and finally the 
//  third triplet. If we are setting the values of the triplet we do a check
//  with the lattice momentum.

//  It turns out I was stack-smashing pretty badly here, I think it is fixed.
//  Actually it looks like the scoping in the openmp'd bit was killing us
static int
recurse_p2( size_t *__restrict *__restrict triplet ,
	    int *__restrict *__restrict momentum , 
	    int *p1 , 
	    int *p2 , 
	    const size_t nn , 
	    const size_t r_nn , 
	    const size_t nmom ,
	    const size_t count , 
	    const size_t check , 
	    const int FIRST_PASS )
{ 
  size_t mu , sum = 0 ; 
  for( mu = 0 ; mu < ND ; mu++ ) {
    sum += p2[mu] * p2[mu] ;
  }

  // if the sum is in the orbit perform one last check
  if( sum == nn ) {

    int test = 0 ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      test += ( p1[mu] + p2[mu] ) * ( p1[mu] + p2[mu] ) ;
    }

    // if we have a final triplet we equate all the triplets
    if( test == nn ) {

	// if we are not doing a first pass, we actually compute something
	if( FIRST_PASS != GLU_TRUE ) {

	  // compute the asymmetry-adjusted momenta ...
	  int *adp3 = malloc( ND * sizeof( int ) ) ;
	  int *adp2 = malloc( ND * sizeof( int ) ) ;
	  int *adp1 = malloc( ND * sizeof( int ) ) ;
	  
	  for( mu = 0 ; mu < ND ; mu++ ) {	 
	    adp1[mu] = (int)rats[mu] * ( p1[mu] ) ;
	    adp2[mu] = (int)rats[mu] * ( p2[mu] ) ;
	    adp3[mu] = -(int)rats[mu] * ( p1[mu] + p2[mu] ) ;
	  }   
	  
	  const int trip_idx = count + check ;

	  triplet[ trip_idx ][0] = get_posit( adp1 , (const int**)momentum , nmom ) ;
	  triplet[ trip_idx ][1] = get_posit( adp2 , (const int**)momentum , nmom ) ;
	  triplet[ trip_idx ][2] = get_posit( adp3 , (const int**)momentum , nmom ) ;

	  double check0 , check1 , check2 ;
	  check0 = get_psq( momentum[ triplet[ trip_idx ][0] ] ) ;
	  check1 = get_psq( momentum[ triplet[ trip_idx ][1] ] ) ;
	  check2 = get_psq( momentum[ triplet[ trip_idx ][2] ] ) ;    

	  // and free the asymmetry corrections ...
	  free( adp3 ) ; free( adp2 ) ; free( adp1 ) ;

	  if( fabs( check0 - check1 ) > PREC_TOL || 
	      fabs( check0 - check2 ) > PREC_TOL || 
	      fabs( check1 - check2 ) > PREC_TOL ) {
	    printf( "Non-conservation of momenta !! \n" ) ;
	    return GLU_FAILURE ;
	  }
	}
	// increment our check
	check ++ ;
      }
  }

  // loop the smallest first
  for( mu = ND - 1 ; mu != 0 ; mu-- ) {
    if( p2[mu] < r_nn ) {
      p2[mu]++ ;
      return recurse_p2( triplet , momentum , p1 , p2 , 
			 nn , r_nn , nmom , 
			 count , check , FIRST_PASS ) ;
    } else {
      p2[mu] = -r_nn ;
    }
  }

  return check ;
}

// recurse the first level and look for equivalent psq's
static int
recurse_p( size_t *__restrict *__restrict triplet ,
	   int *__restrict *__restrict momentum , 
	   int *p1 , 
	   const size_t nn , 
	   const size_t r_nn ,
	   const size_t nmom , 
	   const size_t count ,
	   const int FIRST_PASS )
{ 
  size_t mu , sum = 0 ; 
  for( mu = 0 ; mu < ND ; mu++ ) {
    sum += p1[mu] * p1[mu] ;
  }
  size_t check = 0 ;
  // if we have a hit we recurse down a level
  if( sum == nn ) {
    //int p2[ ND ] ;
    int *p2 = malloc( ND * sizeof( int ) ) ; 
    for( mu = 0 ; mu < ND ; mu++ ) {
      p2[ mu ] = -r_nn ; 
    }  
    count += recurse_p2( triplet , momentum , p1 , p2 , 
			 nn , r_nn , nmom ,  
			 count , check , FIRST_PASS ) ;
    free( p2 ) ;
  }
  // loop the smallest first
  for( mu = ND-1 ; mu != 0 ; mu-- ) {
    if( p1[mu] < r_nn ) {
      p1[mu]++ ;
      return recurse_p( triplet , momentum , p1 , 
			nn , r_nn , nmom , 
			count , FIRST_PASS ) ;
    } else {
      p1[mu] = -r_nn ;
    }
  }
  return count ;
}

/**
   @fn static void get_triplet( int *__restrict *__restrict triplet , int *__restrict *__restrict momentum , const int nnmax , const int nmom )
   \brief computes the triplet and checks for momentum conservation
   @param triplet :: is of the form triplet[ triplet_index ][ 0,1,2 ] as there are three of them
   @param momentum :: a list of the momentum in integer rep of the form momentum[ index ][ 0,1,2 .. ND ]
   @param nnmax :: maximum "p^2"
   for weird non-integer anisotropies this means we can only live on a very small subset
   of the lattice. This is a massive problem for us.

   @return #GLU_SUCCESS or #GLU_FAILURE
 **/
static void
get_triplet( size_t *__restrict *__restrict triplet , 
	     int *__restrict *__restrict momentum , 
	     const size_t nnmax , 
	     const size_t nmom )
{
  size_t nn ;
  simorb_ratios( ND ) ;
  for( nn = 0 ; nn < ND ; nn++ ) { rats[ nn ] = 1. / rats[nn] ; }

  size_t count = 0 ;
  // loop the nn's
  for( nn = 0 ; nn < nnmax ; nn += 2 ) {
    const int root_nn = (int)sqrt( nn ) ;

    //int p[ ND ] , mu ;
    int *p = malloc( ND * sizeof( int ) ) , mu ;
    for( mu = 0 ; mu < ND ; mu ++ ) {
      p[mu] = ( fabs( rats[mu] - (int)rats[mu] ) < PREC_TOL ) ? -root_nn : 0.0 ;
    }

    count = recurse_p( triplet , momentum , p , 
		       nn , root_nn , nmom , 
		       count , GLU_FALSE ) ;

    printf( "Checking mom-conservation orbit :: %zu \n" , nn ) ; 
    free( p ) ;
  }
  return ;
}

/**
   @fn static void get_trip( size_t *__restrict trip , const size_t nnmax )
   @brief calculates the size of each triplet and checks once again for conservation of lattice momentum.
   @param trip :: the momentum list after cutting
   @param nnmax :: maximum \f$ p^2 \f$

   computes the size of the triplet
   @return #GLU_SUCCESS or #GLU_FAILURE
 **/
static void
get_trip( size_t *__restrict trip , 
	  const size_t nnmax )
{
  // unused variables ...
  size_t **triplet = malloc( sizeof( int* ) ) ;
  int **momentum = malloc( sizeof( int* ) ) ;
  size_t nmom = 0 ;

  // feel dirty just doing this
  size_t nn ;
  for( nn = 0 ; nn < 1 ; nn++ ) {
    triplet[nn] = (size_t*) malloc( sizeof( size_t ) ) ;
    momentum[nn] = (int*) malloc( sizeof( int ) ) ;
  }

  // compute the asymmetry factors 
  simorb_ratios( ND ) ;
  for( nn = 0 ; nn < ND ; nn++ ) { rats[ nn ] = 1. / rats[nn] ; }

  // finally loop our triples
  size_t count = 0 ;
  // loop the nn's
  for( nn = 0 ; nn < nnmax ; nn += 2 ) {

    const size_t root_nn = (int)sqrt( nn ) ;
    //int p[ ND ] , mu ;
    int *p = malloc( ND * sizeof( int ) ) , mu ;
    for( mu = 0 ; mu < ND ; mu ++ ) {
      p[mu] = ( fabs( rats[mu] - (int)rats[mu] ) < PREC_TOL ) ? -root_nn : 0.0 ;
    }

    size_t check = count ;
    size_t this = recurse_p( triplet , momentum , p , 
			     nn , root_nn , nmom , 
			     count , GLU_TRUE ) ;
    
    free( p ) ;
    count = this ;
    trip[nn/2] = count - check ;
    printf( "[CUTS] %zu :: %zu \n" , nn/2 , trip[nn/2] ) ;
  }

  // free this crap
  for( nn = 0 ; nn < 1 ; nn++ ) {
    free( triplet[nn] ) ;
    free( momentum[nn] ) ; 
  }
  free( triplet ) ; free( momentum ) ;
  return ;
}

// returns the projector for the symmetric def (BOUCAUD)
// written so that whatever mom-def you choose( sin )
// it should work
#ifndef PROJ_GRACEY

static double
get_proj( const double p1[ ND ] , 
	  const double p2[ ND ] , 
	  const double p3[ ND ] ;
	  const size_t mu , 
	  const size_t nu , 
	  const size_t rho )
{
  size_t mup , nup ,rhop ;
  double psq = 0. ;
  for( mup = 0 ; mup < ND ; mup++ )  {
    const double mom = p1[mup] ;
    psq += mom * mom ;
  }

  // control for NaN's in the zero_mat case
  if( psq < PREC_TOL ) {
    psq = 1.0 ;
  } else {
    psq = 1.0 / psq ;
  } 

  double proj = 0. , temp ;
  // sum over the primes
  for( mup = 0 ; mup < ND ; mup ++ ) {
    for( nup = 0 ; nup < ND ; nup ++ ) {
      for( rhop = 0 ; rhop < ND ; rhop ++ ) {
	// implicitly include the delta's
	temp = 0. ;
	if( mup == nup ) {
	  temp = ( p1[rhop] - p2[rhop] ) ;
	} if( mup == rhop ) {
	  temp += ( p3[nup] - p1[nup] ) ;
	} if( nup == rhop ) {
	  temp += ( p2[mup] - p3[mup] ) ;
	}
	// fairly pointless adding the zero matrix all the time
	if( fabs(temp) < PREC_TOL ) {
	} else {
	  temp *= ( ( d[mup][mu] - p1[mup] * p1[mu] * psq ) *	\
		    ( d[nup][nu] - p2[nup] * p2[nu] * psq ) *	\
		    ( d[rhop][rho] - p3[rhop] * p3[rho] * psq ) ) ;
	  proj += temp ;
	}

      } } }

  // add this last little bit in ....
  proj += (p1[rho] - p2[rho]) * ( p3[nu] - p1[nu] ) * ( p2[mu] - p3[mu] ) * 0.5 * psq ;

  return proj / ( 2.0 * ( ND - 1 ) * ( ND - 1 ) ) ;
}

#endif

// These are the 14 projectors defined in appendix 
// of the paper http://arxiv.org/abs/1108.4806 A.1
// care and due dilligence was required here as these projectors
// were in Minkowski space and there was an overall factor of psq missing
// I have really checked these ( using mathematica ) and am happy with them
#define p1 d[mu][nu] * p[sigma] 
#define p2 d[nu][sigma] * p[mu] 
#define p3 d[sigma][mu] * p[nu] 
#define p4 d[mu][nu] * q[sigma] 
#define p5 d[nu][sigma] * q[mu] 
#define p6 d[sigma][mu] * q[nu] 
#define p7 psq * p[mu] * p[nu] * p[sigma] 
#define p8 psq * p[mu] * p[nu] * q[sigma] 
#define p9 psq * p[mu] * q[nu] * p[sigma] 
#define p10 psq * q[mu] * p[nu] * p[sigma] 
#define p11 psq * p[mu] * q[nu] * q[sigma] 
#define p12 psq * q[mu] * p[nu] * q[sigma] 
#define p13 psq * q[mu] * q[nu] * p[sigma] 
#define p14 psq * q[mu] * q[nu] * q[sigma] 
#define DENOM 1.0 / ( 27.0 * ( ND - 2 ) )

/// This is the projector defined from the gracey paper...
static inline double
get_proj_gracey( const double p[ ND ] , 
		 const double q[ ND ] ,
		 const size_t mu , 
		 const size_t nu , 
		 const size_t sigma )
{
  double psq = 0. ;
  size_t mup ;
  for( mup = 0 ; mup < ND ; mup++ ) {
    psq += p[mup] * p[mup] ;
  }
  psq = ( psq < PREC_TOL ) ? 1.0 : 1.0/psq ;

  // Vaguely mathematica-generated output...
#if PROJ_GRACEY == 1
  return DENOM*(-6.*(6.*p1 - 4.*p10 - 2.*p11 - 2.*p12 - 8.*p13 - 4.*p14 + 3.*p4 - 4.*(2.*p7 + p8 + p9))) ;
#elif PROJ_GRACEY == 2
  return DENOM*(6.*(4.*p10 + 8.*p11 + 2.*p12 + 2.*p13 + 4.*p14 - 6.*p2 - 3.*p5 + 4.*(2.*p7 + p8 + p9))) ;
#elif PROJ_GRACEY == 3
  return DENOM*(6.*(4.*p10 + 2.*p11 + 8.*p12 + 2.*p13 + 4.*p14 - 6.*p3 - 3.*p6 + 4.*(2.*p7 + p8 + p9))) ;
#elif PROJ_GRACEY == 4
  return DENOM*(-6.*(3.*p1 - 2.*(p10 + 2.*p11 + 2.*p12 + 2.*p13 + 4.*p14 - 3.*p4 + 2.*p7 + 4.*p8 + p9))) ;
#elif PROJ_GRACEY == 5
  return DENOM*(6.*(8.*p10 + 4.*p11 + 4.*p12 + 4.*p13 + 8.*p14 - 3.*p2 + 2.*(-3.*p5 + 2.*p7 + p8 + p9))) ;
#elif PROJ_GRACEY == 6
  return DENOM*(6.*(2.*p10 + 4.*p11 + 4.*p12 + 4.*p13 + 8.*p14 - 3.*p3 - 6.*p6 + 4.*p7 + 2.*p8 + 8.*p9)) ;
#elif PROJ_GRACEY == 7
  return DENOM*(8.*(6.*p1 - 4.*(1 + ND).*p10 - 8.*p11 - 8.*p12 - 8.*p13 - 10.*p14 + 6.*p2 + 6.*p3 + 3.*p4 + 3.*p5 + 3.*p6 - 4.*(2.*p7 + p8 + p9) - 
	    ND.*(2.*p11 + 2.*p12 + 2.*p13 + p14 + 4.*(2.*p7 + p8 + p9)))) ;
#elif PROJ_GRACEY == 8
  return DENOM*(4.*(6.*p1 - 4.*(1 + ND).*p10 - 2.*p12 - 8.*p13 - 16.*p14 + 6.*p2 + 6.*p3 + 12.*p4 + 3.*p5 + 3.*p6 - 8.*p7 + 8.*p8 - 4.*p9 - 
	     2.*(p11 + 4.*ND.*p11 + ND.*(4.*p12 + p13 + 2.*(p14 + 2.*p7 + 4.*p8 + p9))))) ;
#elif PROJ_GRACEY == 9
  return DENOM*(4.*(6.*p1 - 4.*(1 + ND).*p10 - 8.*p12 - 2.*p13 - 16.*p14 + 6.*p2 + 6.*p3 + 3.*p4 + 3.*p5 + 12.*p6 - 8.*p7 - 4.*p8 + 8.*p9 - 
	    2.*(p11 + 4.*ND.*p11 + ND.*(p12 + 2.*(2.*p13 + p14 + 2.*p7 + p8 + 4.*p9))))) ;
#elif PROJ_GRACEY == 10
  return DENOM*(4.*(6.*p1 + (8 - 16.*ND).*p10 - 2.*p12 - 2.*p13 - 16.*p14 + 6.*p2 + 6.*p3 + 3.*p4 + 12.*p5 + 3.*p6 - 4.*(2.*p7 + p8 + p9) - 
	     2.*((4 + ND).*p11 + 2.*ND.*(2.*p12 + 2.*p13 + p14 + 2.*p7 + p8 + p9)))) ;
#elif PROJ_GRACEY == 11
  return DENOM*(4.*(3.*p1 + 8.*p11 - 4.*p12 - 4.*p13 - 8.*p14 + 12.*p2 + 3.*p3 + 6.*p4 + 6.*p5 + 6.*p6 - 2.*(8.*p7 + p8 + p9) - 
	    2.*((4 + ND).*p10 + 2.*ND.*(4.*p11 + p12 + p13 + 2.*p14 + p7 + 2.*(p8 + p9))))) ;
#elif PROJ_GRACEY == 12
  return DENOM*(4.*(3.*p1 - 4.*p11 + 8.*p12 - 4.*p13 - 8.*p14 + 3.*p2 + 12.*p3 + 6.*p4 + 6.*p5 + 6.*p6 - 16.*p7 - 2.*p8 - 8.*p9 - 
	     2.*(p10 + 4.*ND.*p10 + ND.*(2.*(p11 + 4.*p12 + p13 + 2.*p14 + p7 + 2.*p8) + p9)))) ;
#elif PROJ_GRACEY == 13
  return DENOM*(4.*(12.*p1 - 4.*p11 - 4.*p12 + 8.*p13 - 8.*p14 + 3.*p2 + 3.*p3 + 6.*p4 + 6.*p5 + 6.*p6 - 2.*(8.*p7 + 4.*p8 + p9) - 
	    2.*(p10 + 4.*ND.*p10 + ND.*(2.*p11 + 2.*p12 + 8.*p13 + 4.*p14 + 2.*p7 + p8 + 4.*p9)))) ;
#elif PROJ_GRACEY == 14
  return DENOM*(8.*(3.*p1 - 2.*(4 + ND).*p10 - 4.*(1 + ND).*p11 - 4.*(1 + ND).*p12 - 4.*(1 + ND).*p13 - 
	     8.*(1 + ND).*p14 + 3.*p2 + 3.*p3 + 6.*p4 + 6.*p5 + 6.*p6 - (10 + ND).*p7 - 2.*(4 + ND).*p8 - 2.*(4 + ND).*p9)) ;
#else
  return DENOM*(3.*(-2.*p1 + 4.*p10 + 2.*p11 - 4.*p12 + 2.*p13 - 4.*p14 - 5.*p2 + 4.*p3 + 2.*p4 - 4.*p5 + 5.*p6 + 4.*p7 - 2.*(p8 + p9)))/2. ;
#endif
  return GLU_FAILURE ;
}

// computes the projector is a [i][ND*ND*ND] array
static void
compute_projector( size_t *__restrict *__restrict triplet , 
		   double *__restrict *__restrict proj , 
		   int *__restrict *__restrict momentum , 
		   const size_t count )
{
  // precompute the projector
  size_t i ;
  #pragma omp parallel for private(i)
  for( i = 0 ; i < count ; i++ ) {
    double mom[ 3 ][ ND ] ;
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu ++ ) {
      #ifdef PSQ_MOM
      // PSQ variant
      mom[ 0 ][ mu ] = momentum[ triplet[ i ][ 0 ] ][ mu ] * Latt.twiddles[ mu ] ;
      mom[ 1 ][ mu ] = momentum[ triplet[ i ][ 1 ] ][ mu ] * Latt.twiddles[ mu ] ;
      mom[ 2 ][ mu ] = momentum[ triplet[ i ][ 2 ] ][ mu ] * Latt.twiddles[ mu ] ;
      #else
      // SIN variant
      mom[ 0 ][ mu ] = 2.0 * sin( momentum[ triplet[ i ][ 0 ] ][ mu ] * Latt.twiddles[ mu ] ) ;
      mom[ 1 ][ mu ] = 2.0 * sin( momentum[ triplet[ i ][ 1 ] ][ mu ] * Latt.twiddles[ mu ] ) ;
      mom[ 2 ][ mu ] = 2.0 * sin( momentum[ triplet[ i ][ 2 ] ][ mu ] * Latt.twiddles[ mu ] ) ;
      #endif
    }

    size_t nu , rho ;
    for( mu = 0 ; mu < ND ; mu ++ ) {
      for( nu = 0 ; nu < ND ; nu ++ ) {
	for( rho = 0 ; rho < ND ; rho ++ ) {
	  const size_t z = rho + ND * ( nu + ND * mu ) ; 
	  // counter for our projector
          #ifdef PROJ_GRACEY
	  proj[i][z] = get_proj_gracey( mom[ 0 ] ,\
					mom[ 1 ] ,	\
					mu , nu , rho ) ;
          #else
	  proj[i][z] = get_proj( mom[ 0 ] ,\
				 mom[ 1 ] ,\
				 mom[ 2 ] ,\
				 mu , nu , rho ) ;
          #endif
	} } }
  }
  return ;
}

////// Read the triplet //////////
int
read_trip( size_t *__restrict trip ,
	   const int nnmax )
{
#ifdef NOT_CONDOR_MODE
  const size_t nn = nnmax / 2 ;
  char str[256] ;
  sprintf( str , "%s/Local/Moments/TRIP_%d.config" ,
	   HAVE_PREFIX ,
	   nnmax ) ;
  FILE *tripfile = fopen( str , "rb" ) ;
  int flag = 0 ;

  //force it to open ->create a file if needed
  if( tripfile == NULL ) {
    flag = 1 ;
  } else {
    fclose( tripfile ) ; 
  }
  // calculate and write out to a file
  if( unlikely( flag == 1 ) ) {

    // get the triplet
    get_trip( trip , nnmax ) ;

    printf("[CUTS] Storing Trip list @@@ ...\n%s\n",str) ;
    // write to a file 
    FILE *tripfile2 = fopen( str , "wb" ) ;
    fwrite( trip , sizeof( size_t ) , nn , tripfile2 ) ;
    fclose( tripfile2 ) ;
  }
  // otherwise we read it in ...
  tripfile = fopen( str , "rb" ) ;
  if( fread( trip , sizeof( size_t ) , nn , tripfile ) != nn ) return GLU_FAILURE ; 
  fclose( tripfile) ;

#else
  get_trip( trip , nnmax ) ;
#endif  
  return GLU_SUCCESS ;
}


// Reader of files and computer of the projector ...
int
read_triplet_and_proj( size_t *__restrict *__restrict triplet , 
		       double *__restrict *__restrict proj , 
		       int *__restrict *__restrict momentum , 
		       const size_t nnmax , 
		       const size_t count ,
		       const size_t nmom )
{
#ifdef NOT_CONDOR_MODE

  // look for a configuration file that has the available triplets in ...
  char str[496] ;
  sprintf( str , "%s/Local/Moments/TRIPNP.PROJ%d." ,
	   HAVE_PREFIX ,
	   #ifndef PROJ_GRACEY
	   -1 
	   #else
	   PROJ_GRACEY 
	   #endif 
	   ) ;

  size_t mu ;
  for( mu = 0 ; mu < ND - 1 ; mu++ ) {
    sprintf( str , "%s%zux" , str , Latt.dims[mu] ) ;
  }
  sprintf( str , "%s%zu_%zu.config" , str , Latt.dims[ ND - 1 ] , nnmax ) ;

  FILE *config = fopen( str , "rb" ) ;
  int flag = 0 ;
  //force it to open ->create a file if needed
  if( config == NULL ) {
    flag = 1 ;
  } else {
    fclose( config ) ; 
  }

  // create a file
  size_t i ;
  if( unlikely( flag == 1 ) ) {
    init_delta( ) ;
    // open a temporary file
    FILE *config2 = fopen( str , "wb" ) ;

    printf("[CUTS] Storing Triplet and Proj list @@@ ...\n%s\n" , str ) ;

    // allocate our triplet 
    size_t **triple = ( size_t** )malloc( count * sizeof( size_t* ) ) ;
    #pragma omp parallel for private(i)
    PFOR( i = 0 ; i < count ; i++ ) {
      triple[i] = ( size_t* ) malloc( 3 *  sizeof ( size_t ) ) ;
    }
    printf( "[CUTS] Precomputing triplet \n" ) ;
    // compute the triplet in the stupidest manner possible
    get_triplet( triple , momentum , nnmax , nmom ) ;
    // write out "triple" our triplet
    for( i = 0 ; i < count ; i++ ) {
      fwrite( triple[i] , sizeof( size_t ) , 3 , config2 ) ;
    }
    printf( "Computing projector ... \n" ) ;
    // compute our projector 
    compute_projector( triple , proj , momentum , count ) ;
    // I want to write the projector here too
    for( i = 0 ; i < count ; i++ ) {
      fwrite( proj[i] , sizeof( double ) , ND * ND * ND , config2 ) ;
    }
    // close the temporary file
    fclose( config2 ) ;
    // free the temporary triplet "triple"
    #pragma omp parallel for private(i)
    PFOR( i = 0 ; i < count ; i++ ) {
      free( triple[i] ) ;
    }
    free( triple ) ; 
    free_delta( ) ;
  }

  // (re)open the file to read the triplet and the projector 
  config = fopen( str , "rb" ) ;
  for( i = 0 ; i < count ; i++ ) {
    if( fread( triplet[i] , sizeof( size_t ) , 3 , config ) != 3 ) { 
      return GLU_FAILURE ; 
    }
  }
  // read in the projector too, is directly beneath the triplet-idxs
  for( i = 0 ; i < count ; i++ ) {
    if( fread( proj[i] , sizeof( double ) , ND * ND * ND , config )
	!= ND * ND * ND ) {
      return GLU_FAILURE ; 
    }
  }
  // and close the file 
  fclose( config ) ;
#else
  // initialise the delta
  init_delta( ) ;
  // If the triplet fails we leave
  get_triplet( triplet , momentum , nnmax , nmom ) ;
  // compute the projector here 
  compute_projector( triplet , proj , momentum , count ) ;
  // and free the delta
  free_delta( ) ;
#endif
  return GLU_SUCCESS ;
}

// clean up the projectors
#undef p1
#undef p2
#undef p3
#undef p4 
#undef p5
#undef p6 
#undef p7 
#undef p8 
#undef p9 
#undef p10 
#undef p11 
#undef p12 
#undef p13 
#undef p14 
#undef DENOM 
