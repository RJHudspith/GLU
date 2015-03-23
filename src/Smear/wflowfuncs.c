/*
    Copyright 2013 Renwick James Hudspith

    This file (wflowfuncs.c) is part of GLU.

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
   @file wflowfuncs.c
   @brief helper functions for the wilson flow
 */

#include "Mainfile.h"

#include "GLU_splines.h"   // spline evaluations
#include "plaqs_links.h"   // clover discretisation
#include "projectors.h"    // stout projection
#include "staples.h"       // computes standard staples

// controls for the wilson flow these get externed!
const double MEAS_START = 0.0 ; // start measuring from 1 lattice spacing flow
const double WFLOW_STOP = 0.3 ; // BMW's choice for the W_0 parameter
const double TMEAS_STOP = 12.0 ; // flow time stopper, be careful after ~10

// RK4 parameters
static const double mthreeOfour = -0.75 ; // 3.0/4.0
static const double mseventeenOthsix = -0.47222222222222222222 ; // -17.0/36.0
static const double eightOnine = 0.88888888888888888889 ; // 8.0/9.0

// shortening function needs to be out of alphabetical order because
// it is called by flow directions
static INLINE_VOID
make_short_log( GLU_complex short_staple[ HERMSIZE ] , 
		const GLU_complex staple[ NCNC ] )
{
#if NC == 3
  *( short_staple + 0 ) = staple[ 0 ] ;
  *( short_staple + 1 ) = staple[ 1 ] ;
  *( short_staple + 2 ) = staple[ 2 ] ;
  *( short_staple + 3 ) = staple[ 4 ] ;
  *( short_staple + 4 ) = staple[ 5 ] ;
#elif NC == 2
  *( short_staple + 0 ) = staple[ 0 ] ;
  *( short_staple + 1 ) = staple[ 1 ] ;
#else
  int i, j , idx = 0 ;
  for( i = 0 ; i < NC-1 ; i++ ) {
    for( j = i ; j < NC ; j++ ) {
      short_staple[ idx ] = staple[ j + i*NC ] ;
      idx ++ ;
    }
  }
#endif
  return ;
}

// wilson flow in all ND directions ...
static void
flow_directions( lat2 , lat , Z , multiplier , delta_t , i , it , mu , SM_TYPE )
     struct spt_site *__restrict lat2 ;
     const struct site *__restrict lat ;
     struct spt_site_herm *__restrict Z ;
     const double multiplier ;
     const double delta_t ;
     const int i ;
     const int it ;
     const int mu ;
     const int SM_TYPE ;
{
  GLU_complex staple[ NCNC ] = { } ;
  // first element
  #ifdef IMPROVED_SMEARING
  all_staples_improve( staple , lat , it , mu , ND , SM_TYPE ) ;
  #else
  all_staples( staple , lat , it , mu , ND , SM_TYPE ) ;
  #endif

  GLU_complex short_staple[ HERMSIZE ] ; // does not need to be inited
  // default to STOUT ...
  if( SM_TYPE == SM_LOG ) {
    make_short_log( short_staple , staple ) ;
  } else {
    // multiply by the daggered link
    GLU_complex temp[ NCNC ] ;
    multab_dag( temp , staple , lat[it].O[mu] ) ; 
    // hermitian projection
    Hermitian_proj_short( short_staple , temp ) ;
  }
  #if NC == 3
  Z[it].O[mu][0] += multiplier * creal( short_staple[0] ) +\
    I * multiplier * creal( short_staple[3] ) ;
  Z[it].O[mu][1] += multiplier * short_staple[1] ;
  Z[it].O[mu][2] += multiplier * short_staple[2] ;
  Z[it].O[mu][3] += multiplier * short_staple[4] ;
  #elif NC == 2
  Z[it].O[mu][0] += multiplier * creal( short_staple[0] ) ;
  Z[it].O[mu][1] += multiplier * short_staple[1] ;
  #else
  int elements ;
  for( elements = 0 ; elements < HERMSIZE ; elements++ ) {
    Z[it].O[mu][ elements ] += multiplier * short_staple[elements] ;
  }
  #endif
  // perform the stout projection on the hermitian-shortened links
  switch( SM_TYPE )
    {
    case SM_LOG :
      project_LOG_wflow_short( lat2[i].O[mu] , Z[it].O[mu] , 
			       lat[it].O[mu] , delta_t ) ; 
      break ;
    default :
      project_STOUT_wflow_short( lat2[i].O[mu] , Z[it].O[mu] , 
				 lat[it].O[mu] , delta_t ) ; 
    }
  return ;
}

// memory-expensive rung-kutta step
static void
RK4step( lat , lat2 , Z , factor , RKfactor , SM_TYPE )
     struct site *__restrict lat ;
     struct spt_site *__restrict lat2 ;
     struct spt_site_herm *__restrict Z ;
     const double factor , RKfactor ;
     const int SM_TYPE ;
{
  int i ;
#pragma omp parallel for private(i) SCHED
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    int mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      GLU_complex staple[ NCNC ] = { } ;
      // first element
      #ifdef IMPROVED_SMEARING
      all_staples_improve( staple , lat , i , mu , ND , SM_TYPE ) ;
      #else
      all_staples( staple , lat , i , mu , ND , SM_TYPE ) ;
      #endif
      GLU_complex short_staple[ HERMSIZE ] ;
      // default to STOUT ...
      if( SM_TYPE == SM_LOG ){
	// log smearing just requires an exact log as
	// we have already multiplied through by the daggered
	// link matrix and taken the log
	make_short_log( short_staple , staple ) ;
      } else {
	// multiply by the daggered link
	GLU_complex temp[ NCNC ] ;
	multab_dag( temp , staple , lat[i].O[mu] ) ; 
	// hermitian projection
	Hermitian_proj_short( short_staple , temp ) ;
      }
      #if NC == 3
      Z[i].O[mu][0] += factor * creal( short_staple[0] ) +\
	I * factor * creal( short_staple[3] ) ;
      Z[i].O[mu][1] += factor * short_staple[1] ;
      Z[i].O[mu][2] += factor * short_staple[2] ;
      Z[i].O[mu][3] += factor * short_staple[4] ;
      #elif NC == 2
      Z[i].O[mu][0] += factor * creal( short_staple[0] ) ;
      Z[i].O[mu][1] += factor * short_staple[1] ;
      #else
      int elements ;
      for( elements = 0 ; elements < HERMSIZE ; elements++ ) {
	Z[i].O[mu][ elements ] += factor * short_staple[ elements ] ;
      }
      #endif
      switch( SM_TYPE )
	{
	case SM_LOG :
	  project_LOG_wflow_short( lat2[i].O[mu] , Z[i].O[mu] , 
				   lat[i].O[mu] , RKfactor ) ; 
	  break ;
	default :
	  project_STOUT_wflow_short( lat2[i].O[mu] , Z[i].O[mu] , 
				     lat[i].O[mu] , RKfactor ) ; 
	  break ;
	}
      // end of the mu-loop
    }
  }
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    memcpy( &lat[i] , &lat2[i] , sizeof( struct spt_site ) ) ; 
  }
  return ;
}

// computes one of the RK steps, doesn't matter which one
static void
RK4step_memcheap( lat , Z , lat2 , lat3 , lat4 , multiplier , step , SM_TYPE )
     struct site *__restrict lat ;
     struct spt_site_herm *__restrict Z ;
     struct spt_site *__restrict lat2 ;
     struct spt_site *__restrict lat3 ;
     struct spt_site *__restrict lat4 ;  
     const double multiplier ;
     const double step ;
     const int SM_TYPE ;
{
  int i ;

  #ifdef IMPROVED_SMEARING
  const int back = lat[ lat[0].back[ ND-1 ] ].back[ ND-1 ] ;
  #else
  const int back = lat[0].back[ ND - 1 ] ;
  #endif

  // split volume - wise
  #pragma omp parallel for private(i) SCHED
  #ifdef IMPROVED_SMEARING
  PFOR( i = 0 ; i < 2*LCU ; i++ ) {
  #else
  PFOR( i = 0 ; i < LCU ; i++ ) {
  #endif
    const int bck = back + i ;
    int mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      flow_directions( lat4 , lat , Z , multiplier , 
		       step , i , bck , mu , SM_TYPE ) ; 
    }
  }
  int t ;
  #ifdef IMPROVED_SMEARING
  for( t = 0 ; t < Latt.dims[ ND - 1 ] - 2 ; t++ ) {
  #else
  for( t = 0 ; t < Latt.dims[ ND - 1 ] - 1 ; t++ ) {
  #endif
    const int slice = LCU * t ;

    // perform the wilson flow for a point on the slice
    #pragma omp parallel for private(i) SCHED
    PFOR( i = 0 ; i < LCU ; i++ ) {
      const int it = slice + i ;
      int mu ;
      for( mu = 0 ; mu < ND ; mu++ ) {
	flow_directions( lat2 , lat , Z , multiplier , 
			 step , i , it , mu , SM_TYPE ) ;
      }
    }
    // swap over the temporary lattice fields
    #ifdef IMPROVED_SMEARING
    const int bck = lat[ lat[ slice ].back[ ND-1 ] ].back[ ND-1 ] ;
    #else
    const int bck = lat[ slice ].back[ ND -1 ] ;
    #endif

   #pragma omp parallel for private(i)
    PFOR( i = 0 ; i < LCU ; i++ ) {
      //put temp into the previous time-slice
      #ifdef IMPROVED_SMEARING
      if( likely( t > 1 ) ) { 
	register const int back = bck + i ;
	memcpy( &lat[back] , &lat3[i] , sizeof( struct spt_site ) ) ;
      }
      // put the evaluation two time slices before into the front half of lat3
      // and the evaluation one time slice before into the lower half
      memcpy( &lat3[i] , &lat3[i+LCU] , sizeof( struct spt_site ) ) ;
      memcpy( &lat3[i+LCU] , &lat2[i] , sizeof( struct spt_site ) ) ;
      #else
      if( likely( t != 0 ) ) { 
	register const int back = bck + i ;
	memcpy( &lat[back] , &lat3[i] , sizeof( struct spt_site ) ) ;
      }
      //make temporary lat3 lat2 again and repeat
      memcpy( &lat3[i] , &lat2[i] , sizeof( struct spt_site ) ) ;
      #endif
    }
  }
  // put the last couple back in ....
  const int slice = LCU * t ;
  const int behind = lat[ slice ].back[ ND - 1 ] ;
  #ifdef IMPROVED_SMEARING
  const int behind2 = lat[ behind ].back[ ND-1 ] ;
  #endif
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LCU ; i++ ) {
    register const int back = behind + i ;
    register const int it = slice + i ; 
    #ifdef IMPROVED_SMEARING
    memcpy( &lat[behind2+i] , &lat3[i] , sizeof( struct spt_site ) ) ; 
    memcpy( &lat[back] , &lat3[i+LCU] , sizeof( struct spt_site ) ) ; 
    memcpy( &lat[it] , &lat4[i] , sizeof( struct spt_site ) ) ; 
    memcpy( &lat[it+LCU] , &lat4[i+LCU] , sizeof( struct spt_site ) ) ; 
    #else
    memcpy( &lat[back] , &lat3[i] , sizeof( struct spt_site ) ) ; 
    memcpy( &lat[it] , &lat4[i] , sizeof( struct spt_site ) ) ; 
    #endif
  }
  return ;
}

// evaluate the flow using splines, evaluate a a particular scale
const double
evaluate_scale( double *der , 
		const double *x ,
		const double *meas ,
		const int Nmeas ,
		const double scale ,
		const char *message )
{
  // compute spline for observable meas, x must be sorted smallest to largest
  int i , change_up = 0 ;
  for( i = 0 ; i < Nmeas ; i++ ) {
    spline_derivative( der , x , meas , Nmeas ) ;
    // W_0 scale
    if( i > 0 ) {
      if( ( meas[ i ] ) > scale && ( meas[ i - 1 ] ) < scale ) {
	change_up = i ;
      }
    }
    #ifdef VERBOSE
    printf( "[%s] %g %g \n" , message , x[ i ] , meas[ i ] ) ;
    #endif
  }
  // print out the spline evaluation?
#ifdef VERBOSE
  double t = 0.0 ;
  for( t = 0.0 ; t < x[ Nmeas - 1 ] ; t+= 0.005 ) {
    printf( "[%s-spline] %g %g \n" , message , t , 
	    cubic_eval( x , meas , der , t , Nmeas ) ) ;
  }
#endif
  // evaluate at "scale" error flag is -1
  return solve_spline( x , meas , der , scale , change_up ) ;
}

// print out the general beginning information
void
print_GG_info( const int SM_TYPE , 
	       const wflow_type WFLOW_TYPE ) 
{
  printf( "[WFLOW] Taking ({W},{GG} and {Qtop}) measurements from t >= %g \n" , 
	  MEAS_START ) ; 
  printf( "[WFLOW] Stopping flow integration at w0 >= %g \n" , WFLOW_STOP ) ; 
  printf( "[WFLOW] OR, Stopping flow integration at t >= %g \n\n" , 
	  TMEAS_STOP ) ; 
  return ;
}

// use the flow of G and W for scale setting
int
scaleset( struct wfmeas *curr , 
	  const double T_0 ,
	  const double W_0 ,
	  const int count ) 
{
  // now we have the number of measurements in count
  double *GT = malloc( count * sizeof( double ) ) ;
  double *time = malloc( count * sizeof( double ) ) ;
  double *der = malloc( count * sizeof( double ) ) ;
  // traverse the list backwards setting the time and GT
  int i ;
  for( i = 0 ; i < count ; i++ ) {
    time[ count - i - 1 ] = curr -> time ;
    GT[ count - i - 1 ] = curr -> Gt ;
    curr = curr -> next ;
  }
  const double t0 = evaluate_scale( der , time , GT , count , T_0 , "GT" ) ;
  if( t0 == -1.0 ) {
    printf( "[WFLOW] cubic solve failure \n" ) ;
    printf( "[WFLOW] solve needs to bound the value you are looking for\n" ) ;
    free( der ) ; free( time ) ; free( GT ) ;
    return GLU_FAILURE ;
  }
  printf( "[T0-scale] %f \n" , sqrt( t0 ) ) ;
  // W(t) = t ( dG(t) / dt )
  for( i = 0 ; i < count ; i++ ) {
    GT[ i ] = time[ i ] * der[ i ] ;
  }
  const double w0 = evaluate_scale( der , time , GT , count , W_0 , "GT" ) ;
  if( w0 == -1.0 ) {
    printf( "[WFLOW] cubic solve failure \n" ) ;
    printf( "[WFLOW] solve needs to bound the value you are looking for\n" ) ;
    free( der ) ; free( time ) ; free( GT ) ;
    return GLU_FAILURE ;
  }
  printf( "[W0-scale] %f \n" , sqrt( w0 ) ) ;

  free( der ) ;
  free( time ) ;
  free( GT ) ;
  return GLU_SUCCESS ;
}

// driver function for the more memory expensive smearing
void
step_distance( struct site *__restrict lat ,
	       struct spt_site *__restrict lat2 ,
	       struct spt_site_herm *__restrict Z ,
	       const double rk1 ,
	       const double rk2 , 
	       const double rk3 , 
	       const int SM_TYPE )
{
  int i ;
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    memset( &Z[i] , 0.0 , sizeof( struct spt_site_herm ) ) ;
  }
  RK4step( lat , lat2 , Z , mseventeenOthsix , rk1 , SM_TYPE ) ;
  RK4step( lat , lat2 , Z , eightOnine , rk2 , SM_TYPE ) ;
  RK4step( lat , lat2 , Z , mthreeOfour , rk3 , SM_TYPE ) ;
  return ;
}

// perform a usual step
void
step_distance_memcheap( struct site *__restrict lat ,
			struct spt_site *__restrict lat2 ,
			struct spt_site *__restrict lat3 ,
			struct spt_site *__restrict lat4 ,
			struct spt_site_herm *__restrict Z ,
			const double rk1 ,
			const double rk2 , 
			const double rk3 , 
			const int SM_TYPE )
{
  // initial zero ...
  int i ;
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    memset( &Z[i] , 0.0 , sizeof( struct spt_site_herm ) ) ;
  }
  RK4step_memcheap( lat , Z , lat2 , lat3 , lat4 , mseventeenOthsix , rk1 , SM_TYPE ) ;
  RK4step_memcheap( lat , Z , lat2 , lat3 , lat4 , eightOnine , rk2 , SM_TYPE ) ;
  RK4step_memcheap( lat , Z , lat2 , lat3 , lat4 , mthreeOfour , rk3 , SM_TYPE ) ;
  return ;
}


