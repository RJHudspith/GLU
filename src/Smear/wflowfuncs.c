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
const double W0_STOP    = 0.3 ; // BMW's choice for the W_0 parameter
const double T0_STOP    = 0.3 ; // Martin's choice for the t_0 scale
double TMEAS_STOP = 20 ;

// set tmeas using some c0 if that is your thing
void 
set_TMEAS_STOP( const double c0 ) 
{ 
  TMEAS_STOP = ( c0 * Latt.dims[0] ) * ( c0 * Latt.dims[0] ) / 8 ; 
  fprintf( stdout , "[WFLOW] TMEAS_STOP changed to %f \n" , TMEAS_STOP ) ;
  return ;
}

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
  size_t i, j , idx = 0 ;
  for( i = 0 ; i < NC-1 ; i++ ) {
    for( j = i ; j < NC ; j++ ) {
      short_staple[ idx ] = staple[ j + i*NC ] ;
      idx ++ ;
    }
  }
#endif
  return ;
}

// accumulate the Z-matrix 
static inline void
set_zmatrix( struct spt_site_herm *__restrict Z ,
	     const GLU_complex short_staple[ HERMSIZE ] ,
	     const double multiplier ,
	     const size_t i ,
	     const size_t mu ) 
{
#if NC == 3
  Z[i].O[mu][0] += multiplier * creal( short_staple[0] ) +\
    I * multiplier * creal( short_staple[3] ) ;
  Z[i].O[mu][1] += multiplier * short_staple[1] ;
  Z[i].O[mu][2] += multiplier * short_staple[2] ;
  Z[i].O[mu][3] += multiplier * short_staple[4] ;
#elif NC == 2
  Z[i].O[mu][0] += multiplier * creal( short_staple[0] ) ;
  Z[i].O[mu][1] += multiplier * short_staple[1] ;
#else
  size_t elements ;
  for( elements = 0 ; elements < HERMSIZE ; elements++ ) {
    Z[i].O[mu][ elements ] += multiplier * short_staple[elements] ;
  }
#endif
}

// wilson flow in all ND directions ...
static void
flow_directions( struct spt_site *__restrict lat2 ,
		 struct spt_site_herm *__restrict Z ,
		 const struct site *__restrict lat ,
		 const double multiplier ,
		 const double delta_t ,
		 const size_t i ,
		 const size_t it ,
		 const size_t mu ,
		 const smearing_types SM_TYPE ,
		 void (*project)( GLU_complex log[ NCNC ] , 
				  GLU_complex *__restrict staple , 
				  const GLU_complex link[ NCNC ] , 
				  const double smear_alpha ) )
{
  GLU_complex staple[ NCNC ] GLUalign , temp[ NCNC ] GLUalign ;
  GLU_complex short_staple[ HERMSIZE ] GLUalign ; // does not need to be inited
  zero_mat( staple ) ;
  // first element
  #ifdef IMPROVED_SMEARING
  all_staples_improve( staple , lat , it , mu , ND , SM_TYPE ) ;
  #else
  all_staples( staple , lat , it , mu , ND , SM_TYPE ) ;
  #endif
  // default to STOUT ...
  switch( SM_TYPE ) {
  case SM_LOG :
    make_short_log( short_staple , staple ) ;
    break ;
  default :
    multab_dag( temp , staple , lat[it].O[mu] ) ; 
    Hermitian_proj_short( short_staple , temp ) ;
    break ;
  }
  // poke back into the accumulated Z-matrix
  set_zmatrix( Z , short_staple , multiplier , it , mu ) ;
  // perform the stout projection on the hermitian-shortened links
  project( lat2[i].O[mu] , Z[it].O[mu] , lat[it].O[mu] , delta_t ) ; 
  return ;
}

// memory-expensive runge-kutta step
static void
RK3step( struct spt_site_herm *__restrict Z ,
	 struct spt_site *__restrict lat2 ,
	 struct site *__restrict lat ,
	 const double multiplier , 
	 const double delta_t ,
	 const smearing_types SM_TYPE ,
	 void (*project)( GLU_complex log[ NCNC ] , 
			  GLU_complex *__restrict staple , 
			  const GLU_complex link[ NCNC ] , 
			  const double smear_alpha ) )
{
  size_t i ;
#pragma omp parallel for private(i) SCHED
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    GLU_complex staple[ NCNC ] GLUalign , temp[ NCNC ] GLUalign ;
    GLU_complex short_staple[ HERMSIZE ] GLUalign ;
    size_t mu ;
    switch( SM_TYPE ) {
    case SM_LOG :
      for( mu = 0 ; mu < ND ; mu++ ) {
	zero_mat( staple ) ;
	// first element
        #ifdef IMPROVED_SMEARING
	all_staples_improve( staple , lat , i , mu , ND , SM_TYPE ) ;
        #else
	all_staples( staple , lat , i , mu , ND , SM_TYPE ) ;
        #endif
	// log smearing just requires an set the log
	make_short_log( short_staple , staple ) ;
	set_zmatrix( Z , short_staple , multiplier , i , mu ) ;
	project( lat2[i].O[mu] , Z[i].O[mu] , lat[i].O[mu] , delta_t ) ; 
      } break ;
    default :
      for( mu = 0 ; mu < ND ; mu++ ) {
	zero_mat( staple ) ;
	// first element
        #ifdef IMPROVED_SMEARING
	all_staples_improve( staple , lat , i , mu , ND , SM_TYPE ) ;
        #else
	all_staples( staple , lat , i , mu , ND , SM_TYPE ) ;
        #endif
	// log smearing just requires an set the log
	// hermitian projection
	multab_dag( temp , staple , lat[i].O[mu] ) ; 
	Hermitian_proj_short( short_staple , temp ) ;
	set_zmatrix( Z , short_staple , multiplier , i , mu ) ;
	project( lat2[i].O[mu] , Z[i].O[mu] , lat[i].O[mu] , delta_t ) ;
      } break ;
    }
  }
  // copy into lat
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    memcpy( &lat[i] , &lat2[i] , sizeof( struct spt_site ) ) ; 
  }
  return ;
}

// computes one of the RK steps, doesn't matter which one
static void
RK3step_memcheap( struct spt_site_herm *__restrict Z ,
		  struct spt_site *__restrict lat2 ,
		  struct spt_site *__restrict lat3 ,
		  struct spt_site *__restrict lat4 ,  
		  struct site *__restrict lat ,
		  const double multiplier ,
		  const double step ,
		  const smearing_types SM_TYPE ,
		  void (*project)( GLU_complex log[ NCNC ] , 
				   GLU_complex *__restrict staple , 
				   const GLU_complex link[ NCNC ] , 
				   const double smear_alpha ) )
{
  size_t i ;
#ifdef IMPROVED_SMEARING
  const size_t back = lat[ lat[0].back[ ND-1 ] ].back[ ND-1 ] ;
#else
  const size_t back = lat[0].back[ ND - 1 ] ;
#endif
  // split volume - wise
  #pragma omp parallel for private(i) SCHED
  #ifdef IMPROVED_SMEARING
  PFOR( i = 0 ; i < 2*LCU ; i++ ) {
  #else
  PFOR( i = 0 ; i < LCU ; i++ ) {
  #endif
    const size_t bck = back + i ;
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      flow_directions( lat4 , Z , lat , multiplier , 
		       step , i , bck , mu , SM_TYPE , project ) ; 
    }
  }
  size_t t ;
  #ifdef IMPROVED_SMEARING
  for( t = 0 ; t < Latt.dims[ ND - 1 ] - 2 ; t++ ) {
  #else
  for( t = 0 ; t < Latt.dims[ ND - 1 ] - 1 ; t++ ) {
  #endif
    const size_t slice = LCU * t ;
    // perform the wilson flow for a point on the slice
    #pragma omp parallel for private(i) SCHED
    PFOR( i = 0 ; i < LCU ; i++ ) {
      const size_t it = slice + i ;
      size_t mu ;
      for( mu = 0 ; mu < ND ; mu++ ) {
	flow_directions( lat2 , Z , lat , multiplier , 
			 step , i , it , mu , SM_TYPE , project ) ;
      }
    }
    // swap over the temporary lattice fields
#ifdef IMPROVED_SMEARING
    const size_t bck = lat[ lat[ slice ].back[ ND-1 ] ].back[ ND-1 ] ;
#else
    const size_t bck = lat[ slice ].back[ ND -1 ] ;
#endif
#pragma omp parallel for private(i)
    PFOR( i = 0 ; i < LCU ; i++ ) {
      //put temp into the previous time-slice
      #ifdef IMPROVED_SMEARING
      if( likely( t > 1 ) ) { 
	register const size_t back = bck + i ;
	memcpy( &lat[back] , &lat3[i] , sizeof( struct spt_site ) ) ;
      }
      // put the evaluation two time slices before into the front half of lat3
      // and the evaluation one time slice before into the lower half
      memcpy( &lat3[i] , &lat3[i+LCU] , sizeof( struct spt_site ) ) ;
      memcpy( &lat3[i+LCU] , &lat2[i] , sizeof( struct spt_site ) ) ;
      #else
      if( likely( t != 0 ) ) { 
	register const size_t back = bck + i ;
	memcpy( &lat[back] , &lat3[i] , sizeof( struct spt_site ) ) ;
      }
      //make temporary lat3 lat2 again and repeat
      memcpy( &lat3[i] , &lat2[i] , sizeof( struct spt_site ) ) ;
      #endif
    }
  }
  // put the last couple back in ....
  const size_t slice = LCU * t ;
  const size_t behind = lat[ slice ].back[ ND - 1 ] ;
#ifdef IMPROVED_SMEARING
  const size_t behind2 = lat[ behind ].back[ ND-1 ] ;
#endif
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LCU ; i++ ) {
    register const size_t back = behind + i ;
    register const size_t it = slice + i ; 
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
static double
evaluate_scale( double *der , 
		const double *x ,
		const double *meas ,
		const size_t Nmeas ,
		const double scale ,
		const char *message )
{
  // compute spline for observable meas, x must be sorted smallest to largest
  size_t i , change_up = 0 ;
  for( i = 0 ; i < Nmeas ; i++ ) {
    spline_derivative( der , x , meas , Nmeas ) ;
    // W_0 scale
    if( i > 0 && ( meas[ i ] > scale ) && ( meas[ i - 1 ] < scale ) ) {
      change_up = i ;
    }
    #ifdef verbose
    fprintf( stdout , "[%s] %g %g \n" , message , x[ i ] , meas[ i ] ) ;
    #endif
  }
  // print out the spline evaluation?
  if( Nmeas > 0 ) {
    double t = 0.0 ;
    for( t = 0.0 ; t < x[ Nmeas - 1 ] ; t+= 0.005 ) {
      fprintf( stdout , "[%s-spline] %g %g \n" , message , t , 
	       cubic_eval( x , meas , der , t , Nmeas ) ) ;
    }
  }
  // evaluate at "scale" error flag is -1
  return solve_spline( x , meas , der , scale , change_up ) ;
}

// print out the general beginning information
void
print_GG_info( void ) 
{
  fprintf( stdout , "[WFLOW] Taking ({W},{GG} and {Qtop}) measurements"
	   "from t >= %g \n" , MEAS_START ) ; 
  fprintf( stdout , "[WFLOW] fine measurements at t_0 >= %g \n" , T0_STOP ) ; 
  fprintf( stdout , "[WFLOW] fine measurements at w_0 >= %g \n" , W0_STOP ) ; 
  fprintf( stdout , "[WFLOW] OR, Stopping flow integration at t >= %g \n\n" , 
	  TMEAS_STOP ) ; 
  return ;
}

// use the flow of G and W for scale setting
int
scaleset( struct wfmeas *curr , 
	  const double T_0 ,
	  const double W_0 ,
	  const size_t count ) 
{
  // now we have the number of measurements in count
  double *GT   = malloc( ( count + 1 ) * sizeof( double ) ) ;
  double *time = malloc( ( count + 1 ) * sizeof( double ) ) ;
  double *der  = malloc( ( count + 1 ) * sizeof( double ) ) ;
  // error handling
  double t0 = -1 , w0 = -1 ;
  int flag = GLU_SUCCESS ;
  // traverse the list backwards setting the time and GT
  size_t i ;
  // traverse back down the linked list
  for( i = 0 ; i < ( count + 1 ) ; i++ ) {
    time[ count - i ] = curr -> time ;
    GT[ count - i ] = curr -> Gt ;
    curr = curr -> next ;
  }
  if( count > 0 ) {
    t0 = evaluate_scale( der , time , GT , count , T_0 , "GT" ) ;
  }
  if( t0 == -1 ) {
    fprintf( stderr , "[WFLOW] cubic solve failure (gt)\n" ) ;
    fprintf( stderr , "[WFLOW] solve needs to bound the value you "
	     "are looking for\n" ) ;
    flag = GLU_FAILURE ;
    goto free ;
  }
  fprintf( stdout , "[GT-scale] G(%g) %f \n" , T_0 , sqrt( t0 ) ) ;
  // W(t) = t ( dG(t) / dt )
  if( count > 0 ) {
    for( i = 0 ; i < ( count + 1 ) ; i++ ) {
      GT[ i ] = time[ i ] * der[ i ] ;
    }
  }
  w0 = evaluate_scale( der , time , GT , count , W_0 , "WT" ) ;
  if( w0 == -1 ) {
    fprintf( stderr , "[WFLOW] cubic solve failure (wt) \n" ) ;
    fprintf( stderr , "[WFLOW] solve needs to bound the value you "
	     "are looking for\n" ) ;
    flag = GLU_FAILURE ;
    goto free ;
  }
  fprintf( stdout , "[WT-scale] W(%g) %g \n" , W_0 , sqrt( w0 ) ) ;

 free :
  free( der ) ;
  free( time ) ;
  free( GT ) ;

  return flag ;
}

// driver function for the more memory expensive smearing
void
step_distance( struct site *__restrict lat ,
	       struct spt_site *__restrict lat2 ,
	       struct spt_site_herm *__restrict Z ,
	       const double rk1 ,
	       const double rk2 , 
	       const double rk3 , 
	       const smearing_types SM_TYPE ,
	       void (*project)( GLU_complex log[ NCNC ] , 
				GLU_complex *__restrict staple , 
				const GLU_complex link[ NCNC ] , 
				const double smear_alpha ) )
{
  // RK4 parameters
  const double mthreeOfour = -3.0/4.0 ;
  const double mseventeenOthsix = -17.0/36.0 ;
  const double eightOnine = 8.0/9.0 ;
  // set z to zero
  size_t i ;
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    memset( &Z[i] , 0.0 , sizeof( struct spt_site_herm ) ) ;
  }
  // flow forwards one fictitious timestep
  RK3step( Z , lat2 , lat , mseventeenOthsix , rk1 , SM_TYPE , project ) ;
  RK3step( Z , lat2 , lat , eightOnine , rk2 , SM_TYPE , project ) ;
  RK3step( Z , lat2 , lat , mthreeOfour , rk3 , SM_TYPE , project ) ;
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
			const smearing_types SM_TYPE ,
			void (*project)( GLU_complex log[ NCNC ] , 
					 GLU_complex *__restrict staple , 
					 const GLU_complex link[ NCNC ] , 
					 const double smear_alpha ) )
{
  // RK4 parameters
  const double mthreeOfour = -3.0/4.0 ;
  const double mseventeenOthsix = -17.0/36.0 ;
  const double eightOnine = 8.0/9.0 ;
  // set z to zero
  size_t i ;
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    memset( &Z[i] , 0.0 , sizeof( struct spt_site_herm ) ) ;
  }
  // flow forwards one timestep
  RK3step_memcheap( Z , lat2 , lat3 , lat4 , 
		    lat , mseventeenOthsix , rk1 , SM_TYPE , project ) ;
  RK3step_memcheap( Z , lat2 , lat3 , lat4 , 
		    lat , eightOnine , rk2 , SM_TYPE , project ) ;
  RK3step_memcheap( Z , lat2 , lat3 , lat4 , 
		    lat , mthreeOfour , rk3 , SM_TYPE , project ) ;
  return ;
}


