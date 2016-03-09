/*
    Copyright 2013 Renwick James Hudspith

    This file (adaptive_flow.c) is part of GLU.

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
   @file adaptive_flow.c
   @brief the (two step) adaptive rk4 wilson flow routine

   Slows down, performing fine measurements at ~t_0 and ~w_0
   W0_STOP and T0_STOP are defined in wflowfuncs.h
 */

#include "Mainfile.h"

#include "geometry.h"     // init_navig is called for the temporary
#include "GLU_splines.h"  // cubic_eval and derivative calculation
#include "plaqs_links.h"  // clover and plaquette measurements
#include "projectors.h"   // smearing projections
#include "wflowfuncs.h"   // wilson flow general routines

// enable this if we are doing a straight shot for a specific time
// to avoid measuring the topological charge and clover and all that,
// be careful when T > 10, will need to increase ADAPTIVE_EPS
//#define TIME_ONLY

/**
   @fn static inline double adaptfmax( const double a , const double b )
   @brief the maximum of two numbers
   Is probably unsafe because Infs and Nans are not considered, although
   what would you do with them?
 */
static inline double
adaptfmax( a , b )
     const double a ;
     const double b ;
{
  return ( b < a ? a : b ) ;
}

/**
   @fn static inline double adaptfmin( const double a , const double b )
   @brief the minimum of two numbers
   Again probably unsafe
 */
static inline double
adaptfmin( a , b )
     const double a ;
     const double b ;
{
  return ( a < b ? a : b ) ;
}

/**
   @brief print the flow progress
 */
static void
print_flow( const struct wfmeas *curr ,
	    const double err ,
	    const double delta_t)
{
  if( delta_t < 0 ) {
    printf( "[WFLOW-TSTOP] {err} %1.3e {t} %f {dt} %g " , err , curr -> time , delta_t ) ;
  } else {
    printf( "[WFLOW] {err} %1.3e {t} %f {dt} %g " , err , curr -> time , delta_t ) ;
  }
  printf( "{p} %g {q} %g {Gt} %g \n" , curr -> avplaq , curr -> qtop , curr -> Gt ) ;
}

/**
   @enum adaptive_control
   @brief when to break our adaptive algorithm if we have done this many halvings and still have no result
 */
enum adaptive_control{ ADAPTIVE_BIG_NUMBER = 20 } ;

/**
   @brief fine measurements
 */
static struct wfmeas *
fine_measurement( struct site *lat , 
		  struct spt_site *lat2 , 
		  struct spt_site *lat3 , 
		  struct spt_site *lat4 , 
		  struct spt_site_herm *Z , 
		  double *flow_next , 
		  double *t , 
		  const double delta_t ,
		  const double preverr , 
		  const smearing_types SM_TYPE  , 
		  void (*project)( GLU_complex log[ NCNC ] , 
				   GLU_complex *__restrict staple , 
				   const GLU_complex link[ NCNC ] , 
				   const double smear_alpha )  )
{
  struct wfmeas *curr = (struct wfmeas*)malloc( sizeof( struct wfmeas ) ) ;
  *t += delta_t ;
  curr -> time = *t ;
  const double rk1 =-0.52941176470588235294 * delta_t ;
  const double rk2 = delta_t ;
  const double rk3 = ( -delta_t ) ;  
  step_distance_memcheap( lat , lat2 , lat3 , lat4 , 
			  Z , rk1 , rk2 , rk3 , SM_TYPE , project ) ;
  // update the linked list
  curr -> Gt = curr -> time * curr -> time * 
    lattice_gmunu( lat , &(curr -> qtop) , &( curr->avplaq ) ) ;
  // set the flow
  *flow_next = curr -> Gt ;
  // the error in a fine measurement will be at least that of the previous
  print_flow( curr , preverr , delta_t ) ;
  return curr ;
}

// Adaptive stepsize version 
int 
flow4d_adaptive_RK( struct site *__restrict lat , 
		    const size_t smiters ,
		    const size_t DIR ,
		    const int SIGN ,
		    const smearing_types SM_TYPE )
{  
  // set this for the coupling measurement
  //set_TMEAS_STOP( 0.4 ) ;

  ////// USUAL STARTUP INFORMATION /////////
  print_GG_info( SM_TYPE , RK4_ADAPTIVE ) ;

  // the error between the two plaquettes
  const double ADAPTIVE_EPS = 1E-7 ;
  // Standard shrink and factor from NRC RK4
  const double ADAPTIVE_SHRINK = -0.25 ; // 0.33?
  // Standard growth and factor from NRC RK4
  const double ADAPTIVE_GROWTH = -0.20 ; // 0.25?
  // define adaptive safe
  const double ADAPTIVE_SAFE = 0.9 ;
  // adaptive error conserving
  const double ADAPTIVE_ERRCON = powl( 5./ADAPTIVE_SAFE , 1./ADAPTIVE_GROWTH ) ;
  // percentage to value we want for performing fine measurements
  const double FINETWIDDLE = 0.05 ;
  // fine measurement step
  const double FINESTEP = 0.02 ;

  // adaptive factors for RK4, we are RK3 could be more lenient?
  fprintf( stdout , "[WFLOW] Adaptive Error :: %e \n" , ADAPTIVE_EPS ) ;
  fprintf( stdout , "[WFLOW] Adaptive ErrCon :: %f \n" , ADAPTIVE_ERRCON ) ;
  fprintf( stdout , "[WFLOW] Safety Factor :: %g \n" , ADAPTIVE_SAFE ) ;
  fprintf( stdout , "[WFLOW] Growth factor :: %g \n" , ADAPTIVE_GROWTH ) ;
  fprintf( stdout , "[WFLOW] Shrink factor :: %g \n\n" , ADAPTIVE_SHRINK ) ; 
  fprintf( stdout , "[WFLOW] Fine measurement %% :: %g \n" , FINETWIDDLE ) ;
  fprintf( stdout ,  "[WFLOW] Fine step :: %g \n\n" , FINESTEP ) ;

  //////////////////////////////////////////
  struct spt_site_herm *Z = NULL ;
  struct spt_site *lat2 = NULL , *lat3 = NULL , *lat4 = NULL ;
  struct site *lat_two = NULL ;

  void (*project) ( GLU_complex log[ NCNC ] , 
		    GLU_complex *__restrict staple , 
		    const GLU_complex link[ NCNC ] , 
		    const double smear_alpha ) ;
  
  // stout is the usual
  project = project_STOUT_wflow_short ;
  switch( SM_TYPE ) {
  case SM_STOUT : break ;
  case SM_LOG :
    project = project_LOG_wflow_short ;
    break ;
  default :
    fprintf( stderr , "[SMEARING] unrecognised smearing projection \n" ) ;
    return GLU_FAILURE ;
  }

  if( GLU_malloc( (void**)&Z , ALIGNMENT , LVOLUME * sizeof( struct spt_site_herm ) ) != 0 ||
      GLU_malloc( (void**)&lat2 , ALIGNMENT , LCU * sizeof( struct spt_site ) )       != 0 ||
      GLU_malloc( (void**)&lat_two , ALIGNMENT , LVOLUME * sizeof( struct site ) )    != 0 ) {
    fprintf( stderr , "[SMEARING] allocation failure \n" ) ;
    return GLU_FAILURE ;
  }
#ifdef IMPROVED_SMEARING
  if( GLU_malloc( (void**)&lat3 , ALIGNMENT , 2 * LCU * sizeof( struct spt_site ) )   != 0 ||
      GLU_malloc( (void**)&lat4 , ALIGNMENT , 2 * LCU * sizeof( struct spt_site ) )   != 0 ) {
    fprintf( stderr , "[SMEARING] allocation failure \n" ) ;
    return GLU_FAILURE ;
  }
#else
  if( GLU_malloc( (void**)&lat3 , ALIGNMENT , LCU * sizeof( struct spt_site ) )   != 0 ||
      GLU_malloc( (void**)&lat4 , ALIGNMENT , LCU * sizeof( struct spt_site ) )   != 0 ) {
    fprintf( stderr , "[SMEARING] allocation failure \n" ) ;
    return GLU_FAILURE ;
  }
#endif
  init_navig( lat_two ) ;

  // set up the step sizes ...
  double delta_t = SIGN * Latt.sm_alpha[0] ;

  struct wfmeas *head = NULL , *curr ;
  curr = (struct wfmeas*)malloc( sizeof( struct wfmeas ) ) ;
  curr -> time = 0.0 ;
  curr -> Gt = curr -> time * curr -> time * 
    lattice_gmunu( lat , &( curr -> qtop ) , &( curr -> avplaq ) ) ;
  print_flow( curr , 0.0 , 0.0 ) ;
  curr -> next = head ;
  head = curr ;

  // counters for the derivative ...
  double yscal = curr -> avplaq , t = 0.0 ;
  double flow = 0. , flow_next = 0. ;
  size_t count = 0 , OK_STEPS = 0 , NOTOK_STEPS = 0 ; 

  // have a flag for whether we mess up
  int FLAG = GLU_FAILURE ;

  // loop up to smiters
  for( count = 1 ; count <= smiters ; count++ ) { 
    curr = (struct wfmeas*)malloc( sizeof( struct wfmeas ) ) ;
    size_t counter = 0 ;
    double errmax = 10. ;
    double new_plaq = 0. ;
    size_t i ;
    // adaptive loop, shrinking or growing the stepsize accordingly
    while( ( errmax > 1.0 ) && ( counter < ADAPTIVE_BIG_NUMBER ) ) {
      #pragma omp parallel for private(i)
      PFOR( i = 0 ; i < LVOLUME ; i++ ) {
	memcpy( &lat_two[i] , &lat[i] , sizeof( struct site ) ) ; 
      }

      // Step forward in two halves ...
      const double rk1 = -0.52941176470588235294 * delta_t ;
      const double rk2 = delta_t ;
      const double rk3 = ( -delta_t ) ;  

      // step forward once and write into lat_two
      step_distance_memcheap( lat_two , lat2 , lat3 , lat4 , Z , 
			      rk1 , rk2 , rk3 , SM_TYPE , project ) ;
      
      // compute the one-step first comparison
      const double old_plaq = av_plaquette( lat_two ) ;

      #pragma omp parallel for private(i)
      PFOR( i = 0 ; i < LVOLUME ; i++ ) {
	memcpy( &lat_two[i] , &lat[i] , sizeof( struct site ) ) ; 
      } 

      // and step forward twice and write into lat_two
      const double half_rk1 = 0.5 * rk1 ;
      const double half_rk2 = 0.5 * rk2 ;
      const double half_rk3 = 0.5 * rk3 ;

      step_distance_memcheap( lat_two , lat2 , lat3 , lat4 , Z , 
			      half_rk1 , half_rk2 , half_rk3 , 
			      SM_TYPE , project ) ;

      step_distance_memcheap( lat_two , lat2 , lat3 , lat4 , Z , 
			      half_rk1 , half_rk2 , half_rk3 , 
			      SM_TYPE , project ) ;
	  
      // compute the error I will use the average plaquette ...
      new_plaq = (double)av_plaquette( lat_two ) ;

      errmax = fabsl( ( new_plaq - old_plaq ) / ( yscal ) ) ;
      errmax /= ADAPTIVE_EPS ;

      // Break the while loop if conditions are satisfied
      if( errmax < 1.0 ) {
	if( counter < 1 ) { OK_STEPS ++ ; } 
	break ;
      }
      // Increment the counter for not adequate steps
      NOTOK_STEPS ++ ;

      // shorten the delta_t ...
      const double del_temp =  ADAPTIVE_SAFE * delta_t * pow( errmax , ADAPTIVE_SHRINK ) ; 
      // set up a tolerance s.t del_temp doesn't go too crazy, only really used when starting guess is bad
      const double tol = 0.1 * delta_t ; 
      delta_t = ( 0. < del_temp ?			\
		  adaptfmax( del_temp , tol ) :		\
		  adaptfmin( del_temp , tol ) ) ;

      // Print if we are in trouble
      if( counter == ADAPTIVE_BIG_NUMBER - 1 ) {
	fprintf( stderr , "[WFLOW] Not stepping to required accuracy"
		 "after < %zu > attempts! \n" , counter ) ;
	goto memfree ;
      }

      // Increment our counter ...
      counter ++ ;
    }

    // set up a scaling parameter to control the adaptation uses a first order finite difference def ...
    const double yscal_new = new_plaq ;
    yscal = 2.0 * yscal_new - yscal ;

    // overwrite lat .. 
    struct site *ptr = lat ;
    lat = lat_two ;
    lat_two = ptr ;

    t += delta_t ; // add one time step to the overall time 

    curr -> time = t ;
    // update the linked list
    curr -> Gt = curr -> time * curr -> time * 
      lattice_gmunu( lat , &(curr -> qtop) , &( curr->avplaq ) ) ;

    // set the flow
    flow_next = curr -> Gt ;

    print_flow( curr , errmax * ADAPTIVE_EPS , delta_t ) ;

    // update the linked list
    curr -> next = head ;
    head = curr ;

    // approximate the derivative
    double wapprox = ( flow_next - flow ) * curr -> time / delta_t ;
    flow = flow_next ;

    // If we get stuck in a rut of updating by zero we leave
    if( fabs( delta_t ) < DBL_MIN ) {
      fprintf( stderr , "[WFLOW] No update made delta_t :: %1.5e \n"
	       "Leaving ... \n" , delta_t ) ; 
      goto memfree ;
    }

    // perform fine measurements around T0_STOP for t_0
    if( fabs( T0_STOP - flow ) <= ( T0_STOP * FINETWIDDLE ) ) {
      while( fabs( T0_STOP - flow ) <= ( T0_STOP * FINETWIDDLE ) ) {
	curr = fine_measurement( lat , lat2 , lat3 , lat4 , Z , 
				 &flow_next , &t , FINESTEP , 
				 errmax * ADAPTIVE_EPS , SM_TYPE , project ) ;
	wapprox = ( flow_next - flow ) * curr -> time / FINESTEP ;
	flow = flow_next ;
	curr -> next = head ;
	head = curr ;
	count++ ;
      }
    }

    // perform some fine measurements around W0_STOP
    if( fabs( W0_STOP - wapprox ) <= ( W0_STOP * FINETWIDDLE ) ) {
      while( fabs( W0_STOP - wapprox ) <= ( W0_STOP * FINETWIDDLE ) ) {
	curr = fine_measurement( lat , lat2 , lat3 , lat4 , Z , 
				 &flow_next , &t , FINESTEP , 
				 errmax * ADAPTIVE_EPS , SM_TYPE , project ) ;
	wapprox = ( flow_next - flow ) * curr -> time / FINESTEP ;
	flow = flow_next ;
	curr -> next = head ;
	head = curr ;
	count++ ;
      }
    }

    // use a poor approximation of the derivative of the flow as a guide to stop
    if( wapprox > ( W0_STOP * 1.5 ) ) {
      break ;
    }

    // stop if we are above the max time
    if( ( curr -> time ) > TMEAS_STOP ) {
      delta_t = TMEAS_STOP - t ;
      curr = fine_measurement( lat , lat2 , lat3 , lat4 , Z , 
			       &flow_next , &t , delta_t , 
			       errmax * ADAPTIVE_EPS , SM_TYPE , project ) ;
      count++ ;
      curr -> next = head ;
      head = curr ;
      break ;
    }
 
    // Increase the step size ...
    if( errmax > ADAPTIVE_ERRCON ) {
      delta_t = ADAPTIVE_SAFE * delta_t * pow( errmax , ADAPTIVE_GROWTH ) ;
    } else {
      delta_t = ADAPTIVE_SAFE * 5.0 * delta_t ;
    }
  }
  curr = head ;

  // Print out the stepping information
  fprintf( stdout , "\n[WFLOW] Inadequate steps :: %zu \n" , NOTOK_STEPS ) ;
  fprintf( stdout , "[WFLOW] Adequate steps :: %zu \n" , OK_STEPS ) ;

  // compute the t_0 and w_0 scales from the measurement
  if( ( fabs( curr -> time - TMEAS_STOP ) > PREC_TOL ) && 
      count <= smiters ) {
    scaleset( curr , T0_STOP , W0_STOP , count ) ;
  }

  // we are successful
  FLAG = GLU_SUCCESS ;

 memfree :

  // reswap pointers if we took a step, lat is freed elsewhere
  if( count > 1 ) {
    struct site *ptr = lat_two ;
    lat_two = lat ;
    lat = ptr ;
  }

  // and free the list
  while( ( curr = head ) != NULL ) {
    head = head -> next ; 
    free( curr ) ;
  }
  
  // free our fields
  free( Z ) ;
  free( lat2 ) ;
  free( lat3 ) ;
  free( lat4 ) ;
  free( lat_two ) ;

  return FLAG ;
}

// if we have set the code in "TIME_ONLY" mode we make sure we clean it up
#ifdef TIME_ONLY
  #undef TIME_ONLY
#endif
