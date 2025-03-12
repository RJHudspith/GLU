/*
Copyright 2013-2025 Renwick James Hudspith

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
   @brief the (two step) adaptive rk3 wilson flow routine

   Slows down, performing fine measurements at ~t_0 and ~w_0
   W0_STOP and T0_STOP are defined in wflowfuncs.h
 */
#include "Mainfile.h"

#include "init.h"         // init_navig is called for the temporary
#include "plaqs_links.h"  // av_plaquette()
#include "projectors.h"   // smearing projections
#include "wflowfuncs.h"   // wilson flow general routines

#include "GLU_timer.h"

/**
   @fn static inline double adaptfmax( const double a , const double b )
   @brief the maximum of two numbers
   Is probably unsafe because Infs and Nans are not considered, although
   what would you do with them?
 */
static inline double
adaptfmax( const double a , const double b )
{
  return ( b < a ? a : b ) ;
}

/**
   @fn static inline double adaptfmin( const double a , const double b )
   @brief the minimum of two numbers
   Again probably unsafe
 */
static inline double
adaptfmin( const double a , const double b )
{
  return ( a < b ? a : b ) ;
}

/**
   @enum adaptive_control
   @brief when to break our adaptive algorithm if we have done this many halvings and still have no result
 */
enum adaptive_control{ ADAPTIVE_BIG_NUMBER = 20 } ;

// the error between the two plaquettes
static const double ADAPTIVE_EPS = 2E-7 ;

// adaptive parameters have been tuned to be reasonably well behaved
#define ADAPTIVE_SHRINK (-0.33)
#define ADAPTIVE_GROWTH (-0.25)
#define ADAPTIVE_SAFE (0.9)

static inline void
setred( struct wflow_temps *WF )
{
  size_t i ;
#pragma omp for private(i)
  for( i = 0 ; i < CLINE*Latt.Nthreads ; i++ ) {
    WF -> red[i] = 0.0 ;
  }
}

// little plaquett utility
static inline double
get_plaq_th( struct site *lat ,
	     struct wflow_temps WF )
{
  av_plaquette_th( WF.red , lat ) ;
  double new_plaq = 0.0 ;
  for( size_t k = 0 ; k < Latt.Nthreads ; k++ ) {
    new_plaq += WF.red[ 3 + CLINE*k ] ;
  }
  return new_plaq/( NC * ( ND-1 ) * ( ND - 2 ) * LVOLUME ) ;
}

// two step adaptive routine
static int
twostep_adaptive( struct site *lat ,
		  struct wflow_temps WF ,
		  double *dt ,
		  double *errmax ,
		  double *new_plaq ,
		  const double t ,
		  const double yscal ,
		  const smearing_types SM_TYPE  , 
		  void (*project)( GLU_complex log[ NCNC ] , 
				   GLU_complex *__restrict staple , 
				   const GLU_complex link[ NCNC ] , 
				   const double smear_alpha )  )
{
  *errmax = 10. ;
  size_t counter = 0 , i , mu ;

  // if we hit the 1/8 we just set it to that and do a normal step
  if( *dt >= 0.125 && t > 1.0 ) {
    *dt = 0.125 ;
    #pragma omp for private(i) collapse(2)
    for( i = 0 ; i < LVOLUME ; i++ ) {
      for( mu = 0 ; mu < ND ; mu++ ) {
	equiv( WF.lat_two[i].O[mu] , lat[i].O[mu] ) ;
      }
    }
    setred( &WF ) ;
    step_distance( WF.lat_two , WF , *dt , SM_TYPE , project ) ;
    *new_plaq = get_plaq_th( WF.lat_two , WF ) ;
    return GLU_SUCCESS ;
  }
  
  // adaptive loop, shrinking or growing the stepsize accordingly
 top :
  
  {
    #pragma omp barrier
  }
  
#pragma omp for private(i) collapse(2)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    for( mu = 0 ; mu < ND ; mu++ ) {
      equiv( WF.lat_two[i].O[mu] , lat[i].O[mu] ) ;
    }
  }
  setred( &WF ) ;
  
  // step forward once and write into lat_two
  step_distance( WF.lat_two , WF , *dt , SM_TYPE , project ) ;
  
  // compute the one-step first comparison
  const double old_plaq = get_plaq_th( WF.lat_two , WF ) ;
  
#pragma omp for private(i) collapse(2)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    for( mu = 0 ; mu < ND ; mu++ ) {
      equiv( WF.lat_two[i].O[mu] , lat[i].O[mu] ) ;
    }
  }
  setred( &WF ) ;
    
  // and step forward twice and write into lat_two
  step_distance( WF.lat_two , WF , 0.5*(*dt) , SM_TYPE , project ) ;
  step_distance( WF.lat_two , WF , 0.5*(*dt) , SM_TYPE , project ) ;
    
  // compute the error I will use the average plaquette ...
  *new_plaq = get_plaq_th( WF.lat_two , WF ) ;
    
  // so the problem is that at large t the plaquettes become very close so the 
  // step gets pretty wild. My way of compensating this is just to multiply by t^2
  // this is in some sense tuning the stepsize for the quantity t^2 E^2
  *errmax = fabsl( (t+*dt)*(t+*dt) * ( *new_plaq - old_plaq ) / ( yscal ) ) ;
  *errmax /= ADAPTIVE_EPS ;
  
  // Break the while loop if conditions are satisfied
  if( *errmax > 1.0 && counter < 20 ) {
    // shorten the delta_t ...
    const double del_temp =  ADAPTIVE_SAFE * (*dt) * pow( *errmax , ADAPTIVE_SHRINK ) ; 
    // set up a tolerance s.t del_temp doesn't go too crazy, only really used when starting guess is bad
    const double tol = 0.1 * ( *dt );
    
    *dt = ( 0. < del_temp ?				\
	    adaptfmax( del_temp , tol ) :		\
	    adaptfmin( del_temp , tol ) ) ;
            
    // Increment our counter ...
    counter ++ ;
      
    goto top ;
  }

  if( counter >= 20 ) return GLU_FAILURE ;

  return GLU_SUCCESS ;
}

// dirty macro
#define stept()							\
  setred( &WF );						\
  step_distance( lat , WF , delta_t , SM_TYPE , project ) ;	\
  new_plaq = get_plaq_th( lat , WF ) ;				\
  curr = (struct wfmeas*)malloc( sizeof( struct wfmeas ) ) ;	\
  update_meas_list( head , curr , WF.red ,			\
		    new_plaq , t = t+delta_t , delta_t ,	\
		    errmax*ADAPTIVE_EPS , lat ) ;		\
  flow_next = curr -> Gt ;					\
  wapprox = ( flow_next - flow ) * curr -> time / delta_t ;	\
  flow = flow_next ; head = curr ; count++ ; meas_count++ ;

// Adaptive stepsize version 
int 
flow_adaptive_RK3( struct site *__restrict lat , 
		   const size_t smiters ,
		   const int SIGN ,
		   const smearing_types SM_TYPE )
{  
  ////// USUAL STARTUP INFORMATION /////////
  print_GG_info( ) ;

  // adaptive error conserving
  const double ADAPTIVE_ERRCON = powl( 5./ADAPTIVE_SAFE , 1./ADAPTIVE_GROWTH ) ;
  // percentage to value we want for performing fine measurements
  const double FINETWIDDLE = 0.05 ;

  // adaptive factors for RK4, we are RK3 could be more lenient?
  fprintf( stdout , "[WFLOW] Adaptive Error :: %e \n" , ADAPTIVE_EPS ) ;
  fprintf( stdout , "[WFLOW] Adaptive ErrCon :: %f \n" , ADAPTIVE_ERRCON ) ;
  fprintf( stdout , "[WFLOW] Safety Factor :: %g \n" , ADAPTIVE_SAFE ) ;
  fprintf( stdout , "[WFLOW] Growth factor :: %g \n" , ADAPTIVE_GROWTH ) ;
  fprintf( stdout , "[WFLOW] Shrink factor :: %g \n\n" , ADAPTIVE_SHRINK ) ; 
  fprintf( stdout , "[WFLOW] Fine measurement %% :: %g \n" , FINETWIDDLE ) ;
  fprintf( stdout , "[WFLOW] Fine step is 0.01 of current flow time\n\n" ) ;

  //////////////////////////////////////////
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
  
  // allocate temps
  struct wflow_temps WF ;
  if( allocate_WF( &WF , GLU_FALSE , GLU_TRUE ) == GLU_FAILURE ) {
    goto memfree ;
  }

  const double inplaq = av_plaquette( lat ) ;
  int FLAG = GLU_SUCCESS ;
  
#pragma omp parallel
  {
    struct wfmeas *head = NULL , *curr ;

    double new_plaq = inplaq ;
    double t = 0.0 ;
    size_t i ;
    double delta_t = SIGN * Latt.sm_alpha[0] ;

    curr = (struct wfmeas*)malloc( sizeof( struct wfmeas ) ) ;
    update_meas_list( head , curr , WF.red ,
		      new_plaq , t , delta_t , 0.0 , lat ) ;
    head = curr ;
    double flow = 0.0 , flow_next = curr -> Gt ;
    double yscal = curr -> avplaq , wapprox = 0.0 ;
    double errmax = 10. ;
    size_t count = 0 , meas_count = 0 ;
    
    // loop up to smiters
  top :
    {
       #pragma omp barrier
    }
    errmax = 10. ;
 
    // perform two-step adaptive if delta_t is small and at some point we just switch to stout smearing at maximal rho?
    if( twostep_adaptive( lat , WF , &delta_t , &errmax , &new_plaq ,
			  t , yscal , SM_TYPE , project ) == GLU_FAILURE ) {
      FLAG = GLU_FAILURE ;
      count = smiters ;
    }
    
    // set up a scaling parameter to control the adaptation uses a first order finite difference def ...
    const double yscal_new = new_plaq ;
    yscal = 2.0 * yscal_new - yscal ;
    
    // overwrite lat .. 
    #pragma omp for private(i) collapse(2)
    for( i = 0 ; i < LVOLUME ; i++ ) {
      for( size_t mu = 0 ; mu < ND ; mu++ ) {
	equiv( lat[i].O[mu] , WF.lat_two[i].O[mu] ) ;
      }
    }

    t += delta_t ; // add one time step to the overall time
    
    // If we get stuck in a rut of updating by zero we leave
    if( fabs( delta_t ) < DBL_MIN ) {
      #pragma omp master
      {
	fprintf( stderr , "[WFLOW] No update made delta_t :: %1.5e \n"
		 "Leaving ... \n" , delta_t ) ;
      }
      FLAG = GLU_FAILURE ;
    }
    
    if( t > WFLOW_MEAS_START ) {
      curr = (struct wfmeas*)malloc( sizeof( struct wfmeas ) ) ;
      update_meas_list( head , curr , WF.red ,
			new_plaq , t , delta_t ,
			errmax*ADAPTIVE_EPS , lat ) ;
      flow_next = curr -> Gt ;
      head = curr ;
      #ifndef WFLOW_TIME_ONLY
      wapprox = flow != 0.0 ? ( flow_next - flow ) * curr -> time / delta_t : 0.0 ;
      #endif
      flow = flow_next ;
      meas_count++ ;
    }

    // perform fine measurements around T0_STOP for t_0 and W0_STOP for w_0
    #ifndef WFLOW_TIME_ONLY
    if( fabs( T0_STOP - flow ) <= ( T0_STOP * FINETWIDDLE ) ) {
      int past_the_post = 0 ;
      // step backwards 1/2 step I don't love this tbH
      if( flow > T0_STOP ) {
	delta_t = -delta_t/2 ;
	stept() ;
      }
      // step forwards
      delta_t = t/100. ;
    t0_top :
      {
         #pragma omp barrier
      }
      stept() ;      
      if( flow > T0_STOP ) { past_the_post++ ; }
      if( fabs( T0_STOP - flow ) <= ( T0_STOP * FINETWIDDLE ) && past_the_post < 3 ) goto t0_top ;
    }

    // perform some fine measurements around W0_STOP
    if( fabs( W0_STOP - wapprox ) <= ( W0_STOP * FINETWIDDLE ) ) {
      int past_the_post = 0 ;
      delta_t = t/100. ;
    w0_top :
      {
         #pragma omp barrier
      }
      stept() ;
      if( wapprox > W0_STOP ) { past_the_post++ ; }
      if( fabs( W0_STOP - wapprox ) <= ( W0_STOP * FINETWIDDLE ) && past_the_post < 3 ) goto w0_top ;
    }
    #endif
    
    // Increase the step size ...
    if( delta_t < 0.125 ) {
      if( errmax > ADAPTIVE_ERRCON ) {
	delta_t = ADAPTIVE_SAFE * delta_t * pow( errmax , ADAPTIVE_GROWTH ) ;
      } else {
	//delta_t = ADAPTIVE_SAFE * 5.0 * delta_t ;
	delta_t = ADAPTIVE_SAFE * 4.0 * delta_t ;
      }
    }
    count++ ;

    // while loop
    if( ( wapprox < W0_STOP*1.5 ) &&
	( count < smiters ) &&
	( t < WFLOW_TIME_STOP ) &&
	FLAG == GLU_SUCCESS ) goto top ; 

    // if we have flowed past the place that we are meant to stop at we
    // step back a little way
    if( t > WFLOW_TIME_STOP ) {
      delta_t = WFLOW_TIME_STOP - t ;
      stept() ;
    }

    // finally point to the right place in the list again
    curr = head ; 
    
    // compute the t_0 and w_0 scales from the measurements
#ifndef WFLOW_TIME_ONLY
    if( ( fabs( curr -> time - WFLOW_TIME_STOP ) > PREC_TOL ) && 
	meas_count > 0 ) {
      scaleset( curr , T0_STOP , W0_STOP , meas_count , GLU_TRUE ) ;
      scaleset( curr , T0_STOP , W0_STOP , meas_count , GLU_FALSE ) ;
    }
#endif
    
    // and free the list
    while( ( curr = head ) != NULL ) {
      head = head -> next ; 
      free( curr ) ;
    }
  }
  
 memfree :

  free_WF( &WF , GLU_FALSE , GLU_TRUE ) ;

  return FLAG ;
}
