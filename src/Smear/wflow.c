/*
    Copyright 2013 Renwick James Hudspith

    This file (wflow.c) is part of GLU.

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
   @file wflow.c
   @brief this includes the new code for the wilson flow routines.

   These have been checked against BMW's code.

   I offer two codes, like the HYP it
   includes memory wasteful and a memory
   cheaper version. The adaptive version, which I recommend you use
   can be found in adaptive.c
 */
#include "Mainfile.h"

#include "plaqs_links.h" // for the clover
#include "projectors.h"
#include "staples.h"
#include "wflowfuncs.h"
#ifdef FAST_SMEAR
  #include "random_config.h"
#endif

//#define WFLOW_TIME_ONLY

// 4D wilson flow algorithm ...
int
flow4d_RK_fast( struct site *__restrict lat , 
		const size_t smiters ,
		const size_t DIR ,
		const int SIGN ,
		const smearing_types SM_TYPE )
{
  /////// Initial information //////////
  print_GG_info( SM_TYPE , RK4_FAST ) ;
  //////////////////////////////////////

  // Set up the temporary fields ...
  struct spt_site_herm *Z = NULL ;
  struct spt_site *lat2 = NULL ;
  if( GLU_malloc( (void**)&Z , ALIGNMENT , LVOLUME * sizeof( struct spt_site_herm ) ) != 0 ||
      GLU_malloc( (void**)&lat2 , ALIGNMENT , LVOLUME * sizeof( struct spt_site ) ) != 0 ) {
    fprintf( stderr , "[WFLOW] temporary field allocation failure\n" ) ;
    return GLU_FAILURE ;
  }

  const double delta_t = SIGN * Latt.sm_alpha[0] ;

  struct wfmeas *head = NULL , *curr ;
  curr = (struct wfmeas*)malloc( sizeof( struct wfmeas ) ) ;

  // initial output information
  lattice_gmunu( lat , &(curr -> qtop) , &(curr -> avplaq) ) ;
  curr -> time = 0.0 ;
  fprintf( stdout , "[WFLOW] {t} %1.12f {dt} %g {p} %1.12f {q} %g {Gt} 0.\n" ,
	   curr -> time , delta_t , curr -> avplaq , curr -> qtop ) ;
  curr -> next = head ;
  head = curr ;

  // counters for the derivative ...
  double flow = 0. , flow_next = 0. ;
  size_t count = 0 ; 

  // forward or backward ?
  const double rk1 = -0.52941176470588235294 * delta_t ;
  const double rk2 = delta_t ;
  const double rk3 = ( -delta_t ) ;

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

  // perform the loop over smearing iterations ... We break at the stopping 
  // point anyway.
  for( count = 1 ; count <= smiters ; count++ ) {

    curr = (struct wfmeas*)malloc( sizeof( struct wfmeas ) ) ;

    curr -> time = count * delta_t ; // add one time step to the overall time 

    // flow forwards using the RK4
    step_distance( lat , lat2 , Z , rk1 , rk2 , rk3 , SM_TYPE , project ) ;

    // update the linked list
    curr -> Gt = curr -> time * curr -> time * 
      lattice_gmunu( lat , &(curr -> qtop) , &( curr->avplaq ) ) ;

    // set the flow
    flow_next = curr -> Gt ;
    fprintf( stdout , "[WFLOW] {t} %1.12f {dt} %g {p} %1.12f {q} %g {Gt} %1.12f\n" , 
	     curr -> time , delta_t  , curr -> avplaq , curr -> qtop , curr -> Gt ) ;

    // use a poor approximation of the derivative of the flow as a guide to stop
    if( ( ( flow_next - flow ) * curr -> time / delta_t ) > ( W0_STOP * 1.25 ) ) {
      curr -> next = head ;
      head = curr ;
      break ;
    }
    flow = flow_next ;

    curr -> next = head ;
    head = curr ;

    // If we are wilson-flowing to a specific time, we compute a negative 
    // flow-time correction
    if( curr -> time > TMEAS_STOP ) {
      curr = (struct wfmeas*)malloc( sizeof( struct wfmeas ) ) ;
      const double delta_tcorr = TMEAS_STOP - count * delta_t ; 
      curr -> time = TMEAS_STOP ;
      const double rk1 = -0.52941176470588235294 * delta_tcorr ;
      const double rk2 = delta_tcorr ;
      const double rk3 = ( -delta_tcorr ) ;  
      step_distance( lat , lat2 , Z , rk1 , rk2 , rk3 , SM_TYPE , project ) ;
      // update the linked list
      curr -> Gt = curr -> time * curr -> time * 
	lattice_gmunu( lat , &(curr -> qtop) , &( curr->avplaq ) ) ;
      // output details
      fprintf( stdout , "[WFLOW] {t} %1.12f {dt} %g {p} %1.12f {q} %g {Gt} %1.12f \n" , 
	       curr -> time , delta_tcorr , curr -> avplaq , curr -> qtop , curr -> Gt ) ;
      curr -> next = head ;
      head = curr ;
      break ;
    }
    // end of the integration step ...
  }
  curr = head ;  

  if( ( fabs( curr -> time - TMEAS_STOP ) > PREC_TOL ) &&
      count <= smiters ) {
    // set the scale at t_0 and w_0
    scaleset( curr , T0_STOP , W0_STOP , count ) ;
  }

  // free the list
  while( ( curr = head ) != NULL ) {
    head = head -> next ; 
    free( curr ) ;
  }

  // free our fields
  free( Z ) ;
  free( lat2 ) ;
   
  return GLU_SUCCESS ;
}

// smaller memory footprint one, still quite cheap to use too
int
flow4d_RK_slow( struct site *__restrict lat , 
		const size_t smiters ,
		const size_t DIR ,
		const int SIGN ,
		const smearing_types SM_TYPE )
{
  ///// USUAL STARTUP INFORMATION //////
  print_GG_info( SM_TYPE , RK4_SLOW ) ;
  //////////////////////////////////////

  struct spt_site_herm *Z = NULL ;
  struct spt_site *lat2 = NULL , *lat3 = NULL , *lat4 = NULL ;
  if( GLU_malloc( (void**)&Z , ALIGNMENT , LVOLUME * sizeof( struct spt_site_herm ) ) != 0 ||
      GLU_malloc( (void**)&lat2 , ALIGNMENT , LVOLUME * sizeof( struct spt_site ) ) != 0 ) {
    fprintf( stderr , "[WFLOW] temporary field allocation failure\n" ) ;
    return GLU_FAILURE ;
  }
#ifdef IMPROVED_SMEARING
  if( GLU_malloc( (void**)&lat3 , ALIGNMENT , 2 * LCU * sizeof( struct spt_site ) ) != 0 ||
      GLU_malloc( (void**)&lat4 , ALIGNMENT , 2 * LCU * sizeof( struct spt_site ) ) != 0 ) {
    fprintf( stderr , "[WFLOW] temporary field allocation failure\n" ) ;
    return GLU_FAILURE ;
  }
#else
  if( GLU_malloc( (void**)&lat3 , ALIGNMENT , LCU * sizeof( struct spt_site ) ) != 0 ||
      GLU_malloc( (void**)&lat4 ,  ALIGNMENT, LCU * sizeof( struct spt_site ) ) != 0 ) {
    fprintf( stderr , "[WFLOW] temporary field allocation failure\n" ) ;
    return GLU_FAILURE ;
  }
#endif

  const double delta_t = SIGN * Latt.sm_alpha[0] ;

  struct wfmeas *head = NULL , *curr ;
  curr = (struct wfmeas*)malloc( sizeof( struct wfmeas ) ) ;

  // initial output information
  lattice_gmunu( lat , &(curr -> qtop) , &(curr -> avplaq) ) ;
  curr -> time = 0.0 ;
  fprintf( stdout , "[WFLOW] {t} %1.12f {dt} %g {p} %1.12f {q} %g {ttGG} 0.\n" ,
	   curr -> time , delta_t , curr -> avplaq , curr -> qtop ) ;
  curr -> next = head ;
  head = curr ;

  // counters for the derivative ...
  double flow = 0. , flow_next = 0. ;
  size_t count = 0 ; 

  // forward or backward ?
  const double rk1 = -0.52941176470588235294 * delta_t;
  const double rk2 = delta_t ;
  const double rk3 = ( -delta_t ) ;  

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

  // loop iterations
  for( count = 1 ; count <= smiters ; count++ )  {
    curr = (struct wfmeas*)malloc( sizeof( struct wfmeas ) ) ;

    curr -> time = count * delta_t ; // add one time step to the overall time 

    step_distance_memcheap( lat , lat2 , lat3 , lat4 , 
			    Z , rk1 , rk2 , rk3 , SM_TYPE , project ) ;

    curr -> time = count * delta_t ; // add one time step to the overall time 

    // update the linked list
    curr -> Gt = curr -> time * curr -> time * 
      lattice_gmunu( lat , &(curr -> qtop) , &( curr->avplaq ) ) ;

    // set the flow
    flow_next = curr -> Gt ;
    fprintf( stdout , "[WFLOW] {t} %1.12f {dt} %g {p} %1.12f {q} %g {Gt} %1.12f \n" , 
	     curr -> time , delta_t , curr -> avplaq , curr -> qtop , curr -> Gt ) ;

    // use a poor approximation of the derivative of the flow as a guide to stop
    if( ( ( flow_next - flow ) * curr -> time / delta_t ) > ( W0_STOP * 1.25 ) ) {
      curr -> next = head ;
      head = curr ;
      break ;
    }
    flow = flow_next ;

    // update the list
    curr -> next = head ;
    head = curr ;

    // If we are wilson-flowing to a specific time, we compute a negative 
    // flow-time correction
    if( curr -> time > TMEAS_STOP ) {
      curr = (struct wfmeas*)malloc( sizeof( struct wfmeas ) ) ;
      const double delta_tcorr = TMEAS_STOP - count * delta_t ; 
      curr -> time = TMEAS_STOP ;
      const double rk1 = -0.52941176470588235294 * delta_tcorr ;
      const double rk2 = delta_tcorr ;
      const double rk3 = ( -delta_tcorr ) ;  
      step_distance_memcheap( lat , lat2 , lat3 , lat4 , 
			      Z , rk1 , rk2 , rk3 , SM_TYPE , project ) ;
      // update the linked list
      curr -> Gt = curr -> time * curr -> time * 
	lattice_gmunu( lat , &(curr -> qtop) , &( curr->avplaq ) ) ;
      // output details
      fprintf( stdout , "[WFLOW] {t} %1.12f {dt} %g {p} %1.12f {q} %g {Gt} %1.12f \n" , 
	       curr -> time , delta_tcorr , curr -> avplaq , curr -> qtop , curr -> Gt ) ;
      curr -> next = head ;
      head = curr ;
      break ;
    }
    // end of the integration step ...
  }
  curr = head ;  

  // set the scales t0 and w0 evaluated by our splines
  if( ( fabs( curr -> time - TMEAS_STOP ) > PREC_TOL ) && 
      count <= smiters ) {
    // set the scale at t_0 and w_0
    scaleset( curr , T0_STOP , W0_STOP , count ) ;
  }

  // free the list
  while( ( curr = head ) != NULL ) {
    head = head -> next ; 
    free( curr ) ;
  }

  // free our fields
  free( Z ) ;
  free( lat2 ) ;
  free( lat3 ) ;
  free( lat4 ) ;
   
  return GLU_SUCCESS ;
}

// make sure we clean this up
#ifdef WFLOW_TIME_ONLY
  #undef WFLOW_TIME_ONLY
#endif
