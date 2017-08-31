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

// 4D wilson flow algorithm ...
int
flow4d_RK_fast( struct site *__restrict lat , 
		const size_t smiters ,
		const int SIGN ,
		const smearing_types SM_TYPE )
{
  /////// Initial information //////////
  print_GG_info( ) ;
  //////////////////////////////////////

  // Set up the temporary fields ...
  struct s_site *Z = NULL , *lat2 = NULL ;
  if( ( lat2 = allocate_s_site( LVOLUME , ND , NCNC ) ) == NULL ||
      ( Z    = allocate_s_site( LVOLUME , ND , TRUE_HERM ) ) == NULL ) {
    fprintf( stderr , "[WFLOW] temporary field allocation failure\n" ) ;
    return GLU_FAILURE ;
  }

  const double delta_t = SIGN * Latt.sm_alpha[0] ;
  const double err_est = delta_t*delta_t*delta_t ;

  struct wfmeas *head = NULL , *curr ;
  curr = (struct wfmeas*)malloc( sizeof( struct wfmeas ) ) ;

  // initial output information
  lattice_gmunu( lat , &(curr -> qtop) , &(curr -> avplaq) ) ;
  curr -> time = 0.0 ;
  print_flow( curr , err_est , delta_t ) ;
  curr -> next = head ;
  head = curr ;

  // counters for the derivative ...
  double flow = 0. , flow_next = 0. ;
  size_t count , meas_count = 1 ; 

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

    // flow forwards using the RK3
    step_distance( lat , lat2 , Z , rk1 , rk2 , rk3 , SM_TYPE , project ) ;

    // update the linked list
    if( (count * delta_t) >= WFLOW_MEAS_START ) {

      curr = (struct wfmeas*)malloc( sizeof( struct wfmeas ) ) ;
      
      curr -> time = count * delta_t ; // add one time step to the overall time 
      
      curr -> Gt = curr -> time * curr -> time * 
	lattice_gmunu( lat , &(curr -> qtop) , &( curr->avplaq ) ) ;

      // set the flow
      flow_next = curr -> Gt ;

      print_flow( curr , err_est , delta_t ) ;

      // use a poor approximation of the derivative of the flow as a guide to stop
      #ifndef WFLOW_TIME_ONLY
      if( flow != 0.0 &&
	  ( ( flow_next - flow ) * curr -> time / delta_t ) > ( W0_STOP * 1.25 ) ) {
	curr -> next = head ;
	head = curr ;
	break ;
      }
      #endif
      flow = flow_next ;

      curr -> next = head ;
      head = curr ;

      meas_count++ ;
    }

    // If we are wilson-flowing to a specific time, we compute a negative 
    // flow-time correction
    if( (count * delta_t) > WFLOW_TIME_STOP ) {
      curr = (struct wfmeas*)malloc( sizeof( struct wfmeas ) ) ;
      const double delta_tcorr = WFLOW_TIME_STOP - count * delta_t ; 
      curr -> time = WFLOW_TIME_STOP ;
      const double rk1 = -0.52941176470588235294 * delta_tcorr ;
      const double rk2 = delta_tcorr ;
      const double rk3 = ( -delta_tcorr ) ;  
      step_distance( lat , lat2 , Z , rk1 , rk2 , rk3 , SM_TYPE , project ) ;
      // update the linked list
      curr -> Gt = curr -> time * curr -> time * 
	lattice_gmunu( lat , &(curr -> qtop) , &( curr->avplaq ) ) ;
      // output details
      print_flow( curr , err_est , delta_t ) ;
      curr -> next = head ;
      head = curr ;
      count++ ;
      meas_count++ ;
      break ;
    }
    // end of the integration step ...
  }
  curr = head ;

  fprintf( stdout , "[WFLOW] Iterations :: %zu || Measurements :: %zu\n" ,
	   count , meas_count ) ;

#ifndef WFLOW_TIME_ONLY
  if( ( fabs( curr -> time - WFLOW_TIME_STOP ) > PREC_TOL ) &&
      count <= smiters &&
      meas_count > 0 ) {
    // set the scale at t_0 and w_0
    scaleset( curr , T0_STOP , W0_STOP , meas_count ) ;
  }
#endif

  // free the list
  while( ( curr = head ) != NULL ) {
    head = head -> next ; 
    free( curr ) ;
  }

  // free our fields
  free_s_site( lat2 , LVOLUME , ND , NCNC ) ;
  free_s_site( Z , LVOLUME , ND , TRUE_HERM ) ;
   
  return GLU_SUCCESS ;
}

// smaller memory footprint one, still quite cheap to use too
int
flow4d_RK_slow( struct site *__restrict lat , 
		const size_t smiters ,
		const int SIGN ,
		const smearing_types SM_TYPE )
{
  ///// USUAL STARTUP INFORMATION //////
  print_GG_info( ) ;
  //////////////////////////////////////
  struct s_site *Z = NULL , *lat2 = NULL , *lat3 = NULL , *lat4 = NULL ;
  if( ( lat2 = allocate_s_site( LVOLUME , ND , NCNC ) ) == NULL || 
      #ifdef IMPROVED_SMEARING
      ( lat3 = allocate_s_site( 2*LCU , ND , NCNC ) ) == NULL || 
      ( lat4 = allocate_s_site( 2*LCU , ND , NCNC ) ) == NULL || 
      #else
      ( lat3 = allocate_s_site( LCU , ND , NCNC ) ) == NULL || 
      ( lat4 = allocate_s_site( LCU , ND , NCNC ) ) == NULL || 
      #endif
      ( Z = allocate_s_site( LVOLUME , ND , TRUE_HERM ) ) == NULL ) {
    fprintf( stderr , "[WFLOW] temporary field allocation failure\n" ) ;
    return GLU_FAILURE ;
  }

  const double delta_t = SIGN * Latt.sm_alpha[0] ;
  const double err_est = delta_t * delta_t * delta_t ;

  struct wfmeas *head = NULL , *curr ;
  curr = (struct wfmeas*)malloc( sizeof( struct wfmeas ) ) ;

  // initial output information
  lattice_gmunu( lat , &(curr -> qtop) , &(curr -> avplaq) ) ;
  curr -> time = 0.0 ;
  print_flow( curr , err_est , delta_t ) ;
  curr -> next = head ;
  head = curr ;

  // counters for the derivative ...
  double flow = 0. , flow_next = 0. ;
  size_t count = 0 , meas_count = 1 ; 

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

    step_distance_memcheap( lat , lat2 , lat3 , lat4 , 
			    Z , rk1 , rk2 , rk3 , SM_TYPE , project ) ;

    if( (count * delta_t) >= WFLOW_MEAS_START ) {

      curr = (struct wfmeas*)malloc( sizeof( struct wfmeas ) ) ;
      curr -> time = count * delta_t ; // add one time step to the overall time
    
      // update the linked list
      curr -> Gt = curr -> time * curr -> time * 
	lattice_gmunu( lat , &(curr -> qtop) , &( curr->avplaq ) ) ;

      // set the flow
      flow_next = curr -> Gt ;
      print_flow( curr , err_est , delta_t ) ;
 
      // use a poor approximation of the derivative of the flow as a guide to stop
      #ifndef WFLOW_TIME_ONLY
      if( flow != 0.0 &&
	  ( ( flow_next - flow ) * curr -> time / delta_t ) > ( W0_STOP * 1.25 ) ) {
	curr -> next = head ;
	head = curr ;
	break ;
      }
      #endif
      flow = flow_next ;

      // update the list
      curr -> next = head ;
      head = curr ;

      meas_count++ ;
    }
    
    // If we are wilson-flowing to a specific time, we compute a negative 
    // flow-time correction
    if( (count * delta_t) > WFLOW_TIME_STOP ) {
      curr = (struct wfmeas*)malloc( sizeof( struct wfmeas ) ) ;
      const double delta_tcorr = WFLOW_TIME_STOP - count * delta_t ; 
      curr -> time = WFLOW_TIME_STOP ;
      const double rk1 = -0.52941176470588235294 * delta_tcorr ;
      const double rk2 = delta_tcorr ;
      const double rk3 = ( -delta_tcorr ) ;  
      step_distance_memcheap( lat , lat2 , lat3 , lat4 , 
			      Z , rk1 , rk2 , rk3 , SM_TYPE , project ) ;
      // update the linked list
      curr -> Gt = curr -> time * curr -> time * 
	lattice_gmunu( lat , &(curr -> qtop) , &( curr->avplaq ) ) ;
      // output details
      print_flow( curr , err_est , delta_t ) ;
      curr -> next = head ;
      head = curr ;
      meas_count++ ;
      count++ ;
      break ;
    }
    // end of the integration step ...
  }
  curr = head ;

  fprintf( stdout , "[WFLOW] Iterations :: %zu || Measurements :: %zu\n" ,
	   count , meas_count ) ;

  // set the scales t0 and w0 evaluated by our splines
#ifndef WFLOW_TIME_ONLY
  if( ( fabs( curr -> time - WFLOW_TIME_STOP ) > PREC_TOL ) && 
      count <= smiters &&
      meas_count > 0 ) {
    // set the scale at t_0 and w_0
    scaleset( curr , T0_STOP , W0_STOP , meas_count ) ;
  }
#endif

  // free the list
  while( ( curr = head ) != NULL ) {
    head = head -> next ; 
    free( curr ) ;
  }

  // free our fields
  free_s_site( lat2 , LVOLUME , ND , NCNC ) ;
  free_s_site( Z , LVOLUME , ND , TRUE_HERM ) ;
#ifdef IMPROVED_SMEARING
  free_s_site( lat3 , 2*LCU , ND , NCNC ) ;
  free_s_site( lat4 , 2*LCU , ND , NCNC ) ;
#else
  free_s_site( lat3 , LCU , ND , NCNC ) ;
  free_s_site( lat4 , LCU , ND , NCNC ) ;
#endif
   
  return GLU_SUCCESS ;
}
