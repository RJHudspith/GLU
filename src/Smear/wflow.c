/*
    Copyright 2013-2018 Renwick James Hudspith

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

// fixed stepsize wilson flow algorithm ...
int
flow_RK3( struct site *__restrict lat , 
	  const size_t smiters ,
	  const int SIGN ,
	  const smearing_types SM_TYPE ,
	  const GLU_bool memcheap )
{
  /////// Initial information //////////
  print_GG_info( ) ;
  //////////////////////////////////////

  // Set up the temporary fields ...
  struct wflow_temps WF ;
  if( allocate_WF( &WF , memcheap , GLU_FALSE ) == GLU_FAILURE ) goto memfree ;

  void (*f) ( struct site *__restrict lat ,
	      struct wflow_temps WF ,
	      const double delta_t ,
	      const smearing_types SM_TYPE ,
	      void (*project)( GLU_complex log[ NCNC ] , 
			       GLU_complex *__restrict staple , 
			       const GLU_complex link[ NCNC ] , 
			       const double smear_alpha ) ) ;
  if( memcheap == GLU_TRUE ) {
    f = step_distance_memcheap ;
  } else {
    f = step_distance ;
  }

  const double delta_t = SIGN * Latt.sm_alpha[0] ;
  const double err_est = delta_t*delta_t*delta_t ;

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

  const double inplaq = av_plaquette( lat ) ;

#pragma omp parallel
  {
    // counters for the derivative ...
    double flow = 0. , flow_next = 0. , t = 0.0 , wapprox = 0.0 , new_plaq = inplaq ;
    size_t count = 0 , meas_count = 0 ;
    
    struct wfmeas *head = NULL , *curr ;
    curr = (struct wfmeas*)malloc( sizeof( struct wfmeas ) ) ;
    update_meas_list( head , curr , WF.red ,
		      new_plaq , t , delta_t ,
		      err_est , lat ) ;
    flow_next = curr -> Gt ;
    head = curr ;

    // perform the loop over smearing iterations ... We break at the stopping 
    // point anyway.
  top :
    #pragma omp barrier

    count++ ;
    
    t = count * delta_t ;
    
    // flow forwards using the RK3
    f( lat , WF , delta_t , SM_TYPE , project ) ;
    
    // update the linked list
    if( t >= WFLOW_MEAS_START ) {
      curr = (struct wfmeas*)malloc( sizeof( struct wfmeas ) ) ;
      // compute plaquette
      av_plaquette_th( WF.red , lat ) ;
      size_t k ;
      new_plaq = 0.0 ;
      for( k = 0 ; k < Latt.Nthreads ; k++ ) {
	new_plaq += WF.red[ 3 + CLINE*k ] ;
      }
      new_plaq /= ( NC * ( ND-1 ) * ( ND - 2 ) * LVOLUME ) ;
      update_meas_list( head , curr , WF.red ,
			new_plaq , t , delta_t ,
			err_est , lat ) ;
      flow_next = curr -> Gt ;
      #ifndef WFLOW_TIME_ONLY
      wapprox = ( flow_next - flow ) * curr -> time / delta_t ;
      #endif
      flow = flow_next ;
      head = curr ;
      meas_count++ ;
    }

    if( ( wapprox < W0_STOP*1.5 ) &&
	( count < smiters ) &&
	( t < WFLOW_TIME_STOP ) ) goto top ; 
    
    // If we are wilson-flowing to a specific time, we compute a negative 
    // flow-time correction
    if( t > WFLOW_TIME_STOP ) {
      const double delta_tcorr = WFLOW_TIME_STOP - count * delta_t ; 
      f( lat , WF , delta_tcorr , SM_TYPE , project ) ;
      // update the linked list
      curr = (struct wfmeas*)malloc( sizeof( struct wfmeas ) ) ;
      update_meas_list( head , curr , WF.red ,
			new_plaq , t , delta_t ,
			err_est , lat ) ;
      head = curr ;
      count++ ;
      meas_count++ ;
    }
    
    // end of the integration step(s) ...
    curr = head ;

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
  }

 memfree : 
  
  free_WF( &WF , memcheap , GLU_FALSE ) ;
   
  return GLU_SUCCESS ;
}
