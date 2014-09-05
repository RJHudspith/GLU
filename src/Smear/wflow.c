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

// multiplier for the RK4
static const double mnineOseventeen = -0.52941176470588235294 ; // -9.0/17.0

//#define WFLOW_TIME_ONLY

// 4D wilson flow algorithm this is no longer the one we use
// it is the simple Eulerian STOUT smearing method
void 
flow4d( struct site *__restrict lat , 
        const int smiters ,
	const int DIR ,
	const int SIGN ,
	const int SM_TYPE )
{
  /////////// Initial blurb ////////
  print_GG_info( SM_TYPE , EULER ) ;
  //////////////////////////////////

  struct spt_site *lat2 = malloc( LVOLUME * sizeof( struct spt_site ) ) ;

  // some numerical constants used in the integration ....
  const double delta_t = SIGN * Latt.sm_alpha[0] ;

  // initial measurement of GG and topology
  double t = 0. ; // this is the flow time ...
  double qtop_out = 0. , avplaq_out = 0. ;
  lattice_gmunu( lat , &qtop_out , &avplaq_out ) ;
  printf("[WFLOW] {t} %g {dt} %g {p} %1.10f {q} %1.10f {ttGG} 0. {W} 0. \n" ,
	 t , delta_t , avplaq_out , qtop_out ) ;

  // counters for the derivative ...
  double flow = 0. , flow_prev = 0. , flow_next = 0. ;
  int count = 0 ; 
  int leapfrog = GLU_FALSE ;

  for( count = 1 ; count <= smiters ; count++ ) {
    t = count * delta_t ;

    int i ;
    #pragma omp parallel for private(i) SCHED
    PFOR( i = 0 ; i < LVOLUME; i++ ) {
      int mu ;
      for( mu = 0 ; mu < DIR ; mu++ ) {
	GLU_complex staple[ NCNC ] = { } ;
        #ifdef IMPROVED_SMEARING
	all_staples_improve( staple , lat , i , mu , DIR , SM_TYPE ) ;
        #else
	all_staples( staple , lat , i , mu , DIR , SM_TYPE ) ;
        #endif
	switch( SM_TYPE )
	  {
	  case SM_LOG :
	    project_LOG_short( lat2[i].O[mu] , staple , 
			       lat[i].O[mu] , delta_t ) ;
	    break ;
	  default :
	    project_STOUT_short( lat2[i].O[mu] , staple , 
				 lat[i].O[mu] , delta_t ) ;
	    break ;
	  } 
      }
    }
    #pragma omp parallel for  private(i) 
    PFOR( i = 0 ; i < LVOLUME; i++ ) {
      memcpy( &lat[i] , &lat2[i] , sizeof( struct spt_site ) ) ; 
    }

    // calculate our "f" this is the wilson flow ...
    // we skip measurements all the way up to t = 1 as a flow time of less 
    // than the lattice spacing is meaningless.
    #ifndef verbose
    if( ( t >= MEAS_START ) ) {
    #endif
      if( leapfrog == GLU_TRUE )  {
	if( deriv_leapfrog( lat , &flow_prev , &flow , &flow_next , 
			    t , delta_t ) > WFLOW_STOP ) {
	  break ;
	}
      } else {
	deriv_euler( lat , &flow , &flow_next , t , delta_t ) ;
        #ifndef verbose
	flow = flow_next ;
        #endif
	leapfrog = GLU_TRUE ;
      }
    #ifndef verbose
    }
    #endif
  }

  #ifdef FAST_SMEAR
  latt_reunitU( lat ) ;
  printf( "\n[WFLOW] A final reunitarisation step to clean things up\n" ) ;
  printf( "[WFLOW] PLAQS :: {w} %1.15f  {s} %1.15f  {t} %1.15f \n" , 
	  av_plaquette( lat ) , s_plaq( lat ) , t_plaq( lat ) ) ; 
  #endif
  printf("[WFLOW] RES:: W0^2 %1.15f , W0 :: %1.15f \n" , t , sqrt( t ) ) ;
  free( lat2 ) ;
  return ;
}

// 4D wilson flow algorithm ...
void 
flow4d_RK_fast( struct site *__restrict lat , 
		const int smiters ,
		const int DIR ,
		const int SIGN ,
		const int SM_TYPE )
{
  /////// Initial information //////////
  print_GG_info( SM_TYPE , RK4_FAST ) ;
  //////////////////////////////////////

  // Set up the temporary fields ...
  struct spt_site_herm *Z = malloc( LVOLUME * sizeof ( struct spt_site_herm ) ) ;
  struct spt_site *lat2 = malloc( LVOLUME * sizeof ( struct spt_site ) ) ;

  const double delta_t = SIGN * Latt.sm_alpha[0] ;

  double qtop_out = 0. , avplaq_out = 0. ;
  double t = 0. ; // this is the flow time ...

  // initial output information
  lattice_gmunu( lat , &qtop_out , &avplaq_out ) ;
  printf("[WFLOW] {t} %g {dt} %g {p} %1.10f {q} %1.10f {ttGG} 0. {W} 0. \n" ,
	 t , delta_t , avplaq_out , qtop_out ) ;

  // counters for the derivative ...
  double flow = 0. , flow_prev = 0. , flow_next = 0. ;
  int count = 0 ; 
  int leapfrog = GLU_FALSE ;

  // forward or backward ?
  const double rk1 = mnineOseventeen * delta_t ;
  const double rk2 = delta_t ;
  const double rk3 = ( -delta_t ) ;  

  // perform the loop over smearing iterations ... We break at the stopping 
  // point anyway.
  for( count = 1 ; count <= smiters ; count++ ) {

    t = count * delta_t ; // add one time step to the overall time 

    step_distance( lat , lat2 , Z , rk1 , rk2 , rk3 , SM_TYPE ) ;

    // calculate our "f" this is the wilson flow ...
    // we skip measurements all the way up to t = 1 as a flow time of less 
    // than the lattice spacing is meaningless.
    #ifndef verbose
    if( ( t >= MEAS_START ) ) {
    #endif
      printf( "[WFLOW] " ) ;
      if( leapfrog == GLU_TRUE ) {
	if( deriv_leapfrog( lat , &flow_prev , &flow , &flow_next , 
			    t , delta_t ) > WFLOW_STOP ) {
	  break ;
	}
      } else  {
	deriv_euler( lat , &flow , &flow_next , t , delta_t ) ;
        #ifndef verbose
	flow = flow_next ;
        #endif
	leapfrog = GLU_TRUE ;
      }
    #ifndef verbose
    }
    #endif

    // If we are wilson-flowing to a specific time, we compute a negative flow-time correction
    if( t > TMEAS_STOP ) {
      const double delta_tcorr = TMEAS_STOP - t ; 
      t = TMEAS_STOP ;
      const double rk1 = mnineOseventeen * delta_tcorr ;
      const double rk2 = delta_tcorr ;
      const double rk3 = ( -delta_tcorr ) ;  
      step_distance( lat , lat2 , Z , rk1 , rk2 , rk3 , SM_TYPE ) ;
      printf( "[WFLOW] " ) ;
      deriv_euler( lat , &flow , &flow_next , t , delta_tcorr ) ;
      break ;
    }
    // end of the integration step ...
  }
  
  printf("[WFLOW] FLOW %d RES:: W0^2 %1.15f , W0 :: %1.15f \n" , 
	 Latt.flow , t , sqrt( t ) ) ;

  // free our fields
  free( Z ) ;
  free( lat2 ) ;
   
  return ;
}

// smaller memory footprint one, still quite cheap to use too
void 
flow4d_RK_slow( struct site *__restrict lat , 
		const int smiters ,
		const int DIR ,
		const int SIGN ,
		const int SM_TYPE )
{
  ///// USUAL STARTUP INFORMATION //////
  print_GG_info( SM_TYPE , RK4_SLOW ) ;
  //////////////////////////////////////

  struct spt_site_herm *Z = malloc( LVOLUME * sizeof ( struct spt_site_herm ) ) ;
  struct spt_site *lat2 = malloc( LCU * sizeof ( struct spt_site ) ) ;
#ifdef IMPROVED_SMEARING
  struct spt_site *lat3 = malloc( 2*LCU * sizeof ( struct spt_site ) ) ;
  struct spt_site *lat4 = malloc( 2*LCU * sizeof ( struct spt_site ) ) ;
#else
  struct spt_site *lat3 = malloc( LCU * sizeof ( struct spt_site ) ) ;
  struct spt_site *lat4 = malloc( LCU * sizeof ( struct spt_site ) ) ;
#endif

  const double delta_t = SIGN * Latt.sm_alpha[0] ;

  double t = 0. ; // this is the flow time ...
  double qtop_out = 0. , avplaq_out = 0. ;
  lattice_gmunu( lat , &qtop_out , &avplaq_out ) ;
  printf("[WFLOW] {t} %g {dt} %g {p} %1.10f {q} %1.10f {ttGG} 0. {W} 0. \n" ,
	 t , delta_t , avplaq_out , qtop_out ) ;
 
  // counters for the derivative ...
  double flow = 0. , flow_prev = 0. , flow_next = 0. ;
  int count = 0 ; 
  int leapfrog = GLU_FALSE ;

  // forward or backward ?
  const double rk1 = mnineOseventeen * delta_t;
  const double rk2 = delta_t ;
  const double rk3 = ( -delta_t ) ;  

  for( count = 1 ; count <= smiters ; count++ )  {
    t = count * delta_t ; // add one time step to the overall time 

    step_distance_memcheap( lat , lat2 , lat3 , lat4 , 
			    Z , rk1 , rk2 , rk3 , SM_TYPE ) ;

#ifndef WFLOW_TIME_ONLY
    // for the W0 scale ...
    #ifndef verbose
    if( ( t >= MEAS_START ) ) {
    #endif
      printf( "[WFLOW] " ) ;
      if( leapfrog == GLU_TRUE ) {
	if( deriv_leapfrog( lat , &flow_prev , &flow , &flow_next , 
			    t , delta_t ) > WFLOW_STOP )
	  break ;
      } else {
	deriv_euler( lat , &flow , &flow_next , t , delta_t ) ;
        #ifndef verbose
	flow = flow_next ;
        #endif
	leapfrog = GLU_TRUE ;
      }
    #ifndef verbose
    }
    #endif
#endif
    // If we are wilson-flowing to a specific time, we compute a negative 
    // flow-time correction
    if( t > TMEAS_STOP ) {
      const double delta_tcorr = TMEAS_STOP - t ; 
      t = TMEAS_STOP ;
      const double rk1 = mnineOseventeen * delta_tcorr ;
      const double rk2 = delta_tcorr ;
      const double rk3 = ( -delta_tcorr ) ;  
      step_distance_memcheap( lat , lat2 , lat3 , lat4 , 
			      Z , rk1 , rk2 , rk3 , SM_TYPE ) ;
      printf( "[WFLOW] " ) ;
      deriv_euler( lat , &flow , &flow_next , t , delta_tcorr ) ;
      break ;
    }
    // end of the integration step ...
  }
  
  printf("[WFLOW] Configuration :: %d || RES:: W0^2 %1.15f || W0 :: %1.15f \n" , 
	 Latt.flow , t , sqrt( t ) ) ;

  // free our fields
  free( Z ) ;
  free( lat2 ) ;
  free( lat3 ) ;
  free( lat4 ) ;
   
  return ;
}

// make sure we clean this up
#ifdef WFLOW_TIME_ONLY
  #undef WFLOW_TIME_ONLY
#endif
