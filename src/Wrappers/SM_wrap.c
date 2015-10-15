/*
    Copyright 2013 Renwick James Hudspith

    This file (SM_wrap.c) is part of GLU.

    GLU is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GLU is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more depsilonils.

    You should have received a copy of the GNU General Public License
    along with GLU.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   @file SM_wrap.c
   @brief wrapper for the smearing methods available

  Wrapper for the smearing methods available,
  APE , STOUT , LOG , HYP , HEX , HYL and the
  Wilson flow routines ...

   Which are enumerated in #smearing_types
 */
#include "Mainfile.h"

// The smearing headers ...
#include "4D_fast.h"        // excessively uses temporaries for precomuputation
#include "adaptive_flow.h"  // adaptive wilson flow routine
#include "GLU_memcheck.h"   // memory checking
#include "GLU_timer.h"      // timing functions
#include "HYP.h"            // 3D (spatial) smearing functions
#include "HYP_4D.h"         // 4D (all directional) smearing functions
#include "ND_generic_HYP.h" // recursive ND-dimensional blocked smearing
#include "smear.h"          // standard smearing
#include "wflow.h"          // fixed-epsilon wilson flow

// Selects the correct code to run ...
static void
hyp_chooser( lat , SMINFO , meminfo )
     struct site *__restrict lat ;
     const struct sm_info SMINFO ;
     const int meminfo ;
{
  // ADAPTIVE RK4 method of integrating the flow equation
  if( SMINFO.type == SM_ADAPTWFLOW_LOG || 
      SMINFO.type == SM_ADAPTWFLOW_STOUT ) {
    int GENTYPE = SM_STOUT ;
    if( SMINFO.type == SM_ADAPTWFLOW_LOG ) {
      GENTYPE = SM_LOG ;
    }
    flow4d_adaptive_RK( lat , SMINFO.smiters , ND , 
			GLU_FORWARD , GENTYPE ) ;
    return ;
  }
  // RK4 method of integrating the flow equation
  if( SMINFO.type == SM_WFLOW_LOG || 
      SMINFO.type == SM_WFLOW_STOUT )
    {
      int GENTYPE = SM_STOUT ;
      if( SMINFO.type == SM_WFLOW_LOG ) {
	GENTYPE = SM_LOG ;
      }
      switch( meminfo )
	{
	case FAST :
	  flow4d_RK_fast( lat , SMINFO.smiters , ND , 
			  GLU_FORWARD , GENTYPE ) ;
	  return ;
	case MODERATE :
	  flow4d_RK_slow( lat , SMINFO.smiters , ND , 
			  GLU_FORWARD , GENTYPE ) ;
	  return ;
	}
    }
  // The original decision for the smearing methods
  const int GENTYPE = ( SMINFO.type - 1 ) % 3 + 1 ;

  if( SMINFO.type == SM_APE || SMINFO.type == SM_STOUT 
      || SMINFO.type == SM_LOG ) {
    if( SMINFO.dir == SPATIAL_LINKS_ONLY ) {
      smear3D( lat , SMINFO.smiters , GENTYPE ) ;
    } else {
      smear4D( lat , SMINFO.smiters , GENTYPE ) ;
    }
  } else {
    #if ND != 4
    HYsmearND( lat , SMINFO.smiters , GENTYPE , SMINFO.dir ) ;
    #else
    if( SMINFO.dir == SPATIAL_LINKS_ONLY ) {
      HYPSLsmear3D( lat , SMINFO.smiters , GENTYPE ) ; 
      return ;
    } else {
      switch( meminfo )
	{
	case FAST :
	  HYPSLsmear4D_expensive( lat , SMINFO.smiters , GENTYPE ) ;
	  return ;
	case MODERATE :
	  HYPSLsmear4D( lat , SMINFO.smiters , GENTYPE ) ;  
	  return ;
	case SLOW :
	  HYsmearND( lat , SMINFO.smiters , GENTYPE , SMINFO.dir ) ;
	  return ;
	}
    }
    #endif
  }
  return ;
}

// is it smear or wilson flow?
static GLU_bool
smear_or_wflow( int type )
{
  if( type == SM_WFLOW_LOG || 
      type == SM_ADAPTWFLOW_LOG ||
      type == SM_WFLOW_STOUT || 
      type == SM_ADAPTWFLOW_STOUT )  {
    printf( "[WFLOW] " ) ;
    return GLU_FALSE ;
  } else { 
    printf( "[SMEAR] " ) ; 
    return GLU_TRUE ;
  }
}

// write to stdout the smearing information
static void
print_smearinfo( const struct sm_info SMINFO ) 
{
  smear_or_wflow( SMINFO.type ) ;
  if( SMINFO.dir == SPATIAL_LINKS_ONLY ) {
    printf( "Spatial SU(%d) smearing\n" , NC ) ;
  } else {
    printf( "All directional SU(%d) smearing\n" , NC ) ;
  }
  #ifdef FAST_SMEAR
  smear_or_wflow( SMINFO.type ) ;
  printf( "Numerically unstable, but fast, routines being used\n" ) ;
  #endif
  smear_or_wflow( SMINFO.type ) ;
  switch( SMINFO.type ) 
    {
    case SM_APE : 
      #ifdef N_APE
      printf( "APE smearing (n-APE determinant-rescaled projection)\n" ) ; 
      #else
      printf( "APE smearing (Cabbibo-Marinari updates) \n" ) ; 
      #endif
      break ;
    case SM_LOG : printf( "Log smearing\n" ) ; break ;
    case SM_STOUT : printf( "STOUT smearing\n" ) ; break ;
    case SM_HYP : 
      #ifdef N_APE
      printf( "HYP smearing (n-APE det rescaled projection)\n" ) ; 
      #else
      printf( "HYP smearing (Cabbibo-Marinari updates) \n" ) ; 
      #endif
      break ;
    case SM_HEX : printf( "HEX smearing\n" ) ; break ;
    case SM_HYL : printf( "HYL smearing\n" ) ; break ;
    case SM_WFLOW_LOG : printf( "RK3-Wilson flow Log\n" ) ; break ;
    case SM_WFLOW_STOUT : printf( "RK3-Wilson flow STOUT\n" ) ; break ;
    case SM_ADAPTWFLOW_LOG : printf( "Adaptive RK3-Wilson flow Log\n" ) ; break ;
    case SM_ADAPTWFLOW_STOUT : printf( "Adaptive R3-Wilson flow STOUT\n" ) ; break ;
    default : printf( "No smearing\n" ) ;
    }
  // wilson flow information about the clover terms and what have you
  if( smear_or_wflow( SMINFO.type ) == GLU_FALSE ) {
  #ifdef CLOVER_IMPROVE
    printf( "O(a^4) Clover Term being used ... \n" ) ;
    #ifndef NK5
    printf( "k5 value :: %f \n" , k5 ) ; 
    #endif
  #elif defined PLAQUETTE_FMUNU
    printf( "Plaquette O(a) Clover Term being used ... \n" ) ;
  #else
    printf( "O(a^2) Clover Term being used ... \n" ) ;
  #endif
    smear_or_wflow( SMINFO.type ) ;
  #ifdef ANTIHERMITIAN
    printf( "Traced, Antihermitian (BMW) field strength tensor def ... \n" ) ;
  #elif defined CLOVER_LOG_DEF
    printf( "Exact Logarithm used in field strength tensor def ... \n" ) ;
  #else
    printf( "Traceless hermitian field strength tensor def ... \n" ) ;
  #endif
  } else {
    printf( "Smearing Alpha(s) used:\n" ) ;
    if( SMINFO.type == SM_HYP || SMINFO.type == SM_HEX || SMINFO.type == SM_HYL ) {
      int i ;
      for( i = 0 ; i < ND-1 ; i++ ) {
	printf( "             (alpha%d) %f\n" , i+1 , Latt.sm_alpha[i]/( ( ND-i-1 ) * ( ND-2 ) ) ) ;
      }
    } else { 
      printf( "             (alpha) %f \n" , Latt.sm_alpha[0]/(2.*(ND-1)) ) ; 
    }
  }
  // rectangle improvement is not implemented for hypercubic smearing ...
#ifdef IMPROVED_SMEARING
  smear_or_wflow( SMINFO.type ) ;
  #ifdef SYMANZIK
    printf( "SYMANZIK Over-improvement factor :: %f \n" , epsilon ) ;
  #elif defined IWASAKI
    printf( "IWASAKI Over-improvement factor :: %f \n" , epsilon ) ;
  #elif defined DBW2
    printf( "DBW2 Over-improvement factor :: %f \n" , epsilon ) ;
  #elif defined SYMANZIK_ONE_LOOP
    printf( "SYMANZIK ONE LOOP Over-improvement factor :: %f \n" , epsilon ) ;
  #else
    printf( "Rectangle coefficients ( %f , %f ) :: Over improvement %f\n" ,
	    IWA_WEIGHT1 , IWA_WEIGHT2 , epsilon ) ;
  #endif
#endif
  smear_or_wflow( SMINFO.type ) ;
  printf( "With a MAXIMUM of %d iterations \n\n" , SMINFO.smiters ) ;
  return ;
}

// the wrapper that decides the routine 5to be called ...
short int
SM_wrap_struct( struct site *__restrict lat ,
		const struct sm_info SMINFO )
{
  if( SMINFO.smiters == 0 ) { return GLU_SUCCESS ; }
  if( SMINFO.type == SM_NOSMEARING ) { return GLU_SUCCESS ; }
  // are we looking for an integer value for the topological charge?
#ifdef TOP_VALUE
  printf( "[SMEAR] Iterating a search for the topological charge commencing after %d iterations \n\n" , TOP_VALUE ) ;
#endif
  start_timer( ) ;
  print_smearinfo( SMINFO ) ;
  if( SMINFO.type == SM_HYP || 
      SMINFO.type == SM_HEX || 
      SMINFO.type == SM_HYL )  {
    const int meminfo = have_memory_hyp( SMINFO ) ;
    hyp_chooser( lat , SMINFO , meminfo ) ;
  } else if( SMINFO.type == SM_WFLOW_LOG || 
	     SMINFO.type == SM_ADAPTWFLOW_LOG ||
	     SMINFO.type == SM_WFLOW_STOUT || 
	     SMINFO.type == SM_ADAPTWFLOW_STOUT )  {
    const int meminfo = have_memory_wf( SMINFO ) ;
    if( meminfo == GLU_FAILURE ) { return GLU_FAILURE ; }
    hyp_chooser( lat , SMINFO , meminfo ) ;
  } else {
    const int meminfo = 0 ;
    hyp_chooser( lat , SMINFO , meminfo ) ;
  }
  print_time( ) ;
  return GLU_SUCCESS ;
}
