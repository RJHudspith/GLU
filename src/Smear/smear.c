/*
    Copyright 2013 Renwick James Hudspith

    This file (smear.c) is part of GLU.

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
   @file smear.c 
   @brief wanted to include both spatial and all dimensional smearing types and we do
   #ND- generic
 */

#include "Mainfile.h"
#include "plaqs_links.h"
#include "projectors.h"
#include "staples.h"

// If we are using the dangerous smearing routines ...
#if ( defined FAST_SMEAR ) || NC > 6
  #include "random_config.h"
#endif

// if we have the symanzik one loop we must include a symanzik tadpole
#ifdef SYMANZIK_ONE_LOOP
  // tadpole improvement only used in symanzik one loop
  GLU_real improve = 1.0 ;
#endif

//////////////////////////////////////////////////////////
int
smear3D( struct site *__restrict lat , 
	 const size_t smiters , 
	 const smearing_types type )
{
#if ND < 3
  return GLU_SUCCESS ;
#else
  // successfully do nothing
  if( unlikely( smiters == 0 ) ) { return GLU_SUCCESS ; }

#if ND != 4 
  const GLU_real alpha1 = Latt.sm_alpha[0] / ( ( ND-1 ) * ( ND -2 ) ) ;
  const GLU_real one_min_a1 = ( 1.0 - Latt.sm_alpha[0] ) ;
#endif

  #ifdef TOP_VALUE
  double qtop_new , qtop_old = 0. ;
  #endif

  // allocate temporary lattice field
  struct sp_site *lat2 = NULL ;
  if( GLU_malloc( (void**)&lat2 , 16 , LCU * sizeof( struct sp_site ) ) != 0 ) {
    fprintf( stderr , "[SMEARING] field allocation failure\n" ) ;
    return GLU_FAILURE ;
  }

  // callback for the projections
  void (*project) ( GLU_complex smeared_link[ NCNC ] , 
		    GLU_complex staple[ NCNC ] , 
		    const GLU_complex link[ NCNC ] , 
		    const double smear_alpha , 	     
		    const double al ) ;

  switch( type ) {
  case SM_APE :
    project = project_APE ;
    break ;
  case SM_STOUT :
    project = project_STOUT ;
    break ;
  case SM_LOG :
    project = project_LOG ;
    break ;
  default :
    fprintf( stderr , "[SMEAR] Unrecognised type [ %d ] ... Leaving \n" , 
	     type ) ; 
    return GLU_FAILURE ;
  }

  size_t count = 0 ; 
  for( count = 1 ; count <= smiters ; count++ ) {
    #ifdef SYMANZIK_ONE_LOOP
    improve = av_plaquette( lat ) ;
    #endif
    /////////////////////
    //loop time slices
    /////////////////////
    size_t i , t ;
    for( t = 0 ; t < Latt.dims[ ND -1 ] ; t++ ) {
      const size_t slice = LCU * t ;
      #pragma omp parallel for private(i) SCHED
      PFOR( i = 0 ; i < LCU ; i++ ) {
	const size_t it = slice + i ; 
	size_t mu ;
	for( mu = 0 ; mu < ND - 1 ; mu++ ) {
	  GLU_complex stap[ NCNC ] ;
	  zero_mat( stap ) ;
          #ifdef IMPROVED_SMEARING
	  all_staples_improve( stap , lat , it , mu , ND - 1 , type ) ; 
          #else
	  all_staples( stap , lat , it , mu , ND - 1 , type ) ; 
          #endif
	  project( lat2[ i ].O[ mu ] , stap , lat[ it ].O[ mu ] , alpha1 , 
		   one_min_a1 ) ; 
	}
      }
	  
      #pragma omp parallel for private(i) 
      PFOR( i = 0 ; i < LCU ; i++ )  {
	const size_t it = slice + i ;
	memcpy( &lat[ it ] , &lat2[ i ] , sizeof( struct sp_site ) ) ;
      }
    }

    #ifdef TOP_VALUE
    if( count > TOP_VALUE ) {
      if( gauge_topological_meas( lat , &qtop_new , &qtop_old , count-1 ) 
	  == GLU_SUCCESS ) { break ; }
    }
    #endif
      
    #ifdef verbose
    print_smearing_obs( lat , type , count , GLU_TRUE ) ;
    #endif
  }

#ifndef verbose
  // -- the counter here as we stop at smiters not simters+1
  print_smearing_obs( lat , type , count-1 , GLU_TRUE ) ;
#endif

  // free our temporary lattice
  free( lat2 ) ;
  return GLU_SUCCESS ;
#endif
}

// General ALL-dimensional smearing routines ...
int
smear4D( struct site *__restrict lat ,
	 const size_t smiters , 
	 const smearing_types type )
{
  // successfully do nothing
  if( unlikely( smiters == 0 ) ) { return GLU_SUCCESS ; }

#if ND != 4 
  const GLU_real alpha1 = ND > 2 ? Latt.sm_alpha[0] / ( ( ND - 1 ) * ( ND - 2 ) ) : Latt.sm_alpha[0] ;
  const GLU_real one_min_a1 = ( 1.0 - Latt.sm_alpha[0] ) ;
#endif	
				     
  #ifdef TOP_VALUE
  double qtop_new , qtop_old = 0. ;
  #endif

  // callback for the projections
  void (*project) ( GLU_complex smeared_link[ NCNC ] , 
		    GLU_complex staple[ NCNC ] , 
		    const GLU_complex link[ NCNC ] , 
		    const double smear_alpha , 	     
		    const double al ) ;

  switch( type ) {
  case SM_APE :
    project = project_APE ;
    break ;
  case SM_STOUT :
    project = project_STOUT ;
    break ;
  case SM_LOG :
    project = project_LOG ;
    break ;
  default :
    fprintf( stderr , "[SMEAR] Unrecognised type [ %d ] ... Leaving \n" , 
	     type ) ; 
    return GLU_FAILURE ;
  }

  struct spt_site *lat2 = NULL , *lat3 = NULL , *lat4 = NULL ;
  if( GLU_malloc( (void**)&lat2 , 16 , LCU * sizeof( struct spt_site ) ) != 0 ) {
    printf( "[SMEARING] field allocation failure\n" ) ;
  }
  #ifdef IMPROVED_SMEARING
  if( GLU_malloc( (void**)&lat3 , 16 , 2 * LCU * sizeof( struct spt_site ) ) != 0 ||
      GLU_malloc( (void**)&lat4 , 16 , 2 * LCU * sizeof( struct spt_site ) ) != 0) {
    printf( "[SMEARING] field allocation failure\n" ) ;
  }
  #else
  if( GLU_malloc( (void**)&lat3 , 16 , LCU * sizeof( struct spt_site ) ) != 0 ||
      GLU_malloc( (void**)&lat4 , 16 , LCU * sizeof( struct spt_site ) ) != 0) {
    printf( "[SMEARING] field allocation failure\n" ) ;
  }
  #endif
  
  size_t count = 0 ;
  for( count = 1 ; count <= smiters ; count++ ) {

    #ifdef SYMANZIK_ONE_LOOP
    improve = av_plaquette( lat ) ;
    #endif

    //this bit initialises the calculation by working out the staples for the last time slice first
    #ifdef IMPROVED_SMEARING
    const size_t back = lat[ lat[0].back[ ND-1 ] ].back[ ND-1 ] ;
    #else
    const size_t back = lat[ 0 ].back[ ND - 1 ] ;
    #endif

    size_t i , t ;
    #pragma omp parallel for private(i) SCHED
    #ifdef IMPROVED_SMEARING
    PFOR( i = 0 ; i < 2*LCU ; i++ ) {
    #else
    PFOR( i = 0 ; i < LCU ; i++ ) {
    #endif
      const size_t bck = back + i ;
      size_t mu ;
      for( mu = 0 ; mu < ND ; mu++ ) { 
	GLU_complex stap[ NCNC ] ;
	zero_mat( stap ) ;
        #ifdef IMPROVED_SMEARING
	all_staples_improve( stap , lat , bck , mu , ND , type ) ;
        #else
	all_staples( stap , lat , bck , mu , ND , type ) ;
        #endif
	project( lat4[ i ].O[ mu ] , stap , lat[ bck ].O[ mu ] , alpha1 ,
		 one_min_a1 ) ; 
      }
    }
    #ifdef IMPROVED_SMEARING
    for( t = 0 ; t < Latt.dims[ ND - 1 ] - 2 ; t++ ) {
    #else
    for( t = 0 ; t < Latt.dims[ ND - 1 ] - 1 ; t++ ) {
    #endif
      const size_t slice = LCU * t ; 
      #pragma omp parallel for private(i) SCHED
      PFOR( i = 0 ; i < LCU ; i++ ) {
	const size_t it = slice + i ;
	size_t mu ;
	for( mu = 0 ; mu < ND ; mu++ ) {
	  GLU_complex stap[ NCNC ] ;
	  zero_mat( stap ) ;
          #ifdef IMPROVED_SMEARING
	  all_staples_improve( stap , lat , it , mu , ND , type ) ;
          #else
	  all_staples( stap , lat , it , mu , ND , type ) ;
          #endif
	  project( lat2[ i ].O[ mu ] , stap , lat[ it ].O[ mu ] , alpha1 ,
		   one_min_a1 ) ; 
	}
      }

      #ifdef IMPROVED_SMEARING
      const size_t bck = lat[ lat[ slice ].back[ ND-1 ] ].back[ ND-1 ] ;
      #else
      const size_t bck = lat[ slice ].back[ ND -1 ] ;
      #endif

      #pragma omp parallel for private(i)
      PFOR( i = 0 ; i < LCU ; i++ ) {
        #ifdef IMPROVED_SMEARING
	if( likely( t > 1 ) ) { 
	  register const size_t back = bck + i ;
	  memcpy( &lat[back] , &lat3[i] , sizeof( struct spt_site ) ) ;
	}
	#else
	if( likely( t != 0 ) ) { 
	  register const size_t back = bck + i ;
	  memcpy( &lat[back] , &lat3[i] , sizeof( struct spt_site ) ) ;
	}
	#endif

        #ifdef IMPROVED_SMEARING
	// put the evaluation two time slices before into the front half of lat3
	// and the evaluation one time slice before into the lower half
	memcpy( &lat3[i] , &lat3[i+LCU] , sizeof( struct spt_site ) ) ;
	memcpy( &lat3[i+LCU] , &lat2[i] , sizeof( struct spt_site ) ) ;
        #else
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

    #ifdef TOP_VALUE
    if( count > TOP_VALUE ) {
      if( gauge_topological_meas( lat , &qtop_new , &qtop_old , count-1 ) 
	  == GLU_SUCCESS ) { break ; }
    }
    #endif
 
    #ifdef verbose
    print_smearing_obs( lat , type , count , GLU_TRUE ) ;
    #endif
    // end of iterations loop
  }

#ifndef verbose
 print_smearing_obs( lat , type , count-1 , GLU_TRUE ) ;
#endif

  // free stuff !
  free( lat2 ) ; 
  free( lat3 ) ; 
  free( lat4 ) ;

  return GLU_SUCCESS ;
}
