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
  struct s_site *lat2 = NULL ;
  if( ( lat2 = allocate_s_site( LCU , (ND-1) , NCNC ) ) == NULL ) {
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
      #pragma omp parallel for private(i)
      PFOR( i = 0 ; i < LCU ; i++ ) {
	const size_t it = slice + i ; 
	size_t mu ;
	for( mu = 0 ; mu < ND - 1 ; mu++ ) {
	  GLU_complex stap[ NCNC ] GLUalign ;
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
	size_t mu ;
	for( mu = 0 ; mu < (ND-1) ; mu++ ) {
	  equiv( lat[it].O[mu] , lat2[i].O[mu] ) ;
	}
      }
    }

    #ifdef TOP_VALUE
    if( count > TOP_VALUE ) {
      if( gauge_topological_meas( lat , &qtop_new , &qtop_old , count-1 ) 
	  == GLU_SUCCESS ) { break ; }
    }
    #endif
      
    #ifdef verbose
    print_smearing_obs( lat , count ) ;
    #endif
  }

#ifndef verbose
  // -- the counter here as we stop at smiters not simters+1
  print_smearing_obs( lat , count-1 ) ;
#endif

  // free our temporary lattice
  free_s_site( lat2 , LCU , (ND-1) , NCNC ) ;

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

  struct s_site *lat2 = NULL , *lat3 = NULL , *lat4 = NULL ;
  if( ( lat2 = allocate_s_site( LCU , ND , NCNC ) ) == NULL ) {
    fprintf( stderr , "[SMEARING] field allocation failure\n" ) ;
    return GLU_FAILURE ;
  }
#ifdef IMPROVED_SMEARING
  if( ( lat3 = allocate_s_site( 2*LCU , ND , NCNC ) ) == NULL || 
      ( lat4 = allocate_s_site( 2*LCU , ND , NCNC ) ) == NULL ) {
    fprintf( stderr , "[SMEARING] field allocation failure\n" ) ;
    return GLU_FAILURE ;
  }
#else
  if( ( lat3 = allocate_s_site( LCU , ND , NCNC ) ) == NULL || 
      ( lat4 = allocate_s_site( LCU , ND , NCNC ) ) == NULL ) {
    fprintf( stderr , "[SMEARING] field allocation failure\n" ) ;
    return GLU_FAILURE ;
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
      register size_t mu ;
      for( mu = 0 ; mu < ND ; mu++ ) { 
	GLU_complex stap[ NCNC ] GLUalign ;
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
	register size_t mu ;
	for( mu = 0 ; mu < ND ; mu++ ) {
	  #ifdef IMPROVED_SMEARING
	  if( t > 1 ) {
	    register const size_t back = bck + i ;
	    equiv( lat[back].O[mu] , lat3[i].O[mu] ) ;
	  }
	  equiv( lat3[i].O[mu]     , lat3[i+LCU].O[mu] ) ;
	  equiv( lat3[i+LCU].O[mu] , lat2[i].O[mu] ) ;
	  #else
	  if( t != 0 ) {
	    register const size_t back = bck + i ;
	    equiv( lat[back].O[mu] , lat3[i].O[mu] ) ;
	  }
	  equiv( lat3[i].O[mu] , lat2[i].O[mu] ) ;
          #endif
	}
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
      size_t mu ;
      for( mu = 0 ; mu < ND ; mu++ ) {
        #ifdef IMPROVED_SMEARING
	equiv( lat[behind2+i].O[mu] , lat3[i].O[mu] ) ;
	equiv( lat[back].O[mu]      , lat3[i].O[mu] ) ;
	equiv( lat[it].O[mu]        , lat4[i].O[mu] ) ;
	equiv( lat[it+LCU].O[mu]    , lat4[i+LCU].O[mu] ) ;
        #else
	equiv( lat[back].O[mu] , lat3[i].O[mu] ) ;
	equiv( lat[it].O[mu]   , lat4[i].O[mu] ) ;
        #endif
      }
    }

    #ifdef TOP_VALUE
    if( count > TOP_VALUE ) {
      if( gauge_topological_meas( lat , &qtop_new , &qtop_old , count-1 ) 
	  == GLU_SUCCESS ) { break ; }
    }
    #endif
 
    #ifdef verbose
    print_smearing_obs( lat , count ) ;
    #endif
    // end of iterations loop
  }

#ifndef verbose
 print_smearing_obs( lat , count-1 ) ;
#endif

  // free stuff !
 free_s_site( lat2 , LCU , ND , NCNC ) ;
#ifdef IMPROVED_SMEARING
 free_s_site( lat3 , 2*LCU , ND , NCNC ) ;
 free_s_site( lat4 , 2*LCU , ND , NCNC ) ;
#else
 free_s_site( lat3 , LCU , ND , NCNC ) ;
 free_s_site( lat4 , LCU , ND , NCNC ) ;
#endif

  return GLU_SUCCESS ;
}
