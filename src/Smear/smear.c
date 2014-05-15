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
   
   ND- generic
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
void 
smear3D( struct site *__restrict lat , 
	 const int smiters , 
	 const smearing_types type )
{
#if ND < 3
  return ;
#else
  if( unlikely( smiters == 0 ) ) { return ; }

  // check the type
  if( ( type != SM_APE ) && ( type != SM_STOUT ) && ( type != SM_LOG ) ) {
    printf( "[SMEAR] Unrecognised type [ %d ] ... Leaving \n" , type ) ; 
    return ; 
  }

#if ND != 4 
  const GLU_real alpha1 = Latt.sm_alpha[0] / ( ( ND-1 ) * ( ND -2 ) ) ;
  const GLU_real one_min_a1 = ( 1.0 - Latt.sm_alpha[0] ) ;
#endif

  #ifdef TOP_VALUE
  double qtop_new , qtop_old = 0. ;
  #endif

  struct sp_site *lat2 = malloc( LCU * sizeof( struct sp_site ) ) ; 
 
  int count = 0 ; 
  for( count = 1 ; count <= smiters ; count++ ) {
    #ifdef SYMANZIK_ONE_LOOP
    improve = av_plaquette( lat ) ;
    #endif
    /////////////////////
    //loop time slices
    /////////////////////
    int t ;
    for( t = 0 ; t < Latt.dims[ ND -1 ] ; t++ ) {
      const int slice = LCU * t ;
      int i ;
	  
      #pragma omp parallel for private(i) SCHED
      PFOR( i = 0 ; i < LCU ; i++ ) {
	const int it = slice + i ; 
	int mu ;

	for( mu = 0 ; mu < ND - 1 ; mu++ ) {
	  GLU_complex stap[ NCNC ] = { } ;
          #ifdef IMPROVED_SMEARING
	  all_staples_improve( stap , lat , it , mu , ND - 1 , type ) ; 
          #else
	  all_staples( stap , lat , it , mu , ND - 1 , type ) ; 
          #endif
		  
	  switch( type )
	    {
	    case SM_APE:
	      project_APE( lat2[ i ].O[ mu ] , stap , 
			   lat[ it ].O[ mu ] , alpha1 , 
			   one_min_a1 ) ; 
	      break ; 
	    case SM_STOUT:
	      project_STOUT_short( lat2[ i ].O[ mu ] , stap , 
				   lat[ it ].O[ mu ] , alpha1 ) ; 
	      break ; 
	    case SM_LOG:
	      project_LOG_short( lat2[ i ].O[ mu ] , stap , 
				 lat[ it ].O[ mu ] , alpha1 ) ; 
	      break ; 
	    default : break ;
	    }
	}
      }
	  
      #pragma omp parallel for private(i) 
      PFOR( i = 0 ; i < LCU ; i++ )  {
	const int it = slice + i ;
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

// need to clean up the lattice fields
  #ifdef FAST_SMEAR
  if( type == SM_STOUT ) {
    latt_reunitU( lat ) ;
    printf( "\n[SMEAR] A final reunitarisation step to clean things up\n" ) ;
    print_smearing_obs( lat , type , count-1 , GLU_TRUE ) ;
  }
  #endif

  // free our temporary lattice
  free( lat2 ) ;
  return ;
#endif
}

// General ALL-dimensional smearing routines ...
void 
smear4D( struct site *__restrict lat ,
	 const int smiters , 
	 const smearing_types type )
{
  if( unlikely( smiters == 0 ) ) { return ; }
  // check our smearing type ...
  if( ( type != SM_APE ) && ( type != SM_STOUT ) && ( type != SM_LOG ) ) {
    printf( "[SMEAR] Unrecognised type [ %d ] ... Leaving \n" , type ) ; 
    return ; 
  }

#if ND != 4 
  // oh god, so I wimped out and allow you to set the alpha in the input file
  const GLU_real alpha1 = ND > 2 ? Latt.sm_alpha[0] / ( ( ND - 1 ) * ( ND - 2 ) ) : Latt.sm_alpha[0] ;
  const GLU_real one_min_a1 = ( 1.0 - Latt.sm_alpha[0] ) ;
#endif	
				     
  #ifdef TOP_VALUE
  double qtop_new , qtop_old = 0. ;
  #endif

  struct spt_site *lat2 = malloc( LCU * sizeof( struct spt_site ) ) ; 
  #ifdef IMPROVED_SMEARING
  struct spt_site *lat3 = malloc( 2*LCU * sizeof( struct spt_site ) ) ; 
  struct spt_site *lat4 = malloc( 2*LCU * sizeof( struct spt_site ) ) ; 
  #else
  struct spt_site *lat3 = malloc( LCU * sizeof( struct spt_site ) ) ; 
  struct spt_site *lat4 = malloc( LCU * sizeof( struct spt_site ) ) ; 
  #endif
 
  int count = 0 ; 

  for( count = 1 ; count <= smiters ; count++ ) {

    #ifdef SYMANZIK_ONE_LOOP
    improve = av_plaquette( lat ) ;
    #endif

    int i ;
    //this bit initialises the calculation by working out the staples for the last time slice first
    #ifdef IMPROVED_SMEARING
    const int back = lat[ lat[0].back[ ND-1 ] ].back[ ND-1 ] ;
    #else
    const int back = lat[0].back[ ND - 1 ] ;
    #endif

    #pragma omp parallel for private(i) SCHED
    #ifdef IMPROVED_SMEARING
    PFOR( i = 0 ; i < 2*LCU ; i++ ) {
    #else
    PFOR( i = 0 ; i < LCU ; i++ ) {
    #endif
      const int bck = back + i ;
      int mu ;

      for( mu = 0 ; mu < ND ; mu++ )   { 
	GLU_complex stap[ NCNC ] = { } ;
        #ifdef IMPROVED_SMEARING
	all_staples_improve( stap , lat , bck , mu , ND , type ) ;
        #else
	all_staples( stap , lat , bck , mu , ND , type ) ;
        #endif

	switch( type )
	  {
	  case SM_APE:
	    project_APE( lat4[ i ].O[ mu ] , stap , 
			 lat[ bck ].O[ mu ] , alpha1 ,
			 one_min_a1 ) ; 
	    break ; 
	  case SM_STOUT:
	    project_STOUT_short( lat4[ i ].O[ mu ] , stap ,
				 lat[ bck ].O[ mu ] , alpha1 ) ; 
	    break ; 
	  case SM_LOG:
	    project_LOG_short( lat4[ i ].O[ mu ] , stap ,
			       lat[ bck ].O[ mu ] , alpha1 ) ; 
	    break ; 
	  default : break ;
	  }
      }
    }
      
    int t ;
    #ifdef IMPROVED_SMEARING
    for( t = 0 ; t < Latt.dims[ ND - 1 ] - 2 ; t++ ) {
    #else
    for( t = 0 ; t < Latt.dims[ ND - 1 ] - 1 ; t++ ) {
    #endif
      const int slice = LCU * t ; 

      #pragma omp parallel for private(i) SCHED
      PFOR( i = 0 ; i < LCU ; i++ ) {
	const int it = slice + i ;
	int mu ;

	for( mu = 0 ; mu < ND ; mu++ ) {
	  GLU_complex stap[ NCNC ] = { } ;
          #ifdef IMPROVED_SMEARING
	  all_staples_improve( stap , lat , it , mu , ND , type ) ;
          #else
	  all_staples( stap , lat , it , mu , ND , type ) ;
          #endif
		
	  switch( type )
	    {
	    case SM_APE:
	      project_APE( lat2[ i ].O[ mu ] , stap ,
			   lat[ it ].O[ mu ] , alpha1 ,
			   one_min_a1 ) ; 
	      break ; 
	    case SM_STOUT:
	      project_STOUT_short( lat2[ i ].O[ mu ] , stap ,
				   lat[ it ].O[ mu ] , alpha1 ) ; 
	      break ; 
	    case SM_LOG:
	      project_LOG_short( lat2[ i ].O[ mu ] , stap , 
				 lat[ it ].O[ mu ] , alpha1 ) ; 
	      break ; 
	    default : break ;
	    }
	}
      }

      #ifdef IMPROVED_SMEARING
      const int bck = lat[ lat[ slice ].back[ ND-1 ] ].back[ ND-1 ] ;
      #else
      const int bck = lat[ slice ].back[ ND -1 ] ;
      #endif

      #pragma omp parallel for private(i)
      PFOR( i = 0 ; i < LCU ; i++ ) {
        #ifdef IMPROVED_SMEARING
	if( likely( t > 1 ) ) { 
	  register const int back = bck + i ;
	  memcpy( &lat[back] , &lat3[i] , sizeof( struct spt_site ) ) ;
	}
	#else
	if( likely( t != 0 ) ) { 
	  register const int back = bck + i ;
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
    const int slice = LCU * t ;
    const int behind = lat[ slice ].back[ ND - 1 ] ;
    #ifdef IMPROVED_SMEARING
    const int behind2 = lat[ behind ].back[ ND-1 ] ;
    #endif

    #pragma omp parallel for private(i)
    PFOR( i = 0 ; i < LCU ; i++ ) {
      register const int back = behind + i ;
      register const int it = slice + i ; 
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

  // need to clean up the lattice fields
#ifdef FAST_SMEAR
  if( type == SM_STOUT ) {
    latt_reunitU( lat ) ;
    printf( "\n[SMEAR] A final reunitarisation step to clean things up\n" ) ;
    print_smearing_obs( lat , type , count , GLU_TRUE ) ;
  }
#endif

  // free stuff !
  free( lat2 ) ; 
  free( lat3 ) ; 
  free( lat4 ) ;

  return ;
}
