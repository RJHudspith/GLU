/*
    Copyright 2013 Renwick James Hudspith

    This file (HYP.c) is part of GLU.

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
   @file HYP.c
   @brief 3D HYP,HEX and HYL

   HEX smearing alphas are related to HYP's
   via HYP::HYP (a1,a2,a3) = (3a1,2a2,a3)
   <br>
   is used for spatial smearing the 4D links

   @warning 4D only   
*/

#include "Mainfile.h"
#include "plaqs_links.h"
#include "projectors.h"

// If we are using the dangerous smearing routines ...
#ifdef FAST_SMEAR
  #include "random_config.h"
#endif

#if ND == 4
// 3D level-1 staples for timeslice t
static void 
get_spatial_lv1( lev1 , lat , t , type )
     struct spatial_lv1 *__restrict lev1 ; 
     const struct site *__restrict lat ; 
     const int t , type ; 
{
  int it ; 
  const int slice = LCU * t ; 
  //do a slice
#pragma omp parallel for private(it) SCHED
  PFOR( it = 0  ;  it < LCU  ;  it++ )  {

    GLU_complex a[ NCNC ] ;
    GLU_complex c[ NCNC ] ;
    const int i = slice + it ;
    int j = -1 ; 
    int mu ;
    int nu ;

    //calculate the level1 staples
    for( mu = 0 ; mu < ND - 1 ; mu++ ) {
      for( nu = 0 ; nu < ND - 1 ; nu++ ) {
	if( likely( nu != mu ) ) {
	  // b is our staple
	  GLU_complex b[ NCNC ] = { } ;
	  // j is our staple counter
	  j ++ ; 
	    
	  int temp = lat[i].neighbor[nu] ; 
	  multab_suNC( a , lat[i].O[nu] , lat[temp].O[mu] ) ; 
	  temp = lat[i].neighbor[mu] ; 
	  multab_dag_suNC( b , a , lat[temp].O[nu] ) ; 
	  if( type == SM_LOG ) {
	    multab_dag_suNC( a , b , lat[i].O[mu] ) ; 
            #ifdef SLOW_SMEAR
	    exact_log_slow( b , a ) ; 
            #else
	    exact_log_fast( b , a ) ; 
            #endif
	  }

	  // put the bottom staple in "c"
	  //bottom staple
	  temp = lat[i].back[nu] ; 
	  multabdag_suNC( a , lat[temp].O[nu] , lat[temp].O[mu] ) ; 
	  temp = lat[temp].neighbor[mu] ; 
	  multab_suNC( c , a , lat[temp].O[nu] ) ; 
	  if( type == SM_LOG ) {
	    multab_dag_suNC( a , c , lat[i].O[mu] ) ; 
            #ifdef SLOW_SMEAR
	    exact_log_slow( c , a ) ; 
            #else
	    exact_log_fast( c , a ) ; 
            #endif	      
	  }
	  a_plus_b( b , c ) ; 
	  
	  switch( type )
	    {
	    case SM_APE :
	      project_APE( lev1[it].O[j] , b , 
			   lat[i].O[mu] , alpha2 , 
			   one_min_a2 ) ; 
	      break ; 
	    case SM_STOUT :
	      project_STOUT_short( lev1[it].O[j] , b , lat[i].O[mu] , alpha2 ) ;  
	      break ; 
	    case SM_LOG :
	      project_LOG_short( lev1[it].O[j] , b , lat[i].O[mu] , alpha2 ) ; 
	      break ; 
	    }
	}
      }
    }
  }
  return ;
}

static void 
staples3D( stap , lat , lev1 , i , mu , t , type )
     struct site *__restrict lat ; 
     const struct spatial_lv1 *__restrict lev1 ; 
     GLU_complex stap[ NCNC ] ; 
     const int i , mu , type ; 
{ 
  GLU_complex a[ NCNC ] , b[ NCNC ] ;
  int nu ;
  const int it = LCU * t + i ; 

  //calculate the staples using the dressed links
  for( nu = 0 ;  nu < ND - 1 ;  ++nu ) {
    if ( likely( nu != mu ) ) {		   
      int jj = -1 , rho = -1 ;
      /* 3rd orthogonal direction: rho */
      for( jj = 0 ;  jj < ND - 1 ;  ++jj ) {
	if( jj != mu && jj != nu ) {
	  rho = jj ; 
	}
      }
	
      jj = ( ND - 2 ) * mu + rho ; 
      if( rho > mu  ) { jj-- ; } 

      int kk = ( ND - 2 ) * nu + rho ; 
      if( rho > nu  ) { kk-- ; } 
	
      //kk , jj , kk are the correct steps for the staples
      int temp = lat[i].neighbor[nu] ; 
      multab_suNC( a , lev1[i].O[kk] , lev1[temp].O[jj] ) ; 
      temp = lat[i].neighbor[mu] ; 
      multab_dag_suNC( b , a , lev1[temp].O[kk] ) ; 

      if( type == SM_LOG ) {
	multab_dag_suNC( a , b , lat[it].O[mu] ) ; 
        #ifdef SLOW_SMEAR
	exact_log_slow( b , a ) ; 
        #else 
	exact_log_fast( b , a ) ; 
        #endif
      }
      a_plus_b( stap , b ) ; 
	
      //bottom staple
      temp = lat[i].back[nu] ; 
      multabdag_suNC( a , lev1[temp].O[kk] , lev1[temp].O[jj] ) ;
      temp = lat[temp].neighbor[mu] ; 
      multab_suNC( b , a , lev1[temp].O[kk] ) ; 
      
      if( type == SM_LOG ) {
	multab_dag_suNC( a , b , lat[it].O[mu] ) ; 
        #ifdef SLOW_SMEAR
	exact_log_slow( b , a ) ; 
        #else
	exact_log_fast( b , a ) ; 
        #endif
      }
      a_plus_b( stap , b ) ; 
    }
  }
  return ;
}
#endif

// spatial only smearing
void 
HYPSLsmear3D( struct site *__restrict lat , 
	      const int smiters , 
	      const int type ) 
{
  if( unlikely( smiters == 0 ) ) { return ; }
#if ND != 4
  printf( "[SMEAR] Sorry, HYP/HEX/HYL not supported for ND == %d \n" , ND ) ;
  return ;
#else

  // check the type
  if( ( type != SM_APE ) && ( type != SM_STOUT ) && ( type != SM_LOG ) ) {
    printf( "[SMEAR] Unrecognised type [ %d ] ... Leaving \n" , type ) ; 
    return ; 
  }

  struct spatial_lv1 *lev1 = NULL ;
  struct sp_site *lat2 = NULL ;
  if( GLU_malloc( (void**)&lev1 , 16 , LCU * sizeof( struct spatial_lv1 ) ) != 0 || 
      GLU_malloc( (void**)&lat2 , 16 , LCU * sizeof( struct sp_site ) ) != 0 ) {
    printf( "[SMAERING] field allocation failure\n" ) ;
    return ;
  }
 
  int count = 0 ; 

  for( count = 1 ; count <= smiters ; count++ )   {
    int t ;
    //loop time slices
    for( t = 0 ; t < Latt.dims[ ND - 1 ] ; t++ ) {
      // get the level 1 links for this slice
      get_spatial_lv1(  lev1  ,  lat  ,  t  ,  type ) ; 
  
      const int slice = LCU * t ; 
      int i ;
      #pragma omp parallel for private(i) SCHED
      PFOR( i = 0 ; i < LCU ; i++ )  {
	const int it = slice + i ;
	int mu ;
	for( mu = 0 ; mu < ND - 1 ; mu++ ) {
	  GLU_complex stap[ NCNC ] = { } ;
	  staples3D( stap , lat , lev1 , i , mu , t , type ) ; 

	  switch( type )
	    {
	    case SM_APE :
	      project_APE( lat2[ i ].O[ mu ] , stap , 
			   lat[ it ].O[ mu ] , alpha1 , 
			   one_min_a1 ) ; 
	      break ; 
	    case SM_STOUT :
	      project_STOUT_short( lat2[ i ].O[ mu ] , stap , 
				   lat[ it ].O[ mu ] , alpha1 ) ; 
	      break ; 
	    case SM_LOG :
	      project_LOG_short( lat2[ i ].O[ mu ] , stap , 
				 lat[ it ].O[ mu ] , alpha1 ) ; 
	      break ; 
	    }
	}
	// swap these round using memcpy
	memcpy( &lat[it] , &lat2[i] , sizeof( struct sp_site ) ) ;
      }
    }

    // only write these out if we are not doing this in parallel ...
    #ifdef verbose
    print_smearing_obs( lat , type , count , GLU_TRUE ) ;
    #endif
  }

#ifndef verbose
  print_smearing_obs( lat , type , count , GLU_TRUE ) ;
#endif

#ifdef FAST_SMEAR
  if( type == SM_STOUT ) {
    latt_reunitU( lat ) ;
    printf( "\n[SMEAR] A final reunitarisation step to clean things up\n" ) ;
    print_smearing_obs( lat , type , count , GLU_TRUE ) ;
  }
#endif

  free( lat2 ) ; 
  free( lev1 ) ; 

  return ;
#endif
}
