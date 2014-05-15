/*
    Copyright 2013 Renwick James Hudspith

    This file (HYP_4D.c) is part of GLU.

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
   @file HYP_4D.c
   @brief 4D HYP,HEX and HYL temporaries stored using shorten() and rebuild()

   HEX smearing alphas are related to HYP's
   via HYP::HYP (a1,a2,a3) = (3a1,2a2,a3)

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
static void 
get_lv1( lev1 , lat , type )
     const struct site *__restrict lat ; 
     struct smallest_lv1 *__restrict lev1 ; 
     const int type ; 
{
  int i ; 
  //do the whole lattice
#pragma omp parallel for private(i) SCHED
  PFOR( i = 0  ;  i < LVOLUME  ;  i++ ) {
    GLU_complex a[ NCNC ] ;
    GLU_complex b[ NCNC ] ;

    int j = -1 ; 
    int mu ;
    //calculate the level1 staples
    for( mu = 0  ;  mu < ND  ;  mu++  ) {
      int nu ; 
      for( nu = 0  ;  nu < ND  ;  nu++ ) {
	if( likely( nu != mu ) ) {
	  GLU_complex stap[ NCNC ] = { } ;

	  // staple counter
	  j++ ; 

	  int temp = lat[i].neighbor[nu] ; 
	  multab_suNC( a , lat[i].O[nu] , lat[temp].O[mu] ) ; 
	  temp = lat[i].neighbor[mu] ; 
	  multab_dag_suNC( b , a , lat[temp].O[nu] ) ; 

	  if( type == SM_LOG ) {
	    multab_dag_suNC( a , b , lat[i].O[mu] ) ; 
            #ifdef FAST_SMEAR
	    exact_log_fast( b , a ) ; 
            #else
	    exact_log_slow( b , a ) ; 
            #endif
	  }
	  a_plus_b( stap , b ) ; 

	  //bottom staple
	  temp = lat[i].back[nu] ; 
	  multabdag_suNC( a , lat[temp].O[nu] , lat[temp].O[mu] ) ; 
	  temp = lat[temp].neighbor[mu] ; 
	  multab_suNC( b , a , lat[temp].O[nu] ) ; 

	  if( type == SM_LOG ) {
	    multab_dag_suNC( a , b , lat[i].O[mu] ) ; 
            #ifdef FAST_SMEAR
	    exact_log_fast( b , a ) ; 
            #else
	    exact_log_slow( b , a ) ; 
            #endif
	  }
	  a_plus_b( stap , b ) ; 

	  switch( type )
	    {
	    case SM_APE :
	      project_APE( b , stap , 
			   lat[i].O[mu] , alpha3 , 
			   one_min_a3 ) ; 
	      break ; 
	    case SM_STOUT :
	      project_STOUT_short( b , stap , lat[i].O[mu] , alpha3 ) ; 
	      break ; 
	    case SM_LOG :
	      project_LOG_short( b , stap , lat[i].O[mu] , alpha3 ) ; 
	      break ; 
	    }
	  // maybe have a shortened SU(2) version????
	  shorten( lev1[i].O[j] , b ) ; 
	}
      }
    }
  }

  return ;
}


//calculate the 4D level2 staples 
static void 
get_lv22( lev2 , lev1 , lat , type , t )
     struct smallest_lv1 *__restrict lev2 ;
     const struct site *__restrict lat ; 
     const struct smallest_lv1 *__restrict lev1 ;
     const int type , t ; 
{
  int i ;

  //could do a slice above and below?
#pragma omp parallel for private(i) 
  PFOR( i = 0  ;  i < LCU  ;  i++ ) {
    const int it = LCU * t + i ; 
    int mu ;
    int ii = -1 ; 

    GLU_complex enlarge[ NCNC ] ;
    GLU_complex enlarge2[ NCNC ] ; 
    GLU_complex a[ NCNC ] ;
    GLU_complex b[ NCNC ] ;

    for( mu = 0  ;  mu < ND  ;  mu++  ) {
      int nu ; 
      //calculate the staples using the dressed links
      for( nu = 0 ;  nu < ND ;  nu++ ) {
	if( unlikely( nu == mu ) ) { continue ; } 
	    
	ii++ ; 
	GLU_complex stap[ NCNC ] = { } ;
	int rho ;

	for( rho = 0 ;  rho < ND ;  rho++ ) {
	  if( rho == mu || rho == nu ) { continue ; } 

	  // 4-th orthogonal direction: sigma 
	  int jj ;
	  int sigma = 0 ;
	  for( jj = 0 ; jj < ND ; jj++ )
	    if( jj != mu && jj != nu && jj != rho ) 
	      sigma = jj ; 
		  
	  jj = ( ND - 1 )*mu + sigma ; 
	  if( sigma > mu  ) jj-- ; 
	  int kk = ( ND - 1 )*rho + sigma ; 
	  if( sigma > rho  ) kk-- ; 
		
	  //kk , jj , kk are the correct steps for the staples rho-mu plane
	  int temp = lat[it].neighbor[rho] ; 
	  rebuild( enlarge , lev1[it].O[kk] ) ; 
	  rebuild( enlarge2 , lev1[temp].O[jj] ) ; 
	  multab_suNC( a , enlarge , enlarge2 ) ; 
	  temp = lat[it].neighbor[mu] ; 
		
	  rebuild( enlarge , lev1[temp].O[kk] ) ; 
	  multab_dag_suNC( b , a , enlarge ) ; 
		
	  if( type == SM_LOG ) {
	    multab_dag_suNC( a , b , lat[it].O[mu] ) ; 
            #ifdef FAST_SMEAR
	    exact_log_fast( b , a ) ; 
            #else
	    exact_log_slow( b , a ) ; 
            #endif
	  }
	  a_plus_b( stap , b ) ; 
		  
	  //bottom staple
	  temp = lat[it].back[rho] ; 
	  rebuild( enlarge , lev1[temp].O[kk] ) ; 
	  rebuild( enlarge2 , lev1[temp].O[jj] ) ; 
	  multabdag_suNC( a , enlarge , enlarge2 ) ; 
	  temp = lat[temp].neighbor[mu] ; 
	  rebuild( enlarge , lev1[temp].O[kk] ) ; 
	  multab_suNC( b , a , enlarge ) ; 
		
	  if( type == SM_LOG ) {
	    multab_dag_suNC( a , b , lat[it].O[mu] ) ; 
            #ifdef FAST
	    exact_log_fast( b , a ) ; 
            #else
	    exact_log_slow( b , a ) ; 
            #endif
	  }
	  a_plus_b( stap , b ) ; 
	}
	// here are the projections; SM_APE, SM_STOUT and SM_LOG
	GLU_complex b[ NCNC ] = { } ;
	switch( type )
	  {
	  case SM_APE :
	    project_APE( b , stap , 
			 lat[it].O[mu] , alpha2 , 
			 one_min_a2 ) ; 
	    break ; 
	  case SM_STOUT :
	    //project_STOUT( b , stap , lat[it].O[mu] , alpha2 ) ; 
	    project_STOUT_short( b , stap , lat[it].O[mu] , alpha2 ) ; 
	    break ; 
	  case SM_LOG :
	    //project_LOG( b , stap , lat[it].O[mu] , alpha2 ) ; 
	    project_LOG_short( b , stap , lat[it].O[mu] , alpha2 ) ;
	    break ; 
	  }
	shorten( lev2[i].O[ii] , b ) ; 
      }
    } 
  }
  return ;
}

//calculate the 4D staples 
inline static void 
gen_staples_4D2( stap , lat , lev2 , lev2_up , lev2_down , i , mu , t , type )
     struct site *__restrict lat ; 
     const struct smallest_lv1 *__restrict lev2 ;
     const struct smallest_lv1 *__restrict lev2_up ;
     const struct smallest_lv1 *__restrict lev2_down ; 
     GLU_complex stap[ NCNC ] ; 
     const int i , mu , type , t ; 
{
  GLU_complex a[ NCNC ] , b[ NCNC ] ;
  GLU_complex enlarge[ NCNC ] , enlarge2[ NCNC ] ; 
  const int it = LCU * t + i ; 

  //calculate the staples using the dressed links
  int nu ;
  for( nu = 0 ;  nu < ND ;  nu++ ) {
    if( likely( mu != nu ) ) {
      int jj = ( ND - 1 ) * mu + nu ; 
      if( nu > mu  ) { jj-- ; } 

      int kk = ( ND - 1 ) * nu + mu ; 
      if( mu > nu  ) { kk-- ; } 
	
      int temp ;
      //kk , jj , kk are the correct steps for the staples nu-mu plane
      rebuild( enlarge , lev2[i].O[kk] ) ; 
      if( nu == ( ND - 1 ) ) {
	rebuild( enlarge2 , lev2_up[i].O[jj] ) ; 
      } else {
	temp = lat[i].neighbor[nu] ; 
	rebuild( enlarge2 , lev2[temp].O[jj] ) ; 
      }
      multab_suNC( a , enlarge , enlarge2 ) ; 

      if( mu == ( ND - 1 ) ) {
	rebuild( enlarge , lev2_up[i].O[kk] ) ; 
      } else {
	temp = lat[i].neighbor[mu] ; 
	rebuild( enlarge , lev2[temp].O[kk] ) ; 
      }
      multab_dag_suNC( b , a , enlarge ) ; 

      if( type == SM_LOG ) {
	multab_dag_suNC( a , b , lat[it].O[mu] ) ; 
        #ifdef FAST_SMEAR
	exact_log_fast( b , a ) ; 
        #else
	exact_log_slow( b , a ) ; 
        #endif
      }
      a_plus_b( stap , b ) ; 

      if( nu == ( ND - 1 ) ) {
	rebuild( enlarge , lev2_down[i].O[kk] ) ; 
	rebuild( enlarge2 , lev2_down[i].O[jj] ) ; 
	multabdag_suNC( a , enlarge , enlarge2 ) ; 
	temp = lat[i].neighbor[mu] ; 
	rebuild( enlarge , lev2_down[temp].O[kk] ) ; 
      } else {
	temp = lat[i].back[nu] ; 
	rebuild( enlarge , lev2[temp].O[kk] ) ; 
	rebuild( enlarge2 , lev2[temp].O[jj] ) ; 
	multabdag_suNC( a , enlarge , enlarge2 ) ; 
	
	if( mu == ( ND - 1 ) ) {
	  rebuild( enlarge , lev2_up[temp].O[kk] ) ; 
	} else {
	  temp = lat[temp].neighbor[mu] ; 
	  rebuild( enlarge , lev2[temp].O[kk] ) ; 
	}
      }

      multab_suNC( b , a , enlarge ) ; 
      if( type == SM_LOG ) {
	multab_dag_suNC( a , b , lat[it].O[mu] ) ; 
        #ifdef FAST_SMEAR
	exact_log_fast( b , a ) ; 
        #else
	exact_log_slow( b , a ) ; 
        #endif
      }
      a_plus_b( stap , b ) ; 
    }
  }
  return ;
}
#endif 

// this code performs the smearing ...
void 
HYPSLsmear4D( struct site *__restrict lat , 
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

  // allocate temporaries ...
  struct spt_site *lat2 = malloc( LCU * sizeof( struct spt_site ) ) ; 
  struct spt_site *lat3 = malloc( LCU * sizeof( struct spt_site ) ) ; 
  struct spt_site *lat4 = malloc( LCU * sizeof( struct spt_site ) ) ; 

  struct smallest_lv1 *lev1 = malloc( LVOLUME * sizeof( struct smallest_lv1 ) ) ; 

  struct smallest_lv1 *lev2 = malloc( LCU * sizeof( struct smallest_lv1 ) ) ; 
  struct smallest_lv1 *lev2_up = malloc( LCU * sizeof( struct smallest_lv1 ) ) ; 
  struct smallest_lv1 *lev2_down = malloc( LCU * sizeof( struct smallest_lv1 ) ) ; 

#ifdef TOP_VALUE
  double qtop_new , qtop_old = 0. ;
#endif
  
  int count ;
  for( count = 1 ; count <= smiters ; count++ ) {
    get_lv1( lev1 , lat , type ) ; 
      
    //init level2's
    get_lv22( lev2_down , lev1 , lat , type , Latt.dims[ND-1] - 2 ) ; 
    get_lv22( lev2 , lev1 , lat , type , Latt.dims[ND-1] - 1 ) ; 
    get_lv22( lev2_up , lev1 , lat , type , 0 ) ; 
      
    const int back = lat[ 0 ].back[ ND - 1 ] ; 
    int i ;
    #pragma omp parallel for private(i) 
    PFOR( i = 0 ; i < LCU ; i++ ) {
      const int bck = back + i ;
      int mu ;
      for( mu = 0 ; mu < ND ; mu++ ) {
	GLU_complex stap[ NCNC ] = { } ;
	gen_staples_4D2( stap , lat , 
			 lev2 , lev2_up , lev2_down , 
			 i , mu , Latt.dims[ ND - 1 ] - 1 , type ) ; 

	switch( type )
	  {
	  case SM_APE :
	    project_APE( lat4[ i ].O[ mu ] , stap , 
			 lat[ bck ].O[ mu ] , alpha1 , 
			 one_min_a1 ) ; 
	    break ; 
	  case SM_STOUT :
	    project_STOUT_short( lat4[ i ].O[ mu ] , stap , 
				 lat[ bck ].O[ mu ] , alpha1 ) ; 
	    break ; 
	  case SM_LOG :
	    project_LOG_short( lat4[ i ].O[ mu ] , stap , 
			       lat[ bck ].O[ mu ] , alpha1 ) ; 
	    break ; 
	  }
      }
    }
    ///exchange the level2's
    #pragma omp parallel for private(i) 
    PFOR( i = 0 ; i < LCU ; i++ ) {
      memcpy( &lev2_down[i] , &lev2[i] , sizeof( struct smallest_lv1 ) ) ;
      memcpy( &lev2[i] , &lev2_up[i] , sizeof( struct smallest_lv1 ) ) ; 
    }

    // Loop time slices
    int t ;
    for( t = 0 ; t < Latt.dims[ ND - 1 ] - 1 ; t++ ) {
      int i ;
      //calculate the next level2 up
      get_lv22( lev2_up , lev1 , lat , type , t + 1 ) ; 
	  
      const int slice = LCU * t ;
      #pragma omp parallel for private(i) 
      PFOR( i = 0 ; i < LCU ; i++ )  {
	const int it = slice + i ; 
	int mu ;

	for( mu = 0 ; mu < ND ; mu++ ) {
	  GLU_complex stap[ NCNC ] = { } ;
	  gen_staples_4D2( stap , lat , 
			   lev2 , lev2_up , lev2_down , 
			   i , mu , t , type ) ; 
		  
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
      }
	  
      const int bck = lat[ slice ].back[ ND - 1 ] ;
      #pragma omp parallel for private(i) 
      PFOR( i = 0 ; i < LCU ; i++ ) {
	//put temp into the previous time-slice
	if( likely( t != 0 ) ) { 
	  register const int back = bck + i ;
	  memcpy( &lat[back] , &lat3[i] , sizeof( struct spt_site ) ) ;
	}
	//make temporary lat3 lat2 again and repeat
	memcpy( &lat3[i] , &lat2[i] , sizeof( struct spt_site ) ) ;

	// swap the level-2's around as usual 
	memcpy( &lev2_down[i] , &lev2[i] , sizeof( struct smallest_lv1 ) ) ;
	memcpy( &lev2[i] , &lev2_up[i] , sizeof( struct smallest_lv1 ) ) ;
      }
    }
    //put last and last but one time slice in
    ////////////////////////////////////////////
    //put last and last but one time slice in
    ////////////////////////////////////////////
    const int slice = LCU * t ;
    const int behind = lat[ slice ].back[ ND - 1 ] ;
    #pragma omp parallel for private(i) 
    PFOR( i = 0 ; i < LCU ; i++ ) {
      register const int back = behind + i ; 
      memcpy( &lat[back] , &lat3[i] , sizeof( struct spt_site ) ) ; 
      register const int it = slice + i ; 
      memcpy( &lat[it] , &lat4[i] , sizeof( struct spt_site ) ) ; 
    }

    // Are we looking for the topological charge? this is the routine for you
    // as in the wilson flow routine YOU provide the maximum iterations available
    // breaks on convergence
    #ifdef TOP_VALUE
    if( count > TOP_VALUE ) {
      if( gauge_topological_meas( lat , &qtop_new , &qtop_old , count-1 ) 
	  == GLU_SUCCESS ) { break ; }
    }
    #endif

    // print out the info ::
    #ifdef verbose
    print_smearing_obs( lat , type , count , GLU_TRUE ) ;
    #endif
  }

#ifndef verbose
  // print out the final details ...
  print_smearing_obs( lat , type , count , GLU_TRUE ) ;
#endif

  // need to clean up the lattice fields
  #ifdef FAST_SMEAR
  if( type == SM_STOUT ) {
    latt_reunitU( lat ) ;
    printf( "\n[SMEAR] A final reunitarisation step to clean things up\n" ) ;
    print_smearing_obs( lat , type , count , GLU_TRUE ) ;
  }
  #endif

  // free that memory //
  free( lev1 ) ; 
  free( lev2 ) ; 
  free( lev2_up ) ;
  free( lev2_down ) ;
  free( lat2 ) ; 
  free( lat3 ) ; 
  free( lat4 ) ; 

  return ;
#endif
}
