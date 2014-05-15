/*
    Copyright 2013 Renwick James Hudspith

    This file (4D_fast.c) is part of GLU.

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
   @file 4D_fast.c 
   @brief 4D HYP,HEX and HYL

   HEX smearing alphas are related to HYP's
   via HYP::HYP (a1,a2,a3) = (3a1,2a2,a3)
   
   @warning only implemented for ND=4
 */

#include "Mainfile.h"
#include "plaqs_links.h"
#include "projectors.h"

// If we are using the dangerous smearing routines ...
#ifdef FAST_SMEAR
  #include "random_config.h"
#endif

#if ND == 4
// precompute the level 1 dressed links ....

static void 
get_lv1( lev1 , lat , type )
     const struct site *__restrict lat ; 
     struct lv1 *__restrict lev1 ; 
     const int type ; 
{
  int i ; 
  //do the whole lattice
#pragma omp parallel for private(i) SCHED
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    GLU_complex a[ NCNC ] ;
    GLU_complex b[ NCNC ] ;
    int j = -1 ;
    int mu ;
    //calculate the level1 staples
    for( mu = 0  ;  mu < ND  ;  mu++  ) {
      int nu ;
      for( nu = 0  ;  nu < ND  ;  nu++ ) {
	if( likely( nu != mu ) ) {
	  // j is our staple counter
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
		
	  //bottom staple
	  GLU_complex c[ NCNC ] ;
	  temp = lat[i].back[nu] ; 
	  multabdag_suNC( a , lat[temp].O[nu] , lat[temp].O[mu] ) ; 
	  temp = lat[temp].neighbor[mu] ; 
	  multab_suNC( c , a , lat[temp].O[nu] ) ; 

	  if( type == SM_LOG ) {
	    multab_dag_suNC( a , c , lat[i].O[mu] ) ; 
            #ifdef FAST_SMEAR
	    exact_log_fast( c , a ) ; 
            #else
	    exact_log_slow( c , a ) ; 
            #endif
	  }
	  a_plus_b( b , c ) ; 
		
	  switch( type )
	    {
	    case SM_APE :
	      project_APE( lev1[i].O[j] , b , lat[i].O[mu] , 
			   alpha3 , one_min_a3 ) ; 
	      break ; 
	    case SM_STOUT :
	      project_STOUT_short( lev1[i].O[j] , b , lat[i].O[mu] , alpha3 ) ; 
	      break ; 
	    case SM_LOG :
	      project_LOG_short( lev1[i].O[j] , b , lat[i].O[mu] , alpha3 ) ;
	      break ; 
	    }

	}
      }
    }
  }
  return ;
}

//calculate the 4D level2 staples 
static void 
get_lv2( lev2 , lev1 , lat , type )
     const struct site *__restrict lat ; 
     const struct lv1 *__restrict lev1 ;
     struct lv1 *__restrict lev2 ; 
     const int type ; 
{
  int i ;
#pragma omp parallel for private(i) SCHED
  PFOR( i = 0  ;  i < LVOLUME  ;  i++ ) {
    int rho  = -1 , sigma = 0 ;
    int ii = -1 ;
    int mu , nu ;
    GLU_complex b[ NCNC ] ;

    for( mu = 0  ;  mu < ND  ;  mu++  ) {
      //calculate the staples using the dressed links
      for( nu = 0 ;  nu < ND ;  ++nu )  {
	if( unlikely( nu == mu ) ) { continue ; } 
	    	    
	ii++ ; 

	// int rho  = -1 , sigma = 0 ;
	GLU_complex stap[ NCNC ] = { } ;

	for( rho = 0 ;  rho < ND ;  ++rho ) {
	  int jj = -1 ;

	  if( rho == mu || rho == nu ) { continue ; } 

	  for( jj = 0 ;  jj < ND ;  ++jj ) {
	    if( jj != mu && jj != nu && jj != rho ) {
	      sigma = jj ; 
	    }
	  }

	  GLU_complex a[ NCNC ] ;
		
	  jj = ( ND - 1 )*mu + sigma ; 
	  if( sigma > mu  ) { jj-- ; } 

	  int kk = ( ND - 1 )*rho + sigma ; 
	  if( sigma > rho  ) { kk-- ; } 
		
	  //kk , jj , kk are the correct steps for the staples rho-mu plane	
	  int temp = lat[i].neighbor[rho] ; 
	  multab_suNC( a , lev1[i].O[kk] , lev1[temp].O[jj] ) ; 
	  temp = lat[i].neighbor[mu] ; 
	  multab_dag_suNC( b , a , lev1[temp].O[kk] ) ; 
		
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
	  temp = lat[i].back[rho] ; 
	  multabdag_suNC( a , lev1[temp].O[kk] , lev1[temp].O[jj]  ) ; 
	  temp = lat[temp].neighbor[mu] ;  
	  multab_suNC( b , a , lev1[temp].O[kk] ) ; 
		
	  if( type == SM_LOG ) {
	    multab_dag_suNC( a , b , lat[i].O[mu] ) ; 
            #ifdef FAST_SMEAR
	    exact_log_fast( b , a ) ; 
            #else
	    exact_log_slow( b , a ) ; 
            #endif
	  }
	  a_plus_b( stap , b ) ; 
	}
	// here are the projections; SM_APE, SM_STOUT and SM_LOG
	switch( type )
	  {
	  case SM_APE :
	    project_APE( lev2[i].O[ii] , stap , lat[i].O[mu] , 
			 alpha2 , one_min_a2 ) ; 
	    break ; 
	  case SM_STOUT :
	    project_STOUT( lev2[i].O[ii] , stap , lat[i].O[mu] , alpha2 ) ; 
	    break ; 
	  case SM_LOG :
	    project_LOG( lev2[i].O[ii] , stap , lat[i].O[mu] , alpha2 ) ; 
	    break ; 
	  }
	// end of projections
      }
    } 
  }
  return ;
}

// complete the staples ....
static void 
gen_staples_4D( stap , lat , lev2 , i , mu , type )
     struct site *__restrict lat ; 
     const struct lv1 *__restrict lev2 ; 
     GLU_complex *__restrict stap ; 
     const int i , mu , type ; 
{
  GLU_complex a[ NCNC ] ;
  GLU_complex b[ NCNC ] ; 
  int nu ;
  //calculate the staples using the dressed links
  for( nu = 0 ;  nu < ND ;  ++nu ) {
    if( likely( mu != nu ) ) {
      int jj = ( ND - 1 ) * mu + nu ; 
      if( nu > mu  ) { jj-- ; } 

      int kk = ( ND - 1 )*nu + mu ; 
      if( mu > nu  ) { kk-- ; } 

      //kk , jj , kk are the correct steps for the staples nu-mu plane
      int temp = lat[i].neighbor[nu] ; 
      multab_suNC( a , lev2[i].O[kk] , lev2[temp].O[jj]  ) ; 
      temp = lat[i].neighbor[mu] ; 
      multab_dag_suNC( b , a , lev2[temp].O[kk] ) ; 

      if( type == SM_LOG ) {
	multab_dag_suNC( a , b , lat[i].O[mu] ) ; 
        #ifdef FAST_SMEAR
	exact_log_fast( b , a ) ; 
        #else
	exact_log_slow( b , a ) ; 
        #endif
      }
      a_plus_b( stap , b ) ; 

      temp = lat[i].back[nu] ; 
      multabdag_suNC( a , lev2[temp].O[kk] , lev2[temp].O[jj] ) ; 
      temp = lat[temp].neighbor[mu] ; 
      multab_suNC( b , a , lev2[temp].O[kk] ) ; 

      if( type == SM_LOG ) {
	multab_dag_suNC( a , b , lat[i].O[mu] ) ; 
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

////////////////////////////////////////////////////

// this code performs the smearing ...
void 
HYPSLsmear4D_expensive( struct site *__restrict lat , 
			const int smiters , 
			const int type )
{
#if ND != 4
  return ;
#else
  if( smiters == 0 ) { return ; }

  // check the type
  if( ( type != SM_APE ) && ( type != SM_STOUT ) && ( type != SM_LOG ) ) {
    printf( "[SMEAR] Unrecognised type [ %d ] ... Leaving \n" , type ) ; 
    return ; 
  }

#ifdef TOP_VALUE
  double qtop_new , qtop_old = 0. ;
#endif

  // allocate temporary lattices ...
  struct spt_site *lat2 = malloc( LCU * sizeof( struct spt_site ) ) ; 
  struct spt_site *lat3 = malloc( LCU * sizeof( struct spt_site ) ) ; 
  struct spt_site *lat4 = malloc( LCU * sizeof( struct spt_site ) ) ; 

  struct lv1 *lev1 = malloc( LVOLUME * sizeof( struct lv1 ) ) ; 
  struct lv1 *lev2 = malloc( LVOLUME * sizeof( struct lv1 ) ) ; 

  int count = 0 ; 
  for( count = 1 ; count <= smiters ; count++ ) {
    int i ;

    get_lv1( lev1 , lat , type ) ; 
    get_lv2( lev2 , lev1 , lat , type ) ;
      
    const int bck = lat[0].back[ ND - 1 ] ; 
    #pragma omp parallel for private(i) SCHED
    PFOR( i = 0 ; i < LCU ; i++ ) {
      const int back = bck + i ;
      int mu ;

      for( mu = 0 ; mu < ND ; mu++ ) {
	GLU_complex stap[ NCNC ] = { } ;

	gen_staples_4D( stap , lat , lev2 , back , mu , type ) ; 
	      
	switch( type )
	  {
	  case SM_APE :
	    project_APE( lat4[i].O[mu] , stap , 
			 lat[back].O[mu] , alpha1 , 
			 one_min_a1 ) ; 
	    break ; 
	  case SM_STOUT :
	    project_STOUT_short( lat4[i].O[mu] , stap , 
				 lat[back].O[mu] , alpha1 ) ; 
	    break ; 
	  case SM_LOG :
	    project_LOG_short( lat4[i].O[mu] , stap , 
			       lat[back].O[mu] , alpha1 ) ; 
	    break ; 
	  }

      }
    }
 
    //loop time slices
    ///////////////////////////////////////
    int t ;
    for( t = 0 ; t < Latt.dims[ ND - 1 ] - 1 ; t++ ) {
      const int slice = LCU * t ; 
      const int bck = lat[ slice ].back[ ND - 1 ] ;
      #pragma omp parallel for private(i) SCHED
      PFOR( i = 0 ; i < LCU ; i++ )  {
	const int it = slice + i ; 
	int mu ;
	      
	for( mu = 0 ; mu < ND ; mu++ ) {
	  GLU_complex stap[ NCNC ] = { } ;
		  
	  gen_staples_4D( stap , lat , lev2 , it , mu , type ) ; 
		  
	  switch( type )
	    {
	    case SM_APE :
	      project_APE( lat2[i].O[mu] , stap , 
			   lat[it].O[mu] , alpha1 , 
			   one_min_a1 ) ; 
	      break ; 
	    case SM_STOUT :
	      project_STOUT_short( lat2[i].O[mu] , stap , lat[it].O[mu] , 
				   alpha1 ) ; 
	      break ; 
	    case SM_LOG :
	      project_LOG_short( lat2[i].O[mu] , stap , lat[it].O[mu] , 
				 alpha1 ) ; 
	      break ; 
	    }
	}
	// this is only a legal maneuver for this method
	//put temp into the previous time-slice
	if( likely( t != 0 ) ) { 
	  register const int back = bck + i ;
	  memcpy( &lat[back] , &lat3[i] , sizeof( struct spt_site ) ) ;
	}
	//make temporary lat3 lat2 again and repeat
	memcpy( &lat3[i] , &lat2[i] , sizeof( struct spt_site ) ) ;
      }
    }
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
  
  // free that memory //
  free( lev1 ) ; 
  free( lev2 ) ; 
  free( lat2 ) ; 
  free( lat3 ) ; 
  free( lat4 ) ; 

  return ;
#endif
}

