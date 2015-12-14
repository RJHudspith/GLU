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

#if ND == 4
static void 
get_lv1( struct smallest_lv1 *__restrict lev1 ,
	 const struct site *__restrict lat ,
	 const int type ,
	 void (*project) ( GLU_complex smeared_link[ NCNC ] ,
			   const GLU_complex staple[ NCNC ] ,
			   const GLU_complex link[ NCNC ] ,
			   const double smear_alpha ,
			   const double al ) )
{
  size_t i ; 
  //do the whole lattice
#pragma omp parallel for private(i) SCHED
  PFOR( i = 0  ;  i < LVOLUME  ;  i++ ) {
    GLU_complex a[ NCNC ] ;
    GLU_complex b[ NCNC ] ;

    size_t j = 0 , mu , nu ;
    //calculate the level1 staples
    for( mu = 0  ;  mu < ND  ;  mu++  ) {
      for( nu = 0  ;  nu < ND  ;  nu++ ) {
	if( likely( nu != mu ) ) {
	  GLU_complex stap[ NCNC ] ;
	  zero_mat( stap ) ;
	  // staple counter
	  size_t temp = lat[i].neighbor[nu] ; 
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
	  project( b , stap , lat[i].O[mu] , alpha3 , one_min_a3 ) ; 
	  // maybe have a shortened SU(2) version????
	  shorten( lev1[i].O[j] , b ) ;
	  j++ ; 
	}
      }
    }
  }
  return ;
}


//calculate the 4D level2 staples 
static void 
get_lv22( struct smallest_lv1 *__restrict lev2 ,
	  const struct smallest_lv1 *__restrict lev1 ,
	  const struct site *__restrict lat , 
	  const size_t t ,
	  const int type ,
	  void (*project) ( GLU_complex smeared_link[ NCNC ] ,
			    const GLU_complex staple[ NCNC ] ,
			    const GLU_complex link[ NCNC ] ,
			    const double smear_alpha ,
			    const double al ) )
{
  size_t i ;
  //could do a slice above and below?
#pragma omp parallel for private(i) 
  PFOR( i = 0 ; i < LCU  ; i++ ) {
    const size_t it = LCU * t + i ; 
    int ii = 0 , mu , nu ; 
    GLU_complex enlarge[ NCNC ] ;
    GLU_complex enlarge2[ NCNC ] ; 
    GLU_complex a[ NCNC ] ;
    GLU_complex b[ NCNC ] ;
    for( mu = 0  ;  mu < ND  ;  mu++  ) {
      //calculate the staples using the dressed links
      for( nu = 0 ;  nu < ND ;  nu++ ) {
	if( unlikely( nu == mu ) ) { continue ; } 
	GLU_complex stap[ NCNC ] ;
	zero_mat( stap ) ;
	size_t rho ;
	for( rho = 0 ;  rho < ND ;  rho++ ) {
	  if( rho == mu || rho == nu ) { continue ; } 
	  // 4-th orthogonal direction: sigma 
	  size_t jj , sigma = 0 ;
	  for( jj = 0 ; jj < ND ; jj++ ) {
	    if( jj != mu && jj != nu && jj != rho ) {
	      sigma = jj ; 
	    }
	  }		  
	  jj = ( ND - 1 )*mu + sigma ; 
	  if( sigma > mu  ) jj-- ; 
	  size_t kk = ( ND - 1 )*rho + sigma ; 
	  if( sigma > rho  ) kk-- ; 	
	  //kk , jj , kk are the correct steps for the staples rho-mu plane
	  size_t temp = lat[it].neighbor[rho] ; 
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
	GLU_complex b[ NCNC ] ;
	zero_mat( b ) ;
	project( b , stap , lat[it].O[mu] , alpha2 , one_min_a2 ) ; 
	shorten( lev2[i].O[ii] , b ) ; 
      }
    } 
  }
  return ;
}

//calculate the 4D staples 
inline static void 
gen_staples_4D2( GLU_complex stap[ NCNC ] ,
		 const struct site *__restrict lat ,
		 const struct smallest_lv1 *__restrict lev2 ,
		 const struct smallest_lv1 *__restrict lev2_up ,
		 const struct smallest_lv1 *__restrict lev2_down ,
		 const size_t i , 
		 const size_t mu , 
		 const size_t t , 
		 const int type )
{
  GLU_complex a[ NCNC ] , b[ NCNC ] ;
  GLU_complex enlarge[ NCNC ] , enlarge2[ NCNC ] ; 
  const size_t it = LCU * t + i ; 

  //calculate the staples using the dressed links
  size_t nu ;
  for( nu = 0 ;  nu < ND ;  nu++ ) {
    if( likely( mu != nu ) ) {
      size_t jj = ( ND - 1 ) * mu + nu ; 
      if( nu > mu  ) { jj-- ; } 

      size_t kk = ( ND - 1 ) * nu + mu ; 
      if( mu > nu  ) { kk-- ; } 
	
      size_t temp ;
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
	      const size_t smiters , 
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

  // callback for the projections
  void (*project) ( GLU_complex smeared_link[ NCNC ] , 
		    const GLU_complex staple[ NCNC ] , 
		    const GLU_complex link[ NCNC ] , 
		    const double smear_alpha , 	     
		    const double al ) ;

  project = project_APE ;
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
    return ;
  }

  // allocate temporaries ...
  struct spt_site *lat2 = NULL , *lat3 = NULL , *lat4 = NULL ;
  struct smallest_lv1 *lev1 = NULL , *lev2 = NULL , 
    *lev2_up = NULL , *lev2_down = NULL ;

  if( GLU_malloc( (void**)&lat2 , 16 , LCU * sizeof( struct spt_site ) ) != 0 ||
      GLU_malloc( (void**)&lat3 , 16 , LCU * sizeof( struct spt_site ) ) != 0 || 
      GLU_malloc( (void**)&lat4 , 16 , LCU * sizeof( struct spt_site ) ) != 0 || 
      GLU_malloc( (void**)&lev1 , 16 , LVOLUME * sizeof( struct smallest_lv1 ) ) != 0 || 
      GLU_malloc( (void**)&lev2 , 16 , LCU * sizeof( struct smallest_lv1 ) ) != 0 || 
      GLU_malloc( (void**)&lev2_up , 16 , LCU * sizeof( struct smallest_lv1 ) ) != 0 || 
      GLU_malloc( (void**)&lev2_down , 16 , LCU * sizeof( struct smallest_lv1 ) ) != 0 ) {
    printf( "[SMEARING] field allocation error\n" ) ;
    return ;
  }

#ifdef TOP_VALUE
  double qtop_new , qtop_old = 0. ;
#endif
  
  size_t count ;
  for( count = 1 ; count <= smiters ; count++ ) {
    get_lv1( lev1 , lat , type , project ) ; 
      
    //init level2's
    get_lv22( lev2_down , lev1 , lat , Latt.dims[ND-1] - 2 , type , project ) ; 
    get_lv22( lev2 , lev1 , lat , Latt.dims[ND-1] - 1 , type , project ) ; 
    get_lv22( lev2_up , lev1 , lat , 0 , type , project ) ; 
      
    const size_t back = lat[ 0 ].back[ ND - 1 ] ; 
    size_t i ;
    #pragma omp parallel for private(i) 
    PFOR( i = 0 ; i < LCU ; i++ ) {
      const size_t bck = back + i ;
      size_t mu ;
      for( mu = 0 ; mu < ND ; mu++ ) {
	GLU_complex stap[ NCNC ] ;
	zero_mat( stap ) ;
	gen_staples_4D2( stap , lat , 
			 lev2 , lev2_up , lev2_down , 
			 i , mu , Latt.dims[ ND - 1 ] - 1 , type ) ; 
	project( lat4[ i ].O[ mu ] , stap , lat[ bck ].O[ mu ] , alpha1 , 
		 one_min_a1 ) ; 
      }
    }
    ///exchange the level2's
    #pragma omp parallel for private(i) 
    PFOR( i = 0 ; i < LCU ; i++ ) {
      memcpy( &lev2_down[i] , &lev2[i] , sizeof( struct smallest_lv1 ) ) ;
      memcpy( &lev2[i] , &lev2_up[i] , sizeof( struct smallest_lv1 ) ) ; 
    }

    // Loop time slices
    size_t t ;
    for( t = 0 ; t < Latt.dims[ ND - 1 ] - 1 ; t++ ) {
      size_t i ;
      //calculate the next level2 up
      get_lv22( lev2_up , lev1 , lat , type , t + 1 , project ) ; 
	  
      const size_t slice = LCU * t ;
      #pragma omp parallel for private(i) 
      PFOR( i = 0 ; i < LCU ; i++ )  {
	const size_t it = slice + i ; 
	size_t mu ;

	for( mu = 0 ; mu < ND ; mu++ ) {
	  GLU_complex stap[ NCNC ] ;
	  zero_mat( stap ) ;
	  gen_staples_4D2( stap , lat , 
			   lev2 , lev2_up , lev2_down , 
			   i , mu , t , type ) ; 
	  project( lat2[ i ].O[ mu ] , stap , lat[ it ].O[ mu ] , 
		   alpha1 , one_min_a1 ) ; 
	}
      }
	  
      const size_t bck = lat[ slice ].back[ ND - 1 ] ;
      #pragma omp parallel for private(i) 
      PFOR( i = 0 ; i < LCU ; i++ ) {
	//put temp into the previous time-slice
	if( likely( t != 0 ) ) { 
	  register const size_t back = bck + i ;
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
    const size_t slice = LCU * t ;
    const size_t behind = lat[ slice ].back[ ND - 1 ] ;
    #pragma omp parallel for private(i) 
    PFOR( i = 0 ; i < LCU ; i++ ) {
      register const size_t back = behind + i ; 
      memcpy( &lat[back] , &lat3[i] , sizeof( struct spt_site ) ) ; 
      register const size_t it = slice + i ; 
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
