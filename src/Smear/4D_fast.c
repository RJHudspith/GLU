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

#if ND == 4
// precompute the level 1 dressed links ....

static void 
get_lv1( struct lv1 *__restrict lev1 ,
	 const struct site *__restrict lat ,
	 const int type ,
	 void (*project) ( GLU_complex smeared_link[ NCNC ] , 
			   GLU_complex staple[ NCNC ] , 
			   const GLU_complex link[ NCNC ] , 
			   const double smear_alpha , 	     
			   const double al ) )
{
  size_t i ; 
  //do the whole lattice
#pragma omp parallel for private(i) SCHED
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    GLU_complex a[ NCNC ] GLUalign , b[ NCNC ] GLUalign ;
    GLU_complex c[ NCNC ] GLUalign ;
    size_t j = 0 , mu , nu ;
    //calculate the level1 staples
    for( mu = 0  ;  mu < ND  ;  mu++  ) {
      for( nu = 0  ;  nu < ND  ;  nu++ ) {
	if( likely( nu != mu ) ) {	
	  size_t temp = lat[i].neighbor[nu] ; 
	  multab_suNC( a , lat[i].O[nu] , lat[temp].O[mu] ) ; 
	  temp = lat[i].neighbor[mu] ; 
	  multab_dag_suNC( b , a , lat[temp].O[nu] ) ;	
	  if( type == SM_LOG ) {
	    multab_dag_suNC( a , b , lat[i].O[mu] ) ; 
	    exact_log_slow( b , a ) ; 
	  }
	  //bottom staple
	  temp = lat[i].back[nu] ; 
	  multabdag_suNC( a , lat[temp].O[nu] , lat[temp].O[mu] ) ; 
	  temp = lat[temp].neighbor[mu] ; 
	  multab_suNC( c , a , lat[temp].O[nu] ) ; 
	  if( type == SM_LOG ) {
	    multab_dag_suNC( a , c , lat[i].O[mu] ) ; 
	    exact_log_slow( c , a ) ; 
	  }
	  a_plus_b( b , c ) ; 
	  project( lev1[i].O[j] , b , lat[i].O[mu] , 
		   alpha3 , one_min_a3 ) ; 
	  // j is our staple counter
	  j++ ; 
	}
      }
    }
  }
  return ;
}

//calculate the 4D level2 staples 
static void 
get_lv2( struct lv1 *__restrict lev2 ,
	 const struct lv1 *__restrict lev1 ,
	 const struct site *__restrict lat ,
	 const int type ,
	 void (*project) ( GLU_complex smeared_link[ NCNC ] , 
			   GLU_complex staple[ NCNC ] , 
			   const GLU_complex link[ NCNC ] , 
			   const double smear_alpha , 	     
			   const double al ) )
{
  size_t i ;
#pragma omp parallel for private(i) SCHED
  PFOR( i = 0  ;  i < LVOLUME  ;  i++ ) {
    size_t rho = 0 , sigma = 0 ;
    size_t ii = 0 , mu , nu ;
    GLU_complex b[ NCNC ] GLUalign , a[ NCNC ] GLUalign ;
    GLU_complex stap[ NCNC ] GLUalign ;

    for( mu = 0  ;  mu < ND  ;  mu++  ) {
      //calculate the staples using the dressed links
      for( nu = 0 ;  nu < ND ;  ++nu )  {
	if( unlikely( nu == mu ) ) { continue ; } 
	zero_mat( stap ) ;

	for( rho = 0 ; rho < ND ; rho++ ) {

	  size_t jj = 0 ;
	  if( rho == mu || rho == nu ) { continue ; } 

	  for( jj = 0 ;  jj < ND ;  jj++ ) {
	    if( jj != mu && jj != nu && jj != rho ) {
	      sigma = jj ; 
	    }
	  }
		
	  jj = ( ND - 1 )*mu + sigma ; 
	  if( sigma > mu  ) { jj-- ; } 

	  size_t kk = ( ND - 1 )*rho + sigma ; 
	  if( sigma > rho ) { kk-- ; } 
		
	  //kk , jj , kk are the correct steps for the staples rho-mu plane	
	  size_t temp = lat[i].neighbor[rho] ; 
	  multab_suNC( a , lev1[i].O[kk] , lev1[temp].O[jj] ) ; 
	  temp = lat[i].neighbor[mu] ; 
	  multab_dag_suNC( b , a , lev1[temp].O[kk] ) ; 
		
	  if( type == SM_LOG ) {
	    multab_dag_suNC( a , b , lat[i].O[mu] ) ; 
	    exact_log_slow( b , a ) ; 
	  }
	  a_plus_b( stap , b ) ; 

	  //bottom staple
	  temp = lat[i].back[rho] ; 
	  multabdag_suNC( a , lev1[temp].O[kk] , lev1[temp].O[jj]  ) ; 
	  temp = lat[temp].neighbor[mu] ;  
	  multab_suNC( b , a , lev1[temp].O[kk] ) ; 
		
	  if( type == SM_LOG ) {
	    multab_dag_suNC( a , b , lat[i].O[mu] ) ; 
	    exact_log_slow( b , a ) ; 
	  }
	  a_plus_b( stap , b ) ; 
	}
	// here are the projections; SM_APE, SM_STOUT and SM_LOG
	project( lev2[i].O[ii] , stap , lat[i].O[mu] , 
		 alpha2 , one_min_a2 ) ; 
	// end of projections
	ii++ ; 
      }
    } 
  }
  return ;
}

// complete the staples ....
static void 
gen_staples_4D( GLU_complex *__restrict stap , 
		const struct lv1 *__restrict lev2 ,
		const struct site *__restrict lat ,
		const size_t i , 
		const size_t mu , 
		const int type )
{
  GLU_complex a[ NCNC ] GLUalign , b[ NCNC ] GLUalign ;
  size_t nu ;
  //calculate the staples using the dressed links
  for( nu = 0 ;  nu < ND ;  ++nu ) {
    if( likely( mu != nu ) ) {
      size_t jj = ( ND - 1 ) * mu + nu ; 
      if( nu > mu  ) { jj-- ; } 

      size_t kk = ( ND - 1 )*nu + mu ; 
      if( mu > nu  ) { kk-- ; } 

      //kk , jj , kk are the correct steps for the staples nu-mu plane
      size_t temp = lat[i].neighbor[nu] ; 
      multab_suNC( a , lev2[i].O[kk] , lev2[temp].O[jj]  ) ; 
      temp = lat[i].neighbor[mu] ; 
      multab_dag_suNC( b , a , lev2[temp].O[kk] ) ; 

      if( type == SM_LOG ) {
	multab_dag_suNC( a , b , lat[i].O[mu] ) ; 
	exact_log_slow( b , a ) ; 
      }
      a_plus_b( stap , b ) ; 

      temp = lat[i].back[nu] ; 
      multabdag_suNC( a , lev2[temp].O[kk] , lev2[temp].O[jj] ) ; 
      temp = lat[temp].neighbor[mu] ; 
      multab_suNC( b , a , lev2[temp].O[kk] ) ; 

      if( type == SM_LOG ) {
	multab_dag_suNC( a , b , lat[i].O[mu] ) ; 
	exact_log_slow( b , a ) ; 
      }
      a_plus_b( stap , b ) ; 
    }
  }
  return ;
}
#endif

////////////////////////////////////////////////////

// this code performs the smearing ...
int
HYPSLsmear4D_expensive( struct site *__restrict lat , 
			const size_t smiters , 
			const int type )
{
#if ND != 4
  return GLU_FAILURE ;
#else
  // if we want zero iterations leaving now is a success
  if( smiters < 1 ) { return GLU_SUCCESS ; }

  // check the type
  if( ( type != SM_APE ) && ( type != SM_STOUT ) && ( type != SM_LOG ) ) {
    fprintf( stderr , "[SMEAR] Unrecognised type [ %d ] ... Leaving \n" , 
	     type ) ; 
    return GLU_FAILURE ; 
  }

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
    printf( "[SMEAR] HYPSLsmear4D unrecognised type %d\n" , type ) ;
    return GLU_FAILURE ;
  }

  // allocate temporary lattices ...
  struct spt_site *lat2 = NULL , *lat3 = NULL , *lat4 = NULL ;
  if( GLU_malloc( (void**)&lat2 , 16 , LCU * sizeof( struct spt_site ) ) != 0 ||
      GLU_malloc( (void**)&lat3 , 16 , LCU * sizeof( struct spt_site ) ) != 0 ||
      GLU_malloc( (void**)&lat4 , 16 , LCU * sizeof( struct spt_site ) ) != 0 ) {
    fprintf( stderr , "[SMEARING] field allocation failure \n" ) ;
    return GLU_FAILURE ;
  }

  // allocate levels
  struct lv1 *lev1 = NULL , *lev2 = NULL ;
  if( GLU_malloc( (void**)&lev1 , 16 , LVOLUME * sizeof( struct lv1 ) ) != 0 ||
      GLU_malloc( (void**)&lev2 , 16 , LVOLUME * sizeof( struct lv1 ) ) != 0 ) {
    fprintf( stderr , "[SMEARING] field allocation failure \n" ) ;
    return GLU_FAILURE ;
  }

  // do the smearing
  size_t count = 0 ; 
  for( count = 1 ; count <= smiters ; count++ ) {
    size_t i ;
    get_lv1( lev1 , lat , type , project ) ; 
    get_lv2( lev2 , lev1 , lat , type , project ) ;
    const size_t bck = lat[0].back[ ND - 1 ] ; 
    #pragma omp parallel for private(i) SCHED
    PFOR( i = 0 ; i < LCU ; i++ ) {
      const size_t back = bck + i ;
      size_t mu ;
      GLU_complex stap[ NCNC ] GLUalign ;
      for( mu = 0 ; mu < ND ; mu++ ) {
	zero_mat( stap ) ;
	gen_staples_4D( stap , lev2 , lat , back , mu , type ) ; 
	project( lat4[i].O[mu] , stap , lat[back].O[mu] , alpha1 , 
		 one_min_a1 ) ; 
      }
    }
    //loop time slices
    ///////////////////////////////////////
    size_t t ;
    for( t = 0 ; t < Latt.dims[ ND - 1 ] - 1 ; t++ ) {
      const size_t slice = LCU * t ; 
      const size_t bck = lat[ slice ].back[ ND - 1 ] ;
      #pragma omp parallel for private(i) SCHED
      PFOR( i = 0 ; i < LCU ; i++ )  {
	const size_t it = slice + i ; 
	GLU_complex stap[ NCNC ] GLUalign ;
	size_t mu ;
	for( mu = 0 ; mu < ND ; mu++ ) {
	  zero_mat( stap ) ;
	  gen_staples_4D( stap , lev2 , lat , it , mu , type ) ; 
	  project( lat2[i].O[mu] , stap , lat[it].O[mu] , alpha1 , 
		   one_min_a1 ) ; 
	}
	// this is only a legal maneuver for this method
	//put temp into the previous time-slice
	if( likely( t != 0 ) ) { 
	  register const size_t back = bck + i ;
	  memcpy( &lat[back] , &lat3[i] , sizeof( struct spt_site ) ) ;
	}
	//make temporary lat3 lat2 again and repeat
	memcpy( &lat3[i] , &lat2[i] , sizeof( struct spt_site ) ) ;
      }
    }
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
    if( count > TOP_VALUE && count%5 == 0 ) {
      if( gauge_topological_meas( lat , &qtop_new , &qtop_old , count-1 ) 
	  == GLU_SUCCESS ) { break ; }
    }
    #endif

    #ifdef verbose
    print_smearing_obs( lat , count ) ;
    #endif
  }

  #ifndef verbose
  print_smearing_obs( lat , count ) ;
  #endif
  
  // free that memory //
  free( lev1 ) ; 
  free( lev2 ) ; 
  free( lat2 ) ; 
  free( lat3 ) ; 
  free( lat4 ) ; 

  return GLU_SUCCESS ;
#endif
}

