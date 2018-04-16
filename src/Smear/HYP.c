/*
    Copyright (2013-2018) Renwick James Hudspith

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

#if ND == 4

// 3D level-1 staples for timeslice t
static void 
get_spatial_lv1( struct s_site *__restrict lev1 ,
		 const struct site *__restrict lat ,
		 const size_t t ,
		 const int type ,
		 void (*project) ( GLU_complex smeared_link[ NCNC ] , 
				   GLU_complex staple[ NCNC ] , 
				   const GLU_complex link[ NCNC ] , 
				   const double smear_alpha , 	     
				   const double al ) ) 
{
  size_t it ; 
  const size_t slice = LCU * t ; 
  //do a slice
#pragma omp for private(it) SCHED
  PFOR( it = 0  ;  it < LCU  ;  it++ )  {
    GLU_complex a[ NCNC ] GLUalign , b[ NCNC ] GLUalign , 
      c[ NCNC ] GLUalign ;
    const size_t i = slice + it ;
    size_t j = 0 , mu , nu ; 
    //calculate the level1 staples
    for( mu = 0 ; mu < ND - 1 ; mu++ ) {
      for( nu = 0 ; nu < ND - 1 ; nu++ ) {
	if( likely( nu != mu ) ) {
	  // b is our staple
	  // j is our staple counter	    
	  size_t temp = lat[i].neighbor[nu] ; 
	  multab_suNC( a , lat[i].O[nu] , lat[temp].O[mu] ) ; 
	  temp = lat[i].neighbor[mu] ; 
	  multab_dag_suNC( b , a , lat[temp].O[nu] ) ; 
	  if( type == SM_LOG ) {
	    multab_dag_suNC( a , b , lat[i].O[mu] ) ; 
	    exact_log_slow( b , a ) ; 
	  }
	  // put the bottom staple in "c"
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
	  project( lev1[it].O[j] , b , lat[i].O[mu] , alpha2 , one_min_a2 ) ; 
	  j++ ; 
	}
      }
    }
  }
  return ;
}

static void 
staples3D( GLU_complex stap[ NCNC ] ,
	   const struct s_site *__restrict lev1 ,
	   const struct site *__restrict lat , 
	   const size_t i , 
	   const size_t mu , 
	   const size_t t , 
	   const size_t type ) 
{ 
  GLU_complex a[ NCNC ] GLUalign , b[ NCNC ] GLUalign ;
  size_t nu ;
  const size_t it = LCU * t + i ; 

  //calculate the staples using the dressed links
  for( nu = 0 ;  nu < ND - 1 ;  ++nu ) {
    if ( likely( nu != mu ) ) {		   
      size_t jj = 0 , rho = 0 ;
      // 3rd orthogonal direction: rho 
      for( jj = 0 ;  jj < ND - 1 ;  ++jj ) {
	if( jj != mu && jj != nu ) {
	  rho = jj ; 
	}
      }
	
      jj = ( ND - 2 ) * mu + rho ; 
      if( rho > mu  ) { jj-- ; } 

      size_t kk = ( ND - 2 ) * nu + rho ; 
      if( rho > nu  ) { kk-- ; }
      //kk , jj , kk are the correct steps for the staples
      size_t temp = lat[i].neighbor[nu] ; 
      multab_suNC( a , lev1[i].O[kk] , lev1[temp].O[jj] ) ; 
      temp = lat[i].neighbor[mu] ; 
      multab_dag_suNC( b , a , lev1[temp].O[kk] ) ; 
      if( type == SM_LOG ) {
	multab_dag_suNC( a , b , lat[it].O[mu] ) ; 
	exact_log_slow( b , a ) ; 
      }
      a_plus_b( stap , b ) ;
      //bottom staple
      temp = lat[i].back[nu] ; 
      multabdag_suNC( a , lev1[temp].O[kk] , lev1[temp].O[jj] ) ;
      temp = lat[temp].neighbor[mu] ; 
      multab_suNC( b , a , lev1[temp].O[kk] ) ; 
      if( type == SM_LOG ) {
	multab_dag_suNC( a , b , lat[it].O[mu] ) ; 
	exact_log_slow( b , a ) ; 
      }
      a_plus_b( stap , b ) ; 
    }
  }
  return ;
}
#endif

// spatial only smearing
int
HYPSLsmear3D( struct site *__restrict lat , 
	      const size_t smiters , 
	      const int type ) 
{
  // successfully do nothing
  if( smiters == 0 ) { return GLU_SUCCESS ; }
  // should never get here ...
#if ND != 4
  fprintf( stderr , "[SMEAR] Should not get here (HYPSLsmear3D)\n" , ND ) ;
  return GLU_FAILURE ;
#else

  struct s_site *lev1 = NULL , *lat2 = NULL ;
  double *red = NULL ;
  int FLAG = GLU_SUCCESS ;
#ifdef TOP_VALUE
  red = malloc( CLINE*Latt.Nthreads*sizeof( double ) ) ;
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
    FLAG = GLU_FAILURE ; goto memfree ;
  }

  // allocate temporary space
  if( ( lev1 = allocate_s_site( LCU , ND*(ND-1)*(ND-2) , NCNC ) ) == NULL ||
      ( lat2 = allocate_s_site( LCU , (ND-1) , NCNC ) ) == NULL ) {
    fprintf( stderr , "[SMEARING] field allocation failure\n" ) ;
    FLAG = GLU_FAILURE ; goto memfree ;
  }
 
  // iteration counter
  #pragma omp parallel
  {
    #ifdef TOP_VALUE
    double qtop_new , qtop_old = 0.0 ;
    #endif
    size_t count = 0 ;
    GLU_bool top_found = GLU_FALSE ; 
    
    for( count = 1 ; count <= smiters && top_found != GLU_TRUE ; count++ ) {
      size_t t , i ;
      //loop time slices
      for( t = 0 ; t < Latt.dims[ ND - 1 ] ; t++ ) {
	// get the level 1 links for this slice
	get_spatial_lv1( lev1  ,  lat  ,  t  ,  type , project ) ; 
	
	const size_t slice = LCU * t ; 
        #pragma omp for private(i)
	for( i = 0 ; i < LCU ; i++ )  {
	  GLU_complex stap[ NCNC ] GLUalign ;
	  const size_t it = slice + i ;
	  size_t mu ;
	  for( mu = 0 ; mu < ND - 1 ; mu++ ) {
	    zero_mat( stap ) ;
	    staples3D( stap , lev1 , lat , i , mu , t , type ) ; 
	    project( lat2[ i ].O[ mu ] , stap , lat[ it ].O[ mu ] , 
		     alpha1 , one_min_a1 ) ; 
	  }
	  // swap these round
	  for( mu = 0 ; mu < ND-1 ; mu++ ) {
	    equiv( lat[it].O[mu] , lat2[i].O[mu] ) ;
	  }
	  //
	}
      }

      #ifdef TOP_VALUE
      if( count > TOP_VALUE && count%5 == 0 ) {
        #pragma omp for private(i)
	for( i = 0 ; i < CLINE*Latt.Nthreads ; i++ ) {
	  red[i] = 0.0 ;
	}
	if( gauge_topological_meas_th( red , lat , &qtop_new ,
				       &qtop_old , count-1 ) 
	    == GLU_SUCCESS ) {
	  top_found = GLU_TRUE ;
	}
      }
      #endif

      // only write these out if we are not doing this in parallel ...
      #ifdef verbose
      print_smearing_obs( lat , count ) ;
      #endif
    }

#ifndef verbose
  print_smearing_obs( lat , count ) ;
#endif
  }

 memfree :

  if( red != NULL ) {
    free( red ) ;
  }
  free_s_site( lat2 , LCU , ND-1 , NCNC ) ;
  free_s_site( lev1 , LCU , ND*(ND-1)*(ND-2) , NCNC ) ;

  return FLAG ;
#endif
}
