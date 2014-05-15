/*
    Copyright 2013 Renwick James Hudspith

    This file (ND_generic_HYP.c) is part of GLU.

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
   @file ND_generic_HYP.c
   @brief very slow recursion of links
   @ingroup Smear
   HEX smearing alphas are related to HYP's
   via HYP::HYP (a1,a2,a3) = (3a1,2a2,a3)
   <br>
   Computes the links when needed rather than precomputation
 */

#include "Mainfile.h"
#include "plaqs_links.h"
#include "projectors.h"

// If we are using the dangerous smearing routines ...
#ifdef FAST_SMEAR
  #include "random_config.h"
#endif

// global maximum smearing direction ( ND-1 == SPATIAL , ND == ALL_DIRECTIONS )
static int MAXDIR ;
static GLU_real smear_alphas[ ND ] , one_minus_smalpha[ ND ] ;

// just for simplicity
enum{ NOT_ORTHOGONAL , ORTHOGONAL } ;

// initialise the smoothing alphas
static void
init_smearing_alphas( )
{
#if ND < 3
  int mu ;
  for( mu = 0 ; mu < MAXDIR-1 ; mu++ ) {
    smear_alphas[ mu+1 ] = Latt.sm_alpha[mu] ;
    one_minus_smalpha[ mu+1 ] = ( 1.0 - Latt.sm_alpha[mu] ) ;
  }
#else
  int mu ;
  for( mu = 0 ; mu < MAXDIR-1 ; mu++ ) {
    smear_alphas[ mu+1 ] = Latt.sm_alpha[mu] / ( ( ND-mu-1 ) * ( ND-2 ) ) ;
    one_minus_smalpha[ mu+1 ] = ( 1.0 - Latt.sm_alpha[mu] ) ;
  }
#endif
  return ;
}

// is the direction orthogonal?
static int
is_orthogonal( jj , size , list_dirs ) 
     const int jj , size , list_dirs[ size ] ;
{
  int nu ;
  for( nu = 0 ; nu < size ; nu++ ) {
    if( list_dirs[nu] == jj ) {
      return NOT_ORTHOGONAL ;
    } 
  }
  return ORTHOGONAL ;
}

// passes link by reference
static void
recurse_staples( link , lat , i , lev , list_dirs , type )
     GLU_complex *__restrict link ; 
     const struct site *__restrict lat ; 
     int lev ;
     const int list_dirs[ MAXDIR - lev ] ; 
     const int type ;
{
  // the last index is our rho plane
  const int rho = list_dirs[ MAXDIR - lev - 1 ] ;
  // generic storage and stuff
  GLU_complex stap[ NCNC ] = { } ;
  GLU_complex a[ NCNC ] , b[ NCNC ] ;

  if( lev == 1 ) { 
    // this is the final one
    int jj ;
    for( jj = 0 ; jj < MAXDIR ; jj++ ) {
      if( is_orthogonal( jj , MAXDIR-lev , list_dirs ) == ORTHOGONAL ) { 
	//jj is our orthogonal direction
	int temp = lat[i].neighbor[jj] ; 

	multab_suNC( a , lat[i].O[jj] , lat[temp].O[rho] ) ; 
	temp = lat[i].neighbor[rho] ; 
	multab_dag_suNC( b , a , lat[temp].O[jj] ) ; 

	if( type == SM_LOG ) {
	  multab_dag_suNC( a , b , lat[i].O[rho] ) ; 
          #ifdef FAST_SMEAR
	  exact_log_fast( b , a ) ; 
          #else
	  exact_log_slow( b , a ) ; 
          #endif
	}
	a_plus_b( stap , b ) ; 
      
	//bottom staple
	temp = lat[i].back[jj] ; 
	multabdag_suNC( a , lat[temp].O[jj] , lat[temp].O[rho] ) ; 
	temp = lat[temp].neighbor[rho] ; 
	multab_suNC( b , a , lat[temp].O[jj] ) ; 

	if( type == SM_LOG ) {
	  multab_dag_suNC( a , b , lat[i].O[rho] ) ; 
          #ifdef FAST_SMEAR
	  exact_log_fast( b , a ) ; 
          #else
	  exact_log_slow( b , a ) ; 
          #endif
	}
	a_plus_b( stap , b ) ; 
      }
    }
    // and then we traverse back up the list
  } else {

    // allocate these two temporaries
    GLU_complex temp[ NCNC ] , temp2[ NCNC ] ;
    int new_list_dirs[ MAXDIR - lev + 1 ] , orthogonal_dirs[ MAXDIR - lev + 1 ] , jj ;
    for( jj = 0 ; jj < MAXDIR ; jj++ ) {
      if( is_orthogonal( jj , MAXDIR-lev , list_dirs ) == ORTHOGONAL ) { 

 	// add direction to the list and recurse
	int nu ;
	for( nu = 0 ; nu < MAXDIR-lev ; nu++ ) {
	  new_list_dirs[nu] = list_dirs[nu] ;
	  orthogonal_dirs[nu] = list_dirs[nu] ;
	}
	orthogonal_dirs[nu-1] = jj ;
	orthogonal_dirs[nu] = rho ; 
	new_list_dirs[nu] = jj ;

	int dir = lat[i].neighbor[jj] ; 
	// first element of the staple is in the jj - rho direction
	recurse_staples( temp , lat , i , lev-1 , new_list_dirs , type ) ; 
	recurse_staples( temp2 , lat , dir , lev-1 , orthogonal_dirs , type ) ; 
	multab_suNC( a , temp , temp2 ) ; 

	dir = lat[i].neighbor[rho] ;
	recurse_staples( temp , lat , dir , lev-1 , new_list_dirs , type ) ; 
	multab_dag_suNC( b , a , temp ) ; 

	if( type == SM_LOG ) {
	  multab_dag_suNC( a , b , lat[i].O[rho] ) ; 
          #ifdef FAST_SMEAR
	  exact_log_fast( b , a ) ; 
          #else
	  exact_log_slow( b , a ) ; 
          #endif
	}
	a_plus_b( stap , b ) ; 
	// end of top staple ...

	dir = lat[i].back[jj] ; 
	recurse_staples( temp , lat , dir , lev-1 , new_list_dirs , type ) ; 
	recurse_staples( temp2 , lat , dir , lev-1 , orthogonal_dirs , type ) ;
	multabdag_suNC( a , temp , temp2 ) ; 

	dir = lat[dir].neighbor[rho] ; 
	recurse_staples( temp , lat , dir , lev-1 , new_list_dirs , type ) ;
	multab_suNC( b , a , temp ) ; 

	if( type == SM_LOG ) {
	  multab_dag_suNC( a , b , lat[i].O[rho] ) ; 
          #ifdef FAST_SMEAR
	  exact_log_fast( b , a ) ; 
          #else
	  exact_log_slow( b , a ) ; 
          #endif
	}
	a_plus_b( stap , b ) ; 
	// end of the bottom staple ...
      }
    }
    // and so it should recurse inside this
  }
  // OK and we now project the fields 
  switch( type )
    {
    case SM_APE :
      project_APE( link , stap , lat[i].O[rho] , 
		   smear_alphas[ MAXDIR-lev ] , one_minus_smalpha[ MAXDIR-lev ] ) ; 
      break ; 
    case SM_STOUT :
      project_STOUT_short( link , stap , lat[i].O[rho] , 
			   smear_alphas[ MAXDIR-lev ] ) ; 
      break ; 
    case SM_LOG :
      project_LOG_short( link , stap , lat[i].O[rho] , 
			 smear_alphas[ MAXDIR-lev ] ) ; 
      break ; 
    }
  return ;
}

// slow, recursive ND-Generic smearing routine
void 
HYsmearND( struct site *__restrict lat , 
	   const int smiters , 
	   const int type ,
	   const int directions )
{
  if( unlikely( smiters == 0 ) ) { return ; }

  // check the type , just in case
  if( ( type != SM_APE ) && ( type != SM_STOUT ) && ( type != SM_LOG ) ) {
    printf( "[SMEAR] Unrecognised type [ %d ] ... Leaving \n" , type ) ; 
    return ; 
  }

#ifdef TOP_VALUE
  double qtop_new , qtop_old = 0.0 ;
#endif

  // globalise the smearing directions
  MAXDIR = directions ;

  // initialise our smearing parameters and print their values to the screen
  init_smearing_alphas( ) ;

  struct spt_site *lat2 = malloc( LCU * sizeof( struct spt_site ) ) ; 
  struct spt_site *lat3 = malloc( LCU * sizeof( struct spt_site ) ) ; 
  struct spt_site *lat4 = malloc( LCU * sizeof( struct spt_site ) ) ; 

  const int lev = MAXDIR - 1 ;
  int count = 0 ; 
  for( count = 1 ; count <= smiters ; count++ ) {
    const int back = lat[ 0 ].back[ ND - 1 ] ;
    int i ;
    #pragma omp parallel for private(i) SCHED
    PFOR( i = 0 ; i < LCU ; i++ ) {
      const int bck = back + i ;
      int mu , list_dirs[ MAXDIR - lev ] ;
      for( mu = 0 ; mu < MAXDIR ; mu++ ) {
	int d = 0 ;
	for( d = 0 ; d < MAXDIR-lev ; d++ ) { list_dirs[d] = mu ; }
	recurse_staples( lat4[i].O[mu] , lat , bck , lev , list_dirs , type ) ; 
      }
    }

    //loop time slice
    ///////////////////////////////////////
    int t ;
    for( t = 0 ; t < ( Latt.dims[ ND - 1 ] - 1 ) ; t++ ) {
      const int slice = LCU * t ;

      #pragma omp parallel for private(i) SCHED
      PFOR( i = 0 ; i < LCU ; i++ ) {
	const int it = slice + i ; 
	int mu , list_dirs[ MAXDIR - lev ] ;
	for( mu = 0 ; mu < MAXDIR ; mu++ ) {
	  int d = 0 ;
	  for( d = 0 ; d < MAXDIR-lev ; d++ ) { list_dirs[d] = mu ; }
	  recurse_staples( lat2[i].O[mu] , lat , it , lev , list_dirs , type ) ; 
	}
      }
	  
      const int bck = lat[slice].back[ ND - 1 ] ;
      #pragma omp parallel for private(i)
      PFOR( i = 0 ; i < LCU ; i++ ) {
	int mu ;
	//put temp into the previous time-slice 
	if( likely( t != 0 ) ) { 
	  register const int back = bck + i ;
	  for( mu = 0 ; mu < MAXDIR ; mu++ ) {
	    equiv( lat[back].O[mu] , lat3[i].O[mu] ) ;
	  }
	}
	//make temporary lat3 lat2 again and repeat
	for( mu = 0 ; mu < MAXDIR ; mu++ ) {
	  equiv( lat3[i].O[mu] , lat2[i].O[mu] ) ;
	}
      }
    }
    //put last and last but one time slice in
    ////////////////////////////////////////////
    const int slice = LCU * t ;    
    const int behind = lat[ slice ].back[ ND - 1 ] ;
    #pragma omp parallel for private(i) 
    PFOR( i = 0 ; i < LCU ; i++ ) {
      int mu ;
      register const int back = behind + i ; 
      register const int it = slice + i ; 
      for( mu = 0 ; mu < MAXDIR ; mu++ ) {
	equiv( lat[back].O[mu] , lat3[i].O[mu] ) ;
	equiv( lat[it].O[mu] , lat4[i].O[mu] ) ;
      }
    }

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

  // need to clean up the lattice fields
  #ifdef FAST_SMEAR
  if( type == SM_STOUT ) {
    latt_reunitU( lat ) ;
    printf( "\n[SMEAR] A final reunitarisation step to clean things up\n" ) ;
    print_smearing_obs( lat , type , count , GLU_TRUE ) ;
  }
  #endif

  // FREE STUFF !! //
  free( lat2 ) ; 
  free( lat3 ) ; 
  free( lat4 ) ; 

  return ;
}
