/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (staples.c) is part of GLU.

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
   @file staples.c
   @brief elementary and tree level improved staple calculators used by smear.c and wflow.c
 */
#include "Mainfile.h"

// staple routine the positive \nu one goes like
//
//         ---<---
//         |      |
//         v      ^
//         |      |
void
all_staples( GLU_complex stap[ NCNC ] , 
	     const struct site *__restrict lat , 
	     const size_t i , 
	     const size_t mu , 
	     const size_t dir , 
	     const smearing_types type )
{
  GLU_complex a[ NCNC ] GLUalign , b[ NCNC ] GLUalign ; 
  register size_t nu , t1 , t2 , b1 , b2 , j ; 
  switch( type ) {
  case SM_LOG :
    for( j = 0 ; j < dir-1 ; j++ ) {
      nu = j < mu ? j : j + 1 ;
      // directions
      t1 = lat[ i ].neighbor[nu] ; 
      t2 = lat[ i ].neighbor[mu] ; 
      b1 = lat[ i ].back[nu] ; 
      b2 = lat[ b1 ].neighbor[mu] ; 
      // top staple
      multab_suNC( a , lat[ i ].O[nu] , lat[ t1 ].O[mu] ) ; 
      multab_dag_suNC( b , a , lat[ t2 ].O[nu] ) ; 
      multab_dag_suNC( a , b , lat[i].O[mu] ) ; 
      exact_log_slow( b , a ) ; 
      a_plus_b( stap , b ) ; 
      //bottom staple
      multabdag_suNC( a , lat[ b1 ].O[nu] , lat[ b1 ].O[mu] ) ; 
      multab_suNC( b , a , lat[ b2 ].O[nu] ) ; 
      multab_dag_suNC( a , b , lat[i].O[mu] ) ; 
      exact_log_slow( b , a ) ; 
      a_plus_b( stap , b ) ; 
    } break ;
  default :
    for( j = 0 ; j < dir-1 ; j++ ) {
      nu = j < mu ? j : j + 1 ;
      // directions
      t1 = lat[ i ].neighbor[nu] ; 
      t2 = lat[ i ].neighbor[mu] ; 
      b1 = lat[ i ].back[nu] ; 
      b2 = lat[ b1 ].neighbor[mu] ; 
      // top staple
      multab_suNC( a , lat[ i ].O[nu] , lat[ t1 ].O[mu] ) ; 
      multab_dag_suNC( b , a , lat[ t2 ].O[nu] ) ; 
      a_plus_b( stap , b ) ; 
      //bottom staple
      multabdag_suNC( a , lat[ b1 ].O[nu] , lat[ b1 ].O[mu] ) ; 
      multab_suNC( b , a , lat[ b2 ].O[nu] ) ; 
      a_plus_b( stap , b ) ; 
    } break ;
  }
  return ;
}

// (OVER)IMPROVED STAPLES ROUTINE
#ifdef IMPROVED_SMEARING
void
all_staples_improve( GLU_complex stap[ NCNC ] , 
		     const struct site *__restrict lat , 
		     const size_t i , 
		     const size_t mu , 
		     const size_t dir , 
		     const smearing_types type )
{
  GLU_complex a[ NCNC ] GLUalign , b[ NCNC ] GLUalign ;
  GLU_complex c0_stap[ NCNC ] GLUalign , c1_stap[ NCNC ] GLUalign ;
  GLU_complex tempstap[ NCNC ] GLUalign ;

  zero_mat( c0_stap ) ; zero_mat( c1_stap ) ;

#ifdef SYMANZIK_ONE_LOOP
  GLU_complex c2_stap[ NCNC ] ;
  zero_mat( c2_stap ) ;
#endif 

  size_t nu ; 
  for( nu = 0 ; nu < dir ; nu++ ) {
    if( nu == mu ) continue ;
    //top staple
    register const size_t t1 = lat[i].neighbor[nu] ; 
    multab_suNC( a , lat[i].O[nu] , lat[ t1 ].O[mu] ) ; 
    register const size_t t2 = lat[i].neighbor[mu] ; 
    multab_dag_suNC( b , a , lat[ t2 ].O[nu] ) ; 
    if( type == SM_LOG ) {
      multab_dag_suNC( a , b , lat[i].O[mu] ) ; 
      exact_log_slow( b , a ) ; 
    }
    a_plus_b( c0_stap , b ) ;
    
    //bottom staple
    const size_t b1 = lat[i].back[nu] ; 
    multabdag_suNC( a , lat[ b1 ].O[nu] , lat[ b1 ].O[mu] ) ; 
    const size_t b2 = lat[ b1 ].neighbor[mu] ; 
    multab_suNC( b , a , lat[ b2 ].O[nu] ) ; 
    if( type == SM_LOG )  {
      multab_dag_suNC( a , b , lat[i].O[mu] ) ; 
      exact_log_slow( b , a ) ; 
    }
    a_plus_b( c0_stap , b ) ;
    
    // IMPROVEMENT factors
    // this is the ( nu , nu , mu , -nu , -nu )
    //
    // The ==<== is the dagger of the link U_\mu(x)
    //
    // x-->--x          x==<==x
    // |     |          |     |
    // ^     v          v     ^
    // |     |          |     |
    // x     x     +    x     x  
    // |     |          |     |
    // ^     v          v     ^
    // |     |          |     |
    // x==<==x          x-->--x
    register const size_t tv1 = lat[ i ].neighbor[nu] ; 
    multab_suNC( a , lat[ i ].O[nu] , lat[ tv1 ].O[nu] ) ; 
    register const size_t tv2 = lat[ tv1 ].neighbor[nu] ; 
    multab_suNC( b , a , lat[ tv2 ].O[mu] ) ; 
    // halfway there ...
    register const size_t tv3 = lat[ i ].neighbor[ mu ] ; 
    register const size_t tv4 = lat[ tv3 ].neighbor[ nu ] ; 
    multab_dag_suNC( a , b , lat[ tv4 ].O[ nu ] ) ; 
    multab_dag_suNC( b , a , lat[ tv3 ].O[ nu ] ) ; 
    if( type == SM_LOG ) {
      multab_dag_suNC( a , b , lat[i].O[mu] ) ; 
      exact_log_slow( b , a ) ; 
    }
    a_plus_b( c1_stap , b ) ;
    
    // include the lower half contributions .....
    register const size_t bv1 = lat[ i ].back[nu] ; 
    register const size_t bv2 = lat[ bv1 ].back[nu] ; 
    multab_dagdag_suNC( a , lat[ bv1 ].O[nu] , lat[ bv2 ].O[nu] ) ; 
    multab_suNC( b , a , lat[ bv2 ].O[mu] ) ; 
    register const size_t bv3 = lat[ bv2 ].neighbor[mu] ; 
    multab_suNC( a , b , lat[ bv3 ].O[nu] ) ; 
    register const size_t bv4 = lat[ bv3 ].neighbor[nu] ; 
    multab_suNC( b , a , lat[ bv4 ].O[nu] ) ; 
    if( type == SM_LOG ) {
      multab_dag_suNC( a , b , lat[i].O[mu] ) ; 
      exact_log_slow( b , a ) ; 
    }
    a_plus_b( c1_stap , b ) ;

    // Horizontal, forwards contribution ....
    //  x-->--x-->--x        x==<==x--<--x
    //  |           |        |           |
    //  ^           v   +    v           ^
    //  |           |        |           |
    //  x==<==x--<--<        x-->--x-->--x
    // The log is a pain as we usually can eliminate 1 matrix multiply
    register const size_t tf1 = lat[ i ].neighbor[nu] ; 
    multab_suNC( a , lat[ i ].O[nu] , lat[ tf1 ].O[mu] ) ; 
    register const size_t tf2 = lat[ tf1 ].neighbor[mu] ; 
    multab_suNC( b , a , lat[ tf2 ].O[mu] ) ; 
    register const size_t tf3 = lat[ i ].neighbor[ mu ] ; 
    register const size_t tf4 = lat[ tf3 ].neighbor[ mu ] ; 
    multab_dag_suNC( a , b , lat[ tf4 ].O[ nu ] ) ; 
    if( type == SM_LOG ) {
      multab_dag_suNC( b , a , lat[ tf3 ].O[ mu ] ) ; 
      multab_dag_suNC( a , b , lat[i].O[mu] ) ; 
      exact_log_slow( b , a ) ; 
      a_plus_b( c1_stap , b ) ;
    } else {
      equiv( tempstap , a ) ;
    }
    
    // bottom rectangle
    register const size_t bf1 = lat[ i ].back[nu] ; 
    multabdag_suNC( a , lat[ bf1 ].O[nu] , lat[ bf1 ].O[mu] ) ; 
    register const size_t bf2 = lat[ bf1 ].neighbor[mu] ; 
    multab_suNC( b , a , lat[ bf2 ].O[mu] ) ; 
    register const size_t bf3 = lat[ bf2 ].neighbor[mu] ; 
    multab_suNC( a , b , lat[ bf3 ].O[nu] ) ; 
    if( type == SM_LOG ) {
      multab_dag_suNC( b , a , lat[ tf3 ].O[ mu ] ) ; 
      multab_dag_suNC( a , b , lat[i].O[mu] ) ; 
      exact_log_slow( b , a ) ; 
    } else  {
      a_plus_b( tempstap , a ) ;
      multab_dag( b , tempstap , lat[ tf3 ].O[mu] ) ;
    }
    a_plus_b( c1_stap , b ) ;

    //// The horizontal, backwards contribution
    //
    // x-->--x-->--x       x--<--x==<==x
    // |           |       |           |
    // ^           v   +   v           ^
    // |           |       |           |
    // x--<--x==<==x       x-->--x-->--x
    // this one is split
    if( type == SM_LOG ) {
      // do the top staple
      register const size_t tb1 = lat[ i ].back[mu] ; 
      multabdag_suNC( a , lat[ tb1 ].O[mu] , lat[ tb1 ].O[nu] ) ; 
      register const size_t tb2 = lat[ tb1 ].neighbor[nu] ; 
      multab_suNC( b , a , lat[ tb2 ].O[mu] ) ; 
      // halfway there ...
      register const size_t tb3 = lat[ tb2 ].neighbor[mu] ; 
      multab_suNC( a , b , lat[ tb3 ].O[ mu ] ) ; 
      register const size_t tb4 = lat[ i ].neighbor[ mu ] ; 
      multab_dag_suNC( b , a , lat[ tb4 ].O[ nu ] ) ; 
      multab_dag_suNC( a , b , lat[i].O[mu] ) ; 
      exact_log_slow( b , a ) ; 
      a_plus_b( c1_stap , b ) ;
      
      // and the bottom staple ...
      register const size_t bb1 = lat[ i ].back[mu] ; 
      register const size_t bb2 = lat[ bb1 ].back[nu] ; 
      multab_dagdag_suNC( a , lat[ bb1 ].O[mu] , lat[ bb2 ].O[nu] ) ; 
      multab_suNC( b , a , lat[ bb2 ].O[mu] ) ; 
      const register size_t bb3 = lat[ bb2 ].neighbor[mu] ; 
      multab_suNC( a , b , lat[ bb3 ].O[mu] ) ; 
      const register size_t bb4 = lat[ bb3 ].neighbor[mu] ; 
      multab_suNC( b , a , lat[ bb4 ].O[nu] ) ; 
      multab_dag_suNC( a , b , lat[i].O[mu] ) ; 
      exact_log_slow( b , a ) ; 
    } else {
      register const size_t tb1 = lat[i].back[mu] ; 
      register const size_t tb2 = lat[ tb1 ].neighbor[nu] ; 
      multab_suNC( b , lat[ tb1 ].O[nu] , lat[ tb2 ].O[mu] ) ; 
      register const size_t tb3 = lat[ tb2 ].neighbor[mu] ; 
      multab_suNC( a , b , lat[ tb3 ].O[ mu ] ) ; 
      register const size_t tb4 = lat[ i ].neighbor[ mu ] ; 
      multab_dag_suNC( b , a , lat[ tb4 ].O[ nu ] ) ; 
      equiv( tempstap , b ) ;

      // bottom rectangle
      register const size_t bb2 = lat[ tb1 ].back[nu] ; 
      multabdag_suNC( b , lat[ bb2 ].O[nu] , lat[ bb2 ].O[mu] ) ; 
      register const size_t bb3 = lat[ bb2 ].neighbor[mu] ; 
      multab_suNC( a , b , lat[ bb3 ].O[mu] ) ; 
      register const size_t bb4 = lat[ bb3 ].neighbor[mu] ; 
      multab_suNC( b , a , lat[ bb4 ].O[nu] ) ; 
      a_plus_b( tempstap , b ) ; 
      multabdag( b , lat[ tb1 ].O[mu] , tempstap ) ;
    }
    a_plus_b( c1_stap , b ) ;

    // ONE LOOP SYMANZIK EXTRA TERM
    //
    // is of the form 
    //
    //      x-->--x          x==<==x
    //     /     /           |     |
    //    ^     v            v     ^
    //   /     /             |     |
    //  x     x      +       x     x
    //  |     |             /     /
    //  ^     v            v     ^
    //  |     |           /     /
    //  x==<==x          x-->--x
#ifdef SYMANZIK_ONE_LOOP
    size_t rho ;
    for( rho = 0 ; rho < dir ; rho++ ) {
      if( rho != mu && rho != nu ) {
	// include the term (nu,rho,mu,-\rho,-\nu)
	//top staple
	size_t temp = lat[i].neighbor[nu] ; 
	multab_suNC( a , lat[i].O[nu] , lat[temp].O[rho] ) ; 
	temp = lat[temp].neighbor[rho] ; 
	multab_suNC( b , a , lat[temp].O[mu] ) ; 
	// halfway there ...
	size_t temp2 = lat[ i ].neighbor[ mu ] ; 
	temp = lat[ temp2 ].neighbor[ rho ] ; 
	multab_dag_suNC( a , b , lat[ temp ].O[ nu ] ) ; 
	multab_dag_suNC( b , a , lat[ temp2 ].O[ rho ] ) ; 
	
	if( type == SM_LOG ) {
	  multab_dag_suNC( a , b , lat[i].O[mu] ) ; 
	  exact_log_slow( b , a ) ; 
	}
	a_plus_b( c2_stap , b ) ; 
	
	//include the bottom staple
	temp = lat[i].back[nu] ; 
	temp2 = lat[temp].back[rho] ; 
	multab_dagdag_suNC( a , lat[temp].O[nu] , lat[temp2].O[rho] ) ; 
	multab_suNC( b , a , lat[temp2].O[mu] ) ; 
	// agree
	temp = lat[temp2].neighbor[mu] ; 
	multab_suNC( a , b , lat[temp].O[nu] ) ; 
	temp2 = lat[temp].neighbor[nu] ; 
	multab_suNC( b , a , lat[temp2].O[rho] ) ; 
	if( type == SM_LOG ) {
	  multab_dag_suNC( a , b , lat[i].O[mu] ) ; 
	  exact_log_slow( b , a ) ; 
	}
	a_plus_b( c2_stap , b ) ; 
	// thats the end of this term ....
      }
    }
#endif
  }

#ifdef SYMANZIK_ONE_LOOP
  const GLU_real c0_multiplier = SYMONELOOP_WEIGHT1 ; 
  const GLU_real c1_multiplier = SYMONELOOP_WEIGHT2 ; 
  const GLU_real c2_multiplier = SYMONELOOP_WEIGHT3 ; 
#else
  const GLU_real c0_multiplier = IWA_WEIGHT1 ; 
  const GLU_real c1_multiplier = IWA_WEIGHT2 ; 
#endif
#if NC == 3
  if( type == SM_LOG ) {
    stap[0] = c0_multiplier * c0_stap[0] + c1_multiplier * c1_stap[0] ;
    stap[1] = c0_multiplier * c0_stap[1] + c1_multiplier * c1_stap[1] ;
    stap[2] = c0_multiplier * c0_stap[2] + c1_multiplier * c1_stap[2] ;
    stap[3] = conj( stap[1] ) ;
    stap[4] = c0_multiplier * c0_stap[4] + c1_multiplier * c1_stap[4] ;
    stap[5] = c0_multiplier * c0_stap[5] + c1_multiplier * c1_stap[5] ;
    stap[6] = conj( stap[2] ) ;
    stap[7] = conj( stap[5] ) ;
    stap[8] = -( stap[0] + stap[4] ) ;
  } else {
    stap[0] = c0_multiplier * c0_stap[0] + c1_multiplier * c1_stap[0] ;
    stap[1] = c0_multiplier * c0_stap[1] + c1_multiplier * c1_stap[1] ;
    stap[2] = c0_multiplier * c0_stap[2] + c1_multiplier * c1_stap[2] ;
    stap[3] = c0_multiplier * c0_stap[3] + c1_multiplier * c1_stap[3] ;
    stap[4] = c0_multiplier * c0_stap[4] + c1_multiplier * c1_stap[4] ;
    stap[5] = c0_multiplier * c0_stap[5] + c1_multiplier * c1_stap[5] ;
    stap[6] = c0_multiplier * c0_stap[6] + c1_multiplier * c1_stap[6] ;
    stap[7] = c0_multiplier * c0_stap[7] + c1_multiplier * c1_stap[7] ;
    stap[8] = c0_multiplier * c0_stap[8] + c1_multiplier * c1_stap[8] ;
  }
#elif NC == 2
  if( type == SM_LOG ) {
    stap[0] = c0_multiplier * c0_stap[0] + c1_multiplier * c1_stap[0] ;
    stap[1] = c0_multiplier * c0_stap[1] + c1_multiplier * c1_stap[1] ;
    stap[2] = conj( stap[1] ) ;
    stap[3] = -stap[0] ;
  } else {
    stap[0] = c0_multiplier * c0_stap[0] + c1_multiplier * c1_stap[0] ;
    stap[1] = c0_multiplier * c0_stap[1] + c1_multiplier * c1_stap[1] ;
    stap[2] = c0_multiplier * c0_stap[2] + c1_multiplier * c1_stap[2] ;
    stap[3] = c0_multiplier * c0_stap[3] + c1_multiplier * c1_stap[3] ;
  }
#else
  for( nu = 0 ; nu < NCNC ; nu++ ) {
    stap[nu] = c0_multiplier * c0_stap[nu] + c1_multiplier * c1_stap[nu] ;
  }
#endif
  
  // include the extra symanzik one loop term if that is what you fancy ...
#ifdef SYMANZIK_ONE_LOOP
  for( nu = 0 ; nu < NCNC ; nu++ ) {
    stap[nu] = stap[nu] + c2_multiplier * c2_stap[nu] ;
  }
#endif

  return ;
}

#endif
