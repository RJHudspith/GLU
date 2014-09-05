/*
    Copyright 2013 Renwick James Hudspith

    This file (log_derivs.c) is part of GLU.

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
  @file log_derivs.c
  @brief finite difference calculators ...

  I have included the Log and definitions
  as well as the stencil terms
 */

#include "Mainfile.h"

#include "lin_derivs.h" // for trace deriv

// approximate log definition to start off with ....
double
approx_log_deriv( GLU_complex sum[ HERMSIZE ] , 
		  const struct site *__restrict lat , 
		  const int i , 
		  const int MAX_DIR )
{
  GLU_complex A[ HERMSIZE ] , shiftA[ HERMSIZE ] ;
  int mu ; 
  for( mu = 0 ; mu < MAX_DIR ; mu++ ) {
    exact_log_fast_short( A , lat[i].O[mu] ) ; 
    exact_log_fast_short( shiftA , lat[lat[i].back[mu]].O[mu] ) ; 
    a_plus_Sxbminc_short( sum , 1.0 , shiftA , A ) ;
  }
  return trace_deriv( sum ) ; 
}
// approximate log definition of the nn type
double
approx_log_deriv_nn( GLU_complex sum[ HERMSIZE ] , 
		     const struct site *__restrict lat , 
		     const int i , 
		     const int MAX_DIR )
{
  GLU_complex A[ HERMSIZE ] , shiftA[ HERMSIZE ] ;
  int mu ; 
  for( mu = 0 ; mu < MAX_DIR ; mu++ ) {
    exact_log_fast_short( A , lat[i].O[mu] ) ; 
    exact_log_fast_short( shiftA , lat[lat[i].back[mu]].O[mu] ) ; 
    a_plus_Sxbminc_short( sum , nn1 , shiftA , A ) ;

    exact_log_fast_short( A , lat[lat[i].neighbor[mu]].O[mu] ) ; 
    exact_log_fast_short( shiftA , lat[lat[lat[i].back[mu]].back[mu]].O[mu] ) ; 
    a_plus_Sxbminc_short( sum , nn2 , shiftA , A ) ;
  }
  return trace_deriv( sum ) ; 
}

// approximate log definition of the nn type
double
approx_log_deriv_nnn( GLU_complex sum[ HERMSIZE ] , 
		      const struct site *__restrict lat , 
		      const int i , 
		      const int MAX_DIR )
{
  GLU_complex A[ HERMSIZE ] , shiftA[ HERMSIZE ] ;
  int mu ; 
  for( mu = 0 ; mu < MAX_DIR ; mu++ ) {
    // first deriv
    exact_log_fast_short( A , lat[i].O[mu] ) ; 
    exact_log_fast_short( shiftA , lat[lat[i].back[mu]].O[mu] ) ; 
    a_plus_Sxbminc_short( sum , nnn1 , shiftA , A ) ;

    exact_log_fast_short( A , lat[lat[i].neighbor[mu]].O[mu] ) ; 
    exact_log_fast_short( shiftA , lat[lat[lat[i].back[mu]].back[mu]].O[mu] ) ; 
    a_plus_Sxbminc_short( sum , nnn2 , shiftA , A ) ;

    const int fi = lat[lat[i].neighbor[mu]].neighbor[mu] ;
    const int bi = lat[lat[lat[i].back[mu]].back[mu]].back[mu] ;
    exact_log_fast_short( A , lat[fi].O[mu] ) ; 
    exact_log_fast_short( shiftA , lat[bi].O[mu] ) ; 
    a_plus_Sxbminc_short( sum , nnn3 , shiftA , A ) ;
  }
  return trace_deriv( sum ) ; 
}

// log definition of the gauge matrices at the point i
double
log_deriv( GLU_complex sum[ HERMSIZE ] , 
	   double *functional ,
	   const struct site *__restrict lat , 
	   const int i , 
	   const int MAX_DIR )
{
  GLU_complex A[ HERMSIZE ] , shiftA[ HERMSIZE ] ;
  GLU_real trAA ;
  *functional = 0.0 ;

  int mu ; 
  for( mu = 0 ; mu < MAX_DIR ; mu++ ) {
    exact_log_slow_short( A , lat[i].O[mu] ) ; 
    exact_log_slow_short( shiftA , lat[lat[i].back[mu]].O[mu] ) ; 
    a_plus_Sxbminc_short( sum , 1.0 , shiftA , A ) ;

    // compute the functional in-step to remove a log
    trace_ab_herm_short( &trAA , A , A ) ;
    *functional += (double)trAA ;
  }
  return trace_deriv( sum ) ; 
}

// log definition of the outer and inner derivatives
double
log_deriv_nn( GLU_complex sum[ HERMSIZE ] , 
	      const struct site *__restrict lat , 
	      const int i , 
	      const int MAX_DIR )
{
  GLU_complex shiftA[ HERMSIZE ] , A[ HERMSIZE ] ;
  int mu ; 
  for( mu = 0 ; mu < ND ; mu++ ) {
    // first deriv
    exact_log_slow_short( A , lat[i].O[mu] ) ; 
    exact_log_slow_short( shiftA , lat[lat[i].back[mu]].O[mu] ) ; 
    a_plus_Sxbminc_short( sum , nn1 , shiftA , A ) ;

    exact_log_slow_short( A , lat[lat[i].neighbor[mu]].O[mu] ) ; 
    exact_log_slow_short( shiftA , lat[lat[lat[i].back[mu]].back[mu]].O[mu] ) ; 
    a_plus_Sxbminc_short( sum , nn2 , shiftA , A ) ;
  }
  return trace_deriv( sum ) ; 
}

// log definition of the outer and inner derivatives
double
log_deriv_nnn( GLU_complex sum[ HERMSIZE ] , 
	      const struct site *__restrict lat , 
	      const int i , 
	      const int MAX_DIR )
{
  GLU_complex shiftA[ HERMSIZE ] , A[ HERMSIZE ] ;
  int mu ; 
  for( mu = 0 ; mu < ND ; mu++ ) {
    // first deriv
    exact_log_slow_short( A , lat[i].O[mu] ) ; 
    exact_log_slow_short( shiftA , lat[lat[i].back[mu]].O[mu] ) ; 
    a_plus_Sxbminc_short( sum , nnn1 , shiftA , A ) ;

    exact_log_slow_short( A , lat[lat[i].neighbor[mu]].O[mu] ) ; 
    exact_log_slow_short( shiftA , lat[lat[lat[i].back[mu]].back[mu]].O[mu] ) ; 
    a_plus_Sxbminc_short( sum , nnn2 , shiftA , A ) ;

    const int fi = lat[lat[i].neighbor[mu]].neighbor[mu] ;
    const int bi = lat[lat[lat[i].back[mu]].back[mu]].back[mu] ;
    exact_log_slow_short( A , lat[fi].O[mu] ) ; 
    exact_log_slow_short( shiftA , lat[bi].O[mu] ) ; 
    a_plus_Sxbminc_short( sum , nnn3 , shiftA , A ) ;
  }
  return trace_deriv( sum ) ; 
}
