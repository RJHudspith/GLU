/*
Copyright 2013-2025 Renwick James Hudspith

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
 */
#include "Mainfile.h"

#include "lin_derivs.h" // for trace deriv

// log definition of the gauge matrices at the point i
double
log_deriv( GLU_complex sum[ HERMSIZE ] , 
	   double *functional ,
	   const struct site *__restrict lat , 
	   const size_t i , 
	   const size_t MAX_DIR )
{
  GLU_complex A[ HERMSIZE ] , shiftA[ HERMSIZE ] ;
  GLU_real trAA ;
  *functional = 0.0 ;

  size_t mu ; 
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
