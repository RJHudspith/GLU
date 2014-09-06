/*
    Copyright 2013 Renwick James Hudspith

    This file (OrLandau.h) is part of GLU.

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
   @file OrLandau.h
   @brief protype functions for Over-relaxed gauge fixing routines
 */
#ifndef ORLANDAU_H
#define ORLANDAU_H

/**
   @fn int OrCoulomb( struct site *__restrict lat , double *theta , const int MAX_ITERS , const double ACC , const double OrParam )
   @brief Over-relaxation Coulomb gauge fixing code
   @param lat :: lattice gauge fields
   @param theta :: gauge fixing accuracy attained
   @param MAX_ITERS :: maximum number of iterations before complaint
   @param ACC :: accuracy intended to be fixed to
   @param OrParam :: Over-relaxation parameter
 */
int
OrCoulomb( struct site *__restrict lat ,
	   double *theta ,
	   const int MAX_ITERS , 
	   const double ACC ,
	   const double OrParam ) ;

/**
   @fn int OrLandau( struct site *__restrict lat , double *theta , const int MAX_ITERS , const double ACC , const double OrParam )
   @brief Over-relaxation Landau gauge fixing code
   @param lat :: lattice gauge fields
   @param theta :: gauge fixing accuracy attained
   @param MAX_ITERS :: maximum number of iterations before complaint
   @param ACC :: accuracy intended to be fixed to
   @param OrParam :: Over-relaxation parameter
 */
int
OrLandau( struct site *__restrict lat ,
	  double *theta ,
	  const int MAX_ITERS , 
	  const double ACC ,
	  const double OrParam ) ;

#endif