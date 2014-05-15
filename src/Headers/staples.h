/*
    Copyright 2013 Renwick James Hudspith

    This file (staples.h) is part of GLU.

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
   @file staples.h
   @brief staple calculations get their own header
 */

#ifndef GLU_STAPLES_H
#define GLU_STAPLES_H

/**
   @fn void all_staples( GLU_complex stap[ NCNC ] , const struct site *__restrict lat , const int i , const int mu , const int dir , const int type )
   @brief the computation of the basic unimproved staple.
   @param stap :: sum of the contributing mu-nu staples
   @param lat :: lattice gauge field
   @param i :: site index
   @param mu :: direction of the link we are improving
   @param dir :: number of directions we smear in
   @param type :: smearing projection type
 **/
void
all_staples( GLU_complex stap[ NCNC ] , 
	     const struct site *__restrict lat , 
	     const int i , 
	     const int mu , 
	     const int dir , 
	     const int type ) ;

/**
   @fn void all_staples_improve( GLU_complex stap[ NCNC ] , const struct site *__restrict lat , const int i , const int mu , const int dir , const int type )
   @brief he computation of the (over)improved staple(s)
   @param stap :: sum of the contributing mu-nu staples
   @param lat :: lattice gauge field
   @param i :: site index
   @param mu :: direction of the link we are improving
   @param dir :: number of directions we smear in
   @param type :: smearing projection type

   for the improvement strategies <br>
      > SYMANZIK          = Tree level Symanzik improvement <br>
      > IWASAKI/DBW2      = Change in the coefficients c0 and c1 <br>
      > SYMANZIK_ONE_LOOP = Tadpole improved action <br>
   
   @warning requires the #IMPROVED_SMEARING definition to be compiled in
 **/
void
all_staples_improve( GLU_complex stap[ NCNC ] , 
		     const struct site *__restrict lat , 
		     const int i , 
		     const int mu , 
		     const int dir , 
		     const int type ) ;

#endif
