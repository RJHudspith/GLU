/*
    Copyright 2013 Renwick James Hudspith

    This file (New_HYP_4D.h) is part of GLU.

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
   @file New_HYP_4D.h
   @brief prototype funtion for a very slow recursive calculation of hypercubically blocked smearing
 */
#ifndef GLU_NEW_HYP_4D_H
#define GLU_NEW_HYP_4D_H

/**
   @fn void HYsmear4D( struct site *__restrict lat , const int smiters , const int type )
   @brief Very slow but memory cheap hypercubic-blocking smearing routine
   @param lat :: lattice gauge field
   @param smiters :: number of smearing iterations
   @param type :: enumerated, can be (SM_LOG,SM_APE,SM_STOUT)
 */
void 
HYsmear4D( struct site *__restrict lat , 
	   const int smiters , 
	   const int type ) ;

#endif
