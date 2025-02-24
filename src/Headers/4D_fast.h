/*
Copyright 2013-2025 Renwick James Hudspith

    This file (4D_fast.h) is part of GLU.

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
   @file 4D_fast.h
   @brief protoype functions for the most memory-expensive and usually fastest Hypercubic-blocking code.
 */
#ifndef GLU_FOURD_FAST_H
#define GLU_FOURD_FAST_H

/**
   @fn void HYPSLsmear4D_expensive( struct site *__restrict lat , const size_t smiters , const int type ) ;
   @brief  This is the most memory-expensive Hypercubic-blocking code that we have.
   @param lat :: The lattice field.
   @param smiters :: The number of smearing iterations to perform.
   @param type :: Can either be SM_LOG, SM_APE or SM_STOUT

   @warning could not run if we dont have memory
   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int
HYPSLsmear4D_expensive( struct site *__restrict lat , 
			const size_t smiters , 
			const int type ) ;

#endif
