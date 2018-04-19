/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (GLU_bswap.h) is part of GLU.

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
   @file GLU_bswap.h
   @brief prototype functions for byte-swapping of different size objects
   @warning The byte swapping routines use OpenMP be warned
 */
#ifndef GLU_BSWAP_H
#define GLU_BSWAP_H

/**
   @fn void bswap_16( const int n , void *u ) 
   @brief swaps the bytes of a 16 bit array
   @param n :: length of the array
   @param u :: pointer to memory
 */
void
bswap_16( const int n , 
	  void *u ) ;

/**
   @fn void bswap_32( const int n , void *u ) 
   @brief swaps the bytes of a 32 bit array
   @param n :: length of the array
   @param u :: pointer to memory
 */
void
bswap_32( const int n , 
	  void *u ) ;

/**
   @fn void bswap_64( const int n , void *u ) 
   @brief swaps the bytes of a 64 bit array
   @param n :: length of the array
   @param u :: pointer to memory
 */
void
bswap_64( const int n , 
	  void *u ) ;

#endif
