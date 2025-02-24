/*
Copyright 2013-2025 Renwick James Hudspith

    This file (init.h) is part of GLU.

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
   @file init.h
   @brief lattice initialisation prototype declarations
 */
#ifndef GLU_INIT_H
#define GLU_INIT_H

/**
   @fn void free_latt( void )
   @brief free the lattice struct information
 */
void
free_latt( void ) ;

/**
   @fn void init_latt( void )
   @brief initialise vital constants
   Inits the constants #LSQ, #LCU and #LVOLUME
 */
void
init_latt( void ) ;

/**
   @fn void init_navig( struct site *__restrict lat )
   @brief Function for generically initialising the lattice navigation
   
   @param lat :: Gauge field

   packs the struct lat's "neighbor" and
   "back" members which will be used for numerical
   derivatives and alike.
 **/
void 
init_navig( struct site *__restrict lat ) ;

#endif
