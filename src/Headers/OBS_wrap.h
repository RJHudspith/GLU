/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (OBS_wrap.h) is part of GLU.

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
   @file OBS_wrap.h
   @brief simple gauge observable calculation function prototypes
 */
#ifndef GLU_OBS_WRAP_H
#define GLU_OBS_WRAP_H

/**
   @fn void gf_check( const struct site *__restrict lat ) ;
   @brief Checks wheteher we are gauge fixed 
   @param lat :: Lattice fields
 */
void
gf_check( const struct site *__restrict lat ) ;

/**
   @fn void gauge( const struct site *__restrict lat )
   @brief compute gauge observables
   @param lat :: lattice gauge field
   computations of temporal and spatial plaquettes and links and polyakov loops
   and derivatives with different field definitions
 */
void 
gauge( const struct site *__restrict lat ) ;

#endif
