/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (relax.h) is part of GLU.

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
   @file relax.h
   @brief performs over-relaxation over the whole lattice
 */
#ifndef GLU_RELAX_H
#define GLU_RELAX_H

/**
   @fn int OR_lattice( struct site *lat , const struct draughtboard db )
   @brief overrelax our links
   @param lat :: lattice gauge fields
   @param db :: draughtboard indexing
   @return #GLU_SUCCESS or #GLU_FAILURE
   @warning overwrites lat
 */
int
OR_lattice( struct site *lat ,
	    const struct draughtboard db ) ;

#endif
