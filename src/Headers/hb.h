/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (hb.h) is part of GLU.

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
   @file hb.h
   @brief prototype declarations for heatbath algorithm
 */
#ifndef GLU_HB_H
#define GLU_HB_H

/**
   @fn int hb_lattice( struct site *lat , const double invbeta , const struct draughtboard db )
   @brief perform heat-bath updates
   @brief lat :: lattice gauge fields
   @brief invbeta :: the inverse of the coupling beta
   @brief db :: draughtboard indexing
 */
int
hb_lattice( struct site *lat ,
	    const double invbeta ,
	    const struct draughtboard db ) ;

#endif
