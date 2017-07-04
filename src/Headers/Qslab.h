/*
    Copyright 2013-2017 Renwick James Hudspith

    This file (Qslab.h) is part of GLU.

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
   @file Qslab.h
   @brief prototype definition of slab methond topological susceptibility
 */
#ifndef GLU_QSLAB_H
#define GLU_QSLAB_H

/**
   @fn int compute_Qsusc( struct site *lat , const struct cut_info CUTINFO , const size_t measurement )
   @brief compute the topological susceptibility from translationally-averaged slabs
   @param lat :: lattice links
   @param CUTINFO :: formatting hangover from old code
   @param measurement :: measurement index
   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int
compute_Qsusc( struct site *lat ,
	       const struct cut_info CUTINFO ,
	       const size_t measurement ) ;

#endif
