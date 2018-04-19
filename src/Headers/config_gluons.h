/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (config_gluons.h) is part of GLU.

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
   @file config_gluons.h
   @brief configuration space gluon propagators
 */
#ifndef GLU_CONFIG_GLUONS_H
#define GLU_CONFIG_GLUONS_H

/**
   @fn int cuts_struct_configspace( struct site *__restrict A , const struct cut_info CUTINFO , const struct sm_info SMINFO )
   @brief computes the configuration-space gluon props
   @param A :: lattice fields
   @param CUTINFO :: cutting info
   @param SMINFO :: passed so that we can use smearing if we so wish

   @warning A is overwritten with the log fields
 */
int 
cuts_struct_configspace( struct site *__restrict A ,
			 const struct cut_info CUTINFO ,
			 const struct sm_info SMINFO ) ;

#endif
