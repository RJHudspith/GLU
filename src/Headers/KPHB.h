/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (KPHB.h) is part of GLU.

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
   @file KPHB.h
   @brief prototype functions for quenched heat bath updates
 */
#ifndef KPHB_H
#define KPHB_H

/**
   @fn int hb_update( struct site *lat , const struct hb_info HBINFO , const char *traj_name , const GLU_output storage , const char *output_details )
   @brief Heatbath-Overrelaxation algorithm
   @param lat :: lattice gauge field
   @param HBINFO :: heatbath information
   @param traj_name :: trajectory name to write out
   @param storage :: gauge field storage type to write out
   @param output_details :: information for the configuration file
   @return #GLU_SUCCESS or #GLU_FAILURE
   @warning overwrites lat
 */
int
hb_update( struct site *lat ,
	   const struct hb_info HBINFO ,
	   const char *traj_name ,
	   const GLU_output storage , 
	   const char *output_details ) ;

#endif
