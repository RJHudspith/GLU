/*
    Copyright 2013 Renwick James Hudspith

    This file (SUNCxU1_config.h) is part of GLU.

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
   @file SUNCxU1_config.h
   @brief prototye declarations for the U(1)-ing of gauge links
 */

#ifndef GLU_SUNCxU1_CONFIG_H
#define GLU_SUNCxU1_CONFIG_H

/**
   @fn void suNC_cross_u1( struct site *__restrict lat , const struct u1_info U1INFO ) 
   @brief creates links which are multiplied by a non-compact U(1) field

   @warning rewrites the gauge field
   @warning cannot write out configurations in any format but NCxNC #config_size

   calls compute_U1_obs()
   
 */
void 
suNC_cross_u1( struct site *__restrict lat , 
	       const struct u1_info U1INFO ) ;

#endif
