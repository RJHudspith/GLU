/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (smearing_param.h) is part of GLU.

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
   @file smearing_param.h
   @brief computes an unsmeared and a smeared gluon propagator for the computation of the momentum space smearing parameter f(q)
   @ingroup Cuts
 */
#ifndef GLU_SMEARING_PARAM_H
#define GLU_SMEARING_PARAM_H

/**
   @fn int cuts_struct_smeared( struct site *__restrict A , const struct cut_info CUTINFO , const struct sm_info SMINFO )
   @brief computes an unsmeared and a smeared gluon propagator file
   @param A :: momentum space gluon fields
   @param CUTINFO :: cutting information struct
   @param SMINFO :: smearing information struct
 */
int 
cuts_struct_smeared( struct site *__restrict A ,
		     const struct cut_info CUTINFO ,
		     const struct sm_info SMINFO ) ;

/**
   @fn int cuts_struct_smeared( struct site *__restrict A , const struct cut_info CUTINFO , const struct sm_info SMINFO )
   @brief computes the Landau condition for each momentum and outputs to a prop file
   @param A :: momentum space gluon fields
   @param CUTINFO :: cutting information struct
   @param SMINFO :: smearing information struct
 */
int 
cuts_struct_smeared_Lcondition( struct site *__restrict A ,
				const struct cut_info CUTINFO ,
				const struct sm_info SMINFO ) ;

#endif
