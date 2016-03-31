/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (MMULdagdag.h) is part of GLU.

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
   @file MMULdagdag.h
   @brief prototype functions for matrix multiplication of generic matrices and SU(NC) only
 */
#ifndef GLU_MMULDAGDAG_H
#define GLU_MMULDAGDAG_H

#ifndef GLU_BGQ

/**
   @fn void multab_dagdag( GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC ] , const GLU_complex c[ NCNC ] )
   @brief computes \f$ a = b^{\dagger} \times c^{\dagger} \f$

   loop unrolled, standard matrix multiply
 **/
void 
multab_dagdag( GLU_complex a[ NCNC ] , 
	       const GLU_complex b[ NCNC ] , 
	       const GLU_complex c[ NCNC ] ) ;

#else // GLU_BGQ ifdef

// inline macros
#include "BGQ_mmuldagdag.h"

#endif // end of the BGQ inlining loop ... 

#endif
