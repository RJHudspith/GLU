/*
Copyright 2013-2025 Renwick James Hudspith

    This file (MMULdag_SUNC.h) is part of GLU.

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
   @file MMULdag.h
   @brief prototype functions for multiplication for \f$ a = b^{\dagger}\times c\f$ 
 */
#ifndef GLU_MMULDAG_SUNC_H
#define GLU_MMULDAG_SUNC_H

#ifndef GLU_BGQ // dirty, force inlining for the Q

#if NC < 4
/**
   @fn void multabdag_suNC( GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC ] , const GLU_complex c[ NCNC ] )
   @brief Computes general Nc x Nc matrix multiplication with \f$ a = b^{\dagger}\times c \f$

   sped up for suNC by computing the signed minors for the bottom row and hand-expanded complex-complex multiplications

   @warning only use-able for \f$ a,b,c \in SU(NC) \f$
 **/
void 
multabdag_suNC( GLU_complex a[ NCNC ] , 
		const GLU_complex b[ NCNC ] , 
		const GLU_complex c[ NCNC ] ) ;
#else
  #define multabdag_suNC multabdag
#endif

#else

// inline BGQmacros
#include "BGQ_mmuldag.h"

#endif // BGQ inline

#endif
