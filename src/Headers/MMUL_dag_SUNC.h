/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (MMUL_dag_SUNC.h) is part of GLU.

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
   @file MMUL_dag_SUNC.h
   @brief prototype functions for computing the matrix product \f$ a = b \times c^{\dagger} \f$
 */
#ifndef GLU_MMUL_DAG_SUNC_H
#define GLU_MMUL_DAG_SUNC_H

#ifndef GLU_BGQ // inlining matrix multiplies for the Q

   #if NC < 4
/**
   @fn void multab_dag_suNC( GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC ] , const GLU_complex c[ NCNC ] ) 
   @brief SU( Nc ) - tuned  matrix multiply \f$ a = b \times c^{\dagger}\f$

   computes \f$ a = b \times c^{\dagger} \f$  tuned for suNC by computing the signed minors for the bottom row and hand written complex-complex multiplications

   @warning only use-able for \f$ a,b,c \in SU(NC) \f$
 */
void 
multab_dag_suNC( GLU_complex a[ NCNC ] , 
		 const GLU_complex b[ NCNC ] , 
		 const GLU_complex c[ NCNC ] ) ;
  #else
    #define multab_dag_suNC multab_dag
  #endif

#else

// include the macros
#include "BGQ_mmul_dag.h"

#endif // end of the inline loop

#endif
