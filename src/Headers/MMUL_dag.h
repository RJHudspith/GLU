/*
    Copyright 2013 Renwick James Hudspith

    This file (MMUL_dag.h) is part of GLU.

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
   @file MMUL_dag.h
   @brief prototype functions for computing the matrix product \f$ a = b \times c^{\dagger} \f$
 */
#ifndef GLU_MMUL_DAG_H
#define GLU_MMUL_DAG_H

#ifndef GLU_BGQ // inlining matrix multiplies for the Q

/**
   @fn void multab_dag( GLU_complex a[ NCNC ] , const GLU_complex b[ NCNC ] , const GLU_complex c[ NCNC ] ) 
   @brief General  matrix multiply \f$ a = b \times c^{\dagger}\f$

   computes \f$ a = b \times c^{\dagger} \f$ general Nc x Nc matrix multiply does not require specific types of matrices uses the in-built complex multiply in complex.h
 **/
void 
multab_dag( GLU_complex a[ NCNC ] , 
	    const GLU_complex b[ NCNC ] , 
	    const GLU_complex c[ NCNC ] ) ;

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

/* Weird, code is slower using the macro inlines for my desktop, I put 
   them here because this might be strong architecture/compiler dependence */
  #if NC==3
    #define multab_dag( a , b , c )					\
    a[0] = b[0] * conj( c[0] ) + b[1] * conj( c[1] ) + b[2] * conj( c[2] ) ; \
    a[1] = b[0] * conj( c[3] ) + b[1] * conj( c[4] ) + b[2] * conj( c[5] ) ; \
    a[2] = b[0] * conj( c[6] ) + b[1] * conj( c[7] ) + b[2] * conj( c[8] ) ; \
    a[3] = b[3] * conj( c[0] ) + b[4] * conj( c[1] ) + b[5] * conj( c[2] ) ; \
    a[4] = b[3] * conj( c[3] ) + b[4] * conj( c[4] ) + b[5] * conj( c[5] ) ; \
    a[5] = b[3] * conj( c[6] ) + b[4] * conj( c[7] ) + b[5] * conj( c[8] ) ; \
    a[6] = b[6] * conj( c[0] ) + b[7] * conj( c[1] ) + b[8] * conj( c[2] ) ; \
    a[7] = b[6] * conj( c[3] ) + b[7] * conj( c[4] ) + b[8] * conj( c[5] ) ; \
    a[8] = b[6] * conj( c[6] ) + b[7] * conj( c[7] ) + b[8] * conj( c[8] ) ; 
  #elif NC==2
    #define multab_dag( a , b , c )				\
    a[0] = b[0] * conj( c[0] ) + b[1] * conj( c[1] ) ;	\
    a[1] = b[0] * conj( c[2] ) + b[1] * conj( c[3] ) ;	\
    a[2] = b[2] * conj( c[0] ) + b[3] * conj( c[1] ) ;	\
    a[3] = b[2] * conj( c[2] ) + b[3] * conj( c[3] ) ;
  #else // instead of inlining we have a function call
    void 
    multab_dag( GLU_complex a[ NCNC ] , 
		const GLU_complex b[ NCNC ] , 
		const GLU_complex c[ NCNC ] ) ;
  #endif

  #if NC==3 
    #define multab_dag_suNC( a , b , c )				\
    a[0] = b[0] * conj( c[0] ) + b[1] * conj( c[1] ) + b[2] * conj( c[2] ) ; \
    a[1] = b[0] * conj( c[3] ) + b[1] * conj( c[4] ) + b[2] * conj( c[5] ) ; \
    a[2] = b[0] * conj( c[6] ) + b[1] * conj( c[7] ) + b[2] * conj( c[8] ) ; \
    a[3] = b[3] * conj( c[0] ) + b[4] * conj( c[1] ) + b[5] * conj( c[2] ) ; \
    a[4] = b[3] * conj( c[3] ) + b[4] * conj( c[4] ) + b[5] * conj( c[5] ) ; \
    a[5] = b[3] * conj( c[6] ) + b[4] * conj( c[7] ) + b[5] * conj( c[8] ) ; \
    a[6] = conj( a[1] * a[5] - a[2] * a[4] ) ;				\
    a[7] = conj( a[2] * a[3] - a[0] * a[5] ) ;				\
    a[8] = conj( a[0] * a[4] - a[1] * a[3] ) ;
  #elif NC==2
    #define multab_dag_suNC( a , b , c )			\
    a[0] = b[0] * conj( c[0] ) + b[1] * conj( c[2] ) ;	\
    a[1] = b[0] * conj( c[2] ) + b[1] * conj( c[3] ) ;	\
    a[2] = -conj( a[1] ) ;				\
    a[3] = conj( a[0] ) ;
  #else
    #define multab_dag_suNC multab_dag
  #endif

#endif // end of the inline loop

#endif
