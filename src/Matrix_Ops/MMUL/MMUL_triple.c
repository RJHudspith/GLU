/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (MMUL_SSE.c) is part of GLU.

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
   @file MMUL_triple.c
   @brief Common triple product A.B.C^\dagger and A^\dagger.B.C
 */
#include "Mainfile.h"

#if ( defined HAVE_IMMINTRIN_H ) && !( defined SINGLE_PREC )  

#include <immintrin.h>
#include "SSE2_OPS.h"

// simple matrix multiplication a = b.c.d^\dag
void 
multabcdag_suNC( GLU_complex a[ NCNC ] , 
		 const GLU_complex b[ NCNC ] , 
		 const GLU_complex c[ NCNC ] ,
		 const GLU_complex d[ NCNC ] )
{
#if NC==3
  __m128d *A = ( __m128d* )a ;
  const __m128d *B = ( __m128d* )b , *C = ( __m128d* )c , *D = ( __m128d* )d ;
  register const __m128d CD0 = _mm_add_pd( SSE2_MUL_CONJ( *C , *(D+0) ) , _mm_add_pd( SSE2_MUL_CONJ( *(C+1) , *(D+1) ) , SSE2_MUL_CONJ( *(C+2) , *(D+2) ) ) ) ;
  register const __m128d CD1 = _mm_add_pd( SSE2_MUL_CONJ( *C , *(D+3) ) , _mm_add_pd( SSE2_MUL_CONJ( *(C+1) , *(D+4) ) , SSE2_MUL_CONJ( *(C+2) , *(D+5) ) ) ) ;
  register const __m128d CD2 = _mm_add_pd( SSE2_MUL_CONJ( *C , *(D+6) ) , _mm_add_pd( SSE2_MUL_CONJ( *(C+1) , *(D+7) ) , SSE2_MUL_CONJ( *(C+2) , *(D+8) ) ) ) ;
  C += 3 ;
  register const __m128d CD3 = _mm_add_pd( SSE2_MUL_CONJ( *C , *(D+0) ) , _mm_add_pd( SSE2_MUL_CONJ( *(C+1) , *(D+1) ) , SSE2_MUL_CONJ( *(C+2) , *(D+2) ) ) ) ;
  register const __m128d CD4 = _mm_add_pd( SSE2_MUL_CONJ( *C , *(D+3) ) , _mm_add_pd( SSE2_MUL_CONJ( *(C+1) , *(D+4) ) , SSE2_MUL_CONJ( *(C+2) , *(D+5) ) ) ) ;
  register const __m128d CD5 = _mm_add_pd( SSE2_MUL_CONJ( *C , *(D+6) ) , _mm_add_pd( SSE2_MUL_CONJ( *(C+1) , *(D+7) ) , SSE2_MUL_CONJ( *(C+2) , *(D+8) ) ) ) ;
  // completion od CD submatrix
  register const __m128d CD6 = SSE2_CONJ( _mm_sub_pd( SSE2_MUL( CD1 , CD5 ) , SSE2_MUL( CD2 , CD4 ) ) ) ;
  register const __m128d CD7 = SSE2_CONJ( _mm_sub_pd( SSE2_MUL( CD2 , CD3 ) , SSE2_MUL( CD0 , CD5 ) ) ) ;
  register const __m128d CD8 = SSE2_CONJ( _mm_sub_pd( SSE2_MUL( CD0 , CD4 ) , SSE2_MUL( CD1 , CD3 ) ) ) ;
  A[0] = _mm_add_pd( SSE2_MUL( *B , CD0 ) , _mm_add_pd( SSE2_MUL( *(B+1) , CD3 ) , SSE2_MUL( *(B+2) , CD6 ) ) ) ;
  A[1] = _mm_add_pd( SSE2_MUL( *B , CD1 ) , _mm_add_pd( SSE2_MUL( *(B+1) , CD4 ) , SSE2_MUL( *(B+2) , CD7 ) ) ) ;
  A[2] = _mm_add_pd( SSE2_MUL( *B , CD2 ) , _mm_add_pd( SSE2_MUL( *(B+1) , CD5 ) , SSE2_MUL( *(B+2) , CD8 ) ) ) ;
  B += 3 ;
  A[3] = _mm_add_pd( SSE2_MUL( *B , CD0 ) , _mm_add_pd( SSE2_MUL( *(B+1) , CD3 ) , SSE2_MUL( *(B+2) , CD6 ) ) ) ;
  A[4] = _mm_add_pd( SSE2_MUL( *B , CD1 ) , _mm_add_pd( SSE2_MUL( *(B+1) , CD4 ) , SSE2_MUL( *(B+2) , CD7 ) ) ) ;
  A[5] = _mm_add_pd( SSE2_MUL( *B , CD2 ) , _mm_add_pd( SSE2_MUL( *(B+1) , CD5 ) , SSE2_MUL( *(B+2) , CD8 ) ) ) ;
  // completion of A
  *( A + 6 ) = SSE2_CONJ( _mm_sub_pd( SSE2_MUL( *( A + 1 ) , *( A + 5 ) ) , SSE2_MUL( *( A + 2 ) , *( A + 4 ) ) ) ) ;
  *( A + 7 ) = SSE2_CONJ( _mm_sub_pd( SSE2_MUL( *( A + 2 ) , *( A + 3 ) ) , SSE2_MUL( *( A + 0 ) , *( A + 5 ) ) ) ) ;
  *( A + 8 ) = SSE2_CONJ( _mm_sub_pd( SSE2_MUL( *( A + 0 ) , *( A + 4 ) ) , SSE2_MUL( *( A + 1 ) , *( A + 3 ) ) ) ) ;
#elif (NC==2)
  __m128d *A = ( __m128d* )a ;
  const __m128d *B = ( __m128d* )b , *C = ( __m128d* )c , *D = ( __m128d* )d ;
  // CD multiply
  register const __m128d CD0 = _mm_add_pd( SSE2_MUL_CONJ( *( C + 0 ) , *( D + 0 ) ) , SSE2_MUL_CONJ( *( C + 1 ) , *( D + 1 ) ) ) ;
  register const __m128d CD1 = _mm_add_pd( SSE2_MUL_CONJ( *( C + 0 ) , *( D + 2 ) ) , SSE2_MUL_CONJ( *( C + 1 ) , *( D + 3 ) ) ) ;
  register const __m128d CD2 = SSE_FLIP( SSE2_CONJ( CD1 ) ) ;
  register const __m128d CD3 = SSE2_CONJ( CD0 ) ;
  *( A + 0 ) = _mm_add_pd( SSE2_MUL( *( B + 0 ) , CD0 ) , SSE2_MUL( *( B + 1 ) , CD2 ) ) ;
  *( A + 1 ) = _mm_add_pd( SSE2_MUL( *( B + 0 ) , CD1 ) , SSE2_MUL( *( B + 1 ) , CD3 ) ) ;
  *( A + 2 ) = SSE_FLIP( SSE2_CONJ( *( A + 1 ) ) ) ; 
  *( A + 3 ) = SSE2_CONJ( *( A + 0 ) ) ;
#else
  GLU_complex tmp[NCNC] GLUalign ;
  multab_dag_suNC( tmp , c, d ) ;
  multab_suNC( a , b, tmp ) ;
#endif
}

void 
multadagbc_suNC( GLU_complex a[ NCNC ] , 
		 const GLU_complex b[ NCNC ] , 
		 const GLU_complex c[ NCNC ] ,
		 const GLU_complex d[ NCNC ] )
{
#if NC==3
  // a = b^dag.c.d
   __m128d *A = ( __m128d* )a ;
  const __m128d *B = ( __m128d* )b ;
  const __m128d *C = ( __m128d* )c ;
  const __m128d *D = ( __m128d* )d ;
  register __m128d CD0 = _mm_add_pd( SSE2_MUL( *C , *(D+0) ) , _mm_add_pd( SSE2_MUL( *(C+1) , *(D+3) ) , SSE2_MUL( *(C+2) , *(D+6) ) ) ) ;
  register __m128d CD3 = _mm_add_pd( SSE2_MUL( *(C+3) , *(D+0) ) , _mm_add_pd( SSE2_MUL( *(C+4) , *(D+3) ) , SSE2_MUL( *(C+5) , *(D+6) ) ) ) ;
  register __m128d CD1 = _mm_add_pd( SSE2_MUL( *C , *(D+1) ) , _mm_add_pd( SSE2_MUL( *(C+1) , *(D+4) ) , SSE2_MUL( *(C+2) , *(D+7) ) ) ) ;
  register __m128d CD4 = _mm_add_pd( SSE2_MUL( *(C+3) , *(D+1) ) , _mm_add_pd( SSE2_MUL( *(C+4) , *(D+4) ) , SSE2_MUL( *(C+5) , *(D+7) ) ) ) ;
  register __m128d CD2 = _mm_add_pd( SSE2_MUL( *C , *(D+2) ) , _mm_add_pd( SSE2_MUL( *(C+1) , *(D+5) ) , SSE2_MUL( *(C+2) , *(D+8) ) ) ) ;
  register __m128d CD5 = _mm_add_pd( SSE2_MUL( *(C+3) , *(D+2) ) , _mm_add_pd( SSE2_MUL( *(C+4) , *(D+5) ) , SSE2_MUL( *(C+5) , *(D+8) ) ) ) ;
  // completion od CD submatrix
  register __m128d CD6 = SSE2_CONJ( _mm_sub_pd( SSE2_MUL( CD1 , CD5 ) , SSE2_MUL( CD2 , CD4 ) ) ) ;
  register __m128d CD7 = SSE2_CONJ( _mm_sub_pd( SSE2_MUL( CD2 , CD3 ) , SSE2_MUL( CD0 , CD5 ) ) ) ;
  register __m128d CD8 = SSE2_CONJ( _mm_sub_pd( SSE2_MUL( CD0 , CD4 ) , SSE2_MUL( CD1 , CD3 ) ) ) ;
  A[0] = _mm_add_pd( SSE2_MULCONJ( *B , CD0 ) , _mm_add_pd( SSE2_MULCONJ( *(B+3) , CD3 ) , SSE2_MULCONJ( *(B+6) , CD6 ) ) ) ;
  A[1] = _mm_add_pd( SSE2_MULCONJ( *B , CD1 ) , _mm_add_pd( SSE2_MULCONJ( *(B+3) , CD4 ) , SSE2_MULCONJ( *(B+6) , CD7 ) ) ) ;
  A[2] = _mm_add_pd( SSE2_MULCONJ( *B , CD2 ) , _mm_add_pd( SSE2_MULCONJ( *(B+3) , CD5 ) , SSE2_MULCONJ( *(B+6) , CD8 ) ) ) ;
  A[3] = _mm_add_pd( SSE2_MULCONJ( *(B+1) , CD0 ) , _mm_add_pd( SSE2_MULCONJ( *(B+4) , CD3 ) , SSE2_MULCONJ( *(B+7) , CD6 ) ) ) ;
  A[4] = _mm_add_pd( SSE2_MULCONJ( *(B+1) , CD1 ) , _mm_add_pd( SSE2_MULCONJ( *(B+4) , CD4 ) , SSE2_MULCONJ( *(B+7) , CD7 ) ) ) ;
  A[5] = _mm_add_pd( SSE2_MULCONJ( *(B+1) , CD2 ) , _mm_add_pd( SSE2_MULCONJ( *(B+4) , CD5 ) , SSE2_MULCONJ( *(B+7) , CD8 ) ) ) ;
  // completion of A
  *( A + 6 ) = SSE2_CONJ( _mm_sub_pd( SSE2_MUL( *( A + 1 ) , *( A + 5 ) ) , SSE2_MUL( *( A + 2 ) , *( A + 4 ) ) ) ) ;
  *( A + 7 ) = SSE2_CONJ( _mm_sub_pd( SSE2_MUL( *( A + 2 ) , *( A + 3 ) ) , SSE2_MUL( *( A + 0 ) , *( A + 5 ) ) ) ) ;
  *( A + 8 ) = SSE2_CONJ( _mm_sub_pd( SSE2_MUL( *( A + 0 ) , *( A + 4 ) ) , SSE2_MUL( *( A + 1 ) , *( A + 3 ) ) ) ) ;
#elif (NC==2)
  __m128d *A = ( __m128d* )a ;
  const __m128d *B = ( __m128d* )b , *C = ( __m128d* )c , *D = ( __m128d* )d ;
  // CD multiply
  register const __m128d CD0 = _mm_add_pd( SSE2_MUL( *( C + 0 ) , *( D + 0 ) ) , SSE2_MUL( *( C + 1 ) , *( D + 2 ) ) ) ;
  register const __m128d CD1 = _mm_add_pd( SSE2_MUL( *( C + 0 ) , *( D + 1 ) ) , SSE2_MUL( *( C + 1 ) , *( D + 3 ) ) ) ;
  register const __m128d CD2 = SSE_FLIP( SSE2_CONJ( CD1 ) ) ;
  register const __m128d CD3 = SSE2_CONJ( CD0 ) ;
  *( A + 0 ) = _mm_add_pd( SSE2_MULCONJ( *( B + 0 ) , CD0 ) , SSE2_MULCONJ( *( B + 2 ) , CD2 ) ) ;
  *( A + 1 ) = _mm_add_pd( SSE2_MULCONJ( *( B + 0 ) , CD1 ) , SSE2_MULCONJ( *( B + 2 ) , CD3 ) ) ;
  *( A + 2 ) = SSE_FLIP( SSE2_CONJ( *( A + 1 ) ) ) ; 
  *( A + 3 ) = SSE2_CONJ( *( A + 0 ) ) ;
#else
  GLU_complex tmp[NCNC] GLUalign ;
  multab_suNC( tmp , c, d ) ;
  multabdag_suNC( a , b, tmp ) ;
#endif
}

#else

void 
multabcdag_suNC( GLU_complex a[ NCNC ] , 
		 const GLU_complex b[ NCNC ] , 
		 const GLU_complex c[ NCNC ] ,
		 const GLU_complex d[ NCNC ] )
{
  GLU_complex tmp[NCNC] GLUalign ;
  multab_dag_suNC( tmp , c , d ) ;
  multab_suNC( a , b , tmp ) ;
}

void 
multadagbc_suNC( GLU_complex a[ NCNC ] , 
		 const GLU_complex b[ NCNC ] , 
		 const GLU_complex c[ NCNC ] ,
		 const GLU_complex d[ NCNC ] )
{
  GLU_complex tmp[NCNC] GLUalign ;
  multab_suNC( tmp , c , d ) ;
  multabdag_suNC( a , b , tmp ) ;
}

#endif
