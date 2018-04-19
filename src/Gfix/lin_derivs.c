/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (lin_derivs.c) is part of GLU.

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
  @file lin_derivs.c
  @brief finite difference calculators ...

  I have included the Log and the Linear definitions
  as well as the stencil terms
 */
#include "Mainfile.h"

// computes the trace of dA
double
trace_deriv( GLU_complex *__restrict sum )
{
  // this takes the trace of AA, used in our GF stopping code
#if NC == 3
  return 4.0 * ( (double)(creal( sum[0] ) * creal( sum[0] )) +		\
		 (double)(creal( sum[3] ) * creal( sum[3] )) +		\
		 (double)(creal( sum[0] ) * creal( sum[3] )) +		\
		 (double)(creal( sum[1] ) * creal( sum[1] )) +		\
		 (double)(cimag( sum[1] ) * cimag( sum[1] )) +		\
		 (double)(creal( sum[2] ) * creal( sum[2] )) +		\
		 (double)(cimag( sum[2] ) * cimag( sum[2] )) +		\
		 (double)(creal( sum[4] ) * creal( sum[4] )) +		\
		 (double)(cimag( sum[4] ) * cimag( sum[4] ))) ; 
#elif NC == 2
  return 4.0 * ( (double)(creal( sum[0] ) * creal( sum[0] )) +		\
		 (double)(creal( sum[1] ) * creal( sum[1] )) +		\
		 (double)(cimag( sum[1] ) * cimag( sum[1] ))) ; 
#else
  GLU_real tr ;
  trace_ab_herm_short( &tr , sum , sum ) ;  
  return 2.0 * (double)tr ;
#endif
}

//lattice deriv is apparently ( AntiHermitian_proj method )
double
latt_deriv_AntiHermitian_proj( GLU_complex sum[ HERMSIZE ] , 
			       const struct site *__restrict lat , 
			       const size_t i , 
			       const size_t MAX_DIR )
{
  GLU_complex A[ HERMSIZE ] , shiftA[ HERMSIZE ] ;
  size_t mu ; 
  for( mu = 0 ; mu < MAX_DIR ; mu++ ) {
    Hermitian_proj_short( A , lat[i].O[mu] ) ; 
    Hermitian_proj_short( shiftA , lat[lat[i].back[mu]].O[mu] ) ; 
    // and accumulate into sum
    a_plus_Sxbminc_short( sum , 1.0 , shiftA , A ) ;
  }
  return trace_deriv( sum ) ;
}

// will return the local tr( || dA || ) here ! 
#define VAL 0.4444444444444444
double
fast_deriv_AntiHermitian_proj( GLU_complex sum[ HERMSIZE ] ,
			       const struct site *__restrict lat , 
			       const size_t i )
{
#if NC == 3 
  double REsum0 = 0. , REsum1 = 0. , IMsum1 = 0. , REsum2 = 0. , IMsum2 = 0. ;
  double REsum3 = 0. , REsum4 = 0. , IMsum4 = 0. ;
  double Ai0 , Ar1 , Ai1 , Ar2 , Ai2 , Ar3 , Ai3 , Ai4 , Ar5 , Ai5 ;
  double Ar6 , Ai6 , Ar7 , Ai7 , Ai8 ;
  double shAi0 , shAr1 , shAi1 , shAr2 , shAi2 , shAr3 , shAi3 , shAi4 , shAr5 , shAi5 ;
  double shAr6 , shAi6 , shAr7 , shAi7 , shAi8 ;
#elif NC == 2
  double REsum0 = 0. , REsum1 = 0. , IMsum1 = 0. ;
  double Ai0 , Ar1 , Ai1 ;
  double shAi0 , shAr1 , shAi1 ;
#else
  size_t nu ;
  for( nu = 0 ; nu < HERMSIZE ; nu++ ) {
    sum[ nu ] = 0. ;
  }
  return latt_deriv_AntiHermitian_proj( sum , lat , i , ND ) ;
#endif

  //calculate the top diagonal for both shift A and A
  size_t mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {

#if NC == 3 
    // cache some stuff in ....
    const GLU_real *qq = (const GLU_real*)lat[i].O[mu] ;
    Ai0 = *( qq + 1 )  ;
    Ar1 = *( qq + 2 )  ; Ai1 = *( qq + 3 ) ;
    Ar2 = *( qq + 4 )  ; Ai2 = *( qq + 5 ) ;
    Ar3 = *( qq + 6 )  ; Ai3 = *( qq + 7 ) ;
    Ai4 = *( qq + 9 )  ;  
    Ar5 = *( qq + 10 ) ; Ai5 = *( qq + 11 ) ;
    Ar6 = *( qq + 12 ) ; Ai6 = *( qq + 13 ) ;
    Ar7 = *( qq + 14 ) ; Ai7 = *( qq + 15 ) ;
    Ai8 = *( qq + 17 ) ;

    const GLU_real *shqq = (const GLU_real*)lat[ lat[i].back[mu] ].O[mu] ;
    shAi0 = *( shqq + 1 )  ;
    shAr1 = *( shqq + 2 )  ; shAi1 = *( shqq + 3 )  ;
    shAr2 = *( shqq + 4 )  ; shAi2 = *( shqq + 5 )  ;
    shAr3 = *( shqq + 6 )  ; shAi3 = *( shqq + 7 )  ;
    shAi4 = *( shqq + 9 )  ;
    shAr5 = *( shqq + 10 ) ; shAi5 = *( shqq + 11 ) ;
    shAr6 = *( shqq + 12 ) ; shAi6 = *( shqq + 13 ) ;
    shAr7 = *( shqq + 14 ) ; shAi7 = *( shqq + 15 ) ;
    shAi8 = *( shqq + 17 ) ;

    REsum0 += 2 * shAi0 - shAi4 - shAi8 ;
    REsum0 += -2 * Ai0 + Ai4 + Ai8 ;     

    REsum1 += shAr1 - shAr3 - Ar1 + Ar3 ;
    IMsum1 += shAi1 + shAi3 - Ai1 - Ai3 ;

    REsum2 += shAr2 - shAr6 - Ar2 + Ar6 ;
    IMsum2 += shAi2 + shAi6 - Ai2 - Ai6 ;

    REsum3 += 2 * shAi4 - shAi0 - shAi8 ;
    REsum3 += -2 * Ai4 + Ai0 + Ai8 ;

    REsum4 += shAr5 - shAr7 - Ar5 + Ar7 ;
    IMsum4 += shAi5 + shAi7 - Ai5 - Ai7 ;

#elif NC == 2

    const GLU_real *qq = (const GLU_real*)lat[i].O[mu] ;
    Ai0 = *( qq + 1 ) ;
    Ar1 = *( qq + 2 ) ; Ai1 = *( qq + 3 ) ;
 
    // and the backward one
    const GLU_real *shqq = (const GLU_real*)lat[ lat[i].back[mu] ].O[mu] ;
    shAi0 = *( shqq + 1 ) ;
    shAr1 = *( shqq + 2 ) ; shAi1 = *( shqq + 3 ) ;

    #ifdef deriv_MLG

    register const GLU_real tr = 1. / ( 1.0 + Ar0 ) ; 
    register const GLU_real tr_back = 1. / ( 1.0 + shAr0 ) ; 
      
    REsum0 += tr_back * shAi0 - tr * Ai0 ;  
    REsum1 += tr_back * shAr1 - tr * Ar1 ; 
    IMsum1 += tr_back * shAi1 - tr * Ai1 ;
      
    #else

    REsum0 += shAi0 - Ai0 ;  
    REsum1 += shAr1 - Ar1 ; 
    IMsum1 += shAi1 - Ai1 ;

    #endif

#endif
    }

#if NC == 3

  *( sum + 0 ) = REsum0 / 3. ;  
  *( sum + 1 ) = -( REsum1 )*0.5*I + 0.5 * IMsum1 ;
  *( sum + 2 ) = -( REsum2 )*0.5*I + 0.5 * IMsum2 ;
  *( sum + 3 ) = REsum3 / 3. ;  
  *( sum + 4 ) = -( REsum4 )*0.5*I + 0.5 * IMsum4 ;

  return VAL * ( REsum0 * REsum0 + REsum3 * ( REsum3 + REsum0 ) ) +	\
    ( REsum1 * REsum1 + IMsum1 * IMsum1 +				\
      REsum2 * REsum2 + IMsum2 * IMsum2 +				\
      REsum4 * REsum4 + IMsum4 * IMsum4 ) ;
  
#elif NC == 2

  *( sum + 0 ) = REsum0 ;
  *( sum + 1 ) = -I * REsum1 + IMsum1 ;

  return ( REsum0 * REsum0 + REsum1 * REsum1 + IMsum1 * IMsum1 ) ;

#endif
}
#undef VAL
