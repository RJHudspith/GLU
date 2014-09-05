/*
    Copyright 2013 Renwick James Hudspith

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
			       const int i , 
			       const int MAX_DIR )
{
  GLU_complex A[ HERMSIZE ] , shiftA[ HERMSIZE ] ;
  int mu ; 
  for( mu = 0 ; mu < MAX_DIR ; mu++ ) {
    Hermitian_proj_short( A , lat[i].O[mu] ) ; 
    Hermitian_proj_short( shiftA , lat[lat[i].back[mu]].O[mu] ) ; 
    // and accumulate into sum
    a_plus_Sxbminc_short( sum , 1.0 , shiftA , A ) ;
  }
  return trace_deriv( sum ) ;
}

//outer and inner lattice derivatives ( AntiHermitian_proj )
double
latt_derivnn_AntiHermitian_proj( GLU_complex sum[ HERMSIZE ] , 
				 const struct site *__restrict lat , 
				 const int i , 
				 const int MAX_DIR ) 
{
  GLU_complex shiftA[ HERMSIZE ] , A[ HERMSIZE ] ;
  int mu ; 
  for( mu = 0 ; mu < MAX_DIR ; mu++ ) {
      // first deriv
    Hermitian_proj_short( A , lat[i].O[mu] ) ; 
    Hermitian_proj_short( shiftA , lat[lat[i].back[mu]].O[mu] ) ; 
    // multiply deriv by a constant and accumulate into sum
    a_plus_Sxbminc_short( sum , nn1 , shiftA , A ) ;
      
    Hermitian_proj_short( A , lat[lat[i].neighbor[mu]].O[mu] ) ; 
    Hermitian_proj_short( shiftA , lat[lat[lat[i].back[mu]].back[mu]].O[mu] ) ; 
    a_plus_Sxbminc_short( sum , nn2 , shiftA , A ) ;
  }
  return trace_deriv( sum ) ;
}

// will return the local tr( || dA || ) here ! 
#define VAL 0.4444444444444444
double
fast_deriv_AntiHermitian_proj( GLU_complex sum[ HERMSIZE ] ,
			       double *functional ,
			       const struct site *__restrict lat , 
			       const int i )
{
  *functional = 0.0 ;
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
  int nu ;
  for( nu = 0 ; nu < HERMSIZE ; nu++ ) {
    sum[ nu ] = 0. ;
  }
  for( nu = 0 ; nu < ND ; nu++ ) {
    *functional += creal( trace( lat[i].O[nu] ) ) ;
  }
  return latt_deriv_AntiHermitian_proj( sum , lat , i , ND ) ;
#endif

  // local functional sum accumulator
  register double loc_sum = 0.0 ;

  //calculate the top diagonal for both shift A and A
  int mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    // compute the functional
    #if NC == 3
    loc_sum += creal( lat[i].O[mu][0] ) ;
    loc_sum += creal( lat[i].O[mu][4] ) ;
    loc_sum += creal( lat[i].O[mu][8] ) ;
    #elif NC == 2
    loc_sum += creal( lat[i].O[mu][0] ) ;
    loc_sum += creal( lat[i].O[mu][3] ) ;
    #else
    loc_sum += creal( trace( lat[i].O[mu] ) ) ;
    #endif

#if NC == 3 
    // cache some stuff in ....
    GLU_real *qq = ( GLU_real* )lat[i].O[mu] ;
    Ai0 = *( qq + 1 )  ;
    Ar1 = *( qq + 2 )  ; Ai1 = *( qq + 3 ) ;
    Ar2 = *( qq + 4 )  ; Ai2 = *( qq + 5 ) ;
    Ar3 = *( qq + 6 )  ; Ai3 = *( qq + 7 ) ;
    Ai4 = *( qq + 9 )  ;  
    Ar5 = *( qq + 10 ) ; Ai5 = *( qq + 11 ) ;
    Ar6 = *( qq + 12 ) ; Ai6 = *( qq + 13 ) ;
    Ar7 = *( qq + 14 ) ; Ai7 = *( qq + 15 ) ;
    Ai8 = *( qq + 17 ) ;

    GLU_real *shqq = ( GLU_real* )lat[ lat[i].back[mu] ].O[mu] ;
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

    GLU_real *qq = ( GLU_real* )lat[i].O[mu] ;
    Ai0 = *( qq + 1 ) ;
    Ar1 = *( qq + 2 ) ; Ai1 = *( qq + 3 ) ;
 
    // and the backward one
    GLU_real *shqq = ( GLU_real* )lat[ lat[i].back[mu] ].O[mu] ;
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
  
  // set the functional
  *functional = loc_sum ;

#if NC == 3

  *( sum + 0 ) = OneO3 * REsum0 ;  
  *( sum + 1 ) = OneOI2 * ( REsum1 ) + 0.5 * IMsum1 ;
  *( sum + 2 ) = OneOI2 * ( REsum2 ) + 0.5 * IMsum2 ;
  *( sum + 3 ) = OneO3 * REsum3 ;  
  *( sum + 4 ) = OneOI2 * ( REsum4 ) + 0.5 * IMsum4 ;

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

// Fast nearest neighbour derivative
double
fast_derivnn_AntiHermitian_proj( GLU_complex sum[ HERMSIZE ] , 
				 const struct site *__restrict lat , 
				 const int i ) 
{
#if NC == 3 
  double REsum0 = 0. , REsum1 = 0. , IMsum1 = 0. , REsum2 = 0. , IMsum2 = 0. ;
  double REsum3 = 0. , REsum4 = 0. , IMsum4 = 0. ;
  double REsum02 = 0. , REsum12 = 0. , IMsum12 = 0. , REsum22 = 0. , IMsum22 = 0. ;
  double REsum32 = 0. , REsum42 = 0. , IMsum42 = 0. ;
  double Ai0 , Ar1 , Ai1 , Ar2 , Ai2 , Ar3 , Ai3 ; 
  double Ai4 , Ar5 , Ai5 , Ar6 , Ai6 , Ar7 , Ai7 , Ai8 ; 
  // and the shifts sh
  double shAi0 , shAr1 , shAi1 , shAr2 , shAi2 , shAr3 , shAi3 ; 
  double shAi4 , shAr5 , shAi5 , shAr6 , shAi6 , shAr7 , shAi7 , shAi8 ; 
#elif NC == 2
  double REsum0 = 0. , REsum1 = 0. , IMsum1 = 0. ;
  double REsum02 = 0. , REsum12 = 0. , IMsum12 = 0. ;
  double Ai0 , Ar1 , Ai1 ;
  double shAi0 , shAr1 , shAi1 ;
#else
  int nu ;
  for( nu = 0 ; nu < HERMSIZE ; nu++ ) { sum[ nu ] = 0. ; }
  return latt_derivnn_AntiHermitian_proj( sum , lat , i , ND ) ;
#endif

  int mu ; 
  for( mu = 0 ; mu < ND ; mu++ ) {
#if NC == 3
    // cache some stuff in ....
    const GLU_real *qq = ( GLU_real* )lat[i].O[mu] ;
    //register GLU_real t = *qq++ ; //increment useless pointer
    Ai0 = *( qq + 1 )  ;
    Ar1 = *( qq + 2 )  ; Ai1 = *( qq + 3 ) ;
    Ar2 = *( qq + 4 )  ; Ai2 = *( qq + 5 ) ;
    Ar3 = *( qq + 6 )  ; Ai3 = *( qq + 7 ) ;
    Ai4 = *( qq + 9 )  ;
    Ar5 = *( qq + 10 ) ; Ai5 = *( qq + 11 ) ;
    Ar6 = *( qq + 12 ) ; Ai6 = *( qq + 13 ) ;
    Ar7 = *( qq + 14 ) ; Ai7 = *( qq + 15 ) ;
    Ai8 = *( qq + 17 ) ;

    register const int it = lat[i].back[mu] ;
    const GLU_real *shqq = ( GLU_real* )lat[ it ].O[mu] ;
    shAi0 = *( shqq + 1 )  ;
    shAr1 = *( shqq + 2 )  ; shAi1 = *( shqq + 3 ) ;
    shAr2 = *( shqq + 4 )  ; shAi2 = *( shqq + 5 ) ;
    shAr3 = *( shqq + 6 )  ; shAi3 = *( shqq + 7 ) ;
    shAi4 = *( shqq + 9 )  ;
    shAr5 = *( shqq + 10 ) ; shAi5 = *( shqq + 11 ) ;
    shAr6 = *( shqq + 12 ) ; shAi6 = *( shqq + 13 ) ;
    shAr7 = *( shqq + 14 ) ; shAi7 = *( shqq + 15 ) ;
    shAi8 = *( shqq + 17 ) ;

    REsum0 += 2.0 * shAi0 - shAi4 - shAi8 ;
    REsum0 += -2.0 * Ai0 + Ai4 + Ai8 ;     

    REsum1 += shAr1 - shAr3 - Ar1 + Ar3 ;
    IMsum1 += shAi1 + shAi3 - Ai1 - Ai3 ;

    REsum2 += shAr2 - shAr6 - Ar2 + Ar6 ;
    IMsum2 += shAi2 + shAi6 - Ai2 - Ai6 ;

    REsum3 += 2.0 * shAi4 - shAi0 - shAi8 ;
    REsum3 += -2.0 * Ai4 + Ai0 + Ai8 ;

    REsum4 += shAr5 - shAr7 - Ar5 + Ar7 ;
    IMsum4 += shAi5 + shAi7 - Ai5 - Ai7 ;

    // second deriv
    // cache some stuff in ....
    GLU_real *qq2 = ( GLU_real* )lat[ lat[i].neighbor[mu] ].O[mu] ;
    Ai0 = *( qq2 + 1 )  ;
    Ar1 = *( qq2 + 2 )  ; Ai1 = *( qq2 + 3 ) ;
    Ar2 = *( qq2 + 4 )  ; Ai2 = *( qq2 + 5 ) ;
    Ar3 = *( qq2 + 6 )  ; Ai3 = *( qq2 + 7 ) ;
    Ai4 = *( qq2 + 9 )  ;
    Ar5 = *( qq2 + 10 ) ; Ai5 = *( qq2 + 11 ) ;
    Ar6 = *( qq2 + 12 ) ; Ai6 = *( qq2 + 13 ) ;
    Ar7 = *( qq2 + 14 ) ; Ai7 = *( qq2 + 15 ) ;
    Ai8 = *( qq2 + 17 ) ;

    register const int it2b = lat[it].back[mu] ;
    GLU_real *shqq2 = ( GLU_real* )lat[it2b].O[mu] ;
    shAi0 = *( shqq2 + 1 )  ;
    shAr1 = *( shqq2 + 2 )  ; shAi1 = *( shqq2 + 3 ) ;
    shAr2 = *( shqq2 + 4 )  ; shAi2 = *( shqq2 + 5 ) ;
    shAr3 = *( shqq2 + 6 )  ; shAi3 = *( shqq2 + 7 ) ;
    shAi4 = *( shqq2 + 9 )  ;
    shAr5 = *( shqq2 + 10 ) ; shAi5 = *( shqq2 + 11 ) ;
    shAr6 = *( shqq2 + 12 ) ; shAi6 = *( shqq2 + 13 ) ;
    shAr7 = *( shqq2 + 14 ) ; shAi7 = *( shqq2 + 15 ) ;
    shAi8 = *( shqq2 + 17 ) ;

    REsum02 += 2.0 * shAi0 - shAi4 - shAi8 ;
    REsum02 += -2.0 * Ai0 + Ai4 + Ai8 ;     

    REsum12 += shAr1 - shAr3 - Ar1 + Ar3 ;
    IMsum12 += shAi1 + shAi3 - Ai1 - Ai3 ;

    REsum22 += shAr2 - shAr6 - Ar2 + Ar6 ;
    IMsum22 += shAi2 + shAi6 - Ai2 - Ai6 ;

    REsum32 += 2.0 * shAi4 - shAi0 - shAi8 ;
    REsum32 += -2.0 * Ai4 + Ai0 + Ai8 ;

    REsum42 += shAr5 - shAr7 - Ar5 + Ar7 ;
    IMsum42 += shAi5 + shAi7 - Ai5 - Ai7 ;

#elif NC == 2

    const GLU_real *qq = ( GLU_real* )lat[i].O[mu] ;
    Ai0 = *( qq + 1 ) ;
    Ar1 = *( qq + 2 ) ; Ai1 = *( qq + 3 ) ;
 
    // and the backward one
    const int it = lat[i].back[mu] ;
    const GLU_real *shqq = ( GLU_real* )lat[it].O[mu] ;
    shAi0 = *( shqq + 1 ) ;
    shAr1 = *( shqq + 2 ) ; shAi1 = *( shqq + 3 ) ;

    REsum0 += shAi0 - Ai0 ;  
    REsum1 += shAr1 - Ar1 ; 
    IMsum1 += shAi1 - Ai1 ;

    register const int it2f = lat[i].neighbor[mu] ;
    const GLU_real *qq2 = ( GLU_real* )lat[it2f].O[mu] ;
    Ai0 = *( qq2 + 1 ) ;
    Ar1 = *( qq2 + 2 ) ; Ai1 = *( qq2 + 3 ) ;

    register const int it2b = lat[it].back[mu] ;
    const GLU_real *shqq2 = ( GLU_real* )lat[it2b].O[mu] ;
    shAi0 = *( shqq2 + 1 ) ;
    shAr1 = *( shqq2 + 2 ) ; shAi1 = *( shqq2 + 3 ) ;

    REsum02 += shAi0 - Ai0 ;  
    REsum12 += shAr1 - Ar1 ; 
    IMsum12 += shAi1 - Ai1 ;

#endif 
    }
  // compute the sums
#if NC == 3
  *( sum + 0 ) = OneO3 * ( nn1 * REsum0 + nn2 * REsum02 ) ;  
  *( sum + 1 ) = OneOI2 * ( nn1 * REsum1 + nn2 * REsum12 ) +\
    0.5 * ( nn1 * IMsum1 + nn2 * IMsum12 ) ;
  *( sum + 2 ) = OneOI2 * ( nn1 * REsum2 + nn2 * REsum22 ) +\
    0.5 * ( nn1 * IMsum2 + nn2 * IMsum22 ) ;
  *( sum + 3 ) = OneO3 * ( nn1 * REsum3 + nn2 * REsum32 ) ;  
  *( sum + 4 ) = OneOI2 * ( nn1 * REsum4 + nn2 * REsum42 ) +\
    0.5 * ( nn1 * IMsum4 + nn2 * IMsum42 ) ;
#elif NC == 2 
  *( sum + 0 ) = ( nn1 * REsum0 + nn2 * REsum02 ) ;
  *( sum + 1 ) = -I * ( nn1 * REsum1 + nn2 * REsum12 ) +\
    0.5 * ( nn1 * IMsum1 + nn2 * IMsum12 ) ;
#endif

  return trace_deriv( sum ) ;
}

