/*
    Copyright 2013 Renwick James Hudspith

    This file (gtrans.c) is part of GLU.

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
   @file gtrans.c
   @brief gauge transformations overwriting the lattice field
*/

#include "Mainfile.h"

/**
   @fn static void gtransform_local( const GLU_complex *__restrict a , GLU_complex *__restrict b , const GLU_complex *__restrict c )
   @brief localised gauge tranform of the form b = a.b.c^{\dagger}

   I have seen marginal performance increases from this so I figure I will stick to it if anyone wants to
   optimise the Landau and Coulomb code this is a good place to start
 */
void
gtransform_local( const GLU_complex *__restrict a ,
		  GLU_complex *__restrict b ,
		  const GLU_complex *__restrict c )
{
#if NC == 3
  // product of three SU(3) matrices?
  //const GLU_complex c0 = b[0] * conj( c[0] ) + b[1] * conj( c[1] ) + b[2] * conj( c[2] ) ;
  const GLU_real REc0 = + creal( b[0] ) * creal( c[0] ) + cimag( b[0] ) * cimag( c[0] ) 
                        + creal( b[1] ) * creal( c[1] ) + cimag( b[1] ) * cimag( c[1] ) 
                        + creal( b[2] ) * creal( c[2] ) + cimag( b[2] ) * cimag( c[2] ) ;
  const GLU_real IMc0 = - ( creal( b[0] ) * cimag( c[0] ) - creal( c[0] ) * cimag( b[0] ) )
                        - ( creal( b[1] ) * cimag( c[1] ) - creal( c[1] ) * cimag( b[1] ) )
                        - ( creal( b[2] ) * cimag( c[2] ) - creal( c[2] ) * cimag( b[2] ) ) ;

  // b[0] * conj( c[3] ) + b[1] * conj( c[4] ) + b[2] * conj( c[5] ) ;
  const GLU_real REc1 = + creal( b[0] ) * creal( c[3] ) + cimag( b[0] ) * cimag( c[3] ) 
                        + creal( b[1] ) * creal( c[4] ) + cimag( b[1] ) * cimag( c[4] ) 
                        + creal( b[2] ) * creal( c[5] ) + cimag( b[2] ) * cimag( c[5] ) ;
  const GLU_real IMc1 = -( creal( b[0] ) * cimag( c[3] ) - creal( c[3] ) * cimag( b[0] ) )
                        -( creal( b[1] ) * cimag( c[4] ) - creal( c[4] ) * cimag( b[1] ) ) 
                        -( creal( b[2] ) * cimag( c[5] ) - creal( c[5] ) * cimag( b[2] ) ) ;

  // b[0] * conj( c[6] ) + b[1] * conj( c[7] ) + b[2] * conj( c[8] ) ;
  const GLU_real REc2 = + creal( b[0] ) * creal( c[6] ) + cimag( b[0] ) * cimag( c[6] ) 
                        + creal( b[1] ) * creal( c[7] ) + cimag( b[1] ) * cimag( c[7] ) 
                        + creal( b[2] ) * creal( c[8] ) + cimag( b[2] ) * cimag( c[8] ) ;
  const GLU_real IMc2 = -( creal( b[0] ) * cimag( c[6] ) - creal( c[6] ) * cimag( b[0] ) ) 
                        -( creal( b[1] ) * cimag( c[7] ) - creal( c[7] ) * cimag( b[1] ) )
                        -( creal( b[2] ) * cimag( c[8] ) - creal( c[8] ) * cimag( b[2] ) ) ;

  // b[3] * conj( c[0] ) + b[4] * conj( c[1] ) + b[5] * conj( c[2] ) ;
  const GLU_real REc3 = + creal( b[3] ) * creal( c[0] ) + cimag( b[3] ) * cimag( c[0] ) 
                        + creal( b[4] ) * creal( c[1] ) + cimag( b[4] ) * cimag( c[1] ) 
                        + creal( b[5] ) * creal( c[2] ) + cimag( b[5] ) * cimag( c[2] ) ;
  const GLU_real IMc3 = -( creal( b[3] ) * cimag( c[0] ) - creal( c[0] ) * cimag( b[3] ) )
                        -( creal( b[4] ) * cimag( c[1] ) - creal( c[1] ) * cimag( b[4] ) )
                        -( creal( b[5] ) * cimag( c[2] ) - creal( c[2] ) * cimag( b[5] ) ) ;

  // b[3] * conj( c[3] ) + b[4] * conj( c[4] ) + b[5] * conj( c[5] ) ;
  const GLU_real REc4 = + creal( b[3] ) * creal( c[3] ) + cimag( b[3] ) * cimag( c[3] ) 
                        + creal( b[4] ) * creal( c[4] ) + cimag( b[4] ) * cimag( c[4] ) 
                        + creal( b[5] ) * creal( c[5] ) + cimag( b[5] ) * cimag( c[5] ) ;
  const GLU_real IMc4 = -( creal( b[3] ) * cimag( c[3] ) - creal( c[3] ) * cimag( b[3] ) )
                        -( creal( b[4] ) * cimag( c[4] ) - creal( c[4] ) * cimag( b[4] ) )
                        -( creal( b[5] ) * cimag( c[5] ) - creal( c[5] ) * cimag( b[5] ) ) ;

  // b[3] * conj( c[6] ) + b[4] * conj( c[7] ) + b[5] * conj( c[8] ) ;
  const GLU_real REc5 = + creal( b[3] ) * creal( c[6] ) + cimag( b[3] ) * cimag( c[6] ) 
                        + creal( b[4] ) * creal( c[7] ) + cimag( b[4] ) * cimag( c[7] ) 
                        + creal( b[5] ) * creal( c[8] ) + cimag( b[5] ) * cimag( c[8] ) ;
  const GLU_real IMc5 = -( creal( b[3] ) * cimag( c[6] ) - creal( c[6] ) * cimag( b[3] ) )
                        -( creal( b[4] ) * cimag( c[7] ) - creal( c[7] ) * cimag( b[4] ) ) 
                        -( creal( b[5] ) * cimag( c[8] ) - creal( c[8] ) * cimag( b[5] ) ) ;

  // conj( c1 * c5 - c2 * c4 ) ;
  const GLU_real REc6 = REc1 * REc5 - IMc1 * IMc5 - REc2 * REc4 + IMc2 * IMc4 ;
  const GLU_real IMc6 = ( REc2 * IMc4 + REc4 * IMc2 ) - ( REc1 * IMc5 + REc5 * IMc1 ) ; 

  // -conj( c0 * c5 - c2 * c3 ) ;
  const GLU_real REc7 = REc2 * REc3 - IMc2 * IMc3 - REc0 * REc5 + IMc0 * IMc5 ;
  const GLU_real IMc7 = ( REc0 * IMc5 + REc5 * IMc0 ) - ( REc2 * IMc3 + REc3 * IMc2 ) ; 

  // conj( c0 * c4 - c1 * c3 ) ;
  const GLU_real REc8 = REc0 * REc4 - IMc0 * IMc4 - REc1 * REc3 + IMc1 * IMc3 ;
  const GLU_real IMc8 = ( REc1 * IMc3 + REc3 * IMc1 ) - ( REc0 * IMc4 + REc4 * IMc0 ) ;

  //b[0] = a[0] * c0 + a[1] * c3 + a[2] * c6 ;
  b[0] = creal( a[0] ) * REc0 - cimag( a[0] ) * IMc0	
    + I * ( creal( a[0] ) * IMc0 + REc0 * cimag( a[0] ) ) 
    + creal( a[1] ) * REc3 - cimag( a[1] ) * IMc3		
    + I * ( creal( a[1] ) * IMc3 + REc3 * cimag( a[1] ) ) 
    + creal( a[2] ) * REc6 - cimag( a[2] ) * IMc6		
    + I * ( creal( a[2] ) * IMc6 + REc6 * cimag( a[2] ) ) ; 

  //b[1] = a[0] * c1 + a[1] * c4 + a[2] * c7 ;
  b[1] = creal( a[0] ) * REc1 - cimag( a[0] ) * IMc1	
    + I * ( creal( a[0] ) * IMc1 + REc1 * cimag( a[0] ) ) 
    + creal( a[1] ) * REc4 - cimag( a[1] ) * IMc4		
    + I * ( creal( a[1] ) * IMc4 + REc4 * cimag( a[1] ) ) 
    + creal( a[2] ) * REc7 - cimag( a[2] ) * IMc7		
    + I * ( creal( a[2] ) * IMc7 + REc7 * cimag( a[2] ) ) ; 

  //b[2] = a[0] * c2 + a[1] * c5 + a[2] * c8 ;
  b[2] = creal( a[0] ) * REc2 - cimag( a[0] ) * IMc2	
    + I * ( creal( a[0] ) * IMc2 + REc2 * cimag( a[0] ) ) 
    + creal( a[1] ) * REc5 - cimag( a[1] ) * IMc5		
    + I * ( creal( a[1] ) * IMc5 + REc5 * cimag( a[1] ) ) 
    + creal( a[2] ) * REc8 - cimag( a[2] ) * IMc8		
    + I * ( creal( a[2] ) * IMc8 + REc8 * cimag( a[2] ) ) ; 

  //b[3] = a[3] * c0 + a[4] * c3 + a[5] * c6 ;
  b[3] = creal( a[3] ) * REc0 - cimag( a[3] ) * IMc0	
    + I * ( creal( a[3] ) * IMc0 + REc0 * cimag( a[3] ) ) 
    + creal( a[4] ) * REc3 - cimag( a[4] ) * IMc3		
    + I * ( creal( a[4] ) * IMc3 + REc3 * cimag( a[4] ) ) 
    + creal( a[5] ) * REc6 - cimag( a[5] ) * IMc6		
    + I * ( creal( a[5] ) * IMc6 + REc6 * cimag( a[5] ) ) ; 

  //b[4] = a[3] * c1 + a[4] * c4 + a[5] * c7 ;
  b[4] = creal( a[3] ) * REc1 - cimag( a[3] ) * IMc1	
    + I * ( creal( a[3] ) * IMc1 + REc1 * cimag( a[3] ) ) 
    + creal( a[4] ) * REc4 - cimag( a[4] ) * IMc4		
    + I * ( creal( a[4] ) * IMc4 + REc4 * cimag( a[4] ) ) 
    + creal( a[5] ) * REc7 - cimag( a[5] ) * IMc7		
    + I * ( creal( a[5] ) * IMc7 + REc7 * cimag( a[5] ) ) ; 

  //b[5] = a[3] * c2 + a[4] * c5 + a[5] * c8 ;
  b[5] = creal( a[3] ) * REc2 - cimag( a[3] ) * IMc2	
    + I * ( creal( a[3] ) * IMc2 + REc2 * cimag( a[3] ) ) 
    + creal( a[4] ) * REc5 - cimag( a[4] ) * IMc5		
    + I * ( creal( a[4] ) * IMc5 + REc5 * cimag( a[4] ) ) 
    + creal( a[5] ) * REc8 - cimag( a[5] ) * IMc8		
    + I * ( creal( a[5] ) * IMc8 + REc8 * cimag( a[5] ) ) ; 

  //b[6] = conj( b[1]  *  b[5] - b[2]  *  b[4] ) ; 
  b[6] = creal( b[1] ) * creal( b[5] ) - cimag( b[1] ) * cimag( b[5] ) - 
    I * ( creal( b[1] ) * cimag( b[5] ) + creal( b[5] ) * cimag( b[1] ) ) \
       - creal( b[2] ) * creal( b[4] ) + cimag( b[2] ) * cimag( b[4] ) + 
    I * ( creal( b[2] ) * cimag( b[4] ) + creal( b[4] ) * cimag( b[2] ) ) ; 

  // b[7] = conj( b[2]  *  b[3] - b[0]  *  b[5] ) ; 
  b[7] = creal( b[2] ) * creal( b[3] ) - cimag( b[2] ) * cimag( b[3] ) - 
    I * ( creal( b[2] ) * cimag( b[3] ) + creal( b[3] ) * cimag( b[2] ) ) \
       - creal( b[0] ) * creal( b[5] ) + cimag( b[0] ) * cimag( b[5] ) + 
    I * ( creal( b[0] ) * cimag( b[5] ) + creal( b[5] ) * cimag( b[0] ) ) ; 

  // b[8] = conj( b[0]  *  b[4] - b[1]  *  b[3] ) ; 
  b[8] = creal( b[0] ) * creal( b[4] ) - cimag( b[0] ) * cimag( b[4] ) 
    - I * ( creal( b[0] ) * cimag( b[4] ) + creal( b[4] ) * cimag( b[0] ) ) \
       - creal( b[1] ) * creal( b[3] ) + cimag( b[1] ) * cimag( b[3] ) 
    + I * ( creal( b[1] ) * cimag( b[3] ) + creal( b[3] ) * cimag( b[1] ) ) ;

#elif NC == 2 // gives us about a 5% speedup, interesting
  register const GLU_real REcache1 = ( creal( b[0] ) * creal( c[0] ) + cimag( b[0] ) * cimag( c[0] )
				       + creal( b[1] ) * creal( c[1] ) + cimag( b[1] ) * cimag( c[1] ) ) ;
  register const GLU_real IMcache1 = ( cimag( b[0] ) * creal( c[0] ) - creal( b[0] ) * cimag( c[0] ) 
				       + cimag( b[1] ) * creal( c[1] ) - creal( b[1] ) * cimag( c[1] ) ) ;
  register const GLU_real REcache2 = ( creal( b[0] ) * creal( c[2] ) + cimag( b[0] ) * cimag( c[2] )
				       + creal( b[1] ) * creal( c[3] ) + cimag( b[1] ) * cimag( c[3] ) ) ;
  register const GLU_real IMcache2 = ( cimag( b[0] ) * creal( c[2] ) - creal( b[0] ) * cimag( c[2] ) 
				       + cimag( b[1] ) * creal( c[3] ) - creal( b[1] ) * cimag( c[3] ) ) ;
  register const GLU_real REA0 = creal( a[0] ) , IMA0 = cimag( a[0] ) ;
  register const GLU_real REA1 = creal( a[1] ) , IMA1 = cimag( a[1] ) ;
  b[0] = REA0 * REcache1 - IMA0 * IMcache1 + I * ( REA0 * IMcache1 + IMA0 * REcache1 ) 
       - REA1 * REcache2 - IMA1 * IMcache2 + I * ( REA1 * IMcache2 - IMA1 * REcache2 ) ;
  b[1] = REA0 * REcache2 - IMA0 * IMcache2 + I * ( REA0 * IMcache2 + IMA0 * REcache2 ) 
       + REA1 * REcache1 + IMA1 * IMcache1 + I * ( -REA1 * IMcache1 + IMA1 * REcache1 ) ;
  b[2] = -conj( b[1] ) ;
  b[3] =  conj( b[0] ) ;
#else
  // standard gauge transform
  GLU_complex temp[ NCNC ] ;
  multab_dag_suNC( temp , b , c ) ;
  multab_suNC( b , a , temp ) ; 
#endif
  return ;
}

//gauge_transform lattice-wide
void 
gtransform( struct site *__restrict lat ,
	    const GLU_complex *__restrict *__restrict gauge )
{
  int i ;
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {

    #if ND == 4
    gtransform_local( gauge[i] , lat[i].O[0] , gauge[lat[i].neighbor[0]] ) ;
    gtransform_local( gauge[i] , lat[i].O[1] , gauge[lat[i].neighbor[1]] ) ;
    gtransform_local( gauge[i] , lat[i].O[2] , gauge[lat[i].neighbor[2]] ) ;
    gtransform_local( gauge[i] , lat[i].O[3] , gauge[lat[i].neighbor[3]] ) ;
    #else
    int mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      gtransform_local( gauge[i] , lat[i].O[mu] , gauge[lat[i].neighbor[mu]] ) ;
    }
    #endif

  } 
  return ;
}

// gauge_transform for the Coulomb definition 
void
gtransform_slice( const GLU_complex *__restrict *__restrict gauge , 
		  struct site *__restrict lat , 
		  const GLU_complex *__restrict *__restrict gauge_up ,
		  const int t )
{
  int i ;
  const int slice = LCU  *  t ; 
#pragma omp parallel for private(i)
  PFOR(  i = 0  ;  i < LCU  ;  i ++ ) {
    const int j = slice + i ;

#if ND == 4   
    gtransform_local( gauge[i] , lat[j].O[0] , gauge[lat[i].neighbor[0]] ) ;
    gtransform_local( gauge[i] , lat[j].O[1] , gauge[lat[i].neighbor[1]] ) ;
    gtransform_local( gauge[i] , lat[j].O[2] , gauge[lat[i].neighbor[2]] ) ;
#else
    int mu ;
    for( mu = 0 ; mu < ND - 1  ; mu++ ) {
      gtransform_local( gauge[i] , lat[j].O[mu] , gauge[lat[i].neighbor[mu]] ) ;
    }
#endif
    gtransform_local( gauge[i] , lat[j].O[ND-1] , gauge_up[i] ) ;
  }
  return ;
}
