/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (effs.c) is part of GLU.

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
   @file effs.c
   @brief calculation of the f-constants used in the technique of "exact exponentiation"

   The paper for this is the well-cited <a href="linkURL"> http://arxiv.org/abs/hep-lat/0311018 </a>
 */

#include "Mainfile.h"

#include "vandermonde.h" // golub and van loan generic VDM solver is pretty unstable

// cheap and dirty solver -> OK for hermitian eigenvalues
void 
calculate_effs_VDM_herm( double complex *__restrict f , 
			 const double *__restrict z )
{
  size_t i ; 
  for( i = 0 ; i < NC ; i++ ) {
    f[ i ] = cos( z[ i ] ) + I * sin( z[i] ) ;
  }
  solvandermonde_new( z , f ) ; 
  return ;
}

// computes the f's given the exponentiated eigenvalues
void 
calculate_effs_VDM_suNC( double complex *__restrict f , 
			 const double complex *__restrict z )
{
  double x[ NC ] ;
  size_t i ; 
  for( i = 0 ; i < NC ; i++ ) {
    f[ i ] = z[ i ] ; 
    x[ i ] = carg( z[ i ] ) ;
  }
  solvandermonde( f , x ) ; 
  return ;
}

// compute the f's from the Hermitian eigenvalues
void 
#if NC == 2
f_hermitian_log_suNC( double complex f[ NC ] , 
		      const double z )
#else
f_hermitian_log_suNC( double complex f[ NC ] , 
		      const double complex z[ NC ] )
#endif
{
#if NC == 3
  register const double a = 0.5 * creal( z[0] ) ;
  const double u = a ;
  register const double w = creal( z[1] ) + a ;
  const double cw = cos( w ) ;

  register const double cu = cos( u ) ;
  register const double su = sin( u ) ; //sqrt( 1 - cu * cu ) ; //sin( u ) ;
  const double complex one1 = cu - I * su ;
  double complex two1 = conj( one1 ) ; //cu + I * su ;
  two1 *= two1 ;  

  const double uu = u * u ;
  const double ww = w * w ; 

  const double denom = 1.0 / ( 9. * uu - ww ) ;
  // we only allow the taylor expansion at a very low acc to allow the 3D
  // hyl smearing to not get stuck, this is similar to the problem I saw with
  // the su2 log and stout expansions
  const double E0 = fabs( w ) < STOL ? 1 - ww / 6. * ( 1 - ww / 20. * ( 1 - ww / 42. ) ) : sin( w ) / w ; 
  
  //const register double twou = 2. * u ;
  f[0] = ( uu - ww ) * two1 + one1 * ( 8. * uu * cw + I * 2. * u * ( 3. * uu + ww ) * E0 ) ; 
  f[1] = 2 * u * two1 - one1 * ( 2. * u * cw - I * ( 3. * uu - ww ) * E0 ) ; 
  f[2] = two1 - one1 * ( cw + 3. * I * u * E0 ) ; 

  f[0] *= denom ; 
  f[1] *= denom ; 
  f[2] *= denom ;  
#elif NC == 2
  f[0] = cos( z ) ;
  f[1] = ( fabs( z ) < SINTOLSU2 ) ? I * ( 1 - z / 6. * ( 1 - z / 20. * ( 1 - z / 42. ) ) ) : I * sin( z ) / z ; 
#else
  calculate_effs_VDM_suNC( f , z ) ;
#endif
  return ;
}

