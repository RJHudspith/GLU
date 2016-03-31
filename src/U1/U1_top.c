/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (U1_top.c) is part of GLU.

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
   @file U1_top.c
   @brief topological measurements using the non-compact U(1) code
 */

#include "Mainfile.h"
#include "geometry.h"

// compute the dirac current
static int 
current( const int *__restrict *__restrict M )
{
  size_t mu , monopole = 0 ;
#pragma omp parallel for private(mu) reduction(+:monopole)
  for( mu = 0 ; mu < ND ; mu++ ) {
    size_t loc_monopole = 0 , i ;
    for( i = 0 ; i < LVOLUME ; i++ )  {
      loc_monopole += abs( M[ mu ][ i ] ) ;
    }
    monopole = monopole + (int)loc_monopole ;
  }
  return monopole;
}

// the mod 2pi bit to get the number of windings
static inline int 
mod_2pi( const GLU_real theta )
{
  return (int)theta / TWOPI ; 
}

// the noncompact face of the plaquette
static GLU_real
noncompact_face( const GLU_real *__restrict *__restrict O , 
		 const size_t i , 
		 const size_t mu , 
		 const size_t nu )
{
  const GLU_real a = O[mu][i];
  size_t temp = gen_shift( i , mu ) ;
  const GLU_real b = O[nu][temp] ;
  const GLU_real d = O[nu][i] ;
  temp = gen_shift( i , nu ) ;
  const GLU_real c = O[mu][temp] ;
  return ( a + b - c - d ) ;
}

// calculates the "s", the winding per face of the plaquette
inline static int 
obtain_S( const GLU_real *__restrict *__restrict O ,
	  const size_t i , 
	  const size_t mu , 
	  const size_t nu )
{
  const GLU_real theta = noncompact_face( O , i , mu , nu ) ;
  return mod_2pi( theta ) ;
}

// compute the "dirac sheet" from the noncompact plaquette
static int 
dirac_sheet( const GLU_real *__restrict *__restrict O )
{
  size_t i ;
  int sheet = 0 ;
#pragma omp parallel for private(i) reduction(+:sheet)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    size_t mu , nu ; 
    int loc_sheet = 0 ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      for( nu = 0 ; nu < mu ; nu++ )	  {
	const int temp = obtain_S( O , i , mu , nu ) ;
	loc_sheet += abs( temp ) ;
      }
    }
    sheet = sheet + (int)loc_sheet ;
  }
  return sheet;
}

// computes the dirac observable "M"
static void 
dirac( int *__restrict *__restrict M ,
       const GLU_real *__restrict *__restrict O )
{
  size_t i ;
#pragma omp parallel for private(i) 
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {// is this OK?
    size_t mu , nu , rho , sigma;
    for( mu = 0 ; mu < ND ; mu++ ) {
      for( nu = 0 ; nu < ND ; nu++ ) {
	if(nu != mu)  {
	  for( rho = 0 ; rho < ND ; rho++ ) {
	    for( sigma = 0 ; sigma < rho ; sigma++ ) {
	      if( sigma == mu || rho == mu || sigma == nu || rho == nu ) {
		continue ;
	      } else {
		int factor = -1 ;
		if( ( nu % 2 ) == 1 ) {
		  factor = 1 ;
		}
		const size_t temp = gen_shift( i , nu ) ;
		M[ mu ][ i ] = M[ mu ][ i ] + factor * ( obtain_S( O , temp , rho , sigma )
							 - obtain_S( O , i , rho , sigma ) ) ;
	      }
	    }
	  }
	}
      }
    }
  }
  return ;
}

// very naive topological charge using the noncompact plaquettes 
// e_\mu\nu\rho\eta G_\mu\nu G_\rho\eta 
static double 
non_Qtop( const GLU_real *__restrict *__restrict O )
{
  double charge = 0.0 ;
  size_t i;
#pragma omp parallel for private(i) reduction(+:charge)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    double qtemp = (double)noncompact_face( O , i , 0 , 1 ) * (double)noncompact_face( O , i , 2 , 3 ) ;
    qtemp += (double)noncompact_face( O , i , 0 , 2 ) * (double)noncompact_face( O , i , 3 , 1 ) ; 
    qtemp += (double)noncompact_face( O , i , 0 , 3 ) * (double)noncompact_face( O , i , 1 , 2 ) ; 
    
    // parallel reduction on charge
    charge = charge + (double)qtemp ;
  }
  return charge * 0.001583143494411527678811 ;
}

// compute all of the topological measures here
void 
U1_topological( int *__restrict monopole , 
		int *__restrict d_sheet , 
		double *__restrict qtop , 
		const GLU_real *__restrict *__restrict O  )
{
  int **M = calloc( ND , sizeof( int * ) ) ;
  size_t mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    M[ mu ] = ( int* )calloc( LVOLUME , sizeof( int ) ) ; 
  }

  // initialise the dirac measurement
  dirac( M , (const GLU_real**)O ) ;
  
  // calculate the observables
  *monopole = current( (const int**)M ) ;
  *d_sheet = dirac_sheet( O ) ;
  *qtop = non_Qtop( O ) ;

  for( mu = 0 ; mu < ND ; mu++ ) {
    free( M[ mu ] ) ;
  }
  free( M ) ;

  return ;
}
