/*
    Copyright 2013 Renwick James Hudspith

    This file (U1_obs.c) is part of GLU.

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
   @file U1_obs.c
   @brief U1 observable measurements 
 */

#include "Mainfile.h"

#include "geometry.h"
#include "U1_top.h"

// calculation of the non-compact U(1) plaquette
static double 
non_plaquette( double *__restrict plaq ,
	       const GLU_real *__restrict *__restrict O )
{
  const double denom = 1. / (double)( LVOLUME * ND * ( ND - 1 ) ) ;
  size_t i ;
  double plaquette = 0. , sum = 0. ;
#pragma omp parallel for private(i) reduction(+:sum) reduction(+:plaquette) 
  for( i = 0 ; i < LVOLUME ; i++ ) {
    register double loc_sum = 0.0 , loc_plaq = 0.0 ;
    size_t mu , nu , s , t ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      s = gen_shift( i , mu ) ;
      for( nu = 0 ; nu < mu ; nu++ ) {
	t = gen_shift( i , nu ) ;
	register const double temp = (double)O[mu][i] + \
	  (double)O[nu][s] -				\
	  (double)O[mu][t] -				\
	  (double)O[nu][i] ;
	loc_sum += cos( temp ) ;
	loc_plaq += ( temp * temp ) ;
      }
    }
    plaquette = plaquette + (double)loc_plaq ;
    sum = sum + (double)loc_sum ;
  }
  *plaq = sum * denom ;
  return plaquette * denom ;
}

// computation of the U1 rectangle
static double
U1_rectangle( double *__restrict U1REC ,
	      const GLU_real *__restrict *__restrict O )
{
  const double denom = 1. / ( LVOLUME * ND * ( ND - 1 ) ) ;
  size_t i ;

  double sum_nc = 0. , sum = 0. ;
#pragma omp parallel for private(i) reduction(+:sum) reduction(+:sum_nc)
  for( i = 0 ; i < LVOLUME ; i++ ) {

    double loc_rec = 0. , loc_ncrec = 0. ;
    size_t mu , nu , s , t , u , v ;
    // first one is the (2x1) rectangle
    for( mu = 0 ; mu < ND ; mu++ ) {
      s = gen_shift( i , mu ) ;
      t = gen_shift( s , mu ) ;
      for( nu = 0 ; nu < mu ; nu++ ) {
	v = gen_shift( i , nu ) ;
	u = gen_shift( v , mu ) ;
	register const double cache = (double)O[mu][i] +		\
	  (double)O[mu][s] +						\
	  (double)O[nu][t] -						\
	  (double)O[mu][u] -						\
	  (double)O[mu][v] -						\
	  (double)O[nu][i] ;
	loc_rec += cos( cache ) ; // taking the cosine is expensive - think!
	loc_ncrec += cache * cache ;
      }
    }
    // second one is the (1x2) rectangle
    for( mu = 0 ; mu < ND ; mu++ ) {
      s = gen_shift( i , mu ) ;
      for( nu = 0 ; nu < mu ; nu++ ) {
	t = gen_shift( s , nu ) ;
	v = gen_shift( i , nu ) ;
	u = gen_shift( v , nu ) ;
	register const double cache = (double)O[mu][i] + (double)O[nu][s] + \
	  (double)O[nu][t] - (double)O[mu][u] -\
	  (double)O[mu][v] - (double)O[nu][i] ;
	loc_rec += cos( cache ) ; // taking the cosine is expensive - think!
	loc_ncrec += cache * cache ;
      }
    }
    // looks pretty sexy
    sum = sum + (double)( loc_rec * 0.5 ) ;
    sum_nc = sum_nc + (double)( loc_ncrec * 0.5 ) ;
  }
  *U1REC = sum * denom ;
  return sum_nc * denom ;
}

// wrapper for the U1 configurations
void
compute_U1_obs( const GLU_real *__restrict *__restrict U , 
		const U1_meas meas )
{
  double plaq = 0. ;
  const double test_noncompact = non_plaquette( &plaq , (const GLU_real**)U ) ;
  const double test_compact = plaq ;
  // should have our U1-ified fields
  printf("\n[U(1)] Plaquettes [NON-COMPACT] %lf  [COMPACT] %lf \n" , test_noncompact , test_compact ) ;

  if( meas == U1_RECTANGLE ) {
    const double test_rectangle = U1_rectangle( &plaq , (const GLU_real**)U ) ;
    printf("\n[U(1)] Rectangle  [NON-COMPACT] %lf  [COMPACT] %lf \n" , test_rectangle , plaq ) ;
  } else if( meas == U1_TOPOLOGICAL ) {
    double qtop ;
    int monopole , dirac_sheet ;
    U1_topological( &monopole , &dirac_sheet , &qtop , (const GLU_real**)U ) ;
    printf("\n[U(1)] Topological :: [Monopole] %d  [Dirac_Sheet] %d [QTOP] %g \n" , monopole , dirac_sheet , qtop ) ;
  }

  return ;
}
