/*
Copyright 2013-2025 Renwick James Hudspith

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
non_plaquette( double *plaq ,
	       const struct site *lat ,
	       const GLU_complex **O )
{
  const double denom = 2. / (double)( LVOLUME * ND * ( ND - 1 ) ) ;
  size_t i ;
  double plaquette = 0. , sum = 0. ;
#pragma omp parallel for private(i) reduction(+:sum) reduction(+:plaquette) 
  for( i = 0 ; i < LVOLUME ; i++ ) {
    register double loc_sum = 0.0 , loc_plaq = 0.0 ;
    size_t mu , nu , s , t ;
    for( mu = 1 ; mu < ND ; mu++ ) {
      s = lat[i].neighbor[mu] ;
      for( nu = 0 ; nu < mu ; nu++ ) {
	t = lat[i].neighbor[nu] ;
	register const double temp =				\
	  creal( (double)O[mu][i] + (double)O[nu][s] -		\
		 (double)O[mu][t] - (double)O[nu][i] ) ;
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
U1_rectangle( double *U1REC ,
	      const struct site *lat ,
	      const GLU_complex **O )
{
  const double denom = 1. / (double)( LVOLUME * ND * ( ND - 1 ) ) ;
  size_t i ;
  double sum_nc = 0. , sum = 0. ;
#pragma omp parallel for private(i) reduction(+:sum) reduction(+:sum_nc)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    double loc_rec = 0. , loc_ncrec = 0. ;
    size_t mu , nu , s , t , t2 , u , v ;
    for( mu = 1 ; mu < ND ; mu++ ) {
      s = lat[i].neighbor[mu] ;
      t = lat[s].neighbor[mu] ;
      for( nu = 0 ; nu < mu ; nu++ ) {
	v = lat[i].neighbor[nu] ;
	u = lat[v].neighbor[mu] ;
	// first one is the (2x1) rectangle
	register double cache =						\
	  creal( (double)O[mu][i] + (double)O[mu][s] + (double)O[nu][t] -
		 (double)O[mu][u] - (double)O[mu][v] - (double)O[nu][i] ) ;
	loc_rec += cos( cache ) ;
	loc_ncrec += cache * cache ;
	// second one is the (1x2) rectangle
	t2 = lat[s].neighbor[nu] ;
	u  = lat[v].neighbor[nu] ;
        cache =								\
	  creal( (double)O[mu][i] + (double)O[nu][s] + (double)O[nu][t2] -
		 (double)O[mu][u] - (double)O[nu][v] - (double)O[nu][i] ) ;
	loc_rec += cos( cache ) ;
	loc_ncrec += cache * cache ;
      }
    }
    sum = sum + (double)( loc_rec ) ;
    sum_nc = sum_nc + (double)( loc_ncrec ) ;
  }
  *U1REC = sum * denom ;
  return sum_nc * denom ;
}

// compute the norm of the fields
static double
normU( const GLU_complex **U )
{
  size_t i ;
  double sum = 0 ;
#pragma omp parallel for private(i) reduction(+:sum)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      sum = sum + U[mu][i];
    }
  }
  return sqrt( sum*sum / LVOLUME ) ;
}

// wrapper for the U1 configurations
void
compute_U1_obs( const GLU_complex **U ,
		const struct site *lat ,
		const U1_meas meas )
{
  double plaq = 0. ;
  const double test_noncompact = \
    non_plaquette( &plaq , lat , (const GLU_complex**)U ) ;
  const double test_compact = plaq ;
  fprintf( stdout , "\n[U(1)] Norm U (should be zero) %e\n" , normU( U ) ) ;
  
  // should have our U1-ified fields
  fprintf( stdout , "\n[U(1)] Plaquettes [NON-COMPACT] %lf [COMPACT] %lf \n" ,
	   test_noncompact , test_compact ) ;

  if( meas == U1_RECTANGLE ) {
    const double test_rectangle = \
      U1_rectangle( &plaq , lat , (const GLU_complex**)U ) ;
    fprintf( stdout , "\n[U(1)] Rectangle  [NON-COMPACT] %lf"
	     " [COMPACT] %lf \n" , test_rectangle , plaq ) ;
  } else if( meas == U1_TOPOLOGICAL ) {
    double qtop ;
    int monopole , dirac_sheet ;
    U1_topological( &monopole , &dirac_sheet , &qtop , (const GLU_complex**)U ) ;
    fprintf( stdout , "\n[U(1)] Topological :: [Monopole] %d"
	     "[Dirac_Sheet] %d [QTOP] %g \n" , monopole , dirac_sheet , qtop ) ;
  }

  return ;
}
