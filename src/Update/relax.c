/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (relax.c) is part of GLU.

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
   @file relax.c
   @brief over relaxation code lives here
 */
#include "Mainfile.h"

#include "staples.h"    // all_staples()
#include "SU2_rotate.h" // rotation
#include "par_rng.h"
#include "gramschmidt.h"

// number of SOR steps
//#define SOR

// microcanonical su(2) update
void
microcanonical( GLU_complex *s0 ,
		GLU_complex *s1 )
{
  // chroma's is better as it has one fewer rotation
  register const double z = ( creal(*s0)*creal(*s0) - cimag(*s0)*cimag(*s0) -
			      creal(*s1)*creal(*s1) - cimag(*s1)*cimag(*s1) ) ;
  *s1 = -2. * creal( *s0 ) * ( *s1 ) ;
  *s0 = z - I * ( 2. * creal(*s0) * cimag(*s0) ) ;
  return ;
}

// overrelaxation algorithm
static void
overrelax( GLU_complex U[ NCNC ] , 
	   const GLU_complex staple[ NCNC ] )
{
  #ifdef SOR
  const uint32_t thread = get_GLU_thread() ;
  #endif
  GLU_complex s0 GLUalign , s1 GLUalign ;
  double scale GLUalign ;
  size_t i ;
  for( i = 0 ; i < NSTOCH ; i++ ) {
    #if (NSTOCH != NSU2SUBGROUPS)
    const size_t stoch = (size_t)( par_rng_dbl( thread ) * NSU2SUBGROUPS ) ;
    #else
    const size_t stoch = i ;
    #endif
    // stochastic-OR?
    #ifdef SOR
    if( par_rng_dbl( thread ) < 0.5 ) continue ;
    #endif
    only_subgroup( &s0 , &s1 , &scale , U , staple , stoch ) ;
    microcanonical( &s0 , &s1 ) ;
    su2_rotate( U , s0 , s1 , stoch ) ;
  }
  return ;
}

// perform a heat-bath over the whole lattice
int
OR_lattice( struct site *lat ,
	    const struct draughtboard db )
{
    // single node until I get the coloring correct
#ifdef IMPROVED_SMEARING
  size_t cmu , i ;
  for( cmu = 0 ; cmu < db.Ncolors ; cmu++ ) {
    const size_t mu = (size_t)cmu/4 ;
    // parallel loop over all sites with this coloring
    #pragma omp for private(i) 
    for( i = 0 ; i < db.Nsquare[cmu] ; i++ ) {
      const size_t it = db.square[cmu][i] ;
      GLU_complex stap[ NCNC ] GLUalign ;
      zero_mat( stap ) ;
      all_staples_improve( stap , lat , it , mu , ND , SM_APE ) ;
      overrelax( lat[ it ].O[mu] , stap ) ;
    }
  }
#else
  size_t cmu , i ;
  for( cmu = 0 ; cmu < ND*db.Ncolors ; cmu++ ) {
    // loop draughtboard coloring
    const size_t mu = cmu/db.Ncolors ;
    const size_t c  = cmu%db.Ncolors ;
#pragma omp for private(i)
    for( i = 0 ; i < db.Nsquare[c] ; i++ ) {
      const size_t it = db.square[c][i] ;
      GLU_complex stap[ NCNC ] GLUalign ;
      zero_mat( stap ) ;
      all_staples( stap , lat , it , mu , ND , SM_APE ) ;
      overrelax( lat[ it ].O[mu] , stap ) ;
    }
    // and that is it
  }
#endif
  return GLU_SUCCESS ;
}

#ifdef SOR
  #undef SOR
#endif
