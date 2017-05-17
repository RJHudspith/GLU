/*
    Copyright 2013-2016 Renwick James Hudspith

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

// chroma's seems the cheapest at the moment
#define CHROMA_RELAX

// number of SOR steps
//#define SOR

// number of stochastic SU(2) updates
//#define NSTOCH (NC*15)

// microcanonical su(2) update
void
microcanonical( GLU_complex *s0 ,
		GLU_complex *s1 )
{
  // chroma's is better as it has one fewer rotation
#ifdef CHROMA_RELAX
  // compute z
  register const double z = ( creal(*s0)*creal(*s0) - cimag(*s0)*cimag(*s0) -
			      creal(*s1)*creal(*s1) - cimag(*s1)*cimag(*s1) ) ;
  *s1 = -2. * creal( *s0 ) * ( *s1 ) ;
  *s0 = z - I * ( 2. * creal(*s0) * cimag(*s0) ) ;
#else
  // milc relax
  *s0 =  creal( *s0 ) - I * cimag( *s0 ) ;
  *s1 = -creal( *s1 ) - I * cimag( *s1 ) ;
#endif
  return ;
}

// overrelaxation algorithm
static void
overrelax( GLU_complex U[ NCNC ] , 
	   const GLU_complex staple[ NCNC ] ,
	   const uint32_t thread )
{
  GLU_complex s0 GLUalign , s1 GLUalign ;
  double scale GLUalign ;
  size_t i ;
#ifdef NSTOCH
  for( i = 0 ; i < NSTOCH ; i++ ) {
    const size_t stoch = (size_t)( par_rng_dbl( thread ) * NSU2SUBGROUPS ) ;
    only_subgroup( &s0 , &s1 , &scale , U , staple , stoch ) ;
    microcanonical( &s0 , &s1 ) ;
    #ifdef CHROMA_RELAX
    su2_rotate( U , s0 , s1 , stoch ) ;
    #else
    su2_rotate( U , s0 , s1 , stoch ) ;  
    su2_rotate( U , s0 , s1 , stoch ) ;
    #endif
  }
#else
  // maybe this is a good idea ...
  size_t k ;
  for( i = 0 ; i < NSU2SUBGROUPS ; i++ ) {
    #ifdef SOR
    // stochastic-OR?
    if( par_rng_dbl( thread ) < 0.5 ) continue ;
    #endif
    only_subgroup( &s0 , &s1 , &scale , U , staple , i ) ;
    microcanonical( &s0 , &s1 ) ;
    #ifdef CHROMA_RELAX
    su2_rotate( U , s0 , s1 , i ) ;
    #else
    su2_rotate( U , s0 , s1 , i ) ;  
    su2_rotate( U , s0 , s1 , i ) ;
    #endif
  }
#endif
  return ;
}

// perform a heat-bath over the whole lattice
int
OR_lattice( struct site *lat ,
	    const struct draughtboard db )
{
  size_t mu , i ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    // single node until I get the coloring correct
#ifdef IMPROVED_SMEARING
    for( i = 0 ; i < LVOLUME ; i++ ) {
      GLU_complex stap[ NCNC ] GLUalign ;
      zero_mat( stap ) ;
      #ifdef IMPROVED_SMEARING
      all_staples_improve( stap , lat , i , mu , ND , SM_APE ) ;
      #else
      all_staples( stap , lat , i , mu , ND , SM_APE ) ;
      #endif
      overrelax( lat[ i ].O[mu] , stap , get_GLU_thread() ) ;
    }
#else
    // loop draughtboard coloring
    size_t c ;
    for( c = 0 ; c < db.Ncolors ; c++ ) {
      #pragma omp parallel for private(i)
      for( i = 0 ; i < db.Nsquare[c] ; i++ ) {
	const size_t it = db.square[c][i] ;
	GLU_complex stap[ NCNC ] GLUalign ;
	zero_mat( stap ) ;
	#ifdef IMPROVED_SMEARING
	all_staples_improve( stap , lat , it , mu , ND , SM_APE ) ;
	#else
	all_staples( stap , lat , it , mu , ND , SM_APE ) ;
	#endif
	overrelax( lat[ it ].O[mu] , stap , get_GLU_thread() ) ;
      }
      // and that is it
    }
#endif
  }
  return GLU_SUCCESS ;
}

// clean up the over-relaxation
#ifdef CHROMA_RELAX
  #undef CHROMA_RELAX
#endif
