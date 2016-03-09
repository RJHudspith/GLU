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
  // standard gauge transform
  GLU_complex temp[ NCNC ] GLUalign ;
  multab_dag_suNC( temp , b , c ) ;
  multab_suNC( b , a , temp ) ; 
  return ;
}

//gauge_transform lattice-wide
void 
gtransform( struct site *__restrict lat ,
	    const GLU_complex *__restrict *__restrict gauge )
{
  size_t i ;
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) {
    #if ND == 4
    gtransform_local( gauge[i] , lat[i].O[0] , gauge[lat[i].neighbor[0]] ) ;
    gtransform_local( gauge[i] , lat[i].O[1] , gauge[lat[i].neighbor[1]] ) ;
    gtransform_local( gauge[i] , lat[i].O[2] , gauge[lat[i].neighbor[2]] ) ;
    gtransform_local( gauge[i] , lat[i].O[3] , gauge[lat[i].neighbor[3]] ) ;
    #else
    size_t mu ;
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
		  const size_t t )
{
  size_t i ;
  const size_t slice = LCU  *  t ; 
#pragma omp parallel for private(i)
  PFOR(  i = 0  ;  i < LCU  ;  i ++ ) {
    const size_t j = slice + i ;
    #if ND == 4   
    gtransform_local( gauge[i] , lat[j].O[0] , gauge[lat[i].neighbor[0]] ) ;
    gtransform_local( gauge[i] , lat[j].O[1] , gauge[lat[i].neighbor[1]] ) ;
    gtransform_local( gauge[i] , lat[j].O[2] , gauge[lat[i].neighbor[2]] ) ;
    #else
    size_t mu ;
    for( mu = 0 ; mu < ND - 1  ; mu++ ) {
      gtransform_local( gauge[i] , lat[j].O[mu] , gauge[lat[i].neighbor[mu]] ) ;
    }
    #endif
    gtransform_local( gauge[i] , lat[j].O[ND-1] , gauge_up[i] ) ;
  }
  return ;
}
