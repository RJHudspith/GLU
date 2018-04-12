/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (Coulomb.c) is part of GLU.

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
  @file Coulomb.c
  @brief coulomb gauge fixing code, fixes each time-slice in turn and restarts with a random gauge transform if there is a problem ...
 */
#include "Mainfile.h"

#include "CFACG.h"        // Coulomb Conjugate Gradient
#include "geometry.h"     // for the p^2 calc
#include "plan_ffts.h"    // fftw plan wrapper
#include "plaqs_links.h"  // plaqutte and link measurements

#ifdef HAVE_FFTW3_H
void
set_FFTW( struct fftw_stuff *FFTW )
{

}
#endif

// Coulomb gauge fixing code
size_t
Coulomb( struct site *__restrict lat , 
	 const double accuracy , 
	 const size_t iter )
{
  double splink , tlink ;
  all_links( lat , &splink , &tlink ) ;

  fprintf( stdout , "[GF] Initial Tlink :: %1.15f || Slink :: %1.15f \n"
	   "[GF] Plaquette :: %1.15f \n", 
	   tlink , splink , av_plaquette( lat ) ) ; 

  struct fftw_stuff FFTW ;
  size_t i ;
  
  // put an ifdef guard here as SD requires none of this ...
#ifdef HAVE_FFTW3_H

  FFTW.psq = malloc( LCU * sizeof( GLU_real ) ) ; 
  
  // we calculate the lattice p-squared here and pass it to the FFT-accelerator
#pragma omp parallel for private(i) 
  PFOR( i = 0 ; i < LCU ; i ++ ) {
    FFTW.psq[i] = MAX_COULOMB / ( gen_p_sq( i , ND - 1 ) ) ;
  }

  // create the fftw plans, or read them if they are stored
  create_plans_DFT( &FFTW , Latt.dims , TRUE_HERM , ND - 1 ) ;

  /////  End of the search for Wisdom  ////
#else
  
  FFTW.in = malloc( ( TRUE_HERM ) * sizeof( GLU_complex* ) ) ; 
  #pragma omp parallel for private(i)
  PFOR(  i = 0 ; i < TRUE_HERM ; i++  ) {
    FFTW.in[i] = ( GLU_complex* )malloc( LCU * sizeof( GLU_complex ) ) ; 
  }
  // these are dummy arrays
  FFTW.out = NULL ;
  FFTW.psq = NULL ;
  FFTW.forward = NULL ;
  FFTW.backward = NULL ;

#endif
  
  const size_t iters =
    #ifdef GLU_GFIX_SD
    Coulomb_FA( lat , &FFTW , accuracy , iter , GLU_FALSE ) ;
    #else
    Coulomb_FA( lat , &FFTW , accuracy , iter , GLU_TRUE ) ;
    #endif
 
#ifdef HAVE_FFTW3_H
    clean_up_fftw( FFTW , TRUE_HERM ) ;
#else
  
  #pragma omp parallel for private(i) 
  PFOR( i = 0 ; i < TRUE_HERM ; i ++  ) {
    free( FFTW.in[i] ) ;
  }
  free( FFTW.in ) ;

#endif

  all_links( lat , &splink , &tlink ) ;
  fprintf( stdout , "[GF] Tuning :: %f || Iterations :: %zu ||"
	   "\n[GF] Final Tlink :: %1.15f "
	   "|| Slink :: %1.15f \n[GF] Plaquette :: %1.15f \n" , 
	   Latt.gf_alpha , iters , tlink , splink , av_plaquette( lat ) ) ; 

  return iters ; 
}

