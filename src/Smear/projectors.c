/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (projectors.c) is part of GLU.

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
   @file projectors.c
   @brief definitions for the projections I use in the smearing and wilson flow rountines ...
 */

#include "Mainfile.h"

#include "plaqs_links.h"  // needed in print_smearing_obs()
#include "givens.h"       // or the standard cabibbo marinari rotations?
#include "taylor_logs.h"  // are we using the unit circle projection?

// this is a common IO pattern for all the smearings
void
print_smearing_obs( const struct site *__restrict lat ,
		    const size_t count )
{
  #pragma omp master
  {
    fprintf( stdout , "[SMEAR] {Iteration} %zu {Link trace} %1.15f \n" , 
	     count , links( lat ) ) ; 
    // print out the info ::
    double splaq , tplaq , plaq = all_plaquettes( lat , &splaq , &tplaq ) ;
    fprintf( stdout , "[SMEAR] {Plaquette} %1.15f {Spatial} %1.15f "
	     "{Temporal} %1.15f \n\n" , plaq , splaq , tplaq ) ; 
  }
  return ;
}

// my version of APE smearing projection full of vitriol, have a good choice
// of projections givens is the trace maximisation, nape_reunit is the projection
// in M-P and durr's paper
void
project_APE( GLU_complex smeared_link[ NCNC ] , 
	     GLU_complex staple[ NCNC ] , 
	     const GLU_complex link[ NCNC ] , 
	     const double smear_alpha , 	     
	     const double al )
{
  // write it all out here ..... ( 1 - alpha )U + lev*alpha/(2*(ND-1))*staple
#if NC == 3
  *( smeared_link + 0 ) = al * link[ 0 ] + smear_alpha * staple[ 0 ] ;
  *( smeared_link + 1 ) = al * link[ 1 ] + smear_alpha * staple[ 1 ] ;
  *( smeared_link + 2 ) = al * link[ 2 ] + smear_alpha * staple[ 2 ] ;
  *( smeared_link + 3 ) = al * link[ 3 ] + smear_alpha * staple[ 3 ] ;
  *( smeared_link + 4 ) = al * link[ 4 ] + smear_alpha * staple[ 4 ] ;
  *( smeared_link + 5 ) = al * link[ 5 ] + smear_alpha * staple[ 5 ] ;
  *( smeared_link + 6 ) = al * link[ 6 ] + smear_alpha * staple[ 6 ] ;
  *( smeared_link + 7 ) = al * link[ 7 ] + smear_alpha * staple[ 7 ] ;
  *( smeared_link + 8 ) = al * link[ 8 ] + smear_alpha * staple[ 8 ] ;
#elif NC == 2
  *( smeared_link + 0 ) = al * link[ 0 ] + smear_alpha * staple[ 0 ] ;
  *( smeared_link + 1 ) = al * link[ 1 ] + smear_alpha * staple[ 1 ] ;
  *( smeared_link + 2 ) = al * link[ 2 ] + smear_alpha * staple[ 2 ] ;
  *( smeared_link + 3 ) = al * link[ 3 ] + smear_alpha * staple[ 3 ] ;
#else
  size_t elem ;
  for( elem = 0 ; elem < NCNC ; elem++ ) {
    smeared_link[ elem ] = al * link[ elem ] + smear_alpha * staple[ elem ] ;  
  }
#endif
#ifdef N_APE
  nape_reunit( smeared_link ) ;
#else
  givens_reunit( smeared_link ) ;
#endif
  return ;
}

// projection for the log smearing .
void
project_LOG( GLU_complex smeared_link[ NCNC ] , 
	     GLU_complex staple[ NCNC ] , 
	     const GLU_complex link[ NCNC ] , 
	     const double smear_alpha ,
	     const double al )
{
  GLU_complex a[ NCNC ] GLUalign ; 
#if NC == 3
  *( staple + 0 ) *= smear_alpha ; 
  *( staple + 1 ) *= smear_alpha ; 
  *( staple + 2 ) *= smear_alpha ;
  *( staple + 3 ) = conj( staple[1] ) ; 
  *( staple + 4 ) *= smear_alpha ; 
  *( staple + 5 ) *= smear_alpha ;
  *( staple + 6 ) = conj( staple[2] ) ; 
  *( staple + 7 ) = conj( staple[5] ) ; 
  *( staple + 8 ) = - staple[0] - staple[4] ;
#elif NC == 2
  *( staple + 0 ) *= smear_alpha ; 
  *( staple + 1 ) *= smear_alpha ;
  *( staple + 2 ) = conj( staple[1] ) ; 
  *( staple + 3 ) = -staple[0] ;
#else
  size_t elem ;
  for( elem = 0 ; elem < NCNC ; elem++ ) {
    *( staple + elem ) *= smear_alpha ;
  }
#endif
  exponentiate( a , staple ) ; 
  multab_suNC( smeared_link , a , link ) ;
  return ;
}

// projection for the Wilson flow with the log
void
project_LOG_wflow_short( GLU_complex log[ NCNC ] , 
			 GLU_complex *__restrict staple , 
			 const GLU_complex link[ NCNC ] , 
			 const double smear_alpha )
{
  GLU_complex a[ HERMSIZE ] , b[ NCNC ] GLUalign ;
  #if NC == 3
  *( a + 0 ) = creal( staple[0] ) * smear_alpha ; 
  *( a + 1 ) = staple[1] * smear_alpha ; 
  *( a + 2 ) = staple[2] * smear_alpha ;
  *( a + 3 ) = cimag( staple[0] ) * smear_alpha ; 
  *( a + 4 ) = staple[3] * smear_alpha ;
  //*( a + 4 ) = staple[4] * smear_alpha ; 
  #elif NC == 2
  *( a + 0 ) = creal( staple[0] ) * smear_alpha ; 
  *( a + 1 ) = staple[1] * smear_alpha ; 
  #else
  size_t elem ;
  for( elem = 0 ; elem < HERMSIZE ; elem++ ) {
    *( a + elem ) = staple[elem] * smear_alpha ;
  }
  #endif
  exponentiate_short( b , a ) ; 
  multab_suNC( log , b , link ) ;
  return ;
}

// projection for stout smearing ..
void
project_STOUT( GLU_complex smeared_link[ NCNC ] , 
	       GLU_complex staple[ NCNC ] , 
	       const GLU_complex link[ NCNC ] , 
	       const double smear_alpha ,
	       const double al )
{
  GLU_complex a[ NCNC ] GLUalign , b[ NCNC ] GLUalign ;  
  multab_dag( b , staple , link ) ;   
  Hermitian_proj( a , b ) ;
#if NC == 3
  *( a + 0 ) *= smear_alpha ; 
  *( a + 1 ) *= smear_alpha ; 
  *( a + 2 ) *= smear_alpha ;
  *( a + 3 ) = conj( a[1] ) ; 
  *( a + 4 ) *= smear_alpha ; 
  *( a + 5 ) *= smear_alpha ;
  *( a + 6 ) = conj( a[2] ) ; 
  *( a + 7 ) = conj( a[5] ) ; 
  *( a + 8 ) = - a[0] - a[4] ;
#elif NC == 2
  *( a + 0 ) *= smear_alpha ; 
  *( a + 1 ) *= smear_alpha ;
  *( a + 2 ) = conj( a[1] ) ; 
  *( a + 3 ) = -a[0] ;
#else
  size_t elem ;
  for( elem = 0 ; elem < NCNC ; elem++ ) {
    *( a + elem ) *= smear_alpha ;
  }
#endif
  exponentiate( b , a ) ; 
  multab_suNC( smeared_link , b , link ) ; 
  return ;
}

// I needed to write a cheaper projection for the wilson flow
void
project_STOUT_wflow_short( GLU_complex stout[ NCNC ] , 
			   GLU_complex *__restrict staple , 
			   const GLU_complex link[ NCNC ] , 
			   const double smear_alpha )
{
  GLU_complex a[ HERMSIZE ] , b[ NCNC ] GLUalign ;
#if NC == 3
  *( a + 0 ) = creal( staple[ 0 ] ) * smear_alpha ; 
  *( a + 1 ) = staple[ 1 ] * smear_alpha ; 
  *( a + 2 ) = staple[ 2 ] * smear_alpha ;
  *( a + 3 ) = cimag( staple[ 0 ] ) * smear_alpha ;
  *( a + 4 ) = staple[ 3 ] * smear_alpha ; 
#elif NC == 2
  *( a + 0 ) = creal( staple[ 0 ] ) * smear_alpha ; 
  *( a + 1 ) = staple[ 1 ] * smear_alpha ; 
#else
  size_t elem ;
  for( elem = 0 ; elem < HERMSIZE ; elem++ ) {
    *( a + elem ) = staple[ elem ] * smear_alpha ;
  }
#endif
  exponentiate_short( b , a ) ; 
  multab_suNC( stout , b , link ) ; 
  return ;
}
