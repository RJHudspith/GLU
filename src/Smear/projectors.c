/*
    Copyright 2013 Renwick James Hudspith

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
#ifdef FAST_SMEAR         // super dangerous, should probably remove this option
 #include "gramschmidt.h" // used for the reunitarisation step in APE
                          // this is so unbearably dangerous
#else
 #include "givens.h"      // or the standard cabibbo marinari rotations?
 #include "taylor_logs.h" // are we using the unit circle projection?
#endif

// this is a common IO pattern for all the smearings
void
print_smearing_obs( const struct site *__restrict lat , 
		    const int type ,
		    const int count ,
		    const GLU_bool hypercubically_blocked )
{
  printf( "[SMEAR] {Iteration} %d {Link trace} %1.15f \n" , count , links( lat ) ) ; 
  // print out the info ::
  double splaq , tplaq , plaq = all_plaquettes( lat , &splaq , &tplaq ) ;
  printf( "[SMEAR] {Plaquette} %1.15f {Spatial} %1.15f {Temporal} %1.15f \n\n" , 
	  plaq , splaq , tplaq ) ; 
  return ;
}

// my version of APE smearing projection full of vitriol, have a good choice
// of projections givens is the trace maximisation, nape_reunit is the projection
// in M-P and durr's paper and simple, stupid reunit if going fast
void
project_APE( GLU_complex ape[ NCNC ] , 
	     const GLU_complex staple[ NCNC ] , 
	     const GLU_complex link[ NCNC ] , 
	     const double smear_alpha , 	     
	     const double al )
{
  // write it all out here ..... ( 1 - alpha )U + lev*alpha/(2*(ND-1))*staple
  #if NC == 3
  *( ape + 0 ) = al * link[ 0 ] + smear_alpha * staple[ 0 ] ;
  *( ape + 1 ) = al * link[ 1 ] + smear_alpha * staple[ 1 ] ;
  *( ape + 2 ) = al * link[ 2 ] + smear_alpha * staple[ 2 ] ;
  *( ape + 3 ) = al * link[ 3 ] + smear_alpha * staple[ 3 ] ;
  *( ape + 4 ) = al * link[ 4 ] + smear_alpha * staple[ 4 ] ;
  *( ape + 5 ) = al * link[ 5 ] + smear_alpha * staple[ 5 ] ;
  *( ape + 6 ) = al * link[ 6 ] + smear_alpha * staple[ 6 ] ;
  *( ape + 7 ) = al * link[ 7 ] + smear_alpha * staple[ 7 ] ;
  *( ape + 8 ) = al * link[ 8 ] + smear_alpha * staple[ 8 ] ;
  #elif NC == 2
  *( ape + 0 ) = al * link[ 0 ] + smear_alpha * staple[ 0 ] ;
  *( ape + 1 ) = al * link[ 1 ] + smear_alpha * staple[ 1 ] ;
  *( ape + 2 ) = al * link[ 2 ] + smear_alpha * staple[ 2 ] ;
  *( ape + 3 ) = al * link[ 3 ] + smear_alpha * staple[ 3 ] ;
  #else
  int elem ;
  for( elem = 0 ; elem < NCNC ; elem++ ) {
    ape[ elem ] = al * link[ elem ] + smear_alpha * staple[ elem ] ;  
  }
  #endif  

  #ifndef FAST_SMEAR
  // ugh, totally included the su3 projection using some sort of givens rotation, feel dirty
    #ifdef N_APE
    nape_reunit( ape ) ;
    #else
    givens_reunit( ape ) ;
    #endif
  #else
  reunit2( ape ) ;
  #endif

  return ;
}

// projection for stout smearing ..
void
project_STOUT( GLU_complex stout[ NCNC ] , 
	       const GLU_complex staple[ NCNC ] , 
	       const GLU_complex link[ NCNC ] , 
	       const double smear_alpha )
{
  GLU_complex b[ NCNC ] ;  
  // staple is no longer in SU(NC) cannot use speedy versions
  multab_dag( b , staple , link ) ;   
  GLU_complex a[ NCNC ] ;
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
  int elem ;
  for( elem = 0 ; elem < NCNC ; elem++ ) {
    *( a + elem ) *= smear_alpha ;
  }
  #endif
  //approx exp is the taylor expanded version of the exponentiate
  #ifdef FAST_SMEAR
  approx_exp( b , a ) ; 
  #else 
  exponentiate( b , a ) ; 
  #endif
  multab_suNC( stout , b , link ) ; 
  return ;
}

// projection for stout smearing ..
void
project_STOUT_short( GLU_complex stout[ NCNC ] , 
		     const GLU_complex staple[ NCNC ] , 
		     const GLU_complex link[ NCNC ] , 
		     const double smear_alpha )
{
#if NC > 3
  return project_STOUT( stout , staple , link , smear_alpha ) ;
#else
  GLU_complex b[ NCNC ] ;  
  // staple is no longer in SU(NC) cannot use speedy versions
  multab_dag( b , staple , link ) ;   
  GLU_complex a[ HERMSIZE ] ;
  Hermitian_proj_short( a , b ) ;
#if NC == 3
  *( a + 0 ) *= smear_alpha ; 
  *( a + 1 ) *= smear_alpha ; 
  *( a + 2 ) *= smear_alpha ;
  *( a + 3 ) *= smear_alpha ; 
  *( a + 4 ) *= smear_alpha ;
#elif NC == 2
  *( a + 0 ) *= smear_alpha ; 
  *( a + 1 ) *= smear_alpha ;
#else
  int mu ;
  for( mu = 0 ; mu < HERMSIZE ; mu++ ) {
    a[ mu ] *= smear_alpha ;
  }
#endif
  exponentiate_short( b , a ) ; 
  multab_suNC( stout , b , link ) ; 
  return ;
#endif
}

// STOUT projection specifically for the Wilson flow
void
project_STOUT_wflow( GLU_complex stout[ NCNC ] , 
		     const GLU_complex *__restrict staple , 
		     const GLU_complex link[ NCNC ] , 
		     const double smear_alpha )
{
  GLU_complex a[ NCNC ] ;
#if NC == 3
  *( a + 0 ) = staple[0] * smear_alpha ; 
  *( a + 1 ) = staple[1] * smear_alpha ; 
  *( a + 2 ) = staple[2] * smear_alpha ;
  *( a + 3 ) = conj( a[1] ) ;
  *( a + 4 ) = staple[3] * smear_alpha ; 
  *( a + 5 ) = staple[4] * smear_alpha ;
  *( a + 6 ) = conj( a[2] ) ; 
  *( a + 7 ) = conj( a[5] ) ; 
  *( a + 8 ) = -a[0] - a[4] ;
#elif NC == 2
  *( a + 0 ) = staple[0] * smear_alpha ; 
  *( a + 1 ) = staple[1] * smear_alpha ; 
  *( a + 2 ) = conj( a[1] ) ; 
  *( a + 3 ) = -a[0] ;
#else
  int elem ;
  rebuild_hermitian( a , staple ) ;
  for( elem = 0 ; elem < NCNC ; elem++ ) {
    *( a + elem ) *= smear_alpha ;
  }
#endif
  //approx exp is the taylor expanded version of the exponentiate
  GLU_complex b[ NCNC ] ;
  #ifdef FAST_SMEAR
     approx_exp( b , a ) ; 
  #else
     exponentiate( b , a ) ; 
  #endif
  multab_suNC( stout , b , link ) ; 
  return ;
}

// I needed to write a cheaper projection for the wilson flow
void
project_STOUT_wflow_short( GLU_complex stout[ NCNC ] , 
			   GLU_complex *__restrict staple , 
			   const GLU_complex link[ NCNC ] , 
			   const double smear_alpha )
{
  GLU_complex a[ HERMSIZE ] ;
#if NC == 3
  *( a + 0 ) = creal( staple[ 0 ] ) * smear_alpha ; 
  *( a + 1 ) = staple[ 1 ] * smear_alpha ; 
  *( a + 2 ) = staple[ 2 ] * smear_alpha ;
  *( a + 3 ) = cimag( staple[ 0 ] ) * smear_alpha ;
  *( a + 4 ) = staple[ 3 ] * smear_alpha ; 
  //*( a + 4 ) = staple[ 4 ] * smear_alpha ;
#elif NC == 2
  *( a + 0 ) = creal( staple[ 0 ] ) * smear_alpha ; 
  *( a + 1 ) = staple[ 1 ] * smear_alpha ; 
#else
  int elem ;
  for( elem = 0 ; elem < HERMSIZE ; elem++ ) {
    *( a + elem ) = staple[ elem ] * smear_alpha ;
  }
#endif
  //approx exp is the taylor expanded version of the exponentiate
  GLU_complex b[ NCNC ] ;
  #ifdef FAST_SMEAR
     approx_exp_short( b , a ) ; 
  #else
     exponentiate_short( b , a ) ; 
  #endif
  multab_suNC( stout , b , link ) ; 
  return ;
}

// projection for the log smearing .
void
project_LOG( GLU_complex log[ NCNC ] , 
	     GLU_complex staple[ NCNC ] , 
	     const GLU_complex link[ NCNC ] , 
	     const double smear_alpha )
{
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
  int elem ;
  for( elem = 0 ; elem < NCNC ; elem++ ) {
    *( staple + elem ) *= smear_alpha ;
  }
#endif
  GLU_complex a[ NCNC ] ; 
  //can't use approx exp 
  exponentiate( a , staple ) ; 
  multab_suNC( log , a , link ) ;
  return ;
}

// projection for the log smearing short def
void
project_LOG_short( GLU_complex log[ NCNC ] , 
		   GLU_complex staple[ NCNC ] , 
		   const GLU_complex link[ NCNC ] , 
		   const double smear_alpha )
{
#if NC > 3
  return project_LOG( log , staple , link , smear_alpha ) ;
#else
  #if NC == 3
  *( staple + 0 ) *= smear_alpha ; 
  *( staple + 1 ) *= smear_alpha ; 
  *( staple + 2 ) *= smear_alpha ;
  *( staple + 3 ) = smear_alpha * staple[ 4 ]; 
  *( staple + 4 ) = smear_alpha * staple[ 5 ] ;
  #elif NC == 2
  *( staple + 0 ) *= smear_alpha ; 
  *( staple + 1 ) *= smear_alpha ;
  #else
  int mu ;
  for( mu = 0 ; mu < HERMSIZE ; mu++ ) {
    *( staple + mu ) *= smear_alpha ;
  }
  #endif
  GLU_complex a[ NCNC ] ; 
  //can't use approx exp 
  exponentiate_short( a , staple ) ; 
  multab_suNC( log , a , link ) ;
  return ;
#endif
}

// LOG smearing wilson flow ....
void
project_LOG_wflow( GLU_complex log[ NCNC ] , 
		   GLU_complex *__restrict staple , 
		   const GLU_complex link[ NCNC ] , 
		   const double smear_alpha )
{
  GLU_complex a[ NCNC ] ;
  #if NC == 3
  *( a + 0 ) = staple[0] * smear_alpha ; 
  *( a + 1 ) = staple[1] * smear_alpha ; 
  *( a + 2 ) = staple[2] * smear_alpha ;
  *( a + 3 ) = conj( a[1] ) ;
  *( a + 4 ) = staple[3] * smear_alpha ; 
  *( a + 5 ) = staple[4] * smear_alpha ;
  *( a + 6 ) = conj( a[2] ) ; 
  *( a + 7 ) = conj( a[5] ) ; 
  *( a + 8 ) = -a[0] - a[4] ;
  #elif NC == 2
  *( a + 0 ) = staple[0] * smear_alpha ; 
  *( a + 1 ) = staple[1] * smear_alpha ; 
  *( a + 2 ) = conj( a[1] ) ; 
  *( a + 3 ) = -a[0] ;
  #else
  int elem ;
  rebuild_hermitian( a , staple ) ;
  for( elem = 0 ; elem < NCNC ; elem++ ) {
    *( a + elem ) *= smear_alpha ;
  }
  #endif
  GLU_complex b[ NCNC ] ;
  //can't use approx exp 
  exponentiate( b , a ) ; 
  multab_suNC( log , b , link ) ;
  return ;
}

// projection for the Wilson flow with the log
void
project_LOG_wflow_short( GLU_complex log[ NCNC ] , 
			 GLU_complex *__restrict staple , 
			 const GLU_complex link[ NCNC ] , 
			 const double smear_alpha )
{
  GLU_complex a[HERMSIZE] ;
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
  int elem ;
  for( elem = 0 ; elem < HERMSIZE ; elem++ ) {
    *( a + elem ) = staple[elem] * smear_alpha ;
  }
  #endif
  GLU_complex b[ NCNC ] ;
  //can't use approx exp 
  exponentiate_short( b , a ) ; 
  multab_suNC( log , b , link ) ;
  return ;
}


