/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (glueprop.c) is part of GLU.

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
   @file glueprop.c
   @brief computes the longitudinal and transverse scalars

   writes the transverse scalar G(p^2) and the Longitudinal F(p^2) to a file

   @warning uses momentum list p,-p symmetry
 */

#include "Mainfile.h"

#include "cut_output.h" // output file formatting
#include "geometry.h"   // lexicographical geometry look-ups

// compute psq using the momentum lattice coordinates
static inline double
psq_calc( double mom[ ND ] ,
	  const size_t posit )
{
  size_t mu ;
  int k[ ND ] ;
  TwoPI_mpipi_momconv( k , posit , ND ) ;
  register double tspsq = 0. ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    #ifdef SIN_MOM
    mom[ mu ] = 2.0 * sin( k[ mu ] * Latt.twiddles[ mu ] * 0.5 ) ;
    #else
    mom[ mu ] = k[ mu ] * Latt.twiddles[ mu ] ;
    #endif
    tspsq += mom[ mu ] * mom[ mu ] ;
  }
  return ( tspsq < PREC_TOL ) ? 1.0 : tspsq ;
}

/**
   @fn static void project_trans_long( double *transverse , double *longitudinal , const struct site *__restrict A , const int posit ) 
  Transverse :
                
                  p_\mu p_\nu
  ( g_{\mu\nu} -  ----------- ) Tr[ G_{\mu\nu}(p^2) ] = G(p^2)
                      p^2

  Longitudinal :

  p_\mu p_\nu
  ----------- Tr[ G_{\mu\nu}(p^2) ] = F(p^2)
      p^2
 */
static void
project_trans_long( double *transverse ,
		    double *longitudinal ,
		    const struct site *__restrict A ,
		    const size_t posit ) 
{
  GLU_complex tr ;
  double mom[ ND ] , common ;
  // make this a constant
  const double spsq = 1.0 / psq_calc( mom , posit ) ;
  size_t mu , nu ;
  *transverse = *longitudinal = 0.0 ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    for( nu = 0 ; nu < ND ; nu++ ) {
      trace_ab_dag( &tr , A[posit].O[mu] , A[posit].O[nu] ) ;
      // common factor
      common = mom[mu] * mom[nu] * spsq * (double)creal(tr) ;
      // compute the longitudinal and transverse scalars
      *transverse += ( mu != nu ) ? -common : (double)creal(tr) - common ;
      *longitudinal += common ;
    }
  }
  return ;
}

// does what it says
int
compute_gluon_prop( FILE *__restrict Ap , 
		    const struct site *__restrict A ,
		    const struct veclist *__restrict list ,
		    int num_mom[1] )
{
  fprintf( stdout , "\n[CUTS] computing the transverse and "
	            "longitudinal gluon propagators \n" ) ;

  // normalisations
  const double g2_norm   = 2.0 / ( ( NCNC - 1 ) * ( ND - 1 ) * LVOLUME ) ;
  const double g0_norm   = 2.0 / ( ( NCNC - 1 ) * ( ND ) * LVOLUME ) ;
  const double long_norm = 2.0 / ( ( NCNC - 1 ) * LVOLUME ) ;

  // allocate the props
  double *transverse = malloc( num_mom[0] * sizeof( double ) ) ;
  double *longitudinal = malloc( num_mom[0] * sizeof( double ) ) ;

  // set up the zero position
  const int zero_pos = ( num_mom[0] ) >> 1 ;
  const int zero_mom = list[ zero_pos ].idx ; // ist das wahr?

  // obviously only need to compute num_mom[0]/2 of these due to symmetry
  size_t i ;
  #pragma omp parallel for private(i)
  PFOR( i = 0 ; i < zero_pos ; i++ ) {
    // uses p of list -p contrib IS the hermitian conjugate.
    const int posit = list[i].idx ;
    
    // projections
    project_trans_long( &transverse[i] , &longitudinal[i] , A , posit ) ;

    // normalisations
    transverse[i] *= g2_norm ;
    longitudinal[i] *= long_norm ;
  }

  // contract the zero
  project_trans_long( &transverse[zero_pos] , 
		      &longitudinal[zero_pos] ,
		      A , zero_mom ) ;
  transverse[zero_pos] *= g0_norm ;
  longitudinal[zero_pos] *= long_norm ;

  // fill up the rest using p, -p symmetry
  #pragma omp parallel for private(i)
  PFOR( i = zero_pos+1 ; i < num_mom[0] ; i++ ) {
    transverse[i] = transverse[ num_mom[0] - i - 1 ] ;
    longitudinal[i] = longitudinal[ num_mom[0] - i - 1 ] ;
  }

  // and write out to a file ...
  write_g2g3_to_list( Ap , transverse , longitudinal , num_mom ) ;

  // free the data
  free( transverse ) ;
  free( longitudinal ) ;

  return GLU_SUCCESS ;
}
