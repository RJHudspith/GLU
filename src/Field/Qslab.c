/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (Qsusc.c) is part of GLU.

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
   @file Qslab.c
   @brief computation of the topological susceptibility on slabs
 */
#include "Mainfile.h"

#include "clover.h"     // computation of the topological charge
#include "cut_output.h" // automatic formatting of our output file
#include "geometry.h"   // for the spacing computation
#include "plan_ffts.h"  // config space correlator is convolution
#include "str_stuff.h"  // append_char()

// local version of the one in geometry
static void
get_mom_2piBZ_slab( int x[ ND ] ,
		    const size_t Dims[ ND ] ,
		    const size_t i ,
		    const size_t DIMS )
{
  int mu , subvol = 1 , sum = 0 ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    if( mu != (int)DIMS ) {
      x[ mu ] = ( ( i - i % subvol ) / subvol ) % Dims[ mu ] ;
      subvol *= Dims[ mu ] ;
      sum += x[mu] ;
    } else {// set it to 0?
      x[ mu ] = 0 ;
    }
  }
  return ;
}

// given a local SLAB site index return the global index taking
// into account for the translational invariance term T0
static size_t
slab_idx( const size_t i ,
	  const size_t slab_dims[ ND ] ,
	  const size_t SLAB_DIR ,
	  const size_t T0 )
{
  int x[ ND ] ;
  get_mom_2piBZ_slab( x , slab_dims , i , ND ) ;
  x[ SLAB_DIR ] = ( x[ SLAB_DIR ] + T0 ) % Latt.dims[ SLAB_DIR ] ;
  return gen_site( x ) ;
}

// compute the slab definition of the topological susceptibility
static int
compute_slabs( const GLU_complex *qtop ,
	       const struct cut_info CUTINFO ,
	       const size_t SLAB_DIR ,
	       const size_t measurement )
{
  // normalisations
  const double NORM = -0.001583143494411527678811 ; // -1.0/(64*Pi*Pi)
  const double NORMSQ = NORM * NORM ;
  
  // various sums and things
  register double sum = 0.0 ;
  size_t T0 , T1 , i , mu ;

  // set up the outputs
  char *str = output_str_struct( CUTINFO ) ;
  char tmp[64] ;
  sprintf( tmp , ".m%zu.tcorr_d%zu" , measurement , SLAB_DIR ) ;
  append_char( &str , tmp ) ;
  FILE *Ap = fopen( str , "wb" ) ;

  // timeslice length
  size_t lt[ 1 ] = { Latt.dims[ SLAB_DIR ] } ;

  // temporal correlator
  double *ct = malloc( Latt.dims[ SLAB_DIR ] * sizeof( double ) ) ;

  // write out the timeslice list ...
  write_tslice_list( Ap , lt ) ;

  // compute the normal topological susceptibility
  for( i = 0 ; i < LVOLUME ; i++ ) {
    sum += creal( qtop[i] ) ;
  }
  fprintf( stdout , "\n[QTOP] Q %zu %1.12e %1.12e \n" ,
	   measurement , sum * NORM , sum * sum * NORMSQ ) ;

  // compute the slab definition
  for( T1 = 1 ; T1 <= Latt.dims[ SLAB_DIR ] ; T1++ ) {

    // set the dimensions of the slab
    size_t dims[ ND ] , subvol = 1 ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      if( mu == SLAB_DIR ) {
	dims[ mu ] = T1 ;
      } else {
	dims[ mu ] = Latt.dims[ mu ] ;
      }
      subvol *= dims[ mu ] ;
    }

#ifdef HAVE_FFTW3_H

    #ifdef verbose
    printf( "DIMS :: %zu %zu %zu %zu -> %zu \n" ,
	    dims[0] , dims[1] , dims[2] , dims[3] , subvol ) ;
    #endif

    struct fftw_small_stuff FFTW ;
    small_create_plans_DFT( &FFTW , dims , ND ) ;
    
#endif

    double tsum = 0.0 ;
    // sum over all possible time cuts
    for( T0 = 0 ; T0 < Latt.dims[ SLAB_DIR ] ; T0++ ) {

      double sum_slab = 0.0 ;
      // perform convolution using FFTW
      #ifdef HAVE_FFTW3_H
      #pragma omp parallel for private(i)
      for( i = 0 ; i < subvol ; i++ ) {
	const size_t src = slab_idx( i , dims , SLAB_DIR , T0 ) ;
	FFTW.in[ i ] = qtop[ src ] ; 
      }
      fftw_execute( FFTW.forward ) ;
      #pragma omp parallel for private(i)
      for( i = 0 ; i < subvol ; i++ ) {
	FFTW.out[ i ] *= conj( FFTW.out[i] ) ;
      }
      fftw_execute( FFTW.backward ) ;

      for( i = 0 ; i < subvol ; i++ ) {
	sum_slab += creal( FFTW.in[ i ] ) ;
      }
      // perform convolution the dumb way
      #else
      #pragma omp parallel for private(i)
      for( i = 0 ; i < subvol ; i++ ) {
	const size_t src = slab_idx( i , dims , SLAB_DIR , T0 ) ;
	const double qsrc = qtop[ src ] ;
	size_t j ;
	for( j = 0 ; j < subvol ; j++ ) {
	  const size_t snk = slab_idx( j , dims , SLAB_DIR , T0 ) ;
	  sum_slab += creal ( qsrc * qtop[ snk ] ) ;
	}
      }
      #endif
      tsum += sum_slab ;
    }

#ifdef HAVE_FFTW3_H
    
    // do the usual convolution norm
    tsum /= ( subvol ) ;

    // cleanup and memory deallocate
    small_clean_up_fftw( FFTW ) ;
#endif

    ct[ T1-1 ] = tsum * NORMSQ / Latt.dims[ SLAB_DIR ] ;
  }

  // write out the list
  for( T1 = 0 ; T1 < Latt.dims[ SLAB_DIR ] ; T1++ ) {
    printf( "[QSUSC] SLAB_%zu %zu %1.12e \n" , SLAB_DIR , T1+1 , ct[T1] ) ;
  }

  // tell us where to go
  fprintf( stdout , "[CUTS] Outputting correlator to %s \n" , str ) ;

  // and write the props ....
  write_g2_to_list( Ap , ct , lt ) ;

  // free allocated memory
  free( str ) ;
  fclose( Ap ) ;

  return GLU_SUCCESS ;
}

// compute the topological susceptibility on a "slab"
int
compute_Qsusc( struct site *lat ,
	       const struct cut_info CUTINFO ,
	       const size_t measurement )
{
  // set up the matrix-valued array of the topological charge
  GLU_complex *qtop = malloc( LVOLUME * sizeof( GLU_complex ) ) ;

  // precompute all of charge densities q(x) over whole lattice
  compute_Gmunu_array( qtop , lat ) ;

  // slabby slab slab
  size_t mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    compute_slabs( qtop , CUTINFO , mu , measurement ) ;
  }
  
  free( qtop ) ;
  
  return GLU_SUCCESS ;
}
