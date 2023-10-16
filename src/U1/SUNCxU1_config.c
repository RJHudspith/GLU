/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (su3xu1_config.c) is part of GLU.

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
   @file SUNCxU1_config.c
   @brief U(1)-ifies the link matrices

   @warning requires FFTW linking to work
 */
#define _GNU_SOURCE // sincos

#include "Mainfile.h"

// and the other headers it uses
#include "geometry.h"
#include "par_rng.h"
#include "GLU_timer.h"
#include "plan_ffts.h"
#include "U1_obs.h"

#ifdef HAVE_FFTW3_H

// just to make it clear what we are doing
enum{ CONJUGATE_NOT_IN_LIST , CONJUGATE_IN_LIST } ;

// little inline for the ( 0 , 0 , .. , 0 ) point in the -Pi -> Pi BZ
static size_t
compute_zeropoint( void )
{
#if ND == 4
  return ( Latt.dims[0] * ( 1 + Latt.dims[1] * ( 1 + Latt.dims[2] * ( 1 + Latt.dims[3] ) ) ) ) >> 1 ;
#else
  size_t check = Latt.dims[ ND-1 ] , mu ;
  for( mu = ND-2 ; mu != 0 ; mu-- ) {
    check = Latt.dims[ mu ] * ( 1 + check ) ;
  }
  return check >> 1 ;
#endif
}

// compute the conjugate site to "i"
static size_t
conjugate_site( const size_t i )
{
  int x[ ND ] ;
  get_mom_pipi( x , i , ND ) ;
  #if ND == 4
  x[0] = -x[0] ;
  x[1] = -x[1] ;
  x[2] = -x[2] ;
  x[3] = -x[3] ;
  #else
  size_t mu ;
  for( mu = 0 ; mu < ND ; mu++ ) { x[mu] = -x[mu] ; }
  #endif
  //  contrary to what its name may suggest get_site_2piBZ translates
  //  the -pi to pi BZ coordinates to a site in the 0->2Pi lattice
  //  which is what the FFT uses
  return get_site_2piBZ( x , ND ) ;
}

// periodic fields using the DFT
static void 
periodic_dft( struct fftw_stuff *FFTW )
{
  // slower version just randomly assigns fields first
  size_t i ;
  #pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    const uint32_t thread = get_GLU_thread( ) ;
    #if ND == 4
    FFTW -> in[0][i] = FFTW -> psq[i] * par_polar_box( thread ) ;
    FFTW -> in[1][i] = FFTW -> psq[i] * par_polar_box( thread ) ;
    FFTW -> in[2][i] = FFTW -> psq[i] * par_polar_box( thread ) ;
    FFTW -> in[3][i] = FFTW -> psq[i] * par_polar_box( thread ) ;
    #else
    for( size_t mu = 0 ; mu < ND ; mu++ ) {
      FFTW -> in[mu][i] = FFTW -> psq[i] * par_polar_box( thread ) ;
    }
    #endif
  }  
  const size_t SYMM_POINT = compute_zeropoint( ) + 1 ;
  for( i = 0 ; i < SYMM_POINT ; i++ ) {
    const size_t b = conjugate_site( i ) ;
    if( i == b ) {
      for( size_t mu = 0 ; mu < ND ; mu++ ) {
	FFTW -> in[mu][i] = creal( FFTW -> in[mu][i] ) + I*0.0 ;
      }
    } else if( b < i ) { // case is already done
      continue ;
    } else {
      for( size_t mu = 0 ; mu < ND ; mu++ ) {
	FFTW -> in[mu][b] = conj( FFTW -> in[mu][i] ) ; 
      }      
    }
  }
  return ;
}

// create the U1 fields and multiplies the gauge field
static int
create_u1( struct site *lat ,
	   struct fftw_stuff *FFTW ,
	   const struct u1_info U1INFO )
{
#ifndef HAVE_FFTW3_H
  fprintf( stderr , "[U(1)] Cannot create U(1) fields, "
	   "requires FFTW ... Leaving\n" ) ;
  return GLU_FAILURE ;
#else
  // alpha = 1/4\pi = 0.0795775387 is beta = 1 , gives noncompact plaq = 0.5
  // note the factor of 2 in the plaquette due to my conventions of giving
  // compact plaquette at \alpha = 0 of 1
  // factor of Volume comes from us not normalizing the FFTs
  const GLU_real Nbeta = LVOLUME / ( 4. * MPI * U1INFO.alpha ) ;
  fprintf( stdout , "\n[U(1)] Beta :: %f \n\n" , 1./( 4*MPI * U1INFO.alpha ) ) ;

  size_t i ;

  fprintf( stdout , "[U(1)] Using the DFT U1 code ... \n\n" ) ;

  // set up the distribution with the lattice mom Feynmann gauge QEDL prescription
  periodic_dft( FFTW ) ;
  
  // fft the fields allow for the parallel omp-ified fftws
  const GLU_real rbeta = 1.0 / (GLU_real)sqrt( Nbeta ) ;
#pragma omp parallel
  {
    #pragma omp for private(i) 
    for( i = 0 ; i < ND ; i++ ) {
      fftw_execute( FFTW -> forward[ i ] ) ;
    }
    #pragma omp for private( i )
    for( i = 0 ; i < LVOLUME*ND ; i++ ) {
      const size_t idx = i/ND ;
      const size_t mu = i%ND ;
      FFTW->out[mu][idx] = creal( FFTW->out[mu][idx] )*rbeta + I*0.0 ;

      // call to cexp is expensive so we just do a simple sincos
      double re , im ;
      sincos( creal( FFTW -> out[mu][idx])*U1INFO.charge , &im , &re ) ;
      // could be a vectorised multiply ... 
      register const GLU_complex U1 = re + I*im ;
      #if NC == 3
      lat[idx].O[mu][0] *= U1 ;
      lat[idx].O[mu][1] *= U1 ;
      lat[idx].O[mu][2] *= U1 ;
      lat[idx].O[mu][3] *= U1 ;
      lat[idx].O[mu][4] *= U1 ;
      lat[idx].O[mu][5] *= U1 ;
      lat[idx].O[mu][6] *= U1 ;
      lat[idx].O[mu][7] *= U1 ;
      lat[idx].O[mu][8] *= U1 ;
      #else
      size_t element ;
      for( element = 0 ; element < NCNC ; element++ ) {
	lat[idx].O[mu][element] *= U1 ;
      }
      #endif
    }
  }
  fprintf( stdout , "[U1] fields set" ) ;
  print_time( ) ;

#endif
  return GLU_SUCCESS ;
}
#endif // HAVE_FFTW3_H

// Possibly a check that the plaquette has changed ..
int
suNC_cross_u1( struct site *lat , 
	       const struct u1_info U1INFO )
{
#ifndef HAVE_FFTW3_H
  fprintf( stderr , "[U1] Require FFTW to be linked to do quenched U(1)\n" ) ;
  return GLU_FAILURE ;
#else

  // begin timing
  start_timer( ) ;

  // initialise the rng
  initialise_par_rng( NULL ) ;

  struct fftw_stuff FFTW ;
  create_plans_DFT( &FFTW , Latt.dims , ND , ND ) ;

  // set psq to the QEDL prescription
  size_t i ; 
  FFTW.psq = malloc( LVOLUME*sizeof( GLU_real ) ) ;
  for( i = 0 ; i < LVOLUME ; i++ ) {
    int flag = 0 ;
    register const GLU_real f = 1.0 / (GLU_real)sqrt( gen_p_sq_feyn( i , &flag ) ) ;
    if( flag == 1 ) {
      FFTW.psq[i] = 0 ;
    } else {
      FFTW.psq[i] = f ;
    }
  }

  // create the U1 field ...
  create_u1( lat , &FFTW, U1INFO ) ;

  fprintf( stdout , "[U1] SUNCxU1 created ....\n" ) ;
  
  // compute some U1 observables why not ?
  start_timer() ;
  compute_U1_obs( (const GLU_complex**)FFTW.out , U1INFO.meas ) ;
  print_time() ;

  fprintf( stdout , "[U1] U1 Memory Cleanup ....\n" ) ;
  clean_up_fftw( FFTW , ND ) ;
#endif
  return GLU_SUCCESS ;
}
