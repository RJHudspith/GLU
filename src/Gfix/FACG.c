/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (FACG.c) is part of GLU.

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
   @file FACG.c
   @brief conjugate gradient Fourier Accelerated Landau gauge fixing routines

   Is slightly different from the routines in sweet_FA.c although I could merge them. At the moment only set up for SU(3)
 */
#include "Mainfile.h"    // for all the definitions

#include "CG.h"          // routines used by both CG codes
#include "CGalloc.h"     // CG temporary allocations
#include "GLU_sums.h"    // round off resistant summations
#include "gftests.h"     // theta_test stopping condition
#include "gtrans.h"      // gauge transformations
#include "lin_derivs.h"  // linear approximation of lie matrices
#include "line_search.h" // evaluate alpha
#include "log_derivs.h"  // log-def of lie matrices 
#include "plaqs_links.h" // plaquette routine might be called

static int
FA_deriv( double *red ,
	  GLU_complex **in , 
	  struct site *lat )
{
  size_t i ;
  #pragma omp for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {

    GLU_complex sum[ HERMSIZE ] GLUalign ;
    memset( &sum , 0 , HERMSIZE * sizeof( GLU_complex ) ) ;

    const size_t th = get_GLU_thread() ;
    red[ LINE_NSTEPS + 2 + th*CLINE ] +=
      fast_deriv_AntiHermitian_proj( sum , lat , i ) ;    

    // make in anti-hermitian here!
    #if NC == 3
    // for SU(3) I pack the first FFT with the two explicitly real diagonal elements
    // to save on a Fourier transform ..
    in[0][i] = I * creal( sum[0] ) - creal( sum[3] ) ;
    in[1][i] = I * creal( sum[1] ) - cimag( sum[1] ) ;
    in[2][i] = I * creal( sum[2] ) - cimag( sum[2] ) ;
    in[3][i] = I * creal( sum[4] ) - cimag( sum[4] ) ;
    #elif NC == 2
    in[0][i] = I * creal( sum[0] ) ;
    in[1][i] = I * creal( sum[1] ) - cimag( sum[1] ) ;
    #else
    size_t mu ;
    for( mu = 0 ; mu < HERMSIZE ; mu++ ) {
      in[mu][i] = I * creal( sum[mu] ) - cimag( sum[mu] ) ;
    }
    #endif
  }
  
  return GLU_SUCCESS ;
}

// is the same for Landau and Coulomb just with different LENGTHS
void
FOURIER_ACCELERATE3( struct fftw_stuff *FFTW ) 
{
#ifdef HAVE_FFTW3_H
  const fftw_plan *forw = ( const fftw_plan* )FFTW -> forward ;
  const fftw_plan *back = ( const fftw_plan* )FFTW -> backward ;
  // single core FFT's
  size_t mu ;
  #pragma omp for private(mu) schedule(dynamic)
  for( mu = 0 ; mu < TRUE_HERM ; mu++ ) {
    fftw_execute( forw[mu] ) ; 
    size_t i ;
    for( i = 0 ; i < LVOLUME ; i++ ) {
      FFTW -> out[ mu ][ i ] *= FFTW -> psq[i] ;
    }
    fftw_execute( back[mu] ) ; 
  }
#endif
  return ;
}

// what does this do
static void
sum_PR3( double *red ,
	 const GLU_complex **in , 
	 const GLU_complex **in_old )
{
  size_t i ;
#pragma omp for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    double loc_sum1 = 0.0 , loc_sum2 = 0.0 ;
    #if NC==3
    loc_sum1 += creal(in[0][i]) * creal(in[0][i]) + cimag(in[0][i])*cimag(in[0][i]) 
      + creal(in[0][i]) * cimag( in[0][i] ) ;
    loc_sum1 += creal(in[1][i]) * creal(in[1][i]) + cimag(in[1][i])*cimag(in[1][i]) ; 
    loc_sum1 += creal(in[2][i]) * creal(in[2][i]) + cimag(in[2][i])*cimag(in[2][i]) ; 
    loc_sum1 += creal(in[3][i]) * creal(in[3][i]) + cimag(in[3][i])*cimag(in[3][i]) ; 
	
    register GLU_complex temp = in[0][i] - in_old[0][i] ;
    loc_sum2 += 2.0 * ( creal( in[0][i] ) * creal( temp ) + cimag( in[0][i] ) * cimag( temp ) ) ;
    loc_sum2 += creal( in[0][i] ) * cimag( temp ) + cimag( in[0][i] ) * creal( temp ) ;
    temp = in[1][i] - in_old[1][i] ;
    loc_sum2 += 2.0 * ( creal( in[1][i] ) * creal( temp ) + cimag( in[1][i] ) * cimag( temp ) ) ;
    temp = in[2][i] - in_old[2][i] ;
    loc_sum2 += 2.0 * ( creal( in[2][i] ) * creal( temp ) + cimag( in[2][i] ) * cimag( temp ) ) ;
    temp = in[3][i] - in_old[3][i] ;
    loc_sum2 += 2.0 * ( creal( in[3][i] ) * creal( temp ) + cimag( in[3][i] ) * cimag( temp ) ) ;
    #elif (NC==2)
    loc_sum1 += creal(in[0][i]) * creal(in[0][i]) + cimag(in[0][i])*cimag(in[0][i]) + creal(in[0][i]) * cimag( in[0][i] ) ;    
    loc_sum1 += creal(in[1][i]) * creal(in[1][i]) + cimag(in[1][i])*cimag(in[1][i]) ;
    register GLU_complex temp = in[0][i] - in_old[0][i] ;
    loc_sum2 += 2.0 * ( creal( in[0][i] ) * creal( temp ) + cimag( in[0][i] ) * cimag( temp ) ) ;
    loc_sum2 += creal( in[0][i] ) * cimag( temp ) + cimag( in[0][i] ) * creal( temp ) ;
    temp = in[1][i] - in_old[1][i] ;
    loc_sum2 += 2.0 * ( creal( in[1][i] ) * creal( temp ) + cimag( in[1][i] ) * cimag( temp ) ) ;    
    #else
    size_t mu ;
    for( mu = 0 ; mu < TRUE_HERM ; mu++ ) {
      test_sum1 += conj(in[mu][i])*in[mu][i] ;
      register GLU_complex temp = in[mu][i] - in_old[mu][i] ;
      test_sum2 += 2*creal( conj(in[mu][i])*temp ) ;
    }
    #endif
    const size_t th = get_GLU_thread() ;
    red[ LINE_NSTEPS + th * CLINE ] += 2*loc_sum1 ;
    red[ LINE_NSTEPS + 1 + th * CLINE ] += loc_sum2 ;
  }
  return ;
}

static void
sum_DER3( double *red ,
	  const GLU_complex **in )
{
  size_t i ;
#pragma omp for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    register double loc_sum = 0.0 ;
#if NC == 3
    loc_sum += creal(in[0][i]) * creal(in[0][i]) + cimag(in[0][i])*cimag(in[0][i]) 
      + creal(in[0][i]) * cimag( in[0][i] ) ;
    loc_sum += creal(in[1][i]) * creal(in[1][i]) + cimag(in[1][i])*cimag(in[1][i]) ; 
    loc_sum += creal(in[2][i]) * creal(in[2][i]) + cimag(in[2][i])*cimag(in[2][i]) ; 
    loc_sum += creal(in[3][i]) * creal(in[3][i]) + cimag(in[3][i])*cimag(in[3][i]) ; 
#elif NC == 2
    loc_sum += creal(in[0][i]) * creal(in[0][i]) + cimag(in[0][i])*cimag(in[0][i]) ;
    loc_sum += creal(in[1][i]) * creal(in[1][i]) + cimag(in[1][i])*cimag(in[1][i]) ;
#else
    GLU_complex temp[ HERMSIZE ] GLUalign ;
    size_t mu ;
    for( mu = 0 ; mu < HERMSIZE ; mu++ ) {
      temp[mu] = in[mu][i] ;
    }
    GLU_real tr ;
    trace_ab_herm_short( &tr , temp , temp ) ;
    loc_sum = 0.5 * (double)tr ;
#endif
    const size_t th = get_GLU_thread() ;
    red[ LINE_NSTEPS + th*CLINE ] += 2*loc_sum ;
  }
  return ;
}

static void
steep_Landau_FA( double *red ,
		 struct site *lat ,
		 struct fftw_stuff *FFTW )
{
  // do a steepest-descents step with the result in "out"
  FA_deriv( red , FFTW -> in , lat ) ;
  
  // and do the fourier acceleration
  FOURIER_ACCELERATE3( FFTW ) ;
  
  // and step length gf_alpha provided from the input file
#if (defined GLU_FR_CG) || (defined GLU_GFIX_SD)
  egauge_Landau( lat , (const GLU_complex**)FFTW -> in , Latt.gf_alpha ) ;
#else
  line_search_Landau( red , lat , (const GLU_complex**)FFTW -> in ) ;
#endif
  
  return ;
}

static int
steep_Landau_FASD( struct site *lat ,
		   struct fftw_stuff *FFTW ,
		   struct CGtemps *CG ,
		   double *tr ,
		   const double acc ,
		   const size_t max_iters )
{
  size_t k , loc_iters = 0 ;
  
 top:
  
  {
    #pragma omp barrier
  }
  
  #pragma omp for private(k)
  for( k = 0 ; k < CLINE*Latt.Nthreads ; k++ ) {
    CG -> red[k] = 0.0 ;
  }
  
  steep_Landau_FA( CG -> red , lat , FFTW ) ;

  double trAA = 0.0 ;
  for( k = 0 ; k < Latt.Nthreads ; k++ ) {
    trAA += CG -> red[ LINE_NSTEPS + 2 + k*CLINE ] ;
  }
  *tr = trAA * GFNORM_LANDAU ;

  loc_iters++ ;
  
  if( ( *tr > acc ) && ( loc_iters < max_iters ) ) goto top ;

  return loc_iters ;
}

static int
steep_Landau_FACG( struct site *lat ,
		   struct fftw_stuff *FFTW ,
		   struct CGtemps *CG ,
		   double *tr ,
		   const double acc ,
		   const size_t max_iters )
{  
  size_t k , loc_iters = 0  ;
  double trAA = 0.0 , inold = 0.0 , insum = 0.0 , sum_conj = 0.0 ;
    
  // perform an SD start
  steep_Landau_FA( CG -> red , lat , FFTW ) ;

  for( k = 0 ; k < Latt.Nthreads ; k++ ) {
    trAA += CG -> red[ LINE_NSTEPS + 2 + k*CLINE ] ;
  }
  *tr = trAA * GFNORM_LANDAU ;
  
  size_t i ;
#pragma omp for private(i)
  for( i = 0 ; i < TRUE_HERM*LVOLUME ; i++ ) {
    const size_t idx = i/LVOLUME , j = i%LVOLUME ;
    CG -> sn[idx][j] = CG -> in_old[idx][j] = FFTW -> in[idx][j] ;
  }

  // compute the quantity Tr( dA dA )
  sum_DER3( CG -> red , (const GLU_complex**)FFTW -> in ) ;
    
  for( k = 0 ; k < Latt.Nthreads ; k++ ) {
    inold += CG -> red[ LINE_NSTEPS + CLINE*k ] ;
  }
  
 top :
  // not totally sure why this needs to be here but result changes
  // if it isn't
  {
    #pragma omp barrier
  }

  // might not need a barrier here... not sure
  trAA = insum = sum_conj = 0.0 ;
  #pragma omp for private(k)
  for( k = 0 ; k < CLINE*Latt.Nthreads ; k++ ) {
    CG -> red[k] = 0.0 ;
  }
  
  // this ONLY works with out and in, make sure that is what we use
  FA_deriv( CG -> red , FFTW -> in , lat ) ;

  // normalise the measure
  for( k = 0 ; k < Latt.Nthreads ; k++ ) {
    trAA += CG -> red[ LINE_NSTEPS + 2 + k*CLINE ] ;
  }
  *tr = trAA * GFNORM_LANDAU ;
  
  // and FA
  FOURIER_ACCELERATE3( FFTW ) ;
  
  // is the sum of the trace of in * in_old
  sum_PR3( CG -> red , (const GLU_complex**)FFTW -> in , 
	   (const GLU_complex**)CG -> in_old ) ; 
  
  // reduction happens for all threads
  for( k = 0 ; k < Latt.Nthreads ; k++ ) {
    insum    += CG -> red[ LINE_NSTEPS + CLINE*k ] ;
    sum_conj += CG -> red[ LINE_NSTEPS + 1 + CLINE*k ] ;
  }
      
  // compute the beta value, who knows what value is best?
  double beta = PRfmax( 0.0 , ( sum_conj ) / inold ) ;
  
  // switch to the fletcher reeves
  if( *tr < CG_TOL ) {
    beta = insum / inold ;
  }
  inold = insum ;
  
#pragma omp for private(i)
  for( i = 0 ; i < TRUE_HERM*LVOLUME ; i++ ) {
    size_t idx = i/LVOLUME , j = i%LVOLUME ;
    GLU_complex *pin = FFTW -> in[idx] ;
    GLU_complex *pin_old = CG -> in_old[idx] ;
    GLU_complex *psn = CG -> sn[idx] ;
    psn[j] = pin[j] + beta * ( psn[j] ) ;
    pin_old[j] = pin[j] ;
  }

  if( *tr > CG_TOL ) {
    line_search_Landau( CG -> red , lat , (const GLU_complex**)CG -> sn ) ;
  } else {
    egauge_Landau( lat , (const GLU_complex**)CG -> sn , Latt.gf_alpha ) ;
  }

  loc_iters++ ;
  
  // increment
  if( ( *tr > acc ) && ( loc_iters < max_iters ) ) goto top ;
  
  return loc_iters ;
}

// overwrites the lattice links in lat to Landau gauge fixed links using CG
size_t
FACG( struct site *lat , 
      struct fftw_stuff *FFTW ,
      double *th ,
      const double acc ,
      const size_t max_iters )
{
  // set up the maximum and what have you
  GLU_real max ; 
  *th = theta_test_lin( lat , &max , ND ) ; 
  size_t iters = 0 ;
  
  // have the option to leave early before allocations
  if( *th < acc ) return iters ;
  
  struct CGtemps CG ;
  if( allocate_temp_lg( &CG , GLU_TRUE ) == GLU_FAILURE ) {
    fprintf( stderr , "[FACG] allocation of temporaries failed\n" ) ;
    goto end ;
  }

#pragma omp parallel
  {
    size_t loc_iters = 1 ;
  top :
    loc_iters = steep_Landau_FACG( lat , FFTW , &CG , th ,
				   acc , max_iters ) ;

    #pragma omp master
    {
      iters = loc_iters ;
    }
    
    if( loc_iters >= max_iters ) {
      if( ( *th < 1E3*acc ) ) {
	#pragma omp master
	{
	  fprintf( stdout , "[GF] Continuation run \n" ) ;
	}
	goto top ;
      } else {
	iters = 123456789 ;
      }
    }
  }

 end :

  free_temp_lg( CG , GLU_TRUE ) ;

  return iters ; 
}

//returns the global gauge transform on lat
size_t
FASD( struct site *lat ,
      struct fftw_stuff *FFTW ,
      double *th ,
      const double acc ,
      const size_t max_iters )
{
  // set up the maximum and what have you
  GLU_real max ; 
  *th = theta_test_lin( lat , &max , ND ) ; 
  size_t iters = 0 ;
  
  // have the option to leave early before allocations
  if( *th < acc ) return iters ;
  
  struct CGtemps CG ;
  if( allocate_temp_lg( &CG , GLU_FALSE ) == GLU_FAILURE ) {
    fprintf( stderr , "[FACG] allocation of temporaries failed\n" ) ;
    goto end ;
  }

#pragma omp parallel
  {
    size_t loc_iters = 1 ;    
  top :
    loc_iters = steep_Landau_FASD( lat , FFTW , &CG , th ,
				   acc , max_iters ) ;
    
    if( loc_iters >= max_iters ) {
      if( ( *th < 1E3*acc ) ) {
	loc_iters = 0 ;
	#pragma omp master
	{
	  iters += loc_iters ;
	  fprintf( stdout , "[GF] Continuation run \n" ) ;
	}
	goto top ;
      } else {
	iters = 123456789 ;
      }
    }
#pragma omp master
    {
      iters += loc_iters ;
    }
  }

 end :

  free_temp_lg( CG , GLU_FALSE ) ;

  return iters ; 
}
