/**
    Copyright 2013-2018 Renwick James Hudspith

    This file (CFACG.c) is part of GLU.

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
   @file CFACG.c
   @brief Coulomb Gauge Fourier Accelerated Conjugate Gradient  
 */
#include "Mainfile.h"      // for all the definitions

#include "CG.h"            // for some of the shared routines
#include "CGalloc.h"       // allocations
#include "gftests.h"       // theta_test stopping condition
#include "gramschmidt.h"   // for reunit2, gram-schmidt reunitarisation
#include "gtrans.h"        // gauge transformations
#include "line_search.h"   // line searches for best alpha
#include "par_rng.h"       // Sunitary_gen()
#include "plaqs_links.h"   // plaquette routine might be called
#include "random_config.h" // for the random transformed gauge slices

#if (defined HAVE_IMMINTRIN_H) && !(defined SINGLE_PREC)
  #include <immintrin.h>
  #include "SSE2_OPS.h"
#endif

// prints out some relevant information ...
static void 
get_info( const size_t t ,
	  const double tr ,
	  const size_t iters , 
	  const int control ,
	  const double accuracy )
{
#pragma omp master
  {
    if( t == 0 ) {
      fprintf( stdout , "\n" ) ;
    } if( tr < accuracy ) {
      fprintf( stdout , "[GF] Slice :: %zu {Stopped by convergence} "
	       "\n[GF] Accuracy :: %1.5e"
	       " || Iterations :: %zu\n[GF] Failures :: %d\n" , 
	       t , tr , iters , control ) ; 
    } else {
      fprintf( stdout , "[GF] Slice :: %zu {Stopped by iterations too high} \n"
	       "[GF] Accuracy :: %1.5e || Iterations :: %zu\n"
	       "[GF] Failures :: %d \n" , 
	       t , tr , iters , control ) ; 
    }
    fprintf( stdout , "\n" ) ;
  }
  return ;
}


// is the same for Landau and Coulomb just with different LENGTHS
void
FOURIER_ACCELERATE2( struct fftw_stuff *FFTW ) 
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
    for( i = 0 ; i < LCU ; i++ ) {
      FFTW -> out[ mu ][ i ] *= FFTW -> psq[i] ;
    }
    fftw_execute( back[mu] ) ; 
  }
#endif
  return ;
}

// bit that calculates the steepest descent with fourier acceleration
static void
steep_deriv_CG( GLU_complex **in ,
		const struct site *lat , 
		const GLU_complex **slice_gauge , 
		const size_t t )
{
  size_t i ;
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i ++ ) {
    // compute gauge rotated derivatives
    GLU_complex sum[ HERMSIZE ] GLUalign ;
    GLU_complex A[ HERMSIZE ] GLUalign ;
    GLU_complex B[ HERMSIZE ] GLUalign ;

    GLU_complex C[ NCNC ] GLUalign ;
    GLU_complex D[ NCNC ] GLUalign ;

    const size_t Uidx = i+LCU*t ;
    size_t mu ;
    for( mu = 0 ; mu < HERMSIZE ; mu++ ) {
      A[ mu ] = B[ mu ] = sum[ mu ] = 0.0 ;
    }
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      const size_t Ubck = lat[Uidx].back[mu] ;
      const size_t bck  = lat[i].back[mu] ;
      const size_t it   = lat[i].neighbor[mu] ;

      memcpy( C , lat[Uidx].O[mu] , NCNC*sizeof(GLU_complex) ) ;
      gtransform_local( slice_gauge[i] , C , slice_gauge[it] ) ;
      
      memcpy( D , lat[Ubck].O[mu] , NCNC*sizeof(GLU_complex) ) ;
      gtransform_local( slice_gauge[bck] , D , slice_gauge[i] ) ;
      
      Hermitian_proj_short( B , C ) ;
      Hermitian_proj_short( A , D ) ;
          
      #if NC == 3
      sum[0] += A[0] - B[0] ;
      sum[1] += A[1] - B[1] ;
      sum[2] += A[2] - B[2] ;
      sum[3] += A[3] - B[3] ;
      sum[4] += A[4] - B[4] ;      
      #elif NC == 2
      sum[0] += A[0] - B[0] ;
      sum[1] += A[1] - B[1] ;
      #else 
      size_t nu ;
      for( nu = 0 ; nu < HERMSIZE ; nu++ ) {
	sum[ nu ] += A[nu] - B[nu] ;
      }
      #endif
    }

    // make it antihermitian
    #if NC == 3
    // for SU(3) I pack the first FFT with the two explicitly
    // real diagonal elements
    // to save on a Fourier transform ..
    in[0][i] = I * creal( sum[0] ) - creal( sum[3] ) ;
    in[1][i] = I * creal( sum[1] ) - cimag( sum[1] ) ;
    in[2][i] = I * creal( sum[2] ) - cimag( sum[2] ) ;
    in[3][i] = I * creal( sum[4] ) - cimag( sum[4] ) ;
    #elif NC == 2
    in[0][i] = I * creal( sum[0] ) ;
    in[1][i] = I * creal( sum[1] ) - cimag( sum[1] ) ;
    #else
    for( mu = 0 ; mu < HERMSIZE ; mu++ ) {
      in[mu][i] = I * creal( sum[mu] ) - cimag( sum[mu] ) ;
    }
    #endif
  }
  return ;
}

static void
sum_DER2( double *red ,
	  const GLU_complex **in )
{
  double c[1] = { 0. } ;
  size_t i ;
  #pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
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
    GLU_complex temp[ HERMSIZE ] ;
    size_t mu ;
    for( mu = 0 ; mu < HERMSIZE ; mu++ ) {
      temp[mu] = in[mu][i] ;
    }
    GLU_real tr ;
    trace_ab_herm_short( &tr , temp , temp ) ;
    loc_sum = 0.5 * (double)tr ;
#endif
    const size_t th = get_GLU_thread() ;

    // kahan summations
    const double y = 2*loc_sum - c[0] ;
    const double t = red[ LINE_NSTEPS + th * CLINE ] + y ;
    c[0] = ( t - red[ LINE_NSTEPS + th * CLINE ] ) - y ;
    red[ LINE_NSTEPS + th * CLINE ] = t ;
  } 
  return ;
}

static void
sum_PR( double *red ,
	const GLU_complex **in , 
	const GLU_complex **in_old )
{
  double c[2] = { 0. , 0. } ;
  size_t i ;
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    double loc_sum1 = 0.0 , loc_sum2 = 0.0 ;
#if NC == 3
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
#else
    GLU_complex temp1[ HERMSIZE ] , temp2[ HERMSIZE ] ;
    size_t mu ;
    for( mu = 0 ; mu < HERMSIZE ; mu++ ) {
      temp1[mu] = in[mu][i] ;
      temp2[mu] = in[mu][i] - in_old[mu][i] ;
    }
    GLU_real tr1 , tr2 ;
    trace_ab_herm_short( &tr1 , temp1 , temp1 ) ;
    trace_ab_herm_short( &tr2 , temp1 , temp2 ) ;
    loc_sum1 = 0.5 * (double)tr1 ;
    loc_sum2 = 0.5 * (double)tr2 ;
#endif
    const size_t th = get_GLU_thread() ;

    // kahan summations
    const double y = 2*loc_sum1 - c[0] ;
    const double t = red[ LINE_NSTEPS + th * CLINE ] + y ;
    c[0] = ( t - red[ LINE_NSTEPS + th * CLINE ] ) - y ;
    red[ LINE_NSTEPS + th * CLINE ] = t ;

    const double y2 = loc_sum2 - c[1] ;
    const double t2 = red[ LINE_NSTEPS + 1 + th * CLINE ] + y2 ;
    c[1] = ( t2 - red[ LINE_NSTEPS + 1 + th * CLINE ] ) - y2 ;
    red[ LINE_NSTEPS + 1 + th * CLINE ] = t2 ;
  }

  return ;
}

// perform a Fourier-accelerated SD step
static void
steep_step_SD( GLU_complex **slice_gauge , 
	       struct CGtemps CG ,
	       struct fftw_stuff *FFTW ,
	       const struct site *lat , 
	       const size_t t )
{
  steep_deriv_CG( FFTW -> in , lat , (const GLU_complex**)slice_gauge , t ) ;

  // if we want to Fourier accelerate, we call this otherwise it is the SD
  FOURIER_ACCELERATE2( FFTW ) ;
  
#if !(defined GLU_GFIX_SD) && !(defined GLU_FR_CG)
  line_search_Coulomb( CG.red , slice_gauge , CG.db , lat , 
		       (const GLU_complex**)FFTW -> in , t ) ;
#else
  exponentiate_gauge_CG( slice_gauge ,
			 (const GLU_complex**)FFTW -> in ,
			 Latt.gf_alpha ) ;
#endif
  return ;
}

// more up to date Polyak-Ribiere CG code
static size_t
steep_step_FACG( GLU_complex **gauge ,
		 struct CGtemps CG ,
		 struct fftw_stuff *FFTW ,
		 double *tr ,
		 const struct site *lat , 
		 const size_t t ,
		 const double accuracy , 
		 const size_t max_iters )
{
  // loop a set number of CG-iterations
  size_t k , i , loc_iters = 0 ;
  double trAA = 0.0 , inold = 0.0 , insum = 0.0 , sum_conj = 0.0 ;
  
  // set everything to zero
  #pragma omp for private(k)
  for( k = 0 ; k < LINE_NSTEPS+CLINE*Latt.Nthreads ; k++ ) {
    CG.red[ k ] = 0.0 ;
  }

  steep_deriv_CG( FFTW -> in , lat , (const GLU_complex**)gauge , t ) ;

  // if we want to Fourier accelerate, we call this otherwise it is the SD
  FOURIER_ACCELERATE2( FFTW ) ;

  sum_DER2( CG.red , (const GLU_complex**)FFTW -> in ) ;
  
  for( k = 0 ; k < Latt.Nthreads ; k++ ) {
    inold += CG.red[ LINE_NSTEPS + CLINE*k ] ;
  }
  *tr = inold * GFNORM_COULOMB ;
  
  line_search_Coulomb( CG.red , gauge , CG.db , lat , 
		       (const GLU_complex**)FFTW -> in , t ) ;
    
#pragma omp for private(i)
  for( i = 0 ; i < TRUE_HERM*LCU ; i++ ) {
    const size_t idx = i/LCU , j = i%LCU ;
    CG.sn[idx][j] = CG.in_old[idx][j] = FFTW -> in[idx][j] ;
  }
    
  // have a goto behaving like a while here
 top :

  // set everything to zero
  insum = sum_conj =  trAA = 0.0 ;
  for( k = 0 ; k < LINE_NSTEPS + CLINE*Latt.Nthreads ; k++ ) {
    CG.red[ k ] = 0.0 ;
  }
  
  steep_deriv_CG( FFTW -> in , lat , (const GLU_complex**)gauge ,  t ) ;
  
  FOURIER_ACCELERATE2( FFTW ) ;
  
  sum_PR( CG.red , (const GLU_complex**)FFTW -> in ,
	  (const GLU_complex**)CG.in_old ) ;

  // reduction happens for all threads
  for( k = 0 ; k < Latt.Nthreads ; k++ ) {
    insum    += CG.red[ LINE_NSTEPS + CLINE*k ] ;
    sum_conj += CG.red[ LINE_NSTEPS + 1 + CLINE*k ] ;
  }
  *tr = insum * GFNORM_COULOMB ;
  
  // compute the beta value, who knows what value is best?
  double beta = PRfmax( 0.0 , sum_conj / inold ) ;
  // switch to the fletcher reeves
  if( *tr < CG_TOL ) {
    beta = insum / inold ;
  }
  inold = insum ;
  
  #pragma omp for private(i)
  for( i = 0 ; i < TRUE_HERM*LCU ; i++ ) {
    const size_t idx = i/LCU , j = i%LCU ;
    const GLU_complex *pin = FFTW -> in[idx] ;
    GLU_complex *pin_old = CG.in_old[idx] ;
    GLU_complex *psn = CG.sn[idx] ;
    psn[j] = pin[j] + beta * ( psn[j] ) ;
    pin_old[j] = pin[j] ;
  }

  if( *tr > CG_TOL) {
    line_search_Coulomb( CG.red , gauge , CG.db , lat ,
			 (const GLU_complex**)CG.sn , t ) ;
  } else {
    exponentiate_gauge_CG( gauge ,
			   (const GLU_complex**)CG.sn , Latt.gf_alpha ) ;
  }

  loc_iters++ ;

  // increment
  if( ( *tr > accuracy ) && ( loc_iters < max_iters ) ) goto top ;
  
  return loc_iters ;
}

// steepest descents step
static size_t
steep_step_FASD( GLU_complex **gauge ,
		 struct CGtemps CG ,
		 struct fftw_stuff *FFTW ,
		 double *tr ,
		 const struct site *lat , 
		 const size_t t ,
		 const double accuracy , 
		 const size_t max_iters )
{
  double trAA ;
  size_t loc_iters = 0 ;
  size_t k ;

 top :

  // need this barrier to sync reduction array
  {
    #pragma omp barrier
  }
  
  *tr = 0.0 ;
  
  // perform a Fourier accelerated step
#pragma omp for private(k)
  for( k = 0 ; k < LINE_NSTEPS + CLINE*Latt.Nthreads ; k++ ) {
    CG.red[ k ] = 0.0 ;
  }
  
  steep_step_SD( gauge , CG , FFTW , lat , t ) ;

  sum_DER2( CG.red , (const GLU_complex**)FFTW -> in ) ;

  // all threads see this
  trAA = 0.0 ;
  for( k = 0 ; k < (size_t)Latt.Nthreads ; k++ ) {
    trAA += CG.red[ LINE_NSTEPS + CLINE*k ] ;
  }  
  *tr = trAA * GFNORM_COULOMB ;
  
  loc_iters++ ;
  
  if( ( *tr > accuracy ) && ( loc_iters < max_iters) ) goto top ;
  
  return loc_iters ;
}

// test the gauge transformation solution
static GLU_bool
adequate_solution( GLU_complex **slice_gauge , 
		   const size_t iters ) 
{
  GLU_bool failure = GLU_TRUE ;
  size_t i ;
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {

    // reunitarise the solution
    gram_reunit( slice_gauge[i] ) ;
    
    const double *p = (const double*)slice_gauge[i] ;
    register double sum = 0.0 ;
    size_t j ;
    for( j = 0 ; j < NCNC ; j++ ) {
      sum += *p * ( *p ) ; p++ ; // re^2
      sum += *p * ( *p ) ; p++ ; // im^2
    }
    if( fabs( sum - NC ) > iters*PREC_TOL ) {
      failure = GLU_FALSE ;
    }
  }
  return failure ;
}

// coulomb gauge fix on slice t
static size_t
steep_fix( GLU_complex **gauge ,
	   struct CGtemps CG ,
	   struct fftw_stuff *FFTW ,
	   const struct site *lat , 
	   const size_t t ,
	   const double accuracy ,
	   const size_t max_iters ,
	   size_t (*f)( GLU_complex **gauge ,
			struct CGtemps CG ,
			struct fftw_stuff *FFTW ,
			double *tr ,
			const struct site *lat , 
			const size_t t ,
			const double accuracy , 
			const size_t max_iters ) )
{
  // perform a Fourier accelerated step, this is where the CG can go ...
    double tr = 1.0 ;
    size_t loc_iters = 0 , tot_iters = 0 , control = 0 ;
    GLU_bool trans_flag = GLU_FALSE ;

 top :
    
    loc_iters = f( gauge , CG , FFTW , &tr ,
		   lat , t , accuracy , max_iters ) ;
    
    // quick test for unitarity of our gauge transformation matrices
    if( tr < accuracy ) {
      if( adequate_solution( gauge , loc_iters > 0 ? loc_iters : 1 ) 
	  == GLU_FALSE ) {
	trans_flag = GLU_TRUE ;
      }
    }

    // attempt a continuation run if we are close to the desired accuracy
    if( loc_iters >= (max_iters) ) {
      if( tr < 1E4*accuracy && !isnan( tr ) && !isinf( tr ) ) {
	tot_iters += loc_iters ;
	goto top ;
      }
      if( tr > accuracy || isnan( tr ) || isinf( tr ) ) {
	trans_flag = GLU_TRUE ;
      }
    }

    // random gauge transformation
    if( ( control < 8 ) && ( trans_flag == GLU_TRUE ) ) {
      #pragma omp master
      {
	fprintf( stdout , "[CG] random gauge transform %zu :: %e (%zu)\n" ,
		 control , tr , loc_iters ) ;
      }
      size_t i ;
      #pragma omp for private(i)
      for( i = 0 ; i < LCU ; i++ ) {
	Sunitary_gen( gauge[i] , get_GLU_thread( ) ) ;
      }
      tot_iters += loc_iters ;
      control++ ;
      trans_flag = GLU_FALSE ;
      goto top ;
    }
    tot_iters += loc_iters ;

    // get the general gauge fixing information
    get_info( t , tr , tot_iters , control , accuracy ) ;

    return tot_iters ;
}

size_t
Coulomb_FA( struct site  *__restrict lat , 
	    struct fftw_stuff *FFTW ,
	    const double accuracy ,
	    const size_t max_iter ,
	    const GLU_bool FACG )
{
  // init the rng
  initialise_par_rng( NULL ) ;
  
  // set the iterations
  static size_t (*f)( GLU_complex **gauge ,
		      struct CGtemps CG ,
		      struct fftw_stuff *FFTW ,
		      double *tr ,
		      const struct site *lat , 
		      const size_t t ,
		      const double accuracy , 
		      const size_t max_iters ) ;
  if( FACG == GLU_TRUE ) {
    f = steep_step_FACG ;
  } else {
    f = steep_step_FASD ;
  }
  
  // CG temporaries
  struct CGtemps CG ;
  CG.sn = NULL ; CG.in_old = NULL ; CG.red = NULL ;

  // gauge allocations 
  struct gauges G ;
  
  // flag for if something goes wrong
  size_t tot_its = 0 ;

  // allocations
  if( allocate_temp_cg( &G , &CG , FACG ) == GLU_FAILURE ) {
    fprintf( stderr , "[CG] problem allocating temporary Coulomb arrays\n" ) ;
    goto memfree ;
  }

  #pragma omp parallel
  {
    size_t t = 0 , i , thread_its = 0 ;
    // OK so we have set up the gauge transformation matrices
    thread_its = steep_fix( G.g_end , CG , FFTW ,
			    lat , t , accuracy , max_iter , f ) ;
    
    // and t+1
    thread_its += steep_fix( G.g , CG , FFTW ,
			     lat , t+1 , accuracy , max_iter , f ) ;
       
    //gauge transform the links at x and set slice_gauge_up to be slice_gauge
    gtransform_slice_th( (const GLU_complex **)G.g_end , lat , 
			 (const GLU_complex **)G.g , 0 ) ;
    
    //now we do the same for all time slices
    for( t = 2 ; t < Latt.dims[ ND - 1 ] ; t++ ) {

      #pragma omp for private(i) 
      for( i = 0 ; i < LCU ; i++ ) { identity( G.g_up[i] ) ; }

      // gauge fix on this slice
      thread_its += steep_fix( G.g_up , CG , FFTW , lat ,
			       t , accuracy , max_iter , f ) ;
            
      //gauge transform the links for this slice "g"
      gtransform_slice_th( (const GLU_complex **)G.g , lat , 
			   (const GLU_complex **)G.g_up , t - 1 ) ;
    
      // and copy "g_up" (the working transformation matrices) into "g"
      #pragma omp single
      {
	GLU_complex **ptr = G.g ;
	G.g = G.g_up ;
	G.g_up = ptr ;
      }
    }

    // no need for a reunitarisation step as it has already been done
    // gauge transform the very final slice
    gtransform_slice_th( (const GLU_complex **)G.g , lat , 
			 (const GLU_complex **)G.g_end , t - 1 ) ;

#pragma omp master
    {
      tot_its = thread_its ;
    }
  }

  // and free all of that memory, especially rotato
 memfree :

  free_par_rng() ;

  free_temp_cg( G , CG , FACG ) ;
 
  // and return the total iterations
  return tot_its ;
}
