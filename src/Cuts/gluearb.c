/**
   @file gluearb.c
   @brief arbitrary gluon prop
 */
#include "Mainfile.h"

#include "geometry.h"

static int
arb_dft( GLU_complex Ap[ND][NCNC] ,
	 const struct site *A ,
	 const double theta_mu[ND] )
{
  double complex fac[ ND ] ;
  double twiddle[ ND ] ;
  size_t i , mu ;
  
  for( mu = 0 ; mu < ND ; mu++ ) {
    zero_mat( Ap[mu] ) ;
    fac[mu] = cos( theta_mu[mu]/2. ) + I * sin( theta_mu[mu]/2. ) ;
  }

  for( i = 0 ; i < LVOLUME ; i++ ) {

    int x[ ND ] ;
    get_mom_2piBZ( x , i , ND ) ;
      
    double q_dot_x = 0.0 ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      if( x[mu] < Latt.dims[mu]/2 ) {
	x[mu] = x[mu] ;
      } else {
	x[mu] = x[mu] - Latt.dims[mu] ;
      }
      q_dot_x += x[mu] * theta_mu[mu] ;
    }
      
    const double complex eiqx = cos( q_dot_x ) + I * sin( q_dot_x ) ;

    // do the fft
    size_t j ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      for( j = 0 ; j < NCNC ; j++ ) {
	Ap[mu][j] += eiqx * A[i].O[mu][j] * fac[mu] ;
      }
    }
    //
  }
  return GLU_SUCCESS ;
}

// do a DFT in the time direction
static int
dft( const struct site *A )
{
  GLU_complex At[Latt.dims[ND-1]][ND][NCNC] ;
  size_t i , j , t , mu , k ;

  // compute t-sum
  for( t = 0 ; t < Latt.dims[ND-1] ; t++ ) {

    for( mu = 0 ; mu < ND ; mu++ ) {
      zero_mat( At[t][mu] ) ;
    }
    
    for( j = LCU*t ; j < LCU*(t+1) ; j++ ) {
      for( mu = 0 ; mu < ND ; mu++ ) {
	for( k = 0 ; k < NCNC ; k++ ) {
	  At[t][mu][k] += A[j].O[mu][k] ;
	}
      }
    }
  }

  // do the DFT
  int p ;
  for( p = -(int)Latt.dims[ND-1]/2 ; p <  (int)Latt.dims[ND-1]/2 ; p++ ) {
    GLU_complex Ap[ND][NCNC] ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      zero_mat( Ap[mu] ) ;
    }
    const double pmu = 2 * sin( p * Latt.twiddles[ND-1] / 2. ) ;
    
    for( t = 0 ; t < Latt.dims[ND-1] ; t++ ) {
      const GLU_complex exp =
	cos( p * Latt.twiddles[ND-1] * t ) +
	I * sin( p * Latt.twiddles[ND-1] * t ) ;
      for( mu = 0 ; mu < ND ; mu++ ) {
	for( k = 0 ; k < NCNC ; k++ ) {
	  Ap[mu][k] += exp * At[t][mu][k] ;
	}
      }
    }

    //
    double complex sum = 0.0 ;
    for( k = 0 ; k < NCNC ; k++ ) {
      sum += Ap[ND-1][k] * pmu ;
    }
    //printf( "sum :: %e \n" , creal( sum ) ) ;

    GLU_complex G , Gsum = 0.0 ;
    trace_ab_dag( &G , Ap[0] , Ap[0] ) ;
    Gsum += G ;
    trace_ab_dag( &G , Ap[1] , Ap[1] ) ;
    Gsum += G ;
    trace_ab_dag( &G , Ap[2] , Ap[2] ) ;
    Gsum += G ;
    trace_ab_dag( &G , Ap[3] , Ap[3] ) ;
    Gsum += G ;
    
    const double g2_norm   = 2.0 / ( ( NCNC - 1 ) * ( ND - 1 ) * LVOLUME ) ;

    printf( "%e %e \n" , pmu*pmu , creal( Gsum ) * g2_norm ) ;
    
  }
  return GLU_SUCCESS ;
}


// compute the log fields, overwrite the link matrices!
static int
Amu_fields( struct site *__restrict A ,
	    const lie_field_def def )
{
  size_t i ;

  // callback for the log definition
  void (*log)( GLU_complex Q[ NCNC ] ,
	       const GLU_complex U[ NCNC ] ) = Hermitian_proj ;
  switch( def ) {
  case LINEAR_DEF :
    break ;
  case LOG_DEF : 
    log = exact_log_slow ; 
    break ;
  }

#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < LVOLUME ; i++ ) { 
    GLU_complex temp[ NCNC ] GLUalign ;
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      log( temp , A[i].O[mu] ) ;
      equiv( A[i].O[mu] , temp ) ;
    }
  }
  return GLU_SUCCESS ;
}

// computes the arbitrary momentum space propagators
int 
cuts_gluearb( struct site *__restrict A ,
	      const struct cut_info CUTINFO ,
	      const struct sm_info SMINFO )
{
  // timelike DFT A_\mu(p) = \sum_x e^{ipx} A_\mu(x)
  Amu_fields( A , CUTINFO.dir ) ;

  const double g2_norm   = 2.0 / ( ( NCNC - 1 ) * ( ND - 1 ) * LVOLUME ) ;

  dft( A ) ;
  
  double theta = 0.0 ;
  for( theta = -6 ; theta < 6 ; theta += 0.1 ) {

    GLU_complex Ap[ND][NCNC] , Ap2[ND][NCNC] ;
    
    const double theta_mu[ND] = { 0,
				  theta*Latt.twiddles[1],
				  theta*Latt.twiddles[2] ,
				  theta*Latt.twiddles[3] } ;
    arb_dft( Ap , A , theta_mu ) ;

    const double theta_mu2[ND] = { 0,
				  -theta*Latt.twiddles[1],
				  -theta*Latt.twiddles[2] ,
				  -theta*Latt.twiddles[3] } ;
    arb_dft( Ap2 , A , theta_mu2 ) ;


    // check p_\mu A_\mu = 0
    double complex sum = 0.0 ;
    double qsq = 0.0 , mom[ND] ;
    size_t mu , nu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      mom[mu] = 2. * sin( theta_mu[mu] / 2. ) ;
      size_t j ;
      for( j = 0 ; j < NCNC ; j++ ) {
	sum += mom[mu] * Ap[mu][j] ;
      }
      qsq += mom[mu]*mom[mu] ;
    }
    printf( "WI %e %e %e\n" , qsq , creal( sum ) , cimag( sum ) ) ;

    GLU_complex Gp[ND][ND] ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      for( nu = 0 ; nu < ND ; nu++ ) {
	trace_ab( &Gp[mu][nu] , Ap[mu] , Ap2[nu] ) ;
      }
    }

    // project
    double trans = 0.0 , longt = 0.0 ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      for( nu = 0 ; nu < ND ; nu++ ) {
	const double common = mom[mu] * mom[nu] * (double)creal( Gp[mu][nu] ) / qsq ;
	trans += ( mu != nu ) ? -common : (double)creal( Gp[mu][nu] ) - common ;
	longt += -common ;
      }
    }

    printf( "LONG %e %e \n" , qsq , longt * g2_norm ) ;
    printf( "TRANS %e %e \n" , qsq , trans * g2_norm ) ;
    printf( "TPL %e %e \n" , qsq , ( trans + longt ) * g2_norm ) ;
  }

  return GLU_SUCCESS ;
}
