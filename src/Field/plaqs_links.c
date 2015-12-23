/*
    Copyright 2013 Renwick James Hudspith

    This file (plaqs_links.c) is part of GLU.

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
   @file plaqs_links.c
   @brief has all of the plaquette based measurements including gauge field strength tensor
 */

#include "Mainfile.h"
#include "clover.h"

// the one with the openmp locks, just an excuse to try this out to be honest
#if ( defined HAVE_OMP_H ) && ( defined _OPENMP )
  #define GLU_OMP_MEAS
  #include <omp.h>
#endif

//// My plaquette code, uses expressions for the trace of four matrices 
static double
complete_plaquette( const GLU_complex *__restrict a ,
		    const GLU_complex *__restrict b ,
		    const GLU_complex *__restrict c , 
		    const GLU_complex *__restrict d )
{
#if NC == 3
  register double complex tra1 = (a[0]*b[0]+a[1]*b[3]+a[2]*b[6]) ;
  register double complex tra2 = (a[0]*b[1]+a[1]*b[4]+a[2]*b[7]) ;
  register double complex tra3 = (a[0]*b[2]+a[1]*b[5]+a[2]*b[8]) ;

  double tra = creal( ( tra1 * conj(c[0]) + tra2 * conj( c[1] ) + tra3 * conj( c[2] ) ) * conj( d[0] ) ) ;
  tra += creal( ( tra1 * conj(c[3]) + tra2 * conj( c[4] ) + tra3 * conj( c[5] ) ) * conj( d[1] ) ) ;
  tra += creal( ( tra1 * conj(c[6]) + tra2 * conj( c[7] ) + tra3 * conj( c[8] ) ) * conj( d[2] ) ) ;

  tra1 = (a[3]*b[0]+a[4]*b[3]+a[5]*b[6]) ;
  tra2 = (a[3]*b[1]+a[4]*b[4]+a[5]*b[7]) ;
  tra3 = (a[3]*b[2]+a[4]*b[5]+a[5]*b[8]) ;

  double trb = creal( ( tra1 * conj(c[0]) + tra2 * conj( c[1] ) + tra3 * conj( c[2] ) ) * conj( d[3] ) ) ;
  trb += creal( ( tra1 * conj(c[3]) + tra2 * conj( c[4] ) + tra3 * conj( c[5] ) ) * conj( d[4] ) ) ;
  trb += creal( ( tra1 * conj(c[6]) + tra2 * conj( c[7] ) + tra3 * conj( c[8] ) ) * conj( d[5] ) ) ;

  tra1 = (a[6]*b[0]+a[7]*b[3]+a[8]*b[6]) ;
  tra2 = (a[6]*b[1]+a[7]*b[4]+a[8]*b[7]) ;
  tra3 = (a[6]*b[2]+a[7]*b[5]+a[8]*b[8]) ;

  double trc = creal( ( tra1 * conj(c[0]) + tra2 * conj( c[1] ) + tra3 * conj( c[2] ) ) * conj( d[6] ) ) ;
  trc += creal( ( tra1 * conj(c[3]) + tra2 * conj( c[4] ) + tra3 * conj( c[5] ) ) * conj( d[7] ) ) ;
  trc += creal( ( tra1 * conj(c[6]) + tra2 * conj( c[7] ) + tra3 * conj( c[8] ) ) * conj( d[8] ) ) ;

  return tra + trb + trc ;

#elif NC == 2
  
  register const double complex tra1 = (a[0]*b[0]+a[1]*b[2]) ;
  register const double complex tra2 = (a[0]*b[1]+a[1]*b[3]) ;
  double tra = (double)creal( ( tra1 * conj(c[0]) + tra2 * conj(c[1]) ) *conj(d[0]) ) ;
  tra += (double)creal( ( tra1 * conj(c[2]) + tra2 * conj(c[3]) )*conj(d[1]) ) ;
  register const double complex tra3 = (a[2]*b[0]+a[3]*b[2]) ;
  register const double complex tra4 = (a[2]*b[1]+a[3]*b[3]) ;
  double trb = (double)creal( ( tra3 * conj(c[0]) + tra4 * conj(c[1]))*conj(d[2]) ) ;
  trb += (double)creal( ( tra3 * conj(c[2]) + tra4 * conj(c[3]))*conj(d[3]) ) ;

  return tra + trb ;
#else
  GLU_complex temp1[NCNC] , temp2[NCNC] ; 
  multab( temp1 , a , b ) ; 
  multab_dag( temp2 , temp1 , c ) ; 
  multab_dag( temp1 , temp2 , d ) ; 
  double tr ;
  speed_trace_Re( &tr , temp1 ) ;
  return tr ;
#endif
}

// all of the plaquettes
double
all_plaquettes( const struct site *__restrict lat ,
		double *__restrict sp_plaq ,
		double *__restrict t_plaq )
{
  double spplaq = 0. , tplaq = 0.0 ;
  size_t i ; 
#pragma omp parallel for private(i) reduction(+:spplaq) reduction(+:tplaq) 
  for( i = 0 ; i < LVOLUME ; i++ ) {
    double p = 0. , face ;
    size_t mu , nu , s , t ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      t = lat[i].neighbor[mu] ; 
      for( nu = 0 ; nu < mu ; nu++ ) {
        s = lat[i].neighbor[nu] ;
	face = complete_plaquette( lat[ i ].O[mu] , lat[ t ].O[nu] , 
				   lat[ s ].O[mu] , lat[ i ].O[nu] ) ; 
	p = p + (double)face ;
      }
    }
    spplaq = spplaq + (double)p ;
    // reinitialise p for the temporal plaquette ...
    p = 0.0 ;
    t = lat[i].neighbor[mu] ; 
    for( nu = 0 ; nu < mu ; nu++ ) {
      s = lat[i].neighbor[nu] ;
      face = complete_plaquette( lat[ i ].O[mu] , lat[ t ].O[nu] , 
				 lat[ s ].O[mu] , lat[ i ].O[nu] ) ; 
      p = p + (double)face ;
    }
    tplaq = tplaq + (double)p ;
  }
  *sp_plaq = 2.0 * spplaq / (double)( ( ND - 1 ) * ( ND - 2 ) * NC * LVOLUME ) ;
  *t_plaq  = 2.0 * tplaq /  (double)( ( ND - 1 ) * ( ND - 2 ) * NC * LVOLUME ) ;
  return 0.5 * ( *sp_plaq + *t_plaq ) ; 
}

// average plaquette
double
av_plaquette( const struct site *__restrict lat )
{
  double plaq = 0. ;
  size_t i ; 
#pragma omp parallel for private(i) reduction(+:plaq) 
  for( i = 0 ; i < LVOLUME ; i++ ) {
    register double p = 0. , face ;
    size_t mu , nu , t , s ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      t = lat[i].neighbor[mu] ; 
      for( nu = 0 ; nu < mu ; nu++ ) {
	s = lat[i].neighbor[nu] ;
	face = complete_plaquette( lat[ i ].O[mu] , 
				   lat[ t ].O[nu] , 
				   lat[ s ].O[mu] , 
				   lat[ i ].O[nu] ) ; 
	p = p + (double)face ;
      }
    }
    plaq = plaq + (double)p ;
  }
  return 2.0 * plaq /(double)( NC * ND * ( ND - 1 ) * LVOLUME ) ; 
}

// spatial plaquettes ..
double 
s_plaq( const struct site *__restrict lat )
{
  double plaq = 0. ;
  size_t i ; 
#pragma omp parallel for private(i) reduction(+:plaq) 
  for( i = 0 ; i < LVOLUME ; i++ ) {
    register double p = 0. , face ;
    size_t mu , nu , t , s ;
    for( mu = 0 ; mu < ND - 1 ; mu++ ) {
      t = lat[i].neighbor[mu] ; 
      for( nu = 0 ; nu < mu ; nu++ ) {
	s = lat[i].neighbor[nu] ;
	face = complete_plaquette( lat[ i ].O[mu] , 
				   lat[ t ].O[nu] , 
				   lat[ s ].O[mu] , 
				   lat[ i ].O[nu] ) ; 
	p = p + (double)face ;
      }
    }
    plaq = plaq + (double)p ;
  }
  return 2.0 * plaq / (double)( ( ND - 1 ) * ( ND - 2 ) * NC * LVOLUME ) ; 
}

// just the temporal plaquettes ..
double 
t_plaq( const struct site *__restrict lat )
{
  double plaq = 0. ;
  size_t i ; 
  const size_t mu = ND - 1 ;
#pragma omp parallel for private(i) reduction(+:plaq) 
  for( i = 0 ; i < LVOLUME ; i++ ) {
    register double p = 0. , face ;
    size_t nu , t = lat[i].neighbor[mu] , s ; 
    for( nu = 0 ; nu < mu ; nu++ ) {
      s = lat[i].neighbor[nu] ;
      face = complete_plaquette( lat[ i ].O[mu] , 
				 lat[ t ].O[nu] , 
				 lat[ s ].O[mu] , 
				 lat[ i ].O[nu] ) ; 
      p = p + (double)face ;
    }
    plaq = plaq + (double)p ;
  }
  return 2.0 * ( plaq ) / (double)( ( ND - 1 ) * ( ND - 2 ) * NC * LVOLUME ) ; 
}

/// Compute the lattice F_\mu\nu term here.
double
lattice_gmunu( const struct site *__restrict lat ,
	       double *__restrict qtop ,
	       double *__restrict avplaq )
{
  *avplaq = av_plaquette( lat ) ;
#if ND == 4
  const double QTOP_DENOM = -0.001583143494411527678811 ; // 1.0/(64*Pi*Pi)
  double GG = 0. ;
  // compute G_munu in clover.c
  compute_Gmunu( &GG , qtop , lat ) ;
  // return the usual suspects, topological charge average plaq and F^2
  *qtop *= QTOP_DENOM ;
  return ( GG ) / ( 16. * LVOLUME ) ;
#else
  // don't have anything for this
  *qtop = -1.0 ;
  return -1.0 ;
#endif
}

// sexy wrapper for the topological charge measurements
int
gauge_topological_meas( const struct site *__restrict lat , 
			double *qtop_new , 
			double *qtop_old , 
			const int iter )
{
  // set up a tolerance for how close to an integer we wish to be
  const double QTOP_TOL = 5E-3 ;
  double avplaq ;
  const double plaq = lattice_gmunu( lat , qtop_new , &avplaq ) ;
  fprintf( stdout , "{iter} %d {w} %g {GG} %g {q} %g {diff} %g \n" ,
	   iter , avplaq , plaq , *qtop_new , fabs( *qtop_old - *qtop_new ) ) ;
  // rounded integer value for the topological charge
  const int qint = (int)( *qtop_new > 0. ? *qtop_new + 0.5 : *qtop_new - 0.5 ) ;
  // The two tests are the difference between the newest topological charge and
  // the previous, suggesting convergence.
  // And the more stringent test of "closeness to an integer" which is allowed to be
  // off by a little bit more. This is certainly a heuristic measure!
  if( fabs( *qtop_old - *qtop_new ) < QTOP_TOL &&
      fabs( qint - *qtop_new ) < ( 3 * QTOP_TOL ) ) {
    fprintf( stdout , "\n{CONFIG} %zu {QTOP} %d \n\n" , Latt.flow , qint ) ;
    return GLU_SUCCESS ;
  }
  *qtop_old = *qtop_new ;
  return GLU_FAILURE ;
}

//////////// LINK TRACE MEASUREMENTS ////////////

// all links
double
all_links( const struct site *__restrict lat ,
	   double *__restrict sp_link ,
	   double *__restrict t_link )
{
  double splink = 0.0 , tlink = 0.0 ; 
  size_t i ; 
#pragma omp parallel for private(i) reduction(+:splink) reduction(+:tlink)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    double p = 0. , res ;
    size_t mu ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      speed_trace_Re( &res , lat[i].O[mu] ) ; 
      p = p + (double)res ; 
    }
    splink = splink + (double)p ;
    // and the temporal
    speed_trace_Re( &res , lat[i].O[ND-1] ) ; 
    tlink = tlink + (double)res ; 
  }
  *sp_link = splink / ( ( ND - 1 ) * NC * LVOLUME ) ;
  *t_link = tlink / ( NC * LVOLUME ) ;
  return ( splink + tlink ) / ( ND * NC * LVOLUME ) ;
}

// General links with the maximum passed by reference, uses OMP locks
double 
indivlinks( const struct site *__restrict lat , GLU_real *max )
{
  double link = 0.0 ;
  *max = 0.0 ;
  size_t i ; 
  #ifdef GLU_OMP_MEAS
  omp_lock_t writelock ;
  omp_init_lock( &writelock ) ;
  #endif

#pragma omp parallel for private(i) reduction(+:link)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    double loc_link = 0. , res ;
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      speed_trace_Re( &res , lat[i].O[mu] ) ;
      loc_link += res ; 
      #ifdef GLU_OMP_MEAS
      if( unlikely( (GLU_real)res > *max ) ) { 
	omp_set_lock( &writelock ) ;
	*max = (GLU_real)res ;
	omp_unset_lock( &writelock ) ;
      } 
      #else
      if( unlikely( res > *max ) ) { *max = (GLU_real)res ; } 
      #endif
    } 
    // reduction 
    link = link + (double)loc_link ;
  }
  #ifdef GLU_OMP_MEAS
  omp_destroy_lock( &writelock ) ;
  #endif

  return ( link ) / (double)( ND * NC * LVOLUME ) ; 
}

// General links
double
links( const struct site *__restrict lat )
{
  double link = 0.0 ; 
  size_t i ; 
#pragma omp parallel for private(i) reduction(+:link)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    double p = 0. , res ;
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      speed_trace_Re( &res , lat[i].O[mu] ) ; 
      p = p + (double)res ; 
    }
    link = link + (double)p ;
  }
  return link / ( ND * NC * LVOLUME ) ; 
}

//spatial links
double 
s_links( const struct site *__restrict lat )
{
  double link = 0.0 ; 
  int i ; 
#pragma omp parallel for private(i) reduction(+:link)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    double p = 0. , res ;
    int mu ;
    for( mu = 0 ; mu < ND - 1 ; mu++ ) {
      speed_trace_Re( &res , lat[i].O[mu] ) ; 
      p += (double)res ; 
    }
    link = link + (double)p ;
  }
  return link/(double)( LVOLUME * ( ND - 1 ) * NC ) ; 
}

//temporal links
double 
t_links( const struct site *__restrict lat )
{    
  double link = 0.0 ; 
  size_t i ; 
#pragma omp parallel for private(i) reduction(+:link)
  for( i = 0 ; i < LVOLUME ; i++ ) { 
    double res = 0. ; 
    speed_trace_Re( &res , lat[i].O[ ND - 1 ] ) ; 
    link = link + (double)res ; 
  }
  return link /(double)( LVOLUME * NC ) ; 
}

// undefine this if it has been set
#ifdef GLU_OMP_MEAS
  #undef GLU_OMP_MEAS
#endif
