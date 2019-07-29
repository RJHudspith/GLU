/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (line_search.c) is part of GLU.

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
   @file line_search.c
   @brief Line search routines
 */
#include "Mainfile.h"

#include "CG.h"          // set gauge matrix
#include "GLU_splines.h" // GLUbic spline interpolation code
#include "gtrans.h"      // gauge transformations

// line search probes for the Coulomb CG-gauge fixing
#define PC1 (0.17)
#define PC2 (0.34)

// line search probes for the Landau gauge fixing
#define PL1 (0.15)
#define PL2 (0.30)

// could have several different searches here
double
approx_minimum( const size_t nmeas , 
		const double alphas[ nmeas ] ,
		const double functional[ nmeas ] )
{
  // find the minimum from the lagrange interpolating polynomial
  const double A = functional[0]/((alphas[0]-alphas[1])*(alphas[0]-alphas[2])) ;
  const double B = functional[1]/((alphas[1]-alphas[0])*(alphas[1]-alphas[2])) ;
  const double C = functional[2]/((alphas[2]-alphas[0])*(alphas[2]-alphas[1])) ;
  
  const double Sol = (alphas[0]*(B+C)+alphas[1]*(A+C)+alphas[2]*(A+B))/(2*(A+B+C)) ;

  if( isnan( Sol ) || isinf( Sol ) || Sol > alphas[2] || Sol < alphas[0] ) {
    return 0.08 ;
  }
  
  return Sol ;
}

void
egauge_Landau( struct site *lat ,
	       const GLU_complex **in ,
	       const double alpha )
{
  size_t i ;
  // the derivative in this form is antihermitian i.e -> i.dA
#pragma omp for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    GLU_complex A[ NCNC ] GLUalign ;
    GLU_complex B[ NCNC ] GLUalign ;
    set_gauge_matrix( A , in , alpha , i ) ;
    #if ND==4
    size_t it = lat[i].neighbor[0] ;
    set_gauge_matrix( B , in , alpha , it ) ;
    gtransform_local( A , lat[i].O[0] , B ) ;

    it = lat[i].neighbor[1] ;
    set_gauge_matrix( B , in , alpha , it ) ;
    gtransform_local( A , lat[i].O[1] , B ) ;
    
    it = lat[i].neighbor[2] ;
    set_gauge_matrix( B , in , alpha , it ) ;
    gtransform_local( A , lat[i].O[2] , B ) ;

    it = lat[i].neighbor[3] ;
    set_gauge_matrix( B , in , alpha , it ) ;
    gtransform_local( A , lat[i].O[3] , B ) ;
    #else
    size_t mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      size_t it = lat[i].neighbor[mu] ;
      set_gauge_matrix( B , in , alpha , it ) ;
      gtransform_local( A , lat[i].O[mu] , B ) ;
    }
    #endif
  }
  return ;
}

void
exponentiate_gauge_CG( GLU_complex **gauge , 
		       const GLU_complex **in ,
		       const double alpha )
{
  size_t i ;
  // the derivative in this form is antihermitian i.e -> i.dA
#pragma omp for private(i)
  for( i = 0 ; i < LCU ; i++ ) {
    GLU_complex temp[ NCNC ] GLUalign , temp2[ NCNC ] GLUalign ;
    set_gauge_matrix( temp , in , alpha , i ) ;
    memcpy( temp2 , gauge[i] , NCNC*sizeof( GLU_complex ) ) ;
    multab_suNC( gauge[i] , temp , temp2 ) ;
  }
  return ;
}

// perform a line search using GLUbic splines for approximately the best alpha
void
line_search_Coulomb( double *red ,
		     GLU_complex **gauge ,
		     const struct draughtboard db ,
		     const struct site *lat ,
		     const GLU_complex **in ,
		     const size_t t )
{
  double c[ LINE_NSTEPS ] = { 0. , 0. , 0. } ;
  size_t i ;
#pragma omp for private(i)
  for( i = 0 ; i < db.Nsquare[0] ; i++ ) {
    GLU_complex A1[ NCNC ] GLUalign ;
    GLU_complex A2[ NCNC ] GLUalign ;
    GLU_complex B[ NCNC ] GLUalign ;
    GLU_complex C[ NCNC ] GLUalign ;
    GLU_complex D[ NCNC ] GLUalign ;
    
    double loc_v[ LINE_NSTEPS ] ;
    size_t mu , n ;
    for( mu = 0 ; mu < LINE_NSTEPS ; mu++ ) {
      loc_v[ mu ] = 0.0 ;
    }
    const size_t idx = db.square[0][i] ;
    size_t bck , fwd ;
    // evaluate for probes -> We are still calling set_gauge_matrix
    // wayyy too often!
    set_gauge_matrix( A1 , in , PC1 , idx ) ;
    set_gauge_matrix( A2 , in , PC2 , idx ) ;
    for( mu = 0 ; mu < ND-1 ; mu++ ) {
      fwd = lat[idx].neighbor[mu] ;
      bck = lat[idx].back[mu] ;

      memcpy( C , lat[idx+LCU*t].O[mu] , NCNC*sizeof( GLU_complex ) ) ;
      gtransform_local( gauge[idx] , C , gauge[fwd] ) ;

      memcpy( D , lat[bck+LCU*t].O[mu] , NCNC*sizeof( GLU_complex ) ) ;
      gtransform_local( gauge[bck] , D , gauge[idx] ) ;

      loc_v[0] += creal( trace( C ) ) ;
      loc_v[0] += creal( trace( D ) ) ;

      // positive ones
      set_gauge_matrix( B , in , PC1 , fwd ) ;
      loc_v[1] += Re_trace_abc_dag_suNC( A1 , C , B ) ;
      // negative ones
      set_gauge_matrix( B , in , PC1 , bck ) ;
      loc_v[1] += Re_trace_abc_dag_suNC( B , D , A1 ) ;
      // positive ones
      set_gauge_matrix( B , in , PC2 , fwd ) ;
      loc_v[2] += Re_trace_abc_dag_suNC( A2 , C , B ) ;
      // negative ones
      set_gauge_matrix( B , in , PC2 , bck ) ;
      loc_v[2] += Re_trace_abc_dag_suNC( B , D , A2 ) ;
    }
    const size_t th = get_GLU_thread() ;
    // reductions
    for( n = 0 ; n < LINE_NSTEPS ; n++ ) {
      const double y = loc_v[n] - c[n] ;
      const double t = red[ n + th * CLINE ] + y ;
      c[n] = ( t - red[ n + th * CLINE ] ) - y ;
      red[ n + th * CLINE ] = t ;
    }
  }
  // compute the vals
  double val[ LINE_NSTEPS ] ;
  size_t j ;
  for( j = 0 ; j < LINE_NSTEPS ; j++ ) {
    val[j] = 0.0 ;
    for( i = 0 ; i < Latt.Nthreads ; i++ ) {
      val[ j ] -= red[ j + CLINE*i ] ;
    }
  }

  const double Calcg[ LINE_NSTEPS ] = { 0.0 , PC1 , PC2 } ;
  const double min = approx_minimum( LINE_NSTEPS , Calcg , val ) ;
 
  exponentiate_gauge_CG( gauge , in , min ) ;

  return ;
}

void
line_search_Landau( double *red ,
		    struct site *lat ,
		    const GLU_complex **in )
{
  size_t i ;
  double c[ LINE_NSTEPS ] = { 0. , 0. , 0. } ;
  
  #pragma omp for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {

    GLU_complex A1[ NCNC ] GLUalign ;
    GLU_complex A2[ NCNC ] GLUalign ;
    GLU_complex B[ NCNC ] GLUalign ;
    
    double loc_v[ LINE_NSTEPS ] ;
    size_t mu , n ;
    for( mu = 0 ; mu < LINE_NSTEPS ; mu++ ) {
      loc_v[ mu ] = 0.0 ;
    }

    // these two don't care about mu so can be set out here
    set_gauge_matrix( A1 , in , PL1 , i ) ;
    set_gauge_matrix( A2 , in , PL2 , i ) ;
    
    for( mu = 0 ; mu < ND ; mu++ ) {

      // first one is just a trace
      #if NC == 3
      loc_v[0] += creal( lat[i].O[mu][0] )
	+ creal( lat[i].O[mu][4] )
	+ creal( lat[i].O[mu][8] ) ;
      #else
      loc_v[0] += creal( trace( lat[i].O[mu] ) ) ;
      #endif

      set_gauge_matrix( B , in , PL1 , lat[i].neighbor[mu] ) ;
      loc_v[1] += Re_trace_abc_dag_suNC( A1 , lat[i].O[mu] , B ) ;

      set_gauge_matrix( B , in , PL2 , lat[i].neighbor[mu] ) ;
      loc_v[2] += Re_trace_abc_dag_suNC( A2 , lat[i].O[mu] , B ) ;
    }
    
    const size_t th = get_GLU_thread() ;
    // reductions
    for( n = 0 ; n < LINE_NSTEPS ; n++ ) {
      const double y = loc_v[n] - c[n] ;
      const double t = red[ n + th * CLINE ] + y ;
      c[n] = ( t - red[ n + th * CLINE ] ) - y ;
      red[ n + th * CLINE ] = t ;
    }
  }
  // compute the vals
  double val[ LINE_NSTEPS ] ;
  size_t j ;
  for( j = 0 ; j < LINE_NSTEPS ; j++ ) {
    val[j] = 0.0 ;
    for( i = 0 ; i < Latt.Nthreads ; i++ ) {
      val[ j ] -= red[ j + CLINE*i ] ;
    }
  }

  const double Lalcg[ LINE_NSTEPS ] = { 0.0 , PL1 , PL2 } ;
  const double min = approx_minimum( LINE_NSTEPS , Lalcg , val ) ;

  egauge_Landau( lat , in , min ) ;
  
  return ; 
}

#undef PC1
#undef PC2
#undef PL1
#undef PL2
