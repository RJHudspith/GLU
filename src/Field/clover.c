/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (clover.c) is part of GLU.

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
   @file clover.c
   @brief Calculation of the \f$ G_{\mu\nu} \f$ terms for the gauge force term from the symmetric lattice def

  Includes the naive gauge force term, using the clover definition.

  Includes the gauge field strength tensor using the simple, single, right-handed plaquette.

  I have also included the improved  gauge force term of Bilson-Thompson and Leinweber.

  @warning NOTE :: the improved clover term is turned on by the CLOVER_IMPROVE define
*/

/*                          
  Ascii picture :: G_{\mu\nu}  - - > ( mu ) ::  | ( -nu ) 
                                                v         
 · ----<----- · · -----<-----  
 |            | |            | 
 |            | |            | 
 v     4      ^ v     1      ^ 
 |            | |            | 
 |            | |            | 
 · ---->----- X X ----->---- · 
 · ----<----- X X -----<---- ·   + 0.5(1x2+2x1) + 0.5(1x3+3x1) + 
 |            | |            |                    2x2 + 3x3 contribs 
 |            | |            | 
 v    3       ^ v     2      ^ 
 |            | |            | 
 |            | |            | 
 · ---->----- · · ----->---- · 

 */

// can test the versions of the projection here ....
//#define TRF_ANTIHERMITIAN
//#define ANTIHERMITIAN
//#define CLOVER_LOG_DEF

#include "Mainfile.h"

// improvement factors from the bilson-thompson paper ...
#ifndef NK5
   GLU_real fClover_k1 ;
   GLU_real fClover_k2 ;
   GLU_real fClover_k3 ;
   GLU_real fClover_k4 ;
   GLU_real fClover_k5 ;
#else
   const GLU_real fClover_k1 = 19./9. ;
   const GLU_real fClover_k2 = 1./36. ;
   const GLU_real fClover_k3 = 32./90. ;
   const GLU_real fClover_k4 = 1./30. ;
   const GLU_real fClover_k5 = 0.0 ;
#endif

#if ND == 4 /// this definition is only correct for ND = 4
// top right
static void
compute_s1( GLU_complex *__restrict sum ,
	    GLU_complex u[ NCNC ] ,
	    GLU_complex v[ NCNC ] ,
	    const struct site *__restrict lat ,
	    const size_t i ,
	    const size_t mu ,
	    const size_t nu ,
	    const GLU_real mul )
{
  const size_t s1 = lat[i].neighbor[mu] ;
  const size_t s2 = lat[i].neighbor[nu] ;
  multab_suNC( u , lat[i].O[mu] , lat[s1].O[nu] ) ;
  multab_dag_suNC( v , u , lat[s2].O[mu] ) ;
  multab_dag_suNC( u , v , lat[i].O[nu] ) ;
#ifdef CLOVER_LOG_DEF
  exact_log_slow( v , u ) ;
  a_plus_Sxb( sum , v , mul ) ;
#else
  a_plus_Sxb( sum , u , mul ) ;
#endif
  return ;
}
#ifndef PLAQUETTE_FMUNU
// bottom right
static void
compute_s2( GLU_complex *__restrict sum ,
	    GLU_complex u[ NCNC ] ,
	    GLU_complex v[ NCNC ] ,
	    const struct site *__restrict lat ,
	    const size_t i ,
	    const size_t mu ,
	    const size_t nu ,
	    const GLU_real mul )
{
  const size_t s1 = lat[i].back[nu] ;
  const size_t s2 = lat[s1].neighbor[mu] ;
  multabdag_suNC( u , lat[s1].O[nu] , lat[s1].O[mu] ) ;
  multab_suNC( v , u , lat[s2].O[nu] ) ;
  multab_dag_suNC( u , v , lat[i].O[mu] ) ; 
#ifdef CLOVER_LOG_DEF
  exact_log_slow( v , u ) ;
  a_plus_Sxb( sum , v , mul ) ;
#else
  a_plus_Sxb( sum , u , mul ) ;
#endif
  return ;
}
// bottom left
static void
compute_s3( GLU_complex *__restrict sum ,
	    GLU_complex u[ NCNC ] ,
	    GLU_complex v[ NCNC ] ,
	    const struct site *__restrict lat ,
	    const size_t i ,
	    const size_t mu ,
	    const size_t nu ,
	    const GLU_real mul )
{
  const size_t s1 = lat[i].back[mu] ;
  const size_t s2 = lat[s1].back[nu] ;
  const size_t s3 = lat[i].back[nu] ;
  multab_dagdag_suNC( u , lat[s1].O[mu] , lat[s2].O[nu] ) ;
  multab_suNC( v , u , lat[s2].O[mu] ) ;
  multab_suNC( u , v , lat[s3].O[nu] ) ;
#ifdef CLOVER_LOG_DEF
  exact_log_slow( v , u ) ;
  a_plus_Sxb( sum , v , mul ) ;
#else
  a_plus_Sxb( sum , u , mul ) ;
#endif
  return ;
}
// top left
static void
compute_s4( GLU_complex *__restrict sum ,
	    GLU_complex u[ NCNC ] ,
	    GLU_complex v[ NCNC ] ,
	    const struct site *__restrict lat ,
	    const size_t i ,
	    const size_t mu ,
	    const size_t nu ,
	    const GLU_real mul )
{
  const size_t s2 = lat[i].back[mu] ;
  const size_t s1 = lat[s2].neighbor[nu] ;
  multab_dag_suNC( u , lat[i].O[nu] , lat[s1].O[mu] ) ;
  multab_dag_suNC( v , u , lat[s2].O[nu] ) ;
  multab_suNC( u , v , lat[s2].O[mu] ) ;
#ifdef CLOVER_LOG_DEF
  exact_log_slow( v , u ) ;
  a_plus_Sxb( sum , v , mul ) ;
#else
  a_plus_Sxb( sum , u , mul ) ;
#endif
  return ;
}
#endif

#ifdef CLOVER_IMPROVE 
// compute sector 1 with improvements 
// return the plaquette. Sector 1 is the
// top right of the clover term.
static void
compute_clover_s1( GLU_complex *__restrict sum ,
		   GLU_complex u[ NCNC ] ,
		   GLU_complex v[ NCNC ] ,
		   const struct site *__restrict lat ,
		   const size_t i ,
		   const size_t mu ,
		   const size_t nu )
{
  // 1x1 contribution
  if( likely( fClover_k1 > CL_TOL ) ) {
    compute_s1( sum , u , v , lat , i , mu , nu , Clover_k1 ) ;
  }
  // 2x2 contribution
  if( likely( fClover_k2 > CL_TOL ) ) {
    const size_t s1 = lat[i].neighbor[mu] ;
    const size_t s2 = lat[s1].neighbor[mu] ;
    const size_t s3 = lat[s2].neighbor[nu] ;
    multab_suNC( u , lat[i].O[mu] , lat[s1].O[mu] ) ;
    multab_suNC( v , u , lat[s2].O[nu] ) ;
    multab_suNC( u , v , lat[s3].O[nu] ) ;
    // the final link posize_ts to the top right corner
    // work backwards
    const size_t s6 = lat[i].neighbor[nu] ;
    const size_t s5 = lat[s6].neighbor[nu] ;    
    const size_t s4 = lat[s5].neighbor[mu] ; 
    multab_dag_suNC( v , u , lat[s4].O[mu] ) ;
    multab_dag_suNC( u , v , lat[s5].O[mu] ) ;
    multab_dag_suNC( v , u , lat[s6].O[nu] ) ;
    multab_dag_suNC( u , v , lat[i].O[nu] ) ;
    #ifdef CLOVER_LOG_DEF
    exact_log_slow( v , u ) ;
    a_plus_Sxb( sum , v , Clover_k2 ) ;
    #else
    a_plus_Sxb( sum , u , Clover_k2 ) ; 
    #endif 
  }
  // 1x2 contribution(s)
  if( likely( fClover_k3 > CL_TOL ) ) {
    // the 1x2 contrib
    const size_t s1 = lat[i].neighbor[mu] ;
    const size_t s2 = lat[s1].neighbor[nu] ;
    const size_t s3 = lat[i].neighbor[nu] ;
    const size_t s4 = lat[s3].neighbor[nu] ;
    multab_suNC( u , lat[i].O[mu] , lat[s1].O[nu] ) ;
    multab_suNC( v , u , lat[s2].O[nu] ) ;
    multab_dag_suNC( u , v , lat[s4].O[mu] ) ;
    multab_dag_suNC( v , u , lat[s3].O[nu] ) ;
    multab_dag_suNC( u , v , lat[i].O[nu] ) ;
    #ifdef CLOVER_LOG_DEF
    exact_log_slow( v , u ) ;
    a_plus_Sxb( sum , v , Clover_k3 ) ;
    #else
    a_plus_Sxb( sum , u , Clover_k3 ) ;
    #endif
    // and the 2x1
    const size_t t1 = lat[i].neighbor[mu] ;
    const size_t t2 = lat[t1].neighbor[mu] ;
    multab_suNC( u , lat[i].O[mu] , lat[t1].O[mu] ) ;
    multab_suNC( v , u , lat[t2].O[nu] ) ;
    const size_t t3 = lat[i].neighbor[nu] ;
    const size_t t4 = lat[t3].neighbor[mu] ;
    multab_dag_suNC( u , v , lat[t4].O[mu] ) ;
    multab_dag_suNC( v , u , lat[t3].O[mu] ) ;
    multab_dag_suNC( u , v , lat[i].O[nu] ) ;
    #ifdef CLOVER_LOG_DEF
    exact_log_slow( v , u ) ;
    a_plus_Sxb( sum , v , Clover_k3 ) ;
    #else
    a_plus_Sxb( sum , u , Clover_k3 ) ;
    #endif
  }
  // 1x3 contribution(s)
  if( likely( fClover_k4 > CL_TOL ) ) {
    // 1x3
    const size_t s1 = lat[i].neighbor[mu] ;
    const size_t s2 = lat[s1].neighbor[nu] ;
    const size_t s3 = lat[s2].neighbor[nu] ;
    multab_suNC( u , lat[i].O[mu] , lat[s1].O[nu] ) ;
    multab_suNC( v , u , lat[s2].O[nu] ) ;
    multab_suNC( u , v , lat[s3].O[nu] ) ;
    const size_t s4 = lat[s3].neighbor[nu] ;
    const size_t s5 = lat[s4].back[mu] ;
    const size_t s6 = lat[s5].back[nu] ;
    const size_t s7 = lat[s6].back[nu] ;
    multab_dag_suNC( v , u , lat[s5].O[mu] ) ;
    multab_dag_suNC( u , v , lat[s6].O[nu] ) ;
    multab_dag_suNC( v , u , lat[s7].O[nu] ) ;
    multab_dag_suNC( u , v , lat[i].O[nu] ) ;
    #ifdef CLOVER_LOG_DEF
    exact_log_slow( v , u ) ;
    a_plus_Sxb( sum , v , Clover_k4 ) ;
    #else
    a_plus_Sxb( sum , u , Clover_k4 ) ;
    #endif
    // 3x1
    const size_t t1 = lat[i].neighbor[mu] ;
    const size_t t2 = lat[t1].neighbor[mu] ;
    const size_t t3 = lat[t2].neighbor[mu] ;
    multab_suNC( u , lat[i].O[mu] , lat[t1].O[mu] ) ;
    multab_suNC( v , u , lat[t2].O[mu] ) ;
    multab_suNC( u , v , lat[t3].O[nu] ) ;
    const size_t t4 = lat[t3].neighbor[nu] ;
    const size_t t5 = lat[t4].back[mu] ;
    const size_t t6 = lat[t5].back[mu] ;
    const size_t t7 = lat[t6].back[mu] ;
    multab_dag_suNC( v , u , lat[t5].O[mu] ) ;
    multab_dag_suNC( u , v , lat[t6].O[mu] ) ;
    multab_dag_suNC( v , u , lat[t7].O[mu] ) ;
    multab_dag_suNC( u , v , lat[i].O[nu] ) ;
    #ifdef CLOVER_LOG_DEF
    exact_log_slow( v , u ) ;
    a_plus_Sxb( sum , v , Clover_k4 ) ;
    #else
    a_plus_Sxb( sum , u , Clover_k4 ) ;
    #endif
  }
#ifndef NK5
  // 3x3 contribution
  if( fClover_k5 > CL_TOL ) {
    GLU_complex x[ NCNC ] , y[ NCNC ] ;
    GLU_complex s[ NCNC ] , t[ NCNC ] ;
    const size_t s1 = lat[i].neighbor[mu] ;
    const size_t s2 = lat[s1].neighbor[mu] ;
    const size_t s3 = lat[s2].neighbor[mu] ;
    const size_t s4 = lat[s3].neighbor[nu] ;
    const size_t s5 = lat[s4].neighbor[nu] ;
    multab_suNC( s , lat[i].O[mu] , lat[s1].O[mu] ) ;
    multab_suNC( v , lat[s2].O[mu] , lat[s3].O[nu] ) ;
    multab_suNC( x , lat[s4].O[nu] , lat[s5].O[nu] ) ;
    multab_suNC( y , s , v ) ;
    multab_suNC( s , y , x ) ;
    // go to s6
    const size_t s7 = lat[i].neighbor[nu] ;
    const size_t s8 = lat[s7].neighbor[nu] ;
    const size_t s9 = lat[s8].neighbor[nu] ;
    const size_t s10 = lat[s9].neighbor[mu] ;
    const size_t s11 = lat[s10].neighbor[mu] ;
    multab_suNC( t , lat[i].O[nu] , lat[s7].O[nu] ) ;
    multab_suNC( v , lat[s8].O[nu] , lat[s9].O[mu] ) ;
    multab_suNC( x , lat[s10].O[mu] , lat[s11].O[mu] ) ;
    multab_suNC( y , t , v ) ;
    multab_suNC( t , y , x ) ;
    multab_dag_suNC( u , s , t ) ;
    #ifdef CLOVER_LOG_DEF
    exact_log_slow( v , u ) ;
    a_plus_Sxb( sum , v , k5 ) ;
    #else
    a_plus_Sxb( sum , u , k5 ) ;
    #endif
  }
#endif
  return ;
}

// Sector 2 is the
// bottom right of the clover term.
static void
compute_clover_s2( GLU_complex *__restrict sum ,
		   GLU_complex u[ NCNC ] ,
		   GLU_complex v[ NCNC ] ,
		   const struct site *__restrict lat ,
		   const size_t i ,
		   const size_t mu ,
		   const size_t nu )
{
  // 1x1 contribution
  if( likely( fClover_k1 > CL_TOL ) ) {
    compute_s2( sum , u , v , lat , i , mu , nu , Clover_k1 ) ;
  }
  // 2x2 contrib
  if( likely( fClover_k2 > CL_TOL ) ) {
    const size_t s1 = lat[i].back[nu] ;
    const size_t s2 = lat[s1].back[nu] ;
    const size_t s3 = lat[s2].neighbor[mu] ;
    multab_dagdag_suNC( u , lat[s1].O[nu] , lat[s2].O[nu] ) ;
    multab_suNC( v , u , lat[s2].O[mu] ) ;
    multab_suNC( u , v , lat[s3].O[mu] ) ;
    // posize_ted to halfway
    const size_t s4 = lat[s3].neighbor[mu] ;
    const size_t s5 = lat[s4].neighbor[nu] ;
    const size_t s6 = lat[i].neighbor[mu] ;
    multab_suNC( v , u , lat[s4].O[nu] ) ;
    multab_suNC( u , v , lat[s5].O[nu] ) ;
    multab_dag_suNC( v , u , lat[s6].O[mu] ) ;
    multab_dag_suNC( u , v , lat[i].O[mu] ) ;
#ifdef CLOVER_LOG_DEF
    exact_log_slow( v , u ) ;
    a_plus_Sxb( sum , v , Clover_k2 ) ;
#else
    a_plus_Sxb( sum , u , Clover_k2 ) ;
#endif
  }
  // 1x2 contrib
  if( likely( fClover_k3 > CL_TOL ) ) {
    // 1x2
    const size_t s1 = lat[i].back[nu] ;
    const size_t s2 = lat[s1].back[nu] ;
    multab_dagdag_suNC( u , lat[s1].O[nu] , lat[s2].O[nu] ) ;
    multab_suNC( v , u , lat[s2].O[mu] ) ;
    const size_t s3 = lat[s2].neighbor[mu] ;
    const size_t s4 = lat[s3].neighbor[nu] ;
    multab_suNC( u , v , lat[s3].O[nu] ) ; 
    multab_suNC( v , u , lat[s4].O[nu] ) ; 
    multab_dag_suNC( u , v , lat[i].O[mu] ) ; 
#ifdef CLOVER_LOG_DEF
    exact_log_slow( v , u ) ;
    a_plus_Sxb( sum , v , Clover_k3 ) ;
#else
    a_plus_Sxb( sum , u , Clover_k3 ) ;
#endif
    // 2x1 
    const size_t t1 = lat[i].back[nu] ;
    const size_t t2 = lat[t1].neighbor[mu] ;
    multabdag_suNC( u , lat[t1].O[nu] , lat[t1].O[mu] ) ;
    multab_suNC( v , u , lat[t2].O[mu] ) ;
    const size_t t3 = lat[t2].neighbor[mu] ;
    const size_t t4 = lat[t3].neighbor[nu] ;
    const size_t t5 = lat[t4].back[mu] ;
    multab_suNC( u , v , lat[t3].O[nu] ) ; 
    multab_dag_suNC( v , u , lat[t5].O[mu] ) ; 
    multab_dag_suNC( u , v , lat[i].O[mu] ) ; 
#ifdef CLOVER_LOG_DEF
    exact_log_slow( v , u ) ;
    a_plus_Sxb( sum , v , Clover_k3 ) ;
#else
    a_plus_Sxb( sum , u , Clover_k3 ) ;
#endif
    }
  // 1x3 contribution
  if( likely( fClover_k4 > CL_TOL ) ) {
    // 1x3
    const size_t s1 = lat[i].back[nu] ;
    const size_t s2 = lat[s1].back[nu] ;
    const size_t s3 = lat[s2].back[nu] ;
    multab_dagdag_suNC( u , lat[s1].O[nu] , lat[s2].O[nu] ) ;
    multab_dag_suNC( v , u , lat[s3].O[nu] ) ;
    multab_suNC( u , v , lat[s3].O[mu] ) ;
    const size_t s4 = lat[s3].neighbor[mu] ;
    const size_t s5 = lat[s4].neighbor[nu] ;
    const size_t s6 = lat[s5].neighbor[nu] ;
    multab_suNC( v , u , lat[s4].O[nu] ) ; 
    multab_suNC( u , v , lat[s5].O[nu] ) ; 
    multab_suNC( v , u , lat[s6].O[nu] ) ; 
    multab_dag_suNC( u , v , lat[i].O[mu] ) ; 
#ifdef CLOVER_LOG_DEF
    exact_log_slow( v , u ) ;
    a_plus_Sxb( sum , v , Clover_k4 ) ;
#else
    a_plus_Sxb( sum , u , Clover_k4 ) ;
#endif
    // 3x1
    const size_t t1 = lat[i].back[nu] ;
    const size_t t2 = lat[t1].neighbor[mu] ;
    const size_t t3 = lat[t2].neighbor[mu] ;
    multabdag_suNC( u , lat[t1].O[nu] , lat[t1].O[mu] ) ;
    multab_suNC( v , u , lat[t2].O[mu] ) ;
    multab_suNC( u , v , lat[t3].O[mu] ) ;
    const size_t t4 = lat[t3].neighbor[mu] ;
    const size_t t5 = lat[t4].neighbor[nu] ;
    const size_t t6 = lat[t5].back[mu] ;
    const size_t t7 = lat[t6].back[mu] ;
    multab_suNC( v , u , lat[t4].O[nu] ) ; 
    multab_dag_suNC( u , v , lat[t6].O[mu] ) ; 
    multab_dag_suNC( v , u , lat[t7].O[mu] ) ; 
    multab_dag_suNC( u , v , lat[i].O[mu] ) ; 
#ifdef CLOVER_LOG_DEF
    exact_log_slow( v , u ) ;
    a_plus_Sxb( sum , v , Clover_k4 ) ;
#else
    a_plus_Sxb( sum , u , Clover_k4 ) ;
#endif
  }
  // 3x3 contribution
#ifndef NK5
  if( fClover_k5 > CL_TOL ) {
    GLU_complex x[ NCNC ] , y[ NCNC ] ;
    GLU_complex s[ NCNC ] , t[ NCNC ] ;
    const size_t s1 = lat[i].back[nu] ;
    const size_t s2 = lat[s1].back[nu] ;
    const size_t s3 = lat[s2].back[nu] ;
    const size_t s4 = lat[s3].neighbor[mu] ;
    const size_t s5 = lat[s4].neighbor[mu] ;
    multab_dagdag_suNC( s , lat[s1].O[nu] , lat[s2].O[nu] ) ;
    multabdag_suNC( v , lat[s3].O[nu] , lat[s3].O[mu] ) ;
    multab_suNC( x , lat[s4].O[mu] , lat[s5].O[mu] ) ;
    multab_suNC( y , s , v ) ;
    multab_suNC( s , y , x ) ;
    // go to s6
    const size_t s6 = lat[i].neighbor[mu] ;
    const size_t s7 = lat[s6].neighbor[mu] ;
    const size_t s8 = lat[s7].neighbor[mu] ;
    const size_t s9 = lat[s8].back[nu] ;
    const size_t s10 = lat[s9].back[nu] ;
    const size_t s11 = lat[s10].back[nu] ;
    multab_suNC( t , lat[i].O[mu] , lat[s6].O[mu] ) ;
    multab_dag_suNC( v , lat[s7].O[mu] , lat[s9].O[nu] ) ;
    multab_dagdag_suNC( x , lat[s10].O[nu] , lat[s11].O[nu] ) ;
    multab_suNC( y , t , v ) ;
    multab_suNC( t , y , x ) ;
    multab_dag_suNC( u , s , t ) ;
#ifdef CLOVER_LOG_DEF
    exact_log_slow( v , u ) ;
    a_plus_Sxb( sum , v , k5 ) ;
#else
    a_plus_Sxb( sum , u , k5 ) ;
#endif
  }
#endif
  return ;
}

// Sector 3 is the
// bottom left of the clover term.
static void
compute_clover_s3( GLU_complex *__restrict sum ,
		   GLU_complex u[ NCNC ] ,
		   GLU_complex v[ NCNC ] ,
		   const struct site *__restrict lat ,
		   const size_t i ,
		   const size_t mu ,
		   const size_t nu )
{
  // 1x1 contribution
  if( likely( fClover_k1 > CL_TOL ) ) {
    compute_s3( sum , u , v , lat , i , mu , nu , Clover_k1 ) ;
  }
  // 2x2 contrib
  if( likely( fClover_k2 > CL_TOL ) ) {
    const size_t s1 = lat[i].back[mu] ;
    const size_t s2 = lat[s1].back[mu] ;
    const size_t s3 = lat[s2].back[nu] ;
    const size_t s4 = lat[s3].back[nu] ;
    multab_dagdag_suNC( u , lat[s1].O[mu] , lat[s2].O[mu] ) ;
    multab_dag_suNC( v , u , lat[s3].O[nu] ) ;
    multab_dag_suNC( u , v , lat[s4].O[nu] ) ;
    const size_t s5 = lat[s4].neighbor[mu] ;
    const size_t s6 = lat[s5].neighbor[mu] ;
    const size_t s7 = lat[s6].neighbor[nu] ;
    multab_suNC( v , u , lat[s4].O[mu] ) ;
    multab_suNC( u , v , lat[s5].O[mu] ) ;
    multab_suNC( v , u , lat[s6].O[nu] ) ;
    multab_suNC( u , v , lat[s7].O[nu] ) ;
#ifdef CLOVER_LOG_DEF
    exact_log_slow( v , u ) ;
    a_plus_Sxb( sum , v , Clover_k2 ) ;
#else
    a_plus_Sxb( sum , u , Clover_k2 ) ;
#endif
  }
  // 1x2 contrib
  if( likely( fClover_k3 > CL_TOL ) ) {
    const size_t s1 = lat[i].back[mu] ;
    const size_t s2 = lat[s1].back[nu] ;
    const size_t s3 = lat[s2].back[nu] ;
    multab_dagdag_suNC( u , lat[s1].O[mu] , lat[s2].O[nu] ) ;
    multab_dag_suNC( v , u , lat[s3].O[nu] ) ;
    const size_t s4 = lat[s3].neighbor[mu] ;
    const size_t s5 = lat[s4].neighbor[nu] ;
    multab_suNC( u , v , lat[s3].O[mu] ) ;
    multab_suNC( v , u , lat[s4].O[nu] ) ;
    multab_suNC( u , v , lat[s5].O[nu] ) ;
#ifdef CLOVER_LOG_DEF
    exact_log_slow( v , u ) ;
    a_plus_Sxb( sum , v , Clover_k3 ) ;
#else
    a_plus_Sxb( sum , u , Clover_k3 ) ;
#endif
    // 2x1
    const size_t t1 = lat[i].back[mu] ;
    const size_t t2 = lat[t1].back[mu] ;
    const size_t t3 = lat[t2].back[nu] ;
    multab_dagdag_suNC( u , lat[t1].O[mu] , lat[t2].O[mu] ) ;
    multab_dag_suNC( v , u , lat[t3].O[nu] ) ;
    const size_t t4 = lat[t3].neighbor[mu] ;
    const size_t t5 = lat[t4].neighbor[mu] ;
    multab_suNC( u , v , lat[t3].O[mu] ) ;
    multab_suNC( v , u , lat[t4].O[mu] ) ;
    multab_suNC( u , v , lat[t5].O[nu] ) ;
#ifdef CLOVER_LOG_DEF
    exact_log_slow( v , u ) ;
    a_plus_Sxb( sum , v , Clover_k3 ) ;
#else
    a_plus_Sxb( sum , u , Clover_k3 ) ;
#endif
  }
  // 1x3 contrib
  if( likely( fClover_k4 > CL_TOL ) ) {
    // 1x3
    const size_t s1 = lat[i].back[mu] ;
    const size_t s2 = lat[s1].back[nu] ;
    const size_t s3 = lat[s2].back[nu] ;
    const size_t s4 = lat[s3].back[nu] ;
    multab_dagdag_suNC( u , lat[s1].O[mu] , lat[s2].O[nu] ) ;
    multab_dag_suNC( v , u , lat[s3].O[nu] ) ;
    multab_dag_suNC( u , v , lat[s4].O[nu] ) ;
    const size_t s5 = lat[s4].neighbor[mu] ;
    const size_t s6 = lat[s5].neighbor[nu] ;
    const size_t s7 = lat[s6].neighbor[nu] ;
    multab_suNC( v , u , lat[s4].O[mu] ) ;
    multab_suNC( u , v , lat[s5].O[nu] ) ;
    multab_suNC( v , u , lat[s6].O[nu] ) ;
    multab_suNC( u , v , lat[s7].O[nu] ) ;
#ifdef CLOVER_LOG_DEF
    exact_log_slow( v , u ) ;
    a_plus_Sxb( sum , v , Clover_k4 ) ;
#else
    a_plus_Sxb( sum , u , Clover_k4 ) ;
#endif
    // 3x1
    const size_t t1 = lat[i].back[mu] ;
    const size_t t2 = lat[t1].back[mu] ;
    const size_t t3 = lat[t2].back[mu] ;
    const size_t t4 = lat[t3].back[nu] ;
    multab_dagdag_suNC( u , lat[t1].O[mu] , lat[t2].O[mu] ) ;
    multab_dag_suNC( v , u , lat[t3].O[mu] ) ;
    multab_dag_suNC( u , v , lat[t4].O[nu] ) ;
    const size_t t5 = lat[t4].neighbor[mu] ;
    const size_t t6 = lat[t5].neighbor[mu] ;
    const size_t t7 = lat[t6].neighbor[mu] ;
    multab_suNC( v , u , lat[t4].O[mu] ) ;
    multab_suNC( u , v , lat[t5].O[mu] ) ;
    multab_suNC( v , u , lat[t6].O[mu] ) ;
    multab_suNC( u , v , lat[t7].O[nu] ) ;
#ifdef CLOVER_LOG_DEF
    exact_log_slow( v , u ) ;
    a_plus_Sxb( sum , v , Clover_k4 ) ;
#else
    a_plus_Sxb( sum , u , Clover_k4 ) ;
#endif
  }
#ifndef NK5
  // 3x3 contribution
  if( fClover_k5 > CL_TOL ) {
    GLU_complex x[ NCNC ] , y[ NCNC ] ;
    GLU_complex s[ NCNC ] , t[ NCNC ] ;
    const size_t s1 = lat[i].back[mu] ;
    const size_t s2 = lat[s1].back[mu] ;
    const size_t s3 = lat[s2].back[mu] ;
    const size_t s4 = lat[s3].back[nu] ;
    const size_t s5 = lat[s4].back[nu] ;
    const size_t s6 = lat[s5].back[nu] ;
    multab_dagdag_suNC( s , lat[s1].O[mu] , lat[s2].O[mu] ) ;
    multab_dagdag_suNC( v , lat[s3].O[mu] , lat[s4].O[nu] ) ;
    multab_dagdag_suNC( x , lat[s5].O[nu] , lat[s6].O[nu] ) ;
    multab_suNC( y , s , v ) ;
    multab_suNC( s , y , x ) ;
    const size_t s7 = lat[s6].neighbor[mu] ;
    const size_t s8 = lat[s7].neighbor[mu] ;
    const size_t s9 = lat[s8].neighbor[mu] ;
    const size_t s10 = lat[s9].neighbor[nu] ;
    const size_t s11 = lat[s10].neighbor[nu] ;
    multab_suNC( t , lat[s6].O[mu] , lat[s7].O[mu] ) ;
    multab_suNC( v , lat[s8].O[mu] , lat[s9].O[nu] ) ;
    multab_suNC( x , lat[s10].O[nu] , lat[s11].O[nu] ) ;
    multab_suNC( y , t , v ) ;
    multab_suNC( t , y , x ) ;
    multab_suNC( u , s , t ) ;
#ifdef CLOVER_LOG_DEF
    exact_log_slow( v , u ) ;
    a_plus_Sxb( sum , v , k5 ) ;
#else
    a_plus_Sxb( sum , u , k5 ) ;
#endif
  }
#endif
  return ;
}

// Sector 4 is the
// top left of the clover term.
static void
compute_clover_s4( GLU_complex *__restrict sum ,
		   GLU_complex u[ NCNC ] ,
		   GLU_complex v[ NCNC ] ,
		   const struct site *__restrict lat ,
		   const size_t i ,
		   const size_t mu ,
		   const size_t nu )
{
  // 1x1 contribution
  if( likely( fClover_k1 > CL_TOL ) ) {
    compute_s4( sum , u , v , lat , i , mu , nu , Clover_k1 ) ;
  }
  // 2x2 contrib
  if( likely( fClover_k2 > CL_TOL ) ) {
    const size_t s1 = lat[i].neighbor[nu] ;
    const size_t s2 = lat[s1].neighbor[nu] ;
    const size_t s3 = lat[s2].back[mu] ;
    const size_t s4 = lat[s3].back[mu] ;
    multab_suNC( u , lat[i].O[nu] , lat[s1].O[nu] ) ;
    multab_dag_suNC( v , u , lat[s3].O[mu] ) ;
    multab_dag_suNC( u , v , lat[s4].O[mu] ) ;
    const size_t s5= lat[s4].back[nu] ;
    const size_t s6 = lat[s5].back[nu] ;
    const size_t s7 = lat[s6].neighbor[mu] ;
    multab_dag_suNC( v , u , lat[s5].O[nu] ) ;
    multab_dag_suNC( u , v , lat[s6].O[nu] ) ;
    multab_suNC( v , u , lat[s6].O[mu] ) ;
    multab_suNC( u , v , lat[s7].O[mu] ) ;
#ifdef CLOVER_LOG_DEF
    exact_log_slow( v , u ) ;
    a_plus_Sxb( sum , v , Clover_k2 ) ;
#else
    a_plus_Sxb( sum , u , Clover_k2 ) ;
#endif
  }
  // 1x2 contrib
  if( likely( fClover_k3 > CL_TOL ) ) {
    const size_t s1 = lat[i].neighbor[nu] ;
    const size_t s2 = lat[s1].neighbor[nu] ;
    const size_t s3 = lat[s2].back[mu] ;
    multab_suNC( u , lat[i].O[nu] , lat[s1].O[nu] ) ;
    multab_dag_suNC( v , u , lat[s3].O[mu] ) ;
    const size_t s4 = lat[s3].back[nu] ;
    const size_t s5 = lat[s4].back[nu] ;
    multab_dag_suNC( u , v , lat[s4].O[nu] ) ;
    multab_dag_suNC( v , u , lat[s5].O[nu] ) ;
    multab_suNC( u , v , lat[s5].O[mu] ) ;
#ifdef CLOVER_LOG_DEF
    exact_log_slow( v , u ) ;
    a_plus_Sxb( sum , v , Clover_k3 ) ;
#else
    a_plus_Sxb( sum , u , Clover_k3 ) ;
#endif
    // 2x1 contrib
    const size_t t1 = lat[i].neighbor[nu] ;
    const size_t t2 = lat[t1].back[mu] ;
    const size_t t3 = lat[t2].back[mu] ;
    multab_dag_suNC( u , lat[i].O[nu] , lat[t2].O[mu] ) ;
    multab_dag_suNC( v , u , lat[t3].O[mu] ) ;
    const size_t t4 = lat[t3].back[nu] ;
    const size_t t5 = lat[t4].neighbor[mu] ;
    multab_dag_suNC( u , v , lat[t4].O[nu] ) ;
    multab_suNC( v , u , lat[t4].O[mu] ) ;
    multab_suNC( u , v , lat[t5].O[mu] ) ;
#ifdef CLOVER_LOG_DEF
    exact_log_slow( v , u ) ;
    a_plus_Sxb( sum , v , Clover_k3 ) ;
#else
    a_plus_Sxb( sum , u , Clover_k3 ) ;
#endif
  }
  // 1x3 contrib
  if( likely( fClover_k4 > CL_TOL ) ) {
    // 1x3
    const size_t s1 = lat[i].neighbor[nu] ;
    const size_t s2 = lat[s1].neighbor[nu] ;
    const size_t s3 = lat[s2].neighbor[nu] ;
    const size_t s4 = lat[s3].back[mu] ;
    multab_suNC( u , lat[i].O[nu] , lat[s1].O[nu] ) ;
    multab_suNC( v , u , lat[s2].O[nu] ) ;
    multab_dag_suNC( u , v , lat[s4].O[mu] ) ;
    const size_t s5 = lat[s4].back[nu] ;
    const size_t s6 = lat[s5].back[nu] ;
    const size_t s7 = lat[s6].back[nu] ;
    multab_dag_suNC( v , u , lat[s5].O[nu] ) ;
    multab_dag_suNC( u , v , lat[s6].O[nu] ) ;
    multab_dag_suNC( v , u , lat[s7].O[nu] ) ;
    multab_suNC( u , v , lat[s7].O[mu] ) ;
#ifdef CLOVER_LOG_DEF
    exact_log_slow( v , u ) ;
    a_plus_Sxb( sum , v , Clover_k4 ) ;
#else
    a_plus_Sxb( sum , u , Clover_k4 ) ;
#endif
    // 3x1
    const size_t t1 = lat[i].neighbor[nu] ;
    const size_t t2 = lat[t1].back[mu] ;
    const size_t t3 = lat[t2].back[mu] ;
    const size_t t4 = lat[t3].back[mu] ;
    multab_dag_suNC( u , lat[i].O[nu] , lat[t2].O[mu] ) ;
    multab_dag_suNC( v , u , lat[t3].O[mu] ) ;
    multab_dag_suNC( u , v , lat[t4].O[mu] ) ;
    const size_t t5 = lat[t4].back[nu] ;
    const size_t t6 = lat[t5].neighbor[mu] ;
    const size_t t7 = lat[t6].neighbor[mu] ;
    multab_dag_suNC( v , u , lat[t5].O[nu] ) ;
    multab_suNC( u , v , lat[t5].O[mu] ) ;
    multab_suNC( v , u , lat[t6].O[mu] ) ;
    multab_suNC( u , v , lat[t7].O[mu] ) ;
#ifdef CLOVER_LOG_DEF
    exact_log_slow( v , u ) ;
    a_plus_Sxb( sum , v , Clover_k4 ) ;
#else
    a_plus_Sxb( sum , u , Clover_k4 ) ;
#endif
  }
#ifndef NK5
  // 3x3 contribution
  if( fClover_k5 > CL_TOL ) {
    GLU_complex x[ NCNC ] , y[ NCNC ] ;
    GLU_complex s[ NCNC ] , t[ NCNC ] ;
    const size_t s1 = lat[i].neighbor[nu] ;
    const size_t s2 = lat[s1].neighbor[nu] ;
    const size_t s3 = lat[s2].neighbor[nu] ;
    const size_t s4 = lat[s3].back[mu] ;
    const size_t s5 = lat[s4].back[mu] ;
    const size_t s6 = lat[s5].back[mu] ;
    multab_suNC( s , lat[i].O[nu] , lat[s1].O[nu] ) ;
    multab_dag_suNC( v , lat[s2].O[nu] , lat[s4].O[mu] ) ;
    multab_dagdag_suNC( x , lat[s5].O[mu] , lat[s6].O[mu] ) ;
    multab_suNC( y , s , v ) ;
    multab_suNC( s , y , x ) ;
    const size_t s7 = lat[s6].back[nu] ;
    const size_t s8 = lat[s7].back[nu] ;
    const size_t s9 = lat[s8].back[nu] ;
    const size_t s10 = lat[s9].neighbor[mu] ;
    const size_t s11 = lat[s10].neighbor[mu] ;
    multab_dagdag_suNC( t , lat[s7].O[nu] , lat[s8].O[nu] ) ;
    multabdag_suNC( v , lat[s9].O[nu] , lat[s9].O[mu] ) ;
    multab_suNC( x , lat[s10].O[mu] , lat[s11].O[mu] ) ;
    multab_suNC( y , t , v ) ;
    multab_suNC( t , y , x ) ;
    multab_suNC( u , s , t ) ;
#ifdef CLOVER_LOG_DEF
    exact_log_slow( v , u ) ;
    a_plus_Sxb( sum , v , k5 ) ;
#else
    a_plus_Sxb( sum , u , k5 ) ;
#endif
  }
#endif
  return ;
}
#endif

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//          After the definitions of the clovers we just                 //
//      need to add them and take the traceless hermitian projection.    //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

// compute the field strength tensor and the topological charge
static void
compute_GG_q( const GLU_complex sum_1[ NCNC ] ,
	      const GLU_complex sum_2[ NCNC ] ,
	      double *__restrict plaq_t , 
	      double *__restrict plaq_sp , 
	      double *__restrict Q )
{
  // take the antihermitian projection, could perform the exact log here?
  GLU_complex GMUNU_1[ NCNC ] GLUalign ;
  GLU_complex GMUNU_2[ NCNC ] GLUalign ;
  #ifdef CLOVER_LOG_DEF
  equiv( GMUNU_1 , sum_1 ) ;
  equiv( GMUNU_2 , sum_2 ) ;
  #elif defined ANTIHERMITIAN
  AntiHermitian_proj( GMUNU_1 , sum_1 ) ;
  AntiHermitian_proj( GMUNU_2 , sum_2 ) ;
  #elif defined TRF_ANTIHERMITIAN
  trf_AntiHermitian_proj( GMUNU_1 , sum_1 ) ;
  trf_AntiHermitian_proj( GMUNU_2 , sum_2 ) ; 
  #else
  Hermitian_proj( GMUNU_1 , sum_1 ) ;
  Hermitian_proj( GMUNU_2 , sum_2 ) ;
  #endif

  GLU_complex tr_sp = 0. ;
  trace_ab( &tr_sp , GMUNU_1 , GMUNU_1 ) ;
  #ifdef ANTIHERMITIAN
  *plaq_sp = -(double)creal( tr_sp ) ;
  #elif defined TRF_ANTIHERMITIAN
  *plaq_sp = -(double)creal( tr_sp ) ; 
  #else
  *plaq_sp = (double)creal( tr_sp ) ;
  #endif

  GLU_complex tr_t = 0. ;
  trace_ab( &tr_t , GMUNU_2 , GMUNU_2 ) ;
  #ifdef ANTIHERMITIAN
  *plaq_t = -(double)creal( tr_t ) ;
  #elif defined TRF_ANTIHERMITIAN
  *plaq_t = -(double)creal( tr_t ) ;
  #else
  *plaq_t = (double)creal( tr_t ) ;
  #endif

  GLU_complex tr_q = 0. ;
  trace_ab( &tr_q , GMUNU_1 , GMUNU_2 ) ;
  #ifdef ANTIHERMITIAN
  *Q = -(double)creal( tr_q ) ;
  #elif defined TRF_ANTIHERMITIAN
  *Q = -(double)creal( tr_q ) ;
  #else
  *Q = (double)creal( tr_q ) ;
  #endif
  return ;
}

#endif // endif for ND = 4

// computes the topological charge at site i for directions
// mu,nu,rho and delta
static void
compute_q( GLU_complex q[ NCNC ] ,
	   const struct site *__restrict lat ,
	   const size_t i ,
	   const size_t mu , 
	   const size_t nu , 
	   const size_t rho , 
	   const size_t delta )
{
#if ND == 4 
  // accumulated mu-nu in sum_1 and rho-delta in sum_2
  GLU_complex sum_1[ NCNC ] GLUalign , sum_2[ NCNC ] GLUalign ;
  zero_mat( sum_1 ) ; zero_mat( sum_2 ) ;
  // temp matrices in u and v
  GLU_complex u[ NCNC ] GLUalign , v[ NCNC ] GLUalign ;
#ifdef CLOVER_IMPROVE
  // highly improved clover definition
  compute_clover_s1( sum_1 , u , v , lat , i , mu , nu ) ;
  compute_clover_s2( sum_1 , u , v , lat , i , mu , nu ) ;
  compute_clover_s3( sum_1 , u , v , lat , i , mu , nu ) ;
  compute_clover_s4( sum_1 , u , v , lat , i , mu , nu ) ;
  compute_clover_s1( sum_2 , u , v , lat , i , rho , delta ) ;
  compute_clover_s2( sum_2 , u , v , lat , i , rho , delta ) ;
  compute_clover_s3( sum_2 , u , v , lat , i , rho , delta ) ;
  compute_clover_s4( sum_2 , u , v , lat , i , rho , delta ) ;
#elif defined PLAQUETTE_FMUNU
  // plaquette definition
  compute_s1( sum_1 , u , v , lat , i , mu , nu , 1.0 ) ;
  compute_s1( sum_2 , u , v , lat , i , rho , delta , 1.0 ) ;
#else
  // standard clover plaquette definition
  compute_s1( sum_1 , u , v , lat , i , mu , nu , 1.0 ) ;
  compute_s2( sum_1 , u , v , lat , i , mu , nu , 1.0 ) ;
  compute_s3( sum_1 , u , v , lat , i , mu , nu , 1.0 ) ;
  compute_s4( sum_1 , u , v , lat , i , mu , nu , 1.0 ) ;
  compute_s1( sum_2 , u , v , lat , i , rho , delta , 1.0 ) ;
  compute_s2( sum_2 , u , v , lat , i , rho , delta , 1.0 ) ;
  compute_s3( sum_2 , u , v , lat , i , rho , delta , 1.0 ) ;
  compute_s4( sum_2 , u , v , lat , i , rho , delta , 1.0 ) ;
#endif
  // take the antihermitian projection, could perform the exact log here?
  GLU_complex GMUNU_1[ NCNC ] GLUalign ;
  GLU_complex GMUNU_2[ NCNC ] GLUalign ;
  #ifdef CLOVER_LOG_DEF
  equiv( GMUNU_1 , sum_1 ) ;
  equiv( GMUNU_2 , sum_2 ) ;
  #elif defined ANTIHERMITIAN
  AntiHermitian_proj( GMUNU_1 , sum_1 ) ;
  AntiHermitian_proj( GMUNU_2 , sum_2 ) ;
  #elif defined TRF_ANTIHERMITIAN
  trf_AntiHermitian_proj( GMUNU_1 , sum_1 ) ;
  trf_AntiHermitian_proj( GMUNU_2 , sum_2 ) ; 
  #else
  Hermitian_proj( GMUNU_1 , sum_1 ) ;
  Hermitian_proj( GMUNU_2 , sum_2 ) ;
  #endif
  // compute the product of G_MUNU_1 with G_MUNU_2
  multab( q , GMUNU_1 , GMUNU_2 ) ;
#endif
  return ;
}

// kernel code 
static void
compute_Gmunu_kernel( double *__restrict plaq_t ,
		      double *__restrict plaq_sp ,
		      double *__restrict qtop ,
		      const struct site *__restrict lat ,
		      const size_t i ,
		      const size_t mu , 
		      const size_t nu , 
		      const size_t rho , 
		      const size_t delta )
{
#if ND == 4
  // accumulated mu-nu in sum_1 and rho-delta in sum_2
  GLU_complex sum_1[ NCNC ] GLUalign , sum_2[ NCNC ] GLUalign ;
  zero_mat( sum_1 ) ; zero_mat( sum_2 ) ;
  // temp matrices in u and v
  GLU_complex u[ NCNC ] GLUalign , v[ NCNC ] GLUalign ;
#ifdef CLOVER_IMPROVE
  // highly improved clover definition
  compute_clover_s1( sum_1 , u , v , lat , i , mu , nu ) ;
  compute_clover_s2( sum_1 , u , v , lat , i , mu , nu ) ;
  compute_clover_s3( sum_1 , u , v , lat , i , mu , nu ) ;
  compute_clover_s4( sum_1 , u , v , lat , i , mu , nu ) ;
  compute_clover_s1( sum_2 , u , v , lat , i , rho , delta ) ;
  compute_clover_s2( sum_2 , u , v , lat , i , rho , delta ) ;
  compute_clover_s3( sum_2 , u , v , lat , i , rho , delta ) ;
  compute_clover_s4( sum_2 , u , v , lat , i , rho , delta ) ;
#elif defined PLAQUETTE_FMUNU
  // plaquette definition
  compute_s1( sum_1 , u , v , lat , i , mu , nu , 1.0 ) ;
  compute_s1( sum_2 , u , v , lat , i , rho , delta , 1.0 ) ;
#else
  // standard clover plaquette definition
  compute_s1( sum_1 , u , v , lat , i , mu , nu , 1.0 ) ;
  compute_s2( sum_1 , u , v , lat , i , mu , nu , 1.0 ) ;
  compute_s3( sum_1 , u , v , lat , i , mu , nu , 1.0 ) ;
  compute_s4( sum_1 , u , v , lat , i , mu , nu , 1.0 ) ;
  compute_s1( sum_2 , u , v , lat , i , rho , delta , 1.0 ) ;
  compute_s2( sum_2 , u , v , lat , i , rho , delta , 1.0 ) ;
  compute_s3( sum_2 , u , v , lat , i , rho , delta , 1.0 ) ;
  compute_s4( sum_2 , u , v , lat , i , rho , delta , 1.0 ) ;
#endif
  // and we are done, perform the measurement
  compute_GG_q( sum_1 , sum_2 , plaq_t , plaq_sp , qtop ) ;
#ifdef PLAQUETTE_FMUNU
  *plaq_sp = 2.0 * ( (double)NC - (double)trace( sum_1 ) ) ;
  *plaq_t = 2.0 * ( (double)NC - (double)trace( sum_2 ) ) ;
#endif
#endif
  return ;
}

// this is the driving code for the computation
void
compute_Gmunu( double *__restrict GG , 
	       double *__restrict qtop ,
	       const struct site *__restrict lat )
{
  // initialise the temporal and spatial plaq and the qtop
  double plaq_sp = 0. , plaq_t = 0. , Q = 0. ;
  size_t i ;
  // control factors
#ifndef NK5
  fClover_k1 = fabs( Clover_k1 ) ;
  fClover_k2 = fabs( Clover_k2 ) ;
  fClover_k3 = fabs( Clover_k3 ) ;
  fClover_k4 = fabs( Clover_k4 ) ;
  fClover_k5 = fabs( k5 ) ;
#endif

#pragma omp parallel for private(i) reduction(+:plaq_sp) reduction(+:plaq_t) reduction(+:Q)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    double plsp , plt , q ;
    compute_Gmunu_kernel( &plt , &plsp , &q , lat ,
			  i , 3 , 0 , 1 , 2  ) ;
    plaq_sp = plaq_sp + plsp ;
    plaq_t = plaq_t + plt ;
    Q = Q + q ;
    //
    compute_Gmunu_kernel( &plt , &plsp , &q , lat ,
			  i , 3 , 1 , 2 , 0  ) ;
    plaq_sp = plaq_sp + plsp ;
    plaq_t = plaq_t + plt ;
    Q = Q + q ;
    //
    compute_Gmunu_kernel( &plt , &plsp , &q , lat ,
			  i , 3 , 2 , 0 , 1  ) ;
    plaq_sp = plaq_sp + plsp ;
    plaq_t = plaq_t + plt ;
    Q = Q + q ;
  }
  *GG = ( plaq_sp + plaq_t ) ;
  *qtop = Q ;
  // just to accommodate for the sum over 4 of the others
#ifdef PLAQUETTE_FMUNU
  *GG *= 16.0 ;
  *qtop *= Q ;
#endif

  return ;
}

// this is the driving code for the computation
void
compute_Gmunu_array( GLU_complex *__restrict qtop , // an LVOLUME array for the qtop
		     const struct site *__restrict lat )
{
  size_t i ;
  // control factors
#ifndef NK5
  fClover_k1 = fabs( Clover_k1 ) ;
  fClover_k2 = fabs( Clover_k2 ) ;
  fClover_k3 = fabs( Clover_k3 ) ;
  fClover_k4 = fabs( Clover_k4 ) ;
  fClover_k5 = fabs( k5 ) ;
#endif

#pragma omp parallel for private(i)
  for( i = 0 ; i < LVOLUME ; i++ ) {
    GLU_complex Qmat[ NCNC ] GLUalign ;
    compute_q( Qmat , lat , i , 3 , 0 , 1 , 2  ) ;
    qtop[i] = trace( Qmat ) ;
    //
    compute_q( Qmat , lat , i , 3 , 1 , 2 , 0  ) ;
    qtop[i] += trace( Qmat ) ;
    //
    compute_q( Qmat , lat , i , 3 , 2 , 0 , 1  ) ;
    qtop[i] += trace( Qmat ) ;
  }
  return ;
}

// undefine all of the field defs
#ifdef TRF_ANTIHERMITIAN
  #undef TRF_ANTIHERMITIAN
#endif
#ifdef ANTIHERMITIAN
  #undef ANTIHERMITIAN
#endif
#ifdef CLOVER_LOG_DEF
  #undef CLOVER_LOG_DEF
#endif
