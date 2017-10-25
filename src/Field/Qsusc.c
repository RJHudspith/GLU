/*
    Copyright 2013-2016 Renwick James Hudspith

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
   @file Qsusc.c
   @brief computation of topological things via a wrapper
   @warning calls the smearing wrapper, which overwrites the links
 */
#include "Mainfile.h"

#include "cut_output.h" // write_qmoments()
#include "Qslab.h"      // compute_Qsusc()
#include "Qcorr.h"      // compute_Qcorr()
#include "Qmoments.h"   // compute_Qmoments()
#include "SM_wrap.h"    // smearing wrapper

// step through smears
int
compute_Qsusc_step( struct site *__restrict lat ,
		    const struct cut_info CUTINFO ,
		    const struct sm_info SMINFO )
{
  // measurement counter
  size_t measurement ;
  int flag = GLU_SUCCESS ;
  
  fprintf( stdout , "[QSUSC] performing %zu measurements at %zu" 
	   " smearing steps\n" , CUTINFO.max_t , SMINFO.smiters ) ;

  // allocations for the moments
  struct Qmoments *Qmom = NULL ;
  if( CUTINFO.dir == TOPOLOGICAL_MOMENTS ) {
    Qmom = malloc( ( CUTINFO.max_t - 1 ) * sizeof( struct Qmoments ) ) ;
    size_t t ;
    for( t = 0 ; t < CUTINFO.max_t-1 ; t++ ) {
      Qmom[t].Q  = malloc( ND * NQMOMENTS * sizeof( double ) ) ;
      Qmom[t].Q2 = malloc( ND * NQMOMENTS * sizeof( double ) ) ;
    }
  }

  // compute the topological correlator in r and the temporal correlator
  if( CUTINFO.dir == TOPOLOGICAL_SUSCEPTIBILITY ) {
    if( compute_Qsusc( lat , CUTINFO , 0 ) == GLU_FAILURE ) {
      flag = GLU_FAILURE ; goto memfree ;
    }
  }

  // perform some measurements each stepping smiters ahead
  for( measurement = 1 ; measurement < CUTINFO.max_t ; measurement++ ) {

    // do some smearing ...
    SM_wrap_struct( lat , SMINFO ) ;
    
    // compute the topological correlator in r and the temporal correlator
    if( CUTINFO.dir == TOPOLOGICAL_SUSCEPTIBILITY ) {
      if( compute_Qsusc( lat , CUTINFO , measurement ) == GLU_FAILURE ) {
        flag = GLU_FAILURE ; break ;
      }
      // compute the correlator
    } else if( CUTINFO.dir == TOPOLOGICAL_CORRELATOR ) {
      if( compute_Qcorr( lat , CUTINFO , measurement ) == GLU_FAILURE ) {
	flag = GLU_FAILURE ; break ;
      }
    } else if( CUTINFO.dir == TOPOLOGICAL_MOMENTS ) {
      if( compute_Qmoments( lat , &Qmom[ measurement-1 ] , 
			    CUTINFO , measurement ) == GLU_FAILURE ) {
        flag = GLU_FAILURE ; break ;
      }
    }
    //
  }

 memfree :
  
  if( CUTINFO.dir == TOPOLOGICAL_MOMENTS ) {

    // write to a file
    write_moments( Qmom , CUTINFO.max_t-1 ) ;
    
    size_t t ;
    for( t = 0 ; t < CUTINFO.max_t-1 ; t++ ) {
      free( Qmom[t].Q ) ; free( Qmom[t].Q2 ) ;
    }
    free( Qmom ) ; 
  }
  
  return flag ;
}
