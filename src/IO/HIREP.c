/*
    Copyright 2013 Renwick James Hudspith

    This file (HIREP.c) is part of GLU.

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
   @file HIREP.c
   @brief code to read and write HIREP configurations

  bizarrely their geometry is
  t,x,y,z which makes little sense to me. I imagine there is a very good reason
  for this.... Probably.
 */

#include "Mainfile.h"

#include "geometry.h"    // for changing to our geometry 
#include "GLU_bswap.h"   // byte swapping arrays
#include "plaqs_links.h" // compute the plaquette

// takes a HiRep index "i" and converts to the Nersc or my index
static size_t
translate_to_GLU_geom( const size_t i )
{
  int x[ ND ] ;
  get_mom_2piBZ_hirep2( x , i ) ;
  return gen_site( x ) ;
}

// generic (ish) reader for files in the HiRep format
int
read_gauge_field( struct site *__restrict lat , 
		  FILE *__restrict in , 
		  uint32_t *chksum )
{
  double Nelements = NCNC ;
#if NC==2
  Nelements = NC ;
#endif
  const size_t stride = 2 * Nelements * ND ;
  // these guys also seem to save only in double, making things easy
  double *uind = malloc( stride * sizeof( double ) ) ; 

  uint32_t k = 0 ;
  size_t i ;
  for( i = 0 ; i < Latt.Volume ; i++ ) {
    const size_t idx = translate_to_GLU_geom( i ) ;

    if( fread( uind , sizeof( double ) , stride , in ) != stride ) {
      fprintf( stderr , "File read error.. Leaving \n " ) ;
      free( uind ) ;
      return GLU_FAILURE ;
    }
      
    if ( !WORDS_BIGENDIAN ) {
      bswap_64( stride , uind ) ; 
    }

#if NC==2
    size_t mu , a = 0 ;
    // t-direction is first
    lat[idx].O[ND-1][0] = (GLU_real)uind[a]   + I * (GLU_real)uind[a+1] ;
    lat[idx].O[ND-1][1] = (GLU_real)uind[a+2] + I * (GLU_real)uind[a+3] ;
    lat[idx].O[ND-1][2] = -conj( lat[idx].O[ND-1][1] ) ;
    lat[idx].O[ND-1][3] =  conj( lat[idx].O[ND-1][0] ) ;
    a += 2*NC ;
    // then the others
    for( mu = 0 ;  mu < ND - 1 ; mu++ ) {
      lat[idx].O[mu][0] = (GLU_real)uind[a] + I * (GLU_real)uind[a+1] ;
      lat[idx].O[mu][1] = (GLU_real)uind[a+2] + I * (GLU_real)uind[a+3] ;
      lat[idx].O[mu][2] = -conj( lat[idx].O[mu][1] ) ;
      lat[idx].O[mu][3] =  conj( lat[idx].O[mu][0] ) ;
      a += 2*NC ;
    }
#else
    size_t mu , j , a = 0 ;
    // t first
    for( j = 0 ; j < NCNC ; j++ ) {
      lat[idx].O[ND-1][j] = (GLU_real)uind[a] + I * (GLU_real)uind[a + 1] ;
      a += 2 ;
    }
    // then the others
    for( mu = 0 ;  mu < ND - 1 ; mu++ ) {
      for( j = 0 ; j < NCNC ; j++ ) {
	lat[idx].O[mu][j] = (GLU_real)uind[a] + I * (GLU_real)uind[a + 1] ;
	a += 2 ;
      }
    }
#endif
  }
  
  *chksum = k ;
  free( uind ) ;
  return GLU_SUCCESS ; 
}

// writes my gauge field out in the HiRep format 
void
write_gauge_field( const struct site *__restrict lat ,
		   FILE *__restrict outfile )
{
#if NC==2
  const size_t mat_stride = NC ;
#else
  const size_t mat_stride = NCNC ;
#endif
  const size_t stride = mat_stride * ND * 2 ; // the size of the link matrices 

  int NAV[ ND + 1 ] ; 
  NAV[ 0 ] = NC ;
  NAV[ 1 ] = Latt.dims[ ND - 1 ] ;

  size_t mu ;
  for( mu = 0 ; mu < ND - 1 ; mu++ ) {
    NAV[ 2 + mu ] = Latt.dims[ mu ] ;
  }

  // first of all set calculate the plaquette checksum in header
  double plaq[1] = { 0 } ;

  plaq[0] = av_plaquette( lat ) ;

  // write the dimensions first and then the plaquette
  // fprintf or fwrite? -- fwrite 
  #ifdef OUT_BIG
  if( !WORDS_BIGENDIAN ) {
    bswap_32( ND + 1 , NAV ) ; 
    bswap_64( 1 , plaq ) ;
  }
  #else
  if( WORDS_BIGENDIAN ) {
    bswap_32( ND + 1 , NAV ) ; 
    bswap_64( 1 , plaq ) ;
  }
  #endif

  // write the checksums
  fwrite( NAV , ( ND + 1 ) * sizeof( int ) , 1 , outfile ) ;
  fwrite( plaq , sizeof( double ) , 1 , outfile ) ;
  // now to put the gauge field in, in HiRep's weird order!

  double *uoutd = ( double* )malloc( stride *sizeof( double ) ) ; 

  size_t i ;
  for( i = 0 ; i < LVOLUME ; i++ ) {
    // loop HiRep indices
    const size_t idx = translate_to_GLU_geom( i ) ;
    size_t mu , j , a = 0 ;
    // t first
    for( j = 0 ; j < mat_stride ; j++ ) {
      uoutd[ a + 0 ] = ( double )creal( lat[idx].O[ ND - 1 ][ j ] ) ; 
      uoutd[ a + 1 ] = ( double )cimag( lat[idx].O[ ND - 1 ][ j ] ) ; 
      a += 2 ;
    }
    // then xyz
    for( mu = 0 ;  mu < ND - 1 ; mu++ ) {
      for( j = 0 ; j < mat_stride ; j++ ) {
	uoutd[ a + 0 ] = ( double )creal( lat[idx].O[mu][j] ) ; 
	uoutd[ a + 1 ] = ( double )cimag( lat[idx].O[mu][j] ) ; 
	a += 2 ;
      }
    } 
    #ifdef OUT_BIG
    if( !WORDS_BIGENDIAN ) {
      bswap_64( stride , uoutd ) ; 
    }
    #else
    if( WORDS_BIGENDIAN ) {
      bswap_64( stride , uoutd ) ;
    }
    #endif
    fwrite( uoutd , sizeof( double ) , stride , outfile ) ; 
  }

  free( uoutd ) ;

  return ;
}
