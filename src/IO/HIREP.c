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
inline static int
translate_to_GLU_geom( i )
     const int i ;
{
  int x[ ND ] ;
  get_mom_2piBZ_hirep2( x , i ) ;
  return gen_site( x ) ;
}

// generic (ish) reader for files in the HiRep format
void
read_gauge_field( struct site *__restrict lat , 
		  FILE *__restrict in , 
		  uint32_t *chksum )
{
  const int stride = 2 * NCNC * ND ;
  // these guys also seem to save only in double, making things easy
  double *uind = malloc( stride * sizeof( double ) ) ; 

  uint32_t k = 0 ;
  int i ;
  for( i = 0 ; i < Latt.Volume ; i++ ) {
    const int idx = translate_to_GLU_geom( i ) ;

    if( fread( uind , sizeof( double ) , stride , in ) != stride ) {
      printf( "File read error.. Leaving \n " ) ;
      free( uind ) ;
      return ;
    }
      
    if ( !WORDS_BIGENDIAN ) {
      bswap_64( stride , uind ) ; 
    }

    int mu , j , a = 0 ;
    for( j = 0 ; j < NCNC ; j++ ) {
      lat[idx].O[ND-1][j] = (GLU_real)uind[a] + I * (GLU_real)uind[a + 1] ;
      a += 2 ;
    }
      
    for( mu = 0 ;  mu < ND - 1 ; mu++ ) {
      for( j = 0 ; j < NCNC ; j++ ) {
	lat[idx].O[mu][j] = (GLU_real)uind[a] + I * (GLU_real)uind[a + 1] ;
	a += 2 ;
      }
    }
  }
  
  *chksum = k ;
  free( uind ) ;
  return ; 
}

// writes my gauge field out in the HiRep format 
void
write_gauge_field( const struct site *__restrict lat ,
		   FILE *__restrict outfile )
{
  const int stride = NCNC * ND * 2 ; // the size of the link matrices at a site

  int NAV[ ND + 1 ] ; 
  NAV[ 0 ] = NC ;
  NAV[ 1 ] = Latt.dims[ ND - 1 ] ;

  int mu ;
  for( mu = 0 ; mu < ND - 1 ; mu++ ) {
    NAV[ 2 + mu ] = Latt.dims[ mu ] ;
  }

  // first of all set calculate the plaquette checksum in header
  double plaq[1] ;

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

  double *uoutd = ( double* ) malloc( stride *sizeof( double ) ) ; 

  int i ;
  for( i = 0 ; i < LVOLUME ; i++ ) {
    // loop HiRep indices
    const int idx = translate_to_GLU_geom( i ) ;
    int mu , j , a = 0 ;
    for( j = 0 ; j < NCNC ; j++ ) {
      uoutd[ a ] = ( double )creal( lat[idx].O[ ND - 1 ][ j ] ) ; 
      uoutd[ a + 1 ] = ( double )cimag( lat[idx].O[ ND - 1 ][ j ] ) ; 
      a += 2 ;
    }
      
    for( mu = 0 ;  mu < ND - 1 ; mu++ ) {
      for( j = 0 ; j < NCNC ; j++ ) {
	uoutd[ a ] = ( double )creal( lat[idx].O[mu][j] ) ; 
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
