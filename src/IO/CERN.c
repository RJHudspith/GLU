/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (CERN.c) is part of GLU.

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
   @file CERN.c
   @brief code to read and write CERN configs generated by openQCD

   bizarrely their geometry is
   t,x,y,z which makes little sense to me. I imagine there is a very good reason
   for this.... Probably.
   
   @warning these aren't #NC or #ND generic. Just for SU3
 */
#include "Mainfile.h"

#include "geometry.h"      // gen_site()
#include "GLU_bswap.h"     // byte swapping arrays
#include "plaqs_links.h"   // compute the plaquette
#include "random_config.h" // latt_reunitU()

// read a CERN gauge field
int
read_CLS_field( struct site *__restrict lat , 
		FILE *__restrict in , 
		uint32_t *chksum )
{
  const size_t Nelements = NCNC ;
  const size_t PM = 2 ; // +/-
  const size_t Complex = 2 ; // double complex

  const size_t stride = PM * Complex * Nelements * ND ;
  // these guys also seem to save only in double, making things easy
  double *uind = malloc( stride * sizeof( double ) ) ; 

  uint32_t k = 0 ;
  size_t x, y,z,t ;
  for( t = 0 ; t < Latt.dims[ND-1] ; t++ ) {
    for( x = 0 ; x < Latt.dims[0] ; x++ ) {
      for( y = 0 ; y < Latt.dims[1] ; y++ ) {
	for( z = 0 ; z < Latt.dims[2] ; z++ ) {
    
	  if( (x+y+z+t)%2 ) {
	    
	    int X[4] = { (int)x , (int)y , (int)z , (int)t } ;

	    const size_t idx = gen_site( X ) ;
	    
	    if( fread( uind , sizeof( double ) , stride , in ) != stride ) {
	      fprintf( stderr , "File read error.. Leaving \n " ) ;
	      free( uind ) ;
	      return GLU_FAILURE ;
	    }
	    if( WORDS_BIGENDIAN ) {
	      bswap_64( stride , uind ) ; 
	    }
      
	    size_t mu , j , a = 0 ;
	    size_t shift = lat[idx].back[ND-1] ;
	    // t first
	    for( j = 0 ; j < NCNC ; j++ ) {
	      lat[idx].O[ND-1][j] = (GLU_real)uind[a] + I * (GLU_real)uind[a + 1] ;
	      a += 2 ;
	    }
	    for( j = 0 ; j < NCNC ; j++ ) {
	      lat[shift].O[ND-1][j] = (GLU_real)uind[a] + I * (GLU_real)uind[a + 1] ;
	      a += 2 ;
	    }
	    // then the others (xyz)?
	    for( mu = 0 ;  mu < ND - 1 ; mu++ ) {
	      for( j = 0 ; j < NCNC ; j++ ) {
		lat[idx].O[mu][j] = (GLU_real)uind[a] + I * (GLU_real)uind[a + 1] ;
		a += 2 ;
	      }
	      size_t shift = lat[idx].back[mu] ;
	      for( j = 0 ; j < NCNC ; j++ ) {
		lat[shift].O[mu][j] = (GLU_real)uind[a] + I * (GLU_real)uind[a + 1] ;
		a += 2 ;
	      }
	    }
	  }
	  
	}}}}
    
  *chksum = k ;
  free( uind ) ;
  return GLU_SUCCESS ; 
}

// writes my gauge field out in the CERN format 
void
write_CLS_field( const struct site *__restrict lat ,
		 FILE *__restrict outfile )
{
  // the size of the link matrices
  const size_t PM = 2 ;
  const size_t Complex = 2 ;
  // chunk we will read
  const size_t stride = PM * Complex * NCNC * ND ; 
  // modifier size
  const size_t Mod[ 4 ] = { 1 , 1 , 1 , 1 } ;
  
  uint32_t NAV[ ND ] ; 
  size_t mu ;
  NAV[ 0 ] = Mod[ ND-1 ] * Latt.dims[ ND - 1 ] ;
  for( mu = 0 ; mu < ND - 1 ; mu++ ) {
    NAV[ 1 + mu ] = Mod[mu] * Latt.dims[ mu ] ;
  }
  if( WORDS_BIGENDIAN ) {
    bswap_32( ND , NAV ) ;
  }
  
  // first of all set calculate the plaquette checksum in header
  double plaq[1] = { 0 } ;

  plaq[0] = NC * av_plaquette( lat ) ;
  if( WORDS_BIGENDIAN ) {
    bswap_64( 1 , plaq ) ;
  }
  
  // write the checksums
  fwrite( NAV , ( ND ) * sizeof( uint32_t ) , 1 , outfile ) ;
  fwrite( plaq , sizeof( double ) , 1 , outfile ) ;

  double *uoutd = calloc( stride , sizeof( double ) ) ; 

  size_t t , x , y , z ;
  for( t = 0 ; t < Mod[3]*Latt.dims[ ND-1 ] ; t++ ) {
    for( x = 0 ; x < Mod[0]*Latt.dims[ 0 ] ; x++ ) {
      for( y = 0 ; y < Mod[1]*Latt.dims[ 1 ] ; y++ ) {
	for( z = 0 ; z < Mod[2]*Latt.dims[ 2 ] ; z++ ) {

	  // convert it into local GLU coordinates
	  if( ( x + y + z + t )%2 ) {

	    // get the local, actual sites
	    int xloc[4] ;
	    xloc[0] = (int)(x%Latt.dims[0]) ;
	    xloc[1] = (int)(y%Latt.dims[1]) ;
	    xloc[2] = (int)(z%Latt.dims[2]) ;
	    xloc[3] = (int)(t%Latt.dims[3]) ;
	    
	    // config idx
	    size_t idx = gen_site( xloc ) ;
	    size_t shift = lat[idx].back[ ND-1 ] ;
	    
	    size_t j , a = 0 ;
	    // t first
	    for( j = 0 ; j < NCNC ; j++ ) {
	      uoutd[ a + 0 ] = ( double )creal( lat[idx].O[ ND - 1 ][ j ] ) ; 
	      uoutd[ a + 1 ] = ( double )cimag( lat[idx].O[ ND - 1 ][ j ] ) ;
	      a += 2 ;
	    }

	    for( j = 0 ; j < NCNC ; j++ ) {
	      uoutd[ a + 0 ] = ( double )creal( lat[shift].O[ ND - 1 ][ j ] ) ; 
	      uoutd[ a + 1 ] = ( double )cimag( lat[shift].O[ ND - 1 ][ j ] ) ;
	      a += 2 ;
	    }
	  
	    // then xyz
	    for( mu = 0 ;  mu < ND-1 ; mu++ ) {
	      for( j = 0 ; j < NCNC ; j++ ) {
		uoutd[ a + 0 ] = ( double )creal( lat[idx].O[mu][j] ) ; 
		uoutd[ a + 1 ] = ( double )cimag( lat[idx].O[mu][j] ) ; 
		a += 2 ;
	      }
	      shift = lat[idx].back[mu] ;
	      for( j = 0 ; j < NCNC ; j++ ) {
		uoutd[ a + 0 ] = ( double )creal( lat[shift].O[mu][j] ) ; 
		uoutd[ a + 1 ] = ( double )cimag( lat[shift].O[mu][j] ) ; 
		a += 2 ;
	      }
	    }

	    // write it out
	    if( WORDS_BIGENDIAN ) {
	      bswap_64( stride , uoutd ) ;
	    }
	    fwrite( uoutd , sizeof( double ) , stride , outfile ) ; 
	  }

	  // tzyx
	}}}}

  free( uoutd ) ;

  return ;
}
