/*
    Copyright 2013 Renwick James Hudspith

    This file (read_config.c) is part of GLU.

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
   @file read_config.c 
   @brief gets the information about our configs from the header

   @warning only HiRep and NERSC and MILC? supported atm.
 */

#include "Mainfile.h"

#include "BPST_config.h"   // generate a BPST config
#include "GLU_memcheck.h"  // do we have enough memory for the fast reader?
#include "HIREP.h"         // read a HIREP file
#include "plaqs_links.h"   // average plaquette, link traces
#include "random_config.h" // random SU(NC) configuration
#include "readers.h"       // read config files

// comparison between checksums
static int 
check_sums( plaq , tr , chksum , HEAD_DATA )
     const double plaq , tr ;
     const uint32_t chksum ; 
     const struct head_data HEAD_DATA ;
{
  // these headers I use the same check, the %29 one
  if( Latt.head == MILC_HEADER || Latt.head == ILDG_SCIDAC_HEADER || 
      Latt.head == SCIDAC_HEADER ) {
    // do a check only on the sum29
    if( chksum != HEAD_DATA.checksum )  {
      printf( "\n[IO] Unequal checksums [Calc] %x [Read] %x \n\n" , 
	      chksum , HEAD_DATA.checksum ) ; 
      return GLU_FAILURE ;
    }
  } else if( Latt.head == ILDG_BQCD_HEADER ) {
    // BQCD provides a CRC and a plaquette
    double TTOL = 0.0 ;
    if( HEAD_DATA.precision == FLOAT_PREC ) {
      TTOL = 1E-6 ;
    } else {
      TTOL = 1E-14 ;
    }
    if( fabs( HEAD_DATA.plaquette - plaq ) > TTOL ) {
      printf( "\n[IO] Unequal Plaquettes %e %e \n\n" , 
	      plaq , HEAD_DATA.plaquette ) ; 
      return GLU_FAILURE ;
    }
    // and the CRC of the binary data
    if( chksum != HEAD_DATA.checksum )  {
      printf( "\n[IO] Unequal checksums Calc %x || Read %x \n\n" , 
	      chksum , HEAD_DATA.checksum ) ; 
      return GLU_FAILURE ;
    }
  } else if( Latt.head == HIREP_HEADER ) {
    // only check available is the plaquette
    if( fabs( plaq - HEAD_DATA.plaquette ) > PREC_TOL ) {
      printf("[IO] HIREP header Plaquette Mismatch %e vs %e ... Leaving \n" ,
	     plaq , HEAD_DATA.plaquette ) ;
      return GLU_FAILURE ;
    }
  } else if( Latt.head == NERSC_HEADER ) {
    enum{ DISASTER = 3 } ;
    // it is disastrous if all three checks fail ...
    int error = 0 ; 
    if( chksum != HEAD_DATA.checksum )  {
      printf( "\n[IO] Unequal checksums Calc %x || Read %x \n\n" , 
	      chksum , HEAD_DATA.checksum ) ; 
      error ++ ; 
      // TOL is defined as 10^-6
    }  if( fabs( plaq - HEAD_DATA.plaquette ) > PLAQ_AND_TRACE_TOL ) {
      printf( "\n[IO] Unequal Plaquettes Calc %f || Read %f \n\n" , 
	      plaq , HEAD_DATA.plaquette ) ; 
      error ++ ; 
    } if( fabs( tr - HEAD_DATA.trace) > PLAQ_AND_TRACE_TOL ) {
      printf( "\n[IO] Unequal Link_Traces Calc %1.8f || Read %1.8f \n\n" , 
	      tr , HEAD_DATA.trace ) ; 
      error ++ ; 
    }
    // pretty printing
    printf( "[IO] Header     Trace :: %f           || Plaq :: %f \n" , 
	    HEAD_DATA.trace , HEAD_DATA.plaquette ) ; 
    // if everything is wrong we leave
    if( unlikely( error == DISASTER ) ) {
      printf("[IO] NONE of the NERSC Checksums match, this is a problem .. Leaving \n") ; 
      return GLU_FAILURE ;
    }
  } 
 printf( "[IO] Calculated Trace :: %1.15f  || Plaq :: %1.15f \n" , 
	 tr , plaq ) ; 
  return GLU_SUCCESS ; // may only be partially successful but I am optimistic
}

// things that check the checksums //
short int
checks( struct site *__restrict lat , 
	uint32_t chksum ,
	struct head_data HEAD_DATA )
{
  const double plaq = av_plaquette( lat ) ; 
  const double tr = links( lat ) ; 
  return check_sums( plaq , tr , chksum , HEAD_DATA ) ;
}

// we wrap this one to our new reader ... 
uint32_t
get_config_SUNC( FILE *__restrict CONFIG , 
		 struct site *__restrict lat ,
		 const struct head_data HEAD_DATA )
{
  uint32_t chksum ;
  switch( Latt.head ) {
  case LIME_HEADER : // is the same but doesn't care about the checksums
  case ILDG_BQCD_HEADER : // basically all the same NERSC NCxNC
  case ILDG_SCIDAC_HEADER : // ILDG and SCIDAC
  case SCIDAC_HEADER : // Scidac's binary data is compliant
  case MILC_HEADER : // MILC's binary data is the same
  case NERSC_HEADER :
    if( HEAD_DATA.config_type == OUTPUT_SMALL ||
	HEAD_DATA.config_type == OUTPUT_GAUGE ||
	HEAD_DATA.config_type == OUTPUT_NCxNC ) {
      const int SAFETY = have_memory_readers_writers( HEAD_DATA.config_type ) ;
      if( SAFETY != FAST ) {
	return lattice_reader_suNC_cheaper( lat , CONFIG , HEAD_DATA ) ;
      } else {
	return lattice_reader_suNC( lat , CONFIG , HEAD_DATA ) ;
      }
    } break ;
  case HIREP_HEADER : // HIREP uses a weird geometry compared to everyone else
    read_gauge_field( lat , CONFIG , &chksum ) ;
    return chksum ;
  case RANDOM_CONFIG :
    random_suNC( lat ) ;
    return GLU_SUCCESS ;
  case UNIT_GAUGE :
    trivial( lat ) ;
    return GLU_SUCCESS ;
  case INSTANTON :
    instanton_config( lat ) ;
    return GLU_SUCCESS ;
  default : 
    printf( "[IO] Unrecognised HEADER type .. Leaving \n" ) ;
    return GLU_FAILURE ;
  }
  return GLU_FAILURE ;
}
