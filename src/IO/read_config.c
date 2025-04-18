/*
Copyright 2013-2025 Renwick James Hudspith

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
 */
#include "Mainfile.h"

#include "CERN.h"          // openQCD config reader/writer
#include "BPST_config.h"   // generate a BPST config
#include "GLU_memcheck.h"  // do we have enough memory for the fast reader?
#include "HIREP.h"         // read a HIREP file
#include "plaqs_links.h"   // average plaquette, link traces
#include "random_config.h" // random SU(NC) configuration
#include "readers.h"       // read config files

// comparison between checksums
static int 
check_sums( const double plaq , 
	    const double tr ,
	    const uint32_t chksum ,
	    const struct head_data HEAD_DATA )
{
  double TTOL = 0.0 ;
  int error = 0 ;
  enum{ DISASTER = 3 } ;
  
  switch( Latt.head ) {
  case MILC_HEADER :
  case ILDG_SCIDAC_HEADER :
  case SCIDAC_HEADER :
    // do a check only on the sum29
    if( chksum != HEAD_DATA.checksum )  {
      fprintf( stderr , "\n[IO] Unequal checksums [Calc] %x [Read] %x \n\n" , 
	       chksum , HEAD_DATA.checksum ) ; 
      return GLU_FAILURE ;
    }
    break ;
  case ILDG_BQCD_HEADER :
    if( HEAD_DATA.precision == FLOAT_PREC ) {
      TTOL = 1E-6 ;
    } else {
      TTOL = PLAQ_AND_TRACE_TOL ;
    }
    if( fabs( HEAD_DATA.plaquette - plaq ) > TTOL ) {
      fprintf( stderr , "\n[IO] Unequal Plaquettes %e %e \n\n" , 
	       plaq , HEAD_DATA.plaquette ) ; 
      return GLU_FAILURE ;
    }
    // and the CRC of the binary data
    if( chksum != HEAD_DATA.checksum )  {
      fprintf( stderr , "\n[IO] Unequal checksums Calc %x || Read %x \n\n" , 
	       chksum , HEAD_DATA.checksum ) ; 
      return GLU_FAILURE ;
    }
    break ;
  case HIREP_HEADER :
  case CERN_HEADER :
    // only check available is the plaquette
    if( fabs( plaq - HEAD_DATA.plaquette ) > 10*PREC_TOL ) {
      fprintf( stderr , "[IO] Config Plaquette Mismatch %e vs %e "
	       " < diff > %e ... Leaving \n" , plaq , HEAD_DATA.plaquette , 
	       fabs( plaq - HEAD_DATA.plaquette ) ) ;
      return GLU_FAILURE ;
    }
    break ;
  case NERSC_HEADER :
    // it is disastrous if all three checks fail ...
    if( chksum != HEAD_DATA.checksum )  {
      fprintf( stderr , "\n[IO] Unequal checksums Calc %x || Read %x \n\n" , 
	       chksum , HEAD_DATA.checksum ) ; 
      error ++ ; 
      // TOL is defined as 10^-6
    } if( fabs( plaq - HEAD_DATA.plaquette ) > PLAQ_AND_TRACE_TOL ) {
      fprintf( stderr , "\n[IO] Unequal Plaquettes Calc %f || Read %f " 
	       " < diff > %e \n\n" , plaq , HEAD_DATA.plaquette , 
	       fabs( plaq - HEAD_DATA.plaquette ) ) ; 
      error ++ ; 
    } if( fabs( tr - HEAD_DATA.trace) > PLAQ_AND_TRACE_TOL ) {
      fprintf( stderr , "\n[IO] Unequal Link_Traces Calc %1.8f || "
	       "Read %1.8f \n\n" , tr , HEAD_DATA.trace ) ; 
      error ++ ; 
    }
    // pretty printing
    fprintf( stdout , "[IO] Header     Trace :: %f           || Plaq :: %f \n" 
	     , HEAD_DATA.trace , HEAD_DATA.plaquette ) ; 
    // if everything is wrong we leave
    if( error == DISASTER ) {
      fprintf( stdout , "[IO] NONE of the NERSC Checksums match, "
	       "this is a problem .. Leaving \n") ; 
      return GLU_FAILURE ;
    }
    break ;
  case LIME_HEADER :
  case RANDOM_CONFIG :
  case UNIT_GAUGE :
  case INSTANTON :
  case JLQCD_HEADER :
    break ;
  case UNSUPPORTED :
    return GLU_FAILURE ;
  }
  fprintf( stdout , "[IO] Calculated Trace :: %1.15f  || Plaq :: %1.15f \n" , 
	   tr , plaq ) ; 
  return GLU_SUCCESS ; // may only be partially successful but I am optimistic
}

// things that check the checksums //
int
checks( struct site *__restrict lat , 
	uint32_t chksum ,
	struct head_data HEAD_DATA )
{
  const double plaq = av_plaquette( lat ) ; 
  const double tr = links( lat ) ; 
  return check_sums( plaq , tr , chksum , HEAD_DATA ) ;
}

// we wrap this one to our new reader ... 
int
get_config_SUNC( FILE *__restrict CONFIG , 
		 struct site *__restrict lat ,
		 const struct head_data HEAD_DATA ,
		 const char *config_in )
{
  uint32_t chksum ;
  switch( Latt.head ) {
  case LIME_HEADER : // is the same but doesn't care about the checksums
  case ILDG_BQCD_HEADER : // basically all the same NERSC NCxNC
  case ILDG_SCIDAC_HEADER : // ILDG and SCIDAC
  case SCIDAC_HEADER : // Scidac's binary data is compliant
  case MILC_HEADER :   // MILC's binary data is the same
  case NERSC_HEADER :
  case JLQCD_HEADER :
    return lattice_reader_suNC_posix( lat , config_in , CONFIG , HEAD_DATA ) ;
  case HIREP_HEADER : // HIREP uses a weird geometry compared to everyone else
    read_gauge_field( lat , CONFIG , &chksum ) ;
    return chksum ;
  case CERN_HEADER :
    read_CLS_field( lat , config_in , &chksum ) ;
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
    fprintf( stderr , "[IO] Unrecognised HEADER type .. Leaving \n" ) ;
    return GLU_FAILURE ;
  }
  return GLU_FAILURE ;
}
