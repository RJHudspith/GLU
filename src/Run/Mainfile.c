/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (Mainfile.c) is part of GLU.

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
   @file Mainfile.c
   @brief the driver for the rest of the code
 */

// Usual dependencies
#include "Mainfile.h"

#include "GLUlib_wrap.h"    // wrappers for general code functionality
#include "input_help.h"     // help about the input file
#include "input_reader.h"   // reading the input file

// Lattice geometry is an important global, used everywhere
struct latt_info Latt ;

// main file
int 
main( const int argc ,
      const char *argv[] )
{
  // for the input file help
  if( argc == 2 ) {
    return GLU_helps_those_who_help_themselves( argv[1] ) ;
  }

  // bit of control here 
  if( argc < READ || !(argc%2) ) { return GLUsage() ; }

  // set this
  Latt.argc = argc ;

  // initialise GLU, its so sticky
  attach_GLU( ) ;

  // parse the input file into a nice big struct
  struct infile_data infile_struct ;
  if( get_input_data( &infile_struct , argv[ INPUT_FILE ] ) == GLU_FAILURE ) {
    return GLU_FAILURE ;
  }

  // and make the input data constant
  const struct infile_data INFILE = infile_struct ;

  // and the options
  switch( INFILE.mode ) {
  case MODE_GF :
    read_and_fix( argv[ READ ] , INFILE.rtrans , INFILE.GFINFO , INFILE.SMINFO ,
		  argv[ WRITE ] , INFILE.storage , INFILE.output_details ) ;
    break ;
  case MODE_CUTS :
    read_and_cut( argv[ READ ] , INFILE.CUTINFO , INFILE.SMINFO ) ;
    break ;
  case MODE_SMEARING :
    read_and_smear( argv[ READ ] , INFILE.rtrans , INFILE.SMINFO ,
		    argv[ WRITE ] , INFILE.storage , INFILE.output_details ) ;
    break ;
  case MODE_CROSS_U1 :
    read_and_U1( argv[ READ ] , INFILE.rtrans , INFILE.U1INFO ,
		 argv[ WRITE ] , INFILE.storage , INFILE.output_details ) ;
    break ;
  case MODE_HEATBATH :
    heatbath( argv[ READ ] , INFILE.HBINFO , INFILE.head , argv[ WRITE ] , 
	      INFILE.storage , INFILE.output_details ) ;
    break ;
  default :
    read_and_check( argv[ READ ] , INFILE.rtrans , argv[ WRITE ] , 
		    INFILE.storage , INFILE.output_details ) ;
    break ;
  }

  // un-initialise GLU
  unstick_GLU( ) ;

  return GLU_SUCCESS ;
}
