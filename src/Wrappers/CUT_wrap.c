/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (CUT_wrap.c) is part of GLU.

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
   @file CUT_wrap.c
   @brief this wraps the momentum cutting routines and what have you.
 */
#include "Mainfile.h"
#include "3Dcuts.h" // instantaneous spatial and temporal props
#include "config_gluons.h" // configuration space gluon propagators
#include "cuts.h" // alpha_s measurements
#include "GLU_timer.h" // for the timer
#include "POLY.h" // for the static potential
#include "Qsusc.h" // for the topological susceptibility correlator
#include "smearing_param.h" // smeared gluon propagator

// wrapper for the variousfa momentum-cutting and other measurement options
void
cuts_wrap_struct( struct site *__restrict lat , 
		  const struct cut_info CUTINFO ,
		  const struct sm_info SMINFO )
{
  start_timer( ) ;

  ///////////// CREATE CUTS ////////////
  switch( CUTINFO.dir ) {
  case INSTANTANEOUS_GLUONS :
    cuts_spatial( lat , CUTINFO ) ;
    break ;
  case CONFIGSPACE_GLUONS :
    cuts_struct_configspace( lat , CUTINFO , SMINFO ) ;
    break ;
  case SMEARED_GLUONS :
    cuts_struct_smeared( lat , CUTINFO , SMINFO ) ;
    break ;
  case STATIC_POTENTIAL :
    Coul_staticpot( lat , CUTINFO , SMINFO ) ;
    break ;
  case TOPOLOGICAL_CORRELATOR :
  case TOPOLOGICAL_SUSCEPTIBILITY :
    compute_Qsusc_step( lat , CUTINFO , SMINFO ) ;
    break ;
  case GLUON_PROPS :
  case EXCEPTIONAL :
  case NONEXCEPTIONAL :
  case FIELDS :
    // all of the other momentum space routines are in here
    cuts_struct( lat , CUTINFO ) ;
    break ;
  default :
    return ;
  }

  print_time( ) ;

  return ;
}
