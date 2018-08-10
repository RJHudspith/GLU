/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (GLU_enum.h) is part of GLU.

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
   @file GLU_enum.h
   @brief I have decided as well to include all of the enumerated types here as well
 **/
#ifndef GLU_ENUM_H
#define GLU_ENUM_H

/**
   @enum cline_arg
   @brief command line argument lookup
 */
typedef enum{ 
  HELP = 1 , 
  INPUT_FILE = 2 , 
  READ = 4 , 
  WRITE = 6 } cline_arg ;

/**
   @enum file_prec
   @brief defs for the readers
*/
typedef enum 
  { FLOAT_PREC , 
    DOUBLE_PREC } file_prec ;

/**
   @enum GLU_speed
   @brief defs for the speed 
*/
typedef enum
  { FAST ,
    MODERATE ,
    SLOW } GLU_speed ;

/**
   @enum GLU_smeardir
   @brief directions for the smearing
   requires definition of ND
 */
typedef enum
  { SPATIAL_LINKS_ONLY = ND - 1 ,
    ALL_DIRECTIONS = ND } GLU_smeardir ;

/**
   @enum GLU_endian
   @brief enumeration of the endianess
*/
typedef enum 
  { L_ENDIAN , 
    B_ENDIAN } GLU_endian ;  

/**
   @enum GLU_direction
   @brief wilson flow direction
 */
enum GLU_direction 
  { GLU_BACKWARD = -1 ,
    GLU_FORWARD = 1 } ;

/**
   @enum GLU_bool
   @brief standard boolean type deal
   In C FALSE = 0 , TRUE = {ANYTHING_ELSE}
 */
typedef enum 
  { GLU_FALSE , 
    GLU_TRUE } GLU_bool ;

/**
   @enum header_mode
   @brief enums for the available headers my code can read
 */
typedef enum 
  { UNSUPPORTED ,
    NERSC_HEADER ,
    HIREP_HEADER ,
    MILC_HEADER ,
    SCIDAC_HEADER ,
    ILDG_SCIDAC_HEADER ,
    ILDG_BQCD_HEADER ,
    LIME_HEADER ,
    CERN_HEADER ,
    RANDOM_CONFIG ,
    UNIT_GAUGE ,
    INSTANTON } header_mode ;

/**
   @enum GLU_fixing
   @brief enums for the gauge fixing type
*/
typedef enum  
  { DEFAULT_NOFIX ,
    GLU_AXIALT_FIX ,
    GLU_COULOMB_FIX ,
    GLU_COULOMB_RESIDUAL_FIX ,
    GLU_LANDAU_FIX ,
    GLU_MAG_FIX } GLU_fixing ;

/**
   @enum GF_improvements
   @brief gauge fixing improvements
 */
typedef enum 
  { NO_IMPROVE , 
    MAG_IMPROVE , 
    SMPREC_IMPROVE ,
    RESIDUAL_IMPROVE } GF_improvements ;

/**
   @enum GLU_mode
   @brief enums for the mode
*/
typedef enum
  { MODE_REWRITE , 
    MODE_GF ,
    MODE_CUTS ,
    MODE_SMEARING ,
    MODE_CROSS_U1 ,
    MODE_HEATBATH } GLU_mode ;

/**
   @enum smearing_types
   @brief available smearing types 
   prepended with SM to avoid clashes
 */
typedef enum 
  { SM_NOSMEARING , 
    SM_APE , 
    SM_STOUT , 
    SM_LOG ,
    SM_HYP ,
    SM_HEX ,
    SM_HYL ,
    SM_WFLOW_STOUT ,
    SM_WFLOW_LOG ,
    SM_ADAPTWFLOW_STOUT , 
    SM_ADAPTWFLOW_LOG } smearing_types ;

/**
   @enum cut_mode
   @brief available cutting modes
 */
typedef enum 
  { EXCEPTIONAL ,
    NONEXCEPTIONAL ,
    FIELDS ,
    SMEARED_GLUONS ,
    INSTANTANEOUS_GLUONS ,
    CONFIGSPACE_GLUONS ,
    GLUON_PROPS ,
    STATIC_POTENTIAL ,
    TOPOLOGICAL_CORRELATOR ,
    TOPOLOGICAL_MOMENTS ,
    TOPOLOGICAL_SUSCEPTIBILITY } cut_mode ;

/**
   @enum lie_field_def
   @brief We offer the approximate (LINEAR) or the exact (LOGARITHMIC) gluon fields
 */
typedef enum 
  { LINEAR_DEF ,
    LOG_DEF } lie_field_def ;

/**
   @enum momentum_cut_def
   @brief enumerate the cutting types
 */
typedef enum
  { HYPERCUBIC_CUT , 
    PSQ_CUT ,
    CYLINDER_CUT ,
    CYLINDER_AND_CONICAL_CUT } momentum_cut_def ;

/**
   @enum GLU_output
   @brief enumerate the output types
 */
typedef enum
  { NO_OUTPUT ,
    OUTPUT_SMALL ,  // logically the fewest number of params saved 
    OUTPUT_GAUGE ,  // NERSC's gauge , the top NC-1 rows and then complete
    OUTPUT_NCxNC ,  // the whole matrix, wasteful
    OUTPUT_MILC ,   // the whole matrix, wasteful
    OUTPUT_ILDG ,   // the whole matrix, wasteful
    OUTPUT_SCIDAC , // the whole matrix, wasteful
    OUTPUT_HIREP ,  // the whole matrix in an odd geometry
    OUTPUT_CERN ,   // the whole matrix in an odd +/- geometry
  } GLU_output ; // the whole matrix in HIREP's order  , wasteful

/**
   @enum config_size
   @brief for writing out files
   loops number of spaces for this format uses NC and NCNC
 */
enum config_size
  { LOOP_SMALL = NCNC - 1 ,
    LOOP_GAUGE = 2 * ( NC - 1 ) * NC ,
    LOOP_NCxNC = 2 * NCNC } ;

/**
   @enum wflow_type
   @brief for wilson flow information
 */
typedef enum 
  { EULER ,
    RK3_FAST ,
    RK3_SLOW , 
    RK3_ADAPTIVE } wflow_type ;

/**
   @enum U1_meas
   @brief type of U1 measurement requested
 */
typedef enum
  { U1_PLAQUETTE ,
    U1_RECTANGLE ,
    U1_TOPOLOGICAL } U1_meas ;

/**
   @enum MOMTYPE
   @brief momentum definition being used
 */
typedef enum
  { MODES ,
    TWOSIN } MOMTYPE ;

#endif
