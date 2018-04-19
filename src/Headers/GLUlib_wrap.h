/*
    Copyright 2013-2018 Renwick James Hudspith

    This file (GLUlib_wrap.h) is part of GLU.

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
#ifdef __cplusplus
extern "C"{
#endif

/**
   @file GLUlib_wrap.h
   @brief prototype wrapper functions for the general library workings
 */
#ifndef GLULIB_WRAP_H
#define GLULIB_WRAP_H

/**
   @fn void attach_GLU( void )
   @brief initialises and precomputes some stuff
 */
void
attach_GLU( void ) ;

/**
   @fn struct site* read_file( struct head_data *HEAD_DATA , const char *config_in )
   @brief read a configuration file that GLU understands and allocates lat
   @param HEAD_DATA :: information from the config header (passed by reference)
   @param config_in :: string of the configuration file's location
   @warning allocates the gauge field
   @return the gauge field, or if it fails: NULL
 */
struct site*
read_file( struct head_data *HEAD_DATA , 
	   const char *config_in ) ;

/**
   @fn int write_configuration( struct site *lat , const char *outfile , const GLU_output storage , const char *output_details )
   @brief write out a gauge configuration
   @param lat :: lattice gauge fields
   @param outfile :: output file configuration
   @param storage :: storage type (NO_STORAGE means we don't write it out)
   @param output_details :: often a header has a string describing it, this is that
   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int
write_configuration( struct site *lat , 
		     const char *outfile , 
		     const GLU_output storage , 
		     const char *output_details ) ;

/**
   @fn int heatbath( const char *infile , const struct hb_info HBINFO , const char *outfile , const GLU_output storage , const char *output_details )
   @brief read and check unitarity and gauge invariant quantities
   @param infile :: input configuration name
   @param outfile :: output file configuration
   @param storage :: storage type (NO_STORAGE means we don't write it out)
   @param output_details :: often a header has a string describing it, this is that
   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int
heatbath( const char *infile ,
	  const struct hb_info HBINFO ,
	  const char *outfile , 
	  const GLU_output storage , 
	  const char *output_details ) ;

/**
   @fn int read_and_check( const char *infile , const GLU_bool rtrans , const char *outfile , const GLU_output storage , const char *output_details )
   @brief read and check unitarity and gauge invariant quantities
   @param infile :: input configuration name
   @param rtrans :: do we want to randomly transform the initial configuration
   @param outfile :: output file configuration
   @param storage :: storage type (NO_STORAGE means we don't write it out)
   @param output_details :: often a header has a string describing it, this is that
   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int
read_and_check( const char *infile ,
		const GLU_bool rtrans ,
		const char *outfile , 
		const GLU_output storage , 
		const char *output_details ) ;

/**
   @fn int read_and_cut( const char *infile , const struct cut_info CUTINFO , const struct sm_info SMINFO )
   @brief read in a configuration file and perform cut-measurements
   @param infile :: input configuration name
   @param CUTINFO :: cutting information of the form #cut_info
   @param SMINFO :: smearing information of the form #sm_info
   @param output_details :: often a header has a string describing it, this is that
   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int
read_and_cut( const char *infile , 
	      const struct cut_info CUTINFO , 
	      const struct sm_info SMINFO ) ;

/**
   @fn int read_and_fix( const char *infile , const GLU_bool rtrans , const struct gf_info GFINFO , const struct sm_info SMINFO , const char *outfile , const GLU_output storage , const char *output_details ) 
   @brief reads a configuration file, gauge fixes and perhaps writes out a configuration too
   @param infile :: input configuration name
   @param rtrans :: do we want to randomly transform the initial configuration
   @param GFINGO :: gauge information of the form #gf_info
   @param SMINFO :: smearing transformations can be used in smeared gauge fixing
   @param outfile :: output file configuration
   @param storage :: storage type (NO_STORAGE means we don't write it out)
   @param output_details :: often a header has a string describing it, this is that
   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int
read_and_fix( const char *infile , 
	      const GLU_bool rtrans , 
	      const struct gf_info GFINFO , 
	      const struct sm_info SMINFO ,
	      const char *outfile , 
	      const GLU_output storage , 
	      const char *output_details ) ;

/**
   @param infile :: input configuration name
   @param outfile :: output file configuration
   @param storage :: storage type (NO_STORAGE means we don't write it out)
   @param output_details :: often a header has a string describing it, this is that
   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int
read_and_smear( const char *infile , 
		const GLU_bool rtrans , 
		const struct sm_info SMINFO ,
		const char *outfile , 
		const GLU_output storage , 
		const char *output_details ) ;

/**
   @fn int read_and_U1( const char *infile , const GLU_bool rtrans , const struct u1_info U1INFO , const char *outfile , const GLU_output storage , const char *output_details )
   @brief reads in a configuration file, creates a quenched SU(N)xU(1) configuration and perhaps writes it out
   @param infile :: input configuration name
   @param rtrans :: perform a random transformation on the initial configuration
   @param U1INFO :: U1 information of the form #u1_info
   @param outfile :: output file configuration
   @param storage :: storage type (NO_STORAGE means we don't write it out)
   @param output_details :: often a header has a string describing it, this is that
   @return #GLU_SUCCESS or #GLU_FAILURE
 */
int
read_and_U1( const char *infile , 
	     const GLU_bool rtrans ,
	     const struct u1_info U1INFO ,
	     const char *outfile , 
	     const GLU_output storage , 
	     const char *output_details ) ;

/**
   @fn void unstick_GLU( void )
   @brief cleans up and frees some allocated precomputations
 */
void
unstick_GLU( void ) ;

#endif

// c++ guard
#ifdef __cplusplus
}
#endif

