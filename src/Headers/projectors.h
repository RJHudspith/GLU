/*
    Copyright 2013 Renwick James Hudspith

    This file (projectors.h) is part of GLU.

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
   @file projectors.h
   @brief function prototypes for smearing projections used, unrolled heavily
 */

#ifndef GLU_PROJECTORS_H
#define GLU_PROJECTORS_H

/**
   @fn void print_smearing_obs( const struct site *__restrict lat , const int type , const size_t count , const GLU_bool hypercubically_blocked )
   @brief prints to stdout the links and plaquettes
 */
void
print_smearing_obs( const struct site *__restrict lat , 
		    const int type ,
		    const size_t count ,
		    const GLU_bool hypercubically_blocked ) ;

/**
   @fn void project_APE( GLU_complex smeared_link[ NCNC ] , const GLU_complex staple[ NCNC ] , const GLU_complex link[ NCNC ] , const double smear_alpha , const double al )
   @brief d and sped up APE-projection
   @param ape :: the ape-projected link
   @param staple :: the computed staple
   @param link :: the link connecting the staples
   @param smear_alpha :: the smearing tuning parameter
   @param al :: ( 1 - alpha ) where alpha is "unnormalised"
   d and sped up APE-projection, uses the
   general reunitarisation rather than anything
   fancy..
 **/
void
project_APE( GLU_complex smeared_link[ NCNC ] , 
	     const GLU_complex staple[ NCNC ] , 
	     const GLU_complex link[ NCNC ] , 
	     const double smear_alpha , 	     
	     const double al ) ;

/**
    @fn void project_LOG( GLU_complex smeared_link[ NCNC ] , const GLU_complex staple[ NCNC ] , const GLU_complex link[ NCNC ] , const double smear_alpha )
   @brief LOG-smearing projection
   @param log :: the value of the loged link
   @param staple :: the computed staple
   @param link :: the link connecting the staples
   @param smear_alpha :: the normalised smearing alpha
   General LOG projection, the staple must already be in hermitian form
 **/
void
project_LOG( GLU_complex smeared_link[ NCNC ] , 
	     GLU_complex staple[ NCNC ] , 
	     const GLU_complex link[ NCNC ] , 
	     const double smear_alpha ,
	     const double al ) ;

/**
    @fn void project_LOG_short( GLU_complex smeared_link[ NCNC ] , const GLU_complex staple[ NCNC ] , const GLU_complex link[ NCNC ] , const double smear_alpha )
   @brief LOG-smearing projection
   @param log :: the value of the loged link
   @param staple :: the computed staple
   @param link :: the link connecting the staples
   @param smear_alpha :: the normalised smearing alpha
   General LOG projection, the staple must already be in hermitian form, uses eponentiate_short()
 **/
void
project_LOG_short( GLU_complex smeared_link[ NCNC ] , 
		   GLU_complex staple[ NCNC ] , 
		   const GLU_complex link[ NCNC ] , 
		   const double smear_alpha ,
		   const double al ) ;

/**
   @fn void project_LOG_wflow( GLU_complex smeared_link[ NCNC ] , GLU_complex *__restrict staple , const GLU_complex link[ NCNC ] , const double smear_alpha )
   @brief log smearing projection for the log-wilson-flow
   @param log :: output link from log smearing transform
   @param staple :: HERMSIZE staple
   @param link :: link being smeared
   @param smear_alpha :: the normalised smearing alpha
 **/
void
project_LOG_wflow( GLU_complex smeared_link[ NCNC ] , 
		   GLU_complex *__restrict staple , 
		   const GLU_complex link[ NCNC ] , 
		   const double smear_alpha ) ;

/**
    @fn void project_LOG_wflow_short( GLU_complex log[ NCNC ] , GLU_complex *__restrict staple , const GLU_complex link[ NCNC ] , const double smear_alpha )
   @brief log smearing projection for the log-wilson-flow
   @param log :: output link from log smearing transform
   @param staple :: #HERMSIZE staple
   @param link :: link being smeared
   @param smear_alpha :: the normalised smearing alpha
   uses the exponentiate_short() routine
 **/
void
project_LOG_wflow_short( GLU_complex log[ NCNC ] , 
			 GLU_complex *__restrict staple , 
			 const GLU_complex link[ NCNC ] , 
			 const double smear_alpha ) ;

/**
   @fn void project_STOUT( GLU_complex stout[ NCNC ] , const GLU_complex staple[ NCNC ] , const GLU_complex link[ NCNC ] , const double smear_alpha )
   @brief STOUT-smearing projection
   @param stout :: the value of the stouted link
   @param staple :: the computed staple
   @param link :: the link connecting the staples
   @param smear_alpha :: the normalised smearing alpha
   General STOUT projection ...
 **/
void
project_STOUT( GLU_complex stout[ NCNC ] , 
	       const GLU_complex staple[ NCNC ] , 
	       const GLU_complex link[ NCNC ] , 
	       const double smear_alpha ,
	       const double al ) ;

/**
   @fn void project_STOUT_short( GLU_complex smeared_link[ NCNC ] , const GLU_complex staple[ NCNC ] , const GLU_complex link[ NCNC ] , const double smear_alpha )
   @brief STOUT-smearing projection
   @param stout :: the value of the stouted link
   @param staple :: the computed staple
   @param link :: the link connecting the staples
   @param smear_alpha :: the normalised smearing alpha
   General STOUT projection uses the shortened versions of links, justifiably taking fewer operations and relying on conjugacy.
 **/
void
project_STOUT_short( GLU_complex smeared_link[ NCNC ] , 
		     const GLU_complex staple[ NCNC ] , 
		     const GLU_complex link[ NCNC ] , 
		     const double smear_alpha ,
		     const double al ) ;

/**
   @fn void project_STOUT_wflow( GLU_complex smeared_link[ NCNC ] , const GLU_complex *__restrict staple , const GLU_complex link[ NCNC ] , const double smear_alpha )
   @brief Had to make a shorter version for large matrices
   @param stout :: the value of the stouted link
   @param staple :: the computed staple shortened to HERMSIZE
   @param link :: the link connecting the staples
   @param smear_alpha :: the normalised smearing alpha
 **/
void
project_STOUT_wflow( GLU_complex smeared_link[ NCNC ] , 
		     const GLU_complex *__restrict staple , 
		     const GLU_complex link[ NCNC ] , 
		     const double smear_alpha ) ;

/**
   @fn void project_STOUT_wflow_short( GLU_complex stout[ NCNC ] , GLU_complex *__restrict staple , const GLU_complex link[ NCNC ] , const double smear_alpha )
   @brief Had to make a shorter version for large matrices
   @param stout :: the value of the stouted link
   @param staple :: the computed staple shortened to HERMSIZE
   @param link :: the link connecting the staples
   @param smear_alpha :: the normalised smearing alpha
   Same as the project_STOUT_wflow() function, but uses the exponentiate_short() function from exactQ.h
 **/
void
project_STOUT_wflow_short( GLU_complex stout[ NCNC ] , 
			   GLU_complex *__restrict staple , 
			   const GLU_complex link[ NCNC ] , 
			   const double smear_alpha ) ;

#endif
