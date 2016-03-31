/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (cut_output.h) is part of GLU.

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
   @file cut_output.h
   @brief Function defs for the various output types used
   @ingroup Cuts
 */

#ifndef GLU_CUT_OUTPUT_H
#define GLU_CUT_OUTPUT_H

/**
   @fn int check_psq( const struct cut_info CUTINFO )
   @brief checks whether the psq is allowed

   @param CUTINFO :: cutting information
   Just a simple check for the cutting type
 **/
int
check_psq( const struct cut_info CUTINFO ) ;

/**
   @fn char* output_str_struct( const struct cut_info CUTINFO )
   @brief From the cutting information, creates the output string necessary
   @param CUTINFO :: general cutting information.
 **/
char*
output_str_struct( const struct cut_info CUTINFO ) ;

/**
   @fn void write_complex_g2g3_to_list( FILE *__restrict Ap , double complex *__restrict g2 , double complex *__restrict g3 , int num_mom[ 1 ] ) ;

   @brief This function writes out the real part of the trace of the gluonic two and three point functions.
   
   @param Ap :: File that we output to.
   @param g2 :: The gluon propagator
   @param g3 :: The gluon three point function
   @param num_mom :: number of momenta in the list after cutting.

   This function writes first the two point function
   and then the three point function corresponding to the
   momenta written at the top of the file.
 **/
void
write_complex_g2g3_to_list( FILE *__restrict Ap , 
			    double complex *__restrict g2 , 
			    double complex *__restrict g3 , 
			    int num_mom[ 1 ] ) ;

/**
   @fn void write_g2_to_list( FILE *__restrict Ap , double *__restrict g2 , int num_mom[ 1 ] ) ;

   @brief This function writes out the real part of the trace of the gluonic two and three point functions.
   
   @param Ap :: File that we output to.
   @param alpha :: The renormalised strong coupling
   @param num_mom :: number of momenta in the list after cutting.

   This function writes first the strong coupling to a file
 **/
void
write_g2_to_list( FILE *__restrict Ap , 
		  double *__restrict g2 , 
		  int num_mom[ 1 ] ) ;

/**
   @fn void write_g2g3_to_list( FILE *__restrict Ap , double *__restrict g2 , double *__restrict g3 , int num_mom[ 1 ] ) ;

   @brief This function writes out the real part of the trace of the gluonic two and three point functions.
   
   @param Ap :: File that we output to.
   @param g2 :: The gluon propagator
   @param g3 :: The gluon three point function
   @param num_mom :: number of momenta in the list after cutting.

   This function writes first the two point function
   and then the three point function corresponding to the
   momenta written at the top of the file.
 **/
void
write_g2g3_to_list( FILE *__restrict Ap , 
		    double *__restrict g2 , 
		    double *__restrict g3 , 
		    int num_mom[ 1 ] ) ;

/**
   @fn void write_lattice_fields( FILE *__restrict Ap , const struct site *__restrict A , const struct veclist *__restrict list , const int num_mom[1] )
   @brief This function writes out the full matrices.  
   
   @param Ap :: File that we output to.
   @param A :: gauge fields in momentum space
   @param list :: momenta in the \f$-\pi \rightarrow \pi\f$ BZ, kept after cutting.
   @param num_mom :: number of momenta in the list after cutting.
 **/
void
write_lattice_fields( FILE *__restrict Ap , 
		      const struct site *__restrict A , 
		      const struct veclist *__restrict list , 
		      int num_mom[1] ) ;

/**
   @fn void write_mom_veclist( FILE *__restrict Ap , int *__restrict num_mom , const struct veclist *__restrict list , const int DIR )
   @brief This writes the momentum list. DIR is the number of momenta to
   write usually ND or ND-1.
   @param Ap :: File that we output to.
   @param num_mom :: number of momenta in the list after cutting.
   @param list :: momenta stored in some definition of the BZ
   @param DIR :: number of directions of momenta used
 **/
void
write_mom_veclist( FILE *__restrict Ap , 
		   int *__restrict num_mom , 
		   const struct veclist *__restrict list ,
		   const int DIR ) ;

/**
   @fn void write_triplet_mom_list( FILE *__restrict Ap , size_t *__restrict num_mom , int *__restrict *__restrict momentum , int *__restrict *__restrict triplet , const size_t DIR )
   @brief write out the momentum list for the triplets

   @warning I only write out one of the three momenta because after projection the vertex is of p^2  
 **/
void
write_triplet_mom_list( FILE *__restrict Ap , 
			int *__restrict num_mom , 
			int *__restrict *__restrict momentum ,
			int *__restrict *__restrict triplet ,
			const size_t DIR ) ;

/**
   @fn void write_rr_values( FILE *__restrict Ap , int size[1] , const int *__restrict rsq , const size_t max_r2 , const size_t ARR_SIZE )
   @brief writer for the topological susceptibility and static potential
   @param Ap :: file we are writing out to
   @param size :: length of the array
   @param rsq :: rsq array
   @param max_r2 :: maximum r^2 we are writin out
   @param ARR_SIZE :: array length of rsq
 */
void
write_rr_values( FILE *__restrict Ap ,
		 int size[1] ,
		 const int *__restrict rsq ,
		 const size_t max_r2 ,
		 const size_t ARR_SIZE ) ;

/**
   @fn void write_tslice_list( FILE *__restrict Ap , int *__restrict LT )
   @brief write out the timeslice index for the configuration space gluon correlator measurement
   @param Ap :: the file we are writing out to
   @param LT :: the length of the time direction (1 element array)
 **/
void
write_tslice_list( FILE *__restrict Ap , 
		   int *__restrict LT ) ;

#endif
