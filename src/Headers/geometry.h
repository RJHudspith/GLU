/*
    Copyright 2013 Renwick James Hudspith

    This file (geometry.h) is part of GLU.

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
   @file geometry.h
   @brief function definitions concerning the lattice geometry
 */

#ifndef GLU_GEOMETRY_H
#define GLU_GEOMETRY_H

/**
   @fn int gen_site( const int x[ ND ] )
   @brief Generic return of the lattice index from the ND-vectors x. 
   @param x :: Converts the ND vector x == ( x , y , z , ... , t ) to our lexicographical index.
   Index vector to periodic lexicographical index value converter. The temporal direction is always last in this ordering.
 **/
int
gen_site( const int x[ ND ] ) ;

/**
   @fn int get_site_2piBZ( int x[ ND ] , const int DIMS )
   @brief Returns our location when given momenta in the \f$ 0 \rightarrow 2\pi \f$ BZ.
   @param x :: Momentum vector
   @param DIMS :: maximum dimension of the problem
 **/
int
get_site_2piBZ( int x[ ND ] , 
		const int DIMS ) ;

/**
   @fn int get_site_pipiBZ( int x[ ND ] , const int DIMS )
   @brief Returns our location when given momenta in the  \f$ -\frac{L_\mu}{2} \rightarrow ( \frac{L_\mu}{2} - 1 ) \f$ BZ.
   @param x :: Momenta
   @param DIMS :: Maximum number of directions of momenta
 **/
int
get_site_pipiBZ( int x[ ND ] , 
		 const int DIMS ) ;

/**
   @fn void get_mom_2piBZ( int x[ ND ] , const int i , const int DIMS )
   @brief Gives the ND-momentum "x" in the \f$ 0 \rightarrow L_\mu \f$ BZ definition when given a lattice index. 
   @param x :: Momenta
   @param i :: Lattice index
   @param DIMS :: Maximum number of directions of momenta
   Useful for translating FFTW's results.
 **/
void 
get_mom_2piBZ( int x[ ND ] , 
	       const int i , 
	       const int DIMS ) ;

/**
   @fn void get_mom_pipi( int x[ ND ] , const int i , const int DIMS )
   @brief Gives the ND-momentum "x" in the \f$ -\frac{L_\mu}{2} \rightarrow ( \frac{L_\mu}{2} - 1 )  \f$ BZ definition. 
   @param x :: Momenta
   @param i :: Lattice index
   @param DIMS :: Maximum number of directions of momenta
   These are often referred to as \f$ \nu's \f$. The integer form of the Fourier modes.
 **/
void
get_mom_pipi( int x[ ND ] , 
	      const int i , 
	      const int DIMS ) ;

/**
   @fn void get_mom_2piBZ_hirep2( int x[ ND ] , const int i )
   @brief ND - generic code to calculate the lattice position vector in the HiRep geometry t,x,y,z. 
   @param x :: ND-momentum
   @param i :: Lattice index
   It's a strange choice of directions, Minkowskian.
 **/
void 
get_mom_2piBZ_hirep2( int x[ ND ] , 
		      const int i ) ;

/**
   @fn void compute_p( GLU_real p[ ND ] , const int n[ ND ] , const int DIMS )
   @brief computes the lattice momentum from the integer n[mu]'s
   @param p :: lattice momentum
   @param n :: Fourier mode
   @param DIMS :: number of dimensions we operate with
 */
void
compute_p( GLU_real p[ ND ] ,
	   const int n[ ND ] ,
	   const int DIMS ) ;

/**
   @fn GLU_real gen_p_sq( const int i , const int DIMS )
   @brief Returns the momentum-squared of the site using the momentum definition for links:
   @param i :: lattice index
   @param DIMS :: Maximum number of dimensions
   \f$ p^2 = 2 \left( DIMS - \sum_mu \cos( \frac{ x_\mu 2\pi }{ L_\mu } ) \right) \f$
   this uses the momenta in the \f$ 0 \rightarrow L_\mu \f$ BZ, but as it is a p^2
   it is \f$\pi\f$ - periodic so we could use either. This one has fewer
   integer operations.
 **/
GLU_real 
gen_p_sq( const int i , 
	  const int DIMS ) ;

/**
   @fn GLU_real gen_p_sq_imp( const int i , const int DIMS )
   @brief Returns the momentum-squared of the site using the momentum definition for links:
   @param i :: lattice index
   @param DIMS :: Maximum number of dimensions
   \f$ p_\mu = 2 \left( \frac{9}{8}\sin\left(\frac{p_\mu}{2}\right) - \frac{1}{24}\sin\left(\frac{3p_\mu}{2}\right)\right)  \f$

   this uses the momenta in the \f$ 0 -> L_\mu \f$ BZ, but as it is a \f$p^2\f$
   it is Pi - periodic so we could use either. This one has fewer
   integer operations.
 **/
GLU_real 
gen_p_sq_imp( const int i , 
	      const int DIMS ) ;

/**
   @fn GLU_real gen_p_sq_feyn( const int i , int *flag )
   @brief Returns the \f$ p^2 \f$ of a site "i"

   @param i :: lattice index
   @param flag :: Flag if we are at a zero momentum.

   With the caveat that the momenta also satisfy the Feynman gauge condition in momentum space. Used in the U(1)-ification.
 **/
GLU_real 
gen_p_sq_feyn( const int i , 
	       int *flag ) ;

/**
   @fn void gen_get_p( GLU_real p[ ND ] , const int i , const int DIMS )
   @brief writes into "p" the lattice momentum

   @param p :: momentum vector
   @param i :: lattice index
   @param DIMS :: maximum number of directions, zero ones above it
   \f$ p_\mu = 2 \sin \left( \frac{ \pi * x_\mu }{ L_\mu} \right)  \f$ 
**/
void 
gen_get_p( GLU_real p[ ND ] , 
	   const int i , 
	   const int DIMS ) ;

/**
   @fn void get_vec_from_origin( int n[ ND ] , const int i , const int DIMS )
   @brief gives the vector from the lattice origin (0,0,0,0)
   @param n :: Fourier mode vector
   @param i :: lexicographical site index
   @param DIMS :: dimensions of the problem
 */
void
get_vec_from_origin( int n[ ND ] , 
		     const int i , 
		     const int DIMS ) ;

/**
   @fn int compute_rsq( const int site , const int dir )
   @brief given a lattice site, returns distance from origin (0,0,0,..,0)
   @param site :: lexicographical site index
   @param DIMS :: maximum number of dimensions, i.e. ND-1 would be in the spatial subcube
   @return the integer valued r-squared
 */
int
compute_rsq( const int site , 
	     const int dir ) ;

/**
   @fn int gen_shift( const int i , const int dir )
   @brief Returns a site index shifted by one in the direction "dir"
   @param i :: lattice index
   @param dir :: direction of the shift
 **/
int
gen_shift( const int i , 
	   const int dir ) ;

/**
   @fn void TwoPI_mpipi_momconv( int MOM[ ND ] , const int i , const int DIR )
   @brief short function to convert momenta from the \f$ 0\rightarrow 2\pi \f$ BZ into the conventional (at least for lattice) \f$-\pi \rightarrow \pi \f$ BZ
   @param MOM :: vector of integer represented momentum that will be changed
   @param i :: lattice index for the \f$ 0\rightarrow 2\pi \f$ BZ
   @param DIR :: maximum number of applicable directions
*/
void
TwoPI_mpipi_momconv( int MOM[ ND ] , 
		     const int i , 
		     const int DIR ) ;

/**
   @fn void init_navig( struct site *__restrict lat )
   @brief Function for generically initialising the lattice navigation
   
   @param lat :: Gauge field

   packs the struct lat's "neighbor" and
   "back" members which will be used for numerical
   derivatives and alike.
 **/
void 
init_navig( struct site *__restrict lat ) ;

/**
   @fn static void init_geom( void )
   @brief initialise vital constants
   Inits the constants #LSQ, #LCU and #LVOLUME
 */
void
init_geom( void ) ;

#endif
