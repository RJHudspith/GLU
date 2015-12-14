/*
    Copyright 2013 Renwick James Hudspith

    This file (GLU_types.h) is part of GLU.

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
   @file GLU_types.h
   @brief includes all of the struct definitions dotted around the code
   I put all of the struct definitions in here ...
 **/

#ifndef GLU_TYPES_H
#define GLU_TYPES_H

#include "Mainfile.h"

/**
   @struct latt_info
   @brief (useful?) lattice information
   @param dims[mu] :: lattice dimensions in c-order, x moves quickest
   @param Volume :: lattice volume
   @param Lcu :: volume with one (slowest moving) dimension removed
   @param Lsq :: volume with two (slowest moving) dimensions removed
   @param flow :: the configuration number
   @param gf_alpha :: the tuning parameter of our gauge fixer. 1.0/12.0 is good.
   @param sm_alpha :: the smearing parameters from the input file
   @param head :: what header type we use
   @param Seed :: the seed we use for our RNG 
 */
struct latt_info {
  size_t dims[ ND ] ; // dimensions in x,y,z,t order opposite to FFTW
  size_t Lsq ; // dims[0] * dims[1] x,y plane
  size_t Lcu ; // dims[2] * Lsq x,y,z cubic volume
  size_t Volume ; // lattice Volume
  size_t flow ; // config number , gets passed around a bit
  double gf_alpha ; //gauge fixing alpha
  double sm_alpha[ ND -1 ] ; // smearing alphas
  header_mode head ;// Which header type are we using
  uint32_t Seed[ 1 ] ; // rng seed, inbuilt KISS uses four of these
  double twiddles[ ND ] ; // fourier transform twiddles
} ;

/**
   @struct QCDheader
   @brief contains the NERSC header information
   This is only used in chklat_stuff.c and should probably be moved there
 */
struct QCDheader {
  int ntoken ; 
  char **token ; 
  char **value ; 
}; 

/**
   @struct site
   @brief the gauge field format
 */
struct site {
  GLU_complex O[ ND ][ NCNC ] ;
  int neighbor[ ND ] ;
  int back[ ND ] ;
} ;

/**
   @struct sp_site
   @brief spatial site format
 */
struct sp_site {
  GLU_complex O[ ND - 1 ][ NCNC ] ;
} ;

/**
   @struct spt_site
   @brief gauge link format
 */
struct spt_site {
  GLU_complex O[ ND ][ NCNC ] ;
} ;

/**
   @struct spt_site_herm
   @brief specifically for the Wilson flow
 */
struct spt_site_herm {
  #if NC == 3
  GLU_complex O[ ND ][ HERMSIZE - 1 ] ;
  #else
  GLU_complex O[ ND ][ HERMSIZE ] ;
  #endif
} ;

/**
   @struct sp_site_herm
   @brief specifically for the CGF
 */
struct sp_site_herm {
  GLU_complex O[ ND-1 ][ HERMSIZE ] ;
} ;

/**
   @struct lv1
   @param the "level1" dressed links for HYP smearing
 */
struct lv1 {
  GLU_complex O[ ND * ( ND - 1 ) ][ NCNC ] ;
} ; 

/**
   @struct smallest_lv1
   @param the "level1" dressed links for HYP smearing shortened using the 8 parameter link definition
 */
struct smallest_lv1 {
  GLU_real O[ ND * ( ND - 1 ) ][ NCNC - 1 ] ;
} ; 

/**
   @struct spatial_lv1
   @param the "level1" dressed links for HYP smearing in spatial directions only
 */
struct spatial_lv1 {
  GLU_complex O[ ( ND - 1 ) * ( ND - 2 ) ][ NCNC ] ;
} ; 

/**
   @struct veclist
   @brief storage for the momenta
 */
struct veclist {
  size_t idx ;
  int MOM[ ND ] ;
} ;

/**
   @struct su2_subgroups
   @brief storage for the indices of su(2) subgroups
*/
struct su2_subgroups {
  size_t idx_a ;
  size_t idx_b ;
  size_t idx_c ;
  size_t idx_d ;
} ;

/**
   @struct wfmeas
   @brief linked list for wilson flow measurements
   @param avplaq :: average plaquette
   @param Gt :: \f$ t^2 \langle G_{\mu\nu} G_{\mu\nu}\rangle \f$ (clover)
   @param qtop :: gauge field topological charge
   @param time :: flow time (has units t^2/a^2 )
   @param next :: pointer to the next element in the list
 */
struct wfmeas {
  double avplaq ;
  double Gt ;
  double qtop ;
  double time ;
  struct wfmeas *next ;
} ;

/**
   @struct gf_info
   @brief gauge fixing information storage
   @param improve :: enumerated, (NONE,MAG,SMEARED_PRECONDITIONED)
   @param max_iters :: maximum iterations before we restart
   @param accuracy :: gauge fixing accuracy we try to meet 1E-20 is good
   @param type :: enumerated (LANDAU,COULOMB)
 */
struct gf_info {
  GF_improvements improve ; // improvement = MAG , SMEARED_PRECONDITIONED
  size_t max_iters ; // maximum iterations of the gauge fixing routine
  double accuracy ; // average accuracy to be used 
  GLU_fixing type ; // type of gauge fixing used, coulomb or landau?
} ;

/**
   @struct cut_info
   @brief cutting information storage
   @param dir :: either spatial or temporal cuts allowed for now
   @param type :: psq,hypercubic,cylinder or conical
   @param max_mom :: maximum allowed p^2 for the vector of ints definition
   @param max_t :: maximum T allowed in static potential
   @param where :: where is our file outputted to?
   @param definition :: are our gauge fields logarithmic or AntiHermitian_projly defined?
   @param angle :: conical angle from the p=0.
   @param cyl_width :: width of the cylinder in lattice units
 */
struct cut_info {
  cut_mode dir ; // overwrite to something more sensible
  momentum_cut_def type ; // enumerated cutting type
  size_t max_mom ; // maximum momentum allowed for the cut
  size_t max_t ; // maximum T for the static potential
  char where[ 256 ] ; // where do we output to ?
  lie_field_def definition ; // LOG or LINEAR definition of gluons
  size_t angle ; // conical cut angle
  double cyl_width ; // cylinder with
} ;

/**
   @struct sm_info
   @brief smearing information storage
   @param dir :: spatial or temporal smearing is allowed
   @param smiters :: number of smearing iterations
   @param type :: can either be (SM_APE,SM_LOG,SM_STOUT)
 */
struct sm_info {
  GLU_smeardir dir ; // direction ND or ND - 1
  size_t smiters ; // number of smearing iterations
  smearing_types type ; // type of smearing
} ;

/**
   @struct head_data
   @brief information taken from the header
   @param endianess :: tells us which end is up
   @param config_type :: NERSC or HiRep?
   @param precision :: what precision the configs were saved in
   @param plaquette :: the read value of the average plaquette
   @param trace :: the read value of the average link trace
   @param checksum :: the read value of the checksum
 */
struct head_data {
  GLU_endian endianess ;
  GLU_output config_type ;
  file_prec precision ;
  double plaquette ;
  double trace ;
  uint32_t checksum ;
  uint32_t checksumb ;
} ;

/**
   @struct u1_info
   @brief information on the U1-ification of the gauge fields

   @param alpha :: the noncompact QED coupling strength
   @param charge :: the sign of the QED charge
   @param meas :: the U(1) measurement being made
 */
struct u1_info {
  double alpha ;
  double charge ;
  U1_meas meas ;    
} ;

/**
   @struct infile_data 
   @brief one struct to rule them all 
 */
struct infile_data {
  struct gf_info GFINFO ;
  struct cut_info CUTINFO ;
  struct sm_info SMINFO ;
  struct u1_info U1INFO ;
  header_mode head ;
  GLU_mode mode ;
  GLU_bool rtrans ;
  char output_details[ 64 ] ;
  GLU_output storage ;
} ;

#endif
