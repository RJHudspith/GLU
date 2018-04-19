/*
    Copyright 2013-2018 Renwick James Hudspith

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
 **/
#ifndef GLU_TYPES_H
#define GLU_TYPES_H

#include "Mainfile.h"

/**
   @struct draughtboard
   @param square :: list of rg(b) checkerboard
   @param Nsquare :: number of sites in each square index
   @param Ncolors :: number of colors we have on our draughtboard
 */
struct draughtboard {
  size_t **square ;
  size_t *Nsquare ;
  size_t Ncolors ;
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
   @struct CGtemps
   @brief coulomb gauge fixing temporaries
   @param red :: reduction array
   @param sn :: conjugate directions
   @param in_old :: old derivative direction
   @param rotato :: gauge rotated links
   @param db :: draughtboard
 */
struct CGtemps {
  double *red ;
  GLU_complex **sn ;
  GLU_complex **in_old ;
  struct s_site *rotato ;
  struct draughtboard db ;
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
   @struct fftw_small_stuff
   @brief small fft plan storage
 */
struct fftw_small_stuff {
  GLU_complex *out ;
  GLU_complex *in ;
#ifdef HAVE_FFTW3_H
  fftw_plan forward ;
  fftw_plan backward ;
#else
  int forward ;
  int backward ;
#endif
  GLU_real *psq ;
} ;

/**
   @struct fftw_stuff
   @brief FFT storage
 */
struct fftw_stuff {
  GLU_complex **out ;
  GLU_complex **in ;
#ifdef HAVE_FFTW3_H
  fftw_plan *forward ;
  fftw_plan *backward ;
#else
  int *forward ;
  int *backward ;
#endif
  GLU_real *psq ;
} ;

/**
   @struct gauges
   @brief temporary gauge rotation matrices for Coulomb gauge fixing
   @param g :: this slice's gauge transformation matrices
   @param g_up :: the next slice's gauge transformation matrices
   @param g_end :: the final slice's gauge transformation matrices
 */
struct gauges {
  GLU_complex **g ;
  GLU_complex **g_up ;
  GLU_complex **g_end ;
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
   @struct hb_info
   @brief heatbath information storage
 */
struct hb_info {
  double beta ;
  size_t iterations ;
  size_t Nmeasure ;
  size_t Nor ;
  size_t Nsave ;
  size_t therm ;
  GLU_bool continuation ;
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
   @struct infile_data 
   @brief one struct to rule them all 
 */
struct infile_data {
  struct cut_info CUTINFO ;
  struct gf_info GFINFO ;
  struct hb_info HBINFO ;
  struct sm_info SMINFO ;
  struct u1_info U1INFO ;
  header_mode head ;
  GLU_mode mode ;
  GLU_bool rtrans ;
  char output_details[ 64 ] ;
  GLU_output storage ;
} ;

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
  uint32_t Nthreads ; // number of threads
  struct su2_subgroups *su2_data ; // su2 subgroups
  cline_arg argc ; // command line arguments
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
} ; 

/**
   @struct Qmoments
   @brief storage of topological moments
 */
struct Qmoments {
  double *Q ;
  double *Q2 ;
} ;

/**
   @struct site
   @brief the gauge field format
 */
struct site {
  GLU_complex **O ;
  size_t neighbor[ ND ] ;
  size_t back[ ND ] ;
} ;

/**
   @struct s_site
   @brief new format
 */
struct s_site {
  GLU_complex **O ;
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
   @struct veclist
   @brief storage for the momenta
 */
struct veclist {
  size_t idx ;
  int MOM[ ND ] ;
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
   @struct wflow_temps
   @brief struct for storing wflow temporaries
   @param lat2 :: time-slice wide temporary
   @param lat3 :: time-slice wide temporary
   @param lat3 :: time-slice wide temporary
   @param Z :: generating functional
   @param lat_two :: copy of the gauge links
   @param red :: reduction array
 */
struct wflow_temps {
  struct s_site *lat2 ;
  struct s_site *lat3 ;
  struct s_site *lat4 ;
  struct s_site *Z ;
  struct site *lat_two ;
  double *red ;
} ;

#endif
