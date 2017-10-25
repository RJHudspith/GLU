/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (input_help.c) is part of GLU.

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
   @file input_help.c
   @brief provides some help about the input file
 */

#include "Mainfile.h"

// shortcut for equivalent strings
static int
are_equal( const char *str_1 , const char *str_2 ) { 
  return !strcmp( str_1 , str_2 ) ; }

// function for generating example input files
static void
create_input_file( const char *mode_str ,
		   const char *gf_str ,
		   const char *cut_str )
{
  fprintf( stdout , "MODE = %s \n" , mode_str ) ;
  fprintf( stdout , "HEADER = NERSC\n" ) ;
  size_t mu ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    fprintf( stdout , "    DIM_%zu = 16\n" , mu ) ;
  }
  fprintf( stdout , "CONFNO = 0\n" 
	   "RANDOM_TRANSFORM = NO\n" 
	   "SEED = 0\n" ) ;
  fprintf( stdout , "GFTYPE = %s\n" , gf_str ) ;
  fprintf( stdout , "    GF_TUNE = 0.08\n" 
	   "    IMPROVEMENTS = NONE\n" 
	   "    ACCURACY = 14\n"
	   "    MAX_ITERS = 1000\n" ) ;
  fprintf( stdout , "CUTTYPE = %s\n"
	   "    FIELD_DEFINITION = LINEAR\n"
	   "    MOM_CUT = SPHERICAL_CUT\n" 
	   "    MAX_T = 7\n"
	   "    MAXMOM = 80\n" 
	   "    CYL_WIDTH = 2.0\n"
	   "    ANGLE = 60\n"
	   "    OUTPUT = ./\n" , cut_str ) ;
  fprintf( stdout , "SMEARTYPE = WFLOW_STOUT\n" 
	   "    DIRECTION = ALL\n"
	   "    SMITERS = 100\n" ) ;
  fprintf( stdout , "    ALPHA1 = 0.02\n" ) ;
  for( mu = 2 ; mu < ND ; mu++ ) {
    fprintf( stdout , "    ALPHA%zu = 0.0\n" , mu ) ;
  }
  fprintf( stdout , "U1_MEAS = U1_RECTANGLE\n" 
	   "    U1_ALPHA = 0.07957753876221914\n" 
	   "    U1_CHARGE = -1.0\n" ) ;
  fprintf( stdout , "CONFIG_INFO = GLU_config\n" ) ;
  fprintf( stdout , "    STORAGE = NERSC_NCxNC\n" ) ;
  fprintf( stdout , "BETA = 6.0\n" ) ;
  fprintf( stdout , "    ITERS = 1500\n" ) ;
  fprintf( stdout , "    MEASURE = 5\n" ) ;
  fprintf( stdout , "    OVER_ITERS = 4\n" ) ;
  fprintf( stdout , "    SAVE = 25\n" ) ;
  fprintf( stdout , "    THERM = 150\n" ) ;
  return ;
}

// possible cut types
static void
cuttype_types( void )
{
  fprintf( stdout , "CUTTYPE = CONFIGSPACE_GLUE           - "
	   "Computes the configuration space gluon propagator,"
	   "must be gauge fixed \n" ) ;
  fprintf( stdout , "CUTTYPE = EXCEPTIONAL                - "
	   "Computes the gluon propagator g(2) and MOMgg-projected"
	   "gluon 3 point function g(3) to a\n" 
	  "                                       list which "
	  "has the format (all in Big Endian binary format)\n"
	  "                                       NMOM\n"
	  "                                       ND ( p_0 , p_1 , "
	  "... , p_ND-1 )\n" 
	  "                                       NMOM\n"
	  "                                       g(2)( p )\n"
	  "                                       NMOM\n"
	  "                                       g(3)( p )\n" 
	  "                                       where NMOM is the number "
	  "of the lists of momenta, g(2) and g(3)\n" ) ;
  fprintf( stdout , "CUTTYPE = FIELDS                     - Writes out full "
	   " NCxNC momentum-space Lie matrices\n" ) ;
  fprintf( stdout , "CUTTYPE = GLUON_PROPS                - Writes the "
	   "transverse and longitudinal Landau gauge fixed gluon propagators "
	   "in the output form of \n" 
	   "                                       EXCEPTIONAL\n") ;
  fprintf( stdout , "CUTTYPE = INSTANTANEOUS_GLUONS       - Writes the "
	   "momentum-space spatial and temporal gluon propagators of "
	   "Coulomb gauge, the output\n"
	   "                                       format is like that of "
	   "EXCEPTIONAL but with momenta list ND-1 ( p_0 , ... , p_{ND-2}\n" ) ;
  fprintf( stdout , "CUTTYPE = NON_EXCEPTIONAL            - Computes the "
	   "momentum space gluon propagator and 3 point function for gluons "
	   "with momenta \n"
	   "                                       p + r + q = 0, p!=r!=q, "
	   "and p^2=r^2=q^2. Writes to an output file in the format of "
	   "EXCEPTIONAL\n" ) ;
  fprintf( stdout , "CUTTYPE = SMEARED_GLUONS             - Computest the "
	   "unsmeared and smeared Landau gauge fixed propagators, with the "
	   "same output as\n"
	   "                                       that of EXCEPTIONAL but "
	   "with the unsmeared propagator and the smeared propagator written.\n"
	   "                                       This requires the smearing "
	   "arguments to be set in the input file\n" ) ;
  fprintf( stdout , "CUTTYPE = STATIC_POTENTIAL           - Computes the "
	   "correlator of Polyakov loops C(r,T) = < L(0,T) L(r,T) > , \n"
	   "                                       it computes the singlet "
	   "C_1(r,t) = < Tr( L(0,T) L(r,T) ) > and the qq C_qq(r,t)\n"
	   "                                       qq : C_qq(r,t) = "
	   "< Tr( L(0,T) ) Tr( L(r,T) ) >. Storage is Big Endian format\n"
	   "                                       binary data, of the form \n"
	   "                                       NR\n"
	   "                                       ND-1 ( r_0 , r_1 "
	   ", ... , r_ND-2 )\n" 
	   "                                       NT\n"
	   "                                       NMOM\n"
	   "                                       C_1( r , T )\n"
	   "                                       NMOM\n"
	   "                                       C_qq( r , T )\n"
	   "                                       .... NT pairs of V_1 and "
	   "Vqq ...\n"
	   "                                       where NR is the number of "
	   " R-vectors from the origin. This code allows for smearing before\n"
	   "                                       the computation of the "
	   "polyakov loops, specified by the smearing options in the "
	   "input file\n" ) ;
  fprintf( stdout , "CUTTYPE = TOPOLOGICAL_CORRELATOR     - Computes the "
	   "configuration space correlator of topological charges C(r) "
	   " = < q(0) q(r) > \n" 
	   "                                       storage is similar to "
	   "EXCEPTIONAL but instead of momenta there are ND r-vectors "
	   "from the origin\n"
	   "                                       and again singlet and qq"
	   "correlators. Also performs smearing of the gauge field, "
	   " specified by\n"
	   "                                       the smearing options in "
	   "the input file\n" ) ;
  fprintf( stdout , "CUTTYPE = TOPOLOGICAL_SUSCEPTIBILITY - Computes the "
	   "topological susceptibility using the slab method with slabs\n"
	   "in multiple directions.\n"
	   "                                       Also performs smearing "
	   "specified by the smearing options in the input file\n" ) ;
  fprintf( stdout , "CUTTYPE = {ALL ELSE} - Do Nothing \n" ) ;
  fprintf( stdout , "\nCONFIGSPACE_GLUE, EXCEPTIONAL, GLUON_PROPS and "
	   "NON_EXCEPTIONAL require the gauge field to be fixed to "
	   "Landau gauge\n" ) ;
  fprintf( stdout , "\nINSTANTANEOUS_GLUONS and STATIC_POTENTIAL "
	   "require the gauge field to be fixed to Coulomb gauge\n" ) ;
  return ;
}

// available gauge fixing types
static void
gftype_types( void )
{
  fprintf( stdout , "GFTYPE = COULOMB    - Coulomb gauge fixing using "
	   "a time-slice by time-slice Cornell method\n" ) ;
  fprintf( stdout , "GFTYPE = LANDAU     - Landau gauge fixing using "
	   "the Cornell method/or OR depnding on compile flags\n" ) ;
  fprintf( stdout , "GFTYPE = {All ELSE} - Do nothing\n" ) ;
  fprintf( stdout , "\nDepending on what configure flags you have "
	   "specified these will either be Fourier accelerated or not\n" ) ;
  return ;
}

// available headers
static void
header_types( void )
{
  fprintf( stdout , "HEADER = HIREP       - Attempts to read a file "
	   "in HiREP format \n" ) ; 
  fprintf( stdout , "       = ILDG_SCIDAC - Attempts to read an ILDG "
	   "configuration in SCIDAC form and compares checksums \n" ) ;
  fprintf( stdout , "       = ILDG_BQCD   - Attempts to read an ILDG "
	   "configuration in BQCD's format and compares checksums \n" ) ;
  fprintf( stdout , "       = INSTANTON   - Generates a BPST instanton field "
	   "with dimensions specified by DIM_ array\n" ) ;
  fprintf( stdout , "       = LIME        - Attempts to read a lime "
	   "configuration *caution* does not bother with checksums!\n" ) ;
  fprintf( stdout , "       = MILC        - Attempts to read a MILC "
	   "configuration\n" ) ;
  fprintf( stdout , "       = NERSC       - Attempts to read a NERSC "
	   "configuration\n" ) ;
  fprintf( stdout , "       = RANDOM      - Generates a random configuration "
	   "with dimensions specified by DIM_ array\n" ) ;
  fprintf( stdout , "       = UNIT        - Generates an identity-matrix "
	   "configuration with dimensions specified by DIM_ array\n" ) ;
  fprintf( stdout , "\nAll of these configuration types have been "
	   "generalised to their obvious NC and ND extensions\n" ) ;
  return ;
}

// typical usage
static void
help_usage( void )
{
  fprintf( stdout , "\nHelp usage : \n\n./GLU --help={options}\n\n" ) ;
  fprintf( stdout , "The available options are:\n" ) ;
  fprintf( stdout , "MODE, HEADER, DIM , CONFNO, RANDOM_TRANSFORM, \n" ) ;
  fprintf( stdout , "GFTYPE, GF_TUNE, IMPROVEMENTS, ACCURACY, MAX_ITERS, \n" ) ;
  fprintf( stdout , "CUTTYPE, FIELD_DEFINITION, MOM_CUT, MAX_T, MAXMOM, "
	            "CYL_WIDTH, ANGLE, OUTPUT,\n" ) ;
  fprintf( stdout , "SMEARTYPE, DIRECTION, SMITERS, ALPHA,\n" ) ;
  fprintf( stdout , "U1_MEAS, U1_ALPHA, U1_CHARGE,\n" ) ;
  fprintf( stdout , "CONFIG_INFO, STORAGE,\n" ) ;
  fprintf( stdout , "*caution* in the {input_file} each one of these have "
	   "to be specified\n"
	  "          ONCE AND ONLY ONCE!\n" ) ;
  fprintf( stdout , "\nIf you would like an example input file, try:\n"
	   "\n./GLU --autoin={options}\n\n"
	   "Where options can be COULOMB, LANDAU, STATIC_POTENTIAL, "
	   "SUNCxU1, TOPOLOGICAL_SUSCEPTIBILITY or WFLOW\n\n" ) ;
  return ;
}

// possible gauge fixing improvement types
static void
improvement_types( void )
{
  fprintf( stdout , "IMPOVEMENTS = MAG        - Fixes to the so-called "
	   "maximal axial gauge before fixing to Landau or Coulomb\n" ) ;
  fprintf( stdout , "            = SMEAR      - Uses the smeared "
	   "preconditioning of Hettrick and de forcrand using the smearing\n"
	  "                           routines specified in the input file\n") ;
  fprintf( stdout , "            = RESIDUAL   - Fixes the remaining temporal "
	   "degree of freedom after Coulomb gauge fixing\n" ) ;
  fprintf( stdout , "            = {ALL ELSE} - Do nothing\n" ) ;
  return ;
}

// standard modes available
static void
mode_types( void ) 
{
  fprintf( stdout , "MODE = CUTTING      - Static Potential, Qsusc correlator "
	   "and momentum-space gluon correlators\n" ) ;
  fprintf( stdout , "     = GAUGE_FIXING - Landau and Coulomb gauge fixing "
	   "mode selection\n" ) ;
  fprintf( stdout , "     = HEATBATH     - Heatbath algorithm for config "
	   "generation \n" ) ;
  fprintf( stdout , "     = SMEARING     - Link smearing mode, overwrites the "
	   "lattice field with the chosen smearing\n" ) ;
  fprintf( stdout , "     = SUNCxU1      - Quenched SU(N)xU(1) configuration "
	   "generation\n" ) ;
  fprintf( stdout , "     = {ALL ELSE}   - Default behaviour, computes "
	   "plaquettes and Polyakov loops\n" ) ;
  return ;
}

// possible momentum cuts available
static void
momcut_types( void )
{
  fprintf( stdout , "MOM_CUT = CONICAL_CUT    - Retain only momenta within a "
	   "cylinder (that have n^{2} < MAX_MOM) of width CYL_WIDTH that lie\n"
	   "                           along one of the principle diagonals of "
	   "the momentum-space lattice, also within a cone\n"
	   "                           whose angle at the apex is specified "
	   "by ANGLE\n" ) ;
  fprintf( stdout , "        = SPHERICAL_CUT  - Retain only momenta that have "
	   "n^{2} < MAX_MOM\n" ) ;
  fprintf( stdout , "        = CYLINDER_CUT   - Retain only momenta within a "
	   "cylinder (that have n^{2} < MAX_MOM) of width CYL_WIDTH that lie\n"
	  "                           along one of the principle diagonals of "
	   "the momentum-space lattice\n" ) ;
  fprintf( stdout , "        = HYPERCUBIC_CUT - Include individual momenta "
	   "that satisfy |p| < sqrt(MAX_MOM) \n" ) ;
  fprintf( stdout , "        = {ALL ELSE}     - Do nothing \n" ) ;
  fprintf( stdout , "\nn^{2} is (anisotropy corrected) n_{mu}n_{mu} Fourier "
	   "mode product, for {mu} being polarisations on the lattice\n" ) ;
  return ;
}

// smearing possibilities
static void
smeartype_types( void )
{
  fprintf( stdout , "SMEARTYPE = ADAPTWFLOW_LOG   - Adaptive RK3 integration "
	   "of the flow equation using LOG links\n" ) ;
  fprintf( stdout , "          = ADAPTWFLOW_STOUT - Adaptive RK3 integration "
	   "of the flow equation using STOUT links\n" ) ;
  fprintf( stdout , "          = APE              - APE smearing U' = ( 1 - "
	   "\\alpha ) U + \\sum staples \n" ) ;
  fprintf( stdout , "          = HEX              - Hypercubically blocked, "
	   "STOUT smeared links\n" ) ;
  fprintf( stdout , "          = HYL              - Hypercubically blocked, "
	   "LOG smeared links\n" ) ;
  fprintf( stdout , "          = HYP              - Hypercubically blocked, "
	   "APE smeared links\n" ) ;
  fprintf( stdout , "          = LOG              - Logarithmic link smearing "
	   "uses the exact log definition of gauge fields\n" ) ;
  fprintf( stdout , "          = STOUT            - STOUT link smearing uses "
	   "the Hermitian projection definition of gauge fields\n" ) ;
  fprintf( stdout , "          = WFLOW_LOG        - RK3 integration of the "
	   "flow equation using LOG links\n" ) ;
  fprintf( stdout , "          = WFLOW_STOUT      - RK3 integration of the "
	   "flow equation using STOUT links\n" ) ;
  fprintf( stdout , "          = {ALL ELSE}       - Do nothing\n" ) ;
  fprintf( stdout , "\n*caution* Hypercubically blocked variants for ND > 4 "
	   "use a very slow recursive definition\n" ) ;
  fprintf( stdout , "\n*caution* Log smearing for NC > 3 requires slow series "
	   "evaluations and can become unstable\n" ) ;
  return ;
}

// gauge configuration storage possibilities
static void
storage_types( void )
{
  fprintf( stdout , "STORAGE = HIREP       - Writes out a HiREP "
	   "format file\n" ) ;
  fprintf( stdout , "        = ILDG        - Writes a SCIDAC_ILDG "
	   "configuration file\n" ) ;
  fprintf( stdout , "        = MILC        - Writes a MILC "
	   "configuration file\n" ) ;
  fprintf( stdout , "        = NERSC_SMALL - (Only for NC < 4! ) writes the"
	   "smallest possible configuration -- less stable than others\n" ) ;
  fprintf( stdout , "        = NERSC_GAUGE - Writes out the top NC-1 rows "
	   "of the link matrices\n" ) ;
  fprintf( stdout , "        = NERSC_NCxNC - Writes out the full "
	   "link matrices\n" ) ;
  fprintf( stdout , "        = SCIDAC      - Writes a SCIDAC format "
	   "configuration\n" ) ;
  fprintf( stdout , "        = {ALL ELSE}  - Defaults to the NERSC_NCxNC "
	   "format\n" ) ;
  fprintf( stdout , "\n*caution* our MILC configurations are not *actually* "
	   "MILC standard, as they don't have a double precision\n"
	   "           output at the moment. We will see if that changes.\n" ) ;
  fprintf( stdout , "\n*caution* if SUNCxU1 is the chosen MODE then "
	   "NERSC_GAUGE and NERSC_SMALL will not be used\n" ) ;
  return ;
}

// U(1) possible measurements
static void
U1meas_types( void ) 
{
  fprintf( stdout , "U1_MEAS = U1_RECTANGLE   - measures the non-compact "
	   "and compact U(1) average trace of the 2x1 Wilson loop\n" ) ;
  fprintf( stdout , "        = U1_TOPOLOGICAL - measures some non-compact "
	   "topological things such as the Dirac sheet and Dirac string\n" ) ;
  fprintf( stdout , "        = {ALL ELSE}     - measures the non-compact "
	   "and compact U(1) average trace of the plaquette\n" ) ;
  fprintf( stdout , "\nThe SU(NC) gauge fields will be overwritten by "
	   "SU(NC)xU(1), this means some configuration output types "
	   "become unavailable -- only full matrices allowed!\n" ) ;
  return ;
}

// yep, that is what I chose to call it : Help functions
int
GLU_helps_those_who_help_themselves( const char *help_str )
{
  if( are_equal( help_str , "--help=MODE" ) ) {
    mode_types( ) ;
  } else if( are_equal( help_str , "--help=HEADER" ) ) {
    header_types( ) ;
  } else if( are_equal( help_str , "--help=DIM" ) ) {
    fprintf( stdout , "DIM_%%d = %%d - Specified lattice dimensions for the "
	     "UNIT, RANDOM and INSTANTON\n"
	     "              types of HEADER. If we are not using these, "
	     "lattice\n"
	     "              dimensions are taken from the configuration "
	     "header\n" ) ;
  } else if( are_equal( help_str , "--help=CONFNO" ) ) {
    fprintf( stdout , "CONFNO = %%d - User-specified configuration number. "
	     "All configuration types except NERSC do not directly provide\n"
	     "              a configuration number. So you should supply one "
	     "for the output configuration\n" ) ;
  } else if( are_equal( help_str , "--help=RANDOM_TRANSFORM" ) ) {
    fprintf( stdout , "RANDOM_TRANSFORM = YES    - Perform a lattice-wide "
	     "random gauge transformation of the fields\n"
	     "                 = {!YES} - Do nothing\n" ) ;
  } else if( are_equal( help_str , "--help=SEED" ) ) {
    fprintf( stdout , "SEED = %%d - User specified RNG seed \n" 
	     "        0 - Generates a SEED from urandom if it is available\n") ;
  } else if( are_equal( help_str , "--help=GFTYPE" ) ) {
    gftype_types( ) ;
  } else if( are_equal( help_str , "--help=GFTUNE" ) ) {
    fprintf( stdout , "GFTUNE = %%lf - User specified tuning parameter for "
	     "the gauge fixing this should generally be < 0.1, less "
	     "important for the default CG routines\n" ) ;
  } else if( are_equal( help_str , "--help=IMPROVEMENTS" ) ) {
    improvement_types( ) ;
  } else if( are_equal( help_str , "--help=ACCURACY" ) ) {
    fprintf( stdout , "ACCURACY = %%d - Stops the gauge fixing after an "
	     "average accuracy of better than 10^{-ACCURACY} has been "
	     "achieved\n" ) ;
  } else if( are_equal( help_str , "--help=MAX_ITERS" ) ) {
    fprintf( stdout , "MAX_ITERS = %%d - After this many iterations the "
	     "routine restarts from the beginning with a random transformation "
	     "of the\n"
	     "                 fields and after %d restarts it complains about "
	     "changing the tuning\n" , GF_GLU_FAILURES ) ;
  } else if( are_equal( help_str , "--help=CUTTYPE" ) ) {
    cuttype_types( ) ;
  } else if( are_equal( help_str , "--help=FIELD_DEFINTION" ) ) {
    fprintf( stdout , "FIELD_DEFINITION = LOG        - Exact logarithm "
	     "definition of the gauge fields\n" ) ;
    fprintf( stdout , "FIELD_DEFINITION = {ALL ELSE} - Hermitian projection "
	     "definition of the gauge fields (denoted LINEAR)\n" ) ;
  } else if( are_equal( help_str , "--help=MOM_CUT" ) ) {
    momcut_types( ) ;
  } else if( are_equal( help_str , "--help=MAX_T" ) ) {
    fprintf( stdout , "MAX_T = %%d - serves as the maximum T separation for "
	     "the Polyakov lines\n"
	     "             in the STATIC_POTENTIAL code. Computes "
	     "T=1,2,...,MAX_T separations\n" ) ;
  } else if( are_equal( help_str , "--help=MAXMOM" ) ) {
    fprintf( stdout , "MAXMOM = %%d - maximum allowed n_{mu}n_{mu} "
	     "where n_{mu} is a Fourier mode\n" ) ;
  } else if( are_equal( help_str , "--help=CYL_WIDTH" ) ) {
    fprintf( stdout , "CYL_WIDTH = %%lf - momenta are allowed within the "
	     "body-diagonal cylinder of width CYL_WIDTH * 2 \\pi / L \n"
	     "                  where L is the smallest lattice direction's "
	     "length\n" ) ;
  } else if( are_equal( help_str , "--help=ANGLE" ) ) {
    fprintf( stdout , "ANGLE = %%lf - conical cut's angle of the apex of "
	     "the cone\n" ) ;
  } else if( are_equal( help_str , "--help=OUTPUT" ) ) {
    fprintf( stdout , "OUTPUT = %%s - prefix destination for where the "
	     "output file for the cut procedure will be written. The actual "
	     "output file\n"
	     "              is generated within the code, the output format "
	     "requires CONFNO to be set \n" ) ;
  } else if( are_equal( help_str , "--help=SMEARTYPE" ) ) {
    smeartype_types( ) ;
  } else if( are_equal( help_str , "--help=DIRECTION" ) ) {
    fprintf( stdout , "DIRECTION = SPATIAL    - smears the fields for each "
	     "time-slice only in the ND-1 polarisation directions\n" ) ;
    fprintf( stdout , "          = {ALL ELSE} - fully smears all the links "
	     "in all the directions\n" ) ;
  } else if( are_equal( help_str , "--help=SMITERS" ) ) {
    fprintf( stdout , "SMITERS = %%d - maximum number of smearing iterations "
	     "to perform, some routines such as the Wilson flow will "
	     "often finish\n"
	     "               before this number is reached\n" ) ;
  } else if( are_equal( help_str , "--help=ALPHA" ) ) {
    fprintf( stdout , "ALPHA_%%d = %%lf - Smearing parameters for each "
	     "level of smearing, expects ND-1 of these\n"
	     "                 e.g. ALPHA_1, ALPHA_2 .. ALPHA_ND-1 "
	     "for Hypercubically blocked smearing in\n"
	     "                 all directions. For the Wilson flow and for "
	     "APE, LOG and STOUT smearings ALPHA_1 is the only relevant\n"
	     "                 parameter and all others are ignored\n" ) ;
  } else if( are_equal( help_str , "--help=U1_MEAS" ) ) {
    U1meas_types( ) ;
  } else if( are_equal( help_str , "--help=U1_ALPHA" ) ) {
    fprintf( stdout , "U1_ALPHA = %%lf - The non-compact bare coupling "
	     "g^{2} / 4\\pi, beta is 1 / ( \\pi ND U1_ALPHA ) \n" ) ;
    fprintf( stdout , "\n*caution* FFTW must be linked\n" ) ; 
  } else if( are_equal( help_str , "--help=U1_CHARGE" ) ) {
    fprintf( stdout , "U1_CHARGE = %%lf - Charges the U(1) fields upon "
	     "compactification \n" 
	     "                  U(1)_{mu} = exp( i U1_CHARGE "
	     "sqrt( 4 \\pi U1_ALPHA ) A_{mu} ) \n" ) ;
  } else if( are_equal( help_str , "--help=CONFIG_INFO" ) ) {
    fprintf( stdout , "CONFIG_INFO = %%s - Provide some small details "
	     "for when we write out the configuration file\n" ) ;
  } else if( are_equal( help_str , "--help=STORAGE" ) ) {
    storage_types( ) ;
  } else if( are_equal( help_str , "--help=BETA" ) ) {
    fprintf( stdout , "BETA = %%f - the parameter 2N/g_0^2 with which "
	     "we weight the ensembles in the heatbath\n" ) ;
  } else if( are_equal( help_str , "--help=ITERS" ) ) {
    fprintf( stdout , "ITERS = %%d - the total number of iterations after "
	     "thermalisation that the HB-OR does\n" ) ;
  } else if( are_equal( help_str , "--help=MEASURE" ) ) {
    fprintf( stdout , "MEASURE = %%d - the number of combined HB-OR iters "
	     "before a measurement of the plaquette and polyakov loops\n" ) ;
  } else if( are_equal( help_str , "--help=OVER_ITERS" ) ) {
    fprintf( stdout , "OVER_ITERS = %%d - the number of overrelaxation "
	     "iterations in the combined HB-OR\n" ) ;
  } else if( are_equal( help_str , "--help=SAVE" ) ) {
    fprintf( stdout , "SAVE = %%d - the iteration count at which we save "
	     "a configuration in the HB-OR\n" ) ;
  } else if( are_equal( help_str , "--help=THERMALISATION" ) ) {
    fprintf( stdout , "THERMALISATION = %%d - the number of HB-OR iterations "
	     "of thermalisation before measurement and saving\n" ) ;    
  } else if( are_equal( help_str , "--autoin=LANDAU" ) ) {
    create_input_file( "GAUGE_FIXING" , "LANDAU" ,
		       "TOPOLOGICAL_SUSCEPTIBILITY" ) ;
  } else if( are_equal( help_str , "--autoin=COULOMB" ) ) {
    create_input_file( "GAUGE_FIXING" , "COULOMB" ,
		       "TOPOLOGICAL_SUSCEPTIBILITY" ) ;
  } else if( are_equal( help_str , "--autoin=HEATBATH" ) ) {
    create_input_file( "HEATBATH" , "COULOMB" ,
		       "TOPOLOGICAL_SUSCEPTIBILITY") ;   
  } else if( are_equal( help_str , "--autoin=STATIC_POTENTIAL" ) ) {
    create_input_file( "CUTTING" , "LANDAU" ,
		       "STATIC_POTENTIAL" ) ;
  } else if( are_equal( help_str , "--autoin=SUNCxU1" ) ) {
    create_input_file( "SUNCxU1" , "LANDAU" ,
		       "TOPOLOGICAL_SUSCEPTIBILITY" ) ;
  } else if( are_equal( help_str , "--autoin=WFLOW" ) ) {
    create_input_file( "SMEARING" , "LANDAU" ,
		       "TOPOLOGICAL_SUSCEPTIBILITY" ) ;
  } else {
    fprintf( stdout , "[IO] Unrecognised {input_file} query \"%s\" for\n" , 
	     help_str ) ;
    help_usage() ;
  }
  return GLU_SUCCESS ;
}

// simple usage information
int
GLUsage( void )
{
  fprintf( stdout , "\nTo run the code use (the last command is optional):\n\n"
	   "./GLU -i {input_file} -c {config_file} "
	   "-o {output_config_name} \n" ) ;
  fprintf( stdout , "\nFor help on various {input_file} options use:\n\n"
	   "./GLU --help={input_file option}\n\n" ) ;
  fprintf( stdout , "To automatically generate a standard input file use:\n\n"
	   "./GLU --autoin={options}\n"
	   "\nWhere {options} can be COULOMB, HEATBATH , LANDAU, "
	   "STATIC_POTENTIAL, SUNCxU1, WFLOW\n\n" ) ;
  return fprintf( stdout , "If using the CG gauge fixing, "
		  "please cite my paper\n"
		  "\"Fourier Accelerated Conjugate Gradient Lattice "
		  "Gauge Fixing\"\n\n"
		  "@article{Hudspith:2014oja,\n"
		  "author         = \"Hudspith, R.J.\",\n"
		  "title          = \"{Fourier Accelerated Conjugate "
		  "Gradient Lattice Gauge Fixing}\",\n"
		  "collaboration  = \"RBC, UKQCD\",\n"
		  "journal        = \"Comput.Phys.Commun.\",\n"
		  "volume         = \"187\",\n"
		  "pages          = \"115-119\",\n"
		  "doi            = \"10.1016/j.cpc.2014.10.017\",\n"
		  "year           = \"2014\",\n"
		  "eprint         = \"1405.5812\",\n"
		  "archivePrefix  = \"arXiv\",\n"
		  "primaryClass   = \"hep-lat\",\n"
		  "SLACcitation   = \"%%%%CITATION = ARXIV:1405.5812;"
		  "%%%%\",\n}\n" ) ;
}
