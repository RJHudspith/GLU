/*
    Copyright 2013-2016 Renwick James Hudspith

    This file (Doxy_mainpage.h) is part of GLU.

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
   @mainpage The GLU Library

   @author Renwick James Hudspith (Jamie)

   @section intro_sec Introduction
   
   This library was written due to the need for a fast Landau gauge fixing
   code for Lattice gauge theories.

   As such it uses the technique of Fourier Acceleration for both its Coulomb
   and Landau gauge fixing routines.

   There are also many other measurements and methods written in this code.
   There is an extensive link smearing library and I have incorporated 
   the U(1)-ing of generic SU(N) links.

   @subsection representation of matrices and geometry

   Lattice fields are often special unitary matrices, I represent these matrices
   as a 1D row-major array. i.e. element "4" of a 3x3 matrix is in row,column
   order the element (2,2), the centre of the matrix. One accesses individual
   elements by computing idx = row_index + Nrows * column_index.

   I represent the gauge field in a similar manner, as a 1D array whose indices
   are in lexicographical order with the direction "x" moving fastest then "y"
   then "z" then .... "t". I always call my "time" direction the slowest moving
   index.

   Lattice fields are stored in struct site. This is an object of the order 
   [ \f$mu\f$ ][ elements ] where \f$mu\f$ is the direction the link is pointing in and
   "elements" are the elements of the SU(N) matrix.
*/

/**
   \defgroup CUTS Cuts
   \ingroup CUTS
   @{
   3Dcuts.h <br>
   cut_output.h <br>
   cut_routines.h <br>
   cuts.h <br>
   glueprop.h <br>
   pspace_landau.h <br>
   MOMgg.h <br>
   MOMggg.h <br>
   smearing_param.h <br>
   triplet_gen.h <br>

   3Dcuts.c <br>  
   cut_output.c <br>
   cut_routines.c <br>
   cuts.c <br>
   glueprop.c <br>
   MOMgg.c <br>
   MOMggg.c <br>
   pspace_landau.c <br>
   smearing_param.c <br>
   triplet_gen.c <br>
   @}

   \defgroup GEOMETRY Geometry
   \ingroup GEOMETRY
   @{
   BPST_config.h <br>
   geometry.h <br>
   plan_ffts.h <br>
   random_config.h <br>

   BPST_config.c <br>
   geometry.c <br>
   plan_ffts.c <br>
   random_config.c <br>
   @}

   \defgroup Io IO
   \ingroup Io
   @{
   chklat_stuff.h <br>
   crc.h <br>
   HIREP.h <br>
   input_reader.h <br>
   input_help.h <br>
   read_config.h <br>
   read_headers.h <br>
   readers.h <br>
   Scidac.h <br>
   write_headers.h <br>
   writers.h <br>
   XML_info.h <br>

   chklat_stuff.c <br>
   crc.c <br>
   HIREP.c <br>
   input_reader.c <br>
   input_help.c <br>
   read_config.c <br>
   read_headers.h <br>
   readers.c <br>
   Scidac.c <br>
   write_headers.c <br>
   writers.c <br>
   XML_info.c <br>
   @}

   \defgroup LANDAU Landau
   \ingroup LANDAU
   @{
   CFACG.h <br>
   Coulomb.h <br>
   FACG.h <br>
   gftests.h <br>
   gtrans.h <br>
   Landau.h <br>
   lin_derivs.h <br>
   log_derivs.h <br>
   MAG.h <br>

   CFACG.c <br>
   Coulomb.c <br>
   FACG.c <br>
   gftests.c <br>
   gtrans.c <br>
   Landau.c <br>
   lin_derivs.c <br>
   log_derivs.c <br>
   MAG.c <br>
   @}

   \defgroup MATRIX_OPS Matrix_OPs
   \ingroup MATRIX_OPS
   @{
   effs.h <br>
   exactQ.h <br>
   gramschmidt.h <br>
   invert.h <br>
   lie_mats.h <br>
   MMUL.h <br>
   MMULdag.h <br>
   MMUL_dag.h <br>
   MMULdagdag.h <br>
   solver.h <br>
   taylor_logs.h <br>
   U_Nops.h <br>
   vandermonde.h  <br>

   effs.c <br>
   exactQ.c <br>
   gramschmidt.c <br>
   invert.c <br>
   lie_mats.c <br>
   MMUL.c <br>
   MMULdag.c <br>
   MMUL_dag.c <br>
   MMULdagdag.c <br>
   solver.c <br>
   taylor_logs.c <br>
   U_Nops.c <br>
   vandermonde.c  <br>
   @}

   \defgroup RUN Run
   \ingroup RUN
   @{
   Mainfile.h  <br>

   Mainfile.c  <br>
   @}

   \defgroup SMEAR Smear
   \ingroup SMEAR
   @{
   4D_fast.h  <br>
   adaptive_flow.h  <br>
   HYP_4D.h  <br>
   HYP.h  <br>
   ND_generic_HYP.h <br>
   projectors.h  <br>
   smear.h  <br>
   staples.h  <br>
   wflow.h  <br>
   wflowfuncs.h  <br>

   4D_fast.c  <br>
   adaptive_flow.c  <br>
   HYP_4D.c  <br>
   HYP.c  <br>
   ND_generic_HYP.c <br>
   projectors.c  <br>
   smear.c  <br>
   staples.c  <br>
   wflow.c  <br>
   wflowfuncs.c  <br>
   @}

   \defgroup UONE U1
   \ingroup UONE
   @{
   su3xu1_config.h  <br>
   U1_obs.h     <br>
   U1_top.h     <br>

   su3xu1_config.c  <br>
   U1_obs.c     <br>
   U1_top.c     <br>
   @}

   \defgroup UTILS Utils
   \ingroup UTILS
   @{
   GLU_bswap.h <br>
   GLU_memcheck.h  <br>
   GLU_rng.h <br>
   GLU_splines.h <br>
   GLU_timer.h  <br>
   KISS_rng.h <br>
   MWC_1038.h <br>
   MWC_4096.h <br>
   well_rng_19937a.h <br>

   GLU_bswap.c <br>
   GLU_memcheck.c  <br>
   GLU_rng.c <br>
   GLU_splines.c <br>
   GLU_timer.c  <br>
   KISS_rng.c <br>
   MWC_1038.c <br>
   MWC_4096.c <br>
   well_rng_19937a.c <br>
   @}

   \defgroup WRAPPERS Wrappers
   \ingroup WRAPPERS
   @{
   CUT_wrap.h  <br>
   GF_wrap.h  <br>
   OBS_wrap.h  <br>
   SM_wrap.h  <br>

   CUT_wrap.c  <br>
   GF_wrap.c  <br>
   OBS_wrap.c  <br>
   SM_wrap.c  <br>
   @}
 */
