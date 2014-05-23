TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

lib_LIBRARIES = libGLU.a

## just include the header directory
AM_CFLAGS=-I${TOPDIR}/src/Headers/

include_HEADERS = \
	 ./Headers/3Dcuts.h ./Headers/4D_fast.h ./Headers/adaptive_flow.h ./Headers/BPST_config.h ./Headers/CFACG.h ./Headers/CG.h ./Headers/chklat_stuff.h ./Headers/clover.h ./Headers/config_gluons.h ./Headers/Coulomb.h ./Headers/crc.h ./Headers/cut_output.h ./Headers/cut_routines.h ./Headers/cuts.h ./Headers/CUT_wrap.h ./Headers/Doxy_mainpage.h ./Headers/effs.h ./Headers/exactQ.h ./Headers/expMat.h ./Headers/FACG.h ./Headers/geometry.h ./Headers/gftests.h ./Headers/GF_wrap.h ./Headers/givens.h ./Headers/GLU_bswap.h ./Headers/GLU_definitions.h ./Headers/GLU_enum.h ./Headers/glueprop.h ./Headers/GLU.h ./Headers/GLUlib_wrap.h ./Headers/GLU_memcheck.h ./Headers/GLU_rng.h ./Headers/GLU_splines.h ./Headers/GLU_timer.h ./Headers/GLU_types.h ./Headers/gramschmidt.h ./Headers/gtrans.h ./Headers/HIREP.h ./Headers/HYP_4D.h ./Headers/HYP.h ./Headers/input_help.h ./Headers/input_reader.h ./Headers/invert.h ./Headers/KISS_rng.h ./Headers/Landau.h ./Headers/lie_mats.h ./Headers/lin_derivs.h ./Headers/log_derivs.h ./Headers/MAG.h ./Headers/Mainfile.h ./Headers/MMULdagdag.h ./Headers/MMUL_dag.h ./Headers/MMULdag.h ./Headers/MMUL.h ./Headers/MOMggg.h ./Headers/MOMgg.h ./Headers/MWC_1038.h ./Headers/MWC_4096.h ./Headers/ND_generic_HYP.h ./Headers/New_HYP_4D.h ./Headers/OBS_wrap.h ./Headers/plan_ffts.h ./Headers/plaqs_links.h ./Headers/POLY.h ./Headers/projectors.h ./Headers/pspace_landau.h ./Headers/Qsusc.h ./Headers/random_config.h ./Headers/read_config.h ./Headers/readers.h ./Headers/read_headers.h ./Headers/Scidac.h ./Headers/smear.h ./Headers/smearing_param.h ./Headers/SM_wrap.h ./Headers/solver.h ./Headers/staples.h ./Headers/SUNCxU1_config.h ./Headers/taylor_logs.h ./Headers/triplet_gen.h ./Headers/U1_obs.h ./Headers/U1_top.h ./Headers/U_Nops.h ./Headers/vandermonde.h ./Headers/well_rng_19937a.h ./Headers/wflowfuncs.h ./Headers/wflow.h ./Headers/write_headers.h ./Headers/writers.h ./Headers/XML_info.h \
	../config.h

libGLU_a_SOURCES = \
	 ./Cuts/3Dcuts.c ./Cuts/config_gluons.c ./Cuts/cut_output.c ./Cuts/cut_routines.c ./Cuts/cuts.c ./Cuts/glueprop.c ./Cuts/MOMgg.c ./Cuts/MOMggg.c ./Cuts/pspace_landau.c ./Cuts/smearing_param.c ./Cuts/triplet_gen.c ./Field/clover.c ./Field/plaqs_links.c ./Field/POLY.c ./Field/Qsusc.c ./Geometry/BPST_config.c ./Geometry/geometry.c ./Geometry/plan_ffts.c ./Geometry/random_config.c ./IO/chklat_stuff.c ./IO/crc.c ./IO/HIREP.c ./IO/input_help.c ./IO/input_reader.c ./IO/read_config.c ./IO/readers.c ./IO/read_headers.c ./IO/Scidac.c ./IO/write_headers.c ./IO/writers.c ./IO/XML_info.c ./Landau/CFACG.c ./Landau/CG.c ./Landau/Coulomb.c ./Landau/FACG.c ./Landau/gftests.c ./Landau/gtrans.c ./Landau/Landau.c ./Landau/lin_derivs.c ./Landau/log_derivs.c ./Landau/MAG.c ./Matrix_Ops/effs.c ./Matrix_Ops/exactQ.c ./Matrix_Ops/expMat.c ./Matrix_Ops/givens.c ./Matrix_Ops/gramschmidt.c ./Matrix_Ops/invert.c ./Matrix_Ops/lie_mats.c ./Matrix_Ops/MMUL.c ./Matrix_Ops/MMUL_dag.c ./Matrix_Ops/MMULdag.c ./Matrix_Ops/MMULdagdag.c ./Matrix_Ops/solver.c ./Matrix_Ops/taylor_logs.c ./Matrix_Ops/U_Nops.c ./Matrix_Ops/vandermonde.c ./Smear/4D_fast.c ./Smear/adaptive_flow.c ./Smear/HYP_4D.c ./Smear/HYP.c ./Smear/ND_generic_HYP.c ./Smear/projectors.c ./Smear/smear.c ./Smear/staples.c ./Smear/wflow.c ./Smear/wflowfuncs.c ./U1/SUNCxU1_config.c ./U1/U1_obs.c ./U1/U1_top.c ./Utils/GLU_bswap.c ./Utils/GLU_memcheck.c ./Utils/GLU_rng.c ./Utils/GLU_splines.c ./Utils/GLU_timer.c ./Utils/KISS_rng.c ./Utils/MWC_1038.c ./Utils/MWC_4096.c ./Utils/well_rng_19937a.c ./Wrappers/CUT_wrap.c ./Wrappers/GF_wrap.c ./Wrappers/GLUlib_wrap.c ./Wrappers/OBS_wrap.c ./Wrappers/SM_wrap.c

if !PREF

Bindir = "${prefix}"/bin

bin_PROGRAMS = GLU

GLU_SOURCES = ./Run/Mainfile.c
GLU_CFLAGS = ${CFLAGS} -I${TOPDIR}/src/Headers/
GLU_LDADD = libGLU.a ${LDFLAGS}

endif