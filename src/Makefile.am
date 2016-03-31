TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

lib_LIBRARIES = libGLU.a

## just include the header directory
AM_CFLAGS=-I${TOPDIR}/src/Headers/

## momentum and config space gluon propagator codes
CUTSFILES=./Cuts/3Dcuts.c ./Cuts/config_gluons.c ./Cuts/cut_output.c \
	./Cuts/cut_routines.c ./Cuts/cuts.c ./Cuts/glueprop.c \
	./Cuts/MOMgg.c ./Cuts/MOMggg.c ./Cuts/pspace_landau.c \
	./Cuts/smearing_param.c ./Cuts/triplet_gen.c 

## plaquettes, traces, field strength tensor
FIELDFILES=./Field/clover.c ./Field/plaqs_links.c ./Field/POLY.c \
	./Field/Qsusc.c

## geometry, fft planning
GEOMETRYFILES=./Geometry/BPST_config.c ./Geometry/geometry.c \
	./Geometry/plan_ffts.c ./Geometry/random_config.c

## file readers and writers
IOFILES=./IO/chklat_stuff.c ./IO/crc.c ./IO/HIREP.c ./IO/input_help.c \
	./IO/input_reader.c ./IO/read_config.c ./IO/readers.c \
	./IO/read_headers.c ./IO/Scidac.c ./IO/write_headers.c \
	./IO/writers.c ./IO/XML_info.c

## Gauge-fixing codes
LANDAUFILES=./Landau/CFACG.c ./Landau/CG.c ./Landau/Coulomb.c ./Landau/FACG.c \
	 ./Landau/gftests.c ./Landau/gtrans.c ./Landau/Landau.c \
	./Landau/lin_derivs.c ./Landau/log_derivs.c ./Landau/MAG.c \
	./Landau/Or.c

## home-spun matrix operations
MATRIXOPSFILES=./Matrix_Ops/effs.c ./Matrix_Ops/evalues.c \
	./Matrix_Ops/exactQ.c ./Matrix_Ops/expMat.c \
	./Matrix_Ops/givens.c \
	./Matrix_Ops/gramschmidt.c ./Matrix_Ops/gramschmidt_SSE.c \
	./Matrix_Ops/invert.c ./Matrix_Ops/lie_mats.c \
	./Matrix_Ops/LU.c ./Matrix_Ops/LU_SSE.c \
	./Matrix_Ops/taylor_logs.c \
	./Matrix_Ops/trace_abc.c ./Matrix_Ops/trace_abc_SSE.c \
	./Matrix_Ops/U_Nops.c ./Matrix_Ops/vandermonde.c

## matrix multiply routines
MMULFILES=./Matrix_Ops/MMUL/MMUL.c ./Matrix_Ops/MMUL/MMUL_SSE.c \
	./Matrix_Ops/MMUL/MMUL_SUNC.c ./Matrix_Ops/MMUL/MMUL_SUNC_SSE.c \
	./Matrix_Ops/MMUL/MMUL_dag.c ./Matrix_Ops/MMUL/MMUL_dag_SSE.c \
	./Matrix_Ops/MMUL/MMUL_dag_SUNC.c \
	./Matrix_Ops/MMUL/MMUL_dag_SUNC_SSE.c \
	./Matrix_Ops/MMUL/MMULdag.c ./Matrix_Ops/MMUL/MMULdag_SSE.c \
	./Matrix_Ops/MMUL/MMULdag_SUNC.c ./Matrix_Ops/MMUL/MMULdag_SUNC_SSE.c \
	./Matrix_Ops/MMUL/MMULdagdag.c ./Matrix_Ops/MMUL/MMULdagdag_SSE.c \
	./Matrix_Ops/MMUL/MMULdagdag_SUNC.c \
	./Matrix_Ops/MMUL/MMULdagdag_SUNC_SSE.c

## smearing codes, generic and targeted blockings, flows etc
SMEARFILES=./Smear/4D_fast.c ./Smear/adaptive_flow.c ./Smear/HYP.c \
	./Smear/ND_generic_HYP.c ./Smear/projectors.c \
	./Smear/smear.c ./Smear/staples.c \
	./Smear/wflow.c ./Smear/wflowfuncs.c

## quenched U1
U1FILES=./U1/SUNCxU1_config.c ./U1/U1_obs.c ./U1/U1_top.c

## utilities: memory checking, rngs splines byte swapping
UTILSFILES=./Utils/GLU_bswap.c ./Utils/GLU_memcheck.c ./Utils/GLU_rng.c \
	./Utils/GLU_splines.c ./Utils/GLU_timer.c ./Utils/KISS_rng.c \
	./Utils/MWC_1038.c ./Utils/MWC_4096.c \
	./Utils/well_rng_19937a.c ./Utils/GLU_malloc.c ./Utils/GLU_sums.c

## wrapper files
WRAPPERSFILES=./Wrappers/CUT_wrap.c ./Wrappers/GF_wrap.c \
	./Wrappers/GLUlib_wrap.c ./Wrappers/OBS_wrap.c ./Wrappers/SM_wrap.c

## all the source files apart from ./Run/Mainfile.c
libGLU_a_SOURCES = \
	${CUTSFILES} \
	${FIELDFILES} \
	${GEOMETRYFILES} \
	${IOFILES} \
	${LANDAUFILES} \
	${MATRIXOPSFILES} \
	${MMULFILES} \
	${SMEARFILES} \
	${U1FILES} \
	${UTILSFILES} \
	${WRAPPERSFILES}

## and their respective headers
include_HEADERS = 
	$(CUTSFILES:./Cuts/%.c=./Headers/%.h) \
	$(FIELDFILES:./Fild/%.c=./Headers/%.h) \
	$(GEOMETRYFILES:./Geometry/%.c=./Headers/%.h) \
	$(IOFILES:./IO/%.c=./Headers/%.h) \
	$(LANDAUFILES:./Landau/%.c=./Headers/%.h) \
	$(MATRIXOPSFILES:./Matrix_Ops/%.c=./Headers/%.h) \
	$(MMULFILES:./Matrix_Ops/MMUL/%.c=./Headers/%.h) \
	$(U1FILES:./U1/%.c=./Headers/%.h) \
	$(UTILSFILES:./Utils/%.c=./Headers/%.h) \
	$(WRAPPERSFILES:./Wrappers/%.c=./Headers/%.h) \
	./Headers/Mainfile.h \
	../config.h

if !PREF

Bindir = "${prefix}"/bin

bin_PROGRAMS = GLU

GLU_SOURCES = ./Run/Mainfile.c
GLU_CFLAGS = ${CFLAGS} -I${TOPDIR}/src/Headers/
GLU_LDADD = libGLU.a ${LDFLAGS}

endif