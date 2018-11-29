## simple makefile for the unit tests

TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

bin_PROGRAMS = unit
unit_SOURCES = units.c geom_test.c U_Nops_test.c MMULs_test.c

## compile against GLUlib installed somewhere
unit_CFLAGS = -Wall ${CFLAGS} -I${TOPDIR}/src/Headers/
unit_LDADD = ../src/libGLU.a ${LDFLAGS}
