## simple makefile for the unit tests

TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

bin_PROGRAMS = unit
unit_SOURCES = units.c

## compile against GLUlib installed somewhere
unit_CFLAGS = -Wall -I${TOPDIR}/src/Headers/ ${CFLAGS}
unit_LDADD = -L${TOPDIR}/src/ -lGLU ${LDFLAGS}
