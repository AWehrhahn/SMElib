AUTOMAKE_OPTIONS = foreign

ACLOCAL_AMFLAGS = -I m4

CXXLD = $(F77)

lib_LTLIBRARIES = libsme.la
libsme_la_SOURCES = \
	src/sme/platform.h \
	src/sme/sme_synth_faster.cpp \
	src/eos/eos.f \
	src/eos/eos_eqns.f \
	src/eos/eos_math_special.f \
	src/sme/hlinop.f \
	src/sme/hlinprof.f
libsme_la_LDFLAGS = -fPIC -version-info 5:22:0 -shared -no-undefined
libsme_la_CPPFLAGS = -fPIC -shared
libsme_la_FFLAGS = -fPIC -shared -fexpensive-optimizations
libsme_la_LIBADD = -lc -lm -lstdc++ -lgfortran

dist_pkgdata_DATA = src/data/bpo_self.grid.INTEL \
	src/data/Fe1_Bautista2017.dat.INTEL \
	src/data/NH_Stancil2018.dat.INTEL \
	src/data/stehle_long.dat.INTEL \
	src/data/vcsbalmer.dat

BUILT_SOURCES = src/sme/platform.h

CLEANFILES = src/sme/platform.h

PLATFORM = $(shell uname -s | tr '[:upper:]' '[:lower:]')
MACHINE = $(shell uname -m)
BITS = 64
SHL_EXT = "$(PLATFORM).$(MACHINE).$(BITS)"

DATA_DIR = @datadir@/libsme/

src/sme/platform.h:
	echo '#define PLATFORM "$(PLATFORM)"' > src/sme/platform.h
	echo '#define DATA_DIR "$(DATA_DIR)"' >> src/sme/platform.h

# Create a symlink to the traditional library name
$(shell mkdir @libdir@; ln -s libsme.so @libdir@/sme_synth.so.$(SHL_EXT))
