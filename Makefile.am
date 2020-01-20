AUTOMAKE_OPTIONS = foreign

ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES = libsme.la
libsme_la_SOURCES = src/sme/hlinop.f \
	src/sme/hlinprof.f \
	src/eos/eos_eqns.f \
	src/eos/eos_math_special.f \
	src/eos/eos_math.f \
	src/eos/eos.f \
	src/sme/sme_synth_faster.cpp \
	src/sme/platform.h
libsme_la_LDFLAGS = -version-info 5:22:0 -shared -static -fexpensive-optimizations

dist_pkgdata_DATA = src/data/bpo_self.grid.INTEL \
	src/data/Fe1_Bautista2017.dat.INTEL \
	src/data/NH_Stancil2018.dat.INTEL \
	src/data/stehle_long.dat.INTEL \
	src/data/vcsbalmer.dat

# EXTRA_DIST = src/data/bpo_self.grid.INTEL \
# 	src/data/Fe1_Bautista2017.dat.INTEL \
# 	src/data/NH_Stancil2018.dat.INTEL \
# 	src/data/stehle_long.dat.INTEL \
# 	vcsbalmer.dat

BUILT_SOURCES = src/sme/platform.h

CLEANFILES = src/sme/platform.h

PLATFORM = $(shell uname -s)
DATA_DIR = @datadir@/smelib/

src/sme/platform.h:
	echo '#define PLATFORM "$(PLATFORM)"' > src/sme/platform.h
	echo '#define DATA_DIR "$(DATA_DIR)"' >> src/sme/platform.h

