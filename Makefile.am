AUTOMAKE_OPTIONS = foreign

ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES = libsme.la
libsme_la_SOURCES = src/sme/hlinop.f \
    src/sme/hlinprof.f \
    src/eos/eos_eqns.f \
    src/eos/eos_math_special.f \
    src/eos/eos_math.f \
    src/eos/eos.f \
    src/sme/sme_synth_faster.cpp
libsme_la_LDFLAGS = -version-info 0:0:0