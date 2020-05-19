# Metadata
AC_INIT([libsme], [5.22], [nikolai.piskunov@physics.uu.se])
AC_CANONICAL_TARGET

# Use Automake
AM_INIT_AUTOMAKE([subdir-objects])
LT_INIT
AC_CONFIG_MACRO_DIRS([m4])

case $target_os in
    linux-gnu)
    target=linux 
    extension=.so
    outdir=libdir
    ;;
    darwin)
    target=osx
    extension=.dylib
    outdir=libdir
    ;;
    win32|mingw32)
    target=win32
    extension=-5.dll
    outdir=bindir
    ;;
    *) 
    target=$target_os
    extension=.so
    outdir=libdir
    ;;
esac

AC_SUBST(TARGET_CPU, $target_cpu)
AC_SUBST(TARGET, $target)
AC_SUBST(EXTENSION, $extension)
AC_SUBST(OUTPUT_DIR, $outdir)

# Find Dependecies
AC_PROG_CXX
AC_PROG_F77



# Create Makefile from Makefile.in
AC_CONFIG_FILES([Makefile])

# Create Output files
AC_OUTPUT
