[![Build Status](https://travis-ci.com/AWehrhahn/SMElib.svg?branch=master)](https://travis-ci.com/AWehrhahn/SMElib)

# SMElib
Spectroscopy Made Easy Source Library

This is just the C and Fortran part of SME. The complete package is available at [download](https://github.com/AWehrhahn/SME). The classic IDL version of SME is available for [download](http://www.stsci.edu/~valenti/sme.html).

Spectroscopy Made Easy (SME) is a software tool that fits an observed
spectrum of a star with a model spectrum. Since its initial release in
[1996](http://adsabs.harvard.edu/abs/1996A%26AS..118..595V).

## Download
You can find compiled versions of the library for Unix, Mac OS, and Windows under [Releases](https://github.com/AWehrhahn/SMElib/releases).
Note that depending on your system you might have to install libgfortran (version 3) as well.

## Notes
 - SMElib requires libgfortran3 on most OSs to be installed (e.g. using `sudo apt install libgfortran3` on Ubuntu)
 - SMELib needs the datafiles to be present, or it will fail silently. It is therefore recommended to use the included `setLibraryPath(path-to-the-datafiles)` function. While SMElib comes with a default location when it is compiled, the location is dependant on the machine it is run on. You can check the currently set path with `getLibraryPath()` and the names of the required datafiles with `GetDataFiles()`.

## Build
It is also possible to build the library yourself. This requires a C and a Fortran 77 compiler.
SMELib can be build using the GNU Autotools using the following commands
```
# Clone the repository
git clone https://github.com/AWehrhahn/SMElib.git
# Move into the new directory
cd SMElib

# Set up the Autotools
./bootstrap
# This creates a Makefile, that will compile SMElib in the local directory
# If you want to compile it for this machine, remove the '--prefix=$PWD' parameter
./configure --prefix=$PWD
# Compile the library
make install
```
The compiled library should now found in "lib" (or in "bin" on Windows), while the datafiles are in "share/libsme".

Common issues with this compilation are:
  - The compiler can't find libgfortran. Find libgfortran on your machine and set `LDFLAGS="-LFortranPath"`, where FortranPath is the path to the directory containing libgfortran. The path might be located using `gfortran --print-file-name=`.
  - The compiler uses the wrong compilers. Set them exlicitly using `CXX=CCompiler` and `F77=FCompiler`.
