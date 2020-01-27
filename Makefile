#
#  Platform independant makefile for the SME Call_External's.
#
SHELL=/bin/sh
CFLAGS=
FFLAGS=
FF1=
FF2=
F_LD_POST=
F_LD_BLAS=
SME_DIR = $(PWD)/src/sme
EOS_DIR = $(PWD)/src/eos
TARGET_DIR = $(PWD)/lib

# The following is the default entry point. This section will determine 
# what system we are on, set the correct flags and call this same makefile
# again with the correct flags.

all : 
	@echo "OS type detected: "`uname`
	PLATFORM=`uname`
	@case `uname` in \
           "SunOS32") if [  -d /proc ]; then \
                        make libs \
                        "SHL_EXT=so.SUN32" \
                        "CC=CC" \
                        "C_LD=cc" \
                        "F_LD=f77" \
                        "F77=f77" \
                        "CFLAGS=-Kpic -DSPARC -G -O3" \
                        "FFLAGS=-pic -G -O3 -ext_names=underscores" \
                        "FF1=-pic -G -O3 -ext_names=underscores" \
                        "FF2=-pic -G -O3 -ext_names=underscores" \
                        "C_LD_FLAGS=-G -DSOLARIS" \
                        "F_LD_FLAGS=-G" \
                        "F_LD_POST= -lsunperf -lsunmath -lF77 -lM77 -lm -lc"; \
                    fi  \
                ;; \
           "SunOS") if [  -d /proc ]; then \
                        make libs \
                        "SHL_EXT=so.SUN64" \
                        "CC=CC" \
                        "C_LD=cc" \
                        "F_LD=f77" \
                        "F77=f77" \
                        "CFLAGS=-xtarget=ultra -xarch=v9 -Kpic -G -O3" \
                        "FFLAGS=-xtarget=ultra -xarch=v9 -pic -G -O3 -ext_names=underscores" \
                        "FF1=-xtarget=ultra -xarch=v9 -pic -G -O3 -ext_names=underscores" \
                        "FF2=-xtarget=ultra -xarch=v9 -pic -G -O3 -ext_names=underscores" \
                        "F_LD_FLAGS=-xtarget=ultra -xarch=v9 -G" \
                        "F_LD_POST= -lsunperf -lsunmath -lF77 -lM77 -lm -lc"; \
                    fi  \
                ;; \
	   "AIX") make libs \
			"SHL_EXT=a.AIX" \
			"CC=cc" \
			"F77=f77" \
			"C_LD=cc" \
			"F_LD=f77" \
			"C_LD_FLAGS=-bM:SRE -bnoentry" \
			"F_LD_FLAGS=-bM:SRE -bnoentry" \
			"SME_EXPORT=-bE:sme_synth.exp.AIX" \
			"F_LD_POST= -lxlf -lxlf90";; \
	   "HP-UX32") make libs \
			"SHL_EXT=so.HPUX32" \
			"CC=aCC" \
			"F77=f90" \
			"C_LD=ld" \
			"F_LD=ld" \
			"CFLAGS= +z +O2 +DD32" \
			"FFLAGS= +z +O2 +DD32 +ppu" \
			"FF1= +z +O2 +DD32 +ppu" \
			"FF2= +z +O2 +DD32 +ppu" \
			"C_LD_FLAGS=-b" \
			"F_LD_FLAGS=-b" \
			"F_LD_POST= /opt/fortran90/lib/libF90.a -lcl -lm" \
                        "F_LD_BLAS= /opt/mlib/lib/pa2.0/liblapack.a" ;; \
	   "HP-UX") make libs \
			"SHL_EXT=so.HPUX64" \
			"CC=aCC" \
			"F77=f90" \
			"C_LD=ld" \
			"F_LD=ld" \
			"CFLAGS= +Z +O2 +DA2.0W +DD64" \
			"FFLAGS= +Z +O2 +DA2.0W +DD64" \
			"FF1= +Z +O2 +DA2.0W +DD64" \
			"FF2= +Z +O2 +DA2.0W +DD64" \
			"C_LD_FLAGS=-b" \
			"F_LD_FLAGS=-b" \
			"F_LD_POST= /opt/fortran90/lib/pa20_64/libF90.a -lcl -lm" \
                        "F_LD_BLAS= /opt/mlib/lib/pa20_64/liblapack.a" ;; \
           "Linux 32g" ) make libs \
                        "SHL_EXT=so.linux.x86.32g" \
                        "CC=g++ -m32" \
                        "F77=gfortran -m32" \
                        "C_LD=g++ -m32" \
                        "F_LD=gfortran -m32" \
                        "CFLAGS= -O3 -fPIC -shared -static" \
                        "FFLAGS=-O -fPIC -shared -static -fexpensive-optimizations" \
                        "FF1=-O0 -fPIC -shared -static -static -fexpensive-optimizations"  \
                        "FF2=-O3 -fPIC -shared -static -fexpensive-optimizations"   \
                        "C_LD_FLAGS=-O -fPIC -shared -static-libgcc" \
                        "F_LD_FLAGS=-O3 -fPIC -shared -fexpensive-optimizations -static-libgcc"  \
                        "F_LD_POST= " \
                        "PLATFORM=Linux_32_Gnu";; \
	   "Linux" ) make libs \
                        "SHL_EXT=so.linux.x86_64.64" \
                        "CC=g++ -m64" \
                        "F77=gfortran -m64" \
                        "C_LD=g++" \
                        "F_LD=gfortran -m64" \
                        "CFLAGS= -fPIC -shared -static" \
                        "FFLAGS= -fPIC -shared -static -fexpensive-optimizations" \
                        "FF1= -fPIC -shared -static -fexpensive-optimizations" \
                        "FF2= -fPIC -shared -static -fexpensive-optimizations" \
                        "C_LD_FLAGS= -fPIC -shared" \
                        "F_LD_FLAGS= -fPIC -shared -fexpensive-optimizations \
                                           -static-libgcc" \
                        "F_LD_POST= -lc -lm -lstdc++" \
                        "CPUPROFILE=cpu_profile " \
                        "HEAPPROFILE=heap_profile " \
                        "PLATFORM=Linux_64_Gnu";; \
           "Linux 32i" ) make libs \
                        "SHL_EXT=so.linux.x86.32i" \
                        "CC=g++" \
                        "F77=ifort" \
                        "C_LD=g++" \
                        "F_LD=g++" \
                        "CFLAGS=-m32 -O2 -fPIC -shared -static " \
                        "FFLAGS=-m32 -O3 -fPIC -shared -static" \
                        "FF1= -m32 -O2 -fPIC -shared -static" \
                        "FF2= -m32 -O2 -fPIC -shared -static" \
                        "F_LD_FLAGS= -m32 -O2 -shared -static-intel" \
                        "F_LD_BLAS= " \
                        "F_LD_POST= -Bstatic" \
                        "PLATFORM=Linux_32_Intel";; \
           "Linux 64i" ) make libs \
                        "SHL_EXT=so.linux.x86_64.64i" \
                        "CC=icc" \
                        "F77=ifort" \
                        "C_LD=icc" \
                        "F_LD=ifort" \
                        "CFLAGS=-m64 -x c++ -O3 -shared -static -fPIC -fno-common -D_REENTRANT" \
                        "FFLAGS=-m64 -O3 -fPIC -shared -static" \
                        "FF1= -m64 -O -fPIC -shared -static" \
                        "FF2= -m64 -O -fPIC -shared -static" \
                        "F_LD_FLAGS= -m64 -O -shared -static" \
                        "F_LD_BLAS= " \
                        "F_LD_POST= -lm -lc -lifport -lstdc++" \
                        "PLATFORM=Linux_64_Intel";; \
           "Darwin 32g" ) make libs \
                        "SHL_EXT=so.darwin.i386.32g" \
                        "CC=g++" \
                        "F77=gfortran -m32" \
                        "C_LD=g++ -m32" \
                        "F_LD=gfortran -m32" \
                        "CFLAGS=-x c++ -O3 -no-cpp-precomp -dynamic -fPIC -fno-common -D_REENTRANT" \
                        "FFLAGS= -O3 -f77rtl -dynamiclib" \
                        "FF1= -O3 -f77rtl -fPIC -shared -static" \
                        "FF2= -O3 -f77rtl -fPIC -shared -static" \
                        "C_LD_FLAGS= -O -fPIC -shared" \
                        "F_LD_FLAGS= -O -fPIC -shared" \
                        "F_LD_BLAS= -llapack -lblas" \
                        "F_LD_POST= -lm -lc -lgcc -lstdc++" \
                        "PLATFORM=Darwin_32_Gnu";; \
           "Darwin 64g" ) make libs \
                        "SHL_EXT=so.darwin.x86_64.64g" \
                       "CC=g++" \
                        "F77=gfortran -m64" \
                        "C_LD=g++ -m64" \
                        "F_LD=gfortran -m64" \
                        "CFLAGS=-x c++ -O3 -no-cpp-precomp -dynamic -fPIC -fno-common -D_REENTRANT" \
                        "FFLAGS= -O3 -dynamiclib" \
                        "FF1= -O3 -fPIC -shared -static" \
                        "FF2= -O3 -fPIC -shared -static" \
                        "F_LD_FLAGS= -O -fPIC -shared -static -static-libgcc" \
                        "F_LD_POST= /usr/local/gfortran/lib/libgfortran.a \
                                    /usr/local/gfortran/lib/libquadmath.a -lm -lc -lgcc -lstdc++" \
                        "PLATFORM=Darwin_64_Gnu";; \
           "Darwin 32i" ) make libs \
                        "SHL_EXT=so.darwin.i386.32" \
                        "CC=icc -m32" \
                        "F77=ifort -m32" \
                        "C_LD=icc -m32" \
                        "F_LD=ifort -m32" \
                        "CFLAGS= -x c++ -O3 -dynamic -fPIC -fno-common -D_REENTRANT" \
                        "FFLAGS= -O3 -f77rtl -dynamiclib" \
                        "FF1= -O3 -f77rtl -fpe0 -fPIC -dynamiclib" \
                        "FF2= -O3 -f77rtl -fpe0 -fPIC -dynamiclib" \
                        "F_LD_FLAGS= -save -O -f77rtl -fpe0 -dynamiclib" \
                        "F_LD_POST= -lm -lc -lstdc++" \
                        "PLATFORM=Darwin_32_Intel";; \
           "Darwin" ) make libs \
                        "SHL_EXT=so.darwin.x86_64.64" \
                        "CC=icc" \
                        "F77=ifort" \
                        "C_LD=icc" \
                        "F_LD=ifort" \
                        "CFLAGS=-m64 -O3 -x c++ -dynamiclib -fPIC -fno-common -D_REENTRANT" \
                        "FFLAGS=-m64 -O3 -f77rtl -fpe0 -dynamiclib" \
                        "FF1= -m64 -O3 -f77rtl -fpe0 -fPIC -dynamiclib" \
                        "FF2= -m64 -O3 -f77rtl -fpe0 -fPIC -dynamiclib" \
                        "F_LD_FLAGS= -save -O -f77rtl -fpe0 -m64 -dynamiclib -static-libstdc++ -static-intel " \
                        "F_LD_POST= -lm -lc -lstdc++" \
                        "PLATFORM=Darwin_64_Intel";; \
 	   *) echo "This system is not supported" ;; \
       esac

# C Only libs

all : sme_synth_faster.o eos.o eos_eqns.o eos_math.o hlinop.o \
                  hlinprof.o smelib.so

install : all

smelib.so:
	$(shell mkdir $(TARGET_DIR))
	$(F_LD) $(F_LD_FLAGS) -o $(TARGET_DIR)/smelib.so $(SME_DIR)/sme_synth_faster.o \
                $(EOS_DIR)/eos.o $(EOS_DIR)/eos_eqns.o $(SME_DIR)/hlinop.o $(SME_DIR)/hlinprof.o $(EOS_DIR)/eos_math.o \
                $(F_LD_POST) $(F_LD_BLAS)
	$(shell ln -s smelib.so $(TARGET_DIR)/sme_synth.so.${SHL_EXT})

sme_synth_faster.o:
	echo '#define PLATFORM "$(PLATFORM)"' > $(SME_DIR)/platform.h
	echo '#define DATA_DIR "$(TARGET_DIR)"' >> $(SME_DIR)/platform.h
	$(CC) $(CFLAGS) -o $(SME_DIR)/sme_synth_faster.o -c $(SME_DIR)/sme_synth_faster.cpp
        
eos.o:
	$(F77) -c $(FFLAGS) -o $(EOS_DIR)/eos.o $(EOS_DIR)/eos.f

eos_eqns.o:
	$(F77) -c $(FFLAGS) -o $(EOS_DIR)/eos_eqns.o $(EOS_DIR)/eos_eqns.f

eos_math.o:
	$(F77) $(FF1) -o $(EOS_DIR)/eos_math.o -c $(EOS_DIR)/eos_math_special.f

hlinop.o:
	$(F77) $(FF2) -o $(SME_DIR)/hlinop.o -c $(SME_DIR)/hlinop.f

hlinprof.o:
	$(F77) $(FF2) -o $(SME_DIR)/hlinprof.o -c $(SME_DIR)/hlinprof.f

# Cleanup

tidy :
	rm -f $(EOS_DIR)/*.o $(SME_DIR)/*.o *.heap cpu_profile* $(TARGET_DIR)/

clean :
	rm -f *.heap cpu_profile*
	rm -f $(SME_DIR)/*.o $(SME_DIR)/*.so $(SME_DIR)/*.sl $(SME_DIR)/*.a $(SME_DIR)/platform.h
	rm -f $(EOS_DIR)/*.o $(EOS_DIR)/*.so $(EOS_DIR)/*.sl $(EOS_DIR)/*.a
