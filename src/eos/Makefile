F77 = gfortran
F77_FLAGS = -fPIC
LNK_FLAGS = -fPIC



eos.o: eos.f DEFAULT.EOS SIZES.EOS Makefile
	$(F77) $(F77_FLAGS) -c eos.f

eos_eqns.o: eos_eqns.f SIZES.EOS Makefile
	$(F77) $(F77_FLAGS) -c eos_eqns.f

eos_math_special.o: eos_math_special.f Makefile
	$(F77) $(F77_FLAGS) -c eos_math_special.f

clean:
	rm *.o
