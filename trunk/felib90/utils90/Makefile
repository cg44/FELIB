
# Makefile for dynamic memory and utilities

# Compiler/Linker options
include ../Makefile.compilers

ifeq  (,$(findstring $(FC),$(COMPILERS)))
# defaults
include ../Makefile.defaults
else
# selected compiler
include ../Makefile.$(FC)
endif


UTILS90 = space.o setopt.o

.SUFFIXES: .f90 .o

.f90.o:
	$(FC) $(FFLAGS) $(MODULES) -c $<
	$(MOD_MV)
	
default: $(UTILS90) felib_utils.o
	ar cr libutils90.a $(UTILS90)

clean:
	-rm *.o

veryclean:
	-rm *.o *.mod *.a
