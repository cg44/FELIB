#
# Global parameters and values
#

# Compiler/Linker options
include ../Makefile.compilers

ifeq  (,$(findstring $(FC),$(COMPILERS)))
# defaults
include ../Makefile.defaults
else
# selected compiler
include ../Makefile.$(FC)
endif

.SUFFIXES: .f90 .o

.f90.o:
	$(FC) $(FFLAGS) -c $<
	$(MOD_MV)
	
default: felib_globals.o

veryclean:
	-rm *.o *.a *.mod
clean:
	-rm *.o
