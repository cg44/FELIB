
#Machine/compile independent

# Compiler/Linker options
include ../Makefile.compilers

ifeq  (,$(findstring $(FC),$(COMPILERS)))
# defaults
include ../Makefile.defaults
else
# selected compiler
include ../Makefile.$(FC)
endif


ROUTINES = adunit.o erunit.o maxint.o ovflo.o unflo.o veps.o vtol.o

.SUFFIXES: .f90 .o

.f90.o:
	$(FC) $(FFLAGS) $(MODULES) -c $<
	$(MOD_MV)

default:  $(ROUTINES) machine90.o
	ar cr libmachine90.a $(ROUTINES) machine90.o
	
machine90.o: $(ROUTINES)

veryclean:
	-rm *.o *.mod *.a

clean:
	-rm *.o

