#
# Makefile for PROG90
#

# FELIB90 directory
include ../Makefile.directories

# Compiler/Linker options
include ../Makefile.compilers

ifeq  (,$(findstring $(FC),$(COMPILERS)))
# defaults
include ../Makefile.defaults
else
# selected compiler
include ../Makefile.$(FC)
endif

# Libraries
LIBRARIES=-L../lib -lfelib90

PROGS=seg3p1

.SUFFIXES: .f90 .o

.f90.o:
	$(FC) -c $(FFLAGS) $(MODULES) $<
#$(MOD_MV)

default: $(PROGS)

seg3p1: def3p1.o seg3p1.o 
	$(FC) $(LFLAGS) $(MODULES)  -o seg3p1 seg3p1.o def3p1.o $(OBJS) $(LIBRARIES)

clean:
	-rm *.o 
	-rm $(PROGS)

veryclean:
	-rm *.o 
	-rm *.mod 
	-rm $(PROGS) 
	-rm \#*
