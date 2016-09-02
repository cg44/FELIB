#
# The basic routine library
#

#default directories
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


GLOBALS90= ../$(GLOBALS_DIR)/felib_globals.o

LIB90 = errmes.o asful.o asrhs.o assym.o bndwth.o chosol.o \
direct.o elgeom.o  getgeo.o gettop.o \
matinv.o matmul.o matran.o matvec.o prtgeo.o \
prttop.o prtval.o qqua4.o quam4.o scaprd.o matnul.o \
vecnul.o quam8.o vecadd.o matadd.o \
matmul_intrinsic.o

MACHINE90 = adunit.o erunit.o maxint.o ovflo.o unflo.o veps.o vtol.o

UTILS90 = space.o setopt.o


.SUFFIXES: .f90 .o

.f90.o:
	$(FC) $(FFLAGS) $(MODULES) -c $<
	$(MOD_MV)

default: $(LIB90) felib90.o 
	ar cr libfelib90.a $(GLOBALS90) $(LIB90) $(addprefix ../$(MACHINE_DIR)/,$(MACHINE90)) $(addprefix ../$(UTILS_DIR)/,$(UTILS90)) felib90.o
	cp libfelib90.a ../$(LIB)

veryclean:
	-rm *.o *.a *.mod 

clean:
	-rm *.o 
