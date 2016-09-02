#
# Make file for FELIB90
#
# Make sure that FC and FFLAGS are set
# sensibly for your system

#default directories
include Makefile.directories

# Compiler/Linker options
include Makefile.compilers
ifeq  (,$(findstring $(FC),$(COMPILERS)))
# defaults
include Makefile.defaults
else
# selected compiler
include Makefile.$(FC)
endif

default:
	@echo "********************************************"
	@echo "Finite Element Library 90 - FELIB90 Makefile"
	@echo "********************************************"
	@echo "Options are:"
	@echo " "
	@echo "all - all libraries and progs are compiled and linked"
	@echo "globals - global parameters and vaules"
	@echo "utils - utility routines - dynamic memory"
	@echo "machine - machine dependent routines"
	@echo "library - main library"
	@echo "progs - examples programs"
	@echo "clean - all *.o are removed"
	@echo "veryclean - all *.a, *.o and *.mod are removed"
	@echo " "
	@echo "The default compiler is $(FC) with the following options:"
	@echo "	- compiler: $(FFLAGS)" 
	@echo "	- linker  : $(LFLAGS)"
	@echo "which can be found in Makefile.defaults."
	@echo " "
	@echo "You should set up Makefile.compilers but try just 'make all' "
	@echo "and see what happens!"
	@echo " "
	
all: globals utils machine library progs

globals:
	cd $(GLOBALS_DIR);make FC=$(FC) FFLAGS="$(FFLAGS)" LFLAGS="$(LFLAGS)" MODULES="$(MODULES)" MOD_MV="$(MOD_MV)"

utils:
	cd $(UTILS_DIR);make FC=$(FC) FFLAGS="$(FFLAGS)" LFLAGS="$(LFLAGS)" MODULES="$(MODULES)" MOD_MV="$(MOD_MV)"

machine:
	cd $(MACHINE_DIR);make FC=$(FC) FFLAGS="$(FFLAGS)" LFLAGS="$(LFLAGS)"  MODULES="$(MODULES)" MOD_MV="$(MOD_MV)"

library:
	cd $(LIB_DIR);make FC=$(FC) FFLAGS="$(FFLAGS)" LFLAGS="$(LFLAGS)"  MODULES="$(MODULES)" MOD_MV="$(MOD_MV)"

progs:
	cd $(PROG_DIR);make FC=$(FC) FFLAGS="$(FFLAGS)" LFLAGS="$(LFLAGS)"  MODULES="$(MODULES)"

clean:
	cd $(GLOBALS_DIR);make clean
	cd $(UTILS_DIR);make clean
	cd $(MACHINE_DIR);make clean
	cd $(LIB_DIR);make clean
	cd $(MODULES_DIR);make clean
	cd $(PROG_DIR); make clean

veryclean:
	cd $(GLOBALS_DIR);make veryclean
	cd $(UTILS_DIR);make veryclean
	cd $(MACHINE_DIR);make veryclean
	cd $(LIB_DIR);make veryclean
	cd $(MODULES_DIR);make veryclean
	cd $(PROG_DIR); make veryclean

