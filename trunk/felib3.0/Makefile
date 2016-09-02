#
# Make file for FELIB
#
# Make sure that FC and FFLAGS are set
# sensibly for your system
# The Makefile will generate a felib-$(FC).a for
# your particular compiler.
#

FC=f77
FFLAGS=
LFLAGS=

default:
	@echo " FELIB 3.0 - Finite Element Libray"
	@echo " "
	@echo " make all should build everything"
	@echo " make clean to clean up after"
	@echo " "
	@echo " Current defaults:"
	@echo "	FC	= "$(FC)
	@echo "	FFLAGS	= "$(FFLAGS)
	@echo "	LFLAGS	= "$(LFLAGS)
	@echo " "

all:
	cd machine;make FC=$(FC)
	cd library;make FC=$(FC)
	cd programs;make FC=$(FC)

clean:
	cd machine;make clean
	cd library;make clean
	cd programs;make clean

