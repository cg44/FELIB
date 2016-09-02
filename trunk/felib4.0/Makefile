#
# Make file for FELIB
#
# Make sure that FC and FFLAGS are set
# sensibly for your system

FC=g77
FFLAGS= 
LFLAGS=

default:
	@echo " FELIB 4.0 - Finite Element Libray"
	@echo " "
	@echo " make all should build everything"
	@echo " make clean to clean up after"
	@echo " make progs will make the programs only"
	@echo " make tests will run some tests"
	@echo " "
	@echo " Current defaults:"
	@echo "	FC	= "$(FC)
	@echo "	FFLAGS	= "$(FFLAGS)
	@echo "	LFLAGS	= "$(LFLAGS)
	@echo
all:
	cd machine; make machine FC=$(FC) FFLAGS="$(FFLAGS)" LFLAGS="$(LFLAGS)"
	cd library; make felib FC=$(FC) FFLAGS="$(FFLAGS)" LFLAGS="$(LFLAGS)"
	cd programs;make programs FC=$(FC) FFLAGS="$(FFLAGS)" LFLAGS="$(LFLAGS)"

tests: 
	cd programs; test-run.sh

progs:
	cd programs; make programs FC=$(FC) FFLAGS=$(FFLAGS) LFLAGS=$(LFLAGS)
clean:
	cd machine; make clean
	cd library; make clean
	cd programs;make clean

