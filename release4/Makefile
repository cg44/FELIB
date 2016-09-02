#
# Make file for FELIB
#
# Make sure that FC and FFLAGS are set
# sensibly for your system

FC=g77
FFLAGS=
LFLAGS=

all:
	cd machine;make FC=$(FC)
	cd library;make FC=$(FC)
	cd programs;make FC=$(FC)

clean:
	cd machine;make clean
	cd library;make clean
	cd programs;make clean

