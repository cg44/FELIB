#f95 - NAG Fortran 95

# FELIB90 Root
FELIB90=/home/mathsoft/felib90

# FELIB90 Directories
include $(FELIB90)/Makefile.directories

# f95 options
FFLAGS=  -g -pg -mdir ../$(MODULES_DIR) -strict95
MODULES=-I../$(MODULES_DIR)
LFLAGS= -pg -g
MOD_MV=
