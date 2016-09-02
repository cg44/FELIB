#
# Makefile to build Machine Dependent Routine Library for FELIB Release 3
#
# (C) 1997 Council for the Central Laboratory of the Research Councils
#          Rutherford Appleton Laboratory
#
# Please address all bugs or problems to author : Chris Greenough 
# Tel. Abingdom (01235) 821900 (ext 5307)
# 

#FFLAGS =  -u

ALL = adunit.o maxint.o unflo.o	veps.o erunit.o ovflo.o

machine.a:	$(ALL)
	ar cr machine.a *.o
	ranlib machine.a

clean:
	-rm *.a *.o
