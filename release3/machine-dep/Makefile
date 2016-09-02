#
# Makefile to build Machine Dependent Routine Library for FELIB Release 3
#
# (C) 1997 Council for the Central Laboratory of the Research Councils
#          Rutherford Appleton Laboratory
#
# Please address all bugs or problems to author : Chris Greenough 
# Tel. Abingdom (01235) 821900 (ext 5307)
# 

FLAGS =  -g

ALL = adunit.o maxint.o unflo.o	veps.o erunit.o ovflo.o

machine-dep.a:	$(ALL)
	ar cr machine-dep.a *.o
	ranlib machine-dep.a
	/bin/cp *.o ../library

adunit.o: adunit.f
	f77 -c $(FLAGS) adunit.f

maxint.o: maxint.f
	f77 -c $(FLAGS) maxint.f

unflo.o: unflo.f
	f77 -c $(FLAGS) unflo.f

veps.o: veps.f
	f77 -c $(FLAGS) veps.f

erunit.o: erunit.f
	f77 -c $(FLAGS) erunit.f

ovflo.o: ovflo.f
	f77 -c $(FLAGS) ovflo.f
