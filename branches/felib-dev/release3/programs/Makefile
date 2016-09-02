#
# Makefile to Run Example Programs from the Finite Element Library
#
# (C) 1997 Council for the Central Laboratory of the Research Councils
#          Rutherford Appleton Laboratory
#
# Please address all bugs or problems to author : Chris Greenough
# Tel. Abingdon (01235) 821900 (ext 5307)
#

#FFLAGS = -g -f68881
FELIB = ../library/felib.a 

ALL = seg1p1 seg1p2 seg2p1 seg2p2 seg2p3 seg3p1 seg3p2\
seg3p3 seg3p4 seg4p1 seg4p2 seg5p1 seg5p2

Temp : $(ALL)

# Segment 1.1
seg1p1 : seg1p1.o
	f77 -o seg1p1 $(FFLAGS) seg1p1.o $(FELIB)
# Segment 1.2
seg1p2 : seg1p2.o
	f77 -o seg1p2  $(FFLAGS) seg1p2.o $(FELIB)
# Segment 2.1
seg2p1 : seg2p1.o
	f77 -o seg2p1  seg2p1.o $(FELIB)
# Segment 2.2
seg2p2 : seg2p2.o
	f77 -o seg2p2  $(FFLAGS) seg2p2.o $(FELIB)
# Segment 2.3
seg2p3 : seg2p3.o
	f77 -o seg2p3  $(FFLAGS) seg2p3.o $(FELIB)
# Segment 3.1
seg3p1 : seg3p1.o
	f77 -o seg3p1  $(FFLAGS) seg3p1.o $(FELIB)
# Segment 3.2
seg3p2 : seg3p2.o
	f77 -o seg3p2  $(FFLAGS) seg3p2.o $(FELIB)
# Segment 3.3
seg3p3 : seg3p3.o
	f77 -o seg3p3  $(FFLAGS) seg3p3.o $(FELIB)
# Segment 3.4
seg3p4 : seg3p4.o
	f77 -o seg3p4  $(FFLAGS) seg3p4.o $(FELIB)
# Segment 4.1
seg4p1 : seg4p1.o
	f77 -o seg4p1  $(FFLAGS) seg4p1.o $(FELIB)
# Segment 4.2
seg4p2 : seg4p2.o
	f77 -o seg4p2  $(FFLAGS) seg4p2.o $(FELIB)
# Segment 5.1
seg5p1 : seg5p1.o
	f77 -o seg5p1  $(FFLAGS) seg5p1.o $(FELIB)
# Segment 5.2
seg5p2 : seg5p2.o
	f77 -o seg5p2  $(FFLAGS) seg5p2.o $(FELIB)

