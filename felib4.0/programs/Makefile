#
# Makefile to Run Example Programs from the Finite Element Library
#
# (C) 1988 Science and Engineering Research Council
#          Rutherford Appleton Laboratory
#
# Please address all bugs or problems to author : Chris Greenough
# Tel. Abingdon (0235) 21900 (ext 5307)
#

FELIB = ../library/felib-$(FC).a


ALL = seg1p1 seg1p2 seg2p1 seg2p2 seg2p3 seg3p1 seg3p2\
seg3p3 seg3p4 seg4p1 seg4p2 seg5p1 seg5p2

default: 
	@echo "*** Do not run this Makefile from here - options will be unknown"

programs: $(ALL)

# Segment 1.1
seg1p1 : seg1p1.o
	$(FC) -o seg1p1 $(LFLAGS) seg1p1.o $(FELIB)
# Segment 1.2
seg1p2 : seg1p2.o
	$(FC) -o seg1p2  $(LFLAGS) seg1p2.o $(FELIB)
# Segment 2.1
seg2p1 : seg2p1.o
	$(FC) -o seg2p1  $(LFLAGS) seg2p1.o $(FELIB)
# Segment 2.2
seg2p2 : seg2p2.o
	$(FC) -o seg2p2  $(LFLAGS) seg2p2.o $(FELIB)
# Segment 2.3
seg2p3 : seg2p3.o
	$(FC) -o seg2p3  $(LFLAGS) seg2p3.o $(FELIB)
# Segment 3.1
seg3p1 : seg3p1.o
	$(FC) -o seg3p1  $(LFLAGS) seg3p1.o $(FELIB)
# Segment 3.2
seg3p2 : seg3p2.o
	$(FC) -o seg3p2  $(LFLAGS) seg3p2.o $(FELIB)
# Segment 3.3
seg3p3 : seg3p3.o
	$(FC) -o seg3p3  $(LFLAGS) seg3p3.o $(FELIB)
# Segment 3.4
seg3p4 : seg3p4.o
	$(FC) -o seg3p4  $(LFLAGS) seg3p4.o $(FELIB)
# Segment 4.1
seg4p1 : seg4p1.o
	$(FC) -o seg4p1  $(LFLAGS) seg4p1.o $(FELIB)
# Segment 4.2
seg4p2 : seg4p2.o
	$(FC) -o seg4p2  $(LFLAGS) seg4p2.o $(FELIB)
# Segment 5.1
seg5p1 : seg5p1.o
	$(FC) -o seg5p1  $(LFLAGS) seg5p1.o $(FELIB)
# Segment 5.2
seg5p2 : seg5p2.o
	$(FC) -o seg5p2  $(LFLAGS) seg5p2.o $(FELIB)

clean:
	-rm *.o
	-rm $(ALL)
	-rm \#diff* \#res* res*.tmp fort.*

