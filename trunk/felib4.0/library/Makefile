#
# Makefile to build the Finite Element Library - Release 4
#
# (C) 1996 Council for the Central Laboratory of the 
#                  Research Councils (CCLRC) 
#          Rutherford Appleton Laboratory
#
# Please address all bugs or problems to author : Dr Chris Greenough
# Tel. Abingdon (01235) 445307
# Fax. Abingdon (01235) 445945
# Email: C.Greenough@rutherford.ac.uk
#

MACHINE=../machine
MACHINE_LIB=..$(MACHINE)/machine-$(FC).a


BINARIES = \
asful.o   brk20.o   csysub.o  fredif.o  matmul.o  qlin2.o   rassym.o  trim6.o \
asfulg.o  brk32.o   cvcnul.o  gaurdn.o  matnul.o  qlin3.o   rasusm.o  update.o \
aslms.o   brk8.o    daxi.o    gausol.o  matran.o  qlval.o   rod2.o    vecadd.o \
asrhs.o   chobak.o  dcsbrk.o  gausub.o  matsub.o  qlvec.o   rod3.o    veccop.o \
assym.o   chofwd.o  dcsqua.o  house.o   matvec.o  qqua4.o   rod4.o    vecmat.o \
assymg.o  chordn.o  dcstri.o  iasrhs.o  mvsyb.o   qqua9.o   scaprd.o  vecnul.o \
asusm.o   chosol.o  direct.o  iassym.o  mvusb.o   qtet4.o   select.o  vecsub.o \
asusmg.o  chosub.o  diso.o    iasusm.o  norm.o    qtri4.o   shapfn.o  vmsyb.o \
b2c2.o    cmtnul.o  dplt.o    jaco.o    prtgeo.o  qtri7.o   sideno.o  vmusb.o \
b2p2.o    cprtmt.o  dpsn.o    linqua.o  prtmat.o  quam12.o  surbrk.o  wdg15.o \
b3c3.o    cprtvc.o  dpss.o    lintri.o  prttop.o  quam4.o   tet10.o   wdg6.o \
bndwth.o  cprtvl.o  dyad.o    matadd.o  prtval.o  quam8.o   tet20.o \
bqbrk.o   csybak.o  elgeom.o  matcop.o  prtvec.o  qwdg6.o   tet4.o \
bqqua.o   csyrdn.o  errmes.o  matidn.o  qbrk6.o   qwdg8.o   trim10.o \
bqtri.o   csysol.o  formnf.o  matinv.o  qbrk8.o   rasrhs.o  trim3.o \
vtol.o

default:
	@echo "*** Do not run this Makefile from here - options will be unknown"

felib:	$(BINARIES) 
	ar cr felib-$(FC).a *.o $(MACHINE)/*.o

clean:
	-rm *.o *.a

