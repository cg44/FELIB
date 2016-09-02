#
# Makefile to build the Finite Element Library - Release 3
#
# (C) 1997 Council for the Central Laboratory of the Research Councils
#          Rutherford Appleton Laboratory
#
# Please address all bugs or problems to author : Chris Greenough
# Tel. Abingdom (01235) 821900 (ext 5307)
#


ALL = cprtmt.o gausub.o prttop.o scaprd.o asful.o cprtvc.o house.o prtval.o select.o asfulg.o cprtvl.o iasrhs.o prtvec.o shapfn.o aslms.o csybak.o iassym.o qbrk6.o sideno.o asrhs.o csyrdn.o iasusm.o qbrk8.o surbrk.o assym.o csysol.o jaco.o qlin2.o tet10.o assymg.o csysub.o linqua.o qlin3.o tet20.o asusm.o cvcnul.o lintri.o qlval.o tet4.o asusmg.o daxi.o matadd.o qlvec.o trim10.o b2c2.o dcsbrk.o matcop.o qqua4.o trim3.o b2p2.o dcsqua.o matidn.o qqua9.o trim6.o b3c3.o dcstri.o matinv.o qtet4.o bndwth.o direct.o matmul.o qtri4.o update.o bqbrk.o diso.o matnul.o qtri7.o vecadd.o bqqua.o dplt.o matran.o quam12.o veccop.o bqtri.o dpsn.o matsub.o quam4.o vecmat.o brk20.o dpss.o matvec.o quam8.o vecnul.o brk32.o dyad.o qwdg6.o vecsub.o brk8.o elgeom.o mvsyb.o qwdg8.o chobak.o errmes.o rasrhs.o vmsyb.o chofwd.o mvusb.o rassym.o vmusb.o chordn.o formnf.o norm.o rasusm.o chosol.o fredif.o rod2.o wdg15.o chosub.o gaurdn.o prtgeo.o rod3.o wdg6.o cmtnul.o gausol.o prtmat.o rod4.o vtol.o

MACHINE=../machine/machine.a

felib-$(FC).a: $(ALL) $(MACHINE)
	ar cr felib-$(FC).a $(ALL) ../machine/*.o
	ranlib felib-$(FC).a 

$(MACHINE):
	cd ../machine;make FC=$(FC)

clean:
	-rm *.o *.a
