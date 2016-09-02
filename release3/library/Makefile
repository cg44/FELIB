#
# Makefile to build the Finite Element Library - Release 3
#
# (C) 1997 Council for the Central Laboratory of the Research Councils
#          Rutherford Appleton Laboratory
#
# Please address all bugs or problems to author : Chris Greenough
# Tel. Abingdom (01235) 821900 (ext 5307)
#

FLAGS = 

ALL = cprtmt.o gausub.o prttop.o scaprd.o asful.o cprtvc.o house.o prtval.o select.o asfulg.o cprtvl.o iasrhs.o prtvec.o shapfn.o aslms.o csybak.o iassym.o qbrk6.o sideno.o asrhs.o csyrdn.o iasusm.o qbrk8.o surbrk.o assym.o csysol.o jaco.o qlin2.o tet10.o assymg.o csysub.o linqua.o qlin3.o tet20.o asusm.o cvcnul.o lintri.o qlval.o tet4.o asusmg.o daxi.o matadd.o qlvec.o trim10.o b2c2.o dcsbrk.o matcop.o qqua4.o trim3.o b2p2.o dcsqua.o matidn.o qqua9.o trim6.o b3c3.o dcstri.o matinv.o qtet4.o bndwth.o direct.o matmul.o qtri4.o update.o bqbrk.o diso.o matnul.o qtri7.o vecadd.o bqqua.o dplt.o matran.o quam12.o veccop.o bqtri.o dpsn.o matsub.o quam4.o vecmat.o brk20.o dpss.o matvec.o quam8.o vecnul.o brk32.o dyad.o qwdg6.o vecsub.o brk8.o elgeom.o mvsyb.o qwdg8.o chobak.o errmes.o rasrhs.o vmsyb.o chofwd.o mvusb.o rassym.o vmusb.o chordn.o formnf.o norm.o rasusm.o chosol.o fredif.o rod2.o wdg15.o chosub.o gaurdn.o prtgeo.o rod3.o wdg6.o cmtnul.o gausol.o prtmat.o rod4.o vtol.o

MACHINE = adunit.o maxint.o unflo.o	veps.o erunit.o ovflo.o

felib:	$(ALL) 
	ar cr felib.a $(ALL) $(MACHINE)
	ranlib felib.a 

asful.o:    asful.f
	f77 -c $(FLAGS) asful.f
asfulg.o:   asfulg.f
	f77 -c $(FLAGS) asfulg.f
aslms.o:    aslms.f
	f77 -c $(FLAGS) aslms.f
asrhs.o:    asrhs.f
	f77 -c $(FLAGS) asrhs.f
assym.o:    assym.f
	f77 -c $(FLAGS) assym.f
assymg.o:   assymg.f
	f77 -c $(FLAGS) assymg.f
asusm.o:    asusm.f
	f77 -c $(FLAGS) asusm.f
asusmg.o:   asusmg.f
	f77 -c $(FLAGS) asusmg.f
b2c2.o:	    b2c2.f
	f77 -c $(FLAGS) b2c2.f
b2p2.o:	    b2p2.f
	f77 -c $(FLAGS) b2p2.f
b3c3.o:	    b3c3.f
	f77 -c $(FLAGS) b3c3.f
bndwth.o:   bndwth.f
	f77 -c $(FLAGS) bndwth.f
bqbrk.o:    bqbrk.f
	f77 -c $(FLAGS) bqbrk.f
bqqua.o:    bqqua.f
	f77 -c $(FLAGS) bqqua.f
bqtri.o:    bqtri.f
	f77 -c $(FLAGS) bqtri.f
brk20.o:    brk20.f
	f77 -c $(FLAGS) brk20.f
brk32.o:    brk32.f
	f77 -c $(FLAGS) brk32.f
brk8.o:	    brk8.f
	f77 -c $(FLAGS) brk8.f
chobak.o:   chobak.f
	f77 -c $(FLAGS) chobak.f
chofwd.o:   chofwd.f
	f77 -c $(FLAGS) chofwd.f
chordn.o:   chordn.f
	f77 -c $(FLAGS) chordn.f
chosol.o:   chosol.f
	f77 -c $(FLAGS) chosol.f
chosub.o:   chosub.f
	f77 -c $(FLAGS) chosub.f
cmtnul.o:   cmtnul.f
	f77 -c $(FLAGS) cmtnul.f
cprtmt.o:   cprtmt.f
	f77 -c $(FLAGS) cprtmt.f
cprtvc.o:   cprtvc.f
	f77 -c $(FLAGS) cprtvc.f
cprtvl.o:   cprtvl.f
	f77 -c $(FLAGS) cprtvl.f
csybak.o:   csybak.f
	f77 -c $(FLAGS) csybak.f
csyrdn.o:   csyrdn.f
	f77 -c $(FLAGS) csyrdn.f
csysol.o:   csysol.f
	f77 -c $(FLAGS) csysol.f
csysub.o:   csysub.f
	f77 -c $(FLAGS) csysub.f
cvcnul.o:   cvcnul.f
	f77 -c $(FLAGS) cvcnul.f
daxi.o:	    daxi.f
	f77 -c $(FLAGS) daxi.f
dcsbrk.o:   dcsbrk.f
	f77 -c $(FLAGS) dcsbrk.f
dcsqua.o:   dcsqua.f
	f77 -c $(FLAGS) dcsqua.f
dcstri.o:   dcstri.f
	f77 -c $(FLAGS) dcstri.f
direct.o:   direct.f
	f77 -c $(FLAGS) direct.f
diso.o:	    diso.f
	f77 -c $(FLAGS) diso.f
dplt.o:	    dplt.f
	f77 -c $(FLAGS) dplt.f
dpsn.o:	    dpsn.f
	f77 -c $(FLAGS) dpsn.f
dpss.o:	    dpss.f
	f77 -c $(FLAGS) dpss.f
dyad.o:	    dyad.f
	f77 -c $(FLAGS) dyad.f
elgeom.o:   elgeom.f
	f77 -c $(FLAGS) elgeom.f
errmes.o:   errmes.f
	f77 -c $(FLAGS) errmes.f
formnf.o:   formnf.f
	f77 -c $(FLAGS) formnf.f
fredif.o:   fredif.f
	f77 -c $(FLAGS) fredif.f
gaurdn.o:   gaurdn.f
	f77 -c $(FLAGS) gaurdn.f
gausol.o:   gausol.f
	f77 -c $(FLAGS) gausol.f
gausub.o:   gausub.f
	f77 -c $(FLAGS) gausub.f
house.o:    house.f
	f77 -c $(FLAGS) house.f
iasrhs.o:   iasrhs.f
	f77 -c $(FLAGS) iasrhs.f
iassym.o:   iassym.f
	f77 -c $(FLAGS) iassym.f
iasusm.o:   iasusm.f
	f77 -c $(FLAGS) iasusm.f
jaco.o:	    jaco.f
	f77 -c $(FLAGS) jaco.f
linqua.o:   linqua.f
	f77 -c $(FLAGS) linqua.f
lintri.o:   lintri.f
	f77 -c $(FLAGS) lintri.f
matadd.o:   matadd.f
	f77 -c $(FLAGS) matadd.f
matcop.o:   matcop.f
	f77 -c $(FLAGS) matcop.f
matidn.o:   matidn.f
	f77 -c $(FLAGS) matidn.f
matinv.o:   matinv.f
	f77 -c $(FLAGS) matinv.f
matmul.o:   matmul.f
	f77 -c $(FLAGS) matmul.f
matnul.o:   matnul.f
	f77 -c $(FLAGS) matnul.f
matran.o:   matran.f
	f77 -c $(FLAGS) matran.f
matsub.o:   matsub.f
	f77 -c $(FLAGS) matsub.f
matvec.o:   matvec.f
	f77 -c $(FLAGS) matvec.f
mvsyb.o:    mvsyb.f
	f77 -c $(FLAGS) mvsyb.f
mvsysp.o:   mvsysp.f
	f77 -c $(FLAGS) mvsysp.f
mvusb.o:    mvusb.f
	f77 -c $(FLAGS) mvusb.f
norm.o:	    norm.f
	f77 -c $(FLAGS) norm.f
prtgeo.o:   prtgeo.f
	f77 -c $(FLAGS) prtgeo.f
prtmat.o:   prtmat.f
	f77 -c $(FLAGS) prtmat.f
prttop.o:   prttop.f
	f77 -c $(FLAGS) prttop.f
prtval.o:   prtval.f
	f77 -c $(FLAGS) prtval.f
prtvec.o:   prtvec.f
	f77 -c $(FLAGS) prtvec.f
qbrk6.o:    qbrk6.f
	f77 -c $(FLAGS) qbrk6.f
qbrk8.o:    qbrk8.f
	f77 -c $(FLAGS) qbrk8.f
qlin2.o:    qlin2.f
	f77 -c $(FLAGS) qlin2.f
qlin3.o:    qlin3.f
	f77 -c $(FLAGS) qlin3.f
qlval.o:    qlval.f
	f77 -c $(FLAGS) qlval.f
qlvec.o:    qlvec.f
	f77 -c $(FLAGS) qlvec.f
qqua4.o:    qqua4.f
	f77 -c $(FLAGS) qqua4.f
qqua9.o:    qqua9.f
	f77 -c $(FLAGS) qqua9.f
qtet4.o:    qtet4.f
	f77 -c $(FLAGS) qtet4.f
qtri4.o:    qtri4.f
	f77 -c $(FLAGS) qtri4.f
qtri7.o:    qtri7.f
	f77 -c $(FLAGS) qtri7.f
quam12.o:   quam12.f
	f77 -c $(FLAGS) quam12.f
quam4.o:    quam4.f
	f77 -c $(FLAGS) quam4.f
quam8.o:    quam8.f
	f77 -c $(FLAGS) quam8.f
qwdg6.o:    qwdg6.f
	f77 -c $(FLAGS) qwdg6.f
qwdg8.o:    qwdg8.f
	f77 -c $(FLAGS) qwdg8.f
rasrhs.o:   rasrhs.f
	f77 -c $(FLAGS) rasrhs.f
rassym.o:   rassym.f
	f77 -c $(FLAGS) rassym.f
rasusm.o:   rasusm.f
	f77 -c $(FLAGS) rasusm.f
rod2.o:	    rod2.f
	f77 -c $(FLAGS) rod2.f
rod3.o:	    rod3.f
	f77 -c $(FLAGS) rod3.f
rod4.o:	    rod4.f
	f77 -c $(FLAGS) rod4.f
scaprd.o:   scaprd.f
	f77 -c $(FLAGS) scaprd.f
select.o:   select.f
	f77 -c $(FLAGS) select.f
shapfn.o:   shapfn.f
	f77 -c $(FLAGS) shapfn.f
sideno.o:   sideno.f
	f77 -c $(FLAGS) sideno.f
surbrk.o:   surbrk.f
	f77 -c $(FLAGS) surbrk.f
tet10.o:    tet10.f
	f77 -c $(FLAGS) tet10.f
tet20.o:    tet20.f
	f77 -c $(FLAGS) tet20.f
tet4.o:	    tet4.f
	f77 -c $(FLAGS) tet4.f
trim10.o:   trim10.f
	f77 -c $(FLAGS) trim10.f
trim3.o:    trim3.f
	f77 -c $(FLAGS) trim3.f
trim6.o:    trim6.f
	f77 -c $(FLAGS) trim6.f
update.o:   update.f
	f77 -c $(FLAGS) update.f
vecadd.o:   vecadd.f
	f77 -c $(FLAGS) vecadd.f
veccop.o:   veccop.f
	f77 -c $(FLAGS) veccop.f
vecmat.o:   vecmat.f
	f77 -c $(FLAGS) vecmat.f
vecnul.o:   vecnul.f
	f77 -c $(FLAGS) vecnul.f
vecsub.o:   vecsub.f
	f77 -c $(FLAGS) vecsub.f
vtol.o:	vtol.f
	f77 -c $(FLAGS) vtol.f
vmsyb.o:    vmsyb.f
	f77 -c $(FLAGS) vmsyb.f
vmusb.o:    vmusb.f
	f77 -c $(FLAGS) vmusb.f
wdg15.o:    wdg15.f
	f77 -c $(FLAGS) wdg15.f
wdg6.o:	    wdg6.f
	f77 -c $(FLAGS) wdg6.f
