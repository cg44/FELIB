C
C
C     Copyright (C) 1997 : clrc, Rutherford Appleton Laboratory
C     Chilton, Didcot, Oxfordshire OX11 0QX
C
C
C
      INTEGER BNDNOD,BNODE,DIMEN,DOFEL,DOFNOD,ELNUM,ELTOP,ELTYP,FRENOD,
     *        HBAND,I,IABSS,IAMOVE,IBNODE,ICOORD,IDTPD,IELK,IELTOP,
     *        IFRNOD,IFUN,IGDER,IGDERT,IGEOM,IJAC,IJACIN,ILDER,INF,IP,
     *        IPD,IQUAD,IRHS,ISPNOD,ISTEER,ISYSK,ITER,ITEST,IWGHT,J,
     *        JABSS,JCOORD,JDTPD,JELK,JELTOP,JGDER,JGDERT,JGEOM,JJAC,
     *        JJACIN,JLDER,JNF,JP,JPD,JSYSK,MAXIT,NELE,NF,NIN,NODE,
     *        NODEL,NODNUM,NOUT,NQP,NSIDE,NUMFRE,NUMSEP,SEPBND,SEPNOD,
     *        STEER,TOTDOF,TOTELS,TOTNOD
      DOUBLE PRECISION ABSS,AMOVE,BVAL,CFACTR,COORD,DET,DTPD,ELK,ETA,
     *                 FSDIFF,FUN,GDER,GDERT,GEOM,JAC,JACIN,LDER,P,PD,
     *                 PI,PREVHT,QUOT,R,RHS,SCALE,SYSK,WGHT,XI
      DIMENSION ABSS(3,9),AMOVE(3),DTPD(8,8),ELK(8,8),FUN(8),GDER(3,8),
     *          GDERT(8,3),GEOM(8,3),JAC(3,3),JACIN(3,3),LDER(3,8),
     *          P(3,3),PD(3,8),STEER(8),WGHT(9)
C
C     Problem size dependent arrays
C
      DIMENSION BNODE(50),BVAL(50),COORD(250,3),ELTOP(350,10),
     *          FRENOD(50),NF(250,1),PREVHT(50),RHS(250),SEPBND(50),
     *          SEPNOD(50),SYSK(250,25)
C
      DATA IABSS/3/,IAMOVE/3/,IDTPD/8/,IELK/8/,IFUN/8/,IGDER/3/,
     *     IGDERT/8/,IGEOM/8/,IJAC/3/,IJACIN/3/,ILDER/3/,IP/3/,IPD/3/,
     *     ISTEER/8/,IWGHT/9/,JABSS/9/,JCOORD/3/,JDTPD/8/,JELK/8/,
     *     JGDER/8/,JGDERT/3/,JGEOM/3/,JJAC/3/,JJACIN/3/,JLDER/8/,
     *     JNF/1/,JP/3/,JPD/8/,SCALE/1.0D+10/
C
C     Problem size dependent data statements
C
      DATA IBNODE/50/,ICOORD/250/,IELTOP/350/,IFRNOD/50/,INF/250/,
     *     IRHS/250/,ISPNOD/50/,ISYSK/250/,JELTOP/10/,JSYSK/25/,NSIDE/2/
C
      DATA DOFNOD/1/,NIN/5/,NOUT/6/
C
C     Set ITEST for full checking
C
      ITEST = 0
C
      PI = ATAN(1.0D0)*4.0D0
C
C     *                           **********************
C     *                           *                    *
C     *                           * Input Data Section *
C     *                           *                    *
C     *                           **********************
C
C     Input of nodal geometry
C
      WRITE (NOUT,FMT=9950)
      READ (NIN,FMT=9990) TOTNOD,DIMEN
      WRITE (NOUT,FMT=9940) TOTNOD,DIMEN
C
      DO 1000 I = 1,TOTNOD
          READ (NIN,FMT=9980) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
          WRITE (NOUT,FMT=9930) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
 1000 CONTINUE
C
C     Input of element topology
C
      WRITE (NOUT,FMT=9920)
      READ (NIN,FMT=9990) TOTELS
      WRITE (NOUT,FMT=9940) TOTELS
C
      DO 1010 I = 1,TOTELS
          READ (NIN,FMT=9990) ELNUM,ELTYP,NODEL,
     *      (ELTOP(ELNUM,J+2),J=1,NODEL)
          WRITE (NOUT,FMT=9940) ELNUM,ELTYP,NODEL,
     *      (ELTOP(ELNUM,J+2),J=1,NODEL)
          ELTOP(ELNUM,1) = ELTYP
          ELTOP(ELNUM,2) = NODEL
 1010 CONTINUE
C
C     Input of permeabilities, construction of permeability matrix P
C
      WRITE (NOUT,FMT=9910)
      CALL MATNUL(P,IP,JP,DIMEN,DIMEN,ITEST)
      READ (NIN,FMT=9970) (P(I,I),I=1,DIMEN)
      WRITE (NOUT,FMT=9900) (P(I,I),I=1,DIMEN)
C
C
C     Input of boundary conditions and construction of nodal freedom
C     array NF
C
      WRITE (NOUT,FMT=9890)
      READ (NIN,FMT=9990) BNDNOD
      WRITE (NOUT,FMT=9940) BNDNOD
C
      DO 1020 I = 1,BNDNOD
          READ (NIN,FMT=9960) BNODE(I),BVAL(I)
          WRITE (NOUT,FMT=9930) BNODE(I),BVAL(I)
 1020 CONTINUE
C
C     Input and initialisation of free surface data. Surface
C
      WRITE (NOUT,FMT=9880)
      READ (NIN,FMT=9990) NUMFRE
      WRITE (NOUT,FMT=9940) NUMFRE
C
      READ (NIN,FMT=9990) (FRENOD(I),I=1,NUMFRE)
      WRITE (NOUT,FMT=9940) (FRENOD(I),I=1,NUMFRE)
      DO 1030 I = 1,NUMFRE
          NODE = FRENOD(I)
          PREVHT(I) = COORD(NODE,2)
 1030 CONTINUE
C
C     Nodes on the seepage boundary
C
      WRITE (NOUT,FMT=9870)
      READ (NIN,FMT=9990) NUMSEP
      WRITE (NOUT,FMT=9940) NUMSEP
C
      READ (NIN,FMT=9990) (SEPNOD(I),I=1,NUMSEP)
      WRITE (NOUT,FMT=9940) (SEPNOD(I),I=1,NUMSEP)
C
C     Modify the seepage boundary conditions
C
      DO 1060 I = 1,NUMSEP
C
C     Find corresponding boundary NODE
C
          DO 1040 J = 1,BNDNOD
              IF (SEPNOD(I).EQ.BNODE(J)) GO TO 1050
 1040     CONTINUE
C
C     UPDATE boundary conditions
C
 1050     CONTINUE
          SEPBND(I) = J
          NODE = SEPNOD(I)
          BVAL(J) = COORD(NODE,2)
 1060 CONTINUE
C
C     Input iteration control data
C
      READ (NIN,FMT=9960) MAXIT,CFACTR
      WRITE (NOUT,FMT=9860) MAXIT,CFACTR
C
C     Form nodal freedom array
C
      TOTDOF = 0
      DO 1080 I = 1,TOTNOD
          DO 1070 J = 1,DOFNOD
              TOTDOF = TOTDOF + 1
              NF(I,J) = TOTDOF
 1070     CONTINUE
 1080 CONTINUE
C
C     Calculation of semi-bandwidth
C
      CALL BNDWTH(ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,TOTELS,HBAND,
     *            ITEST)
      WRITE (NOUT,FMT=9850) HBAND
C
C     *                           *************************
C     *                           *                       *
C     *                           * Free Surface Solution *
C     *                           *                       *
C     *                           *************************
C
      CALL QQUA4(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST)
C
C     Free surface iteration loop
C
      ITER = 0
 1090 CONTINUE
      ITER = ITER + 1
C
C     System matrix assembly
C
      CALL MATNUL(SYSK,ISYSK,JSYSK,TOTDOF,HBAND,ITEST)
      CALL VECNUL(RHS,IRHS,TOTDOF,ITEST)
      DO 1130 NELE = 1,TOTELS
          CALL ELGEOM(NELE,ELTOP,IELTOP,JELTOP,COORD,ICOORD,JCOORD,GEOM,
     *                IGEOM,JGEOM,DIMEN,ITEST)
          NODEL = ELTOP(NELE,2)
          DOFEL = NODEL*DOFNOD
C
C     Integration loop for element stiffness using NQP quadrature
C     points
C
          CALL MATNUL(ELK,IELK,JELK,DOFEL,DOFEL,ITEST)
          DO 1120 IQUAD = 1,NQP
C
C     Form linear shape function and space derivatives in the local
C     corrdinates.
C
              XI = ABSS(1,IQUAD)
              ETA = ABSS(2,IQUAD)
              CALL QUAM4(FUN,IFUN,LDER,ILDER,JLDER,XI,ETA,ITEST)
C
C     Transform local derivatives to global coordinate system
C
              CALL MATMUL(LDER,ILDER,JLDER,GEOM,IGEOM,JGEOM,JAC,IJAC,
     *                    JJAC,DIMEN,NODEL,DIMEN,ITEST)
              CALL MATINV(JAC,IJAC,JJAC,JACIN,IJACIN,JJACIN,DIMEN,DET,
     *                    ITEST)
              CALL MATMUL(JACIN,IJACIN,JJACIN,LDER,ILDER,JLDER,GDER,
     *                    IGDER,JGDER,DIMEN,DIMEN,NODEL,ITEST)
C
C     Formation of element stiffness ELK
C
              CALL MATMUL(P,IP,JP,GDER,IGDER,JGDER,PD,IPD,JPD,DIMEN,
     *                    DIMEN,DOFEL,ITEST)
              CALL MATRAN(GDER,IGDER,JGDER,GDERT,IGDERT,JGDERT,DIMEN,
     *                    DOFEL,ITEST)
              CALL MATMUL(GDERT,IGDERT,JGDERT,PD,IPD,JPD,DTPD,IDTPD,
     *                    JDTPD,DOFEL,DIMEN,DOFEL,ITEST)
              CALL SCAPRD(GEOM(1,1),IGEOM,FUN,IFUN,NODEL,R,ITEST)
              QUOT = 2.D0*PI*R*ABS(DET)*WGHT(IQUAD)
              DO 1110 I = 1,DOFEL
                  DO 1100 J = 1,DOFEL
                      DTPD(I,J) = DTPD(I,J)*QUOT
 1100             CONTINUE
 1110         CONTINUE
              CALL MATADD(ELK,IELK,JELK,DTPD,IDTPD,JDTPD,DOFEL,DOFEL,
     *                    ITEST)
 1120     CONTINUE
C
C     Assembly of system stiffness matrix
C
          CALL DIRECT(NELE,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,STEER,
     *                ISTEER,ITEST)
          CALL ASSYM(SYSK,ISYSK,JSYSK,ELK,IELK,JELK,STEER,ISTEER,HBAND,
     *               DOFEL,ITEST)
 1130 CONTINUE
C
C     Modification of stiffness matrix and right-hand side to
C     implement boundary conditions
C
      DO 1140 I = 1,BNDNOD
          J = BNODE(I)
          SYSK(J,HBAND) = SYSK(J,HBAND)*SCALE
          RHS(J) = SYSK(J,HBAND)*BVAL(I)
 1140 CONTINUE
C
C     Solution of system matrix for the nodal values of the
C     potential
C
      CALL CHOSOL(SYSK,ISYSK,JSYSK,RHS,IRHS,TOTDOF,HBAND,ITEST)
C
C     Free surface UPDATE
C
      FSDIFF = 0.D0
      DO 1150 I = 1,NUMFRE
          NODE = FRENOD(I)
          AMOVE(1) = 0.D0
          AMOVE(2) = COORD(NODE,2) - RHS(NODE)
          CALL ADJUST(ELTOP,IELTOP,JELTOP,COORD,ICOORD,JCOORD,AMOVE,
     *                IAMOVE,TOTELS,DIMEN,FRENOD,IFRNOD,NUMFRE,I,NSIDE,
     *                ITEST)
          FSDIFF = FSDIFF + ABS(RHS(NODE)-PREVHT(I))
          PREVHT(I) = RHS(NODE)
 1150 CONTINUE
      FSDIFF = FSDIFF/DBLE(REAL(NUMFRE-1))
C
C     UPDATE sepage boundary
C
      DO 1160 I = 1,NUMSEP
          J = SEPBND(I)
          NODE = SEPNOD(I)
          BVAL(J) = COORD(NODE,2)
 1160 CONTINUE
C
      WRITE (NOUT,FMT=9830) ITER,FSDIFF
      CALL PRTVAL(RHS,IRHS,NF,INF,JNF,DOFNOD,TOTNOD,NOUT,ITEST)
C
C     Check for converged results
C
      IF ((FSDIFF.GT.CFACTR) .AND. (ITER.LT.MAXIT)) GO TO 1090
C
C     Iteration converged.
C
      IF ((ITER.GE.MAXIT) .AND. (FSDIFF.GT.CFACTR)) WRITE (NOUT,
     *    FMT=9820) MAXIT
C
C     Output final results
C
      WRITE (NOUT,FMT=9810)
      DO 1170 I = 1,TOTNOD
          WRITE (NOUT,FMT=9930) I, (COORD(I,J),J=1,DIMEN)
 1170 CONTINUE
C
C     Program end
C
      STOP
C
 9990 FORMAT (16I5)
 9980 FORMAT (I5,3F10.5)
 9970 FORMAT (2F10.5)
 9960 FORMAT (I5,6F10.5)
 9950 FORMAT (/,/,' **** NODAL GEOMETRY ****',/,' ')
 9940 FORMAT (' ',16I5)
 9930 FORMAT (' ',I5,6F10.5)
 9920 FORMAT (/,/,' **** ELEMENT TOPOLOGY ****',/,' ')
 9910 FORMAT (/,/,' **** PERMEABILITIES ****',/,' ')
 9900 FORMAT (' ',2F10.5)
 9890 FORMAT (/,/,' **** BOUNDARY CONDITIONS ****',/,' ')
 9880 FORMAT (/,/,' **** FREE SURFACE NODES ****',/,' ')
 9870 FORMAT (/,/,' **** SEEPAGE BOUNDARY NODES ****',/,' ')
 9860 FORMAT (/,' MAXIMUM NUMBER OF ITERATIONS = ',I5,/,' CONVERGE',
     *       'NCE CRITERION = ',F10.5)
 9850 FORMAT (/,' SEMI-BANDWIDTH = ',I5)
 9840 FORMAT (/,' *** ERROR : HBAND.GT.JSYSK')
 9830 FORMAT (/,' FREE SURFACE ITERATION NO.',I5,/,' AVERAGE SURF',
     *       'ACE MOVEMENT =',F10.5)
 9820 FORMAT (/,' *** WARNING : NO CONVERGENCE AFTER ',I5,' ITERATIONS')
 9810 FORMAT (/,/,' **** FINAL NODAL GEOMETRY ****',/,' ')
C
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ADJUST(ELTOP,IELTOP,JELTOP,COORD,ICOORD,JCOORD,AMOVE,
     *                  IAMOVE,TOTELS,DIMEN,FRENOD,IFRNOD,NUMFRE,FSPOS,
     *                  NSIDE,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      The ADJUST routine redefines a mesh between a 'top' surface
C      which has been moved and a constant 'base'. the routine
C      takes the node on the top surface, finds the nodes
C      between it and the constant base and moves them (if
C      necessary) to create a well-defined mesh.
C
C HISTORY
C
C      Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
C                           Chilton, Didcot, Oxfordshire OX11 0QX
C
C      Contact: Prof Chris Greenough
C      Email: c.greenough@rl.ac.uk
C      Tel: (44) 1235 445307
C
C      Release 3.0  Oct 1985
C      Release 4.0  Dec 2003 (CG)
C
C ARGUMENTS in
C      ELTOP   the element topology array containing element
C              type number, number of nodes on element and
C              topology list.
C      IELTOP  first dimension of ELTOP
C      JELTOP  second dimension of ELTOP
C      COORD   the global coordinate array
C      ICOORD  first dimension of COORD
C      JCOORD  second dimension of COORD
C      AMOVE   real array containing the amounts by which the
C              top node is to move in the coordinate directions
C      IAMOVE  dimension of AMOVE
C      TOTELS  the total number of elements in the mesh
C      DIMEN   number of space dimensions
C      FRENOD  array containing the nodes which lie along the top
C              surface
C      IFRNOD  dimension of FRENOD
C      NUMFRE  number of nodes on the free surface
C      FSPOS   position of current node in array FRENOD
C      NSIDE   number of nodes on each SIDE of element
C      ITEST   error checking option
C
C ARGUMENTS out
C      COORD   contains adjusted values of coordinates
C
C ROUTINES called
C      ERRMES, WHELEM
C
C      SUBROUTINE ADJUST(ELTOP,IELTOP,JELTOP,COORD,ICOORD,JCOORD,AMOVE,
C     *                  IAMOVE,TOTELS,DIMEN,FRENOD,IFRNOD,NUMFRE,FSPOS,
C     *                  NSIDE,ITEST)
C-----------------------------------------------------------------------
C
      INTEGER ANODE,BNODE,CNODE,DIMEN,DNODE,ELEM,ELEM1,ELEM2,ELTOP,
     *        ENODE,FNODE,FORWD,FRENOD,FSPOS,IAMOVE,ICOORD,IELTOP,
     *        IERROR,IFRNOD,ITEST,JCOORD,JELTOP,JTEST,K,NODEL,NSIDE,
     *        NSIDE1,NUMFRE,OPPSID,SIDE,TOTELS
      INTEGER ERRMES
      DOUBLE PRECISION AMOVE,COORD,ELSID,TOTMOV
      DIMENSION AMOVE(IAMOVE),COORD(ICOORD,JCOORD),ELTOP(IELTOP,JELTOP),
     *          FRENOD(IFRNOD)
      CHARACTER*6 SRNAME
C
      DATA SRNAME/'ADJUST '/
C
C     Parameter checking
C
      JTEST = ITEST
      IF (JTEST.NE.-1) THEN
          IERROR = 0
          IF (TOTELS.LE.0) IERROR = 1
          IF (IELTOP.LT.TOTELS) IERROR = 2
          IF (JCOORD.LT.DIMEN) IERROR = 3
          IF (FSPOS.LE.0 .OR. FSPOS.GT.IFRNOD) IERROR = 4
          IF (NSIDE.LE.1) IERROR = 5
          IF (IAMOVE.LT.DIMEN) IERROR = 6
          ITEST = ERRMES(JTEST,IERROR,SRNAME)
          IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body of SUBROUTINE
C
      ANODE = FRENOD(FSPOS)
C
C     Determine DNODE
C
      NSIDE1 = NSIDE - 1
      IF (FSPOS.LT.NUMFRE) DNODE = FRENOD(FSPOS+NSIDE1)
      IF (FSPOS.EQ.NUMFRE) DNODE = FRENOD(FSPOS-NSIDE1)
C
C     Determine BNODE and ENODE
C
      CALL WHELEM(ELTOP,IELTOP,JELTOP,TOTELS,ANODE,DNODE,0,NSIDE1,ELEM,
     *            SIDE,OPPSID,FORWD,ITEST)
      NODEL = ELTOP(ELEM,2)
      IF (JTEST.NE.-1) THEN
          IERROR = 0
          IF (NODEL+2.GT.JELTOP .OR. NODEL.LT.4) IERROR = 7
          IF (IERROR.NE.0) ITEST = ERRMES(JTEST,IERROR,SRNAME)
          IF (ITEST.NE.0) RETURN
      END IF
      IF (FORWD.EQ.0) THEN
          BNODE = ELTOP(ELEM,OPPSID+2)
          ENODE = ELTOP(ELEM,MOD(OPPSID,NODEL)+3)
      ELSE
          BNODE = ELTOP(ELEM,MOD(OPPSID,NODEL)+3)
          ENODE = ELTOP(ELEM,OPPSID+2)
      END IF
 1000 CONTINUE
C
C     Calculate the length of the element SIDE and determine whether
C     to terminate the procedure.
C
      ELSID = 0.0D0
      TOTMOV = 0.D0
      DO 1010 K = 1,DIMEN
          ELSID = ELSID + (COORD(ANODE,K)-COORD(BNODE,K))*
     *            (COORD(ANODE,K)-COORD(BNODE,K))
          COORD(ANODE,K) = COORD(ANODE,K) - AMOVE(K)
          TOTMOV = TOTMOV + AMOVE(K)*AMOVE(K)
 1010 CONTINUE
      ELSID = SQRT(ELSID)
      TOTMOV = SQRT(TOTMOV)
      IF (TOTMOV.GE.ELSID/3.D0) THEN
C
C     Determine CNODE and FNODE
C
          CALL WHELEM(ELTOP,IELTOP,JELTOP,TOTELS,BNODE,ENODE,ELEM,
     *                NSIDE1,ELEM1,SIDE,OPPSID,FORWD,ITEST)
          IF (ELEM1.NE.0) THEN
              NODEL = ELTOP(ELEM1,2)
              IF (JTEST.NE.-1) THEN
                  IERROR = 0
                  IF (NODEL+2.GT.JELTOP .OR. NODEL.LT.4) IERROR = 7
                  IF (IERROR.NE.0) ITEST = ERRMES(JTEST,IERROR,SRNAME)
                  IF (ITEST.NE.0) RETURN
              END IF
              IF (FORWD.EQ.0) THEN
                  CNODE = ELTOP(ELEM1,OPPSID+2)
                  FNODE = ELTOP(ELEM1,MOD(OPPSID,NODEL)+3)
              ELSE
                  CNODE = ELTOP(ELEM1,MOD(OPPSID,NODEL)+3)
                  FNODE = ELTOP(ELEM1,OPPSID+2)
              END IF
              IF (FSPOS.NE.1 .AND. FSPOS.NE.NUMFRE) THEN
                  CALL WHELEM(ELTOP,IELTOP,JELTOP,TOTELS,BNODE,CNODE,
     *                        ELEM1,NSIDE1,ELEM2,SIDE,OPPSID,FORWD,
     *                        ITEST)
                  IF (ELEM2.EQ.0) RETURN
              END IF
C
C     Set up information for next layer of elements
C
              DO 1020 K = 1,DIMEN
                  AMOVE(K) = COORD(BNODE,K) -
     *                       0.5D0* (COORD(CNODE,K)+COORD(ANODE,K))
 1020         CONTINUE
              ANODE = BNODE
              DNODE = ENODE
              BNODE = CNODE
              ENODE = FNODE
              ELEM = ELEM1
              GO TO 1000
          END IF
      END IF
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE WHELEM(ELTOP,IELTOP,JELTOP,TOTELS,NODE1,NODE2,NELEM,
     *                  NSIDE1,ELEM,SIDE,OPPSID,FORWD,ITEST)
C
C-----------------------------------------------------------------------
C PURPOSE
C      The WHELEM routine takes two nodes and determines in which
C      element they lie on the same SIDE, which SIDE they lie on,
C      the SIDE opposite and whether frenod is defined in the same
C      sense as the element topology.
C
C HISTORY
C
C      Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
C                           Chilton, Didcot, Oxfordshire OX11 0QX
C
C      Contact: Prof Chris Greenough
C      Email: c.greenough@rl.ac.uk
C      Tel: (44) 1235 445307
C
C      Release 3.0  Oct 1985
C      Release 4.0  Dec 2003 (CG)
C
C ARGUMENTS in
C      ELTOP   the element topology array containing element
C              type number, number of nodes on element and
C              topology list.
C      IELTOP  first dimension of ELTOP
C      JELTOP  second dimension of ELTOP
C      TOTELS  total number of elements in the mesh
C      NODE1   first node of interest
C      NODE2   second node of interest
C      NELEM   (optional) element number which contains NODE1
C              and node2, but is not the element required
C      NSIDE1  (number of nodes on an element SIDE) - 1
C      ITEST   error checking option
C
C ARGUMENTS out
C      ELEM    required element number in which NODE1 and NODE2
C              lie on the same SIDE
C      SIDE    the local SIDE number in ELEM on which NODE1 and
C              NODE2 lie
C      OPPSID  the SIDE opposite to SIDE in element ELEM
C      FORWD   FORWD = 1 if free surface is defined in the same sense
C              as ELTOP, FORWD = 0 if free surface is defined in the
C              opposite sense to ELTOP
C
C ROUTINES called
C      ERRMES
C
C      SUBROUTINE WHELEM(ELTOP,IELTOP,JELTOP,TOTELS,NODE1,NODE2,NELEM,
C     *                  NSIDE1,ELEM,SIDE,OPPSID,FORWD,ITEST)
C-----------------------------------------------------------------------
C
      INTEGER ELEM,ELTOP,FORWD,IE,IELTOP,IERROR,ITEST,JELTOP,JTEST,L1,
     *        L2,NELEM,NODE1,NODE2,NODEL,NSIDE1,OPPSID,SIDE,TOTELS
      INTEGER ERRMES
      CHARACTER*6 SRNAME
      DIMENSION ELTOP(IELTOP,JELTOP)
C
      DATA SRNAME/'WHELEM'/
C
C     Parameter checking
C
      JTEST = ITEST
      IF (JTEST.NE.-1) THEN
          IERROR = 0
          IF (NELEM.LT.0) IERROR = 1
          IF (NODE1.LE.0 .OR. NODE2.LE.0) IERROR = 2
          ITEST = ERRMES(JTEST,IERROR,SRNAME)
          IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body of SUBROUTINE
C
      ELEM = 0
C
C     Find element(s) which contain NODE1 and NODE2 on the same SIDE
C
      DO 1030 IE = 1,TOTELS
C
C     Select next element if IE = NELEM
C
          IF (IE.NE.NELEM) THEN
              NODEL = ELTOP(IE,2)
              DO 1000 L1 = 1,NODEL
                  IF (NODE1.EQ.ELTOP(IE,L1+2)) GO TO 1010
 1000         CONTINUE
              GO TO 1030
 1010         CONTINUE
              DO 1020 L2 = 1,NODEL
                  IF (NODE2.EQ.ELTOP(IE,L2+2)) GO TO 1040
 1020         CONTINUE
          END IF
 1030 CONTINUE
C
C     No element found (other than NELEM), so return
C
      RETURN
C
 1040 CONTINUE
      ELEM = IE
C
C     Calculate the SIDE number (SIDE), number of the opposite SIDE
C     (OPPSID) and variable FORWD (indication of the direction of
C     the free surface nodes).
C
      IF (L2.EQ.MOD(L1,NODEL)+NSIDE1) THEN
          FORWD = 1
          SIDE = (L1-1)/NSIDE1 + 1
          OPPSID = MOD(SIDE+1,4) + 1
      ELSE
          FORWD = 0
          SIDE = (L2-1)/NSIDE1 + 1
          OPPSID = MOD(SIDE+1,4) + 1
      END IF
      END
