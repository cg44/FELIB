C $Id: seg3p4.f,v 1.3 2009/02/26 15:05:33 cg44 Exp $
      PROGRAM SEG3P4
C***********************************************************************
C
C    COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                         CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C***********************************************************************
C
      INTEGER BNDNOD, BNODE, DIMEN, DOFEL, DOFNOD, ELNUM, ELTOP,
     *     ELTYP, FRENOD, HBAND, I, IABSS, IAMOVE, IBNODE,
     *     ICOORD, IDTPD, IELK, IELTOP, IFRNOD, IFUN, IGDER,
     *     IGDERT, IGEOM, IJAC, IJACIN, ILDER, INF, IP, IPD,
     *     IQUAD, IRHS, ISPNOD, ISTEER, ISYSK, ITER, ITEST,
     *     IWGHT, J, JABSS, JCOORD, JDTPD, JELK, JELTOP, JGDER,
     *     JGDERT, JGEOM, JJAC, JJACIN, JLDER, JNF, JP, JPD,
     *     JSYSK, MAXIT, NELE, NF, NIN, NODE, NODEL, NODNUM,
     *     NOUT, NQP, NSIDE, NUMFRE, NUMSEP, SEPBND, SEPNOD,
     *     STEER, TOTDOF, TOTELS, TOTNOD
      DOUBLE PRECISION ABSS, AMOVE, BVAL, CFACTR, COORD, DET,
     *     DTPD, ELK, ETA, FSDIFF, FUN, GDER, GDERT, GEOM,
     *     JAC, JACIN, LDER, P, PD, PI, PREVHT, QUOT, R, RHS,
     *     SCALE, SYSK, WGHT, XI
      DIMENSION ABSS(3,9), AMOVE(3), DTPD(8,8), ELK(8,8),
     *     FUN(8), GDER(3,8), GDERT(8,3), GEOM(8,3), JAC(3,3),
     *     JACIN(3,3), LDER(3,8), P(3,3), PD(3,8), STEER(8),
     *     WGHT(9)
C
C                            PROBLEM SIZE DEPENDENT ARRAYS
C
      DIMENSION BNODE(50), BVAL(50), COORD(250,3),
     *     ELTOP(350,10), FRENOD(50), NF(250,1), PREVHT(50),
     *     RHS(250), SEPBND(50), SEPNOD(50), SYSK(250,25)
C
      DATA IABSS /3/, IAMOVE /3/, IDTPD /8/, IELK /8/, IFUN /8/,
     *     IGDER /3/, IGDERT /8/, IGEOM /8/, IJAC /3/, IJACIN /3/,
     *     ILDER /3/, IP /3/, IPD /3/, ISTEER /8/, IWGHT /9/,
     *     JABSS /9/, JCOORD /3/, JDTPD /8/, JELK /8/, JGDER /8/,
     *     JGDERT /3/, JGEOM /3/, JJAC /3/, JJACIN /3/,
     *     JLDER /8/, JNF /1/, JP /3/, JPD /8/, SCALE /1.0D+10/
C
C                            PROBLEM SIZE DEPENDENT DATA STATEMENTS
C
      DATA IBNODE /50/, ICOORD /250/, IELTOP /350/, IFRNOD /50/,
     *     INF /250/, IRHS /250/, ISPNOD /50/, ISYSK /250/,
     *     JELTOP /10/, JSYSK /25/, NSIDE /2/
C
      DATA DOFNOD /1/, NIN /5/, NOUT /6/, ADOUT /10/
C
C                            SET ITEST FOR FULL CHECKING
C
      ITEST = 0
C
      PI = DATAN(1.0D0)*4.0D0
C
C*                           **********************
C*                           *                    *
C*                           * INPUT DATA SECTION *
C*                           *                    *
C*                           **********************
C
C                            INPUT OF NODAL GEOMETRY
C
      WRITE (NOUT,9010)
      READ (NIN,8010) TOTNOD, DIMEN
      WRITE (NOUT,9020) TOTNOD, DIMEN
      IF ((TOTNOD.GT.0) .AND. (TOTNOD.LE.ICOORD)) GO TO 1010
      WRITE (NOUT,9030)
      STOP
C
 1010 DO 1020 I=1,TOTNOD
      READ (NIN,8020) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
      WRITE (NOUT,9040) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
 1020 CONTINUE
C
C                            INPUT OF ELEMENT TOPOLOGY
C
      WRITE (NOUT,9050)
      READ (NIN,8010) TOTELS
      WRITE (NOUT,9020) TOTELS
      IF ((TOTELS.GT.0) .AND. (TOTELS.LE.IELTOP)) GO TO 1030
      WRITE (NOUT,9060)
      STOP
C
 1030 DO 1040 I=1,TOTELS
      READ (NIN,8010) ELNUM, ELTYP, NODEL, (ELTOP(ELNUM,J+2),J=1,
     *     NODEL)
      WRITE (NOUT,9020) ELNUM, ELTYP, NODEL, (ELTOP(ELNUM,J+2),J=1,
     *     NODEL)
      ELTOP(ELNUM,1) = ELTYP
      ELTOP(ELNUM,2) = NODEL
 1040 CONTINUE
C
C                            INPUT OF PERMEABILITIES, CON-
C                            STRUCTION OF PERMEABILITY MATRIX P
C
      WRITE (NOUT,9070)
      CALL MATNUL(P, IP, JP, DIMEN, DIMEN, ITEST)
      READ (NIN,8030) (P(I,I),I=1,DIMEN)
      WRITE (NOUT,9080) (P(I,I),I=1,DIMEN)
C
C
C                            INPUT OF BOUNDARY CONDITIONS AND
C                            CONSTRUCTION OF NODAL FREEDOM ARRAY NF
C
      WRITE (NOUT,9090)
      READ (NIN,8010) BNDNOD
      WRITE (NOUT,9020) BNDNOD
      IF ((BNDNOD.GT.0) .AND. (BNDNOD.LE.IBNODE)) GO TO 1050
      WRITE (NOUT,9100)
      STOP
C
 1050 DO 1060 I=1,BNDNOD
      READ (NIN,8040) BNODE(I), BVAL(I)
      WRITE (NOUT,9040) BNODE(I), BVAL(I)
 1060 CONTINUE
C
C                            INPUT AND INITIALISATION OF FREE
C                            SURFACE DATA.
C                            SURFACE
C
      WRITE (NOUT,9110)
      READ (NIN,8010) NUMFRE
      WRITE (NOUT,9020) NUMFRE
      IF ((NUMFRE.GT.0) .AND. (NUMFRE.LE.IFRNOD)) GO TO 1070
      WRITE (NOUT,9120)
      STOP
C
 1070 READ (NIN,8010) (FRENOD(I),I=1,NUMFRE)
      WRITE (NOUT,9020) (FRENOD(I),I=1,NUMFRE)
      DO 1080 I=1,NUMFRE
      NODE = FRENOD(I)
      PREVHT(I) = COORD(NODE,2)
 1080 CONTINUE
C
C                            NODES ON THE SEEPAGE BOUNDARY
C
      WRITE (NOUT,9130)
      READ (NIN,8010) NUMSEP
      WRITE (NOUT,9020) NUMSEP
      IF ((NUMSEP.GT.0) .AND. (NUMSEP.LE.ISPNOD)) GO TO 1090
      WRITE (NOUT,9140)
      STOP
C
 1090 READ (NIN,8010) (SEPNOD(I),I=1,NUMSEP)
      WRITE (NOUT,9020) (SEPNOD(I),I=1,NUMSEP)
C
C                          CHECK DATA
C
      DO 1120 I=1,NUMSEP
      DO 1100 J=1,BNDNOD
      IF (SEPNOD(I).EQ.BNODE(J)) GO TO 1110
 1100 CONTINUE
      WRITE (NOUT,9150) SEPNOD(I)
      STOP
C
 1110 SEPBND(I) = J
      NODE = SEPNOD(I)
      BVAL(J) = COORD(NODE,2)
 1120 CONTINUE
C
C                            INPUT ITERATION CONTROL DATA
C
      READ (NIN,8040) MAXIT, CFACTR
      WRITE (NOUT,9160) MAXIT, CFACTR
C
C                            FORM NODAL FREEDOM ARRAY
C
      TOTDOF = 0
      DO 1140 I=1,TOTNOD
      DO 1130 J=1,DOFNOD
      TOTDOF = TOTDOF + 1
      NF(I,J) = TOTDOF
 1130 CONTINUE
 1140 CONTINUE
C
C                            CALCULATION OF SEMI-BANDWIDTH
C
      CALL BNDWTH(ELTOP, IELTOP, JELTOP, NF, INF, JNF, DOFNOD,
     *     TOTELS, HBAND, ITEST)
      WRITE (NOUT,9170) HBAND
      IF (HBAND.LE.JSYSK) GO TO 1150
      WRITE (NOUT,9180)
      STOP
C
C*                           *************************
C*                           *                       *
C*                           * FREE SURFACE SOLUTION *
C*                           *                       *
C*                           *************************
C
 1150 CALL QQUA4(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)
C
C                            FREE SURFACE ITERATION LOOP
C
      ITER = 0
 1160 ITER = ITER + 1
C
C                            SYSTEM MATRIX ASSEMBLY
C
      CALL MATNUL(SYSK, ISYSK, JSYSK, TOTDOF, HBAND, ITEST)
      CALL VECNUL(RHS, IRHS, TOTDOF, ITEST)
      DO 1200 NELE=1,TOTELS
      CALL ELGEOM(NELE, ELTOP, IELTOP, JELTOP, COORD, ICOORD,
     *     JCOORD, GEOM, IGEOM, JGEOM, DIMEN, ITEST)
      NODEL = ELTOP(NELE,2)
      DOFEL = NODEL*DOFNOD
C
C                            INTEGRATION LOOP FOR ELEMENT STIFFNESS
C                            USING NQP QUADRATURE POINTS
C
      CALL MATNUL(ELK, IELK, JELK, DOFEL, DOFEL, ITEST)
      DO 1190 IQUAD=1,NQP
C
C                            FORM LINEAR SHAPE FUNCTION AND SPACE
C                            DERIVATIVES IN THE LOCAL CORRDINATES.
C
      XI = ABSS(1,IQUAD)
      ETA = ABSS(2,IQUAD)
      CALL QUAM4(FUN, IFUN, LDER, ILDER, JLDER, XI, ETA, ITEST)
C
C                            TRANSFORM LOCAL DERIVATIVES TO GLOBAL
C                            COORDINATE SYSTEM
C
      CALL MATMUL(LDER, ILDER, JLDER, GEOM, IGEOM, JGEOM, JAC,
     *     IJAC, JJAC, DIMEN, NODEL, DIMEN, ITEST)
      CALL MATINV(JAC, IJAC, JJAC, JACIN, IJACIN, JJACIN, DIMEN,
     *     DET, ITEST)
      CALL MATMUL(JACIN, IJACIN, JJACIN, LDER, ILDER, JLDER, GDER,
     *     IGDER, JGDER, DIMEN, DIMEN, NODEL, ITEST)
C
C                            FORMATION OF ELEMENT STIFFNESS ELK
C
      CALL MATMUL(P, IP, JP, GDER, IGDER, JGDER, PD, IPD, JPD,
     *     DIMEN, DIMEN, DOFEL, ITEST)
      CALL MATRAN(GDER, IGDER, JGDER, GDERT, IGDERT, JGDERT,
     *     DIMEN, DOFEL, ITEST)
      CALL MATMUL(GDERT, IGDERT, JGDERT, PD, IPD, JPD, DTPD,
     *     IDTPD, JDTPD, DOFEL, DIMEN, DOFEL, ITEST)
      CALL SCAPRD(GEOM(1,1), IGEOM, FUN, IFUN, NODEL, R, ITEST)
      QUOT = 2.D0*PI*R*DABS(DET)*WGHT(IQUAD)
      DO 1180 I=1,DOFEL
      DO 1170 J=1,DOFEL
      DTPD(I,J) = DTPD(I,J)*QUOT
 1170 CONTINUE
 1180 CONTINUE
      CALL MATADD(ELK, IELK, JELK, DTPD, IDTPD, JDTPD, DOFEL,
     *     DOFEL, ITEST)
 1190 CONTINUE
C
C                            ASSEMBLY OF SYSTEM STIFFNESS MATRIX
C
      CALL DIRECT(NELE, ELTOP, IELTOP, JELTOP, NF, INF, JNF,
     *     DOFNOD, STEER, ISTEER, ITEST)
      CALL ASSYM(SYSK, ISYSK, JSYSK, ELK, IELK, JELK, STEER,
     *     ISTEER, HBAND, DOFEL, ITEST)
 1200 CONTINUE
C
C                            MODIFICATION OF STIFFNESS MATRIX AND
C                            RIGHT-HAND SIDE TO IMPLEMENT BOUNDARY
C                            CONDITIONS
C
      DO 1210 I=1,BNDNOD
      J = BNODE(I)
      SYSK(J,HBAND) = SYSK(J,HBAND)*SCALE
      RHS(J) = SYSK(J,HBAND)*BVAL(I)
 1210 CONTINUE
C
C                            SOLUTION OF SYSTEM MATRIX FOR THE
C                            NODAL VALUES OF THE POTENTIAL
C
      CALL CHOSOL(SYSK, ISYSK, JSYSK, RHS, IRHS, TOTDOF, HBAND,
     *     ITEST)
C
C                            FREE SURFACE UPDATE
C
      FSDIFF = 0.D0
      DO 1220 I=1,NUMFRE
      NODE = FRENOD(I)
      AMOVE(1) = 0.D0
      AMOVE(2) = COORD(NODE,2) - RHS(NODE)
      CALL ADJUST(ELTOP, IELTOP, JELTOP, COORD, ICOORD, JCOORD,
     *     AMOVE, IAMOVE, TOTELS, DIMEN, FRENOD, IFRNOD, NUMFRE,
     *     I, NSIDE, ITEST)
      FSDIFF = FSDIFF + DABS(RHS(NODE)-PREVHT(I))
      PREVHT(I) = RHS(NODE)
 1220 CONTINUE
      FSDIFF = FSDIFF/DBLE(FLOAT(NUMFRE-1))
C
C                            UPDATE SEPAGE BOUNDARY
C
      DO 1230 I=1,NUMSEP
      J = SEPBND(I)
      NODE = SEPNOD(I)
      BVAL(J) = COORD(NODE,2)
 1230 CONTINUE
C
      WRITE (NOUT,9190) ITER, FSDIFF
      CALL PRTVAL(RHS, IRHS, NF, INF, JNF, DOFNOD, TOTNOD, NOUT,
     *     ITEST)
C
C                            CHECK FOR CONVERGED RESULTS
C
      IF ((FSDIFF.GT.CFACTR) .AND. (ITER.LT.MAXIT)) GO TO 1160
C
C                            ITERATION CONVERGED.
C
      IF ((ITER.GE.MAXIT) .AND. (FSDIFF.GT.CFACTR)) WRITE
     *     (NOUT,9200) MAXIT
C
C                            OUTPUT FINAL RESULTS
C
      WRITE (NOUT,9210)
      DO 1240 I=1,TOTNOD
      WRITE (NOUT,9040) I, (COORD(I,J),J=1,DIMEN)
 1240 CONTINUE
C
C                            PROGRAM END
C
      STOP
C
 8010 FORMAT (16I5)
 8020 FORMAT (I5, 3F10.5)
 8030 FORMAT (2F10.5)
 8040 FORMAT (I5, 6F10.5)
 9010 FORMAT (//25H **** NODAL GEOMETRY ****/1H )
 9020 FORMAT (1H , 16I5)
 9030 FORMAT (/44H *** ERROR : TOTNOD.LE.0 OR TOTNOD.GT.ICOORD)
 9040 FORMAT (1H , I5, 6F10.5)
 9050 FORMAT (//27H **** ELEMENT TOPOLOGY ****/1H )
 9060 FORMAT (/44H *** ERROR : TOTELS.LE.0 OR TOTELS.GT.IELTOP)
 9070 FORMAT (//25H **** PERMEABILITIES ****/1H )
 9080 FORMAT (1H , 2F10.5)
 9090 FORMAT (//30H **** BOUNDARY CONDITIONS ****/1H )
 9100 FORMAT (/44H *** ERROR : BNDNOD.LE.0 OR BNDNOD.GT.IBNODE)
 9110 FORMAT (//29H **** FREE SURFACE NODES ****/1H )
 9120 FORMAT (/44H *** ERROR : NUMFRE.LE.0 OR NUMFRE.GT.IFRNOD)
 9130 FORMAT (//33H **** SEEPAGE BOUNDARY NODES ****/1H )
 9140 FORMAT (/44H *** ERROR : NUMSEP.LE.0 OR NUMSEP.GT.ISPNOD)
 9150 FORMAT (/18H *** ERROR : NODE , I5, 21H IS NOT IN ARRAY BNOD,
     *     1HE)
 9160 FORMAT (/32H MAXIMUM NUMBER OF ITERATIONS = , I5/9H CONVERGE,
     *     16HNCE CRITERION = , F10.5)
 9170 FORMAT (/18H SEMI-BANDWIDTH = , I5)
 9180 FORMAT (/27H *** ERROR : HBAND.GT.JSYSK)
 9190 FORMAT (/27H FREE SURFACE ITERATION NO., I5/13H AVERAGE SURF,
     *     14HACE MOVEMENT =, F10.5)
 9200 FORMAT (/36H *** WARNING : NO CONVERGENCE AFTER , I5,
     *     11H ITERATIONS)
 9210 FORMAT (//31H **** FINAL NODAL GEOMETRY ****/1H )
      END
C
C******************************************************************
C
      SUBROUTINE ADJUST(ELTOP, IELTOP, JELTOP, COORD, ICOORD,
     *     JCOORD, AMOVE, IAMOVE, TOTELS, DIMEN, FRENOD, IFRNOD,
     *     NUMFRE, FSPOS, NSIDE, ITEST)
C------------------------------------------------------------------
C PURPOSE
C      THE ROUTINE REDEFINES A MESH BETWEEN A 'TOP' SURFACE
C      WHICH HAS BEEN MOVED AND A CONSTANT 'BASE'. THE ROUTINE
C      TAKES THE NODE ON THE TOP SURFACE, FINDS THE NODES
C      BETWEEN IT AND THE CONSTANT BASE AND MOVES THEM (IF
C      NECESSARY) TO CREATE A WELL-DEFINED MESH.
C
C HISTORY
C      RELEASE 3.0  OCT 1985
C
C ARGUMENTS IN
C      ELTOP   THE ELEMENT TOPOLOGY ARRAY CONTAINING ELEMENT
C              TYPE NUMBER, NUMBER OF NODES ON ELEMENT AND
C              TOPOLOGY LIST.
C      IELTOP  FIRST DIMENSION OF ELTOP
C      JELTOP  SECOND DIMENSION OF ELTOP
C      COORD   THE GLOBAL COORDINATE ARRAY
C      ICOORD  FIRST DIMENSION OF COORD
C      JCOORD  SECOND DIMENSION OF COORD
C      AMOVE   REAL ARRAY CONTAINING THE AMOUNTS BY WHICH THE
C              TOP NODE IS TO MOVE IN THE COORDINATE DIRECTIONS
C      IAMOVE  DIMENSION OF AMOVE
C      TOTELS  THE TOTAL NUMBER OF ELEMENTS IN THE MESH
C      DIMEN   NUMBER OF SPACE DIMENSIONS
C      FRENOD  ARRAY CONTAINING THE NODES WHICH LIE ALONG THE TOP
C              SURFACE
C      IFRNOD  DIMENSION OF FRENOD
C      NUMFRE  NUMBER OF NODES ON THE FREE SURFACE
C      FSPOS   POSITION OF CURRENT NODE IN ARRAY FRENOD
C      NSIDE   NUMBER OF NODES ON EACH SIDE OF ELEMENT
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      COORD   CONTAINS ADJUSTED VALUES OF COORDINATES
C
C ROUTINES CALLED
C      ERRMES, WHELEM
C
C      SUBROUTINE ADJUST(ELTOP, IELTOP, JELTOP, COORD, ICOORD,
C      *     JCOORD, AMOVE, IAMOVE, TOTELS, DIMEN, FRENOD, IFRNOD,
C      *     NUMFRE, FSPOS, NSIDE, ITEST)
C*******************************************************************
C
      INTEGER ANODE, BNODE, CNODE, DIMEN, DNODE, ELEM, ELEM1,
     *     ELEM2, ELTOP, ENODE, FNODE, FORWD, FRENOD, FSPOS,
     *     IAMOVE, IERROR, IFRNOD, JTEST, K, NODEL, NSIDE,
     *     NSIDE1, NUMFRE, OPPSID, SIDE, TOTELS
      INTEGER ERRMES
      DOUBLE PRECISION AMOVE, COORD, ELSID, SRNAME, TOTMOV
      DIMENSION AMOVE(IAMOVE), COORD(ICOORD,JCOORD),
     *     ELTOP(IELTOP,JELTOP), FRENOD(IFRNOD)
C
      DATA SRNAME /8H ADJUST /
C
C     PARAMETER CHECKING
C
      JTEST = ITEST
      IF (JTEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (TOTELS.LE.0) IERROR = 1
      IF (IELTOP.LT.TOTELS) IERROR = 2
      IF (JCOORD.LT.DIMEN) IERROR = 3
      IF (FSPOS.LE.0 .OR. FSPOS.GT.IFRNOD) IERROR = 4
      IF (NSIDE.LE.1) IERROR = 5
      IF (IAMOVE.LT.DIMEN) IERROR = 6
      ITEST = ERRMES(JTEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
C MAIN BODY OF SUBROUTINE
C
 1010 ANODE = FRENOD(FSPOS)
C
C DETERMINE DNODE
C
      NSIDE1 = NSIDE - 1
      IF (FSPOS.LT.NUMFRE) DNODE = FRENOD(FSPOS+NSIDE1)
      IF (FSPOS.EQ.NUMFRE) DNODE = FRENOD(FSPOS-NSIDE1)
C
C DETERMINE BNODE AND ENODE
C
      CALL WHELEM(ELTOP, IELTOP, JELTOP, TOTELS, ANODE, DNODE, 0,
     *     NSIDE1, ELEM, SIDE, OPPSID, FORWD, ITEST)
      NODEL = ELTOP(ELEM,2)
      IF (JTEST.EQ.-1) GO TO 1020
      IERROR = 0
      IF (NODEL+2.GT.JELTOP .OR. NODEL.LT.4) IERROR = 7
      IF (IERROR.NE.0) ITEST = ERRMES(JTEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
 1020 CONTINUE
      IF (FORWD.EQ.0) GO TO 1030
      BNODE = ELTOP(ELEM,MOD(OPPSID,NODEL)+3)
      ENODE = ELTOP(ELEM,OPPSID+2)
      GO TO 1040
 1030 BNODE = ELTOP(ELEM,OPPSID+2)
      ENODE = ELTOP(ELEM,MOD(OPPSID,NODEL)+3)
 1040 CONTINUE
 1050 CONTINUE
C
C CALCULATE THE LENGTH OF THE ELEMENT SIDE AND DETERMINE
C WHETHER TO TERMINATE THE PROCEDURE.
C
      ELSID = 0.0D0
      TOTMOV = 0.D0
      DO 1060 K=1,DIMEN
         ELSID = ELSID + (COORD(ANODE,K)-COORD(BNODE,K))*
     *        (COORD(ANODE,K)-COORD(BNODE,K))
         COORD(ANODE,K) = COORD(ANODE,K) - AMOVE(K)
         TOTMOV = TOTMOV + AMOVE(K)*AMOVE(K)
 1060 CONTINUE
      ELSID = DSQRT(ELSID)
      TOTMOV = DSQRT(TOTMOV)
      IF (TOTMOV.LT.ELSID/3.D0) RETURN
C
C DETERMINE CNODE AND FNODE
C
      CALL WHELEM(ELTOP, IELTOP, JELTOP, TOTELS, BNODE, ENODE,
     *     ELEM, NSIDE1, ELEM1, SIDE, OPPSID, FORWD, ITEST)
      IF (ELEM1.EQ.0) RETURN
      NODEL = ELTOP(ELEM1,2)
      IF (JTEST.EQ.-1) GO TO 1070
      IERROR = 0
      IF (NODEL+2.GT.JELTOP .OR. NODEL.LT.4) IERROR = 7
      IF (IERROR.NE.0) ITEST = ERRMES(JTEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
 1070 CONTINUE
      IF (FORWD.EQ.0) GO TO 1080
      CNODE = ELTOP(ELEM1,MOD(OPPSID,NODEL)+3)
      FNODE = ELTOP(ELEM1,OPPSID+2)
      GO TO 1090
 1080 CNODE = ELTOP(ELEM1,OPPSID+2)
      FNODE = ELTOP(ELEM1,MOD(OPPSID,NODEL)+3)
 1090 CONTINUE
      IF (FSPOS.EQ.1 .OR. FSPOS.EQ.NUMFRE) GO TO 1100
      CALL WHELEM(ELTOP, IELTOP, JELTOP, TOTELS, BNODE, CNODE,
     *     ELEM1, NSIDE1, ELEM2, SIDE, OPPSID, FORWD, ITEST)
      IF (ELEM2.EQ.0) RETURN
 1100 CONTINUE
C
C SET UP INFORMATION FOR NEXT LAYER OF ELEMENTS
C
      DO 1110 K=1,DIMEN
         AMOVE(K) = COORD(BNODE,K) - 0.5D0*(COORD(CNODE,K)
     *        +COORD(ANODE,K))
 1110 CONTINUE
      ANODE = BNODE
      DNODE = ENODE
      BNODE = CNODE
      ENODE = FNODE
      ELEM = ELEM1
      GO TO 1050
      END
C
C******************************************************************
C
      SUBROUTINE WHELEM(ELTOP, IELTOP, JELTOP, TOTELS, NODE1,
     *     NODE2, NELEM, NSIDE1, ELEM, SIDE, OPPSID, FORWD, ITEST)
C------------------------------------------------------------------
C PURPOSE
C      THE ROUTINE TAKES TWO NODES AND DETERMINES IN WHICH ELEMENT
C      THEY LIE ON THE SAME SIDE, WHICH SIDE THEY LIE ON, THE SIDE
C      OPPOSITE AND WHETHER FRENOD IS DEFINED IN THE SAME
C      SENSE AS THE ELEMENT TOPOLOGY.
C
C HISTORY
C      RELEASE 3.0  OCT 1985
C
C ARGUMENTS IN
C      ELTOP   THE ELEMENT TOPOLOGY ARRAY CONTAINING ELEMENT
C              TYPE NUMBER, NUMBER OF NODES ON ELEMENT AND
C              TOPOLOGY LIST.
C      IELTOP  FIRST DIMENSION OF ELTOP
C      JELTOP  SECOND DIMENSION OF ELTOP
C      TOTELS  TOTAL NUMBER OF ELEMENTS IN THE MESH
C      NODE1   FIRST NODE OF INTEREST
C      NODE2   SECOND NODE OF INTEREST
C      NELEM   (OPTIONAL) ELEMENT NUMBER WHICH CONTAINS NODE1
C              AND NODE2, BUT IS NOT THE ELEMENT REQUIRED
C      NSIDE1  (NUMBER OF NODES ON AN ELEMENT SIDE) - 1
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      ELEM    REQUIRED ELEMENT NUMBER IN WHICH NODE1 AND NODE2
C              LIE ON THE SAME SIDE
C      SIDE    THE LOCAL SIDE NUMBER IN ELEM ON WHICH NODE1 AND
C              NODE2 LIE
C      OPPSID  THE SIDE OPPOSITE TO SIDE IN ELEMENT ELEM
C      FORWD   FORWD = 1 IF FREE SURFACE IS DEFINED IN THE SAME SENSE
C              AS ELTOP, FORWD = 0 IF FREE SURFACE IS DEFINED IN THE
C              OPPOSITE SENSE TO ELTOP
C
C ROUTINES CALLED
C      ERRMES
C
C      SUBROUTINE WHELEM(ELTOP, IELTOP, JELTOP, TOTELS, NODE1,
C     *     NODE2, NELEM, NSIDE1, ELEM, SIDE, OPPSID, FORWD, ITEST)
C*******************************************************************
C
      INTEGER ELEM, ELTOP, FORWD, IE, IELTOP, IERROR,
     *     ITEST, JELTOP, JTEST, L1, L2, NODE1, NODE2, NODEL,
     *     OPPSID, SIDE, TOTELS
      INTEGER ERRMES
      DOUBLE PRECISION SRNAME
      DIMENSION ELTOP(IELTOP,JELTOP)
C
      DATA SRNAME /8H WHELEM /
C
C     PARAMETER CHECKING
C
      JTEST = ITEST
      IF (JTEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (NELEM.LT.0) IERROR = 1
      IF (NODE1.LE.0 .OR. NODE2.LE.0) IERROR = 2
      ITEST = ERRMES(JTEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
C MAIN BODY OF SUBROUTINE
C
 1010 ELEM = 0
C
C FIND ELEMENT(S) WHICH CONTAIN NODE1 AND NODE2 ON THE
C SAME SIDE
C
      DO 1050 IE=1,TOTELS
C
C SELECT NEXT ELEMENT IF IE = NELEM
C
         IF (IE.EQ.NELEM) GO TO 1050
         NODEL = ELTOP(IE,2)
         DO 1020 L1=1,NODEL
            IF (NODE1.EQ.ELTOP(IE,L1+2)) GO TO 1030
 1020    CONTINUE
         GO TO 1050
 1030    DO 1040 L2=1,NODEL
            IF (NODE2.EQ.ELTOP(IE,L2+2)) GO TO 1060
 1040    CONTINUE
 1050 CONTINUE
C
C NO ELEMENT FOUND (OTHER THAN NELEM), SO RETURN
C
      RETURN
C
 1060 ELEM = IE
C
C CALCULATE THE SIDE NUMBER (SIDE), NUMBER OF THE
C OPPOSITE SIDE (OPPSID) AND VARIABLE FORWD (INDICATION
C OF THE DIRECTION OF THE FREE SURFACE NODES).
C
      IF (L2.EQ.MOD(L1,NODEL)+NSIDE1) GO TO 1070
      FORWD = 0
      SIDE = (L2-1)/NSIDE1 + 1
      OPPSID = MOD(SIDE+1,4) + 1
      RETURN
 1070 FORWD = 1
      SIDE = (L1-1)/NSIDE1 + 1
      OPPSID = MOD(SIDE+1,4) + 1
      RETURN
      END
