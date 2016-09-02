C***********************************************************************
C
C    COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                         CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C***********************************************************************
C
      INTEGER BNDNOD, BNODE, DIF, DIMEN, DOFEL, DOFNOD, ELNUM,
     *     ELTOP, ELTYP, HBAND, I, IABSS, ICOORD, IDTPD, IELK,
     *     IELM, IELTOP, IFTF, IFUN, IGDER, IGDERT, IGEOM,
     *     IJAC, IJACIN, ILDER, INF, INTNOD, IP, IPD, IQUAD,
     *     IRHS, ISTEER, ISYSK, ISYSM, ITEST, IWGHT, IWORK1,
     *     IWORK2, J, JABSS, JCOORD, JDTPD, JELK, JELM, JELTOP,
     *     JFTF, JGDER, JGDERT, JGEOM, JJAC, JJACIN, JLDER,
     *     JNF, JP, JPD, JSYSK, JSYSM, K, NELE, NF, NIN, NODEL,
     *     NODNUM, NOUT, NQP, NSTEPS, STEER, TOTDOF, TOTELS,
     *     TOTNOD
      DOUBLE PRECISION ABSS, BVAL, COORD, DET, DTIM, DTPD,
     *     ELK, ELM, ETA, FTF, FUN, GDER, GDERT, GEOM, JAC,
     *     JACIN, LDER, P, PD, PERM, QUOT, RHS, SCALE, SYSK,
     *     SYSM, THETA, WGHT, WORK1, WORK2, XI
      LOGICAL FIRST
      DIMENSION ABSS(3,9), BNODE(30), BVAL(30), DTPD(8,8),
     *     ELK(8,8), ELM(8,8), FTF(8,8), FUN(8), GDER(3,8),
     *     GDERT(8,3), GEOM(8,3), JAC(3,3), JACIN(3,3),
     *     LDER(3,8), P(3,3), PD(3,8), PERM(3), STEER(8),
     *     WGHT(9)
C
C                            PROBLEM SIZE DEPENDENT ARRAYS
C
      DIMENSION COORD(100,3), ELTOP(100,10), NF(100,1),
     *     RHS(100), SYSK(100,25), SYSM(100,25), WORK1(100),
     *     WORK2(100)
C
      DATA IABSS /3/, IDTPD /8/, IELK /8/, IELM /8/, IFTF /8/,
     *     IFUN /8/, IGDER /3/, IGDERT /8/, IGEOM /8/, IJAC /3/,
     *     IJACIN /3/, ILDER /3/, IP /3/, IPD /3/, ISTEER /8/,
     *     IWGHT /9/, JABSS /9/, JCOORD /3/, JDTPD /8/,
     *     JELK /8/, JELM /8/, JFTF /8/, JGDER /8/, JGDERT /3/,
     *     JGEOM /3/, JJAC /3/, JJACIN /3/, JLDER /8/, JNF /1/,
     *     JP /3/, JPD /8/, SCALE /1.0D+10/
C
C                            PROBLEM SIZE DEPENDENT DATA STATEMENTS
C
      DATA ICOORD /100/, IELTOP /100/, INF /100/, IRHS /100/,
     *     ISYSK /100/, ISYSM /100/, IWORK1 /100/, IWORK2 /100/,
     *     JELTOP /10/, JSYSK /25/, JSYSM /25/
C
      DATA NIN /5/, NOUT /6/
C
C                            SET ITEST FOR FULL CHECKING
C
      ITEST = 0
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
      DO 1010 I=1,TOTNOD
      READ (NIN,8020) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
      WRITE (NOUT,9030) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
 1010 CONTINUE
C
C                            INPUT OF ELEMENT TOPOLOGY
C
      WRITE (NOUT,9040)
      READ (NIN,8010) ELTYP, TOTELS, NODEL
      WRITE (NOUT,9020) ELTYP, TOTELS, NODEL
      DO 1020 I=1,TOTELS
      READ (NIN,8010) ELNUM, (ELTOP(ELNUM,J+2),J=1,NODEL)
      WRITE (NOUT,9020) ELNUM, (ELTOP(ELNUM,J+2),J=1,NODEL)
      ELTOP(ELNUM,1) = ELTYP
      ELTOP(ELNUM,2) = NODEL
 1020 CONTINUE
C
C                            INPUT OF PERMEABILITIES AND CON-
C                            STRUCTION OF PERMEABILITY MATRIX P
C
      WRITE (NOUT,9050)
      READ (NIN,8030) (PERM(I),I=1,DIMEN)
      WRITE (NOUT,9060) (PERM(I),I=1,DIMEN)
      CALL MATNUL(P, IP, JP, DIMEN, DIMEN, ITEST)
      DO 1030 I=1,DIMEN
      P(I,I) = PERM(I)
 1030 CONTINUE
C
C                            INPUT OF NUMBER OF DEGREES OF FREEDOM
C                            PER NODE, INPUT OF BOUNDARY CON-
C                            DITIONS AND CONSTRUCTION OF NODAL
C                            FREEDOM ARRAY NF
C
      WRITE (NOUT,9070)
      READ (NIN,8010) DOFNOD
      WRITE (NOUT,9020) DOFNOD
      READ (NIN,8010) BNDNOD
      WRITE (NOUT,9020) BNDNOD
      DO 1040 I=1,BNDNOD
      READ (NIN,8040) BNODE(I), BVAL(I)
      WRITE (NOUT,9030) BNODE(I), BVAL(I)
 1040 CONTINUE
      TOTDOF = 0
      DO 1060 I=1,TOTNOD
      DO 1050 J=1,DOFNOD
      TOTDOF = TOTDOF + 1
      NF(I,J) = TOTDOF
 1050 CONTINUE
 1060 CONTINUE
C
C                            INPUT OF INITIAL CONDINTIONS
C
      WRITE (NOUT,9080)
      CALL VECNUL(RHS, IRHS, TOTDOF, ITEST)
      READ (NIN,8010) INTNOD
      WRITE (NOUT,9090) INTNOD
      DO 1070 I=1,INTNOD
      READ (NIN,8040) NODNUM, WORK1(1)
      WRITE (NOUT,9030) NODNUM, WORK1(1)
      J = NF(NODNUM,1)
      RHS(J) = WORK1(1)
 1070 CONTINUE
C
C                            INPUT TIME STEPPING DATA
C
      WRITE (NOUT,9100)
      READ (NIN,8020) NSTEPS, DTIM, THETA
      WRITE (NOUT,9110) NSTEPS, DTIM, THETA
C
C                            CALCULATION OF SEMI-BANDWIDTH
C
      FIRST = .TRUE.
      DO 1080 NELE=1,TOTELS
      CALL FREDIF(NELE, ELTOP, IELTOP, JELTOP, NF, INF, JNF,
     *     DOFNOD, FIRST, DIF, ITEST)
 1080 CONTINUE
      HBAND = DIF + 1
C
C*                           ************************************
C*                           *                                  *
C*                           * SYSTEM STIFFNESS MATRIX ASSEMBLY *
C*                           *                                  *
C*                           ************************************
C
      CALL MATNUL(SYSK, ISYSK, JSYSK, TOTDOF, HBAND, ITEST)
      CALL MATNUL(SYSM, ISYSM, JSYSM, TOTDOF, HBAND, ITEST)
      DOFEL = NODEL*DOFNOD
      CALL QQUA4(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)
      DO 1120 NELE=1,TOTELS
      CALL ELGEOM(NELE, ELTOP, IELTOP, JELTOP, COORD, ICOORD,
     *     JCOORD, GEOM, IGEOM, JGEOM, DIMEN, ITEST)
C
C                            INTEGRATION LOOP FOR ELEMENT STIFFNESS
C                            USING NQP QUADRATURE POINTS
C
      CALL MATNUL(ELK, IELK, JELK, DOFEL, DOFEL, ITEST)
      CALL MATNUL(ELM, IELM, JELM, DOFEL, DOFEL, ITEST)
      DO 1110 IQUAD=1,NQP
C
C                            FORM LINEAR SHAPE FUNCTION AND SPACE
C                            DERIVATIVES IN THE LOCAL CORRDINATES.
C                            TRANSFORM LOCAL DERIVATIVES TO GLOBAL
C                            COORDINATE SYSTEM
C
      XI = ABSS(1,IQUAD)
      ETA = ABSS(2,IQUAD)
      CALL QUAM4(FUN, IFUN, LDER, ILDER, JLDER, XI, ETA, ITEST)
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
      QUOT = DABS(DET)*WGHT(IQUAD)
      DO 1100 I=1,DOFEL
      DO 1090 J=1,DOFEL
      FTF(I,J) = FUN(I)*FUN(J)*QUOT
      DTPD(I,J) = DTPD(I,J)*QUOT
 1090 CONTINUE
 1100 CONTINUE
      CALL MATADD(ELM, IELM, JELM, FTF, IFTF, JFTF, DOFEL, DOFEL,
     *     ITEST)
      CALL MATADD(ELK, IELK, JELK, DTPD, IDTPD, JDTPD, DOFEL,
     *     DOFEL, ITEST)
 1110 CONTINUE
C
C                            ASSEMBLY OF SYSTEM STIFFNESS MATRIX
C                            SYSK
C
      CALL DIRECT(NELE, ELTOP, IELTOP, JELTOP, NF, INF, JNF,
     *     DOFNOD, STEER, ISTEER, ITEST)
      CALL ASSYM(SYSM, ISYSM, JSYSM, ELM, IELM, JELM, STEER,
     *     ISTEER, HBAND, DOFEL, ITEST)
      CALL ASSYM(SYSK, ISYSK, JSYSK, ELK, IELK, JELK, STEER,
     *     ISTEER, HBAND, DOFEL, ITEST)
 1120 CONTINUE
C
C*                           *********************
C*                           *                   *
C*                           * EQUATION SOLUTION *
C*                           *                   *
C*                           *********************
C
C                            MODIFICATION OF STIFFNESS MATRIX AND
C                            RIGHT-HAND SIDE TO IMPLEMENT BOUNDARY
C                            CONDITIONS
C
      DO 1140 I=1,TOTDOF
      DO 1130 J=1,HBAND
      WORK1(1) = THETA*SYSK(I,J) + SYSM(I,J)/DTIM
      WORK1(2) = (1.0D0-THETA)*SYSK(I,J) - SYSM(I,J)/DTIM
      SYSK(I,J) = WORK1(1)
      SYSM(I,J) = -WORK1(2)
 1130 CONTINUE
 1140 CONTINUE
      DO 1150 I=1,BNDNOD
      J = BNODE(I)
      SYSK(J,HBAND) = SYSK(J,HBAND)*SCALE
      WORK1(I) = SYSK(J,HBAND)*BVAL(I)
 1150 CONTINUE
C
C                            SOLUTION OF SYSTEM MATRIX FOR THE
C                            NODAL VALUES OF THE POTENTIAL
C
      CALL CHORDN(SYSK, ISYSK, JSYSK, TOTDOF, HBAND, ITEST)
      DO 1170 I=1,NSTEPS
      WRITE (NOUT,9120) I
      CALL MVSYB(SYSM, ISYSM, JSYSM, RHS, IRHS, WORK2, IWORK2,
     *     TOTDOF, HBAND, ITEST)
C
C                            INSERTION OF BOUNDARY CONDITIONS
C
      DO 1160 J=1,BNDNOD
      K = BNODE(J)
      WORK2(K) = WORK1(J)
 1160 CONTINUE
C
C                            FORWARD AND BACKWARD SUBSTITUTION
C
      CALL CHOSUB(SYSK, ISYSK, JSYSK, WORK2, IWORK2, TOTDOF,
     *     HBAND, ITEST)
      CALL VECCOP(WORK2, IWORK2, RHS, IRHS, TOTDOF, ITEST)
      CALL PRTVAL(RHS, IRHS, NF, INF, JNF, DOFNOD, TOTNOD, NOUT,
     *     ITEST)
 1170 CONTINUE
      STOP
 8010 FORMAT (16I5)
 8020 FORMAT (I5, 6F10.0)
 8030 FORMAT (2F10.0)
 8040 FORMAT (I5, 6F10.0)
 9010 FORMAT (//25H **** NODAL GEOMETRY ****//1H )
 9020 FORMAT (1H , 16I5)
 9030 FORMAT (1H , I5, 6F10.5)
 9040 FORMAT (//27H **** ELEMENT TOPOLOGY ****//1H )
 9050 FORMAT (//25H **** PERMEABILITIES ****//1H )
 9060 FORMAT (1H , 2F10.5)
 9070 FORMAT (//30H **** BOUNDARY CONDITIONS ****//1H )
 9080 FORMAT (//29H **** INITIAL CONDITIONS ****//1H )
 9090 FORMAT (1H , 16I5)
 9100 FORMAT (//29H **** TIME STEPPING DATA ****//1H )
 9110 FORMAT (1H , I5, 6F10.5)
 9120 FORMAT (//32H **** NODAL POTENTIALS FOR STEP , I3, 6H  ****//
     *     1H )
      END
C***********************************************************************
