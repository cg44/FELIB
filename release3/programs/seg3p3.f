C***********************************************************************
C
C    COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                         CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C***********************************************************************
C
      INTEGER BNDNOD, BNODE, DIMEN, DOFEL, DOFNOD, ELNUM, ELTOP,
     *     ELTYP, HBAND, I, IABSS, IBNODE, ICOORD, IDTPD, IELJ,
     *     IELK, IELM, IELPOT, IELTOP, IFIELD, IFTF, IFUN,
     *     IGDER, IGDERT, IGEOM, IJAC, IJACIN, IJVAL, ILDER,
     *     INF, IPD, IQUAD, IRHS, IS, ISCVEC, ISTEER, ISYSK,
     *     ITEST, IWGHT, J, JABSS, JCOORD, JDTPD, JELK, JELM,
     *     JELTOP, JFTF, JGDER, JGDERT, JGEOM, JJAC, JJACIN,
     *     JLDER, JNF, JPD, JS, JSTEER, JSYSK, NELE, NF, NIN,
     *     NJV, NODEL, NODNUM, NOUT, NQP, NTEMP, STEER, TOTDOF,
     *     TOTELS, TOTNOD
      DOUBLE PRECISION ABSS, BMOD, BVAL, COND, CONST, COORD,
     *     DET, DTPD, ELJ, ELK, ELM, ELPOT, ETA, FIELD, FREQ,
     *     FTF, FUN, GDER, GDERT, GEOM, JAC, JACIN, JVAL, LDER,
     *     MU0, OMEGA, PD, PERM, PI, QUOT, RHS, S, SCALE, SCVEC,
     *     SYSK, WGHT, XI
      DIMENSION ABSS(3,9), DTPD(8,8), ELJ(8), ELK(8,8),
     *     ELM(8,8), ELPOT(8), FIELD(8), FTF(8,8), FUN(8),
     *     GDER(3,8), GDERT(8,3), GEOM(8,3), JAC(3,3),
     *     JACIN(3,3), LDER(3,8), PD(3,8), S(3,3), SCVEC(8),
     *     STEER(8), WGHT(9)
C
C                            PROBLEM SIZE DEPENDENT ARRAYS
C
      DIMENSION BNODE(30), BVAL(2,30), COORD(100,3),
     *     ELTOP(100,10), JVAL(100), NF(100,1), RHS(2,100),
     *     SYSK(2,100,25)
C
      DATA IABSS /3/, IDTPD /8/, IELJ /8/, IELK /8/, IELM /8/,
     *     IELPOT /8/, IFIELD /8/, IFTF /8/, IFUN /8/, IGDER /3/,
     *     IGDERT /8/, IGEOM /8/, IJAC /3/, IJACIN /3/,
     *     ILDER /3/, IPD /3/, IS /3/, ISCVEC /8/, ISTEER /8/,
     *     IWGHT /9/, JABSS /9/, JCOORD /3/, JDTPD /8/,
     *     JELK /8/, JELM /8/, JFTF /8/, JGDER /8/, JGDERT /3/,
     *     JGEOM /3/, JJAC /3/, JJACIN /3/, JLDER /8/, JNF /1/,
     *     JPD /8/, JS /3/, SCALE /1.0D+10/
C
C                            PROBLEM SIZE DEPENDENT DATA STATEMENTS
C
      DATA IBNODE /30/, ICOORD /100/,  IELTOP /100/, IJVAL /100/,
     *     INF /100/, IRHS /100/, ISYSK /100/, JELTOP /10/, JSYSK /25/
C
      DATA NIN /5/, NOUT /6/, NTEMP /7/
C
C                            SET ITEST FOR FULL CHECKING
C
      ITEST = 0
C
C                            SET VALUE OF PI
C
      PI = DATAN(1.0D0)*4.0D0
      MU0 = 4.0D0*PI*1.0D-07
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
      IF (TOTNOD.GT.0 .AND. TOTNOD.LE.ICOORD) GO TO 1010
      WRITE (NOUT,9030)
      STOP
C
 1010 DO 1020 I=1,TOTNOD
      READ (NIN,8020) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
      WRITE (NOUT,9040) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
 1020 CONTINUE
C
C                            INPUT OF ELEMENT TOPOLOGY ELTYP=2 FOR
C                            CONDUCTOR, ELTYP=1 OTHERWISE
C
      WRITE (NOUT,9050)
      READ (NIN,8010) TOTELS
      WRITE (NOUT,9020) TOTELS
      IF (TOTELS.GT.0 .AND. TOTELS.LE.IELTOP) GO TO 1030
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
C                            INPUT OF PERMEABILITIES AND CON-
C                            STRUCTION OF SUSCEPTIBILITY MATRIX S
C
C
      WRITE (NOUT,9070)
      READ (NIN,8030) FREQ, COND, PERM
      WRITE (NOUT,9080) FREQ, COND, PERM
      OMEGA = 2.0D0*PI*FREQ
C
C                            INPUT OF NUMBER OF DEGREES OF FREEDOM
C                            PER NODE, INPUT OF BOUNDARY CONDITIONS
C
      WRITE (NOUT,9090)
      READ (NIN,8010) DOFNOD
      WRITE (NOUT,9020) DOFNOD
C
      READ (NIN,8010) BNDNOD
      WRITE (NOUT,9020) BNDNOD
      IF (BNDNOD.GT.0 .AND. BNDNOD.LE.IBNODE) GO TO 1050
      WRITE (NOUT,9100)
      STOP
C
 1050 DO 1060 I=1,BNDNOD
      READ (NIN,8020) BNODE(I), BVAL(1,I), BVAL(2,I)
      WRITE (NOUT,9040) BNODE(I), BVAL(1,I), BVAL(2,I)
 1060 CONTINUE
C
C                            INPUT NUMBER OF CURRENT CARRYING
C                            ELEMENTS INPUT ELEMENT NUMBER
C                            AND CURRENT VALUE
C
      WRITE (NOUT,9110)
      CALL VECNUL(JVAL, IJVAL, TOTELS, ITEST)
      READ (NIN,8010) NJV
      WRITE (NOUT,9020) NJV
      IF (NJV.EQ.0) GO TO 1090
      IF (NJV.LE.IJVAL) GO TO 1070
      WRITE (NOUT,9120)
      STOP
C
 1070 DO 1080 I=1,NJV
      READ (NIN,8020) ELNUM, JVAL(ELNUM)
      WRITE (NOUT,9040) ELNUM, JVAL(ELNUM)
 1080 CONTINUE
C
C                            CONSTRUCTION OF NODAL FREEDOM ARRAY NF
C
 1090 TOTDOF = 0
      DO 1110 I=1,TOTNOD
      DO 1100 J=1,DOFNOD
      TOTDOF = TOTDOF + 1
      NF(I,J) = TOTDOF
 1100 CONTINUE
 1110 CONTINUE
C
C                            CALCULATION OF SEMI-BANDWIDTH
C
      CALL BNDWTH(ELTOP, IELTOP, JELTOP, NF, INF, JNF, DOFNOD,
     *     TOTELS, HBAND, ITEST)
      IF (HBAND.LE.JSYSK) GO TO 1120
      WRITE (NOUT,9130) HBAND, JSYSK
      STOP
C
C*                           ************************************
C*                           *                                  *
C*                           * SYSTEM STIFFNESS MATRIX ASSEMBLY *
C*                           *                                  *
C*                           ************************************
C
 1120 CALL CVCNUL(RHS, IRHS, TOTDOF, ITEST)
      CALL CMTNUL(SYSK, ISYSK, JSYSK, TOTDOF, HBAND, ITEST)
      DOFEL = NODEL*DOFNOD
      CALL QQUA4(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)
C
C                             MAIN MATRIX ASSEMBLY
C
      DO 1170 NELE=1,TOTELS
      CALL ELGEOM(NELE, ELTOP, IELTOP, JELTOP, COORD, ICOORD,
     *     JCOORD, GEOM, IGEOM, JGEOM, DIMEN, ITEST)
C
      ELTYP = ELTOP(NELE,1)
      CONST = 0.0D0
      IF (ELTYP.EQ.2) CONST = OMEGA*COND
      CALL MATNUL(S, IS, JS, DIMEN, DIMEN, ITEST)
      DO 1130 I=1,DIMEN
      S(I,I) = 1.0D0/MU0
      IF (ELTYP.EQ.2) S(I,I) = 1.0D0/(PERM*MU0)
 1130 CONTINUE
C
C                            INTEGRATION LOOP FOR ELEMENT STIFFNESS
C                            USING NQP QUADRATURE POINTS
C
      CALL MATNUL(ELK, IELK, JELK, DOFEL, DOFEL, ITEST)
      CALL MATNUL(ELM, IELM, JELM, DOFEL, DOFEL, ITEST)
      CALL VECNUL(ELJ, IELJ, DOFEL, ITEST)
      DO 1160 IQUAD=1,NQP
C
C                            FORM LINEAR SHAPE FUNCTION AND SPACE
C                            DERIVATIVES IN THE LOCAL CORRDINATES.
C
      XI = ABSS(1,IQUAD)
      ETA = ABSS(2,IQUAD)
      CALL QUAM4(FUN, IFUN, LDER, ILDER, JLDER, XI, ETA, ITEST)
      CALL MATMUL(LDER, ILDER, JLDER, GEOM, IGEOM, JGEOM, JAC,
     *     IJAC, JJAC, DIMEN, NODEL, DIMEN, ITEST)
C
C                            TRANSFORM LOCAL DERIVATIVES TO GLOBAL
C                            COORDINATE SYSTEM AND OUTPUT TO WORK
C                            FILE FOR LATER USE.
C
      CALL MATINV(JAC, IJAC, JJAC, JACIN, IJACIN, JJACIN, DIMEN,
     *     DET, ITEST)
      CALL MATMUL(JACIN, IJACIN, JJACIN, LDER, ILDER, JLDER, GDER,
     *     IGDER, JGDER, DIMEN, DIMEN, NODEL, ITEST)
      WRITE (NTEMP) ((GDER(I,J),I=1,DIMEN),J=1,DOFEL)
C
C                            FORMATION OF ELEMENT STIFFNESS ELK AND ELM
C
      CALL MATMUL(S, IS, JS, GDER, IGDER, JGDER, PD, IPD, JPD,
     *     DIMEN, DIMEN, DOFEL, ITEST)
      CALL MATRAN(GDER, IGDER, JGDER, GDERT, IGDERT, JGDERT,
     *     DIMEN, DOFEL, ITEST)
      CALL MATMUL(GDERT, IGDERT, JGDERT, PD, IPD, JPD, DTPD,
     *     IDTPD, JDTPD, DOFEL, DIMEN, DOFEL, ITEST)
      CALL DYAD(FUN, IFUN, FUN, IFUN, FTF, IFTF, JFTF, DOFEL,
     *     ITEST)
C
      QUOT = DABS(DET)*WGHT(IQUAD)
      DO 1150 I=1,DOFEL
      DO 1140 J=1,DOFEL
      DTPD(I,J) = DTPD(I,J)*QUOT
      FTF(I,J) = FTF(I,J)*CONST*QUOT
 1140 CONTINUE
      SCVEC(I) = FUN(I)*JVAL(NELE)*QUOT
 1150 CONTINUE
      CALL MATADD(ELK, IELK, JELK, DTPD, IDTPD, JDTPD, DOFEL,
     *     DOFEL, ITEST)
      CALL MATADD(ELM, IELM, JELM, FTF, IFTF, JFTF, DOFEL, DOFEL,
     *     ITEST)
      CALL VECADD(ELJ, IELJ, SCVEC, ISCVEC, DOFEL, ITEST)
 1160 CONTINUE
C
C                            ASSEMBLY OF SYSTEM STIFFNESS
C                            MATRICES ELK AND ELM
C
      CALL DIRECT(NELE, ELTOP, IELTOP, JELTOP, NF, INF, JNF,
     *     DOFNOD, STEER, ISTEER, ITEST)
      CALL RASSYM(SYSK, ISYSK, JSYSK, ELK, IELK, JELK, STEER,
     *     ISTEER, HBAND, DOFEL, ITEST)
      CALL IASSYM(SYSK, ISYSK, JSYSK, ELM, IELM, JELM, STEER,
     *     ISTEER, HBAND, DOFEL, ITEST)
      CALL RASRHS(RHS, IRHS, ELJ, IELJ, STEER, ISTEER, DOFEL,
     *     ITEST)
 1170 CONTINUE
C
C*                           *********************
C*                           *                   *
C*                           * EQUATION SOLUTION *
C*                           *                   *
C*                           *********************
C
C                            IMPLEMENT BOUNDARY CONDITIONS
C
      DO 1180 I=1,BNDNOD
      J = BNODE(I)
      SYSK(1,J,HBAND) = SYSK(1,J,HBAND)*SCALE
      SYSK(2,J,HBAND) = SYSK(2,J,HBAND)*SCALE
      RHS(1,J) = SYSK(1,J,HBAND)*BVAL(1,I) - SYSK(2,J,HBAND)*
     *     BVAL(2,I)
      RHS(2,J) = SYSK(2,J,HBAND)*BVAL(1,I) + SYSK(1,J,HBAND)*
     *     BVAL(2,I)
 1180 CONTINUE
C
C                            SOLUTION OF SYSTEM MATRIX FOR THE
C                            NODAL VALUES OF THE POTENTIAL
C
      CALL CSYSOL(SYSK, ISYSK, JSYSK, RHS, IRHS, TOTDOF, HBAND,
     *     ITEST)
      WRITE (NOUT,9140)
      CALL CPRTVL(RHS, IRHS, NF, INF, JNF, DOFNOD, TOTNOD, NOUT,
     *     ITEST)
C
C*                           ****************************
C*                           *                          *
C*                           *  RECOVER FIELD VALUES    *
C*                           *  FROM VECTOR POTENTIAL   *
C*                           *                          *
C*                           ****************************
C
      REWIND NTEMP
      DO 1210 NELE=1,TOTELS
      WRITE (NOUT,9150) NELE
C
C                            SELECT NODAL POTENTIALS FOR ELEMENT
C                            NELE USING THE STEERING VECTOR
C
      CALL DIRECT(NELE, ELTOP, IELTOP, JELTOP, NF, INF, JNF,
     *     DOFNOD, STEER, ISTEER, ITEST)
      DO 1200 IQUAD=1,NQP
      CALL VECNUL(ELPOT, IELPOT, DOFEL, ITEST)
      DO 1190 I=1,DOFEL
      JSTEER = STEER(I)
      IF (JSTEER.NE.0) ELPOT(I) = DSQRT(RHS(1,JSTEER)**2+RHS(2,
     *     JSTEER)**2)
 1190 CONTINUE
C
C                            RECOVER GDER FOR POINT IQUAD
C
      READ (NTEMP) ((GDER(I,J),I=1,DIMEN),J=1,DOFEL)
      CALL MATVEC(GDER, IGDER, JGDER, ELPOT, IELPOT, DIMEN, DOFEL,
     *     FIELD, IFIELD, ITEST)
      BMOD = DSQRT(FIELD(1)**2+FIELD(2)**2)
      WRITE (NOUT,9160) IQUAD
      WRITE (NOUT,9170) (FIELD(I),I=1,2), BMOD
 1200 CONTINUE
 1210 CONTINUE
C
C                             END OF PROGRAM
C
      STOP
C
 8010 FORMAT (16I5)
 8020 FORMAT (I5, 6F10.0)
 8030 FORMAT (3F10.0)
 9010 FORMAT (//25H **** NODAL GEOMETRY ****/1H )
 9020 FORMAT (1H , 16I5)
 9030 FORMAT (44H *** ERROR : TOTNOD.LE.0 OR TOTNOD.GT.ICOORD)
 9040 FORMAT (1H , I5, 6F10.5)
 9050 FORMAT (//27H **** ELEMENT TOPOLOGY ****/1H )
 9060 FORMAT (44H *** ERROR : TOTELS.LE.0 OR TOTELS.GT.IELTOP)
 9070 FORMAT (//31H ****  FREQ.  COND.  PERM. ****/1H )
 9080 FORMAT (1H , 3F10.3)
 9090 FORMAT (//30H **** BOUNDARY CONDITIONS ****/1H )
 9100 FORMAT (44H *** ERROR : BNDNOD.LE.0 OR BNDNOD.GT.IBNODE)
 9110 FORMAT (//27H **** CURRENT SOURCES **** /1H )
 9120 FORMAT (25H *** ERROR : NJV.GT.IJVAL)
 9130 FORMAT (22H*** ERROR - HBAND WAS , I5, 16H, JSYSK IS ONLY ,
     *     I5)
 9140 FORMAT (//27H **** NODAL POTENTIALS ****/1H )
 9150 FORMAT (//31H **** FIELD VALUES FOR ELEMENT , I3, 5H ****/
     *     1H , 5X, 2HBX, 12X, 2HBY, 14X, 4HMODB)
 9160 FORMAT (/18H QUADRATURE POINT , I3/1H )
 9170 FORMAT (1H , 2(D12.5, 2X), 3X, D12.5)
      END
C***********************************************************************
