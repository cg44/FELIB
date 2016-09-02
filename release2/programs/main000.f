C***********************************************************************
C$SPLIT$SEG1P1$*********************************************************
C***********************************************************************
      INTEGER DIF, DIMEN, DOFEL, DOFNOD, ELNUM, ELTOP, ELTYP,
     *     HBAND, I, IABSS, IB, IBT, IBTDB, ICOORD, ID, IDB,
     *     IELDIS, IELK, IELTOP, IFUN, IGDER, IGEOM, IJAC,
     *     IJACIN, ILDER, ILOADS, INF, IQUAD, IRESTR, ISTEER,
     *     ISTRN, ISTRS, ISYSK, ITEST, IWGHT, J, JABSS, JB,
     *     JBT, JBTDB, JCOORD, JD, JDB, JELK, JELTOP, JGDER,
     *     JGEOM, JJAC, JJACIN, JLDER, JNF, JRESTR, JSTEER,
     *     JSYSK, K, LODNOD, NELE, NF, NIN, NODEL, NODNUM,
     *     NOUT, NQP, NTEMP, NUMSS, RESNOD, RESTR, STEER, TOTDOF,
     *     TOTELS, TOTNOD
      DOUBLE PRECISION ABSS, B, BT, BTDB, COORD, D, DB, DET,
     *     E, ELDIS, ELK, ETA, FUN, GDER, GEOM, JAC, JACIN,
     *     LDER, LOADS, NU, QUOT, STRN, STRS, SYSK, WGHT, WORK,
     *     XI
      LOGICAL FIRST
      DIMENSION ABSS(3,9), B(6,24), BT(24,6), BTDB(24,24),
     *     D(6,6), DB(6,24), ELDIS(24), ELK(24,24), FUN(8),
     *     GDER(3,8), GEOM(8,3), JAC(3,3), JACIN(3,3),
     *     LDER(3,8), STEER(24), STRN(6), STRS(6), WGHT(9),
     *     WORK(3)
C
C                            PROBLEM SIZE DEPENDENT ARRAYS
C
      DIMENSION COORD(100,3), ELTOP(100,10), LOADS(100),
     *     NF(100,3), RESTR(100,4), SYSK(100,25)
C
      DATA IABSS /3/, IB /6/, IBT /24/, IBTDB /24/, ID /6/,
     *     IDB /6/, IELDIS /24/, IELK /24/, IFUN /8/, IGDER /3/,
     *     IGEOM /8/, IJAC /3/, IJACIN /3/, ILDER /3/,
     *     ISTEER /24/, ISTRN /6/, ISTRS /6/, IWGHT /9/,
     *     JABSS /9/, JB /24/, JBT /6/, JBTDB /24/, JCOORD /3/,
     *     JD /6/, JDB /24/, JELK /24/, JGDER /8/, JGEOM /3/,
     *     JJAC /3/, JJACIN /3/, JLDER /8/, JNF /3/, JRESTR /4/
C
C                            PROBLEM SIZE DEPENDENT DATA STATEMENTS
C
      DATA ICOORD /100/, IELTOP /100/, ILOADS /100/, INF /100/,
     *     IRESTR /100/, ISYSK /100/, JELTOP /10/, JSYSK /25/
C
      DATA NIN /5/, NOUT /6/, NTEMP /7/
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
C                            INPUT OF MATERIAL PROPERTIES AND
C                            CONSTRUCTION OF STRESS-STRAIN MATRIX D
C                            FOR PLANE STRAIN
C
      WRITE (NOUT,9050)
      READ (NIN,8030) NU, E
      WRITE (NOUT,9060) NU, E
      CALL DPSN(D, ID, JD, E, NU, NUMSS, ITEST)
C
C                            INPUT OF NUMBER OF DEGREES OF FREEDOM
C                            PER NODE, INPUT OF RESTRAINED NODE
C                            DATA AND CONSTRUCTION OF NODAL FREEDOM
C                            ARRAY NF
C
      WRITE (NOUT,9070)
      READ (NIN,8010) DOFNOD
      WRITE (NOUT,9020) DOFNOD
      READ (NIN,8010) RESNOD
      WRITE (NOUT,9020) RESNOD
      K = DOFNOD + 1
      DO 1030 I=1,RESNOD
      READ (NIN,8010) (RESTR(I,J),J=1,K)
      WRITE (NOUT,9020) (RESTR(I,J),J=1,K)
 1030 CONTINUE
      CALL FORMNF(RESTR, IRESTR, JRESTR, RESNOD, TOTNOD, DOFNOD,
     *     NF, INF, JNF, TOTDOF, ITEST)
C
C                            LOADING DATA INPUT
C
      WRITE (NOUT,9080)
      CALL VECNUL(LOADS, ILOADS, TOTDOF, ITEST)
      READ (NIN,8010) LODNOD
      WRITE (NOUT,9020) LODNOD
      DO 1050 I=1,LODNOD
      READ (NIN,8020) NODNUM, (WORK(J),J=1,DOFNOD)
      WRITE (NOUT,9030) NODNUM, (WORK(J),J=1,DOFNOD)
      DO 1040 J=1,DOFNOD
      K = NF(NODNUM,J)
      IF (K.EQ.0) GO TO 1040
      LOADS(K) = WORK(J)
 1040 CONTINUE
 1050 CONTINUE
C
C                            CALCULATION OF SEMI-BANDWIDTH
C
      FIRST = .TRUE.
      DO 1060 NELE=1,TOTELS
      CALL FREDIF(NELE, ELTOP, IELTOP, JELTOP, NF, INF, JNF,
     *     DOFNOD, FIRST, DIF, ITEST)
 1060 CONTINUE
      HBAND = DIF + 1
C
C*                           ************************************
C*                           *                                  *
C*                           * GLOBAL STIFFNESS MATRIX ASSEMBLY *
C*                           *                                  *
C*                           ************************************
C
      CALL MATNUL(SYSK, ISYSK, JSYSK, TOTDOF, HBAND, ITEST)
      DOFEL = NODEL*DOFNOD
      CALL QQUA4(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)
      DO 1100 NELE=1,TOTELS
      CALL ELGEOM(NELE, ELTOP, IELTOP, JELTOP, COORD, ICOORD,
     *     JCOORD, GEOM, IGEOM, JGEOM, DIMEN, ITEST)
C
C                            INTEGRATION LOOP FOR ELEMENT STIFFNESS
C                            USING NQP QUADRATURE POINTS
C
      CALL MATNUL(ELK, IELK, JELK, DOFEL, DOFEL, ITEST)
      DO 1090 IQUAD=1,NQP
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
C                            FORMATION OF STRAIN-DISPLACEMENT
C                            MATRIX B AND OUTPUT TO WORK FILE FOR
C                            LATER RECOVERY PROCESS
C
      CALL B2C2(B, IB, JB, GDER, IGDER, JGDER, NODEL, ITEST)
      WRITE (NTEMP) ((B(I,J),I=1,NUMSS),J=1,DOFEL)
C
C                            FORMATION OF ELEMENT STIFFNESS ELK
C
      CALL MATMUL(D, ID, JD, B, IB, JB, DB, IDB, JDB, NUMSS,
     *     NUMSS, DOFEL, ITEST)
      CALL MATRAN(B, IB, JB, BT, IBT, JBT, NUMSS, DOFEL, ITEST)
      CALL MATMUL(BT, IBT, JBT, DB, IDB, JDB, BTDB, IBTDB, JBTDB,
     *     DOFEL, NUMSS, DOFEL, ITEST)
      QUOT = DABS(DET)*WGHT(IQUAD)
      DO 1080 I=1,DOFEL
      DO 1070 J=1,DOFEL
      BTDB(I,J) = BTDB(I,J)*QUOT
 1070 CONTINUE
 1080 CONTINUE
      CALL MATADD(ELK, IELK, JELK, BTDB, IBTDB, JBTDB, DOFEL,
     *     DOFEL, ITEST)
 1090 CONTINUE
C
C                            ASSEMBLY OF SYSTEM STIFFNESS MATRIX
C                            SYSK
C
      CALL DIRECT(NELE, ELTOP, IELTOP, JELTOP, NF, INF, JNF,
     *     DOFNOD, STEER, ISTEER, ITEST)
      CALL ASSYM(SYSK, ISYSK, JSYSK, ELK, IELK, JELK, STEER,
     *     ISTEER, HBAND, DOFEL, ITEST)
 1100 CONTINUE
C
C*                           *********************
C*                           *                   *
C*                           * EQUATION SOLUTION *
C*                           *                   *
C*                           *********************
C
C                            SOLUTION OF SYSTEM MATRIX FOR THE
C                            NODAL DISPLACEMENTS
C
      WRITE (NOUT,9090)
      CALL PRTVAL(LOADS, ILOADS, NF, INF, JNF, DOFNOD, TOTNOD,
     *     NOUT, ITEST)
      CALL CHOSOL(SYSK, ISYSK, JSYSK, LOADS, ILOADS, TOTDOF,
     *     HBAND, ITEST)
      WRITE (NOUT,9100)
      CALL PRTVAL(LOADS, ILOADS, NF, INF, JNF, DOFNOD, TOTNOD,
     *     NOUT, ITEST)
C
C*                           **************************
C*                           *                        *
C*                           * STRESS-STRAIN RECOVERY *
C*                           *                        *
C*                           **************************
C
C                            RECOVERY OF STRESSES AND STRAINS AT
C                            THE QUADRATURE SAMPLING POINTS USING
C                            THE ELEMENT STRAIN-DISPLACEMENT MATRIX
C                            B FROM THE WORK FILE
C
      REWIND NTEMP
      DO 1130 NELE=1,TOTELS
      WRITE (NOUT,9110) NELE
C
C                            SELECT NODAL DISPLACEMENTS FOR ELEMENT
C                            NELE USING THE STEERING VECTOR
C
      CALL DIRECT(NELE, ELTOP, IELTOP, JELTOP, NF, INF, JNF,
     *     DOFNOD, STEER, ISTEER, ITEST)
      DO 1120 IQUAD=1,NQP
      CALL VECNUL(ELDIS, IELDIS, DOFEL, ITEST)
      DO 1110 I=1,DOFEL
      JSTEER = STEER(I)
      IF (JSTEER.NE.0) ELDIS(I) = LOADS(JSTEER)
 1110 CONTINUE
C
C                            RECOVER B MATRIX FOR POINT IQUAD AND
C                            CALCULATE THE STRESSES AND STRAINS
C
      READ (NTEMP) ((B(I,J),I=1,NUMSS),J=1,DOFEL)
      CALL MATVEC(B, IB, JB, ELDIS, IELDIS, NUMSS, DOFEL, STRN,
     *     ISTRN, ITEST)
      CALL MATVEC(D, ID, JD, STRN, ISTRN, NUMSS, NUMSS, STRS,
     *     ISTRS, ITEST)
      WRITE (NOUT,9120) IQUAD
      CALL PRTVEC(STRN, ISTRN, NUMSS, NOUT, ITEST)
      CALL PRTVEC(STRS, ISTRS, NUMSS, NOUT, ITEST)
 1120 CONTINUE
 1130 CONTINUE
      STOP
 8010 FORMAT (16I5)
 8020 FORMAT (I5, 6F10.0)
 8030 FORMAT (2F10.0)
 9010 FORMAT (//25H **** NODAL GEOMETRY ****//1H )
 9020 FORMAT (1H , 16I5)
 9030 FORMAT (1H , I5, 6F10.5)
 9040 FORMAT (//27H **** ELEMENT TOPOLOGY ****//1H )
 9050 FORMAT (//30H **** MATERIAL PROPERTIES ****//1H )
 9060 FORMAT (1H , F10.5, E10.3)
 9070 FORMAT (//25H **** RESTRAINT DATA ****//1H )
 9080 FORMAT (//23H **** LOADING DATA ****//1H )
 9090 FORMAT (//26H **** VECTOR OF LOADS ****//1H )
 9100 FORMAT (//42H **** EQUILIBRIUM DISPLACEMENTS (U,V) ****//1H )
 9110 FORMAT (//39H **** STRAINS AND STRESSES FOR ELEMENT , I3,
     *     5H ****//1H )
 9120 FORMAT (/18H QUADRATURE POINT , I3/1H )
      END
