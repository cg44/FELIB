C***********************************************************************
C$SPLIT$SEG2P1$*********************************************************
C***********************************************************************
      INTEGER DIMEN, DOFEL, DOFNOD, ELNUM, ELTOP, ELTYP, I,
     *     IABSS, IB, IBT, IBTDB, ICOORD, ID, IDB, IDIAG, IELK,
     *     IELM, IELTOP, IFTF, IFUN, IGDER, IGEOM, IJAC, IJACIN,
     *     ILDER, INF, IQUAD, IRESTR, ISHP, ISTEER, ISUB, ISYSK,
     *     ISYSM, ITEST, ITSHP, IWGHT, J, JABSS, JB, JBT, JBTDB,
     *     JCOORD, JD, JDB, JELK, JELM, JELTOP, JFTF, JGDER,
     *     JGEOM, JJAC, JJACIN, JLDER, JNF, JRESTR, JSHP, JSYSK,
     *     JTSHP, K, NELE, NF, NIN, NMODES, NODEL, NODNUM,
     *     NOUT, NQP, NUMSS, RESNOD, RESTR, STEER, TOTDOF,
     *     TOTELS, TOTNOD
      DOUBLE PRECISION ABSS, AREA, B, BT, BTDB, COORD, D, DB,
     *     DET, DIAG, E, ELK, ELM, EPS, ETA, FTF, FUN, GDER,
     *     GEOM, JAC, JACIN, LDER, NU, QUOT, RHO, SHP, SUB,
     *     SYSK, SYSM, TOL, TSHP, VEPS, VTOL, WGHT, X, XI
      DIMENSION ABSS(3,9), B(6,24), BT(24,6), BTDB(24,24),
     *     D(6,6), DB(6,24), ELK(24,24), ELM(24,24),
     *     FTF(24,24), FUN(8), GDER(3,8), GEOM(8,3), JAC(3,3),
     *     JACIN(3,3), LDER(3,8), SHP(3,24), STEER(24),
     *     TSHP(24,3), WGHT(9)
C
C                            PROBLEM SIZE DEPENDENT ARRAYS
C
      DIMENSION COORD(100,3), DIAG(100), ELTOP(100,10),
     *     NF(100,3), RESTR(50,4), SUB(100), SYSK(100,100),
     *     SYSM(100)
C
      DATA IABSS /3/, IB /6/, IBT /24/, IBTDB /24/, ID /6/,
     *     IDB /6/, IELK /24/, IELM /24/, IFTF /24/, IFUN /8/,
     *     IGDER /3/, IGEOM /8/, IJAC /3/, IJACIN /3/, ILDER /3/,
     *     ISHP /3/, ISTEER /24/, ITSHP /24/, IWGHT /9/,
     *     JABSS /9/, JB /24/, JBT /6/, JBTDB /8/, JCOORD /3/,
     *     JD /6/, JDB /24/, JELK /24/, JELM /24/, JFTF /24/,
     *     JGDER /8/, JGEOM /3/, JJAC /3/, JJACIN /3/, JLDER /8/,
     *     JNF /3/, JRESTR /4/, JSHP /24/, JTSHP /3/
C
C                            PROBLEM SIZE DEPENDENT DATA STATEMENTS
C
      DATA ICOORD /100/, IDIAG /100/, IELTOP /100/, INF /100/,
     *     IRESTR /50/, ISUB /100/, ISYSK /100/, ISYSM /100/,
     *     JELTOP /10/, JSYSK /100/
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
C                            INPUT OF MATERIAL PROPERTIES AND
C                            CONSTRUCTION OF STRESS-STRAIN MATRIX D
C                            FOR PLANE STRAIN
C
      WRITE (NOUT,9050)
      READ (NIN,8030) NU, E, RHO
      WRITE (NOUT,9060) NU, E, RHO
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
      READ (NIN,8010) NMODES
C
C*                           *******************************
C*                           *                             *
C*                           * SYSTEM STIFFNESS MATRIX AND *
C*                           *     MASS MATRIX ASSEMBLY    *
C*                           *                             *
C*                           *******************************
C
      CALL MATNUL(SYSK, ISYSK, JSYSK, TOTDOF, TOTDOF, ITEST)
      CALL VECNUL(SYSM, ISYSM, TOTDOF, ITEST)
      DOFEL = NODEL*DOFNOD
      CALL QQUA4(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)
      DO 1070 NELE=1,TOTELS
      AREA = 0.0D0
      CALL ELGEOM(NELE, ELTOP, IELTOP, JELTOP, COORD, ICOORD,
     *     JCOORD, GEOM, IGEOM, JGEOM, DIMEN, ITEST)
C
C                            INTEGRATION LOOP FOR ELEMENT STIFFNESS
C                            AND ELEMENT LUMPED MASS USING NQP
C                            QUADRATURE POINTS
C
      CALL MATNUL(ELK, IELK, JELK, DOFEL, DOFEL, ITEST)
      CALL MATNUL(ELM, IELM, JELM, DOFEL, DOFEL, ITEST)
      DO 1060 IQUAD=1,NQP
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
C                            MATRIX B AND FORMATION OF INTEGRAND
C                            FOR ELEMENT STIFFNESS MATRIX ELK
C
      CALL B2C2(B, IB, JB, GDER, IGDER, JGDER, NODEL, ITEST)
      CALL MATMUL(D, ID, JD, B, IB, JB, DB, IDB, JDB, NUMSS,
     *     NUMSS, DOFEL, ITEST)
      CALL MATRAN(B, IB, JB, BT, IBT, JBT, NUMSS, DOFEL, ITEST)
      CALL MATMUL(BT, IBT, JBT, DB, IDB, JDB, BTDB, IBTDB, JBTDB,
     *     DOFEL, NUMSS, DOFEL, ITEST)
C
C                            FORMATION OF INTEGRAND FOR ELEMENT
C                            MASS MATRIX ELM
C
      CALL SHAPFN(SHP, ISHP, JSHP, FUN, IFUN, NODEL, DOFNOD, ITEST)
      CALL MATRAN(SHP, ISHP, JSHP, TSHP, ITSHP, JTSHP, DOFNOD,
     *     DOFEL, ITEST)
      CALL MATMUL(TSHP, ITSHP, JTSHP, SHP, ISHP, JSHP, FTF, IFTF,
     *     JFTF, DOFEL, DOFNOD, DOFEL, ITEST)
C
      QUOT = DABS(DET)*WGHT(IQUAD)
      AREA = AREA + QUOT
      DO 1050 I=1,DOFEL
      DO 1040 J=1,DOFEL
      BTDB(I,J) = BTDB(I,J)*QUOT
      FTF(I,J) = FTF(I,J)*QUOT*RHO
 1040 CONTINUE
 1050 CONTINUE
      CALL MATADD(ELK, IELK, JELK, BTDB, IBTDB, JBTDB, DOFEL,
     *     DOFEL, ITEST)
      CALL MATADD(ELM, IELM, JELM, FTF, IFTF, JFTF, DOFEL, DOFEL,
     *     ITEST)
 1060 CONTINUE
C
C                            ASSEMBLY OF SYSTEM STIFFNESS MATRIX
C                            SYSK, SYSTEM LUMPED MASS MATRIX SYSM
C
      CALL DIRECT(NELE, ELTOP, IELTOP, JELTOP, NF, INF, JNF,
     *     DOFNOD, STEER, ISTEER, ITEST)
      CALL ASFUL(SYSK, ISYSK, JSYSK, ELK, IELK, JELK, STEER,
     *     ISTEER, DOFEL, ITEST)
      CALL ASLMS(SYSM, ISYSM, ELM, IELM, JELM, STEER, ISTEER,
     *     DOFEL, DOFNOD, AREA, ITEST)
 1070 CONTINUE
C
C*                           *********************
C*                           *                   *
C*                           * EQUATION SOLUTION *
C*                           *                   *
C*                           *********************
C
C
C                            REDUCTION OF SYSTEM MASS MATRIX BY
C                            CHOLESKI REDUCTION AND CONSTRUCTION OF
C                            MODIFIED SYSTEM STIFFNESS MATRIX
C
      DO 1080 I=1,TOTDOF
      SYSM(I) = 1.0D0/DSQRT(SYSM(I))
 1080 CONTINUE
      DO 1100 I=1,TOTDOF
      DO 1090 J=1,TOTDOF
      SYSK(I,J) = SYSM(I)*SYSK(I,J)*SYSM(J)
 1090 CONTINUE
 1100 CONTINUE
C
C                            HOUSEHOLDER TRIDIAGONALISATION OF SYSK
C                            AND EIGENVALUE AND EIGENVECTOR
C                            RECOVERY USING QL TRANSFORMATIONS
C
      TOL = VTOL(X)
      CALL HOUSE(SYSK, ISYSK, JSYSK, SYSK, ISYSK, JSYSK, DIAG,
     *     IDIAG, SUB, ISUB, TOTDOF, TOL, ITEST)
      EPS = VEPS(X)
      CALL QLVEC(DIAG, IDIAG, SUB, ISUB, SYSK, ISYSK, JSYSK,
     *     TOTDOF, EPS, ITEST)
      WRITE (NOUT,9080)
      CALL PRTVEC(DIAG, IDIAG, TOTDOF, NOUT, ITEST)
      DO 1120 J=1,NMODES
      DO 1110 I=1,TOTDOF
      DIAG(I) = SYSK(I,J)*SYSM(I)
 1110 CONTINUE
      WRITE (NOUT,9090) J
      CALL PRTVEC(DIAG, IDIAG, TOTDOF, NOUT, ITEST)
 1120 CONTINUE
      STOP
 8010 FORMAT (16I5)
 8020 FORMAT (I5, 6F10.0)
 8030 FORMAT (6F10.0)
 9010 FORMAT (//25H **** NODAL GEOMETRY ****//1H )
 9020 FORMAT (1H , 16I5)
 9030 FORMAT (1H , I5, 6F10.5)
 9040 FORMAT (//27H **** ELEMENT TOPOLOGY ****//1H )
 9050 FORMAT (//30H **** MATERIAL PROPERTIES ****//1H )
 9060 FORMAT (1H , F10.5, E10.3, F10.5)
 9070 FORMAT (//25H **** RESTRAINT DATA ****//1H )
 9080 FORMAT (//22H **** FREQUENCIES ****/1H )
 9090 FORMAT (/21H MODE SHAPE FOR MODE , I5/1H )
      END
