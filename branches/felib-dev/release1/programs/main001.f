C***********************************************************************
C$SPLIT$SEG2P1$*********************************************************
C***********************************************************************
      INTEGER DIMEN, DOFEL, DOFNOD, ELNUM, ELTOP, ELTYP, I,                    1
     *     IABSS, IB, IBT, IBTDB, ICOORD, ID, IDB, IDER, IDERIV,               1
     *     IDIAG, IELK, IELM, IELTOP, IFTF, IFUN, IGEOM, IJAC,                 1
     *     IJACIN, INF, IQUAD, IRESTR, ISHP, ISTEER, ISUB,                     1
     *     ISYSK, ISYSM, ITEST, ITSHP, IWGHT, J, JABSS, JB,                    1
     *     JBT, JBTDB, JCOORD, JD, JDB, JDER, JDERIV, JELK,                    1
     *     JELM, JELTOP, JFTF, JGEOM, JJAC, JJACIN, JNF, JRESTR,               1
     *     JSHP, JSYSK, JTSHP, K, NELE, NF, NIN, NMODES, NODEL,                1
     *     NODNUM, NOUT, NQP, NUMSS, NWORK, RESNOD, RESTR,                     1
     *     STEER, TOTDOF, TOTELS, TOTNOD                                       1
      DOUBLE PRECISION ABSS, AREA, B, BT, BTDB, COORD, D, DB,                  2
     *     DER, DERIV, DET, DIAG, E, ELK, ELM, EPS, ETA, FTF,                  2
     *     FUN, GEOM, JAC, JACIN, NU, QUOT, RHO, SHP, SUB,                     2
     *     SYSK, SYSM, TOL, TSHP, VEPS, VTOL, WGHT, WORK, X,                   2
     *     XI                                                                  2
      DIMENSION ABSS(3,9), B(6,24), BT(24,6), BTDB(24,24),                     3
     *     COORD(100,3), D(6,6), DB(6,24), DER(3,8), DERIV(3,8),               3
     *     DIAG(100), ELK(24,24), ELM(24,24), ELTOP(100,10),                   3
     *     FTF(24,24), FUN(8), GEOM(8,3), JAC(3,3), JACIN(3,3),                3
     *     NF(100,3), NWORK(8), RESTR(50,4), SHP(3,24),                        3
     *     STEER(24), SUB(100), SYSK(100,100), SYSM(100),                      3
     *     TSHP(24,3), WGHT(9), WORK(3)                                        3
      DATA IABSS /3/, IB /6/, IBT /24/, IBTDB /24/, ICOORD /100/,              4
     *     ID /6/, IDB /6/, IDER /3/, IDERIV /3/, IDIAG /100/,                 4
     *     IELK /24/, IELM /24/, IELTOP /100/, IFTF /24/,                      4
     *     IFUN /8/, IGEOM /8/, IJAC /3/, IJACIN /3/,                          4
     *     INF /100/, IRESTR /50/, ISHP /3/, ISTEER /24/,                      4
     *     ISUB /100/, ISYSK /100/, ISYSM /100/, ITSHP /24/,                   4
     *     IWGHT /9/, JABSS /9/, JB /24/, JBT /6/, JBTDB /8/,                  4
     *     JCOORD /3/, JD /6/, JDB /24/, JDER /8/, JDERIV /8/,                 4
     *     JELK /24/, JELM /24/, JELTOP /10/, JFTF /24/,                       4
     *     JGEOM /3/, JJAC /3/, JJACIN /3/, JNF /3/, JRESTR /4/,               4
     *     JSHP /24/, JSYSK /100/, JTSHP /3/                                   4
C
      DATA NIN /5/, NOUT /6/                                                   5
C
C                            SET ITEST FOR FULL CHECKING
C
      ITEST = 0                                                                6
C
C*                           **********************
C*                           *                    *
C*                           * INPUT DATA SECTION *
C*                           *                    *
C*                           **********************
C
C                            INPUT OF NODAL GEOMETRY
C
      WRITE (NOUT,9010)                                                        7
      READ (NIN,8010) TOTNOD, DIMEN                                            8
      WRITE (NOUT,9020) TOTNOD, DIMEN                                          9
      DO 1020 I=1,TOTNOD                                                      10
      READ (NIN,8020) NODNUM, (WORK(J),J=1,DIMEN)                             11
      WRITE (NOUT,9030) NODNUM, (WORK(J),J=1,DIMEN)                           12
      DO 1010 J=1,DIMEN                                                       13
      COORD(NODNUM,J) = WORK(J)                                               14
 1010 CONTINUE                                                                15
 1020 CONTINUE                                                                16
C
C                            INPUT OF ELEMENT TOPOLOGY
C
      WRITE (NOUT,9040)                                                       17
      READ (NIN,8010) ELTYP, TOTELS, NODEL                                    18
      WRITE (NOUT,9020) ELTYP, TOTELS, NODEL                                  19
      DO 1040 I=1,TOTELS                                                      20
      READ (NIN,8010) ELNUM, (NWORK(J),J=1,NODEL)                             21
      WRITE (NOUT,9020) ELNUM, (NWORK(J),J=1,NODEL)                           22
      ELTOP(ELNUM,1) = ELTYP                                                  23
      ELTOP(ELNUM,2) = NODEL                                                  24
      DO 1030 J=1,NODEL                                                       25
      ELTOP(ELNUM,J+2) = NWORK(J)                                             26
 1030 CONTINUE                                                                27
 1040 CONTINUE                                                                28
C
C                            INPUT OF MATERIAL PROPERTIES AND
C                            CONSTRUCTION OF STRESS-STRAIN MATRIX D
C                            FOR PLANE STRAIN
C
      WRITE (NOUT,9050)                                                       29
      READ (NIN,8030) NU, E, RHO                                              30
      WRITE (NOUT,9060) NU, E, RHO                                            31
      CALL DPSN(D, ID, JD, E, NU, NUMSS, ITEST)                               32
C
C                            INPUT OF NUMBER OF DEGREES OF FREEDOM
C                            PER NODE, INPUT OF RESTRAINED NODE
C                            DATA AND CONSTRUCTION OF NODAL FREEDOM
C                            ARRAY NF
C
      WRITE (NOUT,9070)                                                       33
      READ (NIN,8010) DOFNOD                                                  34
      WRITE (NOUT,9020) DOFNOD                                                35
      READ (NIN,8010) RESNOD                                                  36
      WRITE (NOUT,9020) RESNOD                                                37
      K = DOFNOD + 1                                                          38
      DO 1050 I=1,RESNOD                                                      39
      READ (NIN,8010) (RESTR(I,J),J=1,K)                                      40
      WRITE (NOUT,9020) (RESTR(I,J),J=1,K)                                    41
 1050 CONTINUE                                                                42
      CALL FORMNF(RESTR, IRESTR, JRESTR, RESNOD, TOTNOD, DOFNOD,              43
     *     NF, INF, JNF, TOTDOF, ITEST)                                       43
      READ (NIN,8010) NMODES                                                  44
C
C*                           *******************************
C*                           *                             *
C*                           * SYSTEM STIFFNESS MATRIX AND *
C*                           *     MASS MATRIX ASSEMBLY    *
C*                           *                             *
C*                           *******************************
C
      CALL MATNUL(SYSK, ISYSK, JSYSK, TOTDOF, TOTDOF, ITEST)                  45
      CALL VECNUL(SYSM, ISYSM, TOTDOF, ITEST)                                 46
      DOFEL = NODEL*DOFNOD                                                    47
      CALL QQUA4(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)                 48
      DO 1090 NELE=1,TOTELS                                                   49
      AREA = 0.0D0                                                            50
      CALL ELGEOM(NELE, ELTOP, IELTOP, JELTOP, COORD, ICOORD,                 51
     *     JCOORD, GEOM, IGEOM, JGEOM, DIMEN, ITEST)                          51
C
C                            INTEGRATION LOOP FOR ELEMENT STIFFNESS
C                            AND ELEMENT LUMPED MASS USING NQP
C                            QUADRATURE POINTS
C
      CALL MATNUL(ELK, IELK, JELK, DOFEL, DOFEL, ITEST)                       52
      CALL MATNUL(ELM, IELM, JELM, DOFEL, DOFEL, ITEST)                       53
      DO 1080 IQUAD=1,NQP                                                     54
C
C                            FORM LINEAR SHAPE FUNCTION AND SPACE
C                            DERIVATIVES IN THE LOCAL CORRDINATES.
C                            TRANSFORM LOCAL DERIVATIVES TO GLOBAL
C                            COORDINATE SYSTEM
C
      XI = ABSS(1,IQUAD)                                                      55
      ETA = ABSS(2,IQUAD)                                                     56
      CALL QUAM4(FUN, IFUN, DER, IDER, JDER, XI, ETA, ITEST)                  57
      CALL MATMUL(DER, IDER, JDER, GEOM, IGEOM, JGEOM, JAC, IJAC,             58
     *     JJAC, DIMEN, NODEL, DIMEN, ITEST)                                  58
      CALL MATINV(JAC, IJAC, JJAC, JACIN, IJACIN, JJACIN, DIMEN,              59
     *     DET, ITEST)                                                        59
      CALL MATMUL(JACIN, IJACIN, JJACIN, DER, IDER, JDER, DERIV,              60
     *     IDERIV, JDERIV, DIMEN, DIMEN, NODEL, ITEST)                        60
C
C                            FORMATION OF STRAIN-DISPLACEMENT
C                            MATRIX B AND FORMATION OF INTEGRAND
C                            FOR ELEMENT STIFFNESS MATRIX ELK
C
      CALL B2C2(B, IB, JB, DERIV, IDERIV, JDERIV, NODEL, ITEST)               61
      CALL MATMUL(D, ID, JD, B, IB, JB, DB, IDB, JDB, NUMSS, NUMSS,           62
     *     DOFEL, ITEST)                                                      62
      CALL MATRAN(B, IB, JB, BT, IBT, JBT, NUMSS, DOFEL, ITEST)               63
      CALL MATMUL(BT, IBT, JBT, DB, IDB, JDB, BTDB, IBTDB, JBTDB,             64
     *     DOFEL, NUMSS, DOFEL, ITEST)                                        64
C
C                            FORMATION OF INTEGRAND FOR ELEMENT
C                            MASS MATRIX ELM
C
      CALL SHAPFN(SHP, ISHP, JSHP, FUN, IFUN, NODEL, DOFNOD, ITEST)           65
      CALL MATRAN(SHP, ISHP, JSHP, TSHP, ITSHP, JTSHP, DOFNOD,                66
     *     DOFEL, ITEST)                                                      66
      CALL MATMUL(TSHP, ITSHP, JTSHP, SHP, ISHP, JSHP, FTF, IFTF,             67
     *     JFTF, DOFEL, DOFNOD, DOFEL, ITEST)                                 67
C
      QUOT = DET*WGHT(IQUAD)                                                  68
      AREA = AREA + QUOT                                                      69
      DO 1070 I=1,DOFEL                                                       70
      DO 1060 J=1,DOFEL                                                       71
      BTDB(I,J) = BTDB(I,J)*QUOT                                              72
      FTF(I,J) = FTF(I,J)*QUOT                                                73
 1060 CONTINUE                                                                74
 1070 CONTINUE                                                                75
      CALL MATADD(ELK, IELK, JELK, BTDB, IBTDB, JBTDB, DOFEL,                 76
     *     DOFEL, ITEST)                                                      76
      CALL MATADD(ELM, IELM, JELM, FTF, IFTF, JFTF, DOFEL, DOFEL,             77
     *     ITEST)                                                             77
 1080 CONTINUE                                                                78
C
C                            ASSEMBLY OF SYSTEM STIFFNESS MATRIX
C                            SYSK, SYSTEM LUMPED MASS MATRIX SYSM
C
      CALL DIRECT(NELE, ELTOP, IELTOP, JELTOP, NF, INF, JNF,                  79
     *     DOFNOD, STEER, ISTEER, ITEST)                                      79
      CALL ASFUL(SYSK, ISYSK, JSYSK, ELK, IELK, JELK, STEER,                  80
     *     ISTEER, DOFEL, ITEST)                                              80
      CALL ASLMS(SYSM, ISYSM, ELM, IELM, JELM, STEER, ISTEER,                 81
     *     DOFEL, DOFNOD, AREA, ITEST)                                        81
 1090 CONTINUE                                                                82
C
C*                           *********************
C*                           *                   *
C*                           * EQUATION SOLUTION *
C*                           *                   *
C*                           *********************
C
      CALL PRTVEC(SYSM, ISYSM, TOTDOF, ITEST)                                 83
C
C                            REDUCTION OF SYSTEM MASS MATRIX BY
C                            CHOLESKI REDUCTION AND CONSTRUCTION OF
C                            MODIFIED SYSTEM STIFFNESS MATRIX
C
      DO 1100 I=1,TOTDOF                                                      84
      SYSM(I) = 1.0D0/DSQRT(SYSM(I)*RHO)                                      85
 1100 CONTINUE                                                                86
      DO 1120 I=1,TOTDOF                                                      87
      DO 1110 J=1,TOTDOF                                                      88
      SYSK(I,J) = SYSM(I)*SYSK(I,J)*SYSM(J)                                   89
 1110 CONTINUE                                                                90
 1120 CONTINUE                                                                91
C
C                            HOUSEHOLDER TRIDIAGONALISATION OF SYSK
C                            AND EIGENVALUE AND EIGENVECTOR
C                            RECOVERY USING QL TRANSFORMATIONS
C
      TOL = VTOL(X)                                                           92
      CALL HOUSE(SYSK, ISYSK, JSYSK, SYSK, ISYSK, JSYSK, DIAG,                93
     *     IDIAG, SUB, ISUB, TOTDOF, TOL, ITEST)                              93
      EPS = VEPS(X)                                                           94
      CALL QLVEC(DIAG, IDIAG, SUB, ISUB, SYSK, ISYSK, JSYSK,                  95
     *     TOTDOF, EPS, ITEST)                                                95
      CALL PRTVEC(DIAG, IDIAG, TOTDOF, ITEST)                                 96
      DO 1140 J=1,NMODES                                                      97
      DO 1130 I=1,TOTDOF                                                      98
      DIAG(I) = SYSK(I,J)*SYSM(I)                                             99
 1130 CONTINUE                                                               100
      CALL PRTVEC(DIAG, IDIAG, TOTDOF, ITEST)                                101
 1140 CONTINUE                                                               102
      STOP                                                                   103
 8010 FORMAT (16I5)                                                          104
 8020 FORMAT (I5, 6F10.0)                                                    105
 8030 FORMAT (6F10.0)                                                        106
 9010 FORMAT (//24H**** NODAL GEOMETRY ****//)                               107
 9020 FORMAT (16I5)                                                          108
 9030 FORMAT (I5, 6F10.5)                                                    109
 9040 FORMAT (//26H**** ELEMENT TOPOLOGY ****//)                             110
 9050 FORMAT (//29H**** MATERIAL PROPERTIES ****//)                          111
 9060 FORMAT (2F10.5)                                                        112
 9070 FORMAT (//24H**** RESTRAINT DATA ****//)                               113
      END                                                                    114
