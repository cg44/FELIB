C***********************************************************************
C$SPLIT$SEG1P1$*********************************************************
C***********************************************************************
      INTEGER DIF, DIMEN, DOFEL, DOFNOD, ELNUM, ELTOP, ELTYP,                  1
     *     HBAND, I, IABSS, IB, IBT, IBTDB, ICOORD, ID, IDB,                   1
     *     IDER, IDERIV, IELDIS, IELK, IELTOP, IFUN, IGEOM,                    1
     *     IJAC, IJACIN, ILOADS, INF, IQUAD, IRESTR, ISTEER,                   1
     *     ISTRN, ISTRS, ISYSK, ITEST, IWGHT, J, JABSS, JB,                    1
     *     JBT, JBTDB, JCOORD, JD, JDB, JDER, JDERIV, JELK,                    1
     *     JELTOP, JGEOM, JJAC, JJACIN, JNF, JRESTR, JSTEER,                   1
     *     JSYSK, K, LODNOD, NELE, NF, NIN, NODEL, NODNUM,                     1
     *     NOUT, NQP, NTEMP, NUMSS, NWORK, RESNOD, RESTR, STEER,               1
     *     TOTDOF, TOTELS, TOTNOD                                              1
      DOUBLE PRECISION ABSS, B, BT, BTDB, COORD, D, DB, DER,                   2
     *     DERIV, DET, E, ELDIS, ELK, ETA, FUN, GEOM, JAC,                     2
     *     JACIN, LOADS, NU, QUOT, STRN, STRS, SYSK, WGHT,                     2
     *     WORK, XI                                                            2
      LOGICAL FIRST                                                            3
      DIMENSION ABSS(3,9), B(6,24), BT(24,6), BTDB(24,24),                     4
     *     COORD(100,3), D(6,6), DB(6,24), DER(3,8), DERIV(3,8),               4
     *     ELDIS(24), ELK(24,24), ELTOP(100,10), FUN(8),                       4
     *     GEOM(8,3), JAC(3,3), JACIN(3,3), LOADS(100),                        4
     *     NF(100,3), NWORK(8), RESTR(100,4), STEER(24),                       4
     *     STRN(6), STRS(6), SYSK(100,25), WGHT(9), WORK(3)                    4
      DATA IABSS /3/, IB /6/, IBT /24/, IBTDB /24/, ICOORD /100/,              5
     *     ID /6/, IDB /6/, IDER /3/, IDERIV /3/, IELDIS /24/,                 5
     *     IELK /24/, IELTOP /100/, IFUN /8/, IGEOM /8/,                       5
     *     IJAC /3/, IJACIN /3/, ILOADS /100/, INF /100/,                      5
     *     IRESTR /100/, ISTEER /24/, ISTRN /6/, ISTRS /6/,                    5
     *     ISYSK /100/, IWGHT /9/, JABSS /9/, JB /24/, JBT /6/,                5
     *     JBTDB /24/, JCOORD /3/, JD /6/, JDB /24/, JDER /8/,                 5
     *     JDERIV /8/, JELK /24/, JELTOP /10/, JGEOM /3/,                      5
     *     JJAC /3/, JJACIN /3/, JNF /3/, JRESTR /4/,                          5
     *     JSYSK /25/                                                          5
C
      DATA NIN /5/, NOUT /6/, NTEMP /7/                                        6
C
C                            SET ITEST FOR FULL CHECKING
C
      ITEST = 0                                                                7
C
C*                           **********************
C*                           *                    *
C*                           * INPUT DATA SECTION *
C*                           *                    *
C*                           **********************
C
C                            INPUT OF NODAL GEOMETRY
C
      WRITE (NOUT,9010)                                                        8
      READ (NIN,8010) TOTNOD, DIMEN                                            9
      WRITE (NOUT,9020) TOTNOD, DIMEN                                         10
      DO 1020 I=1,TOTNOD                                                      11
      READ (NIN,8020) NODNUM, (WORK(J),J=1,DIMEN)                             12
      WRITE (NOUT,9030) NODNUM, (WORK(J),J=1,DIMEN)                           13
      DO 1010 J=1,DIMEN                                                       14
      COORD(NODNUM,J) = WORK(J)                                               15
 1010 CONTINUE                                                                16
 1020 CONTINUE                                                                17
C
C                            INPUT OF ELEMENT TOPOLOGY
C
      WRITE (NOUT,9040)                                                       18
      READ (NIN,8010) ELTYP, TOTELS, NODEL                                    19
      WRITE (NOUT,9020) ELTYP, TOTELS, NODEL                                  20
      DO 1040 I=1,TOTELS                                                      21
      READ (NIN,8010) ELNUM, (NWORK(J),J=1,NODEL)                             22
      WRITE (NOUT,9020) ELNUM, (NWORK(J),J=1,NODEL)                           23
      ELTOP(ELNUM,1) = ELTYP                                                  24
      ELTOP(ELNUM,2) = NODEL                                                  25
      DO 1030 J=1,NODEL                                                       26
      ELTOP(ELNUM,J+2) = NWORK(J)                                             27
 1030 CONTINUE                                                                28
 1040 CONTINUE                                                                29
C
C                            INPUT OF MATERIAL PROPERTIES AND
C                            CONSTRUCTION OF STRESS-STRAIN MATRIX D
C                            FOR PLANE STRAIN
C
      WRITE (NOUT,9050)                                                       30
      READ (NIN,8030) NU, E                                                   31
      WRITE (NOUT,9060) NU, E                                                 32
      CALL DPSN(D, ID, JD, E, NU, NUMSS, ITEST)                               33
C
C                            INPUT OF NUMBER OF DEGREES OF FREEDOM
C                            PER NODE, INPUT OF RESTRAINED NODE
C                            DATA AND CONSTRUCTION OF NODAL FREEDOM
C                            ARRAY NF
C
      WRITE (NOUT,9070)                                                       34
      READ (NIN,8010) DOFNOD                                                  35
      WRITE (NOUT,9020) DOFNOD                                                36
      READ (NIN,8010) RESNOD                                                  37
      WRITE (NOUT,9020) RESNOD                                                38
      K = DOFNOD + 1                                                          39
      DO 1050 I=1,RESNOD                                                      40
      READ (NIN,8010) (RESTR(I,J),J=1,K)                                      41
      WRITE (NOUT,9020) (RESTR(I,J),J=1,K)                                    42
 1050 CONTINUE                                                                43
      CALL FORMNF(RESTR, IRESTR, JRESTR, RESNOD, TOTNOD, DOFNOD,              44
     *     NF, INF, JNF, TOTDOF, ITEST)                                       44
C
C                            LOADING DATA INPUT
C
      WRITE (NOUT,9080)                                                       45
      CALL VECNUL(LOADS, ILOADS, TOTDOF, ITEST)                               46
      READ (NIN,8010) LODNOD                                                  47
      WRITE (NOUT,9020) LODNOD                                                48
      DO 1070 I=1,LODNOD                                                      49
      READ (NIN,8020) NODNUM, (WORK(J),J=1,DOFNOD)                            50
      WRITE (NOUT,9030) NODNUM, (WORK(J),J=1,DOFNOD)                          51
      DO 1060 J=1,DOFNOD                                                      52
      K = NF(NODNUM,J)                                                        53
      IF (K.EQ.0) GO TO 1060                                                  54
      LOADS(K) = WORK(J)                                                      55
 1060 CONTINUE                                                                56
 1070 CONTINUE                                                                57
C
C                            CALCULATION OF SEMI-BANDWIDTH
C
      FIRST = .TRUE.                                                          58
      DO 1080 NELE=1,TOTELS                                                   59
      CALL FREDIF(NELE, ELTOP, IELTOP, JELTOP, NF, INF, JNF,                  60
     *     DOFNOD, FIRST, DIF, ITEST)                                         60
 1080 CONTINUE                                                                61
      HBAND = DIF + 1                                                         62
C
C*                           ************************************
C*                           *                                  *
C*                           * GLOBAL STIFFNESS MATRIX ASSEMBLY *
C*                           *                                  *
C*                           ************************************
C
      CALL MATNUL(SYSK, ISYSK, JSYSK, TOTDOF, HBAND, ITEST)                   63
      DOFEL = NODEL*DOFNOD                                                    64
      CALL QQUA4(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)                 65
      DO 1120 NELE=1,TOTELS                                                   66
      CALL ELGEOM(NELE, ELTOP, IELTOP, JELTOP, COORD, ICOORD,                 67
     *     JCOORD, GEOM, IGEOM, JGEOM, DIMEN, ITEST)                          67
C
C                            INTEGRATION LOOP FOR ELEMENT STIFFNESS
C                            USING NQP QUADRATURE POINTS
C
      CALL MATNUL(ELK, IELK, JELK, DOFEL, DOFEL, ITEST)                       68
      DO 1110 IQUAD=1,NQP                                                     69
C
C                            FORM LINEAR SHAPE FUNCTION AND SPACE
C                            DERIVATIVES IN THE LOCAL CORRDINATES.
C                            TRANSFORM LOCAL DERIVATIVES TO GLOBAL
C                            COORDINATE SYSTEM
C
      XI = ABSS(1,IQUAD)                                                      70
      ETA = ABSS(2,IQUAD)                                                     71
      CALL QUAM4(FUN, IFUN, DER, IDER, JDER, XI, ETA, ITEST)                  72
      CALL MATMUL(DER, IDER, JDER, GEOM, IGEOM, JGEOM, JAC, IJAC,             73
     *     JJAC, DIMEN, NODEL, DIMEN, ITEST)                                  73
      CALL MATINV(JAC, IJAC, JJAC, JACIN, IJACIN, JJACIN, DIMEN,              74
     *     DET, ITEST)                                                        74
      CALL MATMUL(JACIN, IJACIN, JJACIN, DER, IDER, JDER, DERIV,              75
     *     IDERIV, JDERIV, DIMEN, DIMEN, NODEL, ITEST)                        75
C
C                            FORMATION OF STRAIN-DISPLACEMENT
C                            MATRIX B AND OUTPUT TO WORK FILE FOR
C                            LATER RECOVERY PROCESS
C
      CALL B2C2(B, IB, JB, DERIV, IDERIV, JDERIV, NODEL, ITEST)               76
      WRITE (NTEMP) ((B(I,J),I=1,NUMSS),J=1,DOFEL)                            77
C
C                            FORMATION OF ELEMENT STIFFNESS ELK
C
      CALL MATMUL(D, ID, JD, B, IB, JB, DB, IDB, JDB, NUMSS, NUMSS,           78
     *     DOFEL, ITEST)                                                      78
      CALL MATRAN(B, IB, JB, BT, IBT, JBT, NUMSS, DOFEL, ITEST)               79
      CALL MATMUL(BT, IBT, JBT, DB, IDB, JDB, BTDB, IBTDB, JBTDB,             80
     *     DOFEL, NUMSS, DOFEL, ITEST)                                        80
      QUOT = DET*WGHT(IQUAD)                                                  81
      DO 1100 I=1,DOFEL                                                       82
      DO 1090 J=1,DOFEL                                                       83
      BTDB(I,J) = BTDB(I,J)*QUOT                                              84
 1090 CONTINUE                                                                85
 1100 CONTINUE                                                                86
      CALL MATADD(ELK, IELK, JELK, BTDB, IBTDB, JBTDB, DOFEL,                 87
     *     DOFEL, ITEST)                                                      87
 1110 CONTINUE                                                                88
C
C                            ASSEMBLY OF SYSTEM STIFFNESS MATRIX
C                            SYSK
C
      CALL DIRECT(NELE, ELTOP, IELTOP, JELTOP, NF, INF, JNF,                  89
     *     DOFNOD, STEER, ISTEER, ITEST)                                      89
      CALL ASSYM(SYSK, ISYSK, JSYSK, ELK, IELK, JELK, STEER,                  90
     *     ISTEER, HBAND, DOFEL, ITEST)                                       90
 1120 CONTINUE                                                                91
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
      WRITE (NOUT,9090)                                                       92
      CALL PRTVEC(LOADS, ILOADS, TOTDOF, ITEST)                               93
      CALL CHOSOL(SYSK, ISYSK, JSYSK, LOADS, ILOADS, TOTDOF, HBAND,           94
     *     ITEST)                                                             94
      WRITE (NOUT,9100)                                                       95
      CALL PRTVEC(LOADS, ILOADS, TOTDOF, ITEST)                               96
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
      REWIND NTEMP                                                            97
      DO 1150 NELE=1,TOTELS                                                   98
      WRITE (NOUT,9110) NELE                                                  99
C
C                            SELECT NODAL DISPLACEMENTS FOR ELEMENT
C                            NELE USING THE STEERING VECTOR
C
      CALL DIRECT(NELE, ELTOP, IELTOP, JELTOP, NF, INF, JNF,                 100
     *     DOFNOD, STEER, ISTEER, ITEST)                                     100
      DO 1140 IQUAD=1,NQP                                                    101
      CALL VECNUL(ELDIS, IELDIS, DOFEL, ITEST)                               102
      DO 1130 I=1,DOFEL                                                      103
      JSTEER = STEER(I)                                                      104
      IF (JSTEER.NE.0) ELDIS(I) = LOADS(JSTEER)                              105
 1130 CONTINUE                                                               106
C
C                            RECOVER B MATRIX FOR POINT IQUAD AND
C                            CALCULATE THE STRESSES AND STRAINS
C
      READ (NTEMP) ((B(I,J),I=1,NUMSS),J=1,DOFEL)                            107
      CALL MATVEC(B, IB, JB, ELDIS, IELDIS, NUMSS, DOFEL, STRN,              108
     *     ISTRN, ITEST)                                                     108
      CALL MATVEC(D, ID, JD, STRN, ISTRN, NUMSS, NUMSS, STRS,                109
     *     ISTRS, ITEST)                                                     109
      WRITE (NOUT,9120) IQUAD                                                110
      CALL PRTVEC(STRN, ISTRN, NUMSS, ITEST)                                 111
      CALL PRTVEC(STRS, ISTRS, NUMSS, ITEST)                                 112
 1140 CONTINUE                                                               113
 1150 CONTINUE                                                               114
      STOP                                                                   115
 8010 FORMAT (16I5)                                                          116
 8020 FORMAT (I5, 6F10.0)                                                    117
 8030 FORMAT (2F10.0)                                                        118
 9010 FORMAT (1H //24H**** NODAL GEOMETRY ****//)                            119
 9020 FORMAT (16I5)                                                          120
 9030 FORMAT (I5, 6F10.5)                                                    121
 9040 FORMAT (1H //26H**** ELEMENT TOPOLOGY ****//)                          122
 9050 FORMAT (1H //29H**** MATERIAL PROPERTIES ****//)                       123
 9060 FORMAT (2F10.5)                                                        124
 9070 FORMAT (1H //24H**** RESTRAINT DATA ****//)                            125
 9080 FORMAT (1H //22H**** LOADING DATA ****//)                              126
 9090 FORMAT (1H //26H **** VECTOR OF LOADS ****//)                          127
 9100 FORMAT (1H //42H **** EQUILIBRIUM DISPLACEMENTS (U,V) ****//)          128
 9110 FORMAT (1H //39H **** STRAINS AND STRESSES FOR ELEMENT , I3,           129
     *     5H ****//)                                                        129
 9120 FORMAT (1H /18H QUADRATURE POINT , I3/)                                130
      END                                                                    131
