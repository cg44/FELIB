C***********************************************************************
C$SPLIT$SEG3P1$*********************************************************
C***********************************************************************
      INTEGER BNDNOD, BNODE, DIF, DIMEN, DOFEL, DOFNOD, ELNUM,                 1
     *     ELTOP, ELTYP, HBAND, I, IABSS, ICOORD, ID, IDL,                     1
     *     IDT, IDTPD, IELK, IELTOP, IFUN, IJAC, IJACIN, ILOADS,               1
     *     IGEOM, INF, IP, IPD, IQUAD, ISTEER, ISYSK, ITEST,                   1
     *     IWGHT, J, JABSS, JCOORD, JD, JDL, JDT, JDTPD,                       1
     *     JELK, JELTOP, JJAC, JJACIN, JGEOM, JNF, JP, JPD,                    1
     *     JSYSK, NELE, NF, NIN, NODEL, NODNUM, NOUT, NQP,                     1
     *     NWORK, STEER, TOTDOF, TOTELS, TOTNOD                                1
      DOUBLE PRECISION ABSS, BVAL, COORD, D, DET, DL, DT, DTPD,                2
     *     ELK, ETA, FUN, JAC, JACIN, LOADS, GEOM, P, PD,                      2
     *     PERM, QUOT, SCALE, SYSK, WGHT, WORK, XI                             2
      LOGICAL FIRST                                                            3
      DIMENSION ABSS(3,9), BNODE(30), BVAL(30), COORD(100,3),                  4
     *     D(3,8), DL(3,8), DT(8,3), DTPD(8,8), ELK(8,8),                      4
     *     ELTOP(100,10), FUN(8), JAC(3,3), JACIN(3,3),                        4
     *     LOADS(100), GEOM(8,3), NF(100,1), NWORK(8),                         4
     *     P(3,3), PD(3,8), PERM(3), STEER(8), SYSK(100,25),                   4
     *     WGHT(9), WORK(3)                                                    4
      DATA IABSS /3/, ICOORD /100/, ID /3/, IDL /3/, IDT /8/,                  5
     *     IDTPD /8/, IELK /8/, IELTOP /100/, IFUN /8/,                        5
     *     IJAC /3/, IJACIN /3/, ILOADS /100/, IGEOM /8/,                      5
     *     INF /100/, IP /3/, IPD /3/, ISTEER /8/, ISYSK /100/,                5
     *     IWGHT /9/, JABSS /9/, JCOORD /3/, JD /8/, JDL /8/,                  5
     *     JDT /3/, JDTPD /8/, JELK /8/, JELTOP /10/, JJAC /3/,                5
     *     JJACIN /3/, JGEOM /3/, JNF /1/, JP /3/, JPD /8/,                    5
     *     JSYSK /25/, SCALE /1.0D+10/                                         5
C
      DATA NIN /5/, NOUT /6/                                                   6
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
C                            INPUT OF PERMEABILITIES AND CON-
C                            STRUCTION OF PERMEABILITY MATRIX P
C
      WRITE (NOUT,9050)                                                       30
      READ (NIN,8030) (PERM(I),I=1,DIMEN)                                     31
      WRITE (NOUT,9060) (PERM(I),I=1,DIMEN)                                   32
      CALL MATNUL(P, IP, JP, DIMEN, DIMEN, ITEST)                             33
      DO 1050 I=1,DIMEN                                                       34
      P(I,I) = PERM(I)                                                        35
 1050 CONTINUE                                                                36
C
C                            INPUT OF NUMBER OF DEGREES OF FREEDOM
C                            PER NODE, INPUT OF BOUNDARY CON-
C                            DITIONS AND CONSTRUCTION OF NODAL
C                            FREEDOM ARRAY NF
C
      WRITE (NOUT,9070)                                                       37
      READ (NIN,8010) DOFNOD                                                  38
      WRITE (NOUT,9020) DOFNOD                                                39
      READ (NIN,8010) BNDNOD                                                  40
      WRITE (NOUT,9020) BNDNOD                                                41
      DO 1060 I=1,BNDNOD                                                      42
      READ (NIN,8040) BNODE(I), BVAL(I)                                       43
      WRITE (NOUT,9030) BNODE(I), BVAL(I)                                     44
 1060 CONTINUE                                                                45
      TOTDOF = 0                                                              46
      DO 1080 I=1,TOTNOD                                                      47
      DO 1070 J=1,DOFNOD                                                      48
      TOTDOF = TOTDOF + 1                                                     49
      NF(I,J) = TOTDOF                                                        50
 1070 CONTINUE                                                                51
 1080 CONTINUE                                                                52
C
C                            CALCULATION OF SEMI-BANDWIDTH
C
      FIRST = .TRUE.                                                          53
      DO 1090 NELE=1,TOTELS                                                   54
      CALL FREDIF(NELE, ELTOP, IELTOP, JELTOP, NF, INF, JNF,                  55
     *     DOFNOD, FIRST, DIF, ITEST)                                         55
 1090 CONTINUE                                                                56
      HBAND = DIF + 1                                                         57
C
C*                           ************************************
C*                           *                                  *
C*                           * SYSTEM STIFFNESS MATRIX ASSEMBLY *
C*                           *                                  *
C*                           ************************************
C
      CALL MATNUL(SYSK, ISYSK, JSYSK, TOTDOF, HBAND, ITEST)                   58
      DOFEL = NODEL*DOFNOD                                                    59
      CALL QQUA4(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)                 60
      DO 1130 NELE=1,TOTELS                                                   61
      CALL ELGEOM(NELE, ELTOP, IELTOP, JELTOP, COORD, ICOORD,                 62
     *     JCOORD, GEOM, IGEOM, JGEOM, DIMEN, ITEST)                          62
C
C                            INTEGRATION LOOP FOR ELEMENT STIFFNESS
C                            USING NQP QUADRATURE POINTS
C
      CALL MATNUL(ELK, IELK, JELK, DOFEL, DOFEL, ITEST)                       63
      DO 1120 IQUAD=1,NQP                                                     64
C
C                            FORM LINEAR SHAPE FUNCTION AND SPACE
C                            DERIVATIVES IN THE LOCAL CORRDINATES.
C                            TRANSFORM LOCAL DERIVATIVES TO GLOBAL
C                            COORDINATE SYSTEM
C
      XI = ABSS(1,IQUAD)                                                      65
      ETA = ABSS(2,IQUAD)                                                     66
      CALL QUAM4(FUN, IFUN, DL, IDL, JDL, XI, ETA, ITEST)                     67
      CALL MATMUL(DL, IDL, JDL, GEOM, IGEOM, JGEOM, JAC, IJAC,                68
     *     JJAC, DIMEN, NODEL, DIMEN, ITEST)                                  68
      CALL MATINV(JAC, IJAC, JJAC, JACIN, IJACIN, JJACIN, DIMEN,              69
     *     DET, ITEST)                                                        69
      CALL MATMUL(JACIN, IJACIN, JJACIN, DL, IDL, JDL, D, ID, JD,             70
     *     DIMEN, DIMEN, NODEL, ITEST)                                        70
C
C                            FORMATION OF ELEMENT STIFFNESS ELK
C
      CALL MATMUL(P, IP, JP, D, ID, JD, PD, IPD, JPD, DIMEN, DIMEN,           71
     *     DOFEL, ITEST)                                                      71
      CALL MATRAN(D, ID, JD, DT, IDT, JDT, DIMEN, DOFEL, ITEST)               72
      CALL MATMUL(DT, IDT, JDT, PD, IPD, JPD, DTPD, IDTPD, JDTPD,             73
     *     DOFEL, DIMEN, DOFEL, ITEST)                                        73
      QUOT = DET*WGHT(IQUAD)                                                  74
      DO 1110 I=1,DOFEL                                                       75
      DO 1100 J=1,DOFEL                                                       76
      DTPD(I,J) = DTPD(I,J)*QUOT                                              77
 1100 CONTINUE                                                                78
 1110 CONTINUE                                                                79
      CALL MATADD(ELK, IELK, JELK, DTPD, IDTPD, JDTPD, DOFEL,                 80
     *     DOFEL, ITEST)                                                      80
 1120 CONTINUE                                                                81
C
C                            ASSEMBLY OF SYSTEM STIFFNESS MATRIX
C                            SYSK
C
      CALL DIRECT(NELE, ELTOP, IELTOP, JELTOP, NF, INF, JNF,                  82
     *     DOFNOD, STEER, ISTEER, ITEST)                                      82
      CALL ASSYM(SYSK, ISYSK, JSYSK, ELK, IELK, JELK, STEER,                  83
     *     ISTEER, HBAND, DOFEL, ITEST)                                       83
 1130 CONTINUE                                                                84
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
      CALL VECNUL(LOADS, ILOADS, TOTDOF, ITEST)                               85
      DO 1140 I=1,BNDNOD                                                      86
      J = BNODE(I)                                                            87
      SYSK(J,HBAND) = SYSK(J,HBAND)*SCALE                                     88
      LOADS(J) = SYSK(J,HBAND)*BVAL(I)                                        89
 1140 CONTINUE                                                                90
C
C                            SOLUTION OF SYSTEM MATRIX FOR THE
C                            NODAL VALUES OF THE POTENTIAL
C
      CALL CHOSOL(SYSK, ISYSK, JSYSK, LOADS, ILOADS, TOTDOF, HBAND,           91
     *     ITEST)                                                             91
      WRITE (NOUT,9080)                                                       92
      CALL PRTVEC(LOADS, ILOADS, TOTDOF, ITEST)                               93
      STOP                                                                    94
 8010 FORMAT (16I5)                                                           95
 8020 FORMAT (I5, 6F10.0)                                                     96
 8030 FORMAT (2F10.0)                                                         97
 8040 FORMAT (I5, 6F10.0)                                                     98
 9010 FORMAT (//24H**** NODAL GEOMETRY ****//)                                99
 9020 FORMAT (16I5)                                                          100
 9030 FORMAT(I5,6F10.5)                                                      101
 9040 FORMAT (//26H**** ELEMENT TOPOLOGY ****//)                             102
 9050 FORMAT (//24H**** PERMEABILITIES ****//)                               103
 9060 FORMAT (2F10.5)                                                        104
 9070 FORMAT (//29H**** BOUNDARY CONDITIONS ****//)                          105
 9080 FORMAT (//26H**** NODAL POTENTIALS ****//)                             106
      END                                                                    107
