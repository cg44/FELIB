C***********************************************************************
C$SPLIT$SEG4P1$*********************************************************
C***********************************************************************
      INTEGER BNDNOD, BNODE, DIF, DIMEN, DOFEL, DOFNOD, ELNUM,                 1
     *     ELTOP, ELTYP, HBAND, I, IABSS, ICOORD, ID, IDL,                     1
     *     IDT, IDTPD, IELK, IELM, IELTOP, IFTF, IFUN, IJAC,                   1
     *     IJACIN, ILOADS, IGEOM, INF, INTNOD, IP, IPD, IQUAD,                 1
     *     ISTEER, ISYSK, ISYSM, ITEST, IWGHT, IWORK1, IWORK2,                 1
     *     J, JABSS, JCOORD, JD, JDL, JDT, JDTPD, JELK,                        1
     *     JELM, JELTOP, JFTF, JJAC, JJACIN, JGEOM, JNF, JP,                   1
     *     JPD, JSYSK, JSYSM, NELE, NF, NIN, NODEL,                            1
     *     NODNUM, NOUT, NQP, NSTEPS, NWORK, STEER,                            1
     *     TOTDOF, TOTELS, TOTNOD                                              1
      DOUBLE PRECISION ABSS, BVAL, COORD, D, DET, DL, DT, DTIM,                2
     *     DTPD, ELK, ELM, ETA, FTF, FUN, JAC, JACIN, LOADS,                   2
     *     GEOM, P, PD, PERM, QUOT, SCALE, SYSK, SYSM, THETA,                  2
     *     WGHT, WORK1, WORK2, XI                                              2
      LOGICAL FIRST                                                            3
      DIMENSION ABSS(3,9), BNODE(30), BVAL(30), COORD(100,3),                  4
     *     D(3,8), DL(3,8), DT(8,3), DTPD(8,8), ELK(8,8),                      4
     *     ELM(8,8), ELTOP(100,10), FTF(8,8), FUN(8), JAC(3,3),                4
     *     JACIN(3,3), LOADS(100), GEOM(8,3), NF(100,1),                       4
     *     NWORK(8), P(3,3), PD(3,8), PERM(3), STEER(8),                       4
     *     SYSK(100,25), SYSM(100,25), WGHT(9), WORK1(100),                    4
     *     WORK2(100)                                                          4
      DATA IABSS /3/, ICOORD /100/, ID /3/, IDL /3/, IDT /8/,                  6
     *     IDTPD /8/, IELK /8/, IELM /8/, IELTOP /100/,                        6
     *     IFTF /8/, IFUN /8/, IJAC /3/, IJACIN /3/, ILOADS /100/,             6
     *     IGEOM /8/, INF /100/, IP /3/, IPD /3/, ISTEER /8/,                  6
     *     ISYSK /100/, ISYSM /100/, IWGHT /9/, IWORK1 /100/,                  6
     *     IWORK2 /100/, JABSS /9/, JCOORD /3/, JD /8/,                        6
     *     JDL /8/, JDT /3/, JDTPD /8/, JELK /8/, JELM /8/,                    6
     *     JELTOP /10/, JFTF /8/, JJAC /3/, JJACIN /3/, JGEOM /3/,             6
     *     JNF /1/, JP /3/, JPD /8/, JSYSK /25/, JSYSM /25/,                   6
     *     SCALE /1.0D+20/                                                     6
C
      DATA NIN /5/, NOUT /6/                                                   7
C
C                            SET ITEST FOR FULL CHECKING
C
      ITEST = 0                                                                8
C
C*                           **********************
C*                           *                    *
C*                           * INPUT DATA SECTION *
C*                           *                    *
C*                           **********************
C
C                            INPUT OF NODAL GEOMETRY
C
      WRITE (NOUT,9010)                                                        9
      READ (NIN,8010) TOTNOD, DIMEN                                           10
      WRITE (NOUT,9020) TOTNOD, DIMEN                                         11
      DO 1020 I=1,TOTNOD                                                      12
      READ (NIN,8020) NODNUM, (WORK1(J),J=1,DIMEN)                            13
      WRITE (NOUT,9030) NODNUM, (WORK1(J),J=1,DIMEN)                          14
      DO 1010 J=1,DIMEN                                                       15
      COORD(NODNUM,J) = WORK1(J)                                              16
 1010 CONTINUE                                                                17
 1020 CONTINUE                                                                18
C
C                            INPUT OF ELEMENT TOPOLOGY
C
      WRITE (NOUT,9040)                                                       19
      READ (NIN,8010) ELTYP, TOTELS, NODEL                                    20
      WRITE (NOUT,9020) ELTYP, TOTELS, NODEL                                  21
      DO 1040 I=1,TOTELS                                                      22
      READ (NIN,8010) ELNUM, (NWORK(J),J=1,NODEL)                             23
      WRITE (NOUT,9020) ELNUM, (NWORK(J),J=1,NODEL)                           24
      ELTOP(ELNUM,1) = ELTYP                                                  25
      ELTOP(ELNUM,2) = NODEL                                                  26
      DO 1030 J=1,NODEL                                                       27
      ELTOP(ELNUM,J+2) = NWORK(J)                                             28
 1030 CONTINUE                                                                29
 1040 CONTINUE                                                                30
C
C                            INPUT OF PERMEABILITIES AND CON-
C                            STRUCTION OF PERMEABILITY MATRIX P
C
      WRITE (NOUT,9050)                                                       31
      READ (NIN,8030) (PERM(I),I=1,DIMEN)                                     32
      WRITE (NOUT,9060) (PERM(I),I=1,DIMEN)                                   33
      CALL MATNUL(P, IP, JP, DIMEN, DIMEN, ITEST)                             34
      DO 1050 I=1,DIMEN                                                       35
      P(I,I) = PERM(I)                                                        36
 1050 CONTINUE                                                                37
C
C                            INPUT OF NUMBER OF DEGREES OF FREEDOM
C                            PER NODE, INPUT OF BOUNDARY CON-
C                            DITIONS AND CONSTRUCTION OF NODAL
C                            FREEDOM ARRAY NF
C
      WRITE (NOUT,9070)                                                       38
      READ (NIN,8010) DOFNOD                                                  39
      WRITE (NOUT,9020) DOFNOD                                                40
      READ (NIN,8010) BNDNOD                                                  41
      WRITE (NOUT,9020) BNDNOD                                                42
      DO 1060 I=1,BNDNOD                                                      43
      READ (NIN,8040) BNODE(I), BVAL(I)                                       44
      WRITE (NOUT,9030) BNODE(I), BVAL(I)                                     45
 1060 CONTINUE                                                                46
      TOTDOF = 0                                                              47
      DO 1080 I=1,TOTNOD                                                      48
      DO 1070 J=1,DOFNOD                                                      49
      TOTDOF = TOTDOF + 1                                                     50
      NF(I,J) = TOTDOF                                                        51
 1070 CONTINUE                                                                52
 1080 CONTINUE                                                                53
C
C                            INPUT OF INITIAL CONDINTIONS
C
      WRITE (NOUT,9080)                                                       54
      CALL VECNUL(LOADS, ILOADS, TOTDOF, ITEST)                               55
      READ (NIN,8010) INTNOD                                                  56
      WRITE (NOUT,9090) INTNOD                                                57
      DO 1090 I=1,INTNOD                                                      58
      READ (NIN,8040) NODNUM, WORK1(1)                                        59
      WRITE (NOUT,9030) NODNUM, WORK1(1)                                      60
      J = NF(NODNUM,1)                                                        61
      LOADS(J) = WORK1(1)                                                     62
 1090 CONTINUE                                                                63
C
C                            INPUT TIME STEPPING DATA
C
      WRITE (NOUT,9100)                                                       64
      READ (NIN,8020) NSTEPS, DTIM, THETA                                     65
      WRITE (NOUT,9110) NSTEPS, DTIM, THETA                                   66
C
C                            CALCULATION OF SEMI-BANDWIDTH
C
      FIRST = .TRUE.                                                          67
      DO 1100 NELE=1,TOTELS                                                   68
      CALL FREDIF(NELE, ELTOP, IELTOP, JELTOP, NF, INF, JNF,                  69
     *     DOFNOD, FIRST, DIF, ITEST)                                         69
 1100 CONTINUE                                                                70
      HBAND = DIF + 1                                                         71
C
C*                           ************************************
C*                           *                                  *
C*                           * SYSTEM STIFFNESS MATRIX ASSEMBLY *
C*                           *                                  *
C*                           ************************************
C
      CALL MATNUL(SYSK, ISYSK, JSYSK, TOTDOF, HBAND, ITEST)                   72
      CALL MATNUL(SYSM, ISYSM, JSYSM, TOTDOF, HBAND, ITEST)                   73
      DOFEL = NODEL*DOFNOD                                                    74
      CALL QQUA4(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)                 75
      DO 1140 NELE=1,TOTELS                                                   76
      CALL ELGEOM(NELE, ELTOP, IELTOP, JELTOP, COORD, ICOORD,                 77
     *     JCOORD, GEOM, IGEOM, JGEOM, DIMEN, ITEST)                          77
C
C                            INTEGRATION LOOP FOR ELEMENT STIFFNESS
C                            USING NQP QUADRATURE POINTS
C
      CALL MATNUL(ELK, IELK, JELK, DOFEL, DOFEL, ITEST)                       78
      CALL MATNUL(ELM, IELM, JELM, DOFEL, DOFEL, ITEST)                       79
      DO 1130 IQUAD=1,NQP                                                     80
C
C                            FORM LINEAR SHAPE FUNCTION AND SPACE
C                            DERIVATIVES IN THE LOCAL CORRDINATES.
C                            TRANSFORM LOCAL DERIVATIVES TO GLOBAL
C                            COORDINATE SYSTEM
C
      XI = ABSS(1,IQUAD)                                                      81
      ETA = ABSS(2,IQUAD)                                                     82
      CALL QUAM4(FUN, IFUN, DL, IDL, JDL, XI, ETA, ITEST)                     83
      CALL MATMUL(DL, IDL, JDL, GEOM, IGEOM, JGEOM, JAC, IJAC,                84
     *     JJAC, DIMEN, NODEL, DIMEN, ITEST)                                  84
      CALL MATINV(JAC, IJAC, JJAC, JACIN, IJACIN, JJACIN, DIMEN,              85
     *     DET, ITEST)                                                        85
      CALL MATMUL(JACIN, IJACIN, JJACIN, DL, IDL, JDL, D, ID, JD,             86
     *     DIMEN, DIMEN, NODEL, ITEST)                                        86
C
C                            FORMATION OF ELEMENT STIFFNESS ELK
C
      CALL MATMUL(P, IP, JP, D, ID, JD, PD, IPD, JPD, DIMEN, DIMEN,           87
     *     DOFEL, ITEST)                                                      87
      CALL MATRAN(D, ID, JD, DT, IDT, JDT, DIMEN, DOFEL, ITEST)               88
      CALL MATMUL(DT, IDT, JDT, PD, IPD, JPD, DTPD, IDTPD, JDTPD,             89
     *     DOFEL, DIMEN, DOFEL, ITEST)                                        89
      QUOT = DET*WGHT(IQUAD)                                                  90
      DO 1120 I=1,DOFEL                                                       91
      DO 1110 J=1,DOFEL                                                       92
      FTF(I,J) = FUN(I)*FUN(J)*QUOT                                           93
      DTPD(I,J) = DTPD(I,J)*QUOT                                              94
 1110 CONTINUE                                                                95
 1120 CONTINUE                                                                96
      CALL MATADD(ELM, IELM, JELM, FTF, IFTF, JFTF, DOFEL, DOFEL,             97
     *     ITEST)                                                             97
      CALL MATADD(ELK, IELK, JELK, DTPD, IDTPD, JDTPD, DOFEL,                 98
     *     DOFEL, ITEST)                                                      98
 1130 CONTINUE                                                                99
C
C                            ASSEMBLY OF SYSTEM STIFFNESS MATRIX
C                            SYSK
C
      CALL DIRECT(NELE, ELTOP, IELTOP, JELTOP, NF, INF, JNF,                 100
     *     DOFNOD, STEER, ISTEER, ITEST)                                     100
      CALL ASSYM(SYSM, ISYSM, JSYSM, ELM, IELM, JELM, STEER,                 101
     *     ISTEER, HBAND, DOFEL, ITEST)                                      101
      CALL ASSYM(SYSK, ISYSK, JSYSK, ELK, IELK, JELK, STEER,                 102
     *     ISTEER, HBAND, DOFEL, ITEST)                                      102
 1140 CONTINUE                                                               103
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
      DO 1160 I=1,TOTDOF                                                     104
      DO 1150 J=1,HBAND                                                      105
      WORK1(1) = THETA*SYSK(I,J) + SYSM(I,J)/DTIM                            106
      WORK1(2) = (1.0D0-THETA)*SYSK(I,J) - SYSM(I,J)/DTIM                    107
      SYSK(I,J) = WORK1(1)                                                   108
      SYSM(I,J) = -WORK1(2)                                                  109
 1150 CONTINUE                                                               110
 1160 CONTINUE                                                               111
      CALL VECNUL(WORK1, IWORK1, TOTDOF, ITEST)                              112
      DO 1170 I=1,BNDNOD                                                     113
      J = BNODE(I)                                                           114
      SYSK(J,HBAND) = SYSK(J,HBAND)*SCALE                                    115
      WORK1(J) = SYSK(J,HBAND)*BVAL(I)                                       116
 1170 CONTINUE                                                               117
C
C                            SOLUTION OF SYSTEM MATRIX FOR THE
C                            NODAL VALUES OF THE POTENTIAL
C
      CALL CHORDN(SYSK, ISYSK, JSYSK, TOTDOF, HBAND, ITEST)                  118
      DO 1180 I=1,NSTEPS                                                     119
      WRITE (NOUT,9120) I                                                    120
      CALL MVSLT(SYSM, ISYSM, JSYSM, LOADS, ILOADS, WORK2, IWORK2,           121
     *     TOTDOF, HBAND, ITEST)                                             121
      CALL VECADD(WORK2, IWORK2, WORK1, IWORK1, TOTDOF, ITEST)               122
      CALL CHOFWD(SYSK, ISYSK, JSYSK, WORK2, IWORK2, TOTDOF, HBAND,          123
     *     ITEST)                                                            123
      CALL CHOBAK(SYSK, ISYSK, JSYSK, WORK2, IWORK2, TOTDOF, HBAND,          124
     *     ITEST)                                                            124
      CALL VECCOP(WORK2, IWORK2, LOADS, ILOADS, TOTDOF, ITEST)               125
      CALL PRTVEC(LOADS, ILOADS, TOTDOF, ITEST)                              126
 1180 CONTINUE                                                               127
      STOP                                                                   128
 8010 FORMAT (16I5)                                                          129
 8020 FORMAT (I5, 6F10.0)                                                    130
 8030 FORMAT (2F10.0)                                                        131
 8040 FORMAT (I5, 6F10.0)                                                    132
 9010 FORMAT (//24H**** NODAL GEOMETRY ****//)                               133
 9020 FORMAT (16I5)                                                          134
 9030 FORMAT (I5, 6F10.5)                                                    135
 9040 FORMAT (//26H**** ELEMENT TOPOLOGY ****//)                             136
 9050 FORMAT (//24H**** PERMEABILITIES ****//)                               137
 9060 FORMAT (2F10.5)                                                        138
 9070 FORMAT (//29H**** BOUNDARY CONDITIONS ****//)                          139
 9080 FORMAT (//28H**** INITIAL CONDITIONS ****//)                           140
 9090 FORMAT (16I5)                                                          141
 9100 FORMAT (//28H**** TIME STEPPING DATA ****//)                           142
 9110 FORMAT (I5, 6F10.5)                                                    143
 9120 FORMAT (//31H**** NODAL POTENTIALS FOR STEP , I3, 6H  ****//)          144
      END                                                                    145
