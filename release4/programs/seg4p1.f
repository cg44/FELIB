C***********************************************************************
C
C     Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
C                          Chilton, Didcot, Oxfordshire OX11 0QX
C
C     Contact: Prof Chris Greenough
C     Email: c.greenough@rl.ac.uk
C     Tel: (44) 1235 445307
C
C     Finite Library Program: Seg4p1 - Potential Flow (Diffusion Eqn)
C
C***********************************************************************
C
      PROGRAM SEG4P1
C
      INTEGER BNDNOD,BNODE,DIF,DIMEN,DOFEL,DOFNOD,ELNUM,ELTOP,ELTYP,
     *        HBAND,I,IABSS,ICOORD,IDTPD,IELK,IELM,IELTOP,IFTF,IFUN,
     *        IGDER,IGDERT,IGEOM,IJAC,IJACIN,ILDER,INF,INTNOD,IP,IPD,
     *        IQUAD,IRHS,ISTEER,ISYSK,ISYSM,ITEST,IWGHT,IWORK1,IWORK2,J,
     *        JABSS,JCOORD,JDTPD,JELK,JELM,JELTOP,JFTF,JGDER,JGDERT,
     *        JGEOM,JJAC,JJACIN,JLDER,JNF,JP,JPD,JSYSK,JSYSM,K,NELE,NF,
     *        NIN,NODEL,NODNUM,NOUT,NQP,NSTEPS,STEER,TOTDOF,TOTELS,
     *        TOTNOD
      DOUBLE PRECISION ABSS,BVAL,COORD,DET,DTIM,DTPD,ELK,ELM,ETA,FTF,
     *                 FUN,GDER,GDERT,GEOM,JAC,JACIN,LDER,P,PD,PERM,
     *                 QUOT,RHS,SCALE,SYSK,SYSM,THETA,WGHT,WORK1,WORK2,
     *                 XI
      LOGICAL FIRST
      DIMENSION ABSS(3,9),BNODE(30),BVAL(30),DTPD(8,8),ELK(8,8),
     *          ELM(8,8),FTF(8,8),FUN(8),GDER(3,8),GDERT(8,3),GEOM(8,3),
     *          JAC(3,3),JACIN(3,3),LDER(3,8),P(3,3),PD(3,8),PERM(3),
     *          STEER(8),WGHT(9)
C
C     Problem size dependent arrays
C
      DIMENSION COORD(100,3),ELTOP(100,10),NF(100,1),RHS(100),
     *          SYSK(100,25),SYSM(100,25),WORK1(100),WORK2(100)
C
      DATA IABSS/3/,IDTPD/8/,IELK/8/,IELM/8/,IFTF/8/,IFUN/8/,IGDER/3/,
     *     IGDERT/8/,IGEOM/8/,IJAC/3/,IJACIN/3/,ILDER/3/,IP/3/,IPD/3/,
     *     ISTEER/8/,IWGHT/9/,JABSS/9/,JCOORD/3/,JDTPD/8/,JELK/8/,
     *     JELM/8/,JFTF/8/,JGDER/8/,JGDERT/3/,JGEOM/3/,JJAC/3/,
     *     JJACIN/3/,JLDER/8/,JNF/1/,JP/3/,JPD/8/,SCALE/1.0D+10/
C
C     Problem size dependent data statements
C
      DATA ICOORD/100/,IELTOP/100/,INF/100/,IRHS/100/,ISYSK/100/,
     *     ISYSM/100/,IWORK1/100/,IWORK2/100/,JELTOP/10/,JSYSK/25/,
     *     JSYSM/25/
C
      DATA NIN/5/,NOUT/6/
C
C     Set ITEST for full checking
C
      ITEST = 0
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
      DO 1000 I = 1,TOTNOD
          READ (NIN,FMT=9980) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
          WRITE (NOUT,FMT=9930) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
 1000 CONTINUE
C
C     Input of element topology
C
      WRITE (NOUT,FMT=9920)
      READ (NIN,FMT=9990) ELTYP,TOTELS,NODEL
      WRITE (NOUT,FMT=9940) ELTYP,TOTELS,NODEL
      DO 1010 I = 1,TOTELS
          READ (NIN,FMT=9990) ELNUM, (ELTOP(ELNUM,J+2),J=1,NODEL)
          WRITE (NOUT,FMT=9940) ELNUM, (ELTOP(ELNUM,J+2),J=1,NODEL)
          ELTOP(ELNUM,1) = ELTYP
          ELTOP(ELNUM,2) = NODEL
 1010 CONTINUE
C
C     Input of permeabilities and construction of permeability
C     matrix P
C
      WRITE (NOUT,FMT=9910)
      READ (NIN,FMT=9970) (PERM(I),I=1,DIMEN)
      WRITE (NOUT,FMT=9900) (PERM(I),I=1,DIMEN)
      CALL MATNUL(P,IP,JP,DIMEN,DIMEN,ITEST)
      DO 1020 I = 1,DIMEN
          P(I,I) = PERM(I)
 1020 CONTINUE
C
C     Input of number of degrees of freedom per node, input of
C     boundary conditions and construction of nodal freedom array NF
C
      WRITE (NOUT,FMT=9890)
      READ (NIN,FMT=9990) DOFNOD
      WRITE (NOUT,FMT=9940) DOFNOD
      READ (NIN,FMT=9990) BNDNOD
      WRITE (NOUT,FMT=9940) BNDNOD
      DO 1030 I = 1,BNDNOD
          READ (NIN,FMT=9960) BNODE(I),BVAL(I)
          WRITE (NOUT,FMT=9930) BNODE(I),BVAL(I)
 1030 CONTINUE
      TOTDOF = 0
      DO 1050 I = 1,TOTNOD
          DO 1040 J = 1,DOFNOD
              TOTDOF = TOTDOF + 1
              NF(I,J) = TOTDOF
 1040     CONTINUE
 1050 CONTINUE
C
C     Input of initial condintions
C
      WRITE (NOUT,FMT=9880)
      CALL VECNUL(RHS,IRHS,TOTDOF,ITEST)
      READ (NIN,FMT=9990) INTNOD
      WRITE (NOUT,FMT=9870) INTNOD
      DO 1060 I = 1,INTNOD
          READ (NIN,FMT=9960) NODNUM,WORK1(1)
          WRITE (NOUT,FMT=9930) NODNUM,WORK1(1)
          J = NF(NODNUM,1)
          RHS(J) = WORK1(1)
 1060 CONTINUE
C
C     Input time stepping data
C
      WRITE (NOUT,FMT=9860)
      READ (NIN,FMT=9980) NSTEPS,DTIM,THETA
      WRITE (NOUT,FMT=9850) NSTEPS,DTIM,THETA
C
C     Calculation of semi-bandwidth
C
      FIRST = .TRUE.
      DO 1070 NELE = 1,TOTELS
          CALL FREDIF(NELE,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,FIRST,
     *                DIF,ITEST)
 1070 CONTINUE
      HBAND = DIF + 1
C
C     *                           ************************************
C     *                           *                                  *
C     *                           * System Stiffness Matrix Assembly *
C     *                           *                                  *
C     *                           ************************************
C
      CALL MATNUL(SYSK,ISYSK,JSYSK,TOTDOF,HBAND,ITEST)
      CALL MATNUL(SYSM,ISYSM,JSYSM,TOTDOF,HBAND,ITEST)
      DOFEL = NODEL*DOFNOD
      CALL QQUA4(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST)
      DO 1110 NELE = 1,TOTELS
          CALL ELGEOM(NELE,ELTOP,IELTOP,JELTOP,COORD,ICOORD,JCOORD,GEOM,
     *                IGEOM,JGEOM,DIMEN,ITEST)
C
C     Integration loop for element stiffness using NQP quadrature
C     points
C
          CALL MATNUL(ELK,IELK,JELK,DOFEL,DOFEL,ITEST)
          CALL MATNUL(ELM,IELM,JELM,DOFEL,DOFEL,ITEST)
          DO 1100 IQUAD = 1,NQP
C
C     Form linear shape function and space derivatives in the local
C     corrdinates. Transform local derivatives to global coordinate
C     system
C
              XI = ABSS(1,IQUAD)
              ETA = ABSS(2,IQUAD)
              CALL QUAM4(FUN,IFUN,LDER,ILDER,JLDER,XI,ETA,ITEST)
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
              QUOT = ABS(DET)*WGHT(IQUAD)
              DO 1090 I = 1,DOFEL
                  DO 1080 J = 1,DOFEL
                      FTF(I,J) = FUN(I)*FUN(J)*QUOT
                      DTPD(I,J) = DTPD(I,J)*QUOT
 1080             CONTINUE
 1090         CONTINUE
              CALL MATADD(ELM,IELM,JELM,FTF,IFTF,JFTF,DOFEL,DOFEL,ITEST)
              CALL MATADD(ELK,IELK,JELK,DTPD,IDTPD,JDTPD,DOFEL,DOFEL,
     *                    ITEST)
 1100     CONTINUE
C
C     Assembly of system stiffness matrix SYSK
C
          CALL DIRECT(NELE,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,STEER,
     *                ISTEER,ITEST)
          CALL ASSYM(SYSM,ISYSM,JSYSM,ELM,IELM,JELM,STEER,ISTEER,HBAND,
     *               DOFEL,ITEST)
          CALL ASSYM(SYSK,ISYSK,JSYSK,ELK,IELK,JELK,STEER,ISTEER,HBAND,
     *               DOFEL,ITEST)
 1110 CONTINUE
C
C     *                           *********************
C     *                           *                   *
C     *                           * Equation Solution *
C     *                           *                   *
C     *                           *********************
C
C     Modification of stiffness matrix and right-hand side to
C     implement boundary conditions
C
      DO 1130 I = 1,TOTDOF
          DO 1120 J = 1,HBAND
              WORK1(1) = THETA*SYSK(I,J) + SYSM(I,J)/DTIM
              WORK1(2) = (1.0D0-THETA)*SYSK(I,J) - SYSM(I,J)/DTIM
              SYSK(I,J) = WORK1(1)
              SYSM(I,J) = -WORK1(2)
 1120     CONTINUE
 1130 CONTINUE
      DO 1140 I = 1,BNDNOD
          J = BNODE(I)
          SYSK(J,HBAND) = SYSK(J,HBAND)*SCALE
          WORK1(I) = SYSK(J,HBAND)*BVAL(I)
 1140 CONTINUE
C
C     Solution of system matrix for the nodal values of the
C     potential
C
      CALL CHORDN(SYSK,ISYSK,JSYSK,TOTDOF,HBAND,ITEST)
      DO 1160 I = 1,NSTEPS
          WRITE (NOUT,FMT=9840) I
          CALL MVSYB(SYSM,ISYSM,JSYSM,RHS,IRHS,WORK2,IWORK2,TOTDOF,
     *               HBAND,ITEST)
C
C     Insertion of boundary conditions
C
          DO 1150 J = 1,BNDNOD
              K = BNODE(J)
              WORK2(K) = WORK1(J)
 1150     CONTINUE
C
C     Forward and backward substitution
C
          CALL CHOSUB(SYSK,ISYSK,JSYSK,WORK2,IWORK2,TOTDOF,HBAND,ITEST)
          CALL VECCOP(WORK2,IWORK2,RHS,IRHS,TOTDOF,ITEST)
          CALL PRTVAL(RHS,IRHS,NF,INF,JNF,DOFNOD,TOTNOD,NOUT,ITEST)
 1160 CONTINUE
C
C     End of Program
C
      STOP
C
 9990 FORMAT (16I5)
 9980 FORMAT (I5,6F10.0)
 9970 FORMAT (2F10.0)
 9960 FORMAT (I5,6F10.0)
 9950 FORMAT (/,/,' **** NODAL GEOMETRY ****',/,/,' ')
 9940 FORMAT (' ',16I5)
 9930 FORMAT (' ',I5,6F10.5)
 9920 FORMAT (/,/,' **** ELEMENT TOPOLOGY ****',/,/,' ')
 9910 FORMAT (/,/,' **** PERMEABILITIES ****',/,/,' ')
 9900 FORMAT (' ',2F10.5)
 9890 FORMAT (/,/,' **** BOUNDARY CONDITIONS ****',/,/,' ')
 9880 FORMAT (/,/,' **** INITIAL CONDITIONS ****',/,/,' ')
 9870 FORMAT (' ',16I5)
 9860 FORMAT (/,/,' **** TIME STEPPING DATA ****',/,/,' ')
 9850 FORMAT (' ',I5,6F10.5)
 9840 FORMAT (/,/,' **** NODAL POTENTIALS FOR STEP ',I3,'  ****',/,/,
     *       ' ')
C
      END
