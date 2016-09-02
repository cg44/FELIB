C***********************************************************************
C     $Id: seg3p1.f,v 1.2 2009/02/27 11:44:50 cg44 Exp $
C
C     Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
C                          Chilton, Didcot, Oxfordshire OX11 0QX
C
C     Contact: Prof Chris Greenough
C     Email: c.greenough@rl.ac.uk
C     Tel: (44) 1235 445307
C
C     Finite Library Program: Seg3p1 - Steady State Potential Flow
C
C***********************************************************************
C
      PROGRAM SEG3P1
C
      INTEGER BNDNOD,BNODE,DIF,DIMEN,DOFEL,DOFNOD,ELNUM,ELTOP,ELTYP,
     *        HBAND,I,IABSS,ICOORD,IDTPD,IELK,IELTOP,IELQ,IFUN,IGDER,
     *        IGDERT,IGEOM,IJAC,IJACIN,ILDER,INF,IP,IPD,IQUAD,IRHS,
     *        ISCVEC,ISTEER,ISYSK,ITEST,IWGHT,J,JABSS,JCOORD,JDTPD,JELK,
     *        JELTOP,JGDER,JGDERT,JGEOM,JJAC,JJACIN,JLDER,JNF,JP,JPD,
     *        JSYSK,NELE,NF,NIN,NODEL,NODNUM,NOUT,NQP,STEER,TOTDOF,
     *        TOTELS,TOTNOD
      DOUBLE PRECISION ABSS,BVAL,COORD,DET,DTPD,ELK,ELQ,ETA,FUN,GDER,
     *                 GDERT,GEOM,JAC,JACIN,LDER,P,PD,QUOT,RHS,SCALE,
     *                 SCVEC,SOURCE,STRGTH,SYSK,WGHT,X,XI,Y
      LOGICAL FIRST
      DIMENSION ABSS(3,9),DTPD(8,8),ELK(8,8),ELQ(8),FUN(8),GDER(3,8),
     *          GDERT(8,3),GEOM(8,3),JAC(3,3),JACIN(3,3),LDER(3,8),
     *          P(3,3),PD(3,8),SCVEC(8),STEER(8),WGHT(9)
C
C     Problem size dependent arrays
C
      DIMENSION BNODE(30),BVAL(30),COORD(100,3),ELTOP(100,10),NF(100,1),
     *          RHS(100),SYSK(100,25)
C
      DATA IABSS/3/,IDTPD/8/,IELK/8/,IELQ/8/,IFUN/8/,IGDER/3/,IGDERT/8/,
     *     IGEOM/8/,IJAC/3/,IJACIN/3/,ILDER/3/,IP/3/,IPD/3/,ISCVEC/8/,
     *     ISTEER/8/,IWGHT/9/,JABSS/9/,JCOORD/3/,JDTPD/8/,JELK/8/,
     *     JGDER/8/,JGDERT/3/,JGEOM/3/,JJAC/3/,JJACIN/3/,JLDER/8/,
     *     JNF/1/,JP/3/,JPD/8/,SCALE/1.0D+10/
C
C     Problem size dependent data statements
C
      DATA ICOORD/100/,IELTOP/100/,INF/100/,IRHS/100/,ISYSK/100/,
     *     JELTOP/10/,JSYSK/25/
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
      WRITE (NOUT,9950)
      READ (NIN,9990) TOTNOD,DIMEN
      WRITE (NOUT,9940) TOTNOD,DIMEN
      DO 1000 I = 1,TOTNOD
          READ (NIN,9980) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
          WRITE (NOUT,9930) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
 1000 CONTINUE
C
C     Input of element topology
C
      WRITE (NOUT,9920)
      READ (NIN,9990) ELTYP,TOTELS,NODEL
      WRITE (NOUT,9940) ELTYP,TOTELS,NODEL
      DO 1010 I = 1,TOTELS
          READ (NIN,9990) ELNUM, (ELTOP(ELNUM,J+2),J=1,NODEL)
          WRITE (NOUT,9940) ELNUM, (ELTOP(ELNUM,J+2),J=1,NODEL)
          ELTOP(ELNUM,1) = ELTYP
          ELTOP(ELNUM,2) = NODEL
 1010 CONTINUE
C
C     Input of permeabilities, construction of permeability matrix P
C     and SOURCE strength
C
      WRITE (NOUT,9910)
      CALL MATNUL(P,IP,JP,DIMEN,DIMEN,ITEST)
      READ (NIN,9970) (P(I,I),I=1,DIMEN)
      WRITE (NOUT,9900) (P(I,I),I=1,DIMEN)
C
      WRITE (NOUT,9890)
      READ (NIN,9970) STRGTH
      WRITE (NOUT,9900) STRGTH
C
C     Input of number of degrees of freedom per node, input of
C     boundary conditions and construction of nodal freedom array NF
C
      WRITE (NOUT,9880)
      READ (NIN,9990) DOFNOD
      WRITE (NOUT,9940) DOFNOD
      READ (NIN,9990) BNDNOD
      WRITE (NOUT,9940) BNDNOD
      DO 1020 I = 1,BNDNOD
          READ (NIN,9960) BNODE(I),BVAL(I)
          WRITE (NOUT,9930) BNODE(I),BVAL(I)
 1020 CONTINUE
      TOTDOF = 0
      DO 1040 I = 1,TOTNOD
          DO 1030 J = 1,DOFNOD
              TOTDOF = TOTDOF + 1
              NF(I,J) = TOTDOF
 1030     CONTINUE
 1040 CONTINUE
C
C     Calculation of semi-bandwidth
C
      FIRST = .TRUE.
      DO 1050 NELE = 1,TOTELS
          CALL FREDIF(NELE,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,FIRST,
     *                DIF,ITEST)
 1050 CONTINUE
      HBAND = DIF + 1
C
C     *                           ************************************
C     *                           *                                  *
C     *                           * System Stiffness Matrix Assembly *
C     *                           *                                  *
C     *                           ************************************
C
      CALL VECNUL(RHS,IRHS,TOTDOF,ITEST)
      CALL MATNUL(SYSK,ISYSK,JSYSK,TOTDOF,HBAND,ITEST)
      DOFEL = NODEL*DOFNOD
      CALL QQUA4(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST)
      DO 1090 NELE = 1,TOTELS
          CALL ELGEOM(NELE,ELTOP,IELTOP,JELTOP,COORD,ICOORD,JCOORD,GEOM,
     *                IGEOM,JGEOM,DIMEN,ITEST)
C
C     Integration loop for element stiffness using NQP quadrature
C     points
C
          CALL MATNUL(ELK,IELK,JELK,DOFEL,DOFEL,ITEST)
          CALL VECNUL(ELQ,IELQ,DOFEL,ITEST)
          DO 1080 IQUAD = 1,NQP
C
C     Form linear shape function and space derivatives in the local
C     corrdinates. Transform local derivatives to global coordinate
C     system
C
              XI = ABSS(1,IQUAD)
              ETA = ABSS(2,IQUAD)
              CALL QUAM4(FUN,IFUN,LDER,ILDER,JLDER,XI,ETA,ITEST)
C
              CALL SCAPRD(GEOM(1,1),IGEOM,FUN,IFUN,NODEL,X,ITEST)
              CALL SCAPRD(GEOM(1,2),IGEOM,FUN,IFUN,NODEL,Y,ITEST)
C
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
              DO 1070 I = 1,DOFEL
                  DO 1060 J = 1,DOFEL
                      DTPD(I,J) = DTPD(I,J)*QUOT
 1060             CONTINUE
                  SCVEC(I) = FUN(I)*SOURCE(X,Y,STRGTH)*QUOT
 1070         CONTINUE
              CALL MATADD(ELK,IELK,JELK,DTPD,IDTPD,JDTPD,DOFEL,DOFEL,
     *                    ITEST)
              CALL VECADD(ELQ,IELQ,SCVEC,ISCVEC,DOFEL,ITEST)
 1080     CONTINUE
C
C     Assembly of system stiffness matrix
C
          CALL DIRECT(NELE,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,STEER,
     *                ISTEER,ITEST)
          CALL ASSYM(SYSK,ISYSK,JSYSK,ELK,IELK,JELK,STEER,ISTEER,HBAND,
     *               DOFEL,ITEST)
          CALL ASRHS(RHS,IRHS,ELQ,IELQ,STEER,ISTEER,DOFEL,ITEST)
 1090 CONTINUE
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
      DO 1100 I = 1,BNDNOD
          J = BNODE(I)
          SYSK(J,HBAND) = SYSK(J,HBAND)*SCALE
          RHS(J) = SYSK(J,HBAND)*BVAL(I)
 1100 CONTINUE
C
C     Solution of system matrix for the nodal values of the
C     potential
C
      CALL CHOSOL(SYSK,ISYSK,JSYSK,RHS,IRHS,TOTDOF,HBAND,ITEST)
      WRITE (NOUT,9870)
      CALL PRTVAL(RHS,IRHS,NF,INF,JNF,DOFNOD,TOTNOD,NOUT,ITEST)
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
 9890 FORMAT (/,/,' **** SOURCE STRENGTH ****',/,/,' ')
 9880 FORMAT (/,/,' **** BOUNDARY CONDITIONS ****',/,/,' ')
 9870 FORMAT (/,/,' **** NODAL POTENTIALS ****',/,/,' ')
C
      END
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION SOURCE(X,Y,STRGTH)
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION STRGTH,X,Y
C
      SOURCE = 0.0D0
      IF ((X.GE.1.0D0) .AND. (X.LE.2.0D0) .AND. (Y.GE.1.0D0) .AND.
     *    (Y.LE.2.0D0)) SOURCE = STRGTH
C
      END
