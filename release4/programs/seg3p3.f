C***********************************************************************
C
C     Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
C                          Chilton, Didcot, Oxfordshire OX11 0QX
C
C     Contact: Prof Chris Greenough
C     Email: c.greenough@rl.ac.uk
C     Tel: (44) 1235 445307
C
C     Finite Library Program: Seg1p1 - Two Dimensional Eddy Current
C                                      Solution
C
C***********************************************************************
C
      PROGRAM SEG3P3
C
      INTEGER BNDNOD,BNODE,DIMEN,DOFEL,DOFNOD,ELNUM,ELTOP,ELTYP,HBAND,I,
     *        IABSS,IBNODE,ICOORD,IDTPD,IELJ,IELK,IELM,IELPOT,IELTOP,
     *        IFIELD,IFTF,IFUN,IGDER,IGDERT,IGEOM,IJAC,IJACIN,IJVAL,
     *        ILDER,INF,IPD,IQUAD,IRHS,IS,ISCVEC,ISTEER,ISYSK,ITEST,
     *        IWGHT,J,JABSS,JCOORD,JDTPD,JELK,JELM,JELTOP,JFTF,JGDER,
     *        JGDERT,JGEOM,JJAC,JJACIN,JLDER,JNF,JPD,JS,JSTEER,JSYSK,
     *        NELE,NF,NIN,NJV,NODEL,NODNUM,NOUT,NQP,NTEMP,STEER,TOTDOF,
     *        TOTELS,TOTNOD
      DOUBLE PRECISION ABSS,BMOD,BVAL,COND,CONST,COORD,DET,DTPD,ELJ,ELK,
     *                 ELM,ELPOT,ETA,FIELD,FREQ,FTF,FUN,GDER,GDERT,GEOM,
     *                 JAC,JACIN,JVAL,LDER,MU0,OMEGA,PD,PERM,PI,QUOT,
     *                 RHS,S,SCALE,SCVEC,SYSK,WGHT,XI
      DIMENSION ABSS(3,9),DTPD(8,8),ELJ(8),ELK(8,8),ELM(8,8),ELPOT(8),
     *          FIELD(8),FTF(8,8),FUN(8),GDER(3,8),GDERT(8,3),GEOM(8,3),
     *          JAC(3,3),JACIN(3,3),LDER(3,8),PD(3,8),S(3,3),SCVEC(8),
     *          STEER(8),WGHT(9)
C
C     Problem size dependent arrays
C
      DIMENSION BNODE(30),BVAL(2,30),COORD(100,3),ELTOP(100,10),
     *          JVAL(100),NF(100,1),RHS(2,100),SYSK(2,100,25)
C
      DATA IABSS/3/,IDTPD/8/,IELJ/8/,IELK/8/,IELM/8/,IELPOT/8/,
     *     IFIELD/8/,IFTF/8/,IFUN/8/,IGDER/3/,IGDERT/8/,IGEOM/8/,
     *     IJAC/3/,IJACIN/3/,ILDER/3/,IPD/3/,IS/3/,ISCVEC/8/,ISTEER/8/,
     *     IWGHT/9/,JABSS/9/,JCOORD/3/,JDTPD/8/,JELK/8/,JELM/8/,JFTF/8/,
     *     JGDER/8/,JGDERT/3/,JGEOM/3/,JJAC/3/,JJACIN/3/,JLDER/8/,
     *     JNF/1/,JPD/8/,JS/3/,SCALE/1.0D+10/
C
C     Problem size dependent data statements
C
      DATA IBNODE/30/,ICOORD/100/,IELTOP/100/,IJVAL/100/,INF/100/,
     *     IRHS/100/,ISYSK/100/,JELTOP/10/,JSYSK/25/
C
      DATA NIN/5/,NOUT/6/,NTEMP/7/
C
C     Set ITEST for full checking
C
      ITEST = 0
C
C     Set value of PI
C
      PI = ATAN(1.0D0)*4.0D0
      MU0 = 4.0D0*PI*1.0D-07
C
C     *                           **********************
C     *                           *                    *
C     *                           * Input Data Section *
C     *                           *                    *
C     *                           **********************
C
C     Input of nodal geometry
C
      WRITE (NOUT,FMT=9960)
      READ (NIN,FMT=9990) TOTNOD,DIMEN
      WRITE (NOUT,FMT=9950) TOTNOD,DIMEN
C
      DO 1000 I = 1,TOTNOD
          READ (NIN,FMT=9980) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
          WRITE (NOUT,FMT=9940) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
 1000 CONTINUE
C
C     Input of element topology ELTYP=2 for conductor, ELTYP=1
C     otherwise
C
      WRITE (NOUT,FMT=9920)
      READ (NIN,FMT=9990) TOTELS
      WRITE (NOUT,FMT=9950) TOTELS
C
      DO 1010 I = 1,TOTELS
          READ (NIN,FMT=9990) ELNUM,ELTYP,NODEL,
     *      (ELTOP(ELNUM,J+2),J=1,NODEL)
          WRITE (NOUT,FMT=9950) ELNUM,ELTYP,NODEL,
     *      (ELTOP(ELNUM,J+2),J=1,NODEL)
          ELTOP(ELNUM,1) = ELTYP
          ELTOP(ELNUM,2) = NODEL
 1010 CONTINUE
C
C     Input of permeabilities and construction of susceptibility
C     matrix S
C
C
      WRITE (NOUT,FMT=9910)
      READ (NIN,FMT=9970) FREQ,COND,PERM
      WRITE (NOUT,FMT=9900) FREQ,COND,PERM
      OMEGA = 2.0D0*PI*FREQ
C
C     Input of number of degrees of freedom per node, input of
C     boundary conditions
C
      WRITE (NOUT,FMT=9890)
      READ (NIN,FMT=9990) DOFNOD
      WRITE (NOUT,FMT=9950) DOFNOD
C
      READ (NIN,FMT=9990) BNDNOD
      WRITE (NOUT,FMT=9950) BNDNOD
C
C
      DO 1020 I = 1,BNDNOD
          READ (NIN,FMT=9980) BNODE(I),BVAL(1,I),BVAL(2,I)
          WRITE (NOUT,FMT=9940) BNODE(I),BVAL(1,I),BVAL(2,I)
 1020 CONTINUE
C
C     Input number of current carrying elements input element number
C     and current value
C
      WRITE (NOUT,FMT=9880)
      CALL VECNUL(JVAL,IJVAL,TOTELS,ITEST)
      READ (NIN,FMT=9990) NJV
      WRITE (NOUT,FMT=9950) NJV
      IF (NJV.NE.0) THEN
C
          DO 1030 I = 1,NJV
              READ (NIN,FMT=9980) ELNUM,JVAL(ELNUM)
              WRITE (NOUT,FMT=9930) ELNUM,JVAL(ELNUM)
 1030     CONTINUE
      END IF
C
C     Construction of nodal freedom array NF
C
      TOTDOF = 0
      DO 1050 I = 1,TOTNOD
          DO 1040 J = 1,DOFNOD
              TOTDOF = TOTDOF + 1
              NF(I,J) = TOTDOF
 1040     CONTINUE
 1050 CONTINUE
C
C     Calculation of semi-bandwidth
C
      CALL BNDWTH(ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,TOTELS,HBAND,
     *            ITEST)
C
C     *                           ************************************
C     *                           *                                  *
C     *                           * System Stiffness Matrix Assembly *
C     *                           *                                  *
C     *                           ************************************
C
      CALL CVCNUL(RHS,IRHS,TOTDOF,ITEST)
      CALL CMTNUL(SYSK,ISYSK,JSYSK,TOTDOF,HBAND,ITEST)
      DOFEL = NODEL*DOFNOD
      CALL QQUA4(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST)
C
C     Main matrix assembly
C
      DO 1100 NELE = 1,TOTELS
          CALL ELGEOM(NELE,ELTOP,IELTOP,JELTOP,COORD,ICOORD,JCOORD,GEOM,
     *                IGEOM,JGEOM,DIMEN,ITEST)
C
          ELTYP = ELTOP(NELE,1)
          CONST = 0.0D0
          IF (ELTYP.EQ.2) CONST = OMEGA*COND
          CALL MATNUL(S,IS,JS,DIMEN,DIMEN,ITEST)
          DO 1060 I = 1,DIMEN
              S(I,I) = 1.0D0/MU0
              IF (ELTYP.EQ.2) S(I,I) = 1.0D0/ (PERM*MU0)
 1060     CONTINUE
C
C     Integration loop for element stiffness using NQP quadrature
C     points
C
          CALL MATNUL(ELK,IELK,JELK,DOFEL,DOFEL,ITEST)
          CALL MATNUL(ELM,IELM,JELM,DOFEL,DOFEL,ITEST)
          CALL VECNUL(ELJ,IELJ,DOFEL,ITEST)
          DO 1090 IQUAD = 1,NQP
C
C     Form linear shape function and space derivatives in the local
C     corrdinates.
C
              XI = ABSS(1,IQUAD)
              ETA = ABSS(2,IQUAD)
              CALL QUAM4(FUN,IFUN,LDER,ILDER,JLDER,XI,ETA,ITEST)
              CALL MATMUL(LDER,ILDER,JLDER,GEOM,IGEOM,JGEOM,JAC,IJAC,
     *                    JJAC,DIMEN,NODEL,DIMEN,ITEST)
C
C     Transform local derivatives to global coordinate system and
C     output to work file for later use.
C
              CALL MATINV(JAC,IJAC,JJAC,JACIN,IJACIN,JJACIN,DIMEN,DET,
     *                    ITEST)
              CALL MATMUL(JACIN,IJACIN,JJACIN,LDER,ILDER,JLDER,GDER,
     *                    IGDER,JGDER,DIMEN,DIMEN,NODEL,ITEST)
              WRITE (NTEMP) ((GDER(I,J),I=1,DIMEN),J=1,DOFEL)
C
C     Formation of element stiffness ELK and ELM
C
              CALL MATMUL(S,IS,JS,GDER,IGDER,JGDER,PD,IPD,JPD,DIMEN,
     *                    DIMEN,DOFEL,ITEST)
              CALL MATRAN(GDER,IGDER,JGDER,GDERT,IGDERT,JGDERT,DIMEN,
     *                    DOFEL,ITEST)
              CALL MATMUL(GDERT,IGDERT,JGDERT,PD,IPD,JPD,DTPD,IDTPD,
     *                    JDTPD,DOFEL,DIMEN,DOFEL,ITEST)
              CALL DYAD(FUN,IFUN,FUN,IFUN,FTF,IFTF,JFTF,DOFEL,ITEST)
C
              QUOT = ABS(DET)*WGHT(IQUAD)
              DO 1080 I = 1,DOFEL
                  DO 1070 J = 1,DOFEL
                      DTPD(I,J) = DTPD(I,J)*QUOT
                      FTF(I,J) = FTF(I,J)*CONST*QUOT
 1070             CONTINUE
                  SCVEC(I) = FUN(I)*JVAL(NELE)*QUOT
 1080         CONTINUE
              CALL MATADD(ELK,IELK,JELK,DTPD,IDTPD,JDTPD,DOFEL,DOFEL,
     *                    ITEST)
              CALL MATADD(ELM,IELM,JELM,FTF,IFTF,JFTF,DOFEL,DOFEL,ITEST)
              CALL VECADD(ELJ,IELJ,SCVEC,ISCVEC,DOFEL,ITEST)
 1090     CONTINUE
C
C     Assembly of system stiffness matrices ELK and ELM
C
          CALL DIRECT(NELE,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,STEER,
     *                ISTEER,ITEST)
          CALL RASSYM(SYSK,ISYSK,JSYSK,ELK,IELK,JELK,STEER,ISTEER,HBAND,
     *                DOFEL,ITEST)
          CALL IASSYM(SYSK,ISYSK,JSYSK,ELM,IELM,JELM,STEER,ISTEER,HBAND,
     *                DOFEL,ITEST)
          CALL RASRHS(RHS,IRHS,ELJ,IELJ,STEER,ISTEER,DOFEL,ITEST)
 1100 CONTINUE
C
C     *                           *********************
C     *                           *                   *
C     *                           * Equation Solution *
C     *                           *                   *
C     *                           *********************
C
C     Implement boundary conditions
C
      DO 1110 I = 1,BNDNOD
          J = BNODE(I)
          SYSK(1,J,HBAND) = SYSK(1,J,HBAND)*SCALE
          SYSK(2,J,HBAND) = SYSK(2,J,HBAND)*SCALE
          RHS(1,J) = SYSK(1,J,HBAND)*BVAL(1,I) -
     *               SYSK(2,J,HBAND)*BVAL(2,I)
          RHS(2,J) = SYSK(2,J,HBAND)*BVAL(1,I) +
     *               SYSK(1,J,HBAND)*BVAL(2,I)
 1110 CONTINUE
C
C     Solution of system matrix for the nodal values of the
C     potential
C
      CALL CSYSOL(SYSK,ISYSK,JSYSK,RHS,IRHS,TOTDOF,HBAND,ITEST)
      WRITE (NOUT,FMT=9870)
      CALL CPRTVL(RHS,IRHS,NF,INF,JNF,DOFNOD,TOTNOD,NOUT,ITEST)
C
C     *                           ****************************
C     *                           *                          *
C     *                           *  Recover FIELD Values    *
C     *                           *  From Vector Potential   *
C     *                           *                          *
C     *                           ****************************
C
      REWIND NTEMP
      DO 1140 NELE = 1,TOTELS
          WRITE (NOUT,FMT=9860) NELE
C
C     SELECT nodal potentials for element NELE using the steering
C     vector
C
          CALL DIRECT(NELE,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,STEER,
     *                ISTEER,ITEST)
          DO 1130 IQUAD = 1,NQP
              CALL VECNUL(ELPOT,IELPOT,DOFEL,ITEST)
              DO 1120 I = 1,DOFEL
                  JSTEER = STEER(I)
                  IF (JSTEER.NE.0) ELPOT(I) = SQRT(RHS(1,JSTEER)**2+
     *                RHS(2,JSTEER)**2)
 1120         CONTINUE
C
C     Recover GDER for point IQUAD
C
              READ (NTEMP) ((GDER(I,J),I=1,DIMEN),J=1,DOFEL)
              CALL MATVEC(GDER,IGDER,JGDER,ELPOT,IELPOT,DIMEN,DOFEL,
     *                    FIELD,IFIELD,ITEST)
              BMOD = SQRT(FIELD(1)**2+FIELD(2)**2)
              WRITE (NOUT,FMT=9850) IQUAD
              WRITE (NOUT,FMT=9840) (FIELD(I),I=1,2),BMOD
 1130     CONTINUE
 1140 CONTINUE
C
C     End of program
C
      STOP
C
 9990 FORMAT (16I5)
 9980 FORMAT (I5,6F10.0)
 9970 FORMAT (3F10.0)
 9960 FORMAT (/,/,' **** NODAL GEOMETRY ****',/,' ')
 9950 FORMAT (' ',16I5)
 9940 FORMAT (' ',I5,6F10.5)
 9930 FORMAT (' ',I5,6E10.3)
 9920 FORMAT (/,/,' **** ELEMENT TOPOLOGY ****',/,' ')
 9910 FORMAT (/,/,' ****  FREQ.  COND.  PERM. ****',/,' ')
 9900 FORMAT (' ',F10.3,E10.3,F10.3)
 9890 FORMAT (/,/,' **** BOUNDARY CONDITIONS ****',/,' ')
 9880 FORMAT (/,/,' **** CURRENT SOURCES **** ',/,' ')
 9870 FORMAT (/,/,' **** NODAL POTENTIALS ****',/,' ')
 9860 FORMAT (/,/,' **** FIELD VALUES FOR ELEMENT ',I3,' ****',/,' ',5X,
     *       'BX',12X,'BY',14X,'MODB')
 9850 FORMAT (/,' QUADRATURE POINT ',I3,/,' ')
 9840 FORMAT (' ',2 (D12.5,2X),3X,D12.5)
C
      END
