C***********************************************************************
C     $Id: seg5p2.f,v 1.3 2009/02/27 11:44:50 cg44 Exp $
C
C     Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
C                          Chilton, Didcot, Oxfordshire OX11 0QX
C
C     Contact: Prof Chris Greenough
C     Email: c.greenough@rl.ac.uk
C     Tel: (44) 1235 445307
C
C     Finite Library Program: Seg5p2 - Non-Linear Electrostatics 
C                                      Solution
C
C***********************************************************************
C
      PROGRAM SEG5P2
C
      INTEGER BNDNOD,DIMEN,DOFEL,DOFNOD,ELNUM,ELTOP,HBAND,I,IABSS,
     *        ICOORD,IDTD,IDTDP,IELK,IELQ,IELTOP,IFUN,IGDER,IGDERT,
     *        IGEOM,IJAC,IJACIN,ILDER,INF,IPHI,IQUAD,IRESTR,IRHS,ISTEER,
     *        ISYSK,ITER,ITEST,IWGHT,J,JABSS,JCOORD,JDTD,JELK,JELTOP,
     *        JGDER,JGDERT,JGEOM,JJAC,JJACIN,JLDER,JNF,JRESTR,JSYSK,
     *        MAXIT,NELE,NF,NIN,NODE,NODEL,NODNUM,NOUT,NQP,RESTR,STEER,
     *        TOTDOF,TOTELS,TOTNOD
      DOUBLE PRECISION ABSS,AREA,BOLTZ,BVAL,COORD,DEBYL,DELTA,DET,DMAX,
     *                 DMIN,DOPE,DOPING,DSTOP,DSTRT,DTD,DTDP,ECHRG,ELEC,
     *                 ELK,ELQ,EPSIL,ERROR,ETA,FUN,GDER,GDERT,GEOM,HOLE,
     *                 JAC,JACIN,LDER,NORM,NSUBI,PHI,PHIN,PHIP,QUOT,RHS,
     *                 SYSK,TEMP,VMAX,VMIN,VTHERM,WGHT,X,XI,XSTOP,XSTRT
      DIMENSION ABSS(3,7),DTD(6,6),DTDP(6),ELK(6,6),ELQ(6),FUN(6),
     *          GDER(3,6),GDERT(6,3),GEOM(6,3),JAC(3,3),JACIN(3,3),
     *          LDER(3,6),STEER(6),WGHT(7)
C
C     Problem size dependent arrays
C
      DIMENSION BVAL(30),COORD(100,3),ELTOP(100,10),NF(100,1),PHI(100),
     *          RESTR(30,2),RHS(100),SYSK(100,25)
C
      DATA IABSS/3/,IDTD/6/,IDTDP/6/,IELK/6/,IELQ/6/,IFUN/6/,IGDER/3/,
     *     IGDERT/6/,IGEOM/6/,IJAC/3/,IJACIN/3/,ILDER/3/,ISTEER/6/,
     *     IWGHT/7/,JABSS/7/,JCOORD/3/,JDTD/6/,JELK/6/,JGDER/6/,
     *     JGDERT/3/,JGEOM/3/,JJAC/3/,JJACIN/3/,JLDER/6/,JNF/1/
C
C     Problem size dependent data statements
C
      DATA ICOORD/100/,IELTOP/100/,INF/100/,IPHI/100/,IRESTR/30/,
     *     IRHS/100/,ISYSK/100/,JELTOP/10/,JRESTR/2/,JSYSK/25/
C
      DATA DOFNOD/1/
C
      DATA NIN/5/,NOUT/6/
C
C     Physical constants- assumes input lengths in micrometers
C
      DATA BOLTZ/1.38D-23/,ECHRG/1.6D-19/,EPSIL/1.10D-16/,
     *     NSUBI/1.48D-2/,TEMP/3.00D2/
C
C     Normalisation constants
C
      DEBYL = SQRT(EPSIL/ECHRG*BOLTZ/ECHRG*TEMP/NSUBI)
      VTHERM = BOLTZ*TEMP/ECHRG
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
      WRITE (NOUT,9940)
      READ (NIN,9990) TOTNOD,DIMEN
      WRITE (NOUT,9930) TOTNOD,DIMEN
      DO 1010 I = 1,TOTNOD
          READ (NIN,9980) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
          WRITE (NOUT,9920) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
C
C     Normalise COORD
C
          DO 1000 J = 1,DIMEN
              COORD(I,J) = COORD(I,J)/DEBYL
 1000     CONTINUE
 1010 CONTINUE
C
C     Input of element topology
C
      WRITE (NOUT,9910)
      READ (NIN,9990) TOTELS
      WRITE (NOUT,9930) TOTELS
      DO 1020 I = 1,TOTELS
          READ (NIN,9990) ELNUM,ELTOP(ELNUM,1),NODEL,
     *      (ELTOP(ELNUM,J+2),J=1,NODEL)
          WRITE (NOUT,9930) ELNUM,ELTOP(ELNUM,1),NODEL,
     *      (ELTOP(ELNUM,J+2),J=1,NODEL)
          ELTOP(ELNUM,2) = NODEL
 1020 CONTINUE
C
C     Input of DOPING parameters
C
      WRITE (NOUT,9900)
      READ (NIN,9970) XSTRT,XSTOP
      WRITE (NOUT,9890) XSTRT,XSTOP
      READ (NIN,9970) DSTRT,DSTOP
      WRITE (NOUT,9890) DSTRT,DSTOP
C
C
C     Input of number of degrees of freedom per NODE, input of
C     boundary conditions and construction of nodal freedom array NF
C
      WRITE (NOUT,9880)
      READ (NIN,9990) BNDNOD
      WRITE (NOUT,9930) BNDNOD
      VMAX = 0.D0
      VMIN = 0.D0
      DO 1030 I = 1,BNDNOD
          READ (NIN,9960) RESTR(I,1),BVAL(I)
          WRITE (NOUT,9920) RESTR(I,1),BVAL(I)
          RESTR(I,2) = 1
          BVAL(I) = BVAL(I)/VTHERM
          VMAX = MAX(VMAX,BVAL(I))
          VMIN = MIN(VMIN,BVAL(I))
 1030 CONTINUE
C
C     Max number of iterations and ERROR
C
      WRITE (NOUT,9870)
      READ (NIN,9950) MAXIT,ERROR
      WRITE (NOUT,9920) MAXIT,ERROR
C
C     Create the freedom array
C
      CALL FORMNF(RESTR,IRESTR,JRESTR,BNDNOD,TOTNOD,DOFNOD,NF,INF,JNF,
     *            TOTDOF,ITEST)
C
C     Calculation of semi-bandwidth
C
      CALL BNDWTH(ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,TOTELS,HBAND,
     *            ITEST)
C
C     Calculation of minimum and maximum DOPING
C
      DO 1040 I = 1,TOTNOD
          X = COORD(I,1)
          DOPE = DOPING(X,XSTRT,XSTOP,DSTRT,DSTOP,DEBYL,NSUBI)
          DMIN = MIN(DMIN,DOPE)
          DMAX = MAX(DMAX,DOPE)
 1040 CONTINUE
C
C     Calculation of quasi-fermi potentials
C
      PHIP = VMIN + LOG(-0.5D0*DMIN+SQRT(1.D0+ (0.5D0*DMIN)**2))
      PHIN = VMAX - LOG(0.5D0*DMAX+SQRT(1.D0+ (0.5D0*DMAX)**2))
      DO 1050 I = 1,TOTNOD
          X = COORD(I,1)
          DOPE = DOPING(X,XSTRT,XSTOP,DSTRT,DSTOP,DEBYL,NSUBI)
          IF (DOPE.GT.0.D0) PHI(I) = PHIN +
     *                               LOG(0.5D0*DOPE+SQRT(1.D0+ (0.5D0*
     *                               DOPE)**2))
          IF (DOPE.LT.0.D0) PHI(I) = PHIP -
     *                               LOG(-0.5D0*DOPE+SQRT(1.D0+ (0.5D0*
     *                               DOPE)**2))
          IF (DOPE.EQ.0.D0) PHI(I) = (PHIN+PHIP)/2
 1050 CONTINUE
C
C     Assign the dirichlet boundary values
C
      DO 1060 I = 1,BNDNOD
          J = RESTR(I,1)
          PHI(J) = BVAL(I)
 1060 CONTINUE
C
C     Quadrature weights and abssicae
C
      CALL QTRI4(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST)
C
C     Initialise RHS
C
      CALL VECNUL(RHS,IRHS,TOTDOF,ITEST)
C
C     *                           ***********************************
C     *                           *                                 *
C     *                           *         Iteration Loop          *
C     *                           *                                 *
C     *                           ***********************************
      ITER = 0
 1070 CONTINUE
      ITER = ITER + 1
      WRITE (NOUT,9840) ITER
C
C     *                           ***********************************
C     *                           *                                 *
C     *                           * System Stiffness Matrix Assembly*
C     *                           *                                 *
C     *                           ***********************************
C
      CALL MATNUL(SYSK,ISYSK,JSYSK,TOTDOF,HBAND,ITEST)
      DO 1120 NELE = 1,TOTELS
          NODEL = ELTOP(NELE,2)
          DOFEL = NODEL*DOFNOD
          CALL ELGEOM(NELE,ELTOP,IELTOP,JELTOP,COORD,ICOORD,JCOORD,GEOM,
     *                IGEOM,JGEOM,DIMEN,ITEST)
C
C
          CALL MATNUL(ELK,IELK,JELK,DOFEL,DOFEL,ITEST)
          CALL VECNUL(ELQ,IELQ,DOFEL,ITEST)
C
C     Integration loop for element stiffness using NQP quadrature
C     points
C
          AREA = 0.D0
          DO 1100 IQUAD = 1,NQP
C
C     Form linear shape function and space derivatives in the local
C     corrdinates. Transform local derivatives to global coordinate
C     system
C
              XI = ABSS(1,IQUAD)
              ETA = ABSS(2,IQUAD)
              CALL TRIM3(FUN,IFUN,LDER,ILDER,JLDER,XI,ETA,ITEST)
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
              CALL MATRAN(GDER,IGDER,JGDER,GDERT,IGDERT,JGDERT,DIMEN,
     *                    DOFEL,ITEST)
              CALL MATMUL(GDERT,IGDERT,JGDERT,GDER,IGDER,JGDER,DTD,IDTD,
     *                    JDTD,DOFEL,DIMEN,DOFEL,ITEST)
C
              QUOT = ABS(DET)*WGHT(IQUAD)
              AREA = AREA + QUOT
              DO 1090 I = 1,DOFEL
                  DTDP(I) = 0.D0
                  DO 1080 J = 1,DOFEL
                      DTD(I,J) = DTD(I,J)*QUOT
                      NODE = ELTOP(NELE,J+2)
                      DTDP(I) = DTDP(I) - DTD(I,J)*PHI(NODE)
 1080             CONTINUE
 1090         CONTINUE
C
              CALL MATADD(ELK,IELK,JELK,DTD,IDTD,JDTD,DOFEL,DOFEL,ITEST)
              CALL VECADD(ELQ,IELQ,DTDP,IDTDP,DOFEL,ITEST)
 1100     CONTINUE
C
C     Add in the charge terms (NODEL value * AREA weight)
          DO 1110 I = 1,NODEL
              NODE = ELTOP(NELE,I+2)
              HOLE = EXP(PHIP-PHI(NODE))
              ELEC = EXP(PHI(NODE)-PHIN)
              X = COORD(NODE,1)
              DOPE = DOPING(X,XSTRT,XSTOP,DSTRT,DSTOP,DEBYL,NSUBI)
              ELK(I,I) = ELK(I,I) + (HOLE+ELEC)*AREA/3.D0
              ELQ(I) = ELQ(I) + (DOPE+HOLE-ELEC)*AREA/3.D0
 1110     CONTINUE
C
C     Assembly of system stiffness matrix
C
          CALL DIRECT(NELE,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,STEER,
     *                ISTEER,ITEST)
          CALL ASSYM(SYSK,ISYSK,JSYSK,ELK,IELK,JELK,STEER,ISTEER,HBAND,
     *               DOFEL,ITEST)
          CALL ASRHS(RHS,IRHS,ELQ,IELQ,STEER,ISTEER,DOFEL,ITEST)
 1120 CONTINUE
C
C     *                           *********************
C     *                           *                   *
C     *                           * Equation Solution *
C     *                           *                   *
C     *                           *********************
C
C     Solution of system matrix for the nodal values of the
C     potential
C
      CALL CHOSOL(SYSK,ISYSK,JSYSK,RHS,IRHS,TOTDOF,HBAND,ITEST)
C
C     *                           ***************************
C     *                           *                         *
C     *                           * UPDATE Of Potential And *
C     *                           *  Convergence Checking   *
C     *                           *                         *
C     *                           ***************************
C
      CALL UPDATE(PHI,IPHI,RHS,IRHS,TOTNOD,DOFNOD,TOTDOF,NF,INF,JNF,
     *            ITEST)
      DELTA = NORM(RHS,IRHS,TOTDOF,ITEST)*VTHERM
C
C     Check convergence and number of iterations
C
      IF (DELTA.GT.ERROR .AND. ITER.LT.MAXIT) GO TO 1070
C
      WRITE (NOUT,9860) ITER
      DO 1130 I = 1,TOTNOD
          PHI(I) = PHI(I)*VTHERM
 1130 CONTINUE
C
      WRITE (NOUT,9850)
      CALL PRTVEC(PHI,IPHI,TOTNOD,NOUT,ITEST)
C
C     End of Program
C
      STOP
C
 9990 FORMAT (16I5)
 9980 FORMAT (I5,6F10.0)
 9970 FORMAT (2F10.0)
 9960 FORMAT (I5,6F10.0)
 9950 FORMAT (I5,F10.0)
 9940 FORMAT (/,/,' **** NODAL GEOMETRY ****',/,/,' ')
 9930 FORMAT (' ',16I5)
 9920 FORMAT (' ',I5,6F10.5)
 9910 FORMAT (/,/,' **** ELEMENT TOPOLOGY ****',/,/,' ')
 9900 FORMAT (/,/,' **** DOPING ORIGIN AND GRADIENT ****',/,/,' ')
 9890 FORMAT (' ',2E10.3)
 9880 FORMAT (/,/,' **** BOUNDARY CONDITIONS ****',/,/,' ')
 9870 FORMAT (/,/,' **** MAXIMUM ITERATIONS AND ERROR ****',/,/,' ')
 9860 FORMAT (/,/,'TOTAL ITERATIONS',I5)
 9850 FORMAT (/,/,' **** NODAL POTENTIALS ****',/,/,' ')
 9840 FORMAT (' Iteration : ',I4)
C
      END
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION DOPING(X,XSTRT,XSTOP,DSTRT,DSTOP,DEBYL,
     *                 NSUBI)
C
C-----------------------------------------------------------------------
C PURPOSE
C      The function DOPING is to calculate the impurity concentrations
C      at nodes. The vaules are scaled on return.
C
C HISTORY
C      Release 3.0  24 October 1985 (CJH)
C      Release 4.0   2 December 2003 (CG)
C
C     DOUBLE PRECISION FUNCTION DOPING(X,XSTRT,XSTOP,DSTRT,DSTOP,DEBYL,
C    *                 NSUBI)
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION DEBYL,DSTOP,DSTRT,NSUBI,X,XSTOP,XSTRT
C
      X = X*DEBYL
      IF (X.LE.XSTRT) DOPING = DSTRT
      IF (X.GE.XSTOP) DOPING = DSTOP
      IF (X.GT.XSTRT .AND. X.LT.XSTOP) DOPING = DSTRT +
     *    (DSTOP-DSTRT)* (X-XSTRT)/ (XSTOP-XSTRT)
      DOPING = DOPING/NSUBI
      X = X/DEBYL
C
      END
