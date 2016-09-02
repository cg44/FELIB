C***********************************************************************
C     $Id: seg5p1.f,v 1.2 2009/02/27 11:44:50 cg44 Exp $
C
C     Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
C                          Chilton, Didcot, Oxfordshire OX11 0QX
C
C     Contact: Prof Chris Greenough
C     Email: c.greenough@rl.ac.uk
C     Tel: (44) 1235 445307
C
C     Finite Library Program: Seg5p1 - Non-Linear Magnetostatics 
C                                      Solution
C
C***********************************************************************
C
      PROGRAM SEG5P1
C
      INTEGER BNDNOD,BNODE,DIMEN,DOFEL,DOFNOD,ELNUM,ELTOP,ELTYP,HBAND,I,
     *        IABSS,IB,IBTB,ICOORD,IDTSD,IELJ,IELK,IELPOT,IELT,IELTOP,
     *        IFUN,IGDER,IGDERT,IGEOM,IJAC,IJACIN,IJVAL,ILDER,ILOADS,
     *        INF,IQUAD,IRHS,IS,ISCVEC,ISD,ISTEER,ISYSK,ISYST,ITER,
     *        ITEST,IVEC,IWGHT,IWORK,J,JABSS,JBTB,JCOORD,JDTSD,JELK,
     *        JELT,JELTOP,JGDER,JGDERT,JGEOM,JJAC,JJACIN,JLDER,JNF,JS,
     *        JSD,JSYSK,JSYST,NELE,NF,NIN,MAXIT,NJV,NODEL,NODNUM,NOUT,
     *        NQP,STEER,TOTDOF,TOTELS,TOTNOD,IBNODE
      DOUBLE PRECISION ABSS,B,B2,BMOD,BTB,BVAL,CFACTR,COORD,DET,DNU,
     *                 DTSD,DUMMY,ELJ,ELK,ELPOT,ELT,EPS,ETA,FUN,GDER,
     *                 GDERT,GEOM,JAC,JACIN,JVAL,LDER,LOADS,MAX,NU,NU0,
     *                 PI,QUOT,RHS,S,SCALE,SCVEC,SD,SYSK,SYST,VEC,VEPS,
     *                 WGHT,WORK,XI
      DIMENSION ABSS(3,9),B(8),BTB(8,8),DTSD(8,8),ELJ(8),ELK(8,8),
     *          ELPOT(8),ELT(8,8),FUN(8),GDER(3,8),GDERT(8,3),GEOM(8,3),
     *          JAC(3,3),JACIN(3,3),LDER(3,8),S(3,3),SCVEC(8),SD(3,8),
     *          STEER(8),WGHT(9),WORK(8)
C
C     Problem size dependent arrays
C
      DIMENSION BNODE(30),BVAL(30),COORD(100,3),ELTOP(100,10),JVAL(100),
     *          LOADS(100),NF(100,1),RHS(100),SYSK(100,25),SYST(100,25),
     *          VEC(100)
C
      DATA IABSS/3/,IB/8/,IBTB/8/,IDTSD/8/,IELJ/8/,IELK/8/,IELPOT/8/,
     *     IELT/8/,IFUN/8/,IGDER/3/,IGDERT/8/,IGEOM/8/,IJAC/3/,
     *     IJACIN/3/,ILDER/3/,IS/3/,ISCVEC/8/,ISD/3/,ISTEER/8/,IWGHT/9/,
     *     IWORK/8/,JABSS/9/,JBTB/8/,JCOORD/3/,JDTSD/8/,JELK/8/,JELT/8/,
     *     JGDER/8/,JGDERT/3/,JGEOM/3/,JJAC/3/,JJACIN/3/,JLDER/8/,
     *     JNF/1/,JS/3/,JSD/8/,SCALE/1.0D+10/
C
C     Problem size dependent data statements
C
      DATA ICOORD/100/,IELTOP/100/,IJVAL/100/,ILOADS/100/,INF/100/,
     *     IRHS/100/,ISYSK/100/,ISYST/100/,IVEC/100/,JELTOP/10/,
     *     JSYSK/25/,JSYST/25/,IBNODE/30/
C
      DATA NIN/5/,NOUT/6/
C
C     Set ITEST for full checking
C
      ITEST = 0
C
C     Set values of PI and NU0
C
      PI = ATAN(1.0D0)*4.0D0
      NU0 = 1.0D0/ (4.0D0*PI*1.0D-07)
      EPS = VEPS(DUMMY)
C
C     *                           **********************
C     *                           *                    *
C     *                           * Input Data Section *
C     *                           *                    *
C     *                           **********************
C
C     Input of nodal geometry
C
      WRITE (NOUT,9960)
      READ (NIN,9990) TOTNOD,DIMEN
      WRITE (NOUT,9950) TOTNOD,DIMEN
C
      DO 1000 I = 1,TOTNOD
          READ (NIN,9980) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
          WRITE (NOUT,9940) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
 1000 CONTINUE
C
C     Input of element topology
C
      WRITE (NOUT,9930)
      READ (NIN,9990) TOTELS
      WRITE (NOUT,9950) TOTELS
C
      DO 1010 I = 1,TOTELS
          READ (NIN,9990) ELNUM,ELTYP,NODEL,
     *      (ELTOP(ELNUM,J+2),J=1,NODEL)
          WRITE (NOUT,9950) ELNUM,ELTYP,NODEL,
     *      (ELTOP(ELNUM,J+2),J=1,NODEL)
          ELTOP(ELNUM,1) = ELTYP
          ELTOP(ELNUM,2) = NODEL
 1010 CONTINUE
C
C     Input of number of degrees of freedom per node, input of
C     boundary conditions and construction of nodal freedom array NF
C
      WRITE (NOUT,9920)
      READ (NIN,9990) DOFNOD
      WRITE (NOUT,9950) DOFNOD
C
      READ (NIN,9990) BNDNOD
      WRITE (NOUT,9950) BNDNOD
C
      DO 1020 I = 1,BNDNOD
          READ (NIN,9980) BNODE(I),BVAL(I)
          WRITE (NOUT,9940) BNODE(I),BVAL(I)
 1020 CONTINUE
C
C     Input number of current carrying elements, element number, and
C     value of current
C
      WRITE (NOUT,9910)
      READ (NIN,9990) NJV
      WRITE (NOUT,9950) NJV
      CALL VECNUL(JVAL,IJVAL,TOTELS,ITEST)
C
      IF (NJV.NE.0) THEN
          DO 1030 I = 1,NJV
              READ (NIN,9970) ELNUM,JVAL(ELNUM)
              WRITE (NOUT,9900) ELNUM,JVAL(ELNUM)
 1030     CONTINUE
      END IF
C
C     Input iteration control parameters
C
      READ (NIN,9970) MAXIT,CFACTR
      WRITE (NOUT,9890) MAXIT,CFACTR
C
C     Form nodal freedom array
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
C     *                           *  System Stiffness And Tangent    *
C     *                           *          Matrix Assembly         *
C     *                           *                                  *
C     *                           ************************************
C
      CALL VECNUL(RHS,IRHS,TOTDOF,ITEST)
      CALL QTRI4(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST)
      DOFEL = NODEL*DOFNOD
C
C     Non-linear iteration starts here
C
      ITER = 0
 1060 CONTINUE
      ITER = ITER + 1
C
      CALL MATNUL(SYSK,ISYSK,JSYSK,TOTDOF,HBAND,ITEST)
      CALL MATNUL(SYST,ISYST,JSYST,TOTDOF,HBAND,ITEST)
      CALL VECNUL(LOADS,ILOADS,TOTDOF,ITEST)
      DO 1110 NELE = 1,TOTELS
          ELTYP = ELTOP(NELE,1)
          CALL ELGEOM(NELE,ELTOP,IELTOP,JELTOP,COORD,ICOORD,JCOORD,GEOM,
     *                IGEOM,JGEOM,DIMEN,ITEST)
C
C     Initialise ELPOT nodal values of potential
C
          CALL DIRECT(NELE,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,STEER,
     *                ISTEER,ITEST)
          CALL SELECT(RHS,IRHS,STEER,ISTEER,DOFEL,ELPOT,IELPOT,ITEST)
C
C     Integration loop for element stiffness using NQP quadrature
C     points
C
          CALL MATNUL(ELK,IELK,JELK,DOFEL,DOFEL,ITEST)
          CALL MATNUL(ELT,IELT,JELT,DOFEL,DOFEL,ITEST)
          CALL VECNUL(ELJ,IELJ,DOFEL,ITEST)
          CALL MATNUL(S,IS,JS,DIMEN,DIMEN,ITEST)
          DO 1100 IQUAD = 1,NQP
C
C     Form linear shape function and space derivatives in the local
C     corrdinates. Transform local derivatives to global coordinate
C     system
C
              XI = ABSS(1,IQUAD)
              ETA = ABSS(2,IQUAD)
              CALL TRIM3(FUN,IFUN,LDER,ILDER,JLDER,XI,ETA,ITEST)
              CALL MATMUL(LDER,ILDER,JLDER,GEOM,IGEOM,JGEOM,JAC,IJAC,
     *                    JJAC,DIMEN,NODEL,DIMEN,ITEST)
              CALL MATINV(JAC,IJAC,JJAC,JACIN,IJACIN,JJACIN,DIMEN,DET,
     *                    ITEST)
              CALL MATMUL(JACIN,IJACIN,JJACIN,LDER,ILDER,JLDER,GDER,
     *                    IGDER,JGDER,DIMEN,DIMEN,NODEL,ITEST)
C
C     Generate magnetic flux from previous solution, with modulus in
C     BMOD
C
              CALL MATVEC(GDER,IGDER,JGDER,ELPOT,IELPOT,DIMEN,DOFEL,
     *                    WORK,IWORK,ITEST)
              CALL SCAPRD(WORK,IWORK,WORK,IWORK,DIMEN,B2,ITEST)
              BMOD = SQRT(B2)
C
C     Formation of element tangent ELT
C
              CALL VECMAT(WORK,IWORK,GDER,IGDER,JGDER,DIMEN,DOFEL,B,IB,
     *                    ITEST)
              CALL DYAD(B,IB,B,IB,BTB,IBTB,JBTB,DOFEL,ITEST)
C
C     Formation of element stiffness ELK obtain NU and DNU for
C     element
C
              NU = NU0
              DNU = 0.0D0
              IF (ELTYP.EQ.2) CALL BHCURV(BMOD,NU,DNU)
              DO 1070 I = 1,DIMEN
                  S(I,I) = NU
 1070         CONTINUE
              CALL MATMUL(S,IS,JS,GDER,IGDER,JGDER,SD,ISD,JSD,DIMEN,
     *                    DIMEN,DOFEL,ITEST)
              CALL MATRAN(GDER,IGDER,JGDER,GDERT,IGDERT,JGDERT,DIMEN,
     *                    DOFEL,ITEST)
              CALL MATMUL(GDERT,IGDERT,JGDERT,SD,ISD,JSD,DTSD,IDTSD,
     *                    JDTSD,DOFEL,DIMEN,DOFEL,ITEST)
C
              QUOT = ABS(DET)*WGHT(IQUAD)
              DO 1090 I = 1,DOFEL
                  DO 1080 J = 1,DOFEL
                      DTSD(I,J) = DTSD(I,J)*QUOT
                      IF (BMOD.LT.EPS) BTB(I,J) = 0.0D0
                      IF (BMOD.GT.EPS) BTB(I,J) = BTB(I,J)*QUOT*DNU/BMOD
 1080             CONTINUE
                  SCVEC(I) = -FUN(I)*JVAL(NELE)*QUOT
 1090         CONTINUE
C
              CALL MATADD(ELK,IELK,JELK,DTSD,IDTSD,JDTSD,DOFEL,DOFEL,
     *                    ITEST)
              CALL MATADD(ELT,IELT,JELT,BTB,IBTB,JBTB,DOFEL,DOFEL,ITEST)
              CALL VECADD(ELJ,IELJ,SCVEC,ISCVEC,DOFEL,ITEST)
 1100     CONTINUE
C
C     Assembly of system stiffness and tangent matrix, and source
C     vector
C
          CALL DIRECT(NELE,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,STEER,
     *                ISTEER,ITEST)
          CALL ASSYM(SYSK,ISYSK,JSYSK,ELK,IELK,JELK,STEER,ISTEER,HBAND,
     *               DOFEL,ITEST)
          CALL ASSYM(SYST,ISYST,JSYST,ELT,IELT,JELT,STEER,ISTEER,HBAND,
     *               DOFEL,ITEST)
          CALL ASSYM(SYST,ISYST,JSYST,ELK,IELK,JELK,STEER,ISTEER,HBAND,
     *               DOFEL,ITEST)
          CALL ASRHS(LOADS,ILOADS,ELJ,IELJ,STEER,ISTEER,DOFEL,ITEST)
 1110 CONTINUE
C
C     *                           *********************
C     *                           *                   *
C     *                           * Equation Solution *
C     *                           *                   *
C     *                           *********************
C
      CALL MVSYB(SYSK,ISYSK,JSYSK,RHS,IRHS,VEC,IVEC,TOTDOF,HBAND,ITEST)
      CALL VECADD(VEC,IVEC,LOADS,ILOADS,TOTDOF,ITEST)
C
C     Modify tangent matrix and residual to implement boundary
C     conditions
C
      DO 1120 I = 1,BNDNOD
          J = BNODE(I)
          SYST(J,HBAND) = SYST(J,HBAND)*SCALE
          VEC(J) = VEC(J) + SYST(J,HBAND)* (RHS(J)-BVAL(I))
 1120 CONTINUE
C
C     Form LU decomposition and solve for UPDATE
C
      CALL CHOSOL(SYST,ISYST,JSYST,VEC,IVEC,TOTDOF,HBAND,ITEST)
      CALL VECSUB(RHS,IRHS,VEC,IVEC,TOTDOF,ITEST)
      WRITE (NOUT,9880) ITER
      CALL PRTVEC(RHS,IRHS,TOTDOF,NOUT,ITEST)
C
C     Check for convergence and number of iterations performed
C
      MAX = 0.0D0
      DO 1130 I = 1,TOTDOF
          IF (ABS(VEC(I)).GT.MAX) MAX = ABS(VEC(I))
 1130 CONTINUE
C
      IF (MAX.GT.CFACTR .AND. ITER.LT.MAXIT) GO TO 1060
C
C     End of Program
C
      STOP
C
 9990 FORMAT (16I5)
 9980 FORMAT (I5,6F10.0)
 9970 FORMAT (I5,D10.3)
 9960 FORMAT (/,/,' **** NODAL GEOMETRY ****',/,' ')
 9950 FORMAT (' ',16I5)
 9940 FORMAT (' ',I5,6F10.5)
 9930 FORMAT (/,/,' **** ELEMENT TOPOLOGY ****',/,' ')
 9920 FORMAT (/,/,' **** BOUNDARY CONDITIONS ****',/,' ')
 9910 FORMAT (/,/,' **** CURRENT SOURCES **** ',/,' ')
 9900 FORMAT (' ',I5,D15.6)
 9890 FORMAT (/,/,' MAXIMUM NUMBER OF ITERATIONS ',I5,/,
     *       ' CONVERGENCE FACTOR ',D12.3,/)
 9880 FORMAT (/,/,' **** NODAL POTENTIALS ****',10X,' ITERATION',' NO.',
     *       I5,/,' ')
C
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BHCURV(BMOD,NU,DNU)
C
C-----------------------------------------------------------------------
C PURPOSE
C      The subroutine BHCURV returns values of NU and DNU for
C      a given value of BMOD
C
C METHOD
C      NU vs. BMOD is an hyperbolic tangent function,
C      calculated from:
C      1)  relative NU at BMOD=0 IS set to 1.0D-03
C      2)  when BMOD=2, NU=0.99D0 and will tend to 1.0D0
C          as BMOD increases
C
C HISTORY
C      Release 3.0  24 October 1984 (CRIE)
C      Relase  4.0   2 December 2003 (CG)
C
C ARGUMENTS in
C      BMOD    modulus of magnetic flux
C
C ARGUMENTS out
C      NU      inverse permeability
C      DNU     rate of change of NU with respect to BMOD
C
C      SUBROUTINE BHCURV(BMOD, NU, DNU)
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION A,ATANH,BMOD,C,DNU,NU,NU0,PI
C
      INTRINSIC COSH,TANH
C
      PI = ATAN(1.0D0)*4.0D0
      NU0 = 1.0D0/ (4.0D0*PI*1.0D-07)
      C = TANH(1.0D-3)
      A = 0.5D0* (ATANH(0.99D0)-C)
      NU = TANH(A*BMOD+C)*NU0
      DNU = COSH(A*BMOD+C)
      DNU = A/DNU/DNU*NU0
C
      END
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION ATANH(X)
C
C-----------------------------------------------------------------------
C PURPOSE
C      The function ATANH returns the inverse hyperbolic tangent
C
C HISTORY
C      Release 3.0  24 October 1984 (CRIE)
C      Release 4.0   2 December 2003 (CG)
C
C      DOUBLE PRECISION FUNCTION ATANH(X)
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION X
C
      ATANH = 0.5D0*LOG((1.0D0+X)/ (1.0D0-X))
C
      END
