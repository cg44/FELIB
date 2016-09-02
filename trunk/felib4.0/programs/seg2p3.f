C***********************************************************************
C     $Id: seg2p3.f,v 1.2 2009/02/27 11:44:50 cg44 Exp $
C
C     Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
C                          Chilton, Didcot, Oxfordshire OX11 0QX
C
C     Contact: Prof Chris Greenough
C     Email: c.greenough@rl.ac.uk
C     Tel: (44) 1235 445307
C
C     Finite Library Program: Seg2p3 - Forced Vibration of an Elasic 
C                                      Solid by Modal Superposition
C
C***********************************************************************
C
      PROGRAM SEG2P3
C
      INTEGER DIMEN,DOFEL,DOFNOD,ELNUM,ELTOP,ELTYP,I,IABSS,IB,IBT,IBTDB,
     *        ICOORD,ID,IDB,IDIAG,IELK,IELM,IELTOP,IFUN,IGDER,IGEOM,
     *        IJAC,IJACIN,ILDER,INF,INTN,IQUAD,IRESTR,ISHP,ISTEER,ISTEP,
     *        ISUB,ISYSK,ISYSM,ITEST,ITSHP,IWGHT,IWORK,J,JABSS,JB,JBT,
     *        JBTDB,JCOORD,JD,JDB,JELK,JELM,JELTOP,JGDER,JGEOM,JJAC,
     *        JJACIN,JLDER,JNF,JNTN,JRESTR,JSHP,JSYSK,JTSHP,K,L,LODFRE,
     *        NELE,NF,NIN,NMODES,NODEL,NODNUM,NOUT,NQP,NSTEPS,NUMLOD,
     *        NUMSS,NWORK,RESNOD,RESTR,STEER,TOTDOF,TOTELS,TOTLOD,TOTNOD
      DOUBLE PRECISION ABSS,AMP,AREA,B,BT,BTDB,COORD,D,DAMPR,DB,DET,
     *                 DIAG,DTIM,E,ELK,ELM,EPS,ETA,FUN,GDER,GENSOL,GEOM,
     *                 JAC,JACIN,LDER,NTN,NU,OMEGA,P,QUOT,RHO,SHP,SUB,
     *                 SYSK,SYSM,TIME,TOL,TSHP,VEPS,VTOL,WGHT,WORK,X,XI,
     *                 XMODE
      DIMENSION ABSS(3,9),B(6,24),BT(24,6),BTDB(24,24),D(6,6),DB(6,24),
     *          ELK(24,24),ELM(24,24),FUN(8),GDER(3,8),GEOM(8,3),
     *          JAC(3,3),JACIN(3,3),LDER(3,8),LODFRE(20),NTN(24,24),
     *          SHP(3,24),STEER(24),TSHP(24,3),WGHT(9)
C
C     Problem size dependent arrays
C
      DIMENSION COORD(100,3),DIAG(100),ELTOP(100,10),NF(100,3),NWORK(3),
     *          RESTR(20,4),SUB(100),SYSK(100,100),SYSM(100),WORK(100),
     *          XMODE(100)
C
      DATA IABSS/3/,IB/6/,IBT/24/,IBTDB/24/,ID/6/,IDB/6/,IELK/24/,
     *     IELM/24/,IFUN/8/,IGDER/3/,IGEOM/8/,IJAC/3/,IJACIN/3/,
     *     ILDER/3/,INTN/24/,ISHP/3/,ISTEER/24/,ITSHP/24/,IWGHT/9/,
     *     JABSS/9/,JB/24/,JBT/6/,JBTDB/8/,JD/6/,JDB/24/,JELK/24/,
     *     JELM/24/,JGDER/8/,JGEOM/3/,JJAC/3/,JJACIN/3/,JLDER/8/,
     *     JNTN/24/,JSHP/24/,JTSHP/3/
C
C     Problem size dependent data statements
C
C
      DATA ICOORD/100/,IDIAG/100/,IELTOP/100/,INF/100/,IRESTR/20/,
     *     ISUB/100/,ISYSK/100/,ISYSM/100/,IWORK/100/,JCOORD/3/,
     *     JELTOP/10/,JNF/3/,JRESTR/4/,JSYSK/100/
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
C     Input of material properties and construction of stress-strain
C     matrix D for plane strain
C
      WRITE (NOUT,9910)
      READ (NIN,9970) NU,E,RHO
      WRITE (NOUT,9900) NU,E,RHO
      CALL DPSN(D,ID,JD,E,NU,NUMSS,ITEST)
C
C     Input of number of degrees of freedom per node, input of
C     restrained node data and construction of nodal freedom 
C     array NF
C
      WRITE (NOUT,9890)
      READ (NIN,9990) DOFNOD
      WRITE (NOUT,9940) DOFNOD
      READ (NIN,9990) RESNOD
      WRITE (NOUT,9940) RESNOD
      K = DOFNOD + 1
      DO 1020 I = 1,RESNOD
          READ (NIN,9990) (RESTR(I,J),J=1,K)
          WRITE (NOUT,9940) (RESTR(I,J),J=1,K)
 1020 CONTINUE
      CALL FORMNF(RESTR,IRESTR,JRESTR,RESNOD,TOTNOD,DOFNOD,NF,INF,JNF,
     *            TOTDOF,ITEST)
C
C     Input number of modes and damping ratio
C
      WRITE (NOUT,9880)
      READ (NIN,9980) NMODES,DAMPR
      WRITE (NOUT,9930) NMODES,DAMPR
C
C     Input TIME step length and number of steps
C
      WRITE (NOUT,9870)
      READ (NIN,9960) DTIM,NSTEPS
      WRITE (NOUT,9860) DTIM,NSTEPS
C
C     Input loading data
C
      WRITE (NOUT,9850)
      READ (NIN,9990) NUMLOD
      WRITE (NOUT,9940) NUMLOD
      TOTLOD = 0
      DO 1040 I = 1,NUMLOD
          READ (NIN,9990) NODNUM, (NWORK(J),J=1,DOFNOD)
          WRITE (NOUT,9940) NODNUM, (NWORK(J),J=1,DOFNOD)
          DO 1030 J = 1,DOFNOD
              IF (NWORK(J).NE.0) THEN
                  K = NF(NODNUM,J)
                  IF (K.NE.0) THEN
                      TOTLOD = TOTLOD + 1
                      LODFRE(TOTLOD) = K
                  END IF
              END IF
 1030     CONTINUE
 1040 CONTINUE
C
C     Input forcing function data
C
      WRITE (NOUT,9840)
      READ (NIN,9970) AMP,OMEGA
      WRITE (NOUT,9900) AMP,OMEGA
C
C     *                           *******************************
C     *                           *                             *
C     *                           * System Stiffness Matrix And *
C     *                           *     Mass Matrix Assembly    *
C     *                           *                             *
C     *                           *******************************
C
      CALL MATNUL(SYSK,ISYSK,JSYSK,TOTDOF,TOTDOF,ITEST)
      CALL VECNUL(SYSM,ISYSM,TOTDOF,ITEST)
      DOFEL = NODEL*DOFNOD
      CALL QQUA4(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST)
      DO 1080 NELE = 1,TOTELS
          AREA = 0.0D0
          CALL ELGEOM(NELE,ELTOP,IELTOP,JELTOP,COORD,ICOORD,JCOORD,GEOM,
     *                IGEOM,JGEOM,DIMEN,ITEST)
C
C     Integration loop for element stiffness and element lumped mass
C     using NQP quadrature points
C
          CALL MATNUL(ELK,IELK,JELK,DOFEL,DOFEL,ITEST)
          CALL MATNUL(ELM,IELM,JELM,DOFEL,DOFEL,ITEST)
          DO 1070 IQUAD = 1,NQP
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
C     Formation of strain-displacement matrix B and formation of
C     integrand for element stiffness matrix ELK
C
              CALL B2C2(B,IB,JB,GDER,IGDER,JGDER,NODEL,ITEST)
              CALL MATMUL(D,ID,JD,B,IB,JB,DB,IDB,JDB,NUMSS,NUMSS,DOFEL,
     *                    ITEST)
              CALL MATRAN(B,IB,JB,BT,IBT,JBT,NUMSS,DOFEL,ITEST)
              CALL MATMUL(BT,IBT,JBT,DB,IDB,JDB,BTDB,IBTDB,JBTDB,DOFEL,
     *                    NUMSS,DOFEL,ITEST)
C
C     Formation of integrand for element mass matrix ELM
C
              CALL SHAPFN(SHP,ISHP,JSHP,FUN,IFUN,NODEL,DOFNOD,ITEST)
              CALL MATRAN(SHP,ISHP,JSHP,TSHP,ITSHP,JTSHP,DOFNOD,DOFEL,
     *                    ITEST)
              CALL MATMUL(TSHP,ITSHP,JTSHP,SHP,ISHP,JSHP,NTN,INTN,JNTN,
     *                    DOFEL,DOFNOD,DOFEL,ITEST)
C
              QUOT = ABS(DET)*WGHT(IQUAD)
              AREA = AREA + QUOT
              DO 1060 I = 1,DOFEL
                  DO 1050 J = 1,DOFEL
                      BTDB(I,J) = BTDB(I,J)*QUOT
                      NTN(I,J) = NTN(I,J)*QUOT*RHO
 1050             CONTINUE
 1060         CONTINUE
              CALL MATADD(ELK,IELK,JELK,BTDB,IBTDB,JBTDB,DOFEL,DOFEL,
     *                    ITEST)
              CALL MATADD(ELM,IELM,JELM,NTN,INTN,JNTN,DOFEL,DOFEL,ITEST)
 1070     CONTINUE
C
C     Assembly of system stiffness matrix SYSK, system lumped mass
C     matrix SYSM
C
          CALL DIRECT(NELE,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,STEER,
     *                ISTEER,ITEST)
          CALL ASFUL(SYSK,ISYSK,JSYSK,ELK,IELK,JELK,STEER,ISTEER,DOFEL,
     *               ITEST)
          CALL ASLMS(SYSM,ISYSM,ELM,IELM,JELM,STEER,ISTEER,DOFEL,DOFNOD,
     *               AREA,ITEST)
 1080 CONTINUE
C
C     *                           *********************
C     *                           *                   *
C     *                           * Equation Solution *
C     *                           *                   *
C     *                           *********************
C
C     Reduction of system mass matrix by choleski reduction and
C     construction of modified system stiffness matrix
C
      DO 1090 I = 1,TOTDOF
          SYSM(I) = 1.0D0/SQRT(SYSM(I))
 1090 CONTINUE
      DO 1110 I = 1,TOTDOF
          DO 1100 J = 1,TOTDOF
              SYSK(I,J) = SYSM(I)*SYSK(I,J)*SYSM(J)
 1100     CONTINUE
 1110 CONTINUE
C
C     Householder tridiagonalisation of SYSK and eigenvalue and
C     eigenvector recovery using ql transformations
C
      TOL = VTOL(X)
      CALL HOUSE(SYSK,ISYSK,JSYSK,SYSK,ISYSK,JSYSK,DIAG,IDIAG,SUB,ISUB,
     *           TOTDOF,TOL,ITEST)
      EPS = VEPS(X)
      CALL QLVEC(DIAG,IDIAG,SUB,ISUB,SYSK,ISYSK,JSYSK,TOTDOF,EPS,ITEST)
      WRITE (NOUT,9830)
      CALL PRTVEC(DIAG,IDIAG,TOTDOF,NOUT,ITEST)
      TIME = 0.D0
      DO 1170 ISTEP = 1,NSTEPS
          TIME = TIME + DTIM
          WRITE (NOUT,9820) ISTEP,TIME
          DO 1140 J = 1,NMODES
              DO 1120 I = 1,TOTDOF
                  WORK(I) = SYSK(I,J)*SYSM(I)
 1120         CONTINUE
              P = 0.0D0
              DO 1130 K = 1,TOTLOD
                  L = LODFRE(K)
                  P = P + WORK(L)
 1130         CONTINUE
              XMODE(J) = GENSOL(TIME,P,DAMPR,DIAG(J),AMP,OMEGA)
 1140     CONTINUE
          CALL VECNUL(WORK,IWORK,TOTDOF,ITEST)
          DO 1160 I = 1,TOTDOF
              DO 1150 J = 1,NMODES
                  WORK(I) = WORK(I) + SYSK(I,J)*SYSM(I)*XMODE(J)
 1150         CONTINUE
 1160     CONTINUE
          CALL PRTVAL(WORK,IWORK,NF,INF,JNF,DOFNOD,TOTNOD,NOUT,ITEST)
 1170 CONTINUE
C
C     End of Program
C
      STOP
C
 9990 FORMAT (16I5)
 9980 FORMAT (I5,6F10.0)
 9970 FORMAT (6F10.0)
 9960 FORMAT (F10.0,I5)
 9950 FORMAT (/,/,' **** NODAL GEOMETRY ****',/,' ')
 9940 FORMAT (' ',16I5)
 9930 FORMAT (' ',I5,6F10.5)
 9920 FORMAT (/,/,' **** ELEMENT TOPOLOGY ****',/,' ')
 9910 FORMAT (/,/,' **** MATERIAL PROPERTIES ****',/,' ')
 9900 FORMAT (' ',F10.5,E12.5,F10.5)
 9890 FORMAT (/,/,' **** RESTRAINT DATA ****',/,' ')
 9880 FORMAT (/,/,' **** NUMBER OF MODES AND DAMPING RATIO ****',/,' ')
 9870 FORMAT (/,/,' **** TIME STEP AND NUMBER OF STEPS ****',/,' ')
 9860 FORMAT (' ',F10.5,I5)
 9850 FORMAT (/,/,' **** LOADING DATA ****',/,' ')
 9840 FORMAT (/,/,' **** FORCING FUNCTION DATA ****',/,' ')
 9830 FORMAT (/,/,' **** EIGENVALUES ****',/,' ')
 9820 FORMAT (/,/,' **** DISPLACEMENTS AT STEP',I5,2X,'TIME =',F10.5,
     *       ' ****')
C
      END
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION GENSOL(TIME,P,DAMPR,EIGVAL,AMP,OMEGA)
C
C-----------------------------------------------------------------------
C
C
      DOUBLE PRECISION AMP,COEF1,COEF2,DAMPR,DENOM,EIGVAL,OMEGA,P,TIME
C
      DENOM = (EIGVAL-OMEGA**2)**2 + 4.0D0*DAMPR**2*OMEGA**2*EIGVAL
      COEF1 = P* (EIGVAL-OMEGA**2)/DENOM
      COEF2 = P*2.0D0*OMEGA*SQRT(EIGVAL)*DAMPR/DENOM
      GENSOL = COEF1*COS(OMEGA*TIME) + COEF2*SIN(OMEGA*TIME)
C
      END
