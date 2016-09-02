C***********************************************************************
C     $Id: seg2p2.f,v 1.2 2009/02/27 11:44:50 cg44 Exp $
C
C     Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
C                          Chilton, Didcot, Oxfordshire OX11 0QX
C
C     Contact: Prof Chris Greenough
C     Email: c.greenough@rl.ac.uk
C     Tel: (44) 1235 445307
C
C     Finite Library Program: Seg2p2 - Forced Vibtration of an Elastic
C                                      Solid by Direct Integration
C
C***********************************************************************
C
      PROGRAM SEG2P2
C
      INTEGER DIF,DIMEN,DOFEL,DOFNOD,ELNUM,ELTOP,ELTYP,FREDOM,HBAND,I,
     *        IA0,IA1,IABSS,IB,IBT,IBTDB,ICOORD,ID,ID2A0,ID2A1,IDA0,
     *        IDA1,IDB,IELK,IELM,IELTOP,IFUN,IGDER,IGEOM,IJAC,IJACIN,
     *        ILDER,ILOADS,INF,INTN,IQUAD,IRESTR,ISHP,ISTEER,ISYSK,
     *        ISYSM,ISYSW,ITEST,ITSHP,IWGHT,IWORK,J,JABSS,JB,JBT,JBTDB,
     *        JCOORD,JD,JDB,JELK,JELM,JELTOP,JGDER,JGEOM,JJAC,JJACIN,
     *        JLDER,JNF,JNTN,JRESTR,JSHP,JSYSK,JSYSM,JSYSW,JTSHP,K,
     *        LODFRE,NELE,NF,NIN,NODEL,NODNUM,NOUT,NQP,NSTEPS,NUMLOD,
     *        NUMSS,NWORK,OUTNOD,RESNOD,RESTR,STEER,TOTDOF,TOTELS,
     *        TOTLOD,TOTNOD
      DOUBLE PRECISION A0,A1,ABSS,ALPHA,AMP,AREA,B,BETA,BT,BTDB,C1,C2,
     *                 C3,C4,C5,C6,C7,COORD,D,D2A0,D2A1,DA0,DA1,DB,DET,
     *                 DTIM,E,ELK,ELM,ETA,FORCE,FUN,GDER,GEOM,JAC,JACIN,
     *                 LDER,LOAD,LOADS,NTN,NU,OMEGA,PHASE,QUOT,RHO,SHP,
     *                 SYSK,SYSM,SYSW,T,THETA,TIME,TSHP,WGHT,WORK,XI
      LOGICAL FIRST
      DIMENSION ABSS(3,9),B(6,24),BT(24,6),BTDB(24,24),D(6,6),DB(6,24),
     *          ELK(24,24),ELM(24,24),FUN(8),GDER(3,8),GEOM(8,3),
     *          JAC(3,3),JACIN(3,3),LDER(3,8),LODFRE(20),NTN(24,24),
     *          NWORK(8),SHP(3,24),STEER(8),TSHP(24,3),WGHT(9)
C
C     Problem size dependent arrays
C
      DIMENSION A0(100),A1(100),COORD(100,3),D2A0(100),D2A1(100),
     *          DA0(100),DA1(100),ELTOP(100,10),LOADS(100),NF(100,3),
     *          RESTR(50,4),SYSK(100,25),SYSM(100,25),SYSW(100,25),
     *          WORK(100)
C
      DATA IABSS/3/,IB/6/,IBT/24/,IBTDB/24/,ID/6/,IDB/6/,IELK/24/,
     *     IELM/24/,IFUN/8/,IGDER/3/,IGEOM/8/,IJAC/3/,IJACIN/3/,
     *     ILDER/3/,INTN/24/,ISHP/3/,ISTEER/8/,ITSHP/24/,IWGHT/9/,
     *     JABSS/9/,JB/24/,JBT/6/,JBTDB/24/,JD/6/,JDB/24/,JELK/24/,
     *     JELM/24/,JGDER/8/,JGEOM/3/,JJAC/3/,JJACIN/3/,JLDER/8/,
     *     JNTN/24/,JSHP/24/,JTSHP/3/
C
C     Problem size dependent data statements
C
      DATA IA0/100/,IA1/100/,ICOORD/100/,ID2A0/100/,ID2A1/100/,
     *     IDA0/100/,IDA1/100/,IELTOP/100/,ILOADS/100/,INF/100/,
     *     IRESTR/50/,ISYSK/100/,ISYSM/100/,ISYSW/100/,IWORK/100/,
     *     JCOORD/3/,JELTOP/10/,JNF/3/,JRESTR/4/,JSYSK/25/,JSYSM/25/,
     *     JSYSW/25/
C
      DATA NIN/5/,NOUT/6/
C
C     Set up TIME dependent forcing function as statement function
C
      FORCE(T) = AMP*COS(OMEGA*T+PHASE)
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
      WRITE (NOUT,9960)
      READ (NIN,9990) TOTNOD,DIMEN
      WRITE (NOUT,9950) TOTNOD,DIMEN
      DO 1000 I = 1,TOTNOD
          READ (NIN,9980) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
          WRITE (NOUT,9940) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
 1000 CONTINUE
C
C     Input of element topology
C
      WRITE (NOUT,9930)
      READ (NIN,9990) ELTYP,TOTELS,NODEL
      WRITE (NOUT,9950) ELTYP,TOTELS,NODEL
      DO 1010 I = 1,TOTELS
          READ (NIN,9990) ELNUM, (ELTOP(ELNUM,J+2),J=1,NODEL)
          WRITE (NOUT,9950) ELNUM, (ELTOP(ELNUM,J+2),J=1,NODEL)
          ELTOP(ELNUM,1) = ELTYP
          ELTOP(ELNUM,2) = NODEL
 1010 CONTINUE
C
C     Input of material properties and construction of stress-strain
C     matrix D for plane strain
C
      WRITE (NOUT,9920)
      READ (NIN,9970) NU,E,RHO
      WRITE (NOUT,9910) NU,E,RHO
      CALL DPSN(D,ID,JD,E,NU,NUMSS,ITEST)
C
C     Input of TIME stepping parameters, rayleigh damping
C     coefficients
C
      READ (NIN,9980) NSTEPS,DTIM,THETA
      WRITE (NOUT,9940) NSTEPS,DTIM,THETA
      READ (NIN,9970) ALPHA,BETA
      WRITE (NOUT,9910) ALPHA,BETA
C
C     Input of number of degrees of freedom per node, input of
C     restrained node data and construction of nodal freedom array
C     NF
C
      WRITE (NOUT,9900)
      READ (NIN,9990) DOFNOD
      WRITE (NOUT,9950) DOFNOD
      READ (NIN,9990) RESNOD
      WRITE (NOUT,9950) RESNOD
      K = DOFNOD + 1
      DO 1020 I = 1,RESNOD
          READ (NIN,9990) (RESTR(I,J),J=1,K)
          WRITE (NOUT,9950) (RESTR(I,J),J=1,K)
 1020 CONTINUE
      CALL FORMNF(RESTR,IRESTR,JRESTR,RESNOD,TOTNOD,DOFNOD,NF,INF,JNF,
     *            TOTDOF,ITEST)
C
C     Loading data input and forcing function parameters
C
      READ (NIN,9970) AMP,OMEGA,PHASE
      WRITE (NOUT,9910) AMP,OMEGA,PHASE
      READ (NIN,9990) NUMLOD
      WRITE (NOUT,9950) NUMLOD
      TOTLOD = 0
      DO 1040 I = 1,NUMLOD
          READ (NIN,9990) NODNUM, (NWORK(J),J=1,DOFNOD)
          WRITE (NOUT,9950) NODNUM, (NWORK(J),J=1,DOFNOD)
          DO 1030 J = 1,DOFNOD
              IF (NWORK(J).NE.0) THEN
                  K = NF(NODNUM,J)
                  IF (K.NE.0) THEN
                      TOTLOD = TOTLOD + 1
                      LODFRE(I) = K
                  END IF
              END IF
 1030     CONTINUE
 1040 CONTINUE
C
C     Output control data
C
      READ (NIN,9990) OUTNOD,FREDOM
      WRITE (NOUT,9950) OUTNOD,FREDOM
C
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
C     *                           *******************************
C     *                           *                             *
C     *                           * System Stiffness Matrix And *
C     *                           *     Mass Matrix Assembly    *
C     *                           *                             *
C     *                           *******************************
C
      CALL VECNUL(WORK,IWORK,TOTDOF,ITEST)
      CALL MATNUL(SYSK,ISYSK,JSYSK,TOTDOF,HBAND,ITEST)
      CALL MATNUL(SYSM,ISYSM,JSYSM,TOTDOF,HBAND,ITEST)
      DOFEL = NODEL*DOFNOD
      CALL QQUA4(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST)
      DO 1090 NELE = 1,TOTELS
          AREA = 0.0D0
          CALL ELGEOM(NELE,ELTOP,IELTOP,JELTOP,COORD,ICOORD,JCOORD,GEOM,
     *                IGEOM,JGEOM,DIMEN,ITEST)
C
C     Integration loop for element stiffness and element lumped mass
C     using NQP quadrature points
C
          CALL MATNUL(ELK,IELK,JELK,DOFEL,DOFEL,ITEST)
          CALL MATNUL(ELM,IELM,JELM,DOFEL,DOFEL,ITEST)
          DO 1080 IQUAD = 1,NQP
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
              DO 1070 I = 1,DOFEL
                  DO 1060 J = 1,DOFEL
                      BTDB(I,J) = BTDB(I,J)*QUOT
                      NTN(I,J) = NTN(I,J)*QUOT*RHO
 1060             CONTINUE
 1070         CONTINUE
              CALL MATADD(ELK,IELK,JELK,BTDB,IBTDB,JBTDB,DOFEL,DOFEL,
     *                    ITEST)
              CALL MATADD(ELM,IELM,JELM,NTN,INTN,JNTN,DOFEL,DOFEL,ITEST)
 1080     CONTINUE
C
C     Assembly of system stiffness matrix SYSK, system lumped mass
C     matrix SYSM
C
          CALL DIRECT(NELE,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,STEER,
     *                ISTEER,ITEST)
          CALL ASSYM(SYSK,ISYSK,JSYSK,ELK,IELK,JELK,STEER,ISTEER,HBAND,
     *               DOFEL,ITEST)
          CALL ASSYM(SYSM,ISYSM,JSYSM,ELM,IELM,JELM,STEER,ISTEER,HBAND,
     *               DOFEL,ITEST)
 1090 CONTINUE
C
C     *                           *********************
C     *                           *                   *
C     *                           * Equation Solution *
C     *                           *                   *
C     *                           *********************
C
C     Set initial conditions
C
      CALL VECNUL(A0,IA0,TOTDOF,ITEST)
      CALL VECNUL(DA0,IDA0,TOTDOF,ITEST)
      CALL VECNUL(D2A0,ID2A0,TOTDOF,ITEST)
      C1 = ALPHA + 1.D0/ (THETA*DTIM)
      C2 = BETA + THETA*DTIM
      C3 = (1.D0-THETA)*DTIM
      C4 = BETA - (1.D0-THETA)*DTIM
      C5 = 1.D0/ (THETA*DTIM)
      C6 = (1.D0-THETA)/THETA
      C7 = THETA*DTIM
C
C     Construction of modified system matrices and reduction of left
C     hand side matrix SYSK by choleski reduction.
C
      DO 1110 I = 1,TOTDOF
          DO 1100 J = 1,HBAND
              SYSW(I,J) = C1*SYSM(I,J) + C4*SYSK(I,J)
              SYSK(I,J) = C1*SYSM(I,J) + C2*SYSK(I,J)
              SYSM(I,J) = SYSM(I,J)/THETA
 1100     CONTINUE
 1110 CONTINUE
      CALL CHORDN(SYSK,ISYSK,JSYSK,TOTDOF,HBAND,ITEST)
C
C     Solution for displacements, velocities and accelerations at
C     each TIME step
C
      WRITE (NOUT,9890) OUTNOD,FREDOM
      FREDOM = NF(OUTNOD,FREDOM)
      WRITE (NOUT,9880)
      WRITE (NOUT,9870)
      TIME = 0.0D0
      DO 1140 I = 1,NSTEPS
          TIME = TIME + DTIM
          CALL VECNUL(LOADS,ILOADS,TOTDOF,ITEST)
          LOAD = C7*FORCE(TIME) + C3*FORCE(TIME-DTIM)
          DO 1120 J = 1,TOTLOD
              K = LODFRE(J)
              LOADS(K) = LOAD
 1120     CONTINUE
          CALL MVSYB(SYSM,ISYSM,JSYSM,DA0,IDA0,WORK,IWORK,TOTDOF,HBAND,
     *               ITEST)
          CALL VECADD(LOADS,ILOADS,WORK,IWORK,TOTDOF,ITEST)
          CALL MVSYB(SYSW,ISYSW,JSYSW,A0,IA0,WORK,IWORK,TOTDOF,HBAND,
     *               ITEST)
          CALL VECADD(LOADS,ILOADS,WORK,IWORK,TOTDOF,ITEST)
          CALL CHOSUB(SYSK,ISYSK,JSYSK,LOADS,ILOADS,TOTDOF,HBAND,ITEST)
          CALL VECCOP(LOADS,ILOADS,A1,IA1,TOTDOF,ITEST)
          DO 1130 K = 1,TOTDOF
              DA1(K) = C5* (A1(K)-A0(K)) - C6*DA0(K)
              D2A1(K) = C5* (DA1(K)-DA0(K)) - C6*D2A0(K)
 1130     CONTINUE
          CALL VECCOP(A1,IA1,A0,IA0,TOTDOF,ITEST)
          CALL VECCOP(DA1,IDA1,DA0,IDA0,TOTDOF,ITEST)
          CALL VECCOP(D2A1,ID2A1,D2A0,ID2A0,TOTDOF,ITEST)
          WRITE (NOUT,9860) TIME,A1(FREDOM),DA1(FREDOM),D2A1(FREDOM)
 1140 CONTINUE
C
C     End of Program
C
      STOP
C
 9990 FORMAT (16I5)
 9980 FORMAT (I5,6F10.0)
 9970 FORMAT (6F10.0)
 9960 FORMAT (' ',/,/,' **** NODAL GEOMETRY ****',/,/,' ')
 9950 FORMAT (' ',16I5)
 9940 FORMAT (' ',I5,6F10.5)
 9930 FORMAT (' ',/,/,' **** ELEMENT TOPOLOGY ****',/,/,' ')
 9920 FORMAT (' ',/,/,' **** MATERIAL PROPERTIES ****',/,/,' ')
 9910 FORMAT (' ',10D10.3)
 9900 FORMAT (' ',/,/,' **** RESTRAINT DATA ****',/,/,' ')
 9890 FORMAT ('1RESULTS FOR NODE ',I5,2X,'FREEDOM ',I1,/,/,/,' ')
 9880 FORMAT (' ',/,/,/,45X,'''',14X,'"')
 9870 FORMAT (' ',14X,'T',3 (14X,'V'),/,/,' ')
 9860 FORMAT (' ',10X,4 (D10.3,5X))
C
      END
