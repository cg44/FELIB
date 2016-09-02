C***********************************************************************
C
C     Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
C                          Chilton, Didcot, Oxfordshire OX11 0QX
C
C     Contact: Prof Chris Greenough
C     Email: c.greenough@rl.ac.uk
C     Tel: (44) 1235 445307
C
C     Finite Library Program: Seg2p1 - Free Vibration of an Elastic 
C                                      Solid
C
C***********************************************************************
C
      PROGRAM SEG2P1
C
      INTEGER DIMEN,DOFEL,DOFNOD,ELNUM,ELTOP,ELTYP,I,IABSS,IB,IBT,IBTDB,
     *        ICOORD,ID,IDB,IDIAG,IELK,IELM,IELTOP,IFTF,IFUN,IGDER,
     *        IGEOM,IJAC,IJACIN,ILDER,INF,IQUAD,IRESTR,ISHP,ISTEER,ISUB,
     *        ISYSK,ISYSM,ITEST,ITSHP,IWGHT,J,JABSS,JB,JBT,JBTDB,JCOORD,
     *        JD,JDB,JELK,JELM,JELTOP,JFTF,JGDER,JGEOM,JJAC,JJACIN,
     *        JLDER,JNF,JRESTR,JSHP,JSYSK,JTSHP,K,NELE,NF,NIN,NMODES,
     *        NODEL,NODNUM,NOUT,NQP,NUMSS,RESNOD,RESTR,STEER,TOTDOF,
     *        TOTELS,TOTNOD
      DOUBLE PRECISION ABSS,AREA,B,BT,BTDB,COORD,D,DB,DET,DIAG,E,ELK,
     *                 ELM,EPS,ETA,FTF,FUN,GDER,GEOM,JAC,JACIN,LDER,NU,
     *                 QUOT,RHO,SHP,SUB,SYSK,SYSM,TOL,TSHP,VEPS,VTOL,
     *                 WGHT,X,XI
      DIMENSION ABSS(3,9),B(6,24),BT(24,6),BTDB(24,24),D(6,6),DB(6,24),
     *          ELK(24,24),ELM(24,24),FTF(24,24),FUN(8),GDER(3,8),
     *          GEOM(8,3),JAC(3,3),JACIN(3,3),LDER(3,8),SHP(3,24),
     *          STEER(24),TSHP(24,3),WGHT(9)
C
C     Problem size dependent arrays
C
      DIMENSION COORD(100,3),DIAG(100),ELTOP(100,10),NF(100,3),
     *          RESTR(50,4),SUB(100),SYSK(100,100),SYSM(100)
C
      DATA IABSS/3/,IB/6/,IBT/24/,IBTDB/24/,ID/6/,IDB/6/,IELK/24/,
     *     IELM/24/,IFTF/24/,IFUN/8/,IGDER/3/,IGEOM/8/,IJAC/3/,
     *     IJACIN/3/,ILDER/3/,ISHP/3/,ISTEER/24/,ITSHP/24/,IWGHT/9/,
     *     JABSS/9/,JB/24/,JBT/6/,JBTDB/8/,JCOORD/3/,JD/6/,JDB/24/,
     *     JELK/24/,JELM/24/,JFTF/24/,JGDER/8/,JGEOM/3/,JJAC/3/,
     *     JJACIN/3/,JLDER/8/,JNF/3/,JRESTR/4/,JSHP/24/,JTSHP/3/
C
C     Problem size dependent data statements
C
      DATA ICOORD/100/,IDIAG/100/,IELTOP/100/,INF/100/,IRESTR/50/,
     *     ISUB/100/,ISYSK/100/,ISYSM/100/,JELTOP/10/,JSYSK/100/
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
      WRITE (NOUT,FMT=9960)
      READ (NIN,FMT=9990) TOTNOD,DIMEN
      WRITE (NOUT,FMT=9950) TOTNOD,DIMEN
      DO 1000 I = 1,TOTNOD
          READ (NIN,FMT=9980) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
          WRITE (NOUT,FMT=9940) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
 1000 CONTINUE
C
C     Input of element topology
C
      WRITE (NOUT,FMT=9930)
      READ (NIN,FMT=9990) ELTYP,TOTELS,NODEL
      WRITE (NOUT,FMT=9950) ELTYP,TOTELS,NODEL
      DO 1010 I = 1,TOTELS
          READ (NIN,FMT=9990) ELNUM, (ELTOP(ELNUM,J+2),J=1,NODEL)
          WRITE (NOUT,FMT=9950) ELNUM, (ELTOP(ELNUM,J+2),J=1,NODEL)
          ELTOP(ELNUM,1) = ELTYP
          ELTOP(ELNUM,2) = NODEL
 1010 CONTINUE
C
C     Input of material properties and construction of stress-strain
C     matrix D for plane strain
C
      WRITE (NOUT,FMT=9920)
      READ (NIN,FMT=9970) NU,E,RHO
      WRITE (NOUT,FMT=9910) NU,E,RHO
      CALL DPSN(D,ID,JD,E,NU,NUMSS,ITEST)
C
C     Input of number of degrees of freedom per node, input of
C     restrained node data and construction of nodal freedom array NF
C
      WRITE (NOUT,FMT=9900)
      READ (NIN,FMT=9990) DOFNOD
      WRITE (NOUT,FMT=9950) DOFNOD
      READ (NIN,FMT=9990) RESNOD
      WRITE (NOUT,FMT=9950) RESNOD
      K = DOFNOD + 1
      DO 1020 I = 1,RESNOD
          READ (NIN,FMT=9990) (RESTR(I,J),J=1,K)
          WRITE (NOUT,FMT=9950) (RESTR(I,J),J=1,K)
 1020 CONTINUE
      CALL FORMNF(RESTR,IRESTR,JRESTR,RESNOD,TOTNOD,DOFNOD,NF,INF,JNF,
     *            TOTDOF,ITEST)
      READ (NIN,FMT=9990) NMODES
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
      DO 1060 NELE = 1,TOTELS
          AREA = 0.0D0
          CALL ELGEOM(NELE,ELTOP,IELTOP,JELTOP,COORD,ICOORD,JCOORD,GEOM,
     *                IGEOM,JGEOM,DIMEN,ITEST)
C
C     Integration loop for element stiffness and element lumped mass
C     using NQP quadrature points
C
          CALL MATNUL(ELK,IELK,JELK,DOFEL,DOFEL,ITEST)
          CALL MATNUL(ELM,IELM,JELM,DOFEL,DOFEL,ITEST)
          DO 1050 IQUAD = 1,NQP
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
              CALL MATMUL(TSHP,ITSHP,JTSHP,SHP,ISHP,JSHP,FTF,IFTF,JFTF,
     *                    DOFEL,DOFNOD,DOFEL,ITEST)
C
              QUOT = ABS(DET)*WGHT(IQUAD)
              AREA = AREA + QUOT
              DO 1040 I = 1,DOFEL
                  DO 1030 J = 1,DOFEL
                      BTDB(I,J) = BTDB(I,J)*QUOT
                      FTF(I,J) = FTF(I,J)*QUOT*RHO
 1030             CONTINUE
 1040         CONTINUE
              CALL MATADD(ELK,IELK,JELK,BTDB,IBTDB,JBTDB,DOFEL,DOFEL,
     *                    ITEST)
              CALL MATADD(ELM,IELM,JELM,FTF,IFTF,JFTF,DOFEL,DOFEL,ITEST)
 1050     CONTINUE
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
 1060 CONTINUE
C
C     *                           *********************
C     *                           *                   *
C     *                           * Equation Solution *
C     *                           *                   *
C     *                           *********************
C
C
C     Reduction of system mass matrix by choleski reduction and
C     construction of modified system stiffness matrix
C
      DO 1070 I = 1,TOTDOF
          SYSM(I) = 1.0D0/SQRT(SYSM(I))
 1070 CONTINUE
      DO 1090 I = 1,TOTDOF
          DO 1080 J = 1,TOTDOF
              SYSK(I,J) = SYSM(I)*SYSK(I,J)*SYSM(J)
 1080     CONTINUE
 1090 CONTINUE
C
C     Householder tridiagonalisation of SYSK and eigenvalue and
C     eigenvector recovery using ql transformations
C
      TOL = VTOL(X)
      CALL HOUSE(SYSK,ISYSK,JSYSK,SYSK,ISYSK,JSYSK,DIAG,IDIAG,SUB,ISUB,
     *           TOTDOF,TOL,ITEST)
      EPS = VEPS(X)
      CALL QLVEC(DIAG,IDIAG,SUB,ISUB,SYSK,ISYSK,JSYSK,TOTDOF,EPS,ITEST)
      WRITE (NOUT,FMT=9890)
      CALL PRTVEC(DIAG,IDIAG,TOTDOF,NOUT,ITEST)
      DO 1110 J = 1,NMODES
          DO 1100 I = 1,TOTDOF
              DIAG(I) = SYSK(I,J)*SYSM(I)
 1100     CONTINUE
          WRITE (NOUT,FMT=9880) J
          CALL PRTVEC(DIAG,IDIAG,TOTDOF,NOUT,ITEST)
 1110 CONTINUE
C
C     End of Program
C
      STOP
C
 9990 FORMAT (16I5)
 9980 FORMAT (I5,6F10.0)
 9970 FORMAT (6F10.0)
 9960 FORMAT (/,/,' **** NODAL GEOMETRY ****',/,/,' ')
 9950 FORMAT (' ',16I5)
 9940 FORMAT (' ',I5,6F10.5)
 9930 FORMAT (/,/,' **** ELEMENT TOPOLOGY ****',/,/,' ')
 9920 FORMAT (/,/,' **** MATERIAL PROPERTIES ****',/,/,' ')
 9910 FORMAT (' ',F10.5,E10.3,F10.5)
 9900 FORMAT (/,/,' **** RESTRAINT DATA ****',/,/,' ')
 9890 FORMAT (/,/,' **** FREQUENCIES ****',/,' ')
 9880 FORMAT (/,' MODE SHAPE FOR MODE ',I5,/,' ')
C
      END
