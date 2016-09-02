C***********************************************************************
C
C     Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
C                          Chilton, Didcot, Oxfordshire OX11 0QX
C
C     Contact: Prof Chris Greenough
C     Email: c.greenough@rl.ac.uk
C     Tel: (44) 1235 445307
C
C     Finite Library Program: Seg1p1 - Plane Strain of Elastic Solid
C
C***********************************************************************
C
      PROGRAM SEG1P1
C
      INTEGER DIF,DIMEN,DOFEL,DOFNOD,ELNUM,ELTOP,ELTYP,HBAND,I,IABSS,IB,
     *        IBT,IBTDB,ICOORD,ID,IDB,IELDIS,IELK,IELTOP,IFUN,IGDER,
     *        IGEOM,IJAC,IJACIN,ILDER,ILOADS,INF,IQUAD,IRESTR,ISTEER,
     *        ISTRN,ISTRS,ISYSK,ITEST,IWGHT,J,JABSS,JB,JBT,JBTDB,JCOORD,
     *        JD,JDB,JELK,JELTOP,JGDER,JGEOM,JJAC,JJACIN,JLDER,JNF,
     *        JRESTR,JSTEER,JSYSK,K,LODNOD,NELE,NF,NIN,NODEL,NODNUM,
     *        NOUT,NQP,NTEMP,NUMSS,RESNOD,RESTR,STEER,TOTDOF,TOTELS,
     *        TOTNOD
      DOUBLE PRECISION ABSS,B,BT,BTDB,COORD,D,DB,DET,E,ELDIS,ELK,ETA,
     *                 FUN,GDER,GEOM,JAC,JACIN,LDER,LOADS,NU,QUOT,STRN,
     *                 STRS,SYSK,WGHT,WORK,XI
      LOGICAL FIRST
      DIMENSION ABSS(3,9),B(6,24),BT(24,6),BTDB(24,24),D(6,6),DB(6,24),
     *          ELDIS(24),ELK(24,24),FUN(8),GDER(3,8),GEOM(8,3),
     *          JAC(3,3),JACIN(3,3),LDER(3,8),STEER(24),STRN(6),STRS(6),
     *          WGHT(9),WORK(3)
C
C     Problem size dependent arrays
C
      DIMENSION COORD(100,3),ELTOP(100,10),LOADS(100),NF(100,3),
     *          RESTR(100,4),SYSK(100,25)
C
      DATA IABSS/3/,IB/6/,IBT/24/,IBTDB/24/,ID/6/,IDB/6/,IELDIS/24/,
     *     IELK/24/,IFUN/8/,IGDER/3/,IGEOM/8/,IJAC/3/,IJACIN/3/,
     *     ILDER/3/,ISTEER/24/,ISTRN/6/,ISTRS/6/,IWGHT/9/,JABSS/9/,
     *     JB/24/,JBT/6/,JBTDB/24/,JCOORD/3/,JD/6/,JDB/24/,JELK/24/,
     *     JGDER/8/,JGEOM/3/,JJAC/3/,JJACIN/3/,JLDER/8/,JNF/3/,JRESTR/4/
C
C     Problem size dependent data statements
C
      DATA ICOORD/100/,IELTOP/100/,ILOADS/100/,INF/100/,IRESTR/100/,
     *     ISYSK/100/,JELTOP/10/,JSYSK/25/
C
      DATA NIN/5/,NOUT/6/,NTEMP/7/
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
      READ (NIN,FMT=9970) NU,E
      WRITE (NOUT,FMT=9910) NU,E
      CALL DPSN(D,ID,JD,E,NU,NUMSS,ITEST)
C
C     Input of number of degrees of freedom per node, input of
C     restrained node data and construction of nodal freedom array
C     NF
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
C
C     Loading data input
C
      WRITE (NOUT,FMT=9890)
      CALL VECNUL(LOADS,ILOADS,TOTDOF,ITEST)
      READ (NIN,FMT=9990) LODNOD
      WRITE (NOUT,FMT=9950) LODNOD
      DO 1040 I = 1,LODNOD
          READ (NIN,FMT=9980) NODNUM, (WORK(J),J=1,DOFNOD)
          WRITE (NOUT,FMT=9940) NODNUM, (WORK(J),J=1,DOFNOD)
          DO 1030 J = 1,DOFNOD
              K = NF(NODNUM,J)
              IF (K.NE.0) LOADS(K) = WORK(J)
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
C     *                           * Global Stiffness Matrix Assembly *
C     *                           *                                  *
C     *                           ************************************
C
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
C     Formation of strain-displacement matrix B and output to WORK
C     file for later recovery process
C
              CALL B2C2(B,IB,JB,GDER,IGDER,JGDER,NODEL,ITEST)
              WRITE (NTEMP) ((B(I,J),I=1,NUMSS),J=1,DOFEL)
C
C     Formation of element stiffness ELK
C
              CALL MATMUL(D,ID,JD,B,IB,JB,DB,IDB,JDB,NUMSS,NUMSS,DOFEL,
     *                    ITEST)
              CALL MATRAN(B,IB,JB,BT,IBT,JBT,NUMSS,DOFEL,ITEST)
              CALL MATMUL(BT,IBT,JBT,DB,IDB,JDB,BTDB,IBTDB,JBTDB,DOFEL,
     *                    NUMSS,DOFEL,ITEST)
              QUOT = ABS(DET)*WGHT(IQUAD)
              DO 1070 I = 1,DOFEL
                  DO 1060 J = 1,DOFEL
                      BTDB(I,J) = BTDB(I,J)*QUOT
 1060             CONTINUE
 1070         CONTINUE
              CALL MATADD(ELK,IELK,JELK,BTDB,IBTDB,JBTDB,DOFEL,DOFEL,
     *                    ITEST)
 1080     CONTINUE
C
C     Assembly of system stiffness matrix SYSK
C
          CALL DIRECT(NELE,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,STEER,
     *                ISTEER,ITEST)
          CALL ASSYM(SYSK,ISYSK,JSYSK,ELK,IELK,JELK,STEER,ISTEER,HBAND,
     *               DOFEL,ITEST)
 1090 CONTINUE
C
C     *                           *********************
C     *                           *                   *
C     *                           * Equation Solution *
C     *                           *                   *
C     *                           *********************
C
C     Solution of system matrix for the nodal displacements
C
      WRITE (NOUT,FMT=9880)
      CALL PRTVAL(LOADS,ILOADS,NF,INF,JNF,DOFNOD,TOTNOD,NOUT,ITEST)
      CALL CHOSOL(SYSK,ISYSK,JSYSK,LOADS,ILOADS,TOTDOF,HBAND,ITEST)
      WRITE (NOUT,FMT=9870)
      CALL PRTVAL(LOADS,ILOADS,NF,INF,JNF,DOFNOD,TOTNOD,NOUT,ITEST)
C
C     *                           **************************
C     *                           *                        *
C     *                           * Stress-strain Recovery *
C     *                           *                        *
C     *                           **************************
C
C     Recovery of stresses and strains at the quadrature sampling
C     points using the element strain-displacement matrix B from the
C     WORK file
C
      REWIND NTEMP
      DO 1120 NELE = 1,TOTELS
          WRITE (NOUT,FMT=9860) NELE
C
C     Select nodal displacements for element NELE using the steering
C     vector
C
          CALL DIRECT(NELE,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,STEER,
     *                ISTEER,ITEST)
          DO 1110 IQUAD = 1,NQP
              CALL VECNUL(ELDIS,IELDIS,DOFEL,ITEST)
              DO 1100 I = 1,DOFEL
                  JSTEER = STEER(I)
                  IF (JSTEER.NE.0) ELDIS(I) = LOADS(JSTEER)
 1100         CONTINUE
C
C     Recover B matrix for point IQUAD and calculate the stresses
C     and strains
C
              READ (NTEMP) ((B(I,J),I=1,NUMSS),J=1,DOFEL)
              CALL MATVEC(B,IB,JB,ELDIS,IELDIS,NUMSS,DOFEL,STRN,ISTRN,
     *                    ITEST)
              CALL MATVEC(D,ID,JD,STRN,ISTRN,NUMSS,NUMSS,STRS,ISTRS,
     *                    ITEST)
              WRITE (NOUT,FMT=9850) IQUAD
              CALL PRTVEC(STRN,ISTRN,NUMSS,NOUT,ITEST)
              CALL PRTVEC(STRS,ISTRS,NUMSS,NOUT,ITEST)
 1110     CONTINUE
 1120 CONTINUE
C
C     End of Program
C
      STOP
C
 9990 FORMAT (16I5)
 9980 FORMAT (I5,6F10.0)
 9970 FORMAT (2F10.0)
 9960 FORMAT (/,/,' **** NODAL GEOMETRY ****',/,/,' ')
 9950 FORMAT (' ',16I5)
 9940 FORMAT (' ',I5,6F10.5)
 9930 FORMAT (/,/,' **** ELEMENT TOPOLOGY ****',/,/,' ')
 9920 FORMAT (/,/,' **** MATERIAL PROPERTIES ****',/,/,' ')
 9910 FORMAT (' ',F10.5,E10.3)
 9900 FORMAT (/,/,' **** RESTRAINT DATA ****',/,/,' ')
 9890 FORMAT (/,/,' **** LOADING DATA ****',/,/,' ')
 9880 FORMAT (/,/,' **** VECTOR OF LOADS ****',/,/,' ')
 9870 FORMAT (/,/,' **** EQUILIBRIUM DISPLACEMENTS (U,V) ****',/,/,' ')
 9860 FORMAT (/,/,' **** STRAINS AND STRESSES FOR ELEMENT ',I3,' ****',
     *       /,/,' ')
 9850 FORMAT (/,' QUADRATURE POINT ',I3,/,' ')
C
      END
