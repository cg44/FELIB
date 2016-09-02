C***********************************************************************
C
C     Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
C                          Chilton, Didcot, Oxfordshire OX11 0QX
C
C     Contact: Prof Chris Greenough
C     Email: c.greenough@rl.ac.uk
C     Tel: (44) 1235 445307
C
C     Finite Library Program: Seg4p2 - Sediment Transport in a Column
C                                      of Fluid (Convection/Duffusion)
C
C***********************************************************************
C
      PROGRAM SEG4P2
C
      INTEGER BDCND,BLIST,BTYPE,DIMEN,DOFEL,DOFNOD,ELNUM,ELTOP,ELTYP,
     *        HBAND,I,IABSS,IBDCND,IBELM,IBELV,IBLIST,IBN,IBNTN,ICOORD,
     *        ICOSIN,IDTPD,IELC,IELD,IELK,IELM,IELTOP,IFTF,IFUN,IGDER,
     *        IGDERT,IGEOM,IJAC,IJACIN,ILABSS,ILDER,ILOWER,INF,INTNOD,
     *        INTVD,IP,IPD,IQUAD,IRHS,IROPIV,ISTEER,ISYSK,ISYSM,ITEST,
     *        ITYPE,IVEL,IVGDER,IWGHT,IWORK1,IWORK2,J,JABSS,JBDCND,
     *        JBELM,JBLIST,JBNTN,JCOORD,JDTPD,JELC,JELD,JELK,JELM,
     *        JELTOP,JFTF,JGDER,JGDERT,JGEOM,JJAC,JJACIN,JLDER,JLOWER,
     *        JNF,JNTVD,JP,JPD,JSYSK,JSYSM,K,L,M,NELE,NF,NIN,NODEL,
     *        NODNUM,NODSID,NOUT,NQP,NSTEPS,NUMNOD,NUMSID,ROPIV,SIDNUM,
     *        STEER,TOTBND,TOTDOF,TOTELS,TOTNOD,WIDTH
      DOUBLE PRECISION ABSS,BELM,BELV,BN,BNTN,COEFF,CONST1,CONST2,COORD,
     *                 COSIN,DET,DTIM,DTPD,ELC,ELD,ELK,ELM,ETA,FTF,FUN,
     *                 GDER,GDERT,GEOM,JAC,JACIN,LABSS,LDER,LOWER,NTVD,
     *                 P,PD,QUOT,RHS,SCALE,SYSK,SYSM,THETA,ULEN,VALUE,
     *                 VEL,VGDER,WGHT,WORK1,WORK2,XI
      DIMENSION ABSS(3,9),BELM(8,8),BELV(8),BN(8),BNTN(8,8),COSIN(3),
     *          DTPD(8,8),ELC(8,8),ELD(8,8),ELK(8,8),ELM(8,8),FTF(8,8),
     *          FUN(8),GDER(3,8),GDERT(8,3),GEOM(8,3),JAC(3,3),
     *          JACIN(3,3),LABSS(3),LDER(3,8),NTVD(8,8),P(3,3),PD(3,8),
     *          STEER(8),VEL(8),VGDER(8),WGHT(9)
C
C     Problem size dependent arrays
C
      DIMENSION BDCND(5,40),BLIST(20,2),COORD(100,3),ELTOP(100,10),
     *          LOWER(100,25),NF(100,1),RHS(100),ROPIV(100),
     *          SYSK(100,25),SYSM(100,25),WORK1(100),WORK2(100)
C
      DATA IABSS/3/,IBELM/8/,IBELV/8/,IBN/8/,IBNTN/8/,ICOSIN/3/,
     *     IDTPD/8/,IELC/8/,IELD/8/,IELK/8/,IELM/8/,IFTF/8/,IFUN/8/,
     *     IGDER/3/,IGDERT/8/,IGEOM/8/,IJAC/3/,IJACIN/3/,ILABSS/3/,
     *     ILDER/3/,INTVD/8/,IP/3/,IPD/3/,ISTEER/8/,IVEL/8/,IVGDER/8/,
     *     IWGHT/9/,JABSS/9/,JBELM/8/,JBNTN/8/,JCOORD/3/,JDTPD/8/,
     *     JELC/8/,JELD/8/,JELK/8/,JELM/8/,JFTF/8/,JGDER/8/,JGDERT/3/,
     *     JGEOM/3/,JJAC/3/,JJACIN/3/,JLDER/8/,JNF/1/,JNTVD/8/,JP/3/,
     *     JPD/8/,SCALE/1.0D+10/
C
C     Problem size dependent data statements
C
      DATA IBDCND/5/,IBLIST/20/,ICOORD/100/,IELTOP/100/,ILOWER/100/,
     *     INF/100/,IRHS/100/,IROPIV/100/,ISYSK/100/,ISYSM/100/,
     *     IWORK1/100/,IWORK2/100/,JBDCND/40/,JBLIST/2/,JELTOP/10/,
     *     JLOWER/25/,JSYSK/25/,JSYSM/25/
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
C
      DO 1000 I = 1,TOTNOD
          READ (NIN,FMT=9980) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
          WRITE (NOUT,FMT=9940) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
 1000 CONTINUE
C
C     Input of element topology
C
      WRITE (NOUT,FMT=9930)
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
C     Input of permeabilities and construction of permeability
C     matrix P
C
      WRITE (NOUT,FMT=9920)
      CALL MATNUL(P,IP,JP,DIMEN,DIMEN,ITEST)
      READ (NIN,FMT=9970) (P(I,I),I=1,DIMEN)
      WRITE (NOUT,FMT=9910) (P(I,I),I=1,DIMEN)
C
C     Input of convection velocities
C
      WRITE (NOUT,FMT=9900)
      READ (NIN,FMT=9970) (VEL(I),I=1,DIMEN)
      WRITE (NOUT,FMT=9910) (VEL(I),I=1,DIMEN)
C
C     Input of boundary conditions (BTYPE = 1 (neumann) or 2
C     (cauchy))
C
      WRITE (NOUT,FMT=9890)
      READ (NIN,FMT=9970) CONST1,CONST2
      WRITE (NOUT,FMT=9910) CONST1,CONST2
C
      READ (NIN,FMT=9990) TOTBND
      WRITE (NOUT,FMT=9950) TOTBND
C
      IF (TOTBND.NE.0) THEN
C
          DO 1020 I = 1,TOTBND
              READ (NIN,FMT=9990) BTYPE,NUMNOD,NODSID,
     *          (BDCND(I,J+3),J=1,NUMNOD)
              WRITE (NOUT,FMT=9950) BTYPE,NUMNOD,NODSID,
     *          (BDCND(I,J+3),J=1,NUMNOD)
              BDCND(I,1) = BTYPE
              BDCND(I,2) = NUMNOD
              BDCND(I,3) = NODSID
 1020     CONTINUE
      END IF
C
C     Input number of degrees of freedom per node and form nodal
C     freedom array
C
      READ (NIN,FMT=9990) DOFNOD
      WRITE (NOUT,FMT=9950) DOFNOD
      TOTDOF = 0
      DO 1040 I = 1,TOTNOD
          DO 1030 J = 1,DOFNOD
              TOTDOF = TOTDOF + 1
              NF(I,J) = TOTDOF
 1030     CONTINUE
 1040 CONTINUE
C
C     Input of initial conditions
C
      WRITE (NOUT,FMT=9880)
      CALL VECNUL(RHS,IRHS,TOTDOF,ITEST)
      READ (NIN,FMT=9990) INTNOD
      WRITE (NOUT,FMT=9950) INTNOD
C
      IF (INTNOD.NE.0) THEN
          DO 1050 I = 1,INTNOD
              READ (NIN,FMT=9980) NODNUM,VALUE
              WRITE (NOUT,FMT=9870) NODNUM,VALUE
              J = NF(NODNUM,1)
              RHS(J) = VALUE
 1050     CONTINUE
      END IF
C
C     Input of time stepping data
C
      WRITE (NOUT,FMT=9860)
      READ (NIN,FMT=9980) NSTEPS,DTIM,THETA
      WRITE (NOUT,FMT=9870) NSTEPS,DTIM,THETA
C
C     Calculation of semi-bandwidth
C
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
C     Initialisation of arrays
C
      CALL MATNUL(SYSK,ISYSK,JSYSK,TOTDOF,2*HBAND-1,ITEST)
      CALL MATNUL(SYSM,ISYSM,JSYSM,TOTDOF,2*HBAND-1,ITEST)
      CALL QQUA4(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST)
C
      DO 1090 NELE = 1,TOTELS
          NODEL = ELTOP(NELE,2)
          DOFEL = NODEL*DOFNOD
          CALL ELGEOM(NELE,ELTOP,IELTOP,JELTOP,COORD,ICOORD,JCOORD,GEOM,
     *                IGEOM,JGEOM,DIMEN,ITEST)
C
C     Integration loop for element matrices using NQP quadrature
C     points
C
          CALL MATNUL(ELM,IELM,JELM,DOFEL,DOFEL,ITEST)
          CALL MATNUL(ELD,IELD,JELD,DOFEL,DOFEL,ITEST)
          CALL MATNUL(ELC,IELC,JELC,DOFEL,DOFEL,ITEST)
          DO 1080 IQUAD = 1,NQP
C
C     Form shape function and space derivatives in the local
C     coordinates. Transform local derivatives to global coordinate
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
C     Form integrand of element diffusion matrix ELD
C
              CALL MATMUL(P,IP,JP,GDER,IGDER,JGDER,PD,IPD,JPD,DIMEN,
     *                    DIMEN,DOFEL,ITEST)
              CALL MATRAN(GDER,IGDER,JGDER,GDERT,IGDERT,JGDERT,DIMEN,
     *                    DOFEL,ITEST)
              CALL MATMUL(GDERT,IGDERT,JGDERT,PD,IPD,JPD,DTPD,IDTPD,
     *                    JDTPD,DOFEL,DIMEN,DOFEL,ITEST)
C
C     Form integrand of element mass matrix ELM
C
              CALL DYAD(FUN,IFUN,FUN,IFUN,FTF,IFTF,JFTF,NODEL,ITEST)
C
C     Form integrand of element convection matrix ELC
C
              CALL VECMAT(VEL,IVEL,GDER,IGDER,JGDER,DIMEN,NODEL,VGDER,
     *                    IVGDER,ITEST)
              CALL DYAD(FUN,IFUN,VGDER,IVGDER,NTVD,INTVD,JNTVD,NODEL,
     *                  ITEST)
              QUOT = ABS(DET)*WGHT(IQUAD)
              DO 1070 I = 1,DOFEL
                  DO 1060 J = 1,DOFEL
                      DTPD(I,J) = DTPD(I,J)*QUOT
                      FTF(I,J) = FTF(I,J)*QUOT
                      NTVD(I,J) = NTVD(I,J)*QUOT
 1060             CONTINUE
 1070         CONTINUE
              CALL MATADD(ELM,IELM,JELM,FTF,IFTF,JFTF,DOFEL,DOFEL,ITEST)
              CALL MATADD(ELD,IELD,JELD,DTPD,IDTPD,JDTPD,DOFEL,DOFEL,
     *                    ITEST)
              CALL MATADD(ELC,IELC,JELC,NTVD,INTVD,JNTVD,DOFEL,DOFEL,
     *                    ITEST)
 1080     CONTINUE
C
C     Assembly of system matrices
C
          CALL DIRECT(NELE,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,STEER,
     *                ISTEER,ITEST)
          CALL ASUSM(SYSM,ISYSM,JSYSM,ELM,IELM,JELM,STEER,ISTEER,HBAND,
     *               DOFEL,ITEST)
          CALL ASUSM(SYSK,ISYSK,JSYSK,ELD,IELD,JELD,STEER,ISTEER,HBAND,
     *               DOFEL,ITEST)
          CALL ASUSM(SYSK,ISYSK,JSYSK,ELC,IELC,JELC,STEER,ISTEER,HBAND,
     *               DOFEL,ITEST)
 1090 CONTINUE
C
C     *                         **************************************
C     *                         *                                    *
C     *                         * Application Of Boundary Conditions *
C     *                         *                                    *
C     *                         **************************************
C     *
C     *
C
      CALL QLIN2(WGHT,IWGHT,LABSS,ILABSS,NQP,ITEST)
      CALL VECNUL(WORK1,IWORK1,TOTDOF,ITEST)
C
C     Loop around boundaries
C
      DO 1190 ITYPE = 1,TOTBND
          BTYPE = BDCND(ITYPE,1) - 1
          NUMNOD = BDCND(ITYPE,2)
          CALL SIDENO(TOTELS,ELTOP,IELTOP,JELTOP,ITYPE,BDCND,IBDCND,
     *                JBDCND,NUMSID,BLIST,IBLIST,JBLIST,ITEST)
C
          DO 1180 M = 1,NUMSID
              ELNUM = BLIST(M,1)
              SIDNUM = BLIST(M,2)
              CALL VECNUL(BELV,IBELV,DOFEL,ITEST)
              CALL MATNUL(BELM,IBELM,JBELM,DOFEL,DOFEL,ITEST)
C
C     Construct quadrature rule and local geometry
C
              CALL ELGEOM(ELNUM,ELTOP,IELTOP,JELTOP,COORD,ICOORD,JCOORD,
     *                    GEOM,IGEOM,JGEOM,DIMEN,ITEST)
              CALL BQQUA(ABSS,IABSS,JABSS,LABSS,ILABSS,NQP,SIDNUM,COEFF,
     *                   ITEST)
              DO 1150 J = 1,NQP
                  XI = ABSS(1,J)
                  ETA = ABSS(2,J)
                  CALL QUAM4(FUN,IFUN,LDER,ILDER,JLDER,XI,ETA,ITEST)
                  CALL LINQUA(XI,ETA,GEOM,IGEOM,JGEOM,NODEL,SIDNUM,ULEN,
     *                        ITEST)
                  QUOT = ULEN*WGHT(J)*COEFF
                  GO TO (1100,1120) BTYPE
C
C     Neumann conditions (bottom of column)
C
 1100             CONTINUE
                  DO 1110 K = 1,DOFEL
                      BN(K) = CONST1*VEL(2)*FUN(K)*QUOT
 1110             CONTINUE
                  CALL VECADD(BELV,IBELV,BN,IBN,DOFEL,ITEST)
                  GO TO 1150
C
C     Cauchy conditions (top of column)
C
 1120             CONTINUE
                  DO 1140 K = 1,DOFEL
                      DO 1130 L = 1,DOFEL
                          BNTN(K,L) = CONST2*VEL(2)*FUN(K)*FUN(L)*QUOT
 1130                 CONTINUE
 1140             CONTINUE
                  CALL MATADD(BELM,IBELM,JBELM,BNTN,IBNTN,JBNTN,DOFEL,
     *                        DOFEL,ITEST)
 1150         CONTINUE
C
C     Assembly of boundary conditions
C
              CALL DIRECT(ELNUM,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,
     *                    STEER,ISTEER,ITEST)
              GO TO (1160,1170) BTYPE
 1160         CONTINUE
              CALL ASRHS(WORK1,IWORK1,BELV,IBELV,STEER,ISTEER,NODEL,
     *                   ITEST)
              GO TO 1180
 1170         CONTINUE
              CALL ASUSM(SYSK,ISYSK,JSYSK,BELM,IBELM,JBELM,STEER,ISTEER,
     *                   HBAND,DOFEL,ITEST)
 1180     CONTINUE
 1190 CONTINUE
C
C     *                           *********************
C     *                           *                   *
C     *                           * Equation Solution *
C     *                           *                   *
C     *                           *********************
C
C     Stepping scheme
C
      WIDTH = 2*HBAND - 1
      DO 1210 I = 1,TOTDOF
          DO 1200 J = 1,WIDTH
              WORK2(1) = THETA*SYSK(I,J) + SYSM(I,J)/DTIM
              WORK2(2) = (1.0D0-THETA)*SYSK(I,J) - SYSM(I,J)/DTIM
              SYSK(I,J) = WORK2(1)
              SYSM(I,J) = -WORK2(2)
 1200     CONTINUE
 1210 CONTINUE
C
C     Reduce asymmetric system matrix SYSK
C
      CALL GAURDN(SYSK,ISYSK,JSYSK,LOWER,ILOWER,JLOWER,TOTDOF,HBAND,
     *            ROPIV,IROPIV,ITEST)
C
C     Preform time steps
C
      DO 1220 I = 1,NSTEPS
          WRITE (NOUT,FMT=9850) I
          CALL MVUSB(SYSM,ISYSM,JSYSM,RHS,IRHS,WORK2,IWORK2,TOTDOF,
     *               HBAND,ITEST)
          CALL VECADD(WORK2,IWORK2,WORK1,IWORK1,TOTDOF,ITEST)
          CALL GAUSUB(SYSK,ISYSK,JSYSK,LOWER,ILOWER,JLOWER,TOTDOF,HBAND,
     *                ROPIV,IROPIV,WORK2,IWORK2,ITEST)
          CALL VECCOP(WORK2,IWORK2,RHS,IRHS,TOTDOF,ITEST)
          CALL PRTVAL(RHS,IRHS,NF,INF,JNF,DOFNOD,TOTNOD,NOUT,ITEST)
 1220 CONTINUE
C
C     End of Program
C
      STOP
C
 9990 FORMAT (16I5)
 9980 FORMAT (I5,3F10.0)
 9970 FORMAT (3F10.0)
 9960 FORMAT (/,/,' **** NODAL GEOMETRY ****',/,' ')
 9950 FORMAT (' ',16I5)
 9940 FORMAT (' ',I5,3F10.5)
 9930 FORMAT (/,/,' **** ELEMENT TOPOLOGY ****',/,' ')
 9920 FORMAT (/,/,' **** PERMEABILITIES ****',/,' ')
 9910 FORMAT (' ',3D10.3)
 9900 FORMAT (/,/,' **** CONVECTION VELOCITIES',/,' ')
 9890 FORMAT (/,/,' **** BOUNDARY CONDITIONS ****',/,' ')
 9880 FORMAT (/,/,' **** INITIAL CONDITIONS ****',/,' ')
 9870 FORMAT (' ',I5,3F10.5)
 9860 FORMAT (/,/,' **** TIME STEPPING DATA ****',/,' ')
 9850 FORMAT (/,/,' NODAL CONCENTRATIONS FOR TIME STEP ',I5)
C
      END
