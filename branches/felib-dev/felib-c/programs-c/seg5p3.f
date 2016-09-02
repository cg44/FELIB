C***********************************************************************
C***********************************************************************
C
C    COPYRIGHT (C) 1987 : SERC, RUTHERFORD APPLETON LABORATORY
C                         CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
      INTEGER BNDNOD, DIF, DIMEN, DOFEL, DOFNOD, ELNUM, ELTOP,
     *     ELTYP, HBAND, I, IABSS, ICOORD, IDTD, IDTDP, IELK,
     *     IELQ, IELTOP, IFUN, IGDER, IGDERT, IGEOM, IJAC,
     *     IJACIN, ILDER, INF, IPHI, IQUAD, IRESTR, IRHS, ISTEER,
     *     ISYSK, ITER, ITEST, IWGHT, J, JABSS, JCOORD, JDTD,
     *     JELK, JELTOP, JGDER, JGDERT, JGEOM, JJAC, JJACIN,
     *     JLDER, JNF, JRESTR, JSYSK, MAXIT, NELE, NF, NIN,
     *     NODE, NODEL, NODNUM, NOUT, NQP, RESTR, STEER, TOTDOF,
     *     TOTELS, TOTNOD
      DOUBLE PRECISION ABSS, AREA, BOLTZ, BVAL, COORD, DEBYL,
     *     DELTA, DET, DMAX, DMIN, DOPE, DOPING, DSTOP,
     *     DSTRT, DTD, DTDP, ELEC, ELK, ELQ, EPSIL, ERROR,
     *     ETA, FUN, GDER, GDERT, GEOM, HOLE, JAC, JACIN, LDER,
     *     NORM, NSUBI, PHI, PHIN, PHIP, QUOT, RHS, SYSK, TEMP,
     *     VMAX, VMIN, VTHERM, WGHT, X, XI, XSTOP, XSTRT
      DIMENSION ABSS(3,7), DTD(6,6), DTDP(6), ELK(6,6),
     *     ELQ(6), FUN(6), GDER(3,6), GDERT(6,3), GEOM(6,3),
     *     JAC(3,3), JACIN(3,3), LDER(3,6), STEER(6), WGHT(7)
C
C                            PROBLEM SIZE DEPENDENT ARRAYS
C
      DIMENSION BVAL(30), COORD(100,3), ELTOP(100,10),
     *     NF(100,1), PHI(100), RESTR(30,2), RHS(100),
     *     SYSK(100,25)
C
      DATA IABSS /3/, IDTD /6/, IDTDP /6/, IELK /6/, IELQ /6/,
     *     IFUN /6/, IGDER /3/, IGDERT /6/, IGEOM /6/,
     *     IJAC /3/, IJACIN /3/, ILDER /3/, ISTEER /6/,
     *     IWGHT /7/, JABSS /7/, JCOORD /3/, JDTD /6/, JELK /6/,
     *     JGDER /6/, JGDERT /3/, JGEOM /3/, JJAC /3/, JJACIN /3/,
     *     JLDER /6/, JNF /1/
C
C                            PROBLEM SIZE DEPENDENT DATA STATEMENTS
C
      DATA ICOORD /100/, IELTOP /100/, INF /100/, IPHI /100/,
     *     IRESTR /30/, IRHS /100/, ISYSK /100/, JELTOP /10/,
     *     JRESTR /2/, JSYSK /25/
C
      DATA DOFNOD /1/
C
      DATA NIN /5/, NOUT /6/
C
C                            PHYSICAL CONSTANTS
C
C                            ASSUMES INPUT LENGTHS IN MICROMETERS
      DATA BOLTZ /1.38D-23/, ECHRG /1.6D-19/, EPSIL /1.10D-16/,
     *     NSUBI /1.48D-2/, TEMP /3.00D2/
C
C                            NORMALISATION CONSTANTS
C
      DEBYL =  DSQRT(EPSIL/ECHRG*BOLTZ/ECHRG*TEMP/NSUBI)
C      DEBYL = DSQRT(EPSIL*BOLTZ*TEMP/(ECHRG**2*NSUBI))
      VTHERM = BOLTZ*TEMP/ECHRG
C
C                            SET ITEST FOR FULL CHECKING
C
      ITEST = 0
C
C*                           **********************
C*                           *                    *
C*                           * INPUT DATA SECTION *
C*                           *                    *
C*                           **********************
C
C                            INPUT OF NODAL GEOMETRY
C
      WRITE (NOUT,9010)
      READ (NIN,8010) TOTNOD, DIMEN
      WRITE (NOUT,9020) TOTNOD, DIMEN
      DO 1020 I=1,TOTNOD
      READ (NIN,8020) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
      WRITE (NOUT,9030) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
C
C                            NORMALISE COORD
      DO 1010 J=1,DIMEN
      COORD(I,J) = COORD(I,J)/DEBYL
 1010 CONTINUE
 1020 CONTINUE
C
C                            INPUT OF ELEMENT TOPOLOGY
C
      WRITE (NOUT,9040)
      READ (NIN,8010) TOTELS
      WRITE (NOUT,9020) TOTELS
      DO 1030 I=1,TOTELS
      READ (NIN,8010) ELNUM, ELTOP(ELNUM,1), NODEL,
     *     (ELTOP(ELNUM,J+2),J=1,NODEL)
      WRITE (NOUT,9020) ELNUM, ELTOP(ELNUM,1), NODEL,
     *     (ELTOP(ELNUM,J+2),J=1,NODEL)
      ELTOP(ELNUM,2) = NODEL
 1030 CONTINUE
C
C                            INPUT OF DOPING PARAMETERS
C
      WRITE (NOUT,9050)
      READ (NIN,8030) XSTRT, XSTOP
      WRITE (NOUT,9060) XSTRT, XSTOP
      READ (NIN,8030) DSTRT, DSTOP
      WRITE (NOUT,9060) DSTRT, DSTOP
C
C
C                            INPUT OF NUMBER OF DEGREES OF FREEDOM
C                            PER NODE, INPUT OF BOUNDARY CON-
C                            DITIONS AND CONSTRUCTION OF NODAL
C                            FREEDOM ARRAY NF
C
      WRITE (NOUT,9070)
      READ (NIN,8010) BNDNOD
      WRITE (NOUT,9020) BNDNOD
      VMAX = 0.D0
      VMIN = 0.D0
      DO 1040 I=1,BNDNOD
      READ (NIN,8040) RESTR(I,1), BVAL(I)
      WRITE (NOUT,9030) RESTR(I,1), BVAL(I)
      RESTR(I,2) = 1
      BVAL(I) = BVAL(I)/VTHERM
      VMAX = DMAX1(VMAX,BVAL(I))
      VMIN = DMIN1(VMIN,BVAL(I))
 1040 CONTINUE
C
C                            MAX NUMBER OF ITERATIONS AND ERROR
C
      WRITE (NOUT,9080)
      READ (NIN,8050) MAXIT, ERROR
      WRITE (NOUT,9030) MAXIT, ERROR
C
C
C                            CREATE THE FREEDOM ARRAY
C
      CALL FORMNF(RESTR, IRESTR, JRESTR, BNDNOD, TOTNOD, DOFNOD,
     *     NF, INF, JNF, TOTDOF, ITEST)
C
C                            CALCULATION OF SEMI-BANDWIDTH
C
      CALL BNDWTH(ELTOP, IELTOP, JELTOP, NF, INF, JNF, DOFNOD,
     *     TOTELS, HBAND, ITEST)
C
C                            CALCULATION OF MINIMUM AND
C                            MAXIMUM DOPING
C
      DO 1050 I=1,TOTNOD
      X = COORD(I,1)
      DOPE = DOPING(X,XSTRT,XSTOP,DSTRT,DSTOP,DEBYL,NSUBI)
      DMIN = DMIN1(DMIN,DOPE)
      DMAX = DMAX1(DMAX,DOPE)
 1050 CONTINUE
C
C                            CALCULATION OF QUASI-FERMI POTENTIALS
C
      PHIP = VMIN + DLOG(-0.5D0*DMIN+DSQRT(1.D0+(0.5D0*DMIN)**
     *     2))
      PHIN = VMAX - DLOG(0.5D0*DMAX+DSQRT(1.D0+(0.5D0*DMAX)**
     *     2))
      DO 1055 I=1,TOTNOD
      X=COORD(I,1)
      DOPE=DOPING(X,XSTRT,XSTOP,DSTRT,DSTOP,DEBYL,NSUBI)
      IF(DOPE.GT.0.D0)PHI(I)=PHIN+DLOG(0.5D0*DOPE+DSQRT(1.D0+
     *(0.5D0*DOPE)**2))
      IF(DOPE.LT.0.D0)PHI(I)=PHIP-DLOG(-0.5D0*DOPE+DSQRT(1.D0+
     *(0.5D0*DOPE)**2))
      IF(DOPE.EQ.0.D0)PHI(I)=(PHIN+PHIP)/2
 1055 CONTINUE
C
C                            ASSIGN THE DIRICHLET BOUNDARY VALUES
C
      DO 1060 I=1,BNDNOD
      J = RESTR(I,1)
      PHI(J) = BVAL(I)
 1060 CONTINUE
C
C                            QUADRATURE WEIGHTS AND ABSSICAE
C
      CALL QTRI4(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)
C
C                            INITIALISE RHS
C
      CALL VECNUL(RHS, IRHS, TOTDOF, ITEST)
C
C*                           ***********************************
C*                           *                                 *
C*                           *         ITERATION LOOP          *
C*                           *                                 *
C*                           ***********************************
      ITER = 0
 1070 ITER = ITER + 1
      WRITE(NOUT,*)'Iteration :',ITER
C
C*                           ***********************************
C*                           *                                 *
C*                           * SYSTEM STIFFNESS MATRIX ASSEMBLY*
C*                           *                                 *
C*                           ***********************************
C
      CALL MATNUL(SYSK, ISYSK, JSYSK, TOTDOF, HBAND, ITEST)
      DO 1120 NELE=1,TOTELS
      NODEL = ELTOP(NELE,2)
      DOFEL = NODEL*DOFNOD
      CALL ELGEOM(NELE, ELTOP, IELTOP, JELTOP, COORD, ICOORD,
     *     JCOORD, GEOM, IGEOM, JGEOM, DIMEN, ITEST)
C
C
      CALL MATNUL(ELK, IELK, JELK, DOFEL, DOFEL, ITEST)
      CALL VECNUL(ELQ, IELQ, DOFEL, ITEST)
C
C                            INTEGRATION LOOP FOR ELEMENT STIFFNESS
C                            USING NQP QUADRATURE POINTS
C
      AREA = 0.D0
      DO 1100 IQUAD=1,NQP
C
C                            FORM LINEAR SHAPE FUNCTION AND SPACE
C                            DERIVATIVES IN THE LOCAL CORRDINATES.
C                            TRANSFORM LOCAL DERIVATIVES TO GLOBAL
C                            COORDINATE SYSTEM
C
      XI = ABSS(1,IQUAD)
      ETA = ABSS(2,IQUAD)
      CALL TRIM3(FUN, IFUN, LDER, ILDER, JLDER, XI, ETA, ITEST)
C
      CALL MATMUL(LDER, ILDER, JLDER, GEOM, IGEOM, JGEOM, JAC,
     *     IJAC, JJAC, DIMEN, NODEL, DIMEN, ITEST)
      CALL MATINV(JAC, IJAC, JJAC, JACIN, IJACIN, JJACIN, DIMEN,
     *     DET, ITEST)
      CALL MATMUL(JACIN, IJACIN, JJACIN, LDER, ILDER, JLDER, GDER,
     *     IGDER, JGDER, DIMEN, DIMEN, NODEL, ITEST)
C
C                            FORMATION OF ELEMENT STIFFNESS ELK
C
      CALL MATRAN(GDER, IGDER, JGDER, GDERT, IGDERT, JGDERT,
     *     DIMEN, DOFEL, ITEST)
      CALL MATMUL(GDERT, IGDERT, JGDERT, GDER, IGDER, JGDER, DTD,
     *     IDTD, JDTD, DOFEL, DIMEN, DOFEL, ITEST)
C
      QUOT = DABS(DET)*WGHT(IQUAD)
      AREA = AREA + QUOT
      DO 1090 I=1,DOFEL
      DTDP(I) = 0.D0
      DO 1080 J=1,DOFEL
      DTD(I,J) = DTD(I,J)*QUOT
      NODE = ELTOP(NELE,J+2)
      DTDP(I) = DTDP(I) - DTD(I,J)*PHI(NODE)
 1080 CONTINUE
 1090 CONTINUE
C
      CALL MATADD(ELK, IELK, JELK, DTD, IDTD, JDTD, DOFEL, DOFEL,
     *     ITEST)
      CALL VECADD(ELQ, IELQ, DTDP, IDTDP, DOFEL, ITEST)
 1100 CONTINUE
C
C                            ADD IN THE CHARGE TERMS (NODEL VALUE *
C                            AREA WEIGHT)
      DO 1110 I=1,NODEL
      if(iter.eq.4 .and. nele.eq.48 .and. i.eq.3) then
      iii=0
      endif
      NODE = ELTOP(NELE,I+2)
      HOLE = DEXP(PHIP-PHI(NODE))
      ELEC = DEXP(PHI(NODE)-PHIN)
      X = COORD(NODE,1)
      DOPE = DOPING(X,XSTRT,XSTOP,DSTRT,DSTOP,DEBYL,NSUBI)
      ELK(I,I) = ELK(I,I) + (HOLE+ELEC)*AREA/3.D0
      ELQ(I) = ELQ(I) + (DOPE+HOLE-ELEC)*AREA/3.D0
 1110 CONTINUE
C
C                            ASSEMBLY OF SYSTEM STIFFNESS MATRIX
C
      CALL DIRECT(NELE, ELTOP, IELTOP, JELTOP, NF, INF, JNF,
     *     DOFNOD, STEER, ISTEER, ITEST)
      CALL ASSYM(SYSK, ISYSK, JSYSK, ELK, IELK, JELK, STEER,
     *     ISTEER, HBAND, DOFEL, ITEST)
      CALL ASRHS(RHS, IRHS, ELQ, IELQ, STEER, ISTEER, DOFEL, ITEST)
 1120 CONTINUE
C
C*                           *********************
C*                           *                   *
C*                           * EQUATION SOLUTION *
C*                           *                   *
C*                           *********************
C
C                            SOLUTION OF SYSTEM MATRIX FOR THE
C                            NODAL VALUES OF THE POTENTIAL
C
      CALL CHOSOL(SYSK, ISYSK, JSYSK, RHS, IRHS, TOTDOF, HBAND,
     *     ITEST)
C
C*                           ***************************
C*                           *                         *
C*                           * UPDATE OF POTENTIAL AND *
C*                           *  CONVERGENCE CHECKING   *
C*                           *                         *
C*                           ***************************
C
      CALL UPDATE(PHI, IPHI, RHS, IRHS, TOTNOD, DOFNOD, TOTDOF,
     *     NF, INF, JNF, ITEST)
      DELTA = NORM(RHS,IRHS,TOTDOF,ITEST)*VTHERM
C
C                            CHECK CONVERGENCE AND NUMBER OF
C                            ITERATIONS
C
      IF (DELTA.GT.ERROR .AND. ITER.LT.MAXIT) GO TO 1070
C
      WRITE (NOUT,9090) ITER
      DO 1130 I=1,TOTNOD
      PHI(I) = PHI(I)*VTHERM
 1130 CONTINUE
C
      WRITE (NOUT,9100)
      CALL PRTVEC(PHI, IPHI, TOTNOD, NOUT, ITEST)
      STOP
C
 8010 FORMAT (16I5)
 8020 FORMAT (I5, 6F10.0)
 8030 FORMAT (2F10.0)
 8040 FORMAT (I5, 6F10.0)
 8050 FORMAT (I5, F10.0)
 9010 FORMAT (//25H **** NODAL GEOMETRY ****//1H )
 9020 FORMAT (1H , 16I5)
 9030 FORMAT (1H , I5, 6F10.5)
 9040 FORMAT (//27H **** ELEMENT TOPOLOGY ****//1H )
 9050 FORMAT (//37H **** DOPING ORIGIN AND GRADIENT ****//1H )
 9060 FORMAT (1H , 2E10.3)
 9070 FORMAT (//30H **** BOUNDARY CONDITIONS ****//1H )
 9080 FORMAT (//39H **** MAXIMUM ITERATIONS AND ERROR ****//1H )
 9090 FORMAT (//16HTOTAL ITERATIONS, I5)
 9100 FORMAT (//27H **** NODAL POTENTIALS ****//1H )
      END
C
      DOUBLE PRECISION FUNCTION DOPING(X, XSTRT, XSTOP, DSTRT,
     *     DSTOP, DEBYL, NSUBI)
      DOUBLE PRECISION DEBYL, DSTOP, DSTRT, NSUBI, X, XSTOP,
     *     XSTRT
      X = X*DEBYL
      IF (X.LE.XSTRT) DOPING = DSTRT
      IF (X.GE.XSTOP) DOPING = DSTOP
      IF (X.GT.XSTRT .AND. X.LT.XSTOP) DOPING = DSTRT +
     *     (DSTOP-DSTRT)*(X-XSTRT)/(XSTOP-XSTRT)
      DOPING = DOPING/NSUBI
      X = X/DEBYL
      RETURN
      END
C
