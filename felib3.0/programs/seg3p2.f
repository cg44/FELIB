C $Id: seg3p2.f,v 1.3 2009/02/26 15:05:33 cg44 Exp $
      PROGRAM SEG3P2
C***********************************************************************
C
C    COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                         CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C***********************************************************************
C
      INTEGER BDCND, BLIST, BTYPE, DIF, DIMEN, DOFEL, DOFNOD,
     *     ELNUM, ELTOP, ELTYP, HBAND, I, IABSS, IBDCND,
     *     IBELM, IBELV, IBLIST, IBN, IBNTN, ICOORD, ICOSIN,
     *     IDTPD, IELK, IELTOP, IFUN, IGDER, IGDERT, IGEOM,
     *     IJAC, IJACIN, ILABSS, ILDER, INF, IP, IPD, IQUAD,
     *     IRHS, ISTEER, ISYSK, ITEST, ITYPE, IWGHT, J, JABSS,
     *     JBDCND, JBELM, JBLIST, JBNTN, JCOORD, JDTPD, JELK,
     *     JELTOP, JGDER, JGDERT, JGEOM, JJAC, JJACIN, JLDER,
     *     JNF, JP, JPD, JSYSK, K, L, M, NELE, NF, NIN, NODEL,
     *     NODNUM, NODSID, NOUT, NQP, NUMNOD, NUMSID, SIDNUM,
     *     STEER, TOTBND, TOTDOF, TOTELS, TOTNOD
      DOUBLE PRECISION ABSS, BELM, BELV, BN, BNTN, COEFF, COORD,
     *     COSIN, DET, DTPD, ELK, ETA, F1, F2, FUN, G1, G2,
     *     GDER, GDERT, GEOM, H, JAC, JACIN, LABSS, LDER,
     *     P, PD, PX, PY, QUOT, RHS, SCALE, SYSK, ULEN,
     *     WGHT, X, XI, XX, Y, YY
      LOGICAL FIRST
      DIMENSION ABSS(3,9), BELM(8,8), BELV(8), BN(8),
     *     BNTN(8,8), COSIN(3), DTPD(8,8), ELK(8,8), FUN(8),
     *     GDER(3,8), GDERT(8,3), GEOM(8,3), JAC(3,3),
     *     JACIN(3,3), LABSS(3), LDER(3,8), P(3,3), PD(3,8),
     *     STEER(8), WGHT(9)
C
C                            PROBLEM SIZE DEPENDENT ARRAYS
C
      DIMENSION BDCND(5,40), BLIST(50,5), COORD(100,3),
     *     ELTOP(100,10), NF(100,1), RHS(100), SYSK(100,25)
C
      DATA IABSS /3/, IBELM /8/, IBELV /8/, IBN /8/, IBNTN /8/,
     *     ICOSIN /3/, IDTPD /8/, IELK /8/, IFUN /8/, IGDER /3/,
     *     IGDERT /8/, IGEOM /8/, IJAC /3/, IJACIN /3/,
     *     ILABSS /3/, ILDER /3/, IP /3/, IPD /3/, ISTEER /8/,
     *     IWGHT /9/, JABSS /9/, JBELM /8/, JBNTN /8/, JCOORD /3/,
     *     JDTPD /8/, JELK /8/, JGDER /8/, JGDERT /3/, JGEOM /3/,
     *     JJAC /3/, JJACIN /3/, JLDER /8/, JNF /1/, JP /3/,
     *     JPD /8/, SCALE /1.0D+10/
C
C                            PROBLEM SIZE DEPENDENT DATA STATEMENTS
C
      DATA IBDCND /5/, IBLIST /50/, ICOORD /100/, IELTOP /100/,
     *     INF /100/, IRHS /100/, ISYSK /100/, JBDCND /40/,
     *     JBLIST /5/, JELTOP /10/, JSYSK /25/
C
      DATA NIN /5/, NOUT /6/
C
C                            STATEMENT FUNCTIONS FOR
C                            BOUNDARY CONDITIONS
C
      H(XX,YY) = 0.D0
      F1(XX,YY) = -1.D0
      F2(XX,YY) = 0.D0
      G1(XX,YY) = 0.D0
      G2(XX,YY) = 0.D0
C
C                            STATEMENT FUNCTIONS FOR
C                            PX AND PY
C
      PX(XX,YY) = 1.D0
      PY(XX,YY) = 1.D0
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
      DO 1010 I=1,TOTNOD
      READ (NIN,8020) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
      WRITE (NOUT,9030) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
 1010 CONTINUE
C
C                            INPUT OF ELEMENT TOPOLOGY
C
      WRITE (NOUT,9040)
      READ (NIN,8010) ELTYP, TOTELS, NODEL
      WRITE (NOUT,9020) ELTYP, TOTELS, NODEL
      DO 1020 I=1,TOTELS
      READ (NIN,8010) ELNUM, (ELTOP(ELNUM,J+2),J=1,NODEL)
      WRITE (NOUT,9020) ELNUM, (ELTOP(ELNUM,J+2),J=1,NODEL)
      ELTOP(ELNUM,1) = ELTYP
      ELTOP(ELNUM,2) = NODEL
 1020 CONTINUE
C
C                            INPUT OF NUMBER OF DEGREES OF FREEDOM
C                            PER NODE, INPUT OF BOUNDARY CON-
C                            DITIONS AND CONSTRUCTION OF NODAL
C                            FREEDOM ARRAY NF
C
      WRITE (NOUT,9050)
      READ (NIN,8010) DOFNOD
      WRITE (NOUT,9020) DOFNOD
      READ (NIN,8010) TOTBND
      WRITE (NOUT,9020) TOTBND
      DO 1030 I=1,TOTBND
      READ (NIN,8010) BTYPE, NUMNOD, NODSID, (BDCND(I,J+3),J=1,
     *     NUMNOD)
      WRITE (NOUT,9020) BTYPE, NUMNOD, NODSID,
     *     (BDCND(I,J+3),J=1,NUMNOD)
      BDCND(I,1) = BTYPE
      BDCND(I,2) = NUMNOD
      BDCND(I,3) = NODSID
 1030 CONTINUE
C
C                            SETUP NODAL FREEDOM ARRAY
C
      TOTDOF = 0
      DO 1050 I=1,TOTNOD
      DO 1040 J=1,DOFNOD
      TOTDOF = TOTDOF + 1
      NF(I,J) = TOTDOF
 1040 CONTINUE
 1050 CONTINUE
C
C                            CALCULATION OF SEMI-BANDWIDTH
C
      FIRST = .TRUE.
      DO 1060 NELE=1,TOTELS
      CALL FREDIF(NELE, ELTOP, IELTOP, JELTOP, NF, INF, JNF,
     *     DOFNOD, FIRST, DIF, ITEST)
 1060 CONTINUE
      HBAND = DIF + 1
C
C*                           ************************************
C*                           *                                  *
C*                           * SYSTEM STIFFNESS MATRIX ASSEMBLY *
C*                           *                                  *
C*                           ************************************
C
      CALL MATNUL(P, IP, JP, DIMEN, DIMEN, ITEST)
      CALL MATNUL(SYSK, ISYSK, JSYSK, TOTDOF, HBAND, ITEST)
      DOFEL = NODEL*DOFNOD
      CALL QQUA4(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)
      DO 1100 NELE=1,TOTELS
      CALL ELGEOM(NELE, ELTOP, IELTOP, JELTOP, COORD, ICOORD,
     *     JCOORD, GEOM, IGEOM, JGEOM, DIMEN, ITEST)
C
C                            INTEGRATION LOOP FOR ELEMENT MATRICES
C                            USING NQP QUADRATURE POINTS
C
      CALL MATNUL(ELK, IELK, JELK, DOFEL, DOFEL, ITEST)
      DO 1090 IQUAD=1,NQP
C
C                            FORM SHAPE FUNCTION AND SPACE
C                            DERIVATIVES IN THE LOCAL CORRDINATES.
C                            TRANSFORM LOCAL DERIVATIVES TO GLOBAL
C                            COORDINATE SYSTEM
C
      XI = ABSS(1,IQUAD)
      ETA = ABSS(2,IQUAD)
      CALL QUAM4(FUN, IFUN, LDER, ILDER, JLDER, XI, ETA, ITEST)
      CALL MATMUL(LDER, ILDER, JLDER, GEOM, IGEOM, JGEOM, JAC,
     *     IJAC, JJAC, DIMEN, NODEL, DIMEN, ITEST)
      CALL MATINV(JAC, IJAC, JJAC, JACIN, IJACIN, JJACIN, DIMEN,
     *     DET, ITEST)
      CALL MATMUL(JACIN, IJACIN, JJACIN, LDER, ILDER, JLDER, GDER,
     *     IGDER, JGDER, DIMEN, DIMEN, NODEL, ITEST)
C
C                           CALCULATE (XI,ETA) IN GLOBAL (X,Y)
C                           COORDINATES AND FORM P MATRIX
C
      CALL SCAPRD(GEOM(1,1), IGEOM, FUN, IFUN, NODEL, X, ITEST)
      CALL SCAPRD(GEOM(1,2), IGEOM, FUN, IFUN, NODEL, Y, ITEST)
      P(1,1) = PX(X,Y)
      P(2,2) = PY(X,Y)
C
C                            FORM INTEGRAND ELEMENT STIFFNESS ELK
C
      CALL MATMUL(P, IP, JP, GDER, IGDER, JGDER, PD, IPD, JPD,
     *     DIMEN, DIMEN, DOFEL, ITEST)
      CALL MATRAN(GDER, IGDER, JGDER, GDERT, IGDERT, JGDERT,
     *     DIMEN, DOFEL, ITEST)
      CALL MATMUL(GDERT, IGDERT, JGDERT, PD, IPD, JPD, DTPD,
     *     IDTPD, JDTPD, DOFEL, DIMEN, DOFEL, ITEST)
      QUOT = DABS(DET)*WGHT(IQUAD)
      DO 1080 I=1,DOFEL
      DO 1070 J=1,DOFEL
      DTPD(I,J) = DTPD(I,J)*QUOT
 1070 CONTINUE
 1080 CONTINUE
      CALL MATADD(ELK, IELK, JELK, DTPD, IDTPD, JDTPD, DOFEL,
     *     DOFEL, ITEST)
 1090 CONTINUE
C
C                            ASSEMBLY OF SYSTEM STIFFNESS MATRIX
C
      CALL DIRECT(NELE, ELTOP, IELTOP, JELTOP, NF, INF, JNF,
     *     DOFNOD, STEER, ISTEER, ITEST)
      CALL ASSYM(SYSK, ISYSK, JSYSK, ELK, IELK, JELK, STEER,
     *     ISTEER, HBAND, DOFEL, ITEST)
 1100 CONTINUE
C
C*                           ************************************
C*                           *                                  *
C*                           * INSERTION OF BOUNDARY CONDITIONS *
C*                           *                                  *
C*                           ************************************
C
      CALL VECNUL(RHS, IRHS, TOTDOF, ITEST)
      CALL QLIN2(WGHT, IWGHT, LABSS, ILABSS, NQP, ITEST)
      DO 1230 ITYPE=1,TOTBND
      BTYPE = BDCND(ITYPE,1)
      NUMNOD = BDCND(ITYPE,2)
      GO TO (1110, 1130, 1130), BTYPE
C
C                            PRESCRIBED VALUES  (DIRICHLET)
C
 1110 DO 1120 J=1,NUMNOD
      K = BDCND(ITYPE,J+3)
      SYSK(K,HBAND) = SYSK(K,HBAND)*SCALE
      X = COORD(K,1)
      Y = COORD(K,2)
      RHS(K) = SYSK(K,HBAND)*H(X,Y)
 1120 CONTINUE
      GO TO 1230
C
C                            DERIVATIVES (NEUMANN AND CAUCHY)
C
 1130 CALL SIDENO(TOTELS, ELTOP, IELTOP, JELTOP, ITYPE, BDCND,
     *     IBDCND, JBDCND, NUMSID, BLIST, IBLIST, JBLIST, ITEST)
      DO 1220 M=1,NUMSID
      ELNUM = BLIST(M,1)
      SIDNUM = BLIST(M,2)
      CALL VECNUL(BELV, IBELV, DOFEL, ITEST)
      CALL MATNUL(BELM, IBELM, JBELM, DOFEL, DOFEL, ITEST)
C
C                            CONSTRUCT QUADRATURE RULE AND LOCAL
C                            GEOMETRY
C
      CALL ELGEOM(ELNUM, ELTOP, IELTOP, JELTOP, COORD, ICOORD,
     *     JCOORD, GEOM, IGEOM, JGEOM, DIMEN, ITEST)
      CALL BQQUA(ABSS, IABSS, JABSS, LABSS, ILABSS, NQP, SIDNUM,
     *     COEFF, ITEST)
C
C                            PERFORM BOUNDARY INTEGRATION
C
      DO 1190 J=1,NQP
      XI = ABSS(1,J)
      ETA = ABSS(2,J)
      CALL QUAM4(FUN, IFUN, LDER, ILDER, JLDER, XI, ETA, ITEST)
      CALL LINQUA(XI, ETA, GEOM, IGEOM, JGEOM, NODEL, SIDNUM,
     *     ULEN, ITEST)
      QUOT = ULEN*WGHT(J)*COEFF
C
C                            CALCULATE (XI,ETA) IN GLOBAL (X,Y)
C                            COORDINATES
C
      CALL SCAPRD(GEOM(1,1), IGEOM, FUN, IFUN, NODEL, X, ITEST)
      CALL SCAPRD(GEOM(1,2), IGEOM, FUN, IFUN, NODEL, Y, ITEST)
C
C                            CALCULATION OF THE NORMAL DIRECTION
C                            COSINES
C
      CALL MATMUL(LDER, ILDER, JLDER, GEOM, IGEOM, JGEOM, JAC,
     *     IJAC, JJAC, DIMEN, NODEL, DIMEN, ITEST)
      CALL MATINV(JAC, IJAC, JJAC, JACIN, IJACIN, JJACIN, DIMEN,
     *     DET, ITEST)
      CALL DCSQUA(JACIN, IJACIN, JJACIN, SIDNUM, COSIN, ICOSIN,
     *     ITEST)
      GO TO (1220, 1140, 1160), BTYPE
C
C                            NEUMANN CONDITIONS
C
 1140 DO 1150 K=1,DOFEL
      BN(K) = (COSIN(1)*P(1,1)*F1(X,Y)+COSIN(2)*P(2,2)*F2(X,Y))*
     *     FUN(K)*QUOT
 1150 CONTINUE
      CALL VECADD(BELV, IBELV, BN, IBN, DOFEL, ITEST)
      GO TO 1190
C
C                            CAUCHY CONDITIONS
C
 1160 DO 1180 K=1,DOFEL
      DO 1170 L=1,DOFEL
      BNTN(K,L) = -FUN(K)*FUN(L)*QUOT*(COSIN(1)*P(1,1)*G1(X,Y)+
     *     COSIN(2)*P(2,2)*G2(X,Y))
 1170 CONTINUE
 1180 CONTINUE
      CALL MATADD(BELM, IBELM, JBELM, BNTN, IBNTN, JBNTN, DOFEL,
     *     DOFEL, ITEST)
 1190 CONTINUE
C
C                            ASSEMBLY OF BOUNDARY CONDITIONS
C
      CALL DIRECT(ELNUM, ELTOP, IELTOP, JELTOP, NF, INF, JNF,
     *     DOFNOD, STEER, ISTEER, ITEST)
      GO TO (1220, 1200, 1210), BTYPE
 1200 CALL ASRHS(RHS, IRHS, BELV, IBELV, STEER, ISTEER, NODEL,
     *     ITEST)
      GO TO 1220
 1210 CALL ASSYM(SYSK, ISYSK, JSYSK, BELM, IBELM, JBELM, STEER,
     *     ISTEER, HBAND, DOFEL, ITEST)
 1220 CONTINUE
 1230 CONTINUE
C
C*                           *********************
C*                           *                   *
C*                           * EQUATION SOLUTION *
C*                           *                   *
C*                           *********************
C
      CALL CHOSOL(SYSK, ISYSK, JSYSK, RHS, IRHS, TOTDOF, HBAND,
     *     ITEST)
      WRITE (NOUT,9060)
      CALL PRTVAL(RHS, IRHS, NF, INF, JNF, DOFNOD, TOTNOD, NOUT,
     *     ITEST)
      STOP
 8010 FORMAT (16I5)
 8020 FORMAT (I5, 6F10.0)
 9010 FORMAT (//25H **** NODAL GEOMETRY ****//1H )
 9020 FORMAT (1H , 16I5)
 9030 FORMAT (1H , I5, 6F10.5)
 9040 FORMAT (//27H **** ELEMENT TOPOLOGY ****//1H )
 9050 FORMAT (//30H **** BOUNDARY CONDITIONS ****//1H )
 9060 FORMAT (//27H **** NODAL POTENTIALS ****//1H )
      END
