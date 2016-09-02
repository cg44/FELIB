C***********************************************************************
C     $Id: seg3p2.f,v 1.2 2009/02/27 11:44:50 cg44 Exp $
C
C     Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
C                          Chilton, Didcot, Oxfordshire OX11 0QX
C
C     Contact: Prof Chris Greenough
C     Email: c.greenough@rl.ac.uk
C     Tel: (44) 1235 445307
C
C     Finite Library Program: Seg3p2 - Potential Flow Around a Small
C                                      Smooth Cylinder
C
C***********************************************************************
C
      PROGRAM SEG3P2
C
      INTEGER BDCND,BLIST,BTYPE,DIF,DIMEN,DOFEL,DOFNOD,ELNUM,ELTOP,
     *        ELTYP,HBAND,I,IABSS,IBDCND,IBELM,IBELV,IBLIST,IBN,IBNTN,
     *        ICOORD,ICOSIN,IDTPD,IELK,IELTOP,IFUN,IGDER,IGDERT,IGEOM,
     *        IJAC,IJACIN,ILABSS,ILDER,INF,IP,IPD,IQUAD,IRHS,ISTEER,
     *        ISYSK,ITEST,ITYPE,IWGHT,J,JABSS,JBDCND,JBELM,JBLIST,JBNTN,
     *        JCOORD,JDTPD,JELK,JELTOP,JGDER,JGDERT,JGEOM,JJAC,JJACIN,
     *        JLDER,JNF,JP,JPD,JSYSK,K,L,M,NELE,NF,NIN,NODEL,NODNUM,
     *        NODSID,NOUT,NQP,NUMNOD,NUMSID,SIDNUM,STEER,TOTBND,TOTDOF,
     *        TOTELS,TOTNOD
      DOUBLE PRECISION ABSS,BELM,BELV,BN,BNTN,COEFF,COORD,COSIN,DET,
     *                 DTPD,ELK,ETA,F1,F2,FUN,G1,G2,GDER,GDERT,GEOM,H,
     *                 JAC,JACIN,LABSS,LDER,P,PD,PX,PY,QUOT,RHS,SCALE,
     *                 SYSK,ULEN,WGHT,X,XI,XX,Y,YY
      LOGICAL FIRST
      DIMENSION ABSS(3,9),BELM(8,8),BELV(8),BN(8),BNTN(8,8),COSIN(3),
     *          DTPD(8,8),ELK(8,8),FUN(8),GDER(3,8),GDERT(8,3),
     *          GEOM(8,3),JAC(3,3),JACIN(3,3),LABSS(3),LDER(3,8),P(3,3),
     *          PD(3,8),STEER(8),WGHT(9)
C
C     Problem size dependent arrays
C
      DIMENSION BDCND(5,40),BLIST(50,5),COORD(100,3),ELTOP(100,10),
     *          NF(100,1),RHS(100),SYSK(100,25)
C
      DATA IABSS/3/,IBELM/8/,IBELV/8/,IBN/8/,IBNTN/8/,ICOSIN/3/,
     *     IDTPD/8/,IELK/8/,IFUN/8/,IGDER/3/,IGDERT/8/,IGEOM/8/,IJAC/3/,
     *     IJACIN/3/,ILABSS/3/,ILDER/3/,IP/3/,IPD/3/,ISTEER/8/,IWGHT/9/,
     *     JABSS/9/,JBELM/8/,JBNTN/8/,JCOORD/3/,JDTPD/8/,JELK/8/,
     *     JGDER/8/,JGDERT/3/,JGEOM/3/,JJAC/3/,JJACIN/3/,JLDER/8/,
     *     JNF/1/,JP/3/,JPD/8/,SCALE/1.0D+10/
C
C     Problem size dependent data statements
C
      DATA IBDCND/5/,IBLIST/50/,ICOORD/100/,IELTOP/100/,INF/100/,
     *     IRHS/100/,ISYSK/100/,JBDCND/40/,JBLIST/5/,JELTOP/10/,
     *     JSYSK/25/
C
      DATA NIN/5/,NOUT/6/
C
C     Statement functions for boundary conditions
C
      H(XX,YY) = 0.D0
      F1(XX,YY) = -1.D0
      F2(XX,YY) = 0.D0
      G1(XX,YY) = 0.D0
      G2(XX,YY) = 0.D0
C
C     Statement functions for PX and PY
C
      PX(XX,YY) = 1.D0
      PY(XX,YY) = 1.D0
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
      WRITE (NOUT,9970)
      READ (NIN,9990) TOTNOD,DIMEN
      WRITE (NOUT,9960) TOTNOD,DIMEN
      DO 1000 I = 1,TOTNOD
          READ (NIN,9980) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
          WRITE (NOUT,9950) NODNUM, (COORD(NODNUM,J),J=1,DIMEN)
 1000 CONTINUE
C
C     Input of element topology
C
      WRITE (NOUT,9940)
      READ (NIN,9990) ELTYP,TOTELS,NODEL
      WRITE (NOUT,9960) ELTYP,TOTELS,NODEL
      DO 1010 I = 1,TOTELS
          READ (NIN,9990) ELNUM, (ELTOP(ELNUM,J+2),J=1,NODEL)
          WRITE (NOUT,9960) ELNUM, (ELTOP(ELNUM,J+2),J=1,NODEL)
          ELTOP(ELNUM,1) = ELTYP
          ELTOP(ELNUM,2) = NODEL
 1010 CONTINUE
C
C     Input of number of degrees of freedom per node, input of
C     boundary conditions and construction of nodal freedom array NF
C
      WRITE (NOUT,9930)
      READ (NIN,9990) DOFNOD
      WRITE (NOUT,9960) DOFNOD
      READ (NIN,9990) TOTBND
      WRITE (NOUT,9960) TOTBND
      DO 1020 I = 1,TOTBND
          READ (NIN,9990) BTYPE,NUMNOD,NODSID,
     *      (BDCND(I,J+3),J=1,NUMNOD)
          WRITE (NOUT,9960) BTYPE,NUMNOD,NODSID,
     *      (BDCND(I,J+3),J=1,NUMNOD)
          BDCND(I,1) = BTYPE
          BDCND(I,2) = NUMNOD
          BDCND(I,3) = NODSID
 1020 CONTINUE
C
C     Setup nodal freedom array
C
      TOTDOF = 0
      DO 1040 I = 1,TOTNOD
          DO 1030 J = 1,DOFNOD
              TOTDOF = TOTDOF + 1
              NF(I,J) = TOTDOF
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
C     *                           * System Stiffness Matrix Assembly *
C     *                           *                                  *
C     *                           ************************************
C
      CALL MATNUL(P,IP,JP,DIMEN,DIMEN,ITEST)
      CALL MATNUL(SYSK,ISYSK,JSYSK,TOTDOF,HBAND,ITEST)
      DOFEL = NODEL*DOFNOD
      CALL QQUA4(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST)
      DO 1090 NELE = 1,TOTELS
          CALL ELGEOM(NELE,ELTOP,IELTOP,JELTOP,COORD,ICOORD,JCOORD,GEOM,
     *                IGEOM,JGEOM,DIMEN,ITEST)
C
C     Integration loop for element matrices using NQP quadrature
C     points
C
          CALL MATNUL(ELK,IELK,JELK,DOFEL,DOFEL,ITEST)
          DO 1080 IQUAD = 1,NQP
C
C     Form shape function and space derivatives in the local
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
C     Calculate (XI,ETA) in global (X,Y) coordinates and form P
C     matrix
C
              CALL SCAPRD(GEOM(1,1),IGEOM,FUN,IFUN,NODEL,X,ITEST)
              CALL SCAPRD(GEOM(1,2),IGEOM,FUN,IFUN,NODEL,Y,ITEST)
              P(1,1) = PX(X,Y)
              P(2,2) = PY(X,Y)
C
C     Form integrand element stiffness ELK
C
              CALL MATMUL(P,IP,JP,GDER,IGDER,JGDER,PD,IPD,JPD,DIMEN,
     *                    DIMEN,DOFEL,ITEST)
              CALL MATRAN(GDER,IGDER,JGDER,GDERT,IGDERT,JGDERT,DIMEN,
     *                    DOFEL,ITEST)
              CALL MATMUL(GDERT,IGDERT,JGDERT,PD,IPD,JPD,DTPD,IDTPD,
     *                    JDTPD,DOFEL,DIMEN,DOFEL,ITEST)
              QUOT = ABS(DET)*WGHT(IQUAD)
              DO 1070 I = 1,DOFEL
                  DO 1060 J = 1,DOFEL
                      DTPD(I,J) = DTPD(I,J)*QUOT
 1060             CONTINUE
 1070         CONTINUE
              CALL MATADD(ELK,IELK,JELK,DTPD,IDTPD,JDTPD,DOFEL,DOFEL,
     *                    ITEST)
 1080     CONTINUE
C
C     Assembly of system stiffness matrix
C
          CALL DIRECT(NELE,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,STEER,
     *                ISTEER,ITEST)
          CALL ASSYM(SYSK,ISYSK,JSYSK,ELK,IELK,JELK,STEER,ISTEER,HBAND,
     *               DOFEL,ITEST)
 1090 CONTINUE
C
C     *                           ************************************
C     *                           *                                  *
C     *                           * Insertion Of Boundary Conditions *
C     *                           *                                  *
C     *                           ************************************
C
      CALL VECNUL(RHS,IRHS,TOTDOF,ITEST)
      CALL QLIN2(WGHT,IWGHT,LABSS,ILABSS,NQP,ITEST)
      DO 1220 ITYPE = 1,TOTBND
          BTYPE = BDCND(ITYPE,1)
          NUMNOD = BDCND(ITYPE,2)
          GO TO (1100,1120,1120) BTYPE
C
C     Prescribed values (dirichlet)
C
 1100     CONTINUE
          DO 1110 J = 1,NUMNOD
              K = BDCND(ITYPE,J+3)
              SYSK(K,HBAND) = SYSK(K,HBAND)*SCALE
              X = COORD(K,1)
              Y = COORD(K,2)
              RHS(K) = SYSK(K,HBAND)*H(X,Y)
 1110     CONTINUE
          GO TO 1220
C
C     Derivatives (neumann and cauchy)
C
 1120     CONTINUE
          CALL SIDENO(TOTELS,ELTOP,IELTOP,JELTOP,ITYPE,BDCND,IBDCND,
     *                JBDCND,NUMSID,BLIST,IBLIST,JBLIST,ITEST)
          DO 1210 M = 1,NUMSID
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
C
C     Perform boundary integration
C
              DO 1180 J = 1,NQP
                  XI = ABSS(1,J)
                  ETA = ABSS(2,J)
                  CALL QUAM4(FUN,IFUN,LDER,ILDER,JLDER,XI,ETA,ITEST)
                  CALL LINQUA(XI,ETA,GEOM,IGEOM,JGEOM,NODEL,SIDNUM,ULEN,
     *                        ITEST)
                  QUOT = ULEN*WGHT(J)*COEFF
C
C     Calculate (XI,ETA) in global (X,Y) coordinates
C
                  CALL SCAPRD(GEOM(1,1),IGEOM,FUN,IFUN,NODEL,X,ITEST)
                  CALL SCAPRD(GEOM(1,2),IGEOM,FUN,IFUN,NODEL,Y,ITEST)
C
C     Calculation of the normal direction cosines
C
                  CALL MATMUL(LDER,ILDER,JLDER,GEOM,IGEOM,JGEOM,JAC,
     *                        IJAC,JJAC,DIMEN,NODEL,DIMEN,ITEST)
                  CALL MATINV(JAC,IJAC,JJAC,JACIN,IJACIN,JJACIN,DIMEN,
     *                        DET,ITEST)
                  CALL DCSQUA(JACIN,IJACIN,JJACIN,SIDNUM,COSIN,ICOSIN,
     *                        ITEST)
                  GO TO (1210,1130,1150) BTYPE
C
C     Neumann conditions
C
 1130             CONTINUE
                  DO 1140 K = 1,DOFEL
                      BN(K) = (COSIN(1)*P(1,1)*F1(X,Y)+
     *                        COSIN(2)*P(2,2)*F2(X,Y))*FUN(K)*QUOT
 1140             CONTINUE
                  CALL VECADD(BELV,IBELV,BN,IBN,DOFEL,ITEST)
                  GO TO 1180
C
C     Cauchy conditions
C
 1150             CONTINUE
                  DO 1170 K = 1,DOFEL
                      DO 1160 L = 1,DOFEL
                          BNTN(K,L) = -FUN(K)*FUN(L)*QUOT*
     *                                (COSIN(1)*P(1,1)*G1(X,Y)+
     *                                COSIN(2)*P(2,2)*G2(X,Y))
 1160                 CONTINUE
 1170             CONTINUE
                  CALL MATADD(BELM,IBELM,JBELM,BNTN,IBNTN,JBNTN,DOFEL,
     *                        DOFEL,ITEST)
 1180         CONTINUE
C
C     Assembly of boundary conditions
C
              CALL DIRECT(ELNUM,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,
     *                    STEER,ISTEER,ITEST)
              GO TO (1210,1190,1200) BTYPE
 1190         CONTINUE
              CALL ASRHS(RHS,IRHS,BELV,IBELV,STEER,ISTEER,NODEL,ITEST)
              GO TO 1210
 1200         CONTINUE
              CALL ASSYM(SYSK,ISYSK,JSYSK,BELM,IBELM,JBELM,STEER,ISTEER,
     *                   HBAND,DOFEL,ITEST)
 1210     CONTINUE
 1220 CONTINUE
C
C     *                           *********************
C     *                           *                   *
C     *                           * Equation Solution *
C     *                           *                   *
C     *                           *********************
C
      CALL CHOSOL(SYSK,ISYSK,JSYSK,RHS,IRHS,TOTDOF,HBAND,ITEST)
      WRITE (NOUT,9920)
      CALL PRTVAL(RHS,IRHS,NF,INF,JNF,DOFNOD,TOTNOD,NOUT,ITEST)
C
C     End of Program
C
      STOP
C
 9990 FORMAT (16I5)
 9980 FORMAT (I5,6F10.0)
 9970 FORMAT (/,/,' **** NODAL GEOMETRY ****',/,/,' ')
 9960 FORMAT (' ',16I5)
 9950 FORMAT (' ',I5,6F10.5)
 9940 FORMAT (/,/,' **** ELEMENT TOPOLOGY ****',/,/,' ')
 9930 FORMAT (/,/,' **** BOUNDARY CONDITIONS ****',/,/,' ')
 9920 FORMAT (/,/,' **** NODAL POTENTIALS ****',/,/,' ')
C
      END
