C***********************************************************************
C     $Id: seg1p2.f,v 1.2 2009/02/27 11:44:50 cg44 Exp $
C
C     Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
C                          Chilton, Didcot, Oxfordshire OX11 0QX
C
C     Contact: Prof Chris Greenough
C     Email: c.greenough@rl.ac.uk
C     Tel: (44) 1235 445307
C
C     Finite Library Program: Seg1p2 - Axisymmetric Strain of an
C                                      Elastic Solid
C
C***********************************************************************
C
      PROGRAM SEG1P2
C
      INTEGER BDCND,BLIST,BTYPE,DIMEN,DOFEL,DOFNOD,ELNUM,ELTOP,HBAND,I,
     *        IABS1D,IABSS,IB,IBDCND,IBLIST,IBT,IBTDB,ICOORD,ID,IDB,
     *        IELDIS,IELK,IELP,IELTOP,IFUN,IGDER,IGEOM,IJAC,IJACIN,
     *        ILDER,ILOADS,INF,IPOS,IPROP,IPRSSR,IPRVEC,IQUAD,IRESTR,
     *        ISTEER,ISTRN,ISTRS,ISYSK,ITEST,ITYPE,IWGHT,J,JABSS,JB,
     *        JBDCND,JBLIST,JBT,JBTDB,JCOORD,JD,JDB,JELK,JELTOP,JGDER,
     *        JGEOM,JJAC,JJACIN,JLDER,JNF,JRESTR,JSYSK,K,L,LODNOD,
     *        MATTYP,NELE,NF,NIN,NODEL,NODNUM,NODSID,NOUT,NQP,NTEMP,
     *        NUMMAT,NUMNOD,NUMSID,NUMSS,RESNOD,RESTR,SIDNUM,STEER,
     *        TOTDOF,TOTELS,TOTNOD
      DOUBLE PRECISION ABS1D,ABSS,B,BT,BTDB,COEFF,COORD,D,DB,DET,E,
     *                 ELDIS,ELK,ELP,ETA,FUN,GDER,GEOM,JAC,JACIN,LDER,
     *                 LOADS,NU,PI,PROP,PRSSR,PRVEC,QUOT,RBAR,STRN,STRS,
     *                 SYSK,ULEN,WGHT,XI
      DIMENSION ABS1D(3),ABSS(3,9),B(6,24),BT(24,6),BTDB(24,24),D(6,6),
     *          DB(6,24),ELDIS(24),ELK(24,24),ELP(8),FUN(8),GDER(3,8),
     *          GEOM(8,3),JAC(3,3),JACIN(3,3),LDER(3,8),PRVEC(8),
     *          STEER(24),STRN(6),STRS(6),WGHT(9)
C
C     Problem size dependent arrays
C
      DIMENSION BDCND(10,10),BLIST(10,2),COORD(100,3),ELTOP(100,10),
     *          LOADS(100),NF(100,3),PROP(24,2),PRSSR(10,2),
     *          RESTR(100,4),SYSK(100,30)
C
      DATA IABS1D/3/,IABSS/3/,IB/6/,IBT/24/,IBTDB/24/,ID/6/,IDB/6/,
     *     IELDIS/24/,IELK/24/,IELP/8/,IFUN/8/,IGDER/3/,IGEOM/8/,
     *     IJAC/3/,IJACIN/3/,ILDER/3/,IPROP/24/,IPRVEC/8/,ISTEER/24/,
     *     ISTRN/6/,ISTRS/6/,IWGHT/9/,JABSS/9/,JB/24/,JBT/6/,JBTDB/24/,
     *     JCOORD/3/,JD/6/,JDB/24/,JELK/24/,JGDER/8/,JGEOM/3/,JJAC/3/,
     *     JJACIN/3/,JLDER/8/,JNF/3/,JRESTR/4/,NODSID/2/
C
C     Problem size dependent data statements
C
      DATA IBDCND/10/,IBLIST/10/,ICOORD/100/,IELTOP/100/,ILOADS/100/,
     *     INF/100/,IPRSSR/10/,IRESTR/100/,ISYSK/100/,JBDCND/10/,
     *     JBLIST/2/,JELTOP/10/,JSYSK/30/
C
      DATA NIN/5/,NOUT/6/,NTEMP/7/
C
C     Set ITEST for full checking
C
      ITEST = 0
      PI = 4.D0*ATAN(1.D0)
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
          READ (NIN,9990) ELNUM,MATTYP,NODEL,
     *      (ELTOP(ELNUM,J+2),J=1,NODEL)
          WRITE (NOUT,9950) ELNUM,MATTYP,NODEL,
     *      (ELTOP(ELNUM,J+2),J=1,NODEL)
          ELTOP(ELNUM,1) = MATTYP
          ELTOP(ELNUM,2) = NODEL
 1010 CONTINUE
C
C     Input of material properties and construction of stress-strain
C     matrix D for plane strain
C
      WRITE (NOUT,9920)
      READ (NIN,9990) NUMMAT
      WRITE (NOUT,9910) NUMMAT
C
      READ (NIN,9970) (PROP(I,1),I=1,NUMMAT)
      READ (NIN,9970) (PROP(I,2),I=1,NUMMAT)
      WRITE (NOUT,9900) ((PROP(I,J),I=1,NUMMAT),J=1,2)
C
C     Input of number of degrees of freedom per node, input of
C     restrained node data and construction of nodal freedom array
C     NF
C
      WRITE (NOUT,9890)
      READ (NIN,9990) DOFNOD
      WRITE (NOUT,9950) DOFNOD
C
      READ (NIN,9990) RESNOD
      WRITE (NOUT,9950) RESNOD
C
      K = DOFNOD + 1
      DO 1020 I = 1,RESNOD
          READ (NIN,9990) (RESTR(I,J),J=1,K)
          WRITE (NOUT,9950) (RESTR(I,J),J=1,K)
 1020 CONTINUE
      CALL FORMNF(RESTR,IRESTR,JRESTR,RESNOD,TOTNOD,DOFNOD,NF,INF,JNF,
     *            TOTDOF,ITEST)
C
C     Loading data input
C
      WRITE (NOUT,9880)
      READ (NIN,9990) LODNOD
      WRITE (NOUT,9950) LODNOD
C
      DO 1030 ITYPE = 1,LODNOD
          READ (NIN,9990) BTYPE,NODSID,NUMNOD,
     *      (BDCND(ITYPE,J+3),J=1,NUMNOD)
          READ (NIN,9970) (PRSSR(ITYPE,I),I=1,DOFNOD)
          WRITE (NOUT,9950) BTYPE,NODSID,NUMNOD,
     *      (BDCND(ITYPE,J+3),J=1,NUMNOD)
          WRITE (NOUT,9900) (PRSSR(ITYPE,I),I=1,DOFNOD)
          BDCND(ITYPE,1) = BTYPE
          BDCND(ITYPE,2) = NUMNOD
          BDCND(ITYPE,3) = NODSID
 1030 CONTINUE
C
C     Calculation of semi-bandwidth
C
      CALL BNDWTH(ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,TOTELS,HBAND,
     *            ITEST)
C
C     *                           ************************************
C     *                           *                                  *
C     *                           * Global Stiffness Matrix Assembly *
C     *                           *                                  *
C     *                           ************************************
C
      CALL MATNUL(SYSK,ISYSK,JSYSK,TOTDOF,HBAND,ITEST)
      CALL QQUA4(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST)
C
      DO 1070 NELE = 1,TOTELS
          CALL ELGEOM(NELE,ELTOP,IELTOP,JELTOP,COORD,ICOORD,JCOORD,GEOM,
     *                IGEOM,JGEOM,DIMEN,ITEST)
C
C     Strain-displacement matrix
C
          MATTYP = ELTOP(NELE,1)
          E = PROP(MATTYP,1)
          NU = PROP(MATTYP,2)
          NODEL = ELTOP(NELE,2)
          DOFEL = NODEL*DOFNOD
          CALL DAXI(D,ID,JD,E,NU,NUMSS,ITEST)
C
C     Integration loop for element stiffness using NQP quadrature
C     points
C
          CALL MATNUL(ELK,IELK,JELK,DOFEL,DOFEL,ITEST)
          DO 1060 IQUAD = 1,NQP
C
C     Form linear shape function and space derivatives in the local
C     corrdinates.
C
              XI = ABSS(1,IQUAD)
              ETA = ABSS(2,IQUAD)
              CALL QUAM4(FUN,IFUN,LDER,ILDER,JLDER,XI,ETA,ITEST)
C
C     Calculate radius at gauss point
C
              CALL SCAPRD(GEOM(1,1),IGEOM,FUN,IFUN,NODEL,RBAR,ITEST)
C
C     Transform local derivatives to global coordinate system
C
              CALL MATMUL(LDER,ILDER,JLDER,GEOM,IGEOM,JGEOM,JAC,IJAC,
     *                    JJAC,DIMEN,NODEL,DIMEN,ITEST)
              CALL MATINV(JAC,IJAC,JJAC,JACIN,IJACIN,JJACIN,DIMEN,DET,
     *                    ITEST)
              CALL MATMUL(JACIN,IJACIN,JJACIN,LDER,ILDER,JLDER,GDER,
     *                    IGDER,JGDER,DIMEN,DIMEN,NODEL,ITEST)
C
C     Formation of strain-displacement matrix B and output to work
C     file for later recovery process
C
              CALL B2P2(B,IB,JB,GDER,IGDER,JGDER,FUN,IFUN,GEOM,IGEOM,
     *                  JGEOM,NODEL,ITEST)
              WRITE (NTEMP) ((B(I,J),I=1,NUMSS),J=1,DOFEL)
C
C     Formation of element stiffness ELK
C
              CALL MATMUL(D,ID,JD,B,IB,JB,DB,IDB,JDB,NUMSS,NUMSS,DOFEL,
     *                    ITEST)
              CALL MATRAN(B,IB,JB,BT,IBT,JBT,NUMSS,DOFEL,ITEST)
              CALL MATMUL(BT,IBT,JBT,DB,IDB,JDB,BTDB,IBTDB,JBTDB,DOFEL,
     *                    NUMSS,DOFEL,ITEST)
              QUOT = 2.D0*PI*RBAR*ABS(DET)*WGHT(IQUAD)
              DO 1050 I = 1,DOFEL
                  DO 1040 J = 1,DOFEL
                      BTDB(I,J) = BTDB(I,J)*QUOT
 1040             CONTINUE
 1050         CONTINUE
              CALL MATADD(ELK,IELK,JELK,BTDB,IBTDB,JBTDB,DOFEL,DOFEL,
     *                    ITEST)
 1060     CONTINUE
C
C     Form system matrix SYSK
C
          CALL DIRECT(NELE,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,STEER,
     *                ISTEER,ITEST)
          CALL ASSYM(SYSK,ISYSK,JSYSK,ELK,IELK,JELK,STEER,ISTEER,HBAND,
     *               DOFEL,ITEST)
 1070 CONTINUE
C
C     *                          ******************************
C     *                          *                            *
C     *                          * Assembly Of Pressure LOADS *
C     *                          *                            *
C     *                          ******************************
C
      CALL VECNUL(LOADS,ILOADS,TOTDOF,ITEST)
      CALL QLIN2(WGHT,IWGHT,ABS1D,IABS1D,NQP,ITEST)
C
      DO 1120 ITYPE = 1,LODNOD
C
C     Construct boundary list
C
          CALL SIDENO(TOTELS,ELTOP,IELTOP,JELTOP,ITYPE,BDCND,IBDCND,
     *                JBDCND,NUMSID,BLIST,IBLIST,JBLIST,ITEST)
          DO 1110 J = 1,NUMSID
              ELNUM = BLIST(J,1)
              SIDNUM = BLIST(J,2)
              NODEL = ELTOP(ELNUM,2)
              DOFEL = NODEL*DOFNOD
              CALL VECNUL(ELP,IELP,DOFEL,ITEST)
              CALL ELGEOM(ELNUM,ELTOP,IELTOP,JELTOP,COORD,ICOORD,JCOORD,
     *                    GEOM,IGEOM,JGEOM,DIMEN,ITEST)
              CALL BQQUA(ABSS,IABSS,JABSS,ABS1D,IABS1D,NQP,SIDNUM,COEFF,
     *                   ITEST)
C
C     Integrate along boundary
C
              DO 1100 IQUAD = 1,NQP
                  XI = ABSS(1,IQUAD)
                  ETA = ABSS(2,IQUAD)
                  CALL QUAM4(FUN,IFUN,LDER,ILDER,JLDER,XI,ETA,ITEST)
                  CALL LINQUA(XI,ETA,GEOM,IGEOM,JGEOM,NODEL,SIDNUM,ULEN,
     *                        ITEST)
                  CALL SCAPRD(GEOM(1,1),IGEOM,FUN,IFUN,NODEL,RBAR,ITEST)
                  QUOT = 2.D0*PI*RBAR*ULEN*WGHT(IQUAD)*COEFF
                  DO 1090 K = 1,NODEL
                      DO 1080 L = 1,DOFNOD
                          IPOS = K*2 + L - 2
                          PRVEC(IPOS) = FUN(K)*PRSSR(ITYPE,L)*QUOT
 1080                 CONTINUE
 1090             CONTINUE
                  CALL VECADD(ELP,IELP,PRVEC,IPRVEC,DOFEL,ITEST)
 1100         CONTINUE
C
C     Assemble boundary contributions
C
              CALL DIRECT(ELNUM,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,
     *                    STEER,ISTEER,ITEST)
              CALL ASRHS(LOADS,ILOADS,ELP,IELP,STEER,ISTEER,DOFEL,ITEST)
 1110     CONTINUE
 1120 CONTINUE
C
C     *                           *********************
C     *                           *                   *
C     *                           * Equation Solution *
C     *                           *                   *
C     *                           *********************
C
      WRITE (NOUT,9870)
      CALL PRTVAL(LOADS,ILOADS,NF,INF,JNF,DOFNOD,TOTNOD,NOUT,ITEST)
      CALL CHOSOL(SYSK,ISYSK,JSYSK,LOADS,ILOADS,TOTDOF,HBAND,ITEST)
      WRITE (NOUT,9860)
      CALL PRTVAL(LOADS,ILOADS,NF,INF,JNF,DOFNOD,TOTNOD,NOUT,ITEST)
C
C     *                           **************************
C     *                           *                        *
C     *                           * Stress-strain Recovery *
C     *                           *                        *
C     *                           **************************
C
      REWIND NTEMP
      CALL QQUA4(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST)
C
      DO 1140 NELE = 1,TOTELS
          WRITE (NOUT,9850) NELE
C
C     SELECT nodal displacements for element NELE using the steering
C     vector
C
          CALL DIRECT(NELE,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,STEER,
     *                ISTEER,ITEST)
          DO 1130 IQUAD = 1,NQP
              CALL SELECT(LOADS,ILOADS,STEER,ISTEER,DOFEL,ELDIS,IELDIS,
     *                    ITEST)
C
C     Recover B matrix for point IQUAD and calculate the stresses
C     and strains
C
              READ (NTEMP) ((B(I,J),I=1,NUMSS),J=1,DOFEL)
              CALL MATVEC(B,IB,JB,ELDIS,IELDIS,NUMSS,DOFEL,STRN,ISTRN,
     *                    ITEST)
              CALL MATVEC(D,ID,JD,STRN,ISTRN,NUMSS,NUMSS,STRS,ISTRS,
     *                    ITEST)
              WRITE (NOUT,9840) IQUAD
              CALL PRTVEC(STRN,ISTRN,NUMSS,NOUT,ITEST)
              CALL PRTVEC(STRS,ISTRS,NUMSS,NOUT,ITEST)
 1130     CONTINUE
 1140 CONTINUE
C
C     Program end
C
      STOP
C
 9990 FORMAT (16I5)
 9980 FORMAT (I5,6F10.0)
 9970 FORMAT (2F10.0)
 9960 FORMAT (/,/,' **** NODAL GEOMETRY ****',/,' ')
 9950 FORMAT (' ',16I5)
 9940 FORMAT (' ',I5,6F10.5)
 9930 FORMAT (/,/,' **** ELEMENT TOPOLOGY ****',/,' ')
 9920 FORMAT (/,/,' **** MATERIAL PROPERTIES ****',/,' ')
 9910 FORMAT (/,' NUMBER OF MATERIAL TYPES =',I4)
 9900 FORMAT (' ',2F10.5)
 9890 FORMAT (/,/,' **** RESTRAINT DATA ****',/,' ')
 9880 FORMAT (/,/,' **** LOADING DATA ****',/,' ')
 9870 FORMAT (/,/,' **** VECTOR OF LOADS ****',/,' ')
 9860 FORMAT (/,/,' **** EQUILIBRIUM DISPLACEMENTS (U,V) ****',/,' ')
 9850 FORMAT (/,/,' **** STRAINS AND STRESSES FOR ELEMENT ',I3,' ****',
     *       /,' ')
 9840 FORMAT (/,' QUADRATURE POINT ',I3,/,' ')
C
      END
