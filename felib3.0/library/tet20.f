C***********************************************************************
      SUBROUTINE TET20(FUN, IFUN, DER, IDER, JDER, XI, ETA, ZETA,
     *     ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS THE VALUES OF SHAPE FUNCTIONS AND DERIVATIVES
C      AT A SPECIFIED POINT FOR AN 20-NODED TETRAHEDRAL ELEMENT.
C      THE FUNCTION IS CONTINUOUS ACROSS ELEMENT BOUNDARIES.
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 1.1  29 OCT 1979 (CG)
C      COMMENTED    10 OCT 1980 (KR)
C
C ARGUMENTS IN
C      IFUN    DIMENSION OF VECTOR FUN (.GE.20)
C      IDER    FIRST DIMENSION OF ARRAY DER (.GE.3)
C      JDER    SECOND DIMENSION OF ARRAY DER (.GE.20)
C      XI      VALUE OF LOCAL COORDINATE AT WHICH FUNCTION AND
C              DERIVATIVE VALUES REQUIRED
C      ETA     VALUE OF SECOND LOCAL COORDINATE
C      ZETA    VALUE OF THIRD LOCAL COORDINATE
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      FUN     REAL VECTOR OF DIMENSION IFUN.  FUN(I) CONTAINS
C              VALUE OF I'TH SHAPE FUNCTION AT (XI,ETA,ZETA)
C      DER     REAL ARRAY OF DIMENSIONS (IDER,JDER).  DER(I,J)
C              CONTAINS THE DERIVATIVE OF THE J'TH SHAPE
C              FUNCTION WITH RESPECT TO THE I'TH COORDINATE AT
C              THE POINT (XI,ETA,ZETA)
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE TET20(FUN, IFUN, DER, IDER, JDER, XI, ETA, ZETA,
C    *     ITEST)
C***********************************************************************
C
      INTEGER ERRMES, IDER, IERROR, IFUN, ITEST, JDER
      DOUBLE PRECISION DER, DL1X, DL1Y, DL1Z, DL2X, DL2Y, DL2Z,
     *     DL3X, DL3Y, DL3Z, DL4X, DL4Y, DL4Z, DNMFL, DNMSLD,
     *     DNMSLS, DNVL, ETA, FUN, L1, L2, L3, L4, LA, LB,
     *     LC, NMF, NMS, NV, SRNAME, XI, ZETA
      DIMENSION DER(IDER,JDER), FUN(IFUN)
      DATA SRNAME /8H TET20  /
C
C     STATEMENT FUNCTIONS
C
      NV(LA) = 1.D0/2.D0*(3.D0*LA-1.D0)*(3.D0*LA-2.D0)*LA
      NMF(LA,LB,LC) = 27.D0*LA*LB*LC
      NMS(LA,LB) = 9.D0/2.D0*LA*LB*(3.D0*LA-1.D0)
      DNVL(LA) = 1.D0/2.D0*(27.D0*LA*LA-18.D0*LA+2.D0)
      DNMFL(LA,LB,LC) = 27.D0*LB*LC
      DNMSLS(LA,LB) = 9.D0/2.D0*LB*(6.D0*LA-1.D0)
      DNMSLD(LA,LB) = 9.D0/2.D0*LA*(3.D0*LA-1.D0)
C
C     PARAMETER CHECKING
C
      IF (ITEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (IFUN.LT.20) IERROR = 1
      IF (IDER.LT.3 .OR. JDER.LT.20) IERROR = 2
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
C     MAIN BODY
C
 1010 L1 = 1.D0/12.D0*(3.D0+8.D0*XI-2.D0*DSQRT(2.D0)*ZETA)
      L2 = 1.D0/12.D0*(3.D0-4.D0*(XI+DSQRT(3.D0)*ETA)-2.D0*
     *     DSQRT(2.D0)*ZETA)
      L3 = 1.D0/12.D0*(3.D0-4.D0*(XI-DSQRT(3.D0)*ETA)-2.D0*
     *     DSQRT(2.D0)*ZETA)
      L4 = 1.D0/4.D0*(1.D0+2.D0*DSQRT(2.D0)*ZETA)
C
      DL1X = 2.D0/3.D0
      DL1Y = 0.D0
      DL1Z = -1.D0/(3.D0*DSQRT(2.D0))
      DL2X = -1.D0/3.D0
      DL2Y = -1.D0/DSQRT(3.D0)
      DL2Z = -1.D0/(3.D0*DSQRT(2.D0))
      DL3X = -1.D0/3.D0
      DL3Y = 1.D0/DSQRT(3.D0)
      DL3Z = -1.D0/(3.D0*DSQRT(2.D0))
      DL4X = 0.D0
      DL4Y = 0.D0
      DL4Z = 1.D0/DSQRT(2.D0)
C
C     SHAPE FUNCTIONS
C
      FUN(1) = NV(L1)
      FUN(2) = NMS(L1,L2)
      FUN(3) = NMS(L2,L1)
      FUN(4) = NV(L2)
      FUN(5) = NMS(L2,L3)
      FUN(6) = NMS(L3,L2)
      FUN(7) = NV(L3)
      FUN(8) = NMS(L3,L1)
      FUN(9) = NMS(L1,L3)
      FUN(10) = NMF(L1,L2,L3)
      FUN(11) = NMS(L1,L4)
      FUN(12) = NMF(L4,L1,L2)
      FUN(13) = NMS(L2,L4)
      FUN(14) = NMF(L2,L3,L4)
      FUN(15) = NMS(L3,L4)
      FUN(16) = NMF(L3,L4,L1)
      FUN(17) = NMS(L4,L1)
      FUN(18) = NMS(L4,L2)
      FUN(19) = NMS(L4,L3)
      FUN(20) = NV(L4)
C
C     DERIVATIVES
C
      DER(1,1) = DNVL(L1)*DL1X
      DER(2,1) = DNVL(L1)*DL1Y
      DER(3,1) = DNVL(L1)*DL1Z
      DER(1,2) = DNMSLS(L1,L2)*DL1X + DNMSLD(L1,L2)*DL2X
      DER(2,2) = DNMSLS(L1,L2)*DL1Y + DNMSLD(L1,L2)*DL2Y
      DER(3,2) = DNMSLS(L1,L2)*DL1Z + DNMSLD(L1,L2)*DL2Z
      DER(1,3) = DNMSLD(L2,L1)*DL1X + DNMSLS(L2,L1)*DL2X
      DER(2,3) = DNMSLD(L2,L1)*DL1Y + DNMSLS(L2,L1)*DL2Y
      DER(3,3) = DNMSLD(L2,L1)*DL1Z + DNMSLS(L2,L1)*DL2Z
      DER(1,4) = DNVL(L2)*DL2X
      DER(2,4) = DNVL(L2)*DL2Y
      DER(3,4) = DNVL(L2)*DL2Z
      DER(1,5) = DNMSLS(L2,L3)*DL2X + DNMSLD(L2,L3)*DL3X
      DER(2,5) = DNMSLS(L2,L3)*DL2Y + DNMSLD(L2,L3)*DL3Y
      DER(3,5) = DNMSLS(L2,L3)*DL2Z + DNMSLD(L2,L3)*DL3Z
      DER(1,6) = DNMSLD(L3,L2)*DL2X + DNMSLS(L3,L2)*DL3X
      DER(2,6) = DNMSLD(L3,L2)*DL2Y + DNMSLS(L3,L2)*DL3Y
      DER(3,6) = DNMSLD(L3,L2)*DL2Z + DNMSLS(L3,L2)*DL3Z
      DER(1,7) = DNVL(L3)*DL3X
      DER(2,7) = DNVL(L3)*DL3Y
      DER(3,7) = DNVL(L3)*DL3Z
      DER(1,8) = DNMSLS(L3,L1)*DL3X + DNMSLD(L3,L1)*DL1X
      DER(2,8) = DNMSLS(L3,L1)*DL3Y + DNMSLD(L3,L1)*DL1Y
      DER(3,8) = DNMSLS(L3,L1)*DL3Z + DNMSLD(L3,L1)*DL1Z
      DER(1,9) = DNMSLD(L1,L3)*DL3X + DNMSLS(L1,L3)*DL1X
      DER(2,9) = DNMSLD(L1,L3)*DL3Y + DNMSLS(L1,L3)*DL1Y
      DER(3,9) = DNMSLD(L1,L3)*DL3Z + DNMSLS(L1,L3)*DL1Z
      DER(1,10) = DNMFL(L1,L2,L3)*DL1X + DNMFL(L2,L3,L1)*DL2X +
     *     DNMFL(L3,L1,L2)*DL3X
      DER(2,10) = DNMFL(L1,L2,L3)*DL1Y + DNMFL(L2,L3,L1)*DL2Y +
     *     DNMFL(L3,L1,L2)*DL3Y
      DER(3,10) = DNMFL(L1,L2,L3)*DL1Z + DNMFL(L2,L3,L1)*DL2Z +
     *     DNMFL(L3,L1,L2)*DL3Z
      DER(1,11) = DNMSLS(L1,L4)*DL1X + DNMSLD(L1,L4)*DL4X
      DER(2,11) = DNMSLS(L1,L4)*DL1Y + DNMSLD(L1,L4)*DL4Y
      DER(3,11) = DNMSLS(L1,L4)*DL1Z + DNMSLD(L1,L4)*DL4Z
      DER(1,12) = DNMFL(L1,L2,L4)*DL1X + DNMFL(L2,L1,L4)*DL2X +
     *     DNMFL(L4,L1,L2)*DL4X
      DER(2,12) = DNMFL(L1,L2,L4)*DL1Y + DNMFL(L2,L1,L4)*DL2Y +
     *     DNMFL(L4,L1,L2)*DL4Y
      DER(3,12) = DNMFL(L1,L2,L4)*DL1Z + DNMFL(L2,L1,L4)*DL2Z +
     *     DNMFL(L4,L1,L2)*DL4Z
      DER(1,13) = DNMSLS(L2,L4)*DL2X + DNMSLD(L2,L4)*DL4X
      DER(2,13) = DNMSLS(L2,L4)*DL2Y + DNMSLD(L2,L4)*DL4Y
      DER(3,13) = DNMSLS(L2,L4)*DL2Z + DNMSLD(L2,L4)*DL4Z
      DER(1,14) = DNMFL(L2,L3,L4)*DL2X + DNMFL(L3,L2,L4)*DL3X +
     *     DNMFL(L4,L2,L3)*DL4X
      DER(2,14) = DNMFL(L2,L3,L4)*DL2Y + DNMFL(L3,L2,L4)*DL3Y +
     *     DNMFL(L4,L2,L3)*DL4Y
      DER(3,14) = DNMFL(L2,L3,L4)*DL2Z + DNMFL(L3,L2,L4)*DL3Z +
     *     DNMFL(L4,L2,L3)*DL4Z
      DER(1,15) = DNMSLS(L3,L4)*DL3X + DNMSLD(L3,L4)*DL4X
      DER(2,15) = DNMSLS(L3,L4)*DL3Y + DNMSLD(L3,L4)*DL4Y
      DER(3,15) = DNMSLS(L3,L4)*DL3Z + DNMSLD(L3,L4)*DL4Z
      DER(1,16) = DNMFL(L1,L3,L4)*DL1X + DNMFL(L3,L1,L4)*DL3X +
     *     DNMFL(L4,L1,L3)*DL4X
      DER(2,16) = DNMFL(L1,L3,L4)*DL1Y + DNMFL(L3,L1,L4)*DL3Y +
     *     DNMFL(L4,L1,L3)*DL4Y
      DER(3,16) = DNMFL(L1,L3,L4)*DL1Z + DNMFL(L3,L1,L4)*DL3Z +
     *     DNMFL(L4,L1,L3)*DL4Z
      DER(1,17) = DNMSLD(L4,L1)*DL1X + DNMSLS(L4,L1)*DL4X
      DER(2,17) = DNMSLD(L4,L1)*DL1Y + DNMSLS(L4,L1)*DL4Y
      DER(3,17) = DNMSLD(L4,L1)*DL1Z + DNMSLS(L4,L1)*DL4Z
      DER(1,18) = DNMSLD(L4,L2)*DL2X + DNMSLS(L4,L2)*DL4X
      DER(2,18) = DNMSLD(L4,L2)*DL2Y + DNMSLS(L4,L2)*DL4Y
      DER(3,18) = DNMSLD(L4,L2)*DL2Z + DNMSLS(L4,L2)*DL4Z
      DER(1,19) = DNMSLD(L4,L3)*DL3X + DNMSLS(L4,L3)*DL4X
      DER(2,19) = DNMSLD(L4,L3)*DL3Y + DNMSLS(L4,L3)*DL4Y
      DER(3,19) = DNMSLD(L4,L3)*DL3Z + DNMSLS(L4,L3)*DL4Z
      DER(1,20) = DNVL(L4)*DL4X
      DER(2,20) = DNVL(L4)*DL4Y
      DER(3,20) = DNVL(L4)*DL4Z
C
      RETURN
      END
C***********************************************************************
