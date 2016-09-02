C***********************************************************************
C$SPLIT$TET10$*********************************************************
C***********************************************************************
      SUBROUTINE TET10(FUN, IFUN, DER, IDER, JDER, XI, ETA, ZETA,
     *     ITEST)
C
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS THE VALUES OF SHAPE FUNCTIONS AND DERIVATIVES
C      AT A SPECIFIED POINT FOR AN 10-NODED TETRAHEDRAL ELEMENT.
C      THE FUNCTION IS CONTINUOUS ACROSS ELEMENT BOUNDARIES.
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    10 OCT 1980 (KR)
C
C ARGUMENTS IN
C      IFUN    DIMENSION OF VECTOR FUN (.GE.10)
C      IDER    FIRST DIMENSION OF ARRAY DER (.GE.3)
C      JDER    SECOND DIMENSION OF ARRAY DER (.GE.10)
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
C     SUBROUTINE TET10(FUN, IFUN, DER, IDER, JDER, XI, ETA, ZETA,ITEST)
C***********************************************************************
C
      INTEGER ERRMES, IDER, IERROR, IFUN, ITEST, JDER
      DOUBLE PRECISION DER, DL1X, DL1Y, DL1Z, DL2X, DL2Y, DL2Z,
     *     DL3X, DL3Y, DL3Z, DL4X, DL4Y, DL4Z, DNSL, DNVL,
     *     ETA, FUN, L1, L2, L3, L4, LA, LB, NS, NV, SRNAME,
     *     XI, ZETA
      DIMENSION DER(IDER,JDER), FUN(IFUN)
      DATA SRNAME /8H TET10  /
      NV(LA) = (2.D0*LA-1.D0)*LA
      NS(LA,LB) = 4.D0*LA*LB
      DNVL(LA) = 4.D0*LA - 1.D0
C
      DNSL(LA,LB) = 4.D0*LB
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IFUN.LT.10) IERROR = 1
                        IF (IDER.LT.3 .OR. JDER.LT.10) IERROR = 2
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
C
                        IF (ITEST.NE.0) RETURN
 1010 L1 = 1.D0/12.D0*(3.D0+8.D0*XI-2.D0*DSQRT(2.D0)*ZETA)
      L2 = 1.D0/12.D0*(3.D0-4.D0*(XI+DSQRT(3.D0)*ETA)-2.D0*
     *     DSQRT(2.D0)*ZETA)
      L3 = 1.D0/12.D0*(3.D0-4.D0*(XI-DSQRT(3.D0)*ETA)-2.D0*
     *     DSQRT(2.D0)*ZETA)
C
      L4 = 1.D0/4.D0*(1.D0+2.D0*DSQRT(2.D0)*ZETA)
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
C
      DL4Z = 1.D0/DSQRT(2.D0)
      FUN(1) = NV(L1)
      FUN(2) = NS(L1,L2)
      FUN(3) = NV(L2)
      FUN(4) = NS(L2,L3)
      FUN(5) = NV(L3)
      FUN(6) = NS(L3,L1)
      FUN(7) = NS(L1,L4)
      FUN(8) = NS(L2,L4)
      FUN(9) = NS(L3,L4)
C
      FUN(10) = NV(L4)
      DER(1,1) = DNVL(L1)*DL1X
      DER(2,1) = DNVL(L1)*DL1Y
      DER(3,1) = DNVL(L1)*DL1Z
      DER(1,2) = DNSL(L1,L2)*DL1X + DNSL(L2,L1)*DL2X
      DER(2,2) = DNSL(L1,L2)*DL1Y + DNSL(L2,L1)*DL2Y
      DER(3,2) = DNSL(L1,L2)*DL1Z + DNSL(L2,L1)*DL2Z
      DER(1,3) = DNVL(L2)*DL2X
      DER(2,3) = DNVL(L2)*DL2Y
      DER(3,3) = DNVL(L2)*DL2Z
      DER(1,4) = DNSL(L2,L3)*DL2X + DNSL(L3,L2)*DL3X
      DER(2,4) = DNSL(L2,L3)*DL2Y + DNSL(L3,L2)*DL3Y
      DER(3,4) = DNSL(L2,L3)*DL2Z + DNSL(L3,L2)*DL3Z
      DER(1,5) = DNVL(L3)*DL3X
      DER(2,5) = DNVL(L3)*DL3Y
      DER(3,5) = DNVL(L3)*DL3Z
      DER(1,6) = DNSL(L3,L1)*DL3X + DNSL(L1,L3)*DL1X
      DER(2,6) = DNSL(L3,L1)*DL3Y + DNSL(L1,L3)*DL1Y
      DER(3,6) = DNSL(L3,L1)*DL3Z + DNSL(L1,L3)*DL1Z
      DER(1,7) = DNSL(L1,L4)*DL1X + DNSL(L4,L1)*DL4X
      DER(2,7) = DNSL(L1,L4)*DL1Y + DNSL(L4,L1)*DL4Y
      DER(3,7) = DNSL(L1,L4)*DL1Z + DNSL(L4,L1)*DL4Z
      DER(1,8) = DNSL(L2,L4)*DL2X + DNSL(L4,L2)*DL4X
      DER(2,8) = DNSL(L2,L4)*DL2Y + DNSL(L4,L2)*DL4Y
      DER(3,8) = DNSL(L2,L4)*DL2Z + DNSL(L4,L2)*DL4Z
      DER(1,9) = DNSL(L3,L4)*DL3X + DNSL(L4,L3)*DL4X
      DER(2,9) = DNSL(L3,L4)*DL3Y + DNSL(L4,L3)*DL4Y
      DER(3,9) = DNSL(L3,L4)*DL3Z + DNSL(L4,L3)*DL4Z
      DER(1,10) = DNVL(L4)*DL4X
      DER(2,10) = DNVL(L4)*DL4Y
      DER(3,10) = DNVL(L4)*DL4Z
      RETURN
      END
