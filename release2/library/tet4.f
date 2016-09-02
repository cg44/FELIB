C***********************************************************************
C$SPLIT$TET4$*********************************************************
C***********************************************************************
      SUBROUTINE TET4(FUN, IFUN, DER, IDER, JDER, XI, ETA, ZETA,
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
C      IFUN    DIMENSION OF VECTOR FUN (.GE.4)
C      IDER    FIRST DIMENSION OF ARRAY DER (.GE.3)
C      JDER    SECOND DIMENSION OF ARRAY DER (.GE.4)
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
C     SUBROUTINE TET4(FUN, IFUN, DER, IDER, JDER, XI, ETA, ZETA,ITEST)
C***********************************************************************
C
      INTEGER ERRMES, IDER, IERROR, IFUN, ITEST, JDER
      DOUBLE PRECISION DER, ETA, FUN, SRNAME, XI, ZETA
      DIMENSION DER(IDER,JDER), FUN(IFUN)
      DATA SRNAME /8H TET4   /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IFUN.LT.4) IERROR = 1
                        IF (IDER.LT.3 .OR. JDER.LT.4) IERROR = 2
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 FUN(1) = 1.D0/12.D0*(3.D0+8.D0*XI-2.D0*DSQRT(2.D0)*ZETA)
      FUN(2) = 1.D0/12.D0*(3.D0-4.D0*(XI+DSQRT(3.D0)*ETA)-2.D0*
     *     DSQRT(2.D0)*ZETA)
      FUN(3) = 1.D0/12.D0*(3.D0-4.D0*(XI-DSQRT(3.D0)*ETA)-2.D0*
     *     DSQRT(2.D0)*ZETA)
      FUN(4) = 1.D0/4.D0*(1.D0+2.D0*DSQRT(2.D0)*ZETA)
      DER(1,1) = 2.D0/3.D0
C     DER(1,3)=-DSQRT(2.)/6.
      DER(1,2) = 0.D0
      DER(1,3) = -DSQRT(2.D0)/6.D0
C     DER(2,2)=-DSQRT(3.)/3.
      DER(2,1) = -1.D0/3.D0
C     DER(2,3)=-DSQRT(2.)/6.
      DER(2,2) = -DSQRT(3.D0)/3.D0
      DER(2,3) = -DSQRT(2.D0)/6.D0
C     DER(3,2)=+DSQRT(2.)/3.
      DER(3,1) = -1.D0/3.D0
C     DER(3,3)=-DSQRT(2.)/6.
      DER(3,2) = +DSQRT(2.D0)/3.D0
      DER(3,3) = -DSQRT(2.D0)/6.D0
      DER(4,1) = 0.D0
C     DER(4,3)=DSQRT(2.)/2.
      DER(4,2) = 0.D0
      DER(4,3) = DSQRT(2.D0)/2.D0
      RETURN
      END
