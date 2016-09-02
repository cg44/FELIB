C***********************************************************************
C$SPLIT$WDG6$*********************************************************
C***********************************************************************
      SUBROUTINE WDG6(FUN, IFUN, DER, IDER, JDER, XI, ETA, ZETA,
     *     ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS THE VALUES OF SHAPE FUNCTIONS AND DERIVATIVES
C      AT A SPECIFIED POINT FOR AN 6-NODED PENTAHEDRAL ELEMENT.
C      THE FUNCTION IS CONTINUOUS ACROSS ELEMENT BOUNDARIES.
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    10 OCT 1980 (KR)
C
C ARGUMENTS IN
C      IFUN    DIMENSION OF VECTOR FUN (.GE.6)
C      IDER    FIRST DIMENSION OF ARRAY DER (.GE.3)
C      JDER    SECOND DIMENSION OF ARRAY DER (.GE.6)
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
C     SUBROUTINE WDG6(FUN, IFUN, DER, IDER, JDER, XI, ETA, ZETA,ITEST)
C***********************************************************************
C
      INTEGER ERRMES, IDER, IERROR, IFUN, ITEST, JDER
      DOUBLE PRECISION DER, ETA, FUN, L1, L2, L3, SRNAME, XI,
     *     ZETA, ZETAM, ZETAP, XMIN, XMAX, YMIN, YMAX, ZVAL, VEPS,
     *     DUMMY
      DIMENSION DER(IDER,JDER), FUN(IFUN)
      DATA SRNAME /8H WDG6   /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IFUN.LT.6) IERROR = 1
                        IF (IDER.LT.3 .OR. JDER.LT.6) IERROR = 2
                        YMIN = 1.0D0/DSQRT(3.0D0)*(XI-1.0D0)-
     *                       VEPS(DUMMY)
                        YMAX = 1.0D0/DSQRT(3.0D0)*(1.0D0-XI)+
     *                       VEPS(DUMMY)
                        XMIN=-(0.5D0+VEPS(DUMMY))
                        XMAX=1.0D0+VEPS(DUMMY)
                        ZVAL=1.0D0+VEPS(DUMMY)
                        IF ((XI.LT.XMIN .OR. XI.GT.XMAX)
     *                      .OR. (ETA.LT.YMIN .OR. ETA.GT.YMAX)
     *                       .OR. DABS(ZETA) .GT. ZVAL) IERROR = 3
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 L1 = 1.D0/6.D0*(1.D0+2.D0*XI)
      L2 = 1.D0/6.D0*(1.D0-XI-DSQRT(3.D0)*ETA)
      L3 = 1.D0/6.D0*(1.D0-XI+DSQRT(3.D0)*ETA)
      ZETAM = (1.D0-ZETA)
      ZETAP = (1.D0+ZETA)
      FUN(1) = L1*ZETAM
      FUN(2) = L2*ZETAM
      FUN(3) = L3*ZETAM
      FUN(4) = L1*ZETAP
      FUN(5) = L2*ZETAP
      FUN(6) = L3*ZETAP
      DER(1,1) = 1.D0/3.D0*ZETAM
      DER(2,1) = 0.D0
      DER(3,1) = -L1
      DER(1,2) = -1.D0/6.D0*ZETAM
      DER(2,2) = -1.D0/(2.D0*DSQRT(3.D0))*ZETAM
      DER(3,2) = -L2
      DER(1,3) = -1.D0/6.D0*ZETAM
      DER(2,3) = 1.D0/(2.D0*DSQRT(3.D0))*ZETAM
      DER(3,3) = -L3
      DER(1,4) = 1.D0/3.D0*ZETAP
      DER(2,4) = 0.D0
      DER(3,4) = L1
      DER(1,5) = -1.D0/6.D0*ZETAP
      DER(2,5) = -1.D0/(2.D0*DSQRT(3.D0))*ZETAP
      DER(3,5) = L2
      DER(1,6) = -1.D0/6.D0*ZETAP
      DER(2,6) = 1.D0/(2.D0*DSQRT(3.D0))*ZETAP
      DER(3,6) = L3
      RETURN
      END
