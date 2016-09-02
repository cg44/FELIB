C***********************************************************************
C$SPLIT$LI3FN$*********************************************************
C***********************************************************************
      SUBROUTINE LI3FN(FUN, IFUN, DER, IDER, XI, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS THE VALUES OF THE SHAPE FUNCTION AND ITS
C      DERIVATIVE AT A SPECIFIED POINT FOR A 3-NODED C0
C      CONTINUOUS LINE ELEMENT.  FUNCTION CONTINUOUS ACROSS
C      ELEMENT BOUNDARIES
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    11 OCT 1980 (KR)
C
C ARGUMENTS IN
C      IFUN    DIMENSION OF VECTOR FUN (.GE.3)
C      IDER    DIMENSION OF VECTOR DER (.GE.3)
C      XI      VALUE OF LOCAL COORDINATE
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      FUN     VECTOR OF LENGTH IFUN.  FUN(I) CONTAINS THE
C              VALUE OF THE I'TH SHAPE FUNCTION AT THE POINT XI
C      DER     VECTOR OF LENGTH IDER.  DER(I) CONTAINS THE
C              VALUE OF THE DERIVATIVE OF THE I'TH SHAPE
C              FUNCTION AT THE POINT XI
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE LI3FN(FUN, IFUN, DER, IDER, XI, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, IDER, IERROR, IFUN, ITEST
      DOUBLE PRECISION DER, FUN, SRNAME, XI, XMAX, VEPS, DUMMY
      DIMENSION DER(IDER), FUN(IFUN)
      DATA SRNAME /8H LI3FN  /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IFUN.LT.3 .OR. IDER.LT.3) IERROR = 1
                        XMAX = 1.0D0 + VEPS(DUMMY)
                        IF (DABS(XI) .GT. XMAX) IERROR = 2
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 FUN(1) = 0.5D0*XI*(XI-1.0D0)
      FUN(2) = 1.0D0 - XI*XI
      FUN(3) = 0.5D0*XI*(XI+1.0D0)
      DER(1) = 0.5D0*(2.0D0*XI-1.0D0)
      DER(2) = -2.0D0*XI
      DER(3) = 0.5D0*(2.0D0*XI+1.0D0)
      RETURN
      END
