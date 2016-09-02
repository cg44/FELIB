C***********************************************************************
C$SPLIT$ROD4$*********************************************************
C***********************************************************************
      SUBROUTINE ROD4 (FUN, IFUN, DER, IDER, XI, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS THE VALUES OF THE SHAPE FUNCTIONS AND THEIR
C      DERIVATIVES AT A SPECIFIED POINT FOR A 4-NODED C0
C      CONTINUOUS ELEMENT.  FUNCTION CONTINUOUS ACROSS ELEMENT
C      BOUNDARIES
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    11 OCT 1980 (KR)
C
C ARGUMENTS IN
C      IFUN    DIMENSION OF VECTOR FUN (.GE.4)
C      IDER    DIMENSION OF VECTOR DER (.GE.4)
C      XI      VALUE OF LOCAL COORDINATE
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      FUN     VECTOR OF DIMENSION IFUN.  FUN(I) CONTAINS THE
C              VALUE OF THE I'TH SHAPE FUNCTION AT XI
C      DER     VECTOR OF DIMENSION IDER.  DER(I) CONTAINS THE
C              VALUE OF THE DERIVATIVE OF THE I'TH SHAPE
C              FUNCTION AT XI
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE ROD4 (FUN, IFUN, DER, IDER, XI, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, IDER, IERROR, IFUN, ITEST
      DOUBLE PRECISION DER, FUN, SRNAME, XI, XMAX, VEPS, DUMMY
      DIMENSION DER(IDER), FUN(IFUN)
      DATA SRNAME /8H ROD4   /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IFUN.LT.4 .OR. IDER.LT.4) IERROR = 1
                        XMAX = 1.0D0 + VEPS(DUMMY)
                        IF (DABS(XI) .GT. XMAX) IERROR = 2
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 FUN(1) = 1.0D0/16.0D0*(1.0D0-XI)*(9.0D0*XI*XI-1.0D0)
      FUN(2) = 9.0D0/16.0D0*(XI*XI-1.0D0)*(3.0D0*XI-1.0D0)
      FUN(3) = 9.0D0/16.0D0*(1.0D0-XI*XI)*(3.0D0*XI+1.0D0)
      FUN(4) = 1.0D0/16.0D0*(XI+1.0D0)*(9.0D0*XI*XI-1.0D0)
      DER(1) = -1.0D0/16.0D0*(-1.0D0-18.0D0*XI+27.0D0*XI*XI)
      DER(2) = 9.0D0/16.0D0*(9.0D0*XI*XI-2.0D0*XI-3.0D0)
      DER(3) = 9.0D0/16.0D0*(3.0D0-2.0D0*XI-9.0D0*XI*XI)
      DER(4) = 1.0D0/16.0D0*(27.0D0*XI*XI+18.0D0*XI-1.0D0)
      RETURN
      END