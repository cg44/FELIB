C***********************************************************************
      SUBROUTINE ROD2(FUN, IFUN, DER, IDER, XI, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS THE VALUES OF THE SHAPE FUNCTIONS AND THEIR
C      DERIVATIVES AT A SPECIFIED POINT FOR A 2-NODED C0
C      CONTINUOUS LINE ELEMENT.  ONLY THE FUNCTION WILL BE
C      CONTINUOUS ACROSS ELEMENT BOUNDARIES.
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 1.1  29 OCT 1979 (CG)
C      COMMENTED    26 FEB 1980 (KR)
C
C ARGUMENTS IN
C      IFUN    DIMENSION OF THE VECTOR FUN (.GE. 2)
C      IDER    DIMENSION OF THE VECTOR DER (.GE. 2)
C      XI      SPECIFIES THE VALUE OF THE LOCAL COORDINATE AT
C              WHICH THE FUNCTION AND ITS DERIVATIVE ARE
C              REQUIRED
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      FUN     FUN(I) CONTAINS THE VALUE OF THE I'TH SHAPE
C              FUNCTION AT XI
C      DER     DER(I) CONTAINS THE VALUE OF THE DERIVATIVE OF
C              THE I'TH SHAPE FUNCTION AT XI
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE ROD2(FUN, IFUN, DER, IDER, XI, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, IDER, IERROR, IFUN, ITEST
      DOUBLE PRECISION DER, DUMMY, FUN, SRNAME, VEPS, XI, XMAX
      DIMENSION DER(IDER), FUN(IFUN)
      DATA SRNAME /8H ROD2   /
C
C     PARAMETER CHECKING
C
      IF (ITEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (IFUN.LT.2 .OR. IDER.LT.2) IERROR = 1
      XMAX = 1.0D0 + VEPS(DUMMY)
      IF (DABS(XI).GT.XMAX) IERROR = 2
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
C     MAIN BODY
C
 1010 FUN(1) = 0.5D0*(1.0D0-XI)
      FUN(2) = 0.5D0*(1.0D0+XI)
      DER(1) = -0.5D0
      DER(2) = 0.5D0
C
      RETURN
      END
C***********************************************************************
