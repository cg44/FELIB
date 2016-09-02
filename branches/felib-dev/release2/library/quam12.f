C***********************************************************************
C$SPLIT$QUAM12$*********************************************************
C***********************************************************************
      SUBROUTINE QUAM12(FUN, IFUN, DER, IDER, JDER, XI, ETA, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS THE VALUES OF SHAPE FUNCTIONS AND THEIR
C      DERIVATIVES AT A SPECIFIED POINT FOR AN 12-NODED C0
C      CONTINUOUS QUADRILATERAL ELEMENT.  THE APPROXIMATED
C      FUNCTION WILL BE CONTINUOUS ACROSS ELEMENT BOUNDARIES
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    21 OCT 1980 (KR)
C
C ARGUMENTS IN
C      IFUN    LENGTH OF VECTOR FUN (.GE.12)
C      IDER    FIRST DIMENSION OF ARRAY DER (.GE.2)
C      JDER    SECOND DIMENSION OF ARRAY DER (.GE.12)
C      XI      FIRST LOCAL COORDINATE
C      ETA     SECOND LOCAL COORDINATE
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      FUN     VECTOR OF LENGTH IFUN.  FUN(I) CONTAINS THE
C              VALUE OF THE I'TH SHAPE FUNCTION AT (XI,ETA)
C              FOR I=1(1)12
C      DER     ARRAY OF DIMENSION (IDER,JDER).  DER(I,J)
C              CONTAINS THE VALUE OF THE DERIVATIVE OF THE J'TH
C              SHAPE FUNCTION WITH RESPECT TO THE I'TH
C              COORDINATE, FOR I=1(1)2 AND J=1(1)12
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE QUAM12(FUN, IFUN, DER, IDER, JDER, XI, ETA, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, IDER, IERROR, IFUN, ITEST, JDER
      DOUBLE PRECISION DER, ETA, FUN, SRNAME, XI, VAL, VEPS, DUMMY
      DIMENSION DER(IDER,JDER), FUN(IFUN)
      DATA SRNAME /8H QUAM12 /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IFUN.LT.12) IERROR = 1
                        IF (IDER.LT.2 .OR. JDER.LT.12) IERROR = 2
                        VAL = 1.0D0 + VEPS(DUMMY)
                        IF (DABS(XI).GT.VAL .OR. DABS(ETA).GT.VAL)
     *                      IERROR = 3
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 FUN(1) = 1.0D0/32.0D0*(1.0D0-XI)*(1.0D0-ETA)*(-10.0D0+9.0D0*
     *     (XI*XI+ETA*ETA))
      FUN(2) = 9.0D0/32.0D0*(1.0D0-XI)*(1.0D0-ETA*ETA)*
     *     (1.0D0-3.0D0*ETA)
      FUN(3) = 9.0D0/32.0D0*(1.0D0-XI)*(1.0D0-ETA*ETA)*
     *     (1.0D0+3.0D0*ETA)
      FUN(4) = 1.0D0/32.0D0*(1.0D0-XI)*(1.0D0+ETA)*(-10.0D0+9.0D0*
     *     (XI*XI+ETA*ETA))
      FUN(5) = 9.0D0/32.0D0*(1.0D0+ETA)*(1.0D0-XI*XI)*(1.0D0-3.0D0*
     *     XI)
      FUN(6) = 9.0D0/32.0D0*(1.0D0+ETA)*(1.0D0-XI*XI)*(1.0D0+3.0D0*
     *     XI)
      FUN(7) = 1.0D0/32.0D0*(1.0D0+XI)*(1.0D0+ETA)*(-10.0D0+9.0D0*
     *     (XI*XI+ETA*ETA))
      FUN(8) = 9.0D0/32.0D0*(1.0D0+XI)*(1.0D0-ETA*ETA)*
     *     (1.0D0+3.0D0*ETA)
      FUN(9) = 9.0D0/32.0D0*(1.0D0+XI)*(1.0D0-ETA*ETA)*
     *     (1.0D0-3.0D0*ETA)
      FUN(10) = 1.0D0/32.0D0*(1.0D0+XI)*(1.0D0-ETA)*(-10.0D0+9.0D0*
     *     (XI*XI+ETA*ETA))
      FUN(11) = 9.0D0/32.0D0*(1.0D0-ETA)*(1.0D0-XI*XI)*
     *     (1.0D0+3.0D0*XI)
      FUN(12) = 9.0D0/32.0D0*(1.0D0-ETA)*(1.0D0-XI*XI)*
     *     (1.0D0-3.0D0*XI)
      DER(1,1) = 1.0D0/32.0D0*(1.0D0-ETA)*(10.D0-27.D0*XI*XI+18.D0*
     *     XI-9.D0*ETA*ETA)
      DER(2,1) = 1.D0/32.D0*(1.D0-XI)*(10.D0-27.D0*ETA*ETA+18.D0*
     *     ETA-9.D0*XI*XI)
      DER(1,2) = -9.D0/32.D0*(1.D0-ETA*ETA)*(1.D0-3.D0*ETA)
      DER(2,2) = 9.D0/32.D0*(1.D0-XI)*(9.D0*ETA*ETA-2.D0*ETA-3.D0)
      DER(1,3) = -9.D0/32.D0*(1.D0-ETA*ETA)*(1.D0+3.D0*ETA)
      DER(2,3) = 9.D0/32.D0*(1.D0-XI)*(-9.D0*ETA*ETA-2.D0*ETA+3.D0)
      DER(1,4) = 1.D0/32.D0*(1.D0+ETA)*(10.D0-27.D0*XI*XI+18.D0*
     *     XI-9.D0*ETA*ETA)
      DER(2,4) = 1.D0/32.D0*(1.D0-XI)*(-10.D0+27.D0*ETA*ETA+18.D0*
     *     ETA+9.D0*XI*XI)
      DER(1,5) = 9.D0/32.D0*(1.D0+ETA)*(9.D0*XI*XI-2.D0*XI-3.D0)
      DER(2,5) = 9.D0/32.D0*(1.D0-XI*XI)*(1.D0-3.D0*XI)
      DER(1,6) = 9.D0/32.D0*(1.D0+ETA)*(3.D0-2.D0*XI-9.D0*XI*XI)
      DER(2,6) = 9.D0/32.D0*(1.D0-XI*XI)*(1.D0+3.D0*XI)
      DER(1,7) = 1.D0/32.D0*(1.D0+ETA)*(-10.D0+27.D0*XI*XI+18.D0*
     *     XI+9.D0*ETA*ETA)
      DER(2,7) = 1.D0/32.D0*(1.D0+XI)*(-10.D0+27.D0*ETA*ETA+18.D0*
     *     ETA+9.D0*XI*XI)
      DER(1,8) = 9.D0/32.D0*(1.D0-ETA*ETA)*(1.D0+3.D0*ETA)
      DER(2,8) = 9.D0/32.D0*(1.D0+XI)*(3.D0-2.D0*ETA-9.D0*ETA*ETA)
      DER(1,9) = 9.D0/32.D0*(1.D0-ETA*ETA)*(1.D0-3.D0*ETA)
      DER(2,9) = 9.D0/32.D0*(1.D0+XI)*(9.D0*ETA*ETA-2.D0*ETA-3.D0)
      DER(1,10) = 1.D0/32.D0*(1.D0-ETA)*(-10.D0+27.D0*XI*XI+18.D0*
     *     XI+9.D0*ETA*ETA)
      DER(2,10) = 1.D0/32.D0*(1.D0+XI)*(10.D0-27.D0*ETA*ETA+18.D0*
     *     ETA-9.D0*XI*XI)
      DER(1,11) = 9.D0/32.D0*(1.D0-ETA)*(3.D0-2.D0*XI-9.D0*XI*XI)
      DER(2,11) = -9.D0/32.D0*(1.D0-XI*XI)*(1.D0+3.D0*XI)
      DER(1,12) = 9.D0/32.D0*(1.D0-ETA)*(9.D0*XI*XI-2.D0*XI-3.D0)
      DER(2,12) = -9.D0/32.D0*(1.D0-XI*XI)*(1.D0-3.D0*XI)
      RETURN
      END
