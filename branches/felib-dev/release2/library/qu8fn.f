C***********************************************************************
C$SPLIT$QU8FN$*********************************************************
C***********************************************************************
      SUBROUTINE QU8FN(FUN, IFUN, DER, IDER, JDER, XI, ETA, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS THE VALUES OF SHAPE FUNCTIONS AND THEIR
C      DERIVATIVES AT A SPECIFIED POINT FOR AN 8-NODED C0
C      CONTINUOUS QUADRILATERAL ELEMENT.  THE APPROXIMATED
C      FUNCTION WILL BE CONTINUOUS ACROSS ELEMENT BOUNDARIES
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    21 OCT 1980 (KR)
C
C ARGUMENTS IN
C      IFUN    LENGTH OF VECTOR FUN (.GE.8)
C      IDER    FIRST DIMENSION OF ARRAY DER (.GE.2)
C      JDER    SECOND DIMENSION OF ARRAY DER (.GE.8)
C      XI      FIRST LOCAL COORDINATE
C      ETA     SECOND LOCAL COORDINATE
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      FUN     VECTOR OF LENGTH IFUN.  FUN(I) CONTAINS THE
C              VALUE OF THE I'TH SHAPE FUNCTION AT (XI,ETA)
C              FOR I=1(1)8
C      DER     ARRAY OF DIMENSION (IDER,JDER).  DER(I,J)
C              CONTAINS THE VALUE OF THE DERIVATIVE OF THE J'TH
C              SHAPE FUNCTION WITH RESPECT TO THE I'TH
C              COORDINATE, FOR I=1(1)2 AND J=1(1)8
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE QU8FN(FUN, IFUN, DER, IDER, JDER, XI, ETA, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, IDER, IERROR, IFUN, ITEST, JDER
      DOUBLE PRECISION DER, ETA, ETAM, ETAP, FUN, SRNAME, XI,
     *     XIM, XIP, VAL, VEPS, DUMMY
      DIMENSION DER(IDER,JDER), FUN(IFUN)
      DATA SRNAME /8H QU8FN  /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IFUN.LT.8) IERROR = 1
                        IF (IDER.LT.2 .OR. JDER.LT.8) IERROR = 2
                        VAL = 1.0D0 + VEPS(DUMMY)
                        IF (DABS(XI).GT.VAL .OR. DABS(ETA).GT.VAL)
     *                      IERROR = 3
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 ETAM = 0.250D0*(1.0D0-ETA)
      ETAP = 0.250D0*(1.0D0+ETA)
      XIM = 0.250D0*(1.0D0-XI)
      XIP = 0.250D0*(1.0D0+XI)
      DER(1,1) = ETAM*(2.0D0*XI+ETA)
      DER(1,2) = -8.0D0*ETAM*ETAP
      DER(1,3) = ETAP*(2.0D0*XI-ETA)
      DER(1,4) = -4.0D0*ETAP*XI
      DER(1,5) = ETAP*(2.0D0*XI+ETA)
      DER(1,6) = 8.0D0*ETAP*ETAM
      DER(1,7) = ETAM*(2.0D0*XI-ETA)
      DER(1,8) = -4.0D0*ETAM*XI
      DER(2,1) = XIM*(XI+2.0D0*ETA)
      DER(2,2) = -4.0D0*XIM*ETA
      DER(2,3) = XIM*(2.0D0*ETA-XI)
      DER(2,4) = 8.0D0*XIM*XIP
      DER(2,5) = XIP*(XI+2.0D0*ETA)
      DER(2,6) = -4.0D0*XIP*ETA
      DER(2,7) = XIP*(2.0D0*ETA-XI)
      DER(2,8) = -8.0D0*XIM*XIP
      FUN(1) = 4.0D0*ETAM*XIM*(-XI-ETA-1.0D0)
      FUN(2) = 32.0D0*XIM*ETAM*ETAP
      FUN(3) = 4.0D0*ETAP*XIM*(-XI+ETA-1.0D0)
      FUN(4) = 32.0D0*XIM*XIP*ETAP
      FUN(5) = 4.0D0*XIP*ETAP*(XI+ETA-1.0D0)
      FUN(6) = 32.0D0*XIP*ETAP*ETAM
      FUN(7) = 4.0D0*XIP*ETAM*(XI-ETA-1.0D0)
      FUN(8) = 32.0D0*XIM*XIP*ETAM
      RETURN
      END
