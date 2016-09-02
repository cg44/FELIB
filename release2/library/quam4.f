C***********************************************************************
C$SPLIT$QUAM4$*********************************************************
C***********************************************************************
      SUBROUTINE QUAM4(FUN, IFUN, DER, IDER, JDER, XI, ETA, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS THE VALUES OF SHAPE FUNCTIONS AND THEIR
C      DERIVATIVES AT A SPECIFIED POINT FOR A 4-NODED C0
C      CONTINUOUS QUADRILATERAL ELEMENT.  THE APPROXIMATED
C      FUNCTION WILL BE CONTINUOUS ACROSS ELEMENT BOUNDARIES
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    21 OCT 1980 (KR)
C
C ARGUMENTS IN
C      IFUN    LENGTH OF VECTOR FUN (.GE.4)
C      IDER    FIRST DIMENSION OF ARRAY DER (.GE.2)
C      JDER    SECOND DIMENSION OF ARRAY DER (.GE.4)
C      XI      FIRST LOCAL COORDINATE VALUE
C      ETA     SECOND LOCAL COORDINATE VALUE
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      FUN     VECTOR OF LENGTH IFUN.  FUN(I) CONTAINS THE
C              VALUE OF THE I'TH SHAPE FUNCTION AT THE POINT
C              (XI,ETA)
C      DER     ARRAY OF DIMENSION (IDER,JDER).  DER(I,J)
C              CONTAINS THE VALUE OF THE DERIVATIVE OF THE J'TH
C              SHAPE FUNCTION WITH RESPECT TO THE I'TH
C              COORDINATE AT THE POINT (XI,ETA)
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE QUAM4(FUN, IFUN, DER, IDER, JDER, XI, ETA, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, IDER, IERROR, IFUN, ITEST, JDER
      DOUBLE PRECISION DER, ETA, ETAM, ETAP, FUN, SRNAME, XI,
     *     XIM, XIP, VAL, VEPS, DUMMY
      DIMENSION DER(IDER,JDER), FUN(IFUN)
      DATA SRNAME /8H QUAM4  /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IFUN.LT.4) IERROR = 1
                        IF (IDER.LT.2 .OR. JDER.LT.4) IERROR = 2
                        VAL = 1.0D0 + VEPS(DUMMY)
                        IF (DABS(XI).GT.VAL .OR. DABS(ETA).GT.VAL)
     *                      IERROR = 3
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 ETAM = 0.250D0*(1.0D0-ETA)
      ETAP = 0.250D0*(1.0D0+ETA)
      XIM = 0.250D0*(1.0D0-XI)
      XIP = 0.250D0*(1.0D0+XI)
      FUN(1) = 4.0D0*XIM*ETAM
      FUN(2) = 4.0D0*XIM*ETAP
      FUN(3) = 4.0D0*XIP*ETAP
      FUN(4) = 4.0D0*XIP*ETAM
      DER(1,1) = -ETAM
      DER(2,1) = -XIM
      DER(1,2) = -ETAP
      DER(2,2) = XIM
      DER(1,3) = ETAP
      DER(2,3) = XIP
      DER(1,4) = ETAM
      DER(2,4) = -XIP
      RETURN
      END
