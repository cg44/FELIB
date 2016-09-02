C***********************************************************************
C$SPLIT$TR6FN$*********************************************************
C***********************************************************************
      SUBROUTINE TR6FN(FUN, IFUN, DER, IDER, JDER, XI, ETA, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS THE VALUES OF THE SHAPE FUNCTIONS AND THEIR
C      DERIVATIVES AT A SPECIFIED POINT FOR A 6-NODED C0
C      CONTINUOUS TRIANGULAR ELEMENT.  APPROXIMATED FUNCTION
C      CONTINUOUS ACROSS ELEMENT BOUNDARIES
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    22 OCT 1980 (KR)
C
C ARGUMENTS IN
C      IFUN    LENGTH OF VECTOR FUN (.GE.6)
C      IDER    FIRST DIMENSION OF ARRAY DER (.GE.2)
C      JDER    SECOND DIMENSION OF ARRAY DER (.GE.6)
C      XI      FIRST LOCAL COORDINATE
C      ETA     SECOND LOCAL COORDINATE
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      FUN     VECTOR OF LENGTH IFUN.  FUN(I) CONTAINS THE
C              VALUE OF THE I'TH SHAPE FUNCTION AT THE POINT
C              (XI,ETA), FOR I=1(1)6
C      DER     ARRAY OF DIMENSION (IDER,JDER).  DER(I,J)
C              CONTAINS THE VALUE OF THE DERIVATIVE OF THE J'TH
C              SHAPE FUNCTION WITH RESPECT TO THE I'TH LOCAL
C              COORDINATE AT THE POINT (XI,ETA), FOR I=1(1)2
C              AND J=1(1)6
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE TR6FN(FUN, IFUN, DER, IDER, JDER, XI, ETA, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, IDER, IERROR, IFUN, ITEST, JDER
      DOUBLE PRECISION DER, ETA, FUN, L1, L2, L3, SRNAME, XI,
     *     YMIN, YMAX, XMIN, XMAX, VEPS, DUMMY 
      DIMENSION DER(IDER,JDER), FUN(IFUN)
      DATA SRNAME /8H TR6FN  /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IFUN.LT.6) IERROR = 1
                        IF (IDER.LT.2 .OR. JDER.LT.6) IERROR = 2
                        YMIN = 1.0D0/DSQRT(3.0D0)*(XI-1.0D0)-
     *                       VEPS(DUMMY)
                        YMAX = 1.0D0/DSQRT(3.0D0)*(1.0D0-XI)+
     *                       VEPS(DUMMY)
                        XMIN=-(0.5D0+VEPS(DUMMY))
                        XMAX=1.0D0+VEPS(DUMMY)
                        IF ((XI.LT.XMIN .OR. XI.GT.XMAX)
     *                      .OR. (ETA.LT.YMIN .OR.
     *                       ETA.GT.YMAX)) IERROR = 3
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 L1 = 1.0D0/3.0D0*(1.0D0+2.0D0*XI)
      L2 = 1.0D0/3.0D0*(1.0D0-XI-DSQRT(3.0D0)*ETA)
      L3 = 1.0D0/3.0D0*(1.0D0-XI+DSQRT(3.0D0)*ETA)
      FUN(1) = (2.0D0*L1-1.0D0)*L1
      FUN(2) = 4.0D0*L1*L2
      FUN(3) = (2.0D0*L2-1.0D0)*L2
      FUN(4) = 4.0D0*L2*L3
      FUN(5) = (2.0D0*L3-1.0D0)*L3
      FUN(6) = 4.0D0*L3*L1
      DER(1,1) = 2.0D0/3.0D0*(4.0D0*L1-1.0D0)
      DER(1,2) = 4.0D0/3.0D0*(2.0D0*L2-L1)
      DER(1,3) = -1.0D0/3.0D0*(4.0D0*L2-1.0D0)
      DER(1,4) = -4.0D0/3.0D0*(L2+L3)
      DER(1,5) = -1.0D0/3.0D0*(4.0D0*L3-1.0D0)
      DER(1,6) = 4.0D0/3.0D0*(2.0D0*L3-L1)
      DER(2,1) = 0.0D0
      DER(2,2) = -4.0D0/DSQRT(3.0D0)*L1
      DER(2,3) = -1.0D0/DSQRT(3.0D0)*(4.0D0*L2-1.0D0)
      DER(2,4) = 4.0D0/DSQRT(3.0D0)*(L2-L3)
      DER(2,5) = 1.0D0/DSQRT(3.0D0)*(4.0D0*L3-1.0D0)
      DER(2,6) = 4.0D0/DSQRT(3.0D0)*L1
      RETURN
      END
