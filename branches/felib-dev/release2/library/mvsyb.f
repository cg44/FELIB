C***********************************************************************
C$SPLIT$MVSYB$*********************************************************
C***********************************************************************
      SUBROUTINE MVSYB(A, IA, JA, V, IV, W, IW, N, HBAND, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      POST-MULTIPLIES A REAL SYMMETRIC BANDED MATRIX STORED AS
C      A LOWER TRIANGLE BY A VECTOR
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    14 OCT 1980 (KR)
C
C ARGUMENTS IN
C      A       ARRAY OF DIMENSION (IA,JA).  CONTAINS THE
C              ELEMENTS OF THE LOWER HALF OF THE REAL SYMMETRIC
C              BAND MATRIX OF ORDER N AND SEMI-BANDWIDTH HBAND
C      IA      FIRST DIMENSION OF A (.GE.N)
C      JA      SECOND DIMENSION OF A (.GE.HBAND)
C      V       VECTOR OF DIMENSION IV
C      IV      DIMENSION OF VECTOR V (.GE.N)
C      IW      DIMENSION OF VECTOR W (.GE.N)
C      N       ORDER OF THE REAL SYMMETRIC BAND MATRIX
C      HBAND   SEMI-BANDWIDTH OF THE REAL SYMMETRIC BAND MATRIX
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      W       VECTOR OF DIMENSION IW.  CONTAINS THE RESULT OF
C              THE OPERATION W=A*V
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE MVSYB(A, IA, JA, V, IV, W, IW, N, HBAND, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, HBAND, I, IA, IERROR, IJ, ITEST, IV, IW,
     *     J, JA, JI, N                        
      DOUBLE PRECISION A, SRNAME, V, W, X
      DIMENSION A(IA,JA), V(IV), W(IW)
      DATA SRNAME /8H MVSYB  /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IW.LT.N) IERROR = 4
                        IF (IV.LT.N) IERROR = 3
                        IF (IA.LT.N .OR. JA.LT.HBAND) IERROR = 2
                        IF (N.LE.0) IERROR = 1
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 DO 1060 I=1,N
      X = 0.0D0
      J = HBAND
 1020 IF (I+J.LE.HBAND) GO TO 1030
      IJ = I + J - HBAND
      X = X + A(I,J)*V(IJ)
      J = J - 1
      IF (J.NE.0) GO TO 1020
 1030 J = HBAND - 1
 1040 IF (I-J.GE.N-HBAND+1) GO TO 1050
      JI = I - J + HBAND
      X = X + A(JI,J)*V(JI)
      J = J - 1
      IF (J.NE.0) GO TO 1040
 1050 W(I) = X
 1060 CONTINUE
      RETURN
      END
