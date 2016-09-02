C***********************************************************************
      SUBROUTINE CHORDN(A, IA, JA, N, HBAND, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      PERFORMS CHOLESKI REDUCTION ON A REAL SYMMETRIC POSITIVE
C      DEFINITE BANDED MATRIX.  ONLY THE LOWER BAND AND
C      DIAGONAL ARE STORED IN A RECTANGULAR ARRAY A
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 1.1    29 OCT 1979 (IMS)
C      COMMENTED      12 FEB 1980 (CG)
C      CHECKING ADDED 12 FEB 1980 (KR)
C
C ARGUMENTS IN
C      A       ON ENTRY, CONTAINS LOWER BAND AND DIAGONAL OF
C              PD REAL SYMMETRIC MATRIX
C      IA      FIRST DIMENSION OF A (.GE. N)
C      JA      SECOND DIMENSION OF A (.GE. HBAND)
C      N       ORDER OF MATRIX A
C      HBAND   SEMI-BANDWIDTH OF A (INCLUDING DIAGONAL)
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      A       REDUCED MATRIX L, WHERE THE INPUT MATRIX A HAS
C              BEEN REDUCED TO TRIANGULAR MATRICES L AND LT
C              WHERE A=L LT
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE CHORDN(A, IA, JA, N, HBAND, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, HBAND, I, IA, IERROR, IK, ITEST, J, JA,
     *     JTEST, K, L, LA, LB, LK, N, W
      DOUBLE PRECISION A, SRNAME, X
      DIMENSION A(IA,JA)
      DATA SRNAME /8H CHORDN /
C
C     PARAMETER CHECKING
C
      JTEST = ITEST
      IF (JTEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (IA.LT.N .OR. JA.LT.HBAND) IERROR = 2
      IF (N.LE.0 .OR. HBAND.LE.0) IERROR = 1
      ITEST = ERRMES(JTEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
C     MAIN BODY
C
 1010 W = HBAND - 1
      DO 1070 I=1,N
         X = 0.0D0
         DO 1020 J=1,W
            X = X + A(I,J)*A(I,J)
 1020    CONTINUE
C
C     CHECK ON VALUE OF A(I,J)
C
         IF (JTEST.EQ.-1) GO TO 1030
         IERROR = 0
         IF ((A(I,W+1)-X).LE.0.0D0) IERROR = 3
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
C
 1030    A(I,W+1) = DSQRT(A(I,W+1)-X)
C
         DO 1060 K=1,W
            X = 0.0D0
            IF (I+K.GT.N) GO TO 1060
            IF (K.EQ.W) GO TO 1050
            L = W - K
 1040       IK = I + K
            LK = L + K
            X = X + A(IK,L)*A(I,LK)
            L = L - 1
            IF (L.NE.0) GO TO 1040
 1050       LA = I + K
            LB = W - K + 1
            A(LA,LB) = (A(LA,LB)-X)/A(I,W+1)
 1060    CONTINUE
 1070 CONTINUE
C
      RETURN
      END
C***********************************************************************
