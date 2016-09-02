C***********************************************************************
      SUBROUTINE CHOSOL(A, IA, JA, R, IR, N, HBAND, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      SOLVES A SET OF REAL SYMMETRIC POSITIVE DEFINITE BANDED
C      EQUATIONS WITH A SINGLE RIGHT HAND SIDE BY CHOLESKI
C      DECOMPOSITION.  ONLY THE LOWER BAND AND DIAGONAL ARE
C      STORED IN A RECTANGULAR ARRAY A
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 1.1  29 OCT 1979 (IMS)
C      COMMENTED    12 FEB 1980 (KR)
C
C ARGUMENTS IN
C      A       ON ENTRY CONTAINS LOWER HALF OF PD SYMMETRIC
C              BAND MATRIX STORED AS A RECTANGULAR ARRAY
C      IA      FIRST DIMENSION OF A (.GE. N)
C      JA      SECOND DIMENSION OF A (.GE. HBAND)
C      R       CONTAINS ELEMENTS OF RIGHT HAND SIDE
C      IR      DIMENSION OF R (.GE. N)
C      N       ORDER OF MATRIX A
C      HBAND   SEMI-BANDWIDTH OF A (INCLUDES DIAGONAL)
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      A       ON EXIT, CONTAINS LOWER TRIANGULAR REDUCED
C      R       MATRIX SOLUTION VECTOR
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE CHOSOL(A, IA, JA, R, IR, N, HBAND, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, HBAND, I, IA, IERROR, IJ, IK, IR, ITEST,
     *     J, JA, JTEST, K, L, LA, LB, LK, M, N, W
      DOUBLE PRECISION A, R, SRNAME, X
      DIMENSION A(IA,JA), R(IR)
      DATA SRNAME /8H CHOSOL /
C
C     PARAMETER CHECKING
C
      JTEST = ITEST
      IF (ITEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (IR.LT.N) IERROR = 3
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
C     RANGE CHECKING ON A(I,W+1)
C
         IF (JTEST.EQ.-1) GO TO 1030
         IERROR = 0
         IF ((A(I,W+1)-X).LE.0.0D0) IERROR = 4
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
      R(1) = R(1)/A(1,W+1)
      DO 1090 I=2,N
         X = 0.0D0
         K = 1
         IF (I.LE.W+1) K = W - I + 2
         DO 1080 J=K,W
            IJ = I + J - W - 1
            X = X + A(I,J)*R(IJ)
 1080    CONTINUE
         R(I) = (R(I)-X)/A(I,W+1)
 1090 CONTINUE
      R(N) = R(N)/A(N,W+1)
      I = N - 1
 1100 X = 0.0D0
      L = I + W
      IF (I.GT.N-W) L = N
      M = I + 1
      DO 1110 J=M,L
         IJ = W + I - J + 1
         X = X + A(J,IJ)*R(J)
 1110 CONTINUE
      R(I) = (R(I)-X)/A(I,W+1)
      I = I - 1
      IF (I.NE.0) GO TO 1100
C
      RETURN
      END
C***********************************************************************
