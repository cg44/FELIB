C***********************************************************************
C$SPLIT$CHOBAK$*********************************************************
C***********************************************************************
      SUBROUTINE CHOBAK(A, IA, JA, R, IR, N, HBAND, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      PERFORMS BAKWARD SUBSTITUTION ON A MATRIX PROCESSED BY
C      CHORDN AND CHOFWD
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    10 OCT 1980 (KR)
C
C ARGUMENTS IN
C      A       ARRAY OF DIMENSION (IA,JA).  CONTAINS THE
C              ELEMENTS OF THE LOWER HALF OF THE POSITIVE
C              DEFINITE BAND MATRIX OF ORDER N AND WITH SEMI-
C              BANDWIDTH HBAND, REDUCED BY CHORDN OR CHOSOL
C      IA      FIRST DIMENSION OF A (.GE.N)
C      JA      SECOND DIMENSION OF A (.GE.HBAND)
C      R       ON ENTRY, CONTAINS ELEMENTS OF N RHS'S AFTER
C              PROCESSING BY CHOFWD
C      IR      DIMENSION OF VECTOR R (.GE.N)
C      N       ORDER OF MATRIX A
C      HBAND   SEMI-BANDWIDTH OF MATRIX A
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      R       ON EXIT, CONTAINS SOLUTION VECTOR
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE CHOBAK(A, IA, JA, R, IR, N, HBAND, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, HBAND, I, IA, IERROR, IJ, IR, ITEST,
     *     J, JA, L, M, N, W
      DOUBLE PRECISION A, R, SRNAME, X
      DIMENSION A(IA,JA), R(IR)
      DATA SRNAME /8H CHOBAK /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IR.LT.N) IERROR = 3
                        IF (IA.LT.N .OR. JA.LT.HBAND) IERROR = 2
                        IF (N.LE.0 .OR. HBAND.LE.0) IERROR = 1
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 W = HBAND - 1
      R(N) = R(N)/A(N,W+1)
      I = N - 1
 1020 X = 0.0D0
      L = I + W
      IF (I.GT.N-W) L = N
      M = I + 1
      DO 1030 J=M,L
      IJ = W + I - J + 1
      X = X + A(J,IJ)*R(J)
 1030 CONTINUE
      R(I) = (R(I)-X)/A(I,W+1)
      I = I - 1
      IF (I.NE.0) GO TO 1020
      RETURN
      END
