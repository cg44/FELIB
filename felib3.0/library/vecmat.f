C***********************************************************************
      SUBROUTINE VECMAT(V, IV, A, IA, JA, M, N, W, IW, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      PRE-MULTIPLIES THE MATRIX A BY THE VECTOR V, STORING
C      THE RESULT IN THE VECTOR W
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 2.0  1 FEB 1981 (CG)
C      COMMENTED    1 FEB 1981 (CG)
C
C ARGUMENTS IN
C      V       VECTOR OF DIMENSION IV
C      IV      DIMENSION OF V (.GE.M)
C      A       ARRAY OF DIMENSION (IA,JA)
C      IA      FIRST DIMENSION OF A (.GE.M)
C      JA      SECOND DIMENSION OF A (.GE.N)
C      M       NUMBER OF ROWS OF A TO BE USED IN THE
C              MULTIPLICATION
C      N       NUMBER OF COLUMNS OF A AND THE NUMBER OF
C              ELEMENETS OF V TO BE USED IN THE MULTIPLICATION
C      IW      DIMENSION OF VECTOR W (.GE.N)
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      W       VECTOR OF DIMENSION IW; CONTAINS THE RESULT OF
C              THE OPERATION W=V*A
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE VECMAT(A, IA, JA, V, IV, M, N, W, IW, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IA, IERROR, ITEST, IV, IW, J, JA,
     *     M, N
      DOUBLE PRECISION A, SRNAME, V, W, X
      DIMENSION A(IA,JA), V(IV), W(IW)
C
      DATA SRNAME /8H VECMAT /
C
C     PARAMETER CHECKING
C
      IF (ITEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (N.GT.IW) IERROR = 4
      IF (M.GT.IV) IERROR = 3
      IF (M.GT.IA .OR. N.GT.JA) IERROR = 2
      IF (M.LE.0 .OR. N.LE.0) IERROR = 1
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
C     MAIN BODY
C
 1010 DO 1030 I=1,N
         X = 0.0D0
         DO 1020 J=1,M
            X = X + A(J,I)*V(J)
 1020    CONTINUE
         W(I) = X
 1030 CONTINUE
C
      RETURN
      END
C***********************************************************************
