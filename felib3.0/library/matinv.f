C***********************************************************************
      SUBROUTINE MATINV(A, IA, JA, B, IB, JB, N, DET, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      FORMS THE INVERSE OF MATRIX A IN B
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 1.1  29 OCT 1979 (IMS)
C      COMMENTED    12 OCT 1980 (KR)
C
C ARGUMENTS IN
C      A       ARRAY OF DIMENSION (IA,JA) TO BE INVERTED
C      IA      FIRST DIMENSION OF A (.GE.N)
C      JA      SECOND DIMENSION OF A (.GE.N)
C      IB      FIRST DIMENSION OF ARRAY B (.GE.N)
C      JB      SECOND DIMENSION OF ARRAY B (.GE.N)
C      N       ORDER OF MATRIX A (.GT.1 .AND. .LT. 4)
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      B       ARRAY OF DIMENSION (IB,JB) CONTAINING INVERSE OF
C              A
C      DET     CONTAINS VALUE OF DETERMINANT OF MATRIX A
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE MATINV(A, IA, JA, B, IB, JB, N, DET, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, IA, IB, IERROR, ITEST, JA, JB, JTEST, K, L,
     *     M, N
      DOUBLE PRECISION A, B, DET, SRNAME, UNFLO, X
      DIMENSION A(IA,JA), B(IB,JB)
      DATA SRNAME /8H MATINV /
C
C     PARAMETER CHECKING
C
      JTEST = ITEST
      IF (JTEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (N.GT.IB .OR. N.GT.JB) IERROR = 3
      IF (N.GT.IA .OR. N.GT.JA) IERROR = 2
      IF (N.LE.1 .OR. N.GE.4) IERROR = 1
      ITEST = ERRMES(JTEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
C     MAIN BODY
C
 1010 M = N - 1
      GO TO (1020, 1060), M
C
C     CODE FOR 2X2 MATRIX
C
 1020 DET = A(1,1)*A(2,2) - A(1,2)*A(2,1)
C
C     CHECK THAT DETERMINANT IS NOT NEAR ZERO
C
      IF (JTEST.EQ.-1) GO TO 1030
      IERROR = 0
      IF (DABS(DET).LT.UNFLO(X)) IERROR = 4
      ITEST = ERRMES(JTEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
 1030 B(1,1) = A(2,2)
      B(1,2) = -A(1,2)
      B(2,1) = -A(2,1)
      B(2,2) = A(1,1)
      DO 1050 K=1,2
         DO 1040 L=1,2
            B(K,L) = B(K,L)/DET
 1040    CONTINUE
 1050 CONTINUE
      RETURN
C
C     CODE FOR 3X3 MATRIX
C
 1060 DET = A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      DET = DET - A(1,2)*(A(2,1)*A(3,3)-A(3,1)*A(2,3))
      DET = DET + A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
C
C     CHECK ON DETERMINANT NOT NEAR ZERO
C
      IF (JTEST.EQ.-1) GO TO 1070
      IERROR = 0
      IF (DABS(DET).LT.UNFLO(X)) IERROR = 4
      ITEST = ERRMES(JTEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
 1070 B(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
      B(2,1) = -A(2,1)*A(3,3) + A(3,1)*A(2,3)
      B(3,1) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
      B(1,2) = -A(1,2)*A(3,3) + A(3,2)*A(1,3)
      B(2,2) = A(1,1)*A(3,3) - A(3,1)*A(1,3)
      B(3,2) = -A(1,1)*A(3,2) + A(3,1)*A(1,2)
      B(1,3) = A(1,2)*A(2,3) - A(2,2)*A(1,3)
      B(2,3) = -A(1,1)*A(2,3) + A(2,1)*A(1,3)
      B(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)
C
      DO 1090 K=1,3
         DO 1080 L=1,3
            B(K,L) = B(K,L)/DET
 1080    CONTINUE
 1090 CONTINUE
      RETURN
C
      END
C***********************************************************************
