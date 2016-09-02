      SUBROUTINE MATINV(A, IA, JA, B, IB, JB, N, DET, ITEST)
      INTEGER ERRMES, IA, IB, IERROR, ITEST, JA, JB,
     *     K, L, M, N
      DOUBLE PRECISION A, B, DET, SRNAME, UNFLO, X
      DIMENSION A(IA,JA), B(IB,JB)
      DATA SRNAME /8H MATINV /
      IF (ITEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (N.GT.IB .OR. N.GT.JB) IERROR = 3
      IF (N.GT.IA .OR. N.GT.JA) IERROR = 2
      IF (N.LE.1.OR.N.GE.4) IERROR = 1
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF (IERROR.NE.0) RETURN
 1010 M = N - 1
      GO TO (1020, 1060), M
 1020 DET = A(1,1)*A(2,2) - A(1,2)*A(2,1)
      IF (ITEST.EQ.-1) GO TO 1030
      IERROR = 0
      IF (DABS(DET).LT.UNFLO(X)) IERROR = 4
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF (IERROR.NE.0) RETURN
 1030 B(1,1) = A(2,2)
      B(1,2) = -A(1,2)
      B(2,1) = -A(2,1)
      B(2,2) = A(1,1)
      DO 1050 K=1,2
      DO 1040 L=1,2
      B(K,L) = B(K,L)/DET
 1040 CONTINUE
 1050 CONTINUE
      RETURN
 1060 DET = A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      DET = DET - A(1,2)*(A(2,1)*A(3,3)-A(3,1)*A(2,3))
      DET = DET + A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      IF (ITEST.EQ.-1) GO TO 1070
      IERROR = 0
      IF (DABS(DET).LT.UNFLO(X)) IERROR = 4
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF (IERROR.NE.0) RETURN
 1070 B(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
      B(2,1) = -A(2,1)*A(3,3) + A(3,1)*A(2,3)
      B(3,1) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
      B(1,2) = -A(1,2)*A(3,3) + A(3,2)*A(1,3)
      B(2,2) = A(1,1)*A(3,3) - A(3,1)*A(1,3)
      B(3,2) = -A(1,1)*A(3,2) + A(3,1)*A(1,2)
      B(1,3) = A(1,2)*A(2,3) - A(2,2)*A(1,3)
      B(2,3) = -A(1,1)*A(2,3) + A(2,1)*A(1,3)
      B(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)
      DO 1090 K=1,3
      DO 1080 L=1,3
      B(K,L) = B(K,L)/DET
 1080 CONTINUE
 1090 CONTINUE
      RETURN
      END
