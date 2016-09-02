      SUBROUTINE CHOSOL(KB, IKB, JKB, LOADS, ILOADS, N, HBAND, ITEST)
      INTEGER A, B, I, IK, IKB, ITEST, J, JKB, K, L, W, M, IJ,
     *     LK, N, HBAND, IERROR, ERRMES, ILOADS
      DOUBLE PRECISION KB, X, SRNAME, LOADS
      DIMENSION KB(IKB,JKB), LOADS(ILOADS)
      DATA SRNAME /8H CHOSOL /
      IF(ITEST.EQ.-1) GO TO 999
      IERROR=0
      IF(ILOADS.LT.N) IERROR=3
      IF(IKB.LT.N.OR.JKB.LT.HBAND) IERROR=2
      IF(N.LE.0.OR.HBAND.LE.0) IERROR=1
      ITEST=ERRMES(ITEST,IERROR,SRNAME)
      IF(ITEST.NE.0) RETURN
999   W = HBAND - 1
      DO 1050 I=1,N
      X = 0.0D0
      DO 1010 J=1,W
      X = X + KB(I,J)*KB(I,J)
 1010 CONTINUE
      IF(ITEST.EQ.-1) GO TO 998
      IERROR=0
      IF((KB(I,W+1)-X).LE.0.0D0) IERROR=4
      ITEST=ERRMES(ITEST,IERROR,SRNAME)
      IF(ITEST.NE.0) RETURN
998   KB(I,W+1) = DSQRT(KB(I,W+1)-X)
      DO 1040 K=1,W
      X = 0.0D0
      IF (I+K.GT.N) GO TO 1040
      IF (K.EQ.W) GO TO 1030
      L = W - K
 1020 IK = I + K
      LK = L + K
      X = X + KB(IK,L)*KB(I,LK)
      L = L - 1
      IF (L.NE.0) GO TO 1020
 1030 A = I + K
      B = W - K + 1
      KB(A,B) = (KB(A,B)-X)/KB(I,W+1)
 1040 CONTINUE
 1050 CONTINUE
      LOADS(1) = LOADS(1)/KB(1,W+1)
      DO 2020 I=2,N
      X = 0.0D0
      K = 1
      IF (I.LE.W+1) K = W - I + 2
      DO 2010 J=K,W
      IJ = I + J - W - 1
      X = X + KB(I,J)*LOADS(IJ)
 2010 CONTINUE
      LOADS(I) = (LOADS(I)-X)/KB(I,W+1)
 2020 CONTINUE
      LOADS(N) = LOADS(N)/KB(N,W+1)
      I = N - 1
 3010 X = 0.0D0
      L = I + W
      IF (I.GT.N-W) L = N
      M = I + 1
      DO 3020 J=M,L
      IJ = W + I - J + 1
      X = X + KB(J,IJ)*LOADS(J)
 3020 CONTINUE
      LOADS(I) = (LOADS(I)-X)/KB(I,W+1)
      I = I - 1
      IF (I.NE.0) GO TO 3010
      RETURN
      END
