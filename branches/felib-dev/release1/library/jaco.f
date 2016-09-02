      SUBROUTINE JACO(A, IA, JA, D, ID, E, IE, N, HBAND,
     *     ITEST)
C                            TRANSFORMS A REAL SYMMETRIC BAND
C                            MATRIX A, OF ORDER N AND BAND WIDTH M,
C                            TO TRIDIAGONAL FORM BY AN APPROPRIATE
C                            SEQUENCE OF JACOBI ROTATIONS. DURING
C                            THE TRANSFORMATION THE PROPERTY OF THE
C                            BAND MATRIX IS MAINTAINED. THE METHOD
C                            YIELDS A TRIDIAGONAL MATRIX, THE
C                            DIAGONAL ELEMENTS OF WHICH ARE IN D(N)
C                            AND OFF-DIAGONAL ELEMENTS IN E(N).
      INTEGER I, IA, ID, IE, IR, IRR, ITEST, IUGL, J,
     *     J2, JA, JL, JM, K, KR, L, M, MAXL, MAXR, N,
     *     N2,IERROR,ERRMES,HBAND
      DOUBLE PRECISION A, B, C, C2, CS, D, E, G, S, S2,
     *     U, U1, SRNAME, X
      DIMENSION A(IA,JA), D(ID), E(IE)
      DATA SRNAME /8H JACO   /
      IF(ITEST.EQ.-1) GO TO 999
      IERROR=0
      IF(IE.LT.N) IERROR=4
      IF(ID.LT.N) IERROR=3
      IF(IA.LT.N.OR.JA.LT.HBAND) IERROR=2
      IF(N.LE.0) IERROR=1
      ITEST=ERRMES(ITEST,IERROR,SRNAME)
      IF(ITEST.NE.0) RETURN
999   K=1
      DO 901 J=2,HBAND
      K=K+1
      DO 902 I=K,N
      L=N-I+1
      M=L+K-1
      A(M,J)=A(L,J)
      A(L,J)=0.D0
902   CONTINUE
901   CONTINUE
      M=HBAND/2
      DO 904 J=1,M
      DO 903 I=1,N
      X=A(I,J)
      K=HBAND-J+1
      A(I,J)=A(I,K)
      A(I,K)=X
903   CONTINUE
904   CONTINUE
      M = HBAND - 1
      N2 = N - 2
      IF (N2.LT.1) GO TO 1090
      DO 1080 K=1,N2
      MAXR = M
      IF (N-K.LT.M) MAXR = N - K
      DO 1070 IRR=2,MAXR
      IR = 2 + MAXR - IRR
      KR = K + IR
      DO 1060 J=KR,N,M
      IF (J.EQ.KR) GO TO 1010
      IF (G.EQ.0.0D0) GO TO 1070
      JM = J - M
      B = -A(JM-1,M+1)/G
      IUGL = J - M
      GO TO 1020
 1010 IF (A(K,IR+1).EQ.0.0D0) GO TO 1070
      B = -A(K,IR)/A(K,IR+1)
      IUGL = K
 1020 S = 1.0D0/DSQRT(1.0D0+B*B)
      C = B*S
      C2 = C*C
      S2 = S*S
      CS = C*S
      U = C2*A(J-1,1) - 2.0D0*CS*A(J-1,2) + S2*A(J,1)
      U1 = S2*A(J-1,1) + 2.0D0*CS*A(J-1,2) + C2*A(J,1)
      A(J-1,2) = CS*(A(J-1,1)-A(J,1)) + (C2-S2)*A(J-1,2)
      A(J-1,1) = U
      A(J,1) = U1
      J2 = J - 2
      DO 1030 L=IUGL,J2
      JL = J - L
      U = C*A(L,JL) - S*A(L,JL+1)
      A(L,JL+1) = S*A(L,JL) + C*A(L,JL+1)
      A(L,JL) = U
 1030 CONTINUE
      JM = J - M
      IF (J.NE.KR) A(JM-1,M+1) = C*A(JM-1,M+1) - S*G
      MAXL = M - 1
      IF (N-J.LT.M-1) MAXL = N - J
      IF (MAXL.LE.0) GO TO 1050
      DO 1040 L=1,MAXL
      U = C*A(J-1,L+2) - S*A(J,L+1)
      A(J,L+1) = S*A(J-1,L+2) + C*A(J,L+1)
      A(J-1,L+2) = U
 1040 CONTINUE
 1050 IF (J+M.GT.N) GO TO 1060
      G = -S*A(J,M+1)
      A(J,M+1) = C*A(J,M+1)
 1060 CONTINUE
 1070 CONTINUE
 1080 CONTINUE
 1090 E(1) = 0.0D0
      DO 1100 I=1,N
      D(I) = A(I,1)
 1100 CONTINUE
      IF (2.GT.N) GO TO 1120
      DO 1110 I=2,N
      E(I) = A(I-1,2)
 1110 CONTINUE
 1120 RETURN
      END
