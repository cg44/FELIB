      SUBROUTINE QLVEC(D, ID, E, IE, Z, IZ, JZ, N, EPS, ITEST)
C                            THIS ROUTINE FINDS THE EIGENVALUES
C                            AND EIGENVECTORS OF A TRIDIAGONAL
C                            MATRIX, T, GIVEN WITH ITS DIAGONAL
C                            ELEMENTS IN THE ARRAY D(N) AND ITS
C                            SUB-DIAGONAL ELEMENTS IN THE LAST N -
C                            1 STORES OF THE ARRAY E(N), USING QL
C                            TRANSFORMATIONS. THE EIGENVALUES ARE
C                            OVERWRITTEN ON THE DIAGONAL ELEMENTS
C                            IN THE ARRAY D IN ASCENDING ORDER. THE
C                            EIGENVECTORS ARE FORMED IN THE ARRAY
C                            Z(N,N), OVERWRITING THE ACCUMULATED
C                            TRANSFORMATIONS AS SUPPLIED BY THE
C                            ROUTINE F01AJF. THE ROUTINE WILL
C                            FAIL IF ANY ONE EIGENVALUE TAKES MORE
C                            THAN 30 ITERATIONS.
      INTEGER I, I1, ID, IE, II, ITEST, IZ, J, JZ, K, L, M,
     *     M1, N, IERROR, ERRMES
      DOUBLE PRECISION B, C, D, E, EPS, F, G, H, P, R, S,
     *     Z, SRNAME
      DIMENSION D(ID), E(IE), Z(IZ,JZ)
      DATA SRNAME /8H QLVEC  /
      IF(ITEST.EQ.-1) GO TO 999
      IERROR=0
      IF(IZ.LT.N.OR.JZ.LT.N) IERROR=3
      IF(ID.LT.N.OR.IE.LT.N) IERROR=2
      IF(N.LE.0) IERROR=1
      ITEST=ERRMES(ITEST,IERROR,SRNAME)
      IF(ITEST.NE.0) RETURN
999   IF (N.EQ.1) GO TO 1020
      DO 1010 I=2,N
      E(I-1) = E(I)
 1010 CONTINUE
 1020 E(N) = 0.0D0
      B = 0.0D0
      F = 0.0D0
      DO 1150 L=1,N
      J = 0
      H = EPS*(DABS(D(L))+DABS(E(L)))
      IF (B.LT.H) B = H
C                            LOOK FOR SMALL SUB-DIAG ELEMENT
      DO 1030 M=L,N
      IF (DABS(E(M)).LE.B) GO TO 1040
 1030 CONTINUE
 1040 IF (M.EQ.L) GO TO 1140
 1050 IF (J.EQ.30) GO TO 1200
      J = J + 1
C                            FORM SHIFT
      G = D(L)
      H = D(L+1) - G
      IF (DABS(H).GE.DABS(E(L))) GO TO 1060
      P = H*0.5D0/E(L)
      R = DSQRT(P*P+1.0D0)
      H = P + R
      IF (P.LT.0.0D0) H = P - R
      D(L) = E(L)/H
      GO TO 1070
 1060 P = 2.0D0*E(L)/H
      R = DSQRT(P*P+1.0D0)
      D(L) = E(L)*P/(1.0D0+R)
 1070 H = G - D(L)
      I1 = L + 1
      IF (I1.GT.N) GO TO 1090
      DO 1080 I=I1,N
      D(I) = D(I) - H
 1080 CONTINUE
 1090 F = F + H
C                            QL TRANSFORMATION
      P = D(M)
      C = 1.0D0
      S = 0.0D0
      M1 = M - 1
      DO 1130 II=L,M1
      I = M1 - II + L
      G = C*E(I)
      H = C*P
      IF (DABS(P).LT.DABS(E(I))) GO TO 1100
      C = E(I)/P
      R = DSQRT(C*C+1.0D0)
      E(I+1) = S*P*R
      S = C/R
      C = 1.0D0/R
      GO TO 1110
 1100 C = P/E(I)
      R = DSQRT(C*C+1.0D0)
      E(I+1) = S*E(I)*R
      S = 1.0D0/R
      C = C/R
 1110 P = C*D(I) - S*G
      D(I+1) = H + S*(C*G+S*D(I))
C                            FORM VECTOR
      DO 1120 K=1,N
      H = Z(K,I+1)
      Z(K,I+1) = S*Z(K,I) + C*H
      Z(K,I) = C*Z(K,I) - S*H
 1120 CONTINUE
 1130 CONTINUE
      E(L) = S*P
      D(L) = C*P
      IF (DABS(E(L)).GT.B) GO TO 1050
 1140 D(L) = D(L) + F
 1150 CONTINUE
C                            ORDER EIGENVALUES AND EIGENVECTORS
      DO 1190 I=1,N
      K = I
      P = D(I)
      I1 = I + 1
      IF (I1.GT.N) GO TO 1170
      DO 1160 J=I1,N
      IF (D(J).GE.P) GO TO 1160
      K = J
      P = D(J)
 1160 CONTINUE
 1170 IF (K.EQ.I) GO TO 1190
      D(K) = D(I)
      D(I) = P
      DO 1180 J=1,N
      P = Z(J,I)
      Z(J,I) = Z(J,K)
      Z(J,K) = P
 1180 CONTINUE
 1190 CONTINUE
      ITEST = 0
      RETURN
 1200 IF(ITEST.EQ.-1) RETURN
      IERROR = 1
      ITEST=ERRMES(ITEST,IERROR,SRNAME)
      RETURN
      END
