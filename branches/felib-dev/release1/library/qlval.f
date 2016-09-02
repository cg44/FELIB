      SUBROUTINE QLVAL(D, ID, E, IE, N, EPS, ITEST)
C                            THIS ROUTINE FINDS THE EIGENVALUES
C                            OF A TRIDIAGONAL MATRIX, T, GIVEN WITH
C                            ITS DIAGONAL ELEMENTS IN THE ARRAY
C                            D(N) AND ITS SUBDIAGONAL ELEMENTS IN
C                            THE LAST N - 1 STORES OF THE ARRAY
C                            E(N), USING QL TRANSFORMATIONS. THE
C                            EIGENVALUES ARE OVERWRITTEN ON THE
C                            DIAGONAL ELEMENTS IN THE ARRAY D IN
C                            ASCENDING ORDER. THE ROUTINE WILL
C                            FAIL IF ANY ONE EIGENVALUE TAKES MORE
C                            THAN 30 ITERATIONS.
      INTEGER I, I1, ID, IE, II, ITEST, J, L, M, M1, N,
     *     IERROR,ERRMES
      DOUBLE PRECISION B, C, D, E, EPS, F, G, H, P, R, S,
     *     SRNAME
      DIMENSION D(ID), E(IE)
      DATA SRNAME /8H QLVAL  /
      IF(ITEST.EQ.-1) GO TO 999
      IERROR=0
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
      DO 1170 L=1,N
      J = 0
      H = EPS*(DABS(D(L))+DABS(E(L)))
      IF (B.LT.H) B = H
C                            LOOK FOR SMALL SUB DIAGONAL ELEMENT
      DO 1030 M=L,N
      IF (DABS(E(M)).LE.B) GO TO 1040
 1030 CONTINUE
 1040 IF (M.EQ.L) GO TO 1130
 1050 IF (J.EQ.30) GO TO 1180
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
      DO 1120 II=L,M1
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
 1120 CONTINUE
      E(L) = S*P
      D(L) = C*P
      IF (DABS(E(L)).GT.B) GO TO 1050
 1130 P = D(L) + F
C                            ORDER EIGENVALUE
      IF (L.EQ.1) GO TO 1150
      DO 1140 II=2,L
      I = L - II + 2
      IF (P.GE.D(I-1)) GO TO 1160
      D(I) = D(I-1)
 1140 CONTINUE
 1150 I = 1
 1160 D(I) = P
 1170 CONTINUE
      ITEST = 0
      RETURN
 1180       IF(ITEST.EQ.-1) RETURN
      IERROR=3
      ITEST=ERRMES(ITEST,IERROR,SRNAME)
      RETURN
      END
