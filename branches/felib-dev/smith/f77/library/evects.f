      SUBROUTINE EVECTS(N,ACHEPS,D,E,Z,IZ,IFAIL)
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS OF A
C     TRIDIAGONAL MATRIX, T, GIVEN WITH ITS DIAGONAL ELEMENTS IN
C     THE ARRAY D(N) AND ITS SUB-DIAGONAL ELEMENTS IN THE LAST N
C     - 1 STORES OF THE ARRAY E(N), USING QL TRANSFORMATIONS. THE
C     EIGENVALUES ARE OVERWRITTEN ON THE DIAGONAL ELEMENTS IN THE
C     ARRAY D IN ASCENDING ORDER. THE EIGENVECTORS ARE FORMED IN
C     THE ARRAY Z(N,N), OVERWRITING THE ACCUMULATED
C     TRANSFORMATIONS AS SUPPLIED BY THE SUBROUTINE F01AJF. THE
C     SUBROUTINE WILL FAIL IF ANY ONE EIGENVALUE TAKES MORE THAN 30
C     ITERATIONS.
C
      INTEGER P01AAF, ISAVE, IFAIL, N, I, L, J, M, I1, M1, II, K, IZ
      REAL B, F, H, ACHEPS, G, P, R, C, S, D(*), E(*), Z(IZ,*)
      ISAVE = IFAIL
      IF (N.EQ.1) GO TO 40
      DO 20 I=2,N
         E(I-1) = E(I)
   20 CONTINUE
   40 E(N) = 0.0
      B = 0.0
      F = 0.0
      DO 300 L=1,N
         J = 0
         H = ACHEPS*(ABS(D(L))+ABS(E(L)))
         IF (B.LT.H) B = H
C     LOOK FOR SMALL SUB-DIAG ELEMENT
         DO 60 M=L,N
            IF (ABS(E(M)).LE.B) GO TO 80
   60    CONTINUE
   80    IF (M.EQ.L) GO TO 280
  100    IF (J.EQ.30) GO TO 400
         J = J + 1
C     FORM SHIFT
         G = D(L)
         H = D(L+1) - G
         IF (ABS(H).GE.ABS(E(L))) GO TO 120
         P = H*0.5/E(L)
         R = SQRT(P*P+1.0)
         H = P + R
         IF (P.LT.0.0) H = P - R
         D(L) = E(L)/H
         GO TO 140
  120    P = 2.0*E(L)/H
         R = SQRT(P*P+1.0)
         D(L) = E(L)*P/(1.0+R)
  140    H = G - D(L)
         I1 = L + 1
         IF (I1.GT.N) GO TO 180
         DO 160 I=I1,N
            D(I) = D(I) - H
  160    CONTINUE
  180    F = F + H
C     QL TRANSFORMATION
         P = D(M)
         C = 1.0
         S = 0.0
         M1 = M - 1
         DO 260 II=L,M1
            I = M1 - II + L
            G = C*E(I)
            H = C*P
            IF (ABS(P).LT.ABS(E(I))) GO TO 200
            C = E(I)/P
            R = SQRT(C*C+1.0)
            E(I+1) = S*P*R
            S = C/R
            C = 1.0/R
            GO TO 220
  200       C = P/E(I)
            R = SQRT(C*C+1.0)
            E(I+1) = S*E(I)*R
            S = 1.0/R
            C = C/R
  220       P = C*D(I) - S*G
            D(I+1) = H + S*(C*G+S*D(I))
C     FORM VECTOR
            DO 240 K=1,N
               H = Z(K,I+1)
               Z(K,I+1) = S*Z(K,I) + C*H
               Z(K,I) = C*Z(K,I) - S*H
  240       CONTINUE
  260    CONTINUE
         E(L) = S*P
         D(L) = C*P
         IF (ABS(E(L)).GT.B) GO TO 100
  280    D(L) = D(L) + F
  300 CONTINUE
C     ORDER EIGENVALUES AND EIGENVECTORS
      DO 380 I=1,N
         K = I
         P = D(I)
         I1 = I + 1
         IF (I1.GT.N) GO TO 340
         DO 320 J=I1,N
            IF (D(J).GE.P) GO TO 320
            K = J
            P = D(J)
  320    CONTINUE
  340    IF (K.EQ.I) GO TO 380
         D(K) = D(I)
         D(I) = P
         DO 360 J=1,N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  360    CONTINUE
  380 CONTINUE
      IFAIL = 1
      RETURN
  400 IFAIL=0
      RETURN
      END
