      SUBROUTINE TRIDIA(N,ATOL,A,D,E,Z,IZ)
C
C     THIS SUBROUTINE REDUCES THE GIVEN LOWER TRIANGLE OF A
C     SYMMETRIC MATRIX, A, STORED IN THE ARRAY A(N,N), TO
C     TRIDIAGONAL FORM USING HOUSEHOLDERS REDUCTION. THE DIAGONAL
C     OF THE RESULT IS STORED IN THE ARRAY D(N) AND THE
C     SUB-DIAGONAL IN THE LAST N - 1 STORES OF THE ARRAY E(N)
C     (WITH THE ADDITIONAL ELEMENT E(1) = 0). THE TRANSFORMATION
C     MATRICES ARE ACCUMULATED IN THE ARRAY Z(N,N). THE ARRAY
C     A IS LEFT UNALTERED UNLESS THE ACTUAL PARAMETERS
C     CORRESPONDING TO A AND Z ARE IDENTICAL.
C
      INTEGER N, I, L, IM1, IM2, J, K, IA, IZ
      REAL H, S, A(IZ,*), D(*), E(*), Z(IZ,*)
      DO 40 I=1,N
         DO 20 J=1,I
            Z(I,J) = A(I,J)
   20    CONTINUE
   40 CONTINUE
      IF (N.EQ.1) GO TO 280
      DO 260 II=2,N
         I = N - II + 2
         L = I - 2
         F = Z(I,I-1)
         G = 0.0
         IF (L.EQ.0) GO TO 80
         DO 60 K=1,L
            G = G + Z(I,K)*Z(I,K)
   60    CONTINUE
   80    H = G + F*F
C
C     IF G IS TOO SMALL FOR ORTHOGONALITY TO BE
C     GUARANTEED THE TRANSFORMATION IS SKIPPED
C
         IF (G.GT.ATOL) GO TO 100
         E(I) = F
         H = 0.0
         GO TO 240
  100    L = L + 1
         G = SQRT(H)
         IF (F.GE.0.0) G = -G
         E(I) = G
         H = H - F*G
         Z(I,I-1) = F - G
         F = 0.0
         DO 180 J=1,L
            Z(J,I) = Z(I,J)/H
            G = 0.0
C
C     FORM ELEMENT OF A*U
C
            DO 120 K=1,J
               G = G + Z(J,K)*Z(I,K)
  120       CONTINUE
            J1 = J + 1
            IF (J1.GT.L) GO TO 160
            DO 140 K=J1,L
               G = G + Z(K,J)*Z(I,K)
  140       CONTINUE
C
C     FORM ELEMENT OF P
C
  160       E(J) = G/H
            F = F + G*Z(J,I)
  180    CONTINUE
C
C     FORM K
C
         HH = F/(H+H)
C
C     FORM REDUCED A
C
         DO 220 J=1,L
            F = Z(I,J)
            G = E(J) - HH*F
            E(J) = G
            DO 200 K=1,J
               Z(J,K) = Z(J,K) - F*E(K) - G*Z(I,K)
  200       CONTINUE
  220    CONTINUE
  240    D(I) = H
  260 CONTINUE
  280 E(1) = 0.0
      D(1) = 0.0
C
C     ACCUMULATION OF TRANSFORMATION MATRICES
C
      DO 400 I=1,N
         L = I - 1
         IF (D(I).EQ.0.0) GO TO 360
         DO 340 J=1,L
            G = 0.0
            DO 300 K=1,L
               G = G + Z(I,K)*Z(K,J)
  300       CONTINUE
            DO 320 K=1,L
               Z(K,J) = Z(K,J) - G*Z(K,I)
  320       CONTINUE
  340    CONTINUE
  360    D(I) = Z(I,I)
         Z(I,I) = 1.0
         IF (L.EQ.0) GO TO 400
         DO 380 J=1,L
            Z(I,J) = 0.0
            Z(J,I) = 0.0
  380    CONTINUE
  400 CONTINUE
      RETURN
      END
