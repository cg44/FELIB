      SUBROUTINE HOUSE(A, IA, JA, Z, IZ, JZ, D, ID, E, IE,
     *     N, TOL, ITEST)
C                            THIS ROUTINE REDUCES THE GIVEN
C                            LOWER TRIANGLE OF A SYMMETRIC MATRIX,
C                            A, STORED IN THE ARRAY A(N,N), TO
C                            TRIDIAGONAL FORM USING HOUSEHOLDERS
C                            REDUCTION. THE DIAGONAL OF THE RESULT
C                            IS STORED IN THE ARRAY D(N) AND THE
C                            SUB-DIAGONAL IN THE LAST N - 1 STORES
C                            OF THE ARRAY E(N) (WITH THE ADDITIONAL
C                            ELEMENT E(1) = 0). THE TRANSFORMATION
C                            MATRICES ARE ACCUMULATED IN THE ARRAY
C                            Z(N,N). THE ARRAY A IS LEFT UNALTERED
C                            UNLESS THE ACTUAL PARAMETERS
C                            CORRESPONDING TO A AND Z ARE
C                            IDENTICAL.
      INTEGER I, IA, ID, IE, II, ITEST, IZ, J, J1, JZ, K, L,
     *     N, IERROR, ERRMES
      DOUBLE PRECISION A, D, E, F, G, H, HH, TOL, Z, SRNAME
      DIMENSION A(IA,JA), D(ID), E(IE), Z(IZ,JZ)
      DATA SRNAME /8H HOUSE  /
      IF(ITEST.EQ.-1) GO TO 999
      IERROR=0
      IF(ID.LT.N.OR.IE.LT.N) IERROR=2
      IF(IA.LT.N.OR.JA.LT.N) IERROR=1
      ITEST=ERRMES(ITEST,IERROR,SRNAME)
      IF(ITEST.NE.0) RETURN
999   DO 1020 I=1,N
      DO 1010 J=1,I
      Z(I,J) = A(I,J)
 1010 CONTINUE
 1020 CONTINUE
      IF (N.EQ.1) GO TO 1140
      DO 1130 II=2,N
      I = N - II + 2
      L = I - 2
      F = Z(I,I-1)
      G = 0.0D0
      IF (L.EQ.0) GO TO 1040
      DO 1030 K=1,L
      G = G + Z(I,K)*Z(I,K)
 1030 CONTINUE
 1040 H = G + F*F
C                            IF G IS TOO SMALL FOR ORTHOGONALITY TO
C                            BE GUARANTEED THE TRANSFORMATION IS
C                            SKIPPED
      IF (G.GT.TOL) GO TO 1050
      E(I) = F
      H = 0.0D0
      GO TO 1120
 1050 L = L + 1
      G = DSQRT(H)
      IF (F.GE.0.0D0) G = -G
      E(I) = G
      H = H - F*G
      Z(I,I-1) = F - G
      F = 0.0D0
      DO 1090 J=1,L
      Z(J,I) = Z(I,J)/H
      G = 0.0D0
C                            FORM ELEMENT OF A*U
      DO 1060 K=1,J
      G = G + Z(J,K)*Z(I,K)
 1060 CONTINUE
      J1 = J + 1
      IF (J1.GT.L) GO TO 1080
      DO 1070 K=J1,L
      G = G + Z(K,J)*Z(I,K)
 1070 CONTINUE
C                            FORM ELEMENT OF P
 1080 E(J) = G/H
      F = F + G*Z(J,I)
 1090 CONTINUE
C                            FORM K
      HH = F/(H+H)
C                            FORM REDUCED A
      DO 1110 J=1,L
      F = Z(I,J)
      G = E(J) - HH*F
      E(J) = G
      DO 1100 K=1,J
      Z(J,K) = Z(J,K) - F*E(K) - G*Z(I,K)
 1100 CONTINUE
 1110 CONTINUE
 1120 D(I) = H
 1130 CONTINUE
 1140 E(1) = 0.0D0
      D(1) = 0.0D0
C                            ACCUMULATION OF TRANSFORMATION
C                            MATRICES
      DO 1200 I=1,N
      L = I - 1
      IF (D(I).EQ.0.0D0) GO TO 1180
      DO 1170 J=1,L
      G = 0.0D0
      DO 1150 K=1,L
      G = G + Z(I,K)*Z(K,J)
 1150 CONTINUE
      DO 1160 K=1,L
      Z(K,J) = Z(K,J) - G*Z(K,I)
 1160 CONTINUE
 1170 CONTINUE
 1180 D(I) = Z(I,I)
      Z(I,I) = 1.0D0
      IF (L.EQ.0) GO TO 1200
      DO 1190 J=1,L
      Z(I,J) = 0.0D0
      Z(J,I) = 0.0D0
 1190 CONTINUE
 1200 CONTINUE
      RETURN
      END
