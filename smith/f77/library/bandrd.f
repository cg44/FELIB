      SUBROUTINE BANDRD(N,IW,A,IA,D,E,E2)
C
C      THIS SUBROUTINE TRANSFORMS A REAL SYMMETRIC BAND MATRIX A,
C      OF ORDER N AND BAND WIDTH IW,
C      TO TRIDIAGONAL FORM BY AN APPROPRIATE
C      SEQUENCE OF JACOBI ROTATIONS. DURING THE TRANSFORMATION THE
C      PROPERTY OF THE BAND MATRIX IS MAINTAINED. THE METHOD YIELDS
C      A TRIDIAGONAL MATRIX, THE DIAGONAL ELEMENTS OF WHICH ARE IN
C      D(N) AND OFF-DIAGONAL ELEMENTS IN E(N).
C
      INTEGER M, IW, N2, N, K, MAXR, IRR, IR, KR, J, JM, IUGL, J2,
     + L, JL, MAXL, I, IA
      REAL G, B, S, C, C2, S2, CS, U, U1, A(IA,*), D(*), E(*),E2(*)
      N2 = N - 2
      IF (N2.LT.1) GO TO 180
      DO 160 K=1,N2
         MAXR = IW
         IF (N-K.LT.IW) MAXR = N - K
         DO 140 IRR=2,MAXR
            IR = 2 + MAXR - IRR
            KR = K + IR
            DO 120 J=KR,N,IW
               IF (J.EQ.KR) GO TO 20
               IF (G.EQ.0.0) GO TO 140
               JM = J - IW
               B = -A(JM-1,IW+1)/G
               IUGL = J - IW
               GO TO 40
   20          IF (A(K,IR+1).EQ.0.0) GO TO 140
               B = -A(K,IR)/A(K,IR+1)
               IUGL = K
   40          S = 1.0/SQRT(1.0+B*B)
               C = B*S
               C2 = C*C
               S2 = S*S
               CS = C*S
               U = C2*A(J-1,1) - 2.0*CS*A(J-1,2) + S2*A(J,1)
               U1 = S2*A(J-1,1) + 2.0*CS*A(J-1,2) + C2*A(J,1)
               A(J-1,2) = CS*(A(J-1,1)-A(J,1)) + (C2-S2)*A(J-1,2)
               A(J-1,1) = U
               A(J,1) = U1
               J2 = J - 2
               DO 60 L=IUGL,J2
                  JL = J - L
                  U = C*A(L,JL) - S*A(L,JL+1)
                  A(L,JL+1) = S*A(L,JL) + C*A(L,JL+1)
                  A(L,JL) = U
   60          CONTINUE
               JM = J - IW
               IF (J.NE.KR) A(JM-1,IW+1) = C*A(JM-1,IW+1) - S*G
               MAXL = IW - 1
               IF (N-J.LT.IW-1) MAXL = N - J
               IF (MAXL.LE.0) GO TO 100
               DO 80 L=1,MAXL
                  U = C*A(J-1,L+2) - S*A(J,L+1)
                  A(J,L+1) = S*A(J-1,L+2) + C*A(J,L+1)
                  A(J-1,L+2) = U
   80          CONTINUE
  100          IF (J+IW.GT.N) GO TO 120
               G = -S*A(J,IW+1)
               A(J,IW+1) = C*A(J,IW+1)
  120       CONTINUE
  140    CONTINUE
  160 CONTINUE
  180 E(1) = 0.0
      DO 200 I=1,N
         D(I) = A(I,1)
  200 CONTINUE
      IF (2.GT.N) GO TO 240
      DO 220 I=2,N
         E(I) = A(I-1,2)
  220 CONTINUE
  240 DO 260 I=1,N
         E2(I) = E(I)*E(I)
  260 CONTINUE
      RETURN
      END
