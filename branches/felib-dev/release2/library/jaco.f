C***********************************************************************
C$SPLIT$JACO$*********************************************************
C***********************************************************************
      SUBROUTINE JACO(A, IA, JA, DIAG, IDIAG, SUB, ISUB, N, HBAND,
     *     ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      REDUCES A REAL SYMMETRIC BAND MATRIX TO TRIDIAGONAL FORM
C      USING JACOBI ROTATIONS, FOR USE WITH QLVQL OR QLVEC.
C      THE LOWER TRIANGLE OF THE MATRIX IS STORED IN A
C      RECTANGULAR ARRAY.
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (IMS) --- SERC COPYRIGHT
C      COMMENTED    26 FEB 1980 (KR)
C
C ARGUMENTS IN
C      A       CONTAINS THE ELEMENTS OF THE LOWER TRIANGLE OF
C              THE POSITIVE DEFINITE BAND MATRIX
C      IA      FIRST DIMENSION OF A (.GE. N)
C      JA      SECOND DIMENSION OF A (.GE. HBAND)
C      IDIAG   DIMENSION OF VECTOR DIAG (.GE. N)
C      ISUB    DIMENSION OF VECTOR SUB (.GE. N)
C      N       ORDER OF MATRIX A
C      HBAND   SEMI-BANDWIDTH OF MATRIX A (INCLUDES DIAGONAL)
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      A       DESTROYED ON SUCCESFUL EXIT
C      DIAG    CONTAINS THE DIAGONAL ELEMENTS OF THE
C              TRIDIAGONAL MATRIX
C      SUB     CONTAINS THE N-1 OFF-DIAGONAL ELEMENTS OF THE
C              TRIDIAGONAL MATRIX STORED IN SUB(2) TO SUB(N),
C              WITH SUB(1)=0
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE JACO(A, IA, JA, DIAG, IDIAG, SUB, ISUB, N, HBAND,
C    *     ITEST)
C***********************************************************************
C
      INTEGER ERRMES, HBAND, I, IA, IDIAG, IERROR, IR, IRR,
     *     ISUB, ITEST, IUGL, J, J2, JA, JL, JM, K, KR, L,
     *     M, MAXL, MAXR, N, N2
      DOUBLE PRECISION A, B, C, C2, CS, DIAG, G, S, S2, SRNAME,
     *     SUB, U, U1
      DIMENSION A(IA,JA), DIAG(IDIAG), SUB(ISUB)
      DATA SRNAME /8H JACO   /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (ISUB.LT.N) IERROR = 4
                        IF (IDIAG.LT.N) IERROR = 3
                        IF (IA.LT.N .OR. JA.LT.HBAND) IERROR = 2
                        IF (N.LE.0) IERROR = 1
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 DO 1030 I=1,N
      DO 1020 J=1,HBAND
      L = HBAND + 1 - J
      K = I + 1 - L
      IF (K.LE.0) GO TO 1020
      A(K,L) = A(I,J)
 1020 CONTINUE
 1030 CONTINUE
      KR = HBAND - 1
      DO 1050 I=1,KR
      M = HBAND - I
      DO 1040 J=1,M
      K = N + 1 - I
      L = HBAND + 1 - J
      A(K,L) = 0.0
 1040 CONTINUE
 1050 CONTINUE
      M = HBAND - 1
      N2 = N - 2
      IF (N2.LT.1) GO TO 1140
      DO 1130 K=1,N2
      MAXR = M
      IF (N-K.LT.M) MAXR = N - K
      DO 1120 IRR=2,MAXR
      IR = 2 + MAXR - IRR
      KR = K + IR
      DO 1110 J=KR,N,M
      IF (J.EQ.KR) GO TO 1060
      IF (G.EQ.0.0D0) GO TO 1120
      JM = J - M
      B = -A(JM-1,M+1)/G
      IUGL = J - M
      GO TO 1070
 1060 IF (A(K,IR+1).EQ.0.0D0) GO TO 1120
      B = -A(K,IR)/A(K,IR+1)
      IUGL = K
 1070 S = 1.0D0/DSQRT(1.0D0+B*B)
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
      DO 1080 L=IUGL,J2
      JL = J - L
      U = C*A(L,JL) - S*A(L,JL+1)
      A(L,JL+1) = S*A(L,JL) + C*A(L,JL+1)
      A(L,JL) = U
 1080 CONTINUE
      JM = J - M
      IF (J.NE.KR) A(JM-1,M+1) = C*A(JM-1,M+1) - S*G
      MAXL = M - 1
      IF (N-J.LT.M-1) MAXL = N - J
      IF (MAXL.LE.0) GO TO 1100
      DO 1090 L=1,MAXL
      U = C*A(J-1,L+2) - S*A(J,L+1)
      A(J,L+1) = S*A(J-1,L+2) + C*A(J,L+1)
      A(J-1,L+2) = U
 1090 CONTINUE
 1100 IF (J+M.GT.N) GO TO 1110
      G = -S*A(J,M+1)
      A(J,M+1) = C*A(J,M+1)
 1110 CONTINUE
 1120 CONTINUE
 1130 CONTINUE
 1140 SUB(1) = 0.0D0
      DO 1150 I=1,N
      DIAG(I) = A(I,1)
 1150 CONTINUE
      IF (2.GT.N) GO TO 1170
      DO 1160 I=2,N
      SUB(I) = A(I-1,2)
 1160 CONTINUE
 1170 RETURN
      END
