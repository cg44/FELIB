C***********************************************************************
C$SPLIT$HOUSE$*********************************************************
C***********************************************************************
      SUBROUTINE HOUSE(A, IA, JA, T, IT, JT, DIAG, IDIAG, SUB,
     *     ISUB, N, TOL, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      USES HOUSEHOLDER'S METHOD TO REDUCE A REAL SYMMETRIC
C      MATRIX TO TRIDIAGONAL FORM FOR USE WITH QLVAL, QLVEC
C      A IS STORED AS A FULL MATRIX
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (IMS) --- SERC COPYRIGHT
C      COMMENTED    19 FEB 1980 (KR)
C
C ARGUMENTS IN
C      A       ARRAY CONTAINING THE ELEMENTS OF THE SYMMETRIC
C              MATRIX
C      IA      FIRST DIMENSION OF A (.GE. N)
C      JA      SECOND DIMENSION OF A (.GE. N)
C      IT      FIRST DIMENSION OF ARRAY T (.GE. N)
C      JT      SECOND DIMENSION OF T (.GE. N)
C      IDIAG   DIMENSION OF VECTOR DIAG (.GE. N)
C      ISUB    DIMENSION OF VECTOR SUB (.GE. N)
C      N       ORDER OF MATRIX A
C      TOL     VALUE OF RMIN/EPS, WHERE RMIN IS THE SMALLEST
C              POSITIVE NUMBER EXACTLY REPRESENTABLE ON THE
C              COMPUTER, AND EPS IS THE SMALLEST POSITIVE
C              NUMBER SUCH THAT 1.+EPS>1.
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      T       CONTAINS THE ORTHOGONAL MATRIX 'Q', THE PRODUCT
C              OF THE HOUSEHOLDER TRANSFORMATION MATRICES
C      DIAG    CONTAINS DIAGONAL ELEMENTS OF TRIDIAGONAL MATRIX
C      SUB     CONTAINS THE N-1 OFF-DIAGONAL ELEMENTS OF THE
C              TRIDIAGONAL MATRIX STORED IN SUB(2) TO SUB(N).
C              SUB(1)=0.
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE HOUSE(A, IA, JA, T, IT, JT, DIAG, IDIAG,
C    *     SUB, ISUB, N, TOL, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IA, IDIAG, IERROR, II, ISUB, IT, ITEST,
     *     J, JA, J1, JT, K, L, N
      DOUBLE PRECISION A, DIAG, F, G, H, HH, SRNAME, SUB, T,
     *     TOL
      DIMENSION A(IA,JA), DIAG(IDIAG), SUB(ISUB), T(IT,JT)
      DATA SRNAME /8H HOUSE  /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IDIAG.LT.N .OR. ISUB.LT.N) IERROR = 2
                        IF (IA.LT.N .OR. JA.LT.N) IERROR = 1
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 DO 1030 I=1,N
      DO 1020 J=1,I
      T(I,J) = A(I,J)
 1020 CONTINUE
 1030 CONTINUE
      IF (N.EQ.1) GO TO 1150
      DO 1140 II=2,N
      I = N - II + 2
      L = I - 2
      F = T(I,I-1)
      G = 0.0D0
      IF (L.EQ.0) GO TO 1050
      DO 1040 K=1,L
      G = G + T(I,K)*T(I,K)
C+++++
C     IF G IS TOO SMALL FOR ORTHOGONALITY TO BE GUARANTEED THE
C     TRANSFORMATION IS SKIPPED
C
 1040 CONTINUE
 1050 H = G + F*F
      IF (G.GT.TOL) GO TO 1060
      SUB(I) = F
      H = 0.0D0
      GO TO 1130
 1060 L = L + 1
      G = DSQRT(H)
      IF (F.GE.0.0D0) G = -G
      SUB(I) = G
      H = H - F*G
      T(I,I-1) = F - G
      F = 0.0D0
      DO 1100 J=1,L
C+++++
C     FORM ELEMENT OF A*U
C
      T(J,I) = T(I,J)/H
      G = 0.0D0
      DO 1070 K=1,J
      G = G + T(J,K)*T(I,K)
 1070 CONTINUE
      J1 = J + 1
      IF (J1.GT.L) GO TO 1090
      DO 1080 K=J1,L
C+++++
C     FORM ELEMENT OF P
C
      G = G + T(K,J)*T(I,K)
 1080 CONTINUE
 1090 SUB(J) = G/H
C+++++
C     FORM K
C
      F = F + G*T(J,I)
C+++++
C     FORM REDUCED A
C
 1100 CONTINUE
      HH = F/(H+H)
      DO 1120 J=1,L
      F = T(I,J)
      G = SUB(J) - HH*F
      SUB(J) = G
      DO 1110 K=1,J
      T(J,K) = T(J,K) - F*SUB(K) - G*T(I,K)
 1110 CONTINUE
 1120 CONTINUE
 1130 DIAG(I) = H
 1140 CONTINUE
C+++++
C     ACCUMULATION OF TRANSFORMATION MATRICES
C
 1150 SUB(1) = 0.0D0
      DIAG(1) = 0.0D0
      DO 1210 I=1,N
      L = I - 1
      IF (DIAG(I).EQ.0.0D0) GO TO 1190
      DO 1180 J=1,L
      G = 0.0D0
      DO 1160 K=1,L
      G = G + T(I,K)*T(K,J)
 1160 CONTINUE
      DO 1170 K=1,L
      T(K,J) = T(K,J) - G*T(K,I)
 1170 CONTINUE
 1180 CONTINUE
 1190 DIAG(I) = T(I,I)
      T(I,I) = 1.0D0
      IF (L.EQ.0) GO TO 1210
      DO 1200 J=1,L
      T(I,J) = 0.0D0
      T(J,I) = 0.0D0
 1200 CONTINUE
 1210 CONTINUE
      RETURN
      END
