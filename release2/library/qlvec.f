C***********************************************************************
C$SPLIT$QLVEC$*********************************************************
C***********************************************************************
      SUBROUTINE QLVEC(DIAG, IDIAG, SUB, ISUB, T, IT, JT, N, EPS,
     *     ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CALCULATES ALL THE EIGENVALUES AND EIGENVECTORS OF A
C      REAL SYMMETRIC TRIDIAGONAL MATRIX, OR OF A FULL REAL
C      SYMMETRIC MATRIX REDUCED TO TRIDIAGONAL FORM BY
C      SUBROUTINE HOUSE
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (IMS) --- SERC COPYRIGHT
C      COMMENTED    16 OCT 1980 (KR)
C
C ARGUMENTS IN
C      DIAG    VECTOR OF LENGTH IDIAG.  ON ENTRY, CONTAINS
C              DIAGONAL ELEMENTS OF TRIDIAGONAL MATRIX
C      IDIAG   DIMENSION OF VECTOR DIAG (.GE.N)
C      SUB     VECTOR OF LENGTH ISUB.  ON ENTRY, CONTAINS
C              SUB-DIAGONAL ELEMENTS OF TRIDIAGONAL MATRIX IN
C              SUB(2) TO SUB(N).  ON EXIT, CONTENTS OF SUB ARE
C              DESTROYED
C      ISUB    DIMENSION OF VECTOR SUB (.GE.N)
C      T       ARRAY OF DIMENSION (IT,JT).  THERE ARE TWO MODES
C              OF OPERATION:
C                (I)  EIGEN VECTORS OF TRIDIAGONAL MATRIX
C                     ON ENTRY T SHOULD CONTAIN THE IDENTITY
C                     MATRIX OF ORDER N
C                (II) EIGENVECTORS OF FULL SYMMETRIC MATRIX
C                     ON ENTRY T SHOULD CONTAIN THE ORTHOGONAL
C                     MATRIX OBTAINED FROM SUBROUTINE HOUSE
C      IT      FIRST DIMENSION OF ARRAY T (.GE.N)
C      JT      SECOND DIMENSION OF T (.GE.N)
C      N       ORDER OF TRIDIAGONAL MATRIX
C      EPS     SMALLEST POSITIVE NUMBER SUCH THAT 1.+EPS>1.
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      DIAG    VECTOR OF DIMENSION IDIAG.  ON EXIT, CONTAINS
C              THE EIGENVALUES IN ASCENDING ORDER
C      T       ARRAY OF DIMENSION (IT,JT).  ON EXIT, T CONTAINS
C              THE NORMALISED EIGENVECTORS , WITH T(I,J),
C              I=1(1)N CONTAINING THE EIGENVECTOR CORRESPONDING
C              TO THE EIGENVALUE IN DIAG(J)
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE QLVEC(DIAG, IDIAG, SUB, ISUB, T, IT, JT, N,
C    *                 EPS, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, I1, IDIAG, IERROR, II, ISUB, IT, ITEST,
     *     J, JT, K, L, M, M1, N
      DOUBLE PRECISION B, C, DIAG, EPS, F, G, H, P, R, S, SRNAME,
     *     SUB, T
      DIMENSION DIAG(IDIAG), SUB(ISUB), T(IT,JT)
      DATA SRNAME /8H QLVEC  /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IT.LT.N .OR. JT.LT.N) IERROR = 3
                        IF (IDIAG.LT.N .OR. ISUB.LT.N) IERROR = 2
                        IF (N.LE.0) IERROR = 1
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 IF (N.EQ.1) GO TO 1030
      DO 1020 I=2,N
      SUB(I-1) = SUB(I)
 1020 CONTINUE
 1030 SUB(N) = 0.0D0
      B = 0.0D0
      F = 0.0D0
      DO 1160 L=1,N
      J = 0
      H = EPS*(DABS(DIAG(L))+DABS(SUB(L)))
C+++++
C     LOOK FOR SMALL SUB-DIAGONAL ELEMENT
C
      IF (B.LT.H) B = H
      DO 1040 M=L,N
      IF (DABS(SUB(M)).LE.B) GO TO 1050
 1040 CONTINUE
 1050 IF (M.EQ.L) GO TO 1150
 1060 IF (J.EQ.30) GO TO 1210
C+++++
C     FORM SHIFT
C
      J = J + 1
      G = DIAG(L)
      H = DIAG(L+1) - G
      IF (DABS(H).GE.DABS(SUB(L))) GO TO 1070
      P = H*0.5D0/SUB(L)
      R = DSQRT(P*P+1.0D0)
      H = P + R
      IF (P.LT.0.0D0) H = P - R
      DIAG(L) = SUB(L)/H
      GO TO 1080
 1070 P = 2.0D0*SUB(L)/H
      R = DSQRT(P*P+1.0D0)
      DIAG(L) = SUB(L)*P/(1.0D0+R)
 1080 H = G - DIAG(L)
      I1 = L + 1
      IF (I1.GT.N) GO TO 1100
      DO 1090 I=I1,N
      DIAG(I) = DIAG(I) - H
 1090 CONTINUE
C+++++
C     QL TRANSFORMATION
C
 1100 F = F + H
      P = DIAG(M)
      C = 1.0D0
      S = 0.0D0
      M1 = M - 1
      DO 1140 II=L,M1
      I = M1 - II + L
      G = C*SUB(I)
      H = C*P
      IF (DABS(P).LT.DABS(SUB(I))) GO TO 1110
      C = SUB(I)/P
      R = DSQRT(C*C+1.0D0)
      SUB(I+1) = S*P*R
      S = C/R
      C = 1.0D0/R
      GO TO 1120
 1110 C = P/SUB(I)
      R = DSQRT(C*C+1.0D0)
      SUB(I+1) = S*SUB(I)*R
      S = 1.0D0/R
      C = C/R
 1120 P = C*DIAG(I) - S*G
C+++++
C     FORM VECTOR
C
      DIAG(I+1) = H + S*(C*G+S*DIAG(I))
      DO 1130 K=1,N
      H = T(K,I+1)
      T(K,I+1) = S*T(K,I) + C*H
      T(K,I) = C*T(K,I) - S*H
 1130 CONTINUE
 1140 CONTINUE
      SUB(L) = S*P
      DIAG(L) = C*P
      IF (DABS(SUB(L)).GT.B) GO TO 1060
 1150 DIAG(L) = DIAG(L) + F
C+++++
C     ORDER EIGENVALUES AND EIGENVECTORS
C
 1160 CONTINUE
      DO 1200 I=1,N
      K = I
      P = DIAG(I)
      I1 = I + 1
      IF (I1.GT.N) GO TO 1180
      DO 1170 J=I1,N
      IF (DIAG(J).GE.P) GO TO 1170
      K = J
      P = DIAG(J)
 1170 CONTINUE
 1180 IF (K.EQ.I) GO TO 1200
      DIAG(K) = DIAG(I)
      DIAG(I) = P
      DO 1190 J=1,N
      P = T(J,I)
      T(J,I) = T(J,K)
      T(J,K) = P
 1190 CONTINUE
 1200 CONTINUE
      ITEST = 0
      RETURN
 1210 IF (ITEST.EQ.-1) RETURN
      IERROR = 1
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      RETURN
      END
