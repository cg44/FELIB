C***********************************************************************
C$SPLIT$QLVAL$*********************************************************
C***********************************************************************
      SUBROUTINE QLVAL(DIAG, IDIAG, SUB, ISUB, N, EPS, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CALCULATES ALL THE EIGENVALUES OF A REAL SYMMETRIC
C      TRIDIAGONAL MATRIX OR OF A FULL REAL SYMMETRIC MATRIX
C      THAT HAS BEEN REDUCED TO TRIDIAGONAL FORM USING
C      SUBROUTINE HOUSE
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (IMS) --- SERC COPYRIGHT
C      COMMENTED    15 OCT 1980 (KR)
C
C ARGUMENTS IN
C      DIAG    VECTOR OF LENGTH IDIAG.  ON ENTRY, CONTAINS THE
C              DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX
C      IDIAG   DIMENSION OF VECTOR DIAG (.GE.N)
C      SUB     VECTOR OF DIMENSION ISUB.  ON ENTRY, CONTAINS
C              SUB-DIAGONAL ELEMENTS OF TRIDIAGONAL MATRIX IN
C              ELEMENTS SUB(2) TO SUB(N).  SUB(1) MAY BE
C              ARBITRARY.  CONTENTS DESTROYED DURING EXECUTION
C              OF SUBROUTINE
C      ISUB    DIMENSION OF SUB (.GE.N)
C      N       ORDER OF TRIDIAGONAL MATRIX
C      EPS     SMALLEST POSITIVE NUMBER SUCH THAT 1.+EPS>1.
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      DIAG    VECTOR OF DIMENSION IDIAG.  ON EXIT, CONTAINS
C              EIGENVALUES IN ASCENDING ORDER
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE QLVAL(DIAG, IDIAG, SUB, ISUB, N, EPS, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, I1, IDIAG, IERROR, II, ISUB, ITEST,
     *     J, L, M, M1, N
      DOUBLE PRECISION B, C, DIAG, EPS, F, G, H, P, R, S, SRNAME,
     *     SUB
      DIMENSION DIAG(IDIAG), SUB(ISUB)
      DATA SRNAME /8H QLVAL  /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
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
      DO 1180 L=1,N
      J = 0
      H = EPS*(DABS(DIAG(L))+DABS(SUB(L)))
C+++++
C     LOOK FOR SMALL SUB-DIAGONAL ELEMENT
C
      IF (B.LT.H) B = H
      DO 1040 M=L,N
      IF (DABS(SUB(M)).LE.B) GO TO 1050
 1040 CONTINUE
 1050 IF (M.EQ.L) GO TO 1140
 1060 IF (J.EQ.30) GO TO 1190
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
      DO 1130 II=L,M1
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
      DIAG(I+1) = H + S*(C*G+S*DIAG(I))
 1130 CONTINUE
      SUB(L) = S*P
      DIAG(L) = C*P
      IF (DABS(SUB(L)).GT.B) GO TO 1060
C+++++
C     ORDER EIGENVALUE
C
 1140 P = DIAG(L) + F
      IF (L.EQ.1) GO TO 1160
      DO 1150 II=2,L
      I = L - II + 2
      IF (P.GE.DIAG(I-1)) GO TO 1170
      DIAG(I) = DIAG(I-1)
 1150 CONTINUE
 1160 I = 1
 1170 DIAG(I) = P
 1180 CONTINUE
      ITEST = 0
      RETURN
 1190 IF (ITEST.EQ.-1) RETURN
      IERROR = 3
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      RETURN
      END
