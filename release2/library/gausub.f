C***********************************************************************
C$SPLIT$GAUSUB$*********************************************************
C***********************************************************************
      SUBROUTINE GAUSUB(A, IA, JA, AL, IAL, JAL, N, HBAND, ROPIV,
     *     IROPIV, R, IR, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CALCULATES THE SOLUTION OF A SET OF UNSYMMETRIC REAL
C      BANDED LINEAR EQUATIONS WITH A SINGLE RHS.  THE BANDED
C      MATRIX HAS PREVIOUSLY BEEN DECOMPOSED INTO TRIANGULAR
C      MATRICES USING GAURDN
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    11 OCT 1980 (KR)
C
C ARGUMENTS IN
C      A       ARRAY OF DIMENSIONS (IA,JA).  ON ENTRY, CONTAINS
C              THE ELEMENTS OF THE BAND MATRIX IN LU FORM,
C              AFTER PROCESSING BY GAURDN
C      IA      FIRST DIMENSION OF A (.GE.N)
C      JA      SECOND DIMENSION OF A (.GE.MIN(2*HBAND-1,N))
C      IAL     FIRST DIMENSION OF AL (.GE.N)
C      JAL     SECOND DIMENSION OF AL (.GE.HBAND-1)
C      N       ORDER OF BANDED MATRIX A
C      HBAND   SEMI-BANDWIDTH OF A
C      ROPIV   VECTOR OF DIMENSION IROPIV.  CONTAINS DETAILS
C              OF ROW INTERCHANGES PERFORMED BY GAURDN
C      IROPIV  DIMENSION OF ROPIV (.GE.N)
C      R       ON ENTRY, CONTAINS THE VECTOR OF THE RHS,
C              LENGTH IR
C      IR      DIMENSION OF R (.GE.N)
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      R       ON EXIT, CONTAINS SOLUTION VECTOR
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE GAUSUB(A, IA, JA, AL, IAL, JAL, N, HBAND, ROPIV,
C    *     IROPIV, R, IR, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, HBAND, I, IA, IAL, IERROR, II, IK, IR,
     *     IROPIV, ITEST, IW, J, JA, JAL, K, KK, M, N, ROPIV
      DOUBLE PRECISION A, AL, R, SRNAME, X, Y
      DIMENSION A(IA,JA), AL(IAL,JAL), R(IR), ROPIV(IROPIV)
      DATA SRNAME /8H GAUSUB /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IR.LT.N) IERROR = 5
                        IF (IA.LT.N .OR. JA.LT.2*HBAND-1)
     *                      IERROR = 4
                        IF (IAL.LT.N .OR. JAL.LT.HBAND) IERROR = 3
                        IF (IROPIV.LT.N) IERROR = 2
                        IF (N.LE.0 .OR. HBAND.LE.0) IERROR = 1
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 IW = MIN0(N,2*HBAND-1)
      M = HBAND - 1
      DO 1040 K=1,N
      M = MIN0(M+1,N)
      J = ROPIV(K)
      IF (J.EQ.K) GO TO 1020
      X = R(K)
C+++++
C     ROWS K AND J INTERCHANGED
C
      R(K) = R(J)
      R(J) = X
 1020 IK = K + 1
      IF (IK.GT.M) GO TO 1050
      X = R(K)
      DO 1030 I=IK,M
      II = I - K
      R(I) = R(I) - X*AL(K,II)
C+++++
C     FORWARD SUBSTITUTION COMPLETE
C
 1030 CONTINUE
 1040 CONTINUE
 1050 DO 1080 K=1,N
      M = MIN0(K,IW)
      I = N + 1 - K
      II = I - 1
      Y = A(I,1)
      X = R(I)
      IF (M.EQ.1) GO TO 1070
      DO 1060 J=2,M
      KK = J + II
      X = X - A(I,J)*R(KK)
 1060 CONTINUE
C+++++
C     BACKWARD SUBSTITUTION COMPLETE
C
 1070 R(I) = X*Y
 1080 CONTINUE
      RETURN
      END
