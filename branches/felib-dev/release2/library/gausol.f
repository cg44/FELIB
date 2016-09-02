C***********************************************************************
C$SPLIT$GAUSOL$*********************************************************
C***********************************************************************
      SUBROUTINE GAUSOL(A, IA, JA, AL, IAL, JAL, N, HBAND, ROPIV,
     *     IROPIV, R, IR, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CALCULATES THE SOLUTION OF A SET OF UNSYMMETRIC REAL
C      BANDED EQUATIONS WITH A SINGLE RHS USING GAUSSIAN
C      ELIMINATION WITH PARTIAL PIVOTING
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    11 OCT 1980 (KR)
C
C ARGUMENTS IN
C      A       ARRAY OF DIMENSION (IA,JA).  ON ENTRY, CONTAINS
C              THE ELEMENTS OF THE BAND MATRIX
C      IA      FIRST DIMENSION OF A (.GE.N)
C      JA      SECOND DIMENSION OF A (.GE.MIN(2*HBAND-1,N))
C      IAL     FIRST DIMENSION OF ARRAY AL (.GE.N)
C      JAL     SECOND DIMENSION OF ARRAY AL (.GE.HBAND-1)
C      N       ORDER OF BAND MATRIX A
C      HBAND   SEMI-BANDWIDTH OF MATRIX A
C      IROPIV  DIMENSION OF VECTOR ROPIV (.GE.N)
C      R       ON ENTRY, CONTAINS THE VECTOR OF RHS'S
C      IR      DIMENSION OF R (.GE.N)
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      A       ON EXIT, CONTAINS THE UPPER TRIANGULAR MATRIX U,
C              WITH THE DIAGONAL ELEMENTS OF U STORED AS THEIR
C              RECIPROCALS.  THE I'TH ROW OF U IS STORED IN THE
C              I'TH ROW OF A, WITH THE DIAGONAL ELEMENT OF U IN
C              THE LOCATION A(I,1)
C      AL      CONTAINS THE SUB-DIAGONAL ELEMENTS OF L, THE
C              LOWER TRIANGULAR MATRIX.  THE MULTIPLIERS L(I,R)
C              OBTAINED AT THE R'TH MAJOR STEP OF THE
C              ELIMINATION ARE STORED IN A(R,I-R)
C      ROPIV   CONTAINS DETAILS OF THE ROW INTERCHANGES.
C              ROPIV(R)=R IF NO INTERCHANGE OCCURS AT THE R'TH
C              MAJOR STEP; IF ROWS R AND J ARE INTERCHANGED
C              THEN ROPIV(R)=J
C      R       ON EXIT, CONTAINS THE SOLUTION VECTOR
C
C ROUTINES CALLED
C      VEPS    ERRMES
C
C
C     SUBROUTINE GAUSOL(A, IA, JA, AL, IAL, JAL, N, HBAND, ROPIV,
C    *     IROPIV, R, IR, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, HBAND, I, IA, IAL, IERROR, II, IK, IR,
     *     IRO, IROPIV, ITEST, IW, J, JA, JAL, JJ, JR, K, KK,
     *     L, M, N, ROPIV, JTEST
      DOUBLE PRECISION A, AL, EPS, ONE, R, SRNAME, VEPS, X,
     *     Y, ZERO
      DIMENSION A(IA,JA), AL(IAL,JAL), R(IR), ROPIV(IROPIV)
      DATA ONE /1.0D0/, SRNAME /8H GAUSOL /, ZERO /0.0D0/
      JTEST = ITEST
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IA.LT.N .OR. JA.LT.2*HBAND-1)
     *                      IERROR = 5
                        IF (IAL.LT.N .OR. JAL.LT.HBAND) IERROR = 4
                        IF (IROPIV.LT.N) IERROR = 3
                        IF (IR.LT.N) IERROR = 2
                        IF (N.LE.0 .OR. HBAND.LE.0) IERROR = 1
                        ITEST = ERRMES(JTEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 IERROR = 6
      EPS = VEPS(X)
      IW = MIN0(N,2*HBAND-1)
      M = HBAND
      K = IW - HBAND
      IF (K.LE.0) GO TO 1050
      DO 1040 I=1,K
      L = IW - M
      DO 1020 J=1,M
      JJ = J + L
      A(I,J) = A(I,JJ)
 1020 CONTINUE
      M = M + 1
      DO 1030 J=M,IW
      A(I,J) = ZERO
 1030 CONTINUE
 1040 CONTINUE
 1050 M = N - IW + HBAND + 1
      J = IW + 1
      IF (M.GT.N) GO TO 1080
      DO 1070 I=M,N
      J = J - 1
      DO 1060 K=J,IW
      A(I,K) = ZERO
C+++++
C     INSERT ZEROS
C
 1060 CONTINUE
 1070 CONTINUE
 1080 DO 1110 I=1,N
      X = ZERO
      DO 1090 J=1,IW
      X = X + DABS(A(I,J))
 1090 CONTINUE
      IF (X.GT.ZERO) GO TO 1100
      IRO = I
      GO TO 1270
C+++++
C     ROPIV NORMS OF A CALCULATED AND THEIRO RECIPROCALS
C     STORED IN FIROST COLUMN OF AL
C
 1100 AL(I,1) = ONE/X
 1110 CONTINUE
      IERROR = 7
      DO 1190 IRO=1,N
      X = ZERO
      M = MIN0(IRO+HBAND-1,N)
      DO 1120 I=IRO,M
      Y = DABS(A(I,1))*AL(I,1)
      IF (Y.LE.X) GO TO 1120
      X = Y
      J = I
 1120 CONTINUE
C+++++
C     IRO'TH PIVOT ELEMENT SELECTED
C
      IF (X.LT.EPS) GO TO 1270
      ROPIV(IRO) = J
      IF (J.EQ.IRO) GO TO 1140
      DO 1130 I=1,IW
      X = A(IRO,I)
      A(IRO,I) = A(J,I)
      A(J,I) = X
C+++++
C     ROW PIVOTS IRO AND J INTERCHANGED
C
 1130 CONTINUE
      AL(J,1) = AL(IRO,1)
 1140 JR = IRO + 1
      Y = ONE/A(IRO,1)
      IF (JR.GT.M) GO TO 1180
      DO 1170 I=JR,M
      X = A(I,1)*Y
      IF (IW.LT.2) GO TO 1160
      DO 1150 J=2,IW
      A(I,J-1) = A(I,J) - X*A(IRO,J)
 1150 CONTINUE
 1160 IK = I - IRO
      AL(IRO,IK) = X
      A(I,IW) = ZERO
 1170 CONTINUE
C+++++
C     ELIMINATION COMPLETE
C
 1180 A(IRO,1) = Y
 1190 CONTINUE
      M = HBAND - 1
      DO 1220 K=1,N
      M = MIN0(M+1,N)
      J = ROPIV(K)
      IF (J.EQ.K) GO TO 1200
      X = R(K)
C+++++
C     ROW PIVOTS K AND J INTERCHANGED
C
      R(K) = R(J)
      R(J) = X
 1200 IK = K + 1
      IF (IK.GT.M) GO TO 1230
      X = R(K)
      DO 1210 I=IK,M
      II = I - K
      R(I) = R(I) - X*AL(K,II)
C+++++
C     FORWARD SUBSTITUTION COMPLETE
C
 1210 CONTINUE
 1220 CONTINUE
 1230 DO 1260 K=1,N
      M = MIN0(K,IW)
      I = N + 1 - K
      II = I - 1
      Y = A(I,1)
      X = R(I)
      IF (M.EQ.1) GO TO 1250
      DO 1240 J=2,M
      KK = J + II
      X = X - A(I,J)*R(KK)
 1240 CONTINUE
C+++++
C     BACKWARD SUBSTITUTION COMPLETE
C
 1250 R(I) = X*Y
 1260 CONTINUE
      RETURN
 1270 A(IRO,1) = ZERO
      ITEST = ERRMES(JTEST,IERROR,SRNAME)
      RETURN
      END
