C***********************************************************************
C$SPLIT$GAURDN$*********************************************************
C***********************************************************************
      SUBROUTINE GAURDN(A, IA, JA, AL, IAL, JAL, N, HBAND, ROPIV,
     *     IROPIV, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      DECOMPOSES A REAL UNSYMMETRIC MATRIX OF ORDER N INTO
C      TRIANGULAR MATRICES USING GAUSSIAN ELIMINATION WITH
C      PARTIAL PIVOTING
C
C HISTORY
C      RELEASE 1.1 29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED   10 OCT 1980 (KR)
C
C ARGUMENTS IN
C      A       ARRAY OF DIMENSIONS (IA,JA).  ON ENTRY, CONTAINS
C              THE ELEMENTS OF THE BAND MATRIX
C      IA      FIRST DIMENSION OF A (.GE.N)
C      JA      SECOND DIMENSION OF A (.GE.MIN(2*HBAND-1,N))
C      IAL     FIRST DIMENSION OF ARRAY AL (.GE.N)
C      JAL     SECOND DIMENSION OF ARRAY AL (.GE.HBAND-1)
C      N       ORDER OF MATRIX A
C      HBAND   SEMI-BANDWIDTH OF MATRIX A
C      IROPIV  DIMENSION OF VECTOR ROPIV
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      A       ARRAY OF DIMENSIONS (IA,JA).  ON EXIT, CONTAINS
C              THE UPPER TRIANGULAR MATRIX U, WITH THE
C              DIAGONAL ELEMENTS OF U STORED AS THEIR
C              RECIPROCALS.  THE I'TH ROW OF U IS STORED IN THE
C              I'TH ROW OF A, STARTING WITH THE DIAGONAL
C              ELEMENT OF U IN A(I,1)
C      ROPIV   VECTOR OF LENGTH IROPIV, CONTAINING DETAILS OF
C              ROW INTERCHANGES.  IF NO INTERCHANGE OCCURS AT
C              THE R'TH MAJOR STEP THEN ROPIV(R)=R; IF THE R
C              AND J ROWS ARE INTERCHANGED THEN ROPIV(R)=J
C
C ROUTINES CALLED
C      VEPS    ERRMES
C
C
C     SUBROUTINE GAURDN(A, IA, JA, AL, IAL, JAL, N, HBAND, ROPIV,
C    *     IROPIV, ITEST)
C***********************************************************************
      INTEGER ERRMES, HBAND, I, IA, IAL, IERROR, IK, IR, IROPIV,
     *     ITEST, IW, J, JA, JAL, JJ, JR, K, L, M, N, ROPIV, JTEST
      DOUBLE PRECISION A, AL, EPS, ONE, SRNAME, VEPS, X, Y,
     *     ZERO
      DIMENSION A(IA,JA), AL(IAL,JAL), ROPIV(IROPIV)
      DATA ONE /1.0D0/, SRNAME /8H GAURDN /, ZERO /0.0D0/
      JTEST = ITEST
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IA.LT.N .OR. JA.LT.2*HBAND-1)
     *                      IERROR = 4
                        IF (IAL.LT.N .OR. JAL.LT.HBAND) IERROR = 3
                        IF (IROPIV.LT.N) IERROR = 2
                        IF (N.LE.0 .OR. HBAND.LE.0) IERROR = 1
                        ITEST = ERRMES(JTEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 IERROR = 5
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
C     ZEROS INSERTED
C
 1060 CONTINUE
 1070 CONTINUE
 1080 DO 1110 I=1,N
      X = ZERO
      DO 1090 J=1,IW
      X = X + DABS(A(I,J))
 1090 CONTINUE
      IF (X.GT.ZERO) GO TO 1100
      IR = I
      GO TO 1200
C+++++
C     ROW NORMS OF A CALCULATED AND THEIR RECIPROCALS STORED
C     IN FIRST COLUMN OF AL
C
 1100 AL(I,1) = ONE/X
 1110 CONTINUE
      IERROR = 6
      DO 1190 IR=1,N
      X = ZERO
      M = MIN0(IR+HBAND-1,N)
      DO 1120 I=IR,M
      Y = DABS(A(I,1))*AL(I,1)
      IF (Y.LE.X) GO TO 1120
      X = Y
      J = I
 1120 CONTINUE
C+++++
C     IR'TH PIVOT ELEMENT SELECTED.
C
      IF (X.LT.EPS) GO TO 1200
      ROPIV(IR) = J
      IF (J.EQ.IR) GO TO 1140
      DO 1130 I=1,IW
      X = A(IR,I)
      A(IR,I) = A(J,I)
      A(J,I) = X
C+++++
C     ROW IR AND J INTERCHANGED.
C
 1130 CONTINUE
      AL(J,1) = AL(IR,1)
 1140 JR = IR + 1
      Y = ONE/A(IR,1)
      IF (JR.GT.M) GO TO 1180
      DO 1170 I=JR,M
      X = A(I,1)*Y
      IF (IW.LT.2) GO TO 1160
      DO 1150 J=2,IW
      A(I,J-1) = A(I,J) - X*A(IR,J)
 1150 CONTINUE
 1160 IK = I - IR
      AL(IR,IK) = X
      A(I,IW) = ZERO
 1170 CONTINUE
C+++++
C     ELIMINATION COMPLETED
C
 1180 A(IR,1) = Y
 1190 CONTINUE
      RETURN
 1200 A(IR,1) = ZERO
      ITEST = ERRMES(JTEST,IERROR,SRNAME)
      RETURN
      END
