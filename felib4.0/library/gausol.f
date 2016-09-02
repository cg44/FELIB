C
      SUBROUTINE GAUSOL(A,IA,JA,AL,IAL,JAL,N,HBAND,ROPIV,IROPIV,R,IR,
     *                  ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      GAUSOLc alculates the solution of a set of unsymmetric real
C      banded equations with A single rhs using gaussian
C      elimination with partial pivoting
C
C HISTORY
C
C      Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
C                           Chilton, Didcot, Oxfordshire OX11 0QX
C
C      Contact: Prof Chris Greenough
C      Email: c.greenough@rl.ac.uk
C      Tel: (44) 1235 445307
C
C      Release 1.1  29 Oct 1979 (CG)
C      Commented    11 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      A       array of dimension (IA, JA).  On entry, contains
C              the elements of the band matrix
C      IA      first dimension of A (.GE. N)
C      JA      second dimension of A (.GE. MIN(2*HBAND-1, N))
C      IAL     first dimension of array AL (.GE. N)
C      JAL     second dimension of array AL (.GE. HBAND-1)
C      N       order of band matrix A
C      HBAND   semi-bandwidth of matrix A
C      IROPIV  dimension of vector ROPIV (.GE. N)
C      R       on entry, contains the vector of rhs's
C      IR      dimension of R (.GE. N)
C      ITEST   error checking option
C
C ARGUMENTS out
C      A       on exit, contains the upper triangular matrix U,
C              with the diagonal elements of U stored as their
C              reciprocals.  the i'th row of U is stored in the
C              i'th row of A, with the diagonal element of U in
C              the location A(I, 1)
C      AL      contains the sub-diagonal elements of L, the
C              lower triangular matrix.  The multipliers L(I, R)
C              obtained at the r'th major step of the
C              elimination are stored in A(R, I-R)
C      ROPIV   contains details of the row interchanges.
C              ROPIV(R)=R if no interchange occurs at the r'th
C              major step; if rows R and J are interchanged
C              then ROPIV(R)=J
C      R       on exit, contains the solution vector
C
C ROUTINES called
C      VEPS    ERRMES
C
C     SUBROUTINE GAUSOL(A,IA,JA,AL,IAL,JAL,N,HBAND,ROPIV,IROPIV,R,IR,
C    *                  ITEST)
C***********************************************************************
C
      INTEGER ERRMES,HBAND,I,IA,IAL,IERROR,II,IK,IR,IRO,IROPIV,ITEST,IW,
     *        J,JA,JAL,JJ,JR,JTEST,K,KK,L,M,N,ROPIV
      CHARACTER*6 SRNAME
      DOUBLE PRECISION A,AL,EPS,ONE,R,VEPS,X,Y,ZERO
      DIMENSION A(IA,JA),AL(IAL,JAL),R(IR),ROPIV(IROPIV)
C
      INTRINSIC MIN
      EXTERNAL ERRMES, VEPS
C
      DATA ONE/1.0D0/,SRNAME/'GAUSOL'/,ZERO/0.0D0/
C
C     Parameter checking
C
      JTEST = ITEST
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (IA.LT.N .OR. JA.LT.2*HBAND-1) IERROR = 5
         IF (IAL.LT.N .OR. JAL.LT.HBAND) IERROR = 4
         IF (IROPIV.LT.N) IERROR = 3
         IF (IR.LT.N) IERROR = 2
         IF (N.LE.0 .OR. HBAND.LE.0) IERROR = 1
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      IERROR = 6
      EPS = VEPS(X)
      IW = MIN(N,2*HBAND-1)
      M = HBAND
      K = IW - HBAND
      IF (K.GT.0) THEN
         DO 1020 I = 1,K
            L = IW - M
            DO 1000 J = 1,M
               JJ = J + L
               A(I,J) = A(I,JJ)
 1000       CONTINUE
            M = M + 1
            DO 1010 J = M,IW
               A(I,J) = ZERO
 1010       CONTINUE
 1020    CONTINUE
      END IF
      M = N - IW + HBAND + 1
      J = IW + 1
      IF (M.LE.N) THEN
         DO 1040 I = M,N
            J = J - 1
            DO 1030 K = J,IW
               A(I,K) = ZERO
 1030       CONTINUE
 1040    CONTINUE
C
C     Insert zeros
C
      END IF
C
      DO 1060 I = 1,N
         X = ZERO
         DO 1050 J = 1,IW
            X = X + DABS(A(I,J))
 1050    CONTINUE
         IF (X.GT.ZERO) THEN
C
C     ROPIV norms of A calculated and theiro reciprocals stored in
C     firost column of AL
C
            AL(I,1) = ONE/X
         ELSE
            GO TO 1170
         END IF
 1060 CONTINUE
C
      IERROR = 7
      DO 1110 IRO = 1,N
         X = ZERO
         M = MIN(IRO+HBAND-1,N)
         DO 1070 I = IRO,M
            Y = DABS(A(I,1))*AL(I,1)
            IF (Y.GT.X) THEN
               X = Y
               J = I
            END IF
 1070    CONTINUE
C
C     IRO'TH pivot element selected
C
         IF (X.LT.EPS) THEN
            GO TO 1180
         ELSE
            ROPIV(IRO) = J
            IF (J.NE.IRO) THEN
               DO 1080 I = 1,IW
                  X = A(IRO,I)
                  A(IRO,I) = A(J,I)
                  A(J,I) = X
 1080          CONTINUE
C
C     Row pivots IRO and J interchanged
C
               AL(J,1) = AL(IRO,1)
            END IF
            JR = IRO + 1
            Y = ONE/A(IRO,1)
            IF (JR.LE.M) THEN
               DO 1100 I = JR,M
                  X = A(I,1)*Y
                  IF (IW.GE.2) THEN
                     DO 1090 J = 2,IW
                        A(I,J-1) = A(I,J) - X*A(IRO,J)
 1090                CONTINUE
                  END IF
                  IK = I - IRO
                  AL(IRO,IK) = X
                  A(I,IW) = ZERO
 1100          CONTINUE
            END IF
C
C     Elimination complete
C
            A(IRO,1) = Y
         END IF
 1110 CONTINUE
C
      M = HBAND - 1
      DO 1130 K = 1,N
         M = MIN0(M+1,N)
         J = ROPIV(K)
         IF (J.NE.K) THEN
            X = R(K)
C
C     Row pivots K and J interchanged
C
            R(K) = R(J)
            R(J) = X
         END IF
         IK = K + 1
         IF (IK.GT.M) THEN
            GO TO 1140
         ELSE
            X = R(K)
            DO 1120 I = IK,M
               II = I - K
               R(I) = R(I) - X*AL(K,II)
 1120       CONTINUE
C
C     Forward substitution complete
C
         END IF
 1130 CONTINUE
C
 1140 CONTINUE
      DO 1160 K = 1,N
         M = MIN(K,IW)
         I = N + 1 - K
         II = I - 1
         Y = A(I,1)
         X = R(I)
         IF (M.NE.1) THEN
            DO 1150 J = 2,M
               KK = J + II
               X = X - A(I,J)*R(KK)
 1150       CONTINUE
         END IF
C
C     Backward substitution complete
C
         R(I) = X*Y
 1160 CONTINUE
C
      RETURN
 1170 CONTINUE
      IRO = I
C
 1180 CONTINUE
      A(IRO,1) = ZERO
      ITEST = ERRMES(JTEST,IERROR,SRNAME)
C
      END
