C
      SUBROUTINE CHOSOL(A,IA,JA,R,IR,N,HBAND,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CHOSOL solves a set of real symmetric positive definite banded
C      equations with a single right hand side by choleski decomposition. 
C      Only the lower band and diagonal are stored in a rectangular 
C      array A.
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
C      Release 1.1  29 Oct 1979 (IMS)
C      Commented    12 Feb 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      A       on entry contains lower half of pd symmetric
C              band matrix stored as a rectangular array
C      IA      first dimension of A (.GE. N)
C      JA      second dimension of A (.GE. HBAND)
C      R       contains elements of right hand side
C      IR      dimension of R (.GE. N)
C      N       order of matrix A
C      HBAND   semi-bandwidth of A (includes diagonal)
C      ITEST   error checking option
C
C ARGUMENTS out
C      A       on exit, contains lower triangular reduced
C      R       matrix solution vector
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE CHOSOL(A,IA,JA,R,IR,N,HBAND,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,HBAND,I,IA,IERROR,IJ,IK,IR,ITEST,J,JA,JTEST,K,L,LA,
     *        LB,LK,M,N,W
      CHARACTER*6 SRNAME
      DOUBLE PRECISION A,R,X
      DIMENSION A(IA,JA),R(IR)
      DATA SRNAME/'CHOSOL'/
C
C     Parameter checking
C
      JTEST = ITEST
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IR.LT.N) IERROR = 3
         IF (IA.LT.N .OR. JA.LT.HBAND) IERROR = 2
         IF (N.LE.0 .OR. HBAND.LE.0) IERROR = 1
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      W = HBAND - 1
      DO 1030 I = 1,N
         X = 0.0D0
         DO 1000 J = 1,W
            X = X + A(I,J)*A(I,J)
 1000    CONTINUE
C
C     Range checking on A(I,W+1)
C
         IF (JTEST.NE.-1) THEN
            IERROR = 0
            IF ((A(I,W+1)-X).LE.0.0D0) IERROR = 4
            ITEST = ERRMES(JTEST,IERROR,SRNAME)
            IF (ITEST.NE.0) RETURN
         END IF
C
         A(I,W+1) = DSQRT(A(I,W+1)-X)
C
         DO 1020 K = 1,W
            X = 0.0D0
            IF (I+K.LE.N) THEN
               IF (K.NE.W) THEN
                  L = W - K
 1010             CONTINUE
                  IK = I + K
                  LK = L + K
                  X = X + A(IK,L)*A(I,LK)
                  L = L - 1
                  IF (L.NE.0) GO TO 1010
               END IF
               LA = I + K
               LB = W - K + 1
               A(LA,LB) = (A(LA,LB)-X)/A(I,W+1)
            END IF
 1020    CONTINUE
 1030 CONTINUE
      R(1) = R(1)/A(1,W+1)
      DO 1050 I = 2,N
         X = 0.0D0
         K = 1
         IF (I.LE.W+1) K = W - I + 2
         DO 1040 J = K,W
            IJ = I + J - W - 1
            X = X + A(I,J)*R(IJ)
 1040    CONTINUE
         R(I) = (R(I)-X)/A(I,W+1)
 1050 CONTINUE
      R(N) = R(N)/A(N,W+1)
      I = N - 1
 1060 CONTINUE
      X = 0.0D0
      L = I + W
      IF (I.GT.N-W) L = N
      M = I + 1
      DO 1070 J = M,L
         IJ = W + I - J + 1
         X = X + A(J,IJ)*R(J)
 1070 CONTINUE
      R(I) = (R(I)-X)/A(I,W+1)
      I = I - 1
      IF (I.NE.0) GO TO 1060
C
      END
