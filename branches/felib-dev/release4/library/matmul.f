C
      SUBROUTINE MATMUL(A,IA,JA,B,IB,JB,C,IC,JC,L,M,N,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      MATMUL pre-multiplies matrix B by A, storing the result in C
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
C      Commented    12 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      A       array of dimension (IA, JA)
C      IA      first dimension of A (.GE. L)
C      JA      second dimension of A (.GE. M)
C      B       array of dimension (IB, JB)
C      IB      first dimension of B (.GE. M)
C      JB      second dimension of B (.GE. N)
C      IC      first dimension of array C (.GE. L)
C      JC      second dimension of array C (.GE. N)
C      L       number of rows of A to be used in multiplication
C      M       number of columns of A and number of rows of B
C              to be used in multiplication
C      N       number of columns of B to be used in
C              multiplication
C      ITEST   error checking option
C
C ARGUMENTS out
C      C       contains result of matrix multiplication (C=A*B)
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE MATMUL(A,IA,JA,B,IB,JB,C,IC,JC,L,M,N,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IA,IB,IC,IERROR,ITEST,J,JA,JB,JC,K,L,M,N
      CHARACTER*6 SRNAME
      DOUBLE PRECISION A,B,C,X
      DIMENSION A(IA,JA),B(IB,JB),C(IC,JC)
      DATA SRNAME/'MATMUL'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (L.GT.IC .OR. N.GT.JC) IERROR = 4
         IF (M.GT.IB .OR. N.GT.JB) IERROR = 3
         IF (L.GT.IA .OR. M.GT.JA) IERROR = 2
         IF (L.LE.0 .OR. M.LE.0 .OR. N.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      DO 1020 I = 1,L
         DO 1010 J = 1,N
            X = 0.0D0
            DO 1000 K = 1,M
               X = X + A(I,K)*B(K,J)
 1000       CONTINUE
            C(I,J) = X
 1010    CONTINUE
 1020 CONTINUE
C
      END
