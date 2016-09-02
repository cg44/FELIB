C
      SUBROUTINE MATINV(A,IA,JA,B,IB,JB,N,DET,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      MATINV forms the inverse of matrix A in B
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
C      A       array of dimension (IA, JA) to be inverted
C      IA      first dimension of A (.GE. N)
C      JA      second dimension of A (.GE. N)
C      IB      first dimension of array B (.GE. N)
C      JB      second dimension of array B (.GE. N)
C      N       order of matrix A (.GT. 1 .AND. .LT. 4)
C      ITEST   error checking option
C
C ARGUMENTS out
C      B       array of dimension (IB, JB) containing inverse of
C              A
C      DET     contains value of determinant of matrix A
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE MATINV(A,IA,JA,B,IB,JB,N,DET,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,IA,IB,IERROR,ITEST,JA,JB,JTEST,K,L,M,N
      CHARACTER*6 SRNAME
      DOUBLE PRECISION A,B,DET,UNFLO,X
      DIMENSION A(IA,JA),B(IB,JB)
      DATA SRNAME/'MATINV'/
C
C     Parameter checking
C
      JTEST = ITEST
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (N.GT.IB .OR. N.GT.JB) IERROR = 3
         IF (N.GT.IA .OR. N.GT.JA) IERROR = 2
         IF (N.LE.1 .OR. N.GE.4) IERROR = 1
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      M = N - 1
      GO TO (1000,1030) M
C
C     Code for 2x2 matrix
C
 1000 CONTINUE
      DET = A(1,1)*A(2,2) - A(1,2)*A(2,1)
C
C     Check that determinant is not near zero
C
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (DABS(DET).LT.UNFLO(X)) IERROR = 4
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
      B(1,1) = A(2,2)
      B(1,2) = -A(1,2)
      B(2,1) = -A(2,1)
      B(2,2) = A(1,1)
      DO 1020 K = 1,2
         DO 1010 L = 1,2
            B(K,L) = B(K,L)/DET
 1010    CONTINUE
 1020 CONTINUE
      RETURN
C
C     Code for 3x3 matrix
C
 1030 CONTINUE
      DET = A(1,1)* (A(2,2)*A(3,3)-A(3,2)*A(2,3))
      DET = DET - A(1,2)* (A(2,1)*A(3,3)-A(3,1)*A(2,3))
      DET = DET + A(1,3)* (A(2,1)*A(3,2)-A(3,1)*A(2,2))
C
C     Check on determinant not near zero
C
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (DABS(DET).LT.UNFLO(X)) IERROR = 4
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
      B(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
      B(2,1) = -A(2,1)*A(3,3) + A(3,1)*A(2,3)
      B(3,1) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
      B(1,2) = -A(1,2)*A(3,3) + A(3,2)*A(1,3)
      B(2,2) = A(1,1)*A(3,3) - A(3,1)*A(1,3)
      B(3,2) = -A(1,1)*A(3,2) + A(3,1)*A(1,2)
      B(1,3) = A(1,2)*A(2,3) - A(2,2)*A(1,3)
      B(2,3) = -A(1,1)*A(2,3) + A(2,1)*A(1,3)
      B(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)
C
      DO 1050 K = 1,3
         DO 1040 L = 1,3
            B(K,L) = B(K,L)/DET
 1040    CONTINUE
 1050 CONTINUE
C
      END
