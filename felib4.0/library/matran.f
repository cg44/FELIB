C
      SUBROUTINE MATRAN(A,IA,JA,B,IB,JB,M,N,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      MATRAN forms the transpose of the matrix A
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
C      Release 1.1  20 Oct 1979 (IMS)
C      Commented    12 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      A       array of dimension (IA, JA) to be transposed
C      IA      first dimension of A (.GE. M)
C      JA      second dimension of A (.GE. N)
C      IB      first dimension of array B (.GE. N)
C      JB      second dimension of array B (.GE. M)
C      M       number of rows of A to be transposed
C      N       number of columns of A to be transposed
C      ITEST   error checking option
C
C ARGUMENTS out
C      B       array of dimension (IB, JB). B(J, I) = A(I, J) for
C              I=1(1)M and J=1(1)N
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE MATRAN(A,IA,JA,B,IB,JB,M,N,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IA,IB,IERROR,ITEST,J,JA,JB,M,N
      CHARACTER*6 SRNAME
      DOUBLE PRECISION A,B
      DIMENSION A(IA,JA),B(IB,JB)
      DATA SRNAME/'MATRAN'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (N.GT.IB .OR. M.GT.JB) IERROR = 3
         IF (M.GT.IA .OR. N.GT.JA) IERROR = 2
         IF (M.LE.0 .OR. N.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (IERROR.NE.0) RETURN
      END IF
C
C     Main body
C
      DO 1010 I = 1,M
         DO 1000 J = 1,N
            B(J,I) = A(I,J)
 1000    CONTINUE
 1010 CONTINUE
C
      END
