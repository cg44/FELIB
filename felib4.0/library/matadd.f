C
      SUBROUTINE MATADD(A,IA,JA,B,IB,JB,M,N,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      MATADD adds the matrix A to B, storing the result in A
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
C      IA      first dimension of A (.GE. M)
C      JA      second dimension of A (.GE. N)
C      B       array of dimension (IB, JB)
C      IB      first dimension of B (.GE. M)
C      JB      second dimension of B (.GE. N)
C      M       number of rows of A and B to be added
C      N       number of columns of A and B to be added
C      ITEST   error checking option
C
C ARGUMENTS out
C      A       on exit, contains sum of A and B
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE MATADD(A,IA,JA,B,IB,JB,M,N,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IA,IB,IERROR,ITEST,J,JA,JB,M,N
      CHARACTER*6 SRNAME
      DOUBLE PRECISION A,B
      DIMENSION A(IA,JA),B(IB,JB)
      DATA SRNAME/'MATADD'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (M.GT.IB .OR. N.GT.JB) IERROR = 3
         IF (M.GT.IA .OR. N.GT.JA) IERROR = 2
         IF (M.LE.0 .OR. N.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      DO 1010 J = 1,N
         DO 1000 I = 1,M
            A(I,J) = A(I,J) + B(I,J)
 1000    CONTINUE
 1010 CONTINUE
C
      END
