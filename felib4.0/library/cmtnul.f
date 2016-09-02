C
      SUBROUTINE CMTNUL(A,IA,JA,M,N,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CMTNUL sets complex matrix A to the null matrix
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
C      Release 2.0  20 Oct 1984 (CRIE)
C      Comments      1 Nov 1985 (CG)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      IA      first dimension of array A (.GE.M)
C      JA      second dimension of array A (.GE.N)
C      M       number of rows of A to be set to zero
C      N       number of columns of A to be set to zero
C      ITEST   error checking option
C
C ARGUMENTS out
C      A       array of dimension (2, IA, JA). A(1, I, J)=0
C              and A(2, I, J)=0 for I=1(1)M and J=1(1)N
C              
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE CMTNUL(A,IA,JA,M,N,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IA,IERROR,ITEST,J,JA,M,N
      CHARACTER*6 SRNAME
      DOUBLE PRECISION A
      DIMENSION A(2,IA,JA)
      DATA SRNAME/'CMTNUL'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (M.GT.IA .OR. N.GT.JA) IERROR = 2
         IF (M.LE.0 .OR. N.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main loops
C
      DO 1010 I = 1,M
         DO 1000 J = 1,N
            A(1,I,J) = 0.0D0
            A(2,I,J) = 0.0D0
 1000    CONTINUE
 1010 CONTINUE
C
      END
