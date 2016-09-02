C
      SUBROUTINE MATIDN(A,IA,JA,M,N,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      MATIDN sets the matrix A to the identity matrix
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
C      Commented    12 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      IA      first dimension of array A (.GE. M)
C      JA      second dimension of array A (.GE. N)
C      M       number of rows of A to be assigned values
C      N       number of columns of A to be assigned values
C      ITEST   error checking option
C
C ARGUMENTS out
C      A       array of dimension (IA, JA).  A(I, J) is set to 1
C              if I=J, 0 otherwise
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE MATIDN(A,IA,JA,M,N,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IA,IERROR,ITEST,J,JA,L,M,N
      CHARACTER*6 SRNAME
      DOUBLE PRECISION A
      DIMENSION A(IA,JA)
      DATA SRNAME/'MATIDN'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (M.GT.IA .OR. N.GT.JA) IERROR = 2
         IF (N.LE.0 .OR. M.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (IERROR.NE.0) RETURN
      END IF
C
C     Main body
C
      DO 1010 J = 1,N
         DO 1000 I = 1,M
            A(I,J) = 0.0D0
 1000    CONTINUE
 1010 CONTINUE
C
C     Set diagonal entry
C
      L = MIN(M,N)
      DO 1020 I = 1,L
         A(I,I) = 1.0D0
 1020 CONTINUE
C
      END
