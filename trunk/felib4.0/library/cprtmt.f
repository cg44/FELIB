C
      SUBROUTINE CPRTMT(A,IA,JA,M,N,NOUT,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CPRTMT prints A complex matrix in a standard format
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
C      Release 2.0  29 Oct 1984 (CRIE)
C      Comments      1 Nov 1985 (CG)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      A       matrix of dimension (2, IA, JA) containing numbers to be
C              printed
C      IA      first dimension of A (.GE. M)
C      JA      second dimension of A (.GE. N)
C      M       number of rows of A to be printed (.LE. IA)
C      N       number of columes of A to be printed (.LE. JA)
C      NOUT    fortran unit number
C      ITEST   error checking option
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE CPRTMT(A,IA,JA,M,N,NOUT,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IERROR,ITEST,IA,M,N,NOUT,J,JA,K
      CHARACTER*6 SRNAME
      DOUBLE PRECISION A
      DIMENSION A(2,IA,JA)
      DATA SRNAME/'CPRTMT'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IA.LT.M .OR. JA.LT.N) IERROR = 2
         IF (N.LE.0 .OR. M.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main loops
C
      WRITE (NOUT,FMT=9990) (((A(I,J,K),I=1,2),J=1,M),K=1,N)
C
 9990 FORMAT (' ',3 (2X,2D12.4))
      END
