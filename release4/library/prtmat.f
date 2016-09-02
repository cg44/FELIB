C
      SUBROUTINE PRTMAT(A,IA,JA,M,N,NOUT,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      PRTMAT prints a two-dimensional array in a standard format
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
C      Commented    14 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      A       array of diimension (IA, JA) containing numbers
C              to be printed
C      IA      first dimension of A (.GE. M)
C      JA      second dimension of A (.GE. N)
C      M       number of rows of A to be printed
C      N       number of columns of A to be printed
C      NOUT    fortran unit number
C      ITEST   error checking option
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE PRTMAT(A,IA,JA,M,N,NOUT,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IA,IERROR,ITEST,J,JA,M,N,NOUT
      CHARACTER*6 SRNAME
      DOUBLE PRECISION A
      DIMENSION A(IA,JA)
      DATA SRNAME/'PRTMAT'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IA.LT.M .OR. JA.LT.N) IERROR = 2
         IF (M.LE.0 .OR. N.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      DO 1000 I = 1,M
         WRITE (NOUT,FMT=9990) (A(I,J),J=1,N)
 1000 CONTINUE
C
 9990 FORMAT (' ',6D12.4)
      END
