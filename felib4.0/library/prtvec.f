C
      SUBROUTINE PRTVEC(V,IV,N,NOUT,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      PRTVEC prints a vector in a standard format
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
C      V       vector of dimension IV containing numbers to be
C              printed
C      IV      dimension of V (.GE. N)
C      N       number of elements of V to be printed
C      NOUT    fortran unit number
C      ITEST   error checking option
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE PRTVEC(V,IV,N,NOUT,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IERROR,ITEST,IV,N,NOUT
      CHARACTER*6 SRNAME
      DOUBLE PRECISION V
      DIMENSION V(IV)
      DATA SRNAME/'PRTVEC'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IV.LT.N) IERROR = 2
         IF (N.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      WRITE (NOUT,FMT=9990) (V(I),I=1,N)
C
 9990 FORMAT (' ',6D12.4)
      END
