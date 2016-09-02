C
      SUBROUTINE VECSUB(V,IV,W,IW,N,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      VECSUB subtracts the vector W from vector V, storing the
C      result in V.
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
C      Commented    23 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      V       vector of length IV. On entry, V(I), I=1(1)N,
C              contains values from which corresponding values
C              W(I) are to be subtracted
C      IV      length of V (.GE. N)
C      W       vector of length IW. The elements W(I), I=1(1)N
C              are to be subtracted from the corresponding
C              elements of V
C      IW      length of W (.GE. N)
C      N       number of elements of V and W to be operated on
C      ITEST   error checking option
C
C ARGUMENTS out
C      V       vector of length IV. On exit, V(I) is set to
C              V(I)-W(I) for I=1(1)N
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE VECSUB(V,IV,W,IW,N,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IERROR,ITEST,IV,IW,N
      CHARACTER*6 SRNAME
      DOUBLE PRECISION V,W
      DIMENSION V(IV),W(IW)
C
      DATA SRNAME/'VECSUB'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (N.GT.IW) IERROR = 3
         IF (N.GT.IV) IERROR = 2
         IF (N.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      DO 1000 I = 1,N
         V(I) = V(I) - W(I)
 1000 CONTINUE
C
      END
