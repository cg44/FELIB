C
      SUBROUTINE VECADD(V,IV,W,IW,N,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      VECADD adds the vector V to the vector W, storing the result
C      in V
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
C      V       vector of length IV. On entry, contains the
C              values of one of the vectors to be added
C              in elements 1 to N
C      IV      length of vector V (.GE. N)
C      W       vector of length IW. Contains the values of the
C              second vector to be added in elements 1 to N
C      IW      length of vector W (.GE. N)
C      N       number of elements of V and W to be added
C      ITEST   error checking option
C
C ARGUMENTS out
C      V       vector of length IV. On exit, V(I) contains
C              the sum of V(I) and W(I) for I=1(1)N
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE VECADD(V,IV,W,IW,N,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IERROR,ITEST,IV,IW,N
      CHARACTER*6 SRNAME
      DOUBLE PRECISION V,W
      DIMENSION V(IV),W(IW)
C
      DATA SRNAME/'VECADD'/
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
         V(I) = V(I) + W(I)
 1000 CONTINUE
C
      END
