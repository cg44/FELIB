C
      SUBROUTINE DYAD(V,IV,W,IW,A,IA,JA,N,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      DYAD forms the dyad matrix A from two vectors V and W
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
C      Commented    14 Feb 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      V       first vector in multipcation
C      IV      dimension of vector V (.GE. N)
C      W       second vector in multiplication
C      IW      dimension of vector W (.GE. N)
C      IA      first dimension of dyadic array A (.GE. N)
C      JA      second dimension of A (.GE. N)
C      N       number of elements of vectors V and W to be used
C              in forming the dyad
C      ITEST   error checking option
C
C ARGUMENTS out
C      A       dyadic array formed by product of V and W
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE DYAD(V,IV,W,IW,A,IA,JA,N,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IA,IERROR,ITEST,IV,IW,J,JA,N
      CHARACTER*6 SRNAME
      DOUBLE PRECISION A,V,W
      DIMENSION A(IA,JA),V(IV),W(IW)
      DATA SRNAME/'DYAD'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (N.GT.IW) IERROR = 3
         IF (N.GT.IV) IERROR = 2
         IF (N.LE.0) IERROR = 1
         IF (N.GT.IA .OR. N.GT.JA) IERROR = 4
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      DO 1010 J = 1,N
         DO 1000 I = 1,N
            A(I,J) = V(I)*W(J)
 1000    CONTINUE
 1010 CONTINUE
C
      END
