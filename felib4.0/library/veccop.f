C
      SUBROUTINE VECCOP(V,IV,W,IW,N,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      VECCOP copies first N elements of the vector V into the vector W
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
C      V       vector of length IV to be copied
C      IV      length of vector V (.GE. N)
C      IW      length of vector W
C      N       number of elements of V to be copied
C      ITEST   error checking option
C
C ARGUMENTS out
C      W       vector of length IW; W(I)IS set to V(I) for
C              I=1(1)N
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE VECCOP(V,IV,W,IW,N,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IERROR,ITEST,IV,IW,N
      CHARACTER*6 SRNAME
      DOUBLE PRECISION V,W
      DIMENSION V(IV),W(IW)
C
      DATA SRNAME/'VECCOP'/
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
         W(I) = V(I)
 1000 CONTINUE
C
      END
