C
      SUBROUTINE VECNUL(VEC,IVEC,N,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      VECNUL sets the first N elements of a vector to zero.
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
C      IVEC    length of vector VEC (.GE. N)
C      N       number of elements to be set to zero
C      ITEST   error checking option
C
C ARGUMENTS out
C      VEC     vector of length IVEC. VEC(I)=0.0d0 for I=1(1)N
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE VECNUL(VEC,IVEC,N,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IERROR,ITEST,IVEC,N
      CHARACTER*6 SRNAME
      DOUBLE PRECISION VEC
      DIMENSION VEC(IVEC)
C
      DATA SRNAME/'VECNUL'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (N.GT.IVEC) IERROR = 2
         IF (N.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      DO 1000 I = 1,N
         VEC(I) = 0.0D0
 1000 CONTINUE
C
      END
