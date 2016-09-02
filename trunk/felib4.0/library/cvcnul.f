C
      SUBROUTINE CVCNUL(VEC,IVEC,N,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CVCNUL sets the first N elements of a complex vector to zero
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
C      Commented    23 Oct 1985 (CG)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      IVEC    length of vector VEC (.GE.N)
C      N       number of elements to be set to zero
C      ITEST   error checking option
C
C ARGUMENTS out
C      VEC     vector of length IVEC. VEC(1,I)=0.0D0 and VEC(2,I)=0.0D0
C              for I=1(1)N
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE CVCNUL(VEC,IVEC,N,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IERROR,ITEST,IVEC,N
      CHARACTER*6 SRNAME
      DOUBLE PRECISION VEC
      DIMENSION VEC(2,IVEC)
      DATA SRNAME/'CVCNUL'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (N.GT.IVEC) IERROR = 2
         IF (N.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
C
C     Main loops
C
         IF (IERROR.NE.0) RETURN
      END IF
      DO 1000 I = 1,N
         VEC(1,I) = 0.0D0
         VEC(2,I) = 0.0D0
 1000 CONTINUE
C
      END
