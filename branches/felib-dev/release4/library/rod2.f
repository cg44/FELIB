C
      SUBROUTINE ROD2(FUN,IFUN,DER,IDER,XI,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      ROD2 returns the values of the shape functions and their
C      derivatives at a specified point for a 2-noded c0
C      continuous line element.  Only the function will be
C      continuous across element boundaries.
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
C      Commented    26 Feb 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      IFUN    dimension of the vector FUN (.GE. 2)
C      IDER    dimension of the vector DER (.GE. 2)
C      XI      specifies the value of the local coordinate at
C              which the function and its derivative are
C              required
C      ITEST   error checking option
C
C ARGUMENTS out
C      FUN     FUN(I) contains the value of the i'th shape
C              function at XI
C      DER     DER(I) contains the value of the derivative of
C              the i'th shape function at XI
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE ROD2(FUN,IFUN,DER,IDER,XI,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,IDER,IERROR,IFUN,ITEST
      CHARACTER*6 SRNAME
      DOUBLE PRECISION DER,DUMMY,FUN,VEPS,XI,XMAX
      DIMENSION DER(IDER),FUN(IFUN)
      DATA SRNAME/'ROD2'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IFUN.LT.2 .OR. IDER.LT.2) IERROR = 1
         XMAX = 1.0D0 + VEPS(DUMMY)
         IF (DABS(XI).GT.XMAX) IERROR = 2
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      FUN(1) = 0.5D0* (1.0D0-XI)
      FUN(2) = 0.5D0* (1.0D0+XI)
      DER(1) = -0.5D0
      DER(2) = 0.5D0
C
      END
