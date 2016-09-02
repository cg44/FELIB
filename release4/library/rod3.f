C
      SUBROUTINE ROD3(FUN,IFUN,DER,IDER,XI,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      ROD3 returns the values of the shape function and its
C      derivative at a specified point for a 3-noded c0
C      continuous line element. The function continuous across
C      element boundaries
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
C      Commented    11 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      IFUN    dimension of vector FUN (.GE. 3)
C      IDER    dimension of vector DER (.GE. 3)
C      XI      value of local coordinate
C      ITEST   error checking option
C
C ARGUMENTS out
C      FUN     vector of length IFUN.  FUN(I) contains the
C              value of the i'th shape function at the point XI
C      DER     vector of length IDER.  DER(I) contains the
C              value of the derivative of the i'th shape
C              function at the point XI
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE ROD3(FUN,IFUN,DER,IDER,XI,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,IDER,IERROR,IFUN,ITEST
      CHARACTER*6 SRNAME
      DOUBLE PRECISION DER,DUMMY,FUN,VEPS,XI,XMAX
      DIMENSION DER(IDER),FUN(IFUN)
      DATA SRNAME/'ROD3'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IFUN.LT.3 .OR. IDER.LT.3) IERROR = 1
         XMAX = 1.0D0 + VEPS(DUMMY)
         IF (DABS(XI).GT.XMAX) IERROR = 2
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      FUN(1) = 0.5D0*XI* (XI-1.0D0)
      FUN(2) = 1.0D0 - XI*XI
      FUN(3) = 0.5D0*XI* (XI+1.0D0)
C
      DER(1) = 0.5D0* (2.0D0*XI-1.0D0)
      DER(2) = -2.0D0*XI
      DER(3) = 0.5D0* (2.0D0*XI+1.0D0)
C
      END
