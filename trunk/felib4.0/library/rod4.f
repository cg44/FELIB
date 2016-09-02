C
      SUBROUTINE ROD4(FUN,IFUN,DER,IDER,XI,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      ROD4 returns the values of the shape functions and their
C      derivatives at a specified point for a 4-noded c0
C      continuous element. The function continuous across element
C      boundaries
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
C      IFUN    dimension of vector FUN (.GE. 4)
C      IDER    dimension of vector DER (.GE. 4)
C      XI      value of local coordinate
C      ITEST   error checking option
C
C ARGUMENTS out
C      FUN     vector of dimension IFUN. FUN(I) contains the
C              value of the i'th shape function at XI
C      DER     vector of dimension IDER. DER(I) contains the
C              value of the derivative of the i'th shape
C              function at XI
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE ROD4(FUN,IFUN,DER,IDER,XI,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,IDER,IERROR,IFUN,ITEST
      CHARACTER*6 SRNAME
      DOUBLE PRECISION DER,DUMMY,FUN,VEPS,XI,XMAX
      DIMENSION DER(IDER),FUN(IFUN)
      DATA SRNAME/'ROD4'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IFUN.LT.4 .OR. IDER.LT.4) IERROR = 1
         XMAX = 1.0D0 + VEPS(DUMMY)
         IF (DABS(XI).GT.XMAX) IERROR = 2
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      FUN(1) = 1.0D0/16.0D0* (1.0D0-XI)* (9.0D0*XI*XI-1.0D0)
      FUN(2) = 9.0D0/16.0D0* (XI*XI-1.0D0)* (3.0D0*XI-1.0D0)
      FUN(3) = 9.0D0/16.0D0* (1.0D0-XI*XI)* (3.0D0*XI+1.0D0)
      FUN(4) = 1.0D0/16.0D0* (XI+1.0D0)* (9.0D0*XI*XI-1.0D0)
C
      DER(1) = -1.0D0/16.0D0* (-1.0D0-18.0D0*XI+27.0D0*XI*XI)
      DER(2) = 9.0D0/16.0D0* (9.0D0*XI*XI-2.0D0*XI-3.0D0)
      DER(3) = 9.0D0/16.0D0* (3.0D0-2.0D0*XI-9.0D0*XI*XI)
      DER(4) = 1.0D0/16.0D0* (27.0D0*XI*XI+18.0D0*XI-1.0D0)
C
      END
