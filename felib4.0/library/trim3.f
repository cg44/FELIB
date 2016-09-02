C
      SUBROUTINE TRIM3(FUN,IFUN,DER,IDER,JDER,XI,ETA,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      TRIM3 returns the values of the shape functions and their
C      derivatives at a specified point for a three-noded c0
C      continuous triangular element. The approximation function
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
C      Commented    22 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      IFUN    length of vector FUN (.GE. 3)
C      IDER    first dimension of array DER (.GE. 2)
C      JDER    second dimension of array DER (.GE. 3)
C      XI      first local coordinate
C      ETA     second local coordinate
C      ITEST   error checking option
C
C ARGUMENTS out
C      FUN     vector of length IFUN. FUN(I) contains the
C              value of the i'th shape function at the point
C              (XI,ETA), for i=1(1)3
C      DER     array of dimension (IDER, JDER).  DER(I, J)
C              contains the derivative of the j'th shape
C              function with respect to the i'th coordinate for
C              i=1(1)2 and j=1(1)3
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE TRIM3(FUN,IFUN,DER,IDER,JDER,XI,ETA,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,IDER,IERROR,IFUN,ITEST,JDER
      CHARACTER*6 SRNAME
      DOUBLE PRECISION DER,DUMMY,ETA,FUN,VEPS,XI,XMAX,XMIN,YMAX,YMIN
      DIMENSION DER(IDER,JDER),FUN(IFUN)
C
      DATA SRNAME/'TRIM3'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IFUN.LT.3) IERROR = 1
         IF (IDER.LT.2 .OR. JDER.LT.3) IERROR = 2
         YMIN = 1.0D0/DSQRT(3.0D0)* (XI-1.0D0) - VEPS(DUMMY)
         YMAX = 1.0D0/DSQRT(3.0D0)* (1.0D0-XI) + VEPS(DUMMY)
         XMIN = - (0.5D0+VEPS(DUMMY))
         XMAX = 1.0D0 + VEPS(DUMMY)
         IF ((XI.LT.XMIN.OR.XI.GT.XMAX) .OR.
     *       (ETA.LT.YMIN.OR.ETA.GT.YMAX)) IERROR = 3
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      FUN(1) = 1.0D0/3.0D0* (1.0D0+2.0D0*XI)
      FUN(2) = 1.0D0/3.0D0* (1.0D0-XI-DSQRT(3.0D0)*ETA)
      FUN(3) = 1.0D0/3.0D0* (1.0D0-XI+DSQRT(3.0D0)*ETA)
C
C     Derivatives
C
      DER(1,1) = 2.0D0/3.0D0
      DER(1,2) = -1.0D0/3.0D0
      DER(1,3) = -1.0D0/3.0D0
      DER(2,1) = 0.0D0
      DER(2,2) = -1.0D0/DSQRT(3.0D0)
      DER(2,3) = +1.0D0/DSQRT(3.0D0)
C
      END
