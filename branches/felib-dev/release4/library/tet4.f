C
      SUBROUTINE TET4(FUN,IFUN,DER,IDER,JDER,XI,ETA,ZETA,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      TET4 returns the values of shape functions and derivatives
C      at a specified point for an 10-noded tetrahedral element.
C      The function is continuous across element boundaries.
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
C      Commented    10 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      IFUN    dimension of vector FUN (.GE. 4)
C      IDER    first dimension of array DER (.GE. 3)
C      JDER    second dimension of array DER (.GE. 4)
C      XI      value of local coordinate at which function and
C              derivative values required
C      ETA     value of second local coordinate
C      ZETA    value of third local coordinate
C      ITEST   error checking option
C
C ARGUMENTS out
C      FUN     real vector of dimension IFUN.  FUN(I) contains
C              value of i'th shape function at (XI, ETA, ZETA)
C      DER     real array of dimensions (IDER, JDER).  DER(I, J)
C              contains the derivative of the j'th shape
C              function with respect to the i'th coordinate at
C              the point (XI, ETA, ZETA)
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE TET4(FUN,IFUN,DER,IDER,JDER,XI,ETA,ZETA,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,IDER,IERROR,IFUN,ITEST,JDER
      DOUBLE PRECISION DER,ETA,FUN,XI,ZETA
      CHARACTER*6 SRNAME
      DIMENSION DER(IDER,JDER),FUN(IFUN)
      DATA SRNAME/'TET4  '/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IFUN.LT.4) IERROR = 1
         IF (IDER.LT.3 .OR. JDER.LT.4) IERROR = 2
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      FUN(1) = 1.D0/12.D0* (3.D0+8.D0*XI-2.D0*DSQRT(2.D0)*ZETA)
      FUN(2) = 1.D0/12.D0* (3.D0-4.D0* (XI+DSQRT(3.D0)*ETA)-
     *         2.D0*DSQRT(2.D0)*ZETA)
      FUN(3) = 1.D0/12.D0* (3.D0-4.D0* (XI-DSQRT(3.D0)*ETA)-
     *         2.D0*DSQRT(2.D0)*ZETA)
      FUN(4) = 1.D0/4.D0* (1.D0+2.D0*DSQRT(2.D0)*ZETA)
C
      DER(1,1) = 2.D0/3.D0
      DER(2,1) = 0.D0
      DER(3,1) = -DSQRT(2.D0)/6.D0
      DER(1,2) = -1.D0/3.D0
      DER(2,2) = -DSQRT(3.D0)/3.D0
      DER(3,2) = -DSQRT(2.D0)/6.D0
      DER(1,3) = -1.D0/3.D0
      DER(2,3) = DSQRT(3.D0)/3.D0
      DER(3,3) = -DSQRT(2.D0)/6.D0
      DER(1,4) = 0.D0
      DER(2,4) = 0.D0
      DER(3,4) = DSQRT(2.D0)/2.D0
C
      END
