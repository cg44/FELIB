C
      SUBROUTINE WDG6(FUN,IFUN,DER,IDER,JDER,XI,ETA,ZETA,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      WDG6 returns the values of shape functions and derivatives
C      at a specified point for an 6-noded pentahedral element.
C      the function is continuous across element boundaries.
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
C      IFUN    dimension of vector FUN (.GE. 6)
C      IDER    first dimension of array DER (.GE. 3)
C      JDER    second dimension of array DER (.GE. 6)
C      XI      value of local coordinate at which function and
C              derivative values required
C      ETA     value of second local coordinate
C      ZETA    value of third local coordinate
C      ITEST   error checking option
C
C ARGUMENTS out
C      FUN     real vector of dimension IFUN. FUN(I) contains
C              value of i'th shape function at (XI, ETA, ZETA)
C      DER     real array of dimensions (IDER, JDER).  DER(I, J)
C              contains the derivative of the j'th shape
C              function with respect to the i'th coordinate at
C              the point (XI, ETA, ZETA)
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE WDG6(FUN,IFUN,DER,IDER,JDER,XI,ETA,ZETA,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,IDER,IERROR,IFUN,ITEST,JDER
      CHARACTER*6 SRNAME
      DOUBLE PRECISION DER,DUMMY,ETA,FUN,L1,L2,L3,VEPS,XI,XMAX,XMIN,
     *                 YMAX,YMIN,ZETA,ZETAM,ZETAP,ZVAL
      DIMENSION DER(IDER,JDER),FUN(IFUN)
C
      DATA SRNAME/'WDG6'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IFUN.LT.6) IERROR = 1
         IF (IDER.LT.3 .OR. JDER.LT.6) IERROR = 2
         YMIN = 1.0D0/DSQRT(3.0D0)* (XI-1.0D0) - VEPS(DUMMY)
         YMAX = 1.0D0/DSQRT(3.0D0)* (1.0D0-XI) + VEPS(DUMMY)
         XMIN = - (0.5D0+VEPS(DUMMY))
         XMAX = 1.0D0 + VEPS(DUMMY)
         ZVAL = 1.0D0 + VEPS(DUMMY)
         IF ((XI.LT.XMIN.OR.XI.GT.XMAX) .OR.
     *       (ETA.LT.YMIN.OR.ETA.GT.YMAX) .OR.
     *       DABS(ZETA).GT.ZVAL) IERROR = 3
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      L1 = 1.D0/6.D0* (1.D0+2.D0*XI)
      L2 = 1.D0/6.D0* (1.D0-XI-DSQRT(3.D0)*ETA)
      L3 = 1.D0/6.D0* (1.D0-XI+DSQRT(3.D0)*ETA)
      ZETAM = (1.D0-ZETA)
      ZETAP = (1.D0+ZETA)
C
C     Shape functions
C
      FUN(1) = L1*ZETAM
      FUN(2) = L2*ZETAM
      FUN(3) = L3*ZETAM
      FUN(4) = L1*ZETAP
      FUN(5) = L2*ZETAP
      FUN(6) = L3*ZETAP
C
C     Derivatives
C
      DER(1,1) = 1.D0/3.D0*ZETAM
      DER(2,1) = 0.D0
      DER(3,1) = -L1
      DER(1,2) = -1.D0/6.D0*ZETAM
      DER(2,2) = -1.D0/ (2.D0*DSQRT(3.D0))*ZETAM
      DER(3,2) = -L2
      DER(1,3) = -1.D0/6.D0*ZETAM
      DER(2,3) = 1.D0/ (2.D0*DSQRT(3.D0))*ZETAM
      DER(3,3) = -L3
      DER(1,4) = 1.D0/3.D0*ZETAP
      DER(2,4) = 0.D0
      DER(3,4) = L1
      DER(1,5) = -1.D0/6.D0*ZETAP
      DER(2,5) = -1.D0/ (2.D0*DSQRT(3.D0))*ZETAP
      DER(3,5) = L2
      DER(1,6) = -1.D0/6.D0*ZETAP
      DER(2,6) = 1.D0/ (2.D0*DSQRT(3.D0))*ZETAP
      DER(3,6) = L3
C
      END
