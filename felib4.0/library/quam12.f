C
      SUBROUTINE QUAM12(FUN,IFUN,DER,IDER,JDER,XI,ETA,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      QUAM12 returns the values of shape functions and their
C      derivatives at a specified point for an 12-noded c0
C      continuous quadrilateral element.  The approximated
C      function will be continuous across element boundaries
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
C      Commented    21 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      IFUN    length of vector FUN (.GE. 12)
C      IDER    first dimension of array DER (.GE. 2)
C      JDER    second dimension of array DER (.GE. 12)
C      XI      first local coordinate
C      ETA     second local coordinate
C      ITEST   error checking option
C
C ARGUMENTS out
C      FUN     vector of length IFUN. FUN(I) contains the
C              value of the i'th shape function at (XI, ETA)
C              for i=1(1)12
C      DER     array of dimension (IDER, JDER). DER(I, J)
C              contains the value of the derivative of the j'th
C              shape function with respect to the i'th
C              coordinate, for i=1(1)2 and j=1(1)12
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE QUAM12(FUN,IFUN,DER,IDER,JDER,XI,ETA,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,IDER,IERROR,IFUN,ITEST,JDER
      CHARACTER*6 SRNAME
      DOUBLE PRECISION DER,DUMMY,ETA,FUN,VAL,VEPS,XI
      DIMENSION DER(IDER,JDER),FUN(IFUN)
      DATA SRNAME/'QUAM12'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IFUN.LT.12) IERROR = 1
         IF (IDER.LT.2 .OR. JDER.LT.12) IERROR = 2
         VAL = 1.0D0 + VEPS(DUMMY)
         IF (DABS(XI).GT.VAL .OR. DABS(ETA).GT.VAL) IERROR = 3
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      FUN(1) = 1.0D0/32.0D0* (1.0D0-XI)* (1.0D0-ETA)*
     *         (-10.0D0+9.0D0* (XI*XI+ETA*ETA))
      FUN(2) = 9.0D0/32.0D0* (1.0D0-XI)* (1.0D0-ETA*ETA)*
     *         (1.0D0-3.0D0*ETA)
      FUN(3) = 9.0D0/32.0D0* (1.0D0-XI)* (1.0D0-ETA*ETA)*
     *         (1.0D0+3.0D0*ETA)
      FUN(4) = 1.0D0/32.0D0* (1.0D0-XI)* (1.0D0+ETA)*
     *         (-10.0D0+9.0D0* (XI*XI+ETA*ETA))
      FUN(5) = 9.0D0/32.0D0* (1.0D0+ETA)* (1.0D0-XI*XI)*
     *         (1.0D0-3.0D0*XI)
      FUN(6) = 9.0D0/32.0D0* (1.0D0+ETA)* (1.0D0-XI*XI)*
     *         (1.0D0+3.0D0*XI)
      FUN(7) = 1.0D0/32.0D0* (1.0D0+XI)* (1.0D0+ETA)*
     *         (-10.0D0+9.0D0* (XI*XI+ETA*ETA))
      FUN(8) = 9.0D0/32.0D0* (1.0D0+XI)* (1.0D0-ETA*ETA)*
     *         (1.0D0+3.0D0*ETA)
      FUN(9) = 9.0D0/32.0D0* (1.0D0+XI)* (1.0D0-ETA*ETA)*
     *         (1.0D0-3.0D0*ETA)
      FUN(10) = 1.0D0/32.0D0* (1.0D0+XI)* (1.0D0-ETA)*
     *          (-10.0D0+9.0D0* (XI*XI+ETA*ETA))
      FUN(11) = 9.0D0/32.0D0* (1.0D0-ETA)* (1.0D0-XI*XI)*
     *          (1.0D0+3.0D0*XI)
      FUN(12) = 9.0D0/32.0D0* (1.0D0-ETA)* (1.0D0-XI*XI)*
     *          (1.0D0-3.0D0*XI)
C
      DER(1,1) = 1.0D0/32.0D0* (1.0D0-ETA)*
     *           (10.D0-27.D0*XI*XI+18.D0*XI-9.D0*ETA*ETA)
      DER(2,1) = 1.D0/32.D0* (1.D0-XI)*
     *           (10.D0-27.D0*ETA*ETA+18.D0*ETA-9.D0*XI*XI)
      DER(1,2) = -9.D0/32.D0* (1.D0-ETA*ETA)* (1.D0-3.D0*ETA)
      DER(2,2) = 9.D0/32.D0* (1.D0-XI)* (9.D0*ETA*ETA-2.D0*ETA-3.D0)
      DER(1,3) = -9.D0/32.D0* (1.D0-ETA*ETA)* (1.D0+3.D0*ETA)
      DER(2,3) = 9.D0/32.D0* (1.D0-XI)* (-9.D0*ETA*ETA-2.D0*ETA+3.D0)
      DER(1,4) = 1.D0/32.D0* (1.D0+ETA)*
     *           (10.D0-27.D0*XI*XI+18.D0*XI-9.D0*ETA*ETA)
      DER(2,4) = 1.D0/32.D0* (1.D0-XI)*
     *           (-10.D0+27.D0*ETA*ETA+18.D0*ETA+9.D0*XI*XI)
      DER(1,5) = 9.D0/32.D0* (1.D0+ETA)* (9.D0*XI*XI-2.D0*XI-3.D0)
      DER(2,5) = 9.D0/32.D0* (1.D0-XI*XI)* (1.D0-3.D0*XI)
      DER(1,6) = 9.D0/32.D0* (1.D0+ETA)* (3.D0-2.D0*XI-9.D0*XI*XI)
      DER(2,6) = 9.D0/32.D0* (1.D0-XI*XI)* (1.D0+3.D0*XI)
      DER(1,7) = 1.D0/32.D0* (1.D0+ETA)*
     *           (-10.D0+27.D0*XI*XI+18.D0*XI+9.D0*ETA*ETA)
      DER(2,7) = 1.D0/32.D0* (1.D0+XI)*
     *           (-10.D0+27.D0*ETA*ETA+18.D0*ETA+9.D0*XI*XI)
      DER(1,8) = 9.D0/32.D0* (1.D0-ETA*ETA)* (1.D0+3.D0*ETA)
      DER(2,8) = 9.D0/32.D0* (1.D0+XI)* (3.D0-2.D0*ETA-9.D0*ETA*ETA)
      DER(1,9) = 9.D0/32.D0* (1.D0-ETA*ETA)* (1.D0-3.D0*ETA)
      DER(2,9) = 9.D0/32.D0* (1.D0+XI)* (9.D0*ETA*ETA-2.D0*ETA-3.D0)
      DER(1,10) = 1.D0/32.D0* (1.D0-ETA)*
     *            (-10.D0+27.D0*XI*XI+18.D0*XI+9.D0*ETA*ETA)
      DER(2,10) = 1.D0/32.D0* (1.D0+XI)*
     *            (10.D0-27.D0*ETA*ETA+18.D0*ETA-9.D0*XI*XI)
      DER(1,11) = 9.D0/32.D0* (1.D0-ETA)* (3.D0-2.D0*XI-9.D0*XI*XI)
      DER(2,11) = -9.D0/32.D0* (1.D0-XI*XI)* (1.D0+3.D0*XI)
      DER(1,12) = 9.D0/32.D0* (1.D0-ETA)* (9.D0*XI*XI-2.D0*XI-3.D0)
      DER(2,12) = -9.D0/32.D0* (1.D0-XI*XI)* (1.D0-3.D0*XI)
C
      END
