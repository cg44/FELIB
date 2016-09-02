C
      SUBROUTINE BRK8(FUN,IFUN,DER,IDER,JDER,XI,ETA,ZETA,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      BRK8 returns the values of shape functions and derivatives
C      at a specified point for an 8-noded brick element. The shape
C      function is continuous across element boundaries.
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
C      IFUN    dimension of vector FUN (.GE.8)
C      IDER    first dimension of array DER (.GE.3)
C      JDER    second dimension of array DER (.GE.8)
C      XI      value of local coordinate at which function and
C              derivative values required
C      ETA     value of second local coordinate
C      ZETA    value of third local coordinate
C      ITEST   error checking option
C
C ARGUMENTS out
C      FUN     real vector of dimension IFUN. FUN(I) contains
C              value of i'th shape function at (XI, ETA, ZETA)
C      DER     real array of dimensions (IDER, JDER). DER(I, J)
C              contains the derivative of the j'th shape
C              function with respect to the i'th coordinate at
C              the point (XI, ETA, ZETA)
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE BRK8(FUN,IFUN,DER,IDER,JDER,XI,ETA,ZETA,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,IDER,IERROR,IFUN,ITEST,JDER
      CHARACTER*6 SRNAME
      DOUBLE PRECISION DER,DUMMY,ETA,ETAM,ETAP,FUN,HALF,VAL,VEPS,XI,XIM,
     *                 XIP,ZETA,ZETAM,ZETAP
      DIMENSION DER(IDER,JDER),FUN(IFUN)
      DATA SRNAME/'BRK8'/,HALF/0.5D0/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IFUN.LT.8) IERROR = 1
         IF (IDER.LT.3 .OR. JDER.LT.8) IERROR = 2
         VAL = 1.0D0 + VEPS(DUMMY)
         IF (DABS(XI).GT.VAL .OR. DABS(ETA).GT.VAL .OR.
     *       DABS(ZETA).GT.VAL) IERROR = 3
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main code
C
      ETAM = (1.0D0-ETA)*HALF
      ETAP = (1.0D0+ETA)*HALF
      XIM = (1.0D0-XI)*HALF
      XIP = (1.0D0+XI)*HALF
      ZETAM = (1.0D0-ZETA)*HALF
      ZETAP = (1.0D0+ZETA)*HALF
C
C     Set shape functions
C
      FUN(1) = XIM*ETAM*ZETAM
      FUN(2) = XIM*ETAP*ZETAM
      FUN(3) = XIP*ETAP*ZETAM
      FUN(4) = XIP*ETAM*ZETAM
      FUN(5) = XIM*ETAM*ZETAP
      FUN(6) = XIM*ETAP*ZETAP
      FUN(7) = XIP*ETAP*ZETAP
      FUN(8) = XIP*ETAM*ZETAP
C
C     Set derivatves of shape functions
C
      DER(1,1) = -ETAM*ZETAM*HALF
      DER(1,2) = -ETAP*ZETAM*HALF
      DER(1,3) = ETAP*ZETAM*HALF
      DER(1,4) = ETAM*ZETAM*HALF
      DER(1,5) = -ETAM*ZETAP*HALF
      DER(1,6) = -ETAP*ZETAP*HALF
      DER(1,7) = ETAP*ZETAP*HALF
      DER(1,8) = ETAM*ZETAP*HALF
      DER(2,1) = -XIM*ZETAM*HALF
      DER(2,2) = XIM*ZETAM*HALF
      DER(2,3) = XIP*ZETAM*HALF
      DER(2,4) = -XIP*ZETAM*HALF
      DER(2,5) = -XIM*ZETAP*HALF
      DER(2,6) = XIM*ZETAP*HALF
      DER(2,7) = XIP*ZETAP*HALF
      DER(2,8) = -XIP*ZETAP*HALF
      DER(3,1) = -XIM*ETAM*HALF
      DER(3,2) = -XIM*ETAP*HALF
      DER(3,3) = -XIP*ETAP*HALF
      DER(3,4) = -XIP*ETAM*HALF
      DER(3,5) = XIM*ETAM*HALF
      DER(3,6) = XIM*ETAP*HALF
      DER(3,7) = XIP*ETAP*HALF
      DER(3,8) = XIP*ETAM*HALF
C
      END
