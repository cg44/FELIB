C
      SUBROUTINE QUAM8(FUN,IFUN,DER,IDER,JDER,XI,ETA,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      QUAM8 returns the values of shape functions and their
C      derivatives at a specified point for an 8-noded c0
C      continuous quadrilateral element. The approximated
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
C      IFUN    length of vector FUN (.GE. 8)
C      IDER    first dimension of array DER (.GE. 2)
C      JDER    second dimension of array DER (.GE. 8)
C      XI      first local coordinate
C      ETA     second local coordinate
C      ITEST   error checking option
C
C ARGUMENTS out
C      FUN     vector of length IFUN. FUN(I) contains the
C              value of the i'th shape function at (XI, ETA)
C              for i=1(1)8
C      DER     array of dimension (IDER, JDER). DER(I, J)
C              contains the value of the derivative of the j'th
C              shape function with respect to the i'th
C              coordinate, for i=1(1)2 and j=1(1)8
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE QUAM8(FUN,IFUN,DER,IDER,JDER,XI,ETA,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,IDER,IERROR,IFUN,ITEST,JDER
      CHARACTER*6 SRNAME
      DOUBLE PRECISION DER,DUMMY,ETA,ETAM,ETAP,FUN,VAL,VEPS,XI,XIM,XIP
      DIMENSION DER(IDER,JDER),FUN(IFUN)
      DATA SRNAME/'QUAM8'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IFUN.LT.8) IERROR = 1
         IF (IDER.LT.2 .OR. JDER.LT.8) IERROR = 2
         VAL = 1.0D0 + VEPS(DUMMY)
         IF (DABS(XI).GT.VAL .OR. DABS(ETA).GT.VAL) IERROR = 3
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      ETAM = 0.250D0* (1.0D0-ETA)
      ETAP = 0.250D0* (1.0D0+ETA)
      XIM = 0.250D0* (1.0D0-XI)
      XIP = 0.250D0* (1.0D0+XI)
C
      DER(1,1) = ETAM* (2.0D0*XI+ETA)
      DER(1,2) = -8.0D0*ETAM*ETAP
      DER(1,3) = ETAP* (2.0D0*XI-ETA)
      DER(1,4) = -4.0D0*ETAP*XI
      DER(1,5) = ETAP* (2.0D0*XI+ETA)
      DER(1,6) = 8.0D0*ETAP*ETAM
      DER(1,7) = ETAM* (2.0D0*XI-ETA)
      DER(1,8) = -4.0D0*ETAM*XI
      DER(2,1) = XIM* (XI+2.0D0*ETA)
      DER(2,2) = -4.0D0*XIM*ETA
      DER(2,3) = XIM* (2.0D0*ETA-XI)
      DER(2,4) = 8.0D0*XIM*XIP
      DER(2,5) = XIP* (XI+2.0D0*ETA)
      DER(2,6) = -4.0D0*XIP*ETA
      DER(2,7) = XIP* (2.0D0*ETA-XI)
      DER(2,8) = -8.0D0*XIM*XIP
C
      FUN(1) = 4.0D0*ETAM*XIM* (-XI-ETA-1.0D0)
      FUN(2) = 32.0D0*XIM*ETAM*ETAP
      FUN(3) = 4.0D0*ETAP*XIM* (-XI+ETA-1.0D0)
      FUN(4) = 32.0D0*XIM*XIP*ETAP
      FUN(5) = 4.0D0*XIP*ETAP* (XI+ETA-1.0D0)
      FUN(6) = 32.0D0*XIP*ETAP*ETAM
      FUN(7) = 4.0D0*XIP*ETAM* (XI-ETA-1.0D0)
      FUN(8) = 32.0D0*XIM*XIP*ETAM
C
      END
