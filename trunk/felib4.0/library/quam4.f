C
      SUBROUTINE QUAM4(FUN,IFUN,DER,IDER,JDER,XI,ETA,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      QUAM4 returns the values of shape functions and their
C      derivatives at a specified point for a 4-noded c0
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
C      IFUN    length of vector FUN (.GE. 4)
C      IDER    first dimension of array DER (.GE. 2)
C      JDER    second dimension of array DER (.GE. 4)
C      XI      first local coordinate value
C      ETA     second local coordinate value
C      ITEST   error checking option
C
C ARGUMENTS out
C      FUN     vector of length IFUN. FUN(I) contains the
C              value of the i'th shape function at the point
C              (XI, ETA)
C      DER     array of dimension (IDER, JDER). DER(I, J)
C              contains the value of the derivative of the j'th
C              shape function with respect to the i'th
C              coordinate at the point (XI, ETA)
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE QUAM4(FUN,IFUN,DER,IDER,JDER,XI,ETA,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,IDER,IERROR,IFUN,ITEST,JDER
      CHARACTER*6 SRNAME
      DOUBLE PRECISION DER,DUMMY,ETA,ETAM,ETAP,FUN,VAL,VEPS,XI,XIM,XIP
      DIMENSION DER(IDER,JDER),FUN(IFUN)
      DATA SRNAME/'QUAM4'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IFUN.LT.4) IERROR = 1
         IF (IDER.LT.2 .OR. JDER.LT.4) IERROR = 2
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
      FUN(1) = 4.0D0*XIM*ETAM
      FUN(2) = 4.0D0*XIM*ETAP
      FUN(3) = 4.0D0*XIP*ETAP
      FUN(4) = 4.0D0*XIP*ETAM
C
      DER(1,1) = -ETAM
      DER(2,1) = -XIM
      DER(1,2) = -ETAP
      DER(2,2) = XIM
      DER(1,3) = ETAP
      DER(2,3) = XIP
      DER(1,4) = ETAM
      DER(2,4) = -XIP
C
      END
