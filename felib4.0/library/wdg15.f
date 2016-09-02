C
      SUBROUTINE WDG15(FUN,IFUN,DER,IDER,JDER,XI,ETA,ZETA,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      WDG15 returns the values of shape functions and derivatives
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
C      IFUN    dimension of vector FUN (.GE. 15)
C      IDER    first dimension of array DER (.GE. 3)
C      JDER    second dimension of array DER (.GE. 15)
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
C     SUBROUTINE WDG15(FUN,IFUN,DER,IDER,JDER,XI,ETA,ZETA,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,IDER,IERROR,IFUN,ITEST,JDER
      CHARACTER*6 SRNAME
      DOUBLE PRECISION DER,DL1X,DL1Y,DL2X,DL2Y,DL3X,DL3Y,DNCL,DNCZ,
     *                 DNMQL,DNMQZ,DNMTL,DNMTZ,DUMMY,ETA,FUN,L1,L2,L3,
     *                 LA,LB,NC,NMQ,NMT,VEPS,XI,XMAX,XMIN,YMAX,YMIN,Z,
     *                 ZETA,ZM,ZP,ZVAL
      DIMENSION DER(IDER,JDER),FUN(IFUN)
C
      DATA SRNAME/'WDG15'/
C
C     Statement functions
C
      NC(LA,Z) = 1.D0/2.D0*LA* (2.D0*LA-1.D0)* (1.D0+Z*ZETA) -
     *           1.D0/2.D0*LA* (1.D0-ZETA*ZETA)
      NMT(LA,LB,Z) = 2.D0*LA*LB* (1.D0+Z*ZETA)
      NMQ(LA) = LA* (1.D0-ZETA*ZETA)
      DNCL(LA,Z) = 1.D0/2.D0* ((4.D0*LA-1.D0)* (1.D0+Z*ZETA)+ZETA*ZETA-
     *             1.D0)
      DNCZ(LA,Z) = 1.D0/2.D0*LA* (2.D0*ZETA+Z* (2.D0*LA-1.D0))
      DNMQL(LA) = 1.D0 - ZETA*ZETA
      DNMQZ(LA) = -2.D0*LA*ZETA
      DNMTL(LA,LB,Z) = 2.D0*LB* (1.D0+Z*ZETA)
      DNMTZ(LA,LB,Z) = 2.D0*LA*LB*Z
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         YMIN = 1.0D0/DSQRT(3.0D0)* (XI-1.0D0) - VEPS(DUMMY)
         YMAX = 1.0D0/DSQRT(3.0D0)* (1.0D0-XI) + VEPS(DUMMY)
         XMIN = - (0.5D0+VEPS(DUMMY))
         XMAX = 1.0D0 + VEPS(DUMMY)
         ZVAL = 1.0D0 + VEPS(DUMMY)
         IF ((XI.LT.XMIN.OR.XI.GT.XMAX) .OR.
     *       (ETA.LT.YMIN.OR.ETA.GT.YMAX) .OR.
     *       DABS(ZETA).GT.ZVAL) IERROR = 3
         IF (IDER.LT.3 .OR. JDER.LT.15) IERROR = 2
         IF (IFUN.LT.15) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      L1 = 1.D0/3.D0* (1.D0+2.D0*XI)
      L2 = 1.D0/3.D0* (1.D0-XI-DSQRT(3.D0)*ETA)
      L3 = 1.D0/3.D0* (1.D0-XI+DSQRT(3.D0)*ETA)
      DL1X = 2.D0/3.D0
      DL1Y = 0.D0
      DL2X = -1.D0/3.D0
      DL2Y = -1.D0/DSQRT(3.D0)
      DL3X = -1.D0/3.D0
      DL3Y = 1.D0/DSQRT(3.D0)
      ZP = 1.D0
      ZM = -1.D0
C
C     Shape functions
C
      FUN(1) = NC(L1,ZM)
      FUN(2) = NMT(L1,L2,ZM)
      FUN(3) = NC(L2,ZM)
      FUN(4) = NMT(L2,L3,ZM)
      FUN(5) = NC(L3,ZM)
      FUN(6) = NMT(L3,L1,ZM)
      FUN(7) = NMQ(L1)
      FUN(8) = NMQ(L2)
      FUN(9) = NMQ(L3)
      FUN(10) = NC(L1,ZP)
      FUN(11) = NMT(L1,L2,ZP)
      FUN(12) = NC(L2,ZP)
      FUN(13) = NMT(L2,L3,ZP)
      FUN(14) = NC(L3,ZP)
      FUN(15) = NMT(L3,L1,ZP)
C
C     Derivatives
C
      DER(1,1) = DNCL(L1,ZM)*DL1X
      DER(2,1) = DNCL(L1,ZM)*DL1Y
      DER(3,1) = DNCZ(L1,ZM)
      DER(1,2) = DNMTL(L1,L2,ZM)*DL1X + DNMTL(L2,L1,ZM)*DL2X
      DER(2,2) = DNMTL(L1,L2,ZM)*DL1Y + DNMTL(L2,L1,ZM)*DL2Y
      DER(3,2) = DNMTZ(L1,L2,ZM)
      DER(1,3) = DNCL(L2,ZM)*DL2X
      DER(2,3) = DNCL(L2,ZM)*DL2Y
      DER(3,3) = DNCZ(L2,ZM)
      DER(1,4) = DNMTL(L2,L3,ZM)*DL2X + DNMTL(L3,L2,ZM)*DL3X
      DER(2,4) = DNMTL(L2,L3,ZM)*DL2Y + DNMTL(L3,L2,ZM)*DL3Y
      DER(3,4) = DNMTZ(L2,L3,ZM)
      DER(1,5) = DNCL(L3,ZM)*DL3X
      DER(2,5) = DNCL(L3,ZM)*DL3Y
      DER(3,5) = DNCZ(L3,ZM)
      DER(1,6) = DNMTL(L3,L1,ZM)*DL3X + DNMTL(L1,L3,ZM)*DL1X
      DER(2,6) = DNMTL(L3,L1,ZM)*DL3Y + DNMTL(L1,L3,ZM)*DL1Y
      DER(3,6) = DNMTZ(L3,L1,ZM)
      DER(1,7) = DNMQL(L1)*DL1X
      DER(2,7) = DNMQL(L1)*DL1Y
      DER(3,7) = DNMQZ(L1)
      DER(1,8) = DNMQL(L2)*DL2X
      DER(2,8) = DNMQL(L2)*DL2Y
      DER(3,8) = DNMQZ(L2)
      DER(1,9) = DNMQL(L3)*DL3X
      DER(2,9) = DNMQL(L3)*DL3Y
      DER(3,9) = DNMQZ(L3)
      DER(1,10) = DNCL(L1,ZP)*DL1X
      DER(2,10) = DNCL(L1,ZP)*DL1Y
      DER(3,10) = DNCZ(L1,ZP)
      DER(1,11) = DNMTL(L1,L2,ZP)*DL1X + DNMTL(L2,L1,ZP)*DL2X
      DER(2,11) = DNMTL(L1,L2,ZP)*DL1Y + DNMTL(L2,L1,ZP)*DL2Y
      DER(3,11) = DNMTZ(L1,L2,ZP)
      DER(1,12) = DNCL(L2,ZP)*DL2X
      DER(2,12) = DNCL(L2,ZP)*DL2Y
      DER(3,12) = DNCZ(L2,ZP)
      DER(1,13) = DNMTL(L2,L3,ZP)*DL2X + DNMTL(L3,L2,ZP)*DL3X
      DER(2,13) = DNMTL(L2,L3,ZP)*DL2Y + DNMTL(L3,L2,ZP)*DL3Y
      DER(3,13) = DNMTZ(L2,L3,ZP)
      DER(1,14) = DNCL(L3,ZP)*DL3X
      DER(2,14) = DNCL(L3,ZP)*DL3Y
      DER(3,14) = DNCZ(L3,ZP)
      DER(1,15) = DNMTL(L3,L1,ZP)*DL3X + DNMTL(L1,L3,ZP)*DL1X
      DER(2,15) = DNMTL(L3,L1,ZP)*DL3Y + DNMTL(L1,L3,ZP)*DL1Y
      DER(3,15) = DNMTZ(L3,L1,ZP)
C
      END
