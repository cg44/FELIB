C
      SUBROUTINE BRK20(FUN,IFUN,DER,IDER,JDER,XI,ETA,ZETA,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      BRK20 calculates shape functions and derivatives at a
C      specified point for 20-noded brick element. The shape function
C      is continuous across element boundaries.
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
C      IFUN    dimension of vector FUN (.GE. 20)
C      IDER    first dimension of real array DER (.GE. 3)
C      JDER    second dimension of real array DER (.GE. 20)
C      XI      value of first local coordinate
C      ETA     value of second local coordinate
C      ZETA    value of third local coordinate
C      ITEST   error checking option
C
C ARGUMENTS out
C      FUN     vector of shape function values.  FUN(I)
C              contains the value of the i'th shape function at
C              (XI, ETA, ZETA)
C      DER     array of shape function derivative values.
C              DER(I, J) contains the value of the derivative of
C              the j'th shape function with respect to the i'th
C              coordinate at (XI, ETA, ZETA)
C
C ROUTINES called
C      ERRMES    VEPS
C
C     SUBROUTINE BRK20(FUN,IFUN,DER,IDER,JDER,XI,ETA,ZETA,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,IDER,IERROR,IFUN,ITEST,JDER
      CHARACTER*6 SRNAME
      DOUBLE PRECISION DER,DNCX,DNCY,DNCZ,DNSXX,DNSXY,DNSXZ,DNSYX,DNSYY,
     *                 DNSYZ,DNSZX,DNSZY,DNSZZ,DUMMY,ETA,FUN,M1,NC,NSX,
     *                 NSY,NSZ,P1,VAL,VEPS,X,XI,Y,Z,ZETA
      DIMENSION DER(IDER,JDER),FUN(IFUN)
C
      EXTERNAL ERRMES,VEPS
C
      DATA M1/-1.D0/,P1/1.D0/
      DATA SRNAME/'BRK20 '/
C
C     Statement functions
C
      NC(X,Y,Z) = 1.D0/8.D0* (1.D0+XI*X)* (1.D0+ETA*Y)* (1.D0+ZETA*Z)*
     *            (XI*X+ETA*Y+ZETA*Z-2.D0)
      NSX(Y,Z) = 1.D0/4.D0* (1.D0-XI*XI)* (1.D0+ETA*Y)* (1.D0+ZETA*Z)
      NSY(X,Z) = 1.D0/4.D0* (1.D0+XI*X)* (1.D0-ETA*ETA)* (1.D0+ZETA*Z)
      NSZ(X,Y) = 1.D0/4.D0* (1.D0+XI*X)* (1.D0+ETA*Y)* (1.D0-ZETA*ZETA)
      DNCX(X,Y,Z) = 1.D0/8.D0*X* (1.D0+ETA*Y)* (1.D0+ZETA*Z)*
     *              (2.D0*XI*X+ETA*Y+ZETA*Z-1.D0)
      DNCY(X,Y,Z) = 1.D0/8.D0*Y* (1.D0+XI*X)* (1.D0+ZETA*Z)*
     *              (XI*X+2.D0*ETA*Y+ZETA*Z-1.D0)
      DNCZ(X,Y,Z) = 1.D0/8.D0*Z* (1.D0+XI*X)* (1.D0+ETA*Y)*
     *              (XI*X+ETA*Y+2.D0*ZETA*Z-1.D0)
      DNSXX(Y,Z) = -1.D0/2.D0*XI* (1.D0+ETA*Y)* (1.D0+ZETA*Z)
      DNSXY(Y,Z) = 1.D0/4.D0*Y* (1.D0-XI*XI)* (1.D0+ZETA*Z)
      DNSXZ(Y,Z) = 1.D0/4.D0*Z* (1.D0-XI*XI)* (1.D0+ETA*Y)
      DNSYX(X,Z) = 1.D0/4.D0*X* (1.D0-ETA*ETA)* (1.D0+ZETA*Z)
      DNSYY(X,Z) = -1.D0/2.D0*ETA* (1.D0+XI*X)* (1.D0+ZETA*Z)
      DNSYZ(X,Z) = 1.D0/4.D0*Z* (1.D0-ETA*ETA)* (1.D0+XI*X)
      DNSZX(X,Y) = 1.D0/4.D0*X* (1.D0+ETA*Y)* (1.D0-ZETA*ZETA)
      DNSZY(X,Y) = 1.D0/4.D0*Y* (1.D0+XI*X)* (1.D0-ZETA*ZETA)
      DNSZZ(X,Y) = -1.D0/2.D0*ZETA* (1.D0+XI*X)* (1.D0+ETA*Y)
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IFUN.LT.20) IERROR = 1
         IF (IDER.LT.3 .OR. JDER.LT.20) IERROR = 2
         VAL = 1.0D0 + VEPS(DUMMY)
         IF (DABS(XI).GT.VAL .OR. DABS(ETA).GT.VAL .OR.
     *       DABS(ZETA).GT.VAL) IERROR = 3
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Set shape functions
C
      FUN(1) = NC(M1,M1,M1)
      FUN(2) = NSY(M1,M1)
      FUN(3) = NC(M1,P1,M1)
      FUN(4) = NSX(P1,M1)
      FUN(5) = NC(P1,P1,M1)
      FUN(6) = NSY(P1,M1)
      FUN(7) = NC(P1,M1,M1)
      FUN(8) = NSX(M1,M1)
      FUN(9) = NSZ(M1,M1)
      FUN(10) = NSZ(M1,P1)
      FUN(11) = NSZ(P1,P1)
      FUN(12) = NSZ(P1,M1)
      FUN(13) = NC(M1,M1,P1)
      FUN(14) = NSY(M1,P1)
      FUN(15) = NC(M1,P1,P1)
      FUN(16) = NSX(P1,P1)
      FUN(17) = NC(P1,P1,P1)
      FUN(18) = NSY(P1,P1)
      FUN(19) = NC(P1,M1,P1)
      FUN(20) = NSX(M1,P1)
C
C     Set derivatives
C
      DER(1,1) = DNCX(M1,M1,M1)
      DER(2,1) = DNCY(M1,M1,M1)
      DER(3,1) = DNCZ(M1,M1,M1)
      DER(1,2) = DNSYX(M1,M1)
      DER(2,2) = DNSYY(M1,M1)
      DER(3,2) = DNSYZ(M1,M1)
      DER(1,3) = DNCX(M1,P1,M1)
      DER(2,3) = DNCY(M1,P1,M1)
      DER(3,3) = DNCZ(M1,P1,M1)
      DER(1,4) = DNSXX(P1,M1)
      DER(2,4) = DNSXY(P1,M1)
      DER(3,4) = DNSXZ(P1,M1)
      DER(1,5) = DNCX(P1,P1,M1)
      DER(2,5) = DNCY(P1,P1,M1)
      DER(3,5) = DNCZ(P1,P1,M1)
      DER(1,6) = DNSYX(P1,M1)
      DER(2,6) = DNSYY(P1,M1)
      DER(3,6) = DNSYZ(P1,M1)
      DER(1,7) = DNCX(P1,M1,M1)
      DER(2,7) = DNCY(P1,M1,M1)
      DER(3,7) = DNCZ(P1,M1,M1)
      DER(1,8) = DNSXX(M1,M1)
      DER(2,8) = DNSXY(M1,M1)
      DER(3,8) = DNSXZ(M1,M1)
      DER(1,9) = DNSZX(M1,M1)
      DER(2,9) = DNSZY(M1,M1)
      DER(3,9) = DNSZZ(M1,M1)
      DER(1,10) = DNSZX(M1,P1)
      DER(2,10) = DNSZY(M1,P1)
      DER(3,10) = DNSZZ(M1,P1)
      DER(1,11) = DNSZX(P1,P1)
      DER(2,11) = DNSZY(P1,P1)
      DER(3,11) = DNSZZ(P1,P1)
      DER(1,12) = DNSZX(P1,M1)
      DER(2,12) = DNSZY(P1,M1)
      DER(3,12) = DNSZZ(P1,M1)
      DER(1,13) = DNCX(M1,M1,P1)
      DER(2,13) = DNCY(M1,M1,P1)
      DER(3,13) = DNCZ(M1,M1,P1)
      DER(1,14) = DNSYX(M1,P1)
      DER(2,14) = DNSYY(M1,P1)
      DER(3,14) = DNSYZ(M1,P1)
      DER(1,15) = DNCX(M1,P1,P1)
      DER(2,15) = DNCY(M1,P1,P1)
      DER(3,15) = DNCZ(M1,P1,P1)
      DER(1,16) = DNSXX(P1,P1)
      DER(2,16) = DNSXY(P1,P1)
      DER(3,16) = DNSXZ(P1,P1)
      DER(1,17) = DNCX(P1,P1,P1)
      DER(2,17) = DNCY(P1,P1,P1)
      DER(3,17) = DNCZ(P1,P1,P1)
      DER(1,18) = DNSYX(P1,P1)
      DER(2,18) = DNSYY(P1,P1)
      DER(3,18) = DNSYZ(P1,P1)
      DER(1,19) = DNCX(P1,M1,P1)
      DER(2,19) = DNCY(P1,M1,P1)
      DER(3,19) = DNCZ(P1,M1,P1)
      DER(1,20) = DNSXX(M1,P1)
      DER(2,20) = DNSXY(M1,P1)
      DER(3,20) = DNSXZ(M1,P1)
C
      END
