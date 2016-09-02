C
      SUBROUTINE BRK32(FUN,IFUN,DER,IDER,JDER,XI,ETA,ZETA,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      BRK32 calculates the values of the shape functions and their
C      derivatives at a point for a 32-noded brick element. The shape
C      function is continuous across element boundaries
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
C      IFUN    dimension of vector FUN (.GE.32)
C      IDER    first dimension of array DER (.GE.3)
C      JDER    second dimension of array DER (.GE.32)
C      XI      first local coordinate
C      ETA     second local coordinate
C      ZETA    third local coordinate
C      ITEST   error checking option
C
C ARGUMENTS out
C      FUN     vector containing shape functions.  FUN(I)
C              contains the value of the i'th shape function at
C              (XI, ETA, ZETA)
C      DER     array containing the derivatives of the shape
C              functions.  DER(I, J) contains the value of the
C              derivative of the j'th shape function with
C              respect to the i'th coordinate at (XI, ETA, ZETA)
C
C ROUTINES called
C      ERRMES
C
C
C     SUBROUTINE BRK32(FUN,IFUN,DER,IDER,JDER,XI,ETA,ZETA,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,IDER,IERROR,IFUN,ITEST,JDER
      CHARACTER*6 SRNAME
      DOUBLE PRECISION DER,DNCX,DNCY,DNCZ,DNSXX,DNSXY,DNSXZ,DNSYX,DNSYY,
     *                 DNSYZ,DNSZX,DNSZY,DNSZZ,DUMMY,ETA,FUN,M1,MTHRD,
     *                 NC,NSX,NSY,NSZ,P1,PTHRD,VAL,VEPS,X,XI,Y,Z,ZETA
      DIMENSION DER(IDER,JDER),FUN(IFUN)
      DATA SRNAME/'BRK32'/
C
C     Statement functions
C
      NC(X,Y,Z) = 1.D0/64.D0* (1.D0+XI*X)* (1.D0+ETA*Y)* (1.D0+ZETA*Z)*
     *             (9.D0* (XI*XI+ETA*ETA+ZETA*ZETA)-19.D0)
      NSX(X,Y,Z) = 9.D0/64.D0* (1.D0-XI*XI)* (1.D0+9.D0*XI*X)*
     *             (1.D0+ETA*Y)* (1.D0+ZETA*Z)
      NSY(X,Y,Z) = 9.D0/64.D0* (1.D0-ETA*ETA)* (1.D0+9.D0*ETA*Y)*
     *             (1.D0+XI*X)* (1.D0+ZETA*Z)
      NSZ(X,Y,Z) = 9.D0/64.D0* (1.D0-ZETA*ZETA)* (1.D0+9.D0*ZETA*Z)*
     *             (1.D0+XI*X)* (1.D0+ETA*Y)
      DNCX(X,Y,Z) = 1.D0/64.D0* (1.D0+ETA*Y)* (1.D0+ZETA*Z)*
     *              (X* (9.D0* (XI*XI+ETA*ETA+ZETA*ZETA)-19.D0)+
     *              18.D0*XI* (1.D0+XI*X))
      DNCY(X,Y,Z) = 1.D0/64.D0* (1.D0+XI*X)* (1.D0+ZETA*Z)*
     *              (Y* (9.D0* (XI*XI+ETA*ETA+ZETA*ZETA)-19.D0)+
     *              18.D0*ETA* (1.D0+ETA*Y))
      DNCZ(X,Y,Z) = 1.D0/64.D0* (1.D0+XI*X)* (1.D0+ETA*Y)*
     *              (Z* (9.D0* (XI*XI+ETA*ETA+ZETA*ZETA)-19.D0)+
     *              18.D0*ZETA* (1.D0+ZETA*Z))
      DNSXX(X,Y,Z) = 9.D0/64.D0* (1.D0+ETA*Y)* (1.D0+ZETA*Z)*
     *               (9.D0*X* (1.D0-XI*XI)-2.D0*XI* (1.D0+9.D0*XI*X))
      DNSXY(X,Y,Z) = 9.D0/64.D0*Y* (1.D0-XI*XI)* (1.D0+9.D0*XI*X)*
     *               (1.D0+ZETA*Z)
      DNSXZ(X,Y,Z) = 9.D0/64.D0*Z* (1.D0-XI*XI)* (1.D0+9.D0*XI*X)*
     *               (1.D0+ETA*Y)
      DNSYX(X,Y,Z) = 9.D0/64.D0*X* (1.D0-ETA*ETA)* (1.D0+9.D0*ETA*Y)*
     *               (1.D0+ZETA*Z)
      DNSYY(X,Y,Z) = 9.D0/64.D0* (1.D0+XI*X)* (1.D0+ZETA*Z)*
     *               (9.D0*Y* (1.D0-ETA*ETA)-2.D0*ETA*
     *               (1.D0+9.D0*ETA*Y))
      DNSYZ(X,Y,Z) = 9.D0/64.D0*Z* (1.D0-ETA*ETA)* (1.D0+9.D0*ETA*Y)*
     *               (1.D0+XI*X)
      DNSZX(X,Y,Z) = 9.D0/64.D0*X* (1.D0-ZETA*ZETA)* (1.D0+9.D0*ZETA*Z)*
     *                (1.D0+ETA*Y)
      DNSZY(X,Y,Z) = 9.D0/64.D0*Y* (1.D0-ZETA*ZETA)* (1.D0+9.D0*ZETA*Z)*
     *                (1.D0+XI*X)
      DNSZZ(X,Y,Z) = 9.D0/64.D0* (1.D0+XI*X)* (1.D0+ETA*Y)*
     *               (9.D0*Z* (1.D0-ZETA*ZETA)-
     *               2.D0*ZETA* (1.D0+9.D0*ZETA*Z))
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IFUN.LT.32) IERROR = 1
         IF (IDER.LT.3 .OR. JDER.LT.32) IERROR = 2
         VAL = 1.0D0 + VEPS(DUMMY)
         IF (DABS(XI).GT.VAL .OR. DABS(ETA).GT.VAL .OR.
     *       DABS(ZETA).GT.VAL) IERROR = 3
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Body of code
C
      P1 = 1.D0
      M1 = -1.D0
      PTHRD = 1.D0/3.D0
      MTHRD = -1.D0/3.D0
C
C     Set shape functions
C
      FUN(1) = NC(M1,M1,M1)
      FUN(2) = NSY(M1,MTHRD,M1)
      FUN(3) = NSY(M1,PTHRD,M1)
      FUN(4) = NC(M1,P1,M1)
      FUN(5) = NSX(MTHRD,P1,M1)
      FUN(6) = NSX(PTHRD,P1,M1)
      FUN(7) = NC(P1,P1,M1)
      FUN(8) = NSY(P1,PTHRD,M1)
      FUN(9) = NSY(P1,MTHRD,M1)
      FUN(10) = NC(P1,M1,M1)
      FUN(11) = NSX(PTHRD,M1,M1)
      FUN(12) = NSX(MTHRD,M1,M1)
      FUN(13) = NSZ(M1,M1,MTHRD)
      FUN(14) = NSZ(M1,P1,MTHRD)
      FUN(15) = NSZ(P1,P1,MTHRD)
      FUN(16) = NSZ(P1,M1,MTHRD)
      FUN(17) = NSZ(M1,M1,PTHRD)
      FUN(18) = NSZ(M1,P1,PTHRD)
      FUN(19) = NSZ(P1,P1,PTHRD)
      FUN(20) = NSZ(P1,M1,PTHRD)
      FUN(21) = NC(M1,M1,P1)
      FUN(22) = NSY(M1,MTHRD,P1)
      FUN(23) = NSY(M1,PTHRD,P1)
      FUN(24) = NC(M1,P1,P1)
      FUN(25) = NSX(MTHRD,P1,P1)
      FUN(26) = NSX(PTHRD,P1,P1)
      FUN(27) = NC(P1,P1,P1)
      FUN(28) = NSY(P1,PTHRD,P1)
      FUN(29) = NSY(P1,MTHRD,P1)
      FUN(30) = NC(P1,M1,P1)
      FUN(31) = NSX(PTHRD,M1,P1)
      FUN(32) = NSX(MTHRD,M1,P1)
C
C     Set derivatives
C
      DER(1,1) = DNCX(M1,M1,M1)
      DER(2,1) = DNCY(M1,M1,M1)
      DER(3,1) = DNCZ(M1,M1,M1)
      DER(1,2) = DNSYX(M1,MTHRD,M1)
      DER(2,2) = DNSYY(M1,MTHRD,M1)
      DER(3,2) = DNSYZ(M1,MTHRD,M1)
      DER(1,3) = DNSYX(M1,PTHRD,M1)
      DER(2,3) = DNSYY(M1,PTHRD,M1)
      DER(3,3) = DNSYZ(M1,PTHRD,M1)
      DER(1,4) = DNCX(M1,P1,M1)
      DER(2,4) = DNCY(M1,P1,M1)
      DER(3,4) = DNCZ(M1,P1,M1)
      DER(1,5) = DNSXX(MTHRD,P1,M1)
      DER(2,5) = DNSXY(MTHRD,P1,M1)
      DER(3,5) = DNSXZ(MTHRD,P1,M1)
      DER(1,6) = DNSXX(PTHRD,P1,M1)
      DER(2,6) = DNSXY(PTHRD,P1,M1)
      DER(3,6) = DNSXZ(PTHRD,P1,M1)
      DER(1,7) = DNCX(P1,P1,M1)
      DER(2,7) = DNCY(P1,P1,M1)
      DER(3,7) = DNCZ(P1,P1,M1)
      DER(1,8) = DNSYX(P1,PTHRD,M1)
      DER(2,8) = DNSYY(P1,PTHRD,M1)
      DER(3,8) = DNSYZ(P1,PTHRD,M1)
      DER(1,9) = DNSYX(P1,MTHRD,M1)
      DER(2,9) = DNSYY(P1,MTHRD,M1)
      DER(3,9) = DNSYZ(P1,MTHRD,M1)
      DER(1,10) = DNCX(P1,M1,M1)
      DER(2,10) = DNCY(P1,M1,M1)
      DER(3,10) = DNCZ(P1,M1,M1)
      DER(1,11) = DNSXX(PTHRD,M1,M1)
      DER(2,11) = DNSXY(PTHRD,M1,M1)
      DER(3,11) = DNSXZ(PTHRD,M1,M1)
      DER(1,12) = DNSXX(MTHRD,M1,M1)
      DER(2,12) = DNSXY(MTHRD,M1,M1)
      DER(3,12) = DNSXZ(MTHRD,M1,M1)
      DER(1,13) = DNSZX(M1,M1,MTHRD)
      DER(2,13) = DNSZY(M1,M1,MTHRD)
      DER(3,13) = DNSZZ(M1,M1,MTHRD)
      DER(1,14) = DNSZX(M1,P1,MTHRD)
      DER(2,14) = DNSZY(M1,P1,MTHRD)
      DER(3,14) = DNSZZ(M1,P1,MTHRD)
      DER(1,15) = DNSZX(P1,P1,MTHRD)
      DER(2,15) = DNSZY(P1,P1,MTHRD)
      DER(3,15) = DNSZZ(P1,P1,MTHRD)
      DER(1,16) = DNSZX(P1,M1,MTHRD)
      DER(2,16) = DNSZY(P1,M1,MTHRD)
      DER(3,16) = DNSZZ(P1,M1,MTHRD)
      DER(1,17) = DNSZX(M1,M1,PTHRD)
      DER(2,17) = DNSZY(M1,M1,PTHRD)
      DER(3,17) = DNSZZ(M1,M1,PTHRD)
      DER(1,18) = DNSZX(M1,P1,PTHRD)
      DER(2,18) = DNSZY(M1,P1,PTHRD)
      DER(3,18) = DNSZZ(M1,P1,PTHRD)
      DER(1,19) = DNSZX(P1,P1,PTHRD)
      DER(2,19) = DNSZY(P1,P1,PTHRD)
      DER(3,19) = DNSZZ(P1,P1,PTHRD)
      DER(1,20) = DNSZX(P1,M1,PTHRD)
      DER(2,20) = DNSZY(P1,M1,PTHRD)
      DER(3,20) = DNSZZ(P1,M1,PTHRD)
      DER(1,21) = DNCX(M1,M1,P1)
      DER(2,21) = DNCY(M1,M1,P1)
      DER(3,21) = DNCZ(M1,M1,P1)
      DER(1,22) = DNSYX(M1,MTHRD,P1)
      DER(2,22) = DNSYY(M1,MTHRD,P1)
      DER(3,22) = DNSYZ(M1,MTHRD,P1)
      DER(1,23) = DNSYX(M1,PTHRD,P1)
      DER(2,23) = DNSYY(M1,PTHRD,P1)
      DER(3,23) = DNSYZ(M1,PTHRD,P1)
      DER(1,24) = DNCX(M1,P1,P1)
      DER(2,24) = DNCY(M1,P1,P1)
      DER(3,24) = DNCZ(M1,P1,P1)
      DER(1,25) = DNSXX(MTHRD,P1,P1)
      DER(2,25) = DNSXY(MTHRD,P1,P1)
      DER(3,25) = DNSXZ(MTHRD,P1,P1)
      DER(1,26) = DNSXX(PTHRD,P1,P1)
      DER(2,26) = DNSXY(PTHRD,P1,P1)
      DER(3,26) = DNSXZ(PTHRD,P1,P1)
      DER(1,27) = DNCX(P1,P1,P1)
      DER(2,27) = DNCY(P1,P1,P1)
      DER(3,27) = DNCZ(P1,P1,P1)
      DER(1,28) = DNSYX(P1,PTHRD,P1)
      DER(2,28) = DNSYY(P1,PTHRD,P1)
      DER(3,28) = DNSYZ(P1,PTHRD,P1)
      DER(1,29) = DNSYX(P1,MTHRD,P1)
      DER(2,29) = DNSYY(P1,MTHRD,P1)
      DER(3,29) = DNSYZ(P1,MTHRD,P1)
      DER(1,30) = DNCX(P1,M1,P1)
      DER(2,30) = DNCY(P1,M1,P1)
      DER(3,30) = DNCZ(P1,M1,P1)
      DER(1,31) = DNSXX(PTHRD,M1,P1)
      DER(2,31) = DNSXY(PTHRD,M1,P1)
      DER(3,31) = DNSXZ(PTHRD,M1,P1)
      DER(1,32) = DNSXX(MTHRD,M1,P1)
      DER(2,32) = DNSXY(MTHRD,M1,P1)
      DER(3,32) = DNSXZ(MTHRD,M1,P1)
C
      END
