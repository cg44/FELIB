C***********************************************************************
C$SPLIT$BRK20$*********************************************************
C***********************************************************************
      SUBROUTINE BRK20(FUN, IFUN, DER, IDER, JDER, XI, ETA, ZETA,
     *     ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CALCULATES SHAPE FUNCTIONS AND DERIVATIVES AT A
C      SPECIFIED POINT FOR 20-NODED BRICK ELEMENT.  FUNCTION
C      CONTINUOUS ACROSS ELEMENT BOUNDARIES
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    10 OCT 1980 (KR)
C
C ARGUMENTS IN
C      IFUN    DIMENSION OF VECTOR FUN (.GE.20)
C      IDER    FIRST DIMENSION OF REAL ARRAY DER (.GE.3)
C      JDER    SECOND DIMENSION OF REAL ARRAY DER (.GE.20)
C      XI      VALUE OF FIRST LOCAL COORDINATE
C      ETA     VALUE OF SECOND LOCAL COORDINATE
C      ZETA    VALUE OF THIRD LOCAL COORDINATE
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      FUN     VECTOR OF SHAPE FUNCTION VALUES.  FUN(I
C              CONTAINS THE VALUE OF THE I'TH SHAPE FUNCTION AT
C              (XI,ETA,ZETA)
C      DER     ARRAY OF SHAPE FUNCTION DERIVATIVE VALUES.
C              DER(I,J) CONTAINS THE VALUE OF THE DERIVATIVE OF
C              THE J'TH SHAPE FUNCTION WITH RESPECT TO THE I'TH
C              COORDINATE AT (XI,ETA,ZETA)
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE BRK20(FUN, IFUN, DER, IDER, JDER, XI, ETA, ZETA,
C    *     ITEST)
C***********************************************************************
C
      INTEGER ERRMES, IDER, IERROR, IFUN, ITEST, JDER
      DOUBLE PRECISION DER, DNCX, DNCY, DNCZ, DNSXX, DNSXY,
     *     DNSXZ, DNSYX, DNSYY, DNSYZ, DNSZX, DNSZY, DNSZZ,
     *     ETA, FUN, M1, NC, NSX, NSY, NSZ, P1, SRNAME, X,
     *     XI, Y, Z, ZETA, VAL, VEPS, DUMMY
      DIMENSION DER(IDER,JDER), FUN(IFUN)
      DATA M1 /-1.D0/, P1 /1.D0/, SRNAME /8H BRK20  /
      NC(X,Y,Z) = 1.D0/8.D0*(1.D0+XI*X)*(1.D0+ETA*Y)*(1.D0+ZETA*Z)*
     *     (XI*X+ETA*Y+ZETA*Z-2.D0)
      NSX(Y,Z) = 1.D0/4.D0*(1.D0-XI*XI)*(1.D0+ETA*Y)*(1.D0+ZETA*Z)
      NSY(X,Z) = 1.D0/4.D0*(1.D0+XI*X)*(1.D0-ETA*ETA)*(1.D0+ZETA*Z)
      NSZ(X,Y) = 1.D0/4.D0*(1.D0+XI*X)*(1.D0+ETA*Y)*(1.D0-ZETA*
     *     ZETA)
      DNCX(X,Y,Z) = 1.D0/8.D0*X*(1.D0+ETA*Y)*(1.D0+ZETA*Z)*(2.D0*
     *     XI*X+ETA*Y+ZETA*Z-1.D0)
      DNCY(X,Y,Z) = 1.D0/8.D0*Y*(1.D0+XI*X)*(1.D0+ZETA*Z)*(XI*
     *     X+2.D0*ETA*Y+ZETA*Z-1.D0)
      DNCZ(X,Y,Z) = 1.D0/8.D0*Z*(1.D0+XI*X)*(1.D0+ETA*Y)*(XI*X+ETA*
     *     Y+2.D0*ZETA*Z-1.D0)
      DNSXX(Y,Z) = -1.D0/2.D0*XI*(1.D0+ETA*Y)*(1.D0+ZETA*Z)
      DNSXY(Y,Z) = 1.D0/4.D0*Y*(1.D0-XI*XI)*(1.D0+ZETA*Z)
      DNSXZ(Y,Z) = 1.D0/4.D0*Z*(1.D0-XI*XI)*(1.D0+ETA*Y)
      DNSYX(X,Z) = 1.D0/4.D0*X*(1.D0-ETA*ETA)*(1.D0+ZETA*Z)
      DNSYY(X,Z) = -1.D0/2.D0*ETA*(1.D0+XI*X)*(1.D0+ZETA*Z)
      DNSYZ(X,Z) = 1.D0/4.D0*Z*(1.D0-ETA*ETA)*(1.D0+XI*X)
      DNSZX(X,Y) = 1.D0/4.D0*X*(1.D0+ETA*Y)*(1.D0-ZETA*ZETA)
      DNSZY(X,Y) = 1.D0/4.D0*Y*(1.D0+XI*X)*(1.D0-ZETA*ZETA)
      DNSZZ(X,Y) = -1.D0/2.D0*ZETA*(1.D0+XI*X)*(1.D0+ETA*Y)
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IFUN.LT.20) IERROR = 1
                        IF (IDER.LT.3 .OR. JDER.LT.20) IERROR = 2
                        VAL = 1.0D0+VEPS(DUMMY)
                        IF (DABS(XI).GT.VAL .OR. DABS(ETA)
     *                     .GT.VAL .OR. DABS(ZETA).GT.VAL)
     *                     IERROR = 3
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 FUN(1) = NC(M1,M1,M1)
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
      RETURN
      END
