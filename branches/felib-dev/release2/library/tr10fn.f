C***********************************************************************
C$SPLIT$TR10FN$*********************************************************
C***********************************************************************
      SUBROUTINE TR10FN(FUN, IFUN, DER, IDER, JDER, XI, ETA, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS THE VALUES OF THE SHAPE FUNCTIONS AND THEIR
C      DERIVATIVES AT A SPECIFIED POINT FOR A 10-NODED C0
C      CONTINUOUS TRIANGULAR ELEMENT.  APPROXIMATED FUNCTION
C      CONTINUOUS ACROSS ELEMENT BOUNDARIES
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    22 OCT 1980 (KR)
C
C ARGUMENTS IN
C      IFUN    LENGTH OF VECTOR FUN (.GE.10)
C      IDER    FIRST DIMENSION OF ARRAY DER (.GE.2)
C      JDER    SECOND DIMENSION OF ARRAY DER (.GE.10)
C      XI      FIRST LOCAL COORDINATE
C      ETA     SECOND LOCAL COORDINATE
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      FUN     VECTOR OF LENGTH IFUN.  FUN(I) CONTAINS THE
C              VALUE OF THE I'TH SHAPE FUNCTION AT THE POINT
C              (XI,ETA), FOR I=1(1)10
C      DER     ARRAY OF DIMENSION (IDER,JDER).  DER(I,J)
C              CONTAINS THE VALUE OF THE DERIVATIVE OF THE J'TH
C              SHAPE FUNCTION WITH RESPECT TO THE I'TH LOCAL
C              COORDINATE, FOR I=1(1)2 AND J=1(1)10
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE TR10FN(FUN, IFUN, DER, IDER, JDER, XI, ETA, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, IDER, IERROR, IFUN, ITEST, JDER
      DOUBLE PRECISION DER, DLX, DLY, DNS1, DNS2, DNVL, ETA,
     *     FUN, L1, L2, L3, LA, LB, NS, NV, SRNAME, XI, YMIN,
     *     YMAX, XMIN, XMAX, VEPS, DUMMY
      DIMENSION DER(IDER,JDER), DLX(3), DLY(3), FUN(IFUN)
      DATA SRNAME /8H TR10FN /
      NV(LA) = 1.D0/2.D0*(3.D0*LA-1.D0)*(3.D0*LA-2.D0)*LA
      NS(LA,LB) = 9.D0/2.D0*LA*LB*(3.D0*LA-1.D0)
      DNVL(LA) = 1.D0/2.D0*(27.D0*LA*LA-18.D0*LA+2.D0)
      DNS1(LA,LB) = 9.D0/2.D0*LB*(6.D0*LA-1.D0)
      DNS2(LA,LB) = 9.D0/2.D0*LA*(3.D0*LA-1.D0)
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IFUN.LT.10) IERROR = 1
                        IF (IDER.LT.2 .OR. JDER.LT.10) IERROR = 2
                        YMIN = 1.0D0/DSQRT(3.0D0)*(XI-1.0D0)-
     *                       VEPS(DUMMY)
                        YMAX = 1.0D0/DSQRT(3.0D0)*(1.0D0-XI)+
     *                       VEPS(DUMMY)
                        XMIN=-(0.5D0+VEPS(DUMMY))
                        XMAX=1.0D0+VEPS(DUMMY)
                        IF ((XI.LT.XMIN .OR. XI.GT.XMAX)
     *                      .OR. (ETA.LT.YMIN .OR.
     *                       ETA.GT.YMAX)) IERROR = 3
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 L1 = 1.0D0/3.0D0*(1.0D0+2.0D0*XI)
      L2 = 1.0D0/3.0D0*(1.0D0-XI-DSQRT(3.0D0)*ETA)
      L3 = 1.0D0/3.0D0*(1.0D0-XI+DSQRT(3.0D0)*ETA)
      DLX(1) = 2.D0/3.D0
      DLY(1) = 0.D0
      DLX(2) = -1.D0/3.D0
      DLY(2) = -1.D0/DSQRT(3.D0)
      DLX(3) = -1.D0/3.D0
      DLY(3) = 1.D0/DSQRT(3.D0)
      FUN(1) = NV(L1)
      FUN(2) = NS(L1,L2)
      FUN(3) = NS(L2,L1)
      FUN(4) = NV(L2)
      FUN(5) = NS(L2,L3)
      FUN(6) = NS(L3,L2)
      FUN(7) = NV(L3)
      FUN(8) = NS(L3,L1)
      FUN(9) = NS(L1,L3)
      FUN(10) = 27.D0*L1*L2*L3
      DER(1,1) = DNVL(L1)*DLX(1)
      DER(2,1) = DNVL(L1)*DLY(1)
      DER(1,2) = DNS1(L1,L2)*DLX(1) + DNS2(L1,L2)*DLX(2)
      DER(2,2) = DNS1(L1,L2)*DLY(1) + DNS2(L1,L2)*DLY(2)
      DER(1,3) = DNS1(L2,L1)*DLX(2) + DNS2(L2,L1)*DLX(1)
      DER(2,3) = DNS1(L2,L1)*DLY(2) + DNS2(L2,L1)*DLY(1)
      DER(1,4) = DNVL(L2)*DLX(2)
      DER(2,4) = DNVL(L2)*DLY(2)
      DER(1,5) = DNS1(L2,L3)*DLX(2) + DNS2(L2,L3)*DLX(3)
      DER(2,5) = DNS1(L2,L3)*DLY(2) + DNS2(L2,L3)*DLY(3)
      DER(1,6) = DNS1(L3,L2)*DLX(3) + DNS2(L3,L2)*DLX(2)
      DER(2,6) = DNS1(L3,L2)*DLY(3) + DNS2(L3,L2)*DLY(2)
      DER(1,7) = DNVL(L3)*DLX(3)
      DER(2,7) = DNVL(L3)*DLY(3)
      DER(1,8) = DNS1(L3,L1)*DLX(3) + DNS2(L3,L1)*DLX(1)
      DER(2,8) = DNS1(L3,L1)*DLY(3) + DNS2(L3,L1)*DLY(1)
      DER(1,9) = DNS1(L1,L3)*DLX(1) + DNS2(L1,L3)*DLX(3)
      DER(2,9) = DNS1(L1,L3)*DLY(1) + DNS2(L1,L3)*DLY(3)
      DER(1,10) = 27.D0*(DLX(1)*L2*L3+L1*DLX(2)*L3+L1*L2*DLX(3))
      DER(2,10) = 27.D0*(DLY(1)*L2*L3+L1*DLY(2)*L3+L1*L2*DLY(3))
      RETURN
      END
