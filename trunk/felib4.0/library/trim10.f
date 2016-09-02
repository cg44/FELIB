C
      SUBROUTINE TRIM10(FUN,IFUN,DER,IDER,JDER,XI,ETA,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      TRIM10 returns the values of the shape functions and their
C      derivatives at a specified point for a 10-noded c0
C      continuous triangular element. The function is 
C      continuous across element boundaries.
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
C      Commented    22 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      IFUN    length of vector FUN (.GE. 10)
C      IDER    first dimension of array DER (.GE. 2)
C      JDER    second dimension of array DER (.GE. 10)
C      XI      first local coordinate
C      ETA     second local coordinate
C      ITEST   error checking option
C
C ARGUMENTS out
C      FUN     vector of length IFUN. FUN(I) contains the
C              value of the i'th shape function at the point
C              (XI,ETA), for i=1(1)10
C      DER     array of dimension (IDER, JDER). DER(I, J)
C              contains the value of the derivative of the j'th
C              shape function with respect to the i'th local
C              coordinate, for i=1(1)2 and j=1(1)10
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE TRIM10(FUN,IFUN,DER,IDER,JDER,XI,ETA,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,IDER,IERROR,IFUN,ITEST,JDER
      CHARACTER*6 SRNAME
      DOUBLE PRECISION DER,DLX,DLY,DNS1,DNS2,DNVL,DUMMY,ETA,FUN,L1,L2,
     *                 L3,LA,LB,NS,NV,VEPS,XI,XMAX,XMIN,YMAX,YMIN
      DIMENSION DER(IDER,JDER),DLX(3),DLY(3),FUN(IFUN)
C
      DATA SRNAME/'TRIM10'/
C
C     Statement functions
C
      NV(LA) = 1.D0/2.D0* (3.D0*LA-1.D0)* (3.D0*LA-2.D0)*LA
      NS(LA,LB) = 9.D0/2.D0*LA*LB* (3.D0*LA-1.D0)
      DNVL(LA) = 1.D0/2.D0* (27.D0*LA*LA-18.D0*LA+2.D0)
      DNS1(LA,LB) = 9.D0/2.D0*LB* (6.D0*LA-1.D0)
      DNS2(LA,LB) = 9.D0/2.D0*LA* (3.D0*LA-1.D0)
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IFUN.LT.10) IERROR = 1
         IF (IDER.LT.2 .OR. JDER.LT.10) IERROR = 2
         YMIN = 1.0D0/DSQRT(3.0D0)* (XI-1.0D0) - VEPS(DUMMY)
         YMAX = 1.0D0/DSQRT(3.0D0)* (1.0D0-XI) + VEPS(DUMMY)
         XMIN = - (0.5D0+VEPS(DUMMY))
         XMAX = 1.0D0 + VEPS(DUMMY)
         IF ((XI.LT.XMIN.OR.XI.GT.XMAX) .OR.
     *       (ETA.LT.YMIN.OR.ETA.GT.YMAX)) IERROR = 3
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      L1 = 1.0D0/3.0D0* (1.0D0+2.0D0*XI)
      L2 = 1.0D0/3.0D0* (1.0D0-XI-DSQRT(3.0D0)*ETA)
      L3 = 1.0D0/3.0D0* (1.0D0-XI+DSQRT(3.0D0)*ETA)
      DLX(1) = 2.D0/3.D0
      DLY(1) = 0.D0
      DLX(2) = -1.D0/3.D0
      DLY(2) = -1.D0/DSQRT(3.D0)
      DLX(3) = -1.D0/3.D0
      DLY(3) = 1.D0/DSQRT(3.D0)
C
C     Shape functions
C
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
C
C     Derivatives
C
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
      DER(1,10) = 27.D0* (DLX(1)*L2*L3+L1*DLX(2)*L3+L1*L2*DLX(3))
      DER(2,10) = 27.D0* (DLY(1)*L2*L3+L1*DLY(2)*L3+L1*L2*DLY(3))
C
      END
