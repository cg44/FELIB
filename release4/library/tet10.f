C
      SUBROUTINE TET10(FUN,IFUN,DER,IDER,JDER,XI,ETA,ZETA,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      TET10 returns the values of shape functions and derivatives
C      at a specified point for an 10-noded tetrahedral element.
C      The function is continuous across element boundaries.
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
C      IFUN    dimension of vector FUN (.GE. 10)
C      IDER    first dimension of array DER (.GE. 3)
C      JDER    second dimension of array DER (.GE. 10)
C      XI      value of local coordinate at which function and
C              derivative values required
C      ETA     value of second local coordinate
C      ZETA    value of third local coordinate
C      ITEST   error checking option
C
C ARGUMENTS out
C      FUN     real vector of dimension IFUN. FUN(I) contains
C              value of i'th shape function at (XI, ETA, ZETA)
C      DER     real array of dimensions (IDER, JDER). DER(I,J)
C              contains the derivative of the j'th shape
C              function with respect to the i'th coordinate at
C              the point (XI, ETA, ZETA)
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE TET10(FUN,IFUN,DER,IDER,JDER,XI,ETA,ZETA,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,IDER,IERROR,IFUN,ITEST,JDER
      CHARACTER*6 SRNAME
      DOUBLE PRECISION DER,DL1X,DL1Y,DL1Z,DL2X,DL2Y,DL2Z,DL3X,DL3Y,DL3Z,
     *                 DL4X,DL4Y,DL4Z,DNSL,DNVL,ETA,FUN,L1,L2,L3,L4,LA,
     *                 LB,NS,NV,XI,ZETA
      DIMENSION DER(IDER,JDER),FUN(IFUN)
      DATA SRNAME/'TET10'/
C
C     Statement functions
C
      NV(LA) = (2.D0*LA-1.D0)*LA
      NS(LA,LB) = 4.D0*LA*LB
      DNVL(LA) = 4.D0*LA - 1.D0
      DNSL(LA,LB) = 4.D0*LB
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IFUN.LT.10) IERROR = 1
         IF (IDER.LT.3 .OR. JDER.LT.10) IERROR = 2
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      L1 = 1.D0/12.D0* (3.D0+8.D0*XI-2.D0*DSQRT(2.D0)*ZETA)
      L2 = 1.D0/12.D0* (3.D0-4.D0* (XI+DSQRT(3.D0)*ETA)-
     *     2.D0*DSQRT(2.D0)*ZETA)
      L3 = 1.D0/12.D0* (3.D0-4.D0* (XI-DSQRT(3.D0)*ETA)-
     *     2.D0*DSQRT(2.D0)*ZETA)
      L4 = 1.D0/4.D0* (1.D0+2.D0*DSQRT(2.D0)*ZETA)
C
      DL1X = 2.D0/3.D0
      DL1Y = 0.D0
      DL1Z = -1.D0/ (3.D0*DSQRT(2.D0))
      DL2X = -1.D0/3.D0
      DL2Y = -1.D0/DSQRT(3.D0)
      DL2Z = -1.D0/ (3.D0*DSQRT(2.D0))
      DL3X = -1.D0/3.D0
      DL3Y = 1.D0/DSQRT(3.D0)
      DL3Z = -1.D0/ (3.D0*DSQRT(2.D0))
      DL4X = 0.D0
      DL4Y = 0.D0
      DL4Z = 1.D0/DSQRT(2.D0)
C
C     Shape functions
C
      FUN(1) = NV(L1)
      FUN(2) = NS(L1,L2)
      FUN(3) = NV(L2)
      FUN(4) = NS(L2,L3)
      FUN(5) = NV(L3)
      FUN(6) = NS(L3,L1)
      FUN(7) = NS(L1,L4)
      FUN(8) = NS(L2,L4)
      FUN(9) = NS(L3,L4)
      FUN(10) = NV(L4)
C
C     Derivatives
C
      DER(1,1) = DNVL(L1)*DL1X
      DER(2,1) = DNVL(L1)*DL1Y
      DER(3,1) = DNVL(L1)*DL1Z
      DER(1,2) = DNSL(L1,L2)*DL1X + DNSL(L2,L1)*DL2X
      DER(2,2) = DNSL(L1,L2)*DL1Y + DNSL(L2,L1)*DL2Y
      DER(3,2) = DNSL(L1,L2)*DL1Z + DNSL(L2,L1)*DL2Z
      DER(1,3) = DNVL(L2)*DL2X
      DER(2,3) = DNVL(L2)*DL2Y
      DER(3,3) = DNVL(L2)*DL2Z
      DER(1,4) = DNSL(L2,L3)*DL2X + DNSL(L3,L2)*DL3X
      DER(2,4) = DNSL(L2,L3)*DL2Y + DNSL(L3,L2)*DL3Y
      DER(3,4) = DNSL(L2,L3)*DL2Z + DNSL(L3,L2)*DL3Z
      DER(1,5) = DNVL(L3)*DL3X
      DER(2,5) = DNVL(L3)*DL3Y
      DER(3,5) = DNVL(L3)*DL3Z
      DER(1,6) = DNSL(L3,L1)*DL3X + DNSL(L1,L3)*DL1X
      DER(2,6) = DNSL(L3,L1)*DL3Y + DNSL(L1,L3)*DL1Y
      DER(3,6) = DNSL(L3,L1)*DL3Z + DNSL(L1,L3)*DL1Z
      DER(1,7) = DNSL(L1,L4)*DL1X + DNSL(L4,L1)*DL4X
      DER(2,7) = DNSL(L1,L4)*DL1Y + DNSL(L4,L1)*DL4Y
      DER(3,7) = DNSL(L1,L4)*DL1Z + DNSL(L4,L1)*DL4Z
      DER(1,8) = DNSL(L2,L4)*DL2X + DNSL(L4,L2)*DL4X
      DER(2,8) = DNSL(L2,L4)*DL2Y + DNSL(L4,L2)*DL4Y
      DER(3,8) = DNSL(L2,L4)*DL2Z + DNSL(L4,L2)*DL4Z
      DER(1,9) = DNSL(L3,L4)*DL3X + DNSL(L4,L3)*DL4X
      DER(2,9) = DNSL(L3,L4)*DL3Y + DNSL(L4,L3)*DL4Y
      DER(3,9) = DNSL(L3,L4)*DL3Z + DNSL(L4,L3)*DL4Z
      DER(1,10) = DNVL(L4)*DL4X
      DER(2,10) = DNVL(L4)*DL4Y
      DER(3,10) = DNVL(L4)*DL4Z
C
      END
