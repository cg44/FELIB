C
      SUBROUTINE QQUA9(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      QQUA9 returns weights and abscissae of a 9-point gaussian
C      product quadrature rule for evaluating the integral of a
C      2D function over a rectangular region
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
C      Commented    16 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      IWGHT   dimension of vector WGHT (.GE. NQP)
C      IABSS   first dimension of array ABSS (.GE. 2)
C      JABSS   second dimension of array ABSS (.GE. NQP)
C      ITEST   error checking option
C
C ARGUMENTS out
C      WGHT    vector of length IWGHT.  contains weights to be
C              used in quadrature formula
C      ABSS    contains abscissae of points to be used in
C              formula
C      NQP     number of quadrature points (=9)
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE QQUA9(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,IABSS,IERROR,ITEST,IWGHT,JABSS,NQP
      CHARACTER*6 SRNAME
      DOUBLE PRECISION ABSS,AREA,WGHT,XY
      DIMENSION ABSS(IABSS,JABSS),WGHT(IWGHT)
      DATA SRNAME/'QQUA9'/
C
      NQP = 9
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IWGHT.LT.NQP) IERROR = 1
         IF (IABSS.LT.2 .OR. JABSS.LT.NQP) IERROR = 2
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      AREA = 4.0D0
      WGHT(1) = AREA*25.0D0/324.0D0
      WGHT(3) = WGHT(1)
      WGHT(5) = WGHT(1)
      WGHT(7) = WGHT(1)
      WGHT(2) = AREA*10.0D0/81.0D0
      WGHT(4) = WGHT(2)
      WGHT(6) = WGHT(2)
      WGHT(8) = WGHT(2)
      WGHT(9) = AREA*16.0D0/81.0D0
C
      XY = DSQRT(3.0D0/5.0D0)
      ABSS(1,1) = -XY
      ABSS(2,1) = -XY
      ABSS(1,2) = -XY
      ABSS(2,2) = 0.0D0
      ABSS(1,3) = -XY
      ABSS(2,3) = XY
      ABSS(1,4) = 0.0D0
      ABSS(2,4) = XY
      ABSS(1,5) = XY
      ABSS(2,5) = XY
      ABSS(1,6) = XY
      ABSS(2,6) = 0.0D0
      ABSS(1,7) = XY
      ABSS(2,7) = -XY
      ABSS(1,8) = 0.0D0
      ABSS(2,8) = -XY
      ABSS(1,9) = 0.0D0
      ABSS(2,9) = 0.0D0
C
      END
