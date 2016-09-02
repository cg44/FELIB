C
      SUBROUTINE QBRK8(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      QBRK8 returns weights and abscissae of a 8-point gauss type
C      quadrature rule for use in evaluating the integral of a
C      3D function over a cube
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
C      Commented    15 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      IWGHT   dimension of vector WGHT (.GE. NQP(=8))
C      IABSS   first dimension of array ABSS (.GE. 3)
C      JABSS   second dimension of array ABSS (.GE. NQP(=8))
C      ITEST   error checking option
C
C ARGUMENTS out
C      WGHT    vector of dimension IWGHT.  contains weights to
C              be used in the 8-point quadrature formula
C      ABSS    array of dimension (IABSS, JABSS).  contains
C              abscissae of points to be used in quadrature
C              formula
C      NQP     number of quadrature points to be used (=8)
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE QBRK8(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IABSS,IERROR,ITEST,IWGHT,JABSS,NQP
      CHARACTER*6 SRNAME
      DOUBLE PRECISION ABSS,W,WGHT,XY
      DIMENSION ABSS(IABSS,JABSS),WGHT(IWGHT)
      DATA SRNAME/'QBRK8'/
C
      NQP = 8
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IWGHT.LT.NQP) IERROR = 1
         IF (IABSS.LT.3 .OR. JABSS.LT.NQP) IERROR = 2
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      W = 1.0D0
      XY = DSQRT(1.0D0/3.0D0)
      DO 1000 I = 1,8
         WGHT(I) = W
         ABSS(1,I) = XY
         ABSS(2,I) = XY
         ABSS(3,I) = XY
 1000 CONTINUE
C
      DO 1010 I = 1,4
         ABSS(3,I) = -ABSS(1,1)
 1010 CONTINUE
C
      DO 1020 I = 1,2
         ABSS(1,I+1) = -ABSS(1,1)
         ABSS(1,I+5) = -ABSS(1,1)
         ABSS(2,I) = -ABSS(1,1)
         ABSS(2,I+4) = -ABSS(1,1)
 1020 CONTINUE
C
      END
