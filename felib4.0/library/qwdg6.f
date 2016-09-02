C
      SUBROUTINE QWDG6(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      QWDG6 returns weights and abscissae for a 6-point quadrature
C      rule for evaluating the integral of a 3D function over
C      a pentahedral region
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
C      Release 2.0   6 Oct 1980 (CG)
C      Commented    21 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      IWGHT   dimension of vector WGHT (.GE. NQP(=6))
C      IABSS   first dimension array ABSS (.GE. 3)
C      JABSS   second dimension of array ABSS (.GE. NQP(=6))
C      ITEST   error checking option
C
C ARGUMENTS out
C      WGHT    vector of length IWGHT.  contains weights of
C              quadrature rule
C      ABSS    array of dimension (IABSS,JABSS).  contains
C              abscissae of quadrature points
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE QWDG6(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IABSS,IERROR,ITEST,IWGHT,JABSS,NQP
      CHARACTER*6 SRNAME
      DOUBLE PRECISION ABSS,WGHT
      DIMENSION ABSS(IABSS,JABSS),WGHT(IWGHT)
      DATA SRNAME/'QWDG6'/
C
      NQP = 6
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
      DO 1000 I = 1,6
         WGHT(I) = 1.0D0/4.0D0*DSQRT(3.0D0)
 1000 CONTINUE
C
      ABSS(1,1) = 0.0D0
      ABSS(2,1) = 0.0D0
      ABSS(3,1) = -1.0D0
      ABSS(1,2) = -0.5D0
      ABSS(2,2) = -DSQRT(3.0D0)/2.0D0
      ABSS(3,2) = -1.0D0
      ABSS(1,3) = -0.5D0
      ABSS(2,3) = DSQRT(3.0D0)/2.0D0
      ABSS(3,3) = -1.0D0
      ABSS(1,4) = 0.0D0
      ABSS(2,4) = 0.0D0
      ABSS(3,4) = 1.0D0
      ABSS(1,5) = -0.5D0
      ABSS(2,5) = -DSQRT(3.0D0)/2.0D0
      ABSS(3,5) = 1.0D0
      ABSS(1,6) = -0.5D0
      ABSS(2,6) = DSQRT(3.0D0)/2.0D0
      ABSS(3,6) = 1.0D0
C
      END
