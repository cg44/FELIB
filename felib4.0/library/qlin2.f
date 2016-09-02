C
      SUBROUTINE QLIN2(WGHT,IWGHT,ABSS,IABSS,NQP,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      QLIN2 returns the weights and abscissae of a 2-point gauss-
C      legendre quadrature formula for use in evaluating a 1D
C      integral over a finite range
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
C      IWGHT   dimension of vector WGHT (.GE. NQP(=2))
C      IABSS   dimension of vector ABSS (.GE. NQP(=2))
C      ITEST   error checking option
C
C ARGUMENTS out
C      WGHT    vector of dimension IWGHT.  contains weights to
C              be used in the 2-point quadrature formula
C      ABSS    vector of dimension IABSS.  contains abscissae
C              of points to be used in 2-point quadrature
C              formula
C      NQP     number of quadrature points (=2)
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE QLIN2(WGHT,IWGHT,ABSS,IABSS,NQP,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,IABSS,IERROR,ITEST,IWGHT,NQP
      CHARACTER*6 SRNAME
      DOUBLE PRECISION ABSS,WGHT
      DIMENSION ABSS(IABSS),WGHT(IWGHT)
      DATA SRNAME/'QLIN2'/
C
      NQP = 2
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IWGHT.LT.NQP) IERROR = 1
         IF (IABSS.LT.NQP) IERROR = 2
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      ABSS(1) = 1.0D0/DSQRT(3.0D0)
      ABSS(2) = -ABSS(1)
      WGHT(1) = 1.0D0
      WGHT(2) = 1.0D0
C
      END
