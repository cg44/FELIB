C
      SUBROUTINE QLIN3(WGHT,IWGHT,ABSS,IABSS,NQP,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      QLIN3 returns weights and abscissae of a 3-point gauss-
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
C      IWGHT   dimension of vector WGHT (.GE. NQP(=3))
C      IABSS   dimension of vector ABSS (.GE. NQP(=3))
C      ITEST   error checking option
C
C ARGUMENTS out
C      WGHT    vector of dimension IWGHT. Contains weights for
C              3-point formula
C      ABSS    vector of dimension IABSS.  contains abscissae
C              of points for use in 3-point formula
C      NQP     number of quadrature points (=3)
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE QLIN3(WGHT,IWGHT,ABSS,IABSS,NQP,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,IABSS,IERROR,ITEST,IWGHT,NQP
      CHARACTER*6 SRNAME
      DOUBLE PRECISION ABSS,WGHT
      DIMENSION ABSS(IABSS),WGHT(IWGHT)
      DATA SRNAME/'QLIN3'/
C
      NQP = 3
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
      ABSS(1) = 0.2D0*DSQRT(15.0D0)
      ABSS(2) = 0.0D0
      ABSS(3) = -ABSS(1)
      WGHT(1) = 5.0D0/9.0D0
      WGHT(3) = WGHT(1)
      WGHT(2) = 8.0D0/9.0D0
C
      END
