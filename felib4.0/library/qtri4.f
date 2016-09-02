C
      SUBROUTINE QTRI4(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      QTRI4 returns weights and abscissae of a 4-point gauss-type
C      quadrature rule for evaluating the integral of a 2D
C      function over a triangular region
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
C      Commented    21 Oct 1980 (KR)
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
C      WGHT    vector of length IWGHT, containing NQP weights
C              of quadrature formula
C      ABSS    array of dimension (IABSS, JABSS).  contains
C              abscissae of points to be used in quadrature
C              rule
C      NQP     number of quadrature points (=4)
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE QTRI4(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,IABSS,IERROR,ITEST,IWGHT,JABSS,NQP
      CHARACTER*6 SRNAME
      DOUBLE PRECISION ABSS,AREA,WGHT
      DIMENSION ABSS(IABSS,JABSS),WGHT(IWGHT)
      DATA SRNAME/'QTRI4'/
C
      NQP = 4
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
      AREA = 3.0D0*DSQRT(3.0D0)/4.0D0
C
      WGHT(1) = AREA*1.0D0/12.0D0
      WGHT(2) = AREA*1.0D0/12.0D0
      WGHT(3) = AREA*1.0D0/12.0D0
      WGHT(4) = AREA*3.0D0/4.0D0
C
      ABSS(1,1) = 1.0D0
      ABSS(2,1) = 0.0D0
      ABSS(1,2) = -0.5D0
      ABSS(2,2) = -DSQRT(3.0D0)/2.0D0
      ABSS(1,3) = -0.5D0
      ABSS(2,3) = +DSQRT(3.0D0)/2.0D0
      ABSS(1,4) = 0.0D0
      ABSS(2,4) = 0.0D0
C
      END
