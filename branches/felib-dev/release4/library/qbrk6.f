C
      SUBROUTINE QBRK6(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      QBRK6 returns weights and abscissae of a six-point gauss type
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
C      IWGHT   dimension of vector WGHT (.GE. NQP(=6))
C      IABSS   first dimension of array ABSS (.GE. 3)
C      JABSS   second dimension of array ABSS (.GE. NQP(=6))
C      ITEST   error checking option
C
C ARGUMENTS out
C      WGHT    vector of dimension IWGHT.  contains weights to
C              be used in the 6-point quadrature formula
C      ABSS    array of dimension (IABSS, JABSS).  contains
C              abscissae of points to be used in quadrature
C              formula
C      NQP     number of quadrature points to be used (=6)
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE QBRK6(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IABSS,IERROR,ITEST,IWGHT,JABSS,NQP
      CHARACTER*6 SRNAME
      DOUBLE PRECISION ABSS,W,WGHT
      DIMENSION ABSS(IABSS,JABSS),WGHT(IWGHT)
      DATA SRNAME/'QBRK6'/
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
      W = 8.0D0/6.0D0
C
      DO 1000 I = 1,6
         WGHT(I) = W
         ABSS(1,I) = 0.0D0
         ABSS(2,I) = 0.0D0
         ABSS(3,I) = 0.0D0
 1000 CONTINUE
C
      ABSS(1,5) = -1.0D0
      ABSS(1,6) = 1.0D0
      ABSS(2,3) = -1.0D0
      ABSS(2,4) = 1.0D0
      ABSS(3,1) = -1.0D0
      ABSS(3,2) = 1.0D0
C
      END
