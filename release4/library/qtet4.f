C
      SUBROUTINE QTET4(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      QTET4 returns weights and abscissae of a 4-point rule for
C      evaluating the integral of a 3D function over a
C      tetrahedral region
C
C HISTORY
C
C      Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
C                           Chilton, Didcot, Oxfordshire OX11 0QX
C
C      Contact: Prof Chris Greenough
C      Email: c.greenough@rl.ac.uk
C      Tel: (44) 1235 445307C
C
C      Release 2.0   6 Oct 1980 (CG)
C      Commented    16 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      IWGHT   dimension of vector WGHT (.GE. NQP)
C      IABSS   first dimension of array ABSS (.GE. 3)
C      JABSS   second dimension of array ABSS (.GE. NQP)
C      ITEST   error checking option
C
C ARGUMENTS out
C      WGHT    vector of dimension IWGHT.  contains weights for
C              4-point quadrature rule
C      ABSS    array of dimension (IABSS,JABSS).  contains
C              abscissae of points to be used in quadrature
C              rule
C      NQP     number of quadrature points (=4)
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE QTET4(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IABSS,IERROR,ITEST,IWGHT,JABSS,NQP
      CHARACTER*6 SRNAME
      DOUBLE PRECISION ABSS,ALPHA,BETA,GAMA,WGHT
      DIMENSION ABSS(IABSS,JABSS),WGHT(IWGHT)
      DATA SRNAME/'QTET4'/
C
      NQP = 4
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
      ALPHA = 0.223807D0
      BETA = 0.387298D0
      GAMA = 0.158114D0
C
      DO 1000 I = 1,4
         WGHT(I) = 1.0D0/16.0D0*DSQRT(6.0D0)
 1000 CONTINUE
C
      ABSS(1,1) = 2.0D0*ALPHA
      ABSS(2,1) = 0.0D0
      ABSS(3,1) = -GAMA
      ABSS(1,2) = -ALPHA
      ABSS(2,2) = -BETA
      ABSS(3,2) = -GAMA
      ABSS(1,3) = -ALPHA
      ABSS(2,3) = BETA
      ABSS(3,3) = -GAMA
      ABSS(1,4) = 0.0D0
      ABSS(2,4) = 0.0D0
      ABSS(3,4) = 3.0D0*GAMA
C
      END
