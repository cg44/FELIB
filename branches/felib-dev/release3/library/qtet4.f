C***********************************************************************
      SUBROUTINE QTET4(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS WEIGHTS AND ABSCISSAE OF A 4-POINT RULE FOR
C      EVALUATING THE INTEGRAL OF A 3D FUNCTION OVER A
C      TETRAHEDRAL REGION
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 2.0   6 OCT 1980 (CG)
C      COMMENTED    16 OCT 1980 (KR)
C
C ARGUMENTS IN
C      IWGHT   DIMENSION OF VECTOR WGHT (.GE.NQP)
C      IABSS   FIRST DIMENSION OF ARRAY ABSS (.GE.3)
C      JABSS   SECOND DIMENSION OF ARRAY ABSS (.GE.NQP)
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      WGHT    VECTOR OF DIMENSION IWGHT.  CONTAINS WEIGHTS FOR
C              4-POINT QUADRATURE RULE
C      ABSS    ARRAY OF DIMENSION (IABSS,JABSS).  CONTAINS
C              ABSCISSAE OF POINTS TO BE USED IN QUADRATURE
C              RULE
C      NQP     NUMBER OF QUADRATURE POINTS (=4)
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE QTET4(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IABSS, IERROR, ITEST, IWGHT, JABSS,
     *     NQP
      DOUBLE PRECISION ABSS, ALPHA, BETA, GAMA, SRNAME, WGHT
      DIMENSION ABSS(IABSS,JABSS), WGHT(IWGHT)
      DATA SRNAME /8H QTET4  /
C
      NQP = 4
C
C     PARAMETER CHECKING
C
      IF (ITEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (IWGHT.LT.NQP) IERROR = 1
      IF (IABSS.LT.3 .OR. JABSS.LT.NQP) IERROR = 2
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
C     MAIN BODY
C
 1010 ALPHA = 0.223807D0
      BETA = 0.387298D0
      GAMA = 0.158114D0
C
      DO 1020 I=1,4
         WGHT(I) = 1.0D0/16.0D0*DSQRT(6.0D0)
 1020 CONTINUE
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
      RETURN
      END
C***********************************************************************