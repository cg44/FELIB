C***********************************************************************
C$SPLIT$QWDG8$*********************************************************
C***********************************************************************
      SUBROUTINE QWDG8(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS WEIGHTS AND ABSCISSAE FOR A 8-POINT QUADRATURE
C      RULE FOR EVALUATING THE INTEGRAL OF A 3D FUNCTION OVER
C      A PENTAHEDRAL REGION
C
C HISTORY
C      RELEASE 2.0   6 OCT 1980 (CG) --- SERC COPYRIGHT
C      COMMENTED    21 OCT 1980 (KR)
C
C ARGUMENTS IN
C      IWGHT   DIMENSION OF VECTOR WGHT (.GE.NQP(=8))
C      IABSS   FIRST DIMENSION ARRAY ABSS (.GE.3)
C      JABSS   SECOND DIMENSION OF ARRAY ABSS (.GE.NQP(=8))
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      WGHT    VECTOR OF LENGTH IWGHT.  CONTAINS WEIGHTS OF
C              QUADRATURE RULE
C      ABSS    ARRAY OF DIMENSION (IABSS,JABSS).  CONTAINS
C              ABSCISSAE OF QUADRATURE POINTS
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE QWDG8(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IABSS, IERROR, ITEST, IWGHT, JABSS,
     *     NQP
      DOUBLE PRECISION ABSS, SRNAME, W, WGHT
      DIMENSION ABSS(IABSS,JABSS), WGHT(IWGHT)
      DATA SRNAME /8H QWDG8  /
      NQP = 8
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IWGHT.LT.NQP) IERROR = 1
                        IF (IABSS.LT.3 .OR. JABSS.LT.NQP)
     *                      IERROR = 2
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 W = (3.D0*DSQRT(3.D0))/2.D0
      WGHT(1) = (3.0D0/8.0D0)*W
      WGHT(2) = (3.0D0/8.0D0)*W
      DO 1020 I=3,8
      WGHT(I) = (1.0D0/24.0D0)*W
 1020 CONTINUE
      ABSS(1,1) = 0.0D0
      ABSS(2,1) = 0.0D0
      ABSS(3,1) = 1.0D0
      ABSS(1,2) = 0.0D0
      ABSS(2,2) = 0.0D0
      ABSS(3,2) = -1.0D0
      ABSS(1,3) = 1.0D0
      ABSS(2,3) = 0.0D0
      ABSS(3,3) = -1.0D0
      ABSS(1,4) = -0.5D0
      ABSS(2,4) = DSQRT(3.0D0)/2.0D0
      ABSS(3,4) = -1.0D0
      ABSS(1,5) = -0.5D0
      ABSS(2,5) = -ABSS(2,4)
      ABSS(3,5) = -1.0D0
      ABSS(1,6) = 1.0D0
      ABSS(2,6) = 0.0D0
      ABSS(3,6) = 1.0D0
      ABSS(1,7) = -0.5D0
      ABSS(2,7) = ABSS(2,4)
      ABSS(3,7) = 1.0D0
      ABSS(1,8) = -0.5D0
      ABSS(2,8) = ABSS(2,5)
      ABSS(3,8) = 1.0D0
      RETURN
      END
