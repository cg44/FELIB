C***********************************************************************
C$SPLIT$QWDG6$*********************************************************
C***********************************************************************
      SUBROUTINE QWDG6(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS WEIGHTS AND ABSCISSAE FOR A 6-POINT QUADRATURE
C      RULE FOR EVALUATING THE INTEGRAL OF A 3D FUNCTION OVER
C      A PENTAHEDRAL REGION
C
C HISTORY
C      RELEASE 2.0   6 OCT 1980 (CG) --- SERC COPYRIGHT
C      COMMENTED    21 OCT 1980 (KR)
C
C ARGUMENTS IN
C      IWGHT   DIMENSION OF VECTOR WGHT (.GE.NQP(=6))
C      IABSS   FIRST DIMENSION ARRAY ABSS (.GE.3)
C      JABSS   SECOND DIMENSION OF ARRAY ABSS (.GE.NQP(=6))
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
C     SUBROUTINE QWDG6(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP,
C    *     ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IABSS, IERROR, ITEST, IWGHT, JABSS,
     *     NQP
      DOUBLE PRECISION ABSS, SRNAME, WGHT
      DIMENSION ABSS(IABSS,JABSS), WGHT(IWGHT)
      DATA SRNAME /8H QWDG6  /
      NQP = 6
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IWGHT.LT.NQP) IERROR = 1
                        IF (IABSS.LT.3 .OR. JABSS.LT.NQP)
     *                      IERROR = 2
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 DO 1020 I=1,6
      WGHT(I) = 1.0D0/4.0D0*DSQRT(3.0D0)
 1020 CONTINUE
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
      RETURN
      END
