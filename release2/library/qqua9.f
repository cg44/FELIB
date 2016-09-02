C***********************************************************************
C$SPLIT$QQUA9$*********************************************************
C***********************************************************************
      SUBROUTINE QQUA9(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS WEIGHTS AND ABSCISSAE OF A 9-POINT GAUSSIAN
C      PRODUCT QUADRATURE RULE FOR EVALUATING THE INTEGRAL OF A
C      2D FUNCTION OVER A RECTANGULAR REGION
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    16 OCT 1980 (KR)
C
C ARGUMENTS IN
C      IWGHT   DIMENSION OF VECTOR WGHT (.GE.NQP)
C      IABSS   FIRST DIMENSION OF ARRAY ABSS (.GE.2)
C      JABSS   SECOND DIMENSION OF ARRAY ABSS (.GE.NQP)
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      WGHT    VECTOR OF LENGTH IWGHT.  CONTAINS WEIGHTS TO BE
C              USED IN QUADRATURE FORMULA
C      ABSS    CONTAINS ABSCISSAE OF POINTS TO BE USED IN
C              FORMULA
C      NQP     NUMBER OF QUADRATURE POINTS (=9)
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE QQUA9(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, IABSS, IERROR, ITEST, IWGHT, JABSS, NQP
      DOUBLE PRECISION ABSS, AREA, SRNAME, WGHT, XY
      DIMENSION ABSS(IABSS,JABSS), WGHT(IWGHT)
      DATA SRNAME /8H QQUA9  /
      NQP = 9
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IWGHT.LT.NQP) IERROR = 1
                        IF (IABSS.LT.2 .OR. JABSS.LT.NQP)
     *                      IERROR = 2
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 AREA = 4.0D0
      WGHT(1) = AREA*25.0D0/324.0D0
      WGHT(3) = WGHT(1)
      WGHT(5) = WGHT(1)
      WGHT(7) = WGHT(1)
      WGHT(2) = AREA*10.0D0/81.0D0
      WGHT(4) = WGHT(2)
      WGHT(6) = WGHT(2)
      WGHT(8) = WGHT(2)
      WGHT(9) = AREA*16.0D0/81.0D0
      XY = DSQRT(3.0D0/5.0D0)
      ABSS(1,1) = -XY
      ABSS(2,1) = -XY
      ABSS(1,2) = -XY
      ABSS(2,2) = 0.0D0
      ABSS(1,3) = -XY
      ABSS(2,3) = XY
      ABSS(1,4) = 0.0D0
      ABSS(2,4) = XY
      ABSS(1,5) = XY
      ABSS(2,5) = XY
      ABSS(1,6) = XY
      ABSS(2,6) = 0.0D0
      ABSS(1,7) = XY
      ABSS(2,7) = -XY
      ABSS(1,8) = 0.0D0
      ABSS(2,8) = -XY
      ABSS(1,9) = 0.0D0
      ABSS(2,9) = 0.0D0
      RETURN
      END
