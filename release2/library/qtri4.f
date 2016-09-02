C***********************************************************************
C$SPLIT$QTRI4$*********************************************************
C***********************************************************************
      SUBROUTINE QTRI4(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS WEIGHTS AND ABSCISSAE OF A 4-POINT GAUSS-TYPE
C      QUADRATURE RULE FOR EVALUATING THE INTEGRAL OF A 2D
C      FUNCTION OVER A TRIANGULAR REGION
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    21 OCT 1980 (KR)
C
C ARGUMENTS IN
C      IWGHT   DIMENSION OF VECTOR WGHT (.GE.NQP)
C      IABSS   FIRST DIMENSION OF ARRAY ABSS (.GE.2)
C      JABSS   SECOND DIMENSION OF ARRAY ABSS (.GE.NQP)
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      WGHT    VECTOR OF LENGTH IWGHT, CONTAINING NQP WEIGHTS
C              OF QUADRATURE FORMULA
C      ABSS    ARRAY OF DIMENSION (IABSS,JABSS).  CONTAINS
C              ABSCISSAE OF POINTS TO BE USED IN QUADRATURE
C              RULE
C      NQP     NUMBER OF QUADRATURE POINTS (=4)
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE QTRI4(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, IABSS, IERROR, ITEST, IWGHT, JABSS, NQP
      DOUBLE PRECISION ABSS, AREA, SRNAME, WGHT
      DIMENSION ABSS(IABSS,JABSS), WGHT(IWGHT)
      DATA SRNAME /8H QTRI4  /
      NQP = 4
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IWGHT.LT.NQP) IERROR = 1
                        IF (IABSS.LT.2 .OR. JABSS.LT.NQP)
     *                      IERROR = 2
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 AREA = 3.0D0*DSQRT(3.0D0)/4.0D0
      WGHT(1) = AREA*1.0D0/12.0D0
      WGHT(2) = AREA*1.0D0/12.0D0
      WGHT(3) = AREA*1.0D0/12.0D0
      WGHT(4) = AREA*3.0D0/4.0D0
      ABSS(1,1) = 1.0D0
      ABSS(2,1) = 0.0D0
      ABSS(1,2) = -0.5D0
      ABSS(2,2) = -DSQRT(3.0D0)/2.0D0
      ABSS(1,3) = -0.5D0
      ABSS(2,3) = +DSQRT(3.0D0)/2.0D0
      ABSS(1,4) = 0.0D0
      ABSS(2,4) = 0.0D0
      RETURN
      END
