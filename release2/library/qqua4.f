C***********************************************************************
C$SPLIT$QQUA4$*********************************************************
C***********************************************************************
      SUBROUTINE QQUA4(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS WEIGHTS AND ABSCISSAE OF A 4-POINT GAUSSIAN
C      PRODUCT QUADRATURE RULE FOR USE IN EVALUATING THE
C      INTEGRAL OF A 2D FUNCTION OVER A RECTANGULAR REGION
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
C      WGHT    VECTOR OF DIMENSION IWGHT.  CONTAINS WEIGHTS TO
C              BE USED IN 4-POINT FORMULA
C      ABSS    ARRAY OF DIMENSION (IABSS,JABSS).  CONTAINS
C              ABSCISSAE OF POINTS TO BE USED IN FORMULA
C      NQP     NUMBER OF POINTS TO BE USED IN FORMULA (=4)
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE QQUA4(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, IABSS, IERROR, ITEST, IWGHT, JABSS, NQP
      DOUBLE PRECISION ABSS, AREA, SRNAME, WGHT
      DIMENSION ABSS(IABSS,JABSS), WGHT(IWGHT)
      DATA SRNAME /8H QQUA4  /
      NQP = 4
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IWGHT.LT.NQP) IERROR = 1
                        IF (IABSS.LT.2 .OR. JABSS.LT.NQP)
     *                      IERROR = 2
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 AREA = 4.0D0
      WGHT(1) = AREA*0.25D+0
      WGHT(2) = AREA*0.25D+0
      WGHT(3) = AREA*0.25D+0
      WGHT(4) = AREA*0.25D+0
      ABSS(1,1) = -1.0D0/DSQRT(3.0D0)
      ABSS(1,2) = ABSS(1,1)
      ABSS(1,3) = +1.0D0/DSQRT(3.0D0)
      ABSS(1,4) = ABSS(1,3)
      ABSS(2,1) = -1.0D0/DSQRT(3.0D0)
      ABSS(2,2) = +1.0D0/DSQRT(3.0D0)
      ABSS(2,3) = ABSS(2,2)
      ABSS(2,4) = ABSS(2,1)
      RETURN
      END
