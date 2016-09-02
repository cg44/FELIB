C***********************************************************************
C$SPLIT$QBRK8$*********************************************************
C***********************************************************************
      SUBROUTINE QBRK8(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS WEIGHTS AND ABSCISSAE OF A 8-POINT GAUSS TYPE
C      QUADRATURE RULE FOR USE IN EVALUATING THE INTEGRAL OF A
C      3D FUNCTION OVER A CUBE
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    15 OCT 1980 (KR)
C
C ARGUMENTS IN
C      IWGHT   DIMENSION OF VECTOR WGHT (.GE.NQP(=8))
C      IABSS   FIRST DIMENSION OF ARRAY ABSS (.GE.3)
C      JABSS   SECOND DIMENSION OF ARRAY ABSS (.GE.NQP(=8))
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      WGHT    VECTOR OF DIMENSION IWGHT.  CONTAINS WEIGHTS TO
C              BE USED IN THE 8-POINT QUADRATURE FORMULA
C      ABSS    ARRAY OF DIMENSION (IABSS,JABSS).  CONTAINS
C              ABSCISSAE OF POINTS TO BE USED IN QUADRATURE
C              FORMULA
C      NQP     NUMBER OF QUADRATURE POINTS TO BE USED (=8)
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE QBRK8(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IABSS, IERROR, ITEST, IWGHT, JABSS,
     *     NQP
      DOUBLE PRECISION ABSS, SRNAME, W, WGHT, XY
      DIMENSION ABSS(IABSS,JABSS), WGHT(IWGHT)
      DATA SRNAME /8H QBRK8  /
      NQP = 8
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IWGHT.LT.NQP) IERROR = 1
                        IF (IABSS.LT.3 .OR. JABSS.LT.NQP)
     *                      IERROR = 2
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 W = 1.0D0
      XY = DSQRT(1.0D0/3.0D0)
      DO 1020 I=1,8
      WGHT(I) = W
      ABSS(1,I) = XY
      ABSS(2,I) = XY
      ABSS(3,I) = XY
 1020 CONTINUE
      DO 1030 I=1,4
      ABSS(3,I) = -ABSS(1,1)
 1030 CONTINUE
      DO 1040 I=1,2
      ABSS(1,I+1) = -ABSS(1,1)
      ABSS(1,I+5) = -ABSS(1,1)
      ABSS(2,I) = -ABSS(1,1)
      ABSS(2,I+4) = -ABSS(1,1)
 1040 CONTINUE
      RETURN
      END
