      SUBROUTINE QBRK8(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)
      INTEGER I, IABSS, ITEST, IWGHT, JABSS, NQP, IERROR, ERRMES
      DOUBLE PRECISION ABSS, W, WGHT, XY, SRNAME
      DIMENSION ABSS(IABSS,JABSS), WGHT(IWGHT)
      DATA SRNAME /8H QBRK8  /
      NQP = 8
      IF(ITEST.EQ.-1) GO TO 999
      IERROR=0
      IF(IWGHT.LT.NQP) IERROR=1
      IF(IABSS.LT.3.OR.JABSS.LT.NQP) IERROR=2
      ITEST=ERRMES(ITEST,IERROR,SRNAME)
      IF(ITEST.NE.0) RETURN
999   W = 8.0D0*(1.0D0/8.0D0)
      XY = DSQRT(1.0D0/3.0D0)
      DO 1010 I=1,8
      WGHT(I) = W
      ABSS(1,I) = XY
      ABSS(2,I) = XY
      ABSS(3,I) = XY
 1010 CONTINUE
      DO 1020 I=1,4
      ABSS(1,I+4) = -ABSS(1,I)
      ABSS(2,2*I) = -ABSS(1,I)
 1020 CONTINUE
      DO 1030 I=1,2
      ABSS(3,I+2) = -ABSS(1,I)
      ABSS(3,I+6) = -ABSS(1,I)
 1030 CONTINUE
      RETURN
      END
