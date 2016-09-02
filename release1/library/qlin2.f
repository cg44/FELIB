      SUBROUTINE QLIN2(WGHT, IWGHT, ABSS, IABSS, NQP, ITEST)
      INTEGER IABSS, ITEST, IWGHT, NQP,IERROR, ERRMES
      DOUBLE PRECISION ABSS, WGHT, SRNAME
      DIMENSION ABSS(IABSS), WGHT(IWGHT)
      DATA SRNAME /8H QLIN2  /
      NQP = 2
      IF(ITEST.EQ.-1) GO TO 999
      IERROR=0
      IF(IWGHT.LT.NQP) IERROR=1
      IF(IABSS.LT.NQP) IERROR=2
      ITEST=ERRMES(ITEST,IERROR,SRNAME)
      IF(ITEST.NE.0) RETURN
999   ABSS(1) = 1.0D0/DSQRT(3.0D0)
      ABSS(2) = -ABSS(1)
      WGHT(1) = 1.0D0
      WGHT(2) = 1.0D0
      RETURN
      END