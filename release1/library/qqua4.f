      SUBROUTINE QQUA4(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)
      INTEGER IABSS, ITEST, IWGHT, JABSS, NQP, IERROR, ERRMES
      DOUBLE PRECISION ABSS, AREA, WGHT, SRNAME
      DIMENSION ABSS(IABSS,JABSS), WGHT(IWGHT)
      DATA SRNAME /8H QQUA4  /
      NQP = 4
      IF(ITEST.EQ.-1) GO TO 999
      IERROR=0
      IF(IWGHT.LT.NQP) IERROR=1
      IF(IABSS.LT.2.OR.JABSS.LT.NQP) IERROR=2
      ITEST=ERRMES(ITEST,IERROR,SRNAME)
      IF(ITEST.NE.0) RETURN
999   AREA = 4.0D0
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
