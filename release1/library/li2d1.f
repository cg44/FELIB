      SUBROUTINE LI2D1(FUN, IFUN, DER, IDER, XI, ITEST)
      INTEGER IDER, IFUN, ITEST, IERROR, ERRMES
      DOUBLE PRECISION DER, FUN, XI, XI1, XI2, SRNAME
      DIMENSION DER(IDER), FUN(IFUN)
      DATA SRNAME /8H LI2D1  /
      IF(ITEST.EQ.-1) GO TO 999
      IERROR=0
      IF(IFUN.LT.4.OR.IDER.LT.4) IERROR=1
      IF(XI.LT.-1.0D0.OR.XI.GT.1.0D0) IERROR=2
      ITEST=ERRMES(ITEST,IERROR,SRNAME)
      IF(ITEST.NE.0) RETURN
999   XI1 = (XI-1.0D0)*(XI-1.0D0)
      XI2 = (XI+1.0D0)*(XI+1.0D0)
      FUN(1) = 0.25D0*(2.0D0+XI)*XI1
      FUN(2) = 0.25D0*(1.0D0+XI)*XI1
      FUN(3) = 0.25D0*(2.0D0-XI)*XI2
      FUN(4) = 0.25D0*(XI-1.0D0)*XI2
      RETURN
      END
