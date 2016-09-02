      SUBROUTINE QU8FN(FUN, IFUN, DER, IDER, JDER, XI, ETA, ITEST)
      INTEGER IDER, IFUN, ITEST, JDER, IERROR, ERRMES
      DOUBLE PRECISION DER, ETA, ETAM, ETAP, FUN, XI, XIM,SRNAME,
     *     XIP
      DIMENSION DER(IDER,JDER), FUN(IFUN)
      DATA SRNAME /8H QU8FN  /
      IF(ITEST.EQ.-1) GO TO 999
      IERROR=0
      IF(IFUN.LT.8) IERROR=1
      IF(IDER.LT.2.OR.JDER.LT.8) IERROR=2
      IF(DABS(XI).GT.1.0D0.OR.DABS(ETA).GT.1.0D0)
     *    IERROR=3
      ITEST=ERRMES(ITEST,IERROR,SRNAME)
      IF(ITEST.NE.0) RETURN
999   ETAM = 0.250D0*(1.0D0-ETA)
      ETAP = 0.250D0*(1.0D0+ETA)
      XIM = 0.250D0*(1.0D0-XI)
      XIP = 0.250D0*(1.0D0+XI)
      DER(1,1) = ETAM*(2.0D0*XI+ETA)
      DER(1,2) = -8.0D0*ETAM*ETAP
      DER(1,3) = ETAP*(2.0D0*XI-ETA)
      DER(1,4) = -4.0D0*ETAP*XI
      DER(1,5) = ETAP*(2.0D0*XI+ETA)
      DER(1,6) = 8.0D0*ETAP*ETAM
      DER(1,7) = ETAM*(2.0D0*XI-ETA)
      DER(1,8) = -4.0D0*ETAM*XI
      DER(2,1) = XIM*(XI+2.0D0*ETA)
      DER(2,2) = -4.0D0*XIM*ETA
      DER(2,3) = XIM*(2.0D0*ETA-XI)
      DER(2,4) = 8.0D0*XIM*XIP
      DER(2,5) = XIP*(XI+2.0D0*ETA)
      DER(2,6) = -4.0D0*XIP*ETA
      DER(2,7) = XIP*(2.0D0*ETA-XI)
      DER(2,8) = -8.0D0*XIM*XIP
      FUN(1) = 4.0D0*ETAM*XIM*(-XI-ETA-1.0D0)
      FUN(2) = 32.0D0*XIM*ETAM*ETAP
      FUN(3) = 4.0D0*ETAP*XIM*(-XI+ETA-1.0D0)
      FUN(4) = 32.0D0*XIM*XIP*ETAP
      FUN(5) = 4.0D0*XIP*ETAP*(XI+ETA-1.0D0)
      FUN(6) = 32.0D0*XIP*ETAP*ETAM
      FUN(7) = 4.0D0*XIP*ETAM*(XI-ETA-1.0D0)
      FUN(8) = 32.0D0*XIM*XIP*ETAM
      RETURN
      END