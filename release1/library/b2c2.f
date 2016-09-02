      SUBROUTINE B2C2(BEE, IBEE, JBEE, DERIV, IDERIV, JDERIV,
     *     NODEL, ITEST)
      INTEGER IBEE, IDERIV, ITEST, JBEE, JDERIV, K, L,
     *     M, NODEL,IERROR,ERRMES
      DOUBLE PRECISION BEE, DERIV,SRNAME
      DIMENSION BEE(IBEE,JBEE), DERIV(IDERIV,JDERIV)
      DATA SRNAME /8H B2C2   /
      K = 3
      L = 2*NODEL
      IF(ITEST.EQ.-1) GO TO 999
      IERROR=0
      IF(IDERIV.LT.2.OR.JDERIV.LT.NODEL) IERROR=3
      IF(IBEE.LT.3.OR.JBEE.LT.L) IERROR=2
      IF(NODEL.LE.0) IERROR=1
      ITEST=ERRMES(ITEST,IERROR,SRNAME)
      IF(ITEST.NE.0) RETURN
999   CALL MATNUL(BEE, IBEE, JBEE, K, L, ITEST)
      DO 1010 M=1,NODEL
      K = 2*M
      L = K - 1
      BEE(1,L) = DERIV(1,M)
      BEE(3,K) = DERIV(1,M)
      BEE(2,K) = DERIV(2,M)
      BEE(3,L) = DERIV(2,M)
 1010 CONTINUE
      RETURN
      END