C***********************************************************************
C$SPLIT$QLIN2$*********************************************************
C***********************************************************************
      SUBROUTINE QLIN2(WGHT, IWGHT, ABSS, IABSS, NQP, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS THE WEIGHTS AND ABSCISSAE OF A 2-POINT GAUSS-
C      LEGENDRE QUADRATURE FORMULA FOR USE IN EVALUATING A 1D
C      INTEGRAL OVER A FINITE RANGE
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    15 OCT 1980 (KR)
C
C ARGUMENTS IN
C      IWGHT   DIMENSION OF VECTOR WGHT (.GE.NQP(=2))
C      IABSS   DIMENSION OF VECTOR ABSS (.GE.NQP(=2))
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      WGHT    VECTOR OF DIMENSION IWGHT.  CONTAINS WEIGHTS TO
C              BE USED IN THE 2-POINT QUADRATURE FORMULA
C      ABSS    VECTOR OF DIMENSION IABSS.  CONTAINS ABSCISSAE
C              OF POINTS TO BE USED IN 2-POINT QUADRATURE
C              FORMULA
C      NQP     NUMBER OF QUADRATURE POINTS (=2)
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE QLIN2(WGHT, IWGHT, ABSS, IABSS, NQP, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, IABSS, IERROR, ITEST, IWGHT, NQP
      DOUBLE PRECISION ABSS, SRNAME, WGHT
      DIMENSION ABSS(IABSS), WGHT(IWGHT)
      DATA SRNAME /8H QLIN2  /
      NQP = 2
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IWGHT.LT.NQP) IERROR = 1
                        IF (IABSS.LT.NQP) IERROR = 2
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 ABSS(1) = 1.0D0/DSQRT(3.0D0)
      ABSS(2) = -ABSS(1)
      WGHT(1) = 1.0D0
      WGHT(2) = 1.0D0
      RETURN
      END
