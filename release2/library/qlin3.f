C***********************************************************************
C$SPLIT$QLIN3$*********************************************************
C***********************************************************************
      SUBROUTINE QLIN3(WGHT, IWGHT, ABSS, IABSS, NQP, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS WEIGHTS AND ABSCISSAE OF A 3-POINT GAUSS-
C      LEGENDRE QUADRATURE FORMULA FOR USE IN EVALUATING A 1D
C      INTEGRAL OVER A FINITE RANGE
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    15 OCT 1980 (KR)
C
C ARGUMENTS IN
C      IWGHT   DIMENSION OF VECTOR WGHT (.GE.NQP(=3))
C      IABSS   DIMENSION OF VECTOR ABSS (.GE.NQP(=3))
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      WGHT    VECTOR OF DIMENSION IWGHT.  CONTAINS WEIGHTS FOR
C              3-POINT FORMULA
C      ABSS    VECTOR OF DIMENSION IABSS.  CONTAINS ABSCISSAE
C              OF POINTS FOR USE IN 3-POINT FORMULA
C      NQP     NUMBER OF QUADRATURE POINTS (=3)
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE QLIN3(WGHT, IWGHT, ABSS, IABSS, NQP, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, IABSS, IERROR, ITEST, IWGHT, NQP
      DOUBLE PRECISION ABSS, SRNAME, WGHT
      DIMENSION ABSS(IABSS), WGHT(IWGHT)
      DATA SRNAME /8H QLIN3  /
      NQP = 3
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IWGHT.LT.NQP) IERROR = 1
                        IF (IABSS.LT.NQP) IERROR = 2
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 ABSS(1) = 0.2D0*DSQRT(15.0D0)
      ABSS(2) = 0.0D0
      ABSS(3) = -ABSS(1)
      WGHT(1) = 5.0D0/9.0D0
      WGHT(3) = WGHT(1)
      WGHT(2) = 8.0D0/9.0D0
      RETURN
      END
