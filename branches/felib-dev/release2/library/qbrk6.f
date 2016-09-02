C***********************************************************************
C$SPLIT$QBRK6$*********************************************************
C***********************************************************************
      SUBROUTINE QBRK6(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS WEIGHTS AND ABSCISSAE OF A SIX-POINT GAUSS TYPE
C      QUADRATURE RULE FOR USE IN EVALUATING THE INTEGRAL OF A
C      3D FUNCTION OVER A CUBE
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    15 OCT 1980 (KR)
C
C ARGUMENTS IN
C      IWGHT   DIMENSION OF VECTOR WGHT (.GE.NQP(=6))
C      IABSS   FIRST DIMENSION OF ARRAY ABSS (.GE.3)
C      JABSS   SECOND DIMENSION OF ARRAY ABSS (.GE.NQP(=6))
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      WGHT    VECTOR OF DIMENSION IWGHT.  CONTAINS WEIGHTS TO
C              BE USED IN THE 6-POINT QUADRATURE FORMULA
C      ABSS    ARRAY OF DIMENSION (IABSS,JABSS).  CONTAINS
C              ABSCISSAE OF POINTS TO BE USED IN QUADRATURE
C              FORMULA
C      NQP     NUMBER OF QUADRATURE POINTS TO BE USED (=6)
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE QBRK6(WGHT, IWGHT, ABSS, IABSS, JABSS, NQP, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IABSS, IERROR, ITEST, IWGHT, JABSS,
     *     NQP
      DOUBLE PRECISION ABSS, SRNAME, W, WGHT
      DIMENSION ABSS(IABSS,JABSS), WGHT(IWGHT)
      DATA SRNAME /8H QBRK6  /
      NQP = 6
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IWGHT.LT.NQP) IERROR = 1
                        IF (IABSS.LT.3 .OR. JABSS.LT.NQP)
     *                      IERROR = 2
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 W = 8.0D0/6.0D0
      DO 1020 I=1,6
      WGHT(I) = W
      ABSS(1,I) = 0.0D0
      ABSS(2,I) = 0.0D0
      ABSS(3,I) = 0.0D0
 1020 CONTINUE
      ABSS(1,5) = -1.0D0
      ABSS(1,6) = 1.0D0
      ABSS(2,3) = -1.0D0
      ABSS(2,4) = 1.0D0
      ABSS(3,1) = -1.0D0
      ABSS(3,2) = 1.0D0
      RETURN
      END
