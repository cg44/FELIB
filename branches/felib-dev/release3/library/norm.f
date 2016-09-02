C***********************************************************************
      DOUBLE PRECISION FUNCTION NORM(RHS, IRHS, TOTDOF, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      THIS ROUTINE COMPUTES AN L2 NORM OF A VECTOR FOR USE
C      IN TERMINATING NON-LINEAR ITERATIONS
C
C METHOD
C      SQUARE ROOT OF THE SUMS OF SQUARES
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 3.0  29 JUN 1986 (CJH,CG)
C
C ARGUMENTS IN
C      RHS     VECTOR CONTAINING VALUES
C      IRHS    DIMENSION OF VECTOR RHS (IRHS.GE.TOTDOF)
C      TOTDOF  NUMBER OF ENTRIES IN RHS
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      NORM    THE VALUE OF THE NORM (FUNCTION VALUE)
C
C ROUTINES CALLED
C      ERRMES
C
C
C     DOUBLE PRECISION FUNCTION NORM(RHS, IRHS, TOTDOF, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IERROR, IRHS, ITEST, TOTDOF
      DOUBLE PRECISION RHS, SRNAME
      DIMENSION RHS(IRHS)
      DATA SRNAME /8H NORM   /
C
C     PARAMETER CHECKING
C
      IF (ITEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (TOTDOF.GT.IRHS) IERROR = 1
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
C                            COMPUTE NORM
C
 1010 NORM = 0.D0
      DO 1020 I=1,TOTDOF
      NORM = NORM + RHS(I)**2
 1020 CONTINUE
      NORM = DSQRT(NORM)
C
      RETURN
      END
C***********************************************************************
