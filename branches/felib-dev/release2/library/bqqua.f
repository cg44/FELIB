C***********************************************************************
C$SPLIT$BQQUA$*********************************************************
C***********************************************************************
      SUBROUTINE BQQUA(ABSS, IABSS, JABSS, WORK, IWORK, NQP,
     *     SIDNUM, COEF, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      FORMS AN EQUIVALENT ONE-DIMENSIONAL QUADRATURE RULE OF
C      A GIVEN TWO-DIMENSIONAL RULE FOR INTEGRATION ALONG THE
C      SPECIFIED SIDE OF A RECTANGULAR ELEMENT
C
C HISTORY
C
C      COPYRIGHT (C) 1979 : SERC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 2.0  1 FEB 1981 (CG) 
C      COMMENTED   24 JUN 1981 (CG)
C
C ARGUMENTS IN
C      ABSS    ARRAY HOLDING THE ABSCISSAE OF THE TWO-
C              DIMENSIONAL QUADRATURE RULE TO BE USED
C      IABSS   FIRST DIMENSION OF ARRAY ABSS
C      JABSS   SECOND DIMENSION OF ARRAY ABSS
C      IWORK   DIMENSION OF ARRAY WORK
C      NQP     NUMBER OF QUADRATURE POINTS IN THE RULE
C      SIDNUM  THE SIDE NUMBER FOR WHICH THE ONE-
C              DIMENSIONAL RULE IS REQUIRED
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      WORK    ARRAY CONTAINING THE ABSCISSAE OF THE
C              EQUIVALENT ONE-DIMENSIONAL RULE
C      COEF    A MULTIPLIER OF THE RULE
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE BQQUA(ABSS, IABSS, JABSS, WORK, IWORK, NQP,
C    *     SIDNUM, COEF, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IABSS, IERROR, ITEST, IWORK, J, JABSS,
     *     K, NQP, SIDNUM
      DOUBLE PRECISION ABSS, COEF, SRNAME, VAL, WORK
      DIMENSION ABSS(IABSS,JABSS), WORK(IWORK)
      DATA SRNAME /8H BQQUA  /
C
C     PARAMETER CHECKING
C
      IF (ITEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (IABSS.LT.2 .OR. JABSS.LT.NQP) IERROR = 4
      IF (IWORK.LT.NQP) IERROR = 3
      IF (SIDNUM.LE.0 .OR. SIDNUM.GT.4) IERROR = 2
      IF (NQP.LE.0) IERROR = 1
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
C     MAIN LOOPS
C
 1010 COEF = 1.0D0
      GO TO (1020, 1030, 1040, 1050), SIDNUM
C
C     SIDE NUMBER 1 : XI = -1
C
 1020 I = 1
      J = 2
      VAL = -1.0D0
      GO TO 1060
C
C     SIDE NUMBER 2 : XI = +1
C
 1030 I = 2
      J = 1
      VAL = 1.0D0
      GO TO 1060
C
C     SIDE NUMBER 3 : ETA = -1
C
 1040 I = 1
      J = 2
      VAL = 1.0D0
      GO TO 1060
C
C     SIDE NUMBER 4 : ETA = +1
C
 1050 I = 2
      J = 1
      VAL = -1.0D0
C
 1060 DO 1070 K=1,NQP
         ABSS(I,K) = VAL
         ABSS(J,K) = WORK(K)
 1070 CONTINUE
C
      RETURN
      END
