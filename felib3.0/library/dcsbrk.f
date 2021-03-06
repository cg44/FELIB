C***********************************************************************
      SUBROUTINE DCSBRK(JACIN, IJACIN, JJACIN, FACNUM, COSIN,
     *     ICOSIN, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CALCULATES THE DIRECTION COSINES OF THE OUTWARD NORMAL
C      TO THE SPECIFIED FACE OF A HEXAHEDRAL ELEMENT GIVEN THE
C      JACOBIAN INVERSE
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 2.0  1 FEB 1981 (CG)
C      COMMENTED    1 FEB 1981 (CG)
C
C ARGUMENTS IN
C      JACIN   ARRAY OF DIMENSION (IJACIN,JJACIN). CONTAINS THE
C              INVERSE OF THE TRANSFORMATION JACOBIAN.
C      IJACIN  FIRST DIMENSION OF JACIN (.GE.3)
C      JJACIN  SECOND DIMENSION OF JACIN (.GE.3)
C      FACNUM  THE SIDE NUMBER FOR WHICH THE OUTWARD NORMAL
C              IS REQUIRED (.LE.6)
C      ICOSIN  DIMENSION OF VECTOR COSIN (.GE.3)
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      COSIN   VECTOR OF DIMENSION ICOSIN. CONTAINS THE DIR-
C              ECTION COSINES OF THE OUTWARD NORMAL
C
C ROUTINES CALLED
C      ERRMES  MATVEC
C
C
C     SUBROUTINE DCSBRK(JACIN, IJACIN, JJACIN, FACNUM, COSIN,
C    *     ICOSIN, ITEST)
C***********************************************************************
C
      INTEGER DIMEN, ERRMES, FACNUM, ICOSIN, IERROR, IJACIN,
     *     ITEST, IWORK, JJACIN
      DOUBLE PRECISION AMOD, COSIN, JACIN, SRNAME, WORK
      DIMENSION COSIN(ICOSIN), JACIN(IJACIN,JJACIN), WORK(3)
      DATA DIMEN /3/, IWORK /3/, SRNAME /8H DCSBRK /
C
C     PARAMETER CHECKING
C
      IF (ITEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (IJACIN.LT.3 .OR. JJACIN.LT.3) IERROR = 3
      IF (ICOSIN.LT.3) IERROR = 2
      IF (FACNUM.LE.0 .OR. FACNUM.GT.6) IERROR = 1
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
C     MAIN BODY
C
 1010 GO TO (1060, 1070, 1020, 1050, 1030, 1040), FACNUM
C
C     FACE NUMBER 3 : XI = -1
C
 1020 WORK(1) = -1.0D0
      WORK(2) = 0.0D0
      WORK(3) = 0.0D0
      GO TO 1080
C
C     FACE NUMBER 5 : XI = +1
C
 1030 WORK(1) = 1.0D0
      WORK(2) = 0.0D0
      WORK(3) = 0.0D0
      GO TO 1080
C
C     FACE NUMBER 6 : ETA = -1
C
 1040 WORK(1) = 0.0D0
      WORK(2) = -1.0D0
      WORK(3) = 0.0D0
      GO TO 1080
C
C     FACE NUMBER 4 : ETA = +1
C
 1050 WORK(1) = 0.0D0
      WORK(2) = 1.0D0
      WORK(3) = 0.0D0
      GO TO 1080
C
C     FACE NUMBER 1 : ZETA = -1
C
 1060 WORK(1) = 0.0D0
      WORK(2) = 0.0D0
      WORK(3) = -1.0D0
      GO TO 1080
C
C     FACE NUMBER 2 : ZETA = +1
C
 1070 WORK(1) = 0.0D0
      WORK(2) = 0.0D0
      WORK(3) = 1.0D0
C
C     CALCULATE DIRECTION COSINES
C
 1080 CALL MATVEC(JACIN, IJACIN, JJACIN, WORK, IWORK, DIMEN,
     *     DIMEN, COSIN, ICOSIN, ITEST)
      AMOD = DSQRT(COSIN(1)*COSIN(1)+COSIN(2)*COSIN(2)+COSIN(3)*
     *     COSIN(3))
      COSIN(1) = COSIN(1)/AMOD
      COSIN(2) = COSIN(2)/AMOD
      COSIN(3) = COSIN(3)/AMOD
C
      RETURN
      END
C***********************************************************************
