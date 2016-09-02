C***********************************************************************
      SUBROUTINE ASSYMG(SYSK, ISYSK, JSYSK, ELK, IELK, JELK, STR1,
     *     ISTR1, STR2, ISTR2, HBAND, DOFEL1, DOFEL2, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      FOR REAL SYMMETRIC SYSTEM MATRIX, ADDS THE CONTRIBUTION
C      FROM AN GENERAL ELEMENT MATRIX
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 3.0  27 MAR 1983 (CG)
C
C ARGUMENTS IN
C      SYSK    CONTAINS SYSTEM MATRIX PRIOR TO ADDITION OF
C              CURRENT ELEMENT MATRIX CONTRIBUTION
C      ISYSK   FIRST DIMENSION OF SYSK (.GE. TOTAL NUMBER OF
C              UNCONSTRAINED DEGREES OF FREEDOM)
C      JSYSK   SECOND DIMENSION OF SYSK (.GE. HBAND)
C      ELK     ELEMENT MATRIX
C      IELK    FIRST DIMENSION OF ELK (.GE. DOFEL1)
C      JELK    SECOND DIMENSION OF ELK (.GE. DOFEL2)
C      STR1    CONTAINS FREEDOM NUMBERS ASSOCIATED WITH ELEMENT
C              MATRIX CONTRIBUTIONS TO SYSTEM MATRIX
C      ISTR1   FIRST DIMENSION OF STR1 (.GE. DOFEL1)
C      STR2    CONTAINS FREEDOM NUMBERS ASSOCIATED WITH ELEMENT
C              MATRIX CONTRIBUTIONS TO SYSTEM MATRIX
C      ISTR2   FIRST DIMENSION OF STR2 (.GE. DOFEL2)
C      HBAND   SEMI-BANDWIDTH OF SYSTEM MATRIX, INCLUDING
C              DIAGONAL
C      DOFEL1  MAXIMUM DEGREES OF FREEDOM 1
C      DOFEL2  MAXIMUM DEGREES OF FREEDOM 2
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      SYSK    SYSTEM MATRIX
C
C ROUTINES CALLED
C      ERRMES
C
C     SUBROUTINE ASSYMG(SYSK, ISYSK, JSYSK, ELK, IELK, JELK, STR1,
C    *     ISTR1, STR2, ISTR2, HBAND, DOFEL1, DOFEL2, ITEST)
C***********************************************************************
C
      INTEGER DOFEL1, DOFEL2, ERRMES, HBAND, I, IELK, IERROR,
     *     ISTR1, ISTR2, ISYSK, ITEST, J, JELK, JSYSK, JTEST,
     *     OFFSET, STEERI, STEERJ, STR1, STR2
      DOUBLE PRECISION ELK, SRNAME, SYSK
      DIMENSION ELK(IELK,JELK), STR1(ISTR1), STR2(ISTR2),
     *     SYSK(ISYSK,JSYSK)
      DATA SRNAME /8H ASSYMG /
C
C     PARAMETER CHECKING
C
      JTEST = ITEST
      IF (JTEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (ISTR1.LT.DOFEL1 .OR. ISTR2.LT.DOFEL2) IERROR = 4
      IF (IELK.LT.DOFEL1 .OR. JELK.LT.DOFEL2) IERROR = 3
      IF (JSYSK.LT.HBAND) IERROR = 2
      IF (HBAND.LE.0 .OR. DOFEL1.LE.0 .OR. DOFEL2.LE.0) IERROR = 1
      ITEST = ERRMES(JTEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
C     ASSEMBLY LOOPS
C
 1010 DO 1040 I=1,DOFEL1
         STEERI = STR1(I)
         IF (STEERI.LE.0) GO TO 1040
         DO 1030 J=1,DOFEL2
            STEERJ = STR2(J)
            IF (STEERJ.LE.0) GO TO 1030
            OFFSET = STEERJ - STEERI + HBAND
            IF (OFFSET.GT.HBAND) GO TO 1030
C
C           RANGE CHECKING
C
            IF (JTEST.EQ.-1) GO TO 1020
            IERROR = 0
            IF (ISYSK.LT.STEERI) IERROR = 5
            ITEST = ERRMES(JTEST,IERROR,SRNAME)
            IF (ITEST.NE.0) RETURN
C
 1020       SYSK(STEERI,OFFSET) = SYSK(STEERI,OFFSET) + ELK(I,J)
C
 1030    CONTINUE
 1040 CONTINUE
C
      RETURN
      END
C***********************************************************************
