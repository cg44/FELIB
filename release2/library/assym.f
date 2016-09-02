C***********************************************************************
C$SPLIT$ASSYM$*********************************************************
C***********************************************************************
      SUBROUTINE ASSYM(SYSK, ISYSK, JSYSK, ELK, IELK, JELK, STEER,
     *     ISTEER, HBAND, DOFEL, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      FOR REAL SYMMETRIC SYSTEM MATRIX, ADDS THE CONTRIBUTION
C      FROM AN ELEMENT MATRIX
C
C HISTORY
C
C      COPYRIGHT (C) 1979 : SERC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 1.1  29 OCT 1979 (IMS)
C      COMMENTED     7 FEB 1980 (KR)
C
C ARGUMENTS IN
C      SYSK    CONTAINS SYSTEM MATRIX PRIOR TO ADDITION OF
C              CURRENT ELEMENT MATRIX CONTRIBUTION
C      ISYSK   FIRST DIMENSION OF SYSK (.GE. TOTAL NUMBER OF
C              UNCONSTRAINED DEGREES OF FREEDOM)
C      JSYSK   SECOND DIMENSION OF SYSK (.GE. HBAND)
C      ELK     ELEMENT MATRIX
C      IELK    FIRST DIMENSION OF ELK (.GE. DOFEL)
C      JELK    SECOND DIMENSION OF ELK (.GE. DOFEL)
C      STEER   CONTAINS FREEDOM NUMBERS ASSOCIATED WITH ELEMENT
C              MATRIX CONTRIBUTIONS TO SYSTEM MATRIX
C      ISTEER  FIRST DIMENSION OF STEER (.GE. DOFEL)
C      HBAND   SEMI-BANDWIDTH OF SYSTEM MATRIX, INCLUDING
C              DIAGONAL
C      DOFEL   MAXIMUM DEGREES OF FREEDOM ASSOCIATED WITH
C              ELEMENT TYPE
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      SYSK    SYSTEM MATRIX
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE ASSYM(SYSK, ISYSK, JSYSK, ELK, IELK, JELK,
C    *     STEER, ISTEER, HBAND, DOFEL, ITEST)
C***********************************************************************
C
      INTEGER DOFEL, ERRMES, HBAND, I, IELK, IERROR, ISTEER,
     *     ISYSK, ITEST, J, JELK, JSYSK, JTEST, OFFSET, STEER,
     *     STEERI, STEERJ, JTEST
      DOUBLE PRECISION ELK, SRNAME, SYSK
      DIMENSION ELK(IELK,JELK), STEER(ISTEER),
     *     SYSK(ISYSK,JSYSK)
      DATA SRNAME /8H ASSYM  /
C
C     PARAMETER CHECKING
C
      JTEST = ITEST
      IF (JTEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (ISTEER.LT.DOFEL) IERROR = 4
      IF (IELK.LT.DOFEL .OR. JELK.LT.DOFEL) IERROR = 3
      IF (JSYSK.LT.HBAND) IERROR = 2
      IF (HBAND.LE.0 .OR. DOFEL.LE.0) IERROR = 1
      ITEST = ERRMES(JTEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
C     MAIN LOOPS
C
 1010 DO 1040 I=1,DOFEL
         STEERI = STEER(I)
         IF (STEERI.LE.0) GO TO 1040
         DO 1030 J=1,DOFEL
            STEERJ = STEER(J)
            IF (STEERJ.LE.0) GO TO 1030
            OFFSET = STEERJ - STEERI + HBAND
            IF (OFFSET.GT.HBAND) GO TO 1030
C
C     RANGE CHECKING ON STEERI 
C
            IF (JTEST.EQ.-1) GO TO 1020
            IERROR = 0
            IF (ISYSK.LT.STEERI) IERROR = 4
            ITEST = ERRMES(JTEST,IERROR,SRNAME)
            IF (ITEST.NE.0) RETURN
C
 1020       SYSK(STEERI,OFFSET) = SYSK(STEERI,OFFSET) + ELK(I,J)
C
 1030    CONTINUE
 1040 CONTINUE
      RETURN
C
      END
