C***********************************************************************
      SUBROUTINE ERUNIT(NUNIT, IKEY, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS THE VALUE OF THE CURRENT ERROR MESSAGE UNIT
C      NUMBER OR SETS THE CURRENT ERROR MESSAGE UNIT NUMBER
C      TO A NEW VALUE
C
C      *********************************************
C      **********MACHINE DEPENDENT ROUTINE**********
C      *********************************************
C
C HISTORY
C
C     COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                          CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C     RELEASE 3.0  18 JAN 1984 (CG)
C     DEBUGGED AND COMMENTED 11/10/1989 (CJC)
C
C LOCALS
C      UNIT     STANDARD UNIT NUMBER FOR MACHINE
C
C ARGUMENTS IN
C      NUNIT   NEW UNIT NUMBER (IKEY = 0)
C      IKEY    ACTION KEY (=0 SET NUNIT TO UNIT, = 1 SET UNIT TO NUNIT)
C      ITEST   STANDARD ERROR CHECK FLAG
C
C ARGUMENTS OUT
C      NUNIT   OLD UNIT NUMBER
C      ITEST: ERROR CONDITION RETURNED OR PRINTED
C       1: INVALID IKEY. 2: NUNIT NOT IN VALID RANGE.
C
C     SUBROUTINE ERUNIT(NUNIT,IKEY,ITEST)
C***********************************************************************
C
      INTEGER IKEY, ITEST, JTEST, NUNIT, RANGE1, RANGE2,
     *     UNIT
      DOUBLE PRECISION SRNAME
C
      DATA SRNAME /8H ERUNIT /
C
C-----------------------------------------------------------------------
      DATA RANGE1 /6/, RANGE2 /10/, UNIT /6/
C-----------------------------------------------------------------------
C
      IF (IKEY .EQ. 1) GO TO 1
C     ENQUIRY
      NUNIT = UNIT
      IF (ITEST .EQ. -1) RETURN
      IF (IKEY .NE. 0) GOTO 2
      ITEST = 0
      RETURN
C     ENQUIRY
C          INVALID KEY
    2 JTEST = 1
C     EXIT SEQUENCE IF ERROR DETECTED
    5 IF (ITEST .NE. 0) GOTO 4
      WRITE (UNIT, 9010) SRNAME, JTEST
      STOP
    4 ITEST = JTEST
      RETURN
C     CHANGE ERROR UNIT
    1 IF ((NUNIT .LT. RANGE1) .OR. (NUNIT .GT. RANGE2)) GOTO 3
      UNIT = NUNIT
      IF (ITEST .EQ. -1) RETURN
      ITEST = 0
      RETURN
C          NUNIT OUTSIDE VALID RANGE
    3 IF (ITEST .EQ. -1) RETURN
      JTEST = 2
      GO TO 5
C     CHANGE ERROR UNIT
 9010 FORMAT (43H ERROR DETECTED BY LEVEL 0 LIBRARY ROUTINE , A8,
     *     11H - ITEST = , I5//)
C
      END
C***********************************************************************
