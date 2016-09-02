C
      SUBROUTINE ADUNIT(NUNIT,IKEY,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      ADUNIT returns the value of the current advisory message UNIT
C      number or sets the current advisory message UNIT number
C      to a new value
C
C      *********************************************
C      **********MACHINE DEPENDENT ROUTINE**********
C      *********************************************
C
C HISTORY
C
C      Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
C                           Chilton, Didcot, Oxfordshire OX11 0QX
C
C      Release 3.0  18 Jan 1984 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      NUNIT   new UNIT number (IKEY = 0)
C      IKEY    action key (=0 set NUNIT to UNIT, = 1 set UNIT to NUNIT)
C      ITEST   standard error check flag
C
C ARGUMENTS out
C      NUNIT   old UNIT number
C      ITEST   error conditions returned or printed
C              - 1: invalid IKEY.
C              - 2: NUNIT not in valid range.
C
C***********************************************************************
C
      INTEGER IKEY,ITEST,JTEST,NUNIT,RANGE1,RANGE2,UNIT
      CHARACTER*6 SRNAME
C
      SAVE UNIT
C
      DATA SRNAME/'ADUNIT'/
C
C-----------------------------------------------------------------------
      DATA RANGE1/5/,RANGE2/10/,UNIT/6/
C-----------------------------------------------------------------------
C
      IF (IKEY.NE.1) THEN
C
C     Enquiry
C
         NUNIT = UNIT
         IF (ITEST.EQ.-1) THEN
            RETURN
         ELSE IF (IKEY.NE.0) THEN
            JTEST = 1
         ELSE
            ITEST = 0
            RETURN
         END IF
C
C     Change advisory UNIT
C
      ELSE IF ((NUNIT.GE.RANGE1) .AND. (NUNIT.LE.RANGE2)) THEN
         UNIT = NUNIT
         IF (ITEST.NE.-1) ITEST = 0
         RETURN
C
C     NUNIT outside valid range
C
      ELSE IF (ITEST.EQ.-1) THEN
         RETURN
      ELSE
         JTEST = 2
      END IF
C
C     Exit sequence if error is detected
C
      IF (ITEST.NE.0) THEN
         ITEST = JTEST
      ELSE
         WRITE (UNIT,FMT=9990) SRNAME,JTEST
         STOP
      END IF
C
 9990 FORMAT (' ERROR DETECTED BY LEVEL 0 LIBRARY ROUTINE ',A8,
     *       ' - ITEST = ',I5,//)
      END
