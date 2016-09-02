C
      INTEGER FUNCTION ERRMES(ITEST,IERROR,SRNAME)
C-----------------------------------------------------------------------
C PURPOSE
C      ERRMES returns the value of IERROR or terminates the program,
C      printing a failure message
C
C HISTORY
C
C      Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
C                           Chilton, Didcot, Oxfordshire OX11 0QX
C
C      Contact: Prof Chris Greenough
C      Email: c.greenough@rl.ac.uk
C      Tel: (44) 1235 445307
C
C      Release 1.1  29 Oct 1979 (CG)
C      Commented    14 Feb 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 1996 (CG)
C
C ARGUMENTS in
C      ITEST   contains either 0 (hard fail) or 1 (soft fail).
C              any other entry gives hard fail.
C      IERROR  contains the number of the detected error
C      SRNAME  contains up to 8 characters - usually a library
C              routine name
C
C ARGUMENTS out
C      ERRMES  routine name, contains the value of IERROR
C
C ROUTINES called
C      can call auxiliary routine in some versions of library
C      ERUNIT and ADUNIT
C
C     INTEGER FUNCTION ERRMES(ITEST,IERROR,SRNAME)
C***********************************************************************
C
      INTEGER IERROR,ITEST,JTEST,UNIT,GET
      CHARACTER*6 SRNAME
C
      DATA GET/0/,JTEST/1/
C
C     Hard failure
C
      IF (ITEST.EQ.-99) THEN
C
C     To return Release message
C
         CALL ADUNIT(UNIT,GET,JTEST)
         WRITE (UNIT,FMT=9980)
C
      ELSE IF (ITEST.EQ.1 .OR. IERROR.EQ.0) THEN
C
C     Soft failure
C
         ERRMES = IERROR
         IF (ITEST.NE.0 .AND. IERROR.NE.0) THEN
            CALL ADUNIT(UNIT,GET,JTEST)
            WRITE (UNIT,FMT=9990) SRNAME,IERROR
         END IF
C
      ELSE
C
         CALL ERUNIT(UNIT,GET,JTEST)
         WRITE (UNIT,FMT=9990) SRNAME,IERROR
         STOP
C
      END IF
 9990 FORMAT (' ERROR DETECTED BY LEVEL 0 LIBRARY ROUTINE ',A,
     *       ' - ITEST = ',I5,//)
 9980 FORMAT (' RELEASE 3.0  -  1 JAN 87')
      END
