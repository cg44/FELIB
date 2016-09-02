C***********************************************************************
C$SPLIT$FEMONIT$*********************************************************
C***********************************************************************
C********************************************************************
C
C FINITE ELEMENT LIBRARY MONITORING ROUTINES
C
C
      SUBROUTINE FEFIN
      CHARACTER*20 LOGFIL
      INTEGER*4 ERROUT
      INTEGER*4 MONOUT,LENFIL,ISTREC,MNCUR, LOGGED
      COMMON/CFELIB/MONOUT,LOGGED,LENFIL,ISTREC,MNCUR,
     *              LOGFIL
C
      DATA ERROUT /1/
C
C     CHECK IF ONLY TEST
C
      IF(LOGGED.EQ.6666) RETURN
C
C     THE REAL THING - CONTINUE
C
      IF(LOGGED.EQ.9999) THEN
         CALL MNSTOP(MONOUT,ERROUT,LOGFIL,LENFIL,ISTREC,MNCUR)
C
C     SET FLAG
C
         LOGGED=7777
      ELSE
        WRITE(ERROUT,9010)
        STOP
      ENDIF
      RETURN
 9010 FORMAT('ERROR DETECTED - FEINIT MUST BE CALLED FIRST')
      END
