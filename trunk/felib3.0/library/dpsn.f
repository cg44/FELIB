C***********************************************************************
      SUBROUTINE DPSN(D, ID, JD, E, NU, NUMSS, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      FORMS THE 3 BY 3 STRESS-STRAIN MATRIX FOR USE IN PLANE
C      STRAIN PROBLEMS
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 1.1  29 OCT 1980 (IMS)
C      COMMENTED    14 FEB 1980 (KR)
C
C ARGUMENTS IN
C      ID      FIRST DIMENSION OF ARRAY D (.GE. 3)
C      JD      SECOND DIMENSION OF D (.GE. 3)
C      E       YOUNG'S MODULUS
C      NU      POISSON'S RATIO
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      D       STRESS-STRAIN MATRIX FOR PLANE STRAIN
C      NUMSS   ORDER OF STRESS-STRAIN MATRIX (3)
C
C ROUTINES CALLED
C      MATNUL  ERRMES
C
C
C     SUBROUTINE DPSN(D, ID, JD, E, NU, NUMSS, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, ID, IERROR, ITEST, J, JD, NUMSS
      DOUBLE PRECISION D, E, NU, NU1, NUNU, SRNAME
      DIMENSION D(ID,JD)
      DATA SRNAME /8H DPSN   /
C
C     PARAMETER CHECKING
C
      NUMSS = 3
      IF (ITEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (NU.LT.0.0D0 .OR. NU.GE.0.5D0) IERROR = 2
      IF (ID.LT.NUMSS .OR. JD.LT.NUMSS) IERROR = 1
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
C     MAIN BODY
C
C
C     INITIALISE MATRIX D
C
 1010 CALL MATNUL(D, ID, JD, NUMSS, NUMSS, ITEST)
C
      NU1 = NU/(1.0D0-NU)
      NUNU = (0.5D0-NU)/(1.0D0-NU)
C
      D(1,1) = 1.0D0
      D(2,2) = 1.0D0
      D(3,3) = NUNU
      D(1,2) = NU1
      D(2,1) = NU1
C
      DO 1030 I=1,3
         DO 1020 J=1,3
            D(I,J) = D(I,J)*E/(2.0D0*(1.0D0+NU)*NUNU)
 1020    CONTINUE
 1030 CONTINUE
C
      RETURN
      END
C***********************************************************************
