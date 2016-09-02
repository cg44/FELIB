C***********************************************************************
C$SPLIT$DPLT$*********************************************************
C***********************************************************************
      SUBROUTINE DPLT(D, ID, JD, E, NU, T, NUMSS, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      FORMS 3 BY 3 STRESS-STRAIN MATRIX FOR USE IN PLATE
C      BENDING PROBLEMS
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    13 FEB 1980 (KR)
C
C ARGUMENTS IN
C      ID      FIRST DIMENSION OF D (.GE. 3)
C      JD      SECOND DIMENSION OF D (.GE. 3)
C      E       YOUNG'S MODULUS
C      NU      POISSON'S RATIO
C      T       THICKNESS OF PLATE
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      D       STRESS-STRAIN MATRIX
C      NUMSS   ORDER OF STRESS-STRAIN MATRIX (3)
C
C ROUTINES CALLED
C      MATNUL  ERRMES
C
C
C     SUBROUTINE DPLT(D, ID, JD, E, NU, T, NUMSS, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, ID, IERROR, ITEST, J, JD, NUMSS
      DOUBLE PRECISION D, E, NU, SRNAME, T, T3
      DIMENSION D(ID,JD)
      DATA SRNAME /8H DPLT   /
      NUMSS = 3
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (NU.LT.0.0D0 .OR. NU.GE.0.5D0)
     *                      IERROR = 2
                        IF (ID.LT.NUMSS .OR. JD.LT.NUMSS)
     *                      IERROR = 1
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 CALL MATNUL(D, ID, JD, NUMSS, NUMSS, ITEST)
      T3 = T**3
      D(1,1) = 1.0D0
      D(2,2) = 1.0D0
      D(3,3) = (1.0D0-NU)*0.5D0
      D(1,2) = NU
      D(2,1) = NU
      DO 1030 I=1,3
      DO 1020 J=1,3
      D(I,J) = D(I,J)*T3*E/(12.0D0*(1.0D0-NU*NU))
 1020 CONTINUE
 1030 CONTINUE
      RETURN
      END
