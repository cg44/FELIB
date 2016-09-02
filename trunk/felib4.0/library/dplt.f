C
      SUBROUTINE DPLT(D,ID,JD,E,NU,T,NUMSS,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      DPLT forms 3 by 3 stress-strain matrix for use in plate
C      bending problems
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
C      Commented    13 Feb 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      ID      first dimension of D (.GE. 3)
C      JD      second dimension of D (.GE. 3)
C      E       young's modulus
C      NU      poisson's ratio
C      T       thickness of plate
C      ITEST   error checking option
C
C ARGUMENTS out
C      D       stress-strain matrix
C      NUMSS   order of stress-strain matrix (3)
C
C ROUTINES called
C      MATNUL  ERRMES
C
C     SUBROUTINE DPLT(D,ID,JD,E,NU,T,NUMSS,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,ID,IERROR,ITEST,J,JD,NUMSS
      CHARACTER*6 SRNAME
      DOUBLE PRECISION D,E,NU,T,T3
      DIMENSION D(ID,JD)
      DATA SRNAME/'DPLT'/
C
C     Setting NUMSS
C
      NUMSS = 3
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (NU.LT.0.0D0 .OR. NU.GE.0.5D0) IERROR = 2
         IF (ID.LT.NUMSS .OR. JD.LT.NUMSS) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
C     Initialise matrix D
C
      CALL MATNUL(D,ID,JD,NUMSS,NUMSS,ITEST)
C
      T3 = T**3
      D(1,1) = 1.0D0
      D(2,2) = 1.0D0
      D(3,3) = (1.0D0-NU)*0.5D0
      D(1,2) = NU
      D(2,1) = NU
C
      DO 1010 I = 1,3
         DO 1000 J = 1,3
            D(I,J) = D(I,J)*T3*E/ (12.0D0* (1.0D0-NU*NU))
 1000    CONTINUE
 1010 CONTINUE
C
      END
