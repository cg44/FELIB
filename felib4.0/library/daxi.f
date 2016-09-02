C
      SUBROUTINE DAXI(D,ID,JD,E,NU,NUMSS,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      DAXI forms the 4 by 4 stress-strain matrix for use in
C      axisymmetric isotropic problems
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
C      Release 1.1  29 Oct 1979 (IMS)
C      Commented    13 Feb 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      ID      first dimension of matrix D (.GE. 4)
C      JD      second dimension of D (.GE. 4)
C      E       young's modulus
C      NU      poisson's ratio
C      ITEST   error checking option
C
C ARGUMENTS out
C      D       stress-strain matrix
C      NUMSS   order of stress-strain matrix, 4 in this case
C
C ROUTINES called
C      MATNUL  ERRMES
C
C     SUBROUTINE DAXI(D,ID,JD,E,NU,NUMSS,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,ID,IERROR,ITEST,J,JD,NUMSS
      CHARACTER*6 SRNAME
      DOUBLE PRECISION D,E,NU,NU1,NUNU
      DIMENSION D(ID,JD)
      DATA SRNAME/'DAXI'/
C
C     Setting NUMSS
C
      NUMSS = 4
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
C     Main loops
C
C     Initialise D matrix
C
      CALL MATNUL(D,ID,JD,NUMSS,NUMSS,ITEST)
C
      NU1 = NU/ (1.0D0-NU)
      NUNU = 0.5D0* (1.0D0-2.0D0*NU)/ (1.0D0-NU)
C
      D(1,1) = 1.0D0
      D(2,2) = 1.0D0
      D(3,3) = NUNU
      D(4,4) = 1.0D0
      D(1,2) = NU1
      D(2,1) = NU1
      D(1,4) = NU1
      D(4,1) = NU1
      D(2,4) = NU1
      D(4,2) = NU1
C
      DO 1010 I = 1,4
         DO 1000 J = 1,4
            D(I,J) = D(I,J)*E* (1.0D0-NU)/
     *               ((1.0D0-2.0D0*NU)* (1.0D0+NU))
 1000    CONTINUE
 1010 CONTINUE
C
      END
