C
      SUBROUTINE ASFULG(SYSK,ISYSK,JSYSK,ELK,IELK,JELK,STR1,ISTR1,STR2,
     *                  ISTR2,DOFEL1,DOFEL2,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      ASFULG performs the general assembly of an element matrix
C      into a full system matrix
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
C      Release 3.0    1 Jul 1996 (CG)
C      Release 4.0    2 Oct 2003 (CG)
C
C ARGUMENTS in
C      SYSK    contains system matrix prior to addition of
C              current elemnt matrix contribution
C      ISYSK   first dimension of SYSK (.GE. total number of
C              unconstrained degrees of freedom)
C      JSYSK   second dimension of SYSK (.GE. total number of
C              unconstrained degrees of freedom)
C      ELK     element matrix
C      IELK    first dimension of ELK (.GE. DOFEL1)
C      JELK    second dimension of ELK (.GE. DOFEL1)
C      STR1    contains freedom numbers associated with element
C              matrix contributions to system matrix
C      ISTR2   dimension of STR1 (.GE. DOFEL1)
C      STR1    contains freedom numbers associated with element
C              matrix contributions to system matrix
C      ISTR2   dimension of STR2 (.GE. DOFEL2)
C      DOFEL1  maximum number of degrees of freedom
C      DOFEL2  maximum number of degrees of freedom
C      ITEST   error checking option
C
C ARGUMENTS out
C      SYSK    system matrix
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE ASFULG(SYSK,ISYSK,JSYSK,ELK,IELK,JELK,STR1,ISTR1,STR2,
C    *                  ISTR2,DOFEL1,DOFEL2,ITEST)
C***********************************************************************
C
      INTEGER DOFEL1,DOFEL2,IELK,IERROR,ISTR1,ISTR2,ISYSK,
     *        ITEST,JELK,JSYSK
      INTEGER STR1(ISTR1),STR2(ISTR2)
      DOUBLE PRECISION ELK(IELK,JELK),SYSK(ISYSK,JSYSK)
C
      INTEGER I,J,JTEST,STEERI,STEERJ
      CHARACTER*6 SRNAME
C
      INTEGER ERRMES
      EXTERNAL ERRMES
C
      DATA SRNAME/'ASFULG'/
C
C     Parameter checking
C
      JTEST = ITEST
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (ISTR1.LT.DOFEL1 .OR. ISTR2.LT.DOFEL2) IERROR = 3
         IF (IELK.LT.DOFEL1 .OR. JELK.LT.DOFEL2) IERROR = 2
         IF (DOFEL1.LE.0 .OR. DOFEL2.LE.0) IERROR = 1
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Assembly loops
C
      DO 1010 I = 1,DOFEL1
         STEERI = STR1(I)
         IF (STEERI.GT.0) THEN
            DO 1000 J = 1,DOFEL2
               STEERJ = STR2(J)
               IF (STEERJ.GT.0) THEN
C
C     Range checking on STEERI and STEERJ
C
                  IF (JTEST.NE.-1) THEN
                     IERROR = 0
                     IF (ISYSK.LT.STEERI .OR.
     *                   JSYSK.LT.STEERJ) IERROR = 4
                     ITEST = ERRMES(JTEST,IERROR,SRNAME)
                     IF (ITEST.NE.0) RETURN
                  END IF
C
                  SYSK(STEERI,STEERJ) = SYSK(STEERI,STEERJ) + ELK(I,J)
               END IF
 1000       CONTINUE
C
         END IF
 1010 CONTINUE
C
      END
