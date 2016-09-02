C
      SUBROUTINE ASFUL(SYSK,ISYSK,JSYSK,ELK,IELK,JELK,STEER,ISTEER,
     *                 DOFEL,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      ASFUL assembles full real system matrix
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
C      Release 1.1   29 Oct 1978 (CG)
C      Commented      6 Feb 1980 (KR)
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
C      IELK    first dimension of ELK (.GE. DOFEL)
C      JELK    second dimension of ELK (.GE. DOFEL)
C      STEER   contains freedom numbers associated with element
C              matrix contributions to system matrix
C      ISTEER  dimension of STEER (.GE. DOFEL)
C      DOFEL   maximum number of degrees of freedom associated
C              with element type
C      ITEST   error checking option
C
C ARGUMENTS out
C      SYSK    system matrix
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE ASFUL(SYSK,ISYSK,JSYSK,ELK,IELK,JELK,STEER,ISTEER,
C    *                 DOFEL,ITEST)
C***********************************************************************
C
      INTEGER DOFEL,IELK,ISTEER,ISYSK,ITEST,JELK,JSYSK,STEERI,STEERJ
      INTEGER STEER(ISTEER)
      DOUBLE PRECISION ELK(IELK,JELK),SYSK(ISYSK,JSYSK)
C
      INTEGER I,IERROR,J,JTEST
      CHARACTER*6 SRNAME
C
      INTEGER ERRMES
      EXTERNAL ERRMES
C
      DATA SRNAME/'ASFUL'/
C
C     Parameter checking
C
      JTEST = ITEST
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (ISTEER.LT.DOFEL) IERROR = 3
         IF (IELK.LT.DOFEL .OR. JELK.LT.DOFEL) IERROR = 2
         IF (DOFEL.LE.0) IERROR = 1
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main loops
C
      DO 1010 I = 1,DOFEL
         STEERI = STEER(I)
         IF (STEERI.GT.0) THEN
            DO 1000 J = 1,DOFEL
               STEERJ = STEER(J)
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
