C
      SUBROUTINE IASSYM(SYSK,ISYSK,JSYSK,ELK,IELK,JELK,STEER,ISTEER,
     *                  HBAND,DOFEL,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      IASSYM for real symmetric system matrix, adds the contribution
C      from an element matrix into real part of complex SYSK
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
C      Release 2.0  29 Oct 1985 (CRIE)
C      Commented     7 Nov 1985 (CG)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      SYSK    contains system matrix prior to addition of
C              current element matrix contribution
C      ISYSK   first dimension of SYSK (.GE. total number of
C              unconstrained degrees of freedom)
C      JSYSK   second dimension of SYSK (.GE. HBAND)
C      ELK     element matrix
C      IELK    first dimension of ELK (.GE. DOFEL)
C      JELK    second dimension of ELK (.GE. DOFEL)
C      STEER   contains freedom numbers associated with element
C              matrix contributions to system matrix
C      ISTEER  first dimension of STEER (.GE. DOFEL)
C      HBAND   semi-bandwidth of system matrix, including
C              diagonal
C      DOFEL   maximum degrees of freedom associated with
C              element type
C      ITEST   error checking option
C
C ARGUMENTS out
C      SYSK    system matrix -   ordered pairs
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE IASSYM(SYSK,ISYSK,JSYSK,ELK,IELK,JELK,STEER,ISTEER,
C    *                  HBAND,DOFEL,ITEST)
C***********************************************************************
C
      INTEGER CD,DOFEL,ERRMES,HBAND,I,IELK,IERROR,ISTEER,ISYSK,ITEST,J,
     *        JELK,JSYSK,STEER,STEERI,JTEST
      CHARACTER*6 SRNAME
      DOUBLE PRECISION ELK,SYSK
      DIMENSION ELK(IELK,JELK),STEER(ISTEER),SYSK(2,ISYSK,JSYSK)
C
      EXTERNAL ERRMES
C
      DATA SRNAME/'IASSYM'/
C
C     Parameter checking
C
      JTEST = ITEST
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (ISTEER.LT.DOFEL) IERROR = 4
         IF (IELK.LT.DOFEL .OR. JELK.LT.DOFEL) IERROR = 3
         IF (JSYSK.LT.HBAND) IERROR = 2
         IF (HBAND.LE.0 .OR. DOFEL.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
      DO 1010 I = 1,DOFEL
         IF (STEER(I).NE.0) THEN
            DO 1000 J = 1,DOFEL
               IF (STEER(J).NE.0) THEN
                  CD = STEER(J) - STEER(I) + HBAND
                  IF (CD.LE.HBAND) THEN
                     STEERI = STEER(I)
                     IF (JTEST.NE.-1) THEN
                        IERROR = 0
                        IF (ISYSK.LT.STEERI) IERROR = 5
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
                     END IF
                     SYSK(2,STEERI,CD) = SYSK(2,STEERI,CD) + ELK(I,J)
                  END IF
               END IF
 1000       CONTINUE
         END IF
 1010 CONTINUE
C
      END
