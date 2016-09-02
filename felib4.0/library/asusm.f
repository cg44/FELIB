C
      SUBROUTINE ASUSM(SYSK,ISYSK,JSYSK,ELK,IELK,JELK,STEER,ISTEER,
     *                 HBAND,DOFEL,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      ASUSM assembles the contribution from an element matrix into 
C      into a real unsymmetric system matrix
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
C      Commented     7 Feb 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      SYSK    contains system matrix into which contributions
C              from current element matrix are to be assembled
C      ISYSK   first dimension of SYSK (.GE. total number of
C              unconstrained degrees of freedom)
C      JSYSK   second dimension of SYSK (.GE. 2*HBAND-1)
C      ELK     element matrix
C      IELK    first dimension of ELK (.GE. DOFEL)
C      JELK    second dimension of ELK (.GE. DOFEL)
C      STEER   contains freedom numbers associated with element
C              matrix contributions to system matrix
C      ISTEER  first dimension of STEER (.GE. DOFEL)
C      HBAND   semi-bandwidth of system matrix, including
C              diagonal
C      DOFEL   maximum number of degrees of freedom associated
C              with this element type
C      ITEST   error checking option
C
C ARGUMENTS out
C      SYSK    system matrix
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE ASUSM(SYSK,ISYSK,JSYSK,ELK,IELK,JELK,STEER,ISTEER,
C    *                 HBAND,DOFEL,ITEST)
C***********************************************************************
C
      INTEGER DOFEL,HBAND,IELK,ISTEER,ISYSK,ITEST,JELK,JSYSK
      INTEGER STEER(ISTEER)
      DOUBLE PRECISION ELK(IELK,JELK),SYSK(ISYSK,JSYSK)
C
      INTEGER I,IERROR,J,JTEST,OFFSET,STEERI,STEERJ
      CHARACTER*6 SRNAME
C
      INTEGER  ERRMES
      EXTERNAL ERRMES
C
      DATA SRNAME/'ASUSM'/
C
C     Parameter checking
C
      JTEST = ITEST
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (ISTEER.LT.DOFEL) IERROR = 4
         IF (IELK.LT.DOFEL .OR. JELK.LT.DOFEL) IERROR = 3
         IF (JSYSK.LT.2*HBAND-1) IERROR = 2
         IF (HBAND.LE.0 .OR. DOFEL.LE.0) IERROR = 1
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
                  OFFSET = STEERJ - STEERI + HBAND
C
C     Range checking on STEERI and STEERJ
C
                  IF (ITEST.NE.-1) THEN
                     IERROR = 0
                     IF (ISYSK.LT.STEERI .OR.
     *                   JSYSK.LT.OFFSET) IERROR = 5
                     ITEST = ERRMES(JTEST,IERROR,SRNAME)
                     IF (ITEST.NE.0) RETURN
                  END IF
C
                  SYSK(STEERI,OFFSET) = SYSK(STEERI,OFFSET) + ELK(I,J)
               END IF
 1000       CONTINUE
C
         END IF
 1010 CONTINUE
C
      END
