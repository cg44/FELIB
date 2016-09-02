C
      SUBROUTINE ASRHS(RHS,IRHS,VALUE,IVALUE,STEER,ISTEER,DOFEL,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      ASRHS the routine adds into the right-hand side of a system
C      the values contianed in an element vector, thus
C      forming the right-hand side.
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
C      Release 2.0    1 Jul 1981 (CG)
C      Commented      1 Jul 1981 (CG)
C      Release 3.0    1 Jul 1996 (CG)
C      Release 4.0    2 Oct 2003 (CG)
C
C ARGUMENTS in
C      RHS     the right-hand side of the system
C      IRHS    dimension of array RHS
C      VALUE   the element vector of the current element to
C              be added into the right-hand side
C      ivale   dimension of array VALUE
C      STEER   the steering vector containing the freedom
C              numbers of the freedoms associated with the
C              current element in the local order
C      ISTEER  dimension of array STEER
C      DOFEL   the maximum number of degrees of freedom on
C              an element of the current type
C      ITEST   error checking option
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE ASRHS(RHS,IRHS,VALUE,IVALUE,STEER,ISTEER,DOFEL,ITEST)
C***********************************************************************
C
      INTEGER DOFEL,IRHS,ISTEER,ITEST,IVALUE,STEERI
      INTEGER STEER(ISTEER)
      DOUBLE PRECISION RHS(IRHS),VALUE(IVALUE)
C
      INTEGER IERROR,JTEST,K
      CHARACTER*6 SRNAME
C
      INTEGER ERRMES
      EXTERNAL ERRMES
C
      DATA SRNAME/'ASRHS'/
C
C     Parameter checking
C
      JTEST = ITEST
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (ISTEER.LT.DOFEL) IERROR = 3
         IF (IVALUE.LT.DOFEL) IERROR = 2
         IF (DOFEL.LE.0) IERROR = 1
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
      END IF
C
C     Main loops
C
      IF (ITEST.EQ.0) THEN
         DO 1000 K = 1,DOFEL
            STEERI = STEER(K)
            IF (STEERI.GT.0) THEN
C
C     Range checking on STEERI
C
               IF (JTEST.NE.-1) THEN
                  IERROR = 0
                  IF (STEERI.GT.IRHS) IERROR = 4
                  ITEST = ERRMES(JTEST,IERROR,SRNAME)
                  IF (ITEST.NE.0) RETURN
               END IF
C
               RHS(STEERI) = RHS(STEERI) + VALUE(K)
            END IF
 1000    CONTINUE
C
      END IF
      END
