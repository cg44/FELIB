C
      SUBROUTINE ASLMS(SYSM,ISYSM,ELM,IELM,JELM,STEER,ISTEER,DOFEL,
     *                 DOFNOD,SIZE,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      ASLMS assembles the contribution from an element to the
C      (diagonal) system matrix which is stored as a vector.
C      the 'lumped mass' approximation is assumed, the diagonal
C      elements of the element 'consistent mass' matrix being
C      used, suitably biassed to conserve the relevant quantity
C      (eg mass) on the element.
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
C      Release 1.1   29 Oct 1979 (CG)
C      Commented     31 Oct 1980 (KR)
C      Release 3.0    1 Jul 1996 (CG)
C      Release 4.0    2 Oct 2003 (CG)
C
C ARGUMENTS in
C      ISYSM   dimension of vector SYSM (.GE. total number of
C              unconstrained freedoms)
C      ELM     element mass matrix of dimension (IELEM, JELEM)
C              on entry, ELM(I,J) should contain the consistent
C              mass approximations for the element for
C              I=1(1)DOFEL
C      IELM    first dimension of ELM (.GE. DOFEL)
C      JELM    second dimension of ELM (.GE. DOFEL)
C      STEER   integer vector of length ISTEER containing
C              freedom numbers associating element matrix
C              contributions to system freedom numbers
C      ISTEER  length of vector STEER (.GE. DOFEL)
C      DOFEL   maximum number of degrees of freedom associated
C              with the element type
C      DOFNOD  number of degrees of freedom per node on the
C              element
C      SIZE    in two dimensions, area of the element
C              in three dimensions, volume of the element
C      ITEST   error checking option
C
C ARGUMENTS out
C      SYSM    vector of length ISYSM containing the diagonal
C              elements of the (diagonal) system matrix
C      ELM     element 'mass' matrix, of dimension (IELM,JELM).
C              on exit, ELM(I,J) contains ZERO if I .NE. J, and
C              the calculated 'lumped mass' values if I=J
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE ASLMS(SYSM,ISYSM,ELM,IELM,JELM,STEER,ISTEER,DOFEL,
C    *                 DOFNOD,SIZE,ITEST)
C***********************************************************************
C
      INTEGER DOFEL,DOFNOD,IELM,ISTEER,ISYSM,ITEST,JELM
      INTEGER STEER(ISTEER)
      DOUBLE PRECISION ELM(IELM,JELM),SYSM(ISYSM)
C
      INTEGER I,IERROR,J,JTEST
      DOUBLE PRECISION SIZE,X,ZERO
      CHARACTER*6 SRNAME
C
      INTEGER ERRMES
      EXTERNAL ERRMES
C
      DATA ZERO/0.0D00/
      DATA SRNAME/'ASLMS'/
C
C     Parameter checking
C
      JTEST = ITEST
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (IELM.LT.DOFEL .OR. JELM.LT.DOFEL) IERROR = 2
         IF (ISTEER.LT.DOFEL) IERROR = 3
         IF (DOFEL.EQ.0 .OR. DOFNOD.EQ.0) IERROR = 1
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main loops
C
      X = 0.0D0
      DO 1000 I = 1,DOFEL
         X = X + ELM(I,I)
 1000 CONTINUE
      X = SIZE/X* ((((DBLE(FLOAT(DOFNOD))))))
      DO 1020 I = 1,DOFEL
         DO 1010 J = 1,DOFEL
            IF (I.EQ.J) THEN
               ELM(I,J) = ELM(I,J)*X
            ELSE
               ELM(I,J) = ZERO
            END IF
 1010    CONTINUE
 1020 CONTINUE
      DO 1030 I = 1,DOFEL
         J = STEER(I)
         IF (J.NE.0) THEN
C
C     Range checking on I and J
C
            IF (JTEST.NE.-1) THEN
               IERROR = 0
               IF (ISYSM.LT.J) IERROR = 4
               ITEST = ERRMES(JTEST,IERROR,SRNAME)
               IF (ITEST.NE.0) RETURN
            END IF
C
            SYSM(J) = SYSM(J) + ELM(I,I)
         END IF
 1030 CONTINUE
C
      END
