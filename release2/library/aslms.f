C***********************************************************************
C$SPLIT$ASLMS$*********************************************************
C***********************************************************************
      SUBROUTINE ASLMS(SYSM, ISYSM, ELM, IELM, JELM, STEER,
     *     ISTEER, DOFEL, DOFNOD, SIZE, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      ASSEMBLES THE CONTRIBUTION FROM AN ELEMENT TO THE
C      (DIAGONAL) SYSTEM MATRIX WHICH IS STORED AS A VECTOR.
C      THE 'LUMPED MASS' APPROXIMATION IS ASSUMED, THE DIAGONAL
C      ELEMENTS OF THE ELEMENT 'CONSISTENT MASS' MATRIX BEING
C      USED, SUITABLY BIASSED TO CONSERVE THE RELEVANT QUANTITY
C      (EG MASS) ON THE ELEMENT
C
C HISTORY
C
C      COPYRIGHT (C) 1979 : SERC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 1.1  29 OCT 1979 (CG) 
C      COMMENTED    31 OCT 1980 (KR)
C
C ARGUMENTS IN
C      ISYSM   DIMENSION OF VECTOR SYSM (.GE.TOTAL NUMBER OF
C              UNCONSTRAINED FREEDOMS)
C      ELM     ELEMENT MASS MATRIX OF DIMENSION (IELM,JELM).
C              ON ENTRY, ELM(I,I) SHOULD CONTAIN THE CONSISTENT
C              MASS APPROXIMATIONS FOR THE ELEMENT FOR
C              I=1(1)DOFEL
C      IELM    FIRST DIMENSION OF ELM (.GE.DOFEL)
C      JELM    SECOND DIMENSION OF ELM (.GE.DOFEL)
C      STEER   INTEGER VECTOR OF LENGTH ISTEER CONTAINING
C              FREEDOM NUMBERS ASSOCIATING ELEMENT MATRIX
C              CONTRIBUTIONS TO SYSTEM FREEDOM NUMBERS
C      ISTEER  LENGTH OF VECTOR STEER (.GE.DOFEL)
C      DOFEL   MAXIMUM NUMBER OF DEGREES OF FREEDOM ASSOCIATED
C              WITH THE ELEMENT TYPE
C      DOFNOD  NUMBER OF DEGREES OF FREEDOM PER NODE ON THE
C              ELEMENT
C      SIZE    IN 2D, AREA OF THE ELEMENT
C              IN 3D, VOLUME OF THE ELEMENT
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      SYSM    VECTOR OF LENGTH ISYSM CONTAINING THE DIAGONAL
C              ELEMENTS OF THE (DIAGONAL) SYSTEM MATRIX
C      ELM     ELEMENT 'MASS' MATRIX, OF DIMENSION (IELM,JELM).
C              ON EXIT, ELM(I,J) CONTAINS ZERO IF I.NE.J, AND
C              THE CALCULATED 'LUMPED MASS' VALUES IF I=J
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE ASLMS(SYSM,ISYSM,ELM,IELM,JELM,STEER,ISTEER,
C    *     DOFEL,DOFNOD,SIZE,ITEST)
C***********************************************************************
C
      INTEGER DOFEL, DOFNOD, ERRMES, I, IELM, IERROR, ISTEER,
     *     ISYSM, ITEST, J, JELM, JTEST, STEER
      DOUBLE PRECISION ELM, SIZE, SRNAME, SYSM, X, ZERO
      DIMENSION ELM(IELM,JELM), STEER(ISTEER), SYSM(ISYSM)
      DATA SRNAME /8H ASLMS  /, ZERO /0.0D0/
C
C     PARAMETER CHECKING
C
      JTEST = ITEST
      IF (JTEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (IELM.LT.DOFEL .OR. JELM.LT.DOFEL) IERROR = 2
      IF (ISTEER.LT.DOFEL) IERROR = 3
      IF (DOFEL.EQ.0 .OR. DOFNOD.EQ.0) IERROR = 1
      ITEST = ERRMES(JTEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
C     MAIN LOOPS
C
 1010 X = 0.0D0
      DO 1020 I=1,DOFEL
         X = X + ELM(I,I)
 1020 CONTINUE
C     X = SIZE/X*(((DBLE(FLOAT(DOFNOD)))))
      X = SIZE/X*((((DBLE(FLOAT(DOFNOD))))))
      DO 1050 I=1,DOFEL
         DO 1040 J=1,DOFEL
            IF (I.EQ.J) GO TO 1030
            ELM(I,J) = ZERO
            GO TO 1040
 1030       ELM(I,J) = ELM(I,J)*X
 1040    CONTINUE
 1050 CONTINUE
      DO 1070 I=1,DOFEL
         J = STEER(I)
         IF (J.EQ.0) GO TO 1070
C
C     RANGE CHECKING ON STEERI AND STEERJ
C
         IF (JTEST.EQ.-1) GO TO 1060
         IERROR = 0
         IF (ISYSM.LT.J) IERROR = 4
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
C
 1060    SYSM(J) = SYSM(J) + ELM(I,I)
C
 1070 CONTINUE
      RETURN
      END
