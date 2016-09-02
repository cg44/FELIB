C***********************************************************************
C$SPLIT$ASRHS$*********************************************************
C***********************************************************************
      SUBROUTINE ASRHS(RHS, IRHS, VALUE, IVALUE, STEER, ISTEER,
     *     DOFEL, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      THE ROUTINE ADDS INTO THE RIGHT-HAND SIDE OF A SYSTEM
C      THE VALUES CONTIANED IN AN ELEMENT VECTOR, THUS
C      FORMING THE RIGHT-HAND SIDE
C
C HISTORY
C
C      COPYRIGHT (C) 1981 : SERC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 2.0   1 JUL 1981 (CG) 
C      COMMENTED     1 JUL 1981 (CG)
C
C ARGUMENTS IN
C      RHS     THE RIGHT-HAND SIDE OF THE SYSTEM
C      IRHS    DIMENSION OF ARRAY RSH
C      VALUE   THE ELEMENT VECTOR OF THE CURRENT ELEMENT TO
C              BE ADDED INTO THE RIGHT-HAND SIDE
C      IVALE   DIMENSION OF ARRAY VALUE
C      STEER   THE STEERING VECTOR CONTAINING THE FREEDOM
C              NUMBERS OF THE FREEDOMS ASSOCIATED WITH THE
C              CURRENT ELEMENT IN THE LOCAL ORDER
C      ISTEER  DIMENSION OF ARRAY STEER
C      DOFEL   THE MAXIMUM NUMBER OF DEGREES OF FREEDOM ON
C              AN ELEMENT OF THE CURRENT TYPE
C      ITEST   ERROR CHECKING OPTION
C
C ROUTINES CALLED
C      ERRMES
C
C
C
C     SUBROUTINE ASRHS(RHS,IRHS,VALUE,IVALUE,STEER,ISTEER,
C    * DOFEL,ITEST)
C***********************************************************************
C
      INTEGER DOFEL, ERRMES, IERROR, IRHS, ISTEER, ITEST, IVALUE,
     *     JTEST, K, STEER, STEERI
      DOUBLE PRECISION RHS, SRNAME, VALUE
      DIMENSION RHS(IRHS), STEER(ISTEER), VALUE(IVALUE)
      DATA SRNAME /8H ASRHS  /
C
C     PARAMETER CHECKING
C
      JTEST = ITEST
      IF (JTEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (ISTEER.LT.DOFEL) IERROR = 3
      IF (IVALUE.LT.DOFEL) IERROR = 2
      IF (DOFEL.LE.0) IERROR = 1
      ITEST = ERRMES(JTEST,IERROR,SRNAME)
C
C     MAIN LOOPS
C
 1010 IF (ITEST.NE.0) RETURN
      DO 1030 K=1,DOFEL
         STEERI = STEER(K)
         IF (STEERI.LE.0) GO TO 1030
C
C     RANGE CHECKING ON L
C
         IF (JTEST.EQ.-1) GO TO 1020
         IERROR = 0
         IF (STEERI.GT.IRHS) IERROR = 4
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
C
 1020    RHS(STEERI) = RHS(STEERI) + VALUE(K)
 1030 CONTINUE
C
      RETURN
      END
