C***********************************************************************
      SUBROUTINE IASRHS(RHS, IRHS, VALUE, IVALUE, STEER, ISTEER,
     *     DOFEL, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      THE ROUTINE ADDS INTO THE RIGHT-HAND SIDE OF A SYSTEM
C      THE VALUES CONTIANED IN AN ELEMENT VECTOR, THUS
C      FORMING THE RIGHT-HAND SIDE VALUE IS REAL, RHS IS COMPLEX
C      (ORDERED PAIRS), WITH VALUE ADDED INTO IMAGINARY PART OF RHS.
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 3.0   1 JUL 1984 (CRIE)
C      COMMENTED     1 JUL 1985 (CG)
C
C ARGUMENTS IN
C      RHS     THE RIGHT-HAND SIDE OF THE SYSTEM
C              COMPLEX ( ORDERED PAIRS )
C      IRHS    DIMENSION OF ARRAY RSH
C      VALUE   THE ELEMENT VECTOR OF THE CURRENT ELEMENT TO
C              BE ADDED INTO THE RIGHT-HAND SIDE, IMAGINARY COMPONENT
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
C     SUBROUTINE IASRHS(RHS,IRHS,VALUE,IVALUE,STEER,ISTEER,
C    *     DOFEL,ITEST)
C***********************************************************************
C
      INTEGER DOFEL, ERRMES, IERROR, IRHS, ISTEER, ITEST, IVALUE,
     *     K, L, STEER, JTEST
      DOUBLE PRECISION RHS, SRNAME, VALUE
      DIMENSION RHS(2,IRHS), STEER(ISTEER), VALUE(IVALUE)
      DATA SRNAME /8H IASRHS /
C
C     PARAMETER CHECKING
C
      JTEST = ITEST
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (ISTEER.LT.DOFEL) IERROR = 3
                        IF (IVALUE.LT.DOFEL) IERROR = 2
                        IF (DOFEL.LE.0) IERROR = 1
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
C
C     MAIN LOOPS
C
 1010 DO 1030 K=1,DOFEL
      IF (STEER(K).EQ.0) GO TO 1030
      L = STEER(K)
                        IF (JTEST.EQ.-1) GO TO 1020
                        IERROR = 0
                        IF (L.GT.IRHS) IERROR = 4
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1020 RHS(2,L) = RHS(2,L) + VALUE(K)
 1030 CONTINUE
      RETURN
C
      END
C***********************************************************************
