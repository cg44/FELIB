C***********************************************************************
      SUBROUTINE SELECT(V, IV, STEER, ISTEER, N, W, IW, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      TO CONSTRUCT AN ELEMENT VALUE VECTOR FROM A FULL SYSTEM
C      VECTOR USING THE STEERING VECTOR
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 3.0  29 OCT 1985 (CG)
C
C ARGUMENTS IN
C      V       VECTOR OF SYSTEM VALUES
C      IV      FIRST DIMENSION OF VECTOR V
C      STEER   ELEMENT STEERING VECTOR - CONTAINS FREEDOM NOS.
C      ISTEER  FIRST DIMENSION OF VECTOR STEER (.GE.N)
C      N       NUMBER OF FREEDOMS TO BE ASSEMBLED
C      IW      FIRST DIMENSION OF VECTOR W (.GE.N)
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      W       VECTOR OF ELEMENT VALUES
C
C ROUTINES CALLED
C      ERRMES
C
C      SUBROUTINE SELECT(V,IV,STEER,ISTEER,N,W,IW,ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IERROR, ISTEER, ITEST, IV, IW, J,
     *     N, STEER, JTEST
      DOUBLE PRECISION SRNAME, V, W
      DIMENSION STEER(ISTEER), V(IV), W(IW)
C
      DATA SRNAME /8H SELECT /
C
C     PARAMETER CHECKING
C
      JTEST = ITEST
      IF (JTEST.EQ.-1) GO TO 1010
      IERROR=0
      IF(ISTEER.LT.N .OR. IW.LT.N) IERROR=2
      IF(N.LE.0) IERROR=1
      ITEST = ERRMES(JTEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
C     MAIN LOOPS
C
 1010 DO 1020 I=1,N
         W(I) = 0.D0
         J = IABS(STEER(I))
         IF (J.EQ.0) GO TO 1020
C
         IF (JTEST.EQ.-1) GO TO 1021
         IERROR=0
         IF(J.GT.IV) IERROR=3
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
C
 1021 W(I) = V(J)
 1020 CONTINUE
      RETURN
C
      END
C***********************************************************************