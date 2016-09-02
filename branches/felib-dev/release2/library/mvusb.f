C***********************************************************************
C$SPLIT$MVUSB$*********************************************************
C***********************************************************************
      SUBROUTINE MVUSB(A, IA, JA, V, IV, W, IW, N, HBAND, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      POST-MULTIPLIES A REAL UNSYMMETRIC BANDED MATRIX STORED AS
C      A RECTANGULAR ARRAY BY A VECTOR
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    14 OCT 1980 (KR)
C
C ARGUMENTS IN
C      A       ARRAY OF DIMENSION (IA,JA).  CONTAINS THE
C              ELEMENTS OF THE REAL UNSYMMETRIC BAND MATRIX
C              OF ORDER N AND SEMI-BANDWIDTH HBAND
C      IA      FIRST DIMENSION OF A (.GE.N)
C      JA      SECOND DIMENSION OF A (.GE.2*HBAND-1)
C      V       VECTOR OF DIMENSION IV
C      IV      DIMENSION OF VECTOR V (.GE.N)
C      IW      DIMENSION OF VECTOR W (.GE.N)
C      N       ORDER OF THE REAL UNSYMMETRIC BAND MATRIX
C      HBAND   SEMI-BANDWIDTH OF THE REAL UNSYMMETRIC BAND MATRIX
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      W       VECTOR OF DIMENSION IW.  CONTAINS THE RESULT OF
C              THE OPERATION W=A*V
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE MVUSB(A, IA, JA, V, IV, W, IW, N, HBAND, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, HBAND, I, IERROR, IF, IS, ITEST,
     *     IV, IW, J, N, K
      DOUBLE PRECISION A, SRNAME, V, W, X
      DIMENSION A(IA,JA), V(IV), W(IW)
      DATA SRNAME /8H MVUSB  /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (IW.LT.N) IERROR = 4
                        IF (IV.LT.N) IERROR = 3
                        IF (IA.LT.N .OR. JA.LT.2*HBAND-1)
     *                      IERROR = 2
                        IF (N.LE.0) IERROR = 1
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 DO 1030 I=1,N
      X = 0.0D0
      IS = MAX0(HBAND+1-I,1)
      IF = MIN0(N+HBAND-I,2*HBAND-1)
      DO 1020 J=IS,IF
      K = I + J - HBAND
      X = X + A(I,J)*V(K)
 1020 CONTINUE
      W(I) = X
 1030 CONTINUE
      RETURN
      END
