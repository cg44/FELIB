C***********************************************************************
      SUBROUTINE VMUSB(V, IV, A, IA, JA, W, IW, N, HBAND, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      PRE-MULTIPLIES A REAL UNSYMMETRIC BANDED MATRIX STORED AS
C      A RECTANGULAR ARRAY BY A VECTOR
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 1.1  29 OCT 1979 (CG)
C      COMMENTED    14 OCT 1980 (KR)
C
C ARGUMENTS IN
C      V       VECTOR OF DIMENSION IV
C      IV      DIMENSION OF VECTOR V (.GE.N)
C      A       ARRAY OF DIMENSION (IA,JA).  CONTAINS THE
C              ELEMENTS OF THE REAL UNSYMMETRIC BAND MATRIX
C              OF ORDER N AND SEMI-BANDWIDTH HBAND
C      IA      FIRST DIMENSION OF A (.GE.N)
C      JA      SECOND DIMENSION OF A (.GE.2*HBAND-1)
C      IW      DIMENSION OF VECTOR W (.GE.N)
C      N       ORDER OF THE REAL UNSYMMETRIC BAND MATRIX
C      HBAND   SEMI-BANDWIDTH OF THE REAL UNSYMMETRIC BAND MATRIX
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      W       VECTOR OF DIMENSION IW.  CONTAINS THE RESULT OF
C              THE OPERATION W=V*A
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE VMSYB(V, IV, A, IA, JA, W, IW, N, HBAND, ITEST)
C***********************************************************************
C
      INTEGER BAND, ERRMES, HBAND, I, IA, IERROR, IJ, ITEST, IV,
     *     IW, J, JA, JI, N
      DOUBLE PRECISION A, SRNAME, V, W, X
      DIMENSION A(IA,JA), V(IV), W(IW)
C
      DATA SRNAME /8H VMUSB  /
C
C     PARAMETER CHECKING
C
      IF (ITEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (IW.LT.N) IERROR = 4
      IF (IV.LT.N) IERROR = 3
      IF (IA.LT.N .OR. JA.LT.2*HBAND-1) IERROR = 2
      IF (N.LE.0) IERROR = 1
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
 1010 BAND = 2*HBAND - 1
C
C     MAIN BODY
C
      DO 1030 I=1,N
         X = 0.0D0
         DO 1020 J=1,BAND
            IJ = J - HBAND + I
            JI = BAND + 1 - J
            IF ((IJ.LT.1) .OR. (IJ.GT.N)) GO TO 1020
            X = X + A(IJ,JI)*V(IJ)
 1020    CONTINUE
         W(I) = X
 1030 CONTINUE
C
      RETURN
      END
C***********************************************************************
