C***********************************************************************
      SUBROUTINE CHOFWD(A, IA, JA, R, IR, N, HBAND, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      PERFORMS FORWARD SUBSTITUTION ON A MATRIX REDUCED BY
C      CHORDN
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 1.1  29 OCT 1979 (CG)
C      COMMENTED    10 OCT 1980 (KR)
C
C ARGUMENTS IN
C      A       ARRAY OF DIMENSION (IA,JA).  CONTAINS THE
C              ELEMENTS OF THE LOWER HALF OF THE POSITIVE
C              DEFINITE SYMMETRIC BAND MATRIX OF ORDER N AND
C              SEMI-BANDWIDTH HBAND.  A SHOULD PREVIOUSLY HAVE
C              BEEN REDUCED USING CHORDN OR CHOSOL
C      IA      FIRST DIMENSION OF A (.GE.N)
C      JA      SECOND DIMENSION OF A (.GE.HBAND)
C      R       ON ENTRY, CONTAINS THE VECTOR OF RHS'S
C      IR      DIMENSION OF R (.GE.N)
C      N       ORDER OF MATRIX A
C      HBAND   SEMI-BANDWIDTH OF A
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      R       ON EXIT, CONTAINS THE REDUCED VECTOR OF RHS'S
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE CHOFWD(A, IA, JA, R, IR, N, HBAND, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, HBAND, I, IA, IERROR, IJ, IR, ITEST,
     *     J, JA, K, N, W
      DOUBLE PRECISION A, R, SRNAME, X
      DIMENSION A(IA,JA), R(IR)
      DATA SRNAME /8H CHOFWD /
C
C     PARAMETER CHECKING
C
      IF (ITEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (IR.LT.N) IERROR = 3
      IF (IA.LT.N .OR. JA.LT.HBAND) IERROR = 2
      IF (N.LE.0 .OR. HBAND.LE.0) IERROR = 1
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
C     MAIN BODY
C
 1010 W = HBAND - 1
      R(1) = R(1)/A(1,W+1)
      DO 1030 I=2,N
         X = 0.0D0
         K = 1
         IF (I.LE.W+1) K = W - I + 2
         DO 1020 J=K,W
            IJ = I + J - W - 1
            X = X + A(I,J)*R(IJ)
 1020    CONTINUE
         R(I) = (R(I)-X)/A(I,W+1)
 1030 CONTINUE
C
      RETURN
      END
C***********************************************************************
