C***********************************************************************
      SUBROUTINE CSYRDN(KB, IKB, JKB, N, HBAND, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      PERFORMS SYMMETRIC REDUCTION ON A COMPLEX SYMMETRIC BANDED
C      BANDED MATRIX.  ONLY THE LOWER BAND AND DIAGONAL ARE STORED
C      IN A RECTANGULAR ARRAY KB
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 3.0  29 OCT 1984 (CRIE)
C      COMMENTED     1 NOV 1985 (CG)
C
C ARGUMENTS IN
C      KB      ON ENTRY, CONTAINS LOWER BAND AND DIAGONAL OF
C              COMPLEX SYMMETRIC MATRIX
C      IKB     FIRST DIMENSION OF KB (.GE. N)
C      JKB     SECOND DIMENSION OF KB (.GE. HBAND)
C      N       ORDER OF MATRIX KB
C      HBAND   SEMI-BANDWIDTH OF KB  (INCLUDING DIAGONAL)
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      KB      REDUCED MATRIX L, WHERE THE INPUT MATRIX KB  HAS
C              BEEN REDUCED TO TRIANGULAR MATRICES L AND LT
C              WHERE KB =L LT
C
C ROUTINES CALLED
C      ERRMES
C
C     SUBROUTINE CSYRDN(KB, IKB, JKB, N, HBAND, ITEST)
C**********************************************************************
C
      INTEGER A, B, HBAND, I, IK, IKB, ITEST, J, JKB, K, L, W,
     *     LK, N, IERROR, ERRMES
      DOUBLE PRECISION KB, X, Y, AR, AI, XR, XI
      DOUBLE PRECISION SRNAME
      DIMENSION KB(2,IKB,JKB)
      DATA SRNAME /8H CSYRDN /
C
C     PARAMETER CHECKING
C
                IF(ITEST.EQ.-1) GO TO 999
                IERROR=0
                IF(IKB.LT.N.OR.JKB.LT.HBAND) IERROR=2
                IF(N.LE.0.OR.HBAND.LE.0) IERROR=1
                ITEST=ERRMES(ITEST,IERROR,SRNAME)
                IF(ITEST.NE.0) RETURN
C
C     MAIN LOOPS
C
999   W = HBAND - 1
      DO 1050 I=1,N
      XR = 0.0D00
      XI = 0.0D00
      DO 1010 J=1,W
      X=KB(1,I,J)
      Y=KB(2,I,J)
      XR = XR + X*X - Y*Y
      XI = XI + X*Y*2.0D0
1010  CONTINUE
C
                IF(ITEST.EQ.-1) GO TO 998
                IERROR=0
                X=KB(1,I,W+1)-XR
                Y=KB(2,I,W+1)-XI
                IF(DSQRT(X*X+Y*Y) .LE. 0.0D00) IERROR=3
                ITEST=ERRMES(ITEST,IERROR,SRNAME)
                IF(ITEST.NE.0) RETURN
C
998   IF(KB(1,I,W+1)-XR.GE.0.0D00)GOTO 900
      AR = KB(1,I,W+1)-XR
      AI = KB(2,I,W+1)-XI
      Y = DSIGN(1.0D0,AI) * DSQRT( (DABS(AR)+DSQRT(AR*AR + AI*AI))
     *    /2.0D0)
      X=0.0D0
      IF (Y.GT.0.0D0)
     *     X = AI/(2.0D0*Y)
      KB(1,I,W+1) = X
      KB(2,I,W+1) = Y
      GOTO 901
 900  AR = KB(1,I,W+1)-XR
      AI = KB(2,I,W+1)-XI
      X = DSQRT( (AR+DSQRT(AR*AR+AI*AI))/2.0D0 )
      Y=0.0D0
      IF (X.GT.0.0D0)
     *     Y = AI/(2.0D0*X)
      KB(1,I,W+1) = X
      KB(2,I,W+1) = Y
901   CONTINUE
      DO 1040 K=1,W
      XR = 0.0D0
      XI = 0.0D0
      IF (I+K.GT.N) GO TO 1040
      IF (K.EQ.W) GO TO 1030
      L = W - K
1020  IK = I + K
      LK = L + K
      X = KB(1,IK,L)
      Y = KB(2,IK,L)
      AR = KB(1,I,LK)
      AI = KB(2,I,LK)
      XR = XR + X*AR - Y*AI
      XI = XI + X*AI + Y*AR
      L = L - 1
      IF (L.NE.0) GO TO 1020
1030  A = I + K
      B = W - K + 1
      X = KB(1,A,B)-XR
      Y = KB(2,A,B)-XI
      AR = KB(1,I,W+1)
      AI = KB(2,I,W+1)
      IF(DABS(AR).GT.DABS(AI) ) GOTO 902
      KB(1,A,B) = ((AR/AI)*X + Y)/((AR/AI)*AR + AI)
      KB(2,A,B) = ((AR/AI)*Y - X)/((AR/AI)*AR + AI)
      GOTO 903
 902  CONTINUE
      KB(1,A,B) = (X + (AI/AR)*Y)/((AI/AR)*AI + AR)
      KB(2,A,B) = (Y - (AI/AR)*X)/((AI/AR)*AI + AR)
 903  CONTINUE
1040  CONTINUE
1050  CONTINUE
      RETURN
C
      END
C***********************************************************************
