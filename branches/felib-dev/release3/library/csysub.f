C***********************************************************************
      SUBROUTINE CSYSUB(KB, IKB, JKB, LOADS, ILOADS, N, HBAND,
     *     ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      PERFORMS FORWARD AND BACKWARD SUBSTITUTION ON A COMPLEX MATRIX
C      REDUCED BY CSYRDN
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 3.0  29 OCT 1984 (CRIE)
C      COMMENTED    10 OCT 1985 (CG)
C
C ARGUMENTS IN
C      KB      ARRAY OF DIMENSION (IKB,JKB).  CONTAINS THE
C              ELEMENTS OF THE LOWER HALF OF THE COMPLEX SYMMETRIC
C              BAND MATRIX OF ORDER N AND SEMI-BANDWIDTH HBAND. KB
C              SHOULD PREVIOUSLY HAVE BEEN REDUCED USING CSYSOL OR CSYRDN.
C      IKB     FIRST DIMENSION OF KB (.GE.N)
C      JKB     SECOND DIMENSION OF KB (.GE.HBAND)
C      LOADS   ON ENTRY, CONTAINS THE VECTOR OF RHS'S
C      ILOADS  DIMENSION OF LOADS (.GE.N)
C      N       ORDER OF MATRIX KB
C      HBAND   SEMI-BANDWIDTH OF KB
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      LOADS   ON EXIT, CONTAINS THE SOLUTION VECTOR
C
C ROUTINES CALLED
C      ERRMES
C
C     SUBROUTINE CSYSUB(KB, IKB, JKB, LOADS, ILOADS, N, HBAND,
C    *     ITEST)
C**********************************************************************
C
      INTEGER HBAND, I, IKB, ITEST, J, JKB, K, L, W, M, IJ,
     *     N, IERROR, ERRMES, ILOADS
      DOUBLE PRECISION KB, X, LOADS, Y, AR, AI, XR, XI
      DOUBLE PRECISION SRNAME
      DIMENSION KB(2,IKB,JKB), LOADS(2,ILOADS)
      DATA SRNAME /8H CSYSUB /
C
C     PARAMETER CHECKING
C
                IF(ITEST.EQ.-1) GO TO 999
                IERROR=0
                IF(ILOADS.LT.N) IERROR=3
                IF(IKB.LT.N.OR.JKB.LT.HBAND) IERROR=2
                IF(N.LE.0.OR.HBAND.LE.0) IERROR=1
                ITEST=ERRMES(ITEST,IERROR,SRNAME)
                IF(ITEST.NE.0) RETURN
C
C     MAIN BODY
C
999   W = HBAND - 1
      X = LOADS(1,1)
      Y = LOADS(2,1)
      AR = KB(1,1,W+1)
      AI = KB(2,1,W+1)
      IF (DABS(AR).GT.DABS(AI) ) GOTO 904
      LOADS(1,1) = ((AR/AI)*X + Y)/((AR/AI)*AR + AI)
      LOADS(2,1) = ((AR/AI)*Y - X)/((AR/AI)*AR + AI)
      GOTO 905
 904  CONTINUE
      LOADS(1,1) = (X + (AI/AR)*Y)/((AI/AR)*AI + AR)
      LOADS(2,1) = (Y - (AI/AR)*X)/((AI/AR)*AI + AR)
 905  CONTINUE
      DO 2020 I=2,N
      XR = 0.0D0
      XI = 0.0D0
      K = 1
      IF (I.LE.W+1) K = W - I + 2
      DO 2010 J=K,W
      IJ = I + J - W - 1
      X = KB(1,I,J)
      Y = KB(2,I,J)
      AR = LOADS(1,IJ)
      AI = LOADS(2,IJ)
      XR = XR + X*AR - Y*AI
      XI = XI + X*AI + Y*AR
2010  CONTINUE
      X = LOADS(1,I)-XR
      Y = LOADS(2,I)-XI
      AR = KB(1,I,W+1)
      AI = KB(2,I,W+1)
      IF ( DABS(AR).GT.DABS(AI) ) GOTO 906
      LOADS(1,I) = ((AR/AI)*X + Y)/((AR/AI)*AR + AI)
      LOADS(2,I) = ((AR/AI)*Y - X)/((AR/AI)*AR + AI)
      GOTO 907
 906  CONTINUE
      LOADS(1,I) = (X + (AI/AR)*Y)/((AI/AR)*AI + AR)
      LOADS(2,I) = (Y - (AI/AR)*X)/((AI/AR)*AI + AR)
 907  CONTINUE
2020  CONTINUE
      X = LOADS(1,N)
      Y = LOADS(2,N)
      AR = KB(1,N,W+1)
      AI = KB(2,N,W+1)
      IF ( DABS(AR).GT.DABS(AI) ) GOTO 908
      LOADS(1,N) = ((AR/AI)*X + Y)/((AR/AI)*AR + AI)
      LOADS(2,N) = ((AR/AI)*Y - X)/((AR/AI)*AR + AI)
      GOTO 909
 908  CONTINUE
      LOADS(1,N) = (X + (AI/AR)*Y)/((AI/AR)*AI + AR)
      LOADS(2,N) = (Y - (AI/AR)*X)/((AI/AR)*AI + AR)
 909  CONTINUE
      I = N - 1
3010  XR = 0.0D0
      XI = 0.0D0
      L = I + W
      IF (I.GT.N-W) L = N
      M = I + 1
      DO 3020 J=M,L
      IJ = W + I - J + 1
      X = KB(1,J,IJ)
      Y = KB(2,J,IJ)
      AR = LOADS(1,J)
      AI = LOADS(2,J)
      XR = XR + X*AR - Y*AI
      XI = XI + X*AI + Y*AR
3020  CONTINUE
      X = LOADS(1,I)-XR
      Y = LOADS(2,I)-XI
      AR = KB(1,I,W+1)
      AI = KB(2,I,W+1)
      IF ( DABS(AR).GT.DABS(AI) ) GOTO 910
      LOADS(1,I) = ((AR/AI)*X + Y)/((AR/AI)*AR + AI)
      LOADS(2,I) = ((AR/AI)*Y - X)/((AR/AI)*AR + AI)
      GOTO 911
 910  CONTINUE
      LOADS(1,I) = (X + (AI/AR)*Y)/((AI/AR)*AI + AR)
      LOADS(2,I) = (Y - (AI/AR)*X)/((AI/AR)*AI + AR)
 911  CONTINUE
      I = I - 1
      IF (I.NE.0) GO TO 3010
      RETURN
C
      END
C***********************************************************************
