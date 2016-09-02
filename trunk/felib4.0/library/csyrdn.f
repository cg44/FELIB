C
      SUBROUTINE CSYRDN(KB,IKB,JKB,N,HBAND,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CSYRDN performs symmetric reduction on a complex symmetric banded
C      banded matrix.  Only the lower band and diagonal are stored
C      in a rectangular array KB.
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
C      Release 2.0  29 Oct 1984 (CRIE)
C      Commented     1 Nov 1985 (CG)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      KB      on entry, contains lower band and diagonal of
C              complex symmetric matrix
C      IKB     first dimension of KB (.GE. N)
C      JKB     second dimension of KB (.GE. HBAND)
C      N       order of matrix KB
C      HBAND   semi-bandwidth of KB  (including diagonal)
C      ITEST   error checking option
C
C ARGUMENTS out
C      KB      reduced matrix L, where the input matrix KB  has
C              been reduced to triangular matrices L and LT
C              where KB =L*LT
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE CSYRDN(KB,IKB,JKB,N,HBAND,ITEST)
C**********************************************************************
C
      INTEGER A,B,HBAND,I,IK,IKB,ITEST,J,JKB,K,L,W,LK,N,IERROR,ERRMES
      CHARACTER*6 SRNAME
      DOUBLE PRECISION KB,X,Y,AR,AI,XR,XI
      DIMENSION KB(2,IKB,JKB)
      DATA SRNAME/'CSYRDN'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IKB.LT.N .OR. JKB.LT.HBAND) IERROR = 2
         IF (N.LE.0 .OR. HBAND.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main loops
C
      W = HBAND - 1
      DO 1030 I = 1,N
         XR = 0.0D00
         XI = 0.0D00
         DO 1000 J = 1,W
            X = KB(1,I,J)
            Y = KB(2,I,J)
            XR = XR + X*X - Y*Y
            XI = XI + X*Y*2.0D0
 1000    CONTINUE
C
         IF (ITEST.NE.-1) THEN
            IERROR = 0
            X = KB(1,I,W+1) - XR
            Y = KB(2,I,W+1) - XI
            IF (DSQRT(X*X+Y*Y).LE.0.0D00) IERROR = 3
            ITEST = ERRMES(ITEST,IERROR,SRNAME)
            IF (ITEST.NE.0) RETURN
         END IF
C
         IF (KB(1,I,W+1)-XR.GE.0.0D00) THEN
            AR = KB(1,I,W+1) - XR
            AI = KB(2,I,W+1) - XI
            X = DSQRT((AR+DSQRT(AR*AR+AI*AI))/2.0D0)
            Y = 0.0D0
            IF (X.GT.0.0D0) Y = AI/ (2.0D0*X)
            KB(1,I,W+1) = X
            KB(2,I,W+1) = Y
         ELSE
            AR = KB(1,I,W+1) - XR
            AI = KB(2,I,W+1) - XI
            Y = DSIGN(1.0D0,AI)*DSQRT((DABS(AR)+DSQRT(AR*AR+AI*AI))/
     *          2.0D0)
            X = 0.0D0
            IF (Y.GT.0.0D0) X = AI/ (2.0D0*Y)
            KB(1,I,W+1) = X
            KB(2,I,W+1) = Y
         END IF
         DO 1020 K = 1,W
            XR = 0.0D0
            XI = 0.0D0
            IF (I+K.LE.N) THEN
               IF (K.NE.W) THEN
                  L = W - K
 1010             CONTINUE
                  IK = I + K
                  LK = L + K
                  X = KB(1,IK,L)
                  Y = KB(2,IK,L)
                  AR = KB(1,I,LK)
                  AI = KB(2,I,LK)
                  XR = XR + X*AR - Y*AI
                  XI = XI + X*AI + Y*AR
                  L = L - 1
                  IF (L.NE.0) GO TO 1010
               END IF
               A = I + K
               B = W - K + 1
               X = KB(1,A,B) - XR
               Y = KB(2,A,B) - XI
               AR = KB(1,I,W+1)
               AI = KB(2,I,W+1)
               IF (DABS(AR).GT.DABS(AI)) THEN
                  KB(1,A,B) = (X+ (AI/AR)*Y)/ ((AI/AR)*AI+AR)
                  KB(2,A,B) = (Y- (AI/AR)*X)/ ((AI/AR)*AI+AR)
               ELSE
                  KB(1,A,B) = ((AR/AI)*X+Y)/ ((AR/AI)*AR+AI)
                  KB(2,A,B) = ((AR/AI)*Y-X)/ ((AR/AI)*AR+AI)
               END IF
            END IF
 1020    CONTINUE
 1030 CONTINUE
C
      END
