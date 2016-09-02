C
      SUBROUTINE CSYBAK(KB,IKB,JKB,LOADS,ILOADS,N,HBAND,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CSYBAK performs bakward substitution on a matrix processed by
C      CHORDN and CHOFWD
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
C      Release 1.1  29 Oct 1979 (CG)
C      Commented    10 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      KB      array of dimension (IKB, JKB).  contains the
C              elements of the lower half of the positive
C              definite band matrix of order N and with semi-
C              bandwidth HBAND, reduced by CHORDN or CHOSOL
C      IKB     first dimension of KB (.GE. N)
C      JKB     second dimension of KB (.GE. HBAND)
C      LOADS   on entry, contains elements of N rhs's after
C              processing by CHOFWD
C      ILOADS  dimension of vector LOADS (.GE. N)
C      N       order of matrix KB
C      HBAND   semi-bandwidth of matrix a
C      ITEST   error checking option
C
C ARGUMENTS out
C      LOADS   on exit, contains solution vector
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE CSYBAK(KB,IKB,JKB,LOADS,ILOADS,N,HBAND,ITEST)
C***********************************************************************
C
      INTEGER HBAND,I,IKB,ITEST,J,JKB,K,L,W,M,IJ,N,IERROR,ERRMES,ILOADS
      CHARACTER*6 SRNAME
      DOUBLE PRECISION KB,X,LOADS,Y,AR,AI,XR,XI
      DIMENSION KB(2,IKB,JKB),LOADS(2,ILOADS)
      DATA SRNAME/'CSYBAK'/
C
C     Paramemter cecking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (ILOADS.LT.N) IERROR = 3
         IF (IKB.LT.N .OR. JKB.LT.HBAND) IERROR = 2
         IF (N.LE.0 .OR. HBAND.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main statements
C
      W = HBAND - 1
      X = LOADS(1,1)
      Y = LOADS(2,1)
      AR = KB(1,1,W+1)
      AI = KB(2,1,W+1)
      IF (DABS(AR).GT.DABS(AI)) THEN
         LOADS(1,1) = (X+ (AI/AR)*Y)/ ((AI/AR)*AI+AR)
         LOADS(2,1) = (Y- (AI/AR)*X)/ ((AI/AR)*AI+AR)
      ELSE
         LOADS(1,1) = ((AR/AI)*X+Y)/ ((AR/AI)*AR+AI)
         LOADS(2,1) = ((AR/AI)*Y-X)/ ((AR/AI)*AR+AI)
      END IF
      DO 1010 I = 2,N
         XR = 0.0D0
         XI = 0.0D0
         K = 1
         IF (I.LE.W+1) K = W - I + 2
         DO 1000 J = K,W
            IJ = I + J - W - 1
            X = KB(1,I,J)
            Y = KB(2,I,J)
            AR = LOADS(1,IJ)
            AI = LOADS(2,IJ)
            XR = XR + X*AR - Y*AI
            XI = XI + X*AI*Y*AR
 1000    CONTINUE
         X = LOADS(1,I) - XR
         Y = LOADS(2,I) - XI
         AR = KB(1,I,W+1)
         AI = KB(2,I,W+1)
         IF (DABS(AR).GT.DABS(AI)) THEN
            LOADS(1,I) = (X+ (AI/AR)*Y)/ ((AI/AR)*AI+AR)
            LOADS(2,I) = (Y- (AI/AR)*X)/ ((AI/AR)*AI+AR)
         ELSE
            LOADS(1,I) = ((AR/AI)*X+Y)/ ((AR/AI)*AR+AI)
            LOADS(2,I) = ((AR/AI)*Y-X)/ ((AR/AI)*AR+AI)
         END IF
 1010 CONTINUE
      X = LOADS(1,N)
      Y = LOADS(2,N)
      AR = KB(1,N,W+1)
      AI = KB(2,N,W+1)
      IF (DABS(AR).GT.DABS(AI)) THEN
         LOADS(1,N) = (X+ (AI/AR)*Y)/ ((AI/AR)*AI+AR)
         LOADS(2,N) = (Y- (AI/AR)*X)/ ((AI/AR)*AI+AR)
      ELSE
         LOADS(1,N) = ((AR/AI)*X+Y)/ ((AR/AI)*AR+AI)
         LOADS(2,N) = ((AR/AI)*Y-X)/ ((AR/AI)*AR+AI)
      END IF
      I = N - 1
 1020 CONTINUE
      XR = 0.0D0
      XI = 0.0D0
      L = I + W
      IF (I.GT.N-W) L = N
      M = I + 1
      DO 1030 J = M,L
         IJ = W + I - J + 1
         X = KB(1,J,IJ)
         Y = KB(2,J,IJ)
         AR = LOADS(1,J)
         AI = LOADS(2,J)
         XR = XR + X*AR - Y*AI
         XI = XI + X*AI + Y*AR
 1030 CONTINUE
      X = LOADS(1,I) - XR
      Y = LOADS(2,I) - XI
      AR = KB(1,I,W+1)
      AI = KB(2,I,W+1)
      IF (DABS(AR).GT.DABS(AI)) THEN
         LOADS(1,I) = (X+ (AI/AR)*Y)/ ((AI/AR)*AI+AR)
         LOADS(2,I) = (Y- (AI/AR)*X)/ ((AI/AR)*AI+AR)
      ELSE
         LOADS(1,I) = ((AR/AI)*X+Y)/ ((AR/AI)*AR+AI)
         LOADS(2,I) = ((AR/AI)*Y-X)/ ((AR/AI)*AR+AI)
      END IF
      I = I - 1
      IF (I.NE.0) GO TO 1020
      END
