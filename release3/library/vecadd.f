C***********************************************************************
      SUBROUTINE VECADD(V, IV, W, IW, N, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      ADDS THE VECTOR V TO THE VECTOR W, STORING THE RESULT
C      IN V
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 1.1  29 OCT 1979 (IMS)
C      COMMENTED    23 OCT 1980 (KR)
C
C ARGUMENTS IN
C      V       VECTOR OF LENGTH IV.  ON ENTRY, CONTAINS THE
C              VALUES OF ONE OF THE VECTORS TO BE ADDED
C              IN ELEMENTS 1 TO N
C      IV      LENGTH OF VECTOR V (.GE.N)
C      W       VECTOR OF LENGTH IW.  CONTAINS THE VALUES OF THE
C              SECOND VECTOR TO BE ADDED IN ELEMENTS 1 TO N
C      IW      LENGTH OF VECTOR W (.GE.N)
C      N       NUMBER OF ELEMENTS OF V AND W TO BE ADDED
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      V       VECTOR OF LENGTH IV.  ON EXIT, V(I) CONTAINS
C              THE SUM OF V(I) AND W(I) FOR I=1(1)N
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE VECADD(V, IV, W, IW, N, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IERROR, ITEST, IV, IW, N
      DOUBLE PRECISION SRNAME, V, W
      DIMENSION V(IV), W(IW)
C
      DATA SRNAME /8H VECADD /
C
C     PARAMETER CHECKING
C
      IF (ITEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (N.GT.IW) IERROR = 3
      IF (N.GT.IV) IERROR = 2
      IF (N.LE.0) IERROR = 1
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
C     MAIN BODY
C
 1010 DO 1020 I=1,N
         V(I) = V(I) + W(I)
 1020 CONTINUE
C
      RETURN
      END
C***********************************************************************
