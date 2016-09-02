C***********************************************************************
      SUBROUTINE VECCOP(V, IV, W, IW, N, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      COPIES FIRST N ELEMENTS OF THE VECTOR V INTO THE VECTOR W
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
C      V       VECTOR OF LENGTH IV TO BE COPIED
C      IV      LENGTH OF VECTOR V (.GE.N)
C      IW      LENGTH OF VECTOR W
C      N       NUMBER OF ELEMENTS OF V TO BE COPIED
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      W       VECTOR OF LENGTH IW; W(I)IS SET TO V(I) FOR
C              I=1(1)N
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE VECCOP(V, IV, W, IW, N, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IERROR, ITEST, IV, IW, N
      DOUBLE PRECISION SRNAME, V, W
      DIMENSION V(IV), W(IW)
C
      DATA SRNAME /8H VECCOP /
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
         W(I) = V(I)
 1020 CONTINUE
C
      RETURN
      END
C***********************************************************************
