C***********************************************************************
      SUBROUTINE MATNUL(A, IA, JA, M, N, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      SETS MATRIX A TO THE NULL MATRIX
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 1.1  29 OCT 1979 (CG)
C      COMMENTED    12 OCT 1980 (KR)
C
C ARGUMENTS IN
C      IA      FIRST DIMENSION OF ARRAY A (.GE.M)
C      JA      SECOND DIMENSION OF ARRAY A (.GE.N)
C      M       NUMBER OF ROWS OF A TO BE SET TO ZERO
C      N       NUMBER OF COLUMNS OF A TO BE SET TO ZERO
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      A       ARRAY OF DIMENSION (IA,JA).  A(I,J)=0 FOR
C              I=1(1)M AND J=1(1)N
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE MATNUL(A, IA, JA, M, N, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IA, IERROR, ITEST, J, JA, M, N
      DOUBLE PRECISION A, SRNAME
      DIMENSION A(IA,JA)
      DATA SRNAME /8H MATNUL /
C
C     PARAMETER CHECKING
C
      IF (ITEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (M.GT.IA .OR. N.GT.JA) IERROR = 2
      IF (M.LE.0 .OR. N.LE.0) IERROR = 1
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
C     MAIN BODY
C
 1010 DO 1030 J=1,N
         DO 1020 I=1,M
            A(I,J) = 0.0D0
 1020    CONTINUE
 1030 CONTINUE
C
      RETURN
      END
C***********************************************************************
