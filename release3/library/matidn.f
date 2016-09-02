C***********************************************************************
      SUBROUTINE MATIDN(A, IA, JA, M, N, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      SETS THE MATRIX A TO THE IDENTITY MATRIX
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
C      M       NUMBER OF ROWS OF A TO BE ASSIGNED VALUES
C      N       NUMBER OF COLUMNS OF A TO BE ASSIGNED VALUES
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      A       ARRAY OF DIMENSION (IA,JA).  A(I,J) IS SET TO 1
C              IF I=J, 0 OTHERWISE
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE MATIDN(A, IA, JA, M, N, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IA, IERROR, ITEST, J, JA, L, M, N
      DOUBLE PRECISION A, SRNAME
      DIMENSION A(IA,JA)
      DATA SRNAME /8H MATIDN /
C
C     PARAMETER CHECKING
C
      IF (ITEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (M.GT.IA .OR. N.GT.JA) IERROR = 2
      IF (N.LE.0 .OR. M.LE.0) IERROR = 1
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF (IERROR.NE.0) RETURN
C
C     MAIN BODY
C
 1010 DO 1030 J=1,N
         DO 1020 I=1,M
            A(I,J) = 0.0D0
 1020    CONTINUE
 1030 CONTINUE
C
C     SET DIAGONAL ENTRY
C
      L = MIN0(M,N)
      DO 1040 I=1,L
         A(I,I) = 1.0D0
 1040 CONTINUE
      RETURN
C
      END
C***********************************************************************
