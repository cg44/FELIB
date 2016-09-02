C***********************************************************************
      SUBROUTINE CMTNUL(A, IA, JA, M, N, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      SETS COMPLEX MATRIX A TO THE NULL MATRIX
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 3.0  20 OCT 1984 (CRIE)
C      COMMENTS      1 NOV 1985 (CG)
C
C ARGUMENTS IN
C      IA      FIRST DIMENSION OF ARRAY A (.GE.M)
C      JA      SECOND DIMENSION OF ARRAY A (.GE.N)
C      M       NUMBER OF ROWS OF A TO BE SET TO ZERO
C      N       NUMBER OF COLUMNS OF A TO BE SET TO ZERO
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      A       ARRAY OF DIMENSION (2,IA,JA).  A(1,I,J)=0
C              AND A(2,I,J)=0 FOR
C              I=1(1)M AND J=1(1)N
C
C ROUTINES CALLED
C      ERRMES
C
C     SUBROUTINE CMTNUL(A, IA, JA, M, N, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IA, IERROR, ITEST, J, JA, M, N
      DOUBLE PRECISION A, SRNAME
      DIMENSION A(2,IA,JA)
      DATA SRNAME /8H CMTNUL /
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
C     MAIN LOOPS
C
 1010 DO 1030 I=1,M
         DO 1020 J=1,N
            A(1,I,J) = 0.0D0
            A(2,I,J) = 0.0D0
 1020    CONTINUE
 1030 CONTINUE
C
      RETURN
      END
C***********************************************************************
