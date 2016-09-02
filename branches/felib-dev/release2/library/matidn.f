C***********************************************************************
C$SPLIT$MATIDN$*********************************************************
C***********************************************************************
      SUBROUTINE MATIDN(A, IA, JA, M, N, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      SETS THE MATRIX A TO THE IDENTITY MATRIX
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
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
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (M.GT.IA .OR. N.GT.JA) IERROR = 2
                        IF (N.LE.0 .OR. M.LE.0) IERROR = 1
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (IERROR.NE.0) RETURN
 1010                   DO 1030 I=1,M
                        DO 1020 J=1,N
                        A(I,J) = 0.0D0
 1020                   CONTINUE
 1030                   CONTINUE
                        L = MIN0(M,N)
                        DO 1040 I=1,L
                        A(I,I) = 1.0D0
 1040                   CONTINUE
                        RETURN
                        END
