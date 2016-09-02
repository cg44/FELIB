C***********************************************************************
C$SPLIT$MATMUL$*********************************************************
C***********************************************************************
      SUBROUTINE MATMUL(A, IA, JA, B, IB, JB, C, IC, JC, L, M, N,
     *     ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      PRE-MULTIPLIES MATRIX B BY A, STORING THE RESULT IN C
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (IMS) --- SERC COPYRIGHT
C      COMMENTED    12 OCT 1980 (KR)
C
C ARGUMENTS IN
C      A       ARRAY OF DIMENSION (IA,JA)
C      IA      FIRST DIMENSION OF A (.GE.L)
C      JA      SECOND DIMENSION OF A (.GE.M)
C      B       ARRAY OF DIMENSION (IB,JB)
C      IB      FIRST DIMENSION OF B (.GE.M)
C      JB      SECOND DIMENSION OF B (.GE.N)
C      IC      FIRST DIMENSION OF ARRAY C (.GE.L)
C      JC      SECOND DIMENSION OF ARRAY C (.GE.N)
C      L       NUMBER OF ROWS OF A TO BE USED IN MULTIPLICATION
C      M       NUMBER OF COLUMNS OF A AND NUMBER OF ROWS OF B
C              TO BE USED IN MULTIPLICATION
C      N       NUMBER OF COLUMNS OF B TO BE USED IN
C              MULTIPLICATION
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      C       CONTAINS RESULT OF MATRIX MULTIPLICATION (C=A*B)
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE MATMUL(A, IA, JA, B, IB, JB, C, IC, JC, L, M, N,
C    *     ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IA, IB, IC, IERROR, ITEST, J, JA,
     *     JB, JC, K, L, M, N
      DOUBLE PRECISION A, B, C, SRNAME, X
      DIMENSION A(IA,JA), B(IB,JB), C(IC,JC)
      DATA SRNAME /8H MATMUL /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (L.GT.IC .OR. N.GT.JC) IERROR = 4
                        IF (M.GT.IB .OR. N.GT.JB) IERROR = 3
                        IF (L.GT.IA .OR. M.GT.JA) IERROR = 2
                        IF (L.LE.0 .OR. M.LE.0 .OR. N.LE.0)
     *                      IERROR = 1
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010   DO 1040 I=1,L
        DO 1030 J=1,N
        X = 0.0D0
        DO 1020 K=1,M
        X = X + A(I,K)*B(K,J)
 1020   CONTINUE
        C(I,J) = X
 1030   CONTINUE
 1040   CONTINUE
        RETURN
        END
