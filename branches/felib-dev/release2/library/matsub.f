C***********************************************************************
C$SPLIT$MATSUB$*********************************************************
C***********************************************************************
      SUBROUTINE MATSUB(A, IA, JA, B, IB, JB, M, N, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      SUBTRACTS MATRIX B FROM A, STORING THE RESULT IN A
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (IMS) --- SERC COPYRIGHT
C      COMMENTED    12 OCT 1980 (KR)
C
C ARGUMENTS IN
C      A       ARRAY OF DIMENSION (IA,JA)
C      IA      FIRST DIMENSION OF A (.GE.M)
C      JA      SECOND DIMENSION OF A (.GE.N)
C      B       ARRAY OF DIMENSION (IB,JB) TO BE SUBTRACTED FROM
C              A
C      IB      FIRST DIMENSION OF B (.GE.M)
C      JB      SECOND DIMENSION OF B (.GE.N)
C      M       NUMBER OF ROWS OF A AND B TO BE SUBTRACTED
C      N       NUMBER OF COLUMNS OF A AND B TO BE SUBTRACTED
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      A       ARRAY OF DIMENSION (IA,JA).  ON EXIT, A(I,J)
C              IS SET TO A(I,J)-B(I,J) FOR I=1(1)M AND J=1(1)N
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE MATSUB(A, IA, JA, B, IB, JB, M, N, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IA, IB, IERROR, ITEST, J, JA, JB,
     *     M, N
      DOUBLE PRECISION A, B, SRNAME
      DIMENSION A(IA,JA), B(IB,JB)
      DATA SRNAME /8H MATSUB /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (M.GT.IB .OR. N.GT.JB) IERROR = 3
                        IF (M.GT.IA .OR. N.GT.JA) IERROR = 2
                        IF (M.LE.0 .OR. N.LE.0) IERROR = 1
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (IERROR.NE.0) RETURN
 1010                   DO 1030 I=1,M
                        DO 1020 J=1,N
                        A(I,J) = A(I,J) - B(I,J)
 1020                   CONTINUE
 1030                   CONTINUE
                        RETURN
                        END
