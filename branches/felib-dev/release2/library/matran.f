C***********************************************************************
C$SPLIT$MATRAN$*********************************************************
C***********************************************************************
      SUBROUTINE MATRAN(A, IA, JA, B, IB, JB, M, N, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      FORMS THE TRANSPOSE OF THE MATRIX A
C
C HISTORY
C      RELEASE 1.1  20 OCT 1979 (IMS) --- SERC COPYRIGHT
C      COMMENTED    12 OCT 1980 (KR)
C
C ARGUMENTS IN
C      A       ARRAY OF DIMENSION (IA,JA) TO BE TRANSPOSED
C      IA      FIRST DIMENSION OF A (.GE.M)
C      JA      SECOND DIMENSION OF A (.GE.N)
C      IB      FIRST DIMENSION OF ARRAY B (.GE.N)
C      JB      SECOND DIMENSION OF ARRAY B (.GE.M)
C      M       NUMBER OF ROWS OF A TO BE TRANSPOSED
C      N       NUMBER OF COLUMNS OF A TO BE TRANSPOSED
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      B       ARRAY OF DIMENSION (IB,JB).  B(J,I)=A(I,J) FOR
C              I=1(1)M AND J=1(1)N
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE MATRAN(A, IA, JA, B, IB, JB, M, N, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IA, IB, IERROR, ITEST, J, JA, JB,
     *     M, N
      DOUBLE PRECISION A, B, SRNAME
      DIMENSION A(IA,JA), B(IB,JB)
      DATA SRNAME /8H MATRAN /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (N.GT.IB .OR. M.GT.JB) IERROR = 3
                        IF (M.GT.IA .OR. N.GT.JA) IERROR = 2
                        IF (M.LE.0 .OR. N.LE.0) IERROR = 1
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (IERROR.NE.0) RETURN
 1010                   DO 1030 I=1,M
                        DO 1020 J=1,N
                        B(J,I) = A(I,J)
 1020                   CONTINUE
 1030                   CONTINUE
                        RETURN
                        END
