C***********************************************************************
C$SPLIT$MATCOP$*********************************************************
C***********************************************************************
      SUBROUTINE MATCOP(A, IA, JA, B, IB, JB, M, N, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      COPIES MATRIX A INTO MATRIX B
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (IMS) --- SERC COPYRIGHT
C      COMMENTED    12 OCT 1980 (KR)
C
C ARGUMENTS IN
C      A       ARRAY OF DIMENSION (IA,JA) WHICH IS TO BE COPIED
C      IA      FIRST DIMENSION OF A (.GE.M)
C      JA      SECOND DIMENSION OF A (.GE.N)
C      IB      FIRST DIMENSION OF ARRAY B (.GE.M)
C      JB      SECOND DIMENSION OF B (.GE.N)
C      M       NUMBER OF ROWS OF A TO BE COPIED
C      N       NUMBER OF COLUMNS OF A TO BE COPIED
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      B       ARRAY OF DIMENSION (IB,JB) INTO WHICH A IS
C              COPIED
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE MATCOP(A, IA, JA, B, IB, JB, M, N, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IA, IB, IERROR, ITEST, J, JA, JB,
     *     M, N
      DOUBLE PRECISION A, B, SRNAME
      DIMENSION A(IA,JA), B(IB,JB)
      DATA SRNAME /8H MATCOP /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR=0
                        IF (M.GT.IB .OR. N.GT.JB) IERROR = 3
                        IF (M.GT.IA .OR. N.GT.JA) IERROR = 2
                        IF (M.LE.0 .OR. N.LE.0) IERROR = 1
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 DO 1030 I=1,M
      DO 1020 J=1,N
      B(I,J) = A(I,J)
 1020 CONTINUE
 1030 CONTINUE
      RETURN
      END
