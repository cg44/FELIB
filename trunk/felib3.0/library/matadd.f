C***********************************************************************
      SUBROUTINE MATADD(A, IA, JA, B, IB, JB, M, N, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      ADDS THE MATRIX A TO B, STORING THE RESULT IN A
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 1.1  29 OCT 1979 (IMS)
C      COMMENTED    12 OCT 1980 (KR)
C
C ARGUMENTS IN
C      A       ARRAY OF DIMENSION (IA,JA)
C      IA      FIRST DIMENSION OF A (.GE.M)
C      JA      SECOND DIMENSION OF A (.GE.N)
C      B       ARRAY OF DIMENSION (IB,JB)
C      IB      FIRST DIMENSION OF B (.GE.M)
C      JB      SECOND DIMENSION OF B (.GE.N)
C      M       NUMBER OF ROWS OF A AND B TO BE ADDED
C      N       NUMBER OF COLUMNS OF A AND B TO BE ADDED
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      A       ON EXIT, CONTAINS SUM OF A AND B
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE MATADD(A, IA, JA, B, IB, JB, M, N, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IA, IB, IERROR, ITEST, J, JA, JB,
     *     M, N
      DOUBLE PRECISION A, B, SRNAME
      DIMENSION A(IA,JA), B(IB,JB)
      DATA SRNAME /8H MATADD /
C
C     PARAMETER CHECKING
C
      IF (ITEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (M.GT.IB .OR. N.GT.JB) IERROR = 3
      IF (M.GT.IA .OR. N.GT.JA) IERROR = 2
      IF (M.LE.0 .OR. N.LE.0) IERROR = 1
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
C     MAIN BODY
C
 1010 DO 1030 J=1,N
         DO 1020 I=1,M
            A(I,J) = A(I,J) + B(I,J)
 1020    CONTINUE
 1030 CONTINUE
      RETURN
C
      END
C***********************************************************************
