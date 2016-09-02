C***********************************************************************
C$SPLIT$MATVEC$*********************************************************
C***********************************************************************
      SUBROUTINE MATVEC(A, IA, JA, V, IV, M, N, W, IW, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      POST-MULTIPLIES THE MATRIX A BY THE VECTOR V, STORING
C      THE RESULT IN THE VECTOR W
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (IMS) --- SERC COPYRIGHT
C      COMMENTED    14 OCT 1980 (KR)
C
C ARGUMENTS IN
C      A       ARRAY OF DIMENSION (IA,JA)
C      IA      FIRST DIMENSION OF A (.GE.M)
C      JA      SECOND DIMENSION OF A (.GE.N)
C      V       VECTOR OF DIMENSION IV
C      IV      DIMENSION OF V (.GE.N)
C      M       NUMBER OF ROWS OF A TO BE USED IN THE
C              MULTIPLICATION
C      N       NUMBER OF COLUMNS OF A AND THE NUMBER OF
C              ELEMENETS OF V TO BE USED IN THE MULTIPLICATION
C      IW      DIMENSION OF VECTOR W (.GE.M)
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      W       VECTOR OF DIMENSION IW; CONTAINS THE RESULT OF
C              THE OPERATION W=A*V
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE MATVEC(A, IA, JA, V, IV, M, N, W, IW, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IA, IERROR, ITEST, IV, IW, J, JA,
     *     M, N
      DOUBLE PRECISION A, SRNAME, V, W, X
      DIMENSION A(IA,JA), V(IV), W(IW)
      DATA SRNAME /8H MATVEC /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (M.GT.IW) IERROR = 4
                        IF (N.GT.IV) IERROR = 3
                        IF (M.GT.IA .OR. N.GT.JA) IERROR = 2
                        IF (M.LE.0 .OR. N.LE.0) IERROR = 1
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (IERROR.NE.0) RETURN
 1010                   DO 1030 I=1,M
                        X = 0.0D0
                        DO 1020 J=1,N
                        X = X + A(I,J)*V(J)
 1020                   CONTINUE
                        W(I) = X
 1030                   CONTINUE
                        RETURN
                        END
