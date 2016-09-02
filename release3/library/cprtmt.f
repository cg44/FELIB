C***********************************************************************
      SUBROUTINE CPRTMT(A, IA, JA, M, N, NOUT, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      PRINTS A COMPLEX MATRIX IN A STANDARD FORMAT
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 3.0   29 OCT 1984 (CRIE)
C      COMMENTS       1 NOV 1985  (CG)
C
C ARGUMENTS IN
C      A       MATRIX OF DIMENSION 2,IA,JA CONTAINING NUMBERS TO BE
C              PRINTED
C      IA      FIRST DIMENSION OF A (.GE.M)
C      JA      SECOND DIMENSION OF A (.GE.N)
C      M       NUMBER OF ROWS OF A TO BE PRINTED (.LE.IA)
C      N       NUMBER OF COLUMES OF A TO BE PRINTED (.LE.JA)
C      NOUT    FORTRAN UNIT NUMBER
C      ITEST   ERROR CHECKING OPTION
C
C ROUTINES CALLED
C      ERRMES
C
C     SUBROUTINE CPRTMT(A, IA, JA, M, N, NOUT, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IERROR, ITEST, IA, N, NOUT, J, JA, K
      DOUBLE PRECISION SRNAME, A
      DIMENSION A(2,IA,JA)
      DATA SRNAME /8H CPRTMT /
C
C     PARAMETER CHECKING
C
      IF (ITEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (IA.LT.M .OR. JA.LT.N) IERROR = 2
      IF (N.LE.0 .OR. M.LE.0) IERROR = 1
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
C     MAIN LOOPS
C
 1010 WRITE (NOUT,9010) (((A(I,J,K),I=1,2),J=1,M),K=1,N)
      RETURN
C
 9010 FORMAT (1H , 3(2X,2D12.4))
      END
C***********************************************************************
