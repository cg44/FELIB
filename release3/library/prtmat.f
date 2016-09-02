C***********************************************************************
      SUBROUTINE PRTMAT(A, IA, JA, M, N, NOUT, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      PRINTS A TWO-DIMENSIONAL ARRAY IN A STANDARD FORMAT
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 1.1  29 OCT 1979 (CG)
C      COMMENTED    14 OCT 1980 (KR)
C
C ARGUMENTS IN
C      A       ARRAY OF DIIMENSION (IA,JA) CONTAINING NUMBERS
C              TO BE PRINTED
C      IA      FIRST DIMENSION OF A (.GE.M)
C      JA      SECOND DIMENSION OF A (.GE.N)
C      M       NUMBER OF ROWS OF A TO BE PRINTED
C      N       NUMBER OF COLUMNS OF A TO BE PRINTED
C      NOUT    FORTRAN UNIT NUMBER
C      ITEST   ERROR CHECKING OPTION
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE PRTMAT(A, IA, JA, M, N, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IA, IERROR, ITEST, J, JA, M, N, NOUT
      DOUBLE PRECISION A, SRNAME
      DIMENSION A(IA,JA)
      DATA SRNAME /8H PRTMAT /
C
C     PARAMETER CHECKING
C
      IF (ITEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (IA.LT.M .OR. JA.LT.N) IERROR = 2
      IF (M.LE.0 .OR. N.LE.0) IERROR = 1
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
C     MAIN BODY
C
 1010 DO 1020 I=1,M
         WRITE (NOUT,9010) (A(I,J),J=1,N)
 1020 CONTINUE
C
      RETURN
 9010 FORMAT (1H , 6D12.4)
      END
C***********************************************************************
