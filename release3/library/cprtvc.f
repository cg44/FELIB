C***********************************************************************
      SUBROUTINE CPRTVC(V, IV, N, NOUT, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      PRINTS A COMPLEX VECTOR IN A STANDARD FORMAT
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
C      V       VECTOR OF DIMENSION 2,IV CONTAINING NUMBERS TO BE
C              PRINTED
C      IV      DIMENSION OF V (.GE.N)
C      N       NUMBER OF ELEMENTS OF V TO BE PRINTED
C      NOUT    FORTRAN UNIT NUMBER
C      ITEST   ERROR CHECKING OPTION
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE CPRTVC(V, IV, N, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IERROR, ITEST, IV, N, NOUT, J
      DOUBLE PRECISION SRNAME, V
      DIMENSION V(2,IV)
      DATA SRNAME /8H CPRTVC /
C
C     PARAMETER CHECKING
C
      IF (ITEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (IV.LT.N) IERROR = 2
      IF (N.LE.0) IERROR = 1
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
C     MAIN LOOPS
C
 1010 WRITE (NOUT,9010) ((V(I,J),I=1,2),J=1,N)
      RETURN
C
 9010 FORMAT (1H , 3(2X,2D12.4))
      END
C***********************************************************************
