C***********************************************************************
      SUBROUTINE PRTVEC(V, IV, N, NOUT, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      PRINTS A VECTOR IN A STANDARD FORMAT
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
C      V       VECTOR OF DIMENSION IV CONTAINING NUMBERS TO BE
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
C     SUBROUTINE PRTVEC(V, IV, N, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IERROR, ITEST, IV, N, NOUT
      DOUBLE PRECISION SRNAME, V
      DIMENSION V(IV)
      DATA SRNAME /8H PRTVEC /
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
C     MAIN BODY
C
 1010 WRITE (NOUT,9010) (V(I),I=1,N)
      RETURN
C
 9010 FORMAT (1H , 6D12.4)
      END
C***********************************************************************
