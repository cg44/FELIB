C***********************************************************************
      SUBROUTINE DYAD(V, IV, W, IW, A, IA, JA, N, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      FORMS THE DYAD MATRIX A FROM TWO VECTORS V AND W
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 1.1  29 OCT 1979  (CG)
C      COMMENTED    14 FEB 1980  (KR)
C
C ARGUMENTS IN
C      V       FIRST VECTOR IN MULTIPCATION
C      IV      DIMENSION OF VECTOR V (.GE. N)
C      W       SECOND VECTOR IN MULTIPLICATION
C      IW      DIMENSION OF VECTOR W (.GE. N)
C      IA      FIRST DIMENSION OF DYADIC ARRAY A (.GE. N)
C      JA      SECOND DIMENSION OF A (.GE. N)
C      N       NUMBER OF ELEMENTS OF VECTORS V AND W TO BE USED
C              IN FORMING THE DYAD
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      A       DYADIC ARRAY FORMED BY PRODUCT OF V AND W
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE DYAD(V, IV, W, IW, A, IA, JA, N, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IA, IERROR, ITEST, IV, IW, J, JA,
     *     N
      DOUBLE PRECISION A, SRNAME, V, W
      DIMENSION A(IA,JA), V(IV), W(IW)
      DATA SRNAME /8H DYAD   /
C
C     PARAMETER CHECKING
C
      IF (ITEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (N.GT.IW) IERROR = 3
      IF (N.GT.IV) IERROR = 2
      IF (N.LE.0) IERROR = 1
      IF (N.GT.IA .OR. N.GT.JA) IERROR = 4
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
C     MAIN BODY
C
 1010 DO 1030 J=1,N
         DO 1020 I=1,N
            A(I,J) = V(I)*W(J)
 1020    CONTINUE
 1030 CONTINUE
C
      RETURN
      END
C***********************************************************************
