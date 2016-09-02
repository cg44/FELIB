C
      SUBROUTINE CHOBAK(A,IA,JA,R,IR,N,HBAND,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CHOBAK performs bakward substitution on a matrix processed by
C      CHORDN and CHOFWD
C
C HISTORY
C
C      Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
C                           Chilton, Didcot, Oxfordshire OX11 0QX
C
C      Contact: Prof Chris Greenough
C      Email: c.greenough@rl.ac.uk
C      Tel: (44) 1235 445307
C
C      Release 1.1  29 Oct 1979 (CG)
C      Commented    10 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      A       array of dimension (IA, JA).  contains the
C              elements of the lower half of the positive
C              definite band matrix of order N and with semi-
C              bandwidth HBAND, reduced by CHORDN or CHOSOL
C      IA      first dimension of A (.GE. N)
C      JA      second dimension of A (.GE. HBAND)
C      R       on entry, contains elements of N rhs's after
C              processing by CHOFWD
C      IR      dimension of vector R (.GE. N)
C      N       order of matrix A
C      HBAND   semi-bandwidth of matrix A
C      ITEST   error checking option
C
C ARGUMENTS out
C      R       on exit, contains solution vector
C
C ROUTINES called
C      ERRMES
C
C
C     SUBROUTINE CHOBAK(A,IA,JA,R,IR,N,HBAND,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,HBAND,I,IA,IERROR,IJ,IR,ITEST,J,JA,L,M,N,W
      CHARACTER*6 SRNAME
      DOUBLE PRECISION A,R,X
      DIMENSION A(IA,JA),R(IR)
      DATA SRNAME/'CHOBAK'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IR.LT.N) IERROR = 3
         IF (IA.LT.N .OR. JA.LT.HBAND) IERROR = 2
         IF (N.LE.0 .OR. HBAND.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      W = HBAND - 1
      R(N) = R(N)/A(N,W+1)
C
C     If N greater than 1
C
      IF(I.GT.1) THEN
      I = N - 1
 1000 CONTINUE
      X = 0.0D0
      L = I + W
      IF (I.GT.N-W) L = N
      M = I + 1
      DO 1010 J = M,L
         IJ = W + I - J + 1
         X = X + A(J,IJ)*R(J)
 1010 CONTINUE
      R(I) = (R(I)-X)/A(I,W+1)
      I = I - 1
      IF (I.NE.0) GO TO 1000
C
      ENDIF
C
      END
