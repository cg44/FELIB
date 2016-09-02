C
      SUBROUTINE CHOFWD(A,IA,JA,R,IR,N,HBAND,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CHOFWD performs forward substitution on a matrix reduced by
C      CHORDN
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
C              definite symmetric band matrix of order N and
C              semi-bandwidth HBAND.  A should previously have
C              been reduced using CHORDN or CHOSOL
C      IA      first dimension of A (.GE. N)
C      JA      second dimension of A (.GE. HBAND)
C      R       on entry, contains the vector of rhs's
C      IR      dimension of R (.GE. N)
C      N       order of matrix A
C      HBAND   semi-bandwidth of A
C      ITEST   error checking option
C
C ARGUMENTS out
C      R       on exit, contains the reduced vector of rhs's
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE CHOFWD(A,IA,JA,R,IR,N,HBAND,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,HBAND,I,IA,IERROR,IJ,IR,ITEST,J,JA,K,N,W
      CHARACTER*6 SRNAME
      DOUBLE PRECISION A,R,X
      DIMENSION A(IA,JA),R(IR)
      DATA SRNAME/'CHOFWD'/
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
      R(1) = R(1)/A(1,W+1)
      DO 1010 I = 2,N
         X = 0.0D0
         K = 1
         IF (I.LE.W+1) K = W - I + 2
         DO 1000 J = K,W
            IJ = I + J - W - 1
            X = X + A(I,J)*R(IJ)
 1000    CONTINUE
         R(I) = (R(I)-X)/A(I,W+1)
 1010 CONTINUE
C
      END
