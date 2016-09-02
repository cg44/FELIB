C
      SUBROUTINE CHORDN(A,IA,JA,N,HBAND,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CHORDN performs choleski reduction on a real symmetric positive
C      definite banded matrix.  Only the lower band and
C      diagonal are stored in a rectangular array A.
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
C      Release 1.1  29 Oct 1979 (IMS)
C      Commented    12 Feb 1980 (CG)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      A       on entry, contains lower band and diagonal of
C              pd real symmetric matrix
C      IA      first dimension of A (.GE. N)
C      JA      second dimension of A (.GE. HBAND)
C      N       order of matrix A
C      HBAND   semi-bandwidth of A (including diagonal)
C      ITEST   error checking option
C
C ARGUMENTS out
C      A       reduced matrix L, where the input matrix A has
C              been reduced to triangular matrices L and lt
C              where A=L lt
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE CHORDN(A,IA,JA,N,HBAND,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,HBAND,I,IA,IERROR,IK,ITEST,J,JA,JTEST,K,L,LA,LB,LK,
     *        N,W
      CHARACTER*6 SRNAME
      DOUBLE PRECISION A,X
      DIMENSION A(IA,JA)
      DATA SRNAME/'CHORDN'/
C
C     Parameter checking
C
      JTEST = ITEST
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (IA.LT.N .OR. JA.LT.HBAND) IERROR = 2
         IF (N.LE.0 .OR. HBAND.LE.0) IERROR = 1
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      W = HBAND - 1
      DO 1030 I = 1,N
         X = 0.0D0
         DO 1000 J = 1,W
            X = X + A(I,J)*A(I,J)
 1000    CONTINUE
C
C     Check on value of A(I,J)
C
         IF (JTEST.NE.-1) THEN
            IERROR = 0
            IF ((A(I,W+1)-X).LE.0.0D0) IERROR = 3
            ITEST = ERRMES(JTEST,IERROR,SRNAME)
            IF (ITEST.NE.0) RETURN
         END IF
C
         A(I,W+1) = DSQRT(A(I,W+1)-X)
C
         DO 1020 K = 1,W
            X = 0.0D0
            IF (I+K.LE.N) THEN
               IF (K.NE.W) THEN
                  L = W - K
 1010             CONTINUE
                  IK = I + K
                  LK = L + K
                  X = X + A(IK,L)*A(I,LK)
                  L = L - 1
                  IF (L.NE.0) GO TO 1010
               END IF
               LA = I + K
               LB = W - K + 1
               A(LA,LB) = (A(LA,LB)-X)/A(I,W+1)
            END IF
 1020    CONTINUE
 1030 CONTINUE
C
      END
