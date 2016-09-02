C
      SUBROUTINE MATVEC(A,IA,JA,V,IV,M,N,W,IW,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      MATVEC post-multiplies the matrix A by the vector V, storing
C      the result in the vector W
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
C      Commented    14 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      A       array of dimension (IA, JA)
C      IA      first dimension of A (.GE. M)
C      JA      second dimension of A (.GE. N)
C      V       vector of dimension IV
C      IV      dimension of V (.GE. N)
C      M       number of rows of A to be used in the
C              multiplication
C      N       number of columns of A and the number of
C              elemenets of V to be used in the multiplication
C      IW      dimension of vector W (.GE. M)
C      ITEST   error checking option
C
C ARGUMENTS out
C      W       vector of dimension IW; contains the result of
C              the operation W=A*V
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE MATVEC(A,IA,JA,V,IV,M,N,W,IW,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IA,IERROR,ITEST,IV,IW,J,JA,M,N
      CHARACTER*6 SRNAME
      DOUBLE PRECISION A,V,W,X
      DIMENSION A(IA,JA),V(IV),W(IW)
      DATA SRNAME/'MATVEC'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (M.GT.IW) IERROR = 4
         IF (N.GT.IV) IERROR = 3
         IF (M.GT.IA .OR. N.GT.JA) IERROR = 2
         IF (M.LE.0 .OR. N.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      DO 1010 I = 1,M
         X = 0.0D0
         DO 1000 J = 1,N
            X = X + A(I,J)*V(J)
 1000    CONTINUE
         W(I) = X
 1010 CONTINUE
C
      END
