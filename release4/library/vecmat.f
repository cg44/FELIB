C
      SUBROUTINE VECMAT(V,IV,A,IA,JA,M,N,W,IW,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      VECMAT pre-multiplies the matrix A by the vector V, storing
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
C      Release 2.0   1 Feb 1981 (CG)
C      Commented     1 Feb 1981 (CG)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      V       vector of dimension IV
C      IV      dimension of V (.GE. M)
C      A       array of dimension (IA, JA)
C      IA      first dimension of A (.GE. M)
C      JA      second dimension of A (.GE. N)
C      M       number of rows of A to be used in the
C              multiplication
C      N       number of columns of A and the number of
C              elemenets of V to be used in the multiplication
C      IW      dimension of vector W (.GE. N)
C      ITEST   error checking option
C
C ARGUMENTS out
C      W       vector of dimension IW; contains the result of
C              the operation W=V*A
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE VECMAT(V,IV,A,IA,JA,M,N,W,IW,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IA,IERROR,ITEST,IV,IW,J,JA,M,N
      CHARACTER*6 SRNAME
      DOUBLE PRECISION A,V,W,X
      DIMENSION A(IA,JA),V(IV),W(IW)
C
      DATA SRNAME/'VECMAT'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (N.GT.IW) IERROR = 4
         IF (M.GT.IV) IERROR = 3
         IF (M.GT.IA .OR. N.GT.JA) IERROR = 2
         IF (M.LE.0 .OR. N.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      DO 1010 I = 1,N
         X = 0.0D0
         DO 1000 J = 1,M
            X = X + A(J,I)*V(J)
 1000    CONTINUE
         W(I) = X
 1010 CONTINUE
C
      END
