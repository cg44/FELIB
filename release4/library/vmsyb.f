C
      SUBROUTINE VMSYB(V,IV,A,IA,JA,W,IW,N,HBAND,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      VMSYB pre-multiplies a real unsymmetric banded matrix stored as
C      a rectangular array by a vector.
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
C      Commented    14 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      V       vector of dimension IV
C      IV      dimension of vector V (.GE. N)
C      A       array of dimension (IA, JA). Contains the
C              elements of the real unsymmetric band matrix
C              of order N and semi-bandwidth HBAND
C      IA      first dimension of A (.GE. N)
C      JA      second dimension of A (.GE. HBAND)
C      IW      dimension of vector W (.GE. N)
C      N       order of the real unsymmetric band matrix
C      HBAND   semi-bandwidth of the real unsymmetric band matrix
C      ITEST   error checking option
C
C ARGUMENTS out
C      W       vector of dimension IW. Contains the result of
C              the operation W=A*V
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE VMSYB(V,IV,A,IA,JA,W,IW,N,HBAND,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,HBAND,I,IERROR,IJ,ITEST,IV,IW,J,JI,N,IA,JA
      CHARACTER*6 SRNAME
      DOUBLE PRECISION A,V,W,X
      DIMENSION A(IA,JA),V(IV),W(IW)
C
      DATA SRNAME/'VMSYB'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IW.LT.N) IERROR = 4
         IF (IV.LT.N) IERROR = 3
         IF (IA.LT.N .OR. JA.LT.HBAND) IERROR = 2
         IF (N.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      DO 1020 I = 1,N
         X = 0.0D0
         J = HBAND
 1000    CONTINUE
         IF (I+J.GT.HBAND) THEN
            IJ = I + J - HBAND
            X = X + A(I,J)*V(IJ)
         END IF
         J = J - 1
         IF (J.NE.0) GO TO 1000
         J = HBAND - 1
 1010    CONTINUE
         IF (I-J.LT.N-HBAND+1) THEN
            JI = I - J + HBAND
            X = X + A(JI,J)*V(JI)
         END IF
         J = J - 1
         IF (J.NE.0) GO TO 1010
         W(I) = X
 1020 CONTINUE
C
      END
