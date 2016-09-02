C
      SUBROUTINE SELECT(V,IV,STEER,ISTEER,N,W,IW,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      SELECT constructs an element value vector from a full system
C      vector using the steering vector
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
C      Release 2.0  29 Oct 1985 (CG)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      V       vector of system values
C      IV      first dimension of vector V
C      STEER   element steering vector - contains freedom nos.
C      ISTEER  first dimension of vector STEER (.GE. N)
C      N       number of freedoms to be assembled
C      IW      first dimension of vector W (.GE. N)
C      ITEST   error checking option
C
C ARGUMENTS out
C      W       vector of element values
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE SELECT(V,IV,STEER,ISTEER,N,W,IW,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IERROR,ISTEER,ITEST,IV,IW,J,N,STEER,JTEST
      CHARACTER*6 SRNAME
      DOUBLE PRECISION V,W
      DIMENSION STEER(ISTEER),V(IV),W(IW)
C
      DATA SRNAME/'SELECT'/
C
C     Parameter checking
C
      JTEST = ITEST
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (ISTEER.LT.N .OR. IW.LT.N) IERROR = 2
         IF (N.LE.0) IERROR = 1
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main loops
C
      DO 1000 I = 1,N
         W(I) = 0.D0
         J = IABS(STEER(I))
         IF (J.NE.0) THEN
C
            IF (JTEST.NE.-1) THEN
               IERROR = 0
               IF (J.GT.IV) IERROR = 3
               ITEST = ERRMES(JTEST,IERROR,SRNAME)
               IF (ITEST.NE.0) RETURN
            END IF
C
            W(I) = V(J)
         END IF
 1000 CONTINUE
C
      END
