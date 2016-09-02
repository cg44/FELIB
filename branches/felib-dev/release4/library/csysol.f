C
      SUBROUTINE CSYSOL(KB,IKB,JKB,LOADS,ILOADS,N,HBAND,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CSYSOL solves a set of complex symmetric banded equations with a
C      single right hand side using a symmetric decomposition. Only
C      the lower band and diagonal are stored in a rectangular array KB.
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
C      Release 2.0  29 Oct 1979 (CRIE)
C      Commented    12 Feb 1985 (CG)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      KB      on entry contains lower half of complex symmetric
C              band matrix stored as a rectangular array
C      IKB     first dimension of KB (.GE. N)
C      JKB     second dimension of KB (.GE. HBAND)
C      LOADS   contains elements of right hand side
C      ILOADS  dimension of LOADS (.GE. N)
C      N       order of matrix KB
C      HBAND   semi-bandwidth of KB (includes diagonal)
C      ITEST   error checking option
C
C ARGUMENTS out
C      KB      on exit, contains lower triangular reduced
C      LOADS   matrix solution vector
C
C ROUTINES called
C      ERRMES CSYRDN CSYSUB
C
C     SUBROUTINE CSYSOL(KB,IKB,JKB,LOADS,ILOADS,N,HBAND,ITEST)
C*************************************************************************
C
      INTEGER HBAND,IKB,ITEST,JKB,JTEST,IERROR,N,ERRMES,ILOADS
      CHARACTER*6 SRNAME
      DOUBLE PRECISION KB,LOADS
      DIMENSION KB(2,IKB,JKB),LOADS(2,ILOADS)
C
      EXTERNAL ERRMES,CSYRDN,CSYSUB
C
      DATA SRNAME/'CSYSOL'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (ILOADS.LT.N) IERROR = 3
         IF (IKB.LT.N .OR. JKB.LT.HBAND) IERROR = 2
         IF (N.LE.0 .OR. HBAND.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      JTEST = 1
      CALL CSYRDN(KB,IKB,JKB,N,HBAND,JTEST)
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = JTEST
         IF (JTEST.EQ.3) IERROR = 4
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
      CALL CSYSUB(KB,IKB,JKB,LOADS,ILOADS,N,HBAND,JTEST)
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = JTEST
         IF (JTEST.EQ.3) IERROR = 4
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.EQ.0) RETURN
      END IF
C
C
      END
