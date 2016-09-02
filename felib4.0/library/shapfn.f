C
      SUBROUTINE SHAPFN(N,IN,JN,FUN,IFUN,NODEL,NDE,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      SHAPFN constructs the shape function matrix N for a system of coupled
C      differential equations.  This routine is most frequently used
C      in forming the consistent element and system mass matrices
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
C      Commented    22 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS IN
C      IN      first dimension of array N (.GE. NDE)
C      JN      second dimension of N (.GE. NODEL*NDE)
C      FUN     vector of length IFUN. Contains values of the
C              shape functions at the point where shape
C              function matrix required
C      IFUN    length of vector FUN (.GE. NODEL)
C      NODEL   number of shape functions to be used in
C              constructing shape function matrix N (usually
C              number of nodes IN element under consideration)
C      NDE     number of coupled equations being considered
C      ITEST   error checking option
C
C ARGUMENTS out
C      N       array of dimension (IN, JN).  contains the
C              values of the shape function matrix
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE SHAPFN(N,IN,JN,FUN,IFUN,NODEL,NDE,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IERROR,IFUN,IJ,IN,ITEST,J,JN,NDE,NODEL
      CHARACTER*6 SRNAME
      DOUBLE PRECISION FUN,N
      DIMENSION FUN(IFUN),N(IN,JN)
      DATA SRNAME/'SHAPFN'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (NODEL.LE.0 .OR. NDE.LE.0) IERROR = 1
         IF (IN.LT.NDE .OR. JN.LT.NODEL*NDE) IERROR = 2
         IF (IFUN.LT.NODEL) IERROR = 3
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      CALL MATNUL(N,IN,JN,NDE,NODEL*NDE,ITEST)
C
      DO 1010 J = 1,NODEL
         DO 1000 I = 1,NDE
            IJ = (J-1)*NDE + I
            N(I,IJ) = FUN(J)
 1000    CONTINUE
 1010 CONTINUE
C
      END
