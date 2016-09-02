C
      SUBROUTINE B2P2(B,IB,JB,DER,IDER,JDER,FUN,IFUN,COORD,ICOORD,
     *                JCOORD,NODEL,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      B2P2 forms the strain-displacement matrix for axisymmetric
C      elasticity.
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
C      Commented    12 Feb 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      IB      first dimension of array B (.GE. 3)
C      JB      second dimension of B (.GE. 2*NODEL)
C      DER     DER(I,J) contains the derivative of the j'th
C              shape function with respect to the i'th global
C              coordinate
C      IDER    first dimension of DER (.GE. 2)
C      JDER    second dimension of DER (.GE. NODEL)
C      FUN     FUN(I) conatins the value of the i'th shape
C              function at the point under consideration
C      IFUN    first dimension of FUN (.GE. 4)
C      COORD   COORD(I,J) contains the j'th global coordinate
C              of the i'th node
C      ICOORD  first dimension of COORD (.GE. number of nodes
C              in the mesh)
C      JCOORD  second dimension of COORD (.GE. dimensionality
C              of problem)
C      NODEL   number of nodes on element
C      ITEST   error checking option
C
C ARGUMENTS out
C      B       conatins values of strain-displacement matrix
C
C ROUTINES called
C      MATNUL  ERRMES  VTOL
C
C     SUBROUTINE B2P2(B,IB,JB,DER,IDER,JDER,FUN,IFUN,COORD,ICOORD,
C    *                JCOORD,NODEL,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,IB,ICOORD,IDER,IERROR,IFUN,ITEST,JB,JCOORD,JDER,
     *        JTEST,K,L,M,NODEL
      CHARACTER*6 SRNAME
      DOUBLE PRECISION B,COORD,DER,FUN,SUM,TOL,VTOL,X
      DIMENSION B(IB,JB),COORD(ICOORD,JCOORD),DER(IDER,JDER),FUN(IFUN)
C
      EXTERNAL ERRMES,MATNUL,VTOL
C
      DATA SRNAME/'B2P2'/
C
C     Intialisation
C
      K = 4
      L = 2*NODEL
      TOL = VTOL(X)
C
C     Parameter checking
C
      JTEST = ITEST
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (IDER.LT.2 .OR. JDER.LT.NODEL) IERROR = 3
         IF (IB.LT.3 .OR. JB.LT.L) IERROR = 2
         IF (NODEL.LE.0) IERROR = 1
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main loops
C
      CALL MATNUL(B,IB,JB,K,L,ITEST)
C
C     Calculate average radius
C
      SUM = 0.0D0
      DO 1000 K = 1,NODEL
         SUM = SUM + FUN(K)*COORD(K,1)
 1000 CONTINUE
C
C     Check value of SUM rbar
C
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (DABS(SUM).LT.TOL) IERROR = 4
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
      DO 1010 M = 1,NODEL
         K = 2*M
         L = K - 1
         B(1,L) = DER(1,M)
         B(3,K) = DER(1,M)
         B(2,K) = DER(2,M)
         B(3,L) = DER(2,M)
         B(4,L) = FUN(M)/SUM
 1010 CONTINUE
C
      END
