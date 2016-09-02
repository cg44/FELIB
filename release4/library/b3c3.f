C
      SUBROUTINE B3C3(B,IB,JB,DER,IDER,JDER,NODEL,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      B3C3 forms strain-displacement matrix for 3d elasticity.
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
C      IB      first dimension of array B (.GE. 6)
C      JB      second dimension of B (.GE. 3*NODEL)
C      DER     DER(I,J) contains the derivative of the j'th
C              shape function with respect to the i'th global
C              coordinate
C      IDER    first dimension of DER (.GE. 3)
C      JDER    second dimension of DER (.GE. NODEL)
C      NODEL   number of nodes on element
C      ITEST   error checking option
C
C ARGUMENTS out
C      B       contains values of strain-displacement matrix
C
C ROUTINES called
C      MATNUL  ERRMES
C
C     SUBROUTINE B3C3(B,IB,JB,DER,IDER,JDER,NODEL,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,IB,IDER,IERROR,ITEST,JB,JDER,K,L,M,N,NODEL
      CHARACTER*6 SRNAME
      DOUBLE PRECISION B,DER
      DIMENSION B(IB,JB),DER(IDER,JDER)
C
      EXTERNAL ERRMES,MATNUL
C
      DATA SRNAME/'B3C3'/
C
C     Intialisation
C
      K = 6
      L = 3*NODEL
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IDER.LT.3 .OR. JDER.LT.NODEL) IERROR = 3
         IF (IB.LT.K .OR. JB.LT.L) IERROR = 2
         IF (NODEL.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main loops
C
      CALL MATNUL(B,IB,JB,K,L,ITEST)
C
      DO 1000 M = 1,NODEL
         N = 3*M
         K = N - 1
         L = K - 1
         B(1,L) = DER(1,M)
         B(4,K) = DER(1,M)
         B(6,N) = DER(1,M)
         B(2,K) = DER(2,M)
         B(4,L) = DER(2,M)
         B(5,N) = DER(2,M)
         B(3,N) = DER(3,M)
         B(5,K) = DER(3,M)
         B(6,L) = DER(3,M)
 1000 CONTINUE
C
      END
