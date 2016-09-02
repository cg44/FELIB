C***********************************************************************
C$SPLIT$B3C3$*********************************************************
C***********************************************************************
      SUBROUTINE B3C3(B, IB, JB, DER, IDER, JDER, NODEL, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      FORMS STRAIN-DISPLACEMENT MATRIX FOR 3D ELASTICITY
C
C HISTORY
C
C      COPYRIGHT (C) 1979 : SERC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 1.1  29 OCT 1979 (IMS)
C      COMMENTED    12 FEB 1980 (KR)
C
C ARGUMENTS IN
C      IB      FIRST DIMENSION OF ARRAY B (.GE. 6)
C      JB      SECOND DIMENSION OF B (.GE. 3*NODEL)
C      DER     DER(I,J) CONTAINS THE DERIVATIVE OF THE J'TH
C              SHAPE FUNCTION WITH RESPECT TO THE I'TH GLOBAL
C              COORDINATE
C      IDER    FIRST DIMENSION OF DER (.GE. 3)
C      JDER    SECOND DIMENSION OF DER (.GE. NODEL)
C      NODEL   NUMBER OF NODES ON ELEMENT
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      B       CONTAINS VALUES OF STRAIN-DISPLACEMENT MATRIX
C
C ROUTINES CALLED
C      MATNUL  ERRMES
C
C
C     SUBROUTINE B3C3(B, IB, JB, DER, IDER, JDER, NODEL, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, IB, IDER, IERROR, ITEST, JB, JDER, K,
     *     L, M, N, NODEL
      DOUBLE PRECISION B, DER, SRNAME
      DIMENSION B(IB,JB), DER(IDER,JDER)
      DATA SRNAME /8H B3C3   /
C
C     INTIALISATION
C
      K = 6
      L = 3*NODEL
C
C     PARAMETER CHECKING
C
      IF (ITEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (IDER.LT.3 .OR. JDER.LT.NODEL) IERROR = 3
      IF (IB.LT.K .OR. JB.LT.L) IERROR = 2
      IF (NODEL.LE.0) IERROR = 1
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
C     MAIN LOOPS
C
 1010 CALL MATNUL(B, IB, JB, K, L, ITEST)
C
      DO 1020 M=1,NODEL
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
 1020 CONTINUE
C
      RETURN
      END
