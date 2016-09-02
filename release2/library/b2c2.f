C***********************************************************************
C$SPLIT$B2C2$*********************************************************
C***********************************************************************
      SUBROUTINE B2C2(B, IB, JB, DER, IDER, JDER, NODEL, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      FORMS THE STRAIN-DISPLACEMENT MATRIX FOR 2D PLANE
C      ELASTICITY
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
C      IB      FIRST DIMENSION OF ARRAY B (.GE. 3)
C      JB      SECOND DIMENSION OF B (.GE. 2*NODEL)
C      DER     DER(I,J) CONTAINS THE DERIVATIVE OF THE J'TH
C              SHAPE FUNCTION WITH RESPECT TO THE I'TH GLOBAL
C              COORDINATE
C      IDER    FIRST DIMENSION OF DER (.GE. 2)
C      JDER    SECOND DIMENSION OF DER (.GE. NODEL)
C      NODEL   NUMBER OF NODES ON THE ELEMENT
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      B       CONTAINS THE VALUES OF THE STRAIN-DISPLACEMENT
C              MATRIX
C
C ROUTINES CALLED
C      MATNUL  ERRMES
C
C
C     SUBROUTINE B2C2(B, IB, JB, DER, IDER, JDER, NODEL, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, IB, IDER, IERROR, ITEST, JB, JDER, K,
     *     L, M, NODEL
      DOUBLE PRECISION B, DER, SRNAME
      DIMENSION B(IB,JB), DER(IDER,JDER)
      DATA SRNAME /8H B2C2   /
C
C     INTIALISATION
C
      K = 3
      L = 2*NODEL
C
C     PARAMETER CHECKING
C
      IF (ITEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (IDER.LT.2 .OR. JDER.LT.NODEL) IERROR = 3
      IF (IB.LT.3 .OR. JB.LT.L) IERROR = 2
      IF (NODEL.LE.0) IERROR = 1
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
C     INITIALISE MATRIX
C
 1010 CALL MATNUL(B, IB, JB, K, L, ITEST)
C
      DO 1020 M=1,NODEL
         K = 2*M
         L = K - 1
         B(1,L) = DER(1,M)
         B(3,K) = DER(1,M)
         B(2,K) = DER(2,M)
         B(3,L) = DER(2,M)
 1020 CONTINUE
C
      RETURN
      END
