C***********************************************************************
      SUBROUTINE UPDATE(PHI, IPHI, RHS, IRHS, TOTNOD, DOFNOD,
     *     TOTDOF, NF, INF, JNF, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      THIS ROUTINE TAKES A FULL SOLUTION VECTOR AND A SET OF
C      UPDATES AND UPDATES THE SOLUTION VECTOR
C
C METHOD
C      THE ROUTINE ASSUMES THAT THE SOLUTION VECTOR IS FULL
C      AND THAT THE UPDATE VECTOR IS IMCOMPLETE. THE UPDATE
C      IS PERFORMED USING THE NODAL FREEDOM ARRAY TO DIRECT
C      UPDATE TO FREEDOM MAPPING
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 3.0  29 JUN 1986 (CJH,CG)
C
C ARGUMENTS IN
C      PHI     SOLUTION VECTOR
C      IPHI    DIMENSION OF VECTOR PHI (IPHI.GE.TOTNOD*DOFNOD)
C      RHS     VECTOR OF UPDATES
C      IRHS    DIMENSION OF VECTOR RHS (IRHS.GE.TOTDOF)
C      TOTNOD  THE NUMBER OF NODES IN THE PROBLEM
C      DOFNOD  THE MAXIMUM NUMBER OF NODES PER NODE
C      TOTDOF  THE NUMBER OF FREEDOMS IN RHS
C      NF      THE NODAL FREEDOM ARRAY
C      INF     FIRST DIMENSION OF NF (INF.GE.TOTNOD)
C      JNF     SECOND DIMENSION OF NF (JNF.GE.DOFNOD)
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      PHI     THE UPDATED SOLUTION VECTOR
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE UPDATE(PHI, IPHI, RHS, IRHS, TOTNOD, DOFNOD,
C    *     TOTDOF, NF, INF, JNF, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, DOFNOD, I, IERROR, INF, IPHI, IRHS, ITEST, J,
     *     JNF, K, L, NF, TOTDOF, TOTNOD
      DOUBLE PRECISION PHI, RHS, SRNAME
      DIMENSION NF(INF,JNF), PHI(IPHI), RHS(IRHS)
      DATA SRNAME /8H UPDATE /
C
C                            CHECK ARRAY BOUNDS
C
      IF (ITEST.EQ.-1) GO TO 1010
      IERROR = 0
      IF (IPHI.LT.TOTNOD*DOFNOD) IERROR = 1
      IF (IRHS.LT.TOTDOF) IERROR = 2
      IF (INF.LT.TOTNOD) IERROR = 3
      IF (JNF.LT.DOFNOD) IERROR = 4
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF (ITEST.NE.0) RETURN
C
 1010 DO 1030 I=1,TOTNOD
C
      DO 1020 J=1,DOFNOD
C
      K = NF(I,J)
C
C                            DON'T UPDATE RESTRAINED VARIABLES
C
      IF (K.EQ.0) GO TO 1020
C
      L = I + J - 1
C
C                            UPDATE SOLUTION VECTOR
C
      PHI(L) = PHI(L) + RHS(K)
C
 1020 CONTINUE
 1030 CONTINUE
C
      RETURN
      END