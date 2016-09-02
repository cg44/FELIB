C***********************************************************************
C$SPLIT$DIRECT$*********************************************************
C***********************************************************************
      SUBROUTINE DIRECT(NELE, ELTOP, IELTOP, JELTOP, NF, INF, JNF,
     *     DOFNOD, STEER, ISTEER, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CONSTRUCTS THE STEERING VECTOR TO DIRECT ASSEMBLY OF A
C      SYSTEM MATRIX
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    13 FEB 1980 (KR)
C
C ARGUMENTS IN
C      NELE    ELEMENT NUMBER
C      ELTOP   2D ARRAY CONTAINING ELEMENT TYPE, NUMBER OF
C              NODES ON THE ELEMENT, AND THE ELEMENT TOPOLOGIES
C      IELTOP  FIRST DIMENSION OF ELTOP (.GE. NUMBER OF
C              ELEMENTS IN PROBLEM)
C      JELTOP  SECOND DIMENSION OF ELTOP (.GE. NUMBER OF NODES
C              ON ELEMENT + 2)
C      NF      CONTAINS FREEDOM NUMBERS ASSOCIATED WITH EACH
C              NODE
C      INF     FIRST DIMENSION OF NF (.GE. TOTAL NUMBER OF
C              NODES IN PROBLEM)
C      JNF     SECOND DIMENSION OF NF (.GE. DOFNOD)
C      DOFNOD  NUMBER OF DEGREES OF FREEDOM AT EACH NODE
C      ISTEER  DIMENSION OF VECTOR STEER (.GE. TOTAL NUMBER OF
C              DEGREES OF FREEDOM ON ELEMENT)
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      STEER   VECTOR CONTAINING FREEDOM NUMBERS ASSOCIATED
C              WITH ELEMENT NELE, ARRANGED IN LOCAL NODE ORDER
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE DIRECT(NELE, ELTOP, IELTOP, JELTOP, NF, INF, JNF,
C    *     DOFNOD, STEER, ISTEER, ITEST)
C***********************************************************************
C
      INTEGER DOFNOD, ELTOP, ERRMES, IDEG, IELTOP, IERROR,
     *     INF, INOD, ISTEER, ITEST, JELTOP, JNF, JNOD, K,
     *     NELE, NF, NODEL, STEER
      DOUBLE PRECISION SRNAME
      DIMENSION ELTOP(IELTOP,JELTOP), NF(INF,JNF),
     *     STEER(ISTEER)
      DATA SRNAME /8H DIRECT /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (JNF.LT.DOFNOD) IERROR = 3
                        IF (IELTOP.LT.NELE) IERROR = 2
                        IF (NELE.LE.0 .OR. DOFNOD.LE.0) IERROR = 1
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 NODEL = ELTOP(NELE,2)
                        IF (ITEST.EQ.-1) GO TO 1020
                        IERROR = 0
                        IF (JELTOP.LT.NODEL+2) IERROR = 4
                        IF (ISTEER.LT.DOFNOD*NODEL) IERROR = 6
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1020 K = 1
      DO 1050 INOD=1,NODEL
      JNOD = ELTOP(NELE,INOD+2)
                        IF (ITEST.EQ.-1) GO TO 1030
                        IERROR = 0
                        IF (INF.LT.JNOD) IERROR = 5
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1030 DO 1040 IDEG=1,DOFNOD
      STEER(K) = NF(JNOD,IDEG)
      K = K + 1
 1040 CONTINUE
 1050 CONTINUE
      RETURN
      END
