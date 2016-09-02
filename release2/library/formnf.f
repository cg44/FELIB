C***********************************************************************
C$SPLIT$FORMNF$*********************************************************
C***********************************************************************
      SUBROUTINE FORMNF(REST, IREST, JREST, RESNOD, TOTNOD, DOFNOD, NF,
     *     INF, JNF, TOTDOF, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CONSTRUCTS THE NODAL FREEDOM ARRAY FROM THE RESTRAINED
C      FREEDOM DATA
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979
C      COMMENTED    18 FEB 1980 (KR)
C      RECODED      01 NOV 1981 (NB)
C
C ARGUMENTS IN
C      REST    INTEGER ARRAY; REST(I,J) CONTAINS THE I'TH SET
C              OF RESTRAINT INFORMATION - REST(I,1) CONTAINS
C              THE NODE NUMBER, REST(I,J+1) FOR J=1(1)DOFNOD
C              CONTAINS THE LOCAL FREEDOM NUMBERS OF THE
C              FREEDOMS THAT ARE RESTRAINED
C      IREST   FIRST DIMENSION OF ARRAY REST (.GE. RESNOD)
C      JREST   SECOND DIMENSION OF REST (.GE. DOFNOD)
C      RESNOD  NUMBER OF NODES AT WHICH FREEDOMS ARE RESTRAINED
C      TOTNOD  TOTAL NUMBER OF NODES IN MESH
C      INF     FIRST DIMENSION OF ARRAY INF (.GE. TOTNOD)
C      JNF     SECOND DIMENSION OF INF (.GE. DOFNOD)
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      NF      NF(I,J), J=1(1)DOFNOD, CONTAINS THE FREEDOM
C              NUMBERS ASSOCIATED WITH THE I'TH NODE
C      TOTDOF  TOTAL NUMBER OF FREEDOMS IN PROBLEM UNDER
C              CONSIDERATION
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE FORMNF(REST, IREST, JREST, RESNOD, TOTNOD, DOFNOD, NF,
C    *     INF, JNF, TOTDOF, ITEST)
C***********************************************************************
C
      INTEGER DOFNOD, I, INF, J, K, L, NF, REST, TOTNOD, RESNOD, TOTDOF,
     *     IERROR,ERRMES,M
      DOUBLE PRECISION SRNAME
      DIMENSION NF(INF,JNF), REST(IREST,JREST)
      LOGICAL SWITCH
      DATA SRNAME /8H FORMNF /
                        IF(ITEST.EQ.-1) GO TO 999
                        IERROR=0
                        IF(INF.LT.TOTNOD.OR.JNF.LT.DOFNOD) IERROR=3
                        IF(IREST.LT.RESNOD.OR.JREST.LT.DOFNOD+1)
     *                     IERROR=2
                        IF(RESNOD.LT.0.OR.TOTNOD.LE.0.OR.DOFNOD.LE.0)
     *                     IERROR=1
                        ITEST=ERRMES(ITEST,IERROR,SRNAME)
                        IF(ITEST.NE.0) RETURN
      SWITCH=.TRUE.
999   DO 1020 I=1,TOTNOD
      DO 1010 J=1,DOFNOD
      NF(I,J) = 1
 1010 CONTINUE
 1020 CONTINUE
      IF (RESNOD.EQ.0) GO TO 1045
      DO 1040 I=1,RESNOD
      K = REST(I,1)
      DO 1030 J=1,DOFNOD
      L = REST(I,J+1)
      M=IABS(L)
                        IF(ITEST.EQ.-1) GO TO 888
                        IERROR=0
                        IF(K.GT.TOTNOD.OR.M.GT.DOFNOD) IERROR=4
                        ITEST=ERRMES(ITEST,IERROR,SRNAME)
                        IF(ITEST.NE.0) RETURN
888   IF (L.GT.0) NF(K,L) = 0
      IF (L.GE.0) GO TO 1030
      NF(K,M)=L
      SWITCH=.FALSE.
 1030 CONTINUE
 1040 CONTINUE
 1045 CONTINUE
      K = 1
      DO 1060 I=1,TOTNOD
      DO 1050 J=1,DOFNOD
      IF (NF(I,J).LE.0) GO TO 1050
      NF(I,J) = K
      K = K + 1
 1050 CONTINUE
 1060 CONTINUE
      TOTDOF = K - 1
      IF (SWITCH) RETURN
      DO 1080 I=1,TOTNOD
      DO 1070 J=1,DOFNOD
      IF (NF(I,J).GE.0) GO TO 1070
      NF(I,J)=-K
      K=K+1
 1070 CONTINUE
 1080 CONTINUE
      RETURN
      END
