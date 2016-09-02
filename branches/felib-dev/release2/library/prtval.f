C***********************************************************************
C$SPLIT$PRTVAL$*********************************************************
C***********************************************************************
      SUBROUTINE PRTVAL(VAL, IVAL, NF, INF, JNF, DOFNOD, TOTNOD,
     *     NOUT, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      PRINTS OUT THE NODAL VALUES OF THE SOLUTION IN A
C      STANDARD FORMAT
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    14 OCT 1980 (KR)
C      RECODED      01 NOV 1981 (NB)
C
C ARGUMENTS IN
C      VAL     VECTOR OF DIMENSION IVAL CONTAINING SOLUTION
C              VALUES,AND PRESCRIBED BOUNDARY VALUES
C      IVAL    DIMENSION OF VAL (.GE.TOTAL NUMBER OF FREEDOMS
C              IN SYSTEM)
C      NF      INTEGER ARRAY OF DIMENSION (INF,JNF) CONTAINING
C              FREEDOM NUMBERS ASSOCIATED WITH EACH NODE
C      INF     FIRST DIMENSION OF NF (.GE.TOTNOD)
C      JNF     SECOND DIMENSION OF NF (.GE.DOFNOD)
C      DOFNOD  NUMBER OF DEGREES OF FREEDOM AT EACH NODE
C      TOTNOD  TOTAL NUMBER OF NODES IN MESH
C      NOUT    FORTRAN UNIT NUMBER
C      ITEST   ERROR CHECKING OPTION
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE PRTVAL(VAL, IVAL, NF, INF, JNF, DOFNOD,
C    *      TOTNOD, NOUT, ITEST)
C***********************************************************************
C
      INTEGER DOFNOD, ERRMES, I, IERROR, INF ,ITEST, IVAL, 
     *     J, JNF, K, N, NF, NODE, NOUT, TOTNOD
      DOUBLE PRECISION VAL, WORK, SRNAME
      DIMENSION VAL(IVAL), WORK(5), NF(INF,JNF), NODE(5)
      DATA SRNAME /8H PRTVAL /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (INF.LT.TOTNOD .OR. JNF.LT.DOFNOD)
     *                      IERROR = 2
                        IF (DOFNOD.LE.0 .OR. TOTNOD.LE.0)
     *                      IERROR = 1
                        ITEST = ERRMES(IERROR,ITEST,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 IF (DOFNOD.EQ.1) GO TO 1080
      IF (DOFNOD.EQ.2) WRITE (NOUT,9010)
      IF (DOFNOD.EQ.3) WRITE (NOUT,9020)
      IF (DOFNOD.EQ.4) WRITE (NOUT,9030)
      IF (DOFNOD.GE.5) WRITE (NOUT,9040)
      N = 0
      DO 1070 I=1,TOTNOD
      DO 1060 J=1,DOFNOD
      N = N + 1
      NODE(N) = I
      K = NF(I,J)
                        IF (ITEST.EQ.-1) GO TO 1020
                        IERROR = 0
                        IF (IVAL.LT.K) IERROR = 3
                        ITEST = ERRMES(IERROR,ITEST,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1020 WORK(N) = 0.0D0
      IF (K.EQ.0) GO TO 1040
      IF (K.GT.0) GO TO 1030
      K = IABS(K)
 1030 WORK(N) = VAL(K)
 1040 IF (((DOFNOD.EQ.2) .AND. (N.NE.4)) .OR. (DOFNOD.NE.2)) GO TO
     *     1050
      WRITE (NOUT,9080) NODE(1), WORK(1), WORK(2), NODE(3),
     *     WORK(3), WORK(4)
      N = 0
      GO TO 1060
 1050 IF ((DOFNOD.EQ.2) .OR. (N.NE.DOFNOD)) GO TO 1060
      WRITE (NOUT,9060) I, (WORK(K),K=1,DOFNOD)
      N = 0
 1060 CONTINUE
 1070 CONTINUE
      IF (N.EQ.0) RETURN
      IF (DOFNOD.EQ.2) WRITE (NOUT,9080) NODE(1), WORK(1), WORK(2)
      IF (DOFNOD.NE.2) WRITE (NOUT,9060) I, (WORK(K),K=1,N)
      RETURN
 1080 WRITE (NOUT,9050)
      N = 0
      DO 1110 I=1,TOTNOD
      N = N + 1
      NODE(N) = I
      K = NF(I,1)
      WORK(N) = 0.0D0
      IF (K.EQ.0) GO TO 1100
      IF (K.GT.0) GO TO 1090
      K = IABS(K)
 1090 WORK(N) = VAL(K)
 1100 IF (N.NE.4) GO TO 1110
      WRITE (NOUT,9070) (NODE(J),WORK(J),J=1,4)
      N = 0
 1110 CONTINUE
      IF (N.EQ.0) RETURN
      WRITE (NOUT,9070) (NODE(J),WORK(J),J=1,N)
      RETURN
 9010 FORMAT (/1H , 2(4HNODE, 6X, 2(5HVALUE, 9X))/1H )
 9020 FORMAT (/5H NODE, 6X, 3(5HVALUE, 9X)/1H )
 9030 FORMAT (/5H NODE, 6X, 4(5HVALUE, 9X)/1H )
 9040 FORMAT (/5H NODE, 6X, 5(5HVALUE, 9X)/1H )
 9050 FORMAT (/1H , 4(4HNODE, 5X, 5HVALUE, 5X)/1H )
 9060 FORMAT (1H , I4, 1X, 5(D12.5, 2X)/1H , 5X, 5(D12.5,
     *     2X)/1H , 5X, 5(D12.5, 2X))
 9070 FORMAT (1H , 4(I4, 1X, D12.5, 2X))
 9080 FORMAT (1H , 2(I4, 1X, 2(D12.5, 2X), 5X))
      END
