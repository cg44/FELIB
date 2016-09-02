C***********************************************************************
C$SPLIT$FREDIF$*********************************************************
C***********************************************************************
      SUBROUTINE FREDIF(IELE, ELTOP, IELTOP, JELTOP, NF, INF, JNF,
     *     DOFNOD, FIRST, DIF, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CALCULATES THE MAXIMUM FREEDOM NUMBER DIFFERENCE FOR AN
C      ELEMENT
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      RECODED      12 FEB 1980 (CG)
C      COMMENTED    18 FEB 1980 (KR)
C
C ARGUMENTS IN
C      NELE    ELEMENT NUMBER
C      ELTOP   ELTOP(I,1) = ELEMENT TYPE OF ELEMENT I
C              ELTOP(I,2) = NUMBER OF NODES ON ELEMENT I
C              ELTOP(I,J+2), J=1(1)NUMBER OF NODES ON ELEMENT,
C              CONTAINS THE NODES ASSOCIATED WITH ELEMENT I
C      IELTOP  FIRST DIMENSION OF ARRAY ELTOP (.GE. NELE)
C      JELTOP  SECOND DIMENSION OF ELTOP (.GE. NUMBER OF NODES
C              ON ELEMENT)
C      NF      NF(I,J) CONTAINS THE FREEDOM NUMBERS ASSOCIATED
C              WITH NODE I
C      INF     FIRST DIMENSION OF NF (.GE. MAXIMUM NODE NUMBER
C              ON ELEMENT)
C      JNF     SECOND DIMENSION OF NF (.GE. DOFNOD)
C      DOFNOD  NUMBER OF DEGREES OF FREEDOM PER NODE ON THE
C              ELEMENT
C      FIRST   MUST BE SET TO .TRUE. FOR THE FIRST CALL TO
C              FREDIF AND .FALSE. FOR SUBSEQUENT CALLS
C      DIF     MUST BE ZERO FOR FIRST CALL TO FREDIF;
C              SUBSEQUENTLY CONTAINS THE MAXIMUM FREEDOM
C              DIFFERENCE PRIOR TO THE CURRENT CALL
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      FIRST   SET TO .FALSE.
C      DIF     MAXIMUM FREEDOM DIFFERENCE FOR ALL ELEMENTS UP
C              TO AND INCLUDING ELEMENT NELE
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE FREDIF(IELE, ELTOP, IELTOP, JELTOP, NF, INF, JNF,
C    *     DOFNOD, FIRST, DIF, ITEST)
C***********************************************************************
C
      INTEGER DIF, DOFNOD, ELTOP, ERRMES, I, IDEG, IELE, IELTOP,
     *     IERROR, INF, INOD, ITEST, J, JELTOP, JNF, MAX, MAXINT,
     *     MIN, NF, NODEL
      DOUBLE PRECISION SRNAME
      LOGICAL FIRST
      DIMENSION ELTOP(IELTOP,JELTOP), NF(INF,JNF)
      DATA SRNAME /8H FREDIF /
      IF (.NOT.FIRST) GO TO 1010
      FIRST = .FALSE.
      DIF = 0
 1010                   IF (ITEST.EQ.-1) GO TO 1020
                        IERROR = 0
                        IF (IELTOP.LT.IELE) IERROR = 3
                        IF (JNF.LT.DOFNOD) IERROR = 2
                        IF (IELE.LE.0 .OR. DOFNOD.LE.0) IERROR = 1
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1020  NODEL = ELTOP(IELE,2)
       IF (ITEST.EQ.-1) GO TO 1030
       IERROR = 0
       IF (JELTOP.LT.NODEL+2) IERROR = 4
       ITEST = ERRMES(ITEST,IERROR,SRNAME)
       IF (ITEST.NE.0) RETURN
 1030  MAX = 0
       MIN = MAXINT(MAX)
                        DO 1060 I=1,NODEL
                        INOD = ELTOP(IELE,I+2)
                        IF (ITEST.EQ.-1) GO TO 1040
                        IERROR = 0
                        IF (INF.LT.INOD) IERROR = 5
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1040                   DO 1050 J=1,DOFNOD
                        IDEG = NF(INOD,J)
                        IF (IDEG.EQ.0) GO TO 1050
                        MAX = MAX0(IDEG,MAX)
                        MIN = MIN0(IDEG,MIN)
 1050                   CONTINUE
 1060                   CONTINUE
                        DIF = MAX0(DIF,MAX-MIN)
                        RETURN
                        END
