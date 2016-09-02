C***********************************************************************
C$SPLIT$SHAPFN$*********************************************************
C***********************************************************************
      SUBROUTINE SHAPFN(N, IN, JN, FUN, IFUN, NODEL, NDE, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CONSTRUCTS THE SHAPE FUNCTION MATRIX N FOR A SYSTEM OF COUPLED
C      DIFFERENTIAL EQUATIONS.  THIS ROUTINE IS MOST FREQUENTLY USED
C      IN FORMING THE CONSISTENT ELEMENT AND SYSTEM MASS MATRICES
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    22 OCT 1980 (KR)
C
C ARGUMENTS IN
C      IN      FIRST DIMENSION OF ARRAY N (.GE.NDE)
C      JN      SECOND DIMENSION OF N (.GE.NODEL*NDE)
C      FUN     VECTOR OF LENGTH IFUN.  CONTAINS VALUES OF THE
C              SHAPE FUNCTIONS AT THE POINT WHERE SHAPE
C              FUNCTION MATRIX REQUIRED
C      IFUN    LENGTH OF VECTOR FUN (.GE.NODEL)
C      NODEL   NUMBER OF SHAPE FUNCTIONS TO BE USED IN
C              CONSTRUCTING SHAPE FUNCTION MATRIX N (USUALLY
C              NUMBER OF NODES IN ELEMENT UNDER CONSIDERATION)
C      NDE     NUMBER OF COUPLED EQUATIONS BEING CONSIDERED
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      N       ARRAY OF DIMENSION (IN,JN).  CONTAINS THE
C              VALUES OF THE SHAPE FUNCTION MATRIX
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE SHAPFN(N, IN, JN, FUN, IFUN, NODEL, NDE ,ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IERROR, IFUN, IJ, IN, ITEST, J, JN, NDE,
     *     NODEL
      DOUBLE PRECISION FUN, N, SRNAME
      DIMENSION FUN(IFUN), N(IN,JN)
      DATA SRNAME /8H SHAPFN /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (NODEL.LE.0 .OR. NDE.LE.0) IERROR = 1
                        IF (IN.LT.NDE .OR. JN.LT.NODEL*NDE)
     *                      IERROR = 2
                        IF (IFUN.LT.NODEL) IERROR = 3
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 CALL MATNUL(N, IN, JN, NDE, NODEL*NDE, ITEST)
      DO 1030 J=1,NODEL
      DO 1020 I=1,NDE
      IJ = (J-1)*NDE + I
      N(I,IJ) = FUN(J)
 1020 CONTINUE
 1030 CONTINUE
      RETURN
      END