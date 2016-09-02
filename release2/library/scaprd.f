C***********************************************************************
C$SPLIT$SCAPRD$*********************************************************
C***********************************************************************
      SUBROUTINE SCAPRD(V, IV, W, IW, N, PRODCT, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      FORMS THE SCALAR PRODUCT OF TWO VECTORS V AND W, STORING
C      THE RESULT IN PRODCT
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    22 OCT 1980 (KR)
C
C ARGUMENTS IN
C      V       VECTOR OF LENGTH IV
C      IV      DIMENSION OF V (.GE.N)
C      W       VECTOR OF LENGTH IW
C      IW      DIMENSION OF W (.GE.N)
C      N       NUMBER OF ELEMENTS OF V (AND W) TO BE USED IN
C              FORMING SCALAR PRODUCT
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      PRODCT  CONTAINS: SIGMA I = 1 TO N OF (V(I)*W(I))
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE SCAPRD(V, IV, W, IW, N, PRODCT, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IERROR, ITEST, IV, IW, N
      DOUBLE PRECISION PRODCT, SRNAME, V, W
      DIMENSION V(IV), W(IW)
      DATA SRNAME /8H SCAPRD /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (N.GT.IW) IERROR = 3
                        IF (N.GT.IV) IERROR = 2
                        IF (N.LE.0) IERROR = 1
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (IERROR.NE.0) RETURN
 1010                   PRODCT = 0.0D0
                        DO 1020 I=1,N
                        PRODCT = PRODCT + V(I)*W(I)
 1020                   CONTINUE
                        RETURN
                        END
