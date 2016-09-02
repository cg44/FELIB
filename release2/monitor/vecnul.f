C***********************************************************************
C$SPLIT$VECNUL$*********************************************************
C***********************************************************************
      SUBROUTINE VECNUL(VEC, IVEC, N, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      SETS THE FIRST N ELEMENTS OF A VECTOR TO ZERO
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (IMS) --- SERC COPYRIGHT
C      COMMENTED    23 OCT 1980 (KR)
C
C ARGUMENTS IN
C      IV      LENGTH OF VECTOR V (.GE.N)
C      N       NUMBER OF ELEMENTS TO BE SET TO ZERO
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      V       VECTOR OF LENGTH IV.  V(I)=0.0D0 FOR I=1(1)N
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE VECNUL(VEC, IVEC, N, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IERROR, ITEST, IVEC, N
      DOUBLE PRECISION SRNAME, VEC
      DIMENSION VEC(IVEC)
      DATA SRNAME /8H VECNUL /
C
C     CHECK ON MONITORING
C
      IF(ITEST.EQ.-1) CALL FECHK
C
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (N.GT.IVEC) IERROR = 2
                        IF (N.LE.0) IERROR = 1
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (IERROR.NE.0) RETURN
 1010                   DO 1020 I=1,N
                        VEC(I) = 0.0D0
 1020                   CONTINUE
                        RETURN
                        END
