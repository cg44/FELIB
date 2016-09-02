C***********************************************************************
      SUBROUTINE CVCNUL(VEC, IVEC, N, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      SETS THE FIRST N ELEMENTS OF A COMPLEX VECTOR TO ZERO
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 3.0  29 OCT 1984 (CRIE)
C      COMMENTED    23 OCT 1985 (CG)
C
C ARGUMENTS IN
C      IV      LENGTH OF VECTOR V (.GE.N)
C      N       NUMBER OF ELEMENTS TO BE SET TO ZERO
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      V       VECTOR OF LENGTH IV.  V(1,I)=0.0D0 AND V(2,I)=0.0D0
C              FOR I=1(1)N
C
C ROUTINES CALLED
C      ERRMES
C
C     SUBROUTINE CVCNUL(VEC, IVEC, N, ITEST)
C***********************************************************************
C
      INTEGER ERRMES, I, IERROR, ITEST, IVEC, N
      DOUBLE PRECISION SRNAME, VEC
      DIMENSION VEC(2,IVEC)
      DATA SRNAME /8H CVCNUL /
C
C     PARAMETER CHECKING
C
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (N.GT.IVEC) IERROR = 2
                        IF (N.LE.0) IERROR = 1
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
C
C     MAIN LOOPS
C
                        IF (IERROR.NE.0) RETURN
 1010                   DO 1020 I=1,N
                        VEC(1,I) = 0.0D0
                        VEC(2,I) = 0.0D0
 1020                   CONTINUE
                        RETURN
C
                        END
C***********************************************************************
