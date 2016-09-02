C***********************************************************************
      SUBROUTINE CSYSOL(KB, IKB, JKB, LOADS, ILOADS, N, HBAND,
     *     ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      SOLVES A SET OF COMPLEX SYMMETRIC BANDED EQUATIONS WITH A
C      SINGLE RIGHT HAND SIDE USING A SYMMETRIC DECOMPOSITION.  ONLY
C      THE LOWER BAND AND DIAGONAL ARE STORED IN A RECTANGULAR ARRAY KB.
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 3.0  29 OCT 1979 (CRIE)
C      COMMENTED    12 FEB 1985 (CG)
C
C ARGUMENTS IN
C      KB      ON ENTRY CONTAINS LOWER HALF OF COMPLEX SYMMETRIC
C              BAND MATRIX STORED AS A RECTANGULAR ARRAY
C      IKB     FIRST DIMENSION OF KB (.GE. N)
C      JKB     SECOND DIMENSION OF KB (.GE. HBAND)
C      LOADS   CONTAINS ELEMENTS OF RIGHT HAND SIDE
C      ILOADS  DIMENSION OF LOADS (.GE. N)
C      N       ORDER OF MATRIX KB
C      HBAND   SEMI-BANDWIDTH OF KB (INCLUDES DIAGONAL)
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      KB      ON EXIT, CONTAINS LOWER TRIANGULAR REDUCED
C      LOADS   MATRIX SOLUTION VECTOR
C
C ROUTINES CALLED
C      ERRMES, CSYRDN, CSYSUB
C
C     SUBROUTINE CSYSOL(KB, IKB, JKB, LOADS, ILOADS, N, HBAND,
C    *     ITEST)
C*************************************************************************
C
      INTEGER HBAND, IKB, ITEST, JKB, JTEST, IERROR,
     *     N, ERRMES, ILOADS
      DOUBLE PRECISION KB, LOADS
      DOUBLE PRECISION SRNAME
      DIMENSION KB(2,IKB,JKB), LOADS(2,ILOADS)
      DATA SRNAME /8H CSYSOL /
C
C     PARAMETER CHECKING
C
                IF(ITEST.EQ.-1) GO TO 999
                IERROR=0
                IF(ILOADS.LT.N) IERROR=3
                IF(IKB.LT.N.OR.JKB.LT.HBAND) IERROR=2
                IF(N.LE.0.OR.HBAND.LE.0) IERROR=1
                ITEST=ERRMES(ITEST,IERROR,SRNAME)
                IF(ITEST.NE.0) RETURN
C
C     MAIN BODY
C
 999  CONTINUE
      JTEST=1
      CALL CSYRDN(KB, IKB, JKB, N, HBAND, JTEST)
C
C     PARAMETER CHECKING
C
                IF(ITEST.EQ.-1) GO TO 998
                IERROR=JTEST
                IF(JTEST.EQ.3) IERROR = 4
                ITEST=ERRMES(ITEST, IERROR, SRNAME)
                IF(ITEST.NE.0) RETURN
C
 998  CONTINUE
      CALL CSYSUB(KB, IKB, JKB, LOADS, ILOADS, N, HBAND,
     *     JTEST)
C
C     PARAMETER CHECKING
C
                IF(ITEST.EQ.-1) GO TO 997
                IERROR=JTEST
                IF(JTEST.EQ.3) IERROR = 4
                ITEST=ERRMES(ITEST, IERROR, SRNAME)
                IF(ITEST.NE.0) RETURN
C
 997  CONTINUE
      RETURN
C
      END
C***********************************************************************
