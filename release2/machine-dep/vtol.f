C***********************************************************************
C$SPLIT$VTOL$*********************************************************
C***********************************************************************

      DOUBLE PRECISION FUNCTION VTOL(X)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS THE RATIO OF THE SMALLEST NUMBER WHICH DOES NOT
C      CAUSE UNDERFOW TO THE NUMBER VEPS WHICH IS DEFINED BY
C      1.DO+VEPS > 1.D0
C      ***************************************************
C      **********USES MACHINE DEPENDENT ROUTINES**********
C      **********  CODE IS MACHINE INDEPENDENT  **********
C      ***************************************************
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    23 OCT 1980 (KR)
C
C ARGUMENTS IN
C      X       DUMMY ARGUMENT
C
C ROUTINES CALLED
C      UNFLO   VEPS
C
C
C     DOUBLE PRECISION FUNCTION VTOL(X)
C***********************************************************************
C
      DOUBLE PRECISION X, VALUE
C+PRIME VTOL = 2.0D0**(-32342)
      INTEGER *2 L(4)
      EQUIVALENCE (L(1),VALUE)
      DATA L(1),L(2),L(3),L(4)/:040000,:000000,:000000,:101053/
C
      VTOL = VALUE
      RETURN
      END
