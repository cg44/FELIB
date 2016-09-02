C***********************************************************************
C$SPLIT$UNFLO$*********************************************************
C***********************************************************************
      DOUBLE PRECISION FUNCTION UNFLO(X)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS THE VALUE OF THE SMALLEST POSITIVE REAL FLOATING-
C      POINT NUMBER EXACTLY REPRESENTABLE ON THE COMPUTER
C      *********************************************
C      **********MACHINE DEPENDENT ROUTINE**********
C      *********************************************
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    23 OCT 1980 (KR)
C
C ARGUMENTS IN
C      X       DUMMY ARGUMENT
C
C
C     DOUBLE PRECISION FUNCTION UNFLO(X)
C***********************************************************************
C
      DOUBLE PRECISION X, VALUE
C
C      FOR PRIME 400 X02ABF = 2.0D0**(-32386)
C+PRIME  UNFLO = 2.0D0**(-32386)
      INTEGER *2 L(4)
      EQUIVALENCE (L(1),VALUE)
      DATA L(1),L(2),L(3),L(4)/:040000,:000000,:000000,:100777/
C+IBM+GEC
C+    DATA VALUE /Z0010000000000000/
C
      UNFLO = VALUE
      RETURN
      END
