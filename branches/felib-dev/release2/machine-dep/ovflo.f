C***********************************************************************
C$SPLIT$OVFLO$*********************************************************
C***********************************************************************
      DOUBLE PRECISION FUNCTION OVFLO(X)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS THE VALUE OF THE LARGEST POSITIVE REAL FLOATING-
C      POINT NUMBER REPRESENTABLE ON THE COMPUTER
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
C     DOUBLE PRECISION FUNCTION OVFLO(X)
C***********************************************************************
C
      DOUBLE PRECISION X, VALUE
C
C+PRIME OVFLO = 2.0D0 **(32576)
      INTEGER *2 L(4)
      EQUIVALENCE (L(1),VALUE)
      DATA L(1),L(2),L(3),L(4)/:077777,:177777,:177777,:077500/
C+IBM+GEC
C     DATA VALUE /Z7FFFFFFFFFFFFFFF/
C
      OVFLO = VALUE
      RETURN
      END
