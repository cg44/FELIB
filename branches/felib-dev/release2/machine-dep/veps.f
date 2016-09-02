C***********************************************************************
C$SPLIT$VEPS$*********************************************************
C***********************************************************************
      DOUBLE PRECISION FUNCTION VEPS(X)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS THE VALUE OF THE SMALLEST NUMBER VEPS SUCH THAT
C      1.D0+VEPS > 1.DO ON THIS COMPUTER
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
C     DOUBLE PRECISION FUNCTION VEPS(X)
C***********************************************************************
C
      DOUBLE PRECISION X, VALUE
C+PRIME EPS = 2.0D0**(-44)
      INTEGER*2  L(4)
      EQUIVALENCE (VALUE,L(1))
      DATA L(1),L(2),L(3),L(4)/:040000,:000000,:000000,:000125/
C+GEC+IBM
C+GEC DATA VALUE /Z3410000000000000/
C
      VEPS = VALUE
      RETURN
      END
