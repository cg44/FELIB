C***********************************************************************
      DOUBLE PRECISION FUNCTION UNFLO(X)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS THE VALUE OF THE SMALLEST POSITIVE REAL FLOATING-
C      POINT NUMBER EXACTLY REPRESENTABLE ON THE COMPUTER
C
C      *********************************************
C      **********MACHINE DEPENDENT ROUTINE**********
C      *********************************************
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 1.1  29 OCT 1979 (CG)
C      COMMENTED    23 OCT 1980 (KR)
C
C ARGUMENTS IN
C      X       DUMMY ARGUMENT
C
C
C     DOUBLE PRECISION FUNCTION UNFLO(X)
C***********************************************************************
C+UNIX
      DOUBLE PRECISION DFLMIN
      EXTERNAL DFLMIN
C
      DOUBLE PRECISION X, VALUE
C
      DATA VALUE /1.D-70/
C
C      FOR PRIME 400 X02ABF = 2.0D0**(-32386)
C+PRIME  UNFLO = 2.0D0**(-32386)
C      INTEGER *2 L(4)
C      EQUIVALENCE (L(1),VALUE)
C      DATA L(1),L(2),L(3),L(4)/:040000,:000000,:000000,:100777/
C+IBM+GEC
C+    DATA VALUE /Z0010000000000000/
C
      UNFLO = VALUE
C     UNFLO=DFLMIN(VALUE)
      RETURN
      END
C***********************************************************************
