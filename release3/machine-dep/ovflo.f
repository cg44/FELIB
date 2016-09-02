C***********************************************************************
      DOUBLE PRECISION FUNCTION OVFLO(X)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS THE VALUE OF THE LARGEST POSITIVE REAL FLOATING-
C      POINT NUMBER REPRESENTABLE ON THE COMPUTER
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
C     DOUBLE PRECISION FUNCTION OVFLO(X)
C***********************************************************************
C
      DOUBLE PRECISION X, VALUE, DFLMAX
      EXTERNAL DFLMAX
C
C     DATA VALUE /1.D75/
C
C+PRIME OVFLO = 2.0D0 **(32576)
C      INTEGER *2 L(4)
C      EQUIVALENCE (L(1),VALUE)
C      DATA L(1),L(2),L(3),L(4)/:077777,:177777,:177777,:077500/
C+IBM+GEC
C     DATA VALUE /Z7FFFFFFFFFFFFFFF/
C
C     OVFLO = VALUE
C+UNIX
      OVFLO=DFLMAX(VALUE)
      RETURN
      END
C***********************************************************************
