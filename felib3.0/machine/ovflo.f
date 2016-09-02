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
C      LINUX VALUES  1 OCT 2003 (CG)
C
C ARGUMENTS IN
C      X       DUMMY ARGUMENT
C
C
C     DOUBLE PRECISION FUNCTION OVFLO(X)
C***********************************************************************
C
      DOUBLE PRECISION X, VALUE, DFLMAX
C
C FORTRAN 90/95 OVFLO=HUGE(X)
C
C     OVFLO=1.797693134862316E+308
      OVFLO=1E+38
C
      RETURN
      END
C***********************************************************************
