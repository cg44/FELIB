C***********************************************************************
      DOUBLE PRECISION FUNCTION VEPS(X)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS THE VALUE OF THE SMALLEST NUMBER VEPS SUCH THAT
C      1.D0+VEPS > 1.DO ON THIS COMPUTER
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
C     DOUBLE PRECISION FUNCTION VEPS(X)
C***********************************************************************
C
      DOUBLE PRECISION X
C
C FORTRAN 90/95 VEPS=EPSILON(X)
C
      VEPS = 2.220446049250313E-16
C
      RETURN
      END
C***********************************************************************
