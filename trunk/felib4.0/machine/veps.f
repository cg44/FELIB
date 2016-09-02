C
      DOUBLE PRECISION FUNCTION VEPS(X)
C-----------------------------------------------------------------------
C PURPOSE
C      VEPS returns the value of the smallest number VEPS such that
C      1.D0+VEPS > 1.D0 on this computer
C
C      *********************************************
C      **********MACHINE DEPENDENT ROUTINE**********
C      *********************************************
C
C HISTORY
C
C      Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
C                           Chilton, Didcot, Oxfordshire OX11 0QX
C
C      Release 3.0  18 Jan 1984 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      X       dummy argument
C
C***********************************************************************
C
      DOUBLE PRECISION X
C
C Fortran 90/95 VEPS=EPSILON(X)
C
      VEPS = 2.220446049250313E-16
C
      RETURN
      END
