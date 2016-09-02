C
      DOUBLE PRECISION FUNCTION OVFLO(X)
C-----------------------------------------------------------------------
C PURPOSE
C      OVFLO returns the value of the largest positive real floating-
C      point number representable on the computer
C
C      *********************************************
C      *********MACHINE DEPENDENT ROUTINE***********
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
C Fortran 90/95 OVFLO=HUGE(X)
C
      OVFLO = HUGE(X)
C
C      OVFLO=1.797693134862316E+308
C
      RETURN
      END
