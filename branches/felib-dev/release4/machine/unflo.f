C
      DOUBLE PRECISION FUNCTION UNFLO(X)
C-----------------------------------------------------------------------
C PURPOSE
C      UNFLO returns the value of the smallest positive real floating-
C      point number exactly representable on the computer
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
C FORTRAN 90/95 UNFLO=TINY(X)
C
      UNFLO=2.22507385850721D-308
C
      RETURN
      END
