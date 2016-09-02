C
      INTEGER FUNCTION MAXINT(N)
C-----------------------------------------------------------------------
C PURPOSE
C      MAXINT returns the value of the largest integer capable of
C      being stored by the computer
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
C      Release 1.1  29 Oct 1979 (CG)
C      Release 3.0  18 Jan 1984 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      N       dummy argument
C
C***********************************************************************
C
      INTEGER N
C
C Fortran 90/95 MAXINT=HUGE(N)
C
      MAXINT = 2147483647
C
      RETURN
      END
