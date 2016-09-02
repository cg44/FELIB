C***********************************************************************
      INTEGER FUNCTION MAXINT(N)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS THE VALUE OF THE LARGEST INTEGER CAPABLE OF
C      BEING STORED BY THE COMPUTER
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
C      N       DUMMY ARGUMENT
C
C
C     INTEGER FUNCTION MAXINT(N)
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
C***********************************************************************
