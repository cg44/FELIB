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
C
C ARGUMENTS IN
C      N       DUMMY ARGUMENT
C
C
C     INTEGER FUNCTION MAXINT(N)
C***********************************************************************
C
      INTEGER N
C+UNIX
C     INTEGER INMAX
C     EXTERNAL INMAX
C
C+PRIME+IBM+GEC
      MAXINT = 2147483647
C=UNIX
C
C      MAXINT=INMAX(N)
      RETURN
      END
C***********************************************************************
