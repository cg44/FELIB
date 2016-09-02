C***********************************************************************
C$SPLIT$MAXINT$*********************************************************
C***********************************************************************
      INTEGER FUNCTION MAXINT(N)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS THE VALUE OF THE LARGEST INTEGER CAPABLE OF
C      BEING STORED BY THE COMPUTER
C      *********************************************
C      **********MACHINE DEPENDENT ROUTINE**********
C      *********************************************
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
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
C
C+PRIME+IBM+GEC
      MAXINT = 2147483647
      RETURN
      END
