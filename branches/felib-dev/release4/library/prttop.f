C
      SUBROUTINE PRTTOP(TOTELS,ELTOP,IELTOP,JELTOP,NOUT,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      PRTTOP prints element topologies in a standard format
C
C HISTORY
C
C      Copyright (C) 2003 : CCLRC, Rutherford Appleton Laboratory
C                           Chilton, Didcot, Oxfordshire OX11 0QX
C
C      Contact: Prof Chris Greenough
C      Email: c.greenough@rl.ac.uk
C      Tel: (44) 1235 445307
C
C      Release 1.1  29 Oct 1979 (CG)
C      Commented    14 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      TOTELS  total number of elements in the mesh
C      ELTOP   integer array of dimension (IELTOP, JELTOP)
C              containing element topologies, element type, and
C              number of nodes on the element
C      IELTOP  first dimension of ELTOP (.GE. TOTELS)
C      JELTOP  second dimension of ELTOP (.GE. NUMBER of nodes
C              on element + 2)
C      NOUT    fortran unit number
C      ITEST   error checking option
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE PRTTOP(TOTELS,ELTOP,IELTOP,JELTOP,NOUT,ITEST)
C***********************************************************************
C
      INTEGER ELTOP,ERRMES,I,IELTOP,IERROR,ITEST,J,JELTOP,JTEST,K,L,
     *        NOUT,TOTELS
      CHARACTER*6 SRNAME
      DIMENSION ELTOP(IELTOP,JELTOP)
      DATA SRNAME/'PRTTOP'/
C
C     Parameter checking
C
      JTEST = ITEST
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (IELTOP.LT.TOTELS) IERROR = 2
         IF (TOTELS.LE.0) IERROR = 1
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      WRITE (NOUT,FMT=9990)
      WRITE (NOUT,FMT=9980) TOTELS
      WRITE (NOUT,FMT=9970)
      DO 1000 I = 1,TOTELS
         L = ELTOP(I,2)
         K = L + 2
C
C     Range checking on K
C
         IF (JTEST.NE.-1) THEN
            IERROR = 0
            IF (JELTOP.LT.K) IERROR = 3
            ITEST = ERRMES(JTEST,IERROR,SRNAME)
            IF (ITEST.NE.0) RETURN
         END IF
C
         WRITE (NOUT,FMT=9960) I, (ELTOP(I,J),J=1,K)
 1000 CONTINUE
C
 9990 FORMAT (' ',////' **** ELEMENT TOPOLOGY ****',//' ')
 9980 FORMAT (' NUMBER OF ELEMENTS = ',I3)
 9970 FORMAT (/' ',2X,'ELEM',4X,'ELTYP',4X,'NODEL',4X,'NODES',/' ')
 9960 FORMAT (' ',2X,I4,4X,I3,7X,I3,4X,9 (I4,2X),/' ',27X,9 (I4,2X))
      END
