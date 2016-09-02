C
      SUBROUTINE PRTGEO(TOTNOD,DIMEN,COORD,ICOORD,JCOORD,NOUT,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      PRTGEO prints out element geometry in a standard format
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
C      TOTNOD  total number of nodes in the mesh
C      DIMEN   dimensionality of the geometric data
C      COORD   array of dimension (ICOORD, JCOORD) containing
C              global coordinates of the nodes
C      ICOORD  first dimension of COORD (.GE. TOTNOD)
C      JCOORD  second dimension of COORD (.GE. DIMEN)
C      NOUT    fortran unit number
C      ITEST   error checking option
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE PRTGEO(TOTNOD,DIMEN,COORD,ICOORD,JCOORD,NOUT,ITEST)
C***********************************************************************
C
      INTEGER TOTNOD,DIMEN,ICOORD,JCOORD,ITEST,NOUT,ERRMES,IERROR
      INTEGER I,J
      CHARACTER*6 SRNAME
      DOUBLE PRECISION COORD
      DIMENSION COORD(ICOORD,JCOORD)
      DATA SRNAME/'PRTGEO'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (ICOORD.LT.TOTNOD .OR. JCOORD.LT.DIMEN) IERROR = 2
         IF (TOTNOD.LE.0 .OR. DIMEN.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      WRITE (NOUT,FMT=9990)
      WRITE (NOUT,FMT=9980) TOTNOD
      WRITE (NOUT,FMT=9970) DIMEN
      IF (DIMEN.EQ.1) WRITE (NOUT,FMT=9960)
      IF (DIMEN.EQ.2) WRITE (NOUT,FMT=9950)
      IF (DIMEN.EQ.3) WRITE (NOUT,FMT=9940)
      IF (DIMEN.EQ.4) WRITE (NOUT,FMT=9930)
C
      DO 1000 I = 1,TOTNOD
         WRITE (NOUT,FMT=9920) I, (COORD(I,J),J=1,DIMEN)
 1000 CONTINUE
C
 9990 FORMAT (' ',////' **** NODAL GEOMETRY ****',//' ')
 9980 FORMAT (' NUMBER OF NODES = ',I3)
 9970 FORMAT (' DIMENSIONS      = ',I3)
 9960 FORMAT (/' ',2X,'NODE',9X,'X1',/' ')
 9950 FORMAT (/' ',2X,'NODE',9X,'X1',12X,'X2',/' ')
 9940 FORMAT (/' ',2X,'NODE',9X,'X1',12X,'X2',12X,'X3',/' ')
 9930 FORMAT (/' ',2X,'NODE',9X,'X1',12X,'X2',12X,'X3',12X,'X4',/' ')
 9920 FORMAT (' ',2X,I3,5X,4 (D12.4,2X))
      END
