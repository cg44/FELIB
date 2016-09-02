C
      SUBROUTINE SCAPRD(V,IV,W,IW,N,PRODCT,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      SCAPR forms the scalar product of two vectors V and W, storing
C      the result in PRODCT
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
C      Commented    22 Oct 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      V       vector of length IV
C      IV      dimension of V (.GE. N)
C      W       vector of length IW
C      IW      dimension of W (.GE. N)
C      N       number of elements of V (and W) to be used in
C              forming scalar product
C      ITEST   error checking option
C
C ARGUMENTS out
C      PRODCT  contains: sigma I = 1 to N of (V(I)*W(I))
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE SCAPRD(V,IV,W,IW,N,PRODCT,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,I,IERROR,ITEST,IV,IW,N
      CHARACTER*6 SRNAME
      DOUBLE PRECISION PRODCT,V,W
      DIMENSION V(IV),W(IW)
      DATA SRNAME/'SCAPRD'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (N.GT.IW) IERROR = 3
         IF (N.GT.IV) IERROR = 2
         IF (N.LE.0) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (IERROR.NE.0) RETURN
      END IF
C
C     Main body
C
      PRODCT = 0.0D0
      DO 1000 I = 1,N
         PRODCT = PRODCT + V(I)*W(I)
 1000 CONTINUE
C
      END
