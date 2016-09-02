C
      SUBROUTINE LINQUA(XI,ETA,GEOM,IGEOM,JGEOM,NODEL,SIDNUM,ALEN,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      LINQUA calculates a unit of length along the side of a
C      rectangular element (4, 8 or 12 noded)
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
C      Release 2.0   1 Feb 1981 (CG)
C      Commented    24 Jul 1981 (CG)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      XI      local coordinate of point at which length is
C              required
C      ETA     local coordinate of point at which length is
C              required
C      GEOM    local coordinate array containing the global
C              coordiantes of each node on an element in the
C              local order
C      IGEOM   first dimension of array GEOM (.GE. NODEL)
C      JGEOM   second dimension of array GEOM (.GE. 2)
C      NODEL   number of nodes on the element
C      SIDNUM  the side number of the side of the element to
C              be used in calculating the length (.LE. 3)
C      ITEST   error checking option
C
C ARGUMENTS out
C      ALEN    the unit of length at the specified point
C
C ROUTINES called
C      ERRMES  QUAM4  QUAM8  QUAM12  MATMUL
C
C     SUBROUTINE LINQUA(XI,ETA,GEOM,IGEOM,JGEOM,NODEL,SIDNUM,ALEN,ITEST)
C***********************************************************************
C
      INTEGER DIMEN,ERRMES,IDER,IERROR,IFUN,IGEOM,IJAC,ITEST,JDER,JGEOM,
     *        JJAC,JTEST,N,NODEL,SIDNUM
      CHARACTER*6 SRNAME
      DOUBLE PRECISION ALEN,DER,ETA,FUN,GEOM,JAC,XI
      DIMENSION DER(2,12),FUN(12),GEOM(IGEOM,JGEOM),JAC(2,2)
C
      EXTERNAL  ERRMES,QUAM4,QUAM8,QUAM12,MATMUL
C
      DATA DIMEN/2/,IDER/2/,IFUN/12/,IJAC/2/,JDER/12/,JJAC/2/,
     *     SRNAME/'LINQUA'/
C
C     Parameter checking
C
      JTEST = ITEST
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF ((NODEL.LE.0) .OR. (SIDNUM.LE.0)) IERROR = 1
         IF (IGEOM.LT.NODEL .OR. JGEOM.LT.2) IERROR = 2
         IF (SIDNUM.LT.1 .OR. SIDNUM.GT.4) IERROR = 3
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      N = NODEL/4
C
C     Range checking on N (should 1, 2 or 3)
C
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (N.LE.0 .OR. N.GE.4) IERROR = 4
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
      GO TO (1000,1010,1020) N
C
 1000 CONTINUE
      CALL QUAM4(FUN,IFUN,DER,IDER,JDER,XI,ETA,ITEST)
      GO TO 1030
 1010 CONTINUE
      CALL QUAM8(FUN,IFUN,DER,IDER,JDER,XI,ETA,ITEST)
      GO TO 1030
 1020 CONTINUE
      CALL QUAM12(FUN,IFUN,DER,IDER,JDER,XI,ETA,ITEST)
C
 1030 CONTINUE
      CALL MATMUL(DER,IDER,JDER,GEOM,IGEOM,JGEOM,JAC,IJAC,JJAC,DIMEN,
     *            NODEL,DIMEN,ITEST)
C
C     SELECT correct side of element
C
      GO TO (1040,1050,1040,1050) SIDNUM
 1040 CONTINUE
      ALEN = DSQRT(JAC(2,1)**2+JAC(2,2)**2)
      RETURN
 1050 CONTINUE
      ALEN = DSQRT(JAC(1,1)**2+JAC(1,2)**2)
C
      END
