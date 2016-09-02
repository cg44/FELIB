C
      SUBROUTINE SURBRK(XI,ETA,ZETA,GEOM,IGEOM,JGEOM,NODEL,FACNUM,UAREA,
     *                  ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      SURBRK calculates a unit of area on the face of a
C      brick element (8, 20 or 32 noded)
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
C      Release 2.0  1  Feb 1981 (CG)
C      Commented    24 Jul 1981 (CG)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      XI      local coordinate of point at which area is
C              required
C      ETA     local coordinate of point at which area is
C              required
C      ZETA    local coordinate of point at which area is
C              required
C      GEOM    local coordinate array containing the global
C              coordiantes of each node on an element in the
C              local order
C      IGEOM   first dimension of array GEOM (.GE. NODEL)
C      JGEOM   second dimension of array GEOM (.GE. 3)
C      NODEL   number of nodes on the element
C      FACNUM  the face number of the face of the element to
C              be used in calculating the area (.LE. 6)
C      ITEST   error checking option
C
C ARGUMENTS out
C      UAREA    the unit of area at the specified point
C
C ROUTINES called
C      ERRMES  BRK8  BRK20  BRK32  MATMUL
C
C     SUBROUTINE SURBRK(XI,ETA,ZETA,GEOM,IGEOM,JGEOM,NODEL,FACNUM,UAREA,
C    *                  ITEST)
C***********************************************************************
C
      INTEGER DIMEN,ERRMES,FACNUM,I,IDER,IERROR,IFUN,IGEOM,IJAC,ITEST,J,
     *        JDER,JGEOM,JJAC,JTEST,N,NODEL
      CHARACTER*6 SRNAME
      DOUBLE PRECISION DER,ETA,FUN,G1,G2,G3,GEOM,JAC,UAREA,XI,ZETA
      DIMENSION DER(3,32),FUN(32),GEOM(IGEOM,JGEOM),JAC(3,3)
      DATA DIMEN/3/,IDER/3/,IFUN/32/,IJAC/3/,JDER/32/,JJAC/3/
      DATA SRNAME/'SURBRK'/
C
C     Parameter checking
C
      JTEST = ITEST
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF ((NODEL.LE.0) .OR. (FACNUM.LE.0)) IERROR = 1
         IF (IGEOM.LT.NODEL .OR. JGEOM.LT.3) IERROR = 2
         IF (FACNUM.LT.1 .OR. FACNUM.GT.6) IERROR = 3
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      N = NODEL/8
C
C     Range checking on N
C
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (N.LE.0 .OR. N.GE.4) IERROR = 4
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Get shape function derivatives
C
      GO TO (1020,1000,1010) N
      IF (JTEST.EQ.-1) THEN
         GO TO 1020
      ELSE
         IERROR = 4
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         RETURN
      END IF
 1000 CONTINUE
      CALL BRK20(FUN,IFUN,DER,IDER,JDER,XI,ETA,ZETA,ITEST)
      GO TO 1030
 1010 CONTINUE
      CALL BRK32(FUN,IFUN,DER,IDER,JDER,XI,ETA,ZETA,ITEST)
      GO TO 1030
C
 1020 CONTINUE
      CALL BRK8(FUN,IFUN,DER,IDER,JDER,XI,ETA,ZETA,ITEST)
C
C     Calculate global derivatives
C
 1030 CONTINUE
      CALL MATMUL(DER,IDER,JDER,GEOM,IGEOM,JGEOM,JAC,IJAC,JJAC,DIMEN,
     *            NODEL,DIMEN,ITEST)
C
      GO TO (1040,1040,1050,1060,1050,1060) FACNUM
      GO TO 1050
C
C     ZETA = constant
C
 1040 CONTINUE
      I = 1
      J = 2
      GO TO 1070
C
C     XI = constant
C
 1050 CONTINUE
      I = 2
      J = 3
      GO TO 1070
C
C     ETA = constant
C
 1060 CONTINUE
      I = 1
      J = 3
C
 1070 CONTINUE
      G1 = JAC(I,2)*JAC(J,3) - JAC(I,3)*JAC(J,2)
      G2 = JAC(I,3)*JAC(J,1) - JAC(I,1)*JAC(J,3)
      G3 = JAC(I,1)*JAC(J,2) - JAC(I,2)*JAC(J,1)
      UAREA = DSQRT(G1**2+G2**2+G3**2)
C
      END
