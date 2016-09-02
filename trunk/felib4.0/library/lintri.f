C
      SUBROUTINE LINTRI(XI,ETA,GEOM,IGEOM,JGEOM,NODEL,SIDNUM,ALEN,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      LINTRI calculates a unit of length along the side of a
C      triangular element (3, 6 or 10 noded)
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
C              be used in calculating the length
C      ITEST   error checking option
C
C ARGUMENTS out
C      ALEN    the unit of length at the specified point
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE LINTRI(XI,ETA,GEOM,IGEOM,JGEOM,NODEL,SIDNUM,ALEN,ITEST)
C***********************************************************************
C
      INTEGER ERRMES,IERROR,IGEOM,ITEST,JGEOM,JTEST,N,NODEL,SIDNUM
      CHARACTER*6 SRNAME
      DOUBLE PRECISION ALEN,DXDL,DYDL,ETA,GEOM,L1,L2,L3,T1,T2,T3,T4,XI
      DIMENSION GEOM(IGEOM,JGEOM)
      DATA SRNAME/'LINTRI'/
C
C     Parameter checking
C
      JTEST = ITEST
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF ((NODEL.LE.0) .OR. (SIDNUM.LE.0)) IERROR = 1
         IF (IGEOM.LT.NODEL .OR. JGEOM.LT.2) IERROR = 2
         IF (SIDNUM.LT.1 .OR. SIDNUM.GT.3) IERROR = 3
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      N = NODEL/3
C
C     Range checking on N (should be 1, 2 or 3)
C
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (N.LE.0 .OR. N.GE.4) IERROR = 4
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
      GO TO (1000,1040,1080) N
 1000 CONTINUE
      GO TO (1010,1020,1030) SIDNUM
C
C     Code for 3-node triangle
C
 1010 CONTINUE
      DXDL = GEOM(1,1) - GEOM(2,1)
      DYDL = GEOM(1,2) - GEOM(2,2)
      GO TO 1120
 1020 CONTINUE
      DXDL = GEOM(2,1) - GEOM(3,1)
      DYDL = GEOM(2,2) - GEOM(3,2)
      GO TO 1120
 1030 CONTINUE
      DXDL = GEOM(3,1) - GEOM(1,1)
      DYDL = GEOM(3,2) - GEOM(1,2)
      GO TO 1120
C
C     Code for 6-node triangle
C
 1040 CONTINUE
      GO TO (1050,1060,1070) SIDNUM
 1050 CONTINUE
      L1 = 1.D0/3.D0* (1.D0+2.D0*XI)
      DXDL = (4.D0*L1-1.D0)*GEOM(1,1) + 4.D0* (1.D0-2.D0*L1)*GEOM(2,1) +
     *       (4.D0*L1-3.D0)*GEOM(3,1)
      DYDL = (4.D0*L1-1.D0)*GEOM(1,2) + 4.D0* (1.D0-2.D0*L1)*GEOM(2,2) +
     *       (4.D0*L1-3.D0)*GEOM(3,2)
      GO TO 1120
 1060 CONTINUE
      L2 = 1.D0/3.D0* (1.D0-XI-DSQRT(3.D0)*ETA)
      DXDL = (4.D0*L2-1.D0)*GEOM(3,1) + 4.D0* (1.D0-2.D0*L2)*GEOM(4,1) +
     *       (4.D0*L2-3.D0)*GEOM(5,1)
      DYDL = (4.D0*L2-1.D0)*GEOM(3,2) + 4.D0* (1.D0-2.D0*L2)*GEOM(4,2) +
     *       (4.D0*L2-3.D0)*GEOM(5,2)
      GO TO 1120
 1070 CONTINUE
      L3 = 1.D0/3.D0* (1.D0-XI+DSQRT(3.D0)*ETA)
      DXDL = (4.D0*L3-3.D0)*GEOM(1,1) + (4.D0*L3-1.D0)*GEOM(5,1) +
     *       4.D0* (1.D0-2.D0*L3)*GEOM(6,1)
      DYDL = (4.D0*L3-3.D0)*GEOM(1,2) + (4.D0*L3-1.D0)*GEOM(5,2) +
     *       4.D0* (1.D0-2.D0*L3)*GEOM(6,2)
      GO TO 1120
C
C     Code for 10-node triangle
C
 1080 CONTINUE
      GO TO (1090,1100,1110) SIDNUM
 1090 CONTINUE
      L1 = 1.D0/3.D0* (1.D0+2.D0*XI)
      L2 = 1.D0 - L1
      T1 = 1.D0/2.D0* (27.D0*L1*L1-14.D0*L2+2.D0)
      T2 = 9.D0/2.D0* (3.D0*L1*L1+6.D0*L1*L2-L1-L2)
      T3 = 9.D0/2.D0* (3.D0*L2*L2+6.D0*L1*L2-L1-L2)
      T4 = 1.D0/2.D0* (27*L2*L2-14.D0*L2+2.D0)
      DXDL = T1*GEOM(1,1) + T2*GEOM(2,1) + T3*GEOM(4,1) + T4*GEOM(4,1)
      DYDL = T1*GEOM(1,2) + T2*GEOM(2,2) + T3*GEOM(4,2) + T4*GEOM(4,2)
      GO TO 1120
 1100 CONTINUE
      L2 = 1.D0/3.D0* (1.D0-XI-DSQRT(3.D0)*ETA)
      L3 = 1.D0 - L2
      T1 = 1.D0/2.D0* (27.D0*L2*L2-14.D0*L2+2.D0)
      T2 = 9.D0/2.D0* (3.D0*L2*L2+6.D0*L2*L3-L2-L3)
      T3 = 9.D0/2.D0* (3.D0*L3*L3+6.D0*L2*L3-L2-L3)
      T4 = 1.D0/2.D0* (27.D0*L3*L3-14.D0*L3+2.D0)
      DXDL = T1*GEOM(4,1) + T2*GEOM(5,1) + T3*GEOM(6,1) + T4*GEOM(7,1)
      DXDL = T1*GEOM(4,2) + T2*GEOM(5,2) + T3*GEOM(6,2) + T4*GEOM(7,2)
      GO TO 1120
 1110 CONTINUE
      L3 = 1.D0/3.D0* (1.D0-XI+DSQRT(3.D0)*ETA)
      L1 = 1.D0 - L3
      T1 = 1.D0/2.D0* (27.D0*L1*L1-14.D0*L1+2.D0)
      T2 = 1.D0/2.D0* (27.D0*L3*L3-14.D0*L3+2.D0)
      T3 = 9.D0/2.D0* (3.D0*L3*L3+6.D0*L1*L3-L1-L3)
      T4 = 9.D0/2.D0* (3.D0*L1*L1+6.D0*L1*L3-L1-L3)
      DXDL = T1*GEOM(1,1) + T2*GEOM(7,1) + T3*GEOM(8,1) + T4*GEOM(9,1)
      DXDL = T1*GEOM(1,2) + T2*GEOM(7,2) + T3*GEOM(8,2) + T4*GEOM(9,2)
C
C     Do the calculation
C
 1120 CONTINUE
      ALEN = DSQRT(DXDL**2+DYDL**2)
C
      END
