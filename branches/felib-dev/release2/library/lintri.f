C***********************************************************************
C$SPLIT$LINTRI$*********************************************************
C***********************************************************************
      SUBROUTINE LINTRI(XI, ETA, GEOM, IGEOM, JGEOM, NODEL,
     *     SIDNUM, ALEN, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CALCULATES A UNIT OF LENGTH ALONG THE SIDE OF A
C      TRIANGULAR ELEMENT (3, 6 OR 10 NODED)
C
C HISTORY
C      RELEASE 2.0  1  FEB 1981 (CG) --- SERC COPYRIGHT
C      COMMENTED    24 JULY 1981 (CG)
C
C ARGUMENTS IN
C      XI      LOCAL COORDINATE OF POINT AT WHICH LENGTH IS
C              REQUIRED
C      ETA     LOCAL COORDINATE OF POINT AT WHICH LENGTH IS
C              REQUIRED
C      GEOM    LOCAL COORDINATE ARRAY CONTAINING THE GLOBAL
C              COORDIANTES OF EACH NODE ON AN ELEMENT IN THE
C              LOCAL ORDER
C      IGEOM   FIRST DIMENSION OF ARRAY GEOM (.GE.NODEL)
C      JGEOM   SECOND DIMENSION OF ARRAY GEOM (.GE.2)
C      NODEL   NUMBER OF NODES ON THE ELEMENT
C      SIDNUM  THE SIDE NUMBER OF THE SIDE OF THE ELEMENT TO
C              BE USED IN CALCULATING THE LENGTH
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      ALEN    THE UNIT OF LENGTH AT THE SPECIFIED POINT
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE LINTRI(XI,ETA,GEOM,IGEOM,JGEOM,NODEL,SIDNUM,ALEN,ITEST)
C***********************************************************************
C
      INTEGER IGEOM, ITEST, JGEOM, N, NODEL, SIDNUM, IERROR ,ERRMES
      DOUBLE PRECISION ALEN, DXDL, DYDL, ETA, GEOM, 
     *     L1, L2, L3, SRNAME, T1, T2, T3, T4, XI
      DIMENSION GEOM(IGEOM,JGEOM)
      DATA SRNAME /8H LINTRI /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF ((NODEL.LE.0) .OR. (SIDNUM.LE.0))
     *                      IERROR = 1
                        IF (IGEOM.LT.NODEL .OR. JGEOM.LT.2)
     *                      IERROR = 2
                        IF (SIDNUM.LT.1 .OR. SIDNUM.GT.3)
     *                      IERROR = 3
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 N = NODEL/3
      IF(ITEST.EQ.-1) GOTO 1015
      IERROR=0
      IF(N.LE.0.OR.N.GE.4) IERROR=4
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF(ITEST.NE.0) RETURN
 1015 GO TO (1020, 1060, 1100), N
C+++++
C     CODE FOR 3-NODE TRIANGLE
C
 1020 GO TO (1030, 1040, 1050), SIDNUM
 1030 DXDL = GEOM(1,1) - GEOM(2,1)
      DYDL = GEOM(1,2) - GEOM(2,2)
      GO TO 1140
 1040 DXDL = GEOM(2,1) - GEOM(3,1)
      DYDL = GEOM(2,2) - GEOM(3,2)
      GO TO 1140
 1050 DXDL = GEOM(3,1) - GEOM(1,1)
C+++++
C      CODE FOR 6-NODE TRIANGLE
C
      DYDL = GEOM(3,2) - GEOM(1,2)
      GO TO 1140
 1060 GO TO (1070, 1080, 1090), SIDNUM
 1070 L1 = 1.D0/3.D0*(1.D0+2.D0*XI)
      DXDL = (4.D0*L1-1.D0)*GEOM(1,1) + 4.D0*(1.D0-2.D0*L1)*GEOM(2,
     *     1) + (4.D0*L1-3.D0)*GEOM(3,1)
      DYDL = (4.D0*L1-1.D0)*GEOM(1,2) + 4.D0*(1.D0-2.D0*L1)*GEOM(2,
     *     2) + (4.D0*L1-3.D0)*GEOM(3,2)
      GO TO 1140
 1080 L2 = 1.D0/3.D0*(1.D0-XI-DSQRT(3.D0)*ETA)
      DXDL = (4.D0*L2-1.D0)*GEOM(3,1) + 4.D0*(1.D0-2.D0*L2)*GEOM(4,
     *     1) + (4.D0*L2-3.D0)*GEOM(5,1)
      DYDL = (4.D0*L2-1.D0)*GEOM(3,2) + 4.D0*(1.D0-2.D0*L2)*GEOM(4,
     *     2) + (4.D0*L2-3.D0)*GEOM(5,2)
      GO TO 1140
 1090 L3 = 1.D0/3.D0*(1.D0-XI+DSQRT(3.D0)*ETA)
C+++++
C     CODE FOR 10-NODE TRIANGLE
C
      DXDL = (4.D0*L3-3.D0)*GEOM(1,1) + (4.D0*L3-1.D0)*GEOM(5,1) +
     *     4.D0*(1.D0-2.D0*L3)*GEOM(6,1)
      DYDL = (4.D0*L3-3.D0)*GEOM(1,2) + (4.D0*L3-1.D0)*GEOM(5,2) +
     *     4.D0*(1.D0-2.D0*L3)*GEOM(6,2)
      GO TO 1140
 1100 GO TO (1110, 1120, 1130), SIDNUM
 1110 L1 = 1.D0/3.D0*(1.D0+2.D0*XI)
      L2 = 1.D0 - L1
      T1 = 1.D0/2.D0*(27.D0*L1*L1-14.D0*L2+2.D0)
      T2 = 9.D0/2.D0*(3.D0*L1*L1+6.D0*L1*L2-L1-L2)
      T3 = 9.D0/2.D0*(3.D0*L2*L2+6.D0*L1*L2-L1-L2)
      T4 = 1.D0/2.D0*(27*L2*L2-14.D0*L2+2.D0)
      DXDL = T1*GEOM(1,1) + T2*GEOM(2,1) + T3*GEOM(4,1) +
     *     T4*GEOM(4,1)
      DYDL = T1*GEOM(1,2) + T2*GEOM(2,2) + T3*GEOM(4,2) +
     *     T4*GEOM(4,2)
      GO TO 1140
 1120 L2 = 1.D0/3.D0*(1.D0-XI-DSQRT(3.D0)*ETA)
      L3 = 1.D0 - L2
      T1 = 1.D0/2.D0*(27.D0*L2*L2-14.D0*L2+2.D0)
      T2 = 9.D0/2.D0*(3.D0*L2*L2+6.D0*L2*L3-L2-L3)
      T3 = 9.D0/2.D0*(3.D0*L3*L3+6.D0*L2*L3-L2-L3)
      T4 = 1.D0/2.D0*(27.D0*L3*L3-14.D0*L3+2.D0)
      DXDL = T1*GEOM(4,1) + T2*GEOM(5,1) + T3*GEOM(6,1) +
     *     T4*GEOM(7,1)
      DXDL = T1*GEOM(4,2) + T2*GEOM(5,2) + T3*GEOM(6,2) +
     *     T4*GEOM(7,2)
      GO TO 1140
 1130 L3 = 1.D0/3.D0*(1.D0-XI+DSQRT(3.D0)*ETA)
      L1 = 1.D0 - L3
      T1 = 1.D0/2.D0*(27.D0*L1*L1-14.D0*L1+2.D0)
      T2 = 1.D0/2.D0*(27.D0*L3*L3-14.D0*L3+2.D0)
      T3 = 9.D0/2.D0*(3.D0*L3*L3+6.D0*L1*L3-L1-L3)
      T4 = 9.D0/2.D0*(3.D0*L1*L1+6.D0*L1*L3-L1-L3)
C+++++
C     DO THE CALCULATION
C
      DXDL = T1*GEOM(1,1) + T2*GEOM(7,1) + T3*GEOM(8,1) +
     *     T4*GEOM(9,1)
      DXDL = T1*GEOM(1,2) + T2*GEOM(7,2) + T3*GEOM(8,2) +
     *     T4*GEOM(9,2)
 1140 ALEN = DSQRT(DXDL**2+DYDL**2)
      RETURN
      END
