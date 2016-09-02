C***********************************************************************
C$SPLIT$SURBRK$*********************************************************
C***********************************************************************
      SUBROUTINE SURBRK(XI, ETA, ZETA, GEOM, IGEOM, JGEOM, NODEL,
     *     FACNUM, UAREA, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CALCULATES A UNIT OF AREA ON THE FACE OF A
C      BRICK ELEMENT (8, 20 OR 32 NODED)
C
C HISTORY
C      RELEASE 2.0  1  FEB 1981 (CG) --- SERC COPYRIGHT
C      COMMENTED    24 JULY 1981 (CG)
C
C ARGUMENTS IN
C      XI      LOCAL COORDINATE OF POINT AT WHICH AREA IS
C              REQUIRED
C      ETA     LOCAL COORDINATE OF POINT AT WHICH AREA IS
C              REQUIRED
C      ZETA    LOCAL COORDINATE OF POINT AT WHICH AREA IS
C              REQUIRED
C      GEOM    LOCAL COORDINATE ARRAY CONTAINING THE GLOBAL
C              COORDIANTES OF EACH NODE ON AN ELEMENT IN THE
C              LOCAL ORDER
C      IGEOM   FIRST DIMENSION OF ARRAY GEOM (.GE.NODEL)
C      JGEOM   SECOND DIMENSION OF ARRAY GEOM (.GE.3)
C      NODEL   NUMBER OF NODES ON THE ELEMENT
C      FACNUM  THE FACE NUMBER OF THE FACE OF THE ELEMENT TO
C              BE USED IN CALCULATING THE AREA (.LE.6)
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      UAREA    THE UNIT OF AREA AT THE SPECIFIED POINT
C
C ROUTINES CALLED
C      ERRMES  BRK8  BRK20  BRK32  MATMUL
C
C
C     SUBROUTINE SURBRK(XI,ETA,GEOM,IGEOM,JGEOM,NODEL,FACNUM,UAREA,ITEST
C***********************************************************************
C
      INTEGER DIMEN, FACNUM, IDER, IFUN, IGEOM, IJAC, ITEST,
     *     JDER, JGEOM, JJAC, N, NODEL, I, J, IERROR, ERRMES
      DOUBLE PRECISION DER, ETA, FUN, G1, G2, G3, GEOM, 
     *     JAC, SRNAME, UAREA, XI, ZETA
      DIMENSION DER(3,32), FUN(32), GEOM(IGEOM,JGEOM),
     *     JAC(3,3)
      DATA DIMEN /3/, IDER /3/, IFUN /32/, IJAC /3/, JDER /32/,
     *     JJAC /3/, SRNAME /8H SURBRK /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF ((NODEL.LE.0) .OR. (FACNUM.LE.0))
     *                      IERROR = 1
                        IF (IGEOM.LT.NODEL .OR. JGEOM.LT.3)
     *                      IERROR = 2
                        IF (FACNUM.LT.1 .OR. FACNUM.GT.6)
     *                      IERROR = 3
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 N = NODEL/8
      IF(ITEST.EQ.-1) GOTO 1015
      IERROR = 0
      IF(N.LE.0.OR.N.GE.4) IERROR=4
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF(ITEST.NE.0) RETURN
 1015 GO TO (1020, 1030, 1040), N
 1020 CALL BRK8(FUN, IFUN, DER, IDER, JDER, XI, ETA, ZETA, ITEST)
      GO TO 1050
 1030 CALL BRK20(FUN, IFUN, DER, IDER, JDER, XI, ETA, ZETA, ITEST)
      GO TO 1050
 1040 CALL BRK32(FUN, IFUN, DER, IDER, JDER, XI, ETA, ZETA, ITEST)
 1050 CALL MATMUL(DER, IDER, JDER, GEOM, IGEOM, JGEOM, JAC, IJAC,
     *     JJAC, DIMEN, NODEL, DIMEN, ITEST)
      GO TO (1080, 1080, 1060, 1070, 1060, 1070), FACNUM
C++++
C    XI = CONSTANT
C
 1060 I=2
      J=3
      GO TO 1090
C++++
C    ETA = CONSTANT
C
 1070 I=1
      J=3
      GO TO 1090
C++++
C    ZETA = CONSTANT
C
 1080 I=1
      J=2
 1090 G1 = JAC(I,2)*JAC(J,3) - JAC(I,3)*JAC(J,2)
      G2 = JAC(I,3)*JAC(J,1) - JAC(I,1)*JAC(J,3)
      G3 = JAC(I,1)*JAC(J,2) - JAC(I,2)*JAC(J,1)
      UAREA = DSQRT(G1**2+G2**2+G3**2)
      RETURN
      END
