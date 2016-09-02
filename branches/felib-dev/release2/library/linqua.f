C***********************************************************************
C$SPLIT$LINQUA$*********************************************************
C***********************************************************************
      SUBROUTINE LINQUA(XI, ETA, GEOM, IGEOM, JGEOM, NODEL,
     *     SIDNUM, ALEN, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CALCULATES A UNIT OF LENGTH ALONG THE SIDE OF A
C      RECTANGULAR ELEMENT (4, 8 OR 12 NODED)
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
C              BE USED IN CALCULATING THE LENGTH (.LE.3)
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      ALEN    THE UNIT OF LENGTH AT THE SPECIFIED POINT
C
C ROUTINES CALLED
C      ERRMES  QUAM4  QUAM8  QUAM12  MATMUL
C
C
C     SUBROUTINE LINQUA(XI,ETA,GEOM,IGEOM,JGEOM,NODEL,SIDNUM,ALEN,ITEST)
C***********************************************************************
C
      INTEGER DIMEN, ERRMES, IDER, IERROR, IFUN, IGEOM, IJAC, ITEST,
     *     JDER, JGEOM, JJAC, N, NODEL, SIDNUM                            
      DOUBLE PRECISION ALEN, DER, ETA, FUN, GEOM, JAC,
     *     SRNAME, XI
      DIMENSION DER(2,12), FUN(12), GEOM(IGEOM,JGEOM),
     *     JAC(2,2)
      DATA DIMEN /2/, IDER /2/, IFUN /12/, IJAC /2/, JDER /12/,
     *     JJAC /2/, SRNAME /8H LINQUA /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF ((NODEL.LE.0) .OR. (SIDNUM.LE.0))
     *                      IERROR = 1
                        IF (IGEOM.LT.NODEL .OR. JGEOM.LT.2)
     *                      IERROR = 2
                        IF (SIDNUM.LT.1 .OR. SIDNUM.GT.4)
     *                      IERROR = 3
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010 N = NODEL/4
      IF (ITEST.EQ.-1) GO TO 1015
      IERROR = 0
      IF(N.LE.0.OR.N.GE.4) IERROR=4
      ITEST = ERRMES(ITEST,IERROR,SRNAME)
      IF(ITEST.NE.0) RETURN
 1015 GO TO (1020, 1030, 1040), N
 1020 CALL QUAM4(FUN, IFUN, DER, IDER, JDER, XI, ETA, ITEST)
      GO TO 1050
 1030 CALL QUAM8(FUN, IFUN, DER, IDER, JDER, XI, ETA, ITEST)
      GO TO 1050
 1040 CALL QUAM12(FUN, IFUN, DER, IDER, JDER, XI, ETA, ITEST)
 1050 CALL MATMUL(DER, IDER, JDER, GEOM, IGEOM, JGEOM, JAC, IJAC,
     *     JJAC, DIMEN, NODEL, DIMEN, ITEST)
      GO TO (1060, 1070, 1060, 1070), SIDNUM
 1060 ALEN = DSQRT(JAC(2,1)**2+JAC(2,2)**2)
      RETURN
 1070 ALEN = DSQRT(JAC(1,1)**2+JAC(1,2)**2)
      RETURN
      END
