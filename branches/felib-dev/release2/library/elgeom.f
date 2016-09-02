C***********************************************************************
C$SPLIT$ELGEOM$*********************************************************
C***********************************************************************
      SUBROUTINE ELGEOM(NELE, ELTOP, IELTOP, JELTOP, COORD,
     *     ICOORD, JCOORD, GEOM, IGEOM, JGEOM, DIMEN, ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CONSTRUCTS THE ELEMENT GEOMETRY ARRAY FOR THE SPECIFIED
C      ELEMENT IN THE LOCAL NODE ORDER
C
C HISTORY
C      RELEASE 1.1  29 OCT 1979 (CG) --- SERC COPYRIGHT
C      COMMENTED    14 FEB 1980 (KR)
C
C ARGUMENTS IN
C      NELE    ELEMENT NUMBER FOR WHICH THE NODAL COORDINATE
C              ARRAY IS REQUIRED
C      ELTOP   CONTAINS ELEMENT TYPE, NUMBER OF NODES, AND
C              ELEMENT TOPOLOGIES FOR ALL THE ELEMENTS
C      IELTOP  FIRST DIMENSION OF ELTOP (.GE. MAXIMUM NODE
C              NUMBER ON REGION)
C      JELTOP  SECOND DIMENSION OF ELTOP (.GE. NUMBER OF NODES
C              ON ELEMENT + 2)
C      COORD   COORD(I,J) CONTAINS THE J'TH GLOBAL COORDINATE
C              OF THE I'TH NODE
C      ICOORD  FIRST DIMENSION OF COORD (.GE. MAXIMUM NUMBER OF
C              NODES IN THE MESH)
C      JCOORD  SECOND DIMENSION OF COORD (.GE. DIMEN)
C      IGEOM   FIRST DIMENSION OF ARRAY GEOM (.GE. NUMBER OF
C              NODES ON THE ELEMENT)
C      JGEOM   SECOND DIMENSION OF GEOM (.GE. DIMEN)
C      DIMEN   DIMENSIONALITY OF THE MESH
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      GEOM    CONTAINS GLOBAL COORDINATES OF THE NODES
C              ASSOCIATED WITH ELEMENT NELE IN THE LOCAL NODE
C              NUMBERING ORDER
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE ELGEOM(NELE, ELTOP, IELTOP, JELTOP, COORD, ICOORD,
C    *     JCOORD, GEOM, IGEOM, JGEOM, DIMEN, ITEST)
C***********************************************************************
C
      INTEGER DIMEN, ELTOP, ERRMES, ICOORD, IDIM, IELTOP, IERROR,
     *     IGEOM, INOD, ITEST, JCOORD, JELTOP, JGEOM, JNOD,
     *     NELE, NODEL
      DOUBLE PRECISION COORD, GEOM, SRNAME
      DIMENSION COORD(ICOORD,JCOORD), ELTOP(IELTOP,JELTOP),
     *     GEOM(IGEOM,JGEOM)
      DATA SRNAME /8H ELGEOM /
                        IF (ITEST.EQ.-1) GO TO 1010
                        IERROR = 0
                        IF (JGEOM.LT.DIMEN) IERROR = 4
                        IF (JCOORD.LT.DIMEN) IERROR = 3
                        IF (IELTOP.LT.NELE) IERROR = 2
                        IF (NELE.LE.0 .OR. DIMEN.LE.0) IERROR = 1
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1010                   NODEL = ELTOP(NELE,2)
                        IF (ITEST.EQ.-1) GO TO 1020
                        IERROR = 0
                        IF (JELTOP.LT.NODEL+2) IERROR = 5
                        IF (IGEOM.LT.NODEL) IERROR = 7
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1020                   DO 1050 JNOD=1,NODEL
                        INOD = ELTOP(NELE,JNOD+2)
                        IF (ITEST.EQ.-1) GO TO 1030
                        IERROR = 0
                        IF (ICOORD.LT.INOD) IERROR = 6
                        ITEST = ERRMES(ITEST,IERROR,SRNAME)
                        IF (ITEST.NE.0) RETURN
 1030                   DO 1040 IDIM=1,DIMEN
                        GEOM(JNOD,IDIM) = COORD(INOD,IDIM)
 1040                   CONTINUE
 1050                   CONTINUE
                        RETURN
                        END
