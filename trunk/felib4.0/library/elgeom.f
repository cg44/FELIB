C
      SUBROUTINE ELGEOM(NELE,ELTOP,IELTOP,JELTOP,COORD,ICOORD,JCOORD,
     *                  GEOM,IGEOM,JGEOM,DIMEN,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      ELGEOM constructs the element geometry array for the specified
C      element in the local node order
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
C      Commented    14 Feb 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      NELE    element number for which the nodal coordinate
C              array is required
C      ELTOP   contains element type, number of nodes, and
C              element topologies for all the elements
C      IELTOP  first dimension of ELTOP (.GE. maximum node
C              number on region)
C      JELTOP  second dimension of ELTOP (.GE. number of nodes
C              on element + 2)
C      COORD   COORD(I, J) contains the j'th global coordinate
C              of the i'th node
C      ICOORD  first dimension of COORD (.GE. maximum number of
C              nodes in the mesh)
C      JCOORD  second dimension of COORD (.GE. DIMEN)
C      IGEOM   first dimension of array GEOM (.GE. number of
C              nodes on the element)
C      JGEOM   second dimension of GEOM (.GE. DIMEN)
C      DIMEN   dimensionality of the mesh
C      ITEST   error checking option
C
C ARGUMENTS out
C      GEOM    contains global coordinates of the nodes
C              associated with element NELE in the local node
C              numbering order
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE ELGEOM(NELE,ELTOP,IELTOP,JELTOP,COORD,ICOORD,JCOORD,
C    *                  GEOM,IGEOM,JGEOM,DIMEN,ITEST)
C***********************************************************************
C
      INTEGER DIMEN,ELTOP,ERRMES,ICOORD,IDIM,IELTOP,IERROR,IGEOM,INOD,
     *        ITEST,JCOORD,JELTOP,JGEOM,JNOD,JTEST,NELE,NODEL
      CHARACTER*6 SRNAME
      DOUBLE PRECISION COORD,GEOM
      DIMENSION COORD(ICOORD,JCOORD),ELTOP(IELTOP,JELTOP),
     *          GEOM(IGEOM,JGEOM)
      DATA SRNAME/'ELGEOM'/
C
C     Parameter checking
C
      JTEST = ITEST
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (JGEOM.LT.DIMEN) IERROR = 4
         IF (JCOORD.LT.DIMEN) IERROR = 3
         IF (IELTOP.LT.NELE) IERROR = 2
         IF (NELE.LE.0 .OR. DIMEN.LE.0) IERROR = 1
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      NODEL = ELTOP(NELE,2)
C
C     Range checking on NODEL : NODEL+2 <= JELTOP and NODEL <= IGEOM
C
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (JELTOP.LT.NODEL+2) IERROR = 5
         IF (IGEOM.LT.NODEL) IERROR = 7
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
      DO 1010 JNOD = 1,NODEL
         INOD = ELTOP(NELE,JNOD+2)
C
C     Range checking on INOD : INOD <= ICOORD
C
         IF (JTEST.NE.-1) THEN
            IERROR = 0
            IF (ICOORD.LT.INOD) IERROR = 6
            ITEST = ERRMES(JTEST,IERROR,SRNAME)
            IF (ITEST.NE.0) RETURN
         END IF
C
C     Construct geometry array GEOM
C
         DO 1000 IDIM = 1,DIMEN
            GEOM(JNOD,IDIM) = COORD(INOD,IDIM)
 1000    CONTINUE
 1010 CONTINUE
C
      END
