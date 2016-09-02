C
      SUBROUTINE DIRECT(NELE,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,
     *                  STEER,ISTEER,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      DIRECT constructs the steering vector to direct the assembly of a
C      system matrix
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
C      Commented    13 Feb 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      NELE    element number
C      ELTOP   2D array containing element type, number of
C              nodes on the element, and the element topologies
C      IELTOP  first dimension of ELTOP (.GE. number of
C              elements in problem)
C      JELTOP  second dimension of ELTOP (.GE. number of nodes
C              on element + 2)
C      NF      contains freedom numbers associated with each
C              node
C      INF     first dimension of NF (.GE. total number of
C              nodes in problem)
C      JNF     second dimension of NF (.GE. DOFNOD)
C      DOFNOD  number of degrees of freedom at each node
C      ISTEER  dimension of vector STEER (.GE. total number of
C              degrees of freedom on element)
C      ITEST   error checking option
C
C ARGUMENTS out
C      STEER   vector containing freedom numbers associated
C              with element NELE, arranged in local node order
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE DIRECT(NELE,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,
C    *                  STEER,ISTEER,ITEST)
C***********************************************************************
C
      INTEGER DOFNOD,ELTOP,ERRMES,IDEG,IELTOP,IERROR,INF,INOD,ISTEER,
     *        ITEST,JELTOP,JNF,JNOD,JTEST,K,NELE,NF,NODEL,STEER
      CHARACTER*6 SRNAME
      DIMENSION ELTOP(IELTOP,JELTOP),NF(INF,JNF),STEER(ISTEER)
      DATA SRNAME/'DIRECT'/
C
C     Parameter checking
C
      JTEST = ITEST
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (JNF.LT.DOFNOD) IERROR = 3
         IF (IELTOP.LT.NELE) IERROR = 2
         IF (NELE.LE.0 .OR. DOFNOD.LE.0) IERROR = 1
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main loops
C
      NODEL = ELTOP(NELE,2)
C
C     Range check on NODEL : number of node in element
C
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (JELTOP.LT.NODEL+2) IERROR = 4
         IF (ISTEER.LT.DOFNOD*NODEL) IERROR = 6
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
      K = 1
      DO 1010 INOD = 1,NODEL
         JNOD = ELTOP(NELE,INOD+2)
C
C     Range check on JNOD : node number
C
         IF (JTEST.NE.-1) THEN
            IERROR = 0
            IF (INF.LT.JNOD) IERROR = 5
            ITEST = ERRMES(JTEST,IERROR,SRNAME)
            IF (ITEST.NE.0) RETURN
         END IF
C
C     Construct steering vector
C
         DO 1000 IDEG = 1,DOFNOD
            STEER(K) = NF(JNOD,IDEG)
            K = K + 1
 1000    CONTINUE
 1010 CONTINUE
C
      END
