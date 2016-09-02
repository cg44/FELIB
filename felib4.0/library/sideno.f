C
      SUBROUTINE SIDENO(TOTELS,ELTOP,IELTOP,JELTOP,M,BDCND,IBDCND,
     *                  JBDCND,NUMSID,BLIST,IBLIST,JBLIST,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      SIDENO the routine generates a list of element and side
C      numbers from a list of boundary nodes. the node list
C      must be a continuous sequence of node in the same
C      sense as the local node ordering.
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
C      Release 2.0     Jan 1982 (CG)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      TOTELS  the total number of elements in the mesh
C      ELTOP   the element topology array containing element
C              type number, number of nodes on element and
C              topology list
C      IELTOP  first dimension of ELTOP
C      JELTOP  second dimension of ELTOP
C      M       boundary condition list number
C      BDCND   boundary condition array containing boundary
C              condition type, number of nodes on element side
C              number of nodes in list and node list.
C      IBDCND  first dimension of BDCND array
C      JBDCND  second dimension of BDCND array
C      IBLIST  first dimension BLIST
C      JBLIST  second dimension of BLIST
C      ITEST   error checking option
C
C ARGUMENTS out
C      NUMSID  number of element sides generated from the
C              current boundary node list
C      BLIST   array containing element number and side number
C
C ROUTINES called
C      ERRMES
C
C     SUBROUTINE SIDENO(TOTELS,ELTOP,IELTOP,JELTOP,M,BDCND,IBDCND,
C    *                  JBDCND,NUMSID,BLIST,IBLIST,JBLIST,ITEST)
C***********************************************************************
C
      INTEGER BDCND,BLIST,BNUM1,ELTOP,ERRMES,I,IBDCND,IBLIST,IELTOP,
     *        IERROR,IK,ITEST,J,JBDCND,JBLIST,JELTOP,JTEST,K,L,M,N,N1,
     *        NB,NELE,NODEL,NSIDE1,NUMSID,TOTELS
      CHARACTER*6 SRNAME
      DIMENSION BDCND(IBDCND,JBDCND),BLIST(IBLIST,JBLIST),
     *          ELTOP(IELTOP,JELTOP)
      DATA SRNAME/'SIDENO'/
C
C     Parameter checking
C
      JTEST = ITEST
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (JBLIST.LT.2) IERROR = 4
         IF (IBDCND.LT.M) IERROR = 3
         IF (IELTOP.LT.TOTELS) IERROR = 2
         IF (TOTELS.LE.0 .OR. M.LE.0) IERROR = 1
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      NSIDE1 = BDCND(M,3) - 1
      BNUM1 = BDCND(M,2) - NSIDE1
      N = 0
      DO 1050 I = 1,BNUM1,NSIDE1
C
C     Range checking on I
C
         IF (JTEST.NE.-1) THEN
            IERROR = 0
            IF (I.GT.JBDCND) IERROR = 6
            ITEST = ERRMES(JTEST,IERROR,SRNAME)
            IF (ITEST.NE.0) RETURN
         END IF
C
         NB = BDCND(M,I+3)
         DO 1030 J = 1,TOTELS
            NODEL = ELTOP(J,2)
            DO 1000 K = 1,NODEL
C
C     Range checking on K
C
               IF (JTEST.NE.-1) THEN
                  IERROR = 0
                  IF (K+2.GT.JELTOP) IERROR = 5
                  ITEST = ERRMES(JTEST,IERROR,SRNAME)
                  IF (ITEST.NE.0) RETURN
               END IF
C
               IF (NB.EQ.ELTOP(J,K+2)) GO TO 1010
 1000       CONTINUE
            GO TO 1030
C
 1010       CONTINUE
            NELE = J
            IF (NODEL.EQ.10) NODEL = NODEL - 1
            N1 = K
            DO 1020 K = 1,NSIDE1
               L = MOD(N1+K-1,NODEL) + 3
               IK = I + K + 3
C
C     Range checking on L
C
               IF (JTEST.NE.-1) THEN
                  IERROR = 0
                  IF (L.GT.JELTOP) IERROR = 5
                  ITEST = ERRMES(JTEST,IERROR,SRNAME)
                  IF (ITEST.NE.0) RETURN
               END IF
C
               IF (BDCND(M,IK).NE.ELTOP(J,L)) GO TO 1030
 1020       CONTINUE
            GO TO 1040
 1030    CONTINUE
         GO TO 1050
C
 1040    CONTINUE
         N = N + 1
C
C     Range checking on N
C
         IF (JTEST.NE.-1) THEN
            IERROR = 0
            IF (N.GT.IBLIST) IERROR = 7
            ITEST = ERRMES(JTEST,IERROR,SRNAME)
            IF (ITEST.NE.0) RETURN
         END IF
C
         BLIST(N,1) = NELE
         BLIST(N,2) = (N1-1)/NSIDE1 + 1
 1050 CONTINUE
C
      NUMSID = N
C
      END
