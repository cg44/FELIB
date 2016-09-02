C
      SUBROUTINE BNDWTH(ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,TOTELS,
     *                  HBAND,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      BNDWTH calculates the maximum freedom number difference over
C      all the elements in the mesh.
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
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      ELTOP   ELTOP(I,1) = element type of element I
C              ELTOP(I,2) = number of nodes on element I
C              ELTOP(I,J+2), J=1(1)NUMBER of nodes on element,
C              contains the nodes associated with element I
C      IELTOP  first dimension of array ELTOP (.GE. TOTELS)
C      JELTOP  second dimension of ELTOP (.GE. number of nodes
C              on element)
C      NF      NF(I,J) contains the freedom numbers associated
C              with node I
C      INF     first dimension of NF (.GE. maximum node number
C              on element)
C      JNF     second dimension of NF (.GE. DOFNOD)
C      DOFNOD  number of degrees of freedom per node on the
C              element
C      TOTELS  number of elements in the mesh
C      ITEST   error checking option
C
C ARGUMENTS out
C      HBAND   semi bandwidth
C
C ROUTINES called
C      ERRMES    MAXINT
C
C     SUBROUTINE BNDWTH(ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,TOTELS,
C    *                  HBAND,ITEST)
C*********************************************************
C
      INTEGER DOFNOD,ELTOP,ERRMES,HBAND,I,IDEG,IELE,IELTOP,IERROR,INF,
     *        INOD,ITEST,J,JELTOP,JNF,JTEST,MAX,MAXINT,MIN,NF,NODEL,
     *        TOTELS
      CHARACTER*6 SRNAME
      DIMENSION ELTOP(IELTOP,JELTOP),NF(INF,JNF)
C
      EXTERNAL ERRMES,MAXINT
C
      DATA SRNAME/'BNDWTH'/
C
C     Initialisation
C
      HBAND = 0
      JTEST = ITEST
C
C     Parameter checking
C
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (DOFNOD.LE.0) IERROR = 1
         IF (JNF.LT.DOFNOD) IERROR = 2
         IF (TOTELS.GT.IELTOP) IERROR = 3
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main scanning loops
C
      DO 1020 IELE = 1,TOTELS
         MAX = 0
         MIN = MAXINT(MAX)
         NODEL = ELTOP(IELE,2)
C
C     Range checking on NODEL
C
         IF (JTEST.NE.-1) THEN
            IERROR = 0
            IF (JELTOP.LT.NODEL+2) IERROR = 4
            ITEST = ERRMES(JTEST,IERROR,SRNAME)
            IF (ITEST.NE.0) RETURN
         END IF
C
C     Loop over nodes
C
         DO 1010 I = 1,NODEL
            INOD = ELTOP(IELE,I+2)
C
C     Range checking on INOD
C
            IF (JTEST.NE.-1) THEN
               IERROR = 0
               IF (INF.LT.INOD) IERROR = 5
               ITEST = ERRMES(JTEST,IERROR,SRNAME)
               IF (ITEST.NE.0) RETURN
            END IF
C
C     Loop over freedoms
C
            DO 1000 J = 1,DOFNOD
               IDEG = NF(INOD,J)
               IF (IDEG.NE.0) THEN
                  MAX = MAX0(IDEG,MAX)
                  MIN = MIN0(IDEG,MIN)
               END IF
 1000       CONTINUE
 1010    CONTINUE
C
C     Maximum freedom number difference
C
         HBAND = MAX0(HBAND,MAX-MIN)
 1020 CONTINUE
C
C     Semi band width
C
      HBAND = HBAND + 1
C
      END
