C
      SUBROUTINE FREDIF(IELE,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,
     *                  FIRST,DIF,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      FREDIF calculates the maximum freedom number difference for an
C      element
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
C      Commented    18 Feb 1980 (KR)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      IELE    element number
C      ELTOP   ELTOP(I, 1) = element type of element I
C              ELTOP(I, 2) = number of nodes on element I
C              ELTOP(I, J+2), J=1(1)NUMBER of nodes on element,
C              contains the nodes associated with element I
C      IELTOP  FIRST dimension of array ELTOP (.GE. IELE)
C      JELTOP  second dimension of ELTOP (.GE. number of nodes
C              on element)
C      NF      NF(I, J) contains the freedom numbers associated
C              with node I
C      INF     FIRST dimension of NF (.GE. maximum node number
C              on element)
C      JNF     second dimension of NF (.GE. DOFNOD)
C      DOFNOD  number of degrees of freedom per node on the
C              element
C      FIRST   must be set to .true. for the FIRST call to
C              FREDIF and .false. for subsequent calls
C      DIF     must be zero for FIRST call to FREDIF
C              subsequently contains the maximum freedom
C              difference prior to the current call
C      ITEST   error checking option
C
C ARGUMENTS out
C      FIRST   set to .FALSE.
C      DIF     maximum freedom difference for all elements up
C              to and including element nele
C
C ROUTINES called
C      ERRMES MAXINT
C
C     SUBROUTINE FREDIF(IELE,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,
C    *                  FIRST,DIF,ITEST)
C***********************************************************************
C
      INTEGER DIF,DOFNOD,ELTOP,ERRMES,I,IDEG,IELE,IELTOP,IERROR,INF,
     *        INOD,ITEST,J,JELTOP,JNF,JTEST,MAXMUM,MAXINT,MINMUM,NF,
     *        NODEL
      CHARACTER*6 SRNAME
      LOGICAL FIRST
      DIMENSION ELTOP(IELTOP,JELTOP),NF(INF,JNF)
C
      INTRINSIC MAX,MIN
      EXTERNAL ERRMES,MAXINT
C
      DATA SRNAME/'FREDIF'/
C
C     Initialisation
C
      IF (FIRST) THEN
         FIRST = .FALSE.
         DIF = 0
      END IF
C
C     Parameter checking
C
      JTEST = ITEST
      IF (JTEST.NE.-1) THEN
         IERROR = 0
         IF (IELTOP.LT.IELE) IERROR = 3
         IF (JNF.LT.DOFNOD) IERROR = 2
         IF (IELE.LE.0 .OR. DOFNOD.LE.0) IERROR = 1
         ITEST = ERRMES(JTEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
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
      MAXMUM = 0
      MINMUM = MAXINT(MAXMUM)
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
         DO 1000 J = 1,DOFNOD
            IDEG = NF(INOD,J)
            IF (IDEG.GT.0) THEN
               MAXMUM = MAX(IDEG,MAXMUM)
               MINMUM = MIN(IDEG,MINMUM)
            END IF
 1000    CONTINUE
 1010 CONTINUE
C
C     Calculate maximum freedom number difference
C
      DIF = MAX(DIF,MAXMUM-MINMUM)
C
      END
