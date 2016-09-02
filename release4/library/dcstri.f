C
      SUBROUTINE DCSTRI(JACIN,IJACIN,JJACIN,SIDNUM,COSIN,ICOSIN,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      DCSTRI calculates the direction cosines of the outward normal
C      to the specified side of a triangular element given the
C      jacobian inverse
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
C      Commented     1 Feb 1981 (CG)
C      Release 3.0   1 Jul 1996 (CG)
C      Release 4.0   2 Oct 2003 (CG)
C
C ARGUMENTS in
C      JACIN   array of dimension (IJACIN, JJACIN). contains the
C              inverse of the transformation jacobian.
C      IJACIN  first dimension of JACIN (.GE. 2)
C      JJACIN  second dimension of JACIN (.GE. 2)
C      SIDNUM  the side number for which the outward normal
C              is required (.LE. 3)
C      ICOSIN  dimension of vector COSIN (.GE. 2)
C      ITEST   error checking option
C
C ARGUMENTS out
C      COSIN   vector of dimension ICOSIN. contains the dir-
C              ection cosines of the outward normal
C
C ROUTINES called
C      ERRMES  MATVEC
C
C     SUBROUTINE DCSTRI(JACIN,IJACIN,JJACIN,SIDNUM,COSIN,ICOSIN,ITEST)
C***********************************************************************
C
      INTEGER DIMEN,ERRMES,ICOSIN,IERROR,IJACIN,ITEST,IWORK,JJACIN,
     *        SIDNUM
      CHARACTER*6 SRNAME
      DOUBLE PRECISION AMOD,COSIN,JACIN,WORK
      DIMENSION COSIN(ICOSIN),JACIN(IJACIN,JJACIN),WORK(2)
      DATA DIMEN/2/,IWORK/2/,SRNAME/'DCSTRI'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IJACIN.LT.2 .OR. JJACIN.LT.2) IERROR = 3
         IF (ICOSIN.LT.2) IERROR = 2
         IF (SIDNUM.LE.0 .OR. SIDNUM.GT.3) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      GO TO (1000,1010,1020) SIDNUM
C
C     Side number 1 : l3 = 0
C
 1000 CONTINUE
      WORK(1) = 0.5D0
      WORK(2) = -DSQRT(3.0D0)/2.0D0
      GO TO 1030
C
C     Side number 2 : l1 = 0
C
 1010 CONTINUE
      WORK(1) = -1.0D0
      WORK(2) = 0.0D0
      GO TO 1030
C
C     Side number 3 : l2 = 0
C
 1020 CONTINUE
      WORK(1) = 0.5D0
      WORK(2) = DSQRT(3.0D0)/2.0D0
C
C     Calculate direction cosines
C
 1030 CONTINUE
      CALL MATVEC(JACIN,IJACIN,JJACIN,WORK,IWORK,DIMEN,DIMEN,COSIN,
     *            ICOSIN,ITEST)
      AMOD = DSQRT(COSIN(1)*COSIN(1)+COSIN(2)*COSIN(2))
      COSIN(1) = COSIN(1)/AMOD
      COSIN(2) = COSIN(2)/AMOD
C
      END
