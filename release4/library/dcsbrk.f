C
      SUBROUTINE DCSBRK(JACIN,IJACIN,JJACIN,FACNUM,COSIN,ICOSIN,ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      DCSBRK calculates the direction cosines of the outward normal
C      to the specified face of a hexahedral element given the
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
C      IJACIN  first dimension of JACIN (.GE. 3)
C      JJACIN  second dimension of JACIN (.GE. 3)
C      FACNUM  the side number for which the outward normal
C              is required (.LE. 6)
C      ICOSIN  dimension of vector COSIN (.GE. 3)
C      ITEST   error checking option
C
C ARGUMENTS out
C      COSIN   vector of dimension ICOSIN. contains the dir-
C              ection cosines of the outward normal
C
C ROUTINES called
C      ERRMES  MATVEC
C
C     SUBROUTINE DCSBRK(JACIN,IJACIN,JJACIN,FACNUM,COSIN,ICOSIN,ITEST)
C***********************************************************************
C
      INTEGER DIMEN,ERRMES,FACNUM,ICOSIN,IERROR,IJACIN,ITEST,IWORK,
     *        JJACIN
      CHARACTER*6 SRNAME
      DOUBLE PRECISION AMOD,COSIN,JACIN,WORK
      DIMENSION COSIN(ICOSIN),JACIN(IJACIN,JJACIN),WORK(3)
      DATA DIMEN/3/,IWORK/3/,SRNAME/'DCSBRK'/
C
C     Parameter checking
C
      IF (ITEST.NE.-1) THEN
         IERROR = 0
         IF (IJACIN.LT.3 .OR. JJACIN.LT.3) IERROR = 3
         IF (ICOSIN.LT.3) IERROR = 2
         IF (FACNUM.LE.0 .OR. FACNUM.GT.6) IERROR = 1
         ITEST = ERRMES(ITEST,IERROR,SRNAME)
         IF (ITEST.NE.0) RETURN
      END IF
C
C     Main body
C
      GO TO (1000,1010,1020,1030,1040,1050) FACNUM
      GO TO 1020
C
C     Face number 1 : zeta = -1
C
 1000 CONTINUE
      WORK(1) = 0.0D0
      WORK(2) = 0.0D0
      WORK(3) = -1.0D0
      GO TO 1060
C
C     Face number 2 : zeta = +1
C
 1010 CONTINUE
      WORK(1) = 0.0D0
      WORK(2) = 0.0D0
      WORK(3) = 1.0D0
      GO TO 1060
C
C     Face number 3 : xi = -1
C
 1020 CONTINUE
      WORK(1) = -1.0D0
      WORK(2) = 0.0D0
      WORK(3) = 0.0D0
      GO TO 1060
C
C     Face number 4 : eta = +1
C
 1030 CONTINUE
      WORK(1) = 0.0D0
      WORK(2) = 1.0D0
      WORK(3) = 0.0D0
      GO TO 1060
C
C     Face number 5 : xi = +1
C
 1040 CONTINUE
      WORK(1) = 1.0D0
      WORK(2) = 0.0D0
      WORK(3) = 0.0D0
      GO TO 1060
C
C     Face number 6 : eta = -1
C
 1050 CONTINUE
      WORK(1) = 0.0D0
      WORK(2) = -1.0D0
      WORK(3) = 0.0D0
C
C     Calculate direction cosines
C
 1060 CONTINUE
      CALL MATVEC(JACIN,IJACIN,JJACIN,WORK,IWORK,DIMEN,DIMEN,COSIN,
     *            ICOSIN,ITEST)
      AMOD = DSQRT(COSIN(1)*COSIN(1)+COSIN(2)*COSIN(2)+
     *       COSIN(3)*COSIN(3))
      COSIN(1) = COSIN(1)/AMOD
      COSIN(2) = COSIN(2)/AMOD
      COSIN(3) = COSIN(3)/AMOD
C
      END
